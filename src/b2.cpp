#include "apex/apex.hpp"

#include <deque>
#include <fstream>
#include <sstream>

namespace apex {
namespace {

std::vector<EmbeddingImage> dedupe(std::vector<EmbeddingImage> images) {
    std::unordered_set<std::string> seen;
    std::vector<EmbeddingImage> result;
    for (auto& image : images) {
        const std::string key = canonical_image_key(image);
        if (seen.insert(key).second) result.push_back(std::move(image));
    }
    return result;
}

EmbeddingImage compact_with_source_map(Embedding z, const Morphism& source_to_mutated) {
    CompactResult c = compact(z);
    return {std::move(c.graph), compose(c.map, source_to_mutated)};
}

}  // namespace

std::vector<EmbeddingImage> boundary_completions(const Embedding& z, VertexId v) {
    if (!z.vertex_alive[v] || !z.is_boundary(v) || !z.degree_range[v].fixed() ||
        z.degree_range[v].lower != z.degree(v)) {
        throw std::invalid_argument("boundary_completions: v is not a fixed full boundary vertex");
    }
    if (z.kind == EmbeddingKind::PseudoEmbedding) {
        Embedding target = z;
        const DartId first = target.first_dart(v);
        const DartId last = target.last_dart(v);
        target.darts[first].pred = last;
        target.darts[last].succ = first;
        return {{std::move(target), identity_morphism(z)}};
    }

    std::vector<EmbeddingImage> result;
    if (auto a = add_boundary_darts(z, v); a.has_value()) result.push_back(std::move(*a));
    if (auto b = identify_neighbors(z, v); b.has_value()) result.push_back(std::move(*b));
    return dedupe(std::move(result));
}

std::optional<EmbeddingImage> add_boundary_darts(const Embedding& z, VertexId v) {
    const DartId first = z.first_dart(v);
    const DartId last = z.last_dart(v);
    if (first == kNil || last == kNil) return std::nullopt;
    if (z.tail(first) == z.tail(last)) return std::nullopt;
    Embedding target = z;
    target.darts[first].pred = last;
    target.darts[last].succ = first;
    add_boundary_darts_directly(target, first, last);
    std::string error;
    if (!target.validate_faces(&error)) return std::nullopt;
    return EmbeddingImage{std::move(target), identity_morphism(z)};
}

void add_boundary_darts_directly(Embedding& z, DartId first, DartId last) {
    const VertexId u = z.tail(first);
    const VertexId w = z.tail(last);
    if (u == kNil || w == kNil) throw std::invalid_argument("invalid boundary darts");

    const DartId f = static_cast<DartId>(z.darts.size());
    const DartId g = f + 1;
    z.darts.push_back({u, g, kNil, z.darts[first].rev, true});
    z.darts.push_back({w, f, z.darts[last].rev, kNil, true});
    z.darts[z.darts[first].rev].succ = f;
    z.darts[z.darts[last].rev].pred = g;
}

std::optional<EmbeddingImage> identify_neighbors(const Embedding& z, VertexId v) {
    const DartId first = z.first_dart(v);
    const DartId last = z.last_dart(v);
    if (first == kNil || last == kNil) return std::nullopt;
    const VertexId u = z.tail(first);
    const VertexId w = z.tail(last);
    if (!z.degree_range[u].intersects(z.degree_range[w])) return std::nullopt;

    Embedding target = z;
    target.darts[first].pred = last;
    target.darts[last].succ = first;
    Morphism map = link_incidence_list_ends(target, target.darts[last].rev,
                                            target.darts[first].rev);
    if (target.has_loop()) return std::nullopt;
    auto compacted = compact_with_source_map(std::move(target), map);
    std::string error;
    if (!compacted.target.validate_faces(&error)) return std::nullopt;
    return compacted;
}

Morphism link_incidence_list_ends(Embedding& z, DartId u_first, DartId w_last) {
    const VertexId u = z.darts[u_first].head;
    const VertexId w = z.darts[w_last].head;
    Morphism phi = identity_morphism(z);
    z.darts[u_first].pred = w_last;
    z.darts[w_last].succ = u_first;
    if (u == w) return phi;

    const DegreeRange intersection = z.degree_range[u].intersection(z.degree_range[w]);
    if (intersection.empty()) throw std::runtime_error("linkIncidenceListEnds degree mismatch");
    for (DartId e = u_first; e != kNil; e = z.darts[e].succ) z.darts[e].head = w;
    z.degree_range[w] = intersection;
    phi.vertex[u] = w;
    z.vertex_alive[u] = false;
    return phi;
}

Embedding from_vertex_rotations(int n, const std::vector<std::vector<int>>& rotations,
                                const std::vector<std::pair<int, int>>& digons,
                                EmbeddingKind kind) {
    if (static_cast<int>(rotations.size()) != n) {
        throw std::invalid_argument("from_vertex_rotations: wrong rotation count");
    }
    std::vector<std::vector<bool>> is_digon(n, std::vector<bool>(n, false));
    for (auto [a, b] : digons) {
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) {
            throw std::invalid_argument("invalid digon endpoint");
        }
        is_digon[a][b] = is_digon[b][a] = true;
    }
    std::vector<std::vector<DartId>> ordinary(n, std::vector<DartId>(n, kNil));
    std::vector<std::vector<DartId>> extra(n, std::vector<DartId>(n, kNil));
    Embedding z;
    z.kind = kind;
    for (int i = 0; i < n; ++i) z.add_vertex({1, kInfinity});

    for (int a = 0; a < n; ++a) {
        std::unordered_set<int> seen;
        for (int b : rotations[a]) {
            if (b == -1) continue;
            if (b < 0 || b >= n) throw std::invalid_argument("rotation neighbor out of range");
            if (!seen.insert(b).second) {
                throw std::invalid_argument("duplicate neighbor in a rotation; encode it as a digon");
            }
            ordinary[a][b] = z.add_dart({});
            if (is_digon[a][b]) extra[a][b] = z.add_dart({});
        }
    }

    for (int a = 0; a < n; ++a) {
        const int size = static_cast<int>(rotations[a].size());
        if (size == 0) throw std::invalid_argument("empty rotation");
        for (int i = 0; i < size; ++i) {
            const int b = rotations[a][i];
            if (b == -1) continue;
            if (ordinary[b][a] == kNil) {
                throw std::invalid_argument("rotations of two vertices disagree");
            }
            const int s = rotations[a][(i + 1) % size];
            const int p = rotations[a][(i + size - 1) % size];
            const DartId succ = s == -1 ? kNil : ordinary[a][s];
            const DartId pred = p == -1 ? kNil : (is_digon[a][p] ? extra[a][p] : ordinary[a][p]);
            const DartId e = ordinary[a][b];
            if (is_digon[a][b]) {
                const DartId e2 = extra[a][b];
                z.darts[e] = {a, extra[b][a], e2, pred, true};
                z.darts[e2] = {a, ordinary[b][a], succ, e, true};
            } else {
                z.darts[e] = {a, ordinary[b][a], succ, pred, true};
            }
        }
    }
    std::string error;
    if (!z.validate_single_list(&error)) {
        throw std::runtime_error("invalid dart representation from rotations: " + error);
    }
    return z;
}

std::vector<std::tuple<VertexId, VertexId, VertexId>> find_cut_tuples(
    int n, int ring_size, const std::vector<std::vector<int>>& rotations) {
    std::vector<std::tuple<VertexId, VertexId, VertexId>> result;
    for (int i = ring_size; i < n; ++i) {
        std::set<int> ring_neighbors;
        int transitions = 0;
        const int d = static_cast<int>(rotations[i].size());
        for (int j = 0; j < d; ++j) {
            const int k1 = rotations[i][j];
            const int k2 = rotations[i][(j + 1) % d];
            if (0 <= k1 && k1 < ring_size) ring_neighbors.insert(k1);
            if (0 <= k1 && k1 < ring_size && k2 >= ring_size) ++transitions;
        }
        if (transitions >= 2 && ring_neighbors.size() != 2) {
            throw std::runtime_error("configuration is not normal at vertex " + std::to_string(i));
        }
        if (transitions == 2 && ring_neighbors.size() == 2) {
            auto it = ring_neighbors.begin();
            const int a = *it++;
            const int b = *it;
            result.emplace_back(i, a, b);
        }
    }
    return result;
}

std::vector<RootedConfiguration> extend_from_cut_vertices(
    const ConfigurationFile& configuration) {
    const int n = configuration.vertex_count;
    const int r = configuration.ring_size;
    auto cuts = find_cut_tuples(n, r, configuration.rotations);
    if (cuts.size() >= sizeof(std::size_t) * 8) throw std::runtime_error("too many cut vertices");
    std::vector<RootedConfiguration> result;

    const std::size_t choices = std::size_t{1} << cuts.size();
    for (std::size_t mask = 0; mask < choices; ++mask) {
        std::vector<int> adjacent_cutvertex(r, -1);
        for (std::size_t i = 0; i < cuts.size(); ++i) {
            auto [v, a, b] = cuts[i];
            adjacent_cutvertex[((mask >> i) & 1U) ? a : b] = v;
        }

        std::vector<int> old_to_new(n, -1);
        int next = 0;
        for (int i = 0; i < n; ++i) {
            if (i < r && adjacent_cutvertex[i] == -1) continue;
            old_to_new[i] = next++;
        }
        std::vector<std::vector<int>> rotations(next);
        for (int i = 0; i < n; ++i) {
            if (old_to_new[i] == -1) continue;
            for (int j : configuration.rotations[i]) {
                if (j == -1 || (i < r && j != adjacent_cutvertex[i]) ||
                    (j >= 0 && j < r && i != adjacent_cutvertex[j])) {
                    rotations[old_to_new[i]].push_back(-1);
                } else {
                    rotations[old_to_new[i]].push_back(old_to_new[j]);
                }
            }
        }
        std::vector<std::pair<int, int>> digons;
        for (auto [a, b] : configuration.digons) {
            if (old_to_new[a] != -1 && old_to_new[b] != -1) {
                digons.emplace_back(old_to_new[a], old_to_new[b]);
            }
        }
        Embedding z = from_vertex_rotations(next, rotations, digons,
                                            EmbeddingKind::PseudoEmbedding);
        for (int i = 0; i < r; ++i) {
            if (old_to_new[i] == -1) continue;
            const VertexId v = old_to_new[i];
            const int d = z.degree(v);
            z.degree_range[v] = {d + 1, kInfinity};
        }
        for (int i = r; i < n; ++i) {
            z.degree_range[old_to_new[i]] = {configuration.prescribed_degree[i],
                                             configuration.prescribed_degree[i]};
        }
        DartId special = maximum_degree_dart(z);
        if (special == kNil) throw std::runtime_error("no fixed-degree dart in configuration");
        RootedConfiguration original{configuration.name, z, special};
        result.push_back(original);
        Embedding mirrored_graph = mirror(z);
        result.push_back({configuration.name + ":mirror", std::move(mirrored_graph), special});
    }

    std::unordered_set<std::string> seen;
    std::vector<RootedConfiguration> unique;
    for (auto& c : result) {
        const std::string key = canonical_key(c.graph, c.special);
        if (seen.insert(key).second) unique.push_back(std::move(c));
    }
    return unique;
}

std::optional<std::pair<DartId, DartId>> two_digons_incident_with_same_vertex(
    const Embedding& z) {
    for (VertexId v : z.vertices()) {
        DartId first_digon = kNil;
        for (DartId e0 : z.darts_at(v)) {
            if (z.darts[e0].succ == kNil) continue;
            const DartId e1 = z.darts[z.darts[e0].succ].rev;
            if (z.darts[e1].succ == kNil) continue;
            const DartId e2 = z.darts[z.darts[e1].succ].rev;
            if (e2 == e0) {
                if (first_digon == kNil) first_digon = e0;
                else if (z.tail(first_digon) != z.tail(e0) || first_digon != e0) {
                    return std::make_pair(first_digon, e0);
                }
            }
        }
    }
    return std::nullopt;
}

std::vector<EmbeddingImage> enforce_single_digon_incidence(const Embedding& z) {
    std::queue<EmbeddingImage> queue;
    queue.push({z, identity_morphism(z)});
    std::vector<EmbeddingImage> result;
    std::unordered_set<std::string> visited;

    while (!queue.empty()) {
        EmbeddingImage current = std::move(queue.front());
        queue.pop();
        const std::string state_key = canonical_image_key(current);
        if (!visited.insert(state_key).second) continue;
        auto pair = two_digons_incident_with_same_vertex(current.target);
        if (!pair.has_value()) {
            result.push_back(std::move(current));
            continue;
        }
        for (auto image : free_homomorphism(current.target, {*pair})) {
            image.map = compose(image.map, current.map);
            queue.push(std::move(image));
        }
    }
    return dedupe(std::move(result));
}

std::vector<EmbeddingImage> free_homomorphism_and_enforce_single_digon_incidence(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests) {
    std::vector<EmbeddingImage> result;
    for (auto first : free_homomorphism(z, requests)) {
        for (auto second : enforce_single_digon_incidence(first.target)) {
            second.map = compose(second.map, first.map);
            result.push_back(std::move(second));
        }
    }
    return dedupe(std::move(result));
}

}  // namespace apex

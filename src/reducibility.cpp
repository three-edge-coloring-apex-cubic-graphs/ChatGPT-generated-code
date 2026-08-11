#include "apex/apex.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <sstream>

namespace apex {
namespace {

using Matching = std::set<std::pair<int, int>>;
using ComponentMask = std::uint64_t;
using Topology = std::vector<ComponentMask>;

std::string matching_key(const Matching& matching) {
    std::ostringstream ss;
    for (auto [a, b] : matching) ss << a << '-' << b << ';';
    return ss.str();
}

int checked_power_of_three(int exponent) {
    std::int64_t value = 1;
    for (int i = 0; i < exponent; ++i) {
        value *= 3;
        if (value > std::numeric_limits<int>::max()) {
            throw std::runtime_error("the island ring is too large for the integer coloring encoding");
        }
    }
    return static_cast<int>(value);
}

std::vector<int> powers_of_three(int length) {
    std::vector<int> powers(length, 1);
    for (int i = 1; i < length; ++i) {
        const std::int64_t next = static_cast<std::int64_t>(powers[i - 1]) * 3;
        if (next > std::numeric_limits<int>::max()) {
            throw std::runtime_error("the island ring is too large for the integer coloring encoding");
        }
        powers[i] = static_cast<int>(next);
    }
    return powers;
}

std::vector<int> decode_coloring(int code, int length) {
    std::vector<int> coloring(length);
    for (int i = 0; i < length; ++i) {
        coloring[i] = code % 3;
        code /= 3;
    }
    return coloring;
}

int canonical_color_code(int code, int length, const std::vector<int>& powers) {
    std::array<int, 3> renamed{{-1, -1, -1}};
    int next_name = 0;
    int canonical = 0;
    for (int position = 0; position < length; ++position) {
        const int color = code % 3;
        code /= 3;
        if (renamed[color] == -1) renamed[color] = next_name++;
        canonical += renamed[color] * powers[position];
    }
    return canonical;
}

struct IslandGraph {
    int ring_edges = 0;
    int dummy_edges = 0;
    int edge_count = 0;
    std::vector<std::vector<int>> actual_incident;
    std::vector<std::vector<int>> edge_vertices;
    std::vector<int> non_ring_actual_edges;
};

IslandGraph build_island_graph(const Island& island) {
    IslandGraph graph;
    graph.ring_edges = island.ring_edge_count();
    graph.dummy_edges = island.degree_two_vertices;
    graph.edge_count = island.edge_count();
    graph.actual_incident.resize(island.incident_edges.size());
    graph.edge_vertices.resize(graph.edge_count);
    const int dummy_begin = graph.ring_edges;
    const int dummy_end = dummy_begin + graph.dummy_edges;
    for (int v = 0; v < static_cast<int>(island.incident_edges.size()); ++v) {
        for (int e : island.incident_edges[v]) {
            if (e < 0 || e >= graph.edge_count) throw std::runtime_error("island edge out of range");
            if (dummy_begin <= e && e < dummy_end) continue;
            graph.actual_incident[v].push_back(e);
            graph.edge_vertices[e].push_back(v);
        }
    }
    for (int e = dummy_end; e < graph.edge_count; ++e) graph.non_ring_actual_edges.push_back(e);
    return graph;
}

// Enumerate one representative from every global permutation orbit of edge colors.
// The visitor returns true to stop the search early.
bool visit_boundary_colorings(const Island& island, const std::set<int>& deleted,
                              const std::function<bool(int)>& visitor) {
    const IslandGraph graph = build_island_graph(island);
    const std::vector<int> boundary_powers = powers_of_three(graph.ring_edges);
    std::vector<int> edge_color(graph.edge_count, -1);
    std::vector<int> variable_edges;
    for (int e = 0; e < graph.edge_count; ++e) {
        if (graph.ring_edges <= e && e < graph.ring_edges + graph.dummy_edges) continue;
        if (deleted.contains(e)) continue;
        variable_edges.push_back(e);
    }
    std::stable_sort(variable_edges.begin(), variable_edges.end(), [&](int a, int b) {
        return graph.edge_vertices[a].size() > graph.edge_vertices[b].size();
    });

    std::vector<int> deleted_count(graph.actual_incident.size(), 0);
    for (int e : deleted) {
        for (int v : graph.edge_vertices[e]) ++deleted_count[v];
    }

    auto locally_valid = [&](int edge, int color) {
        for (int v : graph.edge_vertices[edge]) {
            const bool equal_required = deleted_count[v] == 1;
            for (int other : graph.actual_incident[v]) {
                if (other == edge || deleted.contains(other) || edge_color[other] == -1) continue;
                if (equal_required) {
                    if (edge_color[other] != color) return false;
                } else {
                    if (edge_color[other] == color) return false;
                }
            }
        }
        return true;
    };

    std::function<bool(std::size_t, int)> dfs = [&](std::size_t index,
                                                     int maximum_color_used) -> bool {
        if (index == variable_edges.size()) {
            int boundary_code = 0;
            for (int e = 0; e < graph.ring_edges; ++e) {
                if (edge_color[e] < 0) return false;
                boundary_code += edge_color[e] * boundary_powers[e];
            }
            boundary_code = canonical_color_code(boundary_code, graph.ring_edges,
                                                  boundary_powers);
            return visitor(boundary_code);
        }
        const int edge = variable_edges[index];
        // Canonical augmentation: after colors 0,...,m have appeared, the next edge
        // may use only 0,...,min(2,m+1).  This keeps exactly one representative of
        // each global S_3 color orbit.
        const int maximum_allowed = std::min(2, maximum_color_used + 1);
        for (int color = 0; color <= maximum_allowed; ++color) {
            if (!locally_valid(edge, color)) continue;
            edge_color[edge] = color;
            if (dfs(index + 1, std::max(maximum_color_used, color))) return true;
            edge_color[edge] = -1;
        }
        return false;
    };
    return dfs(0, -1);
}

std::unordered_set<int> enumerate_boundary_colorings(const Island& island,
                                                      const std::set<int>& deleted = {}) {
    std::unordered_set<int> result;
    visit_boundary_colorings(island, deleted, [&](int code) {
        result.insert(code);
        return false;
    });
    return result;
}

std::vector<Matching> compute_planar_half_kempe_templates(int n);

const std::vector<Matching>& planar_half_kempe_templates(int n) {
    if (n < 0) throw std::invalid_argument("negative boundary size");
    thread_local std::vector<std::optional<std::vector<Matching>>> cache;
    if (static_cast<int>(cache.size()) <= n) cache.resize(n + 1);
    if (!cache[n].has_value()) cache[n] = compute_planar_half_kempe_templates(n);
    return *cache[n];
}

std::vector<Matching> compute_planar_half_kempe_templates(int n) {
    if (n == 0 || n == 1) return {Matching{}};
    std::vector<Matching> returned;
    std::unordered_set<std::string> all_keys;

    for (int q = n / 2; q >= 1; --q) {
        const int endpoint_count = 2 * q;
        std::vector<int> subset;
        std::function<void(int)> choose = [&](int start) {
            if (static_cast<int>(subset.size()) == endpoint_count) {
                std::vector<bool> selected(n, false);
                for (int x : subset) selected[x] = true;
                bool admissible = true;
                for (int i = 0; i < n; ++i) {
                    if (!selected[i] && !selected[(i + 1) % n]) {
                        admissible = false;
                        break;
                    }
                }
                if (!admissible) return;
                std::vector<int> local(endpoint_count);
                std::iota(local.begin(), local.end(), 0);
                for (const Matching& m0 : noncrossing_perfect_matchings(local)) {
                    Matching m;
                    for (auto [a, b] : m0) m.emplace(subset[a], subset[b]);
                    const std::string key = matching_key(m);
                    all_keys.insert(key);
                    std::vector<int> unmatched;
                    std::vector<bool> matched(n, false);
                    for (auto [a, b] : m) matched[a] = matched[b] = true;
                    for (int i = 0; i < n; ++i) if (!matched[i]) unmatched.push_back(i);
                    bool redundant = false;
                    for (std::size_t a = 0; a < unmatched.size() && !redundant; ++a) {
                        for (std::size_t b = a + 1; b < unmatched.size(); ++b) {
                            Matching extension = m;
                            extension.emplace(unmatched[a], unmatched[b]);
                            if (all_keys.contains(matching_key(extension))) {
                                redundant = true;
                                break;
                            }
                        }
                    }
                    if (!redundant) returned.push_back(std::move(m));
                }
                return;
            }
            for (int i = start; i < n; ++i) {
                subset.push_back(i);
                choose(i + 1);
                subset.pop_back();
            }
        };
        choose(0);
    }
    std::unordered_set<std::string> seen;
    std::vector<Matching> unique;
    unique.reserve(returned.size());
    for (auto& matching : returned) {
        if (seen.insert(matching_key(matching)).second) unique.push_back(std::move(matching));
    }
    return unique;
}

class TopologyCache {
public:
    explicit TopologyCache(const Island& island) : island_(island), ring_edges_(island.ring_edge_count()) {
        if (ring_edges_ >= 63) {
            throw std::runtime_error("too many ring edges for Kempe-component bit masks");
        }
    }

    const std::vector<Topology>& get(ComponentMask selected) {
        if (auto it = cache_.find(selected); it != cache_.end()) return it->second;
        std::vector<std::vector<Topology>> ring_options;
        int offset = 0;
        for (int ring_size : island_.ring_sizes) {
            std::vector<int> positions;
            for (int i = 0; i < ring_size; ++i) {
                const int position = offset + i;
                if (((selected >> position) & 1U) != 0) positions.push_back(position);
            }
            std::vector<Topology> options;
            if (positions.empty()) {
                options.push_back({});
            } else {
                const auto& templates = planar_half_kempe_templates(
                    static_cast<int>(positions.size()));
                for (const Matching& local : templates) {
                    Topology components;
                    std::vector<bool> matched(positions.size(), false);
                    for (auto [a, b] : local) {
                        components.push_back((ComponentMask{1} << positions[a]) |
                                             (ComponentMask{1} << positions[b]));
                        matched[a] = matched[b] = true;
                    }
                    for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
                        if (!matched[i]) components.push_back(ComponentMask{1} << positions[i]);
                    }
                    options.push_back(std::move(components));
                }
            }
            ring_options.push_back(std::move(options));
            offset += ring_size;
        }

        std::vector<Topology> products(1);
        for (const auto& options : ring_options) {
            std::vector<Topology> next;
            for (const Topology& prefix : products) {
                for (const Topology& option : options) {
                    Topology combined = prefix;
                    combined.insert(combined.end(), option.begin(), option.end());
                    next.push_back(std::move(combined));
                }
            }
            products = std::move(next);
        }
        return cache_.emplace(selected, std::move(products)).first->second;
    }

private:
    const Island& island_;
    int ring_edges_ = 0;
    std::unordered_map<ComponentMask, std::vector<Topology>> cache_;
};

class SemiconsistentColoringSet {
public:
    SemiconsistentColoringSet(int ring_edges, const std::unordered_set<int>& extendable)
        : ring_edges_(ring_edges), powers_(powers_of_three(ring_edges)),
          total_(checked_power_of_three(ring_edges)), code_to_orbit_(total_, -1) {
        std::vector<int> representative_to_orbit(total_, -1);
        for (int code = 0; code < total_; ++code) {
            const int representative = canonical_color_code(code, ring_edges_, powers_);
            int& orbit = representative_to_orbit[representative];
            if (orbit == -1) {
                orbit = static_cast<int>(representatives_.size());
                representatives_.push_back(representative);
            }
            code_to_orbit_[code] = orbit;
        }
        active_.assign(representatives_.size(), 1);
        for (int code : extendable) {
            if (0 <= code && code < total_) active_[code_to_orbit_[code]] = 0;
        }
    }

    [[nodiscard]] int ring_edges() const noexcept { return ring_edges_; }
    [[nodiscard]] int total() const noexcept { return total_; }
    [[nodiscard]] const std::vector<int>& powers() const noexcept { return powers_; }
    [[nodiscard]] const std::vector<int>& representatives() const noexcept {
        return representatives_;
    }
    [[nodiscard]] bool active_orbit(std::size_t orbit) const noexcept {
        return active_[orbit] != 0;
    }
    void remove_orbit(std::size_t orbit) noexcept { active_[orbit] = 0; }
    [[nodiscard]] bool contains(int code) const noexcept {
        return 0 <= code && code < total_ && active_[code_to_orbit_[code]] != 0;
    }
    [[nodiscard]] bool empty() const noexcept {
        return std::none_of(active_.begin(), active_.end(), [](std::uint8_t value) {
            return value != 0;
        });
    }

private:
    int ring_edges_ = 0;
    std::vector<int> powers_;
    int total_ = 0;
    std::vector<int> code_to_orbit_;
    std::vector<int> representatives_;
    std::vector<std::uint8_t> active_;
};

bool topology_closure_is_active(const SemiconsistentColoringSet& active,
                                int code, const std::vector<int>& coloring,
                                int x, int y, const Topology& topology) {
    if (topology.size() >= 63) return false;
    std::vector<int> deltas;
    deltas.reserve(topology.size());
    for (ComponentMask component : topology) {
        int delta = 0;
        ComponentMask remaining = component;
        while (remaining != 0) {
            const int position = std::countr_zero(remaining);
            remaining &= remaining - 1;
            if (coloring[position] == x) {
                delta += (y - x) * active.powers()[position];
            } else if (coloring[position] == y) {
                delta += (x - y) * active.powers()[position];
            } else {
                throw std::runtime_error("Kempe component contains an unselected color");
            }
        }
        deltas.push_back(delta);
    }

    // Test the single-component switches first.  Most impossible topologies are
    // rejected here, before the complete subcube is traversed.
    for (int delta : deltas) {
        if (!active.contains(code + delta)) return false;
    }

    const std::uint64_t subset_count = std::uint64_t{1} << deltas.size();
    int changed_code = code;
    std::uint64_t previous_gray = 0;
    for (std::uint64_t index = 1; index < subset_count; ++index) {
        const std::uint64_t gray = index ^ (index >> 1);
        const std::uint64_t changed_bit = gray ^ previous_gray;
        const int component = std::countr_zero(changed_bit);
        if (((gray >> component) & 1U) != 0) changed_code += deltas[component];
        else changed_code -= deltas[component];
        previous_gray = gray;
        if (!active.contains(changed_code)) return false;
    }
    return true;
}

SemiconsistentColoringSet maximal_semiconsistent_set(
    const Island& island, const std::unordered_set<int>& extendable) {
    SemiconsistentColoringSet active(island.ring_edge_count(), extendable);
    TopologyCache topologies(island);
    const std::array<std::pair<int, int>, 3> pairs{{{0, 1}, {0, 2}, {1, 2}}};

    bool changed;
    do {
        changed = false;
        const auto& representatives = active.representatives();
        for (std::size_t orbit = 0; orbit < representatives.size(); ++orbit) {
            if (!active.active_orbit(orbit)) continue;
            const int code = representatives[orbit];
            const std::vector<int> coloring = decode_coloring(code, active.ring_edges());
            bool valid_for_all_pairs = true;
            for (auto [x, y] : pairs) {
                ComponentMask selected = 0;
                for (int position = 0; position < active.ring_edges(); ++position) {
                    if (coloring[position] == x || coloring[position] == y) {
                        selected |= ComponentMask{1} << position;
                    }
                }
                bool topology_exists = false;
                for (const Topology& topology : topologies.get(selected)) {
                    if (topology_closure_is_active(active, code, coloring, x, y,
                                                   topology)) {
                        topology_exists = true;
                        break;
                    }
                }
                if (!topology_exists) {
                    valid_for_all_pairs = false;
                    break;
                }
            }
            if (!valid_for_all_pairs) {
                // Chaotic fixed-point iteration is safe here: the operator is
                // monotone and active colorings are only removed.
                active.remove_orbit(orbit);
                changed = true;
            }
        }
    } while (changed);
    return active;
}

struct FaceData {
    std::vector<std::set<int>> face_edges;
    std::vector<std::array<int, 2>> edge_faces;
    std::set<int> boundary_faces;
};

FaceData island_faces(const Island& island) {
    const IslandGraph graph = build_island_graph(island);
    // Map original edge ids to a compact actual-edge id while preserving original ids in face sets.
    std::vector<std::pair<int, int>> endpoints(graph.edge_count, {-1, -1});
    int leaf_count = 0;
    const int internal_count = static_cast<int>(graph.actual_incident.size());
    for (int e = 0; e < graph.edge_count; ++e) {
        if (graph.ring_edges <= e && e < graph.ring_edges + graph.dummy_edges) continue;
        const auto& vs = graph.edge_vertices[e];
        if (vs.size() == 2) endpoints[e] = {vs[0], vs[1]};
        else if (vs.size() == 1) endpoints[e] = {vs[0], internal_count + leaf_count++};
        else throw std::runtime_error("actual island edge does not have one or two endpoints");
    }
    const int vertex_count = internal_count + leaf_count;
    std::vector<std::vector<int>> rotations(vertex_count);
    for (int v = 0; v < internal_count; ++v) rotations[v] = graph.actual_incident[v];
    int leaf = internal_count;
    for (int e = 0; e < graph.ring_edges; ++e) rotations[leaf++] = {e};

    // Create two directed half-edges per actual edge. The head endpoint follows endpoints.second/first.
    struct Half { int head, rev, succ; int edge; };
    std::vector<Half> halves;
    std::unordered_map<long long, int> half_at;
    auto key = [](int vertex, int edge) { return (static_cast<long long>(vertex) << 32) ^ edge; };
    for (int e = 0; e < graph.edge_count; ++e) {
        if (endpoints[e].first == -1) continue;
        int a = endpoints[e].first, b = endpoints[e].second;
        int ab = static_cast<int>(halves.size());
        int ba = ab + 1;
        halves.push_back({b, ba, -1, e});
        halves.push_back({a, ab, -1, e});
        half_at[key(b, e)] = ab;
        half_at[key(a, e)] = ba;
    }
    for (int v = 0; v < vertex_count; ++v) {
        auto& rotation = rotations[v];
        for (int i = 0; i < static_cast<int>(rotation.size()); ++i) {
            int e = rotation[i];
            int next = rotation[(i + 1) % rotation.size()];
            halves[half_at.at(key(v, e))].succ = half_at.at(key(v, next));
        }
    }
    std::vector<bool> visited(halves.size(), false);
    FaceData data;
    data.edge_faces.assign(graph.edge_count, {-1, -1});
    for (int h = 0; h < static_cast<int>(halves.size()); ++h) {
        if (visited[h]) continue;
        int current = h;
        std::set<int> edges;
        do {
            visited[current] = true;
            edges.insert(halves[current].edge);
            current = halves[halves[current].succ].rev;
        } while (current != h);
        int face = static_cast<int>(data.face_edges.size());
        data.face_edges.push_back(edges);
        for (int e : edges) {
            if (data.edge_faces[e][0] == -1) data.edge_faces[e][0] = face;
            else data.edge_faces[e][1] = face;
        }
    }

    int offset = 0;
    for (int size : island.ring_sizes) {
        std::set<int> ring;
        for (int e = offset; e < offset + size; ++e) ring.insert(e);
        int found = -1;
        for (int f = 0; f < static_cast<int>(data.face_edges.size()); ++f) {
            if (std::includes(data.face_edges[f].begin(), data.face_edges[f].end(),
                              ring.begin(), ring.end())) {
                found = f;
                break;
            }
        }
        if (found != -1) data.boundary_faces.insert(found);
        offset += size;
    }
    return data;
}

bool deletable_size_four(const std::set<int>& deleted, const FaceData& faces) {
    for (int f = 0; f < static_cast<int>(faces.face_edges.size()); ++f) {
        if (faces.boundary_faces.contains(f)) continue;
        int count = 0;
        for (int e : deleted) if (faces.face_edges[f].contains(e)) ++count;
        if (count >= 3) return true;
    }
    for (int f1 = 0; f1 < static_cast<int>(faces.face_edges.size()); ++f1) {
        if (faces.boundary_faces.contains(f1)) continue;
        for (int f2 = f1 + 1; f2 < static_cast<int>(faces.face_edges.size()); ++f2) {
            if (faces.boundary_faces.contains(f2)) continue;
            bool share_edge = false;
            for (int e : faces.face_edges[f1]) if (faces.face_edges[f2].contains(e)) {
                share_edge = true;
                break;
            }
            if (!share_edge) continue;
            bool all_incident = true;
            for (int e : deleted) {
                if (!faces.face_edges[f1].contains(e) && !faces.face_edges[f2].contains(e)) {
                    all_incident = false;
                    break;
                }
            }
            if (all_incident) return true;
        }
    }
    return false;
}

std::vector<std::set<int>> deletable_sets(const Island& island) {
    const IslandGraph graph = build_island_graph(island);
    const FaceData faces = island_faces(island);
    std::vector<std::set<int>> result;
    const auto& edges = graph.non_ring_actual_edges;
    for (int size = 1; size <= 4; ++size) {
        std::vector<int> chosen;
        std::function<void(int)> dfs = [&](int start) {
            if (static_cast<int>(chosen.size()) == size) {
                std::set<int> deleted(chosen.begin(), chosen.end());
                for (const auto& incident : graph.actual_incident) {
                    int count = 0;
                    for (int e : incident) if (deleted.contains(e)) ++count;
                    if (count == 2) return;
                }
                if (size == 4 && !deletable_size_four(deleted, faces)) return;
                result.push_back(std::move(deleted));
                return;
            }
            for (int i = start; i < static_cast<int>(edges.size()); ++i) {
                chosen.push_back(edges[i]);
                dfs(i + 1);
                chosen.pop_back();
            }
        };
        dfs(0);
    }
    return result;
}

}  // namespace

std::vector<std::set<std::pair<int, int>>> noncrossing_perfect_matchings(
    const std::vector<int>& ordered_set) {
    if (ordered_set.empty()) return {Matching{}};
    if (ordered_set.size() % 2 != 0) return {};
    const int first = ordered_set.front();
    std::vector<Matching> result;
    for (std::size_t j_index = 1; j_index < ordered_set.size(); ++j_index) {
        const int j = ordered_set[j_index];
        std::vector<int> inside(ordered_set.begin() + 1, ordered_set.begin() + j_index);
        std::vector<int> outside(ordered_set.begin() + j_index + 1, ordered_set.end());
        if (inside.size() % 2 != 0) continue;
        for (const Matching& in : noncrossing_perfect_matchings(inside)) {
            for (const Matching& out : noncrossing_perfect_matchings(outside)) {
                Matching matching = in;
                matching.insert(out.begin(), out.end());
                matching.emplace(std::min(first, j), std::max(first, j));
                result.push_back(std::move(matching));
            }
        }
    }
    std::unordered_set<std::string> seen;
    std::vector<Matching> unique;
    for (auto& matching : result) {
        if (seen.insert(matching_key(matching)).second) unique.push_back(std::move(matching));
    }
    return unique;
}

std::vector<std::set<std::pair<int, int>>> get_planar_half_kempes(int n) {
    const auto& cached = planar_half_kempe_templates(n);
    return {cached.begin(), cached.end()};
}

ReducibilityResult check_semi_reducibility(const Island& island,
                                           bool search_c_reductions) {
    const auto extendable = enumerate_boundary_colorings(island);
    const auto active = maximal_semiconsistent_set(island, extendable);
    if (active.empty()) {
        return {.semi_d_reducible = true,
                .semi_c_reducible = false,
                .deleting_edges = {}};
    }
    if (!search_c_reductions) return {};

    for (const std::set<int>& deleted : deletable_sets(island)) {
        // A set F witnesses semi-C-reducibility precisely when no boundary
        // coloring modulo F lies in the maximal semi-consistent obstruction.
        // Stop as soon as an intersection is found; most candidate F are then
        // rejected after only a small fraction of their colorings are explored.
        const bool intersects = visit_boundary_colorings(
            island, deleted, [&](int code) { return active.contains(code); });
        if (!intersects) {
            return {.semi_d_reducible = false,
                    .semi_c_reducible = true,
                    .deleting_edges = {deleted.begin(), deleted.end()}};
        }
    }
    return {};
}

}  // namespace apex

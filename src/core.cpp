#include "apex/apex.hpp"

#include <deque>
#include <iomanip>
#include <numeric>
#include <sstream>

namespace apex {
namespace {

void set_error(std::string* error, const std::string& message) {
    if (error != nullptr) *error = message;
}

struct TraversalForm {
    std::string text;
    std::vector<int> vertex_to_canonical;
    std::vector<int> dart_to_canonical;
};

struct TraversalContext {
    std::vector<VertexId> vertices;
    std::vector<DartId> darts;
    std::vector<int> degree;
    std::vector<std::uint8_t> boundary;
};

TraversalContext make_traversal_context(const Embedding& z) {
    TraversalContext summary;
    summary.vertices = z.vertices();
    summary.darts = z.dart_ids();
    summary.degree.assign(z.vertex_alive.size(), 0);
    summary.boundary.assign(z.vertex_alive.size(), 0);
    for (DartId e = 0; e < static_cast<DartId>(z.darts.size()); ++e) {
        if (!z.darts[e].alive) continue;
        const VertexId v = z.darts[e].head;
        ++summary.degree[v];
        if (z.darts[e].succ == kNil || z.darts[e].pred == kNil) {
            summary.boundary[v] = 1;
        }
    }
    return summary;
}

TraversalForm traversal_form(const Embedding& z, DartId root,
                             std::span<const DartId> marked_darts,
                             std::optional<VertexId> marked_vertex,
                             const TraversalContext& summary) {
    TraversalForm out;
    out.vertex_to_canonical.assign(z.vertex_alive.size(), kNil);
    out.dart_to_canonical.assign(z.darts.size(), kNil);
    std::vector<VertexId> canonical_vertices;
    std::vector<DartId> canonical_darts;
    std::deque<DartId> queue;

    auto add_dart = [&](DartId e) {
        if (e == kNil || e < 0 || static_cast<std::size_t>(e) >= z.darts.size() ||
            !z.darts[e].alive) {
            return;
        }
        if (out.dart_to_canonical[e] == kNil) {
            out.dart_to_canonical[e] = static_cast<int>(canonical_darts.size());
            canonical_darts.push_back(e);
            queue.push_back(e);
        }
    };
    auto add_vertex = [&](VertexId v) {
        if (v == kNil || v < 0 || static_cast<std::size_t>(v) >= z.vertex_alive.size() ||
            !z.vertex_alive[v]) {
            return;
        }
        if (out.vertex_to_canonical[v] == kNil) {
            out.vertex_to_canonical[v] = static_cast<int>(canonical_vertices.size());
            canonical_vertices.push_back(v);
        }
    };

    add_dart(root);
    while (!queue.empty()) {
        const DartId e = queue.front();
        queue.pop_front();
        const Dart& d = z.darts[e];
        add_vertex(d.head);
        add_dart(d.rev);
        add_dart(d.succ);
        add_dart(d.pred);
    }

    // Valid inputs are connected, but make canonicalization total for diagnostics too.
    for (DartId e : summary.darts) {
        if (out.dart_to_canonical[e] == kNil) add_dart(e);
        while (!queue.empty()) {
            const DartId f = queue.front();
            queue.pop_front();
            const Dart& d = z.darts[f];
            add_vertex(d.head);
            add_dart(d.rev);
            add_dart(d.succ);
            add_dart(d.pred);
        }
    }
    for (VertexId v : summary.vertices) add_vertex(v);

    std::ostringstream ss;
    ss << (z.kind == EmbeddingKind::PseudoEmbedding ? 'E' : 'T') << ';';
    ss << "V" << canonical_vertices.size() << ';';
    for (VertexId v : canonical_vertices) {
        const auto r = z.degree_range[v];
        ss << r.lower << ',' << (r.unbounded() ? -1 : r.upper) << ','
           << summary.degree[v] << ',' << static_cast<int>(summary.boundary[v]) << ';';
    }
    ss << "D" << canonical_darts.size() << ';';
    for (DartId e : canonical_darts) {
        const Dart& d = z.darts[e];
        auto cv = [&](VertexId v) { return v == kNil ? -1 : out.vertex_to_canonical[v]; };
        auto cd = [&](DartId f) { return f == kNil ? -1 : out.dart_to_canonical[f]; };
        ss << cv(d.head) << ',' << cd(d.rev) << ',' << cd(d.succ) << ',' << cd(d.pred) << ';';
    }
    if (marked_vertex.has_value()) {
        ss << "MV=" << out.vertex_to_canonical[*marked_vertex] << ';';
    }
    if (!marked_darts.empty()) {
        ss << "MD=";
        for (DartId e : marked_darts) ss << out.dart_to_canonical[e] << ',';
        ss << ';';
    }
    out.text = ss.str();
    return out;
}

}  // namespace

std::string DegreeRange::str() const {
    std::ostringstream ss;
    ss << '[' << lower << ',';
    if (unbounded()) ss << "inf";
    else ss << upper;
    ss << ']';
    return ss.str();
}

std::vector<VertexId> Embedding::vertices() const {
    std::vector<VertexId> result;
    result.reserve(vertex_alive.size());
    for (VertexId v = 0; v < static_cast<VertexId>(vertex_alive.size()); ++v) {
        if (vertex_alive[v]) result.push_back(v);
    }
    return result;
}

std::vector<DartId> Embedding::dart_ids() const {
    std::vector<DartId> result;
    result.reserve(darts.size());
    for (DartId e = 0; e < static_cast<DartId>(darts.size()); ++e) {
        if (darts[e].alive) result.push_back(e);
    }
    return result;
}

std::vector<DartId> Embedding::darts_at(VertexId v) const {
    std::vector<DartId> result;
    for (DartId e = 0; e < static_cast<DartId>(darts.size()); ++e) {
        if (darts[e].alive && darts[e].head == v) result.push_back(e);
    }
    return result;
}

int Embedding::degree(VertexId v) const {
    return static_cast<int>(darts_at(v).size());
}

VertexId Embedding::tail(DartId e) const {
    if (e == kNil || e < 0 || static_cast<std::size_t>(e) >= darts.size() || !darts[e].alive) {
        return kNil;
    }
    const DartId r = darts[e].rev;
    if (r == kNil || r < 0 || static_cast<std::size_t>(r) >= darts.size() || !darts[r].alive) {
        return kNil;
    }
    return darts[r].head;
}

bool Embedding::is_boundary(VertexId v) const {
    const auto incident = darts_at(v);
    if (incident.empty()) return false;
    return std::any_of(incident.begin(), incident.end(), [&](DartId e) {
        return darts[e].succ == kNil || darts[e].pred == kNil;
    });
}

DartId Embedding::first_dart(VertexId v) const {
    DartId answer = kNil;
    for (DartId e : darts_at(v)) {
        if (darts[e].pred == kNil) {
            if (answer != kNil) throw std::runtime_error("vertex has more than one first dart");
            answer = e;
        }
    }
    return answer;
}

DartId Embedding::last_dart(VertexId v) const {
    DartId answer = kNil;
    for (DartId e : darts_at(v)) {
        if (darts[e].succ == kNil) {
            if (answer != kNil) throw std::runtime_error("vertex has more than one last dart");
            answer = e;
        }
    }
    return answer;
}

bool Embedding::has_loop() const {
    for (DartId e : dart_ids()) {
        if (tail(e) == darts[e].head) return true;
    }
    return false;
}

bool Embedding::connected() const {
    const auto vs = vertices();
    if (vs.empty()) return true;
    std::unordered_set<VertexId> seen;
    std::queue<VertexId> q;
    seen.insert(vs.front());
    q.push(vs.front());
    while (!q.empty()) {
        const VertexId v = q.front();
        q.pop();
        for (DartId e : darts_at(v)) {
            const VertexId w = tail(e);
            if (w != kNil && seen.insert(w).second) q.push(w);
        }
    }
    return seen.size() == vs.size();
}

bool Embedding::validate_basic(std::string* error) const {
    if (degree_range.size() != vertex_alive.size()) {
        set_error(error, "degree_range and vertex_alive sizes differ");
        return false;
    }
    for (VertexId v : vertices()) {
        if (degree_range[v].empty()) {
            set_error(error, "empty degree range at vertex " + std::to_string(v));
            return false;
        }
        if (darts_at(v).empty()) {
            set_error(error, "isolated live vertex " + std::to_string(v));
            return false;
        }
    }
    for (DartId e : dart_ids()) {
        const Dart& d = darts[e];
        if (d.head < 0 || static_cast<std::size_t>(d.head) >= vertex_alive.size() ||
            !vertex_alive[d.head]) {
            set_error(error, "invalid head at dart " + std::to_string(e));
            return false;
        }
        if (d.rev < 0 || static_cast<std::size_t>(d.rev) >= darts.size() ||
            !darts[d.rev].alive || darts[d.rev].rev != e) {
            set_error(error, "rev is not an involution at dart " + std::to_string(e));
            return false;
        }
        if (d.succ != kNil) {
            if (d.succ < 0 || static_cast<std::size_t>(d.succ) >= darts.size() ||
                !darts[d.succ].alive || darts[d.succ].pred != e ||
                darts[d.succ].head != d.head) {
                set_error(error, "bad successor at dart " + std::to_string(e));
                return false;
            }
        }
        if (d.pred != kNil) {
            if (d.pred < 0 || static_cast<std::size_t>(d.pred) >= darts.size() ||
                !darts[d.pred].alive || darts[d.pred].succ != e ||
                darts[d.pred].head != d.head) {
                set_error(error, "bad predecessor at dart " + std::to_string(e));
                return false;
            }
        }
    }
    return true;
}

bool Embedding::validate_single_list(std::string* error) const {
    if (!validate_basic(error)) return false;
    for (VertexId v : vertices()) {
        const auto incident = darts_at(v);
        std::unordered_set<DartId> seen;
        DartId start = kNil;
        if (is_boundary(v)) {
            start = first_dart(v);
            if (start == kNil || last_dart(v) == kNil) {
                set_error(error, "boundary incidence list lacks an end at vertex " + std::to_string(v));
                return false;
            }
            for (DartId e = start; e != kNil; e = darts[e].succ) {
                if (!seen.insert(e).second) {
                    set_error(error, "cycle in an acyclic incidence list");
                    return false;
                }
            }
        } else {
            start = incident.front();
            DartId e = start;
            do {
                if (e == kNil || !seen.insert(e).second) {
                    set_error(error, "broken cyclic incidence list");
                    return false;
                }
                e = darts[e].succ;
            } while (e != start);
        }
        if (seen.size() != incident.size()) {
            set_error(error, "vertex has more than one incidence list at vertex " + std::to_string(v));
            return false;
        }
    }
    return true;
}

bool Embedding::validate_faces(std::string* error) const {
    if (!validate_single_list(error)) return false;
    if (kind != EmbeddingKind::PseudoTriangulationWithDigons) return true;
    for (DartId e0 : dart_ids()) {
        if (darts[e0].succ == kNil) continue;
        const DartId e1 = darts[darts[e0].succ].rev;
        if (e1 == kNil || darts[e1].succ == kNil) {
            set_error(error, "incomplete inner face");
            return false;
        }
        const DartId e2 = darts[darts[e1].succ].rev;
        if (e2 == e0) continue;
        if (e2 == kNil || darts[e2].succ == kNil) {
            set_error(error, "incomplete inner face");
            return false;
        }
        const DartId e3 = darts[darts[e2].succ].rev;
        if (e3 != e0) {
            set_error(error, "inner face is neither a digon nor a triangle");
            return false;
        }
    }
    return true;
}

std::string Embedding::debug_string() const {
    std::ostringstream ss;
    ss << "Embedding(kind=" << (kind == EmbeddingKind::PseudoEmbedding ? "embedding" : "tri-digon")
       << ", |V|=" << vertices().size() << ", |D|=" << dart_ids().size() << ")\n";
    for (VertexId v : vertices()) {
        ss << "  v" << v << " range=" << degree_range[v].str() << " degree=" << degree(v)
           << (is_boundary(v) ? " boundary" : " inner") << '\n';
    }
    for (DartId e : dart_ids()) {
        const auto& d = darts[e];
        ss << "  e" << e << ": head=" << d.head << " tail=" << tail(e) << " rev=" << d.rev
           << " succ=" << d.succ << " pred=" << d.pred << '\n';
    }
    return ss.str();
}

VertexId Embedding::add_vertex(DegreeRange range) {
    const VertexId id = static_cast<VertexId>(vertex_alive.size());
    vertex_alive.push_back(true);
    degree_range.push_back(range);
    return id;
}

DartId Embedding::add_dart(const Dart& dart) {
    const DartId id = static_cast<DartId>(darts.size());
    darts.push_back(dart);
    return id;
}

Morphism identity_morphism(const Embedding& z) {
    Morphism m;
    m.vertex.assign(z.vertex_alive.size(), kNil);
    m.dart.assign(z.darts.size(), kNil);
    for (VertexId v : z.vertices()) m.vertex[v] = v;
    for (DartId e : z.dart_ids()) m.dart[e] = e;
    return m;
}

Morphism compose(const Morphism& second, const Morphism& first) {
    Morphism result;
    result.vertex.assign(first.vertex.size(), kNil);
    result.dart.assign(first.dart.size(), kNil);
    for (std::size_t i = 0; i < first.vertex.size(); ++i) {
        const int mid = first.vertex[i];
        if (mid != kNil && static_cast<std::size_t>(mid) < second.vertex.size()) {
            result.vertex[i] = second.vertex[mid];
        }
    }
    for (std::size_t i = 0; i < first.dart.size(); ++i) {
        const int mid = first.dart[i];
        if (mid != kNil && static_cast<std::size_t>(mid) < second.dart.size()) {
            result.dart[i] = second.dart[mid];
        }
    }
    return result;
}

CompactResult compact(const Embedding& z) {
    CompactResult result;
    result.graph.kind = z.kind;
    result.map.vertex.assign(z.vertex_alive.size(), kNil);
    result.map.dart.assign(z.darts.size(), kNil);
    for (VertexId v : z.vertices()) {
        result.map.vertex[v] = result.graph.add_vertex(z.degree_range[v]);
    }
    for (DartId e : z.dart_ids()) {
        result.map.dart[e] = result.graph.add_dart({});
    }
    for (DartId e : z.dart_ids()) {
        Dart d = z.darts[e];
        d.head = result.map.vertex[d.head];
        d.rev = result.map.dart[d.rev];
        d.succ = d.succ == kNil ? kNil : result.map.dart[d.succ];
        d.pred = d.pred == kNil ? kNil : result.map.dart[d.pred];
        d.alive = true;
        result.graph.darts[result.map.dart[e]] = d;
    }
    return result;
}

DisjointUnionResult disjoint_union(const Embedding& a, const Embedding& b, EmbeddingKind kind) {
    DisjointUnionResult result;
    result.graph.kind = kind;
    result.left.vertex.assign(a.vertex_alive.size(), kNil);
    result.left.dart.assign(a.darts.size(), kNil);
    result.right.vertex.assign(b.vertex_alive.size(), kNil);
    result.right.dart.assign(b.darts.size(), kNil);

    for (VertexId v : a.vertices()) result.left.vertex[v] = result.graph.add_vertex(a.degree_range[v]);
    for (VertexId v : b.vertices()) result.right.vertex[v] = result.graph.add_vertex(b.degree_range[v]);
    for (DartId e : a.dart_ids()) result.left.dart[e] = result.graph.add_dart({});
    for (DartId e : b.dart_ids()) result.right.dart[e] = result.graph.add_dart({});

    auto copy_darts = [&](const Embedding& src, const Morphism& map) {
        for (DartId e : src.dart_ids()) {
            Dart d = src.darts[e];
            d.head = map.vertex[d.head];
            d.rev = map.dart[d.rev];
            d.succ = d.succ == kNil ? kNil : map.dart[d.succ];
            d.pred = d.pred == kNil ? kNil : map.dart[d.pred];
            d.alive = true;
            result.graph.darts[map.dart[e]] = d;
        }
    };
    copy_darts(a, result.left);
    copy_darts(b, result.right);
    return result;
}

Embedding mirror(const Embedding& z) {
    Embedding result = z;
    for (DartId e : result.dart_ids()) std::swap(result.darts[e].succ, result.darts[e].pred);
    return result;
}

std::string canonical_key(const Embedding& z, std::optional<DartId> rooted_dart,
                          std::span<const DartId> marked_darts,
                          std::optional<VertexId> marked_vertex) {
    const auto ds = z.dart_ids();
    if (ds.empty()) {
        std::ostringstream ss;
        ss << (z.kind == EmbeddingKind::PseudoEmbedding ? 'E' : 'T') << ";nodarts;";
        std::vector<DegreeRange> ranges;
        for (VertexId v : z.vertices()) ranges.push_back(z.degree_range[v]);
        std::sort(ranges.begin(), ranges.end());
        for (const auto& r : ranges) ss << r.str() << ';';
        return ss.str();
    }
    const TraversalContext summary = make_traversal_context(z);
    if (rooted_dart.has_value()) {
        return traversal_form(z, *rooted_dart, marked_darts, marked_vertex, summary).text;
    }
    std::string best;
    bool first = true;
    for (DartId e : summary.darts) {
        const std::string candidate =
            traversal_form(z, e, marked_darts, marked_vertex, summary).text;
        if (first || candidate < best) {
            best = candidate;
            first = false;
        }
    }
    return best;
}

std::string canonical_image_key(const EmbeddingImage& image) {
    const auto ds = image.target.dart_ids();
    if (ds.empty()) return canonical_key(image.target);
    const TraversalContext summary = make_traversal_context(image.target);
    std::string best;
    bool first = true;
    for (DartId root : summary.darts) {
        auto form = traversal_form(image.target, root, {}, std::nullopt, summary);
        std::ostringstream ss;
        ss << form.text << "|phiV=";
        for (VertexId v : image.map.vertex) {
            ss << (v == kNil ? -1 : form.vertex_to_canonical[v]) << ',';
        }
        ss << "|phiD=";
        for (DartId e : image.map.dart) {
            ss << (e == kNil ? -1 : form.dart_to_canonical[e]) << ',';
        }
        const std::string candidate = ss.str();
        if (first || candidate < best) {
            best = candidate;
            first = false;
        }
    }
    return best;
}


std::vector<std::uint64_t> canonical_unordered_dart_pair_signatures(
    const Embedding& z) {
    const std::size_t stride = z.darts.size();
    std::vector<std::uint64_t> signatures(
        stride * stride, std::numeric_limits<std::uint64_t>::max());
    const TraversalContext context = make_traversal_context(z);
    if (context.darts.empty()) return signatures;

    std::string best_form;
    bool first = true;
    std::vector<std::vector<int>> canonical_maps;
    for (DartId root : context.darts) {
        TraversalForm form = traversal_form(z, root, {}, std::nullopt, context);
        if (first || form.text < best_form) {
            best_form = form.text;
            canonical_maps.clear();
            canonical_maps.push_back(std::move(form.dart_to_canonical));
            first = false;
        } else if (form.text == best_form) {
            canonical_maps.push_back(std::move(form.dart_to_canonical));
        }
    }

    for (DartId e : context.darts) {
        for (DartId f : context.darts) {
            std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
            for (const auto& map : canonical_maps) {
                std::uint32_t a = static_cast<std::uint32_t>(map[e]);
                std::uint32_t b = static_cast<std::uint32_t>(map[f]);
                if (b < a) std::swap(a, b);
                const std::uint64_t packed =
                    (static_cast<std::uint64_t>(a) << 32) | b;
                best = std::min(best, packed);
            }
            signatures[static_cast<std::size_t>(e) * stride + f] = best;
        }
    }
    return signatures;
}

int Island::ring_edge_count() const {
    return std::accumulate(ring_sizes.begin(), ring_sizes.end(), 0);
}

int Island::edge_count() const {
    int maximum = -1;
    for (const auto& inc : incident_edges) {
        maximum = std::max({maximum, inc[0], inc[1], inc[2]});
    }
    return maximum + 1;
}

std::string Island::canonical_key() const {
    // Edge labels in the format are semantically significant only by class (ring, dummy, other).
    // A small color-refinement canonicalization over the cubic incidence hypergraph is sufficient here.
    std::ostringstream ss;
    ss << "R";
    for (int r : ring_sizes) ss << r << ',';
    ss << "M" << degree_two_vertices << ';';
    std::vector<std::array<int, 3>> rows = incident_edges;
    for (auto& row : rows) {
        // Preserve cyclic order up to rotation, not reflection.
        std::array<std::array<int, 3>, 3> rotations{{
            row,
            {row[1], row[2], row[0]},
            {row[2], row[0], row[1]},
        }};
        row = *std::min_element(rotations.begin(), rotations.end());
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& row : rows) ss << row[0] << ',' << row[1] << ',' << row[2] << ';';
    return ss.str();
}

}  // namespace apex

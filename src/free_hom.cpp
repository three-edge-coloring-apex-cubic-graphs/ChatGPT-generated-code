#include "apex/apex.hpp"

#include <deque>
#include <numeric>

namespace apex {
namespace {

class UnionFind {
public:
    explicit UnionFind(std::size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int root(int x) {
        if (parent_[x] != x) parent_[x] = root(parent_[x]);
        return parent_[x];
    }

    bool same(int a, int b) { return root(a) == root(b); }

    // Match the pseudocode: root(a) becomes a child of root(b).
    int unite_to_second(int a, int b) {
        a = root(a);
        b = root(b);
        if (a == b) return b;
        parent_[a] = b;
        return b;
    }

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

bool degree_constraint_holds(DegreeConstraint constraint, const DegreeRange& source,
                             const DegreeRange& target) {
    switch (constraint) {
        case DegreeConstraint::Intersection:
            return source.intersects(target);
        case DegreeConstraint::Include:
            return source.includes(target);
    }
    return false;
}

bool inner_subdegree_error(const Embedding& z) {
    for (VertexId v : z.vertices()) {
        if (z.is_inner(v) && z.degree(v) < z.degree_range[v].lower) return true;
    }
    return false;
}

std::optional<VertexId> vertex_single_degree_issue(const Embedding& z) {
    for (VertexId v : z.vertices()) {
        const DegreeRange r = z.degree_range[v];
        if (!r.fixed()) continue;
        if (r.lower < z.degree(v)) return v;
        if (z.is_boundary(v) && r.lower == z.degree(v)) return v;
    }
    return std::nullopt;
}

std::optional<std::pair<Embedding, Embedding>> single_out_lower_degree(const Embedding& z) {
    for (VertexId v : z.vertices()) {
        const DegreeRange r = z.degree_range[v];
        if (r.lower < r.upper && r.lower <= z.degree(v)) {
            Embedding a = z;
            Embedding b = z;
            a.degree_range[v].upper = r.lower;
            b.degree_range[v].lower = r.lower + 1;
            return std::make_pair(std::move(a), std::move(b));
        }
    }
    return std::nullopt;
}

std::vector<EmbeddingImage> deduplicate_images(std::vector<EmbeddingImage> images) {
    std::unordered_set<std::string> seen;
    std::vector<EmbeddingImage> result;
    result.reserve(images.size());
    for (auto& image : images) {
        const std::string key = canonical_image_key(image);
        if (seen.insert(key).second) result.push_back(std::move(image));
    }
    return result;
}

}  // namespace

std::optional<Morphism> homomorphism(const Embedding& source, DartId source_root,
                                     const Embedding& target, DartId target_root,
                                     DegreeConstraint constraint) {
    if (source_root < 0 || target_root < 0 ||
        static_cast<std::size_t>(source_root) >= source.darts.size() ||
        static_cast<std::size_t>(target_root) >= target.darts.size() ||
        !source.darts[source_root].alive || !target.darts[target_root].alive) {
        return std::nullopt;
    }

    Morphism phi;
    phi.vertex.assign(source.vertex_alive.size(), kNil);
    phi.dart.assign(source.darts.size(), kNil);
    std::queue<std::pair<DartId, DartId>> queue;
    queue.emplace(source_root, target_root);

    while (!queue.empty()) {
        auto [f, fs] = queue.front();
        queue.pop();
        if (phi.dart[f] != kNil) {
            if (phi.dart[f] != fs) return std::nullopt;
            continue;
        }
        phi.dart[f] = fs;
        const VertexId h = source.darts[f].head;
        const VertexId hs = target.darts[fs].head;
        if (phi.vertex[h] != kNil && phi.vertex[h] != hs) return std::nullopt;
        phi.vertex[h] = hs;
        if (!degree_constraint_holds(constraint, source.degree_range[h], target.degree_range[hs])) {
            return std::nullopt;
        }

        queue.emplace(source.darts[f].rev, target.darts[fs].rev);
        if (source.darts[f].succ != kNil) {
            if (target.darts[fs].succ == kNil) return std::nullopt;
            queue.emplace(source.darts[f].succ, target.darts[fs].succ);
        }
        if (source.darts[f].pred != kNil) {
            if (target.darts[fs].pred == kNil) return std::nullopt;
            queue.emplace(source.darts[f].pred, target.darts[fs].pred);
        }
    }

    // A connected single-list pseudo-embedding should be fully reached. Treat an unreachable source
    // component as failure instead of silently returning a partial map.
    for (VertexId v : source.vertices()) {
        if (phi.vertex[v] == kNil) return std::nullopt;
    }
    for (DartId e : source.dart_ids()) {
        if (phi.dart[e] == kNil) return std::nullopt;
    }
    return phi;
}

std::optional<EmbeddingImage> dart_identification(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests) {
    const int nv = static_cast<int>(z.vertex_alive.size());
    const int nd = static_cast<int>(z.darts.size());
    UnionFind ufv(nv);
    UnionFind ufd(nd);
    std::vector<DegreeRange> merged_degree_range = z.degree_range;
    std::vector<VertexId> head(nd, kNil);
    std::vector<DartId> rev(nd, kNil), succ(nd, kNil), pred(nd, kNil);
    for (DartId e : z.dart_ids()) {
        head[e] = z.darts[e].head;
        rev[e] = z.darts[e].rev;
        succ[e] = z.darts[e].succ;
        pred[e] = z.darts[e].pred;
    }

    std::queue<std::pair<DartId, DartId>> queue;
    for (auto [e, f] : requests) {
        if (e < 0 || f < 0 || e >= nd || f >= nd || !z.darts[e].alive || !z.darts[f].alive) {
            return std::nullopt;
        }
        queue.emplace(e, f);
    }

    while (!queue.empty()) {
        auto [e0, f0] = queue.front();
        queue.pop();
        DartId e = ufd.root(e0);
        DartId f = ufd.root(f0);
        if (e == f) continue;

        VertexId he = ufv.root(head[e]);
        VertexId hf = ufv.root(head[f]);
        if (he != hf) {
            const DegreeRange intersection =
                merged_degree_range[he].intersection(merged_degree_range[hf]);
            if (intersection.empty()) return std::nullopt;
            ufv.unite_to_second(he, hf);
            merged_degree_range[hf] = intersection;
        }

        // Save representative pointers before making f the representative.
        const DartId rev_e = rev[e];
        const DartId rev_f = rev[f];
        // Identifying a dart with the class of its reverse would create a
        // self-reverse dart, hence a loop in every quotient target.
        if (ufd.root(rev_e) == f || ufd.root(rev_f) == e) return std::nullopt;
        const DartId succ_e = succ[e];
        const DartId succ_f = succ[f];
        const DartId pred_e = pred[e];
        const DartId pred_f = pred[f];
        ufd.unite_to_second(e, f);

        queue.emplace(rev_e, rev_f);
        if (succ_e != kNil && succ_f != kNil) queue.emplace(succ_e, succ_f);
        if (pred_e != kNil && pred_f != kNil) queue.emplace(pred_e, pred_f);
        if (succ_e != kNil && succ_f == kNil) succ[f] = succ_e;
        if (pred_e != kNil && pred_f == kNil) pred[f] = pred_e;
    }

    // Collect roots and compact them.
    Embedding target;
    target.kind = z.kind;
    Morphism phi;
    phi.vertex.assign(z.vertex_alive.size(), kNil);
    phi.dart.assign(z.darts.size(), kNil);
    std::unordered_map<int, int> vroot_to_new;
    std::unordered_map<int, int> droot_to_new;

    for (VertexId v : z.vertices()) {
        const int r = ufv.root(v);
        if (!vroot_to_new.contains(r)) {
            vroot_to_new[r] = target.add_vertex({1, kInfinity});
        }
        phi.vertex[v] = vroot_to_new[r];
    }
    for (DartId e : z.dart_ids()) {
        const int r = ufd.root(e);
        if (!droot_to_new.contains(r)) droot_to_new[r] = target.add_dart({});
        phi.dart[e] = droot_to_new[r];
    }

    // Degree range intersections.
    for (VertexId v : z.vertices()) {
        VertexId tv = phi.vertex[v];
        DegreeRange intersection = target.degree_range[tv].intersection(z.degree_range[v]);
        if (intersection.empty()) return std::nullopt;
        target.degree_range[tv] = intersection;
    }

    for (const auto& [root, ne] : droot_to_new) {
        Dart d;
        d.alive = true;
        d.head = vroot_to_new[ufv.root(head[root])];
        d.rev = droot_to_new[ufd.root(rev[root])];
        d.succ = succ[root] == kNil ? kNil : droot_to_new[ufd.root(succ[root])];
        d.pred = pred[root] == kNil ? kNil : droot_to_new[ufd.root(pred[root])];
        target.darts[ne] = d;
    }

    if (target.has_loop()) return std::nullopt;
    std::string error;
    if (!target.validate_single_list(&error)) {
        // The free-homomorphism algorithm is intended for valid pseudo-embeddings. A malformed result
        // means the requested identification cannot be represented in the target class.
        return std::nullopt;
    }
    if (target.kind == EmbeddingKind::PseudoTriangulationWithDigons &&
        !target.validate_faces(&error)) {
        return std::nullopt;
    }
    return EmbeddingImage{std::move(target), std::move(phi)};
}

std::vector<EmbeddingImage> resolve_degree_issues(const Embedding& z) {
    std::vector<EmbeddingImage> result;
    std::queue<EmbeddingImage> queue;
    queue.push({z, identity_morphism(z)});
    std::unordered_set<std::string> queued;

    while (!queue.empty()) {
        EmbeddingImage current = std::move(queue.front());
        queue.pop();
        if (inner_subdegree_error(current.target)) continue;

        if (auto issue = vertex_single_degree_issue(current.target); issue.has_value()) {
            const VertexId v = *issue;
            const DegreeRange r = current.target.degree_range[v];
            if (r.lower < current.target.degree(v)) {
                DartId e = current.target.is_boundary(v) ? current.target.first_dart(v)
                                                        : current.target.darts_at(v).front();
                DartId f = e;
                bool valid = true;
                for (int i = 0; i < r.lower; ++i) {
                    f = current.target.darts[f].succ;
                    if (f == kNil) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    if (auto fixed = dart_identification(current.target, {{e, f}}); fixed.has_value()) {
                        fixed->map = compose(fixed->map, current.map);
                        const std::string key = canonical_image_key(*fixed);
                        if (queued.insert(key).second) queue.push(std::move(*fixed));
                    }
                }
            } else if (current.target.is_boundary(v) && r.lower == current.target.degree(v)) {
                for (auto completion : boundary_completions(current.target, v)) {
                    completion.map = compose(completion.map, current.map);
                    const std::string key = canonical_image_key(completion);
                    if (queued.insert(key).second) queue.push(std::move(completion));
                }
            }
            continue;
        }

        if (auto split = single_out_lower_degree(current.target); split.has_value()) {
            for (Embedding branch : {std::move(split->first), std::move(split->second)}) {
                EmbeddingImage image{std::move(branch), current.map};
                const std::string key = canonical_image_key(image);
                if (queued.insert(key).second) queue.push(std::move(image));
            }
            continue;
        }
        result.push_back(std::move(current));
    }
    return deduplicate_images(std::move(result));
}

std::vector<EmbeddingImage> free_homomorphism(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests) {
    auto identified = dart_identification(z, requests);
    if (!identified.has_value()) return {};
    std::vector<EmbeddingImage> resolved = resolve_degree_issues(identified->target);
    for (auto& image : resolved) image.map = compose(image.map, identified->map);
    return deduplicate_images(std::move(resolved));
}

}  // namespace apex

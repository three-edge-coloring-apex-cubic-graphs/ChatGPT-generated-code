#include "apex/apex.hpp"

#include <deque>
#include <iostream>
#include <numeric>

namespace apex {

std::vector<std::vector<DartId>> get_walks(const Embedding& z) {
    std::vector<bool> visited(z.darts.size(), false);
    std::vector<std::vector<DartId>> walks;
    for (DartId start : z.dart_ids()) {
        if (visited[start]) continue;
        std::vector<DartId> walk;
        DartId current = start;
        do {
            if (current == kNil || visited[current]) {
                if (current != start) throw std::runtime_error("getWalks reached a previously visited dart");
                break;
            }
            walk.push_back(current);
            visited[current] = true;
            if (z.darts[current].succ == kNil) {
                const DartId first = z.first_dart(z.darts[current].head);
                if (first == kNil) throw std::runtime_error("boundary vertex has no first dart");
                current = z.darts[first].rev;
            } else {
                current = z.darts[z.darts[current].succ].rev;
            }
        } while (current != start);
        walks.push_back(std::move(walk));
    }
    return walks;
}

bool is_planar(const Embedding& z) {
    const auto walks = get_walks(z);
    const long long euler = static_cast<long long>(z.vertices().size()) -
                            static_cast<long long>(z.dart_ids().size()) / 2 +
                            static_cast<long long>(walks.size());
    return euler == 2;
}

std::vector<std::vector<DartId>> enum_cycles(const Embedding& z, int max_length) {
    std::vector<bool> visited_vertex(z.vertex_alive.size(), false);
    std::vector<bool> visited_dart(z.darts.size(), false);
    std::vector<std::vector<DartId>> cycles;

    std::function<void(std::vector<DartId>&)> dfs = [&](std::vector<DartId>& path) {
        if (static_cast<int>(path.size()) > max_length) return;
        const DartId e = path.back();
        for (DartId incoming : z.darts_at(z.darts[e].head)) {
            if (incoming == e) continue;
            const DartId f = z.darts[incoming].rev;
            if (f == path.front()) {
                if (path.size() >= 2) cycles.push_back(path);
                continue;
            }
            const VertexId next = z.darts[f].head;
            if (visited_dart[f] || visited_vertex[next]) continue;
            path.push_back(f);
            visited_vertex[next] = true;
            dfs(path);
            visited_vertex[next] = false;
            path.pop_back();
        }
    };

    for (DartId e : z.dart_ids()) {
        std::vector<DartId> path{e};
        visited_vertex[z.darts[e].head] = true;
        dfs(path);
        visited_vertex[z.darts[e].head] = false;
        visited_dart[e] = true;
    }
    return cycles;
}

std::vector<SideLabel> label_darts(const Embedding& z,
                                   const std::vector<DartId>& cycle) {
    std::vector<SideLabel> labels(z.darts.size(), SideLabel::Unknown);
    for (DartId e : cycle) {
        labels[e] = SideLabel::Left;
        labels[z.darts[e].rev] = SideLabel::Right;
    }

    auto propagate = [&](DartId seed, SideLabel label) {
        if (seed == kNil) return;
        std::vector<DartId> stack{seed};
        while (!stack.empty()) {
            DartId e = stack.back();
            stack.pop_back();
            if (e == kNil) continue;
            if (labels[e] != SideLabel::Unknown) {
                // Under the hypothesis used in Claim 7.4 the labels cannot conflict.  For a
                // branch that is not homomorphic to G*, simply stop at an already-labelled dart;
                // this follows the pseudocode's "propagate only to unknown darts" rule.
                continue;
            }
            labels[e] = label;
            for (DartId f : {z.darts[e].succ, z.darts[e].pred, z.darts[e].rev}) {
                if (f != kNil && labels[f] == SideLabel::Unknown) stack.push_back(f);
            }
        }
    };

    const int l = static_cast<int>(cycle.size());
    for (int i = 0; i < l; ++i) {
        const DartId di = cycle[i];
        const DartId next = cycle[(i + 1) % l];
        if (z.darts[di].succ != kNil && z.darts[di].succ != z.darts[next].rev) {
            propagate(z.darts[di].succ, SideLabel::Left);
        }
        const DartId reverse_next = z.darts[next].rev;
        if (z.darts[reverse_next].succ != kNil && z.darts[reverse_next].succ != di) {
            propagate(z.darts[reverse_next].succ, SideLabel::Right);
        }
    }
    return labels;
}

std::pair<int, int> num_separated_vertices(
    const Embedding& z, const std::vector<DartId>& cycle,
    const std::vector<SideLabel>& labels) {
    std::unordered_set<VertexId> on_cycle;
    for (DartId e : cycle) on_cycle.insert(z.darts[e].head);
    std::unordered_set<VertexId> left, right;
    for (DartId e : z.dart_ids()) {
        const VertexId v = z.darts[e].head;
        if (on_cycle.contains(v)) continue;
        if (labels[e] == SideLabel::Left) left.insert(v);
        else if (labels[e] == SideLabel::Right) right.insert(v);
    }
    int left_inner = 0, left_boundary = 0, right_inner = 0, right_boundary = 0;
    for (VertexId v : left) {
        if (z.is_boundary(v)) ++left_boundary;
        else ++left_inner;
    }
    for (VertexId v : right) {
        if (z.is_boundary(v)) ++right_boundary;
        else ++right_inner;
    }
    const int nleft = left_inner + (left_boundary > 0 ? 1 : 0);
    const int nright = right_inner + (right_boundary > 0 ? 1 : 0);
    return {nleft, nright};
}

bool has_separating_cycle(const Embedding& z) {
    for (const auto& cycle : enum_cycles(z, 4)) {
        const auto labels = label_darts(z, cycle);
        const auto [left, right] = num_separated_vertices(z, cycle, labels);
        if (cycle.size() <= 3 && left > 0 && right > 0) return true;
        if (cycle.size() == 4 && left > 2 && right > 2) return true;
    }
    return false;
}

std::optional<std::array<DartId, 4>> find_four_darts(const Embedding& z) {
    for (DartId e0 : z.dart_ids()) {
        if (z.darts[e0].succ == kNil) continue;
        std::array<DartId, 4> e{e0, kNil, kNil, kNil};
        bool found_boundary = false;
        for (int i = 0; i < 3; ++i) {
            if ((i == 1 || i == 2) && z.darts[e[i]].succ == kNil) {
                found_boundary = true;
                break;
            }
            if (z.darts[e[i]].succ == kNil) {
                found_boundary = true;
                break;
            }
            e[i + 1] = z.darts[z.darts[e[i]].succ].rev;
        }
        if (!found_boundary && e[0] != e[2] && e[0] != e[3]) return e;
    }
    return std::nullopt;
}

std::vector<EmbeddingImage> ensure_outer_extension(
    const Embedding& z, const std::array<DartId, 4>& darts) {
    // Algorithm B.4.10 returns the union of the e0=e3 and e0=e2 branches.
    // Keep every generated branch; do not quotient this vector by isomorphism.
    std::vector<EmbeddingImage> result =
        free_homomorphism_and_enforce_single_digon_incidence(z, {{darts[0], darts[3]}});
    auto second = free_homomorphism_and_enforce_single_digon_incidence(z,
                                                                        {{darts[0], darts[2]}});
    result.insert(result.end(), std::make_move_iterator(second.begin()),
                  std::make_move_iterator(second.end()));
    return result;
}

std::vector<EmbeddingImage> make_outer_extension(const Embedding& z) {
    // Algorithm B.4.8 has no visited-state or canonical-state test.  Every
    // generated queue entry is processed independently.
    std::queue<EmbeddingImage> queue;
    queue.push({z, identity_morphism(z)});
    std::vector<EmbeddingImage> result;
    while (!queue.empty()) {
        EmbeddingImage current = std::move(queue.front());
        queue.pop();
        auto four = find_four_darts(current.target);
        if (!four.has_value()) {
            result.push_back(std::move(current));
            continue;
        }
        for (auto image : ensure_outer_extension(current.target, *four)) {
            image.map = compose(image.map, current.map);
            queue.push(std::move(image));
        }
    }
    return result;
}

Embedding free_completion_from_outer_extension(const Embedding& outer_extension) {
    const auto walks = get_walks(outer_extension);
    Embedding z = outer_extension;
    Morphism phi = identity_morphism(outer_extension);
    for (const auto& walk : walks) {
        for (DartId original_e : walk) {
            DartId e0 = phi.dart[original_e];
            if (e0 == kNil || !z.darts[e0].alive) continue;
            const DartId reverse_e0 = z.darts[e0].rev;
            if (z.darts[reverse_e0].pred != kNil) continue;
            if (z.darts[e0].succ == kNil) {
                throw std::runtime_error("outer extension has an unexpected nil successor");
            }
            const DartId e1 = z.darts[z.darts[e0].succ].rev;
            if (z.darts[e1].succ == kNil) {
                add_boundary_darts_directly(z, z.darts[e1].rev, e0);
                continue;
            }
            const DartId e2 = z.darts[z.darts[e1].succ].rev;
            if (z.darts[e2].succ == kNil) {
                Morphism local = link_incidence_list_ends(z, z.darts[e0].rev, e2);
                phi = compose(local, phi);
                continue;
            }
            throw std::runtime_error("unreachable branch in freeCompletionFromOuterExtension");
        }
    }
    z.kind = EmbeddingKind::PseudoTriangulationWithDigons;
    auto compacted = compact(z);
    std::string error;
    if (!compacted.graph.validate_faces(&error)) {
        throw std::runtime_error("invalid free completion: " + error);
    }
    return compacted.graph;
}

Island island_from_free_completion(const Embedding& free_completion) {
    const auto walks = get_walks(free_completion);
    std::vector<int> dart_to_edge(free_completion.darts.size(), -1);
    int next_edge = 0;
    Island island;

    // B.4.13: ring edges first.
    for (const auto& walk : walks) {
        if (free_completion.darts[walk.front()].succ == kNil && walk.size() > 1) {
            for (DartId e : walk) {
                dart_to_edge[e] = next_edge;
                dart_to_edge[free_completion.darts[e].rev] = next_edge;
                ++next_edge;
            }
            island.ring_sizes.push_back(static_cast<int>(walk.size()));
        }
    }

    // B.4.14: dummy/pendant edges next.
    std::unordered_map<int, int> digon_dummy;
    for (int wi = 0; wi < static_cast<int>(walks.size()); ++wi) {
        const auto& walk = walks[wi];
        const DartId e0 = walk.front();
        if (free_completion.darts[e0].succ == kNil && walk.size() == 1) {
            dart_to_edge[e0] = next_edge;
            dart_to_edge[free_completion.darts[e0].rev] = next_edge;
            ++next_edge;
            ++island.degree_two_vertices;
        }
        if (free_completion.darts[e0].succ != kNil && walk.size() == 2) {
            digon_dummy[wi] = next_edge++;
            ++island.degree_two_vertices;
        }
    }

    // B.4.15: all remaining genuine edges.
    for (const auto& walk : walks) {
        if (free_completion.darts[walk.front()].succ == kNil) continue;
        for (DartId e : walk) {
            const DartId r = free_completion.darts[e].rev;
            if (dart_to_edge[e] == -1 && dart_to_edge[r] == -1) {
                dart_to_edge[e] = dart_to_edge[r] = next_edge++;
            }
        }
    }

    // B.4.16: dual vertices.
    for (int wi = 0; wi < static_cast<int>(walks.size()); ++wi) {
        const auto& walk = walks[wi];
        if (free_completion.darts[walk.front()].succ == kNil) continue;
        if (walk.size() == 2) {
            island.incident_edges.push_back(
                {dart_to_edge[walk[0]], dart_to_edge[walk[1]], digon_dummy.at(wi)});
        } else if (walk.size() == 3) {
            island.incident_edges.push_back(
                {dart_to_edge[walk[0]], dart_to_edge[walk[1]], dart_to_edge[walk[2]]});
        } else {
            throw std::runtime_error("free completion has a face of size other than two or three");
        }
    }
    return island;
}

namespace {

struct B4SearchStatistics {
    std::size_t explored_states = 0;
    std::size_t state_memo_hits = 0;
    std::size_t identification_pairs = 0;
    std::size_t impossible_pair_prunes = 0;
    std::size_t equivalent_pair_memo_hits = 0;
    std::size_t free_images = 0;
    std::size_t outer_images = 0;
    std::size_t outer_memo_hits = 0;
};

bool identification_may_have_an_image(const Embedding& z, DartId e, DartId f) {
    if (z.darts[e].rev == f) return false;
    const VertexId e_head = z.darts[e].head;
    const VertexId f_head = z.darts[f].head;
    if (!z.degree_range[e_head].intersects(z.degree_range[f_head])) return false;
    const VertexId e_tail = z.tail(e);
    const VertexId f_tail = z.tail(f);
    if (!z.degree_range[e_tail].intersects(z.degree_range[f_tail])) return false;
    return true;
}

}  // namespace

std::vector<Island> all_hom_images(
    const Embedding& outer_extension,
    const std::vector<RootedConfiguration>& smaller_configurations,
    const B3SearchOptions& options) {
    using IslandList = std::vector<Island>;
    using IslandListPtr = std::shared_ptr<const IslandList>;
    using EmbeddingList = std::vector<Embedding>;
    using EmbeddingListPtr = std::shared_ptr<const EmbeddingList>;

    B4SearchStatistics statistics;
    std::unordered_map<std::string, IslandListPtr> state_cache;
    std::unordered_map<std::string, EmbeddingListPtr> outer_extension_cache;
    std::size_t cached_island_occurrences = 0;
    std::size_t cached_outer_embeddings = 0;
    constexpr std::size_t kMaximumCachedIslandOccurrences = 250'000;
    constexpr std::size_t kMaximumCachedOuterEmbeddings = 250'000;
    constexpr std::size_t kMaximumCachedStates = 50'000;
    constexpr std::size_t kMaximumCachedOuterStates = 50'000;

    // A target-only version of Algorithm B.4.8.  allHomImages never uses the
    // composed morphisms returned by makeOuterExtension, so avoiding their
    // construction is safe.  When a subproblem is repeated, its complete list
    // of terminal targets is replayed; branch multiplicity is not changed.
    std::function<EmbeddingListPtr(const Embedding&)> outer_targets =
        [&](const Embedding& input) -> EmbeddingListPtr {
        std::string cache_key;
        if (options.memoize_outer_extensions) {
            cache_key = canonical_key(input);
            if (auto it = outer_extension_cache.find(cache_key);
                it != outer_extension_cache.end()) {
                ++statistics.outer_memo_hits;
                return it->second;
            }
        }

        auto result = std::make_shared<EmbeddingList>();
        const auto four = find_four_darts(input);
        if (!four.has_value()) {
            result->push_back(input);
        } else {
            for (auto image : ensure_outer_extension(input, *four)) {
                const EmbeddingListPtr child = outer_targets(image.target);
                result->insert(result->end(), child->begin(), child->end());
            }
        }

        if (options.memoize_outer_extensions &&
            outer_extension_cache.size() < kMaximumCachedOuterStates &&
            cached_outer_embeddings + result->size() <= kMaximumCachedOuterEmbeddings) {
            cached_outer_embeddings += result->size();
            outer_extension_cache.emplace(std::move(cache_key), result);
        }
        return result;
    };

    std::function<IslandListPtr(const Embedding&)> recurse =
        [&](const Embedding& input) -> IslandListPtr {
        ++statistics.explored_states;

        std::string state_key;
        if (options.memoize_recursive_states) {
            // Algorithm B.4.1 has no distinguished center.  The recursive
            // state is therefore keyed by the unrooted pseudo-embedding.
            state_key = canonical_key(input);
            if (auto it = state_cache.find(state_key); it != state_cache.end()) {
                ++statistics.state_memo_hits;
                return it->second;
            }
        }

        auto result = std::make_shared<IslandList>();
        if (!blocked_by_reducible_configuration(input, smaller_configurations)) {
            if (is_planar(input) && !has_separating_cycle(input)) {
                Island island = island_from_free_completion(
                    free_completion_from_outer_extension(input));
                if (island.degree_two_vertices <= 3) {
                    result->push_back(std::move(island));
                }
            }

            std::vector<bool> checked(input.darts.size(), false);
            const auto darts = input.dart_ids();
            std::vector<std::uint64_t> pair_signatures;
            if (options.memoize_equivalent_pair_branches) {
                pair_signatures = canonical_unordered_dart_pair_signatures(input);
            }
            std::unordered_map<std::uint64_t, EmbeddingListPtr> pair_branch_cache;
            for (DartId e : darts) {
                for (DartId f : darts) {
                    // The checked-dart tests are exactly Algorithm B.4.1.  The
                    // additional fast rejection below removes only requests
                    // whose free-homomorphic-image set is provably empty.
                    if (e == f || checked[e] || checked[f]) continue;
                    ++statistics.identification_pairs;
                    if (options.prune_impossible_identifications &&
                        !identification_may_have_an_image(input, e, f)) {
                        ++statistics.impossible_pair_prunes;
                        continue;
                    }
                    if (statistics.identification_pairs % 10000 == 0) {
                        std::cerr << "[B4] explored " << statistics.explored_states
                                  << " quotient-state occurrences, "
                                  << statistics.identification_pairs
                                  << " dart-pair requests\n";
                    }

                    EmbeddingListPtr branch_targets;
                    std::uint64_t pair_signature = 0;
                    if (options.memoize_equivalent_pair_branches) {
                        const std::size_t stride = input.darts.size();
                        pair_signature = pair_signatures[
                            static_cast<std::size_t>(e) * stride + f];
                        if (auto it = pair_branch_cache.find(pair_signature);
                            it != pair_branch_cache.end()) {
                            branch_targets = it->second;
                            ++statistics.equivalent_pair_memo_hits;
                        }
                    }

                    if (!branch_targets) {
                        auto computed_targets = std::make_shared<EmbeddingList>();
                        for (auto image :
                             free_homomorphism_and_enforce_single_digon_incidence(
                                 input, {{e, f}})) {
                            ++statistics.free_images;
                            if (options.memoize_outer_extensions) {
                                const EmbeddingListPtr targets =
                                    outer_targets(image.target);
                                computed_targets->insert(computed_targets->end(),
                                                         targets->begin(), targets->end());
                            } else {
                                // Exact literal B.4.8 path, retained for regression
                                // comparisons via --literal-search.
                                for (auto outer : make_outer_extension(image.target)) {
                                    computed_targets->push_back(std::move(outer.target));
                                }
                            }
                        }
                        branch_targets = computed_targets;
                        if (options.memoize_equivalent_pair_branches) {
                            pair_branch_cache.emplace(pair_signature, branch_targets);
                        }
                    }

                    // Every actual pair is processed.  A cache hit replays the
                    // entire branch list, so isomorphic pair occurrences and all
                    // island multiplicities are retained.
                    for (const Embedding& target : *branch_targets) {
                        ++statistics.outer_images;
                        const IslandListPtr child = recurse(target);
                        result->insert(result->end(), child->begin(), child->end());
                    }
                }
                checked[e] = true;
                checked[input.darts[e].rev] = true;
            }
        }

        if (options.memoize_recursive_states &&
            state_cache.size() < kMaximumCachedStates &&
            cached_island_occurrences + result->size() <=
                kMaximumCachedIslandOccurrences) {
            cached_island_occurrences += result->size();
            state_cache.emplace(std::move(state_key), result);
        }
        return result;
    };

    const IslandListPtr final_result = recurse(outer_extension);
    std::cerr << "[B4] completed allHomImages: state_occurrences="
              << statistics.explored_states
              << ", state_memo_hits=" << statistics.state_memo_hits
              << ", pair_requests=" << statistics.identification_pairs
              << ", impossible_pair_prunes="
              << statistics.impossible_pair_prunes
              << ", equivalent_pair_memo_hits="
              << statistics.equivalent_pair_memo_hits
              << ", free_images=" << statistics.free_images
              << ", outer_images=" << statistics.outer_images
              << ", outer_memo_hits=" << statistics.outer_memo_hits
              << ", generated_island_occurrences=" << final_result->size()
              << "\n";
    return *final_result;
}

}  // namespace apex

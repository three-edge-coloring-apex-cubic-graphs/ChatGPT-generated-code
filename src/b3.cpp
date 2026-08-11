#include "apex/apex.hpp"

#include <deque>
#include <numeric>
#include <sstream>

namespace apex {
namespace {

std::string cartwheel_key(const Cartwheel& c) {
    return canonical_key(c.graph, c.spokes.empty() ? std::nullopt
                                                    : std::optional<DartId>(c.spokes.front()),
                         c.spokes, c.center);
}

std::vector<Cartwheel> dedupe_cartwheels(std::vector<Cartwheel> input) {
    std::unordered_set<std::string> seen;
    std::vector<Cartwheel> result;
    result.reserve(input.size());
    for (auto& c : input) {
        const std::string key = cartwheel_key(c);
        if (seen.insert(key).second) result.push_back(std::move(c));
    }
    return result;
}

std::vector<CombinedCartwheel> dedupe_combined_cartwheels(std::vector<CombinedCartwheel> input) {
    std::unordered_set<std::string> seen;
    std::vector<CombinedCartwheel> result;
    for (auto& c : input) {
        std::ostringstream ss;
        ss << cartwheel_key(c.cartwheel) << "|fixed=";
        for (const Rule& r : c.fixed_incoming) {
            for (bool b : r.members) ss << (b ? '1' : '0');
            ss << ':' << r.charge << ';';
        }
        if (seen.insert(ss.str()).second) result.push_back(std::move(c));
    }
    return result;
}

bool is_rotation_lexicographically_smaller(const std::vector<DegreeRange>& degrees) {
    const std::size_t n = degrees.size();
    for (std::size_t shift = 1; shift < n; ++shift) {
        bool smaller = false;
        bool greater = false;
        for (std::size_t i = 0; i < n; ++i) {
            const DegreeRange& a = degrees[(i + shift) % n];
            const DegreeRange& b = degrees[i];
            if (a < b) {
                smaller = true;
                break;
            }
            if (b < a) {
                greater = true;
                break;
            }
        }
        if (smaller && !greater) return true;
    }
    return false;
}

std::vector<int> vertex_distances(const Embedding& z, VertexId source) {
    std::vector<int> distance(z.vertex_alive.size(), std::numeric_limits<int>::max());
    std::queue<VertexId> q;
    distance[source] = 0;
    q.push(source);
    while (!q.empty()) {
        VertexId v = q.front();
        q.pop();
        for (DartId e : z.darts_at(v)) {
            VertexId w = z.tail(e);
            if (distance[w] == std::numeric_limits<int>::max()) {
                distance[w] = distance[v] + 1;
                q.push(w);
            }
        }
    }
    return distance;
}

bool exact(const DegreeRange& r, int d) { return r.lower == d && r.upper == d; }
bool subset(const DegreeRange& r, int a, int b) { return a <= r.lower && r.upper <= b; }

int center_degree(const Cartwheel& c) {
    const auto r = c.graph.degree_range[c.center];
    if (!r.fixed()) throw std::runtime_error("cartwheel center degree is not fixed");
    return r.lower;
}

template <class Callback>
void for_each_ordinary_wheel(int d, Callback&& callback) {
    std::vector<DegreeRange> seed(d, kCartwheelDegrees.front());
    Cartwheel base = generate_cartwheel(d, seed, false);
    std::vector<DegreeRange> degrees(d);
    std::function<void(int, int)> enumerate = [&](int i, int lowest) {
        if (i == d) {
            if (is_rotation_lexicographically_smaller(degrees)) return;
            Cartwheel wheel = base;
            for (int j = 0; j < d; ++j) {
                const VertexId neighbor = wheel.graph.tail(wheel.spokes[j]);
                wheel.graph.degree_range[neighbor] = degrees[j];
            }
            callback(std::move(wheel));
            return;
        }
        for (int j = lowest; j < static_cast<int>(kCartwheelDegrees.size()); ++j) {
            degrees[i] = kCartwheelDegrees[j];
            enumerate(i + 1, lowest);
        }
    };
    for (int j = 0; j < static_cast<int>(kCartwheelDegrees.size()); ++j) {
        degrees[0] = kCartwheelDegrees[j];
        enumerate(1, j);
    }
}

template <class Callback>
void for_each_digon_wheel(int d, Callback&& callback) {
    const int n = d - 1;
    std::vector<DegreeRange> seed(n, kCartwheelDegrees.front());
    Cartwheel base = generate_cartwheel(d, seed, true);
    // A digon gives one repeated spoke. Assign ranges by the distinct neighboring vertices rather than
    // by spoke index.
    std::vector<VertexId> neighbors;
    for (DartId e : base.spokes) {
        VertexId u = base.graph.tail(e);
        if (std::find(neighbors.begin(), neighbors.end(), u) == neighbors.end()) neighbors.push_back(u);
    }
    if (static_cast<int>(neighbors.size()) != n) throw std::runtime_error("bad digon-wheel template");
    std::vector<DegreeRange> degrees(n);
    std::function<void(int)> enumerate = [&](int i) {
        if (i == n) {
            Cartwheel wheel = base;
            for (int j = 0; j < n; ++j) wheel.graph.degree_range[neighbors[j]] = degrees[j];
            callback(std::move(wheel));
            return;
        }
        for (DegreeRange range : kCartwheelDegrees) {
            degrees[i] = range;
            enumerate(i + 1);
        }
    };
    enumerate(0);
}

}  // namespace

Cartwheel generate_cartwheel(int center_degree_value,
                             const std::vector<DegreeRange>& degrees,
                             bool incident_digon) {
    const int n = incident_digon ? center_degree_value - 1 : center_degree_value;
    if (static_cast<int>(degrees.size()) != n) {
        throw std::invalid_argument("generate_cartwheel: wrong degree array length");
    }
    std::vector<std::vector<int>> rotations(n + 1);
    for (int i = 1; i <= n; ++i) rotations[0].push_back(i);
    for (int i = 1; i <= n; ++i) {
        const int next = i < n ? i + 1 : 1;
        const int prev = i > 1 ? i - 1 : n;
        rotations[i] = {next, 0, prev, -1};
    }
    std::vector<std::pair<int, int>> digons;
    if (incident_digon) digons.emplace_back(0, 1);
    Embedding z = from_vertex_rotations(n + 1, rotations, digons,
                                        EmbeddingKind::PseudoTriangulationWithDigons);
    z.degree_range[0] = {center_degree_value, center_degree_value};
    for (int i = 1; i <= n; ++i) z.degree_range[i] = degrees[i - 1];

    std::vector<DartId> spokes;
    DartId start = z.darts_at(0).front();
    DartId e = start;
    do {
        spokes.push_back(e);
        e = z.darts[e].succ;
    } while (e != start);
    if (static_cast<int>(spokes.size()) != center_degree_value) {
        throw std::runtime_error("generated center has wrong degree");
    }
    return {std::move(z), 0, std::move(spokes)};
}

std::vector<Cartwheel> enum_digon_incident_wheels(int center_degree_value) {
    std::vector<Cartwheel> result;
    for_each_digon_wheel(center_degree_value,
                         [&](Cartwheel wheel) { result.push_back(std::move(wheel)); });
    return result;
}

std::vector<Cartwheel> enum_wheels(int center_degree_value) {
    std::vector<Cartwheel> result;
    for_each_ordinary_wheel(center_degree_value,
                            [&](Cartwheel wheel) { result.push_back(std::move(wheel)); });
    return result;
}

std::set<std::pair<VertexId, VertexId>> enum_digons(const Embedding& z) {
    std::set<std::pair<VertexId, VertexId>> result;
    for (DartId e0 : z.dart_ids()) {
        if (z.darts[e0].succ == kNil) continue;
        const DartId e1 = z.darts[z.darts[e0].succ].rev;
        if (z.darts[e1].succ == kNil) continue;
        const DartId e2 = z.darts[z.darts[e1].succ].rev;
        if (e2 == e0) {
            VertexId a = z.darts[e0].head;
            VertexId b = z.tail(e0);
            if (a > b) std::swap(a, b);
            result.emplace(a, b);
        }
    }
    return result;
}

int lower_bound_of_digon_charge(const Cartwheel& cartwheel) {
    const auto distance = vertex_distances(cartwheel.graph, cartwheel.center);
    int charge = 0;
    for (auto [u, w] : enum_digons(cartwheel.graph)) {
        int du = distance[u];
        int dw = distance[w];
        const DegreeRange ru = cartwheel.graph.degree_range[u];
        const DegreeRange rw = cartwheel.graph.degree_range[w];

        if (du == 0 && dw == 1) {
            if (subset(ru, 7, 10)) charge += 4;
            continue;
        }
        if (du == 1 && dw == 0) {
            if (subset(rw, 7, 10)) charge += 4;
            continue;
        }
        if (du == 1 && dw == 1) {
            if (exact(ru, 5) && exact(rw, 5)) charge += 5;
            else if ((exact(ru, 5) && exact(rw, 6)) ||
                     (exact(ru, 6) && exact(rw, 5))) charge += 3;
            else charge += 2;
            continue;
        }
        if ((du == 1 && dw == 2) || (du == 2 && dw == 1)) {
            DegreeRange near = du == 1 ? ru : rw;
            DegreeRange far = du == 1 ? rw : ru;
            if (exact(near, 5) && exact(far, 5)) charge += 4;
            else if ((exact(near, 5) && subset(far, 6, 7)) ||
                     (exact(far, 5) && subset(near, 6, 7))) charge += 2;
        }
    }
    return charge;
}

std::vector<Cartwheel> concrete_degree_except_tail(const Cartwheel& cartwheel) {
    std::vector<Cartwheel> current{cartwheel};
    for (VertexId u : cartwheel.graph.vertices()) {
        const DegreeRange original = cartwheel.graph.degree_range[u];
        if (original.fixed() || original.unbounded()) continue;
        std::vector<Cartwheel> next;
        for (int d = original.lower; d <= original.upper; ++d) {
            for (const Cartwheel& c : current) {
                Cartwheel copy = c;
                copy.graph.degree_range[u] = {d, d};
                next.push_back(std::move(copy));
            }
        }
        current = dedupe_cartwheels(std::move(next));
    }
    return current;
}

std::vector<Cartwheel> update_degree_by_rule(const Cartwheel& cartwheel,
                                             DartId dart, const Rule& rule) {
    DisjointUnionResult joined = disjoint_union(
        cartwheel.graph, rule.graph, EmbeddingKind::PseudoTriangulationWithDigons);
    const DartId e_cartwheel = joined.left.dart[dart];
    const DartId e_rule = joined.right.dart[rule.distinguished];
    std::vector<Cartwheel> result;
    for (auto image : free_homomorphism_and_enforce_single_digon_incidence(
             joined.graph, {{e_cartwheel, e_rule}})) {
        Cartwheel c;
        c.graph = std::move(image.target);
        c.center = image.map.vertex[joined.left.vertex[cartwheel.center]];
        for (DartId e : cartwheel.spokes) {
            c.spokes.push_back(image.map.dart[joined.left.dart[e]]);
        }
        auto concrete = concrete_degree_except_tail(c);
        result.insert(result.end(), std::make_move_iterator(concrete.begin()),
                      std::make_move_iterator(concrete.end()));
    }
    return dedupe_cartwheels(std::move(result));
}

int amount_of_charge_sent(const Cartwheel& cartwheel, DartId dart,
                          const std::vector<Rule>& rules) {
    int amount = 0;
    for (const Rule& rule : rules) {
        if (always_apply(cartwheel.graph, dart, rule)) amount += rule.charge;
    }
    return amount;
}

int amount_of_possible_charge_sent(const Cartwheel& cartwheel, DartId dart,
                                   const std::vector<Rule>& combined_rules) {
    // Algorithm A.9.4 asks for a maximum.  Testing in non-increasing charge order lets us stop at
    // the first rule that is not impossible; this is a large practical reduction for B.3.14.
    thread_local const std::vector<Rule>* cached_rules = nullptr;
    thread_local std::size_t cached_size = 0;
    thread_local std::vector<std::size_t> order;
    if (cached_rules != &combined_rules || cached_size != combined_rules.size()) {
        cached_rules = &combined_rules;
        cached_size = combined_rules.size();
        order.resize(combined_rules.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return combined_rules[a].charge > combined_rules[b].charge;
        });
    }
    for (std::size_t index : order) {
        const Rule& rule = combined_rules[index];
        if (!never_apply(cartwheel.graph, dart, rule)) return rule.charge;
    }
    return 0;
}

int upper_bound_of_charge(const Cartwheel& cartwheel,
                          const std::vector<Rule>& fixed_incoming,
                          const std::vector<Rule>& rules,
                          const std::vector<Rule>& combined_rules) {
    const int d = center_degree(cartwheel);
    int value = 10 * (6 - d);
    const std::size_t fixed = fixed_incoming.size();
    for (std::size_t j = 0; j < cartwheel.spokes.size(); ++j) {
        if (j < fixed) value += fixed_incoming[j].charge;
        else value += amount_of_possible_charge_sent(cartwheel, cartwheel.spokes[j], combined_rules);
        value -= amount_of_charge_sent(cartwheel,
                                       cartwheel.graph.darts[cartwheel.spokes[j]].rev, rules);
    }
    return value;
}

bool prune_by_non_associated_rule(const Cartwheel& cartwheel,
                                  const std::vector<Rule>& fixed_incoming,
                                  const std::vector<Rule>& rules) {
    for (std::size_t j = 0; j < fixed_incoming.size(); ++j) {
        const Rule& combined = fixed_incoming[j];
        for (std::size_t i = 0; i < rules.size(); ++i) {
            const bool associated = i < combined.members.size() && combined.members[i];
            if (!associated && always_apply(cartwheel.graph, cartwheel.spokes[j], rules[i])) {
                return true;
            }
        }
    }
    return false;
}

bool prune_cartwheel(const Cartwheel& cartwheel,
                     const std::vector<Rule>& fixed_incoming,
                     const std::vector<Rule>& rules,
                     const std::vector<Rule>& combined_rules,
                     const std::vector<RootedConfiguration>& configurations) {
    if (prune_by_non_associated_rule(cartwheel, fixed_incoming, rules)) return true;
    if (upper_bound_of_charge(cartwheel, fixed_incoming, rules, combined_rules) <=
        lower_bound_of_digon_charge(cartwheel)) {
        return true;
    }
    if (blocked_by_reducible_configuration_at_center(
            cartwheel.graph, cartwheel.center, configurations)) {
        return true;
    }
    return false;
}

std::vector<CombinedCartwheel> fix_in_rules(
    const Cartwheel& wheel, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations) {
    std::vector<CombinedCartwheel> current{{wheel, {}}};
    for (std::size_t i = 0; i < wheel.spokes.size(); ++i) {
        std::vector<CombinedCartwheel> next;
        for (const CombinedCartwheel& state : current) {
            for (const Rule& combined : combined_rules) {
                for (Cartwheel updated : update_degree_by_rule(
                         state.cartwheel, state.cartwheel.spokes[i], combined)) {
                    std::vector<Rule> fixed = state.fixed_incoming;
                    fixed.push_back(combined);
                    if (prune_cartwheel(updated, fixed, rules, combined_rules, configurations)) {
                        continue;
                    }
                    next.push_back({std::move(updated), std::move(fixed)});
                }
            }
        }
        current = dedupe_combined_cartwheels(std::move(next));
        if (current.empty()) break;
    }
    return current;
}

bool should_refine(const Cartwheel& cartwheel, int spoke_index,
                   const AuxiliaryCover& cover) {
    if (spoke_index < 0 || spoke_index >= static_cast<int>(cartwheel.spokes.size())) {
        throw std::out_of_range("spoke index");
    }
    const DartId outward = cartwheel.graph.darts[cartwheel.spokes[spoke_index]].rev;
    if (!always_apply(cartwheel.graph, outward, cover.base)) return false;
    for (const Rule& rule : cover.cover) {
        if (always_apply(cartwheel.graph, outward, rule)) return false;
    }
    return true;
}

std::vector<Cartwheel> refinement(const Cartwheel& cartwheel, int spoke_index,
                                  const AuxiliaryCover& cover) {
    std::vector<Cartwheel> result;
    const DartId outward = cartwheel.graph.darts[cartwheel.spokes[spoke_index]].rev;
    for (const Rule& rule : cover.cover) {
        auto updated = update_degree_by_rule(cartwheel, outward, rule);
        result.insert(result.end(), std::make_move_iterator(updated.begin()),
                      std::make_move_iterator(updated.end()));
    }
    return dedupe_cartwheels(std::move(result));
}

std::vector<CombinedCartwheel> fix_out_rules(
    const std::vector<CombinedCartwheel>& fixed, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers) {
    std::queue<CombinedCartwheel> queue;
    for (const auto& c : fixed) queue.push(c);
    std::vector<CombinedCartwheel> result;
    std::unordered_set<std::string> processed;

    while (!queue.empty()) {
        CombinedCartwheel current = std::move(queue.front());
        queue.pop();
        std::ostringstream state;
        state << cartwheel_key(current.cartwheel) << '|';
        for (const Rule& r : current.fixed_incoming) {
            for (bool b : r.members) state << (b ? '1' : '0');
            state << ';';
        }
        if (!processed.insert(state.str()).second) continue;

        bool refined = false;
        for (int i = 0; i < static_cast<int>(current.cartwheel.spokes.size()) && !refined; ++i) {
            for (const AuxiliaryCover& cover : auxiliary_covers) {
                if (!should_refine(current.cartwheel, i, cover)) continue;
                refined = true;
                for (Cartwheel updated : refinement(current.cartwheel, i, cover)) {
                    if (prune_cartwheel(updated, current.fixed_incoming, rules,
                                        combined_rules, configurations)) {
                        continue;
                    }
                    queue.push({std::move(updated), current.fixed_incoming});
                }
                break;
            }
        }
        if (!refined) result.push_back(std::move(current));
    }
    return dedupe_combined_cartwheels(std::move(result));
}

std::vector<Cartwheel> enum_possible_bad_wheels(
    int center_degree_value, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations) {
    std::vector<Cartwheel> result;
    auto process = [&](Cartwheel wheel) {
        if (!prune_cartwheel(wheel, {}, rules, combined_rules, configurations)) {
            result.push_back(std::move(wheel));
        }
    };
    // Stream the exponentially large wheel family instead of materializing it.  This is semantically
    // identical to Algorithm B.3.14 and is necessary for degrees 10 and 11.
    for_each_ordinary_wheel(center_degree_value, process);
    for_each_digon_wheel(center_degree_value, process);
    return dedupe_cartwheels(std::move(result));
}

bool verify_no_bad_cartwheels(
    const Cartwheel& wheel, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers) {
    auto fixed = fix_in_rules(wheel, rules, combined_rules, configurations);
    auto final = fix_out_rules(fixed, rules, combined_rules, configurations, auxiliary_covers);
    return final.empty();
}

}  // namespace apex

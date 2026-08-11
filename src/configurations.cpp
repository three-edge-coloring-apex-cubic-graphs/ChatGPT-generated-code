#include "apex/apex.hpp"

#include <atomic>
#include <cstdint>
#include <sstream>

namespace apex {
namespace {

struct ConfigurationIndex {
    int degree_max = 8;
    std::map<std::pair<int, int>, std::vector<std::size_t>> by_endpoint_degrees;
};

std::string configuration_set_token(const std::vector<RootedConfiguration>& configurations) {
    std::ostringstream ss;
    ss << reinterpret_cast<std::uintptr_t>(&configurations) << ':' << configurations.size();
    if (!configurations.empty()) {
        ss << ':' << configurations.front().name << ':' << configurations.back().name;
    }
    return ss.str();
}

const ConfigurationIndex& configuration_index(
    const std::vector<RootedConfiguration>& configurations) {
    thread_local std::unordered_map<std::string, ConfigurationIndex> cache;
    const std::string token = configuration_set_token(configurations);
    if (auto it = cache.find(token); it != cache.end()) return it->second;
    ConfigurationIndex index;
    for (std::size_t i = 0; i < configurations.size(); ++i) {
        const auto& configuration = configurations[i];
        for (VertexId v : configuration.graph.vertices()) {
            const DegreeRange r = configuration.graph.degree_range[v];
            if (r.fixed() && !r.unbounded()) index.degree_max = std::max(index.degree_max, r.lower);
        }
        const DartId f = configuration.special;
        if (f == kNil) continue;
        const VertexId y = configuration.graph.darts[f].head;
        const VertexId x = configuration.graph.tail(f);
        const DegreeRange ry = configuration.graph.degree_range[y];
        const DegreeRange rx = configuration.graph.degree_range[x];
        if (ry.fixed() && rx.fixed() && !ry.unbounded() && !rx.unbounded()) {
            index.by_endpoint_degrees[{ry.lower, rx.lower}].push_back(i);
        }
    }
    if (cache.size() > 1024) cache.clear();
    return cache.emplace(token, std::move(index)).first->second;
}

bool contains_configuration(const Embedding& target,
                            std::optional<VertexId> center,
                            const std::vector<int>& concrete_degree,
                            const std::vector<RootedConfiguration>& configurations,
                            const ConfigurationIndex& index) {
    std::vector<std::vector<std::vector<DartId>>> by_degree(
        index.degree_max + 1, std::vector<std::vector<DartId>>(index.degree_max + 1));
    for (DartId e : target.dart_ids()) {
        const int dy = concrete_degree[target.darts[e].head];
        const int dx = concrete_degree[target.tail(e)];
        if (0 <= dy && dy <= index.degree_max && 0 <= dx && dx <= index.degree_max) {
            by_degree[dy][dx].push_back(e);
        }
    }

    Embedding concrete = target;
    for (VertexId v : concrete.vertices()) {
        concrete.degree_range[v] = {concrete_degree[v], concrete_degree[v]};
    }

    for (const auto& [degrees, configuration_indices] : index.by_endpoint_degrees) {
        const auto [dy, dx] = degrees;
        if (dy > index.degree_max || dx > index.degree_max || by_degree[dy][dx].empty()) continue;
        for (std::size_t configuration_index_value : configuration_indices) {
            const auto& configuration = configurations[configuration_index_value];
            const DartId f = configuration.special;
            for (DartId fs : by_degree[dy][dx]) {
                if (center.has_value() && dy > 8 &&
                    concrete.darts[fs].head != *center) {
                    continue;
                }
                if (homomorphism(configuration.graph, f, concrete, fs,
                                 DegreeConstraint::Include)
                        .has_value()) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool blocked_by_reducible_configuration_impl(
    const Embedding& z, std::optional<VertexId> center,
    const std::vector<RootedConfiguration>& configurations) {
    if (configurations.empty()) return false;
    if (center.has_value() &&
        (*center < 0 || static_cast<std::size_t>(*center) >= z.vertex_alive.size() ||
         !z.vertex_alive[*center])) {
        throw std::invalid_argument("configuration-blocking center is not a live vertex");
    }

    thread_local std::unordered_map<std::string, bool> cache;
    const std::string cache_key =
        configuration_set_token(configurations) + "|" +
        (center.has_value() ? "centered|" : "uncentered|") +
        canonical_key(z, std::nullopt, {}, center);
    if (auto it = cache.find(cache_key); it != cache.end()) return it->second;

    const ConfigurationIndex& index_data = configuration_index(configurations);
    const auto vertices = z.vertices();
    std::vector<int> concrete(z.vertex_alive.size(), 0);

    std::function<bool(std::size_t)> dfs = [&](std::size_t position) -> bool {
        if (position == vertices.size()) {
            return contains_configuration(z, center, concrete, configurations,
                                          index_data);
        }
        const VertexId v = vertices[position];
        const DegreeRange r = z.degree_range[v];
        std::vector<int> choices;

        if (!center.has_value()) {
            // In Appendix B.4 there is no distinguished center.  Therefore a
            // configuration's (possible) degree-above-8 vertex may map to any
            // target vertex.  Every target vertex consequently uses the
            // CONF_DEG_MAX cutoff, rather than the non-center cutoff 8 from
            // Algorithm A.7.2.
            if (r.upper > index_data.degree_max) {
                choices.push_back(r.upper);
            } else {
                if (r.unbounded()) {
                    throw std::runtime_error(
                        "unbounded uncentered representative range was not reduced");
                }
                for (int d = r.lower; d <= r.upper; ++d) choices.push_back(d);
            }
        } else if (v == *center && r.upper > index_data.degree_max) {
            // Algorithm A.7.2 uses the upper endpoint itself as the
            // representative.  In particular, an unbounded range is
            // represented by kInfinity, not by an arbitrary finite degree just
            // above CONF_DEG_MAX.
            choices.push_back(r.upper);
        } else if (v != *center && r.upper > 8) {
            choices.push_back(r.upper);
        } else {
            if (r.unbounded()) {
                throw std::runtime_error(
                    "unbounded centered representative range was not reduced");
            }
            for (int d = r.lower; d <= r.upper; ++d) choices.push_back(d);
        }

        for (int d : choices) {
            concrete[v] = d;
            if (!dfs(position + 1)) return false;
        }
        return true;
    };

    const bool answer = dfs(0);
    if (cache.size() > 200000) cache.clear();
    cache.emplace(cache_key, answer);
    return answer;
}

std::string application_cache_key(const Embedding& z, DartId e, const Rule& rule) {
    return canonical_key(z, e) + "|RULE|" + canonical_key(rule.graph, rule.distinguished);
}

std::string rule_key(const Rule& rule) {
    std::ostringstream ss;
    ss << canonical_key(rule.graph, rule.distinguished) << "|r=" << rule.charge << "|m=";
    for (bool b : rule.members) ss << (b ? '1' : '0');
    return ss.str();
}

std::vector<Rule> dedupe_rules(std::vector<Rule> rules) {
    std::unordered_set<std::string> seen;
    std::vector<Rule> result;
    result.reserve(rules.size());
    for (auto& rule : rules) {
        const std::string key = rule_key(rule);
        if (seen.insert(key).second) result.push_back(std::move(rule));
    }
    return result;
}

}  // namespace

DartId maximum_degree_dart(const Embedding& z) {
    DartId best = kNil;
    std::pair<int, int> best_degree{0, 0};
    for (DartId e : z.dart_ids()) {
        const VertexId y = z.darts[e].head;
        const VertexId x = z.tail(e);
        if (!z.degree_range[y].fixed() || !z.degree_range[x].fixed()) continue;
        const std::pair<int, int> degree{z.degree_range[y].lower, z.degree_range[x].lower};
        if (degree > best_degree) {
            best_degree = degree;
            best = e;
        }
    }
    return best;
}

bool blocked_by_reducible_configuration(
    const Embedding& z,
    const std::vector<RootedConfiguration>& configurations) {
    return blocked_by_reducible_configuration_impl(z, std::nullopt,
                                                    configurations);
}

bool blocked_by_reducible_configuration_at_center(
    const Embedding& z, VertexId center,
    const std::vector<RootedConfiguration>& configurations) {
    return blocked_by_reducible_configuration_impl(z, center, configurations);
}

bool always_apply(const Embedding& z, DartId e, const Rule& rule) {
    const VertexId zh = z.darts[e].head;
    const VertexId zt = z.tail(e);
    const VertexId rh = rule.graph.darts[rule.distinguished].head;
    const VertexId rt = rule.graph.tail(rule.distinguished);
    if (!rule.graph.degree_range[rh].includes(z.degree_range[zh]) ||
        !rule.graph.degree_range[rt].includes(z.degree_range[zt])) {
        return false;
    }
    // Initial wheels are small and each complete rooted degree assignment is
    // normally queried only once.  Avoid constructing and retaining a large
    // canonical-string cache for these one-shot calls.
    if (z.darts.size() <= 44) {
        return homomorphism(rule.graph, rule.distinguished, z, e,
                            DegreeConstraint::Include)
            .has_value();
    }
    thread_local std::unordered_map<std::string, bool> cache;
    const std::string key = "A|" + application_cache_key(z, e, rule);
    if (auto it = cache.find(key); it != cache.end()) return it->second;
    const bool answer =
        homomorphism(rule.graph, rule.distinguished, z, e, DegreeConstraint::Include).has_value();
    if (cache.size() > 200000) cache.clear();
    cache.emplace(key, answer);
    return answer;
}

bool never_apply(const Embedding& z, DartId e, const Rule& rule) {
    const VertexId zh = z.darts[e].head;
    const VertexId zt = z.tail(e);
    const VertexId rh = rule.graph.darts[rule.distinguished].head;
    const VertexId rt = rule.graph.tail(rule.distinguished);
    if (!z.degree_range[zh].intersects(rule.graph.degree_range[rh]) ||
        !z.degree_range[zt].intersects(rule.graph.degree_range[rt])) {
        return true;
    }
    auto evaluate = [&]() {
        DisjointUnionResult joined = disjoint_union(
            z, rule.graph, EmbeddingKind::PseudoTriangulationWithDigons);
        const DartId ze = joined.left.dart[e];
        const DartId re = joined.right.dart[rule.distinguished];
        return free_homomorphism_and_enforce_single_digon_incidence(
                   joined.graph, {{ze, re}})
            .empty();
    };
    if (z.darts.size() <= 44) return evaluate();

    thread_local std::unordered_map<std::string, bool> cache;
    const std::string key = "N|" + application_cache_key(z, e, rule);
    if (auto it = cache.find(key); it != cache.end()) return it->second;
    const bool answer = evaluate();
    if (cache.size() > 200000) cache.clear();
    cache.emplace(key, answer);
    return answer;
}

std::vector<Rule> combine_rules(
    const std::vector<Rule>& rules, const std::vector<RootedConfiguration>& configurations) {
    const int n = static_cast<int>(rules.size());
    std::vector<std::vector<int>> rotations{{1, -1}, {0, -1}};
    Embedding initial_graph = from_vertex_rotations(
        2, rotations, {}, EmbeddingKind::PseudoTriangulationWithDigons);
    initial_graph.degree_range[0] = {1, kInfinity};
    initial_graph.degree_range[1] = {1, kInfinity};
    DartId initial_dart = kNil;
    for (DartId e : initial_graph.dart_ids()) {
        if (initial_graph.tail(e) == 0 && initial_graph.darts[e].head == 1) initial_dart = e;
    }
    Rule initial{"empty", std::move(initial_graph), initial_dart, 0,
                 std::vector<bool>(n, false)};
    std::vector<Rule> combined{std::move(initial)};

    for (int rule_index = 0; rule_index < n; ++rule_index) {
        std::vector<Rule> next = combined;
        const Rule& added = rules[rule_index];
        for (const Rule& base : combined) {
            DisjointUnionResult joined = disjoint_union(
                base.graph, added.graph, EmbeddingKind::PseudoTriangulationWithDigons);
            const DartId ebase = joined.left.dart[base.distinguished];
            const DartId eadded = joined.right.dart[added.distinguished];
            for (auto image : free_homomorphism_and_enforce_single_digon_incidence(
                     joined.graph, {{ebase, eadded}})) {
                Rule new_rule;
                new_rule.name = base.name + "+" + added.name;
                new_rule.graph = std::move(image.target);
                new_rule.distinguished = image.map.dart[ebase];
                new_rule.charge = base.charge + added.charge;
                new_rule.members = base.members;
                if (new_rule.members.size() != static_cast<std::size_t>(n)) {
                    new_rule.members.resize(n, false);
                }
                new_rule.members[rule_index] = true;
                if (!configurations.empty() &&
                    blocked_by_reducible_configuration_at_center(
                        new_rule.graph, new_rule.graph.darts[new_rule.distinguished].head,
                        configurations)) {
                    continue;
                }
                next.push_back(std::move(new_rule));
            }
        }
        combined = dedupe_rules(std::move(next));
    }
    return combined;
}

}  // namespace apex

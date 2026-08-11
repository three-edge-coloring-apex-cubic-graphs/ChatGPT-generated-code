#include "apex/apex.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace apex {

VerificationMetrics verify_lemma_b1(
    const std::vector<Rule>& rules,
    const std::vector<RootedConfiguration>& configurations,
    std::vector<Rule>* combined_output) {
    VerificationMetrics metrics;
    std::vector<Rule> combined = combine_rules(rules, configurations);
    metrics.combined_rule_count = combined.size();
    metrics.maximum_combined_charge = std::numeric_limits<int>::min();
    for (const Rule& rule : combined) {
        metrics.maximum_combined_charge =
            std::max(metrics.maximum_combined_charge, rule.charge);
    }
    if (combined_output != nullptr) *combined_output = std::move(combined);
    return metrics;
}

VerificationMetrics verify_lemma_b2(
    const std::vector<Rule>& rules, const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers,
    int first_degree, int last_degree) {
    VerificationMetrics metrics;
    for (int d = first_degree; d <= last_degree; ++d) {
        std::cerr << "[B2] center degree " << d
                  << ": enumerating possible bad wheels...\n";
        std::vector<Cartwheel> possible =
            enum_possible_bad_wheels(d, rules, combined_rules, configurations);
        metrics.possible_bad_wheels[d] = possible.size();
        std::cerr << "[B2] center degree " << d << ": " << possible.size()
                  << " possible bad wheels; verifying refinements...\n";
        for (std::size_t i = 0; i < possible.size(); ++i) {
            if (i != 0 && i % 100 == 0) {
                std::cerr << "[B2] center degree " << d << ": verified " << i
                          << '/' << possible.size() << " wheels\n";
            }
            if (!verify_no_bad_cartwheels(possible[i], rules, combined_rules,
                                          configurations, auxiliary_covers)) {
                throw std::runtime_error("Lemma B.2 assertion failed for center degree " +
                                         std::to_string(d) + ", wheel index " +
                                         std::to_string(i));
            }
        }
        std::cerr << "[B2] center degree " << d << ": verification complete\n";
    }
    return metrics;
}

VerificationMetrics verify_lemma_b3(
    const std::vector<ConfigurationFile>& configurations,
    std::vector<Island>* islands_output,
    bool check_reducibility,
    std::size_t first_configuration,
    std::size_t end_configuration,
    const B3SearchOptions& options) {
    VerificationMetrics metrics;
    std::vector<RootedConfiguration> smaller;
    std::vector<Island> all_islands;
    std::unordered_map<std::string, ReducibilityResult> reducibility_cache;
    std::size_t reducibility_cache_hits = 0;
    std::size_t generated_occurrences = 0;

    // Keep the four reported Lemma B.3 categories present even when a shard
    // happens to generate no island in one of them.  Counts for four or more
    // rings are accumulated separately by print_metrics.
    for (int rings = 0; rings <= 3; ++rings) {
        metrics.generated_island_occurrences_by_ring_count[rings] = 0;
    }

    const std::size_t end = std::min(end_configuration, configurations.size());
    if (first_configuration > end) {
        throw std::invalid_argument("invalid B.3 configuration range");
    }
    // Process the prefix even when output is sharded: K_smaller must contain
    // every configuration whose filename precedes the first requested one.
    for (std::size_t configuration_index = 0; configuration_index < end;
         ++configuration_index) {
        const ConfigurationFile& configuration = configurations[configuration_index];
        if (configuration_index >= first_configuration) {
            if ((configuration_index - first_configuration) % 10 == 0) {
                std::cerr << "[B3] processing configuration " << (configuration_index + 1)
                          << '/' << configurations.size() << " (" << configuration.name
                          << "), generated island occurrences in this range so far: "
                          << generated_occurrences << "\n";
            }
            Embedding outer = outer_extension_from_configuration(configuration);
            auto islands = all_hom_images(outer, smaller, options);
            for (Island& island : islands) {
                // Check every occurrence returned by Algorithm B.4.1.  A cached
                // verdict only reuses the computation for an identical island;
                // the occurrence is still checked, counted, and optionally written.
                if (check_reducibility) {
                    ReducibilityResult result;
                    if (options.cache_reducibility_results) {
                        const std::string key = island.canonical_key();
                        if (auto it = reducibility_cache.find(key);
                            it != reducibility_cache.end()) {
                            result = it->second;
                            ++reducibility_cache_hits;
                        } else {
                            result = check_semi_reducibility(island, true);
                            reducibility_cache.emplace(key, result);
                        }
                    } else {
                        result = check_semi_reducibility(island, true);
                    }
                    if (!result.semi_d_reducible && !result.semi_c_reducible) {
                        throw std::runtime_error(
                            "Lemma B.3 found a non-semi-reducible island while processing " +
                            configuration.name);
                    }
                }
                ++metrics.generated_island_occurrences_by_ring_count[
                    static_cast<int>(island.ring_sizes.size())];
                ++generated_occurrences;
                if (islands_output != nullptr) all_islands.push_back(std::move(island));
            }
        }
        auto variants = extend_from_cut_vertices(configuration);
        smaller.insert(smaller.end(), std::make_move_iterator(variants.begin()),
                       std::make_move_iterator(variants.end()));
    }
    if (check_reducibility && options.cache_reducibility_results) {
        std::cerr << "[B3] reducibility cache: unique_islands="
                  << reducibility_cache.size()
                  << ", replayed_occurrences=" << reducibility_cache_hits << "\n";
    }
    if (islands_output != nullptr) *islands_output = std::move(all_islands);
    return metrics;
}

}  // namespace apex

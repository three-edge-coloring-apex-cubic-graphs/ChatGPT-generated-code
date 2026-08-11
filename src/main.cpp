#include "apex/apex.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;
using namespace apex;

namespace {

struct Options {
    std::map<std::string, std::string> values;
    std::set<std::string> flags;
};

Options parse_options(int argc, char** argv, int start) {
    Options options;
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];
        if (!arg.starts_with("--")) throw std::runtime_error("unexpected argument: " + arg);
        if (i + 1 < argc && !std::string(argv[i + 1]).starts_with("--")) {
            options.values[arg] = argv[++i];
        } else {
            options.flags.insert(arg);
        }
    }
    return options;
}

fs::path input_path(const Options& options, const std::string& option,
                    const fs::path& data_relative) {
    if (auto it = options.values.find(option); it != options.values.end()) {
        return it->second;
    }
    if (auto it = options.values.find("--data-root"); it != options.values.end()) {
        return fs::path(it->second) / data_relative;
    }
    throw std::runtime_error("missing required option " + option +
                             " (or supply --data-root)");
}

int integer_option(const Options& options, const std::string& name, int fallback) {
    auto it = options.values.find(name);
    if (it == options.values.end()) return fallback;
    std::size_t used = 0;
    const int value = std::stoi(it->second, &used);
    if (used != it->second.size()) throw std::runtime_error("bad integer for " + name);
    return value;
}

std::vector<AuxiliaryCover> load_auxiliary_covers_from_manifest(
    const fs::path& manifest, const fs::path& directory) {
    std::ifstream input(manifest);
    if (!input) throw std::runtime_error("cannot open auxiliary-cover manifest");
    std::vector<AuxiliaryCover> result;
    std::string line;
    while (std::getline(input, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        std::istringstream stream(line);
        std::vector<std::string> names;
        std::string name;
        while (stream >> name) names.push_back(name);
        if (names.empty()) continue;
        if (names.size() < 2) {
            throw std::runtime_error("an auxiliary cover needs a base and a cover");
        }
        AuxiliaryCover cover;
        cover.base = read_rule_file(directory / names.front());
        for (std::size_t i = 1; i < names.size(); ++i) {
            cover.cover.push_back(read_rule_file(directory / names[i]));
        }
        result.push_back(std::move(cover));
    }
    return result;
}

std::vector<AuxiliaryCover> auxiliary_covers(const Options& options) {
    if (auto it = options.values.find("--auxiliary"); it != options.values.end()) {
        return apex::load_auxiliary_covers(it->second);
    }
    if (auto it = options.values.find("--data-root"); it != options.values.end()) {
        return apex::load_auxiliary_covers(fs::path(it->second) /
                                           "discharging-rules/R_auxiliary");
    }
    auto manifest = options.values.find("--aux-manifest");
    auto directory = options.values.find("--aux-dir");
    if (manifest != options.values.end() && directory != options.values.end()) {
        return load_auxiliary_covers_from_manifest(manifest->second, directory->second);
    }
    throw std::runtime_error(
        "missing --auxiliary DIR (or --aux-manifest FILE with --aux-dir DIR)");
}

void print_metrics(const VerificationMetrics& metrics) {
    if (metrics.combined_rule_count != 0) {
        std::cout << "combined_rule_count=" << metrics.combined_rule_count << '\n';
        std::cout << "maximum_combined_charge=" << metrics.maximum_combined_charge << '\n';
    }
    for (auto [degree, count] : metrics.possible_bad_wheels) {
        std::cout << "possible_bad_wheels_degree_" << degree << '=' << count << '\n';
    }
    if (!metrics.generated_island_occurrences_by_ring_count.empty()) {
        for (int rings = 0; rings <= 3; ++rings) {
            const auto it =
                metrics.generated_island_occurrences_by_ring_count.find(rings);
            const std::size_t count =
                it == metrics.generated_island_occurrences_by_ring_count.end()
                    ? 0
                    : it->second;
            std::cout << "generated_island_occurrences_with_" << rings
                      << "_rings=" << count << '\n';
        }
        std::size_t at_least_four = 0;
        for (const auto& [rings, count] :
             metrics.generated_island_occurrences_by_ring_count) {
            if (rings >= 4) at_least_four += count;
        }
        std::cout << "generated_island_occurrences_with_at_least_4_rings="
                  << at_least_four << '\n';
    }
    std::cout.flush();
}

void usage() {
    std::cerr
        << "Usage:\n"
        << "  apex_verify b1 (--data-root DIR | --rules DIR --configurations DIR) [--combined-out DIR]\n"
        << "  apex_verify b2 (--data-root DIR | --rules DIR --configurations DIR --auxiliary DIR) [--first-degree 7 --last-degree 11] [--enumeration-only]\n"
        << "  apex_verify b3 (--data-root DIR | --configurations DIR) [--first-configuration I --last-configuration J] [--island-occurrences FILE] [--islands-out DIR] [--skip-reducibility] [--literal-search] [--no-reducibility-cache]\n"
        << "  apex_verify validate-data (--data-root DIR | --rules DIR --configurations DIR --auxiliary DIR)\n"
        << "  apex_verify all (--data-root DIR | --rules DIR --configurations DIR --auxiliary DIR)\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            usage();
            return 2;
        }
        const std::string command = argv[1];
        const Options options = parse_options(argc, argv, 2);

        if (command == "b1") {
            auto rules = load_rules(input_path(options, "--rules", "discharging-rules/R"));
            auto configurations = load_configurations(
                input_path(options, "--configurations", "configurations/K"));
            auto rooted = build_rooted_configuration_set(configurations);
            std::vector<Rule> combined;
            VerificationMetrics metrics = verify_lemma_b1(rules, rooted, &combined);
            print_metrics(metrics);
            if (auto it = options.values.find("--combined-out"); it != options.values.end()) {
                fs::create_directories(it->second);
                for (std::size_t i = 0; i < combined.size(); ++i) {
                    std::ostringstream filename;
                    filename << "combined_" << std::setw(5) << std::setfill('0') << i
                             << ".rule";
                    write_rule_file(combined[i], fs::path(it->second) / filename.str(), true);
                }
            }
            return 0;
        }

        if (command == "b2") {
            auto rules = load_rules(input_path(options, "--rules", "discharging-rules/R"));
            auto configurations = load_configurations(
                input_path(options, "--configurations", "configurations/K"));
            auto rooted = build_rooted_configuration_set(configurations);
            std::vector<Rule> combined;
            VerificationMetrics b1 = verify_lemma_b1(rules, rooted, &combined);
            print_metrics(b1);

            const int first_degree = integer_option(options, "--first-degree", 7);
            const int last_degree = integer_option(options, "--last-degree", 11);
            if (first_degree < 7 || last_degree > 11 || first_degree > last_degree) {
                throw std::runtime_error("invalid center-degree range");
            }
            if (options.flags.contains("--enumeration-only")) {
                VerificationMetrics enumeration;
                for (int degree = first_degree; degree <= last_degree; ++degree) {
                    std::cerr << "[B2] center degree " << degree
                              << ": enumerating possible bad wheels only...\n";
                    const auto possible = enum_possible_bad_wheels(
                        degree, rules, combined, rooted);
                    enumeration.possible_bad_wheels[degree] = possible.size();
                }
                print_metrics(enumeration);
                return 0;
            }

            auto auxiliary = auxiliary_covers(options);
            VerificationMetrics b2 = verify_lemma_b2(
                rules, combined, rooted, auxiliary, first_degree, last_degree);
            print_metrics(b2);
            return 0;
        }

        if (command == "b3") {
            auto configurations = load_configurations(
                input_path(options, "--configurations", "configurations/K"));
            const int first_one_based = integer_option(options, "--first-configuration", 1);
            const int last_one_based = integer_option(
                options, "--last-configuration", static_cast<int>(configurations.size()));
            if (first_one_based < 1 || last_one_based < first_one_based ||
                last_one_based > static_cast<int>(configurations.size())) {
                throw std::runtime_error("invalid configuration range");
            }

            auto occurrence_file = options.values.find("--island-occurrences");
            // Retain the old spelling as a compatibility alias.  Records are not
            // canonicalized and duplicate occurrences are intentionally preserved.
            if (occurrence_file == options.values.end()) {
                occurrence_file = options.values.find("--island-keys");
            }
            const auto islands_directory = options.values.find("--islands-out");
            const bool retain_islands = occurrence_file != options.values.end() ||
                                        islands_directory != options.values.end();
            std::vector<Island> islands;

            B3SearchOptions search_options;
            if (options.flags.contains("--literal-search")) {
                search_options.prune_impossible_identifications = false;
                search_options.memoize_recursive_states = false;
                search_options.memoize_outer_extensions = false;
                search_options.memoize_equivalent_pair_branches = false;
            }
            if (options.flags.contains("--no-reducibility-cache")) {
                search_options.cache_reducibility_results = false;
            }

            VerificationMetrics metrics = verify_lemma_b3(
                configurations, retain_islands ? &islands : nullptr,
                !options.flags.contains("--skip-reducibility"),
                static_cast<std::size_t>(first_one_based - 1),
                static_cast<std::size_t>(last_one_based), search_options);
            print_metrics(metrics);

            if (occurrence_file != options.values.end()) {
                std::ofstream occurrences(occurrence_file->second);
                if (!occurrences) {
                    throw std::runtime_error("cannot write " + occurrence_file->second);
                }
                for (const Island& island : islands) {
                    occurrences << island.ring_sizes.size() << '\n';
                }
            }
            if (islands_directory != options.values.end()) {
                fs::create_directories(islands_directory->second);
                for (std::size_t i = 0; i < islands.size(); ++i) {
                    std::ostringstream filename;
                    filename << "island_" << std::setw(6) << std::setfill('0') << i
                             << ".island";
                    write_island_file(islands[i],
                                      fs::path(islands_directory->second) / filename.str());
                }
            }
            return 0;
        }

        if (command == "validate-data") {
            auto rules = load_rules(input_path(options, "--rules", "discharging-rules/R"));
            auto configurations = load_configurations(
                input_path(options, "--configurations", "configurations/K"));
            auto auxiliary = auxiliary_covers(options);
            std::size_t outer_darts = 0;
            for (const ConfigurationFile& configuration : configurations) {
                Embedding outer = outer_extension_from_configuration(configuration);
                std::string error;
                if (!outer.validate_single_list(&error)) {
                    throw std::runtime_error(
                        "invalid outer extension for " + configuration.name + ": " + error);
                }
                outer_darts += outer.dart_ids().size();
            }
            std::cout << "configuration_count=" << configurations.size() << '\n'
                      << "rule_count=" << rules.size() << '\n'
                      << "auxiliary_cover_count=" << auxiliary.size() << '\n'
                      << "outer_extension_dart_count=" << outer_darts << '\n';
            return 0;
        }

        if (command == "all") {
            auto rules = load_rules(input_path(options, "--rules", "discharging-rules/R"));
            auto configurations = load_configurations(
                input_path(options, "--configurations", "configurations/K"));
            auto rooted = build_rooted_configuration_set(configurations);
            std::vector<Rule> combined;
            print_metrics(verify_lemma_b1(rules, rooted, &combined));
            auto auxiliary = auxiliary_covers(options);
            print_metrics(verify_lemma_b2(rules, combined, rooted, auxiliary));
            print_metrics(verify_lemma_b3(configurations, nullptr, true));
            return 0;
        }

        usage();
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}

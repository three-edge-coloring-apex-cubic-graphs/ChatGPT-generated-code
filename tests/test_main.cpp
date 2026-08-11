#include "apex/apex.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace apex;
namespace fs = std::filesystem;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void require_throws_containing(Function&& action, const std::string& expected,
                               const char* message) {
    try {
        action();
    } catch (const std::exception& error) {
        if (std::string(error.what()).find(expected) != std::string::npos) return;
        throw std::runtime_error(std::string(message) + ": unexpected diagnostic: " +
                                 error.what());
    }
    throw std::runtime_error(std::string(message) + ": no exception was thrown");
}
}

int main() {
    try {
        {
            std::vector<std::vector<int>> rotations{{1, 2}, {2, 0}, {0, 1}};
            Embedding triangle = from_vertex_rotations(
                3, rotations, {}, EmbeddingKind::PseudoTriangulationWithDigons);
            for (VertexId v : triangle.vertices()) triangle.degree_range[v] = {2, 2};
            require(triangle.validate_faces(), "triangle validation failed");
            require(is_planar(triangle), "triangle should be planar");
            require(get_walks(triangle).size() == 2, "triangle should have two facial walks");
            // Algorithm B.4.5 emits the two orientations; this implementation does
            // not add a canonical cycle-key filter on top of the pseudocode.
            require(enum_cycles(triangle, 4).size() == 2,
                    "literal triangle cycle enumeration failed");

            // Pairs carried to each other by a rotation of the embedded triangle
            // must receive the same memoization signature; genuinely different
            // pair orbits must remain distinct.
            const auto signatures =
                canonical_unordered_dart_pair_signatures(triangle);
            const std::size_t stride = triangle.darts.size();
            require(signatures[0 * stride + 1] == signatures[2 * stride + 3] &&
                        signatures[2 * stride + 3] == signatures[4 * stride + 5],
                    "automorphic dart pairs received different signatures");
            require(signatures[0 * stride + 1] != signatures[0 * stride + 3],
                    "non-automorphic dart pairs were merged by the signature");
        }
        {
            // Algorithm B.4.10 returns Z_3 union Z_2 directly.  Use the same
            // feasible identification request for the two synthetic branches to
            // verify that the wrapper does not remove isomorphic/duplicate images.
            std::vector<std::vector<int>> rotations{{1, -1}, {2, 0, -1}, {1, -1}};
            Embedding path = from_vertex_rotations(3, rotations, {},
                                                    EmbeddingKind::PseudoEmbedding);
            for (VertexId v : path.vertices()) path.degree_range[v] = {1, 8};
            const std::array<DartId, 4> synthetic_darts{0, 1, 3, 3};
            const auto one_branch =
                free_homomorphism_and_enforce_single_digon_incidence(path, {{0, 3}});
            require(!one_branch.empty(), "synthetic B.4.10 branch should be feasible");
            const auto both_branches = ensure_outer_extension(path, synthetic_darts);
            require(both_branches.size() == 2 * one_branch.size(),
                    "B.4.10 branches were canonically deduplicated");
        }
        {
            std::vector<std::vector<int>> rotations{{1, -1}, {0, -1}};
            Embedding edge = from_vertex_rotations(2, rotations, {},
                                                    EmbeddingKind::PseudoEmbedding);
            edge.degree_range[0] = edge.degree_range[1] = {1, 1};
            auto completions = boundary_completions(edge, 0);
            require(completions.size() == 1, "pseudo-embedding boundary completion count");
            require(completions.front().target.is_inner(0), "boundary list was not closed");
        }
        {
            std::vector<int> ordered{0, 1, 2, 3};
            require(noncrossing_perfect_matchings(ordered).size() == 2,
                    "Catalan matching count for four points");
            require(!get_planar_half_kempes(3).empty(), "half Kempe enumeration failed");
            require(get_planar_half_kempes(9).size() == 210,
                    "cached half-Kempe template count for nine positions");
        }
        {
            std::vector<DegreeRange> degrees(4, {5, 5});
            Cartwheel wheel = generate_cartwheel(4, degrees, false);
            require(wheel.spokes.size() == 4, "wheel spoke count");
            require(enum_digons(wheel.graph).empty(), "ordinary wheel has a digon");
            Cartwheel digon = generate_cartwheel(5, degrees, true);
            require(digon.spokes.size() == 5, "digon wheel spoke count");
            require(enum_digons(digon.graph).size() == 1, "digon wheel missing digon");
        }
        {
            // Outer extension of the one-vertex, three-edge configuration.
            std::vector<std::vector<int>> rotations{{3, -1}, {3, -1}, {3, -1}, {0, 1, 2}};
            Embedding outer = from_vertex_rotations(4, rotations, {},
                                                    EmbeddingKind::PseudoEmbedding);
            outer.degree_range[0] = outer.degree_range[1] = outer.degree_range[2] = {5, kInfinity};
            outer.degree_range[3] = {3, 3};
            Embedding completion = free_completion_from_outer_extension(outer);
            require(completion.validate_faces(), "free completion is invalid");
            Island island = island_from_free_completion(completion);
            require(island.ring_sizes == std::vector<int>{3}, "wrong ring size in dual island");
            require(island.incident_edges.size() == 3, "wrong number of island vertices");
        }
        {
            // Algorithm A.7.2 represents an unbounded non-center range by its
            // actual upper endpoint (infinity), not by an arbitrary finite value.
            std::vector<std::vector<int>> rotations{{1, -1}, {2, 0, -1}, {1, -1}};
            Embedding source = from_vertex_rotations(3, rotations, {},
                                                      EmbeddingKind::PseudoEmbedding);
            source.degree_range[0] = {1, 1};
            source.degree_range[1] = {2, 2};
            source.degree_range[2] = {9, 9};
            DartId special = kNil;
            for (DartId e : source.dart_ids()) {
                if (source.tail(e) == 0 && source.darts[e].head == 1) special = e;
            }
            require(special != kNil, "failed to find synthetic special dart");
            RootedConfiguration configuration{"synthetic", source, special};

            Embedding target = source;
            target.degree_range[2] = {9, kInfinity};
            require(!blocked_by_reducible_configuration_at_center(
                        target, 0, {configuration}),
                    "unbounded representative degree was replaced by a finite degree");
        }
        {
            // Corrected Algorithm B.4.1: configuration blocking in allHomImages
            // is uncentered.  A high-degree vertex of the embedded
            // configuration may map to any target vertex; selecting vertex 0 as
            // an arbitrary center would incorrectly miss this identity image.
            std::vector<std::vector<int>> rotations{{1, -1}, {2, 0, -1}, {1, -1}};
            Embedding source = from_vertex_rotations(
                3, rotations, {}, EmbeddingKind::PseudoEmbedding);
            source.degree_range[0] = {5, 5};
            source.degree_range[1] = {6, 6};
            source.degree_range[2] = {9, 9};
            const DartId special = maximum_degree_dart(source);
            require(special != kNil && source.darts[special].head == 2,
                    "failed to root the high-degree synthetic configuration");
            RootedConfiguration configuration{"uncentered-synthetic", source,
                                                special};

            require(blocked_by_reducible_configuration(source, {configuration}),
                    "uncentered blocking missed a configuration away from vertex 0");
            require(!blocked_by_reducible_configuration_at_center(
                        source, 0, {configuration}),
                    "centered regression did not distinguish the wrong center");
            require(blocked_by_reducible_configuration_at_center(
                        source, 2, {configuration}),
                    "centered regression missed the correct high-degree center");

            // With no center, every vertex uses the configuration-set maximum
            // degree as the representative cutoff.  The target range [9,11]
            // must therefore test 9, 10, and 11; testing only its upper endpoint
            // would incorrectly report this degree-11 configuration as a blocker.
            Embedding ranged_target = source;
            ranged_target.degree_range[2] = {9, 11};
            Embedding degree_eleven_configuration = source;
            degree_eleven_configuration.degree_range[2] = {11, 11};
            const DartId degree_eleven_special =
                maximum_degree_dart(degree_eleven_configuration);
            require(!blocked_by_reducible_configuration(
                        ranged_target,
                        {RootedConfiguration{"degree-eleven-only",
                                             degree_eleven_configuration,
                                             degree_eleven_special}}),
                    "uncentered blocking collapsed a finite high-degree range");
        }
        {
            // A regression island generated from the supplied configurations.
            // The literal reducer and the optimized color-orbit/bit-mask reducer
            // both certify it as semi-C-reducible by deleting edge 5.
            Island island;
            island.ring_sizes = {4};
            island.degree_two_vertices = 1;
            island.incident_edges = {
                {5, 6, 0}, {6, 7, 8}, {7, 9, 10},
                {9, 11, 12}, {11, 5, 3}, {10, 13, 4},
                {13, 14, 15}, {14, 12, 2}, {8, 15, 1}};
            const ReducibilityResult result =
                check_semi_reducibility(island, true);
            require(!result.semi_d_reducible && result.semi_c_reducible,
                    "optimized semi-reducibility result changed");
            require(result.deleting_edges == std::vector<int>{5},
                    "optimized semi-C witness changed");
        }
        {
            Island island;
            island.ring_sizes = {};
            island.degree_two_vertices = 0;
            island.incident_edges = {{0, 1, 2}, {0, 2, 1}};
            const fs::path path = fs::temp_directory_path() / "apex_island_roundtrip.txt";
            write_island_file(island, path);
            Island read_back = read_island_file(path);
            fs::remove(path);
            require(read_back.ring_sizes.empty(), "zero-ring island parsing failed");
            require(read_back.incident_edges == island.incident_edges, "island roundtrip failed");
        }
        {
            // FORMAT.md lists the number of adjacent vertices, not delta_K.  The digon
            // 4--5 contributes one additional edge to the fixed degree of both endpoints.
            const fs::path path =
                fs::temp_directory_path() / "apex_configuration_digon_degree.conf";
            {
                std::ofstream output(path);
                output << "\n"
                       << "5 3\n"
                       << "4 4 1 2 3 5\n"
                       << "5 1 4\n"
                       << "1\n"
                       << "4 5\n";
            }
            ConfigurationFile configuration = read_configuration_file(path);
            fs::remove(path);
            require(configuration.prescribed_degree[3] == 5,
                    "configuration digon was not added to delta_K at vertex 4");
            require(configuration.prescribed_degree[4] == 2,
                    "configuration digon was not added to delta_K at vertex 5");

            Embedding outer = outer_extension_from_configuration(configuration);
            require(outer.degree(0) == 5, "digon dart degree does not match parsed delta_K");
            require(outer.degree_range[0] == DegreeRange{5, 5},
                    "outer extension did not retain the corrected fixed degree");
        }
        {
            // The supplied data contain a few legacy "0 0" sentinels before
            // the real header.  Also verify that a ring vertex with two
            // interior incidences is split into two distinct outer endpoints.
            const fs::path path =
                fs::temp_directory_path() / "apex_configuration_legacy_sentinel.conf";
            {
                std::ofstream output(path);
                output << "0 0\n"
                       << "4 2\n"
                       << "3 3 1 2 4\n"
                       << "4 3 1 2 3\n";
            }
            ConfigurationFile configuration = read_configuration_file(path);
            fs::remove(path);
            require(configuration.vertex_count == 4 && configuration.ring_size == 2,
                    "legacy configuration sentinel was not skipped");
            Embedding outer = outer_extension_from_configuration(configuration);
            int boundary_vertices = 0;
            for (VertexId v : outer.vertices()) {
                if (outer.is_boundary(v)) {
                    ++boundary_vertices;
                    require(outer.degree(v) == 1,
                            "outer-extension boundary endpoint is not degree one");
                }
            }
            require(boundary_vertices == 4,
                    "ring incidences were not split into distinct outer endpoints");
        }
        {
            // FORMAT.md: R, then k, then R_1,...,R_k.  Every rule record has
            // an explicit digon-count line, including the value 0.
            const fs::path path =
                fs::temp_directory_path() / "apex_auxiliary_cover.rule_auxiliary";
            const std::string rule =
                "2 1 2 2\n"
                "1 5 5 2 -1\n"
                "2 5 0 1 -1\n"
                "0\n";
            {
                std::ofstream output(path);
                output << "\n" << rule << "\n2\n\n" << rule << "\n" << rule
                       << "\n# trailing comments and blank lines are ignored\n\n";
            }
            AuxiliaryCover cover = read_auxiliary_cover_file(path);
            fs::remove(path);
            require(cover.cover.size() == 2,
                    "auxiliary-rule cover count was not parsed");
            require(cover.base.charge == 2 && cover.cover[0].charge == 2 &&
                        cover.cover[1].charge == 2,
                    "auxiliary-rule records are wrong");
            require(cover.base.name.ends_with("#R") &&
                        cover.cover[0].name.ends_with("#R_1") &&
                        cover.cover[1].name.ends_with("#R_2"),
                    "auxiliary-rule record names are wrong");
        }
        {
            const fs::path path =
                fs::temp_directory_path() / "apex_auxiliary_missing_m.rule_auxiliary";
            const std::string complete_rule =
                "2 1 2 2\n"
                "1 5 5 2 -1\n"
                "2 5 0 1 -1\n"
                "0\n";
            const std::string rule_without_m =
                "2 1 2 2\n"
                "1 5 5 2 -1\n"
                "2 5 0 1 -1\n";
            {
                std::ofstream output(path);
                output << complete_rule << "1\n" << rule_without_m;
            }
            require_throws_containing(
                [&] { static_cast<void>(read_auxiliary_cover_file(path)); },
                "explicitly contain the digon count M",
                "auxiliary covering rule without explicit M was accepted");
            fs::remove(path);
        }
        {
            const fs::path path =
                fs::temp_directory_path() / "apex_auxiliary_missing_base_m.rule_auxiliary";
            const std::string complete_rule =
                "2 1 2 2\n"
                "1 5 5 2 -1\n"
                "2 5 0 1 -1\n"
                "0\n";
            const std::string rule_without_m =
                "2 1 2 2\n"
                "1 5 5 2 -1\n"
                "2 5 0 1 -1\n";
            {
                std::ofstream output(path);
                output << rule_without_m << "1\n" << complete_rule;
            }
            require_throws_containing(
                [&] { static_cast<void>(read_auxiliary_cover_file(path)); },
                "explicit M line",
                "auxiliary base rule without explicit M was accepted");
            fs::remove(path);
        }
        {
            // Compare the default replay/memoization path with the exact literal
            // Appendix B.4 path on a nontrivial supplied configuration.  The
            // default mode may reuse an isomorphic representative internally,
            // but it must replay every occurrence.
            const fs::path configuration_directory =
                fs::path(APEX_SOURCE_DIR) / "data/configurations/K";
            const std::vector<ConfigurationFile> configurations =
                load_configurations(configuration_directory);
            require(configurations.size() == 915,
                    "wrong number of supplied configurations in B.4 regression");
            std::vector<RootedConfiguration> smaller;
            for (std::size_t i = 0; i < 3; ++i) {
                auto variants = extend_from_cut_vertices(configurations[i]);
                smaller.insert(smaller.end(),
                               std::make_move_iterator(variants.begin()),
                               std::make_move_iterator(variants.end()));
            }
            const Embedding outer =
                outer_extension_from_configuration(configurations[3]);

            B3SearchOptions literal;
            literal.prune_impossible_identifications = false;
            literal.memoize_recursive_states = false;
            literal.memoize_outer_extensions = false;
            literal.memoize_equivalent_pair_branches = false;
            literal.cache_reducibility_results = false;

            const std::vector<Island> literal_islands =
                all_hom_images(outer, smaller, literal);
            const std::vector<Island> accelerated_islands =
                all_hom_images(outer, smaller);
            require(accelerated_islands.size() == literal_islands.size(),
                    "B.4 acceleration changed the occurrence count");

            auto coarse_multiset = [](const std::vector<Island>& islands) {
                std::multiset<std::pair<std::vector<int>, int>> result;
                for (const Island& island : islands) {
                    result.emplace(island.ring_sizes,
                                   island.degree_two_vertices);
                }
                return result;
            };
            require(coarse_multiset(accelerated_islands) ==
                        coarse_multiset(literal_islands),
                    "B.4 acceleration changed the ring/degree-two occurrence multiset");
        }
        {
            const fs::path auxiliary_directory =
                fs::path(APEX_SOURCE_DIR) / "data/discharging-rules/R_auxiliary";
            const std::vector<AuxiliaryCover> covers =
                load_auxiliary_covers(auxiliary_directory);
            require(covers.size() == 4,
                    "wrong number of supplied auxiliary-rule files");
            for (const AuxiliaryCover& cover : covers) {
                require(cover.cover.size() == 3,
                        "a supplied auxiliary-rule file has the wrong cover count");
                require(cover.base.graph.validate_faces(),
                        "a supplied auxiliary base rule is invalid");
                for (const Rule& rule : cover.cover) {
                    require(rule.graph.validate_faces(),
                            "a supplied auxiliary covering rule is invalid");
                }
            }
        }
        std::cout << "all tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}

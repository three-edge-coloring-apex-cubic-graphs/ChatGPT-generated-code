#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iosfwd>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace apex {

using VertexId = int;
using DartId = int;
constexpr int kNil = -1;
constexpr int kInfinity = 1'000'000'000;

struct DegreeRange {
    int lower = 1;
    int upper = kInfinity;

    [[nodiscard]] bool empty() const noexcept { return lower > upper; }
    [[nodiscard]] bool fixed() const noexcept { return lower == upper; }
    [[nodiscard]] bool unbounded() const noexcept { return upper >= kInfinity; }
    [[nodiscard]] bool contains(int d) const noexcept { return lower <= d && d <= upper; }
    [[nodiscard]] bool includes(const DegreeRange& other) const noexcept {
        return lower <= other.lower && other.upper <= upper;
    }
    [[nodiscard]] bool intersects(const DegreeRange& other) const noexcept {
        return std::max(lower, other.lower) <= std::min(upper, other.upper);
    }
    [[nodiscard]] DegreeRange intersection(const DegreeRange& other) const noexcept {
        return {std::max(lower, other.lower), std::min(upper, other.upper)};
    }
    [[nodiscard]] std::string str() const;

    auto operator<=>(const DegreeRange&) const = default;
};

enum class EmbeddingKind {
    PseudoEmbedding,
    PseudoTriangulationWithDigons,
};

struct Dart {
    VertexId head = kNil;
    DartId rev = kNil;
    DartId succ = kNil;
    DartId pred = kNil;
    bool alive = true;
};

struct Morphism {
    std::vector<VertexId> vertex;
    std::vector<DartId> dart;

    [[nodiscard]] bool valid_vertex(VertexId v) const noexcept {
        return v >= 0 && static_cast<std::size_t>(v) < vertex.size() && vertex[v] != kNil;
    }
    [[nodiscard]] bool valid_dart(DartId e) const noexcept {
        return e >= 0 && static_cast<std::size_t>(e) < dart.size() && dart[e] != kNil;
    }
};

struct Embedding {
    EmbeddingKind kind = EmbeddingKind::PseudoEmbedding;
    std::vector<bool> vertex_alive;
    std::vector<DegreeRange> degree_range;
    std::vector<Dart> darts;

    [[nodiscard]] std::vector<VertexId> vertices() const;
    [[nodiscard]] std::vector<DartId> dart_ids() const;
    [[nodiscard]] std::vector<DartId> darts_at(VertexId v) const;
    [[nodiscard]] int degree(VertexId v) const;
    [[nodiscard]] VertexId tail(DartId e) const;
    [[nodiscard]] bool is_boundary(VertexId v) const;
    [[nodiscard]] bool is_inner(VertexId v) const { return !is_boundary(v); }
    [[nodiscard]] DartId first_dart(VertexId v) const;
    [[nodiscard]] DartId last_dart(VertexId v) const;
    [[nodiscard]] bool has_loop() const;
    [[nodiscard]] bool connected() const;
    [[nodiscard]] bool validate_basic(std::string* error = nullptr) const;
    [[nodiscard]] bool validate_single_list(std::string* error = nullptr) const;
    [[nodiscard]] bool validate_faces(std::string* error = nullptr) const;
    [[nodiscard]] std::string debug_string() const;

    VertexId add_vertex(DegreeRange range = {1, kInfinity});
    DartId add_dart(const Dart& dart);
};

struct EmbeddingImage {
    Embedding target;
    Morphism map;
};

struct DisjointUnionResult {
    Embedding graph;
    Morphism left;
    Morphism right;
};

struct CompactResult {
    Embedding graph;
    Morphism map;
};

struct Rule {
    std::string name;
    Embedding graph;
    DartId distinguished = kNil;
    int charge = 0;
    std::vector<bool> members;
};

struct RootedConfiguration {
    std::string name;
    Embedding graph;
    DartId special = kNil;
};

struct ConfigurationFile {
    std::string name;
    int vertex_count = 0;
    int ring_size = 0;
    // The fixed configuration degree delta_K.  In the file format this is computed as
    // the listed adjacent-vertex count plus the number of incident digons.
    std::vector<int> prescribed_degree;
    std::vector<std::vector<int>> rotations;
    std::vector<std::pair<int, int>> digons;
};

struct Cartwheel {
    Embedding graph;
    VertexId center = kNil;
    std::vector<DartId> spokes;
};

struct CombinedCartwheel {
    Cartwheel cartwheel;
    std::vector<Rule> fixed_incoming;
};

struct AuxiliaryCover {
    Rule base;
    std::vector<Rule> cover;
};

struct Island {
    std::vector<int> ring_sizes;
    int degree_two_vertices = 0;
    std::vector<std::array<int, 3>> incident_edges;

    [[nodiscard]] int ring_edge_count() const;
    [[nodiscard]] int edge_count() const;
    [[nodiscard]] std::string canonical_key() const;
};

struct ReducibilityResult {
    bool semi_d_reducible = false;
    bool semi_c_reducible = false;
    std::vector<int> deleting_edges;
};

struct VerificationMetrics {
    std::size_t combined_rule_count = 0;
    int maximum_combined_charge = std::numeric_limits<int>::min();
    std::map<int, std::size_t> possible_bad_wheels;
    std::map<int, std::size_t> generated_island_occurrences_by_ring_count;
};

// Optional, correctness-preserving accelerations for the Appendix B.4 / Lemma B.3
// computation.  Cached subproblems are replayed with their full multiplicity: no
// island occurrence and no feasible dart-pair branch is removed.
struct B3SearchOptions {
    bool prune_impossible_identifications = true;
    bool memoize_recursive_states = true;
    bool memoize_outer_extensions = true;
    bool memoize_equivalent_pair_branches = true;
    bool cache_reducibility_results = true;
};

// Generic helpers.
[[nodiscard]] Morphism identity_morphism(const Embedding& z);
[[nodiscard]] Morphism compose(const Morphism& second, const Morphism& first);
[[nodiscard]] CompactResult compact(const Embedding& z);
[[nodiscard]] DisjointUnionResult disjoint_union(const Embedding& a, const Embedding& b,
                                                 EmbeddingKind kind);
[[nodiscard]] Embedding mirror(const Embedding& z);
[[nodiscard]] std::string canonical_key(const Embedding& z,
                                        std::optional<DartId> rooted_dart = std::nullopt,
                                        std::span<const DartId> marked_darts = {},
                                        std::optional<VertexId> marked_vertex = std::nullopt);
[[nodiscard]] std::string canonical_image_key(const EmbeddingImage& image);
// For every pair (e,f), return an isomorphism-invariant signature of the
// unordered marked pair in z.  The result is flattened with stride
// z.darts.size(); dead-dart entries contain UINT64_MAX.
[[nodiscard]] std::vector<std::uint64_t> canonical_unordered_dart_pair_signatures(
    const Embedding& z);

// A.2--A.4 dependencies from arXiv:2603.24880.
enum class DegreeConstraint { Intersection, Include };
[[nodiscard]] std::optional<Morphism> homomorphism(const Embedding& source, DartId source_root,
                                                   const Embedding& target, DartId target_root,
                                                   DegreeConstraint constraint);
[[nodiscard]] std::optional<EmbeddingImage> dart_identification(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests);
[[nodiscard]] std::vector<EmbeddingImage> resolve_degree_issues(const Embedding& z);
[[nodiscard]] std::vector<EmbeddingImage> free_homomorphism(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests);

// Appendix B.2.
[[nodiscard]] std::vector<EmbeddingImage> boundary_completions(const Embedding& z, VertexId v);
[[nodiscard]] std::optional<EmbeddingImage> add_boundary_darts(const Embedding& z, VertexId v);
void add_boundary_darts_directly(Embedding& z, DartId first, DartId last);
[[nodiscard]] std::optional<EmbeddingImage> identify_neighbors(const Embedding& z, VertexId v);
[[nodiscard]] Morphism link_incidence_list_ends(Embedding& z, DartId u_first, DartId w_last);
[[nodiscard]] Embedding from_vertex_rotations(
    int n, const std::vector<std::vector<int>>& rotations,
    const std::vector<std::pair<int, int>>& digons = {},
    EmbeddingKind kind = EmbeddingKind::PseudoEmbedding);
[[nodiscard]] std::vector<std::tuple<VertexId, VertexId, VertexId>> find_cut_tuples(
    int n, int ring_size, const std::vector<std::vector<int>>& rotations);
[[nodiscard]] std::vector<RootedConfiguration> extend_from_cut_vertices(
    const ConfigurationFile& configuration);
[[nodiscard]] std::optional<std::pair<DartId, DartId>> two_digons_incident_with_same_vertex(
    const Embedding& z);
[[nodiscard]] std::vector<EmbeddingImage> enforce_single_digon_incidence(const Embedding& z);
[[nodiscard]] std::vector<EmbeddingImage> free_homomorphism_and_enforce_single_digon_incidence(
    const Embedding& z, const std::vector<std::pair<DartId, DartId>>& requests);

// Configuration/rule/cartwheel/island parsing and writing.
[[nodiscard]] ConfigurationFile read_configuration_file(const std::filesystem::path& path);
[[nodiscard]] Rule read_rule_file(const std::filesystem::path& path,
                                  std::optional<std::size_t> rule_count = std::nullopt);
[[nodiscard]] Rule read_combined_rule_file(const std::filesystem::path& path,
                                           std::size_t rule_count);
// Parse one auxiliary-rule file: R, k, R_1, ..., R_k.  Every embedded
// rule record must include an explicit digon count M, including M=0.
[[nodiscard]] AuxiliaryCover read_auxiliary_cover_file(
    const std::filesystem::path& path);
[[nodiscard]] Cartwheel read_cartwheel_file(const std::filesystem::path& path);
[[nodiscard]] Island read_island_file(const std::filesystem::path& path);
void write_rule_file(const Rule& rule, const std::filesystem::path& path, bool combined = false);
void write_cartwheel_file(const Cartwheel& cartwheel, const std::filesystem::path& path);
void write_island_file(const Island& island, const std::filesystem::path& path);
[[nodiscard]] std::vector<std::filesystem::path> sorted_regular_files(
    const std::filesystem::path& directory, std::string_view extension = {});
[[nodiscard]] std::vector<Rule> load_rules(const std::filesystem::path& directory);
[[nodiscard]] std::vector<AuxiliaryCover> load_auxiliary_covers(
    const std::filesystem::path& directory);
[[nodiscard]] std::vector<ConfigurationFile> load_configurations(
    const std::filesystem::path& directory);
[[nodiscard]] std::vector<RootedConfiguration> build_rooted_configuration_set(
    const std::vector<ConfigurationFile>& configurations);
[[nodiscard]] Embedding outer_extension_from_configuration(const ConfigurationFile& configuration);

// Blocking by configurations and free combinations (A.6--A.8 dependencies).
[[nodiscard]] DartId maximum_degree_dart(const Embedding& z);
// Uncentered blocking, used by Appendix B.4.  A reducible configuration may
// occur anywhere in z; no target vertex is distinguished.
[[nodiscard]] bool blocked_by_reducible_configuration(
    const Embedding& z, const std::vector<RootedConfiguration>& configurations);
// Centered blocking from Algorithms A.6--A.7 of the companion paper.  This is
// still required for combined rules and cartwheels, where a high-degree vertex
// of a reducible configuration must map to the designated center.
[[nodiscard]] bool blocked_by_reducible_configuration_at_center(
    const Embedding& z, VertexId center,
    const std::vector<RootedConfiguration>& configurations);
[[nodiscard]] bool always_apply(const Embedding& z, DartId e, const Rule& rule);
[[nodiscard]] bool never_apply(const Embedding& z, DartId e, const Rule& rule);
[[nodiscard]] std::vector<Rule> combine_rules(
    const std::vector<Rule>& rules, const std::vector<RootedConfiguration>& configurations);

// Appendix B.3.
inline constexpr std::array<DegreeRange, 5> kCartwheelDegrees{{
    {5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, kInfinity},
}};
[[nodiscard]] Cartwheel generate_cartwheel(int center_degree,
                                           const std::vector<DegreeRange>& degrees,
                                           bool incident_digon);
[[nodiscard]] std::vector<Cartwheel> enum_wheels(int center_degree);
[[nodiscard]] std::vector<Cartwheel> enum_digon_incident_wheels(int center_degree);
[[nodiscard]] std::set<std::pair<VertexId, VertexId>> enum_digons(const Embedding& z);
[[nodiscard]] int lower_bound_of_digon_charge(const Cartwheel& cartwheel);
[[nodiscard]] std::vector<Cartwheel> concrete_degree_except_tail(const Cartwheel& cartwheel);
[[nodiscard]] std::vector<Cartwheel> update_degree_by_rule(const Cartwheel& cartwheel,
                                                           DartId dart, const Rule& rule);
[[nodiscard]] int amount_of_charge_sent(const Cartwheel& cartwheel, DartId dart,
                                        const std::vector<Rule>& rules);
[[nodiscard]] int amount_of_possible_charge_sent(const Cartwheel& cartwheel, DartId dart,
                                                 const std::vector<Rule>& combined_rules);
[[nodiscard]] int upper_bound_of_charge(const Cartwheel& cartwheel,
                                        const std::vector<Rule>& fixed_incoming,
                                        const std::vector<Rule>& rules,
                                        const std::vector<Rule>& combined_rules);
[[nodiscard]] bool prune_by_non_associated_rule(const Cartwheel& cartwheel,
                                                const std::vector<Rule>& fixed_incoming,
                                                const std::vector<Rule>& rules);
[[nodiscard]] bool prune_cartwheel(const Cartwheel& cartwheel,
                                   const std::vector<Rule>& fixed_incoming,
                                   const std::vector<Rule>& rules,
                                   const std::vector<Rule>& combined_rules,
                                   const std::vector<RootedConfiguration>& configurations);
[[nodiscard]] std::vector<CombinedCartwheel> fix_in_rules(
    const Cartwheel& wheel, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations);
[[nodiscard]] bool should_refine(const Cartwheel& cartwheel, int spoke_index,
                                 const AuxiliaryCover& cover);
[[nodiscard]] std::vector<Cartwheel> refinement(const Cartwheel& cartwheel, int spoke_index,
                                                const AuxiliaryCover& cover);
[[nodiscard]] std::vector<CombinedCartwheel> fix_out_rules(
    const std::vector<CombinedCartwheel>& fixed, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers);
[[nodiscard]] std::vector<Cartwheel> enum_possible_bad_wheels(
    int center_degree, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations);
[[nodiscard]] bool verify_no_bad_cartwheels(
    const Cartwheel& wheel, const std::vector<Rule>& rules,
    const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers);

// Appendix B.4.
[[nodiscard]] std::vector<std::vector<DartId>> get_walks(const Embedding& z);
[[nodiscard]] bool is_planar(const Embedding& z);
[[nodiscard]] std::vector<std::vector<DartId>> enum_cycles(const Embedding& z, int max_length);
enum class SideLabel : std::uint8_t { Unknown, Left, Right };
[[nodiscard]] std::vector<SideLabel> label_darts(const Embedding& z,
                                                 const std::vector<DartId>& cycle);
[[nodiscard]] std::pair<int, int> num_separated_vertices(
    const Embedding& z, const std::vector<DartId>& cycle,
    const std::vector<SideLabel>& labels);
[[nodiscard]] bool has_separating_cycle(const Embedding& z);
[[nodiscard]] std::optional<std::array<DartId, 4>> find_four_darts(const Embedding& z);
[[nodiscard]] std::vector<EmbeddingImage> ensure_outer_extension(
    const Embedding& z, const std::array<DartId, 4>& darts);
[[nodiscard]] std::vector<EmbeddingImage> make_outer_extension(const Embedding& z);
[[nodiscard]] Embedding free_completion_from_outer_extension(const Embedding& outer_extension);
[[nodiscard]] Island island_from_free_completion(const Embedding& free_completion);
[[nodiscard]] std::vector<Island> all_hom_images(
    const Embedding& outer_extension,
    const std::vector<RootedConfiguration>& smaller_configurations,
    const B3SearchOptions& options = {});

// Appendix B.1 support used to validate Lemma B.3 outputs.
[[nodiscard]] std::vector<std::set<std::pair<int, int>>> noncrossing_perfect_matchings(
    const std::vector<int>& ordered_set);
[[nodiscard]] std::vector<std::set<std::pair<int, int>>> get_planar_half_kempes(int n);
[[nodiscard]] ReducibilityResult check_semi_reducibility(const Island& island,
                                                         bool search_c_reductions = true);

// High-level verification drivers. Full numerical verification requires the external rule/configuration data.
[[nodiscard]] VerificationMetrics verify_lemma_b1(
    const std::vector<Rule>& rules,
    const std::vector<RootedConfiguration>& configurations,
    std::vector<Rule>* combined_output = nullptr);
[[nodiscard]] VerificationMetrics verify_lemma_b2(
    const std::vector<Rule>& rules, const std::vector<Rule>& combined_rules,
    const std::vector<RootedConfiguration>& configurations,
    const std::vector<AuxiliaryCover>& auxiliary_covers,
    int first_degree = 7, int last_degree = 11);
[[nodiscard]] VerificationMetrics verify_lemma_b3(
    const std::vector<ConfigurationFile>& configurations,
    std::vector<Island>* islands_output = nullptr,
    bool check_reducibility = true,
    std::size_t first_configuration = 0,
    std::size_t end_configuration = std::numeric_limits<std::size_t>::max(),
    const B3SearchOptions& options = {});

}  // namespace apex

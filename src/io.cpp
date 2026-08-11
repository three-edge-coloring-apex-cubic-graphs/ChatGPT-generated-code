#include "apex/apex.hpp"

#include <fstream>
#include <sstream>

namespace apex {
namespace {

std::vector<std::string> read_nonempty_text_lines(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path.string());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.resize(comment);
        const auto begin = line.find_first_not_of(" \t\r\n");
        if (begin == std::string::npos) continue;
        const auto end = line.find_last_not_of(" \t\r\n");
        lines.push_back(line.substr(begin, end - begin + 1));
    }
    return lines;
}

int parse_upper(int value) { return value == 0 ? kInfinity : value; }
int print_upper(int value) { return value >= kInfinity ? 0 : value; }

std::vector<int> incidence_rotation(const Embedding& z, VertexId v) {
    std::vector<DartId> darts;
    if (z.is_boundary(v)) {
        for (DartId e = z.first_dart(v); e != kNil; e = z.darts[e].succ) darts.push_back(e);
    } else {
        const auto incident = z.darts_at(v);
        if (incident.empty()) return {};
        // Choose a start outside a consecutive equal-neighbor pair when possible.
        DartId start = incident.front();
        for (DartId e : incident) {
            const DartId p = z.darts[e].pred;
            if (p != kNil && z.tail(p) != z.tail(e)) {
                start = e;
                break;
            }
        }
        DartId e = start;
        do {
            darts.push_back(e);
            e = z.darts[e].succ;
        } while (e != start);
    }
    std::vector<int> rotation;
    for (DartId e : darts) {
        const int neighbor = z.tail(e);
        if (!rotation.empty() && rotation.back() == neighbor) continue;  // collapse a digon pair
        rotation.push_back(neighbor);
    }
    if (!rotation.empty() && !z.is_boundary(v) && rotation.size() > 1 &&
        rotation.front() == rotation.back()) {
        rotation.pop_back();
    }
    if (z.is_boundary(v)) rotation.push_back(-1);
    return rotation;
}

std::vector<std::pair<int, int>> digon_endpoints(const Embedding& z) {
    const auto digons = enum_digons(z);
    return {digons.begin(), digons.end()};
}

DartId directed_edge(const Embedding& z, int tail, int head) {
    for (DartId e : z.dart_ids()) {
        if (z.tail(e) == tail && z.darts[e].head == head) return e;
    }
    return kNil;
}

enum class DigonCountRequirement {
    OptionalAtEndOfFile,
    Required,
};

std::vector<int> parse_integer_line(const std::string& line,
                                    const std::string& context) {
    std::istringstream input(line);
    std::vector<int> values;
    int value = 0;
    while (input >> value) values.push_back(value);
    if (!input.eof()) {
        throw std::runtime_error("non-integer token in " + context);
    }
    return values;
}

Rule parse_rule_record(const std::vector<std::string>& text_lines,
                       std::size_t& pos, std::string name,
                       DigonCountRequirement digon_count_requirement) {
    const auto fail = [&name](const std::string& message) -> std::runtime_error {
        return std::runtime_error(message + " in " + name);
    };

    if (pos >= text_lines.size()) throw fail("missing rule record");
    const std::vector<int> header =
        parse_integer_line(text_lines[pos++], "rule header " + name);
    if (header.size() != 4 || header[0] <= 0) throw fail("bad rule header");

    const int n = header[0];
    const int source = header[1] - 1;
    const int target = header[2] - 1;
    const int charge = header[3];
    if (source < 0 || source >= n || target < 0 || target >= n) {
        throw fail("distinguished rule endpoint out of range");
    }

    std::vector<std::vector<int>> rotations(n);
    std::vector<DegreeRange> ranges(n);
    std::vector<bool> seen(n, false);
    for (int row = 0; row < n; ++row) {
        if (pos >= text_lines.size()) throw fail("truncated rule vertex list");
        const std::vector<int> values =
            parse_integer_line(text_lines[pos++], "rule vertex row " + name);
        if (values.size() < 3) throw fail("bad rule vertex row");

        const int index = values[0] - 1;
        const int lower = values[1];
        const int upper = parse_upper(values[2]);
        if (index < 0 || index >= n || seen[index]) {
            throw fail("duplicate or out-of-range rule vertex index");
        }
        if (lower < 0 || upper < lower) throw fail("invalid rule degree range");
        seen[index] = true;
        ranges[index] = {lower, upper};

        rotations[index].reserve(values.size() - 3);
        for (std::size_t i = 3; i < values.size(); ++i) {
            const int neighbor = values[i];
            if (neighbor != -1 && (neighbor < 1 || neighbor > n)) {
                throw fail("rule neighbor out of range");
            }
            rotations[index].push_back(neighbor == -1 ? -1 : neighbor - 1);
        }
    }

    int digon_count = 0;
    if (pos == text_lines.size()) {
        if (digon_count_requirement == DigonCountRequirement::Required) {
            throw fail("rule record must explicitly contain the digon count M, including M=0");
        }
    } else {
        const std::vector<int> count =
            parse_integer_line(text_lines[pos++], "rule digon count " + name);
        if (count.size() != 1 || count[0] < 0) throw fail("bad rule digon count");
        digon_count = count[0];
    }

    std::vector<std::pair<int, int>> digons;
    digons.reserve(static_cast<std::size_t>(digon_count));
    for (int i = 0; i < digon_count; ++i) {
        if (pos >= text_lines.size()) {
            throw fail("truncated rule digon list; auxiliary-rule records require an explicit M line");
        }
        const std::vector<int> endpoints =
            parse_integer_line(text_lines[pos++], "rule digon row " + name);
        if (endpoints.size() != 2 || endpoints[0] < 1 || endpoints[0] > n ||
            endpoints[1] < 1 || endpoints[1] > n || endpoints[0] == endpoints[1]) {
            throw fail("bad rule digon row; auxiliary-rule records require an explicit M line");
        }
        const std::pair<int, int> digon = std::minmax(endpoints[0] - 1, endpoints[1] - 1);
        if (std::find(digons.begin(), digons.end(), digon) != digons.end()) {
            throw fail("duplicate rule digon");
        }
        digons.push_back(digon);
    }

    Embedding graph;
    try {
        graph = from_vertex_rotations(
            n, rotations, digons, EmbeddingKind::PseudoTriangulationWithDigons);
    } catch (const std::exception& error) {
        throw fail(std::string("invalid rule rotation/digon data: ") + error.what());
    }
    graph.degree_range = ranges;
    const DartId distinguished = directed_edge(graph, source, target);
    if (distinguished == kNil) throw fail("distinguished rule edge is absent");
    return {std::move(name), std::move(graph), distinguished, charge, {}};
}

Rule read_rule_like(const std::filesystem::path& path, bool combined,
                    std::optional<std::size_t> rule_count) {
    const auto text_lines = read_nonempty_text_lines(path);
    if (text_lines.empty()) throw std::runtime_error("empty rule file " + path.string());

    std::size_t pos = 0;
    Rule result = parse_rule_record(
        text_lines, pos, path.filename().string(),
        combined ? DigonCountRequirement::Required
                 : DigonCountRequirement::OptionalAtEndOfFile);

    if (combined) {
        if (!rule_count.has_value()) throw std::invalid_argument("combined rule count is required");
        if (pos >= text_lines.size()) {
            throw std::runtime_error("combined rule lacks membership flag in " + path.string());
        }
        const std::string& flag = text_lines[pos++];
        if (flag.size() != *rule_count) {
            throw std::runtime_error("combined-rule flag length mismatch in " + path.string());
        }
        result.members.resize(flag.size());
        for (std::size_t i = 0; i < flag.size(); ++i) {
            if (flag[i] != '0' && flag[i] != '1') {
                throw std::runtime_error("bad combined-rule membership flag in " + path.string());
            }
            result.members[i] = flag[i] == '1';
        }
    }
    if (pos != text_lines.size()) {
        throw std::runtime_error("extra lines in rule file " + path.string());
    }
    return result;
}

}  // namespace

ConfigurationFile read_configuration_file(const std::filesystem::path& path) {
    const auto lines = read_nonempty_text_lines(path);
    if (lines.empty()) throw std::runtime_error("empty configuration file " + path.string());
    std::size_t pos = 0;
    ConfigurationFile result;
    result.name = path.filename().string();
    // Four files in the supplied data use a legacy empty-record sentinel "0 0"
    // before the actual configuration header.  Ignore any such leading records.
    for (;;) {
        if (pos >= lines.size()) {
            throw std::runtime_error("configuration contains only empty sentinels: " +
                                     path.string());
        }
        std::istringstream header(lines[pos++]);
        std::string extra;
        if (!(header >> result.vertex_count >> result.ring_size) || (header >> extra)) {
            throw std::runtime_error("bad configuration header in " + path.string());
        }
        if (result.vertex_count == 0 && result.ring_size == 0) continue;
        break;
    }
    const int n = result.vertex_count;
    const int r = result.ring_size;
    if (n < 0 || r < 0 || r > n) throw std::runtime_error("bad configuration sizes");
    // In the configuration format, the integer after the vertex index is the number of
    // adjacent vertices represented in the rotation.  It is not necessarily delta_K(v):
    // every incident digon contributes one additional edge to delta_K(v).
    std::vector<int> adjacent_vertex_count(n, 0);
    std::vector<bool> internal_vertex_seen(n, false);
    result.prescribed_degree.assign(n, 0);
    result.rotations.assign(n, {});
    for (int row = r; row < n; ++row) {
        if (pos >= lines.size()) throw std::runtime_error("truncated configuration");
        std::istringstream line(lines[pos++]);
        int index, adjacency_count;
        if (!(line >> index >> adjacency_count)) {
            throw std::runtime_error("bad configuration vertex line");
        }
        --index;
        if (index < r || index >= n) throw std::runtime_error("configuration index out of range");
        if (internal_vertex_seen[index]) {
            throw std::runtime_error("duplicate configuration vertex line for index " +
                                     std::to_string(index + 1));
        }
        if (adjacency_count < 0) {
            throw std::runtime_error("negative adjacent-vertex count for configuration vertex " +
                                     std::to_string(index + 1));
        }
        internal_vertex_seen[index] = true;
        int neighbor;
        while (line >> neighbor) {
            if (neighbor < 1 || neighbor > n) {
                throw std::runtime_error("configuration neighbor index out of range at vertex " +
                                         std::to_string(index + 1));
            }
            result.rotations[index].push_back(neighbor - 1);
        }
        if (static_cast<int>(result.rotations[index].size()) != adjacency_count) {
            throw std::runtime_error(
                "adjacent-vertex count does not match the rotation at configuration vertex " +
                std::to_string(index + 1));
        }
        adjacent_vertex_count[index] = adjacency_count;
    }
    for (int index = r; index < n; ++index) {
        if (!internal_vertex_seen[index]) {
            throw std::runtime_error("missing configuration vertex line for index " +
                                     std::to_string(index + 1));
        }
    }

    // The ring vertices are implicit. In a free completion each has its two ring neighbors and its
    // interior neighbors. The listed internal rotations determine all interior adjacencies.
    if (r > 0) {
        std::vector<std::vector<int>> interior_neighbors(r);
        for (int v = r; v < n; ++v) {
            for (int u : result.rotations[v]) {
                if (0 <= u && u < r) interior_neighbors[u].push_back(v);
            }
        }
        for (int i = 0; i < r; ++i) {
            const int next = (i + 1) % r;
            const int prev = (i + r - 1) % r;
            result.rotations[i].push_back(next);
            for (int v : interior_neighbors[i]) result.rotations[i].push_back(v);
            result.rotations[i].push_back(prev);
            adjacent_vertex_count[i] = static_cast<int>(result.rotations[i].size());
        }
    }

    std::vector<int> incident_digon_count(n, 0);
    std::set<std::pair<int, int>> seen_digons;
    if (pos < lines.size()) {
        std::istringstream count_line(lines[pos++]);
        int m;
        if (!(count_line >> m)) throw std::runtime_error("bad configuration digon count");
        if (m < 0) throw std::runtime_error("negative configuration digon count");
        for (int i = 0; i < m; ++i) {
            if (pos >= lines.size()) throw std::runtime_error("truncated configuration digons");
            std::istringstream line(lines[pos++]);
            int a, b;
            if (!(line >> a >> b)) throw std::runtime_error("bad configuration digon line");
            --a;
            --b;
            if (a < 0 || b < 0 || a >= n || b >= n || a == b) {
                throw std::runtime_error("invalid configuration digon endpoint");
            }
            const auto endpoints = std::minmax(a, b);
            if (!seen_digons.emplace(endpoints.first, endpoints.second).second) {
                throw std::runtime_error("duplicate configuration digon");
            }
            if (std::find(result.rotations[a].begin(), result.rotations[a].end(), b) ==
                    result.rotations[a].end() ||
                std::find(result.rotations[b].begin(), result.rotations[b].end(), a) ==
                    result.rotations[b].end()) {
                throw std::runtime_error(
                    "configuration digon endpoints are not adjacent in both rotations");
            }
            result.digons.emplace_back(a, b);
            ++incident_digon_count[a];
            ++incident_digon_count[b];
        }
    }
    if (pos != lines.size()) throw std::runtime_error("extra lines in configuration file");

    // FORMAT.md defines delta_K(v) as the adjacent-vertex count plus the number of
    // incident digons.  from_vertex_rotations likewise represents each declared digon
    // by one extra dart at each endpoint, so these fixed degrees now agree with the
    // degree of the resulting dart representation.
    for (int v = 0; v < n; ++v) {
        result.prescribed_degree[v] = adjacent_vertex_count[v] + incident_digon_count[v];
    }
    return result;
}

Rule read_rule_file(const std::filesystem::path& path, std::optional<std::size_t>) {
    return read_rule_like(path, false, std::nullopt);
}

Rule read_combined_rule_file(const std::filesystem::path& path, std::size_t rule_count) {
    return read_rule_like(path, true, rule_count);
}

AuxiliaryCover read_auxiliary_cover_file(const std::filesystem::path& path) {
    const auto lines = read_nonempty_text_lines(path);
    if (lines.empty()) throw std::runtime_error("empty auxiliary-rule file " + path.string());

    std::size_t pos = 0;
    AuxiliaryCover result;
    result.base = parse_rule_record(lines, pos, path.filename().string() + "#R",
                                    DigonCountRequirement::Required);

    if (pos >= lines.size()) {
        throw std::runtime_error("auxiliary-rule file lacks the cover count k: " +
                                 path.string());
    }
    const std::vector<int> count =
        parse_integer_line(lines[pos++], "auxiliary-rule cover count " + path.string());
    if (count.size() != 1 || count[0] < 1) {
        throw std::runtime_error("bad auxiliary-rule cover count k in " + path.string());
    }

    const int cover_count = count[0];
    result.cover.reserve(static_cast<std::size_t>(cover_count));
    for (int i = 0; i < cover_count; ++i) {
        result.cover.push_back(parse_rule_record(
            lines, pos,
            path.filename().string() + "#R_" + std::to_string(i + 1),
            DigonCountRequirement::Required));
    }
    if (pos != lines.size()) {
        throw std::runtime_error("extra lines in auxiliary-rule file " + path.string());
    }
    return result;
}

Cartwheel read_cartwheel_file(const std::filesystem::path& path) {
    const auto lines = read_nonempty_text_lines(path);
    if (lines.empty()) throw std::runtime_error("empty cartwheel file");
    std::size_t pos = 0;
    std::istringstream header(lines[pos++]);
    int n, center;
    if (!(header >> n >> center)) throw std::runtime_error("bad cartwheel header");
    --center;
    std::vector<std::vector<int>> rotations(n);
    std::vector<DegreeRange> ranges(n);
    for (int row = 0; row < n; ++row) {
        if (pos >= lines.size()) throw std::runtime_error("truncated cartwheel");
        std::istringstream line(lines[pos++]);
        int index, lo, hi;
        if (!(line >> index >> lo >> hi)) throw std::runtime_error("bad cartwheel line");
        --index;
        ranges[index] = {lo, parse_upper(hi)};
        int neighbor;
        while (line >> neighbor) rotations[index].push_back(neighbor == -1 ? -1 : neighbor - 1);
    }
    int m = 0;
    if (pos < lines.size()) {
        std::istringstream line(lines[pos++]);
        line >> m;
    }
    std::vector<std::pair<int, int>> digons;
    for (int i = 0; i < m; ++i) {
        std::istringstream line(lines[pos++]);
        int a, b;
        line >> a >> b;
        digons.emplace_back(a - 1, b - 1);
    }
    Embedding z = from_vertex_rotations(n, rotations, digons,
                                        EmbeddingKind::PseudoTriangulationWithDigons);
    z.degree_range = ranges;
    std::vector<DartId> spokes;
    if (z.is_boundary(center)) throw std::runtime_error("cartwheel center is boundary");
    DartId start = z.darts_at(center).front();
    DartId e = start;
    do {
        spokes.push_back(e);
        e = z.darts[e].succ;
    } while (e != start);
    return {std::move(z), center, std::move(spokes)};
}

Island read_island_file(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open " + path.string());
    std::string line;
    auto next_raw_line = [&]() -> std::string {
        if (!std::getline(input, line)) throw std::runtime_error("truncated island file");
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.resize(comment);
        return line;
    };
    // Leading blank lines are not semantic; the blank line after N is semantic when there are no rings.
    do {
        if (!std::getline(input, line)) throw std::runtime_error("empty island file");
    } while (line.find_first_not_of(" \t\r\n") == std::string::npos);
    Island island;
    int n;
    { std::istringstream first(line); if (!(first >> n)) throw std::runtime_error("bad island vertex count"); }
    { std::istringstream rings(next_raw_line()); int r; while (rings >> r) island.ring_sizes.push_back(r); }
    { std::istringstream count(next_raw_line()); if (!(count >> island.degree_two_vertices)) throw std::runtime_error("bad island degree-two count"); }
    for (int i = 0; i < n; ++i) {
        std::istringstream row_line(next_raw_line());
        std::array<int, 3> row;
        if (!(row_line >> row[0] >> row[1] >> row[2])) throw std::runtime_error("bad island row");
        island.incident_edges.push_back(row);
    }
    return island;
}

void write_rule_file(const Rule& rule, const std::filesystem::path& path, bool combined) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    const int s = rule.graph.tail(rule.distinguished);
    const int t = rule.graph.darts[rule.distinguished].head;
    out << '\n' << rule.graph.vertices().size() << ' ' << s + 1 << ' ' << t + 1 << ' '
        << rule.charge << '\n';
    for (VertexId v : rule.graph.vertices()) {
        out << v + 1 << ' ' << rule.graph.degree_range[v].lower << ' '
            << print_upper(rule.graph.degree_range[v].upper);
        for (int u : incidence_rotation(rule.graph, v)) out << ' ' << (u == -1 ? -1 : u + 1);
        out << '\n';
    }
    const auto digons = digon_endpoints(rule.graph);
    if (combined || !digons.empty()) out << digons.size() << '\n';
    for (auto [a, b] : digons) out << a + 1 << ' ' << b + 1 << '\n';
    if (combined) {
        for (bool member : rule.members) out << (member ? '1' : '0');
        out << '\n';
    }
}

void write_cartwheel_file(const Cartwheel& cartwheel, const std::filesystem::path& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << '\n' << cartwheel.graph.vertices().size() << ' ' << cartwheel.center + 1 << '\n';
    for (VertexId v : cartwheel.graph.vertices()) {
        out << v + 1 << ' ' << cartwheel.graph.degree_range[v].lower << ' '
            << print_upper(cartwheel.graph.degree_range[v].upper);
        for (int u : incidence_rotation(cartwheel.graph, v)) out << ' ' << (u == -1 ? -1 : u + 1);
        out << '\n';
    }
    const auto digons = digon_endpoints(cartwheel.graph);
    out << digons.size() << '\n';
    for (auto [a, b] : digons) out << a + 1 << ' ' << b + 1 << '\n';
}

void write_island_file(const Island& island, const std::filesystem::path& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << island.incident_edges.size() << '\n';
    for (std::size_t i = 0; i < island.ring_sizes.size(); ++i) {
        if (i) out << ' ';
        out << island.ring_sizes[i];
    }
    out << '\n' << island.degree_two_vertices << '\n';
    for (const auto& row : island.incident_edges) {
        out << row[0] << ' ' << row[1] << ' ' << row[2] << '\n';
    }
}

std::vector<std::filesystem::path> sorted_regular_files(const std::filesystem::path& directory,
                                                        std::string_view extension) {
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (!entry.is_regular_file()) continue;
        if (!extension.empty() && entry.path().extension() != extension) continue;
        files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end(), [](const auto& a, const auto& b) {
        return a.filename().string() < b.filename().string();
    });
    return files;
}

std::vector<Rule> load_rules(const std::filesystem::path& directory) {
    std::vector<Rule> rules;
    for (const auto& path : sorted_regular_files(directory, ".rule")) {
        rules.push_back(read_rule_file(path));
    }
    for (std::size_t i = 0; i < rules.size(); ++i) {
        rules[i].members.assign(rules.size(), false);
        rules[i].members[i] = true;
    }
    return rules;
}

std::vector<AuxiliaryCover> load_auxiliary_covers(
    const std::filesystem::path& directory) {
    std::vector<AuxiliaryCover> covers;
    for (const auto& path : sorted_regular_files(directory, ".rule_auxiliary")) {
        covers.push_back(read_auxiliary_cover_file(path));
    }
    return covers;
}

std::vector<ConfigurationFile> load_configurations(const std::filesystem::path& directory) {
    std::vector<ConfigurationFile> configurations;
    for (const auto& path : sorted_regular_files(directory, ".conf")) {
        configurations.push_back(read_configuration_file(path));
    }
    return configurations;
}

std::vector<RootedConfiguration> build_rooted_configuration_set(
    const std::vector<ConfigurationFile>& configurations) {
    std::vector<RootedConfiguration> result;
    for (const auto& configuration : configurations) {
        auto variants = extend_from_cut_vertices(configuration);
        result.insert(result.end(), std::make_move_iterator(variants.begin()),
                      std::make_move_iterator(variants.end()));
    }
    return result;
}

Embedding outer_extension_from_configuration(const ConfigurationFile& configuration) {
    const int n = configuration.vertex_count;
    const int r = configuration.ring_size;
    if (r < 0 || r > n) throw std::runtime_error("bad ring size in " + configuration.name);

    // The files encode a free completion.  Remove the first r ring vertices and
    // replace every incidence from an internal vertex to a ring vertex by a
    // fresh degree-one outer endpoint.  A ring vertex can have more than one
    // such incidence, and those incidences must become distinct endpoints.
    const int internal_count = n - r;
    std::vector<std::vector<int>> rotations(internal_count);
    for (int old_v = r; old_v < n; ++old_v) {
        const int new_v = old_v - r;
        for (int old_u : configuration.rotations[old_v]) {
            if (old_u < 0 || old_u >= n) {
                throw std::runtime_error("bad rotation entry in " + configuration.name);
            }
            if (old_u < r) {
                const int endpoint = static_cast<int>(rotations.size());
                rotations.push_back({new_v, -1});
                rotations[new_v].push_back(endpoint);
            } else {
                rotations[new_v].push_back(old_u - r);
            }
        }
    }

    std::vector<std::pair<int, int>> digons;
    for (auto [a, b] : configuration.digons) {
        if (a < r || b < r) {
            throw std::runtime_error("a configuration digon is incident with the ring in " +
                                     configuration.name);
        }
        digons.emplace_back(a - r, b - r);
    }

    Embedding z = from_vertex_rotations(
        static_cast<int>(rotations.size()), rotations, digons,
        EmbeddingKind::PseudoEmbedding);
    for (int old_v = r; old_v < n; ++old_v) {
        const int new_v = old_v - r;
        z.degree_range[new_v] = {configuration.prescribed_degree[old_v],
                                 configuration.prescribed_degree[old_v]};
    }
    for (int v = internal_count; v < static_cast<int>(rotations.size()); ++v) {
        z.degree_range[v] = {5, kInfinity};
    }
    return z;
}

}  // namespace apex

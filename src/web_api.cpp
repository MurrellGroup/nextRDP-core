#include "analysis.hpp"
#include "legacy_method_state.hpp"
#include "legacy_optional/maxchi.hpp"
#include "legacy_optional/chimaera.hpp"
#include "legacy_optional/geneconv.hpp"
#include "legacy_optional/threeseq.hpp"
#include "legacy_optional/bootscan.hpp"
#include "legacy_optional/siscan.hpp"
#include "legacy_optional_bridge.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(__EMSCRIPTEN__)
#include <emscripten/emscripten.h>
#define NEXT_RDP_KEEPALIVE EMSCRIPTEN_KEEPALIVE
#else
#define NEXT_RDP_KEEPALIVE
#endif

#if defined(__EMSCRIPTEN__)
// The worker is inside the synchronous source-faithful call while these
// milestones are produced.  Posting directly from the WASM worker lets the
// browser receive them without waiting for the worker's JS queue to resume.
EM_JS(void, next_rdp_web_emit_native_progress,
      (int phase, int round, int processed, int total, int events), {
    if (typeof self !== "undefined" &&
        typeof self.postMessage === "function") {
      self.postMessage({
        type: "native-progress",
        phase: phase,
        scanRound: round,
        processedTriplets: processed,
        totalTriplets: total,
        eventCount: events
      });
    }
});
#else
void next_rdp_web_emit_native_progress(
    int, int, int, int, int) {}
#endif

namespace {

std::string result_json;
std::string error_message;

struct ParsedFasta {
    std::string format = "FASTA";
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};

using AlignmentRecords = std::vector<std::pair<std::string, std::string>>;

std::string trim_text(const std::string_view value) {
    const auto first = value.find_first_not_of(" \t\n\r");
    if (first == std::string_view::npos) return {};
    const auto last = value.find_last_not_of(" \t\n\r");
    return std::string(value.substr(first, last - first + 1));
}

bool starts_with_case_insensitive(
    const std::string_view value, const std::string_view prefix) {
    if (value.size() < prefix.size()) return false;
    for (std::size_t index = 0; index < prefix.size(); ++index) {
        if (std::tolower(static_cast<unsigned char>(value[index])) !=
            std::tolower(static_cast<unsigned char>(prefix[index]))) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> alignment_lines(const std::string_view input) {
    std::vector<std::string> result;
    std::string line;
    line.reserve(256);
    for (const char character : input) {
        if (character == '\r') continue;
        if (character == '\n') {
            result.push_back(std::move(line));
            line.clear();
        } else {
            line.push_back(character);
        }
    }
    result.push_back(std::move(line));
    return result;
}

std::vector<std::string> split_words(const std::string_view line) {
    std::istringstream input{std::string(line)};
    std::vector<std::string> result;
    for (std::string word; input >> word;) result.push_back(std::move(word));
    return result;
}

std::string sequence_fragment(const std::string_view value) {
    std::string result;
    result.reserve(value.size());
    for (const unsigned char character : value) {
        if (std::isalpha(character)) {
            result.push_back(static_cast<char>(std::toupper(character)));
        } else if (character == '-' || character == '.' || character == '?' ||
                   character == '!' || character == '*') {
            result.push_back(static_cast<char>(character));
        }
    }
    return result;
}

void append_alignment_record(
    AlignmentRecords& records,
    std::unordered_map<std::string, std::size_t>& indices,
    std::string name, std::string fragment) {
    if (name.empty() || fragment.empty()) return;
    const auto found = indices.find(name);
    if (found == indices.end()) {
        indices.emplace(name, records.size());
        records.emplace_back(std::move(name), std::move(fragment));
    } else {
        records[found->second].second += fragment;
    }
}

AlignmentRecords parse_fasta_or_gde(
    const std::vector<std::string>& lines) {
    AlignmentRecords records;
    std::string name;
    std::string sequence;
    const auto finish = [&] {
        if (!name.empty()) records.emplace_back(std::move(name), std::move(sequence));
        name.clear();
        sequence.clear();
    };
    for (const auto& raw : lines) {
        const std::string line = trim_text(raw);
        if (line.empty()) continue;
        if (line.front() == '>' || line.front() == '%') {
            finish();
            name = trim_text(std::string_view(line).substr(1));
            const auto whitespace = name.find_first_of(" \t");
            if (whitespace != std::string::npos) name.resize(whitespace);
        } else if (!name.empty() && line.front() != ';') {
            sequence += sequence_fragment(line);
        }
    }
    finish();
    return records;
}

AlignmentRecords parse_clustal(const std::vector<std::string>& lines) {
    AlignmentRecords records;
    std::unordered_map<std::string, std::size_t> indices;
    for (std::size_t index = 1; index < lines.size(); ++index) {
        const auto& raw = lines[index];
        const std::string line = trim_text(raw);
        if (line.empty()) continue;
        if (!raw.empty() && std::isspace(static_cast<unsigned char>(raw.front())) &&
            line.find_first_not_of("*:. \t") == std::string::npos) {
            continue;
        }
        const auto tokens = split_words(line);
        if (tokens.size() < 2) continue;
        append_alignment_record(
            records, indices, tokens[0], sequence_fragment(tokens[1]));
    }
    return records;
}

AlignmentRecords parse_mega(const std::vector<std::string>& lines) {
    AlignmentRecords records;
    std::unordered_map<std::string, std::size_t> indices;
    std::string current_name;
    for (const auto& raw : lines) {
        const std::string line = trim_text(raw);
        if (line.empty() || line.front() == '!') continue;
        if (line.front() == '#') {
            if (starts_with_case_insensitive(line, "#mega") ||
                starts_with_case_insensitive(line, "#title") ||
                starts_with_case_insensitive(line, "#format")) {
                continue;
            }
            const auto space = line.find_first_of(" \t");
            current_name = trim_text(std::string_view(line).substr(
                1, space == std::string::npos ? std::string::npos : space - 1));
            if (space != std::string::npos) {
                append_alignment_record(
                    records, indices, current_name,
                    sequence_fragment(std::string_view(line).substr(space + 1)));
            }
        } else if (!current_name.empty()) {
            append_alignment_record(
                records, indices, current_name, sequence_fragment(line));
        }
    }
    return records;
}

std::pair<std::string, std::string> nexus_record(
    const std::string_view raw) {
    const std::string line = trim_text(raw);
    if (line.empty()) return {};
    if (line.front() == '\'' || line.front() == '"') {
        const char quote = line.front();
        const auto end = line.find(quote, 1);
        if (end == std::string::npos) return {};
        return {line.substr(1, end - 1),
                sequence_fragment(std::string_view(line).substr(end + 1))};
    }
    const auto tokens = split_words(line);
    if (tokens.size() < 2) return {};
    return {tokens[0], sequence_fragment(tokens[1])};
}

AlignmentRecords parse_nexus(const std::vector<std::string>& lines) {
    AlignmentRecords records;
    std::unordered_map<std::string, std::size_t> indices;
    bool in_matrix = false;
    for (const auto& raw : lines) {
        std::string line = trim_text(raw);
        if (line.empty() || line.front() == '[') continue;
        if (!in_matrix) {
            std::string lowercase = line;
            std::transform(
                lowercase.begin(), lowercase.end(), lowercase.begin(),
                [](const unsigned char character) {
                    return static_cast<char>(std::tolower(character));
                });
            const auto matrix = lowercase.find("matrix");
            if (matrix == std::string::npos) continue;
            in_matrix = true;
            line = trim_text(std::string_view(line).substr(matrix + 6));
        }
        const auto semicolon = line.find(';');
        const std::string row = semicolon == std::string::npos
            ? line : line.substr(0, semicolon);
        auto [name, fragment] = nexus_record(row);
        append_alignment_record(
            records, indices, std::move(name), std::move(fragment));
        if (semicolon != std::string::npos) break;
    }
    return records;
}

std::pair<std::string, std::string> phylip_named_row(
    const std::string_view raw) {
    const auto tokens = split_words(trim_text(raw));
    if (tokens.size() >= 2) {
        std::string fragment;
        for (std::size_t index = 1; index < tokens.size(); ++index) {
            fragment += sequence_fragment(tokens[index]);
        }
        return {tokens[0], std::move(fragment)};
    }
    if (raw.size() > 10) {
        const std::string name = trim_text(raw.substr(0, 10));
        const std::string fragment = sequence_fragment(raw.substr(10));
        if (!name.empty() && !fragment.empty()) return {name, fragment};
    }
    return {};
}

bool phylip_complete(
    const AlignmentRecords& records, const std::size_t count,
    const std::size_t length) {
    return records.size() == count &&
        std::all_of(records.begin(), records.end(), [&](const auto& record) {
            return record.second.size() == length;
        });
}

AlignmentRecords parse_phylip_interleaved(
    const std::vector<std::string>& lines, const std::size_t count,
    const std::size_t length) {
    AlignmentRecords records;
    records.reserve(count);
    std::unordered_map<std::string, std::size_t> indices;
    std::size_t line_index = 1;
    while (line_index < lines.size() && records.size() < count) {
        const auto& raw = lines[line_index++];
        if (trim_text(raw).empty()) continue;
        auto [name, fragment] = phylip_named_row(raw);
        if (name.empty() || fragment.empty() || fragment.size() > length ||
            !indices.emplace(name, records.size()).second) {
            return {};
        }
        records.emplace_back(std::move(name), std::move(fragment));
    }
    if (records.size() != count) return {};
    std::size_t continuation = 0;
    for (; line_index < lines.size(); ++line_index) {
        const std::string line = trim_text(lines[line_index]);
        if (line.empty()) {
            continuation = 0;
            continue;
        }
        const auto tokens = split_words(line);
        if (tokens.empty()) continue;
        const auto named = tokens.size() > 1 ? indices.find(tokens[0]) : indices.end();
        std::size_t target = continuation % count;
        std::string fragment;
        if (named != indices.end()) {
            target = named->second;
            for (std::size_t index = 1; index < tokens.size(); ++index) {
                fragment += sequence_fragment(tokens[index]);
            }
        } else {
            fragment = sequence_fragment(line);
        }
        records[target].second += fragment;
        if (records[target].second.size() > length) return {};
        ++continuation;
    }
    return phylip_complete(records, count, length) ? records : AlignmentRecords{};
}

AlignmentRecords parse_phylip_sequential(
    const std::vector<std::string>& lines, const std::size_t count,
    const std::size_t length) {
    AlignmentRecords records;
    records.reserve(count);
    std::size_t line_index = 1;
    while (records.size() < count) {
        while (line_index < lines.size() && trim_text(lines[line_index]).empty()) {
            ++line_index;
        }
        if (line_index >= lines.size()) return {};
        auto [name, sequence] = phylip_named_row(lines[line_index++]);
        if (name.empty() || sequence.empty() || sequence.size() > length) return {};
        while (sequence.size() < length) {
            while (line_index < lines.size() && trim_text(lines[line_index]).empty()) {
                ++line_index;
            }
            if (line_index >= lines.size()) return {};
            sequence += sequence_fragment(lines[line_index++]);
            if (sequence.size() > length) return {};
        }
        records.emplace_back(std::move(name), std::move(sequence));
    }
    return phylip_complete(records, count, length) ? records : AlignmentRecords{};
}

AlignmentRecords parse_phylip(const std::vector<std::string>& lines) {
    if (lines.empty()) return {};
    const auto header = split_words(lines.front());
    if (header.size() < 2) return {};
    std::size_t count = 0;
    std::size_t length = 0;
    try {
        count = static_cast<std::size_t>(std::stoull(header[0]));
        length = static_cast<std::size_t>(std::stoull(header[1]));
    } catch (...) {
        return {};
    }
    if (count == 0 || length == 0) return {};
    auto records = parse_phylip_interleaved(lines, count, length);
    return records.empty()
        ? parse_phylip_sequential(lines, count, length) : records;
}

ParsedFasta parse_alignment(const std::string& text) {
    if (text.empty()) throw std::runtime_error("The selected alignment file is empty.");
    const auto lines = alignment_lines(text);
    std::size_t first_index = 0;
    std::string first;
    for (; first_index < lines.size(); ++first_index) {
        first = trim_text(lines[first_index]);
        if (!first.empty()) break;
    }
    if (first.empty()) {
        throw std::runtime_error("The selected alignment file contains no readable text.");
    }
    const std::vector<std::string> content(
        lines.begin() + static_cast<std::ptrdiff_t>(first_index), lines.end());
    AlignmentRecords records;
    std::string format;
    if (first.front() == '>' || first.front() == '%') {
        format = first.front() == '>' ? "FASTA" : "GDE";
        records = parse_fasta_or_gde(content);
    } else if (starts_with_case_insensitive(first, "clustal") ||
               starts_with_case_insensitive(first, "muscle")) {
        format = "CLUSTAL";
        records = parse_clustal(content);
    } else if (starts_with_case_insensitive(first, "#nexus")) {
        format = "NEXUS";
        records = parse_nexus(content);
    } else if (starts_with_case_insensitive(first, "#mega")) {
        format = "MEGA";
        records = parse_mega(content);
    } else {
        const auto header = split_words(first);
        const auto decimal = [](const std::string& value) {
            return !value.empty() && std::all_of(
                value.begin(), value.end(), [](const unsigned char character) {
                    return std::isdigit(character) != 0;
                });
        };
        if (header.size() >= 2 && decimal(header[0]) && decimal(header[1])) {
            format = "PHYLIP";
            records = parse_phylip(content);
        }
    }
    if (records.empty()) {
        throw std::runtime_error(
            "The alignment format was not recognised. Accepted text formats are FASTA, GDE, CLUSTAL/MUSCLE, PHYLIP, NEXUS, and MEGA.");
    }
    ParsedFasta parsed;
    parsed.format = std::move(format);
    parsed.names.reserve(records.size());
    parsed.sequences.reserve(records.size());
    for (auto& [name, sequence] : records) {
        parsed.names.push_back(std::move(name));
        parsed.sequences.push_back(std::move(sequence));
    }
    if (parsed.names.size() < 3) {
        throw std::runtime_error("The alignment must contain at least three nucleotide sequences.");
    }
    const std::size_t length = parsed.sequences.front().size();
    if (length == 0) throw std::runtime_error("The alignment contains no nucleotide columns.");
    for (std::size_t index = 0; index < parsed.sequences.size(); ++index) {
        if (parsed.sequences[index].size() != length) {
            std::ostringstream message;
            message << "The sequences are not aligned: '" << parsed.names[index]
                    << "' has " << parsed.sequences[index].size()
                    << " columns; expected " << length << '.';
            throw std::runtime_error(message.str());
        }
    }
    return parsed;
}

std::string canonical_fasta(const ParsedFasta& alignment) {
    std::ostringstream output;
    for (std::size_t index = 0; index < alignment.names.size(); ++index) {
        output << '>' << alignment.names[index] << '\n'
               << alignment.sequences[index] << '\n';
    }
    return output.str();
}

struct WebContext {
    std::string fasta;
    ParsedFasta alignment;
    std::string cache;
    std::string error;
    RdpInitialAnalysisOptions options;
    RdpFullAnalysisResult full;
    bool loaded = false;
    bool started = false;
    bool finished = false;
    int progress_phase = 0;
    int progress_round = 1;
    int progress_processed_triplets = 0;
    int progress_total_triplets = 0;
    int progress_event_count = 0;
    std::vector<unsigned char> masked;
    std::vector<unsigned char> disabled;
    std::vector<unsigned int> reference_groups;
    std::vector<int> review_states;
    alignas(4) std::atomic<std::int32_t> cancel_requested{0};
};

bool web_cancellation_callback(void* user) {
    const auto* context = static_cast<const WebContext*>(user);
    return context != nullptr &&
        context->cancel_requested.load(std::memory_order_relaxed) != 0;
}

void web_progress_callback(
    const int phase, const int round, const int processed_triplets,
    const int total_triplets, const int event_count, void* user) {
    auto* context = static_cast<WebContext*>(user);
    if (context == nullptr) return;
    context->progress_phase = phase;
    context->progress_round = round;
    context->progress_processed_triplets = processed_triplets;
    context->progress_total_triplets = total_triplets;
    context->progress_event_count = event_count;
    next_rdp_web_emit_native_progress(
        phase, round, processed_triplets, total_triplets, event_count);
}

std::vector<std::unique_ptr<WebContext>>& web_contexts() {
    static std::vector<std::unique_ptr<WebContext>> contexts;
    return contexts;
}

WebContext* web_context(const std::uint32_t handle) {
    if (handle == 0 || handle > web_contexts().size()) return nullptr;
    return web_contexts()[handle - 1].get();
}

const char* cached(WebContext& context, std::string value) {
    context.cache = std::move(value);
    return context.cache.c_str();
}

void json_string(std::ostringstream& output, const std::string& value) {
    output << '"';
    for (const unsigned char character : value) {
        switch (character) {
        case '"': output << "\\\""; break;
        case '\\': output << "\\\\"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (character < 0x20) output << ' ';
            else output << static_cast<char>(character);
        }
    }
    output << '"';
}

int web_base_state(const char character) {
    switch (static_cast<char>(std::toupper(
        static_cast<unsigned char>(character)))) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 3;
    case 'T':
    case 'U': return 4;
    default: return 0;
    }
}

std::size_t web_choose_three(const std::size_t count) {
    return count < 3 ? 0 : count * (count - 1) * (count - 2) / 6;
}

std::vector<unsigned char> suggested_web_mask(
    const ParsedFasta& alignment, const std::vector<float>& similarity,
    const std::vector<std::size_t>& valid_sites,
    const double minimum_pair_identity) {
    const std::size_t count = alignment.sequences.size();
    std::vector<unsigned char> mask(count, 0);
    if (count <= 3) return mask;
    if (minimum_pair_identity >= 1.0 - 1.0e-7) {
        std::fill(mask.begin() + 3, mask.end(), 1);
        return mask;
    }
    const auto identity = [&](const std::size_t first, const std::size_t second) {
        return similarity[first * count + second];
    };
    std::vector<float> closest_identity(count, -1.0F);
    std::vector<std::size_t> closest_sequence(count, count);
    const auto active = [&](const std::size_t index) { return mask[index] == 0; };
    const auto refresh = [&](const std::size_t first) {
        closest_identity[first] = -1.0F;
        closest_sequence[first] = count;
        if (!active(first)) return;
        for (std::size_t second = 0; second < count; ++second) {
            if (first == second || !active(second)) continue;
            const float pair_identity = identity(first, second);
            if (pair_identity > closest_identity[first]) {
                closest_identity[first] = pair_identity;
                closest_sequence[first] = second;
            }
        }
    };
    for (std::size_t index = 0; index < count; ++index) refresh(index);

    std::size_t active_count = count;
    const auto length = static_cast<double>(alignment.sequences.front().size());
    while (active_count >= 4) {
        const double correction = static_cast<double>(
            web_choose_three(active_count)) * (length / 30.0);
        if (!(correction > 0.05)) break;
        const double minimum_distance =
            std::log(correction / 0.05) / std::log(4.0) / length;
        float best_identity = -1.0F;
        std::size_t first = count;
        std::size_t second = count;
        for (std::size_t index = 0; index < count; ++index) {
            if (!active(index) || closest_identity[index] <= best_identity) continue;
            best_identity = closest_identity[index];
            first = index;
            second = closest_sequence[index];
        }
        if (first == count || second == count || !active(second) ||
            1.0 - static_cast<double>(best_identity) >= minimum_distance) {
            break;
        }
        std::size_t remove = first;
        const double first_length = static_cast<double>(valid_sites[first]);
        const double second_length = static_cast<double>(valid_sites[second]);
        if (second_length > first_length * 0.95 &&
            first_length > second_length * 0.95) {
            double first_distance = 0.0;
            double second_distance = 0.0;
            for (std::size_t other = 0; other < count; ++other) {
                first_distance += 1.0 - identity(first, other);
                second_distance += 1.0 - identity(second, other);
            }
            remove = first_distance > second_distance ? second : first;
        } else {
            remove = first_length > second_length ? second : first;
        }
        mask[remove] = 1;
        --active_count;
        closest_identity[remove] = -1.0F;
        closest_sequence[remove] = count;
        for (std::size_t index = 0; index < count; ++index) {
            if (active(index) && closest_sequence[index] == remove) refresh(index);
        }
    }
    return mask;
}

std::string summary_json(const WebContext& context) {
    const auto count = context.alignment.sequences.size();
    const auto length = context.alignment.sequences.front().size();
    std::size_t variable_sites = 0;
    std::size_t informative_sites = 0;
    std::vector<std::size_t> valid_sites(count, 0);
    for (std::size_t position = 0; position < length; ++position) {
        std::array<int, 5> states{};
        for (std::size_t sequence = 0; sequence < count; ++sequence) {
            const int state = web_base_state(
                context.alignment.sequences[sequence][position]);
            if (state != 0) {
                ++states[state];
                ++valid_sites[sequence];
            }
        }
        int observed = 0;
        int repeated_states = 0;
        for (int state = 1; state <= 4; ++state) {
            if (states[state] > 0) ++observed;
            if (states[state] > 1) ++repeated_states;
        }
        if (observed > 1) ++variable_sites;
        if (repeated_states > 1) ++informative_sites;
    }
    double minimum_identity = 1.0;
    double identity_sum = 0.0;
    std::size_t identity_pairs = 0;
    std::vector<float> pair_similarity(count * count, 1.0F);
    for (std::size_t first = 0; first < count; ++first) {
        for (std::size_t second = first + 1; second < count; ++second) {
            std::size_t compared = 0;
            std::size_t matches = 0;
            for (std::size_t position = 0; position < length; ++position) {
                const int left = web_base_state(
                    context.alignment.sequences[first][position]);
                const int right = web_base_state(
                    context.alignment.sequences[second][position]);
                if (left == 0 || right == 0) continue;
                ++compared;
                if (left == right) ++matches;
            }
            const float identity = compared == 0 ? 0.0F
                : static_cast<float>(matches) / static_cast<float>(compared);
            pair_similarity[first * count + second] = identity;
            pair_similarity[second * count + first] = identity;
            minimum_identity = std::min(
                minimum_identity, static_cast<double>(identity));
            identity_sum += static_cast<double>(identity);
            ++identity_pairs;
        }
    }
    const auto suggested_mask = suggested_web_mask(
        context.alignment, pair_similarity, valid_sites, minimum_identity);
    const std::size_t masked_count = static_cast<std::size_t>(std::count(
        suggested_mask.begin(), suggested_mask.end(), 1));
    const std::size_t active_count = count - masked_count;
    const double recommended_distance = length == 0 ? 0.0
        : 2.0 * std::log(4.0 * static_cast<double>(
            std::max<std::size_t>(1, active_count))) /
            static_cast<double>(length);
    std::ostringstream output;
    output << "{\"format\":";
    json_string(output, context.alignment.format);
    output << ",\"sequenceCount\":" << count
           << ",\"alignmentLength\":" << length
           << ",\"activeSequenceCount\":" << active_count
           << ",\"tripletCount\":" << web_choose_three(active_count)
           << ",\"variableSiteCount\":" << variable_sites
           << ",\"informativeSiteCount\":" << informative_sites
           << ",\"minimumPairIdentity\":" << (identity_pairs == 0 ? 0.0 : minimum_identity)
           << ",\"meanPairIdentity\":" << (identity_pairs == 0 ? 0.0 : identity_sum / identity_pairs)
           << ",\"recommendedMinimumDistance\":" << recommended_distance
           << ",\"partitionBoundaries\":[1," << length << "]"
           << ",\"sequences\":[";
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0) output << ',';
        const std::size_t missing = length - valid_sites[index];
        output << "{\"index\":" << index << ",\"name\":";
        json_string(output, context.alignment.names[index]);
        output << ",\"validSites\":" << valid_sites[index]
               << ",\"missingSites\":" << missing
               << ",\"missingFraction\":" << std::setprecision(17)
               << static_cast<double>(missing) / static_cast<double>(length)
               << ",\"masked\":" << (suggested_mask[index] ? "true" : "false")
               << '}';
    }
    output << "],\"warnings\":[";
    bool first_warning = true;
    auto warning = [&](const std::string& message) {
        if (!first_warning) output << ',';
        first_warning = false;
        json_string(output, message);
    };
    if (minimum_identity < 0.6) {
        warning("At least one sequence pair is below 60% identity; the RDP5 manual warns that such alignments are especially vulnerable to false signals.");
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (valid_sites[index] == 0) {
            warning("Sequence has no unambiguous nucleotide sites: " +
                    context.alignment.names[index]);
        }
    }
    if (masked_count > 0) {
        warning(std::to_string(masked_count) +
                " near-identical sequence(s) were auto-masked using the supplied RDP5 optimisation workflow; you can change the selection before scanning.");
    }
    output << "]}";
    return output.str();
}

bool plot_base(char value) {
    switch (value) {
    case 'A': case 'a': case 'C': case 'c':
    case 'G': case 'g': case 'T': case 't': case 'U': case 'u':
        return true;
    default:
        return false;
    }
}

char plot_base_at(const std::string& sequence, int coordinate, bool circular) {
    if (sequence.empty()) return 0;
    const int length = static_cast<int>(sequence.size());
    if (circular) {
        coordinate = (coordinate - 1) % length;
        if (coordinate < 0) coordinate += length;
    } else {
        coordinate = std::clamp(coordinate - 1, 0, length - 1);
    }
    return sequence[static_cast<std::size_t>(coordinate)];
}

double plot_pair_identity(
    const std::array<const std::string*, 3>& sequences, int center,
    int first, int second, int half_window, bool circular) {
    int matches = 0;
    int compared = 0;
    for (int offset = -half_window; offset <= half_window; ++offset) {
        const char left = plot_base_at(*sequences[first], center + offset, circular);
        const char right = plot_base_at(*sequences[second], center + offset, circular);
        if (!plot_base(left) || !plot_base(right)) continue;
        ++compared;
        if (std::toupper(static_cast<unsigned char>(left)) ==
            std::toupper(static_cast<unsigned char>(right))) ++matches;
    }
    return compared == 0 ? 0.0 : static_cast<double>(matches) / compared;
}

// RdpFinalEvent keeps the three source prefix sequences in their discovery
// order and stores the consensus winner separately.  The browser contracts,
// like the old review panels, are expressed in recombinant/major/minor order;
// centralising this rotation prevents a winner in slot 1 or 2 from silently
// relabelling every lazy plot and breakpoint panel.
std::array<int, 3> event_role_sequences(const RdpFinalEvent& event) {
    const int winner = std::clamp(event.winning_role, 0, 2);
    return {
        event.representative_sequences[winner],
        event.representative_sequences[(winner + 1) % 3],
        event.representative_sequences[(winner + 2) % 3],
    };
}

int event_role_original_index(const RdpFinalEvent& event, const int role) {
    const int winner = std::clamp(event.winning_role, 0, 2);
    return (winner + std::clamp(role, 0, 2)) % 3;
}

int pair_slot(const int first, const int second) {
    const int low = std::min(first, second);
    const int high = std::max(first, second);
    if (low == 0 && high == 1) return 0;
    if (low == 0 && high == 2) return 1;
    return 2;
}

void write_maxchi_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::MaxChiDiscoveryCandidate& discovery);
void write_chimaera_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::ChimaeraDiscoveryCandidate& discovery);
void write_geneconv_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::GeneconvDiscoveryCandidate& discovery);
void write_threeseq_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::ThreeSeqDiscoveryCandidate& discovery);

template <typename Candidate>
const Candidate* nearest_discovery_candidate(
    const std::vector<Candidate>& candidates, const RdpFinalEvent& event) {
    if (candidates.empty()) return nullptr;
    const auto distance = [&event](const Candidate& candidate) {
        const auto begin = static_cast<std::int64_t>(candidate.beginning);
        const auto end = static_cast<std::int64_t>(candidate.ending);
        return std::llabs(begin - static_cast<std::int64_t>(event.beginning)) +
            std::llabs(end - static_cast<std::int64_t>(event.ending));
    };
    return &*std::min_element(
        candidates.begin(), candidates.end(),
        [&](const Candidate& first, const Candidate& second) {
            return distance(first) < distance(second);
        });
}

std::string signal_plot_json(const WebContext& context, std::uint32_t signal_id) {
    const auto& result = context.full;
    if (signal_id >= result.events.size()) return {};
    const auto& event = result.events[signal_id];
    const int length = context.alignment.sequences.empty()
        ? 0 : static_cast<int>(context.alignment.sequences.front().size());
    if (length <= 0) return {};
    const auto role_sequences = event_role_sequences(event);
    // The old review page plots the selected signal's discovery triplet, not
    // the reconciled event's winner-rotated role triplet.  Keep this order
    // deterministic (the source AnalysisList is canonical sequence order)
    // while leaving event/alignment/tree views in recombinant/parent order.
    std::array<int, 3> plot_sequences = role_sequences;
    if (event.profile_sequences_available) {
        plot_sequences = event.profile_sequences;
    }
    std::sort(plot_sequences.begin(), plot_sequences.end());
    std::array<const std::string*, 3> sequences{};
    for (int role = 0; role < 3; ++role) {
        const int index = plot_sequences[role];
        if (index < 0 || index >= static_cast<int>(context.alignment.sequences.size())) return {};
        sequences[role] = &context.alignment.sequences[static_cast<std::size_t>(index)];
    }
    const int program = event.program;
    const char* method = program == 1 ? "GENECONV"
        : program == 3 ? "MAXCHI"
        : program == 4 ? "CHIMAERA"
        : program == 5 ? "BOOTSCAN"
        : program == 6 ? "SISCAN"
        : program == 8 ? "3SEQ" : "RDP";
    const char* metric = program == 1 ? "negative-log10-p-value"
        : program == 8 ? "random-walk-height"
        : program == 5 ? "bootstrap-support"
        : program == 6 ? "sister-scan-z-score"
        : program == 0 ? "pair-identity" : "chi-square";
    int requested_window = program == 3 ? context.options.maxchi_window_sites
        : program == 4 ? context.options.chimaera_window_sites
        : program == 5 ? context.options.bootscan_window_sites
        : program == 6 ? context.options.siscan_window_sites
        : context.options.window_sites;
    requested_window = std::max(2, requested_window);
    const int half_window = std::max(1, requested_window / 2);
    // Each cyclic RDP event now retains the XOverHomologyP trace produced on
    // the live detection alignment immediately before tract erasure. Use that
    // source-shaped profile for every RDP event. Reconstructing later events
    // from raw pair identity on the original alignment changes both the site
    // set and the vertical scale, which is visibly not an RDP DrawPlots curve.
    const bool exact_rdp_profile = program == 0 &&
        event.rdp_profile.exact &&
        event.rdp_profile.positions.size() ==
            event.rdp_profile.counts[0].size() &&
        event.rdp_profile.positions.size() ==
            event.rdp_profile.counts[1].size() &&
        event.rdp_profile.positions.size() ==
            event.rdp_profile.counts[2].size();
    bool exact_bootscan_profile = false;
    bool exact_siscan_profile = false;
    bool exact_geneconv_profile = false;
    bool exact_maxchi_profile = false;
    bool exact_chimaera_profile = false;
    bool exact_threeseq_profile = false;
    std::string maxchi_discovery_json = "null";
    std::string chimaera_discovery_json = "null";
    std::string geneconv_discovery_json = "null";
    std::string threeseq_discovery_json = "null";
    int optional_target_local = -1;
    int optional_window_sites = requested_window;
    std::vector<std::size_t> optional_coordinates;
    std::array<std::vector<double>, 3> optional_values;
    double optional_minimum = std::numeric_limits<double>::infinity();
    double optional_maximum = -std::numeric_limits<double>::infinity();
    std::optional<next_rdp_legacy_optional::Alignment> optional_alignment;
    if (program == 1 || program == 3 || program == 4 || program == 8) {
        try {
            optional_alignment.emplace(
                next_rdp_legacy_optional_bridge::make_alignment(
                    context.alignment.sequences));
        } catch (const std::exception&) {
            optional_alignment.reset();
        }
    }
    const std::vector<std::uint8_t> no_missing(
        static_cast<std::size_t>(length), 0);
    const std::array<std::uint32_t, 3> role_triplet{
        static_cast<std::uint32_t>(plot_sequences[0]),
        static_cast<std::uint32_t>(plot_sequences[1]),
        static_cast<std::uint32_t>(plot_sequences[2])};
    const auto optional_missing = std::vector<std::uint8_t>(
        static_cast<std::size_t>(length), 0);
    const auto optional_similarities = optional_alignment
        ? next_rdp_legacy_optional_bridge::pair_similarity(
            *optional_alignment, role_triplet)
        : std::array<double, 3>{};
    if (program == 1 && optional_alignment) {
        next_rdp_legacy_optional::GeneconvDiscoveryOptions options;
        options.circular = context.options.circular;
        options.bonferroni = context.options.correction_bonferroni;
        options.p_value_cutoff = context.options.p_value_cutoff;
        options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        options.mismatch_scale = static_cast<std::size_t>(std::max(
            1, context.options.geneconv_mismatch_scale));
        options.maximum_overlapping_fragments = static_cast<std::size_t>(std::max(
            1, context.options.geneconv_max_overlaps));
        next_rdp_legacy_optional::GeneconvWorkspace workspace;
        next_rdp_legacy_optional::MaxChiWorkspace variable_workspace;
        next_rdp_legacy_optional::MaxChiRecheckOptions variable_options;
        variable_options.circular = context.options.circular;
        variable_options.bonferroni = context.options.correction_bonferroni;
        variable_options.p_value_cutoff = context.options.p_value_cutoff;
        variable_options.correction_tests = options.correction_tests;
        (void)next_rdp_legacy_optional::maxchi_recheck(
            *optional_alignment, role_triplet, optional_missing,
            variable_options, variable_workspace);
        std::vector<next_rdp_legacy_optional::GeneconvDiscoveryCandidate> candidates;
        (void)next_rdp_legacy_optional::geneconv_discover_prepared(
            variable_workspace, optional_similarities, options, workspace,
            candidates);
        if (const auto* selected = nearest_discovery_candidate(candidates, event);
            selected != nullptr) {
            std::ostringstream discovery;
            write_geneconv_discovery_json(discovery, *selected);
            geneconv_discovery_json = discovery.str();
        }
        const auto profile = next_rdp_legacy_optional::geneconv_plot_profile(
            *optional_alignment, role_triplet, options, workspace);
        if (profile.available) {
            exact_geneconv_profile = true;
            optional_coordinates = profile.coordinates;
            optional_window_sites = 0;
            for (int pair = 0; pair < 3; ++pair) {
                optional_values[pair] = profile.negative_log10_p_value[pair];
                for (const double value : optional_values[pair]) {
                    optional_minimum = std::min(optional_minimum, value);
                    optional_maximum = std::max(optional_maximum, value);
                }
            }
        }
    } else if (program == 3 && optional_alignment) {
        next_rdp_legacy_optional::MaxChiDiscoveryOptions discovery_options;
        discovery_options.circular = context.options.circular;
        discovery_options.bonferroni = context.options.correction_bonferroni;
        discovery_options.p_value_cutoff = context.options.p_value_cutoff;
        discovery_options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        discovery_options.fixed_window_sites = std::max<std::size_t>(
            2, static_cast<std::size_t>(context.options.maxchi_window_sites));
        std::vector<next_rdp_legacy_optional::MaxChiDiscoveryCandidate> candidates;
        next_rdp_legacy_optional::MaxChiWorkspace workspace;
        (void)next_rdp_legacy_optional::maxchi_discover(
            *optional_alignment, role_triplet, optional_missing,
            discovery_options, workspace, candidates);
        if (const auto* selected = nearest_discovery_candidate(candidates, event);
            selected != nullptr) {
            std::ostringstream discovery;
            write_maxchi_discovery_json(discovery, *selected);
            maxchi_discovery_json = discovery.str();
        }
        const auto profile = next_rdp_legacy_optional::maxchi_plot_profile(
            *optional_alignment, role_triplet, no_missing,
            context.options.circular,
            std::max<std::size_t>(2, context.options.maxchi_window_sites),
            context.options.p_value_cutoff, workspace);
        if (profile.available) {
            exact_maxchi_profile = true;
            optional_coordinates = profile.coordinates;
            optional_window_sites = profile.half_window * 2;
            for (int pair = 0; pair < 3; ++pair) {
                optional_values[pair] = profile.chi_square[pair];
                for (const double value : optional_values[pair]) {
                    optional_minimum = std::min(optional_minimum, value);
                    optional_maximum = std::max(optional_maximum, value);
                }
            }
        }
    } else if (program == 4 && optional_alignment) {
        const int source_target = event.method_target_role >= 0 &&
            event.method_target_role < 3
            ? event.representative_sequences[event.method_target_role]
            : role_sequences[0];
        optional_target_local = 0;
        for (int role = 0; role < 3; ++role) {
            if (plot_sequences[role] == source_target) {
                optional_target_local = role;
                break;
            }
        }
        next_rdp_legacy_optional::MaxChiWorkspace variable_workspace;
        next_rdp_legacy_optional::MaxChiWorkspace target_workspace;
        next_rdp_legacy_optional::MaxChiRecheckOptions variable_options;
        variable_options.circular = context.options.circular;
        variable_options.bonferroni = context.options.correction_bonferroni;
        variable_options.p_value_cutoff = context.options.p_value_cutoff;
        variable_options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        (void)next_rdp_legacy_optional::maxchi_recheck(
            *optional_alignment, role_triplet, optional_missing,
            variable_options, variable_workspace);
        next_rdp_legacy_optional::ChimaeraDiscoveryOptions discovery_options;
        discovery_options.circular = context.options.circular;
        discovery_options.bonferroni = context.options.correction_bonferroni;
        discovery_options.p_value_cutoff = context.options.p_value_cutoff;
        discovery_options.correction_tests = variable_options.correction_tests;
        discovery_options.fixed_window_sites = std::max<std::size_t>(
            2, static_cast<std::size_t>(context.options.chimaera_window_sites));
        std::vector<next_rdp_legacy_optional::ChimaeraDiscoveryCandidate> candidates;
        (void)next_rdp_legacy_optional::chimaera_discover_prepared(
            variable_workspace, optional_missing, optional_similarities,
            discovery_options, target_workspace, candidates);
        if (const auto* selected = nearest_discovery_candidate(candidates, event);
            selected != nullptr) {
            std::ostringstream discovery;
            write_chimaera_discovery_json(discovery, *selected);
            chimaera_discovery_json = discovery.str();
        }
        const auto profile = next_rdp_legacy_optional::chimaera_plot_profile(
            *optional_alignment, role_triplet,
            static_cast<std::uint8_t>(optional_target_local), no_missing,
            context.options.circular,
            std::max<std::size_t>(2, context.options.chimaera_window_sites),
            context.options.p_value_cutoff, variable_workspace,
            target_workspace);
        if (profile.available) {
            exact_chimaera_profile = true;
            optional_coordinates = profile.coordinates;
            optional_window_sites = profile.half_window * 2;
            const std::size_t trace_pair = profile.target_local == 0
                ? 0 : profile.target_local == 1 ? 2 : 1;
            optional_values[trace_pair] = profile.chi_square;
            for (const double value : profile.chi_square) {
                optional_minimum = std::min(optional_minimum, value);
                optional_maximum = std::max(optional_maximum, value);
            }
        }
    } else if (program == 8 && optional_alignment) {
        next_rdp_legacy_optional::ThreeSeqWorkspace workspace;
        next_rdp_legacy_optional::MaxChiWorkspace variable_workspace;
        next_rdp_legacy_optional::MaxChiRecheckOptions variable_options;
        variable_options.circular = context.options.circular;
        variable_options.bonferroni = context.options.correction_bonferroni;
        variable_options.p_value_cutoff = context.options.p_value_cutoff;
        variable_options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        (void)next_rdp_legacy_optional::maxchi_recheck(
            *optional_alignment, role_triplet, optional_missing,
            variable_options, variable_workspace);
        next_rdp_legacy_optional::ThreeSeqDiscoveryOptions discovery_options;
        discovery_options.circular = context.options.circular;
        discovery_options.correction_enabled = context.options.correction_bonferroni;
        discovery_options.p_value_cutoff = context.options.p_value_cutoff;
        discovery_options.correction_tests = variable_options.correction_tests;
        std::vector<next_rdp_legacy_optional::ThreeSeqDiscoveryCandidate> candidates;
        (void)next_rdp_legacy_optional::threeseq_discover_prepared(
            variable_workspace, optional_missing, optional_similarities,
            discovery_options, workspace, candidates);
        if (const auto* selected = nearest_discovery_candidate(candidates, event);
            selected != nullptr) {
            std::ostringstream discovery;
            write_threeseq_discovery_json(discovery, *selected);
            threeseq_discovery_json = discovery.str();
        }
        const auto profile = next_rdp_legacy_optional::threeseq_plot_profile(
            *optional_alignment, role_triplet, workspace, variable_workspace);
        if (profile.available) {
            exact_threeseq_profile = true;
            optional_coordinates = profile.coordinates;
            optional_window_sites = 0;
            for (int target = 0; target < 3; ++target) {
                optional_values[target] = profile.target_walks[target];
                for (const double value : optional_values[target]) {
                    optional_minimum = std::min(optional_minimum, value);
                    optional_maximum = std::max(optional_maximum, value);
                }
            }
        }
    }
    if (program == 5 && event.bootscan_available) {
        const auto optional_alignment =
            next_rdp_legacy_optional_bridge::make_alignment(
                context.alignment.sequences);
        const std::array<std::uint32_t, 3> triplet{
            static_cast<std::uint32_t>(role_sequences[0]),
            static_cast<std::uint32_t>(role_sequences[1]),
            static_cast<std::uint32_t>(role_sequences[2])};
        next_rdp_legacy_optional::BootscanDiscoveryOptions options;
        options.circular = context.options.circular;
        options.bonferroni = context.options.correction_bonferroni;
        options.p_value_cutoff = context.options.p_value_cutoff;
        options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        options.window_sites = std::max(5, context.options.bootscan_window_sites);
        options.step_sites = std::max(1, context.options.bootscan_step_sites);
        options.bootstrap_replicates = std::max(
            10, context.options.bootscan_bootstrap_replicates);
        options.support_cutoff = context.options.bootscan_support_cutoff;
        options.random_seed = context.options.bootscan_random_seed;
        next_rdp_legacy_optional::BootscanWorkspace workspace;
        const auto profile = next_rdp_legacy_optional::bootscan_plot_profile(
            optional_alignment, triplet,
            std::vector<std::uint8_t>(optional_alignment.length, 0), options,
            workspace);
        if (profile.available) {
            exact_bootscan_profile = true;
            optional_coordinates = profile.coordinates;
            options.window_sites = profile.window_sites;
            for (int pair = 0; pair < 3; ++pair) {
                optional_values[pair] = profile.pair_support[pair];
                for (const double value : optional_values[pair]) {
                    optional_minimum = std::min(optional_minimum, value);
                    optional_maximum = std::max(optional_maximum, value);
                }
            }
        }
    } else if (program == 6 && event.siscan_available) {
        const auto optional_alignment =
            next_rdp_legacy_optional_bridge::make_alignment(
                context.alignment.sequences);
        const std::array<std::uint32_t, 3> triplet{
            static_cast<std::uint32_t>(role_sequences[0]),
            static_cast<std::uint32_t>(role_sequences[1]),
            static_cast<std::uint32_t>(role_sequences[2])};
        std::vector<std::uint32_t> origins(optional_alignment.sequence_count());
        for (std::size_t index = 0; index < origins.size(); ++index) {
            origins[index] = static_cast<std::uint32_t>(index);
        }
        next_rdp_legacy_optional::SiscanOptions options;
        options.circular = context.options.circular;
        options.bonferroni = context.options.correction_bonferroni;
        options.p_value_cutoff = context.options.p_value_cutoff;
        options.correction_tests = std::max<std::uint64_t>(
            1, static_cast<std::uint64_t>(context.full.triplet_count));
        options.window_sites = std::max(5, context.options.siscan_window_sites);
        options.step_sites = std::max(1, context.options.siscan_step_sites);
        options.scan_permutations = std::max(
            10, context.options.siscan_scan_permutations);
        options.p_value_permutations = std::max<std::size_t>(
            options.scan_permutations,
            static_cast<std::size_t>(std::max(1, context.options.siscan_p_value_permutations)));
        options.random_seed = context.options.siscan_random_seed;
        next_rdp_legacy_optional::SiscanWorkspace workspace;
        const auto profile = next_rdp_legacy_optional::siscan_plot_profile(
            optional_alignment, triplet,
            std::vector<std::uint8_t>(optional_alignment.length, 0), origins,
            context.disabled, options, workspace,
            static_cast<std::int64_t>(event.siscan_discovery.outlier_sequence));
        if (profile.available) {
            exact_siscan_profile = true;
            optional_coordinates = profile.coordinates;
            options.window_sites = profile.window_sites;
            for (int pair = 0; pair < 3; ++pair) {
                optional_values[pair] = profile.pair_z[pair];
                for (const double value : optional_values[pair]) {
                    optional_minimum = std::min(optional_minimum, value);
                    optional_maximum = std::max(optional_maximum, value);
                }
            }
        }
    }
    const bool exact_optional_profile =
        exact_geneconv_profile || exact_maxchi_profile ||
        exact_chimaera_profile || exact_threeseq_profile ||
        exact_bootscan_profile || exact_siscan_profile;
    const bool exact_detection_profile =
        exact_rdp_profile || (exact_optional_profile && signal_id == 0);
    const int stride = std::max(1, (length + 1999) / 2000);
    std::vector<int> coordinates;
    if (exact_rdp_profile) {
        coordinates = event.rdp_profile.positions;
    } else if (exact_optional_profile) {
        coordinates.reserve(optional_coordinates.size());
        for (const auto coordinate : optional_coordinates) {
            coordinates.push_back(static_cast<int>(coordinate));
        }
    } else {
        for (int coordinate = 1; coordinate <= length; coordinate += stride) coordinates.push_back(coordinate);
        if (coordinates.empty() || coordinates.back() != length) coordinates.push_back(length);
    }
    auto add_coordinate = [&](int coordinate) {
        coordinate = std::clamp(coordinate, 1, length);
        coordinates.push_back(coordinate);
    };
    // Source-kernel optional profiles are already sampled on their own exact
    // coordinate grid, even when they are reconstructed from the original
    // alignment for a later cyclic event. Injecting breakpoint coordinates
    // into that grid desynchronizes coordinates from the parallel value
    // arrays (and could read beyond those arrays). Only the generic identity
    // fallback needs explicit breakpoint samples.
    if (!exact_rdp_profile && !exact_optional_profile) {
        add_coordinate(event.beginning);
        add_coordinate(event.ending);
        std::sort(coordinates.begin(), coordinates.end());
        coordinates.erase(std::unique(coordinates.begin(), coordinates.end()), coordinates.end());
    }

    std::vector<std::array<double, 3>> values;
    values.reserve(coordinates.size());
    std::array<double, 3> walks{};
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t coordinate_index = 0;
         coordinate_index < coordinates.size(); ++coordinate_index) {
        const int coordinate = coordinates[coordinate_index];
        std::array<double, 3> point{};
        if (exact_rdp_profile) {
            const double divisor = static_cast<double>(
                std::max(1, event.rdp_profile.divisor));
            const auto profile_order = event.profile_sequences_available
                ? event.profile_sequences : role_sequences;
            for (int pair = 0; pair < 3; ++pair) {
                const int first_sequence = plot_sequences[pair == 2 ? 1 : 0];
                const int second_sequence = plot_sequences[pair == 0 ? 1 : 2];
                int first_profile_role = 0;
                int second_profile_role = 0;
                for (int role = 0; role < 3; ++role) {
                    if (profile_order[role] == first_sequence) {
                        first_profile_role = role;
                    }
                    if (profile_order[role] == second_sequence) {
                        second_profile_role = role;
                    }
                }
                const int original_pair = pair_slot(first_profile_role, second_profile_role);
                point[pair] = static_cast<double>(
                    event.rdp_profile.counts[original_pair][coordinate_index]) /
                    divisor;
            }
        } else if (exact_optional_profile &&
                   coordinate_index < optional_coordinates.size()) {
            for (int pair = 0; pair < 3; ++pair) {
                // CHIMAERA has a single target/parent trace; its other two
                // pair slots are intentionally absent. Keep those display
                // lanes at zero instead of indexing empty vectors.
                point[pair] = coordinate_index < optional_values[pair].size()
                    ? optional_values[pair][coordinate_index] : 0.0;
            }
        } else if (program == 8) {
            for (int target = 0; target < 3; ++target) {
                const int first = (target + 1) % 3;
                const int second = (target + 2) % 3;
                const char target_base = plot_base_at(*sequences[target], coordinate, context.options.circular);
                const char first_base = plot_base_at(*sequences[first], coordinate, context.options.circular);
                const char second_base = plot_base_at(*sequences[second], coordinate, context.options.circular);
                if (plot_base(target_base) && plot_base(first_base) && plot_base(second_base)) {
                    const auto target_upper = static_cast<char>(std::toupper(static_cast<unsigned char>(target_base)));
                    const auto first_upper = static_cast<char>(std::toupper(static_cast<unsigned char>(first_base)));
                    const auto second_upper = static_cast<char>(std::toupper(static_cast<unsigned char>(second_base)));
                    if (target_upper == first_upper && target_upper != second_upper) ++walks[target];
                    else if (target_upper == second_upper && target_upper != first_upper) --walks[target];
                }
                point[target] = walks[target];
            }
        } else {
            std::array<double, 3> identity{
                plot_pair_identity(sequences, coordinate, 0, 1, half_window, context.options.circular),
                plot_pair_identity(sequences, coordinate, 0, 2, half_window, context.options.circular),
                plot_pair_identity(sequences, coordinate, 1, 2, half_window, context.options.circular)};
            if (program == 1) {
                for (int pair = 0; pair < 3; ++pair) {
                    const double tail = 1.0 - identity[pair];
                    point[pair] = -std::log10(tail < 1.0e-6 ? 1.0e-6 : tail);
                }
            } else if (program == 3 || program == 4) {
                for (int pair = 0; pair < 3; ++pair) {
                    const double other = (identity[(pair + 1) % 3] + identity[(pair + 2) % 3]) / 2.0;
                    point[pair] = 4.0 * (identity[pair] - other) * (identity[pair] - other) * requested_window;
                }
            } else {
                point = identity;
            }
        }
        values.push_back(point);
        for (const double value : point) {
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
    }
    if (exact_rdp_profile) {
        // SignalPlot defaults the RDP pair-identity domain to [0, 1], even
        // when the observed profile never reaches one. Preserve that source
        // metadata as well as the visual scale.
        minimum = 0.0;
        maximum = 1.0;
    } else if (exact_optional_profile) {
        minimum = optional_minimum;
        maximum = optional_maximum;
    }
    if (!std::isfinite(minimum)) minimum = 0.0;
    if (!std::isfinite(maximum)) maximum = 0.0;
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"signalId\":" << signal_id
           << ",\"windowSites\":" << (exact_rdp_profile ? context.options.window_sites
               : exact_optional_profile ? optional_window_sites
               : requested_window)
           << ",\"alignmentLength\":" << length
           << ",\"method\":\"" << method << "\",\"metric\":\"" << metric
           << "\",\"targetLocal\":";
    if (optional_target_local < 0) output << "null";
    else output << optional_target_local;
    output
           << ",\"maxChiDiscovery\":" << maxchi_discovery_json
           << ",\"chimaeraDiscovery\":" << chimaera_discovery_json
           << ",\"geneconvDiscovery\":" << geneconv_discovery_json
           << ",\"threeSeqDiscovery\":" << threeseq_discovery_json
           << ",\"profileContext\":\""
           << (exact_detection_profile ? "detection-alignment" :
               "original-alignment-reconstruction") << "\""
           << ",\"detectionProfileExact\":"
           << (exact_detection_profile ? "true" : "false")
           << ",\"minimumValue\":" << minimum
           << ",\"maximumValue\":" << maximum << ",\"points\":[";
    // DrawPlots allocates an x=0 sentinel (and repairs its coordinate), but
    // RedrawPlotAA starts at index 1.  Keep that source bookkeeping out of
    // the browser path; otherwise the first displayed segment would be an
    // artificial vertical drop to zero at XDiffPos(1).
    const std::size_t first_display_index = exact_rdp_profile &&
        coordinates.size() > 1 ? 1 : 0;
    bool first_point = true;
    for (std::size_t index = first_display_index; index < coordinates.size(); ++index) {
        if (!first_point) output << ',';
        first_point = false;
        output << "{\"alignmentPosition\":" << coordinates[index]
               << ",\"pair12\":" << values[index][0]
               << ",\"pair13\":" << values[index][1]
               << ",\"pair23\":" << values[index][2] << '}';
    }
    output << "]}";
    return output.str();
}

std::vector<int> event_review_sequences(
    const WebContext& context, const RdpFinalEvent& event, std::size_t limit) {
    std::vector<int> sequences;
    const auto add = [&](const int sequence) {
        if (sequence < 0 || sequence >= static_cast<int>(context.alignment.sequences.size())) return;
        if (std::find(sequences.begin(), sequences.end(), sequence) == sequences.end()) {
            sequences.push_back(sequence);
        }
    };
    for (const int sequence : event.representative_sequences) add(sequence);
    for (const auto& group : event.sequence_groups) for (const int sequence : group) add(sequence);
    // The source review panel keeps the broader alignment context available
    // even when the compact core result has no serialized correlation lists.
    // Append the remaining input rows as neutral evidence rows; this is a
    // display fallback only and does not turn them into event support.
    for (int sequence = 0;
         sequence < static_cast<int>(context.alignment.sequences.size());
         ++sequence) {
        add(sequence);
    }
    const std::size_t candidate_count = sequences.size();
    if (sequences.size() > limit) sequences.resize(limit);
    (void)candidate_count;
    return sequences;
}

std::vector<int> review_window_coordinates(
    int center, int flank, int length, bool circular) {
    std::vector<int> coordinates;
    if (length <= 0) return coordinates;
    flank = std::clamp(flank, 5, 100);
    coordinates.reserve(static_cast<std::size_t>(flank * 2 + 1));
    for (int offset = -flank; offset <= flank; ++offset) {
        int coordinate = center + offset;
        if (circular) {
            coordinate = (coordinate - 1) % length;
            if (coordinate < 0) coordinate += length;
            coordinate += 1;
        } else {
            if (coordinate < 1 || coordinate > length) continue;
        }
        coordinates.push_back(coordinate);
    }
    return coordinates;
}

void confidence_boundary_json(
    std::ostringstream& output, const RdpFinalEvent& event, int boundary,
    const char* name, int input, int polished) {
    const int offset = boundary == 0 ? 0 : 3;
    const int c99_beginning = event.burt.confidence[static_cast<std::size_t>(offset)];
    const int c99_ending = event.burt.confidence[static_cast<std::size_t>(offset + 1)];
    const int hmm = event.burt.confidence[static_cast<std::size_t>(offset + 2)];
    const int c95_beginning = event.burt.confidence[static_cast<std::size_t>(boundary == 0 ? 6 : 8)];
    const int c95_ending = event.burt.confidence[static_cast<std::size_t>(boundary == 0 ? 7 : 9)];
    const bool available = event.burt_attempted && event.burt_applied &&
        !event.burt.intervals.empty() && c99_beginning != 0 && c99_ending != 0;
    output << "{\"name\":";
    json_string(output, name);
    output << ",\"inputCoordinate\":" << input
           << ",\"polishedCoordinate\":" << polished
           << ",\"intervalAvailable\":" << (available ? "true" : "false")
           << ",\"sourceIntervalContainsInput\":" << (available ? "true" : "false")
           << ",\"confidence99\":{\"beginning\":" << c99_beginning
           << ",\"ending\":" << c99_ending
           << ",\"wrapsOrigin\":" << (c99_beginning > c99_ending ? "true" : "false")
           << "},\"hmmCoordinate\":" << hmm
           << ",\"confidence95\":{\"beginning\":" << c95_beginning
           << ",\"ending\":" << c95_ending
           << ",\"wrapsOrigin\":" << (c95_beginning > c95_ending ? "true" : "false")
           << "},\"repositioned\":" << (input != polished ? "true" : "false")
           << ",\"missingDataAdjusted\":false,\"finalGapAdjusted\":false}";
}

std::string event_alignment_json(const WebContext& context, std::uint32_t event_id,
                                 std::uint32_t requested_flank, std::uint32_t row_limit) {
    if (event_id >= context.full.events.size()) return {};
    const auto& event = context.full.events[event_id];
    const auto role_sequences = event_role_sequences(event);
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const int flank = std::clamp(static_cast<int>(requested_flank), 5, 100);
    const auto all_sequences = event_review_sequences(context, event, 64);
    const std::size_t candidate_count = all_sequences.size();
    const std::size_t limit = std::clamp<std::size_t>(row_limit, 3, 64);
    const std::size_t visible_count = std::min(candidate_count, limit);
    const auto is_group_member = [&](const int sequence) {
        const auto& group = event.sequence_groups[static_cast<std::size_t>(
            std::clamp(event.winning_role, 0, 2))];
        return std::find(group.begin(), group.end(), sequence) != group.end();
    };
    const auto input_role = [&](const int sequence) {
        if (!context.options.query_reference_mode || sequence < 0 ||
            sequence >= static_cast<int>(context.reference_groups.size())) {
            return "not-applied";
        }
        return context.reference_groups[static_cast<std::size_t>(sequence)] > 0
            ? "reference" : "query";
    };
    const auto input_group = [&](const int sequence) {
        if (!context.options.query_reference_mode || sequence < 0 ||
            sequence >= static_cast<int>(context.reference_groups.size()) ||
            context.reference_groups[static_cast<std::size_t>(sequence)] == 0) return 0U;
        return context.reference_groups[static_cast<std::size_t>(sequence)];
    };
    const auto role = [&](const int sequence) {
        if (sequence == role_sequences[0]) return "recombinant";
        if (sequence == role_sequences[1]) return "major-parent";
        if (sequence == role_sequences[2]) return "minor-parent";
        return is_group_member(sequence) ? "co-recombinant" : "evidence";
    };
    const auto make_panel = [&](const char* name, int center, int expected_left, int expected_right) {
        const auto coordinates = review_window_coordinates(center, flank, length, context.options.circular);
        std::ostringstream panel;
        panel << "{\"name\":";
        json_string(panel, name);
        const auto center_index = std::find(coordinates.begin(), coordinates.end(), center);
        panel << ",\"center\":" << center
              << ",\"centerIndex\":" << (center_index == coordinates.end() ? 0 : center_index - coordinates.begin())
              << ",\"statisticalConfidence\":";
        confidence_boundary_json(panel, event, std::strcmp(name, "beginning") == 0 ? 0 : 1, name,
                                 center, center);
        panel << ",\"parentTransition\":{\"expectedLeft\":\""
              << (expected_left == 1 ? "major" : "minor")
              << "\",\"expectedRight\":\"" << (expected_right == 1 ? "major" : "minor")
              << "\",\"leftInformativeCoordinate\":null,\"rightInformativeCoordinate\":null"
              << ",\"spanSites\":null,\"supported\":false}"
              << ",\"adjacentErasureEventIds\":[],\"erasureAdjacent\":false"
              << ",\"uncertainErasureEventIds\":[],\"erasureWithinRdpWindow\":false"
              << ",\"uncertainDueToErasure\":false,\"nativeCheckEndsApplied\":false"
              << ",\"nativeCheckEndsWarning\":false,\"informationProfileAvailable\":true"
              << ",\"inputMissingDataInCheckRange\":false,\"linearEdgeWithinRdpWindow\":false"
              << ",\"nativeCheckRange\":{\"beginning\":0,\"ending\":0,\"wrapsOrigin\":false,\"coordinateCount\":0}"
              << ",\"rdpWindowInformativeSites\":" << context.options.window_sites
              << ",\"nearestErasureInformativeSites\":null,\"coordinates\":[";
        for (std::size_t index = 0; index < coordinates.size(); ++index) {
            if (index != 0) panel << ',';
            panel << coordinates[index];
        }
        panel << "]}";
        return panel.str();
    };
    const auto beginning_panel = make_panel("beginning", event.beginning, 1, 2);
    const auto ending_panel = make_panel("ending", event.ending, 2, 1);
    const auto beginning_coordinates = review_window_coordinates(event.beginning, flank, length, context.options.circular);
    const auto ending_coordinates = review_window_coordinates(event.ending, flank, length, context.options.circular);
    std::ostringstream output;
    output << "{\"eventId\":" << event_id << ",\"alignmentLength\":" << length
           << ",\"circular\":" << (context.options.circular ? "true" : "false")
           << ",\"fragmentAssisted\":false,\"requestedFlankSites\":" << flank
           << ",\"candidateRowCount\":" << candidate_count
           << ",\"omittedRowCount\":" << (candidate_count - visible_count) << ",\"rows\":[";
    for (std::size_t row = 0; row < visible_count; ++row) {
        if (row != 0) output << ',';
        const int sequence = all_sequences[row];
        output << "{\"sequenceIndex\":" << sequence << ",\"sequenceName\":";
        json_string(output, context.alignment.names[static_cast<std::size_t>(sequence)]);
        output << ",\"role\":\"" << role(sequence)
               << "\",\"queryReferenceInputRole\":\"" << input_role(sequence)
               << "\",\"referenceGroup\":" << (input_group(sequence) == 0U ? "null" : std::to_string(input_group(sequence)))
               << ",\"masked\":" << (context.masked[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"disabled\":" << (context.disabled[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"currentGroupMember\":" << (is_group_member(sequence) ? "true" : "false")
               << ",\"automaticGroupMember\":" << (is_group_member(sequence) ? "true" : "false")
               << ",\"trace\":false,\"panels\":[";
        for (int panel = 0; panel < 2; ++panel) {
            if (panel != 0) output << ',';
            const auto& coordinates = panel == 0 ? beginning_coordinates : ending_coordinates;
            std::string bases;
            bases.reserve(coordinates.size());
            for (const int coordinate : coordinates) {
                bases.push_back(plot_base_at(context.alignment.sequences[static_cast<std::size_t>(sequence)], coordinate,
                                             context.options.circular));
            }
            json_string(output, bases);
        }
        output << "]}";
    }
    output << "],\"panels\":[" << beginning_panel << ',' << ending_panel << "]}";
    return output.str();
}

std::string event_trees_json(const WebContext& context, std::uint32_t event_id) {
    if (event_id >= context.full.events.size()) return {};
    const auto& event = context.full.events[event_id];
    const auto role_sequences = event_role_sequences(event);
    const auto sequences = event_review_sequences(context, event, 48);
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const auto in_tract = [&](const int coordinate) {
        if (context.options.circular && event.beginning == event.ending) return true;
        if (event.beginning <= event.ending) return coordinate >= event.beginning && coordinate <= event.ending;
        return coordinate >= event.beginning || coordinate <= event.ending;
    };
    const auto region_coordinates = [&](const int region) {
        std::vector<int> coordinates;
        coordinates.reserve(static_cast<std::size_t>(length));
        for (int coordinate = 1; coordinate <= length; ++coordinate) {
            const bool inside = in_tract(coordinate);
            bool keep = false;
            switch (region) {
            case 0: keep = coordinate < event.beginning && !inside; break;
            case 1: keep = inside; break;
            case 2: keep = coordinate > event.ending && !inside; break;
            case 3: keep = inside; break;
            case 4: keep = !inside; break;
            default: keep = inside; break;
            }
            if (keep) coordinates.push_back(coordinate);
        }
        if (coordinates.empty()) {
            for (int coordinate = 1; coordinate <= std::min(length, 60); ++coordinate) coordinates.push_back(coordinate);
        }
        return coordinates;
    };
    const auto role = [&](const int sequence) {
        if (sequence == role_sequences[0]) return "recombinant";
        if (sequence == role_sequences[1]) return "major-parent";
        if (sequence == role_sequences[2]) return "minor-parent";
        return "evidence";
    };
    const auto is_group_member = [&](const int sequence) {
        const auto& group = event.sequence_groups[static_cast<std::size_t>(
            std::clamp(event.winning_role, 0, 2))];
        return std::find(group.begin(), group.end(), sequence) != group.end();
    };
    const auto input_role = [&](const int sequence) {
        if (!context.options.query_reference_mode || sequence < 0 ||
            sequence >= static_cast<int>(context.reference_groups.size())) {
            return "not-applied";
        }
        return context.reference_groups[static_cast<std::size_t>(sequence)] > 0
            ? "reference" : "query";
    };
    const auto input_group = [&](const int sequence) {
        if (!context.options.query_reference_mode || sequence < 0 ||
            sequence >= static_cast<int>(context.reference_groups.size()) ||
            context.reference_groups[static_cast<std::size_t>(sequence)] == 0) return 0U;
        return context.reference_groups[static_cast<std::size_t>(sequence)];
    };
    const auto distance = [&](const int first, const int second, const std::vector<int>& coordinates) {
        int compared = 0;
        int differences = 0;
        const auto& left = context.alignment.sequences[static_cast<std::size_t>(first)];
        const auto& right = context.alignment.sequences[static_cast<std::size_t>(second)];
        for (const int coordinate : coordinates) {
            const char a = plot_base_at(left, coordinate, false);
            const char b = plot_base_at(right, coordinate, false);
            if (!plot_base(a) || !plot_base(b)) continue;
            ++compared;
            if (std::toupper(static_cast<unsigned char>(a)) != std::toupper(static_cast<unsigned char>(b))) ++differences;
        }
        return compared == 0 ? 0.0 : static_cast<double>(differences) / compared;
    };
    struct TreeEdge { int from = 0; int to = 0; double length = 0.0; };
    const auto make_nj_tree = [&](const std::vector<int>& coordinates) {
        const int leaf_count = static_cast<int>(sequences.size());
        const int node_capacity = std::max(1, 2 * leaf_count - 1);
        std::vector<std::vector<double>> matrix(
            static_cast<std::size_t>(node_capacity),
            std::vector<double>(static_cast<std::size_t>(node_capacity), 0.0));
        for (int first = 0; first < leaf_count; ++first) {
            for (int second = first + 1; second < leaf_count; ++second) {
                const double p = distance(sequences[first], sequences[second], coordinates);
                // Clearcut's JC transform returns a bounded large distance
                // when p >= 0.75; retaining that bound keeps the browser tree
                // finite and lets the display preserve source-style negative
                // branch repair without inventing infinite SVG coordinates.
                const double transformed = p >= 0.75
                    ? 10.0 : p <= 0.0 ? 0.0 : -0.75 * std::log(1.0 - (4.0 * p / 3.0));
                matrix[first][second] = matrix[second][first] =
                    std::isfinite(transformed) ? transformed : 10.0;
            }
        }
        std::vector<int> active;
        active.reserve(static_cast<std::size_t>(leaf_count));
        for (int leaf = 0; leaf < leaf_count; ++leaf) active.push_back(leaf);
        std::vector<TreeEdge> edges;
        edges.reserve(static_cast<std::size_t>(2 * leaf_count));
        int next_node = leaf_count;
        while (active.size() > 2) {
            const double denominator = static_cast<double>(active.size() - 2);
            int best_i = active[0], best_j = active[1];
            double best_q = std::numeric_limits<double>::infinity();
            std::vector<double> row_sum(static_cast<std::size_t>(node_capacity), 0.0);
            for (const int node : active) {
                for (const int other : active) row_sum[node] += matrix[node][other];
            }
            for (std::size_t i = 0; i < active.size(); ++i) {
                for (std::size_t j = i + 1; j < active.size(); ++j) {
                    const int left = active[i], right = active[j];
                    const double q = denominator * matrix[left][right] - row_sum[left] - row_sum[right];
                    if (q < best_q) { best_q = q; best_i = left; best_j = right; }
                }
            }
            const double pair_distance = matrix[best_i][best_j];
            const long double left_unclamped = 0.5L * pair_distance +
                (row_sum[best_i] - row_sum[best_j]) / (2.0L * denominator);
            const double left_length = left_unclamped > 0.0L
                ? static_cast<double>(left_unclamped) : 0.0;
            const double right_length = pair_distance > left_length
                ? pair_distance - left_length : 0.0;
            edges.push_back({next_node, best_i, left_length});
            edges.push_back({next_node, best_j, right_length});
            for (const int node : active) {
                if (node == best_i || node == best_j) continue;
                matrix[next_node][node] = matrix[node][next_node] =
                    0.5 * (matrix[best_i][node] + matrix[best_j][node] - pair_distance);
            }
            active.erase(std::remove(active.begin(), active.end(), best_i), active.end());
            active.erase(std::remove(active.begin(), active.end(), best_j), active.end());
            active.push_back(next_node++);
        }
        const int root = next_node;
        if (active.size() == 2) {
            const double half = matrix[active[0]][active[1]] > 0.0
                ? matrix[active[0]][active[1]] / 2.0 : 0.0;
            edges.push_back({root, active[0], half});
            edges.push_back({root, active[1], half});
        } else if (active.size() == 1) {
            edges.push_back({root, active[0], 0.0});
        }
        return std::pair<int, std::vector<TreeEdge>>(root, std::move(edges));
    };
    static constexpr std::array<const char*, 6> names{
        "5-prime-outside", "5-prime-inside", "3-prime-outside",
        "3-prime-inside", "outside-tract", "inside-tract"};
    std::ostringstream output;
    output << "{\"eventId\":" << event_id
           << ",\"method\":\"neighbour-joining\",\"distance\":\"Jukes-Cantor\""
           << ",\"njKernel\":\"source-shaped-neighbor-joining\",\"distanceEncoding\":\"source-tree2arrayp2-midpoint-ranks\""
           << ",\"bootstrapGenerator\":\"disabled-rdp-5.93-event-path\",\"bootstrapSupport\":\"not-applied\""
           << ",\"negativeBranchPolicy\":\"absolute-five-decimal-serialization\",\"analyticalBranchParsing\":\"four-decimal-clamped-complete-edge-repair\""
           << ",\"treeRooting\":\"source-tree2arrayp2-midpoint\",\"collapseEncoding\":\"unbootstrapped-raw-tree-copy\""
           << ",\"displayRooting\":\"arbitrary-internal-node\",\"bootstrapCollapseCutoff\":null"
           << ",\"bootstrapReplicates\":0,\"randomSeed\":3,\"flankVariableSiteTarget\":60"
           << ",\"subsampled\":false,\"sequenceCap\":48,\"fragmentAssisted\":false,\"leaves\":[";
    for (std::size_t index = 0; index < sequences.size(); ++index) {
        if (index != 0) output << ',';
        const int sequence = sequences[index];
        output << "{\"node\":" << index << ",\"workingSequenceIndex\":" << sequence
               << ",\"sequenceIndex\":" << sequence << ",\"sequenceName\":";
        json_string(output, context.alignment.names[static_cast<std::size_t>(sequence)]);
        output << ",\"fragmentEventId\":null,\"role\":\"" << role(sequence)
               << "\",\"queryReferenceInputRole\":\"" << input_role(sequence)
               << "\",\"referenceGroup\":" << (input_group(sequence) == 0U ? "null" : std::to_string(input_group(sequence)))
               << ",\"masked\":" << (context.masked[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"disabled\":" << (context.disabled[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"currentGroupMember\":" << (is_group_member(sequence) ? "true" : "false")
               << ",\"automaticGroupMember\":" << (is_group_member(sequence) ? "true" : "false")
               << ",\"trace\":false}";
    }
    output << "],\"regions\":[";
    for (int region = 0; region < 6; ++region) {
        if (region != 0) output << ',';
        const auto coordinates = region_coordinates(region);
        const bool usable = sequences.size() >= 3 && coordinates.size() >= 3;
        const auto tree = usable ? make_nj_tree(coordinates)
                                 : std::pair<int, std::vector<TreeEdge>>{0, {}};
        const int tree_root = tree.first;
        const int node_count = usable ? tree_root + 1 : 0;
        const int internal_branches = usable
            ? std::max(0, static_cast<int>(sequences.size()) - 3) : 0;
        output << "{\"name\":\"" << names[region] << "\",\"sites\":" << coordinates.size()
               << ",\"sequences\":" << sequences.size() << ",\"usable\":" << (usable ? "true" : "false")
               << ",\"nodeCount\":" << node_count << ",\"root\":" << (usable ? tree_root : 0)
               << ",\"bootstrapReplicates\":0,\"supportedInternalBranches\":0,\"internalBranches\":" << internal_branches
               << ",\"rawDistanceRankLevels\":" << (usable ? 1 : 0)
               << ",\"collapsedDistanceRankLevels\":" << (usable ? 1 : 0)
               << ",\"negativeBranchesNormalized\":0,\"bootstrapRandomSeed\":3,\"edges\":[";
        if (usable) {
            for (std::size_t index = 0; index < tree.second.size(); ++index) {
                if (index != 0) output << ',';
                const auto& edge = tree.second[index];
                output << "{\"from\":" << edge.from << ",\"to\":" << edge.to
                       << ",\"length\":" << std::setprecision(17) << edge.length
                       << ",\"bootstrapSupport\":null,\"internal\":"
                       << (edge.from >= static_cast<int>(sequences.size()) ? "true" : "false")
                       << ",\"collapsed\":false}";
            }
        }
        output << "]}";
    }
    output << "]}";
    return output.str();
}

std::string event_phylpro_json(const WebContext& context, std::uint32_t event_id,
                               std::uint32_t requested_window, int gap_mode, int include_self) {
    if (event_id >= context.full.events.size()) return {};
    const auto& event = context.full.events[event_id];
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const int window = std::clamp(static_cast<int>(requested_window), 10, 5000);
    if (length < 2) return {};
    // The source calls VB CInt(window / 2), which is round-to-nearest-even.
    const int half_requested = std::max(1, (window / 2) +
        ((window & 1) != 0 && ((window / 2) & 1) != 0 ? 1 : 0));
    std::vector<int> context_sequences;
    for (std::size_t index = 0; index < context.alignment.sequences.size(); ++index) {
        if (index < context.disabled.size() && context.disabled[index]) continue;
        context_sequences.push_back(static_cast<int>(index));
    }
    const std::array<int, 3> targets = event_role_sequences(event);
    std::array<int, 3> sorted_targets = targets;
    std::sort(sorted_targets.begin(), sorted_targets.end());
    if (std::adjacent_find(sorted_targets.begin(), sorted_targets.end()) != sorted_targets.end()) return {};
    for (const int target : targets) {
        if (target < 0 || target >= static_cast<int>(context.alignment.sequences.size()) ||
            (target < static_cast<int>(context.disabled.size()) && context.disabled[target])) return {};
    }
    if (context_sequences.size() < 3) return {};
    const auto state = [](const char value) -> std::uint8_t {
        switch (static_cast<unsigned char>(std::toupper(static_cast<unsigned char>(value)))) {
        case 'A': return 1; case 'C': return 2; case 'G': return 3;
        case 'T': case 'U': return 4; default: return 0;
        }
    };
    std::vector<int> eligible_columns;
    eligible_columns.reserve(static_cast<std::size_t>(length));
    for (int position = 0; position < length; ++position) {
        std::array<bool, 5> seen{};
        bool missing = false;
        int observed = 0;
        for (const int sequence : context_sequences) {
            const auto code = state(context.alignment.sequences[static_cast<std::size_t>(sequence)][static_cast<std::size_t>(position)]);
            if (code == 0) { missing = true; continue; }
            if (!seen[code]) { seen[code] = true; ++observed; }
        }
        if ((gap_mode == 0 || !missing) && observed >= 2) eligible_columns.push_back(position);
    }
    if (eligible_columns.size() < 2) return {};
    const auto vb_round_half = [](const std::size_t value) {
        const std::size_t lower = value / 2;
        return lower + ((value & 1U) != 0U && (lower & 1U) != 0U ? 1U : 0U);
    };
    const std::size_t requested_half = static_cast<std::size_t>(half_requested);
    const std::size_t maximum_half = context.options.circular
        ? std::max<std::size_t>(1, vb_round_half(eligible_columns.size()))
        : std::max<std::size_t>(1, eligible_columns.size() / 2);
    const std::size_t half = std::min(requested_half, maximum_half);
    const auto mismatch = [&](const int first, const int second, const int position) {
        const auto left = state(context.alignment.sequences[static_cast<std::size_t>(first)][static_cast<std::size_t>(position)]);
        const auto right = state(context.alignment.sequences[static_cast<std::size_t>(second)][static_cast<std::size_t>(position)]);
        return left != 0 && right != 0 && left != right;
    };
    std::array<std::size_t, 3> target_context_index{};
    for (int role = 0; role < 3; ++role) {
        const auto found = std::find(context_sequences.begin(), context_sequences.end(), targets[role]);
        if (found == context_sequences.end()) return {};
        target_context_index[role] = static_cast<std::size_t>(found - context_sequences.begin());
    }
    std::array<std::vector<float>, 3> left, right;
    for (int role = 0; role < 3; ++role) {
        left[role].assign(context_sequences.size(), 0.0F);
        right[role].assign(context_sequences.size(), 0.0F);
    }
    const auto add_position = [&](const int target, const int position, const float direction,
                                  std::vector<float>& distances) {
        for (std::size_t index = 0; index < context_sequences.size(); ++index) {
            if (mismatch(target, context_sequences[index], position)) distances[index] += direction;
        }
    };
    const std::size_t column_count = eligible_columns.size();
    const auto populate = [&](const std::size_t left_begin, const std::size_t right_begin) {
        for (std::size_t offset = 0; offset < half; ++offset) {
            const int left_position = eligible_columns[left_begin + offset];
            const int right_position = eligible_columns[right_begin + offset];
            for (int role = 0; role < 3; ++role) {
                add_position(targets[role], left_position, 1.0F, left[role]);
                add_position(targets[role], right_position, 1.0F, right[role]);
            }
        }
    };
    if (context.options.circular) populate(column_count - half, 0);
    else populate(0, half);
    std::vector<int> coordinates;
    std::array<std::vector<double>, 3> profiles;
    const auto source_pearson = [&](const std::vector<float>& x_values,
                                    const std::vector<float>& y_values,
                                    const std::size_t self_index) {
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
        std::size_t observations = 0;
        for (std::size_t index = 0; index < x_values.size(); ++index) {
            if (!include_self && index == self_index) continue;
            const double x = x_values[index], y = y_values[index];
            sum_x += x; sum_y += y; sum_xy += x * y; sum_x2 += x * x; sum_y2 += y * y;
            ++observations;
        }
        if (observations == 0 || sum_x2 <= 0.0 || sum_y2 <= 0.0) return 1.0;
        const double count = static_cast<double>(observations);
        const double numerator = count * sum_xy - sum_x * sum_y;
        const double left_variance = count * sum_x2 - sum_x * sum_x;
        const double right_variance = count * sum_y2 - sum_y * sum_y;
        if (left_variance <= 0.0 || right_variance <= 0.0) return 1.0;
        const double denominator = std::sqrt(left_variance) * std::sqrt(right_variance);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) return 1.0;
        const double correlation = numerator / denominator;
        return std::isfinite(correlation) ? correlation : 1.0;
    };
    const auto emit_point = [&](const std::size_t column_index) {
        coordinates.push_back(eligible_columns[column_index] + 1);
        for (int role = 0; role < 3; ++role) {
            profiles[role].push_back(source_pearson(left[role], right[role], target_context_index[role]));
        }
    };
    std::size_t rolling_updates = 0;
    if (context.options.circular) {
        coordinates.reserve(column_count);
        for (std::size_t center = 0; center < column_count; ++center) {
            emit_point(center);
            if (center + 1 == column_count) break;
            const std::size_t old_left = (center + column_count - half) % column_count;
            const std::size_t boundary = center;
            const std::size_t new_right = (center + half) % column_count;
            for (int role = 0; role < 3; ++role) {
                add_position(targets[role], eligible_columns[old_left], -1.0F, left[role]);
                add_position(targets[role], eligible_columns[boundary], 1.0F, left[role]);
                add_position(targets[role], eligible_columns[boundary], -1.0F, right[role]);
                add_position(targets[role], eligible_columns[new_right], 1.0F, right[role]);
            }
            rolling_updates += 3 * context_sequences.size() * 4;
        }
    } else {
        const std::size_t first_center = half;
        const std::size_t last_center = column_count - half;
        if (last_center < first_center) return {};
        coordinates.reserve(last_center - first_center + 1);
        for (std::size_t center = first_center; center <= last_center; ++center) {
            emit_point(std::min(center, column_count - 1));
            if (center == last_center) break;
            const std::size_t old_left = center - half;
            const std::size_t boundary = center;
            const std::size_t new_right = center + half;
            for (int role = 0; role < 3; ++role) {
                add_position(targets[role], eligible_columns[old_left], -1.0F, left[role]);
                add_position(targets[role], eligible_columns[boundary], 1.0F, left[role]);
                add_position(targets[role], eligible_columns[boundary], -1.0F, right[role]);
                add_position(targets[role], eligible_columns[new_right], 1.0F, right[role]);
            }
            rolling_updates += 3 * context_sequences.size() * 4;
        }
    }
    if (coordinates.empty()) return {};
    std::array<std::size_t, 3> minimum_indices{};
    for (int target = 0; target < 3; ++target) {
        minimum_indices[target] = static_cast<std::size_t>(
            std::min_element(profiles[target].begin(), profiles[target].end()) - profiles[target].begin());
    }
    double minimum = 1.0;
    double maximum = -1.0;
    for (const auto& profile : profiles) {
        for (const double value : profile) {
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
    }
    const auto coordinate_distance = [&](const int first, const int second) {
        const int direct = std::abs(first - second);
        return context.options.circular ? std::min(direct, length - std::min(direct, length)) : direct;
    };
    const auto nearest = [&](int coordinate) {
        std::size_t best = 0;
        int best_distance = std::numeric_limits<int>::max();
        for (std::size_t index = 0; index < coordinates.size(); ++index) {
            const int difference = coordinate_distance(coordinates[index], coordinate);
            if (difference < best_distance) {
                best = index;
                best_distance = difference;
            }
        }
        return best;
    };
    constexpr std::size_t maximum_plot_points = 2048;
    const std::size_t sample_stride = std::max<std::size_t>(
        1, (coordinates.size() + maximum_plot_points - 1) / maximum_plot_points);
    std::vector<std::size_t> retained;
    retained.reserve(std::min(coordinates.size(), maximum_plot_points + 8));
    for (std::size_t index = 0; index < coordinates.size(); index += sample_stride) {
        retained.push_back(index);
    }
    retained.push_back(coordinates.size() - 1);
    retained.push_back(nearest(event.beginning));
    retained.push_back(nearest(event.ending));
    for (int target = 0; target < 3; ++target) retained.push_back(minimum_indices[target]);
    std::sort(retained.begin(), retained.end());
    retained.erase(std::unique(retained.begin(), retained.end()), retained.end());
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"eventId\":" << event_id << ",\"status\":\"source-shaped-active-unvalidated\",\"kernel\":\"FindSubSeqPP-MakePDstMat-UpdatePDstMat-PPRegression\""
           << ",\"roleOrder\":\"recombinant-major-parent-minor-parent\",\"columnSelection\":\"polymorphic-after-gap-policy\",\"distance\":\"pairwise-Hamming-count\",\"correlation\":\"Pearson-source-single-output\",\"significanceTest\":\"not-implemented-in-supplied-rdp5\",\"optimization\":\"three-target-rows-linear-in-context\""
           << ",\"maskedContextIncluded\":true,\"disabledContextExcluded\":true,\"fragmentContextIncluded\":false"
           << ",\"circular\":" << (context.options.circular ? "true" : "false") << ",\"windowSites\":" << window
           << ",\"halfWindowSites\":" << half << ",\"windowCapped\":"
           << (half != requested_half ? "true" : "false") << ",\"gapMode\":\""
           << (gap_mode ? "strip-any-missing-column" : "ignore-missing-pairwise") << "\",\"includeSelf\":" << (include_self ? "true" : "false")
           << ",\"eligibleColumns\":" << eligible_columns.size() << ",\"contextSequences\":" << context_sequences.size()
           << ",\"targetContextComparisons\":" << (3 * context_sequences.size() * 2 * half)
           << ",\"rollingUpdates\":" << rolling_updates << ",\"evaluatedPoints\":" << coordinates.size()
           << ",\"returnedPoints\":" << retained.size() << ",\"minimumValue\":" << minimum << ",\"maximumValue\":" << maximum
           << ",\"sequenceIndices\":[" << targets[0] << ',' << targets[1] << ',' << targets[2] << "],\"sequenceNames\":[";
    for (int target = 0; target < 3; ++target) {
        if (target != 0) output << ',';
        json_string(output, context.alignment.names[static_cast<std::size_t>(targets[target])]);
    }
    output << "],\"minimumBySequence\":[";
    for (int target = 0; target < 3; ++target) {
        if (target != 0) output << ',';
        output << "{\"sequenceIndex\":" << targets[target] << ",\"position\":" << coordinates[minimum_indices[target]]
               << ",\"correlation\":" << profiles[target][minimum_indices[target]] << '}';
    }
    output << "],\"breakpoints\":[";
    for (int boundary = 0; boundary < 2; ++boundary) {
        if (boundary != 0) output << ',';
        const std::size_t point_index = nearest(boundary == 0 ? event.beginning : event.ending);
        output << "{\"name\":\"" << (boundary == 0 ? "beginning" : "ending") << "\",\"eventPosition\":"
               << (boundary == 0 ? event.beginning : event.ending) << ",\"profilePosition\":" << coordinates[point_index] << ",\"correlations\":[";
        for (int target = 0; target < 3; ++target) {
            if (target != 0) output << ',';
            output << profiles[target][point_index];
        }
        output << "]}";
    }
    output << "],\"points\":[";
    for (std::size_t output_index = 0; output_index < retained.size(); ++output_index) {
        const std::size_t index = retained[output_index];
        if (output_index != 0) output << ',';
        output << "{\"alignmentPosition\":" << coordinates[index]
               << ",\"recombinant\":" << profiles[0][index]
               << ",\"majorParent\":" << profiles[1][index]
               << ",\"minorParent\":" << profiles[2][index] << '}';
    }
    output << "]}";
    return output.str();
}

void write_maxchi_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::MaxChiDiscoveryCandidate& discovery) {
    output << std::setprecision(17)
           << "{\"status\":\"source-shaped-active-unvalidated\""
           << ",\"kernel\":\"MCXoverF-multi-peak-destroy-retry\""
           << ",\"peakOrdering\":\"raw-chi-square-lazy-heap\""
           << ",\"smoothingUse\":\"twelve-term-eleven-divisor-source-basin-destruction-only\""
           << ",\"peakAttempt\":" << discovery.peak_attempt
           << ",\"peakPair\":";
    if (discovery.peak_pair < 0) output << "null";
    else output << static_cast<unsigned int>(discovery.peak_pair);
    output << ",\"tractSide\":\"";
    if (discovery.tract_side == next_rdp_legacy_optional::MaxChiTractSide::left) output << "left";
    else if (discovery.tract_side == next_rdp_legacy_optional::MaxChiTractSide::right) output << "right";
    else output << "unavailable";
    output << "\",\"peakAlignmentPosition\":" << discovery.peak_alignment_position
           << ",\"variableSites\":" << discovery.variable_sites
           << ",\"initialHalfWindow\":" << discovery.initial_half_window
           << ",\"grownHalfWindow\":" << discovery.grown_half_window
           << ",\"criticalDifference\":" << discovery.critical_difference
           << ",\"maximumChiSquare\":" << discovery.maximum_chi_square
           << ",\"rawPValue\":" << discovery.raw_p_value
           << ",\"withinTripletPValue\":" << discovery.within_triplet_p_value
           << ",\"correctedPValue\":" << discovery.corrected_p_value
           << ",\"leftFlankChiSquare\":" << discovery.left_flank_chi_square
           << ",\"rightFlankChiSquare\":" << discovery.right_flank_chi_square
           << ",\"missingDataWindowFilterApplied\":"
           << (discovery.missing_data_window_filter_applied ? "true" : "false")
           << ",\"linearEdgeWindowFilterApplied\":"
           << (discovery.linear_edge_window_filter_applied ? "true" : "false") << '}';
}

void write_chimaera_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::ChimaeraDiscoveryCandidate& discovery) {
    output << std::setprecision(17)
           << "{\"status\":\"source-shaped-active-unvalidated\""
           << ",\"kernel\":\"AlistChi-FastRecCheckChim-CXoverA\""
           << ",\"profile\":\"target-specific-information-rich-binary-string\""
           << ",\"peakOrdering\":\"raw-chi-square-lazy-heap-per-target\""
           << ",\"smoothingUse\":\"twelve-term-eleven-divisor-source-basin-destruction-only\""
           << ",\"targetLocal\":" << static_cast<unsigned int>(discovery.target_local)
           << ",\"peakAttempt\":" << discovery.peak_attempt
           << ",\"tractSide\":\"";
    if (discovery.tract_side == next_rdp_legacy_optional::MaxChiTractSide::left) output << "left";
    else if (discovery.tract_side == next_rdp_legacy_optional::MaxChiTractSide::right) output << "right";
    else output << "unavailable";
    output << "\",\"peakAlignmentPosition\":" << discovery.peak_alignment_position
           << ",\"informationRichSites\":" << discovery.information_rich_sites
           << ",\"initialHalfWindow\":" << discovery.initial_half_window
           << ",\"grownHalfWindow\":" << discovery.grown_half_window
           << ",\"criticalDifference\":" << discovery.critical_difference
           << ",\"maximumChiSquare\":" << discovery.maximum_chi_square
           << ",\"rawPValue\":" << discovery.raw_p_value
           << ",\"withinTripletPValue\":" << discovery.within_triplet_p_value
           << ",\"correctedPValue\":" << discovery.corrected_p_value
           << ",\"leftFlankChiSquare\":" << discovery.left_flank_chi_square
           << ",\"rightFlankChiSquare\":" << discovery.right_flank_chi_square
           << ",\"insideParentOneMatchRate\":" << discovery.inside_parent_one_match_rate
           << ",\"outsideParentOneMatchRate\":" << discovery.outside_parent_one_match_rate
           << ",\"missingDataWindowFilterApplied\":"
           << (discovery.missing_data_window_filter_applied ? "true" : "false")
           << ",\"linearEdgeWindowFilterApplied\":"
           << (discovery.linear_edge_window_filter_applied ? "true" : "false") << '}';
}

void write_geneconv_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::GeneconvDiscoveryCandidate& discovery) {
    output << std::setprecision(17)
           << "{\"status\":\"source-shaped-active-unvalidated\""
           << ",\"kernel\":\"FindSubSeqGCAP6-GetFragsP-GetMaxFragScoreP-CalcKMaxP-GCCalcPValP2-GCXoverD\""
           << ",\"probabilityModel\":\"karlin-altschul\""
           << ",\"indelMode\":\"ignored\""
           << ",\"overlapPolicy\":\"stable-lowest-p-configured-coverage\""
           << ",\"minimumFragmentFiltersApplied\":false"
           << ",\"track\":" << static_cast<unsigned int>(discovery.track)
           << ",\"polymorphicSites\":" << discovery.polymorphic_sites
           << ",\"positiveSites\":" << discovery.positive_sites
           << ",\"discordantSites\":" << discovery.discordant_sites
           << ",\"mismatchPenalty\":" << discovery.mismatch_penalty
           << ",\"fragmentScore\":" << discovery.fragment_score
           << ",\"criticalScore\":" << discovery.critical_score
           << ",\"lambda\":" << discovery.lambda
           << ",\"karlinAltschulK\":" << discovery.karlin_altschul_k
           << ",\"rawPValue\":" << discovery.raw_p_value
           << ",\"correctedPValue\":" << discovery.corrected_p_value
           << ",\"karlinAltschulProbability\":"
           << (discovery.karlin_altschul_probability ? "true" : "false")
           << ",\"ignoredIndels\":" << (discovery.ignored_indels ? "true" : "false")
           << ",\"overlapFilterApplied\":"
           << (discovery.overlap_filter_applied ? "true" : "false") << '}';
}

void write_threeseq_discovery_json(
    std::ostringstream& output,
    const next_rdp_legacy_optional::ThreeSeqDiscoveryCandidate& discovery) {
    output << std::setprecision(17)
           << "{\"status\":\"source-shaped-active-unvalidated\""
           << ",\"kernel\":\"FindSubSeqTS-Seq3PVals-CheckwrapC-TSXOver\""
           << ",\"profile\":\"target-specific-information-rich-random-walk\""
           << ",\"probabilityModel\":\"exact-hypergeometric-walk-with-siegmund-fallback\""
           << ",\"correctionModel\":\"dunn-sidak-when-project-correction-enabled\""
           << ",\"targetLocal\":" << static_cast<unsigned int>(discovery.target_local)
           << ",\"walkDirection\":\""
           << (discovery.direction == next_rdp_legacy_optional::ThreeSeqWalkDirection::ascent ? "ascent" : "descent")
           << "\",\"informationRichSites\":" << discovery.information_rich_sites
           << ",\"parentOneMatches\":" << discovery.parent_one_matches
           << ",\"parentTwoMatches\":" << discovery.parent_two_matches
           << ",\"probabilityExcursion\":" << discovery.probability_excursion
           << ",\"maximumExcursion\":" << discovery.maximum_excursion
           << ",\"rawPValue\":" << discovery.raw_p_value
           << ",\"correctedPValue\":" << discovery.corrected_p_value
           << ",\"exactProbability\":" << (discovery.exact_probability ? "true" : "false")
           << ",\"siegmundFallback\":" << (discovery.siegmund_fallback ? "true" : "false")
           << ",\"missingDataSplitApplied\":"
           << (discovery.missing_data_split_applied ? "true" : "false") << '}';
}

std::string full_json(const WebContext& context) {
    const auto& result = context.full;
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"sourceFaithfulCore\":true,\"engineVersion\":\"nextRDP-core 0.1.0\""
           << ",\"sequenceCount\":" << result.sequence_count
           << ",\"sequenceLength\":" << result.sequence_length
           << ",\"tripletCount\":" << result.triplet_count
           << ",\"rawCandidateCount\":" << result.raw_candidate_count
           << ",\"rdpEnabled\":" << (context.options.enable_rdp ? "true" : "false")
           << ",\"analysisMode\":\""
           << (context.options.query_reference_mode ? "query-reference" : "exploratory")
           << "\",\"queryReference\":{\"active\":"
           << (context.options.query_reference_mode ? "true" : "false")
           << ",\"querySequenceCount\":";
    std::size_t query_count = 0;
    std::size_t reference_count = 0;
    std::size_t reference_group_count = 0;
    if (context.options.query_reference_mode) {
        std::vector<unsigned int> groups = context.reference_groups;
        std::sort(groups.begin(), groups.end());
        groups.erase(std::remove(groups.begin(), groups.end(), 0U), groups.end());
        groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
        reference_group_count = groups.size();
        for (std::size_t index = 0; index < context.alignment.names.size(); ++index) {
            if (index < context.masked.size() && context.masked[index]) continue;
            if (index < context.disabled.size() && context.disabled[index]) continue;
            if (index < context.reference_groups.size() && context.reference_groups[index] == 0) ++query_count;
            else if (index < context.reference_groups.size() && context.reference_groups[index] > 0) ++reference_count;
        }
    }
    output << query_count << ",\"referenceSequenceCount\":" << reference_count
           << ",\"referenceGroupCount\":" << reference_group_count
           << ",\"tripletConstraint\":\"one-query-two-different-reference-groups\""
           << ",\"sourceCorrectionRule\":\"reference-group-pairs-times-query-origins\"}"
           << ",\"geneconvEnabled\":" << (context.options.enable_geneconv ? "true" : "false")
           << ",\"maxChiEnabled\":" << (context.options.enable_maxchi ? "true" : "false")
           << ",\"chimaeraEnabled\":" << (context.options.enable_chimaera ? "true" : "false")
           << ",\"threeSeqEnabled\":" << (context.options.enable_three_seq ? "true" : "false")
           << ",\"bootscanPrimaryEnabled\":" << (context.options.enable_bootscan ? "true" : "false")
           << ",\"bootscanSecondaryEnabled\":" << (context.options.enable_bootscan_secondary ? "true" : "false")
           << ",\"siscanPrimaryEnabled\":" << (context.options.enable_siscan ? "true" : "false")
           << ",\"siscanSecondaryEnabled\":" << (context.options.enable_siscan_secondary ? "true" : "false")
           << ",\"enabledMethods\":[";
    bool first_method = true;
    const auto append_method = [&](const char* name, bool enabled) {
        if (!enabled) return;
        if (!first_method) output << ',';
        first_method = false;
        output << '"' << name << '"';
    };
    append_method("RDP", context.options.enable_rdp);
    append_method("GENECONV", context.options.enable_geneconv);
    append_method("MAXCHI", context.options.enable_maxchi);
    append_method("CHIMAERA", context.options.enable_chimaera);
    append_method("3SEQ", context.options.enable_three_seq);
    append_method("BOOTSCAN", context.options.enable_bootscan || context.options.enable_bootscan_secondary);
    append_method("SISCAN", context.options.enable_siscan || context.options.enable_siscan_secondary);
    output << ']'
           << ",\"events\":[";
    for (std::size_t index = 0; index < result.events.size(); ++index) {
        if (index != 0) output << ',';
        const auto& event = result.events[index];
        output << "{\"id\":" << index
               << ",\"program\":" << event.program
               << ",\"winningRole\":" << event.winning_role
               << ",\"methodTargetRole\":";
        if (event.method_target_role < 0) output << "null";
        else output << event.method_target_role;
        output
               << ",\"probability\":" << event.probability
               << ",\"beginning\":" << event.beginning
               << ",\"ending\":" << event.ending << ",\"representativeSequences\":["
               << event.representative_sequences[0] << ','
               << event.representative_sequences[1] << ','
               << event.representative_sequences[2] << "],\"profileSequences\":";
        if (!event.profile_sequences_available) {
            output << "null";
        } else {
            output << '[' << event.profile_sequences[0] << ','
                   << event.profile_sequences[1] << ','
                   << event.profile_sequences[2] << ']';
        }
        output << ",\"sequenceGroups\":[";
        for (int role = 0; role < 3; ++role) {
            if (role != 0) output << ',';
            output << '[';
            for (std::size_t member = 0; member < event.sequence_groups[role].size(); ++member) {
                if (member != 0) output << ',';
                output << event.sequence_groups[role][member];
            }
            output << ']';
        }
        output << "]"
               << ",\"burtAttempted\":" << (event.burt_attempted ? "true" : "false")
               << ",\"burtApplied\":" << (event.burt_applied ? "true" : "false")
               << ",\"burtInformationRichSites\":" << event.burt.information_rich_sites
               << ",\"burtCandidateIntervalCount\":" << event.burt.intervals.size()
               << ",\"burtBestLogLikelihood\":" << event.burt.best_log_likelihood
               << ",\"burtInputBeginning\":" << event.burt.input_beginning
               << ",\"burtInputEnding\":" << event.burt.input_ending
               << ",\"burtConfidence\":[";
        for (std::size_t confidence = 0; confidence < event.burt.confidence.size(); ++confidence) {
            if (confidence != 0) output << ',';
            output << event.burt.confidence[confidence];
        }
        output << "]"
               << ",\"maxChiDiscovery\":";
        if (!event.maxchi_available) {
            output << "null";
        } else {
            write_maxchi_discovery_json(output, event.maxchi_discovery);
        }
        output << ",\"chimaeraDiscovery\":";
        if (!event.chimaera_available) {
            output << "null";
        } else {
            write_chimaera_discovery_json(output, event.chimaera_discovery);
        }
        output << ",\"geneconvDiscovery\":";
        if (!event.geneconv_available) {
            output << "null";
        } else {
            write_geneconv_discovery_json(output, event.geneconv_discovery);
        }
        output << ",\"threeSeqDiscovery\":";
        if (!event.three_seq_available) {
            output << "null";
        } else {
            write_threeseq_discovery_json(output, event.three_seq_discovery);
        }
        output << ",\"bootscanDiscovery\":";
        if (!event.bootscan_available) {
            output << "null";
        } else {
            const auto& discovery = event.bootscan_discovery;
            output << "{\"status\":\"source-shaped-active-unvalidated\",\"kernel\":\"BSXoverR-SEQBOOT2-FastBootDist-GetPltVal-ScanBSPlots-MakeBSEvent\",\"mode\":\"jukes-cantor-distance\",\"probabilityModel\":\"MakeScoresBS-binomial\",\"strictClosestPairVoting\":true"
                   << ",\"supportedPair\":" << static_cast<int>(discovery.supported_pair)
                   << ",\"windowsScored\":" << discovery.windows_scored
                   << ",\"usableWindows\":" << discovery.usable_windows
                   << ",\"informativeSites\":" << discovery.informative_sites
                   << ",\"tractInformativeSites\":" << discovery.tract_informative_sites
                   << ",\"tractPairMatches\":" << discovery.tract_pair_matches
                   << ",\"outsidePairMatches\":" << discovery.outside_pair_matches
                   << ",\"maximumPairSupport\":" << discovery.maximum_pair_support
                   << ",\"meanPairSupport\":" << discovery.mean_pair_support
                   << ",\"bootstrapPValue\":" << discovery.bootstrap_p_value
                   << ",\"rawPValue\":" << discovery.raw_p_value
                   << ",\"correctedPValue\":" << discovery.corrected_p_value
                   << ",\"erasedWindowFilterApplied\":"
                   << (discovery.erased_window_filter_applied ? "true" : "false")
                   << '}';
        }
        output << ",\"siscanDiscovery\":";
        if (!event.siscan_available) {
            output << "null";
        } else {
            const auto& discovery = event.siscan_discovery;
            output << "{\"status\":\"source-shaped-active-unvalidated\",\"kernel\":\"SSXoverC-GetSSOL-Get3Score-GetPScores2-DoPerms3-MakeZValue2-DoSums-FindMaxZ-ShrinkRegionC\",\"outlierMode\":\"nearest-source-wpgma\",\"permutationGenerator\":\"microsoft-crt-flat-prefix\",\"gapMode\":\"strip\",\"variablePatternMode\":\"one-two-three-variable\",\"sourceFastWindowQuirk\":true"
                   << ",\"globalPair\":" << static_cast<int>(discovery.global_pair)
                   << ",\"candidatePair\":" << static_cast<int>(discovery.candidate_pair)
                   << ",\"outlierSequence\":" << discovery.outlier_sequence
                   << ",\"windowsInRegion\":" << discovery.windows_in_region
                   << ",\"informativeSites\":" << discovery.informative_sites
                   << ",\"permutationDraws\":" << discovery.permutation_draws
                   << ",\"selectedScore\":" << static_cast<int>(discovery.selected_score)
                   << ",\"selectedScoreFamily\":\""
                   << (discovery.selected_score_family == next_rdp_legacy_optional::SiscanScoreFamily::partition
                       ? "partition" : discovery.selected_score_family == next_rdp_legacy_optional::SiscanScoreFamily::summed
                       ? "summed" : "unavailable")
                   << "\",\"maximumZ\":" << discovery.maximum_z
                   << ",\"normalTailPValue\":" << discovery.normal_tail_p_value
                   << ",\"regionLengthAdjustedPValue\":" << discovery.region_length_adjusted_p_value
                   << ",\"windowAdjustedPValue\":" << discovery.window_adjusted_p_value
                   << ",\"correctedPValue\":" << discovery.corrected_p_value
                   << '}';
        }
        output << ",\"bootscanRecheck\":";
        if (!event.bootscan_recheck.requested) {
            output << "null";
        } else {
            const auto& recheck = event.bootscan_recheck;
            output << "{\"status\":\"complete-active-unvalidated\",\"kernel\":\"BSXoverM-SEQBOOT2-FastBootDistIP-DrawBSPlotsIII\",\"eventDiscoveryApplied\":false,\"coordinateChanging\":false"
                   << ",\"requested\":true,\"representativeSkipped\":" << (recheck.representative_skipped ? "true" : "false")
                   << ",\"profileAvailable\":" << (recheck.profile_available ? "true" : "false")
                   << ",\"sourceDistanceMode\":true,\"sourceBinomialProbability\":true,\"sourceCircularWindows\":true"
                   << ",\"erasedWindowFilterApplied\":" << (recheck.erased_window_filter_applied ? "true" : "false")
                   << ",\"bonferroniApplied\":" << (recheck.bonferroni_applied ? "true" : "false")
                   << ",\"correctionTests\":" << recheck.correction_tests
                   << ",\"windowSites\":" << recheck.window_sites << ",\"stepSites\":" << recheck.step_sites
                   << ",\"bootstrapReplicates\":" << recheck.bootstrap_replicates << ",\"randomSeed\":" << recheck.random_seed
                   << ",\"supportCutoff\":" << recheck.support_cutoff
                   << ",\"windowsScanned\":" << recheck.windows_scanned << ",\"eventWindowsScored\":" << recheck.event_windows_scored
                   << ",\"usableEventWindows\":" << recheck.usable_event_windows << ",\"informativeSites\":" << recheck.informative_sites
                   << ",\"tractInformativeSites\":" << recheck.tract_informative_sites << ",\"tractPairMatches\":" << recheck.tract_pair_matches
                   << ",\"outsidePairMatches\":" << recheck.outside_pair_matches
                   << ",\"scoredPair\":" << static_cast<int>(recheck.scored_pair)
                   << ",\"maximumPairSupport\":" << recheck.maximum_pair_support
                   << ",\"meanScoredPairSupport\":" << recheck.mean_scored_pair_support
                   << ",\"localPValue\":" << recheck.local_p_value << ",\"correctedPValue\":" << recheck.corrected_p_value
                   << ",\"supportGatePassed\":" << (recheck.support_gate_passed ? "true" : "false")
                   << ",\"sourceRecheckHit\":" << (recheck.source_recheck_hit ? "true" : "false") << '}';
        }
        output << ",\"siscanRecheck\":";
        if (!event.siscan_recheck.requested) {
            output << "null";
        } else {
            const auto& recheck = event.siscan_recheck;
            output << "{\"status\":\"complete-active-unvalidated\",\"kernel\":\"GetSSOL-Get3Score-GetPScores2-DoPerms3P-MakeZValue2-DoSums\",\"eventDiscoveryApplied\":false,\"coordinateChanging\":false"
                   << ",\"requested\":true,\"representativeSkipped\":" << (recheck.representative_skipped ? "true" : "false")
                   << ",\"profileAvailable\":" << (recheck.profile_available ? "true" : "false")
                   << ",\"outlierAvailable\":" << (recheck.outlier_available ? "true" : "false")
                   << ",\"sourceNearestOutlier\":true,\"sourceGapStripping\":true,\"sourceVariablePatterns\":true"
                   << ",\"bonferroniApplied\":" << (recheck.bonferroni_applied ? "true" : "false")
                   << ",\"correctionTests\":" << recheck.correction_tests
                   << ",\"outlierSequence\":" << recheck.outlier_sequence
                   << ",\"informativeSites\":" << recheck.informative_sites << ",\"permutationDraws\":" << recheck.permutation_draws
                   << ",\"globalPair\":" << static_cast<int>(recheck.global_pair)
                   << ",\"scoredPair\":" << static_cast<int>(recheck.scored_pair)
                   << ",\"selectedScore\":" << static_cast<int>(recheck.selected_score)
                   << ",\"selectedScoreFamily\":\""
                   << (recheck.selected_score_family == next_rdp_legacy_optional::SiscanScoreFamily::partition
                       ? "partition" : recheck.selected_score_family == next_rdp_legacy_optional::SiscanScoreFamily::summed
                       ? "summed" : "unavailable")
                   << "\",\"maximumZ\":" << recheck.maximum_z
                   << ",\"normalTailPValue\":" << recheck.normal_tail_p_value
                   << ",\"regionLengthAdjustedPValue\":" << recheck.region_length_adjusted_p_value
                   << ",\"windowAdjustedPValue\":" << recheck.window_adjusted_p_value
                   << ",\"correctedPValue\":" << recheck.corrected_p_value
                   << ",\"sourceRecheckHit\":" << (recheck.source_recheck_hit ? "true" : "false") << '}';
        }
        output << "}";
    }
    output << "]}";
    return output.str();
}

std::string serialize(const RdpInitialAnalysisResult& result) {
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"sequenceCount\":" << result.alignment.next_no + 1
           << ",\"sequenceLength\":" << result.alignment.sequence_length
           << ",\"tripletCount\":"
           << result.alignment.analysis_list_last + 1
           << ",\"scannedTriplets\":" << result.events.scanned_triplets
           << ",\"significantIntervals\":"
           << result.events.significant_candidates << ",\"events\":[";
    bool first = true;
    for (int sequence = 0; sequence <= result.alignment.next_no; ++sequence) {
        const int count = result.events.current_xover[sequence];
        for (int slot = 0; slot < count; ++slot) {
            const auto& event = result.events.xover_list[sequence][slot];
            if (!first) output << ',';
            first = false;
            output << "{\"row\":" << sequence
                   << ",\"slot\":" << slot + 1
                   << ",\"daughter\":" << event.daughter
                   << ",\"minorParent\":" << event.minor_parent
                   << ",\"majorParent\":" << event.major_parent
                   << ",\"beginning\":" << event.beginning
                   << ",\"ending\":" << event.ending
                   << ",\"probability\":" << event.probability
                   << ",\"program\":"
                   << static_cast<int>(event.program_flag) << '}';
        }
    }
    output << "]}";
    return output.str();
}

}  // namespace

extern "C" {

NEXT_RDP_KEEPALIVE int next_rdp_analyze(
    const char* fasta, const std::size_t length, const int circular,
    const double p_value_cutoff, const int window_sites) {
    try {
        if (fasta == nullptr || length == 0) {
            throw std::runtime_error("No FASTA alignment was supplied");
        }
        const RdpInitialAnalysisOptions options{
            circular != 0,
            p_value_cutoff,
            window_sites,
        };
        result_json = serialize(run_rdp_initial_analysis_from_fasta_text(
            std::string(fasta, length), options));
        error_message.clear();
        return 1;
    } catch (const std::exception& error) {
        result_json.clear();
        error_message = error.what();
        return 0;
    }
}

NEXT_RDP_KEEPALIVE const char* next_rdp_result_json() {
    return result_json.c_str();
}

NEXT_RDP_KEEPALIVE const char* next_rdp_error() {
    return error_message.c_str();
}

NEXT_RDP_KEEPALIVE const char* next_rdp_version() {
    return "nextRDP-core 0.1.0";
}

NEXT_RDP_KEEPALIVE std::uint32_t rdp_create() {
    auto& contexts = web_contexts();
    for (std::size_t index = 0; index < contexts.size(); ++index) {
        if (!contexts[index]) {
            contexts[index] = std::make_unique<WebContext>();
            return static_cast<std::uint32_t>(index + 1);
        }
    }
    contexts.push_back(std::make_unique<WebContext>());
    return static_cast<std::uint32_t>(contexts.size());
}

NEXT_RDP_KEEPALIVE void rdp_destroy(const std::uint32_t handle) {
    if (handle != 0 && handle <= web_contexts().size()) web_contexts()[handle - 1].reset();
}

NEXT_RDP_KEEPALIVE const char* rdp_version() {
    return "nextRDP-core 0.1.0";
}

NEXT_RDP_KEEPALIVE std::uint32_t rdp_set_worker_threads(
    const std::uint32_t handle, const std::uint32_t requested) {
    if (web_context(handle) == nullptr) return 0;
    return static_cast<std::uint32_t>(set_rdp_method_worker_threads(
        static_cast<int>(requested == 0 ? 1 : requested)));
}

NEXT_RDP_KEEPALIVE int rdp_load_alignment(
    const std::uint32_t handle, const std::uint8_t* bytes, const std::size_t length) {
    auto* context = web_context(handle);
    if (context == nullptr) return 0;
    context->error.clear();
    try {
        if (bytes == nullptr || length == 0) throw std::runtime_error("The selected alignment file is empty.");
        const std::string input(reinterpret_cast<const char*>(bytes), length);
        context->alignment = parse_alignment(input);
        // The source-faithful command-line core consumes FASTA. Preserve the
        // detected input format for the browser/project metadata while
        // canonicalising the same aligned records at this boundary.
        context->fasta = canonical_fasta(context->alignment);
        context->review_states.assign(context->alignment.names.size(), 0);
        context->loaded = true;
        context->started = false;
        context->finished = false;
        context->full = {};
        return 1;
    } catch (const std::exception& error) {
        context->error = error.what();
        context->loaded = false;
        return 0;
    }
}

NEXT_RDP_KEEPALIVE const char* rdp_get_summary_json(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    if (!context->loaded) {
        context->error = "No alignment is loaded.";
        return "";
    }
    return cached(*context, summary_json(*context));
}

NEXT_RDP_KEEPALIVE int rdp_scan_begin(
    const std::uint32_t handle, const int circular, const int correction_mode,
    const double p_value_cutoff, const std::uint32_t window_sites,
    const int rdp_enabled,
    const int maxchi_enabled, const std::uint32_t maxchi_window_sites,
    const int chimaera_enabled, const std::uint32_t chimaera_window_sites,
    const int geneconv_enabled, const std::uint32_t geneconv_mismatch_scale,
    const std::uint32_t geneconv_max_overlaps, const int threeseq_enabled,
    const int bootscan_primary_enabled, const int bootscan_secondary_enabled,
    const std::uint32_t bootscan_window_sites, const std::uint32_t bootscan_step_sites,
    const std::uint32_t bootscan_bootstrap_replicates, const double bootscan_support_cutoff,
    const std::uint32_t bootscan_random_seed, const int siscan_primary_enabled,
    const int siscan_secondary_enabled, const std::uint32_t siscan_window_sites,
    const std::uint32_t siscan_step_sites, const std::uint32_t siscan_scan_permutations,
    const std::uint32_t siscan_p_value_permutations, const std::uint32_t siscan_random_seed,
    const int polish_breakpoints, const int query_reference_mode,
    const std::uint32_t* reference_groups, const std::size_t reference_group_count,
    const std::uint8_t* masked_sequences, const std::size_t mask_length,
    const std::uint8_t* disabled_sequences, const std::size_t disabled_length) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->loaded) return 0;
    const auto count = context->alignment.names.size();
    if (masked_sequences == nullptr || disabled_sequences == nullptr ||
        mask_length != count || disabled_length != count) {
        context->error = "The scan's sequence-role buffers do not match the loaded alignment.";
        return 0;
    }
    RdpInitialAnalysisOptions options;
    options.circular = circular != 0;
    options.correction_bonferroni = correction_mode == 0;
    options.p_value_cutoff = p_value_cutoff;
    options.window_sites = static_cast<int>(window_sites);
    options.enable_rdp = rdp_enabled != 0;
    options.enable_geneconv = geneconv_enabled != 0;
    options.geneconv_mismatch_scale = static_cast<int>(geneconv_mismatch_scale);
    options.geneconv_max_overlaps = static_cast<int>(geneconv_max_overlaps);
    options.enable_maxchi = maxchi_enabled != 0;
    options.maxchi_window_sites = static_cast<int>(maxchi_window_sites);
    options.enable_chimaera = chimaera_enabled != 0;
    options.chimaera_window_sites = static_cast<int>(chimaera_window_sites);
    options.enable_three_seq = threeseq_enabled != 0;
    options.enable_bootscan = bootscan_primary_enabled != 0;
    options.enable_bootscan_secondary = bootscan_secondary_enabled != 0;
    options.bootscan_window_sites = static_cast<int>(bootscan_window_sites);
    options.bootscan_step_sites = static_cast<int>(bootscan_step_sites);
    options.bootscan_bootstrap_replicates = static_cast<int>(bootscan_bootstrap_replicates);
    options.bootscan_support_cutoff = bootscan_support_cutoff;
    options.bootscan_random_seed = bootscan_random_seed;
    options.enable_siscan = siscan_primary_enabled != 0;
    options.enable_siscan_secondary = siscan_secondary_enabled != 0;
    options.siscan_window_sites = static_cast<int>(siscan_window_sites);
    options.siscan_step_sites = static_cast<int>(siscan_step_sites);
    options.siscan_scan_permutations = static_cast<int>(siscan_scan_permutations);
    options.siscan_p_value_permutations = static_cast<int>(siscan_p_value_permutations);
    options.siscan_random_seed = siscan_random_seed;
    options.polish_breakpoints_with_burt = polish_breakpoints != 0;
    options.query_reference_mode = query_reference_mode != 0;
    if (options.query_reference_mode) {
        if (reference_groups == nullptr || reference_group_count != count) {
            context->error = "The query-vs-reference group buffer does not match the loaded alignment.";
            return 0;
        }
        options.reference_groups.assign(reference_groups, reference_groups + reference_group_count);
    }
    if (!options.enable_rdp && !options.enable_geneconv &&
        !options.enable_maxchi && !options.enable_chimaera &&
        !options.enable_three_seq && !options.enable_bootscan &&
        !options.enable_siscan) {
        context->error = "Select RDP or at least one other discovery method.";
        return 0;
    }
    options.progress_callback = &web_progress_callback;
    options.progress_user = context;
    options.cancellation_callback = &web_cancellation_callback;
    options.cancellation_user = context;
    context->masked.assign(masked_sequences, masked_sequences + mask_length);
    context->disabled.assign(disabled_sequences, disabled_sequences + disabled_length);
    options.masked_sequences = context->masked;
    options.disabled_sequences = context->disabled;
    context->reference_groups = options.reference_groups;
    context->options = std::move(options);
    context->cancel_requested.store(0, std::memory_order_relaxed);
    context->started = true;
    context->finished = false;
    context->progress_phase = 0;
    context->progress_round = 1;
    context->progress_processed_triplets = 0;
    context->progress_total_triplets = count < 3
        ? 0 : static_cast<int>(count * (count - 1) * (count - 2) / 6);
    if (context->options.query_reference_mode) {
        std::size_t references = 0;
        std::size_t queries = 0;
        for (std::size_t index = 0; index < count; ++index) {
            if (context->masked[index] || context->disabled[index]) continue;
            if (context->reference_groups[index] == 0) ++queries;
            else ++references;
        }
        std::size_t pairs = 0;
        for (std::size_t first = 0; first < count; ++first) {
            if (context->masked[first] || context->disabled[first] ||
                context->reference_groups[first] == 0) continue;
            for (std::size_t second = first + 1; second < count; ++second) {
                if (context->masked[second] || context->disabled[second] ||
                    context->reference_groups[second] == 0 ||
                    context->reference_groups[first] == context->reference_groups[second]) continue;
                ++pairs;
            }
        }
        (void)references;
        context->progress_total_triplets = static_cast<int>(
            std::min<std::size_t>(std::numeric_limits<int>::max(), pairs * queries));
    }
    context->progress_event_count = 0;
    context->error.clear();
    return 1;
}

NEXT_RDP_KEEPALIVE int rdp_scan_batch(const std::uint32_t handle, const std::uint32_t /*triplet_budget*/) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->loaded || !context->started) return -1;
    if (context->finished) return 1;
    try {
        context->full = run_rdp_full_analysis_from_fasta_text(context->fasta, context->options);
        context->finished = true;
        return 1;
    } catch (const RdpAnalysisCancelled& error) {
        context->error = error.what();
        context->started = false;
        context->finished = false;
        return 2;
    } catch (const std::exception& error) {
        context->error = error.what();
        return -1;
    }
}

NEXT_RDP_KEEPALIVE int rdp_reconcile(const std::uint32_t handle) {
    auto* context = web_context(handle);
    return context != nullptr && context->finished ? 1 : 0;
}

NEXT_RDP_KEEPALIVE void rdp_cancel(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context != nullptr) {
        context->cancel_requested.store(1, std::memory_order_relaxed);
    }
}

NEXT_RDP_KEEPALIVE std::uintptr_t rdp_get_cancel_flag_address(
    const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr) return 0;
    static_assert(std::atomic<std::int32_t>::is_always_lock_free);
    return reinterpret_cast<std::uintptr_t>(&context->cancel_requested);
}

NEXT_RDP_KEEPALIVE const char* rdp_get_progress_json(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto count = context->alignment.names.size();
    // The query/reference scheme replaces MakeAListP2's exploratory total
    // with reference-group pairs × query origins.  Recomputing choose-three
    // here made the live monitor report the wrong denominator even though
    // rdp_scan_begin had already calculated the constrained work list.
    std::size_t total = count < 3 ? 0 : count * (count - 1) * (count - 2) / 6;
    if (context->options.query_reference_mode) {
        std::size_t queries = 0;
        std::size_t pairs = 0;
        for (std::size_t index = 0; index < count; ++index) {
            if (index < context->masked.size() && context->masked[index]) continue;
            if (index < context->disabled.size() && context->disabled[index]) continue;
            if (index < context->reference_groups.size() &&
                context->reference_groups[index] == 0) {
                ++queries;
            }
        }
        for (std::size_t first = 0; first < count; ++first) {
            if (first >= context->reference_groups.size() ||
                context->reference_groups[first] == 0 ||
                (first < context->masked.size() && context->masked[first]) ||
                (first < context->disabled.size() && context->disabled[first])) {
                continue;
            }
            for (std::size_t second = first + 1; second < count; ++second) {
                if (second >= context->reference_groups.size() ||
                    context->reference_groups[second] == 0 ||
                    context->reference_groups[first] == context->reference_groups[second] ||
                    (second < context->masked.size() && context->masked[second]) ||
                    (second < context->disabled.size() && context->disabled[second])) {
                    continue;
                }
                ++pairs;
            }
        }
        total = pairs * queries;
    }
    const int phase = context->finished ? 2 : context->progress_phase;
    const char* phase_name = phase == 1 ? "cyclic-rescan"
        : phase == 2 ? "complete" : "primary";
    const auto processed = context->finished
        ? static_cast<std::size_t>(total)
        : static_cast<std::size_t>(std::max(0, context->progress_processed_triplets));
    const auto event_count = context->finished
        ? context->full.events.size()
        : static_cast<std::size_t>(std::max(0, context->progress_event_count));
    std::size_t query_count = 0;
    std::size_t reference_count = 0;
    std::size_t reference_group_count = 0;
    if (context->options.query_reference_mode) {
        std::vector<unsigned int> groups;
        for (std::size_t index = 0; index < count; ++index) {
            if (index < context->masked.size() && context->masked[index]) continue;
            if (index < context->disabled.size() && context->disabled[index]) continue;
            const unsigned int group = index < context->reference_groups.size()
                ? context->reference_groups[index] : 0U;
            if (group == 0U) ++query_count;
            else {
                ++reference_count;
                groups.push_back(group);
            }
        }
        std::sort(groups.begin(), groups.end());
        groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
        reference_group_count = groups.size();
    }
    std::ostringstream output;
    output << "{\"state\":\"" << (context->finished ? "done" : context->started ? "running" : "idle")
           << "\",\"phase\":\"" << phase_name
           << "\",\"processedTriplets\":" << processed
           << ",\"totalTriplets\":" << total
           << ",\"correctionTests\":" << total
           << ",\"activeWorkingSequenceCount\":" << count
           << ",\"queryWorkingSequenceCount\":" << query_count
           << ",\"referenceWorkingSequenceCount\":" << reference_count
           << ",\"activeReferenceGroupCount\":" << reference_group_count
           << ",\"cumulativeTriplets\":" << processed
           << ",\"tripletKernelEvaluations\":" << processed
           << ",\"tripletSummariesReused\":0,\"cleanTripletsPruned\":0,\"cachedSignalsReused\":0"
           << ",\"methodScansSkipped\":0,\"invalidScheduleTripletsSkipped\":0,\"pairShortlistTripletsSkipped\":0"
           << ",\"fragmentSequencesPruned\":0,\"scanRound\":" << context->progress_round << ",\"maximumDetectionCycles\":1000"
           << ",\"fixedEventCount\":0,\"signalCount\":0,\"eventCount\":" << event_count
           << ",\"cycleTermination\":\"" << (context->finished ? "no-significant-events" : "running")
           << "\",\"fraction\":" << (context->finished ? 1 : total == 0 ? 0 : static_cast<double>(processed) / static_cast<double>(total))
           << ",\"maxChiProfilesScanned\":0,\"maxChiPeakAttempts\":0,\"maxChiCandidatesFound\":0,\"maxChiPeakLimitTriplets\":0"
           << ",\"chimaeraProfilesScanned\":0,\"chimaeraPeakAttempts\":0,\"chimaeraCandidatesFound\":0,\"chimaeraPeakLimitTargets\":0"
           << ",\"geneconvFragmentsScored\":0,\"geneconvQualifiedFragments\":0,\"geneconvCandidatesFound\":0,\"geneconvOverlapRejections\":0,\"geneconvNumericalFallbackTracks\":0"
           << ",\"threeSeqProfilesScanned\":0,\"threeSeqExactEvaluations\":0,\"threeSeqApproximateEvaluations\":0,\"threeSeqCandidatesFound\":0"
           << ",\"bootscanProfilesScanned\":0,\"bootscanCandidateRegionsScored\":0,\"bootscanCandidatesFound\":0,\"bootscanPairProfilesRequested\":0,\"bootscanPairProfileCacheHits\":0,\"bootscanPairProfileCacheMisses\":0,\"bootscanPairProfileCacheEvictions\":0,\"bootscanPairProfileCacheBytes\":0,\"bootscanPairProfileCachePeakBytes\":0"
           << ",\"siscanProfilesScanned\":0,\"siscanWindowsScored\":0,\"siscanCandidateRegionsScored\":0,\"siscanCandidatesFound\":0,\"siscanPermutationDraws\":0,\"siscanContextBuilds\":0,\"siscanContextPairComparisons\":0,\"siscanContextTreeMerges\":0,\"siscanRandomValuesGenerated\":0}"
           ;
    return cached(*context, output.str());
}

NEXT_RDP_KEEPALIVE const char* rdp_get_results_json(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->finished) return "";
    return cached(*context, full_json(*context));
}

NEXT_RDP_KEEPALIVE const char* rdp_get_error(const std::uint32_t handle) {
    auto* context = web_context(handle);
    return context == nullptr ? "The RDP analysis context is invalid." : context->error.c_str();
}

NEXT_RDP_KEEPALIVE const char* rdp_get_signal_plot_json(const std::uint32_t handle, const std::uint32_t signal_id) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto plot = signal_plot_json(*context, signal_id);
    if (plot.empty()) {
        context->error = "The selected signal does not exist or has invalid representative sequences.";
        return "";
    }
    return cached(*context, plot);
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_alignment_json(const std::uint32_t handle, const std::uint32_t event_id, const std::uint32_t flank_sites, const std::uint32_t row_limit) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto view = event_alignment_json(*context, event_id, flank_sites, row_limit);
    if (view.empty()) {
        context->error = "The selected event does not exist.";
        return "";
    }
    return cached(*context, view);
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_trees_json(const std::uint32_t handle, const std::uint32_t event_id) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto view = event_trees_json(*context, event_id);
    if (view.empty()) {
        context->error = "The selected event does not exist.";
        return "";
    }
    return cached(*context, view);
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_phylpro_json(const std::uint32_t handle, const std::uint32_t event_id, const std::uint32_t window_sites, const int gap_mode, const int include_self) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto view = event_phylpro_json(*context, event_id, window_sites, gap_mode, include_self);
    if (view.empty()) {
        context->error = "The selected event does not exist.";
        return "";
    }
    return cached(*context, view);
}

NEXT_RDP_KEEPALIVE int rdp_set_review_state(const std::uint32_t handle, const std::uint32_t signal_id, const int state) {
    auto* context = web_context(handle);
    if (context == nullptr) return 0;
    if (signal_id >= context->review_states.size()) context->review_states.resize(signal_id + 1, 0);
    context->review_states[signal_id] = state;
    return 1;
}

NEXT_RDP_KEEPALIVE int rdp_set_event_review_state(const std::uint32_t handle, const std::uint32_t event_id, const int state) {
    auto* context = web_context(handle);
    if (context == nullptr) return 0;
    if (event_id >= context->review_states.size()) context->review_states.resize(event_id + 1, 0);
    context->review_states[event_id] = state;
    return 1;
}

NEXT_RDP_KEEPALIVE int rdp_update_event(
    const std::uint32_t handle, const std::uint32_t event_id,
    const std::uint32_t recombinant, const std::uint32_t major_parent,
    const std::uint32_t minor_parent, const std::uint32_t beginning,
    const std::uint32_t ending) {
    auto* context = web_context(handle);
    if (context == nullptr || event_id >= context->full.events.size()) return 0;
    auto& event = context->full.events[event_id];
    const auto count = context->alignment.names.size();
    if (recombinant >= count || major_parent >= count || minor_parent >= count ||
        recombinant == major_parent || recombinant == minor_parent ||
        major_parent == minor_parent || beginning == 0 || ending == 0 ||
        beginning > context->alignment.sequences.front().size() ||
        ending > context->alignment.sequences.front().size()) {
        context->error = "The edited event contains an invalid role or breakpoint.";
        return 0;
    }
    const int winner = std::clamp(event.winning_role, 0, 2);
    event.representative_sequences[winner] = static_cast<int>(recombinant);
    event.representative_sequences[(winner + 1) % 3] = static_cast<int>(major_parent);
    event.representative_sequences[(winner + 2) % 3] = static_cast<int>(minor_parent);
    event.beginning = static_cast<int>(beginning);
    event.ending = static_cast<int>(ending);
    context->cache.clear();
    context->error.clear();
    return 1;
}
NEXT_RDP_KEEPALIVE int rdp_update_event_group(
    const std::uint32_t handle, const std::uint32_t event_id,
    const std::uint32_t* sequence_indices, const std::size_t sequence_count,
    const int manual_override) {
    auto* context = web_context(handle);
    if (context == nullptr || event_id >= context->full.events.size()) return 0;
    if (sequence_count > 0 && sequence_indices == nullptr) return 0;
    const auto count = context->alignment.names.size();
    auto& event = context->full.events[event_id];
    const int winner = std::clamp(event.winning_role, 0, 2);
    auto& group = event.sequence_groups[winner];
    group.clear();
    for (std::size_t index = 0; index < sequence_count; ++index) {
        if (sequence_indices[index] >= count) {
            context->error = "The edited co-recombinant group contains an invalid sequence.";
            return 0;
        }
        group.push_back(static_cast<int>(sequence_indices[index]));
    }
    if (group.empty()) group.push_back(event.representative_sequences[winner]);
    (void)manual_override;
    context->cache.clear();
    context->error.clear();
    return 1;
}
NEXT_RDP_KEEPALIVE int rdp_reconcile_after(const std::uint32_t handle, const std::uint32_t) { return web_context(handle) == nullptr ? 0 : 1; }

NEXT_RDP_KEEPALIVE const char* rdp_export_csv(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->finished) return "";
    std::ostringstream output;
    output << "event,recombinant,major_parent,minor_parent,beginning,ending,p_value\n";
    for (std::size_t index = 0; index < context->full.events.size(); ++index) {
        const auto& event = context->full.events[index];
        output << index + 1 << ',' << event.representative_sequences[event.winning_role]
               << ',' << event.representative_sequences[(event.winning_role + 1) % 3]
               << ',' << event.representative_sequences[(event.winning_role + 2) % 3]
               << ',' << event.beginning << ',' << event.ending << ',' << event.probability << '\n';
    }
    return cached(*context, output.str());
}

void write_fasta_record(
    std::ostringstream& output, const std::string& name,
    const std::string& sequence) {
    output << '>' << name << '\n';
    constexpr std::size_t width = 80;
    for (std::size_t offset = 0; offset < sequence.size(); offset += width) {
        output << sequence.substr(offset, std::min(width, sequence.size() - offset)) << '\n';
    }
}

std::string curated_fasta(
    const WebContext& context,
    const std::vector<std::uint8_t>& masked,
    const std::vector<std::uint8_t>& disabled,
    const bool include_enabled) {
    std::ostringstream output;
    for (std::size_t index = 0; index < context.alignment.names.size(); ++index) {
        const bool is_masked = index < masked.size() && masked[index] != 0;
        const bool is_disabled = index < disabled.size() && disabled[index] != 0;
        if (include_enabled ? (is_masked || is_disabled) : (!is_masked && !is_disabled)) {
            continue;
        }
        write_fasta_record(output, context.alignment.names[index],
                           context.alignment.sequences[index]);
    }
    return output.str();
}

bool coordinate_in_tract(
    const std::size_t coordinate, const std::size_t beginning,
    const std::size_t ending, const bool wraps_origin) {
    if (beginning == 0 || ending == 0) return false;
    // ModSeqNumY's wrapped branch treats BPos == EPos as a full circular
    // tract, rather than a one-column interval.
    if (beginning == ending) return true;
    return wraps_origin ? coordinate >= beginning || coordinate <= ending
                        : coordinate >= beginning && coordinate <= ending;
}

bool event_accepted(const WebContext& context, const std::size_t index) {
    return index < context.review_states.size() &&
        context.review_states[index] == 1;
}

bool final_alignment_ready(const WebContext& context) {
    if (!context.finished) return false;
    for (std::size_t index = 0; index < context.full.events.size(); ++index) {
        if (index >= context.review_states.size() ||
            context.review_states[index] == 0) return false;
    }
    return true;
}

std::vector<std::uint8_t> accepted_sequence_mask(const WebContext& context) {
    std::vector<std::uint8_t> removed(context.alignment.sequences.size(), 0);
    for (std::size_t index = 0; index < context.full.events.size(); ++index) {
        if (!event_accepted(context, index)) continue;
        const auto& event = context.full.events[index];
        if (event.winning_role >= 0 && event.winning_role < 3) {
            const int recombinant = event.representative_sequences[event.winning_role];
            if (recombinant >= 0 && recombinant < static_cast<int>(removed.size())) {
                removed[static_cast<std::size_t>(recombinant)] = 1;
            }
            for (const int sequence : event.sequence_groups[event.winning_role]) {
                if (sequence >= 0 && sequence < static_cast<int>(removed.size())) {
                    removed[static_cast<std::size_t>(sequence)] = 1;
                }
            }
        }
    }
    return removed;
}

std::string export_fasta(const WebContext& context) {
    std::ostringstream output;
    for (std::size_t index = 0; index < context.alignment.names.size(); ++index) {
        output << '>' << context.alignment.names[index] << '\n' << context.alignment.sequences[index] << '\n';
    }
    return output.str();
}

NEXT_RDP_KEEPALIVE const char* rdp_export_enabled_sequences_fasta(
    const std::uint32_t handle, const std::uint8_t* masked_sequences,
    const std::size_t mask_length, const std::uint8_t* disabled_sequences,
    const std::size_t disabled_length) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->loaded) return "";
    if (masked_sequences == nullptr || disabled_sequences == nullptr ||
        mask_length != context->alignment.names.size() ||
        disabled_length != context->alignment.names.size()) {
        context->error = "The sequence curation state does not match the loaded alignment.";
        return "";
    }
    const std::vector<std::uint8_t> masked(masked_sequences, masked_sequences + mask_length);
    const std::vector<std::uint8_t> disabled(disabled_sequences, disabled_sequences + disabled_length);
    return cached(*context, curated_fasta(*context, masked, disabled, true));
}
NEXT_RDP_KEEPALIVE const char* rdp_export_masked_or_disabled_sequences_fasta(
    const std::uint32_t handle, const std::uint8_t* masked_sequences,
    const std::size_t mask_length, const std::uint8_t* disabled_sequences,
    const std::size_t disabled_length) {
    auto* context = web_context(handle);
    if (context == nullptr || !context->loaded) return "";
    if (masked_sequences == nullptr || disabled_sequences == nullptr ||
        mask_length != context->alignment.names.size() ||
        disabled_length != context->alignment.names.size()) {
        context->error = "The sequence curation state does not match the loaded alignment.";
        return "";
    }
    const std::vector<std::uint8_t> masked(masked_sequences, masked_sequences + mask_length);
    const std::vector<std::uint8_t> disabled(disabled_sequences, disabled_sequences + disabled_length);
    return cached(*context, curated_fasta(*context, masked, disabled, false));
}
NEXT_RDP_KEEPALIVE const char* rdp_export_recombinant_sequences_removed_fasta(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !final_alignment_ready(*context)) return "";
    const auto removed = accepted_sequence_mask(*context);
    std::ostringstream output;
    for (std::size_t index = 0; index < context->alignment.names.size(); ++index) {
        if (removed[index] == 0) write_fasta_record(output, context->alignment.names[index], context->alignment.sequences[index]);
    }
    return cached(*context, output.str());
}
NEXT_RDP_KEEPALIVE const char* rdp_export_recombinant_columns_removed_fasta(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !final_alignment_ready(*context)) return "";
    const std::size_t length = context->alignment.sequences.front().size();
    std::vector<std::uint8_t> removed(length + 1, 0);
    for (std::size_t index = 0; index < context->full.events.size(); ++index) {
        if (!event_accepted(*context, index)) continue;
        const auto& event = context->full.events[index];
        const std::size_t beginning = static_cast<std::size_t>(std::clamp(event.beginning, 1, static_cast<int>(length)));
        const std::size_t ending = static_cast<std::size_t>(std::clamp(event.ending, 1, static_cast<int>(length)));
        const bool wraps = beginning > ending;
        for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
            if (coordinate_in_tract(coordinate, beginning, ending, wraps)) removed[coordinate] = 1;
        }
    }
    std::ostringstream output;
    for (std::size_t sequence = 0; sequence < context->alignment.names.size(); ++sequence) {
        std::string retained;
        for (std::size_t coordinate = 1; coordinate <= length; ++coordinate) {
            if (removed[coordinate] == 0) retained.push_back(context->alignment.sequences[sequence][coordinate - 1]);
        }
        if (!retained.empty()) write_fasta_record(output, context->alignment.names[sequence], retained);
    }
    return cached(*context, output.str());
}
NEXT_RDP_KEEPALIVE const char* rdp_export_recombination_free_fasta(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !final_alignment_ready(*context)) return "";
    auto sequences = context->alignment.sequences;
    for (std::size_t index = 0; index < context->full.events.size(); ++index) {
        if (!event_accepted(*context, index)) continue;
        const auto& event = context->full.events[index];
        std::vector<int> affected = event.sequence_groups[event.winning_role];
        affected.push_back(event.representative_sequences[event.winning_role]);
        for (const int sequence : affected) {
            if (sequence < 0 || sequence >= static_cast<int>(sequences.size())) continue;
            const std::size_t beginning = static_cast<std::size_t>(std::clamp(event.beginning, 1, static_cast<int>(sequences[sequence].size())));
            const std::size_t ending = static_cast<std::size_t>(std::clamp(event.ending, 1, static_cast<int>(sequences[sequence].size())));
            const bool wraps = beginning > ending;
            for (std::size_t coordinate = 1; coordinate <= sequences[sequence].size(); ++coordinate) {
                if (coordinate_in_tract(coordinate, beginning, ending, wraps)) sequences[sequence][coordinate - 1] = '-';
            }
        }
    }
    std::ostringstream output;
    for (std::size_t index = 0; index < sequences.size(); ++index) write_fasta_record(output, context->alignment.names[index], sequences[index]);
    return cached(*context, output.str());
}
NEXT_RDP_KEEPALIVE const char* rdp_export_fragmented_fasta(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr || !final_alignment_ready(*context)) return "";
    struct Fragment { std::string name; std::string sequence; };
    auto remainder = context->alignment.sequences;
    std::vector<Fragment> fragments;
    for (std::size_t index = 0; index < context->full.events.size(); ++index) {
        if (!event_accepted(*context, index)) continue;
        const auto& event = context->full.events[index];
        std::vector<int> affected = event.sequence_groups[event.winning_role];
        affected.push_back(event.representative_sequences[event.winning_role]);
        for (const int sequence : affected) {
            if (sequence < 0 || sequence >= static_cast<int>(remainder.size())) continue;
            const std::size_t beginning = static_cast<std::size_t>(std::clamp(event.beginning, 1, static_cast<int>(remainder[sequence].size())));
            const std::size_t ending = static_cast<std::size_t>(std::clamp(event.ending, 1, static_cast<int>(remainder[sequence].size())));
            const bool wraps = beginning > ending;
            std::string fragment(remainder[sequence].size(), '-');
            bool copied = false;
            for (std::size_t coordinate = 1; coordinate <= remainder[sequence].size(); ++coordinate) {
                if (!coordinate_in_tract(coordinate, beginning, ending, wraps)) continue;
                fragment[coordinate - 1] = remainder[sequence][coordinate - 1];
                remainder[sequence][coordinate - 1] = '-';
                copied = true;
            }
            if (!copied) continue;
            std::ostringstream name;
            name << context->alignment.names[sequence] << "|RDP_event_" << index + 1
                 << "|bp_" << event.beginning << '_' << event.ending;
            if (wraps) name << "|circular";
            fragments.push_back({name.str(), std::move(fragment)});
        }
    }
    std::ostringstream output;
    for (std::size_t index = 0; index < remainder.size(); ++index) write_fasta_record(output, context->alignment.names[index], remainder[index]);
    for (const auto& fragment : fragments) write_fasta_record(output, fragment.name, fragment.sequence);
    return cached(*context, output.str());
}

NEXT_RDP_KEEPALIVE const char* rdp_export_project_json(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    std::ostringstream output;
    output << "{\"schema\":\"org.rdp-web.project/v1alpha20\",\"engineVersion\":\"nextRDP-core 0.1.0\",\"dataset\":{\"format\":";
    json_string(output, context->alignment.format);
    output << ",\"sequences\":[";
    for (std::size_t index = 0; index < context->alignment.names.size(); ++index) {
        if (index) output << ',';
        output << "{\"name\":";
        json_string(output, context->alignment.names[index]);
        output << ",\"sequence\":";
        json_string(output, context->alignment.sequences[index]);
        output << '}';
    }
    output << "],\"alignmentLength\":" << context->alignment.sequences.front().size() << "},\"analysis\":";
    if (context->finished) output << full_json(*context); else output << "null";
    output << '}';
    return cached(*context, output.str());
}

}  // extern "C"

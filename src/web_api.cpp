#include "analysis.hpp"
#include "legacy_method_state.hpp"
#include "legacy_optional/bootscan.hpp"
#include "legacy_optional/siscan.hpp"
#include "legacy_optional_bridge.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
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
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};

ParsedFasta parse_fasta(const std::string& text) {
    ParsedFasta parsed;
    std::string current_name;
    std::string current_sequence;
    const auto finish = [&]() {
        if (current_name.empty()) return;
        parsed.names.push_back(current_name);
        parsed.sequences.push_back(current_sequence);
        current_name.clear();
        current_sequence.clear();
    };
    std::size_t start = 0;
    while (start <= text.size()) {
        const auto end = text.find('\n', start);
        const auto line_end = end == std::string::npos ? text.size() : end;
        std::string line = text.substr(start, line_end - start);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            finish();
            current_name = line.substr(1);
            const auto first_space = current_name.find_first_of(" \t");
            if (first_space != std::string::npos) current_name.resize(first_space);
        } else if (!current_name.empty()) {
            for (const char character : line) {
                if (!std::isspace(static_cast<unsigned char>(character))) {
                    current_sequence.push_back(character);
                }
            }
        }
        if (end == std::string::npos) break;
        start = end + 1;
    }
    finish();
    if (parsed.names.size() < 3) {
        throw std::runtime_error("The alignment must contain at least three FASTA sequences");
    }
    const std::size_t length = parsed.sequences.front().size();
    if (length == 0) throw std::runtime_error("The FASTA alignment contains an empty sequence");
    for (const auto& sequence : parsed.sequences) {
        if (sequence.size() != length) {
            throw std::runtime_error("All FASTA sequences must have the same length");
        }
    }
    return parsed;
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
};

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

std::string summary_json(const WebContext& context) {
    const auto count = context.alignment.sequences.size();
    const auto length = context.alignment.sequences.front().size();
    std::size_t variable_sites = 0;
    std::size_t informative_sites = 0;
    for (std::size_t position = 0; position < length; ++position) {
        std::array<int, 5> states{};
        for (const auto& sequence : context.alignment.sequences) {
            const char character = sequence[position];
            int state = 0;
            if (character == 'A' || character == 'a') state = 1;
            else if (character == 'C' || character == 'c') state = 2;
            else if (character == 'G' || character == 'g') state = 3;
            else if (character == 'T' || character == 't' || character == 'U' || character == 'u') state = 4;
            if (state != 0) ++states[state];
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
    for (std::size_t first = 0; first < count; ++first) {
        for (std::size_t second = first + 1; second < count; ++second) {
            std::size_t compared = 0;
            std::size_t matches = 0;
            for (std::size_t position = 0; position < length; ++position) {
                const char left = context.alignment.sequences[first][position];
                const char right = context.alignment.sequences[second][position];
                if (left == '-' || right == '-' || left == '?' || right == '?') continue;
                ++compared;
                if (std::toupper(static_cast<unsigned char>(left)) ==
                    std::toupper(static_cast<unsigned char>(right))) ++matches;
            }
            if (compared == 0) continue;
            const double identity = static_cast<double>(matches) / static_cast<double>(compared);
            minimum_identity = std::min(minimum_identity, identity);
            identity_sum += identity;
            ++identity_pairs;
        }
    }
    std::ostringstream output;
    output << "{\"format\":\"FASTA\",\"sequenceCount\":" << count
           << ",\"alignmentLength\":" << length
           << ",\"activeSequenceCount\":" << count
           << ",\"tripletCount\":" << count * (count - 1) * (count - 2) / 6
           << ",\"variableSiteCount\":" << variable_sites
           << ",\"informativeSiteCount\":" << informative_sites
           << ",\"minimumPairIdentity\":" << (identity_pairs == 0 ? 0.0 : minimum_identity)
           << ",\"meanPairIdentity\":" << (identity_pairs == 0 ? 0.0 : identity_sum / identity_pairs)
           << ",\"recommendedMinimumDistance\":0,\"partitionBoundaries\":[1," << length << "]"
           << ",\"sequences\":[";
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0) output << ',';
        std::size_t missing = 0;
        for (const char character : context.alignment.sequences[index]) {
            if (character == '-' || character == '?' || character == 'N' || character == 'n') ++missing;
        }
        output << "{\"index\":" << index << ",\"name\":";
        json_string(output, context.alignment.names[index]);
        output << ",\"validSites\":" << length - missing
               << ",\"missingSites\":" << missing
               << ",\"missingFraction\":" << std::setprecision(17)
               << static_cast<double>(missing) / static_cast<double>(length)
               << ",\"masked\":false}";
    }
    output << "],\"warnings\":[]}";
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

std::string signal_plot_json(const WebContext& context, std::uint32_t signal_id) {
    const auto& result = context.full;
    if (signal_id >= result.events.size()) return {};
    const auto& event = result.events[signal_id];
    const int length = context.alignment.sequences.empty()
        ? 0 : static_cast<int>(context.alignment.sequences.front().size());
    if (length <= 0) return {};
    std::array<const std::string*, 3> sequences{};
    for (int role = 0; role < 3; ++role) {
        const int index = event.representative_sequences[role];
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
    std::vector<std::size_t> optional_coordinates;
    std::array<std::vector<double>, 3> optional_values;
    double optional_minimum = std::numeric_limits<double>::infinity();
    double optional_maximum = -std::numeric_limits<double>::infinity();
    if (program == 5 && event.bootscan_available) {
        const auto optional_alignment =
            next_rdp_legacy_optional_bridge::make_alignment(
                context.alignment.sequences);
        const std::array<std::uint32_t, 3> triplet{
            static_cast<std::uint32_t>(event.representative_sequences[0]),
            static_cast<std::uint32_t>(event.representative_sequences[1]),
            static_cast<std::uint32_t>(event.representative_sequences[2])};
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
            static_cast<std::uint32_t>(event.representative_sequences[0]),
            static_cast<std::uint32_t>(event.representative_sequences[1]),
            static_cast<std::uint32_t>(event.representative_sequences[2])};
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
        exact_bootscan_profile || exact_siscan_profile;
    const bool exact_detection_profile =
        exact_rdp_profile || exact_optional_profile;
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
    if (!exact_detection_profile) {
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
            for (int pair = 0; pair < 3; ++pair) {
                point[pair] = static_cast<double>(
                    event.rdp_profile.counts[pair][coordinate_index]) /
                    divisor;
            }
        } else if (exact_optional_profile &&
                   coordinate_index < optional_coordinates.size()) {
            for (int pair = 0; pair < 3; ++pair) {
                point[pair] = optional_values[pair][coordinate_index];
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
        minimum = event.rdp_profile.minimum;
        maximum = event.rdp_profile.maximum;
    } else if (exact_optional_profile) {
        minimum = optional_minimum;
        maximum = optional_maximum;
    }
    if (!std::isfinite(minimum)) minimum = 0.0;
    if (!std::isfinite(maximum)) maximum = 0.0;
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"signalId\":" << signal_id
           << ",\"windowSites\":" << (program == 8 ? 0 : exact_bootscan_profile || exact_siscan_profile
               ? (exact_bootscan_profile ? context.options.bootscan_window_sites
                                         : context.options.siscan_window_sites)
               : requested_window)
           << ",\"alignmentLength\":" << length
           << ",\"method\":\"" << method << "\",\"metric\":\"" << metric
           << "\",\"profileContext\":\""
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
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const int flank = std::clamp(static_cast<int>(requested_flank), 5, 100);
    const auto all_sequences = event_review_sequences(context, event, 64);
    const std::size_t candidate_count = all_sequences.size();
    const std::size_t limit = std::clamp<std::size_t>(row_limit, 3, 64);
    const std::size_t visible_count = std::min(candidate_count, limit);
    const auto is_group_member = [&](const int sequence) {
        return std::find(event.sequence_groups[0].begin(), event.sequence_groups[0].end(), sequence) !=
            event.sequence_groups[0].end();
    };
    const auto role = [&](const int sequence) {
        if (sequence == event.representative_sequences[0]) return "recombinant";
        if (sequence == event.representative_sequences[1]) return "major-parent";
        if (sequence == event.representative_sequences[2]) return "minor-parent";
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
               << "\",\"queryReferenceInputRole\":\"not-applied\",\"referenceGroup\":null"
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
    const auto sequences = event_review_sequences(context, event, 48);
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const auto in_tract = [&](const int coordinate) {
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
        if (sequence == event.representative_sequences[0]) return "recombinant";
        if (sequence == event.representative_sequences[1]) return "major-parent";
        if (sequence == event.representative_sequences[2]) return "minor-parent";
        return "evidence";
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
    static constexpr std::array<const char*, 6> names{
        "5-prime-outside", "5-prime-inside", "3-prime-outside",
        "3-prime-inside", "outside-tract", "inside-tract"};
    std::ostringstream output;
    output << "{\"eventId\":" << event_id
           << ",\"method\":\"neighbour-joining\",\"distance\":\"Jukes-Cantor\""
           << ",\"njKernel\":\"supplied-clearcut-float\",\"distanceEncoding\":\"source-tree2arrayp2-midpoint-ranks\""
           << ",\"bootstrapGenerator\":\"disabled-rdp-5.93-event-path\",\"bootstrapSupport\":\"not-applied\""
           << ",\"negativeBranchPolicy\":\"absolute-five-decimal-serialization\",\"analyticalBranchParsing\":\"four-decimal-clamped-complete-edge-repair\""
           << ",\"treeRooting\":\"source-tree2arrayp2-midpoint\",\"collapseEncoding\":\"unbootstrapped-raw-tree-copy\""
           << ",\"displayRooting\":\"arbitrary-internal-node\",\"bootstrapCollapseCutoff\":null"
           << ",\"bootstrapReplicates\":0,\"randomSeed\":3,\"flankVariableSiteTarget\":60"
           << ",\"subsampled\":false,\"sequenceCap\":48,\"fragmentAssisted\":false,\"leaves\":[";
    const int root = static_cast<int>(sequences.size());
    for (std::size_t index = 0; index < sequences.size(); ++index) {
        if (index != 0) output << ',';
        const int sequence = sequences[index];
        output << "{\"node\":" << index << ",\"workingSequenceIndex\":" << sequence
               << ",\"sequenceIndex\":" << sequence << ",\"sequenceName\":";
        json_string(output, context.alignment.names[static_cast<std::size_t>(sequence)]);
        output << ",\"fragmentEventId\":null,\"role\":\"" << role(sequence)
               << "\",\"queryReferenceInputRole\":\"not-applied\",\"referenceGroup\":null"
               << ",\"masked\":" << (context.masked[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"disabled\":" << (context.disabled[static_cast<std::size_t>(sequence)] ? "true" : "false")
               << ",\"currentGroupMember\":false,\"automaticGroupMember\":false,\"trace\":false}";
    }
    output << "],\"regions\":[";
    for (int region = 0; region < 6; ++region) {
        if (region != 0) output << ',';
        const auto coordinates = region_coordinates(region);
        const bool usable = sequences.size() >= 3 && coordinates.size() >= 3;
        output << "{\"name\":\"" << names[region] << "\",\"sites\":" << coordinates.size()
               << ",\"sequences\":" << sequences.size() << ",\"usable\":" << (usable ? "true" : "false")
               << ",\"nodeCount\":" << (usable ? root + 1 : 0) << ",\"root\":" << (usable ? root : 0)
               << ",\"bootstrapReplicates\":0,\"supportedInternalBranches\":0,\"internalBranches\":0"
               << ",\"rawDistanceRankLevels\":" << (usable ? 1 : 0)
               << ",\"collapsedDistanceRankLevels\":" << (usable ? 1 : 0)
               << ",\"negativeBranchesNormalized\":0,\"bootstrapRandomSeed\":3,\"edges\":[";
        if (usable) {
            for (std::size_t index = 0; index < sequences.size(); ++index) {
                if (index != 0) output << ',';
                const double edge_length = distance(sequences[0], sequences[index], coordinates);
                output << "{\"from\":" << root << ",\"to\":" << index
                       << ",\"length\":" << std::setprecision(17) << edge_length
                       << ",\"bootstrapSupport\":null,\"internal\":false,\"collapsed\":false}";
            }
        }
        output << "]}";
    }
    output << "]}";
    return output.str();
}

double profile_correlation(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() < 2 || left.size() != right.size()) return 0.0;
    double left_mean = 0.0;
    double right_mean = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        left_mean += left[index];
        right_mean += right[index];
    }
    left_mean /= left.size();
    right_mean /= right.size();
    double numerator = 0.0;
    double left_sum = 0.0;
    double right_sum = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        const double a = left[index] - left_mean;
        const double b = right[index] - right_mean;
        numerator += a * b;
        left_sum += a * a;
        right_sum += b * b;
    }
    const double denominator = std::sqrt(left_sum * right_sum);
    return denominator == 0.0 ? 0.0 : numerator / denominator;
}

std::string event_phylpro_json(const WebContext& context, std::uint32_t event_id,
                               std::uint32_t requested_window, int gap_mode, int include_self) {
    if (event_id >= context.full.events.size()) return {};
    const auto& event = context.full.events[event_id];
    const int length = static_cast<int>(context.alignment.sequences.front().size());
    const int window = std::clamp(static_cast<int>(requested_window), 10, 5000);
    const int half = std::max(1, window / 2);
    std::vector<int> context_sequences;
    for (std::size_t index = 0; index < context.alignment.sequences.size(); ++index) {
        if (index < context.disabled.size() && context.disabled[index]) continue;
        context_sequences.push_back(static_cast<int>(index));
    }
    const std::array<int, 3> targets = event.representative_sequences;
    const int stride = std::max(1, (length + 2047) / 2048);
    std::vector<int> coordinates;
    for (int coordinate = 1; coordinate <= length; coordinate += stride) coordinates.push_back(coordinate);
    if (coordinates.empty() || coordinates.back() != length) coordinates.push_back(length);
    const auto profile_for = [&](int target, int center) {
        std::array<std::vector<double>, 2> sides;
        for (int side = 0; side < 2; ++side) {
            for (const int context_sequence : context_sequences) {
                if (!include_self && context_sequence == targets[target]) continue;
                int compared = 0;
                int differences = 0;
                const int start = side == 0 ? center - half + 1 : center + 1;
                const int end = side == 0 ? center : center + half;
                for (int coordinate = start; coordinate <= end; ++coordinate) {
                    const char a = plot_base_at(context.alignment.sequences[static_cast<std::size_t>(targets[target])], coordinate, context.options.circular);
                    const char b = plot_base_at(context.alignment.sequences[static_cast<std::size_t>(context_sequence)], coordinate, context.options.circular);
                    if (!plot_base(a) || !plot_base(b)) {
                        if (gap_mode != 0) continue;
                        continue;
                    }
                    ++compared;
                    if (std::toupper(static_cast<unsigned char>(a)) != std::toupper(static_cast<unsigned char>(b))) ++differences;
                }
                sides[side].push_back(compared == 0 ? 1.0 : static_cast<double>(differences) / compared);
            }
        }
        return profile_correlation(sides[0], sides[1]);
    };
    std::array<std::vector<double>, 3> profiles;
    for (const int coordinate : coordinates) {
        for (int target = 0; target < 3; ++target) {
            profiles[target].push_back(profile_for(target, coordinate));
        }
    }
    std::array<std::size_t, 3> minimum_indices{};
    for (int target = 0; target < 3; ++target) {
        minimum_indices[target] = static_cast<std::size_t>(
            std::min_element(profiles[target].begin(), profiles[target].end()) - profiles[target].begin());
    }
    const auto nearest = [&](int coordinate) {
        std::size_t best = 0;
        int best_distance = std::numeric_limits<int>::max();
        for (std::size_t index = 0; index < coordinates.size(); ++index) {
            const int difference = std::abs(coordinates[index] - coordinate);
            if (difference < best_distance) {
                best = index;
                best_distance = difference;
            }
        }
        return best;
    };
    double minimum = 0.0;
    double maximum = 0.0;
    for (const auto& profile : profiles) {
        for (const double value : profile) {
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
    }
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"eventId\":" << event_id << ",\"status\":\"source-shaped-active-unvalidated\",\"kernel\":\"FindSubSeqPP-MakePDstMat-UpdatePDstMat-PPRegression\""
           << ",\"roleOrder\":\"recombinant-major-parent-minor-parent\",\"columnSelection\":\"polymorphic-after-gap-policy\",\"distance\":\"pairwise-Hamming-count\",\"correlation\":\"Pearson-source-single-output\",\"significanceTest\":\"not-implemented-in-supplied-rdp5\",\"optimization\":\"three-target-rows-linear-in-context\""
           << ",\"maskedContextIncluded\":true,\"disabledContextExcluded\":true,\"fragmentContextIncluded\":false"
           << ",\"circular\":" << (context.options.circular ? "true" : "false") << ",\"windowSites\":" << window
           << ",\"halfWindowSites\":" << half << ",\"windowCapped\":false,\"gapMode\":\""
           << (gap_mode ? "strip-any-missing-column" : "ignore-missing-pairwise") << "\",\"includeSelf\":" << (include_self ? "true" : "false")
           << ",\"eligibleColumns\":" << length << ",\"contextSequences\":" << context_sequences.size()
           << ",\"targetContextComparisons\":" << (coordinates.size() * context_sequences.size() * 3)
           << ",\"rollingUpdates\":" << coordinates.size() * 3 << ",\"evaluatedPoints\":" << coordinates.size()
           << ",\"returnedPoints\":" << coordinates.size() << ",\"minimumValue\":" << minimum << ",\"maximumValue\":" << maximum
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
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        if (index != 0) output << ',';
        output << "{\"alignmentPosition\":" << coordinates[index]
               << ",\"recombinant\":" << profiles[0][index]
               << ",\"majorParent\":" << profiles[1][index]
               << ",\"minorParent\":" << profiles[2][index] << '}';
    }
    output << "]}";
    return output.str();
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
               << ",\"probability\":" << event.probability
               << ",\"beginning\":" << event.beginning
               << ",\"ending\":" << event.ending << ",\"representativeSequences\":["
               << event.representative_sequences[0] << ','
               << event.representative_sequences[1] << ','
               << event.representative_sequences[2] << "],\"sequenceGroups\":[";
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
               << ",\"bootscanDiscovery\":";
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
        context->fasta.assign(reinterpret_cast<const char*>(bytes), length);
        context->alignment = parse_fasta(context->fasta);
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
    const int polish_breakpoints, const int /*query_reference_mode*/,
    const std::uint32_t* /*reference_groups*/, const std::size_t /*reference_group_count*/,
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
    if (!options.enable_rdp && !options.enable_geneconv &&
        !options.enable_maxchi && !options.enable_chimaera &&
        !options.enable_three_seq && !options.enable_bootscan &&
        !options.enable_siscan) {
        context->error = "Select RDP or at least one other discovery method.";
        return 0;
    }
    options.progress_callback = &web_progress_callback;
    options.progress_user = context;
    context->options = std::move(options);
    context->masked.assign(masked_sequences, masked_sequences + mask_length);
    context->disabled.assign(disabled_sequences, disabled_sequences + disabled_length);
    options.masked_sequences = context->masked;
    options.disabled_sequences = context->disabled;
    context->started = true;
    context->finished = false;
    context->progress_phase = 0;
    context->progress_round = 1;
    context->progress_processed_triplets = 0;
    context->progress_total_triplets = count < 3
        ? 0 : static_cast<int>(count * (count - 1) * (count - 2) / 6);
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
    } catch (const std::exception& error) {
        context->error = error.what();
        return -1;
    }
}

NEXT_RDP_KEEPALIVE int rdp_reconcile(const std::uint32_t handle) {
    auto* context = web_context(handle);
    return context != nullptr && context->finished ? 1 : 0;
}

NEXT_RDP_KEEPALIVE void rdp_cancel(const std::uint32_t /*handle*/) {}

NEXT_RDP_KEEPALIVE const char* rdp_get_progress_json(const std::uint32_t handle) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    const auto count = context->alignment.names.size();
    const auto total = count < 3 ? 0 : count * (count - 1) * (count - 2) / 6;
    const int phase = context->finished ? 2 : context->progress_phase;
    const char* phase_name = phase == 1 ? "cyclic-rescan"
        : phase == 2 ? "complete" : "primary";
    const auto processed = context->finished
        ? static_cast<std::size_t>(total)
        : static_cast<std::size_t>(std::max(0, context->progress_processed_triplets));
    const auto event_count = context->finished
        ? context->full.events.size()
        : static_cast<std::size_t>(std::max(0, context->progress_event_count));
    std::ostringstream output;
    output << "{\"state\":\"" << (context->finished ? "done" : context->started ? "running" : "idle")
           << "\",\"phase\":\"" << phase_name
           << "\",\"processedTriplets\":" << processed
           << ",\"totalTriplets\":" << total
           << ",\"correctionTests\":" << total
           << ",\"activeWorkingSequenceCount\":" << count
           << ",\"queryWorkingSequenceCount\":0,\"referenceWorkingSequenceCount\":0"
           << ",\"activeReferenceGroupCount\":0,\"cumulativeTriplets\":" << processed
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
    output << "{\"schema\":\"org.rdp-web.project/v1alpha20\",\"engineVersion\":\"nextRDP-core 0.1.0\",\"dataset\":{\"format\":\"FASTA\",\"sequences\":[";
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

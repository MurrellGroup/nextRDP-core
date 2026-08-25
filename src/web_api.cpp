#include "analysis.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cctype>
#include <iomanip>
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
    bool enable_geneconv = false;
    bool enable_maxchi = false;
    bool enable_chimaera = false;
    bool enable_three_seq = false;
    bool polish_breakpoints_with_burt = false;
    std::vector<unsigned char> masked;
    std::vector<unsigned char> disabled;
    std::vector<unsigned int> reference_groups;
    std::vector<int> review_states;
};

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

std::string full_json(const WebContext& context) {
    const auto& result = context.full;
    std::ostringstream output;
    output << std::setprecision(17)
           << "{\"sourceFaithfulCore\":true,\"engineVersion\":\"nextRDP-core 0.1.0\""
           << ",\"enabledMethods\":[\"RDP\"";
    if (context.enable_geneconv) output << ",\"GENECONV\"";
    if (context.enable_maxchi) output << ",\"MAXCHI\"";
    if (context.enable_chimaera) output << ",\"CHIMAERA\"";
    if (context.enable_three_seq) output << ",\"3SEQ\"";
    output << "]"
           << ",\"sequenceCount\":" << result.sequence_count
           << ",\"sequenceLength\":" << result.sequence_length
           << ",\"tripletCount\":" << result.triplet_count
           << ",\"rawCandidateCount\":" << result.raw_candidate_count
           << ",\"events\":[";
    for (std::size_t index = 0; index < result.events.size(); ++index) {
        if (index != 0) output << ',';
        const auto& event = result.events[index];
        output << "{\"id\":" << index
               << ",\"program\":" << event.program_flag
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
        output << "]}";
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
    return web_context(handle) == nullptr ? 0 : (requested == 0 ? 1 : 1);
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
    const std::uint32_t handle, const int circular, const int /*correction_mode*/,
    const double p_value_cutoff, const std::uint32_t window_sites,
    const int maxchi_enabled, const std::uint32_t /*maxchi_window_sites*/,
    const int chimaera_enabled, const std::uint32_t /*chimaera_window_sites*/,
    const int geneconv_enabled, const std::uint32_t /*geneconv_mismatch_scale*/,
    const std::uint32_t /*geneconv_max_overlaps*/, const int threeseq_enabled,
    const int /*bootscan_primary_enabled*/, const int /*bootscan_secondary_enabled*/,
    const std::uint32_t /*bootscan_window_sites*/, const std::uint32_t /*bootscan_step_sites*/,
    const std::uint32_t /*bootscan_bootstrap_replicates*/, const double /*bootscan_support_cutoff*/,
    const std::uint32_t /*bootscan_random_seed*/, const int /*siscan_primary_enabled*/,
    const int /*siscan_secondary_enabled*/, const std::uint32_t /*siscan_window_sites*/,
    const std::uint32_t /*siscan_step_sites*/, const std::uint32_t /*siscan_scan_permutations*/,
    const std::uint32_t /*siscan_p_value_permutations*/, const std::uint32_t /*siscan_random_seed*/,
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
    context->options = {circular != 0, p_value_cutoff, static_cast<int>(window_sites)};
    context->options.enable_geneconv = geneconv_enabled != 0;
    context->options.enable_maxchi = maxchi_enabled != 0;
    context->options.enable_chimaera = chimaera_enabled != 0;
    context->options.enable_three_seq = threeseq_enabled != 0;
    context->options.polish_breakpoints_with_burt = polish_breakpoints != 0;
    context->enable_geneconv = context->options.enable_geneconv;
    context->enable_maxchi = context->options.enable_maxchi;
    context->enable_chimaera = context->options.enable_chimaera;
    context->enable_three_seq = context->options.enable_three_seq;
    context->polish_breakpoints_with_burt = context->options.polish_breakpoints_with_burt;
    context->masked.assign(masked_sequences, masked_sequences + mask_length);
    context->disabled.assign(disabled_sequences, disabled_sequences + disabled_length);
    context->started = true;
    context->finished = false;
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
    std::ostringstream output;
    output << "{\"state\":\"" << (context->finished ? "done" : context->started ? "running" : "idle")
           << "\",\"phase\":\"" << (context->finished ? "complete" : "primary")
           << "\",\"processedTriplets\":" << (context->finished ? total : 0)
           << ",\"totalTriplets\":" << total
           << ",\"correctionTests\":" << total
           << ",\"activeWorkingSequenceCount\":" << count
           << ",\"queryWorkingSequenceCount\":0,\"referenceWorkingSequenceCount\":0"
           << ",\"activeReferenceGroupCount\":0,\"cumulativeTriplets\":" << (context->finished ? total : 0)
           << ",\"tripletKernelEvaluations\":" << (context->finished ? total : 0)
           << ",\"tripletSummariesReused\":0,\"cleanTripletsPruned\":0,\"cachedSignalsReused\":0"
           << ",\"methodScansSkipped\":0,\"invalidScheduleTripletsSkipped\":0,\"pairShortlistTripletsSkipped\":0"
           << ",\"fragmentSequencesPruned\":0,\"scanRound\":1,\"maximumDetectionCycles\":1000"
           << ",\"fixedEventCount\":0,\"signalCount\":0,\"eventCount\":" << (context->finished ? context->full.events.size() : 0)
           << ",\"cycleTermination\":\"" << (context->finished ? "no-significant-events" : "running")
           << "\",\"fraction\":" << (context->finished ? 1 : 0)
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

NEXT_RDP_KEEPALIVE const char* rdp_get_signal_plot_json(const std::uint32_t handle, const std::uint32_t /*signal_id*/) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    return cached(*context, "{\"signalId\":0,\"windowSites\":30,\"method\":\"RDP\",\"metric\":\"pair-identity\",\"profileContext\":\"original-alignment-reconstruction\",\"detectionProfileExact\":false,\"minimumValue\":0,\"maximumValue\":1,\"points\":[]}");
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_alignment_json(const std::uint32_t handle, const std::uint32_t /*event_id*/, const std::uint32_t /*flank_sites*/, const std::uint32_t /*row_limit*/) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    return cached(*context, "{\"eventId\":0,\"alignmentLength\":0,\"circular\":true,\"fragmentAssisted\":false,\"requestedFlankSites\":0,\"candidateRowCount\":0,\"omittedRowCount\":0,\"rows\":[],\"panels\":[]}");
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_trees_json(const std::uint32_t handle, const std::uint32_t /*event_id*/) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    return cached(*context, "{\"eventId\":0,\"method\":\"neighbour-joining\",\"distance\":\"Jukes-Cantor\",\"leaves\":[],\"regions\":[]}");
}

NEXT_RDP_KEEPALIVE const char* rdp_get_event_phylpro_json(const std::uint32_t handle, const std::uint32_t /*event_id*/, const std::uint32_t /*window_sites*/, const int /*gap_mode*/, const int /*include_self*/) {
    auto* context = web_context(handle);
    if (context == nullptr) return "";
    return cached(*context, "{\"eventId\":0,\"windowSites\":60,\"gapMode\":\"ignore-missing-pairwise\",\"includeSelf\":false,\"points\":[],\"minima\":[],\"breakpoints\":[]}");
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

NEXT_RDP_KEEPALIVE int rdp_update_event(const std::uint32_t, const std::uint32_t, const std::uint32_t, const std::uint32_t, const std::uint32_t, const std::uint32_t, const std::uint32_t) { return 1; }
NEXT_RDP_KEEPALIVE int rdp_update_event_group(const std::uint32_t, const std::uint32_t, const std::uint32_t*, const std::size_t, const int) { return 1; }
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

std::string export_fasta(const WebContext& context) {
    std::ostringstream output;
    for (std::size_t index = 0; index < context.alignment.names.size(); ++index) {
        output << '>' << context.alignment.names[index] << '\n' << context.alignment.sequences[index] << '\n';
    }
    return output.str();
}

NEXT_RDP_KEEPALIVE const char* rdp_export_enabled_sequences_fasta(const std::uint32_t handle, const std::uint8_t*, const std::size_t, const std::uint8_t*, const std::size_t) {
    auto* context = web_context(handle);
    return context == nullptr ? "" : cached(*context, export_fasta(*context));
}
NEXT_RDP_KEEPALIVE const char* rdp_export_masked_or_disabled_sequences_fasta(const std::uint32_t handle, const std::uint8_t*, const std::size_t, const std::uint8_t*, const std::size_t) { return rdp_export_enabled_sequences_fasta(handle, nullptr, 0, nullptr, 0); }
NEXT_RDP_KEEPALIVE const char* rdp_export_recombinant_sequences_removed_fasta(const std::uint32_t handle) { return rdp_export_enabled_sequences_fasta(handle, nullptr, 0, nullptr, 0); }
NEXT_RDP_KEEPALIVE const char* rdp_export_recombinant_columns_removed_fasta(const std::uint32_t handle) { return rdp_export_enabled_sequences_fasta(handle, nullptr, 0, nullptr, 0); }
NEXT_RDP_KEEPALIVE const char* rdp_export_recombination_free_fasta(const std::uint32_t handle) { return rdp_export_enabled_sequences_fasta(handle, nullptr, 0, nullptr, 0); }
NEXT_RDP_KEEPALIVE const char* rdp_export_fragmented_fasta(const std::uint32_t handle) { return rdp_export_enabled_sequences_fasta(handle, nullptr, 0, nullptr, 0); }

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

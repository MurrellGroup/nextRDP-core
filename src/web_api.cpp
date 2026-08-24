#include "analysis.hpp"

#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#if defined(__EMSCRIPTEN__)
#include <emscripten/emscripten.h>
#define NEXT_RDP_KEEPALIVE EMSCRIPTEN_KEEPALIVE
#else
#define NEXT_RDP_KEEPALIVE
#endif

namespace {

std::string result_json;
std::string error_message;

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

}  // extern "C"

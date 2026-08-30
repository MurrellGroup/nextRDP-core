#pragma once

#include "burt_state.hpp"
#include "legacy_optional/bootscan.hpp"
#include "legacy_optional/chimaera.hpp"
#include "legacy_optional/geneconv.hpp"
#include "legacy_optional/maxchi.hpp"
#include "legacy_optional/siscan.hpp"
#include "legacy_optional/threeseq.hpp"
#include "scan_state.hpp"
#include "round_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <exception>
#include <string>
#include <vector>

using RdpProgressCallback = void (*)(
    int phase, int round, int processed_triplets, int total_triplets,
    int event_count, void* user);
using RdpCancellationCallback = bool (*)(void* user);

class RdpAnalysisCancelled final : public std::exception {
 public:
    const char* what() const noexcept override {
        return "The scan was cancelled.";
    }
};

struct RdpInitialAnalysisOptions {
    bool circular = true;
    double p_value_cutoff = 0.05;
    int window_sites = 30;
    // RDP is the default primary method, but the source UI permits an
    // optional-method-only scan.  Keep this independent from the optional
    // method switches so callers can request (for example) MaxChi alone.
    bool enable_rdp = true;
    bool enable_geneconv = false;
    int geneconv_mismatch_scale = 1;
    int geneconv_max_overlaps = 1;
    bool enable_maxchi = false;
    int maxchi_window_sites = 70;
    bool enable_chimaera = false;
    int chimaera_window_sites = 60;
    bool enable_three_seq = false;
    bool correction_bonferroni = true;
    // Optional BootScan/SISCAN lanes are entered only when requested. Their
    // implementations are isolated from the source-faithful RDP scheduler,
    // so the default RDP path remains unchanged.
    bool enable_bootscan = false;
    bool enable_bootscan_secondary = false;
    int bootscan_window_sites = 200;
    int bootscan_step_sites = 20;
    int bootscan_bootstrap_replicates = 100;
    double bootscan_support_cutoff = 0.70;
    unsigned int bootscan_random_seed = 3;
    bool enable_siscan = false;
    bool enable_siscan_secondary = false;
    int siscan_window_sites = 200;
    int siscan_step_sites = 20;
    int siscan_scan_permutations = 100;
    int siscan_p_value_permutations = 1000;
    unsigned int siscan_random_seed = 3;
    std::vector<unsigned char> masked_sequences;
    std::vector<unsigned char> disabled_sequences;
    // When enabled, replace the exploratory MakeAListP2 list with the
    // source MakeAnalysisListQvR ordering: one query (group 0) and two
    // references from different positive reference groups.  Keeping this in
    // the core options makes the browser and native entry points share the
    // same list construction rather than merely labelling an exploratory
    // scan as query-vs-reference.
    bool query_reference_mode = false;
    std::vector<unsigned int> reference_groups;
    bool polish_breakpoints_with_burt = false;

    // Optional source call-order captures.  GENECONV now has a direct
    // AlistGC2 scheduler and does not need a capture; these remain useful for
    // the methods whose list schedulers are still being lifted.
    std::string geneconv_call_order_path;
    std::string geneconv_call_count_path;
    std::string maxchi_call_order_path;
    std::string maxchi_call_count_path;
    std::string chimaera_call_order_path;
    std::string chimaera_call_count_path;

    // Optional observer used by the web wrapper. It reports source-stage
    // boundaries and completed rows from the low-level XOver walk; it never
    // fabricates progress for work that has not completed.
    RdpProgressCallback progress_callback = nullptr;
    void* progress_user = nullptr;
    RdpCancellationCallback cancellation_callback = nullptr;
    void* cancellation_user = nullptr;
};

struct RdpInitialAnalysisResult {
    RdpScanState alignment;
    RdpRawEventState events;
    std::vector<double> store_lpv;
    int store_lpv_upper_bound = 8;
};

// The values passed to the original VB DrawPlots routine.  RDP plots are not
// ordinary nucleotide-window identities: they are XOverHomologyP rolling
// agreement counts over information-rich sites, drawn at XDiffPos.
struct RdpSlidingWindowProfile {
    int window_sites = 0;
    int divisor = 0;
    int alignment_length = 0;
    double minimum = 0.0;
    double maximum = 0.0;
    bool exact = false;
    std::vector<int> positions;
    std::array<std::vector<int>, 3> counts;
};

struct RdpFinalEvent {
    int event_number = 0;
    int program = 0;
    int winning_role = 0;
    // For target-specific optional methods (currently CHIMAERA), retain the
    // selected method's candidate role separately from the later consensus
    // winner.  The browser plotter uses this to draw the same one-target
    // trace that the source method inspected.
    int method_target_role = -1;
    double probability = 1.0;
    int beginning = 0;
    int ending = 0;
    std::array<int, 3> representative_sequences{};
    // The source review plot is keyed to the discovery triplet, which is
    // independent from the later consensus winner/role rotation. Keep that
    // triplet so the browser can draw and label pair traces in source order.
    std::array<int, 3> profile_sequences{};
    bool profile_sequences_available = false;
    std::array<std::vector<int>, 3> sequence_groups;
    std::array<double, 3> consensus{};
    bool burt_attempted = false;
    bool burt_applied = false;
    RdpBurtResult burt{};
    RdpSlidingWindowProfile rdp_profile{};
    bool maxchi_available = false;
    next_rdp_legacy_optional::MaxChiDiscoveryCandidate maxchi_discovery{};
    bool chimaera_available = false;
    next_rdp_legacy_optional::ChimaeraDiscoveryCandidate chimaera_discovery{};
    bool geneconv_available = false;
    next_rdp_legacy_optional::GeneconvDiscoveryCandidate geneconv_discovery{};
    bool three_seq_available = false;
    next_rdp_legacy_optional::ThreeSeqDiscoveryCandidate three_seq_discovery{};
    bool bootscan_available = false;
    next_rdp_legacy_optional::BootscanDiscoveryCandidate bootscan_discovery{};
    next_rdp_legacy_optional::BootscanPlotProfile bootscan_profile{};
    next_rdp_legacy_optional::BootscanRecheckEvidence bootscan_recheck{};
    bool siscan_available = false;
    next_rdp_legacy_optional::SiscanDiscoveryCandidate siscan_discovery{};
    next_rdp_legacy_optional::SiscanPlotProfile siscan_profile{};
    next_rdp_legacy_optional::SiscanRecheckEvidence siscan_recheck{};
};

struct RdpFullAnalysisResult {
    int sequence_count = 0;
    int sequence_length = 0;
    int triplet_count = 0;
    int raw_candidate_count = 0;
    std::vector<RdpFinalEvent> events;
};

struct RdpFullAnalysisTrace {
    std::vector<std::array<int, 2>> selected_slots;
    std::vector<std::vector<unsigned char>> done_before_selection;
    std::vector<int> done_row_upper_bounds;
    std::vector<RdpConsensusState> consensus_states;
    std::vector<RdpCompleteRoundState> complete_round_states;
    std::vector<RdpFinalTrimState> consensus_candidate_states;
    std::vector<RdpFinalTrimState> final_candidate_states;
    std::vector<std::vector<std::array<int, 3>>> inner_triplets;
    std::vector<std::vector<std::array<int, 3>>> outer_triplets;
    std::vector<std::vector<unsigned char>> inner_redo;
    std::vector<std::vector<unsigned char>> outer_redo;
    std::vector<std::vector<int>> trace_sub_after_expansion;
    std::vector<std::vector<int>> expanded_actual_sequence_sizes;
    std::vector<std::vector<int>> retained_actual_sequence_sizes;
    std::vector<std::vector<int>> fragment_reference_counts;
    std::vector<std::vector<int>> fragment_reference_counts_before_drop;
    std::vector<std::vector<int>> actual_sequence_sizes_before_drop;
    std::vector<RdpRawEventState> collection_events_before_adjustment;
    std::vector<std::vector<unsigned char>> adjusted_pairs_to_rescan;
    std::vector<RdpRawEventState> adjusted_events;
    std::vector<RdpRawEventState> events_before_drop;
    std::vector<RdpRawEventState> post_rescan_events;
};

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_text(
    const std::string& fasta_text,
    const RdpInitialAnalysisOptions& options = {});

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_file(
    const std::string& fasta_path,
    const RdpInitialAnalysisOptions& options = {});

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_text(
    const std::string& fasta_text,
    const RdpInitialAnalysisOptions& options = {});

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_file(
    const std::string& fasta_path,
    const RdpInitialAnalysisOptions& options = {});

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_file_with_trace(
    const std::string& fasta_path, RdpFullAnalysisTrace& trace,
    const RdpInitialAnalysisOptions& options = {});

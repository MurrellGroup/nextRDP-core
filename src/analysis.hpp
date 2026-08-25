#pragma once

#include "scan_state.hpp"
#include "round_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <string>
#include <vector>

struct RdpInitialAnalysisOptions {
    bool circular = true;
    double p_value_cutoff = 0.05;
    int window_sites = 30;
    bool enable_geneconv = false;
    bool enable_maxchi = false;
    bool enable_chimaera = false;
    bool enable_three_seq = false;
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
};

struct RdpInitialAnalysisResult {
    RdpScanState alignment;
    RdpRawEventState events;
    std::vector<double> store_lpv;
    int store_lpv_upper_bound = 8;
};

struct RdpFinalEvent {
    int event_number = 0;
    int winning_role = 0;
    double probability = 1.0;
    int beginning = 0;
    int ending = 0;
    std::array<int, 3> representative_sequences{};
    std::array<std::vector<int>, 3> sequence_groups;
    std::array<double, 3> consensus{};
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

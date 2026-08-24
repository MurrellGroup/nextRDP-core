#pragma once

#include "scan_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <vector>

struct RdpBreakpointFlanks {
    std::array<int, 4> positions{};
    std::array<double, 4> informative_counts{};
};

RdpBreakpointFlanks make_rdp_breakpoint_flanks(
    const RdpScanState& scan_state, int beginning, int ending,
    const std::array<int, 3>& sequences, int variable_site_target = 60,
    int total_site_target = 0);

struct RdpCorrelationState {
    std::vector<float> correlation;
    std::vector<float> inversion;
    std::vector<float> tested_correlation;
    std::array<double, 2> intermediate{};
    std::array<double, 3> results{};
};

RdpCorrelationState calculate_rdp_correlations(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<std::vector<double>, 3>& regional_matrices,
    double minimum_offset = 1.0e-14);

struct RdpCorrelationDecisionState {
    RdpCorrelationState correlations;
    std::array<unsigned char, 3> warnings{};
};

RdpCorrelationDecisionState finalize_rdp_correlations(
    int next_no, RdpCorrelationState correlations,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<double>& summary_matrix,
    const std::vector<double>& regional_distance_matrix);

std::vector<float> make_rdp_local_distance_panels(
    const RdpScanState& scan_state, const std::array<int, 4>& starts,
    const std::array<int, 4>& ends,
    const std::array<int, 3>& sequences);

void apply_rdp_distance_warnings(
    int next_no, const std::array<int, 3>& sequences,
    const std::vector<float>& local_distance_panels,
    std::array<unsigned char, 3>& warnings);

std::vector<int> make_rdp_good_comparisons(
    const RdpScanState& scan_state,
    const std::array<int, 4>& breakpoint_flanks);

struct RdpRoleLists {
    std::array<unsigned char, 3> inside{};
    std::array<unsigned char, 3> outside{};
};

RdpRoleLists make_rdp_role_lists(
    const std::array<unsigned char, 2>& minimum_pair);

std::vector<unsigned char> make_rdp_acceptable_correlations(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside,
    const std::vector<float>& first_direct,
    const std::vector<float>& second_direct,
    std::vector<float>& first_adjusted,
    std::vector<float>& second_adjusted,
    std::vector<float>& first_collapsed,
    std::vector<float>& second_collapsed);

struct RdpCandidateLists {
    std::vector<int> list_by_threshold;
    std::vector<int> inverse_by_threshold;
    std::array<int, 30> last_by_threshold{};
    std::vector<int> list;
    std::vector<int> inverse;
    std::array<int, 3> last{};
    std::vector<double> probability_scores;
    std::vector<double> probability_values;
    std::vector<double> t_values;
    std::vector<double> totals;
    std::vector<double> list_scores;
    double result{};
};

RdpCandidateLists make_rdp_candidate_lists(
    int next_no, const std::vector<int>& good_comparisons,
    const std::array<int, 3>& sequences,
    RdpCorrelationDecisionState& correlations,
    const std::vector<unsigned char>& dont_redo,
    const std::vector<unsigned char>& acceptable_correlations);

struct RdpActualEventCallState {
    int role{};
    std::array<int, 6> region_sizes_before{};
    std::vector<int> breakpoint_matches_before;
    std::vector<float> best_matches_before;
    std::vector<double> acceptable_sequences_before;
    std::vector<int> found_before;
    std::array<int, 3> candidate_last_before{};
    std::vector<int> candidate_list_before;
    std::vector<unsigned char> inversion_state_before;
    std::array<int, 2> candidate_scratch_before{};
    std::array<int, 3> trace_sequences_before{};
    std::array<double, 2> match_before{};
    std::array<int, 4> sequence_scratch_before{};
    std::array<unsigned char, 6> tried_permutations_before{};
    std::vector<unsigned char> role_membership_before;
    int result{};
    std::array<int, 6> region_sizes_after{};
    std::vector<int> breakpoint_matches_after;
    std::vector<float> best_matches_after;
    std::vector<double> acceptable_sequences_after;
    std::vector<int> found_after;
    std::array<int, 2> candidate_scratch_after{};
    std::array<int, 3> trace_sequences_after{};
    std::array<double, 2> match_after{};
    std::array<int, 4> sequence_scratch_after{};
    std::array<unsigned char, 6> tried_permutations_after{};
};

struct RdpActualEventResolution {
    RdpCandidateLists candidates;
    RdpCorrelationDecisionState correlations;
    std::vector<double> acceptable_sequences;
    std::vector<unsigned char> unfound;
    std::vector<int> breakpoint_matches;
    std::vector<float> best_matches;
    std::vector<int> event_overlap_mask;
    std::vector<int> beginning_overlap_mask;
    std::vector<int> ending_overlap_mask;
    std::array<int, 6> region_sizes{};
    std::array<int, 3> candidate_last_before_strip{};
    std::vector<int> candidate_list_before_strip;
    std::vector<int> candidate_inverse_before_strip;
    std::array<int, 3> inversion_penalty{};
    std::array<RdpActualEventCallState, 3> calls;
};

RdpActualEventResolution resolve_rdp_actual_events(
    int sequence_length, int next_no,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 6>& starts,
    const std::array<int, 6>& ends,
    RdpCorrelationDecisionState correlations,
    RdpCandidateLists candidates,
    const std::vector<unsigned char>& dont_redo,
    const RdpRawEventState& events,
    int permanent_next_no = -1);

struct RdpEventSetState {
    std::array<int, 3> candidate_last{};
    std::vector<int> candidate_list;
    std::array<int, 6> set_totals{};
    std::vector<unsigned char> role_sets;
};

RdpEventSetState find_rdp_event_sets(
    int sequence_length, int next_no, int beginning, int ending,
    const std::array<int, 3>& sequences,
    const RdpRawEventState& events);

struct RdpTreeCompatibilityCallState {
    int role{};
    std::array<int, 3> compatibility_before{};
    std::array<int, 3> reverse_compatibility_before{};
    std::array<int, 3> nonrecombinant_last_before{};
    std::vector<int> done_before;
    std::vector<int> nonrecombinant_list_before;
    std::array<double, 3> list_distances{};
    std::array<int, 3> compatibility_after{};
    std::array<int, 3> reverse_compatibility_after{};
    std::array<int, 3> nonrecombinant_last_after{};
};

struct RdpTreeCompatibilityState {
    std::array<double, 3> background_list_distances{};
    std::array<double, 3> region_list_distances{};
    std::array<int, 3> background_compatibility{};
    std::array<int, 3> background_reverse_compatibility{};
    std::array<int, 3> region_compatibility{};
    std::array<int, 3> region_reverse_compatibility{};
    std::array<RdpTreeCompatibilityCallState, 6> calls;
};

struct RdpTreeCompatibilityFlowState {
    RdpEventSetState event_sets;
    std::array<int, 3> background{};
    std::array<int, 3> background_secondary{};
    std::array<int, 3> background_sets{};
    std::array<int, 3> background_secondary_sets{};
    std::array<int, 3> region{};
    std::array<int, 3> region_secondary{};
    std::array<int, 3> region_sets{};
    std::array<int, 3> region_secondary_sets{};
    std::vector<RdpTreeCompatibilityCallState> calls;
};

RdpTreeCompatibilityFlowState run_rdp_tree_compatibility_flow(
    int sequence_length, int next_no, int beginning, int ending,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& background_ancestor_matrix,
    const std::vector<float>& region_ancestor_matrix,
    const std::vector<float>& background_secondary_matrix,
    const std::vector<float>& region_secondary_matrix,
    const RdpRawEventState& events);

struct RdpPatternState {
    std::vector<double> pattern;
    std::vector<unsigned char> done;
    std::vector<double> acceptable_sequences;
};

RdpPatternState check_rdp_sequence_patterns(
    int sequence_length, int next_no,
    const std::array<int, 3>& sequences,
    const std::array<int, 3>& starts,
    const std::array<int, 3>& ends,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<short>& sequence_data,
    const std::vector<double>& acceptable_sequences);

struct RdpFinalTrimState {
    std::array<int, 3> candidate_last{};
    std::vector<int> candidate_list;
    std::array<int, 3> nonrecombinant_last{};
    std::vector<int> nonrecombinant_list;
    std::vector<double> acceptable_sequences;
    std::vector<std::array<int, 2>> synthetic_event_roles;
    std::vector<unsigned char> strict_removed;
    std::vector<unsigned char> strict_readded;
};

RdpFinalTrimState run_rdp_final_trim_candidate_maintenance(
    int next_no, int rwinpp,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<unsigned char, 2>& minimum_pair,
    const std::array<unsigned char, 3>& inside_roles,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::vector<unsigned char>& unfound,
    const std::vector<float>& correlations,
    const std::vector<float>& inversions,
    const std::vector<float>& local_distance_panels,
    const std::vector<float>& first_direct_small,
    const std::vector<float>& region_direct_small,
    const std::vector<float>& first_ancestor_small,
    const std::vector<float>& region_ancestor_small,
    const std::vector<float>& first_collapsed_small,
    const std::vector<float>& region_collapsed_small,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<double>& acceptable_sequences);

std::vector<double> calculate_rdp_match_evidence(
    int sequence_length, int next_no, int beginning, int ending,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<short>& sequence_data,
    const std::vector<double>& acceptable_sequences,
    bool conservative_grouping);

RdpFinalTrimState make_rdp_consensus_candidates(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::vector<float>& correlations,
    const std::vector<float>& inversions,
    const std::vector<float>& first_direct,
    const std::vector<float>& region_direct,
    const std::vector<float>& first_ancestor,
    const std::vector<float>& region_ancestor,
    const std::vector<float>& first_direct_small,
    const std::vector<float>& region_direct_small,
    const std::vector<float>& first_ancestor_small,
    const std::vector<float>& region_ancestor_small,
    const std::vector<float>& first_collapsed_small,
    const std::vector<float>& region_collapsed_small,
    RdpFinalTrimState state, bool conservative_grouping);

struct RdpMaximumDistanceState {
    int result = 0;
    int included_last = -1;
    std::vector<int> included_sequences;
    std::vector<unsigned char> included_mask;
    std::vector<unsigned char> representative_mask;
    std::vector<int> informative_to_position;
    std::vector<int> position_to_informative;
    std::vector<unsigned char> nucleotide_map;
    std::vector<float> split_scores;
    std::array<float, 3> distance_totals{};
    std::array<int, 3> distance_counts{};
    std::array<float, 3> maximum_distances{};
};

RdpMaximumDistanceState calculate_rdp_maximum_distances(
    int sequence_length, int next_no,
    const std::array<int, 3>& sequences,
    int beginning, int ending,
    const std::vector<short>& sequence_data,
    const std::vector<unsigned char>& masked_sequences);

struct RdpConsensusInputs {
    int next_no{};
    int permanent_next_no{};
    std::array<int, 6> comparison_matrix{};
    std::array<double, 3> list_correlation{};
    std::array<double, 3> simple_distance_strength{};
    std::array<int, 3> simple_distance_score{};
    std::array<double, 3> phylpro{};
    std::array<double, 3> phylpro_secondary{};
    std::array<double, 3> phylpro_collapsed{};
    std::array<double, 3> subtree_score{};
    std::array<double, 3> split_distance{};
    std::array<int, 3> outlier_index{};
    std::array<double, 3> subtree_phylpro{};
    std::array<double, 3> subtree_score_secondary{};
    std::array<double, 3> subtree_phylpro_secondary{};
    std::array<int, 3> compatibility{};
    std::array<int, 3> compatibility_secondary{};
    std::array<int, 3> compatibility_tertiary{};
    std::array<int, 3> compatibility_quaternary{};
    std::array<int, 3> region_compatibility{};
    std::array<int, 3> region_compatibility_secondary{};
    std::array<int, 3> region_compatibility_tertiary{};
    std::array<int, 3> region_compatibility_quaternary{};
    std::array<int, 3> post_trim_compatibility{};
    std::array<int, 3> post_trim_region_compatibility{};
    std::array<double, 3> triplet_score{};
    std::array<double, 3> bad_distances{};
    std::array<int, 3> outside_list{};
    std::array<double, 3> list_correlation_secondary{};
    std::array<double, 3> list_correlation_tertiary{};
    std::array<int, 3> outlier_check{};
    std::array<float, 3> maximum_distance{};
    std::array<std::array<int, 2>, 3> ranks{};
};

struct RdpConsensusState {
    RdpConsensusInputs rounded_inputs;
    std::array<double, 3> consensus{};
    std::array<double, 78> decision_scores{};
    int winning_role{};
};

RdpConsensusState make_rdp_consensus(RdpConsensusInputs inputs);

struct RdpSplitDistanceState {
    std::array<double, 3> distances{};
    std::array<int, 3> outlier_index{};
};

RdpSplitDistanceState calculate_rdp_split_distances(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::vector<float>& ancestor_rows,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix,
    const std::vector<int>& done);

struct RdpSimpleDistanceState {
    std::array<std::array<int, 2>, 3> ranks{};
    std::array<int, 3> scores{};
    std::array<double, 3> strengths{};
};

RdpSimpleDistanceState calculate_rdp_simple_distances(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix);

std::array<int, 3> calculate_rdp_outlier_checks(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::vector<float>& background_ancestor_rows,
    const std::vector<float>& region_ancestor_rows);

std::array<double, 3> calculate_rdp_bad_distances(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<unsigned char>& unfound,
    const std::vector<float>& correlations,
    const std::vector<float>& ancestor_rows,
    const std::vector<float>& local_distance_panels);

struct RdpListCorrelationState {
    std::array<double, 3> mismatches{};
    std::array<double, 3> expected_strength{};
    std::array<double, 3> absent_strength{};
    std::array<std::array<double, 3>, 3> regional_strength{};
};

RdpListCorrelationState calculate_rdp_list_correlations(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<float>& inversions,
    const std::vector<float>& tested_correlations,
    const std::vector<float>& background_ancestor_rows,
    const std::vector<float>& region_ancestor_rows);

RdpFinalTrimState apply_rdp_strict_group_constraints(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix,
    const std::vector<float>& background_rows,
    const std::vector<float>& region_rows,
    const std::vector<float>& background_ancestor_rows,
    const std::vector<float>& region_ancestor_rows,
    RdpFinalTrimState state);

std::vector<int> make_rdp_relevant_sequences(
    int next_no, const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list);

RdpRawEventState prepare_rdp_collection_event_list(
    int next_no, int winning_role,
    const std::array<int, 3>& sequences,
    const std::array<int, 2>& trace,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<double>& acceptable_sequences,
    const RdpRawEventState& events);

struct RdpCollectedEventsState {
    int result{};
    std::vector<RdpRawEvent> events;
};

RdpCollectedEventsState make_rdp_parent_collect_events(
    int sequence_length, int next_no, int role,
    std::array<int, 6> region_sizes,
    const std::vector<int>& overlap_sequence,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    int add_num, const std::array<int, 3>& sequences,
    const std::array<int, 2>& trace,
    const RdpRawEventState& source_events);

RdpTreeCompatibilityCallState make_rdp_tree_compatibility_call(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix, int role,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& ancestor_matrix,
    const std::array<double, 3>& list_distances,
    std::array<int, 3>& compatibility,
    std::array<int, 3>& reverse_compatibility,
    std::array<int, 3>& nonrecombinant_last,
    std::vector<int>& nonrecombinant_list);

RdpTreeCompatibilityState evaluate_rdp_tree_compatibility(
    int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& background_ancestor_matrix,
    const std::vector<float>& region_ancestor_matrix);

struct RdpPhylProScoreState {
    std::vector<int> trace_involved;
    std::array<double, 3> scores{};
    std::array<double, 3> temporary_scores{};
    std::array<double, 3> sub_scores{};
    std::array<double, 3> sub_distance_scores{};
    double result{};
};

RdpPhylProScoreState make_rdp_phylpro_scores(
    int next_no, double minimum_offset,
    const std::vector<int>& done_this,
    const std::array<int, 3>& sequences,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix);

std::vector<int> make_rdp_score_filter(
    int next_no, const std::array<int, 3>& sequences,
    const std::vector<float>& raw_background_rows,
    const std::vector<float>& ancestor_background_rows,
    const std::vector<float>& ancestor_region_rows);

struct RdpTripletGroupState {
    std::vector<int> counts;
    std::vector<int> done;
    std::vector<int> groups;
    std::array<double, 3> minimum_distances{};
    int result{};
};

RdpTripletGroupState make_rdp_triplet_groups(
    int role, int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<float>& ancestor_background_rows,
    std::array<double, 3> minimum_distances = {});

double make_rdp_triplet_tree_score(
    int role, int next_no, const std::array<int, 3>& sequences,
    const std::vector<float>& ancestor_background_rows,
    const std::vector<float>& ancestor_region_rows,
    const RdpTripletGroupState& groups);

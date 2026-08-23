#pragma once

#include "scan_state.hpp"

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

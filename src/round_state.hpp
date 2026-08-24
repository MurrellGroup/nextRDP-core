#pragma once

#include "distance_state.hpp"
#include "identification_state.hpp"

#include <array>
#include <vector>

struct RdpRoundPrefixState {
    std::array<int, 3> sequences{};
    std::vector<float> breakpoint_distance;
    std::vector<float> remainder_distance;
    RdpDistanceState region_distance;
    RdpEventDistanceMatrices matrices;
    std::array<unsigned char, 2> minimum_pair{};
    std::array<unsigned char, 3> sequence_pair{};
    std::vector<float> background_adjusted;
    std::vector<float> region_adjusted;
    RdpBreakpointFlanks breakpoint_flanks;
    std::array<int, 6> starts{};
    std::array<int, 6> ends{};
    std::vector<double> summary_matrix;
    std::vector<double> regional_distance_matrix;
    std::array<std::vector<double>, 3> correlation_matrices;
    RdpCorrelationDecisionState correlation_decisions;
    std::vector<float> local_distance_panels;
    std::vector<int> good_comparisons;
    std::vector<float> first_direct_small;
    std::vector<float> region_direct_small;
    std::vector<float> first_adjusted_small;
    std::vector<float> region_adjusted_small;
    std::vector<float> first_collapsed_small;
    std::vector<float> region_collapsed_small;
    RdpRoleLists role_lists;
    std::vector<unsigned char> acceptable_correlations;
    std::vector<unsigned char> dont_redo;
    RdpActualEventResolution actual_resolution;
};

RdpRoundPrefixState identify_rdp_round_prefix(
    const RdpScanState& scan_state,
    const RdpDistanceState& full_distance,
    const RdpRawEventState& events,
    const RdpRawEvent& selected,
    const std::vector<unsigned char>& missing_data,
    int minimum_sequence_size = 20);

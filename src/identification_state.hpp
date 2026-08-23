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

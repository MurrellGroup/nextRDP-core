#pragma once

#include "scan_state.hpp"

#include <array>

struct RdpBreakpointFlanks {
    std::array<int, 4> positions{};
    std::array<double, 4> informative_counts{};
};

RdpBreakpointFlanks make_rdp_breakpoint_flanks(
    const RdpScanState& scan_state, int beginning, int ending,
    const std::array<int, 3>& sequences, int variable_site_target = 60,
    int total_site_target = 0);


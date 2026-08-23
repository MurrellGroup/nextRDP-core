#pragma once

#include "scan_state.hpp"

#include <vector>

struct RdpDistanceState {
    std::vector<float> differences;
    std::vector<float> valid_sites;
    std::vector<float> distance;
    std::vector<short> redo_distance;
    double average_distance_accumulator = 0.0;
    double upper_distance = 0.0;
};

RdpDistanceState build_rdp_distance_state(
    const RdpScanState& scan_state, int start_position, int end_position);

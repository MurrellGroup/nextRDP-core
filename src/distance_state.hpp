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

struct RdpEventDistanceMatrices {
    std::vector<float> background;
    std::vector<float> event_region;
};

RdpDistanceState build_rdp_distance_state(
    const RdpScanState& scan_state, int start_position, int end_position,
    bool apply_initial_caller_normalization = true);

RdpEventDistanceMatrices finish_rdp_event_distances(
    int next_no, const RdpDistanceState& full_alignment,
    const RdpDistanceState& event_region);

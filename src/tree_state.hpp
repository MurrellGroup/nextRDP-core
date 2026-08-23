#pragma once

#include "distance_state.hpp"

#include <vector>

struct RdpTreeState {
    std::vector<float> tree_distance;
    std::vector<short> tree_x;
    std::vector<short> tree_y;
    std::vector<double> node_length;
};

RdpTreeState build_rdp_upgma_tree_state(
    int next_no, const RdpDistanceState& distance_state);

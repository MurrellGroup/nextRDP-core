#pragma once

#include "xover_state.hpp"

#include <array>
#include <vector>

struct RdpSelectionState {
    int next_no = -1;
    int slot_upper_bound = 0;
    int make_test_result = 0;
    int total_candidates = 0;
    int done_target = 0;
    bool found = false;
    double probability = 0.0;
    std::array<int, 2> trace{0, 0};
    std::vector<unsigned char> done_sequence;
    std::vector<double> test_probabilities;
};

RdpSelectionState select_rdp_best_event(
    const RdpRawEventState& events, int next_no,
    double probability_cutoff = 0.05,
    const std::vector<unsigned char>& existing_done = {},
    int existing_done_row_upper_bound = -1);

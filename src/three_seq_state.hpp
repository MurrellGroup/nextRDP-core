#pragma once

#include "scan_state.hpp"

#include <array>
#include <vector>

struct RdpThreeSeqSide {
    int beginning = 0;
    int ending = 0;
    int first_count = 0;
    int second_count = 0;
    int excursion = 0;
    double probability = 1.0;
    bool significant = false;
};

struct RdpThreeSeqResult {
    int informative_last = -1;
    std::array<RdpThreeSeqSide, 2> sides{};
    // The TSXOver prefix retains these maps for its subsequent CentreBP
    // call.  Keeping them in the result lets the FinalTrim port perform the
    // same alignment-coordinate centering after FindSubSeqTS2/CheckWrap.
    std::vector<int> difference_position;
    std::vector<int> position_difference;
};

std::vector<float> make_rdp_three_seq_probability_table(
    int maximum_permutation_distance = 44);

// Source port of the FindSubSeqTS2/CheckwrapC/GetTSPVal prefix used by
// TSXOver(1).  The supplied RDP 3seqTable is a VB6 column-major Single array;
// table_bound is its inclusive upper bound.
RdpThreeSeqResult evaluate_rdp_three_seq(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    bool circular, int mc_flag, int mc_correction,
    double lowest_probability, const std::vector<float>& probability_table,
    int table_bound, bool use_ts2 = false);

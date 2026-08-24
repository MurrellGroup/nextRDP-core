#pragma once

#include "distance_state.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

struct RdpRescanScreenSettings {
    int circular = 1;
    int correction_tests = 1;
    int correction_flag = 0;
    double probability_cutoff = 0.05;
    int target = 0;
    int short_output = 100;
    int fss_upper_bound = 125;
    int half_window = 15;
    short full_window = 30;
    int probability_file_flag = 0;
    int probability_one_upper_bound = 171;
    int probability_two_upper_bound = 171;
    int factorial_three_upper_bound = 97;
};

std::vector<short> flatten_rdp_triplets(
    const std::vector<std::array<int, 3>>& triplets);

std::vector<unsigned char> screen_rdp_rescan_triplets(
    const std::vector<std::array<int, 3>>& triplets,
    const RdpScanState& scan_state,
    const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state,
    const RdpRescanScreenSettings& settings,
    std::vector<unsigned char>& fss_rdp,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three,
    std::vector<double>& fact);

void propagate_rdp_group_pairs(
    int next_no, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    std::vector<unsigned char>& pairs_to_rescan);

std::vector<std::array<int, 3>> make_rdp_inner_scan_triplets(
    const RdpScanState& original_scan_state,
    const std::vector<unsigned char>& initially_screened_triplets,
    int winning_role, const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    int permanent_next_no, int minimum_sequence_size,
    const std::vector<unsigned char>& pairs_to_rescan);

std::vector<std::array<int, 3>> make_rdp_outer_scan_triplets(
    const RdpScanState& original_scan_state,
    const std::vector<unsigned char>& initially_screened_triplets,
    int current_next_no, int starting_next_no, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    int permanent_next_no, int minimum_sequence_size,
    const std::vector<unsigned char>& pairs_to_rescan,
    int subvalid_upper_bound, const std::vector<float>& subvalid);

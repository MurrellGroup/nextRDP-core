#pragma once

#include "scan_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

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

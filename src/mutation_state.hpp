#pragma once

#include "xover_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

struct RdpErasedTracts {
    std::vector<short> sequence_data;
    std::vector<unsigned char> missing_data;
    std::vector<short> saved_tracts;
    std::vector<int> breakpoints;
};

RdpErasedTracts erase_rdp_recombinant_tracts(
    int sequence_length, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    int beginning, int ending, std::vector<int> breakpoints,
    const std::vector<short>& sequence_data,
    const std::vector<unsigned char>& missing_data);

void rebuild_rdp_recombinant_tracts(
    int sequence_length, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& breakpoints,
    const std::vector<short>& saved_tracts,
    std::vector<short>& sequence_data);

void make_rdp_fragment_rows(
    int sequence_length, int next_no, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    int beginning, int ending, std::vector<int>& breakpoints,
    std::vector<short>& sequence_data,
    std::vector<unsigned char>& missing_data);

void erase_rdp_original_tracts(
    int sequence_length, int next_no, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& selected_candidates,
    int beginning, int ending, std::vector<int>& breakpoints,
    std::vector<short>& sequence_data,
    std::vector<unsigned char>& missing_data);

std::vector<int> calculate_rdp_actual_sequence_sizes(
    int sequence_length, int next_no, int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<short>& sequence_data);

struct RdpRedistributedEvents {
    RdpRawEventState events;
    std::vector<unsigned char> pairs_to_rescan;
    std::array<int, 101> event_counts{};
};

struct RdpAdjustedEvents {
    RdpRawEventState events;
    std::vector<unsigned char> done_sequence;
    int done_row_upper_bound = 0;
    int done_slot_upper_bound = 0;
    std::vector<unsigned char> pairs_to_rescan;
};

RdpAdjustedEvents adjust_rdp_events_exact(
    int next_no, int winning_role, double lowest_probability,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const RdpRawEventState& source_events,
    const std::vector<unsigned char>& source_done,
    int source_done_row_upper_bound,
    int source_done_slot_upper_bound);

struct RdpDroppedFragmentEvents {
    int next_no = -1;
    RdpRawEventState events;
    std::vector<int> reference_counts;
    std::vector<int> trace_sub;
    std::vector<int> actual_sequence_sizes;
};

RdpRedistributedEvents redistribute_rdp_events(
    int next_no, int winning_role, double lowest_probability,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const RdpRawEventState& source_events);

RdpDroppedFragmentEvents drop_rdp_unused_fragment_events(
    int original_next_no, int expanded_next_no, int minimum_sequence_size,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const RdpRawEventState& events_before_outer_scan,
    const RdpRawEventState& outer_scan_events);

#pragma once

#include "distance_state.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"

#include <array>
#include <vector>

#if defined(_WIN32)
#define RDP_XOVER_CALL __stdcall
#else
#define RDP_XOVER_CALL
#endif

struct Dna5XoverApi {
    int(RDP_XOVER_CALL* find_subsequence)(
        int*, int, int, int, int, int, int, int, int, unsigned char*, int,
        char*, unsigned char*);
    int(RDP_XOVER_CALL* calculate_homology)(
        short, int, int, short, char*, int*);
    int(RDP_XOVER_CALL* find_next)(
        int, int, int, int, int, int, int, int*);
};

struct RdpFirstXoverState {
    std::array<int, 3> sequences{};
    std::array<int, 3> agreement_counts{};
    int informative_length = 0;
    int homology_length = 0;
    int homology_sequence_length = 0;
    int initial_high_homology = 0;
    int high_homology = 0;
    int med_homology = 0;
    int low_homology = 0;
    int homology_start = 0;
    int next_position = -1;
    int active_sequence = -1;
    int active_major_parent = -1;
    int active_minor_parent = -1;
    int sequence_daughter = -1;
    int sequence_minor = -1;
    int xover_sequence_ub = 0;
    int homology_ub = 0;
    std::array<double, 3> average_homology{};
    std::vector<char> xover_sequence;
    std::vector<int> homology;
};

RdpFirstXoverState build_rdp_first_xover_state(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, int triplet_index, int fss_ub,
    std::vector<unsigned char>& fss_rdp, int xover_window,
    short xover_window_x, const Dna5XoverApi& api);

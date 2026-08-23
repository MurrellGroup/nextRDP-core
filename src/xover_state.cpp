#include "xover_state.hpp"

#include <cstddef>
#include <stdexcept>

RdpFirstXoverState build_rdp_first_xover_state(
    const RdpScanState& scan_state,
    const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const int triplet_index,
    const int fss_ub, std::vector<unsigned char>& fss_rdp,
    const int xover_window, const short xover_window_x,
    const Dna5XoverApi& api) {
    if (triplet_index < 0 || triplet_index > scan_state.analysis_list_last) {
        throw std::runtime_error("XOver received an invalid triplet index");
    }
    const int sequence_count = scan_state.next_no + 1;
    if (tree_state.tree_distance.size() !=
        static_cast<std::size_t>(sequence_count) * sequence_count) {
        throw std::runtime_error("XOver received a mis-sized tree matrix");
    }
    if (distance_state.distance.size() !=
        static_cast<std::size_t>(sequence_count) * sequence_count) {
        throw std::runtime_error("XOver received a mis-sized distance matrix");
    }

    RdpFirstXoverState state;
    for (int role = 0; role < 3; ++role) {
        state.sequences[role] =
            scan_state.analysis_list[role + triplet_index * 3];
    }
    state.xover_sequence_ub =
        scan_state.sequence_length + (xover_window_x / 2) * 2;
    state.xover_sequence.assign(
        static_cast<std::size_t>(state.xover_sequence_ub + 1) * 3, 0);
    state.informative_length = api.find_subsequence(
        state.agreement_counts.data(), fss_ub, xover_window,
        scan_state.compressed_sequence_ub, scan_state.sequence_length,
        scan_state.next_no, state.sequences[0], state.sequences[1],
        state.sequences[2],
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        state.xover_sequence_ub, state.xover_sequence.data(), fss_rdp.data());

    if (state.informative_length < xover_window * 2 ||
        state.agreement_counts[0] < xover_window / 3 ||
        state.agreement_counts[1] < xover_window / 3 ||
        state.agreement_counts[2] < xover_window / 3) {
        return state;
    }
    state.homology_length = state.informative_length - 1;
    const int first = state.sequences[0];
    const int second = state.sequences[1];
    const int third = state.sequences[2];
    const auto& distance = tree_state.tree_distance;
    const float d12 = distance[first + second * sequence_count];
    const float d13 = distance[first + third * sequence_count];
    const float d23 = distance[second + third * sequence_count];
    if (d12 >= d13 && d12 >= d23) {
        state.high_homology = 1;
    } else if (d13 >= d12 && d13 >= d23) {
        state.high_homology = 2;
    } else if (d23 >= d12 && d23 >= d13) {
        state.high_homology = 3;
    }
    state.initial_high_homology = state.high_homology;

    state.homology_ub = scan_state.sequence_length + xover_window * 2;
    state.homology.assign(
        static_cast<std::size_t>(state.homology_ub + 1) * 3, 0);
    // XOver passes Len(StrainSeq(0)) to FindSubSeqPB3, but passes the global
    // LenStrainSeq (initialized to Len(StrainSeq(0)) + 1) here.  DNA5 uses
    // this argument as its flattened row stride, so the inconsistency is
    // observable and must be preserved.
    state.homology_sequence_length = scan_state.sequence_length + 1;
    state.homology_start = api.calculate_homology(
        static_cast<short>(state.initial_high_homology),
        state.homology_sequence_length, state.homology_length,
        static_cast<short>(xover_window),
        state.xover_sequence.data(), state.homology.data());

    for (int pair = 0; pair < 3; ++pair) {
        state.average_homology[pair] =
            static_cast<double>(state.agreement_counts[pair]) /
            static_cast<double>(state.informative_length);
    }

    const auto& raw_distance = distance_state.distance;
    const float raw_d12 = raw_distance[first + second * sequence_count];
    const float raw_d13 = raw_distance[first + third * sequence_count];
    const float raw_d23 = raw_distance[second + third * sequence_count];
    auto& a1 = state.average_homology[0];
    auto& a2 = state.average_homology[1];
    auto& a3 = state.average_homology[2];
    if (a1 == a2 && a1 == a3) {
        if (raw_d12 >= raw_d13 && raw_d12 >= raw_d23) {
            if (raw_d13 > raw_d23) {
                a2 -= 0.00001;
                a3 -= 0.00002;
            } else {
                a2 -= 0.00002;
                a3 -= 0.00001;
            }
        } else if (raw_d13 >= raw_d12 && raw_d13 >= raw_d23) {
            if (raw_d12 > raw_d23) {
                a1 -= 0.00001;
                a3 -= 0.00002;
            } else {
                a1 -= 0.00002;
                a3 -= 0.00001;
            }
        } else {
            if (raw_d12 > raw_d13) {
                a1 -= 0.00001;
                a2 -= 0.00002;
            } else {
                a1 -= 0.00002;
                a2 -= 0.00001;
            }
        }
    } else if (a1 == a2) {
        if (raw_d12 > raw_d13) {
            a2 -= 0.00001;
        } else {
            a1 -= 0.00001;
        }
    } else if (a1 == a3) {
        if (raw_d12 > raw_d23) {
            a3 -= 0.00001;
        } else {
            a1 -= 0.00001;
        }
    } else if (a2 == a3) {
        if (raw_d13 > raw_d23) {
            a3 -= 0.00001;
        } else {
            a2 -= 0.00001;
        }
    }

    if (a1 >= a2 && a1 >= a3) {
        state.high_homology = 1;
        if (a2 >= a3) {
            state.med_homology = 2;
            state.low_homology = 3;
            state.active_sequence = first;
            state.active_major_parent = second;
            state.active_minor_parent = third;
            state.sequence_daughter = 0;
            state.sequence_minor = 2;
        } else {
            state.med_homology = 3;
            state.low_homology = 2;
            state.active_sequence = second;
            state.active_major_parent = first;
            state.active_minor_parent = third;
            state.sequence_daughter = 1;
            state.sequence_minor = 2;
        }
    } else if (a2 >= a1 && a2 >= a3) {
        state.high_homology = 2;
        if (a1 >= a3) {
            state.med_homology = 1;
            state.low_homology = 3;
            state.active_sequence = first;
            state.active_major_parent = third;
            state.active_minor_parent = second;
            state.sequence_daughter = 0;
            state.sequence_minor = 1;
        } else {
            state.med_homology = 3;
            state.low_homology = 1;
            state.active_sequence = third;
            state.active_major_parent = first;
            state.active_minor_parent = second;
            state.sequence_daughter = 2;
            state.sequence_minor = 1;
        }
    } else {
        state.high_homology = 3;
        if (a1 >= a2) {
            state.med_homology = 1;
            state.low_homology = 2;
            state.active_sequence = second;
            state.active_major_parent = third;
            state.active_minor_parent = first;
            state.sequence_daughter = 1;
            state.sequence_minor = 0;
        } else {
            state.med_homology = 2;
            state.low_homology = 1;
            state.active_sequence = third;
            state.active_major_parent = second;
            state.active_minor_parent = first;
            state.sequence_daughter = 2;
            state.sequence_minor = 0;
        }
    }

    state.next_position = api.find_next(
        state.homology_ub, 1, state.high_homology, state.med_homology,
        state.low_homology, state.homology_length, xover_window,
        state.homology.data());
    return state;
}

#include "xover_state.hpp"

#include <cstddef>
#include <cmath>
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

void define_rdp_first_xover_event(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& settings,
    const Dna5XoverApi& api) {
    if (state.next_position < 0) {
        state.event_position = state.next_position;
        return;
    }
    const int stride = state.homology_ub + 1;
    const int med_offset = (state.med_homology - 1) * stride;
    const int high_offset = (state.high_homology - 1) * stride;
    int search_start = 1;
    int position = state.next_position;
    int old_position = -1;
    while (position > -1 && position != old_position) {
        old_position = position;
        if (settings.circular == 1 && position == 1 &&
            state.homology[position + med_offset] >
                state.homology[position + high_offset]) {
            state.used_find_first = true;
            position = api.find_first(
                position, state.med_homology, state.high_homology,
                state.homology_length, state.homology_ub,
                state.homology.data());
            if (position < state.homology_length + 1 &&
                position > search_start) {
                search_start = position + 1;
                position = api.find_next(
                    state.homology_ub, search_start, state.high_homology,
                    state.med_homology, state.low_homology,
                    state.homology_length, xover_window,
                    state.homology.data());
                continue;
            }
            state.event_position = position;
            return;
        }

        // VB Long locals begin as zero, and XOver explicitly resets NCommon
        // and XOverLength before the DefineEventP2 branch.
        state.number_in_common = 0;
        state.event_length = 0;
        state.define_input_position = position;
        state.xover_sequence_at_define = state.xover_sequence;
        state.homology_at_define = state.homology;
        // Despite DefineEventP2's declaration order, XOver passes SeqMinorP
        // first and SeqDaughter second. Preserve that observable caller bug.
        state.event_position = api.define_event(
            state.homology_ub, settings.short_output, settings.long_winded,
            state.med_homology, state.high_homology, state.low_homology,
            settings.target, settings.circular, position, xover_window,
            scan_state.sequence_length, state.homology_length,
            state.sequence_minor, state.sequence_daughter, &state.end_flag,
            &state.event_begin, &state.event_end, &state.number_in_common,
            &state.event_length, state.xover_sequence.data(),
            state.homology.data());
        return;
    }
    state.event_position = position;
}

void calculate_rdp_first_xover_probability(
    RdpFirstXoverState& state, const RdpProbabilitySettings& settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api) {
    if (state.define_input_position < 0 || state.event_length <= 2 ||
        state.event_end == state.event_begin ||
        (state.event_end <= state.event_begin && settings.circular != 1)) {
        return;
    }

    state.probability_different =
        state.event_length - state.number_in_common;
    if (state.number_in_common <= state.probability_different * 0.8) return;

    state.probability_length = state.event_length;
    state.probability_same = state.number_in_common;
    if (state.event_length >= 170) {
        state.probability_scale =
            static_cast<double>(state.event_length) / 169.0;
        state.probability_different = static_cast<int>(std::nearbyint(
            static_cast<double>(state.probability_different) * 169.0 /
            static_cast<double>(state.event_length)));
        state.probability_length = 169;
        state.probability_same =
            state.probability_length - state.probability_different;
    }
    state.individual_probability =
        state.average_homology[state.med_homology - 1];

    bool proceed = settings.probability_file_flag != 0;
    if (!proceed) {
        const int first_stride = settings.probability_one_ub + 1;
        const int second_stride = settings.probability_two_ub + 1;
        const int category =
            static_cast<int>(state.individual_probability * 50.0);
        const auto index = static_cast<std::size_t>(state.probability_length) +
            static_cast<std::size_t>(state.probability_same) * first_stride +
            static_cast<std::size_t>(category) * first_stride * second_stride;
        if (index >= probability_estimate.size()) {
            throw std::runtime_error("ProbEstimate lookup exceeds its bounds");
        }
        state.probability_prefilter_value = probability_estimate[index];
        proceed =
            state.probability_prefilter_value < settings.lowest_probability;
    }
    if (!proceed) return;

    state.probability_tested = true;
    if (state.event_length <= settings.fact_three_ub) {
        state.used_probability_p2 = true;
        state.event_probability = api.probability_p2(
            fact_three.data(), settings.fact_three_ub,
            state.probability_length, state.probability_same,
            state.individual_probability, state.homology_length);
    } else {
        state.event_probability = api.probability_p(
            fact.data(), state.probability_length, state.probability_same,
            state.individual_probability, state.homology_length);
    }
}

void continue_rdp_xover_to_first_probability(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api) {
    while (!state.probability_tested) {
        int position = state.event_position;
        if (state.end_flag == 1) {
            state.end_flag = 0;
            position = state.homology_length;
        }
        if (position >= state.homology_length + 1 ||
            position <= state.define_input_position) {
            return;
        }
        const int search_start = position + 1;
        position = api.find_next(
            state.homology_ub, search_start, state.high_homology,
            state.med_homology, state.low_homology, state.homology_length,
            xover_window, state.homology.data());
        if (position < 0) return;

        state.number_in_common = 0;
        state.event_length = 0;
        state.define_input_position = position;
        state.event_position = api.define_event(
            state.homology_ub, xover_settings.short_output,
            xover_settings.long_winded, state.med_homology,
            state.high_homology, state.low_homology, xover_settings.target,
            xover_settings.circular, position, xover_window,
            scan_state.sequence_length, state.homology_length,
            state.sequence_minor, state.sequence_daughter, &state.end_flag,
            &state.event_begin, &state.event_end, &state.number_in_common,
            &state.event_length, state.xover_sequence.data(),
            state.homology.data());
        calculate_rdp_first_xover_probability(
            state, probability_settings, probability_estimate, fact_three,
            fact, api);
    }
}

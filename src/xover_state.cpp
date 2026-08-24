#include "xover_state.hpp"

#include <cstddef>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

void assign_active_roles(RdpFirstXoverState& state) {
    const int first = state.sequences[0];
    const int second = state.sequences[1];
    const int third = state.sequences[2];
    if (state.high_homology == 1 && state.med_homology == 2) {
        state.active_sequence = first;
        state.active_major_parent = second;
        state.active_minor_parent = third;
        state.sequence_daughter = 0;
        state.sequence_minor = 2;
    } else if (state.high_homology == 1 && state.med_homology == 3) {
        state.active_sequence = second;
        state.active_major_parent = first;
        state.active_minor_parent = third;
        state.sequence_daughter = 1;
        state.sequence_minor = 2;
    } else if (state.high_homology == 2 && state.med_homology == 1) {
        state.active_sequence = first;
        state.active_major_parent = third;
        state.active_minor_parent = second;
        state.sequence_daughter = 0;
        state.sequence_minor = 1;
    } else if (state.high_homology == 2 && state.med_homology == 3) {
        state.active_sequence = third;
        state.active_major_parent = first;
        state.active_minor_parent = second;
        state.sequence_daughter = 2;
        state.sequence_minor = 1;
    } else if (state.high_homology == 3 && state.med_homology == 1) {
        state.active_sequence = second;
        state.active_major_parent = third;
        state.active_minor_parent = first;
        state.sequence_daughter = 1;
        state.sequence_minor = 0;
    } else if (state.high_homology == 3 && state.med_homology == 2) {
        state.active_sequence = third;
        state.active_major_parent = second;
        state.active_minor_parent = first;
        state.sequence_daughter = 2;
        state.sequence_minor = 0;
    }
}

int vb_clng(const double value) {
    return static_cast<int>(std::nearbyint(value));
}

void map_and_centre_breakpoints(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& settings,
    int& beginning, int& ending, int& beginning_warning,
    int& ending_warning) {
    const int original_beginning = beginning;
    int original_ending = ending;
    if (settings.circular == 0) {
        if (beginning == 1) {
            if (state.xdiffpos[beginning] < settings.target) {
                beginning = 1;
            } else {
                beginning = state.xdiffpos[beginning];
            }
        } else {
            beginning = state.xdiffpos[beginning];
        }
    } else {
        beginning = state.xdiffpos[beginning];
    }

    if (ending == state.homology_length && settings.circular == 0) {
        if (settings.short_output == 0 || settings.short_output == 6 ||
            settings.short_output == 10) {
            ending = scan_state.sequence_length;
        } else {
            ending = state.xdiffpos[state.homology_length];
        }
        original_ending = state.homology_length;
    } else {
        if (ending >= state.homology_length) {
            if (settings.long_winded == 0) {
                state.xdiffpos[ending] = 0;
            } else {
                ending = state.homology_length + 1;
            }
        } else if (ending < 1) {
            ending = state.homology_length + 1;
        }
        original_ending = ending;
        ending = state.xdiffpos[ending];
    }
    if (ending == 0) {
        if (settings.short_output == 0 || settings.short_output == 6 ||
            settings.short_output == 10) {
            ending = scan_state.sequence_length;
        } else {
            ending = state.xdiffpos[state.homology_length];
        }
        original_ending = state.homology_length;
    }

    // The compressed PB3 path sets XPDDone before the scan. Consequently the
    // source calls CentreBP with OBE=oEN=0, not the saved variable-site
    // coordinates. SEventNumber is zero during this initial scan, so the
    // missing-data branches are intentionally absent here.
    const int mapped_beginning = beginning;
    if (state.xposdiff[beginning] - 1 > 0) {
        beginning -= vb_clng(
            ((beginning - state.xdiffpos[state.xposdiff[beginning] - 1]) /
                2.0) -
            0.1);
    } else {
        beginning -= vb_clng(
            ((beginning + scan_state.sequence_length -
                 state.xdiffpos[state.homology_length]) /
                2.0) -
            0.1);
    }
    if (beginning == 0) {
        beginning = 1;
    } else if (beginning < 1) {
        beginning = settings.circular == 0 ?
            1 : scan_state.sequence_length + beginning;
    }
    if (settings.circular == 0 &&
        state.xposdiff[beginning] < xover_window) {
        beginning_warning = 1;
    }

    state.xposdiff[scan_state.sequence_length] = state.homology_length;
    if (state.xposdiff[ending] + 1 <= state.homology_length) {
        ending += vb_clng(
            ((state.xdiffpos[state.xposdiff[ending] + 1] - ending) / 2.0) -
            0.1);
    } else {
        ending += vb_clng(
            ((state.xdiffpos[1] +
                 (scan_state.sequence_length - ending)) /
                2.0) -
            0.1);
    }
    if (ending > scan_state.sequence_length) {
        ending = settings.circular == 0 ?
            scan_state.sequence_length : ending - scan_state.sequence_length;
    }
    if (settings.circular == 0 &&
        state.xposdiff[ending] > state.homology_length - xover_window) {
        ending_warning = 1;
    }

    (void)original_beginning;
    (void)original_ending;
    (void)mapped_beginning;
}

std::array<int, 3> choose_storage_roles(
    const RdpFirstXoverState& state, const RdpRawEventState& events,
    const std::vector<double>& store_lpv, const int store_lpv_ub) {
    const int first = state.sequences[0];
    const int second = state.sequences[1];
    const int third = state.sequences[2];
    const auto count = [&events](const int sequence) {
        return events.current_xover[sequence];
    };
    if (count(first) < count(second) && count(first) < count(third)) {
        return {first, second, third};
    }
    if (count(second) < count(first) && count(second) < count(third)) {
        return {second, first, third};
    }
    if (count(third) < count(first) && count(third) < count(second)) {
        return {third, first, second};
    }
    const auto lpv = [&store_lpv, store_lpv_ub](const int sequence) {
        const auto index = static_cast<std::size_t>(sequence) *
            static_cast<std::size_t>(store_lpv_ub + 1);
        if (index >= store_lpv.size()) {
            throw std::runtime_error("StoreLPV lookup exceeds its bounds");
        }
        return store_lpv[index];
    };
    if (lpv(first) >= lpv(second) && lpv(first) >= lpv(third)) {
        return {first, second, third};
    }
    if (lpv(second) >= lpv(first) && lpv(second) >= lpv(third)) {
        return {second, first, third};
    }
    return {third, first, second};
}

}  // namespace

RdpFirstXoverState build_rdp_first_xover_state(
    const RdpScanState& scan_state,
    const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const int triplet_index,
    const int fss_ub, std::vector<unsigned char>& fss_rdp,
    const int xover_window, const short xover_window_x,
    const Dna5XoverApi& api,
    const std::array<int, 3>* explicit_sequences) {
    if (explicit_sequences == nullptr &&
        (triplet_index < 0 || triplet_index > scan_state.analysis_list_last)) {
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
        state.sequences[role] = explicit_sequences == nullptr
            ? scan_state.analysis_list[role + triplet_index * 3]
            : (*explicit_sequences)[role];
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
        if (position == state.old_find_position) return;
        state.old_find_position = position;
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
    state.probability_tested = false;
    state.used_probability_p2 = false;
    state.probability_length = 0;
    state.probability_same = 0;
    state.probability_different = 0;
    state.probability_scale = 1.0;
    state.individual_probability = 0.0;
    state.probability_prefilter_value = 0.0;
    state.event_probability = 0.0;
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
        if (position == state.old_find_position) return;
        state.old_find_position = position;

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

bool apply_rdp_probability_cutoff(
    RdpFirstXoverState& state, const RdpProbabilitySettings& settings) {
    state.significant_event = false;
    state.adjusted_event_probability = state.event_probability;
    if (!state.probability_tested || state.event_probability >= 0.5) {
        return false;
    }
    if (state.probability_scale != 1.0) {
        if (state.adjusted_event_probability > 0.0) {
            state.adjusted_event_probability = std::pow(
                state.adjusted_event_probability, state.probability_scale);
        } else {
            state.adjusted_event_probability = 0.05;
        }
    }
    const double minimum_probability = std::pow(10.0, -300.0);
    if (state.adjusted_event_probability < minimum_probability) {
        state.adjusted_event_probability = minimum_probability;
    }
    if (settings.mc_flag == 0) {
        state.adjusted_event_probability *= settings.mc_correction;
    }
    state.significant_event = state.adjusted_event_probability <
            settings.lowest_probability &&
        state.adjusted_event_probability > 0.0;
    return state.significant_event;
}

void build_rdp_first_position_maps(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int fss_ub, const int xover_window,
    std::vector<unsigned char>& fss_rdp, const Dna5XoverApi& api) {
    if (!state.significant_event) return;
    state.xdiffpos.assign(
        static_cast<std::size_t>(scan_state.sequence_length + 201), 0);
    state.xposdiff.assign(
        static_cast<std::size_t>(scan_state.sequence_length + 201), 0);
    state.position_map_result = api.find_subsequence_with_positions(
        state.agreement_counts.data(), fss_ub, xover_window,
        scan_state.compressed_sequence_ub, scan_state.sequence_length,
        scan_state.next_no, state.sequences[0], state.sequences[1],
        state.sequences[2],
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        state.xover_sequence_ub, state.xover_sequence.data(),
        state.xdiffpos.data(), state.xposdiff.data(), fss_rdp.data());
}

bool advance_rdp_role_cycle(RdpFirstXoverState& state, const int do_all) {
    if (state.find_cycle == 0) {
        const int temp = state.med_homology;
        state.med_homology = state.low_homology;
        state.low_homology = temp;
    } else if (state.find_cycle == 1) {
        if (state.average_homology[state.high_homology - 1] >= 0.7 &&
            do_all != 1) {
            return false;
        }
        const int temp = state.high_homology;
        state.high_homology = state.low_homology;
        state.low_homology = state.med_homology;
        state.med_homology = temp;
    } else {
        return false;
    }
    assign_active_roles(state);
    ++state.find_cycle;
    return true;
}

void scan_rdp_current_roles_to_first_probability(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api) {
    state.probability_tested = false;
    int search_start = 1;
    int position = api.find_next(
        state.homology_ub, search_start, state.high_homology,
        state.med_homology, state.low_homology, state.homology_length,
        xover_window, state.homology.data());
    while (position > -1 && position != state.old_find_position) {
        state.old_find_position = position;
        const int stride = state.homology_ub + 1;
        if (xover_settings.circular == 1 && position == 1 &&
            state.homology[position + (state.med_homology - 1) * stride] >
                state.homology[position +
                    (state.high_homology - 1) * stride]) {
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
            return;
        }

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
        if (state.probability_tested) return;
        continue_rdp_xover_to_first_probability(
            state, scan_state, xover_window, xover_settings,
            probability_settings, probability_estimate, fact_three, fact,
            api);
        return;
    }
}

RdpRawEventState scan_rdp_redo_triplets(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const std::vector<unsigned char>& redo,
    std::vector<unsigned char>& fss_rdp, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int fss_ub, const int xover_window,
    const short xover_window_x, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api, const int do_all,
    const RdpRawEventState* initial_events,
    const std::array<int, 3>* explicit_sequences) {
    RdpRawEventState events = initial_events == nullptr
        ? RdpRawEventState{} : *initial_events;
    if (initial_events == nullptr) {
        events.current_xover.assign(scan_state.next_no + 1, 0);
        events.xover_list.resize(scan_state.next_no + 1);
    }
    for (int triplet = 0; triplet <= scan_state.analysis_list_last;
         ++triplet) {
        if (static_cast<std::size_t>(triplet) >= redo.size() ||
            redo[triplet] != 1) {
            continue;
        }
        ++events.scanned_triplets;
        auto state = build_rdp_first_xover_state(
            scan_state, distance_state, tree_state, triplet, fss_ub, fss_rdp,
            xover_window, xover_window_x, api, explicit_sequences);
        if (state.informative_length < xover_window * 2 ||
            state.agreement_counts[0] < xover_window / 3 ||
            state.agreement_counts[1] < xover_window / 3 ||
            state.agreement_counts[2] < xover_window / 3) {
            api.clean_xover_sequence(
                state.homology_length + xover_window * 2, xover_window,
                state.xover_sequence_ub, state.xover_sequence.data());
            continue;
        }

        bool position_maps_built = false;
        int old_x = -1;
        for (int role_cycle = 0; role_cycle < 3; ++role_cycle) {
            int next_position = 1;
            while (true) {
                int position = api.find_next(
                    state.homology_ub, next_position, state.high_homology,
                    state.med_homology, state.low_homology,
                    state.homology_length, xover_window,
                    state.homology.data());
                if (position <= -1 || position == old_x) break;
                old_x = position;
                const int stride = state.homology_ub + 1;
                if (xover_settings.circular == 1 && position == 1 &&
                    state.homology[position +
                        (state.med_homology - 1) * stride] >
                        state.homology[position +
                            (state.high_homology - 1) * stride]) {
                    position = api.find_first(
                        position, state.med_homology, state.high_homology,
                        state.homology_length, state.homology_ub,
                        state.homology.data());
                } else {
                    state.number_in_common = 0;
                    state.event_length = 0;
                    state.define_input_position = position;
                    state.event_position = api.define_event(
                        state.homology_ub, xover_settings.short_output,
                        xover_settings.long_winded, state.med_homology,
                        state.high_homology, state.low_homology,
                        xover_settings.target, xover_settings.circular,
                        position, xover_window, scan_state.sequence_length,
                        state.homology_length, state.sequence_minor,
                        state.sequence_daughter, &state.end_flag,
                        &state.event_begin, &state.event_end,
                        &state.number_in_common, &state.event_length,
                        state.xover_sequence.data(), state.homology.data());
                    position = state.event_position;
                    calculate_rdp_first_xover_probability(
                        state, probability_settings, probability_estimate,
                        fact_three, fact, api);
                    if (apply_rdp_probability_cutoff(
                            state, probability_settings)) {
                        ++events.significant_candidates;
                        if (!position_maps_built) {
                            build_rdp_first_position_maps(
                                state, scan_state, fss_ub, xover_window,
                                fss_rdp, api);
                            position_maps_built = true;
                        }
                        int beginning = state.event_begin;
                        int ending = state.event_end;
                        int beginning_warning = 0;
                        int ending_warning = 0;
                        map_and_centre_breakpoints(
                            state, scan_state, xover_window, xover_settings,
                            beginning, ending, beginning_warning,
                            ending_warning);
                        const auto storage = choose_storage_roles(
                            state, events, store_lpv, store_lpv_ub);
                        RdpRawEvent event;
                        event.beginning = beginning;
                        event.ending = ending;
                        event.major_parent =
                            static_cast<std::int16_t>(storage[1]);
                        event.minor_parent =
                            static_cast<std::int16_t>(storage[2]);
                        event.daughter =
                            static_cast<std::int16_t>(storage[0]);
                        event.program_flag = 0;
                        event.probability = state.adjusted_event_probability;
                        if (beginning_warning == 1 && ending_warning == 1) {
                            event.sbp_flag = 3;
                        } else if (beginning_warning == 1) {
                            event.sbp_flag = 1;
                        } else if (ending_warning == 1) {
                            event.sbp_flag = 2;
                        }
                        events.xover_list[storage[0]].push_back(event);
                        events.current_xover[storage[0]] =
                            static_cast<std::int16_t>(
                                events.xover_list[storage[0]].size());
                    }
                }
                if (state.end_flag == 1) {
                    state.end_flag = 0;
                    position = state.homology_length;
                }
                if (position < state.homology_length + 1 &&
                    position > next_position) {
                    next_position = position + 1;
                } else {
                    break;
                }
            }
            if (!advance_rdp_role_cycle(state, do_all)) break;
        }
        api.clean_xover_sequence(
            state.homology_length + xover_window * 2, xover_window,
            state.xover_sequence_ub, state.xover_sequence.data());
    }
    return events;
}

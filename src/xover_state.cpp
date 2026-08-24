#include "xover_state.hpp"

#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>
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

struct MappedBreakpoints {
    int beginning = 0;
    int ending = 0;
    int original_beginning = 0;
    int original_ending = 0;
};

MappedBreakpoints map_breakpoints(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const RdpXoverSettings& settings, int beginning, int ending) {
    MappedBreakpoints mapped;
    mapped.original_beginning = beginning;
    mapped.original_ending = ending;
    if (settings.circular == 0) {
        if (beginning == 1) {
            if (state.xdiffpos[beginning] < settings.target) {
                beginning = 1;
                mapped.original_beginning = 1;
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
        mapped.original_ending = state.homology_length;
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
        mapped.original_ending = ending;
        ending = state.xdiffpos[ending];
    }
    if (ending == 0) {
        if (settings.short_output == 0 || settings.short_output == 6 ||
            settings.short_output == 10) {
            ending = scan_state.sequence_length;
        } else {
            ending = state.xdiffpos[state.homology_length];
        }
        mapped.original_ending = state.homology_length;
    }
    mapped.beginning = beginning;
    mapped.ending = ending;
    return mapped;
}

void centre_mapped_breakpoints(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& settings,
    int& beginning, int& ending, int& beginning_warning,
    int& ending_warning, const int sequence_event_number,
    const std::vector<unsigned char>* missing_data,
    const int original_beginning = 0, const int original_ending = 0,
    const bool use_compress = false) {
    // XOver has two distinct CentreBP call forms. The compressed initial scan
    // has XPDDone set and passes (0, 0), while FinalTrim/rescans use the plain
    // FindSubSeqP path and pass the saved variable-site OBE/oEN values.
    const bool use_original_positions = sequence_event_number > 0 &&
        !use_compress &&
        (original_beginning > 0 || original_ending > 0);
    if (use_original_positions && original_beginning - 1 > 0) {
        beginning -= vb_clng(
            ((beginning - state.xdiffpos[original_beginning - 1]) / 2.0) -
            0.1);
    } else if (use_original_positions) {
        beginning -= vb_clng(
            ((beginning + scan_state.sequence_length -
                 state.xdiffpos[state.homology_length]) / 2.0) - 0.1);
    } else if (state.xposdiff[beginning] - 1 > 0) {
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
    const auto has_missing = [&](const int position) {
        if (missing_data == nullptr || sequence_event_number <= 0) return false;
        const int stride = scan_state.sequence_length + 1;
        for (const int sequence : state.sequences) {
            const auto index = static_cast<std::size_t>(position) +
                static_cast<std::size_t>(sequence) * stride;
            if (index < missing_data->size() && (*missing_data)[index] == 1) {
                return true;
            }
        }
        return false;
    };
    if (has_missing(beginning)) {
        beginning_warning = 1;
        for (int step = 1; step <= scan_state.sequence_length; ++step) {
            ++beginning;
            if (beginning > scan_state.sequence_length) {
                beginning = 1;
                if (settings.circular == 0) break;
            }
            if (!has_missing(beginning)) break;
        }
    }
    if (settings.circular == 0) {
        if (!use_original_positions &&
            state.xposdiff[beginning] < xover_window) {
            beginning_warning = 1;
        } else if (use_original_positions) {
            for (int variable = 1; variable <= xover_window; ++variable) {
                if (variable < static_cast<int>(state.xdiffpos.size()) &&
                    state.xdiffpos[variable] > beginning) {
                    beginning_warning = 1;
                }
            }
        }
    }

    state.xposdiff[scan_state.sequence_length] = state.homology_length;
    if (use_original_positions && original_ending + 1 <=
            state.homology_length) {
        ending += vb_clng(
            ((state.xdiffpos[original_ending + 1] - ending) / 2.0) - 0.1);
    } else if (use_original_positions) {
        ending += vb_clng(
            ((state.xdiffpos[1] +
                 (scan_state.sequence_length - ending)) / 2.0) - 0.1);
    } else if (state.xposdiff[ending] + 1 <= state.homology_length) {
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
    if (has_missing(ending)) {
        ending_warning = 1;
        for (int step = 1; step <= scan_state.sequence_length; ++step) {
            --ending;
            if (ending < 1) {
                ending = scan_state.sequence_length;
                if (settings.circular == 0) break;
            }
            if (!has_missing(ending)) break;
        }
    }
    if (settings.circular == 0) {
        if (!use_original_positions &&
            state.xposdiff[ending] > state.homology_length - xover_window) {
            ending_warning = 1;
        } else if (use_original_positions) {
            for (int variable = state.homology_length;
                 variable >= state.homology_length - xover_window;
                 --variable) {
                if (variable >= 0 &&
                    variable < static_cast<int>(state.xdiffpos.size()) &&
                    state.xdiffpos[variable] < ending) {
                    ending_warning = 1;
                }
            }
        }
    }

}

void map_and_centre_breakpoints(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& settings,
    int& beginning, int& ending, int& beginning_warning,
    int& ending_warning, const int sequence_event_number,
    const std::vector<unsigned char>* missing_data) {
    const auto mapped = map_breakpoints(
        state, scan_state, settings, beginning, ending);
    beginning = mapped.beginning;
    ending = mapped.ending;
    centre_mapped_breakpoints(
        state, scan_state, xover_window, settings, beginning, ending,
        beginning_warning, ending_warning, sequence_event_number,
        missing_data, mapped.original_beginning, mapped.original_ending,
        false);
}

int check_split(
    const int step, const int sequence_length, const int beginning,
    const int ending, const std::array<int, 3>& sequences, int& split,
    const std::vector<unsigned char>& missing_data) {
    const int stride = sequence_length + 1;
    const auto missing = [&](const int position) {
        for (const int sequence : sequences) {
            if (missing_data[position + sequence * stride] == 1) return true;
        }
        return false;
    };
    int position = beginning;
    if (beginning < ending) {
        for (position = beginning; position <= ending; position += step) {
            if (missing(position)) {
                split = 1;
                break;
            }
        }
    } else {
        for (position = beginning; position <= sequence_length;
             position += step) {
            if (missing(position)) {
                split = 1;
                break;
            }
        }
        if (split == 0) {
            for (position = 1; position <= ending; position += step) {
                if (missing(position)) {
                    split = 1;
                    break;
                }
            }
        }
    }
    return position;
}

int find_missing(
    const int sequence_length, const std::array<int, 3>& sequences,
    const int split_position, const int ending,
    const std::vector<unsigned char>& missing_data) {
    const int stride = sequence_length + 1;
    const auto missing = [&](const int position) {
        for (const int sequence : sequences) {
            if (missing_data[position + sequence * stride] == 1) return true;
        }
        return false;
    };
    int position = ending;
    if (split_position < ending) {
        for (position = ending; position >= split_position; --position) {
            if (missing(position)) break;
        }
    } else {
        bool found = false;
        for (position = ending; position >= 1; --position) {
            if (missing(position)) {
                found = true;
                break;
            }
        }
        if (!found) {
            for (position = sequence_length; position >= split_position;
                 --position) {
                if (missing(position)) break;
            }
        }
    }
    return position;
}

int split_event(
    const int xover_window, const int sequence_length,
    const int informative_length, const int sequence_daughter,
    const int sequence_minor, const int beginning, const int ending,
    int& event_length, int& number_in_common,
    const std::vector<int>& xposdiff,
    const std::vector<char>& xover_sequence) {
    int comparison = 0;
    if ((sequence_daughter == 0 && sequence_minor == 2) ||
        (sequence_daughter == 2 && sequence_minor == 0)) {
        comparison = 1;
    } else if ((sequence_daughter == 2 && sequence_minor == 1) ||
               (sequence_daughter == 1 && sequence_minor == 2)) {
        comparison = 2;
    }
    const int offset = xover_window +
        comparison * (sequence_length + 1 + xover_window * 2);
    number_in_common = 0;
    if (ending >= beginning) {
        for (int position = xposdiff[beginning];
             position <= xposdiff[ending]; ++position) {
            if (xover_sequence[position + offset] != 0) ++number_in_common;
        }
        event_length = xposdiff[ending] - xposdiff[beginning] + 1;
    } else {
        for (int position = xposdiff[beginning];
             position <= informative_length; ++position) {
            if (xover_sequence[position + offset] != 0) ++number_in_common;
        }
        for (int position = 1; position <= xposdiff[ending]; ++position) {
            if (xover_sequence[position + offset] != 0) ++number_in_common;
        }
        event_length = informative_length - xposdiff[beginning] + 1 +
            xposdiff[ending];
    }
    return event_length - number_in_common;
}

int check_ends(
    const RdpFirstXoverState& state, const RdpScanState& scan_state,
    const int xover_window, const RdpXoverSettings& settings,
    const int check_ending, const int original_beginning,
    const int original_ending,
    const std::vector<unsigned char>& missing_data,
    const int source_beginning = 0, const int source_ending = 0) {
    // CheckEndsVB has a second branch for the plain XOver path: OBE/oEN are
    // the variable-site coordinates saved before mapping to full sequence
    // positions. Preserve that branch separately from the (0,0) compressed
    // path used by the initial scan and split fragments.
    if (source_beginning > 0 || source_ending > 0) {
        int beginning = original_beginning;
        int ending = original_ending;
        const int saved_beginning = beginning;
        const int saved_ending = ending;
        const int stride = scan_state.sequence_length + 1;
        const auto missing = [&](const int position) {
            if (position < 0 || position > scan_state.sequence_length) {
                return false;
            }
            for (const int sequence : state.sequences) {
                const auto index = static_cast<std::size_t>(position) +
                    static_cast<std::size_t>(sequence) * stride;
                if (index < missing_data.size() && missing_data[index] == 1) {
                    return true;
                }
            }
            return false;
        };
        int warning = 0;
        if (check_ending == 0) {
            int target = 1;
            if (source_beginning - xover_window > 0) {
                target = state.xdiffpos[source_beginning - xover_window];
            } else if (settings.circular == 1 &&
                       source_beginning - xover_window +
                               state.homology_length >= 0) {
                target = state.xdiffpos[source_beginning - xover_window +
                                        state.homology_length];
            }
            if (source_beginning + xover_window <
                    state.homology_length) {
                beginning = state.xdiffpos[source_beginning + xover_window];
            } else if (settings.circular == 1) {
                beginning = state.xdiffpos[source_beginning + xover_window -
                                           state.homology_length];
            } else {
                beginning = scan_state.sequence_length;
            }
            if (target < beginning) {
                for (int position = target; position <= beginning; ++position) {
                    if (missing(position)) { warning = 1; break; }
                }
            } else {
                for (int position = target;
                     position <= scan_state.sequence_length - 1; ++position) {
                    if (missing(position)) { warning = 1; break; }
                }
                if (warning == 0 && settings.circular == 0 &&
                    (missing(scan_state.sequence_length) || missing(1))) {
                    warning = 1;
                }
                if (warning == 0) {
                    for (int position = 2; position <= beginning; ++position) {
                        if (missing(position)) { warning = 1; break; }
                    }
                }
            }
            if (settings.circular == 0 && source_beginning < xover_window) {
                warning = 1;
            }
        } else {
            int target = scan_state.sequence_length;
            if (source_ending + xover_window < state.homology_length) {
                target = state.xdiffpos[source_ending + xover_window];
            } else if (settings.circular == 1) {
                target = state.xdiffpos[source_ending + xover_window -
                                        state.homology_length];
            }
            if (source_ending - xover_window > 0) {
                ending = state.xdiffpos[source_ending - xover_window];
            } else if (settings.circular == 1 &&
                       source_ending - xover_window +
                               state.homology_length >= 0) {
                ending = state.xdiffpos[source_ending - xover_window +
                                        state.homology_length];
            } else {
                ending = 1;
            }
            if (ending < 1) {
                ending = settings.circular == 1
                    ? ending + scan_state.sequence_length
                    : scan_state.sequence_length;
            }
            if (target > ending) {
                for (int position = ending; position <= target; ++position) {
                    if (missing(position)) { warning = 1; break; }
                }
            } else {
                for (int position = ending;
                     position <= scan_state.sequence_length - 1; ++position) {
                    if (missing(position)) { warning = 1; break; }
                }
                if (warning == 0 && settings.circular == 0 &&
                    (missing(scan_state.sequence_length) || missing(1))) {
                    warning = 1;
                }
                if (warning == 0) {
                    for (int position = 2; position <= target; ++position) {
                        if (missing(position)) { warning = 1; break; }
                    }
                }
            }
            if (settings.circular == 0 &&
                source_ending + xover_window > state.homology_length) {
                warning = 1;
            }
        }
        if (std::getenv("RDP_TRACE_SBP") != nullptr &&
            source_beginning <= 1 && source_ending > 0) {
            std::cerr << "checkends-trace ch=" << check_ending
                      << " src=" << source_beginning << ':' << source_ending
                      << " win=" << xover_window
                      << " circ=" << settings.circular
                      << " tmp=" << beginning << ':' << ending
                      << " warning=" << warning;
            int first_missing = -1;
            for (int position = 1;
                 position <= scan_state.sequence_length && first_missing < 0;
                 ++position) {
                if (missing(position)) first_missing = position;
            }
            std::cerr << " firstmiss=" << first_missing << '\n';
        }
        (void)saved_beginning;
        (void)saved_ending;
        return warning;
    }
    int beginning = original_beginning;
    int ending = original_ending;
    int cycle = 0;
    if (state.xposdiff[ending] == 0) {
        while (state.xposdiff[ending] == 0) {
            --ending;
            if (ending < 1) {
                ending = scan_state.sequence_length;
                if (cycle == 1) return 1;
                cycle = 1;
            }
        }
    }
    if (state.xposdiff[beginning] == 0) {
        while (state.xposdiff[beginning] == 0) {
            ++beginning;
            if (beginning > scan_state.sequence_length) beginning = 1;
        }
    }

    const int stride = scan_state.sequence_length + 1;
    const auto missing = [&](const int position) {
        for (const int sequence : state.sequences) {
            if (missing_data[position + sequence * stride] == 1) return true;
        }
        return false;
    };
    int warning = 0;
    if (check_ending == 0) {
        int target = 1;
        if (state.xposdiff[beginning] - xover_window > 0) {
            target = state.xdiffpos[
                state.xposdiff[beginning] - xover_window];
        } else if (settings.circular == 1 &&
                   state.xposdiff[beginning] - xover_window +
                           state.homology_length >=
                       0) {
            target = state.xdiffpos[
                state.xposdiff[beginning] - xover_window +
                state.homology_length];
        }
        if (target < beginning) {
            for (int position = target; position <= beginning; ++position) {
                if (missing(position)) {
                    warning = 1;
                    break;
                }
            }
        } else {
            for (int position = target;
                 position < scan_state.sequence_length; ++position) {
                if (missing(position)) {
                    warning = 1;
                    break;
                }
            }
            if (warning == 0 && settings.circular == 0 &&
                (missing(scan_state.sequence_length) || missing(1))) {
                warning = 1;
            }
            if (warning == 0) {
                for (int position = 2; position <= beginning; ++position) {
                    if (missing(position)) {
                        warning = 1;
                        break;
                    }
                }
            }
        }
    } else {
        int target = scan_state.sequence_length;
        if (state.xposdiff[ending] + xover_window <
            state.homology_length) {
            target = state.xdiffpos[
                state.xposdiff[ending] + xover_window];
        } else if (settings.circular == 1) {
            target = state.xdiffpos[
                state.xposdiff[ending] + xover_window -
                state.homology_length];
        }
        if (state.xposdiff[ending] - xover_window > 0) {
            ending = state.xdiffpos[
                state.xposdiff[ending] - xover_window];
        } else if (settings.circular == 1) {
            if (state.xposdiff[ending] - xover_window +
                    state.homology_length >=
                0) {
                ending = state.xdiffpos[
                    state.xposdiff[ending] - xover_window +
                    state.homology_length];
            } else {
                ending = 1;
            }
        } else {
            ending = 1;
        }
        if (ending < 1) {
            ending = settings.circular == 1
                ? ending + scan_state.sequence_length
                : scan_state.sequence_length;
        }
        if (target > ending) {
            for (int position = ending; position <= target; ++position) {
                if (missing(position)) {
                    warning = 1;
                    break;
                }
            }
        } else {
            for (int position = ending;
                 position <= scan_state.sequence_length - 1; ++position) {
                if (missing(position)) {
                    warning = 1;
                    break;
                }
            }
            if (warning == 0 && settings.circular == 0 &&
                (missing(scan_state.sequence_length) || missing(1))) {
                warning = 1;
            }
            if (warning == 0) {
                for (int position = 2; position <= target; ++position) {
                    if (missing(position)) {
                        warning = 1;
                        break;
                    }
                }
            }
        }
    }
    if (settings.circular == 0) {
        if (check_ending == 0) {
            if (state.xposdiff[original_beginning] < xover_window) {
                warning = 1;
            }
        } else if (state.xposdiff[original_ending] + xover_window >
                   state.homology_length) {
            warning = 1;
        }
    }
    if (std::getenv("RDP_TRACE_SBP") != nullptr &&
        (original_beginning == 1 || original_beginning == 7)) {
        int first_missing = -1;
        int first_missing_from_target = -1;
        for (int position = 1;
             position <= scan_state.sequence_length && first_missing < 0;
             ++position) {
            if (missing(position)) first_missing = position;
        }
        const int mapped_begin = (original_beginning >= 0 &&
                                  original_beginning <
                                      static_cast<int>(state.xposdiff.size()))
            ? state.xposdiff[original_beginning] : -1;
        const int mapped_end = (original_ending >= 0 &&
                                original_ending <
                                    static_cast<int>(state.xposdiff.size()))
            ? state.xposdiff[original_ending] : -1;
        const int normalized_begin = (beginning >= 0 &&
                                     beginning < static_cast<int>(state.xposdiff.size()))
            ? state.xposdiff[beginning] : -1;
        const int normalized_end = (ending >= 0 &&
                                   ending < static_cast<int>(state.xposdiff.size()))
            ? state.xposdiff[ending] : -1;
        const int xpd_value = check_ending == 0 ? normalized_begin : normalized_end;
        int target_debug = -1;
        if (check_ending == 0) {
            if (xpd_value - xover_window > 0) {
                target_debug = state.xdiffpos[xpd_value - xover_window];
            } else if (settings.circular == 1 &&
                       xpd_value - xover_window + state.homology_length >= 0) {
                target_debug = state.xdiffpos[
                    xpd_value - xover_window + state.homology_length];
            } else {
                target_debug = 1;
            }
        } else {
            if (xpd_value + xover_window < state.homology_length) {
                target_debug = state.xdiffpos[xpd_value + xover_window];
            } else if (settings.circular == 1) {
                target_debug = state.xdiffpos[
                    xpd_value + xover_window - state.homology_length];
            } else {
                target_debug = scan_state.sequence_length;
            }
        }
        for (int position = std::max(1, target_debug);
             position <= scan_state.sequence_length &&
             first_missing_from_target < 0; ++position) {
            if (missing(position)) first_missing_from_target = position;
        }
        std::cerr << "checkends-zero ch=" << check_ending
                  << " orig=" << original_beginning << ':'
                  << original_ending << " win=" << xover_window
                  << " circ=" << settings.circular << " tmp="
                  << beginning << ':' << ending << " warning=" << warning
                  << " firstmiss=" << first_missing << " xpd="
                  << mapped_begin << ':' << mapped_end << " norm="
                  << normalized_begin << ':' << normalized_end << " target="
                  << target_debug << " from-target="
                  << first_missing_from_target << " edge="
                  << missing(scan_state.sequence_length) << ':' << missing(1)
                  << " hom=" << state.homology_length << '\n';
    }
    return warning;
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
    const std::array<int, 3>* explicit_sequences,
    const int sequence_event_number, const bool use_compress) {
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
    state.xdiffpos.assign(
        static_cast<std::size_t>(scan_state.sequence_length + 201), 0);
    state.xposdiff.assign(
        static_cast<std::size_t>(scan_state.sequence_length + 201), 0);
    if (sequence_event_number > 0 && !use_compress) {
        // FinalTrim calls XOver with UseCompress=0.  That is the source's
        // plain FindSubSeqP path, which consumes SeqNum directly and produces
        // the uncompressed XOverSeqNumW/XDiffPos/XPosDiff buffers.
        std::vector<short> xover_sequence_num(
            static_cast<std::size_t>(scan_state.sequence_length + 1) * 3, 0);
        std::vector<short> spacer_sequences(1, 0);
        std::vector<short> valid_spacer(1, 0);
        state.informative_length = api.find_subsequence_plain(
            state.agreement_counts.data(), 0, 0, xover_window,
            scan_state.sequence_length + 1, scan_state.next_no,
            state.sequences[0], state.sequences[1], state.sequences[2], 0,
            const_cast<short*>(scan_state.sequence_data.data()),
            xover_sequence_num.data(), state.xover_sequence.data(),
            spacer_sequences.data(), state.xdiffpos.data(),
            state.xposdiff.data(), valid_spacer.data());
        if (std::getenv("RDP_TRACE_XOVER_COMPARE") != nullptr &&
            ((state.sequences[0] == 6 && state.sequences[1] == 12 &&
              state.sequences[2] == 15) ||
             (state.sequences[0] == 15 && state.sequences[1] == 6 &&
              state.sequences[2] == 12))) {
            std::array<int, 3> pb3_ah{};
            std::vector<char> pb3_xover(state.xover_sequence.size(), 0);
            const int pb3_result = api.find_subsequence(
                pb3_ah.data(), fss_ub, xover_window,
                scan_state.compressed_sequence_ub, scan_state.sequence_length,
                scan_state.next_no, state.sequences[0], state.sequences[1],
                state.sequences[2], const_cast<unsigned char*>(
                    scan_state.compressed_sequence.data()),
                state.xover_sequence_ub, pb3_xover.data(), fss_rdp.data());
            std::cerr << "xover-compare " << state.sequences[0] << ':'
                      << state.sequences[1] << ':' << state.sequences[2]
                      << " plain="
                      << state.informative_length << ':'
                      << state.agreement_counts[0] << ','
                      << state.agreement_counts[1] << ','
                      << state.agreement_counts[2] << " pb3="
                      << pb3_result << ':' << pb3_ah[0] << ','
                      << pb3_ah[1] << ',' << pb3_ah[2] << '\n';
        }
    } else {
        state.informative_length = api.find_subsequence(
            state.agreement_counts.data(), fss_ub, xover_window,
            scan_state.compressed_sequence_ub, scan_state.sequence_length,
            scan_state.next_no, state.sequences[0], state.sequences[1],
            state.sequences[2],
            const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
            state.xover_sequence_ub, state.xover_sequence.data(),
            fss_rdp.data());
    }
    if (std::getenv("RDP_TRACE_XOVER_PATH") != nullptr) {
        std::cerr << "xover-path="
                  << (sequence_event_number > 0 && !use_compress
                          ? "plain" : "pb3")
                  << " se=" << sequence_event_number << " triplet="
                  << state.sequences[0] << ':' << state.sequences[1] << ':'
                  << state.sequences[2] << " result="
                  << state.informative_length << " ah="
                  << state.agreement_counts[0] << ','
                  << state.agreement_counts[1] << ','
                  << state.agreement_counts[2] << '\n';
    }

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
    const auto map_size =
        static_cast<std::size_t>(scan_state.sequence_length + 201);
    if (state.xdiffpos.size() != map_size) {
        state.xdiffpos.assign(map_size, 0);
    }
    if (state.xposdiff.size() != map_size) {
        state.xposdiff.assign(map_size, 0);
    }
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
    const std::array<int, 3>* explicit_sequences,
    const int sequence_event_number,
    const std::vector<unsigned char>* missing_data,
    const bool use_compress,
    int* shared_xdiffpos0) {
    RdpRawEventState events = initial_events == nullptr
        ? RdpRawEventState{} : *initial_events;
    if (initial_events == nullptr) {
        events.current_xover.assign(scan_state.next_no + 1, 0);
        events.xover_list.resize(scan_state.next_no + 1);
    }
    events.triplets_with_events.assign(
        static_cast<std::size_t>(scan_state.analysis_list_last + 1), 0);
    for (int triplet = 0; triplet <= scan_state.analysis_list_last;
         ++triplet) {
        if (static_cast<std::size_t>(triplet) >= redo.size() ||
            redo[triplet] != 1) {
            continue;
        }
        ++events.scanned_triplets;
        auto state = build_rdp_first_xover_state(
            scan_state, distance_state, tree_state, triplet, fss_ub, fss_rdp,
            xover_window, xover_window_x, api, explicit_sequences,
            sequence_event_number, use_compress);
        if (shared_xdiffpos0 != nullptr) {
            state.xdiffpos[0] = *shared_xdiffpos0;
        }
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
                        // The compressed initial XOver path obtains its
                        // breakpoint maps lazily through FindSubSeqPB4.  A
                        // FinalTrim XOver uses the plain FindSubSeqP path,
                        // which has already filled the maps for this exact
                        // triplet.  Replacing those per-triplet maps with a
                        // map retained from an earlier triplet changes event
                        // fragmentation and is not what the VB execution
                        // does.
                        if ((sequence_event_number == 0 || use_compress) &&
                            !position_maps_built) {
                            build_rdp_first_position_maps(
                                state, scan_state, fss_ub, xover_window,
                                fss_rdp, api);
                            position_maps_built = true;
                        }
                        const auto store_event = [&] (
                            const int beginning, const int ending,
                            const double probability,
                            const int beginning_warning,
                            const int ending_warning) {
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
                            event.probability = probability;
                            if (beginning_warning == 1 &&
                                ending_warning == 1) {
                                event.sbp_flag = 3;
                            } else if (beginning_warning == 1) {
                                event.sbp_flag = 1;
                            } else if (ending_warning == 1) {
                                event.sbp_flag = 2;
                            }
                            events.xover_list[storage[0]].push_back(event);
                            events.triplets_with_events[triplet] = 1;
                            events.current_xover[storage[0]] =
                                static_cast<std::int16_t>(
                                    events.xover_list[storage[0]].size());
                        };

                        const auto mapped = map_breakpoints(
                            state, scan_state, xover_settings,
                            state.event_begin, state.event_end);
                        int beginning_warning = 0;
                        int ending_warning = 0;
                        int split = 0;
                        int split_position = 0;
                        if (xover_settings.long_winded == 1 &&
                            sequence_event_number > 0 &&
                            missing_data != nullptr) {
                            int prior_beginning = mapped.beginning;
                            if (state.xdiffpos[mapped.beginning] > 0) {
                                prior_beginning = state.xdiffpos[
                                    mapped.original_beginning - 1];
                            }
                            int next_ending = mapped.ending;
                            if (state.xdiffpos[mapped.ending] <
                                scan_state.sequence_length) {
                                next_ending = state.xdiffpos[
                                    mapped.original_ending + 1];
                            }
                            if (next_ending < mapped.ending) {
                                next_ending = mapped.ending;
                            }
                            if (prior_beginning < 1) {
                                prior_beginning = scan_state.sequence_length;
                            }
                            if (next_ending > scan_state.sequence_length) {
                                next_ending = 1;
                            }
                            if (prior_beginning != mapped.beginning) {
                                int warning_split = 0;
                                check_split(
                                    2, scan_state.sequence_length,
                                    prior_beginning, mapped.beginning,
                                    state.sequences, warning_split,
                                    *missing_data);
                                if (warning_split == 1) {
                                    beginning_warning = 1;
                                }
                            }
                            if (mapped.ending != next_ending) {
                                int warning_split = 0;
                                check_split(
                                    2, scan_state.sequence_length,
                                    mapped.ending, next_ending,
                                    state.sequences, warning_split,
                                    *missing_data);
                                if (warning_split == 1) ending_warning = 1;
                            }
                            split_position = check_split(
                                10, scan_state.sequence_length,
                                mapped.beginning, mapped.ending,
                                state.sequences, split, *missing_data);
                        }

                        if (split == 1) {
                            for (int fragment = 0; fragment <= 1;
                                 ++fragment) {
                                int fragment_beginning = 0;
                                int fragment_ending = 0;
                                if (fragment == 0) {
                                    fragment_beginning = mapped.beginning;
                                    if (mapped.beginning == split_position) {
                                        ++split_position;
                                    }
                                    if (state.xposdiff[split_position] > 0) {
                                        fragment_ending = state.xdiffpos[
                                            state.xposdiff[split_position] - 1];
                                    } else {
                                        fragment_ending =
                                            scan_state.sequence_length;
                                    }
                                } else {
                                    int missing_position = find_missing(
                                        scan_state.sequence_length,
                                        state.sequences, split_position,
                                        mapped.ending, *missing_data);
                                    if (mapped.ending == missing_position) {
                                        fragment_beginning =
                                            mapped.ending - 1;
                                        fragment_ending = mapped.ending;
                                    } else {
                                        const int position_index =
                                            state.xposdiff[missing_position];
                                        if (position_index + 1 <
                                                state.homology_length &&
                                            (state.xdiffpos[
                                                position_index + 1] >
                                                    missing_position ||
                                             state.xdiffpos[
                                                position_index + 1] <
                                                    mapped.ending)) {
                                            fragment_beginning =
                                                state.xdiffpos[
                                                    position_index + 1];
                                        } else {
                                            const int stride =
                                                scan_state.sequence_length + 1;
                                            while (true) {
                                                const short first =
                                                    scan_state.sequence_data[
                                                        missing_position +
                                                        state.sequences[0] *
                                                            stride];
                                                const short second =
                                                    scan_state.sequence_data[
                                                        missing_position +
                                                        state.sequences[1] *
                                                            stride];
                                                const short third =
                                                    scan_state.sequence_data[
                                                        missing_position +
                                                        state.sequences[2] *
                                                            stride];
                                                if (first != 46 &&
                                                    second != 46 &&
                                                    third != 46 &&
                                                    (first != second ||
                                                     first != third)) {
                                                    break;
                                                }
                                                ++missing_position;
                                                if (missing_position ==
                                                    mapped.ending) {
                                                    ++missing_position;
                                                    break;
                                                }
                                                if (missing_position >
                                                    scan_state.sequence_length) {
                                                    missing_position = 1;
                                                }
                                            }
                                            if (missing_position >
                                                scan_state.sequence_length) {
                                                missing_position = 1;
                                            }
                                            fragment_beginning =
                                                missing_position;
                                        }
                                        fragment_ending = mapped.ending;
                                    }
                                }

                                int event_length = 0;
                                int number_in_common = 0;
                                const int number_different = split_event(
                                    xover_window, scan_state.sequence_length,
                                    state.homology_length,
                                    state.sequence_daughter,
                                    state.sequence_minor,
                                    fragment_beginning, fragment_ending,
                                    event_length, number_in_common,
                                    state.xposdiff, state.xover_sequence);
                                position =
                                    state.xposdiff[fragment_ending] + 1;
                                double fragment_probability = 1.0;
                                if (event_length > 2) {
                                    int probability_length = event_length;
                                    int probability_different =
                                        number_different;
                                    if (event_length >= 170) {
                                        probability_different = vb_clng(
                                            static_cast<double>(
                                                number_different) * 169.0 /
                                            static_cast<double>(event_length));
                                        probability_length = 169;
                                    }
                                    const int probability_same =
                                        probability_length -
                                        probability_different;
                                    if (event_length <=
                                        probability_settings.fact_three_ub) {
                                        fragment_probability =
                                            api.probability_p2(
                                                fact_three.data(),
                                                probability_settings.
                                                    fact_three_ub,
                                                probability_length,
                                                probability_same,
                                                state.individual_probability,
                                                state.homology_length);
                                    } else {
                                        fragment_probability =
                                            api.probability_p(
                                                fact.data(),
                                                probability_length,
                                                probability_same,
                                                state.individual_probability,
                                                state.homology_length);
                                    }
                                    if (fragment_probability > 1.0) {
                                        fragment_probability = 1.0;
                                    }
                                }
                                const double minimum_probability =
                                    std::pow(10.0, -300.0);
                                if (fragment_probability <
                                    minimum_probability) {
                                    fragment_probability =
                                        minimum_probability;
                                }
                                if (probability_settings.mc_flag == 0) {
                                    fragment_probability *=
                                        probability_settings.mc_correction;
                                }
                                if (fragment_probability <
                                        probability_settings.
                                            lowest_probability &&
                                    fragment_probability > 0.0) {
                                    int centred_beginning =
                                        fragment_beginning;
                                    int centred_ending = fragment_ending;
                                    int fragment_beginning_warning =
                                        beginning_warning;
                                    int fragment_ending_warning =
                                        ending_warning;
                                    centre_mapped_breakpoints(
                                        state, scan_state, xover_window,
                                        xover_settings, centred_beginning,
                                        centred_ending,
                                        fragment_beginning_warning,
                                        fragment_ending_warning,
                                        sequence_event_number,
                                        missing_data, 0, 0, use_compress);
                                    if (sequence_event_number > 0 &&
                                        missing_data != nullptr) {
                                        if (fragment_ending_warning == 0) {
                                            fragment_ending_warning =
                                                check_ends(
                                                    state, scan_state,
                                                    xover_window,
                                        xover_settings, 1,
                                        centred_beginning,
                                        centred_ending,
                                        *missing_data);
                                        }
                                        if (fragment_beginning_warning == 0) {
                                            fragment_beginning_warning =
                                                check_ends(
                                                    state, scan_state,
                                                    xover_window,
                                                    xover_settings, 0,
                                                    centred_beginning,
                                                    centred_ending,
                                                    *missing_data);
                                        }
                                    }
                                    store_event(
                                        centred_beginning, centred_ending,
                                        fragment_probability,
                                        fragment_beginning_warning,
                                        fragment_ending_warning);
                                }
                            }
                        } else {
                            int beginning = mapped.beginning;
                            int ending = mapped.ending;
                            centre_mapped_breakpoints(
                                state, scan_state, xover_window,
                                xover_settings, beginning, ending,
                                        beginning_warning, ending_warning,
                                sequence_event_number, missing_data,
                                mapped.original_beginning,
                                mapped.original_ending, use_compress);
                            if (sequence_event_number > 0 &&
                                missing_data != nullptr) {
                                if (ending_warning == 0) {
                                    ending_warning = check_ends(
                                        state, scan_state, xover_window,
                                        xover_settings, 1, beginning, ending,
                                        *missing_data,
                                        use_compress ? 0 :
                                            mapped.original_beginning,
                                        use_compress ? 0 :
                                            mapped.original_ending);
                                }
                                if (beginning_warning == 0) {
                                    beginning_warning = check_ends(
                                        state, scan_state, xover_window,
                                        xover_settings, 0, beginning, ending,
                                        *missing_data,
                                        use_compress ? 0 :
                                            mapped.original_beginning,
                                        use_compress ? 0 :
                                            mapped.original_ending);
                                }
                            }
                            if (std::getenv("RDP_TRACE_SBP") != nullptr &&
                                sequence_event_number > 0 &&
                                ((beginning == 1 && ending == 1231) ||
                                 (beginning == 7 && ending == 436))) {
                                std::cerr << "sbp-trace seq="
                                          << state.sequences[0] << ':'
                                          << state.sequences[1] << ':'
                                          << state.sequences[2] << " bp="
                                          << beginning << '-' << ending
                                          << " mapped=" << mapped.beginning
                                          << '-' << mapped.ending
                                          << " orig="
                                          << mapped.original_beginning << '-'
                                          << mapped.original_ending
                                          << " circ=" << xover_settings.circular
                                          << " xdp=" << state.xdiffpos[1]
                                          << ',' << state.xdiffpos[2]
                                          << ',' << state.xdiffpos[3]
                                          << ',' << state.xdiffpos[4]
                                          << " miss=";
                                if (missing_data != nullptr) {
                                    const int stride =
                                        scan_state.sequence_length + 1;
                                    int missing_count = 0;
                                    for (int pos = 1; pos <=
                                             scan_state.sequence_length; ++pos) {
                                        for (const int seq : state.sequences) {
                                            const auto index =
                                                static_cast<std::size_t>(pos) +
                                                static_cast<std::size_t>(seq) *
                                                    stride;
                                            if (index < missing_data->size() &&
                                                (*missing_data)[index] == 1) {
                                                ++missing_count;
                                                break;
                                            }
                                        }
                                    }
                                    std::cerr << missing_count;
                                } else {
                                    std::cerr << "null";
                                }
                                if (missing_data != nullptr) {
                                    const int stride =
                                        scan_state.sequence_length + 1;
                                    const auto endpoint_missing =
                                        [&](const int position) {
                                            for (const int seq : state.sequences) {
                                                const auto index =
                                                    static_cast<std::size_t>(position) +
                                                    static_cast<std::size_t>(seq) * stride;
                                                if (index < missing_data->size() &&
                                                    (*missing_data)[index] == 1) {
                                                    return true;
                                                }
                                            }
                                            return false;
                                        };
                                    std::cerr << " endmiss="
                                              << endpoint_missing(beginning)
                                              << ':' << endpoint_missing(ending);
                                }
                                std::cerr << " warn=" << beginning_warning
                                          << ':' << ending_warning << '\n';
                            }
                            store_event(
                                beginning, ending,
                                state.adjusted_event_probability,
                                beginning_warning, ending_warning);
                        }
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

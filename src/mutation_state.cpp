#include "mutation_state.hpp"

#include "MathFuncsDll.h"

#include <algorithm>
#include <stdexcept>

namespace {

XOVERDEFINE to_dna_event(const RdpRawEvent& source) {
    XOVERDEFINE event{};
    event.OutsideFlag = source.outside_flag;
    event.MissIdentifyFlag = source.misidentify_flag;
    event.ProgramFlag = source.program_flag;
    event.SBPFlag = source.sbp_flag;
    event.Accept = source.accept;
    event.MajorP = source.major_parent;
    event.MinorP = source.minor_parent;
    event.Daughter = source.daughter;
    event.Beginning = source.beginning;
    event.Ending = source.ending;
    event.LHolder = source.length_holder;
    event.Eventnumber = source.event_number;
    event.PermPVal = source.permutation_pvalue;
    event.BeginP = source.begin_parent;
    event.EndP = source.end_parent;
    event.Probability = source.probability;
    event.DHolder = source.distance_holder;
    return event;
}

RdpRawEvent from_dna_event(const XOVERDEFINE& source) {
    RdpRawEvent event;
    event.outside_flag = source.OutsideFlag;
    event.misidentify_flag = source.MissIdentifyFlag;
    event.program_flag = source.ProgramFlag;
    event.sbp_flag = source.SBPFlag;
    event.accept = source.Accept;
    event.major_parent = source.MajorP;
    event.minor_parent = source.MinorP;
    event.daughter = source.Daughter;
    event.beginning = source.Beginning;
    event.ending = source.Ending;
    event.length_holder = source.LHolder;
    event.event_number = source.Eventnumber;
    event.permutation_pvalue = source.PermPVal;
    event.begin_parent = source.BeginP;
    event.end_parent = source.EndP;
    event.probability = source.Probability;
    event.distance_holder = source.DHolder;
    return event;
}

std::size_t sequence_cell(const int sequence_length, const int position,
                          const int sequence) {
    return static_cast<std::size_t>(position) +
        static_cast<std::size_t>(sequence) * (sequence_length + 1);
}

int candidate(const std::vector<int>& list, const int role, const int slot) {
    return list[role + slot * 3];
}

void validate_candidate_state(const int next_no, const int winning_role,
                              const std::array<int, 3>& last,
                              const std::vector<int>& list) {
    if (winning_role < 0 || winning_role > 2 ||
        list.size() <= static_cast<std::size_t>(
            winning_role + std::max(last[winning_role], 0) * 3) ||
        last[winning_role] < 0 || last[winning_role] > next_no) {
        throw std::runtime_error("RDP mutation candidate state differs");
    }
}

}  // namespace

RdpAdjustedEvents adjust_rdp_events_exact(
    const int next_no, const int winning_role,
    const double lowest_probability,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const RdpRawEventState& source_events,
    const std::vector<unsigned char>& source_done,
    const int source_done_row_upper_bound,
    const int source_done_slot_upper_bound) {
    const int row_count = next_no + 1;
    int source_slot_upper_bound = 1;
    for (const auto& row : source_events.xover_list) {
        source_slot_upper_bound = std::max<int>(
            source_slot_upper_bound, row.size());
    }
    const int output_slot_upper_bound = source_slot_upper_bound + 4;
    std::vector<XOVERDEFINE> source_matrix(
        static_cast<std::size_t>(row_count) *
        (source_slot_upper_bound + 1));
    auto source_current = source_events.current_xover;
    source_current.resize(row_count, 0);
    for (int row = 0; row <= next_no; ++row) {
        if (row >= static_cast<int>(source_events.xover_list.size())) continue;
        for (std::size_t slot = 0;
             slot < source_events.xover_list[row].size(); ++slot) {
            source_matrix[row + (slot + 1) * row_count] =
                to_dna_event(source_events.xover_list[row][slot]);
        }
    }
    std::vector<unsigned char> aligned_done(
        static_cast<std::size_t>(row_count) *
        (source_slot_upper_bound + 1), 0);
    const int copy_slots = std::min(
        source_slot_upper_bound, source_done_slot_upper_bound);
    const int copy_rows = std::min(
        next_no, source_done_row_upper_bound);
    for (int slot = 0; slot <= copy_slots; ++slot) {
        for (int row = 0; row <= copy_rows; ++row) {
            aligned_done[row + slot * row_count] = source_done[
                row + slot * (source_done_row_upper_bound + 1)];
        }
    }

    RdpAdjustedEvents output;
    output.done_row_upper_bound = next_no;
    output.done_slot_upper_bound = output_slot_upper_bound;
    output.done_sequence.assign(
        static_cast<std::size_t>(row_count) *
        (output_slot_upper_bound + 1), 0);
    output.pairs_to_rescan.assign(
        static_cast<std::size_t>(row_count) * row_count, 0);
    std::vector<std::int16_t> output_current(row_count, 0);
    std::vector<XOVERDEFINE> output_matrix(
        static_cast<std::size_t>(row_count) *
        (output_slot_upper_bound + 1));
    std::array<int, 101> event_counts{};
    auto mutable_candidate_last = candidate_last;
    auto mutable_candidate_list = candidate_list;
    auto mutable_trace_sub = trace_sub;
    MathFuncs::MyMathFuncs::AddjustCXO(
        next_no, winning_role, lowest_probability,
        next_no, source_slot_upper_bound, aligned_done.data(),
        next_no, output_slot_upper_bound, output.done_sequence.data(),
        event_counts.data(), mutable_candidate_last.data(),
        mutable_candidate_list.data(), output.pairs_to_rescan.data(),
        static_cast<int>(mutable_trace_sub.size()) - 1,
        mutable_trace_sub.data(), output_current.data(), next_no,
        output_slot_upper_bound, output_matrix.data(), source_current.data(),
        next_no, source_slot_upper_bound, source_matrix.data());
    output.events.current_xover = output_current;
    output.events.xover_list.resize(row_count);
    for (int row = 0; row <= next_no; ++row) {
        for (int slot = 1; slot <= output_current[row]; ++slot) {
            output.events.xover_list[row].push_back(from_dna_event(
                output_matrix[row + slot * row_count]));
        }
    }
    return output;
}

RdpErasedTracts erase_rdp_recombinant_tracts(
    const int sequence_length, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const int beginning, const int ending, std::vector<int> breakpoints,
    const std::vector<short>& sequence_data,
    const std::vector<unsigned char>& missing_data) {
    const int next_no = static_cast<int>(
        sequence_data.size() / (sequence_length + 1)) - 1;
    validate_candidate_state(
        next_no, winning_role, candidate_last, candidate_list);
    if (breakpoints.size() <
        static_cast<std::size_t>(2 * (candidate_last[winning_role] + 1)) ||
        sequence_data.size() != missing_data.size()) {
        throw std::runtime_error("ModSeqNumY input dimensions differ");
    }

    RdpErasedTracts output;
    output.sequence_data = sequence_data;
    output.missing_data = missing_data;
    output.breakpoints = std::move(breakpoints);
    output.saved_tracts.assign(
        static_cast<std::size_t>(sequence_length + 1) *
            (candidate_last[winning_role] + 1),
        0);

    // Literal DNA.dll ModSeqNumY. Its comments describe an endpoint penalty,
    // but the active loops erase and save both endpoints inclusively.
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int sequence = candidate(candidate_list, winning_role, slot);
        if (output.breakpoints[slot * 2] == 0 &&
            output.breakpoints[slot * 2 + 1] == 0) {
            output.breakpoints[slot * 2] = beginning;
            output.breakpoints[slot * 2 + 1] = ending;
        }
        const int first = output.breakpoints[slot * 2];
        const int last = output.breakpoints[slot * 2 + 1];
        const auto erase_position = [&](const int position) {
            const auto source = sequence_cell(
                sequence_length, position, sequence);
            const auto saved = sequence_cell(sequence_length, position, slot);
            output.saved_tracts[saved] = output.sequence_data[source];
            output.sequence_data[source] = 46;
            output.missing_data[source] = 1;
        };
        if (first < last) {
            for (int position = first; position <= last; ++position) {
                erase_position(position);
            }
        } else {
            for (int position = first; position <= sequence_length; ++position) {
                erase_position(position);
            }
            for (int position = 1; position <= last; ++position) {
                erase_position(position);
            }
        }
    }
    return output;
}

void rebuild_rdp_recombinant_tracts(
    const int sequence_length, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& breakpoints,
    const std::vector<short>& saved_tracts,
    std::vector<short>& sequence_data) {
    const int next_no = static_cast<int>(
        sequence_data.size() / (sequence_length + 1)) - 1;
    validate_candidate_state(
        next_no, winning_role, candidate_last, candidate_list);
    // Literal DNA.dll RebuildSeqNum: unlike ModSeqNumY it restores only the
    // strict interior, so the two endpoints remain erased.
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int sequence = candidate(candidate_list, winning_role, slot);
        const int first = breakpoints[slot * 2];
        const int last = breakpoints[slot * 2 + 1];
        const auto restore = [&](const int position) {
            sequence_data[sequence_cell(sequence_length, position, sequence)] =
                saved_tracts[sequence_cell(sequence_length, position, slot)];
        };
        if (first < last) {
            for (int position = first + 1; position < last; ++position) {
                restore(position);
            }
        } else {
            for (int position = first + 1; position <= sequence_length;
                 ++position) {
                restore(position);
            }
            for (int position = 1; position < last; ++position) {
                restore(position);
            }
        }
    }
}

void make_rdp_fragment_rows(
    const int sequence_length, const int next_no, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const int beginning, const int ending, std::vector<int>& breakpoints,
    std::vector<short>& sequence_data,
    std::vector<unsigned char>& missing_data) {
    validate_candidate_state(
        next_no, winning_role, candidate_last, candidate_list);
    const auto required = static_cast<std::size_t>(next_no + 1) *
        (sequence_length + 1);
    if (sequence_data.size() < required) sequence_data.resize(required, 0);
    if (missing_data.size() < required) missing_data.resize(required, 0);

    // Literal DNA.dll ModSN. The caller has already increased NextNo by the
    // winning list size, so the new rows occupy its final contiguous block.
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int new_sequence = next_no - candidate_last[winning_role] + slot;
        const int source_sequence =
            candidate(candidate_list, winning_role, slot);
        if (breakpoints[slot * 2] == 0 &&
            breakpoints[slot * 2 + 1] == 0) {
            breakpoints[slot * 2] = beginning;
            breakpoints[slot * 2 + 1] = ending;
        }
        const int first = breakpoints[slot * 2];
        const int last = breakpoints[slot * 2 + 1];
        const auto copy = [&](const int position) {
            sequence_data[sequence_cell(
                sequence_length, position, new_sequence)] =
                sequence_data[sequence_cell(
                    sequence_length, position, source_sequence)];
        };
        const auto mask = [&](const int position) {
            const auto cell = sequence_cell(
                sequence_length, position, new_sequence);
            sequence_data[cell] = 46;
            missing_data[cell] = 1;
        };
        if (first < last) {
            for (int position = 1; position < first; ++position) mask(position);
            for (int position = first; position <= last; ++position) copy(position);
            for (int position = last + 1; position <= sequence_length;
                 ++position) mask(position);
        } else {
            for (int position = first; position <= sequence_length; ++position) {
                copy(position);
            }
            for (int position = last + 1; position < first; ++position) {
                mask(position);
            }
            for (int position = 1; position <= last; ++position) copy(position);
        }
    }
}

void erase_rdp_original_tracts(
    const int sequence_length, const int next_no, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& selected_candidates,
    const int beginning, const int ending, std::vector<int>& breakpoints,
    std::vector<short>& sequence_data,
    std::vector<unsigned char>& missing_data) {
    if (winning_role < 0 || winning_role > 2 ||
        candidate_last[winning_role] < 0 ||
        selected_candidates.size() <=
            static_cast<std::size_t>(candidate_last[winning_role])) {
        throw std::runtime_error("ModSeqNumZ input dimensions differ");
    }
    // Literal DNA.dll ModSeqNumZ. The caller passes the compact winning-role
    // list here, rather than the three-row RList layout used by ModSN.
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int sequence = selected_candidates[slot];
        if (breakpoints[slot * 2] == 0 &&
            breakpoints[slot * 2 + 1] == 0) {
            breakpoints[slot * 2] = beginning;
            breakpoints[slot * 2 + 1] = ending;
        }
        const int first = breakpoints[slot * 2];
        const int last = breakpoints[slot * 2 + 1];
        const auto erase = [&](const int position) {
            const auto cell = sequence_cell(sequence_length, position, sequence);
            sequence_data[cell] = 46;
            missing_data[cell] = 1;
        };
        if (first < last) {
            for (int position = first; position <= last; ++position) erase(position);
        } else {
            for (int position = first; position <= sequence_length; ++position) {
                erase(position);
            }
            for (int position = 1; position <= last; ++position) erase(position);
        }
    }
}

std::vector<int> calculate_rdp_actual_sequence_sizes(
    const int sequence_length, const int next_no, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<short>& sequence_data) {
    validate_candidate_state(
        next_no, winning_role, candidate_last, candidate_list);
    std::vector<int> sizes(next_no + 1, 0);
    // Literal DNA.dll MakeActualSeqSize for the affected original/new pairs.
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int new_sequence = next_no - candidate_last[winning_role] + slot;
        const int source_sequence =
            candidate(candidate_list, winning_role, slot);
        for (int position = 1; position <= sequence_length; ++position) {
            if (sequence_data[sequence_cell(
                    sequence_length, position, source_sequence)] != 46) {
                ++sizes[source_sequence];
            }
            if (sequence_data[sequence_cell(
                    sequence_length, position, new_sequence)] != 46) {
                ++sizes[new_sequence];
            }
        }
    }
    return sizes;
}

RdpRedistributedEvents redistribute_rdp_events(
    const int next_no, const int winning_role,
    const double lowest_probability,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const RdpRawEventState& source_events) {
    validate_candidate_state(
        next_no, winning_role, candidate_last, candidate_list);
    if (trace_sub.size() < static_cast<std::size_t>(next_no + 1) ||
        source_events.xover_list.size() < static_cast<std::size_t>(next_no + 1)) {
        throw std::runtime_error("AddjustCXO input dimensions differ");
    }
    RdpRedistributedEvents output;
    output.events.current_xover.assign(next_no + 1, 0);
    output.events.xover_list.resize(next_no + 1);
    output.pairs_to_rescan.assign(
        static_cast<std::size_t>(next_no + 1) * (next_no + 1), 0);
    const auto mark_pairs = [&](const int daughter, const int major,
                                const int minor) {
        for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
            const int member = candidate(candidate_list, winning_role, slot);
            if (member == daughter || member == major || member == minor) {
                const int count = next_no + 1;
                output.pairs_to_rescan[minor + major * count] = 1;
                output.pairs_to_rescan[major + minor * count] = 1;
                output.pairs_to_rescan[daughter + major * count] = 1;
                output.pairs_to_rescan[major + daughter * count] = 1;
                output.pairs_to_rescan[minor + daughter * count] = 1;
                output.pairs_to_rescan[daughter + minor * count] = 1;
                return false;
            }
        }
        return true;
    };

    // Literal DNA5.dll AddjustCXO without its fixed backing-array cap. Only
    // active source records are represented by the vector form.
    for (int source_row = 0; source_row <= next_no; ++source_row) {
        for (const auto& source_event : source_events.xover_list[source_row]) {
            int daughter = source_event.daughter;
            int minor = source_event.minor_parent;
            int major = source_event.major_parent;
            if (daughter > next_no) daughter = trace_sub[daughter];
            if (minor >= static_cast<int>(trace_sub.size())) minor = 0;
            if (major >= static_cast<int>(trace_sub.size())) major = 0;
            if (minor > next_no) minor = trace_sub[minor];
            if (major > next_no) major = trace_sub[major];
            const bool uninvolved = mark_pairs(daughter, major, minor);
            if (source_event.probability <= lowest_probability && uninvolved) {
                int destination = daughter;
                if (output.events.xover_list[daughter].size() <=
                        output.events.xover_list[minor].size() &&
                    output.events.xover_list[daughter].size() <=
                        output.events.xover_list[major].size()) {
                    destination = daughter;
                } else if (output.events.xover_list[minor].size() <=
                               output.events.xover_list[daughter].size() &&
                           output.events.xover_list[minor].size() <=
                               output.events.xover_list[major].size()) {
                    destination = minor;
                } else {
                    destination = major;
                }
                auto event = source_event;
                if (destination == minor) {
                    event.daughter = static_cast<std::int16_t>(minor);
                    event.minor_parent = static_cast<std::int16_t>(daughter);
                } else if (destination == major) {
                    event.daughter = static_cast<std::int16_t>(major);
                    event.major_parent = static_cast<std::int16_t>(daughter);
                }
                output.events.xover_list[destination].push_back(event);
                ++output.event_counts[100];
                if (event.program_flag < 100) ++output.event_counts[event.program_flag];
            }
        }
    }
    for (int row = 0; row <= next_no; ++row) {
        output.events.current_xover[row] = static_cast<std::int16_t>(
            output.events.xover_list[row].size());
    }
    return output;
}

RdpDroppedFragmentEvents drop_rdp_unused_fragment_events(
    const int original_next_no, const int expanded_next_no,
    const int minimum_sequence_size, const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const RdpRawEventState& events_before_outer_scan,
    const RdpRawEventState& outer_scan_events) {
    if (expanded_next_no < original_next_no ||
        trace_sub.size() < static_cast<std::size_t>(expanded_next_no + 1) ||
        actual_sequence_sizes.size() <
            static_cast<std::size_t>(expanded_next_no + 1)) {
        throw std::runtime_error("DropSeqs input dimensions differ");
    }
    RdpDroppedFragmentEvents output;
    output.next_no = expanded_next_no;
    output.events = events_before_outer_scan;
    output.events.xover_list.resize(expanded_next_no + 1);
    output.events.current_xover.resize(expanded_next_no + 1, 0);
    output.reference_counts.assign(expanded_next_no + 1, 0);
    output.trace_sub.assign(
        trace_sub.begin(), trace_sub.begin() + expanded_next_no + 1);
    output.actual_sequence_sizes.assign(
        actual_sequence_sizes.begin(),
        actual_sequence_sizes.begin() + expanded_next_no + 1);

    const auto map_to_original = [&](int sequence) {
        if (sequence > original_next_no) {
            sequence = sequence <= expanded_next_no
                ? output.trace_sub[sequence] : -1;
        }
        return sequence;
    };
    // Literal section 15 bookkeeping: retained and inner-scan records count
    // daughter and minor references, but accidentally omit the major parent.
    for (const auto& row : output.events.xover_list) {
        for (const auto& event : row) {
            const int daughter = map_to_original(event.daughter);
            const int minor = map_to_original(event.minor_parent);
            if (daughter >= 0 && daughter <= original_next_no) {
                ++output.reference_counts[daughter];
            }
            if (minor >= 0 && minor <= original_next_no) {
                ++output.reference_counts[minor];
            }
        }
    }
    // CopyXOListsX appends outer-scan records and counts all three roles.
    for (std::size_t row = 0; row < outer_scan_events.xover_list.size(); ++row) {
        for (const auto& event : outer_scan_events.xover_list[row]) {
            output.events.xover_list[row].push_back(event);
            for (const int sequence : {
                     static_cast<int>(event.daughter),
                     static_cast<int>(event.major_parent),
                     static_cast<int>(event.minor_parent)}) {
                if (sequence >= 0 && sequence <= expanded_next_no) {
                    ++output.reference_counts[sequence];
                }
            }
        }
        output.events.current_xover[row] = static_cast<std::int16_t>(
            output.events.xover_list[row].size());
    }

    const auto erase_references = [&](const int removed, const int last,
                                      const int last_row_to_check) {
        for (int row = 0; row <= last_row_to_check; ++row) {
            if (row == removed && removed < last) continue;
            auto& records = output.events.xover_list[row];
            std::size_t slot = 0;
            while (slot < records.size()) {
                const auto& event = records[slot];
                if (event.major_parent == removed ||
                    event.minor_parent == removed ||
                    event.daughter == removed) {
                    records[slot] = records.back();
                    records.pop_back();
                } else {
                    ++slot;
                }
            }
            output.events.current_xover[row] =
                static_cast<std::int16_t>(records.size());
        }
    };

    int sequence = original_next_no + 1;
    while (sequence <= output.next_no) {
        if (output.actual_sequence_sizes[sequence] >= minimum_sequence_size &&
            output.reference_counts[sequence] != 0) {
            ++sequence;
            continue;
        }
        const int last = output.next_no;
        if (sequence < last) {
            const int old_reference_count = output.reference_counts[sequence];
            output.actual_sequence_sizes[sequence] =
                output.actual_sequence_sizes[last];
            output.trace_sub[sequence] = output.trace_sub[last];
            output.reference_counts[sequence] = output.reference_counts[last];
            output.events.xover_list[sequence] = output.events.xover_list[last];
            for (auto& event : output.events.xover_list[sequence]) {
                event.daughter = static_cast<std::int16_t>(sequence);
            }
            output.events.current_xover[sequence] = static_cast<std::int16_t>(
                output.events.xover_list[sequence].size());
            if (old_reference_count > 0) {
                erase_references(sequence, last, last);
            }
            for (int row = 0; row <= last; ++row) {
                if (row == sequence) continue;
                for (auto& event : output.events.xover_list[row]) {
                    if (event.major_parent == last) {
                        event.major_parent = static_cast<std::int16_t>(sequence);
                    } else if (event.minor_parent == last) {
                        event.minor_parent = static_cast<std::int16_t>(sequence);
                    }
                }
            }
        } else {
            erase_references(last, last, original_next_no);
        }
        --output.next_no;
    }
    output.events.xover_list.resize(output.next_no + 1);
    output.events.current_xover.resize(output.next_no + 1);
    output.reference_counts.resize(output.next_no + 1);
    output.trace_sub.resize(output.next_no + 1);
    output.actual_sequence_sizes.resize(output.next_no + 1);
    return output;
}

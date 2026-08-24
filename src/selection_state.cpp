#include "selection_state.hpp"

#include "MathFuncsDll.h"

#include <algorithm>
#include <cstddef>

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

}  // namespace

RdpSelectionState select_rdp_best_event(
    const RdpRawEventState& events, const int next_no,
    const double probability_cutoff,
    const std::vector<unsigned char>& existing_done,
    const int existing_done_row_upper_bound) {
    RdpSelectionState output;
    output.next_no = next_no;
    output.probability = probability_cutoff;

    int maximum_slots = 0;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (sequence < static_cast<int>(events.xover_list.size())) {
            maximum_slots = std::max(
                maximum_slots,
                static_cast<int>(events.xover_list[sequence].size()));
        }
    }
    output.slot_upper_bound = maximum_slots;
    const int row_count = next_no + 1;
    const int slot_count = maximum_slots + 1;
    std::vector<XOVERDEFINE> xovers(
        static_cast<std::size_t>(row_count) * slot_count);
    std::vector<short> current_xover(row_count, 0);
    output.done_sequence.assign(
        static_cast<std::size_t>(row_count) * slot_count, 0);
    output.test_probabilities.assign(
        static_cast<std::size_t>(row_count) * slot_count, 0.0);

    if (!existing_done.empty() && existing_done_row_upper_bound >= 0) {
        const int existing_rows = existing_done_row_upper_bound + 1;
        const int existing_slots = static_cast<int>(
            existing_done.size() / static_cast<std::size_t>(existing_rows));
        for (int slot = 0; slot < std::min(slot_count, existing_slots); ++slot) {
            for (int row = 0; row < std::min(row_count, existing_rows); ++row) {
                output.done_sequence[row + slot * row_count] =
                    existing_done[row + slot * existing_rows];
            }
        }
    }

    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (sequence >= static_cast<int>(events.xover_list.size())) continue;
        const auto& source_row = events.xover_list[sequence];
        current_xover[sequence] = static_cast<short>(source_row.size());
        for (std::size_t slot = 0; slot < source_row.size(); ++slot) {
            xovers[sequence + (slot + 1) * static_cast<std::size_t>(row_count)] =
                to_dna_event(source_row[slot]);
        }
    }

    output.make_test_result = MathFuncs::MyMathFuncs::MakeTestPVs(
        next_no, output.done_sequence.data(), next_no, next_no,
        maximum_slots, current_xover.data(), xovers.data(),
        output.test_probabilities.data());

    const auto run_selection = [&](const int done_target) {
        output.done_target = done_target;
        output.trace = {0, 0};
        output.probability = probability_cutoff;
        output.total_candidates =
            MathFuncs::MyMathFuncs::FindBestRecSignalP2(
                static_cast<char>(done_target), next_no, maximum_slots,
                next_no, &output.probability,
                reinterpret_cast<char*>(output.done_sequence.data()),
                output.trace.data(), current_xover.data(),
                output.test_probabilities.data());
        output.found = output.probability != probability_cutoff;
    };
    run_selection(0);
    if (!output.found) run_selection(1);
    return output;
}

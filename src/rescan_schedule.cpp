#include "rescan_schedule.hpp"

#include "MathFuncsDll.h"

#include <stdexcept>
#include <cmath>

namespace {

int list_member(const std::vector<int>& list, int role, int slot) {
    return list[role + slot * 3];
}

void validate_pairs(int upper_bound,
                    const std::vector<unsigned char>& pairs) {
    const auto stride = static_cast<std::size_t>(upper_bound + 1);
    if (pairs.size() < stride * stride) {
        throw std::runtime_error("RDP rescan pair matrix differs");
    }
}

int pair_upper_bound(const std::vector<unsigned char>& pairs) {
    const auto stride = static_cast<int>(std::sqrt(
        static_cast<double>(pairs.size())));
    if (stride <= 0 ||
        static_cast<std::size_t>(stride) * stride != pairs.size()) {
        throw std::runtime_error("pair-rescan matrix is not square");
    }
    return stride - 1;
}

}  // namespace

std::vector<short> flatten_rdp_triplets(
    const std::vector<std::array<int, 3>>& triplets) {
    std::vector<short> output;
    output.reserve(triplets.size() * 3);
    for (const auto& triplet : triplets) {
        for (const int sequence : triplet) {
            output.push_back(static_cast<short>(sequence));
        }
    }
    return output;
}

std::vector<unsigned char> screen_rdp_rescan_triplets(
    const std::vector<std::array<int, 3>>& triplets,
    const RdpScanState& scan_state,
    const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state,
    const RdpRescanScreenSettings& settings,
    std::vector<unsigned char>& fss_rdp,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three,
    std::vector<double>& fact) {
    if (triplets.empty()) return {};
    auto analysis_list = flatten_rdp_triplets(triplets);
    std::vector<unsigned char> redo(triplets.size(), 0);
    const double uncorrected_threshold = settings.correction_flag == 0
        ? settings.probability_cutoff / settings.correction_tests
        : settings.probability_cutoff;
    MathFuncs::MyMathFuncs::AlistRDP3(
        analysis_list.data(), static_cast<int>(triplets.size()) - 1, 0,
        static_cast<int>(triplets.size()) - 1, scan_state.next_no,
        uncorrected_threshold, redo.data(), settings.circular,
        settings.correction_tests, settings.correction_flag,
        settings.probability_cutoff, settings.target,
        scan_state.sequence_length, settings.short_output,
        scan_state.next_no,
        const_cast<float*>(distance_state.distance.data()),
        scan_state.next_no,
        const_cast<float*>(tree_state.tree_distance.data()),
        settings.fss_upper_bound, scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        const_cast<short*>(scan_state.sequence_data.data()),
        settings.half_window, settings.full_window, fss_rdp.data(),
        settings.probability_file_flag,
        settings.probability_one_upper_bound,
        settings.probability_two_upper_bound, probability_estimate.data(),
        settings.factorial_three_upper_bound, fact_three.data(), fact.data());
    return redo;
}

void propagate_rdp_group_pairs(
    const int next_no, const int winning_role,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    std::vector<unsigned char>& pairs_to_rescan) {
    validate_pairs(next_no, pairs_to_rescan);
    const int stride = next_no + 1;
    // Literal Module3 section 10 propagation after AddjustCXO.
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
            const int member = list_member(candidate_list, winning_role, slot);
            if (pairs_to_rescan[member + sequence * stride] == 1) {
                for (int other = 0;
                     other <= candidate_last[winning_role]; ++other) {
                    const int group_member =
                        list_member(candidate_list, winning_role, other);
                    pairs_to_rescan[group_member + sequence * stride] = 1;
                    pairs_to_rescan[sequence + group_member * stride] = 1;
                }
                break;
            }
        }
    }
}

std::vector<std::array<int, 3>> make_rdp_inner_scan_triplets(
    const RdpScanState& original_scan_state,
    const std::vector<unsigned char>& initially_screened_triplets,
    const int winning_role, const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const int permanent_next_no, const int minimum_sequence_size,
    const std::vector<unsigned char>& pairs_to_rescan) {
    // AnalysisList remains the permanent-sequence list while DoPairs grows
    // and shrinks with fragment rows. Its upper bound therefore comes from
    // the pair matrix, not the permanent analysis-list state.
    const int next_no = pair_upper_bound(pairs_to_rescan);
    validate_pairs(next_no, pairs_to_rescan);
    if (initially_screened_triplets.size() <
            static_cast<std::size_t>(original_scan_state.analysis_list_last + 1) ||
        trace_sub.size() < static_cast<std::size_t>(next_no + 1) ||
        actual_sequence_sizes.size() < static_cast<std::size_t>(next_no + 1)) {
        throw std::runtime_error("MakeAListISP3 input dimensions differ");
    }
    const int stride = next_no + 1;
    std::vector<std::array<int, 3>> output;
    // Literal DNA5 MakeAListISP3 at ProbDo=1.1. Worthwhile/ProgBinRead
    // collapse to the initial RDP screening bit in an RDP-only run.
    for (int triplet = 0;
         triplet <= original_scan_state.analysis_list_last; ++triplet) {
        if (initially_screened_triplets[triplet] == 0) continue;
        std::array<int, 3> sequences{
            original_scan_state.analysis_list[triplet * 3],
            original_scan_state.analysis_list[triplet * 3 + 1],
            original_scan_state.analysis_list[triplet * 3 + 2]};
        // The source intentionally does not reset Seqs inside this loop;
        // replacements therefore carry into the next winning-list member.
        for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
            const int member = list_member(candidate_list, winning_role, slot);
            const int base = member > permanent_next_no ?
                trace_sub[member] : member;
            for (int role = 0; role < 3; ++role) {
                if (sequences[role] == base) {
                    sequences[role] = member;
                    break;
                }
            }
            if (sequences[0] != member && sequences[1] != member &&
                sequences[2] != member) {
                continue;
            }
            if (actual_sequence_sizes[sequences[0]] <= minimum_sequence_size ||
                actual_sequence_sizes[sequences[1]] <= minimum_sequence_size ||
                actual_sequence_sizes[sequences[2]] <= minimum_sequence_size) {
                continue;
            }
            if (pairs_to_rescan[sequences[0] + sequences[1] * stride] == 1 &&
                pairs_to_rescan[sequences[0] + sequences[2] * stride] == 1 &&
                pairs_to_rescan[sequences[1] + sequences[2] * stride] == 1) {
                output.push_back(sequences);
            }
        }
    }
    return output;
}

std::vector<std::array<int, 3>> make_rdp_outer_scan_triplets(
    const RdpScanState& original_scan_state,
    const std::vector<unsigned char>& initially_screened_triplets,
    const int current_next_no, const int starting_next_no,
    const int winning_role, const std::array<int, 3>& candidate_last,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const int permanent_next_no, const int minimum_sequence_size,
    const std::vector<unsigned char>& pairs_to_rescan,
    const int subvalid_upper_bound, const std::vector<float>& subvalid) {
    (void)winning_role;
    const int pair_ub = pair_upper_bound(pairs_to_rescan);
    validate_pairs(pair_ub, pairs_to_rescan);
    const int pair_stride = pair_ub + 1;
    const int valid_stride = subvalid_upper_bound + 1;
    if (trace_sub.size() < static_cast<std::size_t>(current_next_no + 1) ||
        actual_sequence_sizes.size() <
            static_cast<std::size_t>(current_next_no + 1) ||
        subvalid.size() < static_cast<std::size_t>(valid_stride) * valid_stride) {
        throw std::runtime_error("MakeAListOSP input dimensions differ");
    }
    std::vector<std::array<int, 3>> output;
    // Literal DNA5 MakeAListOSP for BusyWithExcludes=0 and program zero.
    for (int triplet = 0;
         triplet <= original_scan_state.analysis_list_last; ++triplet) {
        if (initially_screened_triplets[triplet] == 0) continue;
        std::array<int, 3> sequences{
            original_scan_state.analysis_list[triplet * 3],
            original_scan_state.analysis_list[triplet * 3 + 1],
            original_scan_state.analysis_list[triplet * 3 + 2]};
        std::array<int, 3> bases{};
        bool go_on = true;
        for (int role = 0; role < 3; ++role) {
            if (sequences[role] >= static_cast<int>(trace_sub.size())) {
                go_on = false;
                break;
            }
            bases[role] = sequences[role] > starting_next_no ?
                trace_sub[sequences[role]] : sequences[role];
        }
        if (!go_on ||
            pairs_to_rescan[bases[0] + bases[1] * pair_stride] != 1 ||
            pairs_to_rescan[bases[0] + bases[2] * pair_stride] != 1 ||
            pairs_to_rescan[bases[1] + bases[2] * pair_stride] != 1) {
            continue;
        }
        // The source likewise carries replacements between fragment members.
        for (int member = current_next_no - candidate_last[winning_role];
             member <= current_next_no; ++member) {
            const int base = member > permanent_next_no ?
                trace_sub[member] : member;
            for (int role = 0; role < 3; ++role) {
                if (sequences[role] == base) {
                    sequences[role] = member;
                    break;
                }
            }
            if (sequences[0] != member && sequences[1] != member &&
                sequences[2] != member) {
                continue;
            }
            if (actual_sequence_sizes[sequences[0]] <= minimum_sequence_size ||
                actual_sequence_sizes[sequences[1]] <= minimum_sequence_size ||
                actual_sequence_sizes[sequences[2]] <= minimum_sequence_size) {
                continue;
            }
            // Native MakeAListOSP checks only pairs 0-2 and 0-1 here.
            if (subvalid[bases[0] + bases[2] * valid_stride] > 20.0F &&
                subvalid[bases[0] + bases[1] * valid_stride] > 20.0F) {
                output.push_back(sequences);
            }
        }
    }
    return output;
}

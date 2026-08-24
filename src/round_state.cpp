#include "round_state.hpp"

#include "MathFuncsDll.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace {

std::vector<float> make_zero_rep_collapsed_matrix(
    const int next_no,
    const int local_last,
    const int legacy_pass,
    const std::vector<int>& trace,
    const std::vector<int>& redo,
    const std::vector<float>& direct,
    std::vector<float>& adjusted,
    const std::vector<float>& temporary) {
    const int stride = next_no + 1;
    const int local_stride = local_last + 1;
    std::vector<float> collapsed(
        static_cast<std::size_t>(stride) * stride, 0.0F);

    // TestMoveInTreeAlt's Reps=0/BootFlag=1 branch maps tSAMat on pass
    // zero and tFAMat on pass one. Preserve its Y = x + 1 loop bound,
    // including the pass-one omission of the first local pair.
    for (int local_first = 0; local_first <= local_last; ++local_first) {
        const int first = trace[1 + local_first * 2];
        for (int local_second = legacy_pass + 1;
             local_second <= local_last; ++local_second) {
            const int second = trace[1 + local_second * 2];
            const float value = temporary[
                local_first + local_second * local_stride];
            collapsed[first + second * stride] = value;
            collapsed[second + first * stride] = value;
        }
    }

    MathFuncs::MyMathFuncs::CleanFCMat2P(
        next_no, next_no, next_no, collapsed.data(),
        const_cast<int*>(redo.data()));

    std::vector<int> best_sequence(stride, 0);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        best_sequence[sequence] = sequence;
        if (redo[sequence] == 0) {
            continue;
        }
        float minimum_distance = 10.0F;
        for (int candidate = 0; candidate <= next_no; ++candidate) {
            if (candidate != sequence && redo[candidate] == 0 &&
                direct[sequence + candidate * stride] < minimum_distance) {
                minimum_distance = direct[sequence + candidate * stride];
                best_sequence[sequence] = candidate;
            }
        }
    }

    // ReAddDistsB mutates both the adjusted and collapsed matrices in place.
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (redo[sequence] == 0) {
            continue;
        }
        const int closest = best_sequence[sequence];
        for (int other = 0; other <= next_no; ++other) {
            const auto forward = static_cast<std::size_t>(
                sequence + other * stride);
            const auto reverse = static_cast<std::size_t>(
                other + sequence * stride);
            if (other == sequence) {
                collapsed[forward] = 0.0F;
                adjusted[forward] = 0.0F;
                continue;
            }
            const auto replacement = static_cast<std::size_t>(
                closest + best_sequence[other] * stride);
            if (adjusted[forward] > adjusted[replacement]) {
                collapsed[forward] = collapsed[replacement];
                collapsed[reverse] = collapsed[forward];
                adjusted[forward] = adjusted[replacement];
                adjusted[reverse] = adjusted[forward];
            }
        }
    }
    return collapsed;
}

}  // namespace

RdpRoundPrefixState identify_rdp_round_prefix(
    const RdpScanState& scan_state,
    const RdpDistanceState& full_distance,
    const RdpRawEventState& events,
    const RdpRawEvent& selected,
    const std::vector<unsigned char>& missing_data,
    const int minimum_sequence_size) {
    RdpRoundPrefixState state;
    state.sequences = {
        selected.daughter, selected.minor_parent, selected.major_parent};
    const int next_no = scan_state.next_no;
    const int sequence_length = scan_state.sequence_length;
    const int stride = next_no + 1;

    state.breakpoint_distance.assign(3, 0.0F);
    state.remainder_distance.assign(3, 0.0F);
    auto mutable_sequences = std::vector<int>(
        state.sequences.begin(), state.sequences.end());
    MathFuncs::MyMathFuncs::UFDist(
        sequence_length, selected.beginning, selected.ending, next_no,
        const_cast<float*>(full_distance.valid_sites.data()),
        const_cast<float*>(full_distance.differences.data()),
        state.breakpoint_distance.data(), state.remainder_distance.data(),
        mutable_sequences.data(), sequence_length,
        const_cast<short*>(scan_state.sequence_data.data()));

    state.region_distance = build_rdp_distance_state(
        scan_state, selected.beginning, selected.ending, false);
    state.matrices = finish_rdp_event_distances(
        next_no, full_distance, state.region_distance);
    std::vector<unsigned char> minimum_pair{3, 3, 0};
    std::vector<unsigned char> sequence_pair(3, 0);
    constexpr std::array<int, 3> role_outlier{2, 1, 0};
    float minimum_background = 1000000.0F;
    float minimum_region = 1000000.0F;
    int pair = 0;
    for (int first = 0; first < 2; ++first) {
        for (int second = first + 1; second < 3; ++second) {
            const auto offset = static_cast<std::size_t>(
                state.sequences[first] + state.sequences[second] * stride);
            if (state.matrices.background[offset] < minimum_background) {
                minimum_background = state.matrices.background[offset];
                minimum_pair[0] = static_cast<unsigned char>(pair);
                sequence_pair = {
                    static_cast<unsigned char>(first),
                    static_cast<unsigned char>(second),
                    static_cast<unsigned char>(role_outlier[pair])};
            }
            if (state.matrices.event_region[offset] < minimum_region) {
                minimum_region = state.matrices.event_region[offset];
                minimum_pair[1] = static_cast<unsigned char>(pair);
            }
            ++pair;
        }
    }
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        state.matrices.background[sequence + sequence * stride] = 0.0F;
    }
    int region_half = selected.beginning < selected.ending
        ? static_cast<int>(std::nearbyint(
              (selected.ending - selected.beginning) / 2.0))
        : static_cast<int>(std::nearbyint(
              (selected.ending + sequence_length - selected.beginning) /
              2.0));
    region_half = std::min(region_half, 20);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        bool erase = false;
        for (const int selected_sequence : state.sequences) {
            if (selected_sequence <= next_no &&
                sequence != selected_sequence) {
                const auto offset = static_cast<std::size_t>(
                    selected_sequence + sequence * stride);
                if (full_distance.valid_sites[offset] -
                            state.region_distance.valid_sites[offset] <
                        minimum_sequence_size ||
                    state.region_distance.valid_sites[offset] < region_half) {
                    erase = true;
                    break;
                }
            }
        }
        if (erase) {
            for (int other = 0; other <= next_no; ++other) {
                state.matrices.background[sequence + other * stride] = 3.0F;
                state.matrices.background[other + sequence * stride] = 3.0F;
                state.matrices.event_region[sequence + other * stride] = 3.0F;
                state.matrices.event_region[other + sequence * stride] = 3.0F;
            }
        }
    }
    std::vector<int> minimums(stride, 0);
    std::vector<unsigned char> missing_pair(
        static_cast<std::size_t>(stride) * stride, 0);
    std::vector<int> background_total(stride, 0);
    std::vector<int> region_total(stride, 0);
    MathFuncs::MyMathFuncs::CheckMatrixP(
        minimums.data(), mutable_sequences.data(), next_no, region_half,
        minimum_sequence_size, next_no, missing_pair.data(), next_no,
        const_cast<float*>(full_distance.valid_sites.data()), next_no,
        state.region_distance.valid_sites.data(), next_no,
        state.matrices.background.data(), state.matrices.event_region.data(),
        background_total.data(), region_total.data());

    std::vector<int> outlier{2, 1, 0};
    std::vector<int> redo(stride, 0);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (state.matrices.background[sequence + sequence * stride] == 3.0F) {
            redo[sequence] = 1;
        }
    }
    int local_last = -1;
    std::vector<int> trace(static_cast<std::size_t>(2) * stride, 0);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (redo[sequence] == 0) {
            ++local_last;
            trace[sequence * 2] = local_last;
            trace[1 + local_last * 2] = sequence;
        }
    }
    state.background_adjusted.assign(
        static_cast<std::size_t>(stride) * stride, 0.0F);
    state.region_adjusted.assign(
        static_cast<std::size_t>(stride) * stride, 0.0F);
    const int local_stride = local_last + 1;
    std::vector<char> background_holder(
        static_cast<std::size_t>(local_stride) * 80 + 1, 0);
    std::vector<char> region_holder(
        static_cast<std::size_t>(local_stride) * 80 + 1, 0);
    std::vector<float> temporary_background(
        static_cast<std::size_t>(local_stride) * local_stride, 0.0F);
    std::vector<float> temporary_region(
        static_cast<std::size_t>(local_stride) * local_stride, 0.0F);
    const int name_length = std::max(
        2, static_cast<int>(std::to_string(local_last).size()));
    MathFuncs::MyMathFuncs::MakeNJTreesP2(
        1, local_last, next_no, mutable_sequences.data(), minimum_pair.data(),
        sequence_pair.data(), 3, name_length, sequence_length, 1,
        outlier.data(), trace.data(), next_no,
        state.matrices.background.data(), next_no,
        state.matrices.event_region.data(), next_no,
        state.background_adjusted.data(), next_no,
        state.region_adjusted.data(), redo.data(), background_holder.data(),
        region_holder.data(), temporary_background.data(),
        temporary_region.data());
    state.region_collapsed = make_zero_rep_collapsed_matrix(
        next_no, local_last, 0, trace, redo,
        state.matrices.event_region, state.region_adjusted,
        temporary_region);
    state.background_collapsed = make_zero_rep_collapsed_matrix(
        next_no, local_last, 1, trace, redo,
        state.matrices.background, state.background_adjusted,
        temporary_background);
    state.minimum_pair = {minimum_pair[0], minimum_pair[1]};
    state.sequence_pair = {
        sequence_pair[0], sequence_pair[1], sequence_pair[2]};

    state.breakpoint_flanks = make_rdp_breakpoint_flanks(
        scan_state, selected.beginning, selected.ending, state.sequences);
    state.starts = {
        state.breakpoint_flanks.positions[0], selected.beginning,
        state.breakpoint_flanks.positions[2], selected.ending + 1,
        selected.beginning, 0};
    state.ends = {
        selected.beginning - 1, state.breakpoint_flanks.positions[1],
        selected.ending, state.breakpoint_flanks.positions[3],
        selected.ending, 0};
    for (int region = 0; region < 4; ++region) {
        if (state.starts[region] > sequence_length) {
            state.starts[region] -= sequence_length;
        } else if (state.starts[region] < 1) {
            state.starts[region] += sequence_length;
        }
        if (state.ends[region] > sequence_length) {
            state.ends[region] -= sequence_length;
        } else if (state.ends[region] < 1) {
            state.ends[region] += sequence_length;
        }
    }
    constexpr std::array<int, 6> comparison{1, 0, 0, 2, 2, 1};
    state.summary_matrix.assign(static_cast<std::size_t>(9) * stride, 0.0);
    state.regional_distance_matrix.assign(
        static_cast<std::size_t>(45) * stride, 0.0);
    auto mutable_missing = missing_data;
    if (mutable_missing.size() <
        static_cast<std::size_t>(stride) * (sequence_length + 1)) {
        mutable_missing.resize(
            static_cast<std::size_t>(stride) * (sequence_length + 1), 0);
    }
    MathFuncs::MyMathFuncs::MakeSDMP2(
        next_no, sequence_length, state.starts.data(), state.ends.data(),
        mutable_sequences.data(), const_cast<int*>(comparison.data()),
        mutable_missing.data(),
        const_cast<short*>(scan_state.sequence_data.data()),
        state.summary_matrix.data(), state.regional_distance_matrix.data());
    std::vector<unsigned char> correlation_positions(2, 0);
    for (int target = 0; target < 3; ++target) {
        auto& matrix = state.correlation_matrices[target];
        matrix.assign(static_cast<std::size_t>(18) * stride, 0.0);
        MathFuncs::MyMathFuncs::FillRmat(
            target, next_no, 2, 5, 2, 4, next_no, matrix.data(),
            state.regional_distance_matrix.data(),
            correlation_positions.data());
    }
    state.correlation_decisions = finalize_rdp_correlations(
        next_no,
        calculate_rdp_correlations(
            next_no, state.sequences, comparison,
            state.correlation_matrices),
        state.sequences, comparison, state.summary_matrix,
        state.regional_distance_matrix);
    if (state.correlation_decisions.warnings[0] != 0 &&
        state.correlation_decisions.warnings[1] != 0) {
        state.correlation_decisions.warnings[2] = 0;
    }
    const std::array<int, 4> local_starts{
        state.starts[0], state.starts[1], state.starts[2], state.starts[3]};
    const std::array<int, 4> local_ends{
        state.ends[0], state.ends[1], state.ends[2], state.ends[3]};
    state.local_distance_panels = make_rdp_local_distance_panels(
        scan_state, local_starts, local_ends, state.sequences);
    apply_rdp_distance_warnings(
        next_no, state.sequences, state.local_distance_panels,
        state.correlation_decisions.warnings);
    state.good_comparisons = make_rdp_good_comparisons(
        scan_state, state.breakpoint_flanks.positions);
    const auto make_small = [&](const std::vector<float>& matrix) {
        std::vector<float> small(static_cast<std::size_t>(3) * stride, 0.0F);
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            for (int role = 0; role < 3; ++role) {
                small[role + sequence * 3] = matrix[
                    state.sequences[role] + sequence * stride];
            }
        }
        return small;
    };
    state.first_direct_small = make_small(state.matrices.background);
    state.region_direct_small = make_small(state.matrices.event_region);
    state.first_adjusted_small = make_small(state.background_adjusted);
    state.region_adjusted_small = make_small(state.region_adjusted);
    std::array<unsigned char, 2> final_pair{};
    std::array<float, 2> final_distance{1000000.0F, 1000000.0F};
    pair = 0;
    for (int first = 0; first < 2; ++first) {
        for (int second = first + 1; second < 3; ++second) {
            const auto offset = static_cast<std::size_t>(
                first + state.sequences[second] * 3);
            if (state.first_adjusted_small[offset] < final_distance[0]) {
                final_distance[0] = state.first_adjusted_small[offset];
                final_pair[0] = static_cast<unsigned char>(pair);
            }
            if (state.region_adjusted_small[offset] < final_distance[1]) {
                final_distance[1] = state.region_adjusted_small[offset];
                final_pair[1] = static_cast<unsigned char>(pair);
            }
            ++pair;
        }
    }
    state.role_lists = make_rdp_role_lists(final_pair);
    state.first_collapsed_small = make_small(state.background_collapsed);
    state.region_collapsed_small = make_small(state.region_collapsed);
    state.acceptable_correlations = make_rdp_acceptable_correlations(
        next_no, state.sequences, state.role_lists.inside,
        state.first_direct_small, state.region_direct_small,
        state.first_adjusted_small, state.region_adjusted_small,
        state.first_collapsed_small, state.region_collapsed_small);
    state.dont_redo.assign(static_cast<std::size_t>(3) * stride, 0);
    auto candidates = make_rdp_candidate_lists(
        next_no, state.good_comparisons, state.sequences,
        state.correlation_decisions, state.dont_redo,
        state.acceptable_correlations);
    state.actual_resolution = resolve_rdp_actual_events(
        sequence_length, next_no, state.sequences, comparison,
        state.starts, state.ends, state.correlation_decisions,
        std::move(candidates), state.dont_redo, events, next_no);
    return state;
}

RdpCompleteRoundState identify_rdp_complete_round(
    const RdpScanState& scan_state,
    const RdpDistanceState& full_distance,
    const RdpRawEventState& events,
    const RdpRawEvent& selected,
    const std::vector<unsigned char>& missing_data,
    int permanent_next_no,
    const int minimum_sequence_size) {
    RdpCompleteRoundState state;
    state.prefix = identify_rdp_round_prefix(
        scan_state, full_distance, events, selected, missing_data,
        minimum_sequence_size);
    const int next_no = scan_state.next_no;
    if (permanent_next_no < 0) permanent_next_no = next_no;
    const auto& prefix = state.prefix;
    const auto& sequences = prefix.sequences;
    constexpr std::array<int, 6> comparison{1, 0, 0, 2, 2, 1};

    const auto initial_compatibility = evaluate_rdp_tree_compatibility(
        next_no, sequences, comparison,
        prefix.actual_resolution.inversion_penalty,
        prefix.actual_resolution.candidates.last,
        prefix.actual_resolution.candidates.list, prefix.good_comparisons,
        prefix.background_adjusted, prefix.region_adjusted);
    const auto tied = [](const std::array<int, 3>& values) {
        return values[0] == values[1] && values[0] == values[2];
    };
    std::vector<float> background_secondary;
    std::vector<float> region_secondary;
    if (tied(initial_compatibility.background_compatibility)) {
        background_secondary = prefix.background_collapsed;
    }
    if (tied(initial_compatibility.region_compatibility)) {
        region_secondary = prefix.region_collapsed;
    }
    state.tree_compatibility = run_rdp_tree_compatibility_flow(
        scan_state.sequence_length, next_no, selected.beginning,
        selected.ending, sequences, comparison,
        prefix.actual_resolution.inversion_penalty,
        prefix.actual_resolution.candidates.last,
        prefix.actual_resolution.candidates.list, prefix.good_comparisons,
        prefix.background_adjusted, prefix.region_adjusted,
        background_secondary, region_secondary, events);

    state.score_filters[0] = make_rdp_score_filter(
        next_no, sequences, prefix.first_direct_small,
        prefix.first_direct_small, prefix.region_direct_small);
    state.phylpro[0] = make_rdp_phylpro_scores(
        next_no, 1.0e-14, state.score_filters[0], sequences,
        prefix.matrices.background, prefix.matrices.event_region);
    state.phylpro[1] = make_rdp_phylpro_scores(
        next_no, 1.0e-14, state.score_filters[0], sequences,
        prefix.background_adjusted, prefix.region_adjusted);
    state.score_filters[1] = make_rdp_score_filter(
        next_no, sequences, prefix.first_direct_small,
        prefix.first_direct_small, prefix.region_direct_small);
    state.phylpro[2] = make_rdp_phylpro_scores(
        next_no, 1.0e-14, state.score_filters[1], sequences,
        prefix.background_collapsed, prefix.region_collapsed);

    for (int role = 0; role < 3; ++role) {
        state.triplet_groups[role] = make_rdp_triplet_groups(
            role, next_no, sequences, comparison,
            prefix.first_adjusted_small);
        state.triplet_scores[role] = make_rdp_triplet_tree_score(
            role, next_no, sequences, prefix.first_adjusted_small,
            prefix.region_adjusted_small, state.triplet_groups[role]);
    }

    const std::array<int, 3> pattern_starts{
        prefix.starts[0], prefix.starts[2], prefix.starts[4]};
    const std::array<int, 3> pattern_ends{
        prefix.ends[1], prefix.ends[3], prefix.ends[4]};
    state.pattern = check_rdp_sequence_patterns(
        scan_state.sequence_length, next_no, sequences, pattern_starts,
        pattern_ends, comparison, scan_state.sequence_data,
        prefix.actual_resolution.acceptable_sequences);
    state.consensus_candidates.candidate_last =
        prefix.actual_resolution.candidates.last;
    state.consensus_candidates.candidate_list =
        prefix.actual_resolution.candidates.list;
    state.consensus_candidates.acceptable_sequences =
        state.pattern.acceptable_sequences;

    // Module3 tests the first FAMat/SAMat MakeRCompat families here. The
    // later fallback/set families feed MakeConsensusC but do not control
    // RetrimFlag.
    state.consensus_retrimmed =
        std::all_of(initial_compatibility.background_compatibility.begin(),
                    initial_compatibility.background_compatibility.end(),
                    [](const int value) { return value > 0; }) &&
        std::all_of(initial_compatibility.region_compatibility.begin(),
                    initial_compatibility.region_compatibility.end(),
                    [](const int value) { return value > 0; });
    if (state.consensus_retrimmed) {
        auto trimmed = run_rdp_final_trim_candidate_maintenance(
            next_no, sequences, comparison, prefix.minimum_pair,
            prefix.role_lists.inside,
            prefix.correlation_decisions.warnings,
            prefix.actual_resolution.unfound,
            prefix.actual_resolution.correlations.correlations.correlation,
            prefix.actual_resolution.correlations.correlations.inversion,
            prefix.local_distance_panels, prefix.first_adjusted_small,
            prefix.region_adjusted_small, prefix.first_collapsed_small,
            prefix.region_collapsed_small,
            prefix.actual_resolution.candidates.last,
            prefix.actual_resolution.candidates.list,
            state.pattern.acceptable_sequences);
        trimmed.acceptable_sequences = calculate_rdp_match_evidence(
            scan_state.sequence_length, next_no, selected.beginning,
            selected.ending, sequences, comparison, scan_state.sequence_data,
            trimmed.acceptable_sequences, true);
        state.consensus_candidates = make_rdp_consensus_candidates(
            next_no, sequences, comparison,
            prefix.correlation_decisions.warnings,
            prefix.actual_resolution.correlations.correlations.correlation,
            prefix.actual_resolution.correlations.correlations.inversion,
            prefix.matrices.background, prefix.matrices.event_region,
            prefix.background_adjusted, prefix.region_adjusted,
            prefix.first_direct_small, prefix.region_direct_small,
            prefix.first_adjusted_small, prefix.region_adjusted_small,
            prefix.first_collapsed_small, prefix.region_collapsed_small,
            std::move(trimmed), true);
    }

    const std::vector<unsigned char> maximum_distance_mask(
        static_cast<std::size_t>(next_no + 1), 0);
    state.maximum_distance = calculate_rdp_maximum_distances(
        scan_state.sequence_length, next_no, sequences, selected.beginning,
        selected.ending, scan_state.sequence_data, maximum_distance_mask);
    state.split_distances = calculate_rdp_split_distances(
        next_no, sequences, prefix.role_lists.inside,
        prefix.first_adjusted_small, prefix.matrices.background,
        prefix.matrices.event_region, state.score_filters[0]);
    state.simple_distances = calculate_rdp_simple_distances(
        next_no, sequences, prefix.role_lists.inside,
        prefix.actual_resolution.candidates.last,
        prefix.actual_resolution.candidates.list, prefix.matrices.background,
        prefix.matrices.event_region);
    if (prefix.minimum_pair[0] != prefix.minimum_pair[1]) {
        state.outlier_checks = calculate_rdp_outlier_checks(
            next_no, sequences, prefix.role_lists.inside,
            prefix.first_adjusted_small, prefix.region_adjusted_small);
        state.list_correlations = calculate_rdp_list_correlations(
            next_no, sequences, prefix.role_lists.inside,
            prefix.correlation_decisions.warnings,
            prefix.actual_resolution.candidates.last,
            prefix.actual_resolution.candidates.list,
            prefix.actual_resolution.correlations.correlations.inversion,
            prefix.actual_resolution.correlations.correlations
                .tested_correlation,
            prefix.first_adjusted_small, prefix.region_adjusted_small);
    }
    state.bad_distances = calculate_rdp_bad_distances(
        next_no, sequences, comparison,
        prefix.actual_resolution.candidates.last,
        prefix.actual_resolution.candidates.list,
        prefix.actual_resolution.unfound,
        prefix.actual_resolution.correlations.correlations.correlation,
        prefix.first_adjusted_small, prefix.local_distance_panels);

    if (state.consensus_retrimmed) {
        const auto run_post_trim = [&](const std::vector<float>& matrix,
                                       const std::array<double, 3>& distances,
                                       std::array<int, 3>& compatibility) {
            std::array<int, 3> reverse{};
            for (int role = 0; role < 3; ++role) {
                std::array<int, 3> nonrecombinant_last{};
                std::vector<int> nonrecombinant_list(
                    static_cast<std::size_t>(3) * (next_no + 1), 0);
                make_rdp_tree_compatibility_call(
                    next_no, sequences, comparison, role,
                    prefix.actual_resolution.inversion_penalty,
                    state.consensus_candidates.candidate_last,
                    state.consensus_candidates.candidate_list,
                    prefix.good_comparisons, matrix, distances,
                    compatibility, reverse, nonrecombinant_last,
                    nonrecombinant_list);
            }
        };
        run_post_trim(prefix.background_adjusted,
                      state.tree_compatibility.calls[0].list_distances,
                      state.post_trim_background);
        run_post_trim(prefix.region_adjusted,
                      state.tree_compatibility.calls[3].list_distances,
                      state.post_trim_region);
    }

    RdpConsensusInputs consensus_inputs;
    consensus_inputs.next_no = next_no;
    consensus_inputs.permanent_next_no = permanent_next_no;
    consensus_inputs.comparison_matrix = comparison;
    consensus_inputs.list_correlation = state.list_correlations.mismatches;
    consensus_inputs.simple_distance_strength =
        state.simple_distances.strengths;
    consensus_inputs.simple_distance_score = state.simple_distances.scores;
    consensus_inputs.phylpro = state.phylpro[0].scores;
    consensus_inputs.phylpro_secondary = state.phylpro[1].scores;
    consensus_inputs.phylpro_collapsed = state.phylpro[2].scores;
    consensus_inputs.subtree_score = state.phylpro[0].sub_distance_scores;
    consensus_inputs.split_distance = state.split_distances.distances;
    consensus_inputs.outlier_index = state.split_distances.outlier_index;
    consensus_inputs.subtree_phylpro = state.phylpro[0].sub_scores;
    consensus_inputs.subtree_score_secondary =
        state.phylpro[1].sub_distance_scores;
    consensus_inputs.subtree_phylpro_secondary = state.phylpro[1].sub_scores;
    consensus_inputs.compatibility = state.tree_compatibility.background;
    consensus_inputs.compatibility_secondary =
        state.tree_compatibility.background_secondary;
    consensus_inputs.compatibility_tertiary =
        state.tree_compatibility.background_sets;
    consensus_inputs.compatibility_quaternary =
        state.tree_compatibility.background_secondary_sets;
    consensus_inputs.region_compatibility = state.tree_compatibility.region;
    consensus_inputs.region_compatibility_secondary =
        state.tree_compatibility.region_secondary;
    consensus_inputs.region_compatibility_tertiary =
        state.tree_compatibility.region_sets;
    consensus_inputs.region_compatibility_quaternary =
        state.tree_compatibility.region_secondary_sets;
    consensus_inputs.post_trim_compatibility = state.post_trim_background;
    consensus_inputs.post_trim_region_compatibility = state.post_trim_region;
    consensus_inputs.triplet_score = state.triplet_scores;
    consensus_inputs.bad_distances = state.bad_distances;
    for (int role = 0; role < 3; ++role) {
        consensus_inputs.outside_list[role] =
            prefix.role_lists.outside[role];
    }
    consensus_inputs.list_correlation_secondary =
        state.list_correlations.expected_strength;
    consensus_inputs.list_correlation_tertiary =
        state.list_correlations.absent_strength;
    consensus_inputs.outlier_check = state.outlier_checks;
    consensus_inputs.maximum_distance =
        state.maximum_distance.maximum_distances;
    consensus_inputs.ranks = state.simple_distances.ranks;
    state.consensus = make_rdp_consensus(std::move(consensus_inputs));
    state.final_candidates = apply_rdp_strict_group_constraints(
        next_no, sequences, comparison, prefix.matrices.background,
        prefix.matrices.event_region, prefix.first_direct_small,
        prefix.region_direct_small, prefix.first_adjusted_small,
        prefix.region_adjusted_small, state.consensus_candidates);
    return state;
}

#include "identification_state.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

double vb_round(const double value, const double scale) {
    return std::nearbyint(value * scale) / scale;
}

double rdp_zip(const double q, const double initial, const int last,
               const double offset) {
    double factor = initial;
    double term = 1.0;
    double result = term;
    for (int index = static_cast<int>(factor); index <= last; index += 2) {
        term *= q * factor / (factor - offset);
        result += term;
        factor += 2.0;
    }
    return result;
}

double rdp_buzz(double value, const int degrees_of_freedom,
                const double half_pi) {
    value = std::abs(value);
    const double angle = std::atan(value / std::sqrt(degrees_of_freedom));
    if (degrees_of_freedom == 1) return 1.0 - angle / half_pi;
    return 1.0 - std::sin(angle) *
        rdp_zip(std::cos(angle) * std::cos(angle), 1.0,
                degrees_of_freedom - 3, -1.0);
}

double rdp_ttest_probability(const double value,
                             const int degrees_of_freedom) {
    const double half_pi = 4.0 * std::atan(1.0) / 2.0;
    return 1.0 - rdp_buzz(value, degrees_of_freedom, half_pi) / 2.0;
}

}  // namespace

RdpBreakpointFlanks make_rdp_breakpoint_flanks(
    const RdpScanState& scan_state, const int beginning, const int ending,
    const std::array<int, 3>& sequences, const int variable_site_target,
    const int total_site_target) {
    if (beginning < 1 || beginning > scan_state.sequence_length ||
        ending < 1 || ending > scan_state.sequence_length) {
        throw std::runtime_error("MakeBPosLR interval is outside alignment");
    }
    for (const int sequence : sequences) {
        if (sequence < 0 || sequence > scan_state.next_no) {
            throw std::runtime_error("MakeBPosLR sequence is outside alignment");
        }
    }

    RdpBreakpointFlanks result;
    const int length = scan_state.sequence_length;
    const int stride = length + 1;
    const std::array<int, 4> increments{-1, 1, -1, 1};
    for (int flank = 0; flank < 4; ++flank) {
        int position;
        int stop;
        if (flank == 0) {
            position = beginning - 1;
            stop = ending + 1;
        } else if (flank == 1) {
            position = beginning;
            stop = ending;
        } else if (flank == 2) {
            position = ending;
            stop = beginning;
        } else {
            position = ending + 1;
            stop = beginning - 1;
        }
        if (position > length) position -= length;
        if (stop > length) stop -= length;

        int variable_count = 0;
        int total_count = 0;
        while (position != stop) {
            const short first = scan_state.sequence_data[
                position + sequences[0] * stride];
            const short second = scan_state.sequence_data[
                position + sequences[1] * stride];
            const short third = scan_state.sequence_data[
                position + sequences[2] * stride];
            if (first != 46 && second != 46 && third != 46) {
                ++total_count;
                if ((first != second || first != third) &&
                    (first == second || first == third || second == third)) {
                    ++variable_count;
                    if (variable_count == variable_site_target &&
                        total_count > total_site_target) {
                        break;
                    }
                }
            }
            position += increments[flank];
            if (position < 1) {
                position = length;
            } else if (position > length) {
                position = 1;
            }
        }
        result.positions[flank] = position;
        result.informative_counts[flank] =
            variable_count == 0 ? 1.0 : variable_count;
    }
    return result;
}

RdpCorrelationState calculate_rdp_correlations(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<std::vector<double>, 3>& regional_matrices,
    const double minimum_offset) {
    const int sequence_count = next_no + 1;
    const auto regional_size = static_cast<std::size_t>(18) * sequence_count;
    for (const auto& matrix : regional_matrices) {
        if (matrix.size() != regional_size) {
            throw std::runtime_error("CalCR regional matrix dimensions differ");
        }
    }
    RdpCorrelationState state;
    state.correlation.assign(static_cast<std::size_t>(9) * sequence_count, 0);
    state.inversion.assign(static_cast<std::size_t>(9) * sequence_count, 0);
    state.tested_correlation.assign(
        static_cast<std::size_t>(45) * sequence_count, 0);

    for (int target = 0; target < 3; ++target) {
        const auto& regional = regional_matrices[target];
        for (int region = 0; region < 3; ++region) {
            const int region_offset = region * 3;
            state.correlation[target + region_offset +
                sequences[target] * 9] = 1.0F;
            state.correlation[target + region_offset +
                sequences[comparison_matrix[target]] * 9] = 0.0F;
            state.correlation[target + region_offset +
                sequences[comparison_matrix[target + 3]] * 9] = 0.0F;
        }

        const int target_offset = sequences[target] * 18;
        for (int region = 0; region < 3; ++region) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                const int sequence_offset = sequence * 9;
                const int regional_sequence_offset = sequence * 18;
                state.correlation[target + region * 3 + sequence_offset] = 0;
                state.inversion[target + region * 3 + sequence_offset] = 0;
                if (sequence == sequences[target]) {
                    state.tested_correlation[
                        target + region * 3 + sequence * 45] = 1;
                    state.correlation[
                        target + region * 3 + sequence_offset] = 1;
                    continue;
                }

                for (int permutation = 0; permutation < 6; ++permutation) {
                    double sum_x = 0;
                    double sum_y = 0;
                    double sum_xy = 0;
                    double sum_x2 = 0;
                    double sum_y2 = 0;
                    const int stored_permutation =
                        permutation == 5 ? 4 : permutation;
                    const auto tested_index = static_cast<std::size_t>(target) +
                        region * 3 + stored_permutation * 9 + sequence * 45;
                    if (permutation <= 4) {
                        state.tested_correlation[tested_index] = 0;
                    }
                    int sample_count = 6;
                    for (int category = 0; category < 6; ++category) {
                        int permuted_category = category;
                        if (permutation == 1) {
                            if (category == 0 || category == 3) {
                                permuted_category = category + 1;
                            } else if (category != 2 && category != 5) {
                                permuted_category = category - 1;
                            }
                        } else if (permutation == 2) {
                            if (category == 2 || category == 5) {
                                permuted_category = category - 2;
                            } else if (category == 0 || category == 3) {
                                permuted_category = category + 2;
                            }
                        } else if (permutation == 3) {
                            if (category == 2 || category == 5) {
                                permuted_category = category - 1;
                            } else if (category != 0 && category != 3) {
                                permuted_category = category + 1;
                            }
                        } else if (permutation == 4) {
                            if (category == 2 || category == 5) {
                                permuted_category = category - 2;
                            } else {
                                permuted_category = category + 1;
                            }
                        } else if (permutation == 5) {
                            if (category == 2 || category == 5) {
                                permuted_category = category - 1;
                            } else if (category == 0 || category == 3) {
                                permuted_category = category + 2;
                            } else {
                                permuted_category = category - 1;
                            }
                        }
                        state.intermediate[0] = regional[
                            region + category * 3 + target_offset];
                        state.intermediate[1] = regional[
                            region + permuted_category * 3 +
                            regional_sequence_offset];
                        if (state.intermediate[0] > 4 ||
                            state.intermediate[1] > 4) {
                            --sample_count;
                            break;
                        }
                        sum_x += state.intermediate[0];
                        sum_y += state.intermediate[1];
                        sum_xy +=
                            state.intermediate[0] * state.intermediate[1];
                        sum_x2 +=
                            state.intermediate[0] * state.intermediate[0];
                        sum_y2 +=
                            state.intermediate[1] * state.intermediate[1];
                    }

                    bool calculated = false;
                    if (sample_count == 6) {
                        if (sum_x2 > 0 && sum_y2 > 0) {
                            double denominator_x =
                                6.0 * sum_x2 - sum_x * sum_x;
                            if (denominator_x > 0.000000001) {
                                denominator_x = std::pow(denominator_x, 0.5);
                                double denominator_y =
                                    6.0 * sum_y2 - sum_y * sum_y;
                                if (denominator_y > 0.000000001) {
                                    calculated = true;
                                    denominator_y =
                                        std::pow(denominator_y, 0.5);
                                    const double denominator =
                                        denominator_x * denominator_y;
                                    const float stored_denominator =
                                        static_cast<float>(denominator);
                                    double extended_correlation;
                                    if (stored_denominator > 0) {
                                        extended_correlation =
                                            (6.0 * sum_xy - sum_x * sum_y) /
                                            stored_denominator;
                                    } else {
                                        extended_correlation = 0;
                                    }
                                    extended_correlation += minimum_offset;
                                    const float correlation =
                                        static_cast<float>(
                                            extended_correlation);
                                    if (permutation == 5) {
                                        if (state.tested_correlation[
                                                tested_index] <
                                            extended_correlation) {
                                            state.tested_correlation[
                                                tested_index] = correlation;
                                        }
                                    } else {
                                        state.tested_correlation[
                                            tested_index] = correlation;
                                    }
                                    const auto result_index =
                                        static_cast<std::size_t>(target) +
                                        region * 3 + sequence_offset;
                                    if (stored_permutation == 0) {
                                        state.correlation[result_index] =
                                            correlation;
                                        state.inversion[result_index] = 0;
                                    } else if (state.tested_correlation[
                                                       target + region * 3 +
                                                       sequence * 45] < 0.83F &&
                                               state.correlation[result_index] <
                                                   extended_correlation) {
                                        state.correlation[result_index] =
                                            correlation;
                                        state.inversion[result_index] =
                                            static_cast<float>(
                                                stored_permutation);
                                    }
                                }
                            }
                        }
                        if (!calculated) {
                            const bool constant =
                                6.0 * sum_y2 == sum_y * sum_y ||
                                6.0 * sum_x2 == sum_x * sum_x;
                            if (stored_permutation == 0) {
                                const auto result_index =
                                    static_cast<std::size_t>(target) +
                                    region * 3 + sequence_offset;
                                state.correlation[result_index] =
                                    constant ? 1.0F : 0.0F;
                                state.tested_correlation[tested_index] =
                                    constant ? 1.0F : 0.0F;
                                if (constant) state.inversion[result_index] = 0;
                            } else {
                                state.tested_correlation[tested_index] =
                                    constant ? 1.0F : 0.0F;
                            }
                        }
                    } else {
                        state.tested_correlation[tested_index] = 0;
                        if (stored_permutation == 0) {
                            state.correlation[target + region * 3 +
                                sequence_offset] = 0;
                        }
                    }
                }
            }
        }
        state.results[target] = 1.0;
    }
    return state;
}

RdpCorrelationDecisionState finalize_rdp_correlations(
    const int next_no, RdpCorrelationState correlations,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<double>& summary_matrix,
    const std::vector<double>& regional_distance_matrix) {
    const int sequence_count = next_no + 1;
    if (summary_matrix.size() !=
            static_cast<std::size_t>(9) * sequence_count ||
        regional_distance_matrix.size() !=
            static_cast<std::size_t>(45) * sequence_count) {
        throw std::runtime_error("MakeProperRCorr matrix dimensions differ");
    }
    const auto correlation_index = [](const int role, const int region,
                                      const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto tested_index = [](const int role, const int region,
                                 const int permutation, const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 +
            permutation * 9 + sequence * 45;
    };
    const auto distance_index = [sequence_count](
                                    const int role, const int region,
                                    const int sequence, const int category) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 15 +
            static_cast<std::size_t>(category) * 15 * sequence_count;
    };

    // Module2.MakeProperRCorr: discard isolated direct correlations whose
    // anchor-to-candidate distance is worse than either comparison role.
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            int direct_count = 0;
            for (int region = 0; region < 2; ++region) {
                const auto index = correlation_index(role, region, sequence);
                if (correlations.correlation[index] > 0.83F &&
                    correlations.inversion[index] == 0.0F) {
                    ++direct_count;
                }
            }
            const auto combined = correlation_index(role, 2, sequence);
            if (direct_count == 1 ||
                (direct_count == 0 &&
                 correlations.correlation[combined] > 0.83F &&
                 correlations.inversion[combined] == 0.0F)) {
                for (int region = 0; region < 3; ++region) {
                    const auto index = correlation_index(role, region, sequence);
                    if (correlations.correlation[index] > 0.83F &&
                        correlations.correlation[index] < 0.99F &&
                        correlations.inversion[index] == 0.0F &&
                        (summary_matrix[index] > summary_matrix[
                            correlation_index(
                                comparison_matrix[role], region, sequence)] ||
                         summary_matrix[index] > summary_matrix[
                            correlation_index(
                                comparison_matrix[role + 3], region,
                                sequence)])) {
                        correlations.correlation[index] = 0.0F;
                        correlations.inversion[index] = 0.0F;
                    }
                }
            }
        }
    }

    RdpCorrelationDecisionState result;
    // The two breakpoint-side warning tests compare the strictly greatest
    // triplet category after VB6 CLng rounding to five decimal places.
    for (int role = 0; role < 3; ++role) {
        for (int side = 0; side < 2; ++side) {
            if (result.warnings[side] != 0) continue;
            std::array<double, 6> values{};
            const int first_region = side == 0 ? 0 : 2;
            const int second_region = side == 0 ? 1 : 3;
            for (int category = 0; category < 3; ++category) {
                values[category] = vb_round(
                    regional_distance_matrix[distance_index(
                        role, first_region, sequences[role], category)],
                    100000.0);
                values[category + 3] = vb_round(
                    regional_distance_matrix[distance_index(
                        role, second_region, sequences[role], category)],
                    100000.0);
            }
            for (int category = 0; category < 3; ++category) {
                const int other1 = (category + 1) % 3;
                const int other2 = (category + 2) % 3;
                if (values[category] > values[other1] &&
                    values[category] > values[other2] &&
                    values[category + 3] > values[other1 + 3] &&
                    values[category + 3] > values[other2 + 3]) {
                    result.warnings[side] = 1;
                }
            }
        }
    }
    if (result.warnings[0] != 0 || result.warnings[1] != 0) {
        result.warnings[2] = 1;
    }

    for (float& value : correlations.tested_correlation) {
        value = static_cast<float>(vb_round(value, 1000000.0));
    }
    for (int role = 0; role < 3; ++role) {
        for (int region = 0; region < 3; ++region) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                const auto index = correlation_index(role, region, sequence);
                if (correlations.inversion[index] > 0.0F) {
                    float winner = -1.0F;
                    for (int permutation = 0; permutation < 5;
                         ++permutation) {
                        const float candidate = correlations.tested_correlation[
                            tested_index(role, region, permutation, sequence)];
                        if (candidate > winner) {
                            correlations.inversion[index] =
                                static_cast<float>(permutation);
                            winner = candidate;
                        }
                    }
                }
                if (correlations.inversion[index] == 3.0F) {
                    correlations.inversion[index] = 2.0F;
                }
            }
        }
    }
    result.correlations = std::move(correlations);
    return result;
}

std::vector<float> make_rdp_local_distance_panels(
    const RdpScanState& scan_state, const std::array<int, 4>& starts,
    const std::array<int, 4>& ends,
    const std::array<int, 3>& sequences) {
    const int sequence_count = scan_state.next_no + 1;
    const int alignment_stride = scan_state.sequence_length + 1;
    std::vector<float> result(
        static_cast<std::size_t>(12) * sequence_count, 0.0F);
    for (int panel = 0; panel < 4; ++panel) {
        for (int role = 0; role < 3; ++role) {
            const int selected_offset = sequences[role] * alignment_stride;
            for (int sequence = 0; sequence <= scan_state.next_no;
                 ++sequence) {
                const int candidate_offset = sequence * alignment_stride;
                int valid = 0;
                int differences = 0;
                int position = starts[panel];
                while (true) {
                    const short selected = scan_state.sequence_data[
                        selected_offset + position];
                    const short candidate = scan_state.sequence_data[
                        candidate_offset + position];
                    if (selected != 46 && candidate != 46) {
                        ++valid;
                        if (selected != candidate) ++differences;
                    }
                    if (position == ends[panel]) break;
                    ++position;
                    if (position > scan_state.sequence_length) position = 1;
                }
                float distance = 10.0F;
                if (valid > 0) {
                    const float identity = static_cast<float>(valid - differences) /
                        static_cast<float>(valid);
                    if (identity > 0.25F) {
                        const float transformed =
                            (4.0F * identity - 1.0F) / 3.0F;
                        distance = -0.75F * std::log(transformed);
                    }
                }
                const auto index = static_cast<std::size_t>(panel) + role * 4 +
                    sequence * 12;
                result[index] = valid > 10
                    ? static_cast<float>(vb_round(distance, 10000.0))
                    : 3.0F;
            }
        }
    }
    return result;
}

void apply_rdp_distance_warnings(
    const int next_no, const std::array<int, 3>& sequences,
    const std::vector<float>& local_distance_panels,
    std::array<unsigned char, 3>& warnings) {
    const int sequence_count = next_no + 1;
    if (local_distance_panels.size() !=
        static_cast<std::size_t>(12) * sequence_count) {
        throw std::runtime_error("MakeDMatS matrix dimensions differ");
    }
    const auto at = [&](const int panel, const int role, const int sequence) {
        return local_distance_panels[panel + role * 4 + sequence * 12];
    };
    for (int panel = 0; panel < 4; ++panel) {
        double total = at(panel, 0, sequences[1]) +
            at(panel, 0, sequences[2]) + at(panel, 2, sequences[1]);
        const int side = panel < 2 ? 0 : 1;
        if (total > 0.0) {
            total = 2.0 / total;
            if (1.0 - at(panel, 0, sequences[1]) * total < 0.4 &&
                1.0 - at(panel, 0, sequences[2]) * total < 0.4 &&
                1.0 - at(panel, 2, sequences[1]) * total < 0.4) {
                warnings[side] = 1;
            }
        } else {
            warnings[side] = 1;
        }
    }
    if (warnings[0] != 0 && warnings[1] != 0) {
        warnings[2] = 0;
    } else if (warnings[0] != 0 || warnings[1] != 0) {
        warnings[2] = 1;
    }
}

std::vector<int> make_rdp_good_comparisons(
    const RdpScanState& scan_state,
    const std::array<int, 4>& breakpoint_flanks) {
    const int sequence_count = scan_state.next_no + 1;
    const int stride = scan_state.sequence_length + 1;
    std::vector<int> result(static_cast<std::size_t>(2) * sequence_count, 0);
    for (int side = 0; side < 2; ++side) {
        const int start = breakpoint_flanks[side * 2];
        const int end = breakpoint_flanks[side * 2 + 1];
        for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
            int valid = 0;
            int position = start;
            while (true) {
                if (scan_state.sequence_data[position + sequence * stride] !=
                    46) {
                    ++valid;
                }
                if (position == end) break;
                ++position;
                if (position > scan_state.sequence_length) position = 1;
            }
            if (valid > 10) result[sequence + side * sequence_count] = 1;
        }
    }
    return result;
}

RdpRoleLists make_rdp_role_lists(
    const std::array<unsigned char, 2>& minimum_pair) {
    RdpRoleLists result;
    const int first = minimum_pair[0];
    const int second = minimum_pair[1];
    if (first == 0 && second == 1) {
        result.inside = {1, 0, 2}; result.outside = {1, 0, 2};
    } else if (first == 0 && second == 2) {
        result.inside = {0, 1, 2}; result.outside = {0, 1, 2};
    } else if (first == 1 && second == 0) {
        result.inside = {2, 0, 1}; result.outside = {1, 2, 0};
    } else if (first == 1 && second == 2) {
        result.inside = {0, 2, 1}; result.outside = {0, 2, 1};
    } else if (first == 2 && second == 0) {
        result.inside = {2, 1, 0}; result.outside = {2, 1, 0};
    } else if (first == 2 && second == 1) {
        result.inside = {1, 2, 0}; result.outside = {2, 0, 1};
    }
    return result;
}

std::vector<unsigned char> make_rdp_acceptable_correlations(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside,
    const std::vector<float>& first_direct,
    const std::vector<float>& second_direct,
    std::vector<float>& first_adjusted,
    std::vector<float>& second_adjusted,
    std::vector<float>& first_collapsed,
    std::vector<float>& second_collapsed) {
    const int sequence_count = next_no + 1;
    const auto expected = static_cast<std::size_t>(3) * sequence_count;
    if (first_direct.size() != expected || second_direct.size() != expected ||
        first_adjusted.size() != expected || second_adjusted.size() != expected ||
        first_collapsed.size() != expected || second_collapsed.size() != expected) {
        throw std::runtime_error("MakeACOR matrix dimensions differ");
    }
    float minimum_first = 1000.0F;
    float minimum_second = 1000.0F;
    int first_role = -1;
    int first_sequence = -1;
    int second_role = -1;
    int second_sequence = -1;
    const auto at = [](const int row, const int column) {
        return static_cast<std::size_t>(row) +
            static_cast<std::size_t>(column) * 3;
    };
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (sequence == sequences[role]) continue;
            const auto index = at(role, sequence);
            if (minimum_first > first_adjusted[index]) {
                minimum_first = first_adjusted[index];
                first_role = role;
                first_sequence = sequence;
            }
            if (minimum_second > second_adjusted[index]) {
                minimum_second = second_adjusted[index];
                second_role = role;
                second_sequence = sequence;
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (sequence == sequences[role]) continue;
            const auto index = at(role, sequence);
            if (first_direct[index] ==
                    first_direct[at(first_role, first_sequence)] &&
                first_direct[index] != 3.0F) {
                first_adjusted[index] = minimum_first;
                first_collapsed[index] = minimum_first;
            }
            if (second_direct[index] ==
                    second_direct[at(second_role, second_sequence)] &&
                second_direct[index] != 3.0F) {
                second_adjusted[index] = minimum_second;
                second_collapsed[index] = minimum_second;
            }
        }
    }
    std::vector<unsigned char> result(
        static_cast<std::size_t>(3) * sequence_count, 0);
    const float first_threshold = first_adjusted[
        at(inside[2], sequences[inside[0]])];
    const float second_threshold = second_adjusted[
        at(inside[1], sequences[inside[0]])];
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (first_adjusted[at(inside[0], sequence)] < first_threshold ||
            second_adjusted[at(inside[1], sequence)] < second_threshold) {
            result[inside[0] + sequence * 3] = 1;
            result[inside[1] + sequence * 3] = 1;
        }
        if (first_adjusted[at(inside[2], sequence)] < first_threshold ||
            second_adjusted[at(inside[2], sequence)] < second_threshold) {
            result[inside[2] + sequence * 3] = 1;
        }
    }
    return result;
}

RdpCandidateLists make_rdp_candidate_lists(
    const int next_no, const std::vector<int>& good_comparisons,
    const std::array<int, 3>& sequences,
    RdpCorrelationDecisionState& decision,
    const std::vector<unsigned char>& dont_redo,
    const std::vector<unsigned char>& acceptable_correlations) {
    const int sequence_count = next_no + 1;
    if (good_comparisons.size() !=
            static_cast<std::size_t>(2) * sequence_count ||
        dont_redo.size() != static_cast<std::size_t>(3) * sequence_count ||
        acceptable_correlations.size() !=
            static_cast<std::size_t>(3) * sequence_count) {
        throw std::runtime_error("MakeRList input dimensions differ");
    }
    auto& rcorr = decision.correlations.correlation;
    const auto& rinv = decision.correlations.inversion;
    RdpCandidateLists result;
    result.list_by_threshold.assign(
        static_cast<std::size_t>(30) * sequence_count, 0);
    result.inverse_by_threshold.assign(
        static_cast<std::size_t>(30) * sequence_count, 0);
    result.list.assign(static_cast<std::size_t>(3) * sequence_count, 0);
    result.inverse.assign(static_cast<std::size_t>(3) * sequence_count, 0);
    result.probability_scores.assign(
        static_cast<std::size_t>(9) * sequence_count, 0.0);
    result.probability_values.assign(
        static_cast<std::size_t>(9) * sequence_count, 0.0);
    result.t_values.assign(static_cast<std::size_t>(9) * sequence_count, 0.0);
    result.totals.assign(static_cast<std::size_t>(6) * sequence_count, 0.0);
    result.list_scores.assign(static_cast<std::size_t>(3) * sequence_count, 0.0);

    int warning_total = 0;
    int unflagged_count = 0;
    for (const auto warning : decision.warnings) {
        warning_total += warning;
        if (warning == 0) ++unflagged_count;
    }
    const double target = warning_total > 0
        ? 0.9 - (0.9 / 3.0) * warning_total : 0.9;
    const auto rc_index = [](const int role, const int region,
                             const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto row_index = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };
    const auto score_index = [sequence_count](
                                 const int role, const int sequence,
                                 const int region) {
        return static_cast<std::size_t>(role) + sequence * 3 +
            static_cast<std::size_t>(region) * 3 * sequence_count;
    };
    const auto append_threshold = [&](const int role, const int threshold,
                                      const int sequence) {
        const int count_index = role + threshold * 3;
        const auto output_index = static_cast<std::size_t>(role) +
            result.last_by_threshold[count_index] * 3 +
            static_cast<std::size_t>(threshold) * 3 * sequence_count;
        result.list_by_threshold[output_index] = sequence;
        ++result.last_by_threshold[count_index];
    };
    const auto append_main = [&](const int role, const int sequence,
                                 const bool inverse) {
        const auto output_index = static_cast<std::size_t>(role) +
            result.last[role] * 3;
        result.list[output_index] = sequence;
        if (inverse) result.inverse[output_index] = 1;
        ++result.last[role];
        result.list_scores[row_index(role, sequence)] =
            result.totals[row_index(role, sequence)];
    };
    const auto append_inverse_if_supported = [&](const int role,
                                                  const int sequence) {
        for (int region = 0; region < 3; ++region) {
            const auto index = rc_index(role, region, sequence);
            if (rinv[index] > 0.0F && rcorr[index] > 0.83F &&
                decision.warnings[region] == 0) {
                append_main(role, sequence, true);
                return true;
            }
        }
        return false;
    };

    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const auto row = row_index(role, sequence);
            if (dont_redo[row] == 0) {
                int correlation_count = 0;
                for (int region = 0; region < 2; ++region) {
                    const float value = rcorr[rc_index(role, region, sequence)];
                    if (value > 0.95F && decision.warnings[region] == 0) {
                        correlation_count += value > 0.98F ? 2 : 1;
                    }
                }
                if (acceptable_correlations[row] == 1 ||
                    sequences[role] == sequence || correlation_count == 2) {
                    for (int region = 0; region < 3; ++region) {
                        const auto correlation = rc_index(role, region, sequence);
                        if (rinv[correlation] == 0.0F &&
                            decision.warnings[region] == 0) {
                            if (rcorr[correlation] >= 1.0F) {
                                rcorr[correlation] = 0.99999999F;
                            }
                            const auto score = score_index(role, sequence, region);
                            if (rcorr[correlation] <= 0.0F) {
                                result.t_values[score] = 0.0;
                            } else {
                                const double squared =
                                    rcorr[correlation] * rcorr[correlation];
                                result.t_values[score] = 1.0 - squared > 0.0
                                    ? rcorr[correlation] *
                                        std::sqrt(4.0 / (1.0 - squared))
                                    : 0.0;
                            }
                            const double probability = rdp_ttest_probability(
                                result.t_values[score], 4);
                            result.probability_values[score] =
                                (1.0 - probability) * 2.0;
                            result.probability_scores[score] =
                                1.0 - result.probability_values[score];
                            const double normalized =
                                result.probability_scores[score] >= 0.95
                                ? (result.probability_scores[score] - 0.95) /
                                    0.05
                                : (result.probability_scores[score] - 0.95) /
                                    0.95;
                            if (normalized > 0.0) {
                                result.totals[row] += normalized;
                            } else {
                                result.totals[row + 3 * sequence_count] =
                                    result.totals[row] - normalized;
                            }
                        }
                    }
                    double threshold = 0.0;
                    if (result.totals[row] > target) {
                        if (good_comparisons[sequence] == 1 ||
                            good_comparisons[sequence + sequence_count] == 1 ||
                            sequence == sequences[role]) {
                            append_main(role, sequence, false);
                            for (int level = 0; level < 10; ++level) {
                                if (result.totals[row] > threshold) {
                                    append_threshold(role, level, sequence);
                                } else {
                                    break;
                                }
                                threshold +=
                                    (unflagged_count - threshold) / 3.0;
                            }
                        }
                    } else {
                        for (int level = 0; level < 10; ++level) {
                            if (result.totals[row] > threshold ||
                                sequence == sequences[role]) {
                                append_threshold(role, level, sequence);
                            } else {
                                break;
                            }
                            threshold += (unflagged_count - threshold) / 3.0;
                        }
                        if (sequence == sequences[role]) {
                            append_main(role, sequence, false);
                            // The source appends the representative to all ten
                            // threshold lists a second time in this branch.
                            for (int level = 0; level < 10; ++level) {
                                append_threshold(role, level, sequence);
                            }
                        } else {
                            append_inverse_if_supported(role, sequence);
                        }
                    }
                } else {
                    append_inverse_if_supported(role, sequence);
                }
            } else if (sequence == sequences[role]) {
                append_main(role, sequence, false);
                for (int level = 0; level < 10; ++level) {
                    append_threshold(role, level, sequence);
                }
            }
        }
        --result.last[role];
        for (int level = 0; level < 10; ++level) {
            --result.last_by_threshold[role + level * 3];
        }
    }
    result.result = 1.0;
    return result;
}

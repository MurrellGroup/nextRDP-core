#include "identification_state.hpp"

#include <algorithm>
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

RdpActualEventResolution resolve_rdp_actual_events(
    const int sequence_length, const int next_no,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 6>& starts,
    const std::array<int, 6>& ends,
    RdpCorrelationDecisionState correlations,
    RdpCandidateLists candidates,
    const std::vector<unsigned char>& dont_redo,
    const RdpRawEventState& events,
    int permanent_next_no) {
    const int count = next_no + 1;
    if (sequence_length < 1 || next_no < 2 ||
        candidates.list.size() != static_cast<std::size_t>(3 * count) ||
        candidates.inverse.size() != static_cast<std::size_t>(3 * count) ||
        candidates.probability_values.size() !=
            static_cast<std::size_t>(9 * count) ||
        candidates.list_scores.size() != static_cast<std::size_t>(3 * count) ||
        correlations.correlations.correlation.size() !=
            static_cast<std::size_t>(9 * count) ||
        correlations.correlations.inversion.size() !=
            static_cast<std::size_t>(9 * count) ||
        dont_redo.size() != static_cast<std::size_t>(3 * count) ||
        events.current_xover.size() != static_cast<std::size_t>(count) ||
        events.xover_list.size() != static_cast<std::size_t>(count)) {
        throw std::runtime_error("FindActualEvents input dimensions differ");
    }
    if (permanent_next_no < 0) permanent_next_no = next_no;

    const auto row = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };
    const auto rc = [](const int role, const int region,
                       const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto ok = [](const int role, const int category,
                       const int sequence) {
        return static_cast<std::size_t>(role) + category * 3 + sequence * 57;
    };
    const auto probability = [count](const int role, const int sequence,
                                     const int region) {
        return static_cast<std::size_t>(role) + sequence * 3 +
            static_cast<std::size_t>(region) * 3 * count;
    };

    RdpActualEventResolution output;
    output.candidates = std::move(candidates);
    output.correlations = std::move(correlations);
    output.acceptable_sequences.assign(
        static_cast<std::size_t>(57) * count, 0.0);
    output.unfound.assign(static_cast<std::size_t>(3) * count, 0);
    output.breakpoint_matches.assign(static_cast<std::size_t>(6) * count, 0);
    output.best_matches.assign(static_cast<std::size_t>(3) * count, 0.0F);
    auto& rlist = output.candidates.list;
    auto& invlist = output.candidates.inverse;
    auto& rnum = output.candidates.last;
    auto& rcorr = output.correlations.correlations.correlation;
    const auto& rinv = output.correlations.correlations.inversion;

    // Module3.Check6: if direct and inverted evidence contradict each other,
    // RDP keeps the direct interpretation and zeroes the inverted regions.
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= rnum[role]; ++slot) {
            const auto list_index = row(role, slot);
            const int sequence = rlist[list_index];
            if (invlist[list_index] == 1) {
                for (int region = 0; region < 3; ++region) {
                    if (rcorr[rc(role, region, sequence)] > 0.83F &&
                        rinv[rc(role, region, sequence)] == 0.0F) {
                        invlist[list_index] = 0;
                        for (int other = 0; other < 3; ++other) {
                            if (rinv[rc(role, other, sequence)] > 0.0F) {
                                rcorr[rc(role, other, sequence)] = 0.0F;
                            }
                        }
                    }
                }
            }
        }
    }

    // Module3 category 4 followed by Module2.AddOK1. The legacy category
    // count is nineteen here even though FindActualEvents itself later uses
    // its older hard-coded eighteen-category stride when writing category 1.
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= rnum[role]; ++slot) {
            output.acceptable_sequences[
                ok(role, 4, rlist[row(role, slot)])] = 1.0;
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (dont_redo[row(role, sequence)] == 0) {
                for (int region = 0; region < 3; ++region) {
                    const double value = output.candidates.probability_values[
                        probability(role, sequence, region)];
                    const auto destination = ok(role, 0, sequence);
                    if (value > 0.0 &&
                        (value * permanent_next_no <
                             output.acceptable_sequences[destination] ||
                         output.acceptable_sequences[destination] == 0.0)) {
                        output.acceptable_sequences[destination] =
                            value * permanent_next_no;
                    }
                }
            } else {
                output.acceptable_sequences[ok(role, 0, sequence)] = 0.049;
            }
        }
    }

    std::vector<unsigned char> inversion_state(
        static_cast<std::size_t>(3) * count, 0);
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= rnum[role]; ++slot) {
            if (invlist[row(role, slot)] == 1) {
                inversion_state[row(role, rlist[row(role, slot)])] = 1;
            }
        }
    }

    const auto make_overlap = [&](const int beginning, const int ending,
                                  const int size_index) {
        std::vector<int> mask(sequence_length + 1, 0);
        if (beginning < ending) {
            output.region_sizes[size_index] = ending - beginning + 1;
            for (int position = beginning; position <= ending; ++position) {
                mask[position] = 1;
            }
        } else {
            output.region_sizes[size_index] =
                ending + sequence_length - beginning + 1;
            for (int position = 1; position <= ending; ++position) {
                mask[position] = 1;
            }
            for (int position = beginning; position <= sequence_length;
                 ++position) {
                mask[position] = 1;
            }
        }
        return mask;
    };
    const auto overlap_begin = make_overlap(starts[0], ends[1], 2);
    const auto overlap_end = make_overlap(starts[2], ends[3], 4);
    const auto overlap_event = make_overlap(starts[4], ends[4], 0);
    output.beginning_overlap_mask = overlap_begin;
    output.ending_overlap_mask = overlap_end;
    output.event_overlap_mask = overlap_event;

    std::vector<unsigned char> don(static_cast<std::size_t>(3) * count, 0);
    std::array<int, 2> candidate_scratch{};
    std::array<int, 3> trace_sequences{};
    std::array<double, 2> match{};
    std::array<int, 4> sequence_scratch{};
    std::array<unsigned char, 6> tried_permutations{};

    const auto overlap_size = [&](const RdpRawEvent& event,
                                  const std::vector<int>& mask,
                                  const int size_index) {
        int overlap = 0;
        if (event.beginning < event.ending) {
            output.region_sizes[size_index] =
                event.ending - event.beginning + 1;
            for (int position = event.beginning; position <= event.ending;
                 ++position) {
                overlap += mask[position];
            }
        } else {
            output.region_sizes[size_index] =
                event.ending + sequence_length - event.beginning + 1;
            for (int position = 1; position <= event.ending; ++position) {
                overlap += mask[position];
            }
            for (int position = event.beginning; position <= sequence_length;
                 ++position) {
                overlap += mask[position];
            }
        }
        return overlap;
    };
    const auto interval_overlap = [&](const int beginning, const int ending,
                                      const std::vector<int>& mask,
                                      const int size_index) {
        int overlap = 0;
        if (beginning < ending) {
            output.region_sizes[size_index] = ending - beginning + 1;
            for (int position = beginning; position <= ending; ++position) {
                overlap += mask[position];
            }
        } else {
            output.region_sizes[size_index] =
                ending + sequence_length - beginning + 1;
            for (int position = 1; position <= ending; ++position) {
                overlap += mask[position];
            }
            for (int position = beginning; position <= sequence_length;
                 ++position) {
                overlap += mask[position];
            }
        }
        return overlap;
    };

    for (int role = 0; role < 3; ++role) {
        auto& call = output.calls[role];
        call.role = role;
        std::vector<unsigned char> role_membership(
            static_cast<std::size_t>(3) * count, 0);
        role_membership[row(comparison_matrix[role],
                            sequences[comparison_matrix[role]])] = 1;
        role_membership[row(comparison_matrix[role + 3],
                            sequences[comparison_matrix[role + 3]])] = 1;
        for (int slot = 0; slot <= rnum[role]; ++slot) {
            role_membership[row(role, rlist[row(role, slot)])] = 1;
        }
        std::vector<int> found(count, 0);

        call.region_sizes_before = output.region_sizes;
        call.breakpoint_matches_before = output.breakpoint_matches;
        call.best_matches_before = output.best_matches;
        call.acceptable_sequences_before = output.acceptable_sequences;
        call.found_before = found;
        call.candidate_last_before = rnum;
        call.candidate_list_before = rlist;
        call.inversion_state_before = inversion_state;
        call.candidate_scratch_before = candidate_scratch;
        call.trace_sequences_before = trace_sequences;
        call.match_before = match;
        call.sequence_scratch_before = sequence_scratch;
        call.tried_permutations_before = tried_permutations;
        call.role_membership_before = role_membership;

        int old_slot = -1;
        int repeat_count = 0;
        for (int daughter = 0; daughter <= next_no; ++daughter) {
            if (role_membership[row(role, daughter)] == 1 ||
                sequences[comparison_matrix[role]] == daughter ||
                sequences[comparison_matrix[role + 3]] == daughter) {
                old_slot = -1;
                for (int slot = 1;
                     slot <= events.current_xover[daughter]; ++slot) {
                    if (old_slot != slot) {
                        tried_permutations.fill(0);
                        old_slot = slot;
                        repeat_count = 0;
                    } else {
                        ++repeat_count;
                        if (repeat_count > 6) {
                            tried_permutations.fill(0);
                            old_slot = slot;
                            ++slot;
                            repeat_count = 0;
                            if (slot > events.current_xover[daughter]) break;
                        }
                        if (slot == old_slot) {
                            int tried = 0;
                            for (const auto value : tried_permutations) {
                                tried += value;
                            }
                            if (tried == 6) {
                                ++slot;
                                tried_permutations.fill(0);
                                old_slot = slot;
                                repeat_count = 0;
                                if (slot > events.current_xover[daughter]) {
                                    break;
                                }
                            }
                        }
                    }
                    const auto& event = events.xover_list[daughter][slot - 1];
                    sequence_scratch[1] = event.major_parent;
                    if (role_membership[row(role, sequence_scratch[1])] == 1 ||
                        sequences[comparison_matrix[role]] ==
                            sequence_scratch[1] ||
                        sequences[comparison_matrix[role + 3]] ==
                            sequence_scratch[1]) {
                        int go_on = 0;
                        sequence_scratch[2] = event.minor_parent;
                        if (role_membership[row(role, sequence_scratch[2])] == 1 ||
                            sequences[comparison_matrix[role]] ==
                                sequence_scratch[2] ||
                            sequences[comparison_matrix[role + 3]] ==
                                sequence_scratch[2]) {
                            sequence_scratch[0] = daughter;
                            match[0] = 0.0;
                            const int s0 = sequence_scratch[0];
                            const int s1 = sequence_scratch[1];
                            const int s2 = sequence_scratch[2];
                            const std::array<std::array<int, 3>, 6> orders{{
                                {{s0, s1, s2}}, {{s0, s2, s1}},
                                {{s1, s2, s0}}, {{s1, s0, s2}},
                                {{s2, s1, s0}}, {{s2, s0, s1}},
                            }};
                            for (int permutation = 0; permutation < 6;
                                 ++permutation) {
                                const auto& order = orders[permutation];
                                if (tried_permutations[permutation] == 0 &&
                                    ((role == 0 && don[row(0, order[0])] == 0) ||
                                     (role == 1 && don[row(1, order[1])] == 0) ||
                                     (role == 2 && don[row(2, order[2])] == 0)) &&
                                    role_membership[row(0, order[0])] == 1 &&
                                    role_membership[row(1, order[1])] == 1 &&
                                    role_membership[row(2, order[2])] == 1) {
                                    tried_permutations[permutation] = 1;
                                    trace_sequences = order;
                                    match[0] = 3.0;
                                    break;
                                }
                            }
                            if (match[0] == 3.0) {
                                if (inversion_state[row(
                                        comparison_matrix[role],
                                        trace_sequences[
                                            comparison_matrix[role]])] == 0 &&
                                    inversion_state[row(
                                        comparison_matrix[role + 3],
                                        trace_sequences[
                                            comparison_matrix[role + 3]])] == 0) {
                                    int candidate = 0;
                                    for (; candidate <= rnum[role]; ++candidate) {
                                        if (rlist[row(role, candidate)] ==
                                            trace_sequences[role]) {
                                            break;
                                        }
                                    }
                                    if (candidate > rnum[role]) {
                                        match[0] = 0.0;
                                    } else {
                                        candidate_scratch[1] = candidate;
                                        go_on = 1;
                                    }
                                } else {
                                    match[0] = 0.0;
                                }
                                if (match[0] == 3.0 && go_on == 1) {
                                    const int overlap = overlap_size(
                                        event, overlap_event, 1);
                                    if (overlap > 0) {
                                        match[1] =
                                            (static_cast<double>(overlap) * 2.0) /
                                            (output.region_sizes[0] +
                                             output.region_sizes[1]);
                                    } else {
                                        match[1] = 0.0;
                                    }
                                    const double original_match = match[1];
                                    if (match[0] * match[1] > 1.0) {
                                        const int candidate_sequence = rlist[row(
                                            role, candidate_scratch[1])];
                                        if (rcorr[rc(role, 2,
                                                    candidate_sequence)] > 0.83F &&
                                            match[1] > 0.6) {
                                            match[0] = 1.0;
                                        }
                                        // Literal threshold.CPP indexing is
                                        // preserved in the second term: it
                                        // omits the normal three-column stride.
                                        const auto legacy_slot =
                                            static_cast<std::size_t>(role) +
                                            candidate_scratch[1];
                                        const int legacy_sequence =
                                            legacy_slot < rlist.size()
                                                ? rlist[legacy_slot] : 0;
                                        if (rcorr[rc(role, 2,
                                                    candidate_sequence)] > 0.83F ||
                                            rcorr[rc(role, 0,
                                                    legacy_sequence)] > 0.83F) {
                                            const int flank_overlap =
                                                interval_overlap(
                                                    starts[0], ends[1],
                                                    overlap_begin, 3);
                                            if (flank_overlap > 0) {
                                                match[1] =
                                                    (static_cast<double>(
                                                         flank_overlap) * 2.0) /
                                                    (output.region_sizes[2] +
                                                     output.region_sizes[3]);
                                            } else {
                                                match[1] = 0.0;
                                            }
                                            if (match[1] > 0.2) {
                                                match[0] += 1.0;
                                            } else if (rcorr[rc(
                                                           role, 0,
                                                           candidate_sequence)] >
                                                       0.83F) {
                                                if (match[1] == 0.0 ||
                                                    flank_overlap ==
                                                        output.region_sizes[2]) {
                                                    match[0] -= 0.5;
                                                }
                                            }
                                        }
                                        if (rcorr[rc(role, 2,
                                                    candidate_sequence)] > 0.83F ||
                                            rcorr[rc(role, 1,
                                                    candidate_sequence)] > 0.83F) {
                                            const int flank_overlap =
                                                interval_overlap(
                                                    starts[2], ends[3],
                                                    overlap_end, 5);
                                            if (flank_overlap > 0) {
                                                match[1] =
                                                    (static_cast<double>(
                                                         flank_overlap) * 2.0) /
                                                    (output.region_sizes[4] +
                                                     output.region_sizes[5]);
                                            } else {
                                                match[1] = 0.0;
                                            }
                                            if (match[1] > 0.2) {
                                                match[0] += 1.0;
                                            } else if (rcorr[rc(
                                                           role, 0,
                                                           candidate_sequence)] >
                                                       0.83F) {
                                                if (match[1] == 0.0 ||
                                                    flank_overlap ==
                                                        output.region_sizes[4]) {
                                                    match[0] -= 0.5;
                                                }
                                            }
                                        }
                                        if (match[0] >= 1.0) {
                                            found[candidate_scratch[1]] = 1;
                                            const auto best =
                                                row(role, candidate_sequence);
                                            if (output.best_matches[best] <
                                                original_match) {
                                                // FindActualEvents was compiled
                                                // for OKSeq(2,17,N), so this
                                                // write intentionally uses 54,
                                                // not the caller's 57 stride.
                                                output.acceptable_sequences[
                                                    static_cast<std::size_t>(role) +
                                                    3 +
                                                    candidate_sequence * 54] =
                                                    original_match;
                                                output.best_matches[best] =
                                                    static_cast<float>(
                                                        original_match);
                                                output.breakpoint_matches[
                                                    static_cast<std::size_t>(role) +
                                                    candidate_sequence * 6] =
                                                    event.beginning;
                                                output.breakpoint_matches[
                                                    static_cast<std::size_t>(role) +
                                                    3 +
                                                    candidate_sequence * 6] =
                                                    event.ending;
                                            }
                                            --slot;
                                        }
                                    } else {
                                        --slot;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        call.result = 0;
        call.region_sizes_after = output.region_sizes;
        call.breakpoint_matches_after = output.breakpoint_matches;
        call.best_matches_after = output.best_matches;
        call.acceptable_sequences_after = output.acceptable_sequences;
        call.found_after = found;
        call.candidate_scratch_after = candidate_scratch;
        call.trace_sequences_after = trace_sequences;
        call.match_after = match;
        call.sequence_scratch_after = sequence_scratch;
        call.tried_permutations_after = tried_permutations;

        for (int slot = 0; slot <= rnum[role]; ++slot) {
            if (found[slot] == 0 || invlist[row(role, slot)] == 1) {
                output.unfound[row(role, rlist[row(role, slot)])] = 1;
            }
        }

        // threshold.CPP::StripUnfound, including its swap-with-last ordering.
        int slot = 0;
        while (slot <= rnum[role]) {
            if (found[slot] == 0) {
                const int sequence = rlist[row(role, slot)];
                if (invlist[row(role, slot)] == 1 ||
                    rcorr[rc(role, 0, sequence)] < 0.95F ||
                    rcorr[rc(role, 1, sequence)] < 0.95F ||
                    rcorr[rc(role, 2, sequence)] < 0.95F) {
                    if (slot < rnum[role]) {
                        output.candidates.list_scores[
                            row(role, sequence)] = 0.0;
                        rlist[row(role, slot)] =
                            rlist[row(role, rnum[role])];
                        invlist[row(role, slot)] =
                            invlist[row(role, rnum[role])];
                        found[slot] = found[rnum[role]];
                    }
                    --rnum[role];
                } else {
                    ++slot;
                }
            } else {
                ++slot;
            }
        }
    }
    output.candidate_last_before_strip = rnum;
    output.candidate_list_before_strip = rlist;
    output.candidate_inverse_before_strip = invlist;

    // threshold.CPP::StripDupInv. Its duplicate-removal block is dead because
    // MScore is initialized to zero and tested against ten, so only the
    // inversion penalty and swap-with-last removal execute.
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= rnum[role]; ++slot) {
            if (invlist[row(role, slot)] > 0) {
                output.inversion_penalty[role] = 1;
                break;
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        int slot = 0;
        while (slot <= rnum[role]) {
            if (invlist[row(role, slot)] == 1) {
                if (slot < rnum[role]) {
                    rlist[row(role, slot)] = rlist[row(role, rnum[role])];
                    invlist[row(role, slot)] =
                        invlist[row(role, rnum[role])];
                }
                --rnum[role];
            } else {
                ++slot;
            }
        }
    }
    return output;
}

RdpEventSetState find_rdp_event_sets(
    const int sequence_length, const int next_no, const int beginning,
    const int ending, const std::array<int, 3>& sequences,
    const RdpRawEventState& events) {
    const int count = next_no + 1;
    if (sequence_length < 1 || beginning < 1 || ending < 1 ||
        beginning > sequence_length || ending > sequence_length ||
        events.current_xover.size() != static_cast<std::size_t>(count) ||
        events.xover_list.size() != static_cast<std::size_t>(count)) {
        throw std::runtime_error("FindSets input dimensions differ");
    }
    const auto cell = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };
    const int target_size = ending > beginning
        ? ending - beginning + 1
        : ending + sequence_length - beginning + 1;
    std::vector<unsigned char> target_mask(
        static_cast<std::size_t>(sequence_length + 1), 0);
    int target_position = beginning;
    while (target_position != ending) {
        target_mask[target_position] = 1;
        ++target_position;
        if (target_position > sequence_length) target_position = 1;
    }
    target_mask[ending] = 1;
    // DNA!DoSetsAP computes the initial overlap by walking the event through
    // MakeOLSeqB's target mask. Preserve its begin==end full-circle behavior.
    const auto initial_event_overlap = [&](const RdpRawEvent& event) {
        int overlap = 0;
        int event_size = 0;
        if (event.beginning < event.ending) {
            event_size = event.ending - event.beginning + 1;
            for (int position = event.beginning;
                 position <= event.ending; ++position) {
                overlap += target_mask[position];
            }
        } else {
            event_size = event.ending + sequence_length -
                event.beginning + 1;
            for (int position = event.beginning;
                 position <= sequence_length; ++position) {
                overlap += target_mask[position];
            }
            for (int position = 1; position <= event.ending; ++position) {
                overlap += target_mask[position];
            }
        }
        return static_cast<float>(overlap) /
            ((static_cast<float>(target_size) + event_size) / 2.0F) > 0.3F;
    };
    // DNA5!FillSetsP3 replaces the mask walk during closure with this
    // hand-expanded branch table. It is intentionally not simplified: its
    // wrapped-interval arithmetic is part of the reference execution path.
    const auto closure_event_overlap = [&](const RdpRawEvent& event) {
        int overlap = 0;
        int event_size = 0;
        if (beginning < ending) {
            if (event.beginning < event.ending) {
                event_size = event.ending - event.beginning + 1;
                if (event.ending >= beginning && event.beginning <= ending) {
                    if (ending >= event.ending &&
                        beginning <= event.beginning) {
                        overlap = event_size;
                    } else if (ending <= event.ending &&
                               beginning >= event.beginning) {
                        overlap = target_size;
                    } else if (ending <= event.ending) {
                        overlap = ending - event.beginning + 1;
                    } else if (beginning >= event.beginning) {
                        overlap = event.ending - beginning + 1;
                    }
                }
            } else {
                event_size = event.ending + sequence_length -
                    event.beginning + 1;
                if (event.ending >= beginning) {
                    overlap = ending >= event.ending
                        ? event.ending - beginning + 1
                        : target_size;
                }
                if (event.beginning <= ending) {
                    overlap = beginning <= event.beginning
                        ? ending - event.beginning + 1
                        : target_size;
                }
            }
        } else {
            if (event.beginning < event.ending) {
                event_size = event.ending - event.beginning + 1;
                if (ending >= event.beginning) {
                    overlap = ending >= event.ending
                        ? event_size : ending - event.beginning + 1;
                }
                if (beginning <= event.ending) {
                    overlap += beginning <= event.beginning
                        ? event_size : event.ending - beginning + 1;
                }
            } else if (event.beginning > event.ending) {
                event_size = event.ending + sequence_length -
                    event.beginning + 1;
                overlap = event.ending <= ending
                    ? event.ending : ending;
                overlap += event.beginning <= beginning
                    ? sequence_length - beginning + 1
                    : sequence_length - event.beginning + 1;
            }
        }
        return static_cast<float>(overlap) /
            ((static_cast<float>(target_size) + event_size) / 2.0F) > 0.3F;
    };
    const auto event_members = [](const RdpRawEvent& event) {
        return std::array<int, 3>{
            event.daughter, event.major_parent, event.minor_parent};
    };

    RdpEventSetState output;
    output.candidate_list.assign(static_cast<std::size_t>(3 * count), 0);
    std::vector<unsigned char> sets(static_cast<std::size_t>(3 * count), 0);
    for (int row = 0; row <= next_no; ++row) {
        for (const auto& event : events.xover_list[row]) {
            const auto members = event_members(event);
            std::array<unsigned char, 3> involved{};
            bool any = false;
            for (int role = 0; role < 3; ++role) {
                if (std::find(members.begin(), members.end(),
                              sequences[role]) != members.end()) {
                    involved[role] = 1;
                    any = true;
                }
            }
            if (any && initial_event_overlap(event)) {
                for (int role = 0; role < 3; ++role) {
                    if (involved[role] == 0) continue;
                    for (const int member : members) {
                        if (member >= 0 && member <= next_no) {
                            sets[cell(role, member)] = 1;
                        }
                    }
                }
            }
        }
    }

    std::vector<unsigned char> membership(
        static_cast<std::size_t>(3 * count), 0);
    for (int role = 0; role < 3; ++role) {
        membership[cell(role, sequences[role])] = 1;
    }
    bool changed = false;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        const int total = sets[cell(0, sequence)] +
            sets[cell(1, sequence)] + sets[cell(2, sequence)];
        if (total == 2) {
            for (int role = 0; role < 3; ++role) {
                if (sets[cell(role, sequence)] == 0 &&
                    membership[cell(role, sequence)] == 0) {
                    membership[cell(role, sequence)] = 1;
                    changed = true;
                }
            }
        }
    }
    const auto rebuild_lists = [&] {
        output.candidate_last = {-1, -1, -1};
        std::fill(output.candidate_list.begin(),
                  output.candidate_list.end(), 0);
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                if (membership[cell(role, sequence)] != 0) {
                    ++output.candidate_last[role];
                    output.candidate_list[cell(
                        role, output.candidate_last[role])] = sequence;
                }
            }
        }
    };
    rebuild_lists();

    std::vector<unsigned char> unique_sets(
        static_cast<std::size_t>(7 * count), 0);
    while (changed) {
        std::fill(sets.begin(), sets.end(), 0);
        for (int row = 0; row <= next_no; ++row) {
            for (const auto& event : events.xover_list[row]) {
                if (!closure_event_overlap(event)) continue;
                const auto members = event_members(event);
                for (int role = 0; role < 3; ++role) {
                    bool touches = false;
                    for (const int member : members) {
                        if (member >= 0 && member <= next_no &&
                            membership[cell(role, member)] != 0) {
                            touches = true;
                            break;
                        }
                    }
                    if (touches) {
                        for (const int member : members) {
                            if (member >= 0 && member <= next_no) {
                                sets[cell(role, member)] = 1;
                            }
                        }
                    }
                }
            }
        }
        changed = false;
        std::fill(unique_sets.begin(), unique_sets.end(), 0);
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const int total = sets[cell(0, sequence)] +
                sets[cell(1, sequence)] + sets[cell(2, sequence)];
            if (total == 2) {
                for (int role = 0; role < 3; ++role) {
                    if (sets[cell(role, sequence)] == 0 &&
                        membership[cell(role, sequence)] == 0) {
                        membership[cell(role, sequence)] = 1;
                        changed = true;
                    }
                }
            } else if (total == 1) {
                for (int role = 0; role < 3; ++role) {
                    if (sets[cell(role, sequence)] != 0) {
                        unique_sets[static_cast<std::size_t>(role + 4) +
                                    sequence * 7] = 1;
                    }
                }
            }
        }
        rebuild_lists();
    }
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            output.set_totals[role * 2] +=
                unique_sets[static_cast<std::size_t>(role + 4) +
                            sequence * 7];
        }
    }
    output.role_sets = std::move(sets);
    return output;
}

RdpTreeCompatibilityCallState make_rdp_tree_compatibility_call(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix, const int role,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& matrix,
    const std::array<double, 3>& list_distances,
    std::array<int, 3>& compatibility,
    std::array<int, 3>& reverse_compatibility,
    std::array<int, 3>& nonrecombinant_last,
    std::vector<int>& nonrecombinant_list) {
    const int count = next_no + 1;
    if (role < 0 || role > 2 ||
        candidate_list.size() != static_cast<std::size_t>(3 * count) ||
        good_comparisons.size() != static_cast<std::size_t>(2 * count) ||
        matrix.size() != static_cast<std::size_t>(count) * count ||
        nonrecombinant_list.size() != static_cast<std::size_t>(3 * count)) {
        throw std::runtime_error("MakeRCompat input dimensions differ");
    }
    const auto cell = [](const int selected_role, const int item) {
        return static_cast<std::size_t>(selected_role) + item * 3;
    };
    RdpTreeCompatibilityCallState call;
    call.role = role;
    call.compatibility_before = compatibility;
    call.reverse_compatibility_before = reverse_compatibility;
    call.nonrecombinant_last_before = nonrecombinant_last;
    call.nonrecombinant_list_before = nonrecombinant_list;
    call.list_distances = list_distances;
    std::vector<int> categories(next_no * 3 + 1, 0);
    std::vector<int> done(count, 0);
    call.done_before = done;

    const int first_selected = sequences[comparison_matrix[role]];
    const int second_selected = sequences[comparison_matrix[role + 3]];
    done[first_selected] = 1;
    done[second_selected] = 1;
    for (int slot = 0; slot <= candidate_last[role]; ++slot) {
        const int recombinant = candidate_list[cell(role, slot)];
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (done[sequence] == 0 &&
                (good_comparisons[sequence] == 1 ||
                 good_comparisons[sequence + count] == 1) &&
                matrix[recombinant +
                       static_cast<std::size_t>(sequence) * count] <
                    list_distances[role]) {
                int comparison = 0;
                for (; comparison <= candidate_last[role]; ++comparison) {
                    if (sequence == candidate_list[cell(role, comparison)]) {
                        break;
                    }
                }
                if (comparison == candidate_last[role] + 1) {
                    done[sequence] = 1;
                    nonrecombinant_list[cell(
                        role, nonrecombinant_last[role])] = sequence;
                    ++nonrecombinant_last[role];
                }
            }
        }
    }
    --nonrecombinant_last[role];

    compatibility[role] = 0;
    for (int slot = 0; slot <= candidate_last[role]; ++slot) {
        std::fill(categories.begin(), categories.end(), 0);
        const int recombinant = candidate_list[cell(role, slot)];
        if (nonrecombinant_last[role] > -1) {
            for (int item = 0; item <= nonrecombinant_last[role]; ++item) {
                const int nonrecombinant =
                    nonrecombinant_list[cell(role, item)];
                const int category = static_cast<int>(
                    matrix[recombinant +
                           static_cast<std::size_t>(nonrecombinant) * count] *
                        1000.0F + 0.0000001F);
                categories[category] = 1;
            }
        }
        for (const int selected : {first_selected, second_selected}) {
            if (matrix[recombinant +
                       static_cast<std::size_t>(selected) * count] <
                list_distances[role]) {
                for (int item = 0; item <= candidate_last[role]; ++item) {
                    const int other = candidate_list[cell(role, item)];
                    const int category = static_cast<int>(
                        matrix[other +
                               static_cast<std::size_t>(selected) * count] *
                            1000.0F + 0.0000001F);
                    categories[category] = 1;
                }
            }
        }
        int category_count = 0;
        for (int category = 0; category <= next_no; ++category) {
            category_count += categories[category];
        }
        compatibility[role] = std::max(compatibility[role], category_count);
    }

    int added_first = 0;
    int added_second = 0;
    for (int slot = 0; slot <= candidate_last[role]; ++slot) {
        const int recombinant = candidate_list[cell(role, slot)];
        if (matrix[recombinant +
                   static_cast<std::size_t>(first_selected) * count] <
                list_distances[role] && added_first == 0) {
            ++nonrecombinant_last[role];
            nonrecombinant_list[cell(role, nonrecombinant_last[role])] =
                first_selected;
            added_first = 1;
        }
        if (matrix[recombinant +
                   static_cast<std::size_t>(second_selected) * count] <
                list_distances[role] && added_second == 0) {
            ++nonrecombinant_last[role];
            nonrecombinant_list[cell(role, nonrecombinant_last[role])] =
                second_selected;
            added_second = 1;
        }
    }
    reverse_compatibility[role] = 0;
    if (nonrecombinant_last[role] > -1) {
        for (int item = 0; item <= nonrecombinant_last[role]; ++item) {
            std::fill(categories.begin(), categories.end(), 0);
            const int nonrecombinant =
                nonrecombinant_list[cell(role, item)];
            for (int slot = 0; slot <= candidate_last[role]; ++slot) {
                const int recombinant = candidate_list[cell(role, slot)];
                const int category = static_cast<int>(
                    matrix[nonrecombinant +
                           static_cast<std::size_t>(recombinant) * count] *
                        1000.0F + 0.0000001F);
                categories[category] = 1;
            }
            int category_count = 0;
            for (int category = 0; category <= next_no; ++category) {
                category_count += categories[category];
            }
            --category_count;
            reverse_compatibility[role] = std::max(
                reverse_compatibility[role], category_count);
        }
    }
    if (nonrecombinant_last[role] > -1 &&
        reverse_compatibility[role] < compatibility[role]) {
        compatibility[role] = reverse_compatibility[role];
    }
    if (compatibility[role] > candidate_last[role]) {
        compatibility[role] = candidate_last[role];
    }
    if (compatibility[role] > 0) {
        compatibility[role] += inversion_penalty[role];
    }
    call.compatibility_after = compatibility;
    call.reverse_compatibility_after = reverse_compatibility;
    call.nonrecombinant_last_after = nonrecombinant_last;
    return call;
}

RdpTreeCompatibilityState evaluate_rdp_tree_compatibility(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& background_ancestor_matrix,
    const std::vector<float>& region_ancestor_matrix) {
    const int count = next_no + 1;
    if (candidate_list.size() != static_cast<std::size_t>(3 * count) ||
        good_comparisons.size() != static_cast<std::size_t>(2 * count) ||
        background_ancestor_matrix.size() !=
            static_cast<std::size_t>(count) * count ||
        region_ancestor_matrix.size() !=
            static_cast<std::size_t>(count) * count) {
        throw std::runtime_error("MakeRCompat input dimensions differ");
    }
    const auto cell = [](const int role, const int item) {
        return static_cast<std::size_t>(role) + item * 3;
    };
    const auto make_list_distances = [&](const std::vector<float>& matrix) {
        std::array<double, 3> distances{};
        for (int role = 0; role < 3; ++role) {
            for (int first = 0; first < candidate_last[role]; ++first) {
                const int first_sequence =
                    candidate_list[cell(role, first)];
                const auto first_offset =
                    static_cast<std::size_t>(first_sequence) * count;
                for (int second = first + 1;
                     second <= candidate_last[role]; ++second) {
                    const float distance = matrix[
                        candidate_list[cell(role, second)] + first_offset];
                    if (distance > distances[role]) {
                        distances[role] = distance;
                    }
                }
            }
        }
        return distances;
    };

    RdpTreeCompatibilityState output;
    output.background_list_distances =
        make_list_distances(background_ancestor_matrix);
    output.region_list_distances =
        make_list_distances(region_ancestor_matrix);

    const auto run_call = [&](const int call_index, const int role,
                              const std::vector<float>& matrix,
                              const std::array<double, 3>& list_distances,
                              std::array<int, 3>& compatibility,
                              std::array<int, 3>& reverse_compatibility,
                              std::array<int, 3>& nonrecombinant_last,
                              std::vector<int>& nonrecombinant_list) {
        auto& call = output.calls[call_index];
        call.role = role;
        std::vector<int> categories(next_no * 3 + 1, 0);
        std::vector<int> done(count, 0);
        call.compatibility_before = compatibility;
        call.reverse_compatibility_before = reverse_compatibility;
        call.nonrecombinant_last_before = nonrecombinant_last;
        call.done_before = done;
        call.nonrecombinant_list_before = nonrecombinant_list;
        call.list_distances = list_distances;

        const int first_selected =
            sequences[comparison_matrix[role]];
        const int second_selected =
            sequences[comparison_matrix[role + 3]];
        done[first_selected] = 1;
        done[second_selected] = 1;
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            const int recombinant = candidate_list[cell(role, slot)];
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                if (done[sequence] == 0 &&
                    (good_comparisons[sequence] == 1 ||
                     good_comparisons[sequence + count] == 1) &&
                    matrix[recombinant +
                           static_cast<std::size_t>(sequence) * count] <
                        list_distances[role]) {
                    int comparison = 0;
                    for (; comparison <= candidate_last[role]; ++comparison) {
                        if (sequence == candidate_list[cell(role, comparison)]) {
                            break;
                        }
                    }
                    if (comparison == candidate_last[role] + 1) {
                        done[sequence] = 1;
                        nonrecombinant_list[cell(
                            role, nonrecombinant_last[role])] = sequence;
                        ++nonrecombinant_last[role];
                    }
                }
            }
        }
        --nonrecombinant_last[role];

        compatibility[role] = 0;
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            std::fill(categories.begin(), categories.end(), 0);
            const int recombinant = candidate_list[cell(role, slot)];
            if (nonrecombinant_last[role] > -1) {
                for (int item = 0; item <= nonrecombinant_last[role]; ++item) {
                    const int nonrecombinant =
                        nonrecombinant_list[cell(role, item)];
                    const int category = static_cast<int>(
                        matrix[recombinant +
                               static_cast<std::size_t>(nonrecombinant) *
                                   count] *
                            1000.0F +
                        0.0000001F);
                    categories[category] = 1;
                }
            }
            if (matrix[recombinant +
                       static_cast<std::size_t>(first_selected) * count] <
                list_distances[role]) {
                for (int item = 0; item <= candidate_last[role]; ++item) {
                    const int other = candidate_list[cell(role, item)];
                    const int category = static_cast<int>(
                        matrix[other +
                               static_cast<std::size_t>(first_selected) *
                                   count] *
                            1000.0F +
                        0.0000001F);
                    categories[category] = 1;
                }
            }
            if (matrix[recombinant +
                       static_cast<std::size_t>(second_selected) * count] <
                list_distances[role]) {
                for (int item = 0; item <= candidate_last[role]; ++item) {
                    const int other = candidate_list[cell(role, item)];
                    const int category = static_cast<int>(
                        matrix[other +
                               static_cast<std::size_t>(second_selected) *
                                   count] *
                            1000.0F +
                        0.0000001F);
                    categories[category] = 1;
                }
            }
            int category_count = 0;
            for (int category = 0; category <= next_no; ++category) {
                category_count += categories[category];
            }
            if (category_count > compatibility[role]) {
                compatibility[role] = category_count;
            }
        }

        int added_first = 0;
        int added_second = 0;
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            const int recombinant = candidate_list[cell(role, slot)];
            if (matrix[recombinant +
                       static_cast<std::size_t>(first_selected) * count] <
                    list_distances[role] &&
                added_first == 0) {
                ++nonrecombinant_last[role];
                nonrecombinant_list[cell(
                    role, nonrecombinant_last[role])] = first_selected;
                added_first = 1;
            }
            if (matrix[recombinant +
                       static_cast<std::size_t>(second_selected) * count] <
                    list_distances[role] &&
                added_second == 0) {
                ++nonrecombinant_last[role];
                nonrecombinant_list[cell(
                    role, nonrecombinant_last[role])] = second_selected;
                added_second = 1;
            }
        }
        reverse_compatibility[role] = 0;
        if (nonrecombinant_last[role] > -1) {
            for (int item = 0; item <= nonrecombinant_last[role]; ++item) {
                std::fill(categories.begin(), categories.end(), 0);
                const int nonrecombinant =
                    nonrecombinant_list[cell(role, item)];
                for (int slot = 0; slot <= candidate_last[role]; ++slot) {
                    const int recombinant = candidate_list[cell(role, slot)];
                    const int category = static_cast<int>(
                        matrix[nonrecombinant +
                               static_cast<std::size_t>(recombinant) * count] *
                            1000.0F +
                        0.0000001F);
                    categories[category] = 1;
                }
                int category_count = 0;
                for (int category = 0; category <= next_no; ++category) {
                    category_count += categories[category];
                }
                --category_count;
                if (category_count > reverse_compatibility[role]) {
                    reverse_compatibility[role] = category_count;
                }
            }
        }
        if (nonrecombinant_last[role] > -1 &&
            reverse_compatibility[role] < compatibility[role]) {
            compatibility[role] = reverse_compatibility[role];
        }
        if (compatibility[role] > candidate_last[role]) {
            compatibility[role] = candidate_last[role];
        }
        if (compatibility[role] > 0) {
            compatibility[role] += inversion_penalty[role];
        }

        call.compatibility_after = compatibility;
        call.reverse_compatibility_after = reverse_compatibility;
        call.nonrecombinant_last_after = nonrecombinant_last;
    };

    std::array<int, 3> background_last{};
    std::vector<int> background_list(static_cast<std::size_t>(3) * count, 0);
    for (int role = 0; role < 3; ++role) {
        run_call(role, role, background_ancestor_matrix,
                 output.background_list_distances,
                 output.background_compatibility,
                 output.background_reverse_compatibility,
                 background_last, background_list);
    }
    for (int role = 0; role < 3; ++role) {
        // The VB caller ReDims NRNum2 and NRList2 inside this loop, so every
        // region-matrix role starts with a completely fresh list family.
        std::array<int, 3> region_last{};
        std::vector<int> region_list(static_cast<std::size_t>(3) * count, 0);
        run_call(3 + role, role, region_ancestor_matrix,
                 output.region_list_distances,
                 output.region_compatibility,
                 output.region_reverse_compatibility,
                 region_last, region_list);
    }
    return output;
}

RdpTreeCompatibilityFlowState run_rdp_tree_compatibility_flow(
    const int sequence_length, const int next_no, const int beginning,
    const int ending, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& inversion_penalty,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<int>& good_comparisons,
    const std::vector<float>& background_ancestor_matrix,
    const std::vector<float>& region_ancestor_matrix,
    const std::vector<float>& background_secondary_matrix,
    const std::vector<float>& region_secondary_matrix,
    const RdpRawEventState& events) {
    const int count = next_no + 1;
    const auto matrix_cells = static_cast<std::size_t>(count) * count;
    if (background_ancestor_matrix.size() != matrix_cells ||
        region_ancestor_matrix.size() != matrix_cells ||
        (!background_secondary_matrix.empty() &&
         background_secondary_matrix.size() != matrix_cells) ||
        (!region_secondary_matrix.empty() &&
         region_secondary_matrix.size() != matrix_cells)) {
        throw std::runtime_error("tree-compatibility flow matrix dimensions differ");
    }
    const auto cell = [](const int role, const int item) {
        return static_cast<std::size_t>(role) + item * 3;
    };
    const auto make_list_distances = [&](const std::array<int, 3>& last,
                                         const std::vector<int>& list,
                                         const std::vector<float>& matrix) {
        std::array<double, 3> distances{};
        for (int role = 0; role < 3; ++role) {
            for (int first = 0; first < last[role]; ++first) {
                const int first_sequence = list[cell(role, first)];
                for (int second = first + 1; second <= last[role]; ++second) {
                    distances[role] = std::max<double>(
                        distances[role],
                        matrix[list[cell(role, second)] +
                               static_cast<std::size_t>(first_sequence) *
                                   count]);
                }
            }
        }
        return distances;
    };
    const auto tied = [](const std::array<int, 3>& values) {
        return values[0] == values[1] && values[0] == values[2];
    };

    RdpTreeCompatibilityFlowState output;
    output.event_sets = find_rdp_event_sets(
        sequence_length, next_no, beginning, ending, sequences, events);

    const auto run_family = [&](const std::array<int, 3>& last,
                                const std::vector<int>& list,
                                const std::vector<float>& matrix,
                                std::array<int, 3>& compatibility,
                                const bool preserve_between_roles) {
        std::array<int, 3> reverse{};
        const auto distances = make_list_distances(last, list, matrix);
        std::array<int, 3> nonrecombinant_last{};
        std::vector<int> nonrecombinant_list(
            static_cast<std::size_t>(3) * count, 0);
        for (int role = 0; role < 3; ++role) {
            if (!preserve_between_roles) {
                nonrecombinant_last = {};
                nonrecombinant_list.assign(
                    static_cast<std::size_t>(3) * count, 0);
            }
            output.calls.push_back(make_rdp_tree_compatibility_call(
                next_no, sequences, comparison_matrix, role,
                inversion_penalty, last, list, good_comparisons, matrix,
                distances, compatibility, reverse, nonrecombinant_last,
                nonrecombinant_list));
        }
    };

    // Module3.ProcessEvent: FAMat and SAMat are always tried first. The
    // bootstrap/alternate-set families are entered only while each preceding
    // family remains tied.
    run_family(candidate_last, candidate_list, background_ancestor_matrix,
               output.background, true);
    run_family(candidate_last, candidate_list, region_ancestor_matrix,
               output.region, false);
    if (tied(output.background) && next_no > 2 &&
        !background_secondary_matrix.empty()) {
        run_family(candidate_last, candidate_list, background_secondary_matrix,
                   output.background_secondary, false);
        if (tied(output.background_secondary)) {
            run_family(output.event_sets.candidate_last,
                       output.event_sets.candidate_list,
                       background_ancestor_matrix, output.background_sets,
                       false);
            if (tied(output.background_sets) && next_no > 2) {
                run_family(output.event_sets.candidate_last,
                           output.event_sets.candidate_list,
                           background_secondary_matrix,
                           output.background_secondary_sets, false);
            }
        }
    }
    if (tied(output.region) && next_no > 2 &&
        !region_secondary_matrix.empty()) {
        run_family(candidate_last, candidate_list, region_secondary_matrix,
                   output.region_secondary, false);
        if (tied(output.region_secondary)) {
            run_family(output.event_sets.candidate_last,
                       output.event_sets.candidate_list,
                       region_ancestor_matrix, output.region_sets, false);
            if (tied(output.region_sets) && next_no > 2) {
                run_family(output.event_sets.candidate_last,
                           output.event_sets.candidate_list,
                           region_secondary_matrix,
                           output.region_secondary_sets, false);
            }
        }
    }
    return output;
}

RdpPatternState check_rdp_sequence_patterns(
    const int sequence_length, const int next_no,
    const std::array<int, 3>& sequences,
    const std::array<int, 3>& starts,
    const std::array<int, 3>& ends,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<short>& sequence_data,
    const std::vector<double>& acceptable_sequences) {
    const int count = next_no + 1;
    const int positions = sequence_length + 1;
    if (sequence_data.size() <
            static_cast<std::size_t>(positions) * count ||
        acceptable_sequences.size() != static_cast<std::size_t>(57) * count) {
        throw std::runtime_error("CheckPattern input dimensions differ");
    }
    const auto nucleotide = [&](const int position, const int sequence) {
        return sequence_data[position +
            static_cast<std::size_t>(sequence) * positions];
    };
    const auto pattern_cell = [](const int role, const int region,
                                 const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto done_cell = [](const int region, const int sequence) {
        return static_cast<std::size_t>(region) + sequence * 3;
    };
    const auto acceptable_cell = [](const int role, const int category,
                                    const int sequence) {
        return static_cast<std::size_t>(role) + category * 3 + sequence * 57;
    };

    RdpPatternState output;
    output.pattern.assign(static_cast<std::size_t>(9) * count, 0.0);
    output.done.assign(static_cast<std::size_t>(3) * count, 0);
    output.acceptable_sequences = acceptable_sequences;

    // Literal DNA!CheckPatternX traversal. G is mostly redundant because
    // DonePattern makes each (region, sequence) cell execute only once, but
    // retaining the loop also retains the selected-sequence skip order.
    const auto score_position = [&](const int position, const int target,
                                    const int region) {
        const short target_nucleotide = nucleotide(position, target);
        const short first = nucleotide(position, sequences[0]);
        const short second = nucleotide(position, sequences[1]);
        const short third = nucleotide(position, sequences[2]);
        if (target_nucleotide == 46 || first == 46 || second == 46 ||
            third == 46 || (first == second && first == third)) {
            return;
        }
        if (target_nucleotide == first) {
            if (target_nucleotide != second) {
                if (second == third) {
                    output.pattern[pattern_cell(0, region, target)] += 1.0;
                } else if (target_nucleotide != third) {
                    output.pattern[pattern_cell(0, region, target)] += 0.5;
                }
            }
        } else if (target_nucleotide == second) {
            if (first == third) {
                output.pattern[pattern_cell(1, region, target)] += 1.0;
            } else if (target_nucleotide != third) {
                output.pattern[pattern_cell(1, region, target)] += 0.5;
            }
        } else if (target_nucleotide == third) {
            if (first == second) {
                output.pattern[pattern_cell(2, region, target)] += 1.0;
            } else if (target_nucleotide != second) {
                output.pattern[pattern_cell(2, region, target)] += 0.5;
            }
        }
    };
    for (int selected_role = 0; selected_role < 3; ++selected_role) {
        for (int target = 0; target <= next_no; ++target) {
            if (target == sequences[selected_role]) continue;
            for (int region = 0; region < 3; ++region) {
                const auto done_index = done_cell(region, target);
                if (output.done[done_index] != 0) continue;
                output.done[done_index] = 1;
                if (starts[region] < ends[region]) {
                    for (int position = starts[region];
                         position <= ends[region]; ++position) {
                        score_position(position, target, region);
                    }
                } else {
                    for (int position = 1; position <= ends[region];
                         ++position) {
                        score_position(position, target, region);
                    }
                    for (int position = starts[region];
                         position <= sequence_length; ++position) {
                        score_position(position, target, region);
                    }
                }
            }
        }
    }

    // The active VB CheckPattern analysis block never removes a candidate;
    // its only surviving state change is this per-role proportion panel.
    for (int role = 0; role < 3; ++role) {
        for (int target = 0; target <= next_no; ++target) {
            const auto destination = acceptable_cell(role, 3, target);
            if (target == sequences[role]) {
                output.acceptable_sequences[destination] = 3.0;
                continue;
            }
            for (int region = 0; region < 3; ++region) {
                const double total =
                    output.pattern[pattern_cell(role, region, target)] +
                    output.pattern[pattern_cell(
                        comparison_matrix[role], region, target)] +
                    output.pattern[pattern_cell(
                        comparison_matrix[role + 3], region, target)];
                if (total > 0.0) {
                    output.acceptable_sequences[destination] +=
                        output.pattern[pattern_cell(role, region, target)] /
                        total;
                }
            }
        }
    }
    return output;
}

RdpPhylProScoreState make_rdp_phylpro_scores(
    const int next_no, const double minimum_offset,
    const std::vector<int>& done_this,
    const std::array<int, 3>& sequences,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix) {
    const int count = next_no + 1;
    const auto matrix_cells = static_cast<std::size_t>(count) * count;
    if (done_this.size() != static_cast<std::size_t>(2 * count) ||
        background_matrix.size() != matrix_cells ||
        region_matrix.size() != matrix_cells) {
        throw std::runtime_error("MakePhPrScore input dimensions differ");
    }

    // Literal port of threshold.CPP::MakePhPrScore. In particular,
    // NumInvolved is a last index, while NS is the number of observations
    // remaining after the current reference sequence is excluded.
    RdpPhylProScoreState output;
    output.trace_involved.assign(count, 0);
    for (int role = 0; role < 3; ++role) {
        output.trace_involved[role] = sequences[role];
    }
    int last_involved = 2;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (done_this[sequence * 2] == 0 ||
            done_this[sequence * 2 + 1] == 0) {
            ++last_involved;
            if (last_involved <= next_no) {
                output.trace_involved[last_involved] = sequence;
            } else {
                --last_involved;
                break;
            }
        }
    }

    for (int omitted = 0; omitted < 4; ++omitted) {
        const double sample_count = omitted == 0
            ? static_cast<double>(last_involved)
            : static_cast<double>(last_involved - 1);
        if (omitted > 0) output.sub_scores[omitted - 1] = 0.0;
        if (sample_count > 1.0) {
            for (int role = 0; role < 3; ++role) {
                if (role == omitted - 1) continue;
                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_xy = 0.0;
                double sum_x2 = 0.0;
                double sum_y2 = 0.0;
                for (int item = 0; item <= last_involved; ++item) {
                    if (item == omitted - 1) continue;
                    const int reference = output.trace_involved[role];
                    const int other = output.trace_involved[item];
                    if (reference == other) continue;
                    const auto index = static_cast<std::size_t>(reference) +
                        static_cast<std::size_t>(other) * count;
                    const double first = background_matrix[index];
                    const double second = region_matrix[index];
                    if (omitted == 0) {
                        output.sub_distance_scores[role] +=
                            std::fabs(first - second);
                    }
                    sum_x += first;
                    sum_y += second;
                    sum_xy += first * second;
                    sum_x2 += first * first;
                    sum_y2 += second * second;
                }

                double score = 1.0;
                if (sum_x2 > 0.0 && sum_y2 > 0.0 &&
                    !((sample_count * sum_y2) * 0.99999 < sum_y * sum_y &&
                      (sample_count * sum_y2) / 0.99999 > sum_y * sum_y) &&
                    !((sample_count * sum_x2) * 0.99999 < sum_x * sum_x &&
                      (sample_count * sum_x2) / 0.99999 > sum_x * sum_x)) {
                    const double denominator =
                        std::sqrt(sample_count * sum_x2 - sum_x * sum_x) *
                        std::sqrt(sample_count * sum_y2 - sum_y * sum_y);
                    score = (sample_count * sum_xy - sum_x * sum_y) /
                        denominator + minimum_offset;
                }
                output.temporary_scores[role] = score;
                if (omitted == 0) {
                    output.scores[role] = score;
                } else {
                    output.sub_scores[omitted - 1] += score;
                }
            }
        } else if (omitted == 0) {
            output.scores[0] = 1.0;
        } else {
            output.sub_scores[omitted - 1] = 1.0;
        }
        if (omitted > 0) output.sub_scores[omitted - 1] /= 2.0;
    }
    output.result = 1.0;
    return output;
}

std::vector<int> make_rdp_score_filter(
    const int next_no, const std::array<int, 3>& sequences,
    const std::vector<float>& raw_background_rows,
    const std::vector<float>& ancestor_background_rows,
    const std::vector<float>& ancestor_region_rows) {
    const int count = next_no + 1;
    const auto cells = static_cast<std::size_t>(3 * count);
    if (raw_background_rows.size() != cells ||
        ancestor_background_rows.size() != cells ||
        ancestor_region_rows.size() != cells) {
        throw std::runtime_error("MakeDoneThis3 input dimensions differ");
    }
    const auto cell = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };

    // MakeDoneThis3 is the selected-row form of MakeDoneThis2: its first
    // matrix index is a triplet role, not an absolute sequence number.
    float upper_background = 0.0F;
    float lower_background = 10000.0F;
    float upper_region = 0.0F;
    float lower_region = 10000.0F;
    for (int first = 0; first < 3; ++first) {
        for (int second = 0; second < 3; ++second) {
            const float background =
                ancestor_background_rows[cell(first, sequences[second])];
            upper_background = std::max(upper_background, background);
            lower_background = std::min(lower_background, background);
            const float region =
                ancestor_region_rows[cell(first, sequences[second])];
            upper_region = std::max(upper_region, region);
            lower_region = std::min(lower_region, region);
        }
    }

    std::vector<int> done(static_cast<std::size_t>(2 * count), 0);
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const float background =
                ancestor_background_rows[cell(role, sequence)];
            if (background < lower_background ||
                background > upper_background) {
                done[sequence * 2] = 1;
            }
            const float region = ancestor_region_rows[cell(role, sequence)];
            if (region < lower_region || region > upper_region) {
                done[sequence * 2 + 1] = 1;
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const float background =
                ancestor_background_rows[cell(role, sequence)];
            if (background > lower_background &&
                background < upper_background) {
                done[sequence * 2] = 0;
            }
            const float region = ancestor_region_rows[cell(role, sequence)];
            if (region > lower_region && region < upper_region) {
                done[sequence * 2 + 1] = 0;
            }
        }
    }
    for (const int sequence : sequences) {
        done[sequence * 2] = 1;
        done[sequence * 2 + 1] = 1;
    }
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (raw_background_rows[cell(0, sequence)] == 3.0F ||
            raw_background_rows[cell(1, sequence)] == 3.0F ||
            raw_background_rows[cell(2, sequence)] == 3.0F) {
            done[sequence * 2] = 1;
            done[sequence * 2 + 1] = 1;
        }
    }
    return done;
}

RdpTripletGroupState make_rdp_triplet_groups(
    const int role, const int next_no,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<float>& ancestor_background_rows,
    std::array<double, 3> minimum_distances) {
    const int count = next_no + 1;
    if (role < 0 || role > 2 ||
        ancestor_background_rows.size() !=
            static_cast<std::size_t>(3 * count)) {
        throw std::runtime_error("MakeTrpGroups2 input dimensions differ");
    }
    const auto cell = [](const int selected_role, const int sequence) {
        return static_cast<std::size_t>(selected_role) + sequence * 3;
    };
    RdpTripletGroupState output;
    output.counts.assign(count, 0);
    output.done.assign(count, 0);
    output.groups.assign(count, 0);
    output.minimum_distances = minimum_distances;
    output.minimum_distances[role] = 0.0;
    const int first_parent = sequences[comparison_matrix[role]];
    const int second_parent = sequences[comparison_matrix[role + 3]];
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        const float distance = ancestor_background_rows[cell(role, sequence)];
        if (distance < ancestor_background_rows[cell(role, first_parent)] &&
            distance < ancestor_background_rows[cell(role, second_parent)]) {
            output.done[sequence] = 1;
            output.groups[sequence] = 0;
            if (distance > output.minimum_distances[role]) {
                output.minimum_distances[role] = distance;
            }
        }
    }
    int group_number = 1;
    while (true) {
        double minimum = 10001.0;
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const float distance =
                ancestor_background_rows[cell(role, sequence)];
            if (output.done[sequence] == 0 &&
                distance < minimum &&
                distance > output.minimum_distances[role]) {
                minimum = distance;
            }
        }
        if (minimum == 10001.0) break;
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (ancestor_background_rows[cell(role, sequence)] == minimum) {
                output.groups[sequence] = group_number;
                output.done[sequence] = 1;
            }
        }
        ++group_number;
    }
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        ++output.counts[output.groups[sequence]];
    }
    output.result = 1;
    return output;
}

double make_rdp_triplet_tree_score(
    const int role, const int next_no,
    const std::array<int, 3>& sequences,
    const std::vector<float>& ancestor_background_rows,
    const std::vector<float>& ancestor_region_rows,
    const RdpTripletGroupState& groups) {
    const int count = next_no + 1;
    if (role < 0 || role > 2 ||
        ancestor_background_rows.size() !=
            static_cast<std::size_t>(3 * count) ||
        ancestor_region_rows.size() !=
            static_cast<std::size_t>(3 * count) ||
        groups.counts.size() != static_cast<std::size_t>(count) ||
        groups.groups.size() != static_cast<std::size_t>(count)) {
        throw std::runtime_error("MakeTrpScore2 input dimensions differ");
    }
    const auto cell = [](const int selected_role, const int sequence) {
        return static_cast<std::size_t>(selected_role) + sequence * 3;
    };
    (void)sequences;
    double score = 0.0;
    for (int first = 0; first <= next_no; ++first) {
        const float first_background =
            ancestor_background_rows[cell(role, first)];
        const float first_region = ancestor_region_rows[cell(role, first)];
        for (int second = first + 1; second <= next_no; ++second) {
            const float second_background =
                ancestor_background_rows[cell(role, second)];
            const float second_region =
                ancestor_region_rows[cell(role, second)];
            const bool changed =
                (second_background > first_background &&
                 second_region < first_region) ||
                (second_background < first_background &&
                 second_region > first_region) ||
                (second_background == first_background &&
                 second_region != first_region) ||
                (second_background != first_background &&
                 second_region == first_region);
            if (changed) {
                score += 1.0 / static_cast<double>(
                    groups.counts[groups.groups[first]] *
                    groups.counts[groups.groups[second]]);
            }
        }
    }
    return score;
}

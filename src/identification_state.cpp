#include "identification_state.hpp"

#include <cmath>
#include <stdexcept>

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

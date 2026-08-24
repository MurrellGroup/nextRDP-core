#include "identification_state.hpp"
#include "MathFuncsDll.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <iostream>
#include <sstream>
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

int make_var_map2(const int next_no, const int sequence_length,
                  const int variable_site_upper_bound,
                  const short* sequence_data, short* variable_site_map,
                  int* variable_region_position, int* variable_position,
                  const int* sequences, const int* comparison_matrix) {
    int variable_sites = 0;
    const int position_stride = sequence_length + 1;
    const int map_sequence_stride = (variable_site_upper_bound + 1) * 3;
    for (int position = 0; position <= sequence_length; ++position) {
        variable_region_position[position] = variable_sites;
        const int first = sequence_data[
            position + sequences[0] * position_stride];
        const int second = sequence_data[
            position + sequences[1] * position_stride];
        const int third = sequence_data[
            position + sequences[2] * position_stride];
        if (first == 46 || second == 46 || third == 46 ||
            (first == second && first == third)) {
            continue;
        }
        ++variable_sites;
        variable_position[variable_sites] = position;
        const int variable_offset = variable_sites * 3;
        for (int role = 0; role < 3; ++role) {
            const int representative = sequence_data[
                position + sequences[role] * position_stride];
            const int first_parent = sequence_data[
                position + sequences[comparison_matrix[role]] *
                    position_stride];
            const int second_parent = sequence_data[
                position + sequences[comparison_matrix[role + 3]] *
                    position_stride];
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                const int value = sequence_data[
                    position + sequence * position_stride];
                auto& destination = variable_site_map[
                    role + variable_offset + sequence * map_sequence_stride];
                if (representative == value) {
                    destination = 2;
                } else if (value != 46) {
                    if ((value == first_parent && representative == second_parent) ||
                        (representative == first_parent && value == second_parent)) {
                        destination = -1;
                    } else if ((value == first_parent && value != second_parent) ||
                               (value != first_parent && value == second_parent)) {
                        destination = 0;
                    } else if (value != first_parent && value != second_parent) {
                        destination = 1;
                    }
                }
            }
        }
    }
    return variable_sites;
}

void make_count_hits(const int beginning, const int ending,
                     const int smoothing_window, const int next_no,
                     const int variable_sites, const int sequence_length,
                     const int variable_site_upper_bound,
                     std::vector<double>& count_hits,
                     const std::vector<short>& variable_site_map,
                     std::vector<float>& variable_site_smooth,
                     const std::vector<int>& variable_region_position) {
    const int smooth_sequence_stride = (variable_sites + 1) * 3;
    const int map_sequence_stride = (variable_site_upper_bound + 1) * 3;
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            int total = 0;
            const int map_base = role + sequence * map_sequence_stride;
            const int smooth_base = role + sequence * smooth_sequence_stride;
            for (int offset = 1 - smoothing_window;
                 offset <= smoothing_window + 1; ++offset) {
                int site = offset;
                if (site < 1) site = variable_sites + site;
                else if (site > variable_sites) site -= variable_sites;
                total += variable_site_map[site * 3 + map_base];
            }
            variable_site_smooth[3 + smooth_base] = static_cast<float>(
                static_cast<double>(total) /
                ((smoothing_window * 2 + 1) * 2));
            for (int site = 2; site <= variable_sites; ++site) {
                int remove = site - smoothing_window - 1;
                if (remove < 1) remove = variable_sites + remove;
                else if (remove > variable_sites) remove -= variable_sites;
                total -= variable_site_map[remove * 3 + map_base];
                int add = site + smoothing_window;
                if (add > variable_sites) add -= variable_sites;
                total += variable_site_map[add * 3 + map_base];
                variable_site_smooth[site * 3 + smooth_base] =
                    static_cast<float>(static_cast<double>(total) /
                        ((smoothing_window * 2 + 1) * 2));
            }
        }
    }
    const auto accumulate = [&](const int first, const int last,
                                const int side, double& samples) {
        for (int site = first; site <= last; ++site) {
            samples += 1.0;
            const int site_offset = site * 3;
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                const int smooth_base =
                    site_offset + sequence * smooth_sequence_stride;
                for (int role = 0; role < 3; ++role) {
                    const float value = variable_site_smooth[
                        role + smooth_base];
                    if (value > 0.6F) {
                        count_hits[role + side * 3 + sequence * 6] +=
                            (value - 0.6) / 0.4;
                    }
                }
            }
        }
    };
    const auto divide = [&](const int side, const double samples) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            for (int role = 0; role < 3; ++role) {
                count_hits[role + side * 3 + sequence * 6] /= samples;
            }
        }
    };
    int stop = beginning > 1 ? beginning - 1 : beginning;
    if (beginning < ending) {
        double samples = 0.0;
        accumulate(1, variable_region_position[stop - 1], 0, samples);
        accumulate(variable_region_position[ending + 1],
                   variable_region_position[sequence_length], 0, samples);
        divide(0, samples);
        samples = 0.0;
        accumulate(variable_region_position[beginning],
                   variable_region_position[ending], 1, samples);
        divide(1, samples);
    } else {
        double samples = 0.0;
        accumulate(1, variable_region_position[ending], 1, samples);
        accumulate(variable_region_position[beginning],
                   variable_region_position[sequence_length], 1, samples);
        divide(1, samples);
        samples = 0.0;
        accumulate(variable_region_position[ending + 1],
                   variable_region_position[beginning - 1], 0, samples);
        divide(0, samples);
    }
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

RdpFinalTrimState run_rdp_final_trim_candidate_maintenance(
    const int next_no, const int rwinpp,
    const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<unsigned char, 2>& minimum_pair,
    const std::array<unsigned char, 3>& inside_roles,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::vector<unsigned char>& unfound,
    const std::vector<float>& correlations,
    const std::vector<float>& inversions,
    const std::vector<float>& local_distance_panels,
    const std::vector<float>& first_direct_small,
    const std::vector<float>& region_direct_small,
    const std::vector<float>& first_ancestor_small,
    const std::vector<float>& region_ancestor_small,
    const std::vector<float>& first_collapsed_small,
    const std::vector<float>& region_collapsed_small,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<double>& acceptable_sequences) {
    const int count = next_no + 1;
    const auto small_cells = static_cast<std::size_t>(3 * count);
    if (next_no < 2 || unfound.size() != small_cells ||
        correlations.size() != static_cast<std::size_t>(9 * count) ||
        inversions.size() != static_cast<std::size_t>(9 * count) ||
        local_distance_panels.size() != static_cast<std::size_t>(12 * count) ||
        first_direct_small.size() != small_cells ||
        region_direct_small.size() != small_cells ||
        first_ancestor_small.size() != small_cells ||
        region_ancestor_small.size() != small_cells ||
        first_collapsed_small.size() != small_cells ||
        region_collapsed_small.size() != small_cells ||
        candidate_list.size() != small_cells ||
        acceptable_sequences.size() != static_cast<std::size_t>(57 * count)) {
        throw std::runtime_error("FinalTrim candidate dimensions differ");
    }
    const auto row = [](const int role, const int value) {
        return static_cast<std::size_t>(role) + value * 3;
    };
    const auto rc = [](const int role, const int region, const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto dm = [](const int panel, const int role, const int sequence) {
        return static_cast<std::size_t>(panel) + role * 4 + sequence * 12;
    };
    const auto ok = [](const int role, const int category, const int sequence) {
        return static_cast<std::size_t>(role) + category * 3 + sequence * 57;
    };

    RdpFinalTrimState output;
    output.candidate_last = candidate_last;
    output.candidate_list = candidate_list;
    output.nonrecombinant_last.fill(0);
    output.nonrecombinant_list.assign(small_cells, 0);
    output.acceptable_sequences = acceptable_sequences;

    // Module2.FinalTrim 23398-23477: preserve the correlation cube, then
    // erase evidence from warned regions and inverted interpretations.
    auto trimmed_correlations = correlations;
    for (int region = 0; region < 3; ++region) {
        if (correlation_warnings[region] == 0) continue;
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                trimmed_correlations[rc(role, region, sequence)] = 0.0F;
            }
        }
    }
    for (std::size_t index = 0; index < inversions.size(); ++index) {
        if (inversions[index] != 0.0F) trimmed_correlations[index] = 0.0F;
    }
    std::vector<unsigned char> duplicates(small_cells, 0);
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
            const int sequence = output.candidate_list[row(role, slot)];
            if (sequence > next_no) continue;
            for (int region = 0; region < 3; ++region) {
                if (trimmed_correlations[rc(role, region, sequence)] > 0.83F) {
                    ++duplicates[row(region, sequence)];
                }
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
            const int sequence = output.candidate_list[row(role, slot)];
            if (sequence > next_no) continue;
            for (int region = 0; region < 3; ++region) {
                if (duplicates[row(region, sequence)] > 1) {
                    trimmed_correlations[rc(role, region, sequence)] = 0.0F;
                }
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
            const int sequence = output.candidate_list[row(role, slot)];
            if (sequence <= next_no &&
                (trimmed_correlations[rc(role, 0, sequence)] >= 0.83F ||
                 trimmed_correlations[rc(role, 1, sequence)] >= 0.83F ||
                 trimmed_correlations[rc(role, 2, sequence)] >= 0.83F)) {
                output.acceptable_sequences[ok(role, 5, sequence)] = 1.0;
            }
            output.acceptable_sequences[ok(role, 6, sequence)] = 1.0;
        }
    }

    const bool have_collapsed = !first_collapsed_small.empty();
    const auto& first_prune = have_collapsed
        ? first_collapsed_small : first_ancestor_small;
    const auto& region_prune = have_collapsed
        ? region_collapsed_small : region_ancestor_small;
    std::array<int, 4> old_last{};
    old_last.fill(0);

    // Module2.FinalTrim 23504-23831. VB deliberately compares only RNum,
    // not list contents, when deciding whether another maintenance pass is
    // needed; the same termination quirk is retained here.
    for (int role = 0; role < 3; ++role) {
        std::array<double, 2> nearest_distances{};
        while (old_last[role] != output.candidate_last[role]) {
            old_last[role] = output.candidate_last[role];
            std::vector<unsigned char> taken(count, 0);
            nearest_distances = {0.0, 0.0};
            for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
                const int sequence = output.candidate_list[row(role, slot)];
                if (sequence > next_no) continue;
                taken[sequence] = 1;
                nearest_distances[0] = std::max<double>(
                    nearest_distances[0], first_ancestor_small[row(role, sequence)]);
                nearest_distances[1] = std::max<double>(
                    nearest_distances[1], region_ancestor_small[row(role, sequence)]);
            }
            output.nonrecombinant_last[role] = -1;
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                if (first_ancestor_small[row(role, sequence)] <=
                        nearest_distances[0] && taken[sequence] == 0) {
                    ++output.nonrecombinant_last[role];
                    output.nonrecombinant_list[row(
                        role, output.nonrecombinant_last[role])] = sequence;
                    taken[sequence] = 1;
                }
                if (region_ancestor_small[row(role, sequence)] <=
                        nearest_distances[1] && taken[sequence] == 0) {
                    ++output.nonrecombinant_last[role];
                    output.nonrecombinant_list[row(
                        role, output.nonrecombinant_last[role])] = sequence;
                    taken[sequence] = 1;
                }
            }
            std::vector<unsigned char> remove(static_cast<std::size_t>(2 * count), 0);
            if (output.nonrecombinant_last[role] > -1) {
                nearest_distances[0] = 1000000.0;
                for (int slot = 0; slot <= output.nonrecombinant_last[role]; ++slot) {
                    const int sequence = output.nonrecombinant_list[row(role, slot)];
                    const float value = first_prune[row(role, sequence)];
                    // The source compares the ancestor value against a
                    // threshold formed using integer arithmetic.
                    if (value < nearest_distances[0] &&
                        first_ancestor_small[row(role, sequence)] <
                            static_cast<float>((next_no * 2) / 1000)) {
                        nearest_distances[0] = value;
                    }
                }
                auto mark_from_panel_pair = [&](const int sequence,
                                                const int first_panel,
                                                const int second_panel) {
                    int go_on = 0;
                    if (correlation_warnings[0] == 0) {
                        if (trimmed_correlations[rc(role, 0, sequence)] > 0.95F) {
                            for (const int panel : {first_panel, second_panel}) {
                                if (local_distance_panels[dm(panel, role, sequence)] >
                                        local_distance_panels[dm(
                                            panel, role,
                                            sequences[comparison_matrix[role]])] ||
                                    local_distance_panels[dm(panel, role, sequence)] >
                                        local_distance_panels[dm(
                                            panel, role,
                                            sequences[comparison_matrix[role + 3]])]) {
                                    if (local_distance_panels[dm(panel, role, sequence)] >
                                            local_distance_panels[dm(
                                                panel, comparison_matrix[role], sequence)] ||
                                        local_distance_panels[dm(panel, role, sequence)] >
                                            local_distance_panels[dm(
                                                panel, comparison_matrix[role + 3], sequence)]) {
                                        go_on += 2;
                                    } else {
                                        go_on += 1;
                                    }
                                } else {
                                    go_on -= 1;
                                }
                            }
                        } else {
                            go_on = 1;
                        }
                    } else {
                        go_on = 1;
                    }
                    return go_on;
                };
                auto prune_pass = [&](const std::vector<float>& matrix,
                                      const int evidence_region,
                                      const int first_panel,
                                      const int second_panel,
                                      const int mark_side,
                                      const double nearest) {
                    int slot = 0;
                    while (slot <= output.candidate_last[role]) {
                        const int sequence = output.candidate_list[row(role, slot)];
                        if (sequence > next_no ||
                            trimmed_correlations[rc(role, 0, sequence)] >= 0.99F ||
                            trimmed_correlations[rc(role, 1, sequence)] >= 0.99F) {
                            ++slot;
                            continue;
                        }
                        const float first_correlation =
                            trimmed_correlations[rc(role, 0, sequence)];
                        const float second_correlation =
                            trimmed_correlations[rc(role, 1, sequence)];
                        if (!((first_correlation < 0.99F && first_correlation > 0.83F) ||
                              (second_correlation < 0.99F && second_correlation > 0.83F))) {
                            ++slot;
                            continue;
                        }
                        const float value = matrix[row(role, sequence)];
                        const float parent0 = matrix[row(
                            role, sequences[comparison_matrix[role]])];
                        const float parent1 = matrix[row(
                            role, sequences[comparison_matrix[role + 3]])];
                        if (!(value > nearest || value > parent0 || value > parent1)) {
                            ++slot;
                            continue;
                        }
                        int go_on = 1;
                        if (first_correlation > 0.95F || second_correlation > 0.95F) {
                            go_on = 0;
                        }
                        if (unfound[row(role, sequence)] == 1 || go_on == 1) {
                            output.acceptable_sequences[ok(role, 6, sequence)] = 0.0;
                            if (slot < output.candidate_last[role]) {
                                output.candidate_list[row(role, slot)] =
                                    output.candidate_list[row(
                                        role, output.candidate_last[role])];
                            }
                            --output.candidate_last[role];
                            continue;
                        }
                        go_on = mark_from_panel_pair(
                            sequence, first_panel, second_panel);
                        if (go_on > 0 && correlation_warnings[evidence_region] == 0) {
                            go_on = 0;
                            if (trimmed_correlations[
                                    rc(role, evidence_region, sequence)] > 0.95F) {
                                for (const int panel : {first_panel, second_panel}) {
                                    if (local_distance_panels[dm(panel, role, sequence)] >
                                            local_distance_panels[dm(
                                                panel, role,
                                                sequences[comparison_matrix[role]])] ||
                                        local_distance_panels[dm(panel, role, sequence)] >
                                            local_distance_panels[dm(
                                                panel, role,
                                                sequences[comparison_matrix[role + 3]])]) {
                                        if (local_distance_panels[dm(panel, role, sequence)] >
                                                local_distance_panels[dm(
                                                    panel, comparison_matrix[role], sequence)] ||
                                            local_distance_panels[dm(panel, role, sequence)] >
                                                local_distance_panels[dm(
                                                    panel, comparison_matrix[role + 3], sequence)]) {
                                            go_on += 2;
                                        } else {
                                            go_on += 1;
                                        }
                                    } else {
                                        go_on -= 1;
                                    }
                                }
                            } else {
                                go_on = 1;
                            }
                        }
                        if (go_on > 0) remove[sequence + mark_side * count] = 1;
                        ++slot;
                    }
                };
                prune_pass(first_prune, 1, 2, 3, 0, nearest_distances[0]);

                nearest_distances[1] = 1000000.0;
                for (int slot = 0; slot <= output.nonrecombinant_last[role]; ++slot) {
                    const int sequence = output.nonrecombinant_list[row(role, slot)];
                    const float value = region_prune[row(role, sequence)];
                    if (value < nearest_distances[1] &&
                        value < static_cast<float>((next_no * 2) / 1000)) {
                        nearest_distances[1] = value;
                    }
                }
                prune_pass(region_prune, 1, 2, 3, 1, nearest_distances[1]);
            } else {
                nearest_distances = {0.0, 0.0};
            }
            int slot = 0;
            while (slot <= output.candidate_last[role]) {
                const int sequence = output.candidate_list[row(role, slot)];
                if (sequence <= next_no && remove[sequence] == 1 &&
                    remove[sequence + count] == 1) {
                    output.acceptable_sequences[ok(role, 6, sequence)] = 0.0;
                    if (slot < output.candidate_last[role]) {
                        output.candidate_list[row(role, slot)] =
                            output.candidate_list[row(
                                role, output.candidate_last[role])];
                    }
                    --output.candidate_last[role];
                } else {
                    ++slot;
                }
            }
        }

        // FinalTrim's two source-faithful near-sequence collection passes.
        std::vector<unsigned char> already_in(count, 0);
        for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
            const int sequence = output.candidate_list[row(role, slot)];
            if (sequence <= next_no) already_in[sequence] = 1;
        }
        for (const int sequence : sequences) already_in[sequence] = 1;
        double nearest_first = 1000000.0;
        double nearest_region = 1000000.0;
        for (int slot = 0; slot <= output.nonrecombinant_last[role]; ++slot) {
            const int sequence = output.nonrecombinant_list[row(role, slot)];
            nearest_first = std::min<double>(
                nearest_first, first_prune[row(role, sequence)]);
            nearest_region = std::min<double>(
                nearest_region, region_prune[row(role, sequence)]);
        }
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (already_in[sequence] == 0 &&
                first_prune[row(role, sequence)] <= nearest_first &&
                region_prune[row(role, sequence)] <= nearest_region &&
                (trimmed_correlations[rc(role, 0, sequence)] > 0.83F ||
                 trimmed_correlations[rc(role, 1, sequence)] > 0.83F ||
                 trimmed_correlations[rc(role, 2, sequence)] > 0.83F)) {
                ++output.candidate_last[role];
                output.candidate_list[row(role, output.candidate_last[role])] =
                    sequence;
                already_in[sequence] = 1;
            }
        }

        // Module2.FinalTrim 23873-24249. These six OKSeq panels are consumed
        // by ConsensusOK; omitting them changes which co-recombinants survive
        // even when all tree and correlation inputs are otherwise identical.
        std::array<double, 2> evidence_weight{1.0, 1.0};
        if (minimum_pair[0] == minimum_pair[1]) {
            const int pair = minimum_pair[0];
            const bool role_in_pair =
                (pair == 0 && (role == 0 || role == 1)) ||
                (pair == 1 && (role == 0 || role == 2)) ||
                (pair == 2 && (role == 1 || role == 2));
            if (!role_in_pair) evidence_weight = {0.5, 0.5};
        }
        const int other0 = comparison_matrix[role];
        const int other1 = comparison_matrix[role + 3];
        const int parent0 = sequences[other0];
        const int parent1 = sequences[other1];
        const auto add_relative_score = [](double& score, const float value,
                                           const float bound0,
                                           const float bound1,
                                           const double weight,
                                           const double tied_bonus) {
            if (value < bound0 && value < bound1) {
                score += 4.0 * weight;
            } else if (value <= bound0 && value <= bound1) {
                if (value < bound0 || value < bound1) {
                    score += tied_bonus * weight;
                }
            } else if (value > bound0 && value > bound1) {
                score -= 10.0 / weight;
            } else if (value > bound0 || value > bound1) {
                score -= 2.0 / weight;
            }
        };
        if (have_collapsed) {
            const float first_parent0 =
                first_collapsed_small[row(role, parent0)];
            const float first_parent1 =
                first_collapsed_small[row(role, parent1)];
            const float region_parent0 =
                region_collapsed_small[row(role, parent0)];
            const float region_parent1 =
                region_collapsed_small[row(role, parent1)];
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                auto& score = output.acceptable_sequences[
                    ok(role, 7, sequence)];
                const float first_value =
                    first_collapsed_small[row(role, sequence)];
                const float region_value =
                    region_collapsed_small[row(role, sequence)];
                add_relative_score(score, first_value, first_parent0,
                                   first_parent1, evidence_weight[0], 2.0);
                add_relative_score(score, region_value, region_parent0,
                                   region_parent1, evidence_weight[1], 2.0);
                add_relative_score(
                    score, first_value,
                    first_collapsed_small[row(other0, sequence)],
                    first_collapsed_small[row(other1, sequence)],
                    evidence_weight[0], 2.0);
                add_relative_score(
                    score, region_value,
                    region_collapsed_small[row(other0, sequence)],
                    region_collapsed_small[row(other1, sequence)],
                    evidence_weight[1], 2.0);
            }
        }
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            auto& score = output.acceptable_sequences[ok(role, 8, sequence)];
            const float first_value = first_ancestor_small[row(role, sequence)];
            const float region_value = region_ancestor_small[row(role, sequence)];
            add_relative_score(
                score, first_value, first_ancestor_small[row(role, parent0)],
                first_ancestor_small[row(role, parent1)], evidence_weight[0],
                2.0);
            // The source uses OKMod(0), not OKMod(1), for this parent-distance
            // SAMat contribution.
            add_relative_score(
                score, region_value, region_ancestor_small[row(role, parent0)],
                region_ancestor_small[row(role, parent1)], evidence_weight[0],
                2.0);
            add_relative_score(
                score, first_value, first_ancestor_small[row(other0, sequence)],
                first_ancestor_small[row(other1, sequence)], evidence_weight[0],
                1.0);
            add_relative_score(
                score, region_value,
                region_ancestor_small[row(other0, sequence)],
                region_ancestor_small[row(other1, sequence)],
                evidence_weight[1], 1.0);

        }

        const auto fill_direct_distance_evidence = [&](const int category,
                                                       const auto& value_at) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                auto& score = output.acceptable_sequences[
                    ok(role, category, sequence)];
                const float value0 = value_at(0, role, sequence);
                const float value1 = value_at(1, role, sequence);
                const int baseline_role = category == 9 ? other1 : other0;
                const int baseline_sequence = category == 9 ? parent0 : parent1;
                const float parent_value0 =
                    value_at(0, baseline_role, baseline_sequence);
                const float parent_value1 =
                    value_at(1, baseline_role, baseline_sequence);
                if (category != 9 && (value0 >= 3.0F || value1 >= 3.0F)) {
                    continue;
                }
                if (value0 < value_at(0, role, parent1) &&
                    value0 < value_at(0, role, parent0) &&
                    value1 < value_at(1, role, parent1) &&
                    value1 < value_at(1, role, parent0)) {
                    if (value0 < parent_value0) {
                        if (value1 < parent_value1) score = category == 9 ? 8.0 : 6.0;
                        else if (role != inside_roles[0]) score = 2.0;
                    } else if (value1 < parent_value1) {
                        if (role != inside_roles[0]) score = 2.0;
                    } else if (role != inside_roles[0]) {
                        score = category == 9 ? 0.5 : 1.0;
                    }
                } else if (value0 > value_at(0, other1, sequence) &&
                           value0 > value_at(0, other0, sequence) &&
                           value1 > value_at(1, other1, sequence) &&
                           value1 > value_at(1, other0, sequence)) {
                    if (value0 > parent_value0) {
                        score = value1 > parent_value1 ? -1.0 : -0.5;
                    } else {
                        score = value1 > parent_value1 ? -0.5 : -0.25;
                    }
                } else if (category != 9 &&
                           value0 > value_at(0, role, parent1) &&
                           value0 > value_at(0, role, parent0) &&
                           value1 > value_at(1, role, parent1) &&
                           value1 > value_at(1, role, parent0)) {
                    if (value0 > parent_value0) {
                        score = value1 > parent_value1 ? -1.0 : -0.5;
                    } else {
                        score = value1 > parent_value1 ? -0.5 : -0.25;
                    }
                }
                if (minimum_pair[0] != minimum_pair[1]) {
                    if ((role == inside_roles[2] || role == inside_roles[0]) &&
                        score > 0.0) score /= 2.0;
                } else if (score > 0.0) {
                    score *= evidence_weight[0];
                }
            }
        };
        fill_direct_distance_evidence(9, [&](const int panel, const int matrix_role,
                                             const int sequence) {
            return panel == 0
                ? first_direct_small[row(matrix_role, sequence)]
                : region_direct_small[row(matrix_role, sequence)];
        });
        if (correlation_warnings[0] == 0) {
            fill_direct_distance_evidence(12, [&](const int panel,
                                                  const int matrix_role,
                                                  const int sequence) {
                return local_distance_panels[
                    dm(panel, matrix_role, sequence)];
            });
        }
        if (correlation_warnings[1] == 0) {
            fill_direct_distance_evidence(13, [&](const int panel,
                                                  const int matrix_role,
                                                  const int sequence) {
                return local_distance_panels[
                    dm(panel + 2, matrix_role, sequence)];
            });
        }

        const float first_parent0 = first_ancestor_small[row(
            role, sequences[comparison_matrix[role]])];
        const float first_parent1 = first_ancestor_small[row(
            role, sequences[comparison_matrix[role + 3]])];
        const float region_parent0 = region_ancestor_small[row(
            role, sequences[comparison_matrix[role]])];
        const float region_parent1 = region_ancestor_small[row(
            role, sequences[comparison_matrix[role + 3]])];
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (already_in[sequence] == 0 &&
                first_ancestor_small[row(role, sequence)] < first_parent0 &&
                first_ancestor_small[row(role, sequence)] < first_parent1 &&
                region_ancestor_small[row(role, sequence)] < region_parent0 &&
                region_ancestor_small[row(role, sequence)] < region_parent1) {
                ++output.candidate_last[role];
                output.candidate_list[row(role, output.candidate_last[role])] =
                    sequence;
                already_in[sequence] = 1;
            }
        }
    }

    // Module2.FinalTrim 24321-24407.  The early RetrimFlag call uses
    // RWinPP=4 and skips this block; the ordinary final RFF=0 call uses the
    // already selected WinPP.  Preserve the source's missing `Y = 0` before
    // the third loop: it starts at the value left by the INList(1) loop.
    if (minimum_pair[0] != minimum_pair[1] && rwinpp < 3) {
        const auto candidate_value = [&](const std::vector<float>& collapsed,
                                         const std::vector<float>& ancestor,
                                         const int role,
                                         const int sequence) {
            return (have_collapsed ? collapsed : ancestor)[
                row(role, sequence)];
        };
        const auto correlations_below_ticket = [&](const int role,
                                                    const int sequence) {
            return correlations[rc(role, 0, sequence)] < 0.99F &&
                correlations[rc(role, 1, sequence)] < 0.99F &&
                correlations[rc(role, 2, sequence)] < 0.99F;
        };
        const auto prune_role = [&](const int role, const int first_parent,
                                    const int second_parent,
                                    const bool require_both,
                                    int slot) {
            while (slot <= output.candidate_last[role]) {
                const int sequence = output.candidate_list[row(role, slot)];
                const float first_value = candidate_value(
                    first_collapsed_small, first_ancestor_small,
                    role, sequence);
                const float first_bound = candidate_value(
                    first_collapsed_small, first_ancestor_small,
                    role, sequences[first_parent]);
                const float second_value = candidate_value(
                    region_collapsed_small, region_ancestor_small,
                    role, sequence);
                // The collapsed-matrix source branch uses the other role as
                // the S-matrix row for INList(0), unlike its ancestor branch.
                const int second_row = have_collapsed && require_both
                    ? first_parent : role;
                const int second_selected_role =
                    have_collapsed && require_both ? role : second_parent;
                const float second_bound = candidate_value(
                    region_collapsed_small, region_ancestor_small,
                    second_row, sequences[second_selected_role]);
                const bool keep = require_both
                    ? first_value <= first_bound && second_value <= second_bound
                    : first_value <= first_bound || second_value <= second_bound;
                if (!keep && correlations_below_ticket(role, sequence)) {
                    if (slot < output.candidate_last[role]) {
                        output.candidate_list[row(role, slot)] =
                            output.candidate_list[row(
                                role, output.candidate_last[role])];
                    }
                    --output.candidate_last[role];
                } else {
                    ++slot;
                }
            }
            return slot;
        };

        int slot = prune_role(inside_roles[0], inside_roles[1],
                              inside_roles[0], true, 0);
        slot = prune_role(inside_roles[1], inside_roles[0],
                          inside_roles[2], false, 0);
        (void)prune_role(inside_roles[2], inside_roles[1],
                         inside_roles[1], false, slot);
    }

    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= output.candidate_last[role]; ++slot) {
            output.acceptable_sequences[ok(
                role, 15, output.candidate_list[row(role, slot)])] = 1.0;
        }
    }
    return output;
}

std::vector<double> calculate_rdp_match_evidence(
    const int sequence_length, const int next_no, const int beginning,
    const int ending, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<short>& sequence_data,
    const std::vector<double>& acceptable_sequences,
    const bool conservative_grouping) {
    const int count = next_no + 1;
    const int position_stride = sequence_length + 1;
    if (sequence_length < 1 || next_no < 2 || beginning < 1 ||
        beginning > sequence_length || ending < 1 ||
        ending > sequence_length ||
        sequence_data.size() <
            static_cast<std::size_t>(position_stride * count) ||
        acceptable_sequences.size() != static_cast<std::size_t>(57 * count)) {
        throw std::runtime_error("CalcMatchY input dimensions differ");
    }
    const auto nucleotide = [&](const int position, const int sequence) {
        return sequence_data[position + sequence * position_stride];
    };
    const auto is_variable = [&](const int position) {
        return nucleotide(position, sequences[0]) != 46 &&
            nucleotide(position, sequences[1]) != 46 &&
            nucleotide(position, sequences[2]) != 46 &&
            (nucleotide(position, sequences[0]) !=
                 nucleotide(position, sequences[1]) ||
             nucleotide(position, sequences[0]) !=
                 nucleotide(position, sequences[2]));
    };
    const auto ok = [](const int role, const int category, const int sequence) {
        return static_cast<std::size_t>(role) + category * 3 + sequence * 57;
    };
    auto output = acceptable_sequences;
    int interval_variable_sites = 0;
    int position = beginning;
    do {
        if (is_variable(position)) ++interval_variable_sites;
        ++position;
        if (position > sequence_length) position = 1;
    } while (position != ending);
    const int target_window = interval_variable_sites >= 30
        ? 15 : static_cast<int>(std::nearbyint(interval_variable_sites / 2.0));
    if (target_window < 2) return output;

    std::array<int, 4> cycles{};
    std::array<int, 4> bounds{};
    int fragment_length = 0;
    int variable_sites = 0;
    position = beginning;
    do {
        ++fragment_length;
        if (is_variable(position) && ++variable_sites == 40) break;
        --position;
        if (position < 1) {
            position = sequence_length;
            ++cycles[0];
            if (cycles[0] > 40 || (cycles[0] > 2 && variable_sites == 0)) {
                break;
            }
        }
    } while (true);
    bounds[0] = position;
    if (fragment_length > sequence_length * 3) {
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                output[ok(role, 17, sequence)] = 0.0;
            }
        }
        return output;
    }
    const int local_beginning = fragment_length;
    position = beginning + 1;
    if (position > sequence_length) {
        position = 1;
        ++cycles[1];
    }
    variable_sites = 0;
    do {
        ++fragment_length;
        if (is_variable(position) && ++variable_sites == 40) break;
        ++position;
        if (position > sequence_length) {
            position = 1;
            ++cycles[1];
            if (cycles[1] > 40 || (cycles[1] > 2 && variable_sites == 0)) {
                break;
            }
        }
    } while (true);
    bounds[1] = position;
    position = ending;
    variable_sites = 0;
    do {
        ++fragment_length;
        if (is_variable(position) && ++variable_sites == 40) break;
        --position;
        if (position < 1) {
            position = sequence_length;
            ++cycles[2];
            if (cycles[2] > 40 || (cycles[2] > 2 && variable_sites == 0)) {
                break;
            }
        }
    } while (true);
    bounds[2] = position;
    if (fragment_length > sequence_length * 3) {
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                output[ok(role, 17, sequence)] = 0.0;
            }
        }
        return output;
    }
    const int local_ending = fragment_length;
    position = ending + 1;
    if (position > sequence_length) {
        position = 1;
        ++cycles[3];
    }
    variable_sites = 0;
    do {
        ++fragment_length;
        if (is_variable(position) && ++variable_sites == 40) break;
        ++position;
        if (position > sequence_length) {
            position = 1;
            ++cycles[3];
            if (cycles[3] > 40 || (cycles[3] > 2 && variable_sites == 0)) {
                break;
            }
        }
    } while (true);
    bounds[3] = position;
    if (fragment_length > sequence_length * 3) {
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                output[ok(role, 17, sequence)] = 0.0;
            }
        }
        return output;
    }

    const int sequence2_upper_bound = fragment_length + 1;
    std::vector<short> sequence2(
        static_cast<std::size_t>(sequence2_upper_bound + 1) * count, 0);
    fragment_length = MathFuncs::MyMathFuncs::MakeLenFrag(
        sequence_length, next_no, beginning, ending, cycles.data(),
        bounds.data(), sequence2_upper_bound, sequence2.data(),
        sequence_length, const_cast<short*>(sequence_data.data()));
    if (fragment_length > sequence_length * 3) {
        for (int role = 0; role < 3; ++role) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                output[ok(role, 17, sequence)] = 0.0;
            }
        }
        return output;
    }

    static unsigned int trace_invocation = 0;
    const unsigned int this_trace_invocation = ++trace_invocation;
    const char* trace_directory = std::getenv("RDP_CALCMATCH_TRACE_DIR");
    std::ofstream trace;
    const auto trace_write = [&](const auto* data, const std::size_t count) {
        if (trace.is_open()) {
            trace.write(reinterpret_cast<const char*>(data),
                        static_cast<std::streamsize>(count * sizeof(*data)));
        }
    };

    constexpr int map_upper_bound = 160;
    std::vector<int> variable_positions(sequence2_upper_bound + 1, 0);
    std::vector<int> variable_region_positions(sequence2_upper_bound + 2, 0);
    variable_region_positions[sequence2_upper_bound + 1] =
        sequence2_upper_bound + 1;
    std::vector<short> variable_site_map(
        static_cast<std::size_t>(3) * (map_upper_bound + 1) * count, 0);
    const int smoothing_window = std::max(target_window, 5);
    if (trace_directory != nullptr && *trace_directory != '\0') {
        std::ostringstream path;
        path << trace_directory << "/calcmatch-call-"
             << this_trace_invocation << ".bin";
        trace.open(path.str(), std::ios::binary | std::ios::trunc);
        const unsigned int header[] = {
            0x594d4343U, this_trace_invocation,
            static_cast<unsigned int>(sequence_length),
            static_cast<unsigned int>(next_no),
            static_cast<unsigned int>(beginning),
            static_cast<unsigned int>(ending),
            static_cast<unsigned int>(sequence2_upper_bound),
            static_cast<unsigned int>(local_beginning),
            static_cast<unsigned int>(local_ending),
            static_cast<unsigned int>(smoothing_window)};
        trace_write(header, 10);
        trace_write(sequences.data(), sequences.size());
        trace_write(comparison_matrix.data(), comparison_matrix.size());
        trace_write(sequence2.data(), sequence2.size());
    }
    variable_sites = make_var_map2(
        next_no, sequence2_upper_bound, map_upper_bound, sequence2.data(),
        variable_site_map.data(), variable_region_positions.data(),
        variable_positions.data(), sequences.data(), comparison_matrix.data());
    trace_write(&variable_sites, 1);
    trace_write(variable_site_map.data(), variable_site_map.size());
    trace_write(variable_region_positions.data(),
                variable_region_positions.size());
    trace_write(variable_positions.data(), variable_positions.size());
    std::vector<float> variable_site_smooth(
        static_cast<std::size_t>(3) * (variable_sites + 1) * count, 0.0F);
    std::vector<double> count_hits(static_cast<std::size_t>(6) * count, 0.0);
    make_count_hits(local_beginning, local_ending, smoothing_window, next_no,
                    variable_sites, sequence2_upper_bound, map_upper_bound,
                    count_hits, variable_site_map, variable_site_smooth,
                    variable_region_positions);
    trace_write(count_hits.data(), count_hits.size());
    trace_write(variable_site_smooth.data(), variable_site_smooth.size());
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            output[ok(role, 17, sequence)] =
                count_hits[role + sequence * 6] *
                count_hits[role + 3 + sequence * 6];
        }
    }

    const auto smooth = [&](const int role, int site, const int sequence) {
        const int maximum = variable_region_positions[sequence2_upper_bound];
        if (site < 0) site += maximum;
        else if (site > maximum) site -= maximum;
        return variable_site_smooth[
            role + site * 3 + sequence * 3 * (variable_sites + 1)];
    };
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const std::array<float, 6> track{
                smooth(role,
                       variable_region_positions[local_beginning] -
                           smoothing_window,
                       sequence),
                smooth(role,
                       variable_region_positions[local_beginning] +
                           smoothing_window,
                       sequence),
                smooth(role,
                       variable_region_positions[local_ending] -
                           smoothing_window,
                       sequence),
                smooth(role,
                       variable_region_positions[local_ending] +
                           smoothing_window,
                       sequence),
                smooth(role, variable_region_positions[local_beginning],
                       sequence),
                smooth(role, variable_region_positions[local_ending],
                       sequence),
            };
            double value = 0.0;
            if (!conservative_grouping) {
                if ((track[0] > 0.8F && track[1] > 0.8F &&
                     track[4] > 0.8F) ||
                    (track[2] > 0.8F && track[3] > 0.8F &&
                     track[5] > 0.8F)) {
                    value = 1.0;
                } else if (track[4] > 0.9F || track[5] > 0.9F) {
                    value = 1.0;
                } else if (track[4] > 0.8F || track[5] > 0.8F) {
                    value = 2.0;
                } else if ((track[4] > 0.75F || track[5] > 0.75F) &&
                           (track[1] > 0.75F || track[2] > 0.75F) &&
                           (track[0] > 0.75F || track[3] > 0.75F)) {
                    if ((track[0] + track[1] + track[4]) / 3.0F > 0.75F ||
                        (track[2] + track[3] + track[5]) / 3.0F > 0.75F ||
                        (track[0] > 0.75F && track[1] > 0.75F &&
                         track[4] > 0.75F) ||
                        (track[2] > 0.75F && track[3] > 0.75F &&
                         track[5] > 0.75F) ||
                        ((track[0] > 0.7F && track[1] > 0.7F &&
                          track[4] > 0.7F) &&
                         (track[2] > 0.7F && track[3] > 0.7F &&
                          track[5] > 0.7F))) {
                        value = 2.0;
                    }
                } else if ((track[4] < 0.65F && track[5] < 0.65F) ||
                           (track[1] < 0.65F && track[2] < 0.65F) ||
                           (track[0] < 0.65F && track[3] < 0.65F) ||
                           ((track[0] < 0.7F && track[1] < 0.7F &&
                             track[4] < 0.7F) &&
                            (track[2] < 0.7F && track[3] < 0.7F &&
                             track[5] < 0.7F)) ||
                           ((track[0] + track[1] + track[4]) / 3.0F < 0.7F &&
                            (track[2] + track[3] + track[5]) / 3.0F < 0.7F)) {
                    value = -1.0;
                }
            } else {
                if ((track[0] > 0.75F && track[1] > 0.75F &&
                     track[4] > 0.75F) ||
                    (track[2] > 0.75F && track[3] > 0.75F &&
                     track[5] > 0.75F)) {
                    value = 1.0;
                } else if (track[4] > 0.85F || track[5] > 0.85F) {
                    value = 1.0;
                } else if (track[4] > 0.75F || track[5] > 0.75F) {
                    value = 2.0;
                } else if ((track[4] > 0.7F || track[5] > 0.7F) &&
                           (track[1] > 0.7F || track[2] > 0.7F) &&
                           (track[0] > 0.7F || track[3] > 0.7F)) {
                    if ((track[0] + track[1] + track[4]) / 3.0F > 0.75F ||
                        (track[2] + track[3] + track[5]) / 3.0F > 0.75F ||
                        (track[0] > 0.7F && track[1] > 0.7F &&
                         track[4] > 0.7F) ||
                        (track[2] > 0.7F && track[3] > 0.7F &&
                         track[5] > 0.7F) ||
                        ((track[0] > 0.65F && track[1] > 0.65F &&
                          track[4] > 0.65F) &&
                         (track[2] > 0.65F && track[3] > 0.65F &&
                          track[5] > 0.65F))) {
                        value = 2.0;
                    }
                } else if ((track[4] < 0.6F && track[5] < 0.6F) ||
                           (track[1] < 0.6F && track[2] < 0.6F) ||
                           (track[0] < 0.6F && track[3] < 0.6F) ||
                           ((track[0] < 0.65F && track[1] < 0.65F &&
                             track[4] < 0.65F) &&
                            (track[2] < 0.65F && track[3] < 0.65F &&
                             track[5] < 0.65F)) ||
                           ((track[0] + track[1] + track[4]) / 3.0F < 0.6F &&
                            (track[2] + track[3] + track[5]) / 3.0F < 0.6F)) {
                    value = -1.0;
                }
            }
            output[ok(role, 18, sequence)] = value;
        }
    }
    return output;
}

RdpFinalTrimState make_rdp_consensus_candidates(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::vector<float>& correlations,
    const std::vector<float>& inversions,
    const std::vector<float>& first_direct,
    const std::vector<float>& region_direct,
    const std::vector<float>& first_ancestor,
    const std::vector<float>& region_ancestor,
    const std::vector<float>& first_direct_small,
    const std::vector<float>& region_direct_small,
    const std::vector<float>& first_ancestor_small,
    const std::vector<float>& region_ancestor_small,
    const std::vector<float>& first_collapsed_small,
    const std::vector<float>& region_collapsed_small,
    RdpFinalTrimState state, const bool conservative_grouping) {
    const int count = next_no + 1;
    const auto full_cells = static_cast<std::size_t>(count) * count;
    const auto small_cells = static_cast<std::size_t>(3 * count);
    if (first_direct.size() != full_cells || region_direct.size() != full_cells ||
        first_ancestor.size() != full_cells || region_ancestor.size() != full_cells ||
        first_direct_small.size() != small_cells ||
        region_direct_small.size() != small_cells ||
        first_ancestor_small.size() != small_cells ||
        region_ancestor_small.size() != small_cells ||
        first_collapsed_small.size() != small_cells ||
        region_collapsed_small.size() != small_cells ||
        correlations.size() != static_cast<std::size_t>(9 * count) ||
        inversions.size() != static_cast<std::size_t>(9 * count) ||
        state.candidate_list.size() != small_cells ||
        state.acceptable_sequences.size() !=
            static_cast<std::size_t>(57 * count)) {
        throw std::runtime_error("ConsensusOK input dimensions differ");
    }
    const auto small = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };
    const auto full = [count](const int first, const int second) {
        return static_cast<std::size_t>(first) + second * count;
    };
    const auto rc = [](const int role, const int region, const int sequence) {
        return static_cast<std::size_t>(role) + region * 3 + sequence * 9;
    };
    const auto ok = [](const int role, const int category, const int sequence) {
        return static_cast<std::size_t>(role) + category * 3 + sequence * 57;
    };
    auto& ok_sequences = state.acceptable_sequences;
    std::vector<float> maximum_correlation(small_cells, 0.0F);
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            for (int region = 0; region < 3; ++region) {
                if (correlation_warnings[region] == 0 &&
                    inversions[rc(role, region, sequence)] == 0.0F) {
                    maximum_correlation[small(role, sequence)] = std::max(
                        maximum_correlation[small(role, sequence)],
                        correlations[rc(role, region, sequence)]);
                }
            }
        }
    }
    std::vector<double> consensus_score(small_cells, 0.0);
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (sequence == sequences[role]) {
                consensus_score[small(role, sequence)] = 1000.0;
                continue;
            }
            if (sequence == sequences[0] || sequence == sequences[1] ||
                sequence == sequences[2]) {
                continue;
            }
            double score = 0.0;
            const double probability = ok_sequences[ok(role, 0, sequence)];
            if (probability > 0.0 && probability < 1.0) {
                score += -std::log10(probability) * 20.0;
            }
            score += ok_sequences[ok(role, 1, sequence)] * 5.0;
            score += ok_sequences[ok(role, 3, sequence)] * 4.0;
            if (ok_sequences[ok(role, 2, sequence)] == 1.0 &&
                ok_sequences[ok(role, 4, sequence)] == 1.0) {
                score *= 1.2;
            } else if (ok_sequences[ok(role, 2, sequence)] == 1.0 ||
                       ok_sequences[ok(role, 4, sequence)] == 1.0) {
                score *= 1.1;
            }
            if (ok_sequences[ok(role, 5, sequence)] == 1.0) score *= 1.1;
            if (ok_sequences[ok(role, 6, sequence)] == 1.0) score *= 1.1;
            if (ok_sequences[ok(role, 15, sequence)] == 1.0) score *= 2.0;
            const double corr = maximum_correlation[small(role, sequence)];
            score = (score + corr) * corr;
            // NS is a VB Long. Preserve conversion after each compound
            // assignment, including the surprising -0.5/2^-NS path.
            long evidence = static_cast<long>(std::nearbyint(
                (ok_sequences[ok(role, 7, sequence)] +
                 ok_sequences[ok(role, 8, sequence)]) / 2.0));
            for (int category = 9; category <= 14; ++category) {
                if (ok_sequences[ok(role, category, sequence)] != 0.0) {
                    evidence = static_cast<long>(std::nearbyint(
                        evidence + ok_sequences[ok(role, category, sequence)]));
                }
            }
            if (evidence < 1) {
                evidence = static_cast<long>(std::nearbyint(-0.5));
                evidence = 1 - evidence;
                evidence = static_cast<long>(std::nearbyint(
                    std::pow(2.0, -evidence)));
            }
            consensus_score[small(role, sequence)] = score * evidence;
        }
    }
    const auto original_last = state.candidate_last;
    const auto original_list = state.candidate_list;
    if (std::getenv("RDP_TRACE_CONSENSUS") != nullptr && next_no == 25) {
        std::cerr << "consensus-input seq=" << sequences[0] << ','
                  << sequences[1] << ',' << sequences[2] << " last=";
        for (const int value : original_last) std::cerr << value << ',';
        std::cerr << " role0=";
        for (int slot = 0; slot <= original_last[0]; ++slot)
            std::cerr << original_list[small(0, slot)] << ',';
        std::cerr << "\n";
        for (int role = 0; role < 3; ++role) {
            for (int sequence : {4, 5}) {
                std::cerr << "consensus-cell role=" << role << " seq="
                          << sequence << " topo="
                          << ok_sequences[ok(role, 18, sequence)]
                          << " match=" << ok_sequences[ok(role, 17, sequence)]
                          << " score=" << consensus_score[small(role, sequence)]
                          << " corr=" << maximum_correlation[small(role, sequence)]
                          << " direct=" << first_direct[full(sequences[role], sequence)]
                          << ',' << region_direct[full(sequences[role], sequence)]
                          << " collapsed=" << first_collapsed_small[small(role, sequence)]
                          << ',' << region_collapsed_small[small(role, sequence)]
                          << "\n";
            }
        }
    }
    state.candidate_last.fill(-1);

    // ConsensusOK's role-topology veto can be expressed uniformly through
    // CompMat even though the VB source spells out all three branches.
    for (int role = 0; role < 3; ++role) {
        const int other0 = comparison_matrix[role];
        const int other1 = comparison_matrix[role + 3];
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            auto& topology = ok_sequences[ok(role, 18, sequence)];
            if (topology <= -1.0) continue;
            if (!((first_ancestor_small[small(role, sequence)] <=
                       first_ancestor_small[small(other1, sequence)] ||
                   region_ancestor_small[small(role, sequence)] <=
                       region_ancestor_small[small(other1, sequence)]) &&
                  (first_ancestor_small[small(role, sequence)] <=
                       first_ancestor_small[small(other0, sequence)] ||
                   region_ancestor_small[small(role, sequence)] <=
                       region_ancestor_small[small(other0, sequence)]))) {
                topology = -1.0;
                continue;
            }
            const auto apply_cross_veto = [&](const std::vector<float>& matrix) {
                const float first = matrix[small(role, sequences[other0])];
                const float second = matrix[small(role, sequences[other1])];
                if ((first < second &&
                     matrix[small(other0, sequence)] >
                         matrix[small(other1, sequence)]) ||
                    (first > second &&
                     matrix[small(other0, sequence)] <
                         matrix[small(other1, sequence)])) {
                    topology = -1.0;
                }
            };
            apply_cross_veto(first_ancestor_small);
            apply_cross_veto(region_ancestor_small);
        }
    }

    const auto append = [&](const int role, const int sequence) {
        ++state.candidate_last[role];
        state.candidate_list[small(role, state.candidate_last[role])] = sequence;
    };
    for (int role = 0; role < 3; ++role) {
        const int other0 = comparison_matrix[role];
        const int other1 = comparison_matrix[role + 3];
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const double topology = ok_sequences[ok(role, 18, sequence)];
            if (topology <= -1.0) continue;
            const double match = ok_sequences[ok(role, 17, sequence)];
            const double score = consensus_score[small(role, sequence)];
            bool added = false;
            if ((((score > 30.0 && match > 0.05) || score * match > 2.0) &&
                 topology == 1.0) || (match > 0.1 && topology == 1.0)) {
                append(role, sequence);
                continue;
            }
            if (((score > 30.0 && match > 0.05) || score * match > 2.0) &&
                topology > 0.0) {
                for (int candidate = 0; candidate <= next_no; ++candidate) {
                    if (ok_sequences[ok(role, 18, candidate)] == 1.0 &&
                        candidate != sequence &&
                        region_ancestor_small[small(role, candidate)] >=
                            region_ancestor_small[small(role, sequence)] &&
                        first_ancestor_small[small(role, candidate)] >=
                            first_ancestor_small[small(role, sequence)]) {
                        append(role, sequence);
                        added = true;
                        break;
                    }
                }
            } else if (topology > 1.0) {
                for (int candidate = 0; candidate <= next_no; ++candidate) {
                    if (ok_sequences[ok(role, 18, candidate)] == 1.0 &&
                        candidate != sequence &&
                        region_ancestor_small[small(role, candidate)] >=
                            region_ancestor_small[small(role, sequence)] &&
                        first_ancestor_small[small(role, candidate)] >=
                            first_ancestor_small[small(role, sequence)]) {
                        append(role, sequence);
                        added = true;
                        break;
                    }
                }
            }
            // ConsensusOK does not stop at topology == 0 here.  Its fallback
            // FCMat/SCMat checks can still admit a candidate with a zero
            // breakpoint-topology score; only an explicitly vetoed -1 is
            // excluded above.
            if (added) continue;
            const bool collapsed_role_is_nearest =
                ((first_collapsed_small[small(role, sequence)] <=
                      first_collapsed_small[small(other0, sequence)] &&
                  first_collapsed_small[small(role, sequence)] <=
                      first_collapsed_small[small(other1, sequence)]) ||
                 (region_collapsed_small[small(role, sequence)] <=
                      region_collapsed_small[small(other0, sequence)] &&
                  region_collapsed_small[small(role, sequence)] <=
                      region_collapsed_small[small(other1, sequence)])) &&
                first_collapsed_small[small(role, sequences[other0])] >=
                    first_collapsed_small[small(other0, sequence)] &&
                first_collapsed_small[small(role, sequences[other1])] >=
                    first_collapsed_small[small(other1, sequence)] &&
                region_collapsed_small[small(role, sequences[other0])] >=
                    region_collapsed_small[small(other0, sequence)] &&
                region_collapsed_small[small(role, sequences[other1])] >=
                    region_collapsed_small[small(other1, sequence)];
            if (collapsed_role_is_nearest) {
                if (!((first_direct_small[small(other0, sequence)] == 0.0F &&
                       region_direct_small[small(other0, sequence)] == 0.0F) ||
                      (first_direct_small[small(other1, sequence)] == 0.0F &&
                       region_direct_small[small(other1, sequence)] == 0.0F))) {
                    append(role, sequence);
                    added = true;
                }
            } else if (first_direct[full(sequences[role], sequence)] == 0.0F &&
                       region_direct[full(sequences[role], sequence)] == 0.0F) {
                append(role, sequence);
                added = true;
            }
            if (added) continue;
            const bool equal_ancestor_pattern =
                first_ancestor_small[small(role, sequences[other0])] ==
                    first_ancestor_small[small(other0, sequence)] &&
                first_ancestor_small[small(role, sequences[other1])] ==
                    first_ancestor_small[small(other1, sequence)] &&
                region_ancestor_small[small(role, sequences[other0])] ==
                    region_ancestor_small[small(other0, sequence)] &&
                region_ancestor_small[small(role, sequences[other1])] ==
                    region_ancestor_small[small(other1, sequence)];
            // The x=0 source branch accidentally tests FMat(0,Y), while the
            // other roles test FMat(ISeqs(x),Y). Preserve that typo.
            const int zero_test_sequence = role == 0 ? 0 : sequences[role];
            if (equal_ancestor_pattern ||
                (first_direct[full(zero_test_sequence, sequence)] == 0.0F &&
                 region_direct[full(zero_test_sequence, sequence)] == 0.0F)) {
                append(role, sequence);
            }
        }

        std::vector<unsigned char> in_list(count, 0);
        for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
            in_list[state.candidate_list[small(role, slot)]] = 1;
        }
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const double topology = ok_sequences[ok(role, 18, sequence)];
            if (topology <= -1.0 || in_list[sequence] != 0 ||
                (topology <= 0.0 && !conservative_grouping)) {
                continue;
            }
            const double match = ok_sequences[ok(role, 17, sequence)];
            const double score = consensus_score[small(role, sequence)];
            bool added = false;
            if ((score > 30.0 && match > 0.05) || score * match > 2.0) {
                for (int candidate = 0; candidate <= next_no; ++candidate) {
                    if (in_list[candidate] != 0 &&
                        ok_sequences[ok(role, 18, candidate)] > 0.0 &&
                        candidate != sequence &&
                        region_ancestor_small[small(2, candidate)] ==
                            region_ancestor_small[small(2, sequence)] &&
                        first_ancestor_small[small(2, candidate)] ==
                            first_ancestor_small[small(2, sequence)] &&
                        region_ancestor_small[small(1, candidate)] ==
                            region_ancestor_small[small(1, sequence)] &&
                        first_ancestor_small[small(1, candidate)] ==
                            first_ancestor_small[small(1, sequence)] &&
                        region_ancestor_small[small(0, candidate)] ==
                            region_ancestor_small[small(0, sequence)] &&
                        first_ancestor_small[small(0, candidate)] ==
                            first_ancestor_small[small(0, sequence)]) {
                        append(role, sequence);
                        added = true;
                        break;
                    }
                }
            } else {
                for (int candidate = 0; candidate <= next_no; ++candidate) {
                    if (in_list[candidate] != 0 && candidate != sequence &&
                        region_direct[full(sequences[role], candidate)] >
                            region_direct[full(candidate, sequence)] &&
                        first_direct[full(sequences[role], candidate)] >
                            first_direct[full(candidate, sequence)]) {
                        append(role, sequence);
                        added = true;
                        break;
                    }
                }
            }
            if (added) in_list[sequence] = 1;
        }
        std::fill(in_list.begin(), in_list.end(), 0);
        for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
            in_list[state.candidate_list[small(role, slot)]] = 1;
        }
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (ok_sequences[ok(role, 18, sequence)] <= -1.0 ||
                in_list[sequence] != 0) {
                continue;
            }
            for (int candidate = 0; candidate <= next_no; ++candidate) {
                if (in_list[candidate] != 0 && candidate != sequence &&
                    region_ancestor_small[small(role, candidate)] >=
                        region_ancestor[full(candidate, sequence)] &&
                    first_ancestor_small[small(role, candidate)] >=
                        first_ancestor[full(candidate, sequence)] &&
                    region_direct_small[small(role, candidate)] >=
                        region_direct[full(candidate, sequence)] &&
                    first_direct_small[small(role, candidate)] >=
                        first_direct[full(candidate, sequence)]) {
                    append(role, sequence);
                    break;
                }
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        if (state.candidate_last[role] != -1) continue;
        state.candidate_last = original_last;
        state.candidate_list = original_list;
        break;
    }
    if (std::getenv("RDP_TRACE_CONSENSUS") != nullptr) {
        std::cerr << "consensus-output";
        for (int role = 0; role < 3; ++role) {
            std::cerr << " role" << role << '[';
            for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                std::cerr << state.candidate_list[small(role, slot)] << ',';
            }
            std::cerr << ']';
        }
        std::cerr << '\n';
    }
    return state;
}

RdpMaximumDistanceState calculate_rdp_maximum_distances(
    const int sequence_length, const int next_no,
    const std::array<int, 3>& sequences,
    const int beginning, const int ending,
    const std::vector<short>& sequence_data,
    const std::vector<unsigned char>& masked_sequences) {
    const int count = next_no + 1;
    const int stride = sequence_length + 1;
    if (sequence_length < 1 || next_no < 3 || beginning < 1 ||
        beginning > sequence_length || ending < 1 ||
        ending > sequence_length ||
        sequence_data.size() < static_cast<std::size_t>(stride * count) ||
        masked_sequences.size() != static_cast<std::size_t>(count)) {
        throw std::runtime_error("CalcMaxD input dimensions differ");
    }

    RdpMaximumDistanceState state;
    state.representative_mask.assign(count, 0);
    state.included_mask.assign(count, 0);
    for (const int sequence : sequences) {
        if (sequence < 0 || sequence > next_no) {
            throw std::runtime_error("CalcMaxD representative is out of range");
        }
        state.representative_mask[sequence] = 1;
        state.included_mask[sequence] = 1;
    }

    // Module5.CalcMaxD reduces StartSize for very large alignments before it
    // samples additional rows. All supplied demo runs remain in the branch
    // that includes every unmasked row, but retain the exact size arithmetic
    // so a future >30-row fixture can add the VB Rnd selection branch cleanly.
    int start_size = 0;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (masked_sequences[sequence] == 0) ++start_size;
    }
    start_size = std::min(start_size, 30);
    const double target_length =
        static_cast<double>(start_size) * (start_size - 1) *
        (start_size - 2) / 6.0 * 40000.0;
    while (static_cast<double>(start_size) * (start_size - 1) *
               (start_size - 2) / 6.0 * sequence_length >= target_length) {
        --start_size;
        if (start_size < 10) {
            start_size = 10;
            break;
        }
    }
    if (start_size > next_no - 3) start_size = next_no - 3;
    if (next_no <= start_size + 3) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (masked_sequences[sequence] == 0) {
                state.included_mask[sequence] = 1;
            }
        }
    } else {
        throw std::runtime_error(
            "CalcMaxD VB Rnd sampling requires a captured >30-row fixture");
    }
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (state.included_mask[sequence] != 0) {
            state.included_sequences.push_back(sequence);
        }
    }
    state.included_last =
        static_cast<int>(state.included_sequences.size()) - 1;

    state.informative_to_position.assign(stride, 0);
    state.position_to_informative.assign(stride, 0);
    int informative = 0;
    for (int position = 1; position <= sequence_length; ++position) {
        std::array<int, 4> counts{};
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const short nucleotide =
                sequence_data[position + sequence * stride];
            if (nucleotide == 66) ++counts[0];
            else if (nucleotide == 68) ++counts[1];
            else if (nucleotide == 72) ++counts[2];
            else if (nucleotide == 85) ++counts[3];
        }
        state.position_to_informative[position] = informative;
        int repeated_states = 0;
        for (const int value : counts) {
            if (value >= 2) ++repeated_states;
        }
        if (repeated_states >= 2) {
            ++informative;
            state.informative_to_position[informative] = position;
            state.position_to_informative[position] = informative;
        }
    }

    state.nucleotide_map.assign(86, 0);
    state.nucleotide_map[66] = 1;
    state.nucleotide_map[68] = 2;
    state.nucleotide_map[72] = 3;
    state.nucleotide_map[85] = 4;
    state.split_scores.assign(1875, 0.0F);
    // threshold.CPP::MakeVScoreMat. The source fills only states 1..4;
    // every tuple containing a missing-state zero remains zero-initialized.
    for (int first = 1; first < 5; ++first) {
        for (int second = 1; second < 5; ++second) {
            for (int third = 1; third < 5; ++third) {
                for (int fourth = 1; fourth < 5; ++fourth) {
                    const int cell = first + second * 5 + third * 25 +
                        fourth * 125;
                    if (first == second) {
                        if (first != third) {
                            if (third == fourth) {
                                state.split_scores[cell] = 1.0F;
                            } else if (fourth != first) {
                                state.split_scores[cell] = 0.5F;
                            }
                        }
                    } else if (first == third) {
                        if (first != fourth) {
                            state.split_scores[cell + 625] =
                                fourth == second ? 1.0F : 0.5F;
                        }
                    } else if (first == fourth) {
                        state.split_scores[cell + 1250] =
                            second == third ? 1.0F : 0.5F;
                    } else if (second == third) {
                        if (second != fourth) {
                            state.split_scores[cell + 1250] = 0.5F;
                        }
                    } else if (third == fourth) {
                        if (second != fourth) state.split_scores[cell] = 0.5F;
                    } else if (second == fourth) {
                        state.split_scores[cell + 625] = 0.5F;
                    }
                }
            }
        }
    }
    std::vector<short> sequence_scratch(
        static_cast<std::size_t>(stride) * count, 0);
    std::array<float, 3> outside{};
    std::array<float, 3> inside{};
    state.result = MathFuncs::MyMathFuncs::CMaxD2P3(
        state.included_last, sequences[0], sequences[1], sequences[2],
        beginning, ending, next_no, sequence_length,
        const_cast<short*>(sequence_data.data()), sequence_scratch.data(),
        state.informative_to_position.data(),
        state.position_to_informative.data(), state.nucleotide_map.data(),
        state.included_sequences.data(), state.included_mask.data(),
        state.representative_mask.data(), outside.data(), inside.data(),
        state.split_scores.data(), state.distance_totals.data(),
        state.distance_counts.data());
    if (state.distance_counts[0] > 0 && state.distance_counts[1] > 0 &&
        state.distance_counts[2] > 0) {
        for (int role = 0; role < 3; ++role) {
            state.maximum_distances[role] = state.distance_totals[role] /
                static_cast<float>(state.distance_counts[role]);
        }
    }
    return state;
}

RdpConsensusState make_rdp_consensus(RdpConsensusInputs inputs) {
    // Direct port of Module2.MakeConsensusC's active decision-tree path.
    // The logistic/NN alternatives are deliberately outside this routine:
    // the command-line oracle is configured with ConsensusStrat = 0.
    RdpConsensusState state;
    auto& v = inputs;
    const auto round = [](const double value, const double scale) {
        return std::nearbyint(value * scale) / scale;
    };
    for (int role = 0; role < 3; ++role) {
        v.triplet_score[role] = round(v.triplet_score[role], 100000.0);
        v.simple_distance_strength[role] =
            round(v.simple_distance_strength[role], 100000.0);
        v.phylpro[role] = round(v.phylpro[role], 100000.0);
        v.phylpro_secondary[role] =
            round(v.phylpro_secondary[role], 100000.0);
        v.phylpro_collapsed[role] =
            round(v.phylpro_collapsed[role], 100000.0);
        if (v.split_distance[role] < 10000.0) {
            v.split_distance[role] = round(v.split_distance[role], 100000.0);
        }
        v.subtree_phylpro[role] =
            round(v.subtree_phylpro[role], 10000.0);
        if (std::abs(v.subtree_score[role]) < 100.0) {
            v.subtree_score[role] = round(v.subtree_score[role], 1000000.0);
        }
        v.subtree_phylpro_secondary[role] =
            round(v.subtree_phylpro_secondary[role], 100000.0);
        if (v.subtree_score_secondary[role] > 100000.0) {
            v.subtree_score_secondary[role] = 100000.0;
        }
        v.subtree_score_secondary[role] =
            round(v.subtree_score_secondary[role], 100000.0);
        v.list_correlation[role] =
            round(v.list_correlation[role], 100000.0);
        if (v.list_correlation_secondary[role] > 0.0) {
            v.list_correlation_secondary[role] =
                round(v.list_correlation_secondary[role], 100000.0);
        }
        if (v.list_correlation_tertiary[role] > 0.0) {
            v.list_correlation_tertiary[role] =
                round(v.list_correlation_tertiary[role], 100000.0);
        }
    }

    int ps1 = 3;
    int ps2 = 3;
    int ps3 = 3;
    for (int role = 0; role < 3; ++role) {
        if (std::abs(v.phylpro[role]) < 0.999999) ++ps1;
        if (std::abs(v.phylpro_secondary[role]) < 0.999999) ++ps2;
        if (std::abs(v.phylpro_collapsed[role]) < 0.999999) ++ps3;
    }
    for (int role = 0; role < 3; ++role) {
        if (std::abs(v.phylpro[role]) == 1.0) ps1 = 0;
        if (std::abs(v.phylpro_secondary[role]) == 1.0) ps2 = 0;
        if (std::abs(v.phylpro_collapsed[role]) == 1.0) ps3 = 0;
    }

    auto score = [&](const int index, const int role) -> double& {
        return state.decision_scores[index + role * 26];
    };
    const auto parent = [&](const int role, const int which) {
        return v.comparison_matrix[role + which * 3];
    };
    const auto all_equal = [](const auto& values) {
        return values[0] == values[1] && values[0] == values[2];
    };
    const auto add_min = [&](const auto& values, const int role,
                             const double first_prize,
                             const double second_prize) {
        const int p0 = parent(role, 0);
        const int p1 = parent(role, 1);
        if (values[role] <= values[p0] && values[role] <= values[p1]) {
            state.consensus[role] += first_prize;
            return first_prize;
        }
        if (values[role] < values[p0] || values[role] < values[p1]) {
            state.consensus[role] += second_prize;
            return second_prize;
        }
        return 0.0;
    };
    const auto add_max = [&](const auto& values, const int role,
                             const double first_prize,
                             const double second_prize) {
        const int p0 = parent(role, 0);
        const int p1 = parent(role, 1);
        if (values[role] >= values[p0] && values[role] >= values[p1]) {
            state.consensus[role] += first_prize;
            return first_prize;
        }
        if (values[role] > values[p0] || values[role] > values[p1]) {
            state.consensus[role] += second_prize;
            return second_prize;
        }
        return 0.0;
    };

    for (int role = 0; role < 3; ++role) {
        const int p0 = parent(role, 0);
        const int p1 = parent(role, 1);
        const double dtotal = v.maximum_distance[p0] +
            v.maximum_distance[role] + v.maximum_distance[p1];
        score(10, role) = dtotal > 0.0
            ? v.maximum_distance[role] / dtotal * 20.0 : 0.0;
        state.consensus[role] = score(10, role);
        if (v.maximum_distance[role] >= v.maximum_distance[p0] * 1.1 &&
            v.maximum_distance[role] >= v.maximum_distance[p1] * 1.1) {
            state.consensus[role] += 30.0;
        } else if (v.maximum_distance[role] >= v.maximum_distance[p0] &&
                   v.maximum_distance[role] >= v.maximum_distance[p1]) {
            state.consensus[role] += 20.0;
        } else if (v.maximum_distance[role] >= v.maximum_distance[p0] * 1.1 ||
                   v.maximum_distance[role] >= v.maximum_distance[p1] * 1.1) {
            state.consensus[role] += 10.0;
        } else if (v.maximum_distance[role] >= v.maximum_distance[p0] ||
                   v.maximum_distance[role] >= v.maximum_distance[p1]) {
            state.consensus[role] += 5.0;
        }
        if (v.outlier_index[role] == 1) {
            state.consensus[role] += 5.0;
            score(16, role) = 5.0;
        }
        if (v.split_distance[role] != v.split_distance[p0] &&
            v.split_distance[role] != v.split_distance[p1]) {
            if (v.split_distance[role] >= v.split_distance[p0] &&
                v.split_distance[role] >= v.split_distance[p1]) {
                state.consensus[role] += 5.0;
                score(15, role) = 5.0;
            } else if (v.split_distance[role] > v.split_distance[p0] ||
                       v.split_distance[role] > v.split_distance[p1]) {
                state.consensus[role] += 2.5;
                score(15, role) = 2.5;
            }
        }
        if (v.outlier_check[role] > 0) {
            state.consensus[role] += 5.0;
            score(5, role) = 10.0;
        } else if (v.outlier_check[role] < 0) {
            state.consensus[role] -= 5.0;
            score(5, role) = 0.0;
        } else {
            score(5, role) = 5.0;
        }
        if (v.bad_distances[role] == v.bad_distances[p0] &&
            v.bad_distances[role] == v.bad_distances[p1]) {
            state.consensus[role] += 1.0;
            score(2, role) = 1.0;
        } else {
            score(2, role) = add_min(v.bad_distances, role, 10.0, 5.0);
        }

        const auto add_compatibility_ladder = [&](const auto& first,
                                                   const auto& second,
                                                   const auto& third,
                                                   const auto& fourth,
                                                   const int score_index) {
            const auto* selected = &first;
            if (all_equal(first)) selected = &second;
            if (all_equal(first) && all_equal(second)) selected = &third;
            if (all_equal(first) && all_equal(second) && all_equal(third)) {
                selected = &fourth;
            }
            if (!all_equal(*selected)) {
                score(score_index, role) =
                    add_min(*selected, role, 20.0, 10.0);
            }
        };
        add_compatibility_ladder(
            v.compatibility, v.compatibility_secondary,
            v.compatibility_tertiary, v.compatibility_quaternary, 1);
        add_compatibility_ladder(
            v.region_compatibility, v.region_compatibility_secondary,
            v.region_compatibility_tertiary,
            v.region_compatibility_quaternary, 14);

        if (v.permanent_next_no > 10) {
            score(9, role) = add_max(v.triplet_score, role, 8.0, 4.0);
        }
        if (ps1 >= 2) {
            if ((v.subtree_phylpro[0] != -1.0 &&
                 v.subtree_phylpro[0] != 1.0) ||
                (v.subtree_phylpro[1] != -1.0 &&
                 v.subtree_phylpro[1] != 1.0) ||
                (v.subtree_phylpro[2] != -1.0 &&
                 v.subtree_phylpro[2] != 1.0)) {
                score(3, role) =
                    add_max(v.subtree_phylpro, role, 10.0, 5.0);
            }
            score(13, role) = add_max(v.subtree_score, role, 2.0, 1.0);
            score(8, role) = add_min(v.phylpro, role, 8.0, 4.0);
        }
        if (ps2 >= 2) {
            score(6, role) = add_max(
                v.subtree_phylpro_secondary, role, 8.0, 4.0);
            score(4, role) = add_max(
                v.subtree_score_secondary, role, 10.0, 5.0);
            score(7, role) = add_min(
                v.phylpro_secondary, role, 18.0, 14.0);
        }
        score(12, role) = add_min(v.list_correlation, role, 4.0, 2.0);
        if (v.list_correlation_secondary[role] > 0.0 && v.next_no > 10) {
            score(17, role) = add_max(
                v.list_correlation_secondary, role, 20.0, 15.0);
        }

        constexpr double super_maximum = 20.0;
        constexpr double maximum = 10.0;
        constexpr double middle = 5.0;
        if (ps1 >= 2) {
            if (v.subtree_phylpro[role] >= v.subtree_phylpro[p0] &&
                v.subtree_phylpro[role] >= v.subtree_phylpro[p1] &&
                ps3 >= 2 &&
                v.phylpro_collapsed[role] <= v.phylpro_collapsed[p0] &&
                v.phylpro_collapsed[role] <= v.phylpro_collapsed[p1]) {
                state.consensus[role] += super_maximum;
                score(21, role) = super_maximum / 2.0;
            }
            if (v.phylpro[role] <= v.phylpro[p0] &&
                v.phylpro[role] <= v.phylpro[p1] &&
                (v.compatibility[role] != v.compatibility[p0] ||
                 v.compatibility[role] != v.compatibility[p1]) &&
                v.compatibility[role] <= v.compatibility[p0] &&
                v.compatibility[role] <= v.compatibility[p1]) {
                state.consensus[role] += maximum;
                score(19, role) = maximum / 2.0;
            }
        }
        if ((v.region_compatibility[role] != v.region_compatibility[p0] ||
             v.region_compatibility[role] != v.region_compatibility[p1]) &&
            v.region_compatibility[role] <= v.region_compatibility[p0] &&
            v.region_compatibility[role] <= v.region_compatibility[p1] &&
            v.triplet_score[role] >= v.triplet_score[p0] &&
            v.triplet_score[role] >= v.triplet_score[p1]) {
            state.consensus[role] += maximum;
            score(25, role) = maximum / 2.0;
        }
        if ((v.outlier_check[role] != v.outlier_check[p0] ||
             v.outlier_check[role] != v.outlier_check[p1]) &&
            v.outlier_check[role] >= v.outlier_check[p0] &&
            v.outlier_check[role] >= v.outlier_check[p1] &&
            v.simple_distance_score[role] == 1) {
            state.consensus[role] += maximum;
            score(22, role) = super_maximum / 2.0;
        }
        if (ps2 >= 2 &&
            v.subtree_phylpro_secondary[role] >=
                v.subtree_phylpro_secondary[p0] &&
            v.subtree_phylpro_secondary[role] >=
                v.subtree_phylpro_secondary[p1] &&
            v.simple_distance_strength[role] >=
                v.simple_distance_strength[p0] &&
            v.simple_distance_strength[role] >=
                v.simple_distance_strength[p1]) {
            state.consensus[role] += middle;
            score(23, role) = middle / 2.0;
        }
        if (ps2 >= 2 &&
            v.subtree_score_secondary[role] >=
                v.subtree_score_secondary[p0] &&
            v.subtree_score_secondary[role] >=
                v.subtree_score_secondary[p1] &&
            (v.bad_distances[role] != v.bad_distances[p0] ||
             v.bad_distances[role] != v.bad_distances[p1]) &&
            v.bad_distances[role] <= v.bad_distances[p0] &&
            v.bad_distances[role] <= v.bad_distances[p1]) {
            state.consensus[role] += middle;
            score(24, role) = middle / 3.0;
        }
        score(11, role) = static_cast<double>(
            v.ranks[role][1] - v.ranks[role][0]) / v.next_no * 10.0;
        if (score(11, role) < 0.0) score(11, role) = 0.0;
    }

    for (int role = 0; role < 3; ++role) {
        const int penalty = std::min(
            v.post_trim_compatibility[role],
            v.post_trim_region_compatibility[role]);
        if (penalty > 1) state.consensus[role] -= (penalty - 1) * 10.0;
    }
    double lowest = 0.0;
    for (const double value : state.consensus) {
        if (value < lowest) lowest = value;
    }
    for (double& value : state.consensus) value -= lowest;
    double total = state.consensus[0] + state.consensus[1] +
        state.consensus[2];
    if (total == 0.0) total = 1.0;
    for (int role = 0; role < 3; ++role) {
        score(0, role) = state.consensus[role] / total * 30.0;
    }

    int winner = 0;
    for (int role = 0; role < 3; ++role) {
        if (state.consensus[role] >= state.consensus[parent(role, 0)] &&
            state.consensus[role] >= state.consensus[parent(role, 1)]) {
            winner = role;
        }
    }
    if (state.consensus[0] == state.consensus[1] &&
        state.consensus[0] == state.consensus[2]) {
        const bool phylpro_informative =
            (v.phylpro[0] != -1.0 && v.phylpro[0] != 1.0) ||
            (v.phylpro[1] != -1.0 && v.phylpro[1] != 1.0) ||
            (v.phylpro[2] != -1.0 && v.phylpro[2] != 1.0);
        if (phylpro_informative) {
            std::array<double, 3> tie_break{};
            for (int role = 0; role < 3; ++role) {
                tie_break[role] = v.phylpro[role] - v.triplet_score[role] +
                    v.compatibility[role] + v.region_compatibility[role];
            }
            if (tie_break[0] < tie_break[1] && tie_break[0] < tie_break[2]) {
                winner = 0;
            } else if (tie_break[1] < tie_break[0] &&
                       tie_break[1] < tie_break[2]) {
                winner = 1;
            } else {
                winner = 2;
            }
        }
    }
    state.winning_role = winner;
    state.rounded_inputs = std::move(inputs);
    return state;
}

RdpSplitDistanceState calculate_rdp_split_distances(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::vector<float>& ancestor_rows,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix,
    const std::vector<int>& done) {
    const int count = next_no + 1;
    if (ancestor_rows.size() != static_cast<std::size_t>(3 * count) ||
        background_matrix.size() != static_cast<std::size_t>(count * count) ||
        region_matrix.size() != static_cast<std::size_t>(count * count) ||
        done.size() != static_cast<std::size_t>(2 * count)) {
        throw std::runtime_error("MakeSSDistB input dimensions differ");
    }
    const auto full = [count](const int row, const int column) {
        return row + column * count;
    };
    const auto small = [](const int role, const int sequence) {
        return role + sequence * 3;
    };
    float background_total = 0.0F;
    float region_total = 0.0F;
    for (int first = 0; first < next_no; ++first) {
        for (int second = first + 1; second <= next_no; ++second) {
            if ((done[first * 2] == 0 || done[first * 2 + 1] == 0) &&
                (done[second * 2] == 0 || done[second * 2 + 1] == 0)) {
                background_total += background_matrix[full(first, second)];
                region_total += region_matrix[full(first, second)];
            }
        }
    }
    RdpSplitDistanceState state;
    if (region_total <= 0.0F) return state;
    const float adjustment = background_total / region_total;
    for (int role = 0; role < 3; ++role) {
        float maximum = 0.0F;
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const float candidate = static_cast<float>(
                ancestor_rows[small(role, sequence)] * 1000.0F);
            if (candidate > maximum) maximum = candidate;
        }
        const int maximum_bin =
            static_cast<int>(std::nearbyint(maximum + 1.0F));
        std::vector<float> sums(maximum_bin + 1, 0.0F);
        std::vector<int> sizes(maximum_bin + 1, 0);
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (done[sequence * 2] != 0 && done[sequence * 2 + 1] != 0) {
                continue;
            }
            const int bin = static_cast<int>(std::nearbyint(
                ancestor_rows[small(role, sequence)] * 1000.0F));
            float distance = std::abs(
                background_matrix[full(sequences[role], sequence)] -
                region_matrix[full(sequences[role], sequence)]) * adjustment;
            distance *= distance;
            if (bin <= maximum_bin) {
                sums[bin] += distance;
                ++sizes[bin];
            }
        }
        for (int bin = 0; bin <= maximum_bin; ++bin) {
            if (sizes[bin] > 0) {
                state.distances[role] += sums[bin] / sizes[bin];
            }
        }
    }
    const int outside = inside_roles[0];
    const int inside = inside_roles[1];
    const int remaining = inside_roles[2];
    if (state.distances[outside] > state.distances[inside] &&
        state.distances[outside] > state.distances[remaining]) {
        state.outlier_index[outside] = 1;
    } else if (state.distances[outside] < state.distances[inside] &&
               state.distances[outside] < state.distances[remaining]) {
        state.outlier_index[inside] = 1;
        state.outlier_index[remaining] = 1;
    }
    return state;
}

RdpSimpleDistanceState calculate_rdp_simple_distances(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix) {
    const int count = next_no + 1;
    const auto full = [count](const int row, const int column) {
        return row + column * count;
    };
    std::vector<unsigned char> dont_use(count, 0);
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            dont_use[candidate_list[role + slot * 3]] = 1;
        }
    }
    float background_total = 0.0F;
    float region_total = 0.0F;
    for (int first = 0; first < next_no; ++first) {
        for (int second = first + 1; second <= next_no; ++second) {
            if (dont_use[first] == 0 && dont_use[second] == 0 &&
                background_matrix[full(first, second)] < 3.0F) {
                background_total += background_matrix[full(first, second)];
                region_total += region_matrix[full(first, second)];
            }
        }
    }
    std::vector<float> move_background(count, 0.0F);
    std::vector<float> move_region(count, 0.0F);
    for (int first = 0; first <= next_no; ++first) {
        for (int second = 0; second <= next_no; ++second) {
            move_background[first] += background_matrix[full(first, second)];
            move_region[first] += region_matrix[full(first, second)];
        }
    }
    RdpSimpleDistanceState state;
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (move_background[sequences[role]] > move_background[sequence]) {
                ++state.ranks[role][0];
            }
            if (move_region[sequences[role]] > move_region[sequence]) {
                ++state.ranks[role][1];
            }
        }
    }
    const double adjustment = background_total > 0.0F && region_total > 0.0F
        ? background_total / region_total : 1.0;
    const auto f = [&](const int first, const int second) {
        return static_cast<double>(background_matrix[
            full(sequences[first], sequences[second])]);
    };
    const auto s = [&](const int first, const int second) {
        return static_cast<double>(region_matrix[
            full(sequences[first], sequences[second])]) * adjustment;
    };
    for (int test = 0; test < 3; ++test) {
        const int target = inside_roles[test];
        double baseline = 0.0;
        double first_difference = 0.0;
        double second_difference = 0.0;
        if (test == 0) {
            baseline = std::abs(f(inside_roles[1], inside_roles[2]) -
                                s(inside_roles[1], inside_roles[2]));
            first_difference = s(inside_roles[0], inside_roles[1]) -
                f(inside_roles[0], inside_roles[1]);
            second_difference = s(inside_roles[0], inside_roles[2]) -
                f(inside_roles[0], inside_roles[2]);
        } else if (test == 1) {
            baseline = std::abs(f(inside_roles[0], inside_roles[2]) -
                                s(inside_roles[0], inside_roles[2]));
            first_difference = s(inside_roles[0], inside_roles[1]) -
                f(inside_roles[0], inside_roles[1]);
            second_difference = f(inside_roles[1], inside_roles[2]) -
                s(inside_roles[1], inside_roles[2]);
        } else {
            baseline = std::abs(f(inside_roles[0], inside_roles[1]) -
                                s(inside_roles[0], inside_roles[1]));
            first_difference = f(inside_roles[1], inside_roles[2]) -
                s(inside_roles[1], inside_roles[2]);
            second_difference = f(inside_roles[0], inside_roles[2]) -
                s(inside_roles[0], inside_roles[2]);
        }
        if (baseline < first_difference && baseline < second_difference) {
            state.scores[target] += 1;
        }
        state.strengths[target] =
            first_difference + second_difference - baseline;
    }
    return state;
}

std::array<int, 3> calculate_rdp_outlier_checks(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& inside_roles,
    const std::vector<float>& background_ancestor_rows,
    const std::vector<float>& region_ancestor_rows) {
    const auto at = [](const int role, const int sequence) {
        return role + sequence * 3;
    };
    std::array<int, 3> result{};
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (background_ancestor_rows[at(inside_roles[0], sequence)] >
                background_ancestor_rows[at(
                    inside_roles[0], sequences[inside_roles[1]])] &&
            background_ancestor_rows[at(inside_roles[0], sequence)] <
                background_ancestor_rows[at(
                    inside_roles[0], sequences[inside_roles[2]])]) {
            if (region_ancestor_rows[at(inside_roles[1], sequence)] <
                region_ancestor_rows[at(inside_roles[0], sequence)]) {
                ++result[inside_roles[0]];
                --result[inside_roles[1]];
                --result[inside_roles[2]];
            } else if (region_ancestor_rows[at(inside_roles[1], sequence)] >
                       region_ancestor_rows[at(inside_roles[0], sequence)]) {
                --result[inside_roles[0]];
                ++result[inside_roles[1]];
                --result[inside_roles[2]];
            } else if (region_ancestor_rows[at(inside_roles[1], sequence)] >
                       region_ancestor_rows[at(
                           inside_roles[0], sequences[inside_roles[1]])]) {
                --result[inside_roles[0]];
                --result[inside_roles[1]];
                ++result[inside_roles[2]];
            }
        }
    }
    return result;
}

std::array<double, 3> calculate_rdp_bad_distances(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<unsigned char>& unfound,
    const std::vector<float>& correlations,
    const std::vector<float>& ancestor_rows,
    const std::vector<float>& local_distance_panels) {
    const auto corr = [&](const int role, const int region, const int sequence) {
        return correlations[role + region * 3 + sequence * 9];
    };
    const auto dist = [&](const int panel, const int role, const int sequence) {
        return local_distance_panels[panel + role * 4 + sequence * 12];
    };
    const auto comp = [&](const int role, const int parent_index) {
        return comparison_matrix[role + parent_index * 3];
    };
    std::array<double, 3> result{};
    for (int role = 0; role < 3; ++role) {
        std::vector<unsigned char> bins(next_no + 1, 0);
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            const int sequence = candidate_list[role + slot * 3];
            if (unfound[role + sequence * 3] != 0) continue;
            bool bad = false;
            if (corr(role, 0, sequence) > 0.83F) {
                bad = dist(0, role, sequence) >
                          dist(0, role, sequences[comp(role, 1)]) ||
                    dist(1, role, sequence) >=
                          dist(1, role, sequences[comp(role, 1)]) ||
                    dist(0, role, sequence) >
                          dist(0, role, sequences[comp(role, 0)]) ||
                    dist(1, role, sequence) >=
                          dist(1, role, sequences[comp(role, 0)]);
            }
            if (corr(role, 1, sequence) > 0.83F) {
                bad = bad || dist(2, role, sequence) >=
                          dist(2, role, sequences[comp(role, 1)]) ||
                    dist(3, role, sequence) >
                          dist(3, role, sequences[comp(role, 1)]) ||
                    dist(2, role, sequence) >=
                          dist(2, role, sequences[comp(role, 0)]) ||
                    dist(3, role, sequence) >
                          dist(3, role, sequences[comp(role, 0)]);
            }
            if (!bad) continue;
            int bin = static_cast<int>(std::nearbyint(
                ancestor_rows[role + sequence * 3] * 1000.0F));
            if (bin >= next_no) bin = next_no;
            bins[bin] = 1;
        }
        result[role] = std::accumulate(bins.begin(), bins.end(), 0.0);
    }
    return result;
}

RdpListCorrelationState calculate_rdp_list_correlations(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<unsigned char, 3>& in,
    const std::array<unsigned char, 3>& correlation_warnings,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<float>& inversions,
    const std::vector<float>& tested_correlations,
    const std::vector<float>& f,
    const std::vector<float>& s) {
    const int count = next_no + 1;
    const auto row = [](const int role, const int sequence) {
        return role + sequence * 3;
    };
    const auto cube = [count](const int hypothesis, const int role,
                              const int sequence) {
        return hypothesis + role * 3 + sequence * 9;
    };
    std::vector<int> expected(static_cast<std::size_t>(9) * count, -1);
    const auto set = [&](const int hypothesis, const int role,
                         const int sequence, const int value) {
        expected[cube(hypothesis, role, sequence)] = value;
    };
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (f[row(in[0], sequence)] < f[row(in[0], sequences[in[1]])] &&
            s[row(in[0], sequence)] < s[row(in[0], sequences[in[1]])]) {
            set(0, in[0], sequence, 0);
        }
        if (f[row(in[0], sequence)] < f[row(in[0], sequences[in[2]])] &&
            f[row(in[0], sequence)] > 0 &&
            s[row(in[1], sequence)] < s[row(in[0], sequences[in[1]])]) {
            set(0, in[1], sequence, 0);
        }
        if (f[row(in[0], sequence)] > f[row(in[0], sequences[in[2]])] &&
            s[row(in[0], sequence)] > s[row(in[2], sequence)]) {
            set(0, in[1], sequence, 2);
        }
        if (f[row(in[0], sequence)] > f[row(in[0], sequences[in[1]])] &&
            s[row(in[0], sequence)] > s[row(in[1], sequence)]) {
            set(0, in[2], sequence, 0);
        }
        if (f[row(in[0], sequence)] < f[row(in[0], sequences[in[1]])] &&
            f[row(in[0], sequence)] > 0 &&
            s[row(in[0], sequence)] > s[row(in[1], sequence)]) {
            set(0, in[2], sequence, 2);
        }

        if (f[row(in[1], sequence)] > 0 &&
            f[row(in[1], sequence)] < f[row(in[1], sequences[in[2]])] &&
            s[row(in[1], sequence)] > s[row(in[0], sequence)]) {
            set(1, in[0], sequence, 0);
        }
        if (f[row(in[2], sequence)] > 0 &&
            f[row(in[2], sequence)] < f[row(in[1], sequences[in[2]])] &&
            s[row(in[2], sequence)] > s[row(in[1], sequences[in[2]])]) {
            set(1, in[0], sequence, 1);
        }
        if (f[row(in[1], sequence)] < f[row(in[1], sequences[in[0]])] &&
            s[row(in[1], sequence)] < s[row(in[1], sequences[in[0]])]) {
            set(1, in[1], sequence, 0);
        }
        // Preserve the source's FAMatSmall typo in the second comparison.
        if (f[row(in[2], sequence)] < f[row(in[0], sequences[in[2]])] &&
            s[row(in[2], sequence)] < f[row(in[0], sequences[in[2]])]) {
            set(1, in[2], sequence, 0);
        }
        if (f[row(in[1], sequence)] > f[row(in[1], sequences[in[0]])] &&
            f[row(in[1], sequence)] < f[row(in[1], sequences[in[2]])]) {
            set(1, in[2], sequence, 1);
        }
        if (f[row(in[1], sequence)] > 0 &&
            f[row(in[1], sequence)] < f[row(in[1], sequences[in[0]])] &&
            s[row(in[1], sequence)] < s[row(in[1], sequences[in[0]])]) {
            set(1, in[2], sequence, 4);
        }

        if (f[row(in[1], sequence)] > 0 &&
            f[row(in[1], sequence)] < f[row(in[1], sequences[in[2]])]) {
            set(2, in[0], sequence, 0);
        }
        if (f[row(in[2], sequence)] > 0 &&
            f[row(in[2], sequence)] < f[row(in[0], sequences[in[2]])] &&
            s[row(in[2], sequence)] < s[row(in[0], sequence)]) {
            set(2, in[0], sequence, 1);
        }
        if (f[row(in[1], sequence)] < f[row(in[1], sequences[in[0]])] &&
            s[row(in[1], sequence)] < s[row(in[1], sequences[in[0]])]) {
            set(2, in[1], sequence, 0);
        }
        if (f[row(in[2], sequence)] > 0 &&
            f[row(in[2], sequence)] < f[row(in[0], sequences[in[2]])] &&
            s[row(in[2], sequence)] < f[row(in[0], sequence)]) {
            set(2, in[1], sequence, 4);
        }
        if (f[row(in[1], sequence)] > f[row(in[1], sequences[in[0]])] &&
            f[row(in[1], sequence)] < f[row(in[1], sequences[in[2]])]) {
            set(2, in[1], sequence, 2);
        }
        if (f[row(in[2], sequence)] < f[row(in[0], sequences[in[2]])] &&
            s[row(in[2], sequence)] < s[row(in[0], sequences[in[2]])]) {
            set(2, in[2], sequence, 0);
        }
    }

    const auto tested = [](const int role, const int region,
                           const int category, const int sequence) {
        return role + region * 3 + category * 9 + sequence * 45;
    };
    std::vector<float> trc(tested_correlations.size(), 0.0F);
    for (int role = 0; role < 3; ++role) {
        for (int region = 0; region < 3; ++region) {
            for (int category = 0; category < 5; ++category) {
                for (int sequence = 0; sequence <= next_no; ++sequence) {
                    const auto index = tested(role, region, category, sequence);
                    const float value = tested_correlations[index];
                    if (value > 0.5F && value < 1.0F) trc[index] = value;
                }
            }
        }
    }
    for (int warned_region = 0; warned_region < 3; ++warned_region) {
        if (correlation_warnings[warned_region] == 0) continue;
        for (int role = 0; role < 3; ++role) {
            for (int category = 0; category < 5; ++category) {
                for (int sequence = 0; sequence <= next_no; ++sequence) {
                    trc[tested(role, warned_region, category, sequence)] = 0;
                }
            }
        }
    }
    std::vector<int> actual(static_cast<std::size_t>(9) * count, -1);
    for (int role = 0; role < 3; ++role) {
        for (int region = 0; region < 3; ++region) {
            for (int slot = 0; slot <= candidate_last[role]; ++slot) {
                const int sequence = candidate_list[role + slot * 3];
                actual[role + region * 3 + sequence * 9] =
                    static_cast<int>(inversions[
                        role + region * 3 + sequence * 9]);
            }
        }
    }
    std::vector<int> combined(static_cast<std::size_t>(3) * count, -1);
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            const int a = actual[role + sequence * 9];
            const int b = actual[role + 3 + sequence * 9];
            const int c = actual[role + 6 + sequence * 9];
            if (a > -1 || b > -1 || c > -1) {
                if (a == 1 || b == 1 || c == 1) combined[row(role, sequence)] = 1;
                else if (a == 2 || b == 2 || c == 2) combined[row(role, sequence)] = 2;
                else if (a == 4 || b == 4 || c == 4) combined[row(role, sequence)] = 4;
                else combined[row(role, sequence)] = 0;
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        for (int region = 0; region < 3; ++region) {
            for (int sequence = 0; sequence <= next_no; ++sequence) {
                const auto two = tested(role, region, 2, sequence);
                const auto three = tested(role, region, 3, sequence);
                const float maximum = std::max(trc[two], trc[three]);
                trc[two] = maximum;
                trc[three] = maximum;
            }
        }
    }
    RdpListCorrelationState state;
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (combined[row(role, sequence)] <= -1) continue;
            for (int hypothesis = 0; hypothesis < 3; ++hypothesis) {
                const int mapped = in[hypothesis];
                const int expected_value = expected[cube(hypothesis, role, sequence)];
                if (combined[row(role, sequence)] != expected_value) {
                    state.mismatches[mapped] += 1.0;
                }
                if (expected_value > -1 &&
                    std::abs(trc[tested(role, 0, expected_value, sequence)]) < 1 &&
                    std::abs(trc[tested(role, 1, expected_value, sequence)]) < 1 &&
                    std::abs(trc[tested(role, 2, expected_value, sequence)]) < 1) {
                    for (int region = 0; region < 3; ++region) {
                        const float value = trc[tested(
                            role, region, expected_value, sequence)];
                        if (value > 0) state.regional_strength[mapped][region] += value;
                    }
                }
            }
        }
    }
    std::array<std::array<int, 2>, 3> counts{};
    for (int role = 0; role < 3; ++role) {
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            for (int hypothesis = 0; hypothesis < 3; ++hypothesis) {
                const int mapped = in[hypothesis];
                const int expected_value = expected[cube(hypothesis, role, sequence)];
                if (expected_value > -1) {
                    for (int region = 0; region < 3; ++region) {
                        const float value = trc[tested(
                            role, region, expected_value, sequence)];
                        if (value > 0) state.expected_strength[mapped] += value;
                    }
                    counts[mapped][0] += 3;
                } else {
                    for (int region = 0; region < 3; ++region) {
                        const float value = trc[tested(role, region, 0, sequence)];
                        if (value > 0) state.absent_strength[mapped] += value;
                    }
                    counts[mapped][1] += 3;
                }
            }
        }
    }
    for (int role = 0; role < 3; ++role) {
        if (counts[role][0] == 0 || counts[role][1] == 0) {
            state.expected_strength = {};
            state.absent_strength = {};
            break;
        }
        state.expected_strength[role] /= counts[role][0];
        state.absent_strength[role] /= counts[role][1];
    }
    return state;
}

RdpFinalTrimState apply_rdp_strict_group_constraints(
    const int next_no, const std::array<int, 3>& sequences,
    const std::array<int, 6>& comparison_matrix,
    const std::vector<float>& background_matrix,
    const std::vector<float>& region_matrix,
    const std::vector<float>& f,
    const std::vector<float>& s,
    const std::vector<float>& fa,
    const std::vector<float>& sa,
    RdpFinalTrimState state) {
    const int count = next_no + 1;
    const auto matrix_cells = static_cast<std::size_t>(count) * count;
    const auto row_cells = static_cast<std::size_t>(3) * count;
    if (next_no < 1 || background_matrix.size() != matrix_cells ||
        region_matrix.size() != matrix_cells || f.size() != row_cells ||
        s.size() != row_cells || fa.size() != row_cells ||
        sa.size() != row_cells || state.candidate_list.size() != row_cells) {
        throw std::runtime_error("FinalTrim strict-group dimensions differ");
    }
    const auto row = [](const int role, const int sequence) {
        return static_cast<std::size_t>(role) + sequence * 3;
    };
    const auto matrix = [count](const std::vector<float>& values,
                                const int first, const int second) {
        return values[first + static_cast<std::size_t>(second) * count];
    };
    const auto comp = [&](const int role, const int parent) {
        return comparison_matrix[role + parent * 3];
    };
    state.strict_removed.assign(row_cells, 0);
    state.strict_readded.assign(row_cells, 0);
    if (std::getenv("RDP_TRACE_STRICT") != nullptr) {
        std::cerr << "strict-input";
        for (int role = 0; role < 3; ++role) {
            std::cerr << " role" << role << '[';
            for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                std::cerr << state.candidate_list[row(role, slot)] << ',';
            }
            std::cerr << ']';
        }
        std::cerr << '\n';
    }

    // Module2.FinalTrim 24534-24944, reached by the second RFF=1 call when
    // ConservativeGroup=0. Keep the arithmetic in Single precision: these
    // row sums feed the source's rank-based outlier constraints.
    std::vector<float> move_f(count, 0.0F);
    std::vector<float> move_s(count, 0.0F);
    for (int first = 0; first <= next_no; ++first) {
        for (int second = 0; second <= next_no; ++second) {
            move_f[first] += matrix(background_matrix, first, second);
            move_s[first] += matrix(region_matrix, first, second);
        }
    }
    std::array<int, 3> rank_f{};
    std::array<int, 3> rank_s{};
    for (int role = 0; role < 3; ++role) {
        const float selected_f = move_f[sequences[role]];
        const float selected_s = move_s[sequences[role]];
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (selected_f > move_f[sequence]) ++rank_f[role];
            if (selected_s > move_s[sequence]) ++rank_s[role];
        }
    }

    for (int role = 0; role < 3; ++role) {
        std::vector<unsigned char> remove(count, 0);
        std::array<int, 2> ilp{};
        const int seq1 = sequences[role];
        const int seq2_role = comp(role, 0);
        const int seq3_role = comp(role, 1);
        const int seq2 = sequences[seq2_role];
        const int seq3 = sequences[seq3_role];
        int go_on = 1;

        if (sa[row(role, seq2)] > sa[row(seq2_role, seq3)] ||
            (s[row(role, seq2)] > s[row(seq2_role, seq3)] &&
             s[row(role, seq3)] > s[row(seq2_role, seq3)])) {
            ilp[0] = 1;
            const float bound = sa[row(seq2_role, seq3)];
            for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                const int seq4 = state.candidate_list[row(role, slot)];
                if (seq4 != seq1 && s[row(role, seq4)] > 0.0F) {
                    if (sa[row(role, seq2)] != sa[row(seq2_role, seq4)] ||
                        sa[row(role, seq3)] != sa[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                    if (sa[row(role, seq4)] > bound) remove[seq4] = 1;
                    // The source repeats the identical SAMatSmall test here.
                    if (sa[row(role, seq4)] > bound) remove[seq4] = 1;
                    if (s[row(role, seq4)] > s[row(seq2_role, seq4)] ||
                        s[row(role, seq4)] > s[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            }
            if (fa[row(role, seq2)] < fa[row(role, seq3)]) {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 != seq1 &&
                        fa[row(role, seq2)] != fa[row(seq2_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            } else {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 != seq1 &&
                        fa[row(role, seq3)] != fa[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            }
            go_on = 0;
        }

        if (fa[row(role, seq2)] > fa[row(seq2_role, seq3)] ||
            (f[row(role, seq2)] > f[row(seq2_role, seq3)] &&
             f[row(role, seq3)] > f[row(seq2_role, seq3)])) {
            ilp[1] = 1;
            const float bound = fa[row(seq2_role, seq3)];
            for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                const int seq4 = state.candidate_list[row(role, slot)];
                if (seq4 != seq1 && f[row(role, seq4)] > 0.0F) {
                    if (fa[row(role, seq2)] != fa[row(seq2_role, seq4)] ||
                        fa[row(role, seq3)] != fa[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                    if (fa[row(role, seq4)] > bound) remove[seq4] = 1;
                    if (f[row(role, seq4)] > f[row(seq2_role, seq4)] ||
                        f[row(role, seq4)] > f[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            }
            if (sa[row(role, seq2)] < sa[row(role, seq3)]) {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 != seq1 &&
                        sa[row(role, seq2)] != sa[row(seq2_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            } else {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 != seq1 &&
                        sa[row(role, seq3)] != sa[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
            }
            go_on = 0;
        }

        if (go_on == 1) {
            if (sa[row(role, seq2)] < sa[row(role, seq3)]) {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 == seq1) continue;
                    if (sa[row(role, seq4)] >= sa[row(seq3_role, seq4)] ||
                        s[row(role, seq4)] > s[row(role, seq2)] * 6.0F) {
                        remove[seq4] = 1;
                    }
                    if (sa[row(role, seq4)] > sa[row(role, seq2)] &&
                        s[row(role, seq2)] < f[row(role, seq2)] &&
                        s[row(seq2_role, seq4)] > f[row(seq2_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
                if ((static_cast<double>(rank_s[seq3_role]) / next_no > 0.95 &&
                     static_cast<double>(rank_f[seq3_role]) / next_no < 0.75) ||
                    static_cast<double>(rank_s[seq3_role] - rank_f[seq3_role]) /
                            next_no > 0.5) {
                    for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                        const int seq4 = state.candidate_list[row(role, slot)];
                        if (seq4 != seq1 &&
                            sa[row(seq2_role, seq4)] > sa[row(role, seq2)]) {
                            remove[seq4] = 1;
                        }
                    }
                }
            } else {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 == seq1) continue;
                    if (sa[row(role, seq4)] >= sa[row(seq2_role, seq4)] ||
                        s[row(role, seq4)] > s[row(role, seq3)] * 6.0F) {
                        remove[seq4] = 1;
                    }
                    if (sa[row(role, seq4)] > sa[row(role, seq3)] &&
                        s[row(role, seq3)] < f[row(role, seq3)] &&
                        s[row(seq3_role, seq4)] > f[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
                if ((static_cast<double>(rank_s[seq2_role]) / next_no > 0.95 &&
                     static_cast<double>(rank_f[seq2_role]) / next_no < 0.75) ||
                    static_cast<double>(rank_s[seq2_role] - rank_f[seq2_role]) /
                            next_no > 0.5) {
                    for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                        const int seq4 = state.candidate_list[row(role, slot)];
                        if (seq4 != seq1 &&
                            sa[row(seq3_role, seq4)] > sa[row(role, seq3)]) {
                            remove[seq4] = 1;
                        }
                    }
                }
            }

            if (fa[row(role, seq2)] < fa[row(role, seq3)]) {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 == seq1) continue;
                    if (fa[row(role, seq4)] >= fa[row(seq3_role, seq4)] ||
                        f[row(role, seq4)] > f[row(role, seq2)] * 6.0F) {
                        remove[seq4] = 1;
                    }
                    if (fa[row(role, seq4)] > fa[row(role, seq2)] &&
                        f[row(role, seq2)] < s[row(role, seq2)] &&
                        f[row(seq2_role, seq4)] > s[row(seq2_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
                if ((static_cast<double>(rank_f[seq3_role]) / next_no > 0.95 &&
                     static_cast<double>(rank_s[seq3_role]) / next_no < 0.75) ||
                    static_cast<double>(rank_f[seq3_role] - rank_s[seq3_role]) /
                            next_no > 0.5) {
                    for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                        const int seq4 = state.candidate_list[row(role, slot)];
                        if (seq4 != seq1 &&
                            fa[row(seq2_role, seq4)] > fa[row(role, seq2)]) {
                            remove[seq4] = 1;
                        }
                    }
                }
            } else {
                for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                    const int seq4 = state.candidate_list[row(role, slot)];
                    if (seq4 == seq1) continue;
                    if (fa[row(role, seq4)] >= fa[row(seq2_role, seq4)] ||
                        f[row(role, seq4)] > f[row(role, seq3)] * 6.0F) {
                        remove[seq4] = 1;
                    }
                    if (fa[row(role, seq4)] > fa[row(role, seq3)] &&
                        f[row(role, seq3)] < s[row(role, seq3)] &&
                        f[row(seq3_role, seq4)] > s[row(seq3_role, seq4)]) {
                        remove[seq4] = 1;
                    }
                }
                if ((static_cast<double>(rank_f[seq2_role]) / next_no > 0.95 &&
                     static_cast<double>(rank_s[seq2_role]) / next_no < 0.75) ||
                    static_cast<double>(rank_f[seq2_role] - rank_s[seq2_role]) /
                            next_no > 0.5) {
                    for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                        const int seq4 = state.candidate_list[row(role, slot)];
                        if (seq4 != seq1 &&
                            fa[row(seq3_role, seq4)] > fa[row(role, seq3)]) {
                            remove[seq4] = 1;
                        }
                    }
                }
            }
        }

        // The VB loop removes in place by replacing a marked entry with the
        // current tail, then immediately retesting that slot.
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            state.strict_removed[row(role, sequence)] = remove[sequence];
        }
        int slot = 0;
        while (slot <= state.candidate_last[role]) {
            const int sequence = state.candidate_list[row(role, slot)];
            if (remove[sequence] != 0) {
                if (slot < state.candidate_last[role]) {
                    state.candidate_list[row(role, slot)] =
                        state.candidate_list[row(role,
                                                  state.candidate_last[role])];
                }
                --state.candidate_last[role];
            } else {
                ++slot;
            }
        }

        float highest_s = 0.0F;
        float highest_f = 0.0F;
        float highest_sa = 0.0F;
        float highest_fa = 0.0F;
        for (int item = 0; item <= state.candidate_last[role]; ++item) {
            const int sequence = state.candidate_list[row(role, item)];
            highest_sa = std::max(highest_sa, sa[row(role, sequence)]);
            highest_fa = std::max(highest_fa, fa[row(role, sequence)]);
            highest_s = std::max(highest_s, s[row(role, sequence)]);
            highest_f = std::max(highest_f, f[row(role, sequence)]);
        }
        const auto contains = [&](const int sequence) {
            for (int item = 0; item <= state.candidate_last[role]; ++item) {
                if (state.candidate_list[row(role, item)] == sequence) {
                    return true;
                }
            }
            return false;
        };
        const auto append = [&](const int sequence) {
            if (!contains(sequence)) {
                ++state.candidate_last[role];
                state.candidate_list[row(role,
                                          state.candidate_last[role])] = sequence;
                state.strict_readded[row(role, sequence)] = 1;
            }
        };
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            go_on = 0;
            // Preserve the source typo: the fourth bound is FAMatSmall <= HDF,
            // not FMatSmall <= HDF.
            if (sa[row(role, sequence)] <= highest_sa &&
                fa[row(role, sequence)] <= highest_fa &&
                s[row(role, sequence)] <= highest_s &&
                fa[row(role, sequence)] <= highest_f) {
                if (sa[row(role, sequence)] < sa[row(seq2_role, sequence)] &&
                    sa[row(role, sequence)] < sa[row(seq3_role, sequence)] &&
                    fa[row(role, sequence)] < fa[row(seq2_role, sequence)] &&
                    fa[row(role, sequence)] < fa[row(seq3_role, sequence)]) {
                    append(sequence);
                    go_on = 1;
                }
            }
            // `... Or x = x` in the original makes the ILP guard tautological.
            if (go_on == 0) {
                if (sa[row(role, sequence)] < sa[row(seq2_role, sequence)] &&
                    sa[row(role, sequence)] < sa[row(seq3_role, sequence)] &&
                    s[row(role, sequence)] < s[row(seq2_role, sequence)] &&
                    s[row(role, sequence)] < s[row(seq3_role, sequence)] &&
                    fa[row(role, sequence)] < fa[row(seq2_role, sequence)] &&
                    fa[row(role, sequence)] < fa[row(seq3_role, sequence)] &&
                    f[row(role, sequence)] < f[row(seq2_role, sequence)] &&
                    f[row(role, sequence)] < f[row(seq3_role, sequence)]) {
                    append(sequence);
                    go_on = 1;
                }
            }
            if (go_on == 1) {
                state.synthetic_event_roles.push_back({role, sequence});
            }
        }

        (void)ilp;
    }
    if (std::getenv("RDP_TRACE_STRICT") != nullptr) {
        std::cerr << "strict-output";
        for (int role = 0; role < 3; ++role) {
            std::cerr << " role" << role << '[';
            for (int slot = 0; slot <= state.candidate_last[role]; ++slot) {
                std::cerr << state.candidate_list[row(role, slot)] << ',';
            }
            std::cerr << ']';
        }
        std::cerr << '\n';
    }
    return state;
}

std::vector<int> make_rdp_relevant_sequences(
    const int next_no, const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list) {
    const int count = next_no + 1;
    if (next_no < 0 || candidate_list.size() !=
            static_cast<std::size_t>(3) * count) {
        throw std::runtime_error("MakeRelevant dimensions differ");
    }
    std::vector<int> relevant(count, 0);
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= candidate_last[role]; ++slot) {
            relevant[candidate_list[role + slot * 3]] = 1;
        }
    }
    return relevant;
}

RdpRawEventState prepare_rdp_collection_event_list(
    const int next_no, const int winning_role,
    const std::array<int, 3>& sequences,
    const std::array<int, 2>& trace,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const std::vector<double>& acceptable_sequences,
    const RdpRawEventState& events) {
    const int count = next_no + 1;
    if (winning_role < 0 || winning_role > 2 ||
        trace[0] < 0 || trace[0] > next_no || trace[1] < 1 ||
        static_cast<std::size_t>(trace[0]) >= events.xover_list.size() ||
        static_cast<std::size_t>(trace[1]) >
            events.xover_list[trace[0]].size() ||
        candidate_list.size() != static_cast<std::size_t>(3) * count ||
        acceptable_sequences.size() !=
            static_cast<std::size_t>(3) * 19 * count) {
        throw std::runtime_error("collection event-list dimensions differ");
    }
    RdpRawEventState prepared = events;
    const RdpRawEvent selected = events.xover_list[trace[0]][trace[1] - 1];
    for (int slot = 0; slot <= candidate_last[winning_role]; ++slot) {
        const int sequence = candidate_list[winning_role + slot * 3];
        const auto ok_index = static_cast<std::size_t>(winning_role) +
            18U * 3U + static_cast<std::size_t>(sequence) * 57U;
        if (sequence == sequences[winning_role] ||
            acceptable_sequences[ok_index] <= 1.0) {
            continue;
        }
        RdpRawEvent copy = selected;
        // MakeCollecteventsC receives the persistent PXOList record.  The
        // candidate copy therefore retains its original probability; the
        // source only assigns 1.0 later when it builds the temporary rescan
        // list, not in the collection-event input itself.
        if (copy.daughter == sequences[winning_role]) {
            copy.daughter = static_cast<std::int16_t>(sequence);
        } else if (copy.major_parent == sequences[winning_role]) {
            copy.major_parent = static_cast<std::int16_t>(sequence);
        } else if (copy.minor_parent == sequences[winning_role]) {
            copy.minor_parent = static_cast<std::int16_t>(sequence);
        }
        copy.distance_holder = std::abs(copy.distance_holder);
        prepared.xover_list[sequence].push_back(copy);
        prepared.current_xover[sequence] = static_cast<std::int16_t>(
            prepared.xover_list[sequence].size());
    }

    return prepared;
}

RdpCollectedEventsState make_rdp_parent_collect_events(
    const int sequence_length, const int next_no, const int role,
    std::array<int, 6> region_sizes,
    const std::vector<int>& overlap_sequence,
    const std::array<int, 6>& comparison_matrix,
    const std::array<int, 3>& candidate_last,
    const std::vector<int>& candidate_list,
    const int add_num, const std::array<int, 3>& sequences,
    const std::array<int, 2>& trace,
    const RdpRawEventState& source_events) {
    const int count = next_no + 1;
    RdpCollectedEventsState state;
    state.events.resize(static_cast<std::size_t>(count) * (add_num + 1));
    if (candidate_last[role] == -1) return state;
    if (role < 0 || role > 2 || overlap_sequence.size() !=
            static_cast<std::size_t>(sequence_length + 1) ||
        candidate_list.size() != static_cast<std::size_t>(3) * count ||
        source_events.current_xover.size() != static_cast<std::size_t>(count) ||
        source_events.xover_list.size() != static_cast<std::size_t>(count)) {
        throw std::runtime_error("MakeCollecteventsC dimensions differ");
    }
    const int role_count = candidate_last[role] + 1;
    std::vector<double> best_probability(
        static_cast<std::size_t>(role_count) * (add_num * 2 + 1), 0.0);
    for (int slot = 0; slot < role_count; ++slot) {
        for (int method = 0; method <= add_num; ++method) {
            best_probability[slot + method * role_count] = 1.0;
        }
    }
    std::vector<unsigned char> role_membership(
        static_cast<std::size_t>(3) * count, 0);
    const int comparison0 = comparison_matrix[role];
    const int comparison1 = comparison_matrix[role + 3];
    for (const int member_role : {role, comparison0, comparison1}) {
        for (int slot = 0; slot <= candidate_last[member_role]; ++slot) {
            const int sequence = candidate_list[member_role + slot * 3];
            role_membership[member_role + sequence * 3] = 1;
        }
    }
    const auto in_any_role = [&](const int sequence) {
        return role_membership[role + sequence * 3] == 1 ||
            role_membership[comparison0 + sequence * 3] == 1 ||
            role_membership[comparison1 + sequence * 3] == 1;
    };
    const auto role_has_event_member = [&](const int member_role,
                                           const RdpRawEvent& event) {
        return role_membership[member_role + event.daughter * 3] == 1 ||
            role_membership[member_role + event.major_parent * 3] == 1 ||
            role_membership[member_role + event.minor_parent * 3] == 1;
    };
    for (int row = 0; row <= next_no; ++row) {
        for (const auto& event : source_events.xover_list[row]) {
            if (!in_any_role(event.daughter) || event.major_parent > next_no ||
                !in_any_role(event.major_parent) ||
                event.minor_parent > next_no ||
                !in_any_role(event.minor_parent)) {
                continue;
            }
            const bool bad_match = !role_has_event_member(role, event) ||
                !role_has_event_member(comparison0, event) ||
                !role_has_event_member(comparison1, event);
            int overlap_size = 0;
            if (event.beginning < event.ending) {
                region_sizes[1] = event.ending - event.beginning + 1;
                for (int position = event.beginning;
                     position <= event.ending; ++position) {
                    overlap_size += overlap_sequence[position];
                }
            } else {
                region_sizes[1] = event.ending + sequence_length -
                    event.beginning + 1;
                for (int position = 1; position <= event.ending; ++position) {
                    overlap_size += overlap_sequence[position];
                }
                for (int position = event.beginning;
                     position <= sequence_length; ++position) {
                    overlap_size += overlap_sequence[position];
                }
            }
            float match = 0.0F;
            if (overlap_size > 0) {
                const float denominator = static_cast<float>(
                    region_sizes[0] + region_sizes[1]);
                match = (static_cast<float>(overlap_size) * 2.0F) /
                    denominator;
            }
            const double rounded_match =
                std::round(static_cast<double>(match) * 100000.0) / 100000.0;
            if (!((rounded_match > 0.1 && !bad_match) ||
                  (rounded_match > 0.5 && bad_match))) {
                continue;
            }
            int collection_slot = -1;
            for (int slot = 0; slot < role_count; ++slot) {
                const int candidate = candidate_list[role + slot * 3];
                if (candidate == event.daughter ||
                    candidate == event.major_parent ||
                    candidate == event.minor_parent) {
                    collection_slot = slot;
                    break;
                }
            }
            const int method = event.program_flag;
            if (collection_slot >= 0 && method >= 0 && method <= add_num &&
                (best_probability[collection_slot + method * role_count] >
                     event.probability ||
                 best_probability[collection_slot + method * role_count] ==
                     0.0)) {
                best_probability[collection_slot + method * role_count] =
                    event.probability;
                state.events[collection_slot + method * count] = event;
            }
        }
    }
    int selected_slot = 0;
    for (int slot = 0; slot < role_count; ++slot) {
        const int candidate = candidate_list[role + slot * 3];
        if (candidate == sequences[0] || candidate == sequences[1] ||
            candidate == sequences[2]) {
            selected_slot = slot;
            break;
        }
    }
    const auto& selected = source_events.xover_list[trace[0]][trace[1] - 1];
    state.events[selected_slot + selected.program_flag * count] = selected;
    state.result = 1;
    return state;
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
    if (std::getenv("RDP_TRACE_PHPR") != nullptr) {
        std::uint64_t hash_f = 1469598103934665603ULL;
        std::uint64_t hash_s = 1469598103934665603ULL;
        const auto mix = [](std::uint64_t& hash, const void* data,
                            const std::size_t bytes) {
            const auto* p = static_cast<const unsigned char*>(data);
            for (std::size_t i = 0; i < bytes; ++i) {
                hash ^= p[i];
                hash *= 1099511628211ULL;
            }
        };
        mix(hash_f, background_matrix.data(),
            background_matrix.size() * sizeof(float));
        mix(hash_s, region_matrix.data(), region_matrix.size() * sizeof(float));
        std::cerr << "phpr next=" << next_no << " seqs=" << sequences[0]
                  << ':' << sequences[1] << ':' << sequences[2]
                  << " done=";
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            if (done_this[sequence * 2] || done_this[sequence * 2 + 1]) {
                std::cerr << sequence << ':' << done_this[sequence * 2] << ':'
                          << done_this[sequence * 2 + 1] << ',';
            }
        }
        std::cerr << " hash=" << hash_f << ':' << hash_s << " samples="
                  << background_matrix[17 + 20 * (next_no + 1)] << ':'
                  << background_matrix[17 + 3 * (next_no + 1)] << ':'
                  << region_matrix[17 + 20 * (next_no + 1)] << ':'
                  << region_matrix[17 + 3 * (next_no + 1)] << " scores="
                  << output.scores[0]
                  << ':' << output.scores[1] << ':' << output.scores[2]
                  << " sub=" << output.sub_scores[0] << ':'
                  << output.sub_scores[1] << ':' << output.sub_scores[2]
                  << '\n';
    }
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

#include "tree_state.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

void make_distance_backup(
    const int next_no, const std::vector<float>& distance,
    std::vector<float>& backup) {
    const int stride = next_no + 1;
    for (int first = 0; first <= next_no; ++first) {
        for (int second = first + 1; second <= next_no; ++second) {
            if (distance[first + second * stride] > 0.0F) {
                backup[first + second * stride] =
                    1.0F - distance[first + second * stride];
                backup[second + first * stride] =
                    backup[first + second * stride];
            } else {
                backup[first + second * stride] = 0.999F;
                backup[second + first * stride] = 0.999F;
            }
        }
    }
}

void make_distance_map(
    const int next_no, const std::vector<float>& distance,
    std::vector<float>& map, std::vector<int>& winning) {
    const int stride = next_no + 1;
    std::fill(map.begin(), map.end(), 100.0F);
    for (int first = 0; first <= next_no; ++first) {
        for (int second = first + 1; second <= next_no; ++second) {
            const float value = distance[first + second * stride];
            if (map[first] > value) {
                map[first] = value;
                winning[first] = second;
            }
            if (map[second] > value) {
                map[second] = value;
                winning[second] = first;
            }
        }
    }
}

float shortest_distance(
    const int next_no, const int sorted_no, const std::vector<float>& map,
    const std::vector<int>& winning, const std::vector<float>& distance,
    std::vector<short>& tree_x, std::vector<short>& tree_y) {
    float minimum = 1.0F;
    int winner = 0;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (map[sequence] < minimum) {
            minimum = map[sequence];
            winner = sequence;
        }
    }
    if (minimum < 1.0F) {
        const float result =
            distance[winning[winner] + winner * (next_no + 1)];
        if (winner < winning[winner]) {
            tree_x[sorted_no] = static_cast<short>(winner);
            tree_y[sorted_no] = static_cast<short>(winning[winner]);
        } else {
            tree_x[sorted_no] = static_cast<short>(winning[winner]);
            tree_y[sorted_no] = static_cast<short>(winner);
        }
        return result;
    }
    return minimum;
}

void add_sequence_to_upgma(
    const int next_no, const int sorted_no, std::vector<short>& instances,
    std::vector<float>& distance, std::vector<short>& tree_x,
    const std::vector<short>& tree_y, std::vector<int>& check,
    std::vector<short>& node_y_position) {
    const int stride = next_no + 1;
    const int first = tree_x[sorted_no];
    const int second = tree_y[sorted_no];
    const int first_offset = first * stride;
    const int second_offset = second * stride;
    int multiple_include = 0;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (sequence != first && sequence != second) {
            if (distance[sequence + first_offset] < 1.0F &&
                distance[sequence + second_offset] < 1.0F) {
                distance[sequence + first_offset] =
                    (distance[sequence + first_offset] +
                     distance[sequence + second_offset]) /
                    2.0F;
            } else if (distance[sequence + second_offset] < 1.0F) {
                distance[sequence + first_offset] =
                    distance[sequence + second_offset];
            } else {
                distance[sequence + first_offset] = 0.999999F;
            }
            distance[first + sequence * stride] =
                distance[sequence + first_offset];
        } else if (sequence == first) {
            if (instances[first] > 0) {
                for (int row = sorted_no - 1; row >= 0; --row) {
                    if (tree_x[row] == first) {
                        multiple_include = row;
                        break;
                    }
                }
                const int count = instances[first];
                for (int member = 0; member <= count; ++member) {
                    tree_x[sorted_no + member * stride] =
                        tree_x[multiple_include + member * stride];
                }
            }
            if (instances[second] > 0) {
                for (int row = sorted_no - 1; row >= 0; --row) {
                    if (tree_x[row] == second) {
                        multiple_include = row;
                        break;
                    }
                }
                const int count = instances[second];
                for (int member = 0; member <= count; ++member) {
                    ++instances[first];
                    if (instances[first] > next_no) {
                        std::fill(check.begin(), check.end(), 0);
                        for (int slot = 0; slot <= next_no; ++slot) {
                            ++check[tree_x[sorted_no + slot * stride]];
                        }
                        int replacement = 0;
                        for (int slot = 0; slot <= next_no; ++slot) {
                            if (check[slot] > 1) replacement = slot;
                        }
                        instances[first] = static_cast<short>(replacement);
                    }
                    tree_x[sorted_no + instances[first] * stride] =
                        tree_x[multiple_include + member * stride];
                }
            } else {
                ++instances[first];
                if (instances[first] > next_no) {
                    std::fill(check.begin(), check.end(), 0);
                    for (int slot = 0; slot <= next_no; ++slot) {
                        ++check[tree_x[sorted_no + slot * stride]];
                    }
                    int replacement = 0;
                    for (int slot = 0; slot <= next_no; ++slot) {
                        if (check[slot] > 1) replacement = slot;
                    }
                    instances[first] = static_cast<short>(replacement);
                }
                tree_x[sorted_no + instances[first] * stride] =
                    static_cast<short>(second);
            }
        }
        node_y_position[sequence] = instances[sequence];
    }
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        if (sequence != second) {
            distance[second + sequence * stride] = 100.0F;
            distance[sequence + second * stride] = 100.0F;
        }
    }
}

void update_distance_map(
    const int first, const int second, const int next_no,
    std::vector<float>& map, const std::vector<float>& distance,
    std::vector<int>& winning) {
    const int stride = next_no + 1;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        const int offset = sequence * stride;
        if (first != sequence && second != sequence) {
            if (map[first] > distance[first + offset]) {
                map[first] = distance[first + offset];
                winning[first] = sequence;
            }
            if (winning[sequence] == first || winning[sequence] == second) {
                if (map[sequence] < 100.0F) {
                    map[sequence] = 100.0F;
                    for (int candidate = 0; candidate <= next_no;
                         ++candidate) {
                        if (candidate != sequence &&
                            map[sequence] > distance[candidate + offset]) {
                            map[sequence] = distance[candidate + offset];
                            winning[sequence] = candidate;
                        }
                    }
                }
            }
        }
    }
}

void finish_tree_distances(
    const int next_no, const std::vector<short>& tree_x,
    const std::vector<short>& tree_y, std::vector<float>& tree_distance) {
    const int stride = next_no + 1;
    int member = 0;
    for (int row = next_no - 1; row >= 0; --row) {
        const int second = tree_y[row];
        int slot = 1;
        for (; slot <= next_no; ++slot) {
            member = tree_x[row + slot * stride];
            if (member > -1) {
                if (member == second) break;
                tree_distance[member + second * stride] =
                    tree_distance[tree_x[row] + second * stride];
                tree_distance[second + member * stride] =
                    tree_distance[member + second * stride];
            } else {
                break;
            }
        }
        if (member > -1) {
            for (int later_slot = slot + 1; later_slot <= next_no;
                 ++later_slot) {
                const int later = tree_x[row + later_slot * stride];
                if (later > -1) {
                    for (int earlier_slot = 0; earlier_slot < slot;
                         ++earlier_slot) {
                        const int earlier =
                            tree_x[row + earlier_slot * stride];
                        tree_distance[later + earlier * stride] =
                            tree_distance[tree_x[row] + second * stride];
                        tree_distance[earlier + later * stride] =
                            tree_distance[later + earlier * stride];
                    }
                } else {
                    break;
                }
            }
        }
    }
}

}  // namespace

RdpTreeState build_rdp_upgma_tree_state(
    const int next_no, const RdpDistanceState& distance_state) {
    const int sequence_count = next_no + 1;
    const auto matrix_size =
        static_cast<std::size_t>(sequence_count) * sequence_count;
    if (next_no < 1 || distance_state.distance.size() != matrix_size) {
        throw std::runtime_error("UPGMA received a mis-sized distance matrix");
    }

    RdpTreeState result;
    result.tree_distance.assign(matrix_size, 0.0F);
    result.tree_x.assign(matrix_size, -1);
    result.tree_y.assign(sequence_count, 0);
    result.node_length.assign(sequence_count, 0.0);
    std::vector<float> distance_backup(matrix_size, 0.0F);
    std::vector<short> instances(sequence_count + 1, 0);
    std::vector<int> check(sequence_count, 0);
    std::vector<short> node_y_position(sequence_count, 0);
    std::vector<float> distance_map(sequence_count, 0.0F);
    std::vector<int> winning(sequence_count, 0);

    make_distance_backup(next_no, distance_state.distance, distance_backup);
    make_distance_map(next_no, distance_backup, distance_map, winning);

    for (int sorted = 0; sorted < next_no; ++sorted) {
        float shortest = shortest_distance(
            next_no, sorted, distance_map, winning, distance_backup,
            result.tree_x, result.tree_y);
        if (shortest == 1.0F) {
            for (int first = 0; first <= next_no; ++first) {
                for (int second = first + 1; second <= next_no; ++second) {
                    if (distance_backup[first + second * sequence_count] <
                        shortest) {
                        shortest =
                            distance_backup[first + second * sequence_count];
                        result.tree_x[sorted] = static_cast<short>(first);
                        result.tree_y[sorted] = static_cast<short>(second);
                    }
                }
            }
        }
        result.node_length[sorted] =
            static_cast<int>((shortest / 2.0F) * 100000.0F) / 100000.0;
        if (shortest < 0.999999F) {
            const int first = result.tree_x[sorted];
            const int second = result.tree_y[sorted];
            result.tree_distance[first + second * sequence_count] =
                1.0F - shortest / 2.0F;
            result.tree_distance[second + first * sequence_count] =
                1.0F - shortest / 2.0F;
        }
        add_sequence_to_upgma(
            next_no, sorted, instances, distance_backup, result.tree_x,
            result.tree_y, check, node_y_position);
        const int first = result.tree_x[sorted];
        const int second = result.tree_y[sorted];
        distance_map[first] = 100.0F;
        distance_map[second] = 100.0F;
        update_distance_map(
            first, second, next_no, distance_map, distance_backup, winning);
    }

    finish_tree_distances(
        next_no, result.tree_x, result.tree_y, result.tree_distance);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        result.tree_distance[sequence + sequence * sequence_count] = 1.0F;
    }
    return result;
}

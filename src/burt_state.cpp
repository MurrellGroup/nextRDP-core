#include "burt_state.hpp"

#include "MathFuncsDll.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kNegativeLattice = -10000000000000000.0;

int source_clng(const double value) {
    // VB6 CLng uses banker's rounding.  The default IEC rounding mode used by
    // nearbyint is ties-to-even, which is the source behavior for positives.
    return static_cast<int>(std::nearbyint(value));
}

int bounded_coordinate(const int value, const int length) {
    if (length <= 0) return 0;
    const long long absolute = value < 0 ? -static_cast<long long>(value) : value;
    if (absolute < 1) return 1;
    if (absolute > length) return length;
    return static_cast<int>(absolute);
}

int xdiff_at(const std::vector<int>& xdiff, int coordinate, int length) {
    coordinate = bounded_coordinate(coordinate, length);
    if (coordinate < 0 || coordinate >= static_cast<int>(xdiff.size())) return 0;
    return xdiff[static_cast<std::size_t>(coordinate)];
}

struct HmmInput {
    std::vector<unsigned char> recode;
    std::vector<int> xdiff;
    std::vector<int> xpos;
    int informative = 0;
    int sequence_length = 0;
    int slen = 0;
    int circular_offset = 0;
};

HmmInput make_hmm_input(
    const RdpScanState& scan_state, std::array<int, 3> representatives,
    bool circular) {
    HmmInput input;
    input.sequence_length = scan_state.sequence_length;
    std::sort(representatives.begin(), representatives.end());
    const int length = scan_state.sequence_length;
    const int stride = length + 1;
    input.recode.assign(static_cast<std::size_t>(length * 2 + 4), 0);
    input.xdiff.assign(static_cast<std::size_t>(length * 2 + 4), 0);
    input.xpos.assign(static_cast<std::size_t>(length * 2 + 4), 0);
    for (int coordinate = 1; coordinate <= length; ++coordinate) {
        const short first = scan_state.sequence_data[
            coordinate + representatives[0] * stride];
        const short second = scan_state.sequence_data[
            coordinate + representatives[1] * stride];
        const short third = scan_state.sequence_data[
            coordinate + representatives[2] * stride];
        if (first != 46 && second != 46 && third != 46) {
            unsigned char symbol = 0;
            bool informative = false;
            if (first == second && first != third) {
                symbol = 0;
                informative = true;
            } else if (first == third && first != second) {
                symbol = 2;
                informative = true;
            } else if (second == third && first != second) {
                symbol = 1;
                informative = true;
            }
            if (informative) {
                input.recode[static_cast<std::size_t>(input.informative)] = symbol;
                ++input.informative;
                input.xdiff[static_cast<std::size_t>(input.informative)] = coordinate;
            }
        }
        input.xpos[static_cast<std::size_t>(coordinate)] = input.informative;
    }
    if (!circular) {
        input.slen = input.informative;
        return input;
    }

    input.circular_offset = source_clng(input.informative / 2.0);
    const int expanded = input.informative + input.circular_offset * 2;
    std::vector<unsigned char> expanded_recode(
        static_cast<std::size_t>(expanded + 1), 0);
    for (int index = 0; index < input.circular_offset; ++index) {
        expanded_recode[static_cast<std::size_t>(
            input.informative + index + input.circular_offset + 1)] =
            input.recode[static_cast<std::size_t>(index)];
    }
    for (int index = input.circular_offset; index < input.informative; ++index) {
        expanded_recode[static_cast<std::size_t>(index - input.circular_offset)] =
            input.recode[static_cast<std::size_t>(index)];
    }
    for (int index = 0; index < input.informative; ++index) {
        expanded_recode[static_cast<std::size_t>(index + input.circular_offset + 1)] =
            input.recode[static_cast<std::size_t>(index)];
    }
    input.recode = std::move(expanded_recode);
    input.slen = expanded;
    return input;
}

struct HmmOutput {
    std::vector<double> posterior;
    std::vector<int> path;
    double best_like = kNegativeLattice;
    bool trained = false;
};

HmmOutput train_hmm(const HmmInput& input, int hmm_cycles) {
    HmmOutput output;
    if (input.informative <= 0 || input.slen <= 0) return output;
    constexpr int states = 3;
    constexpr int symbols = 3;
    const int stride = input.slen + 1;
    std::vector<double> transition(states * states, 0.0);
    std::vector<double> emission(symbols * states, 0.0);
    std::vector<double> initial(states, 0.0);
    std::vector<int> path(static_cast<std::size_t>(input.slen + 4), 0);
    const float best = MathFuncs::MyMathFuncs::DoHMMCyclesSerial(
        3, input.slen, hmm_cycles, input.sequence_length, states, symbols,
        const_cast<unsigned char*>(input.recode.data()), transition.data(),
        emission.data(), initial.data(), path.data());
    output.best_like = static_cast<double>(best);
    output.path = path;
    output.trained = best > -1000000.0;
    if (!output.trained) return output;

    std::vector<double> forward(static_cast<std::size_t>(stride * states), 0.0);
    std::vector<double> reverse(static_cast<std::size_t>(stride * states), 0.0);
    std::vector<double> values(states, 0.0);
    std::vector<double> options(states * states, 0.0);
    for (int state = 0; state < states; ++state) {
        forward[static_cast<std::size_t>(state * stride)] =
            emission[input.recode[0] + state * symbols] + initial[state];
    }
    MathFuncs::MyMathFuncs::ForwardCP(
        input.slen, symbols, states, values.data(), options.data(),
        const_cast<unsigned char*>(input.recode.data()), forward.data(),
        transition.data(), emission.data());
    for (int state = 0; state < states; ++state) {
        reverse[static_cast<std::size_t>(input.slen + state * stride)] = 0.0;
    }
    MathFuncs::MyMathFuncs::ReverseCP(
        input.slen, symbols, states, values.data(), options.data(),
        const_cast<unsigned char*>(input.recode.data()), reverse.data(),
        transition.data(), emission.data());
    output.posterior.assign(static_cast<std::size_t>(stride * states), 0.0);
    for (int position = 0; position <= input.slen; ++position) {
        double maximum = -1.0e22;
        for (int state = 0; state < states; ++state) {
            const double combined =
                forward[position + state * stride] +
                reverse[position + state * stride];
            output.posterior[position * states + state] = combined;
            maximum = std::max(maximum, combined);
        }
        double total = 0.0;
        for (int state = 0; state < states; ++state) {
            double& value = output.posterior[position * states + state];
            value = std::exp(value - maximum);
            total += value;
        }
        for (int state = 0; state < states; ++state) {
            output.posterior[position * states + state] /= total;
        }
    }

    return output;
}

// We keep this as a separate helper so the source's one-based CIs and its
// CurCI=0 first-interval behavior remain visible and testable.
std::array<int, 5> match_bp_to_ci(
    std::vector<RdpBurtInterval>& intervals,
    int cur_ci,
    int breakpoint,
    int length,
    const std::vector<int>& xpos,
    bool circular) {
    std::array<int, 5> result{};
    if (cur_ci < 1) {
        // MatchBPtoCI still falls through to its SmallestM branch when
        // CurCI=0.  Since its one-based loop does no work, M remains zero and
        // the routine returns the negated first interval (rather than an
        // empty result).  This is an observable VB quirk.
        if (!intervals.empty()) {
            result = intervals.front().values;
            for (int& value : result) value = -value;
        }
        return result;
    }
    std::vector<int> match(static_cast<std::size_t>(cur_ci + 1), 0);
    for (int x = 1; x <= cur_ci && x < static_cast<int>(intervals.size()); ++x) {
        auto& ci = intervals[static_cast<std::size_t>(x)].values;
        const int factor = length > 0 ? ci[0] / length : 0;
        if (factor > 0) {
            for (int& value : ci) {
                if (value > length) value -= length * factor;
            }
        }
        std::vector<unsigned char> range(static_cast<std::size_t>(length + 1), 0);
        if (ci[0] <= length && ci[1] <= length && ci[2] <= length) {
            if (ci[0] < ci[1]) {
                for (int y = std::max(0, ci[0]); y <= ci[1] && y <= length; ++y)
                    range[static_cast<std::size_t>(y)] = 1;
            } else {
                for (int y = std::max(0, ci[0]); y <= length; ++y)
                    range[static_cast<std::size_t>(y)] = 1;
                for (int y = 0; y <= ci[1] && y <= length; ++y)
                    range[static_cast<std::size_t>(y)] = 1;
            }
            if (ci[2] >= 0) {
                match[static_cast<std::size_t>(x)] = std::abs(
                    xdiff_at(xpos, breakpoint, length) -
                    xdiff_at(xpos, ci[2], length));
            }
            if (circular) {
                match[static_cast<std::size_t>(x)] = std::min(
                    match[static_cast<std::size_t>(x)],
                    xdiff_at(xpos, length, length) -
                        match[static_cast<std::size_t>(x)]);
            }
            if (breakpoint >= 0 && breakpoint <= length &&
                range[static_cast<std::size_t>(breakpoint)] != 0) {
                const int wrapped = xdiff_at(xpos, length, length) -
                    match[static_cast<std::size_t>(x)];
                match[static_cast<std::size_t>(x)] = -std::min(
                    match[static_cast<std::size_t>(x)], wrapped);
            }
        }
    }
    int smallest = length;
    int selected = 0;
    for (int x = 1; x <= cur_ci && x < static_cast<int>(match.size()); ++x) {
        if (smallest == 0) break;
        if (smallest > 0) {
            if (smallest > match[static_cast<std::size_t>(x)]) {
                smallest = match[static_cast<std::size_t>(x)];
                selected = x;
            }
        } else if (std::abs(smallest) > std::abs(match[static_cast<std::size_t>(x)]) &&
                   match[static_cast<std::size_t>(x)] <= 0) {
            smallest = match[static_cast<std::size_t>(x)];
            selected = x;
        }
    }
    if (smallest <= 0 && selected < static_cast<int>(intervals.size())) {
        result = intervals[static_cast<std::size_t>(selected)].values;
    } else if (selected < static_cast<int>(intervals.size())) {
        result = intervals[static_cast<std::size_t>(selected)].values;
        for (int& value : result) value = -value;
    }
    return result;
}

} // namespace

RdpBurtResult run_rdp_burt(
    const RdpScanState& scan_state,
    const std::array<int, 3>& representatives,
    int beginning,
    int ending,
    bool circular,
    int hmm_cycles,
    int repos_flag) {
    RdpBurtResult result;
    result.attempted = true;
    result.input_beginning = beginning;
    result.input_ending = ending;
    result.polished_beginning = beginning;
    result.polished_ending = ending;
    const int length = scan_state.sequence_length;
    if (length <= 0 || beginning < 1 || ending < 1 || beginning > length ||
        ending > length || hmm_cycles < 0 ||
        std::any_of(representatives.begin(), representatives.end(),
                    [&](const int sequence) {
                        return sequence < 0 || sequence > scan_state.next_no;
                    })) {
        return result;
    }

    auto input = make_hmm_input(scan_state, representatives, circular);
    result.information_rich_sites = input.informative;
    auto hmm = train_hmm(input, hmm_cycles);
    result.best_log_likelihood = hmm.best_like;
    result.trained = hmm.trained;
    if (!hmm.trained) return result;

    int slen = input.slen;
    if (circular && slen > 0) {
        const int selected = source_clng(slen / 2.0);
        for (int x = 0; x <= selected; ++x) {
            const int source = x + input.circular_offset + 1;
            if (source >= static_cast<int>(hmm.path.size())) break;
            hmm.path[static_cast<std::size_t>(x)] =
                hmm.path[static_cast<std::size_t>(source)];
            for (int state = 0; state < 3; ++state) {
                hmm.posterior[static_cast<std::size_t>(x * 3 + state)] =
                    hmm.posterior[static_cast<std::size_t>(source * 3 + state)];
            }
        }
        slen = selected;
    }
    if (slen <= 0) return result;
    if (slen >= 3) {
        hmm.path[static_cast<std::size_t>(slen - 1)] = circular
            ? hmm.path[2] : hmm.path[static_cast<std::size_t>(slen - 2)];
    }
    int cur_ci = -1;
    for (int x = 2; x <= slen - 2; ++x) {
        if (hmm.path[static_cast<std::size_t>(x)] ==
            hmm.path[static_cast<std::size_t>(x + 1)]) continue;
        ++cur_ci;
        RdpBurtInterval interval;
        if (x < slen - 2) {
            interval.values[2] = source_clng(
                input.xdiff[static_cast<std::size_t>(x)] +
                (input.xdiff[static_cast<std::size_t>(x + 1)] -
                 input.xdiff[static_cast<std::size_t>(x)]) / 2.0);
        } else {
            interval.values[2] = source_clng(
                input.xdiff[static_cast<std::size_t>(x)] +
                (length - input.xdiff[static_cast<std::size_t>(x)]) / 2.0);
        }
        for (int y = x; y >= -slen; --y) {
            int a = y;
            if (y <= 0) {
                if (circular) a = y + slen;
                else { interval.values[0] = 1; break; }
            }
            bool stop = false;
            for (int state = 0; state < 3; ++state) {
                const double probability = hmm.posterior[
                    static_cast<std::size_t>(a * 3 + state)];
                if (probability > 0.995 && interval.values[3] == 0) {
                    interval.values[3] = a > 1
                        ? input.xdiff[static_cast<std::size_t>(a - 1)] + 1 : 1;
                }
                if (probability > 0.999) {
                    interval.values[0] = a > 1
                        ? input.xdiff[static_cast<std::size_t>(a - 1)] + 1 : 1;
                    stop = true;
                    break;
                }
            }
            if (stop) break;
        }
        for (int y = x + 1; y <= slen * 2; ++y) {
            int a = y;
            if (y > slen) {
                if (circular) a = y - slen;
                else {
                    interval.values[1] = length;
                    interval.values[4] = length;
                    break;
                }
            }
            bool stop = false;
            for (int state = 0; state < 3; ++state) {
                const double probability = hmm.posterior[
                    static_cast<std::size_t>(a * 3 + state)];
                if (probability > 0.995 && interval.values[4] == 0) {
                    interval.values[4] = a < slen - 1
                        ? input.xdiff[static_cast<std::size_t>(a + 1)] - 1 : length;
                }
                if (probability > 0.999) {
                    interval.values[1] = a < slen - 1
                        ? input.xdiff[static_cast<std::size_t>(a + 1)] - 1 : length;
                    stop = true;
                    break;
                }
            }
            if (stop) break;
        }
        if (circular) {
            for (int field : {0, 1, 3, 4}) {
                if (interval.values[field] != 1 && interval.values[field] != length)
                    continue;
                if (field == 1 || field == 4) {
                    interval.values[field] = input.xdiff[1] - 1;
                    if (interval.values[field] < 1) interval.values[field] = 1;
                } else {
                    interval.values[field] = input.xdiff[static_cast<std::size_t>(slen)] + 1;
                    if (interval.values[field] > length) interval.values[field] = length;
                }
            }
        }
        if (interval.values[3] < interval.values[4]) {
            if (!(interval.values[2] > interval.values[3] &&
                  interval.values[2] < interval.values[4])) {
                interval.values[2] = interval.values[3] + source_clng(
                    (interval.values[4] - interval.values[3]) / 2.0);
            }
        } else if (!(interval.values[2] > interval.values[3] ||
                     interval.values[2] < interval.values[4])) {
            interval.values[2] = source_clng(
                (interval.values[3] + length - interval.values[4]) / 2.0);
            if (interval.values[2] > length) interval.values[2] -= length;
        }
        result.intervals.push_back(interval);
    }

    result.confidence.fill(0);
    if (cur_ci > -1) {
        const auto beginning_ci = match_bp_to_ci(
            result.intervals, cur_ci, beginning, length, input.xpos, circular);
        const auto ending_ci = match_bp_to_ci(
            result.intervals, cur_ci, ending, length, input.xpos, circular);
        result.confidence[0] = beginning_ci[0];
        result.confidence[1] = beginning_ci[1];
        result.confidence[2] = beginning_ci[2];
        result.confidence[6] = beginning_ci[3];
        result.confidence[7] = beginning_ci[4];
        result.confidence[3] = ending_ci[0];
        result.confidence[4] = ending_ci[1];
        result.confidence[5] = ending_ci[2];
        result.confidence[8] = ending_ci[3];
        result.confidence[9] = ending_ci[4];
    }

    const int original_beginning = beginning;
    const int original_ending = ending;
    const int span = beginning < ending ? ending - beginning : length + ending - beginning;
    auto near_enough = [&](const int input_coordinate, const int candidate) {
        const int distance = std::abs(input_coordinate - std::abs(candidate));
        return distance < span / 2.0 ||
            (circular && length - distance < span / 2.0);
    };
    if (std::abs(std::abs(result.confidence[2]) -
                 std::abs(result.confidence[5])) > 2) {
        if (result.confidence[2] >= 0 || repos_flag == 1) {
            if (circular) {
                beginning = std::abs(result.confidence[2]);
            } else if (near_enough(beginning, result.confidence[2])) {
                beginning = bounded_coordinate(result.confidence[2], length);
            }
        } else if (near_enough(beginning, result.confidence[2])) {
            beginning = bounded_coordinate(result.confidence[2], length);
            for (const int index : {0, 1, 2, 6, 7})
                result.confidence[static_cast<std::size_t>(index)] =
                    std::abs(result.confidence[static_cast<std::size_t>(index)]);
        }
        if (result.confidence[5] >= 0 || repos_flag == 1) {
            if (circular) {
                ending = std::abs(result.confidence[5]);
            } else if (near_enough(ending, result.confidence[5])) {
                ending = bounded_coordinate(result.confidence[5], length);
            }
        } else if (near_enough(ending, result.confidence[5])) {
            ending = bounded_coordinate(result.confidence[5], length);
            for (const int index : {3, 4, 5, 8, 9})
                result.confidence[static_cast<std::size_t>(index)] =
                    std::abs(result.confidence[static_cast<std::size_t>(index)]);
        }
    } else {
        result.single_transition_assignment = true;
        const int point = std::abs(result.confidence[2]);
        const int begin_distance = std::abs(beginning - point);
        const int end_distance = std::abs(ending - point);
        const bool begin_closer = begin_distance <= end_distance ||
            (circular && length - begin_distance <= end_distance);
        if (begin_closer) {
            if (near_enough(beginning, result.confidence[2]))
                beginning = bounded_coordinate(result.confidence[2], length);
            else result.confidence[0] = result.confidence[1] = result.confidence[2] = -1;
            result.confidence[3] = result.confidence[4] = result.confidence[5] = -1;
        } else {
            if (near_enough(ending, result.confidence[5]))
                ending = bounded_coordinate(result.confidence[5], length);
            else result.confidence[3] = result.confidence[4] = result.confidence[5] = -1;
            result.confidence[0] = result.confidence[1] = result.confidence[2] = -1;
        }
    }
    if (result.confidence[0] == 0 && result.confidence[1] == 0 && result.confidence[2] == 0) {
        result.confidence[0] = result.confidence[1] = result.confidence[2] = -1;
        beginning = original_beginning;
    }
    if (result.confidence[3] == 0 && result.confidence[4] == 0 && result.confidence[5] == 0) {
        result.confidence[3] = result.confidence[4] = result.confidence[5] = -1;
        ending = original_ending;
    }
    result.polished_beginning = beginning;
    result.polished_ending = ending;
    result.available = result.confidence[2] != -1 || result.confidence[5] != -1;
    return result;
}

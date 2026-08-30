#include "burt_state.hpp"

#include "MathFuncsDll.h"
#include "legacy_method_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <thread>
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

#if defined(NEXT_RDP_METHOD_MANUAL_THREADS) && !defined(NEXT_RDP_USE_REAL_OPENMP)

// One BenHMM restart is independent of every other restart except for two
// source-order choices: the MSVCRT random stream that seeds it, and the
// PathMax equality test carried between restarts.  Capture the random inputs
// in source order, evaluate both possible first-iteration/full-convergence
// outcomes in parallel, then replay those two choices serially.  This retains
// DoHMMCyclesSerial's exact winning model while making the supplied 21
// restarts useful work for the browser pthread pool.
struct HmmCycleSeed {
    double imbalance = 0.0;
    std::array<int, 3> favored_symbols{};
};

struct HmmCycleState {
    double maximum_likelihood = kNegativeLattice;
    std::array<double, 9> transition{};
    std::array<double, 9> emission{};
    std::array<double, 3> initial{};
    std::vector<int> path;
};

struct HmmCycleOutcome {
    std::vector<HmmCycleState> iterations;
    HmmCycleState exhausted_iterations;
    bool exhausted = false;
};

std::uint32_t next_msvcrt_random(std::uint32_t& state) {
    state = state * UINT32_C(214013) + UINT32_C(2531011);
    return (state >> 16U) & UINT32_C(0x7fff);
}

std::vector<HmmCycleSeed> make_hmm_cycle_seeds(const int cycles) {
    std::vector<HmmCycleSeed> seeds(static_cast<std::size_t>(cycles));
    std::uint32_t random_state = 3;
    constexpr double random_maximum = 32767.0;
    for (auto& seed : seeds) {
        const double imbalance_random =
            static_cast<double>(next_msvcrt_random(random_state)) /
            random_maximum;
        seed.imbalance = static_cast<int>(3.0 * imbalance_random + 1.0) /
            10.0;
        std::array<unsigned char, 3> used{};
        for (int state = 0; state < 3; ++state) {
            for (;;) {
                const double choice_random =
                    static_cast<double>(next_msvcrt_random(random_state)) /
                    random_maximum;
                const int choice = static_cast<int>(6.0 * choice_random) - 2;
                if (choice < 0 || choice >= 3 || used[choice] != 0) continue;
                used[choice] = 1;
                seed.favored_symbols[static_cast<std::size_t>(state)] = choice;
                break;
            }
        }
    }
    return seeds;
}

HmmCycleState capture_hmm_cycle_state(
    const double maximum_likelihood,
    const std::array<double, 9>& transition,
    const std::array<double, 9>& emission,
    const std::array<double, 3>& initial,
    const std::vector<int>& path) {
    return {maximum_likelihood, transition, emission, initial, path};
}

HmmCycleOutcome run_hmm_cycle(
    const HmmInput& input, const HmmCycleSeed& seed) {
    constexpr int states = 3;
    constexpr int symbols = 3;
    const int lattice_stride = input.slen + 1;
    std::array<double, states * states> transition{};
    std::array<double, symbols * states> emission{};
    std::array<double, states> initial{};
    std::array<double, states * states> options{};
    std::array<double, states * states> transition_count{};
    std::array<double, symbols * states> state_count{};
    std::vector<double> lattice_backtrack(
        static_cast<std::size_t>(lattice_stride * states), 0.0);
    std::vector<double> lattice_values(
        static_cast<std::size_t>(lattice_stride * states), 0.0);
    std::vector<int> path(static_cast<std::size_t>(input.slen + 4), 0);

    const float initial_change_probability =
        static_cast<float>(5.0 / static_cast<double>(input.sequence_length));
    for (int from = 0; from < states; ++from) {
        for (int to = 0; to < states; ++to) {
            transition[from + to * states] = from == to
                // The supplied body evaluates `1 - iVal` at float precision
                // before casting it to double for log().
                ? std::log(static_cast<double>(
                    1.0F - initial_change_probability))
                : std::log(initial_change_probability / (states - 1.0));
        }
    }
    for (int state = 0; state < states; ++state) {
        for (int symbol = 0; symbol < symbols; ++symbol) {
            const double probability =
                symbol == seed.favored_symbols[static_cast<std::size_t>(state)]
                ? 1.0 / symbols + seed.imbalance * 2.0
                : 1.0 / symbols - seed.imbalance;
            emission[symbol + state * symbols] = std::log(probability);
        }
        initial[state] = std::log(1.0 / states);
    }

    HmmCycleOutcome outcome;
    outcome.iterations.reserve(16);
    double preceding_likelihood = std::numeric_limits<double>::quiet_NaN();
    for (int iteration = 1; iteration <= 100; ++iteration) {
        for (int state = 0; state < states; ++state) {
            lattice_values[state * lattice_stride] =
                emission[input.recode[0] + state * symbols] + initial[state];
        }
        (void)MathFuncs::MyMathFuncs::ViterbiCP(
            input.slen, symbols, states, options.data(),
            const_cast<unsigned char*>(input.recode.data()),
            lattice_values.data(), transition.data(), emission.data(),
            lattice_backtrack.data());
        const double maximum_likelihood =
            MathFuncs::MyMathFuncs::GetLaticePathP(
                input.slen, states, lattice_values.data(),
                lattice_backtrack.data(), path.data());
        outcome.iterations.push_back(capture_hmm_cycle_state(
            maximum_likelihood, transition, emission, initial, path));
        if (preceding_likelihood == maximum_likelihood) {
            return outcome;
        }
        preceding_likelihood = maximum_likelihood;

        transition_count.fill(0.0);
        state_count.fill(0.0);
        (void)MathFuncs::MyMathFuncs::UpdateCountsP(
            input.slen, symbols, states, path.data(),
            const_cast<unsigned char*>(input.recode.data()),
            transition_count.data(), state_count.data());
        constexpr double fudge = 0.01;
        for (int from = 0; from < states; ++from) {
            double total = 0.0;
            for (int to = 0; to < states; ++to) {
                total += transition_count[from + to * states];
            }
            total += fudge * states;
            for (int to = 0; to < states; ++to) {
                const int index = from + to * states;
                transition_count[index] += fudge;
                transition[index] = std::log(transition_count[index] / total);
            }
        }
        for (int state = 0; state < states; ++state) {
            double total = 0.0;
            for (int symbol = 0; symbol < symbols; ++symbol) {
                total += state_count[symbol + state * symbols];
            }
            total += fudge * symbols;
            for (int symbol = 0; symbol < symbols; ++symbol) {
                emission[symbol + state * symbols] = std::log(
                    (state_count[symbol + state * symbols] + fudge) / total);
            }
        }
    }
    // DoHMMCyclesSerial still copies the post-update model when the fixed
    // 100-iteration bound is exhausted, even though MaxL and LaticePath came
    // from the preceding Viterbi pass.
    outcome.exhausted_iterations = capture_hmm_cycle_state(
        preceding_likelihood, transition, emission, initial, path);
    outcome.exhausted = true;
    return outcome;
}

float train_hmm_restarts_parallel(
    const HmmInput& input, const int hmm_cycles,
    double* transition, double* emission, double* initial, int* path) {
    const int cycle_count = hmm_cycles + 1;
    const auto seeds = make_hmm_cycle_seeds(cycle_count);
    std::vector<HmmCycleOutcome> outcomes(static_cast<std::size_t>(cycle_count));
    const int workers = std::min(
        std::max(1, rdp_method_worker_threads()), cycle_count);
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(workers));
    for (int worker = 0; worker < workers; ++worker) {
        const int first = cycle_count * worker / workers;
        const int last = cycle_count * (worker + 1) / workers;
        threads.emplace_back([&, first, last] {
            for (int cycle = first; cycle < last; ++cycle) {
                outcomes[static_cast<std::size_t>(cycle)] = run_hmm_cycle(
                    input, seeds[static_cast<std::size_t>(cycle)]);
            }
        });
    }
    for (auto& thread : threads) thread.join();

    double best_likelihood = -1000000.0;
    double path_maximum = 0.0;
    for (const auto& outcome : outcomes) {
        const HmmCycleState* selected = nullptr;
        for (const auto& iteration : outcome.iterations) {
            if (path_maximum == iteration.maximum_likelihood) {
                selected = &iteration;
                break;
            }
            path_maximum = iteration.maximum_likelihood;
        }
        if (selected == nullptr) {
            if (!outcome.exhausted) {
                selected = &outcome.iterations.back();
            } else {
                selected = &outcome.exhausted_iterations;
            }
        }
        if (selected->maximum_likelihood <= best_likelihood) continue;
        best_likelihood = selected->maximum_likelihood;
        std::copy(selected->transition.begin(), selected->transition.end(), transition);
        std::copy(selected->emission.begin(), selected->emission.end(), emission);
        std::copy(selected->initial.begin(), selected->initial.end(), initial);
        std::copy(selected->path.begin(), selected->path.end(), path);
    }
    return static_cast<float>(best_likelihood);
}

#endif

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
    float best = kNegativeLattice;
#if defined(NEXT_RDP_METHOD_MANUAL_THREADS) && !defined(NEXT_RDP_USE_REAL_OPENMP)
    if (rdp_method_worker_threads() > 1) {
        best = train_hmm_restarts_parallel(
            input, hmm_cycles, transition.data(), emission.data(),
            initial.data(), path.data());
    } else
#endif
    {
        best = MathFuncs::MyMathFuncs::DoHMMCyclesSerial(
            3, input.slen, hmm_cycles, input.sequence_length, states, symbols,
            const_cast<unsigned char*>(input.recode.data()), transition.data(),
            emission.data(), initial.data(), path.data());
    }
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

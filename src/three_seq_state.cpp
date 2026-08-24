#include "three_seq_state.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

double approx_norm_pdf(const double x) {
    constexpr double pi = 3.14159265359;
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * pi);
}

double approx_norm_cdf(const double x) {
    constexpr double b1 = 0.31938153;
    constexpr double b2 = -0.356563782;
    constexpr double b3 = 1.781477937;
    constexpr double b4 = -1.821255978;
    constexpr double b5 = 1.330274429;
    constexpr double p = 0.2316419;
    constexpr double c = 0.39894228;
    if (x > 6.0) return 1.0;
    if (x < -6.0) return 0.0;
    if (x >= 0.0) {
        const double t = 1.0 / (1.0 + p * x);
        return 1.0 - c * std::exp(-x * x / 2.0) * t *
            (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
    }
    const double t = 1.0 / (1.0 - p * x);
    return c * std::exp(-x * x / 2.0) * t *
        (t * (t * (t * (t * b5 + b4) + b3) + b2) + b1);
}

double approx_nu(const double x) {
    if (x <= 0.0) return 0.0;
    return ((approx_norm_cdf(x / 2.0) - 0.5) * 2.0) / x /
        (approx_norm_pdf(x / 2.0) + x * approx_norm_cdf(x / 2.0) / 2.0);
}

double siegmund_discrete(const int m, const int n, const int k) {
    const double total = static_cast<double>(m + n);
    if (total == 0.0) return -1.0;
    const double b = k - 0.5;
    const double xi = n - m;
    double h1 = -2.0 * b * (b - xi) / total;
    if (h1 > 700.0) h1 = 700.0;
    double p1 = std::exp(h1);
    if (p1 > 1.0e200) p1 = 1.0e200;
    const double p2 = p1 *
        (2.0 * (2.0 * b - xi) * (b - xi) / total + 1.0);
    const double an = 2.0 * (2.0 * b - xi) / total;
    double p3 = approx_nu(an);
    p3 = p3 * p3 * p2;
    if (p3 < -1.0) p3 = -1.0;
    return 1.0 - std::exp(-p3);
}

double get_ts_pvalue(const int m, const int n, const int k,
                     const std::vector<float>& table,
                     const int table_bound) {
    if (table_bound < 1) return siegmund_discrete(m, n, k);
    if (m >= table_bound - 1 || n >= table_bound - 1 ||
        k >= table_bound - 1) {
        const double direct = siegmund_discrete(m, n, k);
        if (direct > 0.0 && direct < 1.0) return direct;
    }
    int sm = m;
    int sn = n;
    int sk = k;
    double scale = 1.0;
    if (sm >= table_bound - 1 || sn >= table_bound - 1 ||
        sk >= table_bound - 1) {
        const int maximum = std::max({sm, sn, sk});
        scale = static_cast<double>(maximum) / (table_bound - 2);
        sm = static_cast<int>(sm / scale);
        sn = static_cast<int>(sn / scale);
        sk = static_cast<int>(sk / scale);
        if (sm > 0) {
            scale = static_cast<double>(m) / sm;
        } else if (sn > 0) {
            scale = static_cast<double>(n) / sn;
        }
    }
    const int width = table_bound + 1;
    const int lookup_k = std::max(sk, 0);
    const auto index = static_cast<std::size_t>(sm) +
        static_cast<std::size_t>(width) * sn +
        static_cast<std::size_t>(width) * width * lookup_k;
    if (index >= table.size()) {
        throw std::runtime_error("3seqTable dimensions differ");
    }
    const double original = table[index];
    double result = original;
    if (scale > 1.0) result = std::pow(result, scale);
    if (result == 0.0 && original > 0.0) result = 1.0e-300;
    return result;
}

double correct_probability(const double raw, const int mc_flag,
                           const int correction, double& product) {
    if (mc_flag != 0) {
        product = raw;
        return raw;
    }
    product = raw * correction;
    if (raw >= 1.0 || product >= 1.0) return 1.0;
    if (raw > 1.0e-15) return 1.0 - std::pow(1.0 - raw, correction);
    return product;
}

bool passes(const double probability, const double product,
            const double threshold) {
    return ((probability < 1.0 && probability <= threshold) ||
            (probability == 1.0 && product < threshold)) && product > 0.0;
}

void check_wrap(const int sequence_length, const int informative_last,
                const int sign, const bool circular, int& excursion,
                int& beginning, int& ending,
                const std::vector<int>& difference_position,
                const std::vector<int>& position_difference,
                const std::vector<int>& walk) {
    if (informative_last < 0) return;
    if (beginning == 0) beginning = difference_position[0];
    int maximum = walk[position_difference[beginning]] * sign;
    const int tail = walk[informative_last] * sign;
    const int stop = beginning < ending
        ? position_difference[beginning] : position_difference[ending];
    for (int x = 0; x <= stop; ++x) {
        const int height = tail + walk[x] * sign;
        if (height > maximum) {
            maximum = height;
            beginning = difference_position[x];
        }
        if (maximum - height > excursion) {
            excursion = maximum - height;
            ending = difference_position[x];
        }
    }
    if (position_difference[beginning] < informative_last) {
        beginning = difference_position[position_difference[beginning] + 1];
    } else {
        beginning = difference_position[0];
    }
    if (!circular && beginning > ending) {
        const int new_beginning = position_difference[ending] < informative_last
            ? difference_position[position_difference[ending] + 1] : 1;
        ending = position_difference[beginning] > 0
            ? difference_position[position_difference[beginning] - 1]
            : sequence_length;
        beginning = new_beginning;
    }
}

}  // namespace

RdpThreeSeqResult evaluate_rdp_three_seq(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const bool circular, const int mc_flag, const int mc_correction,
    const double lowest_probability, const std::vector<float>& table,
    const int table_bound) {
    const int length = scan_state.sequence_length;
    const int stride = length + 1;
    for (const int sequence : sequences) {
        if (sequence < 0 || sequence > scan_state.next_no) {
            throw std::runtime_error("three-sequence index differs");
        }
    }
    RdpThreeSeqResult result;
    std::vector<int> position_difference(length + 1, 0);
    std::vector<int> difference_position(length + 1, 0);
    std::vector<int> walk(length + 1, 0);
    int current = 0;
    int maximum_descent = 0;
    int maximum_ascent = 0;
    int informative = 0;
    int maximum = 0;
    int minimum = 0;
    int first_count = 0;
    int second_count = 0;
    int beginning = 0;
    int ending = 0;
    int beginning2 = 0;
    int ending2 = 0;
    const auto& data = scan_state.sequence_data;
    for (int position = 0; position <= length; ++position) {
        position_difference[position] = informative;
        const short first = data[position + sequences[0] * stride];
        const short second = data[position + sequences[1] * stride];
        if (second == first || first == 46 || second == 46) continue;
        const short third = data[position + sequences[2] * stride];
        if (third == 46) continue;
        if (third == first) {
            ++current;
            walk[informative] = current;
            difference_position[informative] = position;
            if (current > maximum) {
                maximum = current;
                beginning = position;
            }
            if (current - minimum > maximum_ascent) {
                maximum_ascent = current - minimum;
                ending2 = position;
            }
            ++first_count;
            ++informative;
            position_difference[position] = informative;
        } else if (third == second) {
            --current;
            walk[informative] = current;
            difference_position[informative] = position;
            ++second_count;
            if (maximum - current > maximum_descent) {
                maximum_descent = maximum - current;
                ending = position;
            }
            if (current < minimum) {
                minimum = current;
                beginning2 = position;
            }
            ++informative;
            position_difference[position] = informative;
        }
    }
    result.informative_last = informative - 1;
    auto& first = result.sides[0];
    auto& second = result.sides[1];
    first = {beginning, ending, first_count, second_count,
             maximum_descent,
             get_ts_pvalue(first_count, second_count, maximum_descent,
                           table, table_bound), false};
    second = {beginning2, ending2, second_count, first_count,
              maximum_ascent,
              get_ts_pvalue(second_count, first_count, maximum_ascent,
                            table, table_bound), false};
    if (result.informative_last < 3) return result;
    check_wrap(length, result.informative_last, 1, circular,
               first.excursion, first.beginning, first.ending,
               difference_position, position_difference, walk);
    check_wrap(length, result.informative_last, -1, circular,
               second.excursion, second.beginning, second.ending,
               difference_position, position_difference, walk);
    if (second.probability < first.probability) std::swap(first, second);
    if ((first.second_count > 0 && first.excursion == 1) ||
        first.second_count - first.first_count == first.excursion) {
        return result;
    }
    double product = 1.0;
    const double corrected = correct_probability(
        first.probability, mc_flag, mc_correction, product);
    first.significant = passes(corrected, product, lowest_probability);
    double second_product = 1.0;
    const double second_corrected = correct_probability(
        second.probability, mc_flag, mc_correction, second_product);
    second.significant = first.significant &&
        passes(second_corrected, second_product, lowest_probability);
    return result;
}


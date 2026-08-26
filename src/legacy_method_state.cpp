#include "legacy_method_state.hpp"

#include "MathFuncsDll.h"
#include "three_seq_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace {

std::uint64_t fnv1a64(const unsigned char* data, const std::size_t bytes) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (std::size_t index = 0; index < bytes; ++index) {
        hash ^= data[index];
        hash *= 1099511628211ULL;
    }
    return hash;
}

double store_probability(const std::vector<double>& values,
                         const int upper_bound, const int program,
                         const int sequence) {
    const auto index = static_cast<std::size_t>(program) +
        static_cast<std::size_t>(sequence) * (upper_bound + 1);
    if (index >= values.size()) {
        throw std::runtime_error("StoreLPV lookup exceeds its bounds");
    }
    return values[index];
}

std::array<int, 3> choose_active(
    const std::array<int, 3>& sequences, const int program,
    const std::vector<double>& store_lpv, const int store_lpv_ub,
    const RdpLegacyEventAllocator& allocator) {
    const int first = sequences[0];
    const int second = sequences[1];
    const int third = sequences[2];
    if (allocator.count(first) < allocator.count(second) &&
        allocator.count(first) < allocator.count(third)) {
        return {first, second, third};
    }
    if (allocator.count(second) < allocator.count(first) &&
        allocator.count(second) < allocator.count(third)) {
        return {second, first, third};
    }
    if (allocator.count(third) < allocator.count(first) &&
        allocator.count(third) < allocator.count(second)) {
        return {third, first, second};
    }
    const double first_p = store_probability(
        store_lpv, store_lpv_ub, program, first);
    const double second_p = store_probability(
        store_lpv, store_lpv_ub, program, second);
    const double third_p = store_probability(
        store_lpv, store_lpv_ub, program, third);
    if (first_p >= second_p && first_p >= third_p) {
        return {first, second, third};
    }
    if (second_p >= first_p && second_p >= third_p) {
        return {second, first, third};
    }
    return {third, first, second};
}

// VB6 CLng uses bankers' rounding.
int vb_clng(const double value) {
    const double floor_value = std::floor(value);
    const double fraction = value - floor_value;
    if (fraction < 0.5) return static_cast<int>(floor_value);
    if (fraction > 0.5) return static_cast<int>(floor_value + 1.0);
    const int lower = static_cast<int>(floor_value);
    return (lower & 1) == 0 ? lower : lower + 1;
}

// Module8.QCTopolChange, used by the GENECONV and MaxChi paths before their
// breakpoints are centred.  Scores is Module30.SubSeq, whose second dimension
// uses the alignment-length stride even though only informative sites are read.
std::uint8_t legacy_qc_topology_change(
    const int beginning, const int ending,
    const std::vector<int>& position_difference,
    const std::vector<char>& scores, const int score_stride,
    const int informative) {
    const int left = position_difference[beginning];
    const int right = position_difference[ending];
    std::array<int, 3> outside{};
    std::array<int, 3> inside{};
    const auto add = [&](std::array<int, 3>& totals, const int position) {
        for (int comparison = 0; comparison < 3; ++comparison) {
            totals[comparison] += static_cast<unsigned char>(
                scores[position + comparison * score_stride]);
        }
    };
    if (left < right) {
        for (int position = 0; position < left; ++position)
            add(outside, position);
        for (int position = left; position <= right; ++position)
            add(inside, position);
        for (int position = right + 1; position <= informative; ++position)
            add(outside, position);
    } else {
        for (int position = 0; position <= right; ++position)
            add(inside, position);
        for (int position = right + 1; position < left; ++position)
            add(outside, position);
        for (int position = left; position <= informative; ++position)
            add(inside, position);
    }
    for (int comparison = 0; comparison < 3; ++comparison) {
        const int other_a = (comparison + 1) % 3;
        const int other_b = (comparison + 2) % 3;
        if (inside[comparison] > inside[other_a] &&
            inside[comparison] > inside[other_b] &&
            outside[comparison] > outside[other_a] &&
            outside[comparison] > outside[other_b]) {
            return 1;
        }
    }
    return 0;
}

// Module3.CentreBP(0, 0, ...), restricted to the initial scan where
// SEventNumber=0.  The omitted missing-data branches are guarded by
// SEventNumber > 0 in the VB source; the circular/linear clamps remain exact.
void legacy_centre_geneconv_breakpoints(
    const int length, const int circular, const int informative,
    std::vector<int>& position_difference,
    const std::vector<int>& difference_position,
    int& beginning, int& ending) {
    if (position_difference[beginning] - 1 > 0) {
        beginning -= vb_clng(
            ((beginning -
              difference_position[position_difference[beginning] - 1]) /
             2.0) - 0.1);
    } else {
        beginning -= vb_clng(
            ((beginning + length - difference_position[informative]) / 2.0) -
            0.1);
    }
    if (beginning == 0) {
        beginning = 1;
    } else if (beginning < 1) {
        beginning = circular == 0 ? 1 : length + beginning;
    }

    position_difference[length] = informative;
    if (position_difference[ending] + 1 <= informative) {
        ending += vb_clng(
            ((difference_position[position_difference[ending] + 1] - ending) /
             2.0) - 0.1);
    } else {
        ending += vb_clng(
            ((difference_position[1] + (length - ending)) / 2.0) - 0.1);
    }
    if (ending > length) {
        ending = circular == 0 ? length : ending - length;
    }
}

int critical_difference(const int window_size, const double lowest_probability) {
    int half_window = vb_clng(static_cast<double>(window_size) / 2.0);
    double low_probability = lowest_probability / 6.0;
    if (low_probability < 0.0001) low_probability = 0.0001;
    const double high_probability = low_probability * 1.000000001;
    low_probability *= 0.999999999;

    double chi = 5.0;
    double probability = 10.0;
    while (probability > low_probability) {
        chi *= 2.0;
        probability = MathFuncs::MyMathFuncs::ChiPVal2P(chi);
        if (probability == 1.0e-20) {
            low_probability = probability;
            break;
        }
    }
    double lower_chi = 0.0;
    double upper_chi = chi * 2.0;
    double last_probability = 0.0;
    while (true) {
        probability = MathFuncs::MyMathFuncs::ChiPVal2P(chi);
        if (probability * 1.000000001 > last_probability &&
            probability * 0.9999999999 < last_probability) break;
        if (probability < low_probability) {
            const double prior = chi;
            chi -= (chi - lower_chi) / 2.0;
            if (chi == prior) break;
            upper_chi = prior;
        } else if (probability > high_probability) {
            const double prior = chi;
            chi += (upper_chi - chi) / 2.0;
            if (chi == prior) break;
            lower_chi = prior;
        } else {
            break;
        }
        last_probability = probability;
    }

    int lower_difference = 0;
    int upper_difference = half_window;
    int difference = vb_clng(static_cast<double>(half_window) / 2.0);
    while (true) {
        if (difference < 2) {
            difference = 2;
            break;
        }
        const double a = 0.0;
        const double b = half_window;
        const double c = difference;
        const double d = half_window - difference;
        double value = std::pow(a * d - b * c, 2.0) * half_window * 2.0;
        value /= a + b;
        value /= c + d;
        value /= a + c;
        value /= b + d;
        if (value > chi) {
            const int prior = difference;
            difference -= vb_clng(
                static_cast<double>(difference - lower_difference) / 2.0);
            upper_difference = prior;
        } else if (value < chi) {
            const int prior = difference;
            difference += vb_clng(
                static_cast<double>(upper_difference - difference) / 2.0);
            lower_difference = prior;
        } else {
            break;
        }
        if (upper_difference == difference || difference == lower_difference)
            break;
    }
    return difference - 1;
}

// These are the non-P DNA.dll routines used by CXoverA during FinalTrim.
// Keep their indexing and tie behavior literal: the P variants are very
// similar, but this path is part of the compatibility target and must not be
// silently substituted with a newer helper.
int legacy_win_score_calc4(
    const int criticaldiff, const int hwindow_width,
    const int len_xover_seq, const int len_seq, const int seq1,
    const int seq2, const int seq3, std::vector<unsigned char>& scores,
    const std::vector<int>& xdiff_pos, const std::vector<short>& seq_num,
    std::vector<int>& win_scores) {
    (void)seq3;
    int goon = 0;
    for (int x = 1; x <= len_xover_seq; ++x) {
        const int first = xdiff_pos[x] + len_seq * seq1;
        const int second = xdiff_pos[x] + len_seq * seq2;
        scores[x] = seq_num[first] == seq_num[second] ? 1 : 0;
        win_scores[x] = 0;
    }
    win_scores[0] = 0;
    for (int x = len_xover_seq - hwindow_width + 1;
         x <= len_xover_seq; ++x) {
        win_scores[0] += scores[x];
    }
    for (int x = 1; x <= hwindow_width; ++x) {
        const int lo = (len_xover_seq - hwindow_width) + x;
        win_scores[x] = win_scores[x - 1] - scores[lo] + scores[x];
    }
    for (int x = hwindow_width + 1; x <= len_xover_seq; ++x) {
        const int lo = x - hwindow_width;
        win_scores[x] = win_scores[x - 1] - scores[lo] + scores[x];
    }
    for (int x = len_xover_seq + 1;
         x < len_xover_seq + hwindow_width; ++x) {
        const int lo = x - hwindow_width;
        const int ro = x - len_xover_seq;
        win_scores[x] = win_scores[x - 1] - scores[lo] + scores[ro];
    }
    for (int x = 1; x < len_xover_seq + hwindow_width; ++x) {
        const int delta = win_scores[x] - win_scores[x + hwindow_width - 1];
        if (delta > criticaldiff || delta < -criticaldiff) {
            goon = 1;
            break;
        }
    }
    return goon;
}

double legacy_calc_chi_vals3(
    const int /*criticaldiff*/, const int hwindow_width,
    const int len_xover_seq, const int /*len_seq*/,
    const std::vector<int>& win_scores, std::vector<double>& chi_values) {
    double maximum = 0.0;
    for (int x = 0; x < len_xover_seq; ++x) {
        const double a = win_scores[x];
        const double c = win_scores[x + hwindow_width];
        const double b = hwindow_width - a;
        const double d = hwindow_width - c;
        if (a + c > 0.0 && b + d > 0.0) {
            const double e = a * d - b * c;
            const double chi = e * e * 2.0 /
                (hwindow_width * (a + c) * (b + d));
            chi_values[x] = chi;
            if (chi > maximum) maximum = chi;
        } else {
            chi_values[x] = 0.0;
        }
    }
    return maximum;
}

void legacy_find_mchi3(
    const int /*len_seq*/, const int len_xover_seq, int& maximum_position,
    short& comparison, double& maximum,
    const std::vector<double>& chi_values) {
    int best_position = -1;
    short best_comparison = -1;
    double best_value = 0.0;
    for (int x = 0; x < len_xover_seq; ++x) {
        const double value = chi_values[x];
        if (value > best_value) {
            best_value = value;
            best_position = x;
            best_comparison = 0;
        }
    }
    maximum = best_value;
    maximum_position = best_position;
    comparison = best_comparison;
}

struct ChiLookupTable {
    static constexpr int max_window = 150;
    std::array<int, max_window + 1> map{};
    std::vector<float> values;

    ChiLookupTable() {
        map.fill(-1);
        std::size_t size = 0;
        for (int window = 5; window <= max_window; ++window) {
            size += static_cast<std::size_t>(window + 1) * (window + 1);
        }
        values.reserve(size);
        for (int window = 5; window <= max_window; ++window) {
            map[window] = static_cast<int>(values.size());
            for (int left = 0; left <= window; ++left) {
                for (int right = 0; right <= window; ++right) {
                    const double b = window - left;
                    const double d = window - right;
                    double chi = 0.0;
                    if (left + right > 0 && b + d > 0) {
                        chi = std::pow(left * d - b * right, 2.0) *
                              window * 2.0;
                        chi /= left + b;
                        chi /= right + d;
                        chi /= left + right;
                        chi /= b + d;
                    }
                    values.push_back(static_cast<float>(chi));
                }
            }
        }
    }
};

const ChiLookupTable& chi_lookup_table() {
    static const ChiLookupTable table;
    return table;
}

// Module4.FillFSSRDP constructs this lookup with nine nested 0..4 loops.
// GCXoverD uses the resulting FSSGC table through FindSubSeqGCAP7 whenever
// UseCompress=1, which is the active initial command-line scan path.
std::vector<unsigned char>& geneconv_fss_table() {
    constexpr int upper_bound = 125;
    constexpr int width = upper_bound + 1;
    static std::vector<unsigned char> table(
        static_cast<std::size_t>(4) * width * width * width, 0);
    return table;
}

void initialize_geneconv_fss_table(std::vector<unsigned char>& table) {
    constexpr int width = 126;
    const auto decode = [](const int code) {
        return std::array<int, 3>{code / 25, (code / 5) % 5, code % 5};
    };
    // FillFSSRDP's VB layout is (site, first packed symbol, second packed
    // symbol, third packed symbol), with the first dimension contiguous.
    // DNA5 receives the same flat layout through the DLL ABI.
    for (int first_code = 0; first_code < 125; ++first_code) {
        for (int second_code = 0; second_code < 125; ++second_code) {
            for (int third_code = 0; third_code < 125; ++third_code) {
                const auto first = decode(first_code);
                const auto second = decode(second_code);
                const auto third = decode(third_code);
                const std::size_t offset =
                    static_cast<std::size_t>(first_code) * 4 +
                    static_cast<std::size_t>(second_code) * 4 * width +
                    static_cast<std::size_t>(third_code) * 4 * width * width;
                int total = 0;
                for (int site = 0; site < 3; ++site) {
                    int action = 0;
                    if (first[site] != 0 && second[site] != 0 &&
                        third[site] != 0 &&
                        (first[site] != third[site] ||
                         first[site] != second[site])) {
                        if (first[site] == second[site]) action = 1;
                        else if (first[site] == third[site]) action = 2;
                        else if (third[site] == second[site]) action = 3;
                        else action = 7;
                    }
                    table[offset + site] =
                        static_cast<unsigned char>(action);
                    total += action;
                }
                table[offset + 3] = static_cast<unsigned char>(total);
            }
        }
    }
}

void ensure_geneconv_fss_table() {
    static std::once_flag initialized;
    std::call_once(initialized, [] {
        initialize_geneconv_fss_table(geneconv_fss_table());
    });
}

std::vector<unsigned char>& maxchi_fss_table() {
    constexpr int width = 126;
    static std::vector<unsigned char> table(
        static_cast<std::size_t>(3) * width * width * width, 0);
    static std::once_flag initialized;
    std::call_once(initialized, [] {
        const auto decode = [](const int code) {
            return std::array<int, 3>{code / 25, (code / 5) % 5,
                                      code % 5};
        };
        for (int first_code = 0; first_code < 125; ++first_code) {
            for (int second_code = 0; second_code < 125; ++second_code) {
                for (int third_code = 0; third_code < 125; ++third_code) {
                    const auto first = decode(first_code);
                    const auto second = decode(second_code);
                    const auto third = decode(third_code);
                    for (int site = 0; site < 3; ++site) {
                        const std::size_t offset =
                            static_cast<std::size_t>(site) +
                            3 * static_cast<std::size_t>(first_code) +
                            3 * static_cast<std::size_t>(second_code) * width +
                            3 * static_cast<std::size_t>(third_code) * width * width;
                        table[offset] =
                            first[site] != 0 && second[site] != 0 &&
                            third[site] != 0 &&
                            (first[site] != third[site] ||
                             first[site] != second[site]) ? 1 : 0;
                    }
                }
            }
        }
    });
    return table;
}

std::vector<unsigned char>& chimaera_fss_table() {
    constexpr int width = 126;
    static std::vector<unsigned char> table(
        static_cast<std::size_t>(4) * width * width * width, 0);
    static std::once_flag initialized;
    std::call_once(initialized, [] {
        const auto decode = [](const int code) {
            return std::array<int, 3>{code / 25, (code / 5) % 5,
                                      code % 5};
        };
        for (int first_code = 0; first_code < 125; ++first_code) {
            for (int second_code = 0; second_code < 125; ++second_code) {
                for (int third_code = 0; third_code < 125; ++third_code) {
                    const auto first = decode(first_code);
                    const auto second = decode(second_code);
                    const auto third = decode(third_code);
                    int total = 0;
                    for (int site = 0; site < 3; ++site) {
                        unsigned char action = 0;
                        if (first[site] != 0 && second[site] != 0 &&
                            third[site] != 0 &&
                            (first[site] != third[site] ||
                             first[site] != second[site])) {
                            if (first[site] == second[site]) action = 1;
                            else if (first[site] == third[site]) action = 2;
                            else if (third[site] == second[site]) action = 3;
                        }
                        const std::size_t offset =
                            static_cast<std::size_t>(site) +
                            4 * static_cast<std::size_t>(first_code) +
                            4 * static_cast<std::size_t>(second_code) * width +
                            4 * static_cast<std::size_t>(third_code) * width * width;
                        table[offset] = action;
                        total += action;
                    }
                    const std::size_t total_offset =
                        3 + 4 * static_cast<std::size_t>(first_code) +
                        4 * static_cast<std::size_t>(second_code) * width +
                        4 * static_cast<std::size_t>(third_code) * width * width;
                    table[total_offset] = static_cast<unsigned char>(total);
                }
            }
        }
    });
    return table;
}

void populate_geneconv_fss_triplet(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    std::vector<unsigned char>& table) {
    ensure_geneconv_fss_table();
    constexpr int width = 126;
    const int stride = scan_state.compressed_sequence_ub + 1;
    const auto decode = [](const int code) {
        return std::array<int, 3>{
            code / 25, (code / 5) % 5, code % 5};
    };
    for (int position = 1;
         position <= scan_state.compressed_sequence_ub; ++position) {
        const int first_code = scan_state.compressed_sequence[
            position + sequences[0] * stride];
        const int second_code = scan_state.compressed_sequence[
            position + sequences[1] * stride];
        const int third_code = scan_state.compressed_sequence[
            position + sequences[2] * stride];
        // FillFSSRDP leaves the unused code-125 plane zero-filled.
        if (first_code >= 125 || second_code >= 125 || third_code >= 125)
            continue;
        const auto first = decode(first_code);
        const auto second = decode(second_code);
        const auto third = decode(third_code);
        const std::size_t offset =
            static_cast<std::size_t>(first_code) * 4 +
            static_cast<std::size_t>(second_code) * 4 * width +
            static_cast<std::size_t>(third_code) * 4 * width * width;
        int total = 0;
        for (int site = 0; site < 3; ++site) {
            int action = 0;
            if (first[site] != 0 && second[site] != 0 && third[site] != 0 &&
                (first[site] != third[site] ||
                 first[site] != second[site])) {
                if (first[site] == second[site]) action = 1;
                else if (first[site] == third[site]) action = 2;
                else if (third[site] == second[site]) action = 3;
                else action = 7;
            }
            table[offset + site] = static_cast<unsigned char>(action);
            total += action;
        }
        table[offset + 3] = static_cast<unsigned char>(total);
    }
}

RdpMethodScreenResult screen_rdp_geneconv_candidates_impl(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int target) {
    RdpMethodScreenResult result;
    const int triplet_count = scan_state.analysis_list_last + 1;
    if (triplet_count <= 0) return result;
    if (scan_state.analysis_list.size() <
        static_cast<std::size_t>(triplet_count) * 3) {
        throw std::runtime_error("GENECONV analysis list is truncated");
    }
    if (store_lpv_ub < 1 || store_lpv.size() <
        static_cast<std::size_t>(store_lpv_ub + 1) *
            static_cast<std::size_t>(scan_state.next_no + 1)) {
        throw std::runtime_error("GENECONV StoreLPV is truncated");
    }
    ensure_geneconv_fss_table();
    std::vector<unsigned char> redo(static_cast<std::size_t>(triplet_count), 0);
    const int correction = std::max(1, correction_tests);
    const double sub_threshold = mc_flag != 0
        ? lowest_probability
        : lowest_probability / static_cast<double>(correction);
    constexpr int gc_dimension = 2999;
    (void)MathFuncs::MyMathFuncs::AlistGC2(
        store_lpv_ub, const_cast<double*>(store_lpv.data()), 0, 1,
        gc_dimension, const_cast<short*>(scan_state.analysis_list.data()),
        scan_state.analysis_list_last, 0, scan_state.analysis_list_last,
        scan_state.next_no, sub_threshold, redo.data(), circular, correction,
        mc_flag, lowest_probability,
        target != 0 ? target : static_cast<int>(std::nearbyint(
            static_cast<double>(scan_state.sequence_length) / 10.0)),
        scan_state.sequence_length, 171, 125,
        scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        geneconv_fss_table().data());
    result.redo = redo;
    result.candidates.reserve(static_cast<std::size_t>(triplet_count));
    for (int index = 0; index < triplet_count; ++index) {
        // Module8 calls GCXoverD only for RedoL3==1.  A value of 2 is a
        // deferred AddToRedoList entry and must not become an immediate call.
        if (redo[index] != 1) continue;
        result.candidates.push_back({
            scan_state.analysis_list[index * 3],
            scan_state.analysis_list[index * 3 + 1],
            scan_state.analysis_list[index * 3 + 2]});
    }
    return result;
}

RdpMethodScreenResult screen_rdp_maxchi_candidates_impl(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int event_number,
    const std::vector<unsigned char>* supplied_missing_data) {
    RdpMethodScreenResult result;
    const int triplet_count = scan_state.analysis_list_last + 1;
    if (triplet_count <= 0) return result;
    if (scan_state.analysis_list.size() <
        static_cast<std::size_t>(triplet_count) * 3) {
        throw std::runtime_error("MAXCHI analysis list is truncated");
    }
    ensure_geneconv_fss_table();
    auto& fss_mc = maxchi_fss_table();
    std::vector<unsigned char> redo(static_cast<std::size_t>(triplet_count), 0);
    std::vector<unsigned char> worthwhilescan(static_cast<std::size_t>(triplet_count), 1);
    const std::size_t missing_size =
        static_cast<std::size_t>(scan_state.sequence_length + 1) *
        static_cast<std::size_t>(scan_state.next_no + 1);
    std::vector<unsigned char> missing_data(
        supplied_missing_data != nullptr &&
                supplied_missing_data->size() >= missing_size
            ? supplied_missing_data->begin()
            : std::vector<unsigned char>::const_iterator{},
        supplied_missing_data != nullptr &&
                supplied_missing_data->size() >= missing_size
            ? supplied_missing_data->begin() + missing_size
            : std::vector<unsigned char>::const_iterator{});
    if (missing_data.size() != missing_size) missing_data.assign(missing_size, 0);
    const int correction = std::max(1, correction_tests);
    const double sub_threshold = mc_flag != 0
        ? lowest_probability
        : lowest_probability / static_cast<double>(correction);
    constexpr int window_size = 70;
    constexpr int max_window = ChiLookupTable::max_window;
    const int half_window = vb_clng(window_size / 2.0);
    const int critical = critical_difference(window_size, lowest_probability);
    const auto& chi = chi_lookup_table();
    const double window_fraction = 0.1;
    const short* seq_num = scan_state.sequence_data.data();
    const int target_unused = 0;
    (void)target_unused;
    (void)MathFuncs::MyMathFuncs::AlistMC3(
        event_number, worthwhilescan.data(), 0, scan_state.analysis_list_last,
        0, 171, max_window, half_window, half_window, critical, 0,
        scan_state.next_no, store_lpv_ub,
        const_cast<double*>(store_lpv.data()),
        const_cast<short*>(scan_state.analysis_list.data()),
        scan_state.analysis_list_last, redo.data(),
        static_cast<short>(circular), correction, static_cast<short>(mc_flag),
        sub_threshold, lowest_probability, window_fraction, window_size, 0,
        scan_state.sequence_length, scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        125, fss_mc.data(), const_cast<short*>(seq_num), missing_data.data(),
        const_cast<int*>(chi.map.data()),
        const_cast<float*>(chi.values.data()));
    result.redo = redo;
    result.candidates.reserve(static_cast<std::size_t>(triplet_count));
    for (int index = 0; index < triplet_count; ++index) {
        if (redo[index] != 1) continue;
        result.candidates.push_back({
            scan_state.analysis_list[index * 3],
            scan_state.analysis_list[index * 3 + 1],
            scan_state.analysis_list[index * 3 + 2]});
    }
    return result;
}

RdpMethodScreenResult screen_rdp_chimaera_candidates_impl(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int event_number,
    const std::vector<unsigned char>* supplied_missing_data) {
    RdpMethodScreenResult result;
    const int triplet_count = scan_state.analysis_list_last + 1;
    if (triplet_count <= 0) return result;
    if (scan_state.analysis_list.size() <
        static_cast<std::size_t>(triplet_count) * 3) {
        throw std::runtime_error("CHIMAERA analysis list is truncated");
    }
    const int correction = std::max(1, correction_tests);
    const double sub_threshold = mc_flag != 0
        ? lowest_probability
        : lowest_probability / static_cast<double>(correction);
    constexpr int window_size = 70;
    constexpr int max_window = ChiLookupTable::max_window;
    const int half_window = vb_clng(window_size / 2.0);
    const int critical = critical_difference(window_size, lowest_probability);
    const auto& chi = chi_lookup_table();
    std::vector<unsigned char> redo(static_cast<std::size_t>(triplet_count), 0);
    std::vector<unsigned char> worthwhile_scan(
        static_cast<std::size_t>(triplet_count), 127);
    const std::size_t missing_size =
        static_cast<std::size_t>(scan_state.sequence_length + 1) *
        static_cast<std::size_t>(scan_state.next_no + 1);
    std::vector<unsigned char> missing_data(
        supplied_missing_data != nullptr &&
                supplied_missing_data->size() >= missing_size
            ? supplied_missing_data->begin()
            : std::vector<unsigned char>::const_iterator{},
        supplied_missing_data != nullptr &&
                supplied_missing_data->size() >= missing_size
            ? supplied_missing_data->begin() + missing_size
            : std::vector<unsigned char>::const_iterator{});
    if (missing_data.size() != missing_size) missing_data.assign(missing_size, 0);
    const short* seq_num = scan_state.sequence_data.data();
    (void)MathFuncs::MyMathFuncs::AlistChi(
        event_number, missing_data.data(), worthwhile_scan.data(), 0,
        scan_state.analysis_list_last,
        0, 171, max_window, half_window, half_window, critical, 0,
        scan_state.next_no, store_lpv_ub,
        const_cast<double*>(store_lpv.data()),
        const_cast<short*>(scan_state.analysis_list.data()),
        scan_state.analysis_list_last, redo.data(), static_cast<short>(circular),
        correction, static_cast<short>(mc_flag), sub_threshold,
        lowest_probability, 0.1, window_size, 0, scan_state.sequence_length,
        scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        125, chimaera_fss_table().data(), const_cast<short*>(seq_num),
        const_cast<int*>(chi.map.data()), const_cast<float*>(chi.values.data()));
    result.redo = redo;
    result.candidates.reserve(static_cast<std::size_t>(triplet_count));
    for (int index = 0; index < triplet_count; ++index) {
        if (redo[index] != 1) continue;
        result.candidates.push_back({
            scan_state.analysis_list[index * 3],
            scan_state.analysis_list[index * 3 + 1],
            scan_state.analysis_list[index * 3 + 2]});
    }
    return result;
}

int find_subseq_maxchi_compressed(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    std::vector<int>& difference_position,
    std::vector<int>& position_difference) {
    const int stride = scan_state.compressed_sequence_ub + 1;
    const auto decode = [](const int code) {
        return std::array<int, 3>{
            code / 25, (code / 5) % 5, code % 5};
    };
    int informative = 0;
    int site_position = 0;
    for (int packed = 1;
         packed <= scan_state.compressed_sequence_ub; ++packed) {
        const int first_code = scan_state.compressed_sequence[
            packed + sequences[0] * stride];
        const int second_code = scan_state.compressed_sequence[
            packed + sequences[1] * stride];
        const int third_code = scan_state.compressed_sequence[
            packed + sequences[2] * stride];
        const auto first = decode(first_code);
        const auto second = decode(second_code);
        const auto third = decode(third_code);
        for (int site = 0; site < 3; ++site) {
            ++site_position;
            if (first_code < 125 && second_code < 125 && third_code < 125 &&
                first[site] != 0 && second[site] != 0 && third[site] != 0 &&
                (first[site] != third[site] ||
                 first[site] != second[site])) {
                ++informative;
                difference_position[informative] = site_position;
            }
            position_difference[site_position] = informative;
        }
    }
    return informative;
}

std::vector<std::array<int, 3>> make_rdp_inner_method_triplets_impl(
    const RdpScanState& scan_state, const std::array<int, 3>& rnum,
    const std::vector<int>& rlist, const int win_pp,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const std::vector<unsigned char>& do_pairs,
    const int permanent_next_no, const int min_sequence_size,
    const int method_program, const int selected_program_bits,
    const float probability_step) {
    const int pair_stride = static_cast<int>(std::sqrt(
        static_cast<double>(do_pairs.size())));
    if (pair_stride <= 0 ||
        static_cast<std::size_t>(pair_stride) * pair_stride != do_pairs.size()) {
        throw std::runtime_error("MakeAListISP3 pair matrix is not square");
    }
    const int triplet_count = scan_state.analysis_list_last + 1;
    if (triplet_count <= 0) return {};
    if (win_pp < 0 || win_pp > 2 || rnum[win_pp] < -1 ||
        rlist.size() < static_cast<std::size_t>(3) *
            static_cast<std::size_t>(std::max(0, rnum[win_pp] + 1)) ||
        trace_sub.size() < static_cast<std::size_t>(pair_stride) ||
        actual_sequence_sizes.size() < static_cast<std::size_t>(pair_stride)) {
        throw std::runtime_error("MakeAListISP3 input dimensions differ");
    }
    if (method_program < 0 || method_program > 7 ||
        (selected_program_bits & (1 << method_program)) == 0) {
        return {};
    }

    // Feed the vendored DNA5 routine directly.  MakeAListISP3's subtle
    // RestartPos/ProgBinRead/Worthwhilescan/ProbDo behavior is part of the
    // observable method call order and is not safely approximated in a loop.
    constexpr int worthwhile = 127;
    constexpr int prog_bin_upper_bound = 30;
    std::vector<unsigned char> prog_bin_read(
        static_cast<std::size_t>(256) * (prog_bin_upper_bound + 1), 0);
    prog_bin_read[method_program + worthwhile * (prog_bin_upper_bound + 1)] =
        (selected_program_bits & (1 << method_program)) != 0 ? 1 : 0;
    std::vector<unsigned char> worthwhile_scan(
        static_cast<std::size_t>(std::max(0, triplet_count)),
        static_cast<unsigned char>(worthwhile));
    const int max_alist = std::max(0, triplet_count +
        std::max(0, rnum[win_pp]) * 3);
    std::vector<short> alist(
        static_cast<std::size_t>(max_alist + 1) * 3, 0);
    int restart_pos = 0;
    const int alist_last = MathFuncs::MyMathFuncs::MakeAListISP3(
        probability_step, method_program, &restart_pos,
        prog_bin_upper_bound, prog_bin_read.data(),
        const_cast<int*>(trace_sub.data()), win_pp,
        const_cast<int*>(rnum.data()), 2, const_cast<int*>(rlist.data()),
        2, const_cast<short*>(scan_state.analysis_list.data()),
        scan_state.analysis_list_last, worthwhile_scan.data(),
        const_cast<int*>(actual_sequence_sizes.data()), permanent_next_no,
        pair_stride - 1, min_sequence_size, 2, max_alist,
        alist.data(), pair_stride - 1,
        const_cast<unsigned char*>(do_pairs.data()));
    std::vector<std::array<int, 3>> output;
    if (alist_last < 0) return output;
    output.reserve(static_cast<std::size_t>(alist_last + 1));
    for (int index = 0; index <= alist_last; ++index) {
        output.push_back({alist[index * 3], alist[index * 3 + 1],
                          alist[index * 3 + 2]});
    }
    return output;
}

int circular_position(int position, const int length) {
    while (position < 1) position += length;
    while (position > length) position -= length;
    return position;
}

// The 32-bit DNA.dll evaluates this integer-count chi expression in the x87
// register stack before truncating chi*10000 to an int.  Keeping the rational
// arithmetic exact reproduces that observable extended-precision boundary;
// evaluating through a WebAssembly double can round an exact N/10 value just
// below the integer and suppress a strict improvement.
int legacy_chi_score_10000(const int a, const int c, const int window) {
    const int b = window - a;
    const int d = window - c;
    if (a + c <= 0 || b + d <= 0) return 0;
    const std::int64_t cross =
        static_cast<std::int64_t>(a) * d -
        static_cast<std::int64_t>(b) * c;
    const std::int64_t numerator = cross * cross * 20000;
    const std::int64_t denominator =
        static_cast<std::int64_t>(a + c) * (b + d) * window;
    return static_cast<int>(numerator / denominator);
}

// Module3.CentreBP with OBE=OEN=0 and SEventNumber=0, the initial-scan
// branch. It centres each raw informative-site boundary in the flanking
// invariant interval before MCXoverF stores the event.
void centre_initial_breakpoints(
    int& beginning, int& ending, std::vector<int>& position_difference,
    const std::vector<int>& difference_position, const int sequence_length,
    const int informative, const bool circular) {
    const int beginning_index = position_difference[beginning];
    if (beginning_index - 1 > 0) {
        beginning -= vb_clng(
            (static_cast<double>(
                beginning - difference_position[beginning_index - 1]) /
             2.0) - 0.1);
    } else {
        beginning -= vb_clng(
            (static_cast<double>(
                beginning + sequence_length -
                difference_position[informative]) /
             2.0) - 0.1);
    }
    if (beginning == 0) beginning = 1;
    else if (beginning < 1)
        beginning = circular ? sequence_length + beginning : 1;

    position_difference[sequence_length] = informative;
    const int ending_index = position_difference[ending];
    if (ending_index + 1 <= informative) {
        ending += vb_clng(
            (static_cast<double>(
                difference_position[ending_index + 1] - ending) /
             2.0) - 0.1);
    } else {
        ending += vb_clng(
            (static_cast<double>(
                difference_position[1] + sequence_length - ending) /
             2.0) - 0.1);
    }
    if (ending > sequence_length)
        ending = circular ? ending - sequence_length : sequence_length;
}

// Module2.FindDaughter, active PrgF=3/4 endpoint-polishing branch.  Despite
// its name, this branch does not alter the three roles; it only walks the
// informative-site map and conditionally tightens the stored breakpoints.
void legacy_polish_maxchi_breakpoints(
    const int comparison, const bool low_high,
    const std::array<int, 3>& sequences, const int sequence_length,
    const int informative, const bool circular,
    const std::vector<int>& position_difference,
    const std::vector<int>& difference_position,
    const std::vector<unsigned char>& missing_map,
    const std::vector<short>& sequence_data,
    int& beginning, int& ending, MaxchiPeakTrace* polish_trace) {
    int high_first = 0;
    int high_second = 0;
    if (comparison == 0) {
        high_first = sequences[0];
        high_second = sequences[1];
    } else if (comparison == 1) {
        high_first = sequences[0];
        high_second = sequences[2];
    } else {
        high_first = sequences[2];
        high_second = sequences[1];
    }
    const int stride = sequence_length + 1;
    const auto equal_at = [&](const int informative_position) {
        const int site = difference_position[informative_position];
        return sequence_data[site + high_first * stride] ==
            sequence_data[site + high_second * stride];
    };
    const int original_beginning = beginning;
    const int original_ending = ending;
    int region_start = beginning;
    int region_end = ending;
    int polished_beginning = beginning;
    int polished_ending = ending;
    int x = position_difference[region_start];
    if (x == 0) x = 1;
    if (polish_trace) polish_trace->polish_begin_index = x;
    for (int guard = 0; guard <= informative * 2 + 2; ++guard) {
        if (low_high ? !equal_at(x) : equal_at(x)) break;
        --x;
        if (x >= 0 && x < static_cast<int>(missing_map.size()) &&
            missing_map[x] == 1) break;
        if (x < 1) {
            if (!circular) break;
            if (x == position_difference[region_end]) { ++x; break; }
            x = informative;
        }
        if (x == position_difference[region_end]) { ++x; break; }
    }
    if (polish_trace) polish_trace->polish_left_stop_index = x;
    region_start = difference_position[x];
    x = position_difference[region_start];
    for (int guard = 0; guard <= informative * 2 + 2; ++guard) {
        if (low_high ? equal_at(x) : !equal_at(x)) {
            polished_beginning = difference_position[x];
            if (polish_trace) {
                polish_trace->polish_begin_candidate_index = x;
                polish_trace->polish_begin_candidate = polished_beginning;
            }
            break;
        }
        ++x;
        if (x >= 0 && x < static_cast<int>(missing_map.size()) &&
            missing_map[x] == 1) break;
        if (x > informative) {
            if (!circular) break;
            x = 1;
        }
        if (x == position_difference[region_end]) { --x; break; }
    }

    x = position_difference[region_end];
    if (x >= static_cast<int>(missing_map.size()))
        x = static_cast<int>(missing_map.size()) - 1;
    if (missing_map[x] == 0) {
        for (int guard = 0; guard <= informative * 2 + 2; ++guard) {
            if (low_high ? !equal_at(x) : equal_at(x)) break;
            ++x;
            if (x >= 0 && x < static_cast<int>(missing_map.size()) &&
                missing_map[x] == 1) break;
            if (x > informative) {
                if (!circular) break;
                x = 1;
                if (missing_map[x] == 1) break;
            }
            if (x == position_difference[region_start]) { --x; break; }
        }
    }
    region_end = difference_position[x];
    x = position_difference[region_end];
    if (x == 0) x = 1;
    for (int guard = 0; guard <= informative * 2 + 2; ++guard) {
        if (low_high ? equal_at(x) : !equal_at(x)) {
            polished_ending = difference_position[x];
            break;
        }
        if (x > 0 && x < static_cast<int>(missing_map.size()) &&
            missing_map[x] == 1) break;
        --x;
        if (x >= 0 && x < static_cast<int>(missing_map.size()) &&
            missing_map[x] == 1) break;
        if (x < 1) {
            if (!circular) break;
            x = informative;
        }
        if (x == position_difference[region_start]) { ++x; break; }
    }

    const auto near_original = [&](const int candidate, const int original) {
        const int candidate_position = position_difference[candidate];
        const int original_position = position_difference[original];
        return std::abs(candidate_position - original_position) < 5 ||
            std::abs(candidate_position -
                     (original_position - sequence_length)) < 5 ||
            std::abs((candidate_position - sequence_length) -
                     original_position) < 5;
    };
    if (near_original(polished_beginning, original_beginning)) {
        beginning = polished_beginning;
        if (polish_trace) polish_trace->polish_begin_accepted = true;
    }
    if (near_original(polished_ending, original_ending))
        ending = polished_ending;
}

int event_half_window(const int event_beginning, const int event_ending,
                      const int informative, const int critical,
                      const int prior_half_window,
                      const std::vector<int>& position_difference) {
    int half_window = 0;
    if (event_beginning < event_ending) {
        half_window = position_difference[event_ending] -
            position_difference[event_beginning] + 1;
    } else {
        half_window = position_difference[event_ending] +
            (informative - position_difference[event_beginning]) + 1;
    }
    if (half_window * 2 > informative) {
        half_window = vb_clng(informative * 0.75 / 2.0 + 0.00001) - 1;
    }
    if (half_window <= critical) {
        half_window = vb_clng(informative / 2.0 + 0.00001) - 1;
    }
    return half_window < 6 ? -prior_half_window : half_window;
}

// threshold.cpp:FindSide from the DNA.dll source called by MCXoverF/CXoverA.
void legacy_find_side(const int top_left, const int top_right,
                      const int sequence_length, const int left,
                      const int right, const int window,
                      const int informative, const int comparison,
                      const std::vector<unsigned char>& scores,
                      double& high_left, double& high_right) {
    const int stride = sequence_length + 1;
    const int offset = comparison * stride;
    const auto score = [&](const int position) {
        return static_cast<int>(
            scores[circular_position(position, informative) + offset]);
    };
    int outside = 0;
    for (int position = left - window; position < left; ++position)
        outside += score(position);
    auto chi = [&](const int first, const int second) -> double {
        const int b = window - first;
        const int d = window - second;
        if (first + second <= 0 || b + d <= 0) return 0.0;
        const double cross = static_cast<double>(first) * d -
                             static_cast<double>(b) * second;
        return cross * cross * 2.0 /
            (static_cast<double>(window) * (first + second) * (b + d));
    };
    high_left = chi(top_left, outside);
    outside = 0;
    for (int position = right + 1; position <= right + window; ++position)
        outside += score(position);
    high_right = chi(top_right, outside);
}

int legacy_opt_left(const int left, const double high_left, const int top_left,
                    const int maximum, const int comparison, const int window,
                    const int informative, const int sequence_length,
                    const std::vector<unsigned char>& scores,
                    const std::vector<unsigned char>& missing_map,
                    MaxchiPeakTrace* trace) {
    const int offset = comparison * (sequence_length + 1);
    int outside = 0;
    for (int x = left - window; x <= left - 1; ++x) {
        int position = x;
        if (position < 1) position = informative + position;
        else if (position > informative) position -= informative;
        outside += scores[position + offset];
    }
    int p = maximum - 1;
    int m = left - 1;
    int n = left - window;
    if (p < 1) p = informative + p;
    if (n < 1) n = informative + n;
    if (m < 1) m = informative + m;
    double a = outside;
    double c = top_left;
    int best = static_cast<int>(high_left * 10000.0);
    if (trace) trace->opt_initial_score = best;
    int best_position = m;
    int failures = 0;
    while (true) {
        a -= scores[m + offset];
        c -= scores[p + offset];
        c += scores[m + offset];
        --p; --n; --m;
        if (missing_map[p] || missing_map[n] || missing_map[m]) break;
        if (p < 1) p = informative;
        if (n < 1) n = informative;
        if (m < 1) m = informative;
        if (missing_map[p] || missing_map[n] || missing_map[m]) break;
        a += scores[n + offset];
        const double b = window - a;
        const double d = window - c;
        if (a + c > 0 && b + d > 0) {
            const int candidate = legacy_chi_score_10000(
                static_cast<int>(a), static_cast<int>(c), window);
            if (trace && trace->opt_first_position == 0) {
                trace->opt_first_position = m;
                trace->opt_first_score = candidate;
                trace->opt_first_improved = trace->opt_first_score > best;
            }
            if (candidate > best) {
                best = candidate;
                best_position = m;
                failures = 0;
            } else if (++failures > window / 10) break;
        }
    }
    return best_position + 1;
}

int legacy_opt_right(const int right, const double high_right,
                     const int top_right, const int maximum,
                     const int comparison, const int window,
                     const int informative, const int sequence_length,
                     const std::vector<unsigned char>& scores,
                     const std::vector<unsigned char>& missing_map) {
    const int offset = comparison * (sequence_length + 1);
    int outside = 0;
    for (int x = right + 1; x <= right + window; ++x) {
        int position = x;
        if (position < 1) position = informative + position;
        else if (position > informative) position -= informative;
        outside += scores[position + offset];
    }
    int p = right + window;
    int m = right + 1;
    int n = maximum;
    if (p > informative) p -= informative;
    if (n > informative) n -= informative;
    if (m > informative) m -= informative;
    double a = top_right;
    double c = outside;
    int best = static_cast<int>(high_right * 10000.0);
    int best_position = m;
    int failures = 0;
    while (true) {
        a -= scores[n + offset];
        c -= scores[m + offset];
        a += scores[m + offset];
        ++p; ++n; ++m;
        if (missing_map[p] || missing_map[n] || missing_map[m]) break;
        if (p > informative) p = 1;
        if (n > informative) n = 1;
        if (m > informative) m = 1;
        if (missing_map[p] || missing_map[n] || missing_map[m]) break;
        c += scores[p + offset];
        const double b = window - a;
        const double d = window - c;
        if (a + c > 0 && b + d > 0) {
            const int candidate = legacy_chi_score_10000(
                static_cast<int>(a), static_cast<int>(c), window);
            if (candidate > best) {
                best = candidate;
                best_position = m;
                failures = 0;
            } else if (++failures > window / 10) break;
        }
    }
    return best_position - 1;
}

// threshold.cpp:DestroyPeaks. This function's broad peak deletion is what
// makes the enumerating routines differ from the DLL fast-check wrappers.
void legacy_destroy_peaks(const int comparison, const int informative,
                          const int length, int left, int right,
                          const std::vector<double>& smooth,
                          std::vector<double>& chi_values) {
    const int offset = comparison * (length + 1);
    if (left == 0) left = 1;
    if (right == 0) right = 1;
    int start_left = left;
    int circuits = 0;
    bool erase_all = false;
    if (informative > 5) {
        while (smooth[left + offset] >= smooth[left + offset + 1] ||
               smooth[left + offset] >= smooth[left + offset + 2]) {
            --left;
            if (left == right) { erase_all = true; break; }
            if (left < 1) {
                if (++circuits == 2) { erase_all = true; break; }
                left = informative + left;
            }
            if (left == right) { erase_all = true; break; }
            if (left > informative - 2) left = informative - 2;
        }
        if (left == right) erase_all = true;
        if (!erase_all) {
            circuits = 0;
            while ((smooth[left + offset] <= smooth[left + offset + 1] ||
                    smooth[left + offset] <= smooth[left + offset + 2]) &&
                   smooth[left + offset] > 1) {
                --left;
                if (left == right) { erase_all = true; break; }
                if (left < 1) {
                    if (++circuits == 2) { erase_all = true; break; }
                    left = informative + left;
                }
                if (left == right) { erase_all = true; break; }
            }
        }
        if (left == right) erase_all = true;
        if (!erase_all) {
            if (right <= 0) right = informative + right;
            while (true) {
                ++right;
                if (left == right || right == start_left) {
                    erase_all = true; break;
                }
                if (right > informative - 1) {
                    right = 1;
                    if (right > left) { erase_all = true; break; }
                }
                if (left == right || right == start_left) {
                    erase_all = true; break;
                }
                if (smooth[right + offset] < smooth[right - 1 + offset]) break;
            }
        }
        if (!erase_all) {
            while (true) {
                ++right;
                if (left == right || right == start_left) {
                    erase_all = true; break;
                }
                if (right > informative - 1) {
                    right = 1;
                    if (right > left) { erase_all = true; break; }
                }
                if (smooth[right + offset] == 0 ||
                    smooth[right + offset] > smooth[right - 1 + offset]) break;
                if (left == right || right == start_left) {
                    erase_all = true; break;
                }
            }
        }
    } else {
        erase_all = true;
    }
    if (erase_all) {
        for (int x = 0; x <= informative; ++x) chi_values[x + offset] = 0;
    } else if (left < right) {
        for (int x = left; x <= right; ++x) chi_values[x + offset] = 0;
    } else {
        for (int x = 0; x <= right; ++x) chi_values[x + offset] = 0;
        for (int x = left; x <= informative; ++x) chi_values[x + offset] = 0;
    }
}

void fill_legacy_event(RdpLegacyEventAllocator& allocator,
                       const std::array<int, 3>& active, const int program,
                       const double probability, const int beginning,
                       const int ending, const int window_width,
                       const bool companion) {
    const int slot = allocator.allocate(active[0], program, probability);
    if (slot < 1) return;
    auto& event = allocator.event(active[0], slot);
    event.daughter = static_cast<std::int16_t>(active[0]);
    event.minor_parent = static_cast<std::int16_t>(active[1]);
    event.major_parent = static_cast<std::int16_t>(active[2]);
    event.beginning = beginning;
    event.ending = ending;
    event.program_flag = static_cast<std::uint8_t>(program);
    event.probability = probability;
    event.length_holder = window_width * 2;
    if (!companion) return;
    // allocate() may grow the row and invalidate references into its vector.
    // Copy the completed source record before requesting the companion slot.
    auto reverse = event;
    const int reverse_slot = allocator.allocate(active[0], program, probability);
    if (reverse_slot > 0) {
        std::swap(reverse.beginning, reverse.ending);
        allocator.event(active[0], reverse_slot) = reverse;
    }
}

}  // namespace

std::vector<std::array<int, 3>> make_rdp_inner_method_triplets(
    const RdpScanState& scan_state, const std::array<int, 3>& rnum,
    const std::vector<int>& rlist, const int win_pp,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const std::vector<unsigned char>& do_pairs,
    const int permanent_next_no, const int min_sequence_size,
    const int method_program, const int selected_program_bits,
    const float probability_step) {
    return make_rdp_inner_method_triplets_impl(
        scan_state, rnum, rlist, win_pp, trace_sub, actual_sequence_sizes,
        do_pairs, permanent_next_no, min_sequence_size, method_program,
        selected_program_bits, probability_step);
}

RdpMethodScreenResult screen_rdp_geneconv_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int target) {
    return screen_rdp_geneconv_candidates_impl(
        scan_state, store_lpv, store_lpv_ub, correction_tests,
        lowest_probability, circular, mc_flag, target);
}

RdpMethodScreenResult screen_rdp_maxchi_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int event_number,
    const std::vector<unsigned char>* missing_data) {
    return screen_rdp_maxchi_candidates_impl(
        scan_state, store_lpv, store_lpv_ub, correction_tests,
        lowest_probability, circular, mc_flag, event_number, missing_data);
}

RdpMethodScreenResult screen_rdp_chimaera_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    const int store_lpv_ub, const int correction_tests,
    const double lowest_probability, const int circular, const int mc_flag,
    const int event_number,
    const std::vector<unsigned char>* missing_data) {
    return screen_rdp_chimaera_candidates_impl(
        scan_state, store_lpv, store_lpv_ub, correction_tests,
        lowest_probability, circular, mc_flag, event_number, missing_data);
}

RdpLegacyEventAllocator::RdpLegacyEventAllocator(
    RdpRawEventState& events, const int selected_program_count,
    const int max_events)
    : events_(events), selected_program_count_(selected_program_count),
      max_events_(max_events),
      max_xop_(static_cast<std::size_t>(method_count_) *
               events.xover_list.size(), 0) {
    if (events_.current_xover.size() != events_.xover_list.size() ||
        selected_program_count_ < 1) {
        throw std::runtime_error("UpdateXOList3 state dimensions differ");
    }
    // FinalTrim starts with ReDim XOverList(NextNo, 10), then earlier method
    // calls may already have populated it. Reproduce their MaxXOP counters.
    for (std::size_t sequence = 0; sequence < events_.xover_list.size();
         ++sequence) {
        events_.current_xover[sequence] = static_cast<std::int16_t>(
            events_.xover_list[sequence].size());
        for (const auto& event : events_.xover_list[sequence]) {
            if (event.program_flag < method_count_) {
                ++max_xop_[event.program_flag + sequence * method_count_];
            }
        }
    }
}

int RdpLegacyEventAllocator::allocate(
    const int active_sequence, const int program, const double probability) {
    (void)probability;
    if (active_sequence < 0 ||
        active_sequence >= static_cast<int>(events_.xover_list.size()) ||
        program < 0 || program >= method_count_) {
        throw std::runtime_error("UpdateXOList3 index differs");
    }
    const auto max_index = static_cast<std::size_t>(program) +
        static_cast<std::size_t>(active_sequence) * method_count_;
    // PN is the number of selected programs. The RDP-only run has PN=1.
    const int per_program_limit = static_cast<int>(
        static_cast<double>(max_events_) / selected_program_count_ - 0.49 + 0.5);
    if (max_xop_[max_index] >= per_program_limit ||
        count(active_sequence) >= max_events_) {
        return -1;
    }
    const int slot = count(active_sequence) + 1;
    events_.current_xover[active_sequence] =
        static_cast<std::int16_t>(slot);
    if (slot > upper_bound_ && upper_bound_ <= max_events_) {
        // VB6's untyped 1.1 literal is a Single expression here.
        upper_bound_ = static_cast<int>(
            static_cast<float>(upper_bound_) * 1.1F) + 1;
        upper_bound_ = std::min(upper_bound_, max_events_);
        if (upper_bound_ < slot) upper_bound_ = slot + 3;
    }
    auto& row = events_.xover_list[active_sequence];
    if (row.size() < static_cast<std::size_t>(slot)) row.resize(slot);
    ++max_xop_[max_index];
    return slot;
}

bool RdpLegacyEventAllocator::has_strictly_later_slot(const int slot) const {
    return upper_bound_ > slot;
}

RdpRawEvent& RdpLegacyEventAllocator::event(
    const int sequence, const int slot) {
    if (sequence < 0 || sequence >= static_cast<int>(events_.xover_list.size()) ||
        slot < 1 || slot > static_cast<int>(events_.xover_list[sequence].size())) {
        throw std::runtime_error("XOverList slot differs");
    }
    return events_.xover_list[sequence][slot - 1];
}

int RdpLegacyEventAllocator::count(const int sequence) const {
    return events_.current_xover[sequence];
}

void run_rdp_geneconv_recheck(
    const RdpScanState& scan_state, std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, const int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, const bool long_winded,
    std::vector<GeneconvEmissionTrace>* trace) {
    constexpr short mismatch_penalty = 1;
    // The command-line scan resets GCMaxOverlapFrags to its default of one
    // for the exploratory pass used by these fixtures.  DelPVals compares
    // the overlap counter inclusively, so this value is observable in both
    // event order and the final cyclic selection.
    constexpr short max_overlap_fragments = 1;
    const int length = scan_state.sequence_length;
    // GCXoverD receives the global GCDimSize separately from the alignment
    // length and passes that upper bound to every fragment/probability DLL
    // routine.  The pinned command-line oracle initializes it to 2999; using
    // the alignment length here changes the VB column stride and therefore
    // which peaks DelPValsP removes on longer alignments.
    constexpr int dimension = 2999;
    const int stride = dimension + 1;
    const int subsequence_stride = length + 1;
    std::vector<char> subsequence(
        static_cast<std::size_t>(subsequence_stride) * 7, 0);
    // Module2 dimensions both maps as Len(StrainSeq(0)) + 200.  The
    // compressed routine expands each packed byte to three sites and can
    // write through the final padded triplet beyond the alignment endpoint.
    std::vector<int> difference_position(length + 201, 0);
    std::vector<int> position_difference(length + 201, 0);
    std::array<int, 8> differences{};
    auto& fss_gc = geneconv_fss_table();
    populate_geneconv_fss_triplet(scan_state, sequences, fss_gc);
    const int informative = MathFuncs::MyMathFuncs::FindSubSeqGCAP7(
        scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        125, const_cast<unsigned char*>(fss_gc.data()), 0, length,
        sequences[0], sequences[1], sequences[2], subsequence.data(),
        differences.data(), difference_position.data(),
        position_difference.data());
    if (differences[0] == informative || differences[1] == informative ||
        differences[2] == informative || informative < 1) {
        return;
    }
    subsequence[informative + 1] = 0;
    subsequence[informative + 1 + (length + 1)] = 0;
    subsequence[informative + 1 + 2 * (length + 1)] = 0;
    subsequence[informative + 1 + 6 * (length + 1)] = 0;
    differences[3] = differences[0] + differences[1];
    differences[4] = differences[0] + differences[2];
    differences[5] = differences[1] + differences[2];
    differences[0] = informative - differences[0];
    differences[1] = informative - differences[1];
    differences[2] = informative - differences[2];
    const float total = static_cast<float>(informative * mismatch_penalty);
    std::array<double, 8> miss_penalty{};
    for (int comparison = 0; comparison < 3; ++comparison) {
        if (differences[comparison] == 0) return;
        miss_penalty[comparison] =
            static_cast<int>(total / differences[comparison]) + 1;
    }
    int minimum = length;
    int maximum = 0;
    for (int comparison = 0; comparison < 6; ++comparison) {
        minimum = std::min(minimum, differences[comparison]);
        maximum = std::max(maximum, differences[comparison]);
    }
    if (minimum < 3 && maximum > minimum * 10) return;
    for (int comparison = 3; comparison < 6; ++comparison) {
        if (differences[comparison] == 0) differences[comparison] = 1;
        miss_penalty[comparison] =
            static_cast<int>(total / differences[comparison]) + 1;
    }

    std::vector<int> fragment_start(static_cast<std::size_t>(stride) * 7, 0);
    std::vector<int> fragment_end(static_cast<std::size_t>(stride) * 7, 0);
    std::vector<int> fragment_score(static_cast<std::size_t>(stride) * 7, 0);
    std::array<int, 8> fragment_count{};
    if (MathFuncs::MyMathFuncs::GetFragsP(
            static_cast<short>(settings.circular), informative, length,
            dimension, subsequence.data(), fragment_start.data(),
            fragment_end.data(), fragment_score.data(),
            fragment_count.data()) == 0) {
        return;
    }
    std::vector<int> max_score_position(static_cast<std::size_t>(stride) * 6, 0);
    std::vector<int> fragment_max_score(static_cast<std::size_t>(stride) * 6, 0);
    std::array<int, 8> high_score{};
    MathFuncs::MyMathFuncs::GetMaxFragScoreP(
        informative, dimension, static_cast<short>(settings.circular),
        mismatch_penalty, miss_penalty.data(), max_score_position.data(),
        fragment_max_score.data(), fragment_score.data(),
        fragment_count.data(), high_score.data());
    std::array<int, 8> high_enough{};
    for (int comparison = 0; comparison < 6; ++comparison) {
        high_enough[comparison] = high_score[comparison] > 3 ? 1 : 0;
    }
    std::array<double, 8> critical{};
    std::array<double, 8> ll{};
    std::array<double, 8> kmax{};
    double cutoff = 0.0;
    if (MathFuncs::MyMathFuncs::CalcKMaxP(
            mismatch_penalty, informative, static_cast<short>(settings.mc_flag),
            settings.mc_correction, settings.lowest_probability, &cutoff,
            high_score.data(), critical.data(), miss_penalty.data(), ll.data(),
            kmax.data(), differences.data(), high_enough.data()) == 0) {
        return;
    }
    std::vector<double> pvalues(static_cast<std::size_t>(stride) * 6, 0.0);
    const double maximum_score = MathFuncs::MyMathFuncs::GCCalcPValP2(
        dimension, informative, fragment_max_score.data(), pvalues.data(),
        fragment_count.data(), kmax.data(), ll.data(), high_enough.data(),
        critical.data());
    if (maximum_score > cutoff) return;
    // GCXoverD does not apply the RDP StoreLPV pre-filter.  Its only global
    // threshold check is MaxScore > PCO; StoreLPV is used for neither the
    // initial exploratory pass nor the LongWinded role choice.
    std::vector<int> deleted(informative + 2, 0);
    double prior = 0.0;
    int repetitions = 0;
    while (true) {
        int y = -1;
        int comparison = -1;
        const double pvalue = MathFuncs::MyMathFuncs::GCGetHiPValP(
            dimension, informative, fragment_count.data(), pvalues.data(),
            &y, &comparison, high_enough.data());
        ++repetitions;
        if ((pvalue == prior && repetitions > 10) || pvalue > cutoff ||
            pvalue == 1.0 || y < 0 || comparison < 0) {
            return;
        }
        prior = pvalue;
        const int deletion_accepted = MathFuncs::MyMathFuncs::DelPValsP(
                max_overlap_fragments, y, comparison, dimension,
                pvalues.data(), fragment_count.data(), fragment_start.data(),
                fragment_end.data(), max_score_position.data(),
                deleted.data());
        if (deletion_accepted == 1 &&
            pvalues[y + comparison * stride] < cutoff) {
            const std::array<int, 3> counts{
                allocator.count(sequences[0]),
                allocator.count(sequences[1]),
                allocator.count(sequences[2]),
            };
            std::array<int, 3> active{};
            if (long_winded) {
                active = choose_active(
                    sequences, 1, store_lpv, store_lpv_ub, allocator);
            } else if (comparison == 0 || comparison == 3) {
                active = {sequences[0], sequences[1], sequences[2]};
            } else if (comparison == 1) {
                active = {sequences[0], sequences[2], sequences[1]};
            } else if (comparison == 2) {
                active = {sequences[1], sequences[2], sequences[0]};
            } else if (comparison == 4) {
                active = {sequences[1], sequences[0], sequences[2]};
            } else {
                active = {sequences[2], sequences[1], sequences[0]};
            }
            const double adjusted = settings.mc_flag == 0
                ? pvalues[y + comparison * stride] * settings.mc_correction
                : pvalues[y + comparison * stride];
            if (long_winded && adjusted > store_probability(
                    store_lpv, store_lpv_ub, 1, active[0])) {
                return;
            }
            const int slot = allocator.allocate(active[0], 1, adjusted);
            if (slot > 0) {
                auto& event = allocator.event(active[0], slot);
                event.daughter = static_cast<std::int16_t>(active[0]);
                event.minor_parent = static_cast<std::int16_t>(active[1]);
                event.major_parent = static_cast<std::int16_t>(active[2]);
                const int fragment_begin =
                    fragment_start[y + comparison * stride];
                const int fragment_finish = fragment_end[
                    max_score_position[y + comparison * stride] +
                    comparison * stride];
                if (long_winded) {
                    event.beginning = difference_position[fragment_begin];
                    event.ending = difference_position[fragment_finish];
                } else {
                    event.beginning = fragment_begin > 0
                        ? difference_position[fragment_begin - 1] + 1 : 1;
                    event.ending = difference_position[fragment_finish + 1] - 1;
                }
                event.outside_flag = legacy_qc_topology_change(
                    event.beginning, event.ending, position_difference,
                    subsequence, subsequence_stride, informative);
                // CheckEndsVB and FixEnds precede CentreBP in GCXoverD.  At
                // this first MakeTestPVs boundary SEventNumber is zero and
                // MissingData is the zero-filled initial matrix, hence both
                // warnings and FixEnds are no-ops.
                legacy_centre_geneconv_breakpoints(
                    length, settings.circular, informative,
                    position_difference, difference_position,
                    event.beginning, event.ending);
                event.program_flag = 1;
                event.probability = adjusted;
                if (trace != nullptr) {
                    trace->push_back(GeneconvEmissionTrace{
                        sequences,
                        counts,
                        {
                            store_probability(store_lpv, store_lpv_ub, 1,
                                              sequences[0]),
                            store_probability(store_lpv, store_lpv_ub, 1,
                                              sequences[1]),
                            store_probability(store_lpv, store_lpv_ub, 1,
                                              sequences[2]),
                        },
                        active,
                        event,
                    });
                }
                // The current GCXoverD source restores ActiveSeq/Minor/Major
                // after each peak and does not rewrite Seq1/Seq2/Seq3.
            }
            // This call is inside both `If GoOn = 1` and
            // `If PVals(Y, X) < PCO` in GCXoverD.  Rejected overlapping
            // candidates are set to 100 by DelPValsP but must not expand the
            // occupied-region map themselves.
            MathFuncs::MyMathFuncs::MakeDeleteArrayP(
                fragment_start[y + comparison * stride],
                fragment_end[max_score_position[y + comparison * stride] +
                             comparison * stride],
                fragment_count[comparison], deleted.data());
        }
    }
}

void run_rdp_maxchi_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, const int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, const int event_beginning,
    const int event_ending, const bool initial_scan,
    std::vector<MaxchiPeakTrace>* trace, const bool inner_scan,
    const std::vector<unsigned char>* missing_data) {
    constexpr int window_size = 70;
    constexpr int max_window = ChiLookupTable::max_window;
    const int length = scan_state.sequence_length;
    int half_window = vb_clng(window_size / 2.0);
    int critical = critical_difference(window_size, settings.lowest_probability);
    if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
        std::cerr << "temp-maxchi-settings half=" << half_window
                  << " critical=" << critical << " initial=" << initial_scan
                  << " event=" << event_beginning << ':' << event_ending << '\n';
    }
    std::vector<int> difference_position(length + 201, 0);
    std::vector<int> position_difference(length + 201, 0);
    // MCXoverF uses the packed FindSubSeqMCPB path only for the initial
    // SEventNumber=0 scan.  FinalTrim's recheck (SEventNumber>0) switches to
    // the ordinary FindSubSeqCP path; using compressed triplets here silently
    // changes the informative-site map on gapped alignments.
    const int informative = initial_scan
        ? find_subseq_maxchi_compressed(
              scan_state, sequences, difference_position,
              position_difference)
        : MathFuncs::MyMathFuncs::FindSubSeqCP(
              scan_state.sequence_length + 1,
              static_cast<short>(scan_state.next_no),
              static_cast<short>(sequences[0]),
              static_cast<short>(sequences[1]),
              static_cast<short>(sequences[2]),
              const_cast<short*>(scan_state.sequence_data.data()),
              difference_position.data(), position_difference.data());
    const std::uint64_t difference_position_hash = fnv1a64(
        reinterpret_cast<const unsigned char*>(difference_position.data()),
        static_cast<std::size_t>(informative + 1) * sizeof(int));
    const std::uint64_t position_difference_hash = fnv1a64(
        reinterpret_cast<const unsigned char*>(position_difference.data()),
        static_cast<std::size_t>(
            scan_state.compressed_sequence_ub * 3 + 1) * sizeof(int));
    if (informative < critical * 2 || informative < 7) return;
    if (initial_scan) {
        // MakeWindowSize takes the configured fixed width when BEP=ENP=0,
        // including the exploratory FindAllFlag=1 path.
        if (half_window * 2 > informative)
            half_window = vb_clng(informative * 0.75 / 2.0 + 0.00001) - 1;
        if (half_window <= critical)
            half_window = vb_clng(informative / 2.0 + 0.00001) - 1;
        if (half_window < 6) return;
    } else if (inner_scan) {
        // InnerScan2 calls MCXoverF with FindAllFlag=0 and BEP=ENP=0.
        // MakeWindowSize therefore uses the configured fixed MC window; the
        // event-specific window is only used by FinalTrim's FindAllFlag=1
        // recheck.
        const int go_on = MathFuncs::MyMathFuncs::MakeWindowSizeP(
            0, 0, &critical, informative, 0.1, window_size, &half_window,
            half_window, 0);
        if (go_on == 0) return;
    } else {
        const int go_on = MathFuncs::MyMathFuncs::MakeWindowSizeP(
            event_beginning, event_ending, &critical, informative, 0.1,
            window_size, &half_window, half_window, 0);
        if (go_on == 0) return;
    }
    if (half_window < 0) return;
    const int win_upper = length + half_window * 2;
    const int win_stride = win_upper + 1;
    const int sequence_stride = length + 1;
    std::vector<unsigned char> scores(
        static_cast<std::size_t>(sequence_stride) * 3, 0);
    std::vector<int> window_scores(
        static_cast<std::size_t>(win_stride) * 3, 0);
    std::vector<double> chi_values(
        static_cast<std::size_t>(sequence_stride) * 3, 0.0);
    std::vector<double> smooth(
        static_cast<std::size_t>(sequence_stride) * 3, 0.0);
    // MCXoverF allocates BanWin through the full alignment length plus both
    // half-window extensions.  The extension is accessed by CalcChiVals4P3
    // even when only the informative prefix is populated.
    std::vector<int> banned_windows(length + half_window * 2 + 2, 0);
    std::vector<unsigned char> missing_map(length + 3, 0);
    if (inner_scan) {
        // InnerScan2 calls MCXoverF with SEventNumber > 0 and FindAllFlag=0.
        // Its source rebuilds the informative-coordinate missing-data map on
        // every call; this is not equivalent to marking raw alignment sites.
        // MakeBanWinP consumes XPosDiff (informative -> alignment) and
        // XDiffPos (alignment -> informative), exactly as FindSubSeqCP fills
        // them above.
        if (missing_data != nullptr &&
            missing_data->size() >=
                static_cast<std::size_t>(length + 1) *
                    static_cast<std::size_t>(scan_state.next_no + 1)) {
            MathFuncs::MyMathFuncs::MakeBanWinP(
                length + half_window * 2, sequences[0], sequences[1],
                sequences[2], half_window, length, informative,
                banned_windows.data(), missing_map.data(),
                const_cast<unsigned char*>(missing_data->data()),
                position_difference.data(), difference_position.data());
        } else {
            std::fill(banned_windows.begin(), banned_windows.end(), 0);
            std::fill(missing_map.begin(), missing_map.end(), 0);
        }
    }
    if (inner_scan && settings.circular == 0) {
        // MCXoverF also excludes windows crossing the ends of a linear
        // alignment after the missing-data map has been prepared.
        missing_map[1] = 1;
        missing_map[informative] = 1;
        for (int position = informative - half_window + 2;
             position <= informative; ++position) {
            if (position >= 0 && position < static_cast<int>(banned_windows.size()))
                banned_windows[position] = 1;
        }
    }
    MathFuncs::MyMathFuncs::WinScoreCalcP(
        win_upper, critical, half_window, informative, length + 1,
        sequences[0], sequences[1], sequences[2], scores.data(),
        difference_position.data(),
        const_cast<short*>(scan_state.sequence_data.data()),
        window_scores.data());
    const std::uint64_t score_hash = fnv1a64(scores.data(), scores.size());
    std::array<std::uint64_t, 3> score_plane_hash{};
    for (int plane = 0; plane < 3; ++plane) {
        score_plane_hash[plane] = fnv1a64(
            scores.data() + plane * sequence_stride, informative + 1);
    }
    const auto& table = chi_lookup_table();
    const std::uint64_t chi_table_hash = fnv1a64(
        reinterpret_cast<const unsigned char*>(table.values.data()),
        table.values.size() * sizeof(float));
    const std::uint64_t chi_map_hash = fnv1a64(
        reinterpret_cast<const unsigned char*>(table.map.data()),
        table.map.size() * sizeof(int));
    // MCXoverF's FindAllFlag=1 path (the FinalTrim recheck) always uses the
    // formula CalcChiValsP branch, even for a linear alignment.  The lookup
    // table is selected only by the FindAllFlag=0/CircularFlag=0 branch.
    // `initial_scan` denotes that latter first-pass path here.
    double maximum = (initial_scan || inner_scan) && settings.circular == 0
        ? MathFuncs::MyMathFuncs::CalcChiVals4P3(
              win_upper, critical, half_window, informative, length,
              window_scores.data(), chi_values.data(), banned_windows.data(),
              const_cast<float*>(
                  table.values.data() + table.map[half_window]))
        : MathFuncs::MyMathFuncs::CalcChiValsP(
              win_upper, critical, half_window, informative, length,
              window_scores.data(), chi_values.data());
    if (MathFuncs::MyMathFuncs::ChiPVal2P(maximum) *
            (static_cast<double>(informative) / half_window) * 3.0 >
        settings.lowest_probability) return;
    MathFuncs::MyMathFuncs::SmoothChiValsP(
        informative, length, chi_values.data(), smooth.data());

    int wasted_peaks = 0;
    for (int repetition = 1; repetition <= 100; ++repetition) {
        int maximum_position = -1;
        short comparison = -1;
        MathFuncs::MyMathFuncs::FindMChiP(
            length, informative, &maximum_position, &comparison, &maximum,
            chi_values.data());
        if (maximum_position < 0 || comparison < 0) return;
        const double initial_probability =
            MathFuncs::MyMathFuncs::ChiPVal2P(maximum);
        if (initial_probability *
                (static_cast<double>(informative) / half_window) * 3.0 >
                settings.lowest_probability ||
            initial_probability == 1.0) return;
        if (initial_probability >= settings.lowest_probability) return;

        const int original_maximum = maximum_position;
        const double initial_maximum = maximum;
        if (maximum_position == 0) maximum_position = 1;
        int growing_window = half_window;
        MathFuncs::MyMathFuncs::MakeTWinP(
            (initial_scan || inner_scan) ? 0 : 1, half_window,
            &growing_window, informative);
        int maximum_failures = half_window * 2;
        maximum_failures = std::min(
            maximum_failures, (informative - growing_window * 2) / 2);
        if (maximum_failures == 0) maximum_failures = 1;
        int left_count = 0;
        int right_count = 0;
        MathFuncs::MyMathFuncs::GetACP(
            informative, length, comparison, maximum_position,
            growing_window, &left_count, &right_count, scores.data());
        const int peak_position = maximum_position;
        if (maximum_position < 1) {
            maximum_position = settings.circular
                ? informative + maximum_position : 1;
        } else if (maximum_position > informative) {
            maximum_position = settings.circular
                ? maximum_position - informative : informative - 1;
        }
        int top_left = growing_window >= half_window ? left_count : 0;
        int top_right = growing_window >= half_window ? right_count : 0;
        ++growing_window;
        double probability = initial_probability *
            (static_cast<double>(informative) / half_window);
        int best_window = half_window;
        int top_left_position = maximum_position - half_window + 1;
        int top_right_position = maximum_position + half_window;
        int left = maximum_position - growing_window + 1;
        if (left < 0) left = informative + left;
        int right = maximum_position + growing_window;
        if (right >= informative * 2) right -= informative * 2;
        if (right >= informative) right -= informative;
        if (inner_scan) {
            MathFuncs::MyMathFuncs::GrowMChiWin2P2(
                max_window, left, right, informative, half_window,
                growing_window, comparison, length, left_count, right_count,
                maximum_failures, &probability, &best_window, &maximum,
                &top_left, &top_right, &top_left_position,
                &top_right_position, scores.data(),
                missing_map.data(),
                const_cast<float*>(table.values.data()),
                const_cast<int*>(table.map.data()));
        } else {
            MathFuncs::MyMathFuncs::GrowMChiWinP2(
                max_window, left, right, informative, half_window,
                growing_window, comparison, length, left_count, right_count,
                maximum_failures, &probability, &best_window, &maximum,
                &top_left, &top_right, &top_left_position, &top_right_position,
                scores.data(), const_cast<float*>(table.values.data()),
                const_cast<int*>(table.map.data()));
        }
        probability = (maximum < 20000.0)
            ? MathFuncs::MyMathFuncs::ChiPVal2P(maximum) *
                (static_cast<double>(informative) /
                 std::min(best_window, half_window))
            : 1.0e-200;
        probability *= 3.0;
        if (settings.mc_flag == 0) probability *= settings.mc_correction;

        if (trace) {
            trace->push_back({
                repetition, original_maximum, comparison, growing_window,
                left_count, right_count, best_window, top_left, top_right,
                top_left_position,
                top_right_position, initial_maximum, maximum, probability,
                difference_position_hash, position_difference_hash,
                score_hash, score_plane_hash,
                chi_table_hash, chi_map_hash,
                probability < settings.lowest_probability});
        }

        left = maximum_position - best_window;
        right = maximum_position + best_window - 1;
        if (right >= informative * 2) right -= informative * 2;
        if (right > informative) right -= informative;
        if (left < 1) left += informative;

        if (probability < settings.lowest_probability) {
            if (left - best_window + informative < 0)
                left = -informative + best_window;
            double high_left = 0.0;
            double high_right = 0.0;
            legacy_find_side(top_left, top_right, length, left, right,
                             best_window, informative, comparison, scores,
                             high_left, high_right);
            int beginning = 0;
            int ending = 0;
            bool low_high = false;
            if (high_left >= high_right) {
                low_high = top_right <= top_left;
                left = legacy_opt_left(
                    left, high_left, top_left, maximum_position, comparison,
                    best_window, informative, length, scores, missing_map,
                    trace ? &trace->back() : nullptr);
                if (left > informative) left = 1;
                if (left < 1) left += informative;
                if (maximum_position < 1)
                    right = informative + maximum_position;
                else if (maximum_position > informative)
                    right = maximum_position - informative;
                else if (maximum_position == 1) {
                    right = missing_map[1] == 0 && banned_windows[1] == 0
                        ? maximum_position : informative;
                }
                else
                    right = maximum_position;
                if (right > 0) --right;
                while (missing_map[right] != 0) {
                    --right;
                    if (right < 1) right = informative;
                }
                ++right;
                if (right > informative) right = informative;
                beginning = difference_position[
                    left + 1 > informative ? left : left + 1];
                ending = difference_position[right];
            } else {
                low_high = top_right > top_left;
                right = legacy_opt_right(
                    right, high_right, top_right, maximum_position, comparison,
                    best_window, informative, length, scores, missing_map);
                if (right < 1) right = informative;
                if (right >= informative * 2) right -= informative * 2;
                if (right > informative) right -= informative;
                if (maximum_position < 0)
                    left = informative + maximum_position;
                else if (maximum_position > informative)
                    left = maximum_position - informative;
                else
                    left = maximum_position;
                if (left < informative) ++left;
                while (missing_map[left] != 0) {
                    ++left;
                    if (left > informative) left = 1;
                }
                --left;
                if (left < 1) left = 1;
                beginning = difference_position[
                    left + 1 > informative ? left : left + 1];
                ending = difference_position[right];
            }
            if (initial_scan) {
                if (trace) {
                    trace->back().raw_beginning = beginning;
                    trace->back().raw_ending = ending;
                }
                centre_initial_breakpoints(
                    beginning, ending, position_difference,
                    difference_position, length, informative,
                    settings.circular != 0);
                if (trace) {
                    trace->back().first_centered_beginning = beginning;
                    trace->back().first_centered_ending = ending;
                }
                legacy_polish_maxchi_breakpoints(
                    comparison, low_high, sequences, length, informative,
                    settings.circular != 0, position_difference,
                    difference_position, missing_map,
                    scan_state.sequence_data, beginning, ending,
                    trace ? &trace->back() : nullptr);
                if (trace) {
                    trace->back().polished_beginning = beginning;
                    trace->back().polished_ending = ending;
                }
                centre_initial_breakpoints(
                    beginning, ending, position_difference,
                    difference_position, length, informative,
                    settings.circular != 0);
            }
            if (trace) {
                trace->back().centered_beginning = beginning;
                trace->back().centered_ending = ending;
            }
            maximum_position = peak_position;
            if (left < right) {
                if (!(maximum_position >= left && maximum_position <= right)) {
                    if (std::abs(maximum_position - left) >
                        std::abs(maximum_position - right)) right = maximum_position;
                    else left = maximum_position;
                }
            } else if (!(maximum_position > left || maximum_position < right)) {
                if (std::abs(maximum_position - left) >
                    std::abs(maximum_position - right)) right = maximum_position;
                else left = maximum_position;
            }
            if (trace) {
                trace->back().destroy_left = left;
                trace->back().destroy_right = right;
            }
            legacy_destroy_peaks(comparison, informative, length, left, right,
                                 smooth, chi_values);
            chi_values[original_maximum + comparison * sequence_stride] = 0;
            const auto active = choose_active(
                sequences, 3, store_lpv, store_lpv_ub, allocator);
            // MCXoverF stores the first event only after CentreBP.  In the
            // FinalTrim (FindallFlag=1) path it then runs FindDaughter on the
            // stored event and centres the polished endpoints once more,
            // before copying that final event into the reverse companion.
            const int slot = allocator.allocate(active[0], 3, probability);
            if (slot > 0) {
                auto& event = allocator.event(active[0], slot);
                event.daughter = static_cast<std::int16_t>(active[0]);
                event.minor_parent = static_cast<std::int16_t>(active[1]);
                event.major_parent = static_cast<std::int16_t>(active[2]);
                event.beginning = beginning;
                event.ending = ending;
                event.program_flag = 3;
                event.probability = probability;
                event.length_holder = best_window * 2;
                if (!initial_scan) {
                    centre_initial_breakpoints(
                        event.beginning, event.ending, position_difference,
                        difference_position, length, informative,
                        settings.circular != 0);
                    legacy_polish_maxchi_breakpoints(
                        comparison, low_high, sequences, length, informative,
                        settings.circular != 0, position_difference,
                        difference_position, missing_map,
                        scan_state.sequence_data, event.beginning,
                        event.ending, nullptr);
                    centre_initial_breakpoints(
                        event.beginning, event.ending, position_difference,
                        difference_position, length, informative,
                        settings.circular != 0);
                    if (!inner_scan) {
                        // The companion allocation can reallocate this row;
                        // retain the polished primary record before allocating.
                        auto reverse = event;
                        const int reverse_slot = allocator.allocate(
                            active[0], 3, probability);
                        // MCXoverF does not apply TSXOver's rectangular-slot
                        // guard here: it copies whenever UpdateXOList3 returns a
                        // positive SIP, growing the second dimension as needed.
                        if (reverse_slot > 0) {
                            std::swap(reverse.beginning, reverse.ending);
                            allocator.event(active[0], reverse_slot) = reverse;
                        }
                    }
                }
            }
            // MCXoverF resets WasteOfTime after every accepted grown peak.
            wasted_peaks = 0;
        } else {
            // Three consecutive grown peaks above LowestProb terminate the
            // source loop before the third peak is deleted.
            if (++wasted_peaks == 3) return;
            // MCXoverF saves maxz before normalizing a position-zero peak,
            // but saves pMaxX afterwards. DestroyPeakP receives pMaxX while
            // the following scalar clear still addresses maxz.
            int point = peak_position;
            double scratch[2]{};
            MathFuncs::MyMathFuncs::DestroyPeakP(
                comparison, length, point, point, informative, scratch,
                smooth.data(), chi_values.data());
            chi_values[original_maximum + comparison * sequence_stride] = 0;
        }
    }
}

void run_rdp_chimaera_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, const int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, const int event_beginning,
    const int event_ending, const bool initial_scan,
    int* const shared_xdiffpos0, const bool inner_scan,
    const std::vector<unsigned char>* missing_data) {
    constexpr int window_size = 60;
    const int length = scan_state.sequence_length;
    int half_window = vb_clng(window_size / 2.0);
    int critical = critical_difference(window_size, settings.lowest_probability);
    if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
        std::cerr << "temp-chimaera-settings half=" << half_window
                  << " critical=" << critical << " initial=" << initial_scan
                  << " event=" << event_beginning << ':' << event_ending << '\n';
    }
    std::vector<int> difference_position(length + 2, 0);
    std::vector<int> position_difference(length + 2, 0);
    const int informative = MathFuncs::MyMathFuncs::FindSubSeqDP(
        length + 1, static_cast<short>(scan_state.next_no),
        static_cast<short>(sequences[0]), static_cast<short>(sequences[1]),
        static_cast<short>(sequences[2]),
        const_cast<short*>(scan_state.sequence_data.data()),
        difference_position.data(), position_difference.data());
    if (informative < critical * 2 || informative < 7) return;
    if (initial_scan) {
        if (half_window * 2 > informative)
            half_window = vb_clng(informative * 0.75 / 2.0 + 0.00001) - 1;
        if (half_window <= critical)
            half_window = vb_clng(informative / 2.0 + 0.00001) - 1;
    } else {
        half_window = event_half_window(
            event_beginning, event_ending, informative, critical, half_window,
            position_difference);
    }
    if (half_window < 0) return;

    const int win_upper = length + half_window * 2;
    const int sequence_stride = length + 1;
    std::vector<unsigned char> scores(sequence_stride, 0);
    std::vector<int> window_scores(win_upper + 1, 0);
    std::vector<double> chi_values(sequence_stride, 0.0);
    std::vector<double> smooth(sequence_stride, 0.0);
    std::vector<int> banned_windows(length + half_window * 2 + 2, 0);
    std::vector<unsigned char> missing_map(length + 3, 0);
    if (inner_scan) {
        // CXoverA's FindallFlag=0/S event-number>0 path rebuilds the
        // informative-coordinate missing-data map before CalcChiVals5 and
        // GrowMChiWin2P.  Pass the same XPosDiff/XDiffPos orientation as the
        // DNA5 MakeBanWinP routine; raw alignment coordinates are not valid
        // indices for MDMap.
        if (missing_data != nullptr &&
            missing_data->size() >=
                static_cast<std::size_t>(length + 1) *
                    static_cast<std::size_t>(scan_state.next_no + 1)) {
            MathFuncs::MyMathFuncs::MakeBanWinP(
                length + half_window * 2, sequences[0], sequences[1],
                sequences[2], half_window, length, informative,
                banned_windows.data(), missing_map.data(),
                const_cast<unsigned char*>(missing_data->data()),
                position_difference.data(), difference_position.data());
        } else {
            std::fill(banned_windows.begin(), banned_windows.end(), 0);
            std::fill(missing_map.begin(), missing_map.end(), 0);
        }
    }
    if (inner_scan && settings.circular == 0) {
        missing_map[1] = 1;
        missing_map[informative] = 1;
        for (int position = informative - half_window + 2;
             position <= informative; ++position) {
            if (position >= 0 &&
                position < static_cast<int>(banned_windows.size())) {
                banned_windows[position] = 1;
            }
        }
    }
    const int score_goon = initial_scan
        ? MathFuncs::MyMathFuncs::WinScoreCalc4P2(
              critical, half_window, informative, length + 1, sequences[0],
              sequences[1], sequences[2], scores.data(),
              difference_position.data(),
              const_cast<short*>(scan_state.sequence_data.data()),
              window_scores.data())
        : legacy_win_score_calc4(
              critical, half_window, informative, length + 1, sequences[0],
              sequences[1], sequences[2], scores, difference_position,
              scan_state.sequence_data, window_scores);
    if (score_goon == 0) return;
    double maximum = initial_scan
        ? MathFuncs::MyMathFuncs::CalcChiVals3P(
              critical, half_window, informative, length,
              window_scores.data(), chi_values.data())
        : inner_scan
        ? MathFuncs::MyMathFuncs::CalcChiVals5P(
              critical, half_window, informative, length,
              window_scores.data(), chi_values.data(), banned_windows.data())
        : legacy_calc_chi_vals3(
              critical, half_window, informative, length, window_scores,
              chi_values);
    if (MathFuncs::MyMathFuncs::ChiPVal2P(maximum) *
            (static_cast<double>(informative) / half_window) * 3.0 >
        settings.lowest_probability) return;
    MathFuncs::MyMathFuncs::SmoothChiVals3P(
        informative, length, chi_values.data(), smooth.data());

    // Module30.MCXoverG updates the module-level XDiffPos(0) here, after
    // the initial chi-value pass but before it starts consuming peaks.  The
    // FindSubSeqDP call itself leaves element zero untouched.  This is a
    // deliberately unconditional per-call assignment: each successful
    // Chimaera prepass owns a fresh source map, so a later call may replace
    // the value left by an earlier one.
    if (shared_xdiffpos0 != nullptr && informative > 0) {
        *shared_xdiffpos0 = difference_position[1];
    }
    int wasted_peaks = 0;
    double prior_initial_probability = 0.0;
    for (int repetition = 1; repetition <= 100; ++repetition) {
        int maximum_position = -1;
        short comparison = -1;
        if (initial_scan) {
            MathFuncs::MyMathFuncs::FindMChi3P(
                length, informative, &maximum_position, &comparison, &maximum,
                chi_values.data());
        } else {
            legacy_find_mchi3(
                length, informative, maximum_position, comparison, maximum,
                chi_values);
        }
        if (maximum_position < 0 || comparison < 0) return;
        const double initial_probability =
            MathFuncs::MyMathFuncs::ChiPVal2P(maximum);
        if (initial_probability == prior_initial_probability) return;
        prior_initial_probability = initial_probability;
        if (initial_probability *
                (static_cast<double>(informative) / half_window) * 3.0 >
                settings.lowest_probability ||
            initial_probability == 1.0) return;
        if (initial_probability >= settings.lowest_probability) return;

        const int original_maximum = maximum_position;
        if (maximum_position == 0) maximum_position = 1;
        int growing_window = half_window;
        MathFuncs::MyMathFuncs::MakeTWinP(
            initial_scan ? 0 : 2, half_window, &growing_window, informative);
        int maximum_failures = std::min(
            half_window * 2, (informative - growing_window * 2) / 2);
        if (maximum_failures == 0) maximum_failures = 1;
        int left_count = 0;
        int right_count = 0;
        MathFuncs::MyMathFuncs::GetACP(
            informative, length, 0, maximum_position, growing_window,
            &left_count, &right_count, scores.data());
        const int peak_position = maximum_position;
        maximum_position = circular_position(maximum_position, informative);
        int top_left = growing_window >= half_window ? left_count : 0;
        int top_right = growing_window >= half_window ? right_count : 0;
        ++growing_window;
        double probability = initial_probability *
            (static_cast<double>(informative) / half_window);
        int best_window = half_window;
        int top_left_position = maximum_position - half_window + 1;
        int top_right_position = maximum_position + half_window;
        int left = maximum_position - growing_window + 1;
        if (left < 0) left += informative;
        int right = maximum_position + growing_window;
        if (right >= informative * 2) right -= informative * 2;
        if (right >= informative) right -= informative;
        if (inner_scan) {
            MathFuncs::MyMathFuncs::GrowMChiWin2P(
                left, right, informative, half_window, growing_window, 0,
                length, left_count, right_count, maximum_failures,
                &probability, &best_window, &maximum, &top_left, &top_right,
                &top_left_position, &top_right_position, scores.data(),
                missing_map.data());
        } else {
            MathFuncs::MyMathFuncs::GrowMChiWinP(
                left, right, informative, half_window, growing_window, 0,
                length, left_count, right_count, maximum_failures,
                &probability, &best_window, &maximum, &top_left, &top_right,
                &top_left_position, &top_right_position, scores.data());
        }
        if (maximum < 0) maximum = chi_values[original_maximum];
        probability = (maximum < 20000.0)
            ? MathFuncs::MyMathFuncs::ChiPVal2P(maximum) *
                (static_cast<double>(informative) /
                 std::min(best_window, half_window))
            : 1.0e-200;
        probability *= 3.0;
        if (settings.mc_flag == 0) probability *= settings.mc_correction;

        left = maximum_position - best_window;
        right = maximum_position + best_window - 1;
        if (right >= informative * 2) right -= informative * 2;
        if (right > informative) right -= informative;
        if (left < 1) left += informative;
        if (probability < settings.lowest_probability) {
            if (left - best_window + informative < 0)
                left = -informative + best_window;
            double high_left = 0.0;
            double high_right = 0.0;
            legacy_find_side(top_left, top_right, length, left, right,
                             best_window, informative, 0, scores,
                             high_left, high_right);
            if (high_left == 0.0 && high_right == 0.0) return;
            int beginning = 0;
            int ending = 0;
            int polish_comparison = 0;
            if (high_left >= high_right) {
                polish_comparison = top_right > top_left ? 1 : 0;
                left = legacy_opt_left(left, high_left, top_left,
                    maximum_position, 0, best_window, informative,
                    length, scores, missing_map, nullptr);
                if (left > informative) left = 1;
                if (left < 1) left += informative;
                if (maximum_position < 1)
                    right = informative + maximum_position;
                else if (maximum_position > informative)
                    right = maximum_position - informative;
                else if (maximum_position == 1)
                    right = missing_map[1] == 0 ? maximum_position : informative;
                else
                    right = maximum_position;
                if (right > 0) --right;
                while (missing_map[right] != 0) {
                    --right;
                    if (right < 1) right = informative;
                }
                ++right;
                if (right > informative) right = informative;
                beginning = difference_position[
                    left + 1 > informative ? left : left + 1];
                ending = difference_position[right];
            } else {
                polish_comparison = top_right > top_left ? 0 : 1;
                right = legacy_opt_right(right, high_right, top_right,
                    maximum_position, 0, best_window, informative,
                    length, scores, missing_map);
                if (right < 1) right = informative;
                if (right >= informative * 2) right -= informative * 2;
                if (right > informative) right -= informative;
                if (maximum_position < 0)
                    left = informative + maximum_position;
                else if (maximum_position > informative)
                    left = maximum_position - informative;
                else
                    left = maximum_position;
                if (left < informative) ++left;
                while (missing_map[left] != 0) {
                    ++left;
                    if (left > informative) left = 1;
                }
                --left;
                if (left < 1) left = 1;
                beginning = difference_position[
                    left + 1 > informative ? left : left + 1];
                ending = difference_position[right];
            }
            // CXoverA's FindDaughter(PrgF=4) uses the same MCMaxY choice as
            // the endpoint search.  Unlike MCXoverF, CXoverA forces
            // LoHiFlag=1 in both branches.
            const int daughter_comparison = polish_comparison;
            std::uint8_t outside_flag = 0;
            // Source CXoverA calls CentreBP before storing the event on both
            // the initial scan and the FinalTrim recheck.
            centre_initial_breakpoints(
                beginning, ending, position_difference,
                difference_position, length, informative,
                settings.circular != 0);
            if (initial_scan) {
                const int left_index = position_difference[beginning];
                const int right_index = position_difference[ending];
                std::array<int, 2> inside{};
                std::array<int, 2> outside{};
                const auto count_score = [&](std::array<int, 2>& counts,
                                             const int index) {
                    ++counts[scores[index] == 0 ? 0 : 1];
                };
                if (left_index < right_index) {
                    for (int x = 0; x < left_index; ++x)
                        count_score(outside, x);
                    for (int x = left_index; x <= right_index; ++x)
                        count_score(inside, x);
                    for (int x = right_index + 1; x <= informative; ++x)
                        count_score(outside, x);
                } else {
                    for (int x = 0; x <= right_index; ++x)
                        count_score(inside, x);
                    for (int x = right_index + 1; x < left_index; ++x)
                        count_score(outside, x);
                    for (int x = left_index; x <= informative; ++x)
                        count_score(inside, x);
                }
                outside_flag =
                    (inside[0] > inside[1] && outside[0] > outside[1]) ||
                    (inside[1] > inside[0] && outside[1] > outside[0]);
                legacy_polish_maxchi_breakpoints(
                    polish_comparison, true, sequences, length, informative,
                    settings.circular != 0, position_difference,
                    difference_position, missing_map,
                    scan_state.sequence_data, beginning, ending, nullptr);
                centre_initial_breakpoints(
                    beginning, ending, position_difference,
                    difference_position, length, informative,
                    settings.circular != 0);
            }
            if (!initial_scan) {
                // FinalTrim calls FindDaughter followed by a second CentreBP
                // before creating CXoverA's FindallFlag=2 companion.
                legacy_polish_maxchi_breakpoints(
                    daughter_comparison, true, sequences, length, informative,
                    settings.circular != 0, position_difference,
                    difference_position, missing_map,
                    scan_state.sequence_data, beginning, ending, nullptr);
                centre_initial_breakpoints(
                    beginning, ending, position_difference,
                    difference_position, length, informative,
                    settings.circular != 0);
            }
            maximum_position = peak_position;
            if (left < right) {
                if (!(maximum_position >= left && maximum_position <= right)) {
                    if (std::abs(maximum_position - left) >
                        std::abs(maximum_position - right)) right = maximum_position;
                    else left = maximum_position;
                }
            } else if (!(maximum_position > left || maximum_position < right)) {
                if (std::abs(maximum_position - left) >
                    std::abs(maximum_position - right)) right = maximum_position;
                else left = maximum_position;
            }
            legacy_destroy_peaks(0, informative, length, left, right,
                                 smooth, chi_values);
            chi_values[original_maximum] = 0;
            const auto active = choose_active(
                sequences, 4, store_lpv, store_lpv_ub, allocator);
            // FinalTrim sets DontWorryAboutSplitsFlag=1, so CXoverA's
            // FindallFlag>1 missing-data block is skipped, but the reverse
            // slot is still allocated and copied for FindallFlag=2.
            const int slot = allocator.allocate(active[0], 4, probability);
            if (slot > 0) {
                auto& event = allocator.event(active[0], slot);
                event.daughter = static_cast<std::int16_t>(active[0]);
                event.minor_parent = static_cast<std::int16_t>(active[1]);
                event.major_parent = static_cast<std::int16_t>(active[2]);
                event.beginning = beginning;
                event.ending = ending;
                event.program_flag = 4;
                event.probability = probability;
                event.length_holder = best_window * 2;
                event.outside_flag = outside_flag;
                event.distance_holder = sequences[0];
                // The source's CXoverA(2) copy is bounded by the rectangular
                // XOverList upper bound; the native transition fixtures show
                // no additional method-4 slots for this pass. Keep the
                // primary event here until that strict-bound behavior is
                // reproduced independently.
            }
            wasted_peaks = 0;
        } else {
            if (++wasted_peaks == 3) return;
            int point = peak_position;
            double scratch[2]{};
            MathFuncs::MyMathFuncs::DestroyPeakP(
                0, length, point, point, informative, scratch,
                smooth.data(), chi_values.data());
            chi_values[original_maximum] = 0;
        }
    }
}

void run_rdp_three_seq_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, const int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    const std::vector<float>& probability_table, const int table_bound,
    RdpLegacyEventAllocator& allocator) {
    auto result = evaluate_rdp_three_seq(
        scan_state, sequences, settings.circular != 0, settings.mc_flag,
        settings.mc_correction, settings.lowest_probability,
        probability_table, table_bound, true);
    const int informative = result.informative_last + 1;
    std::vector<int> difference_position(
        static_cast<std::size_t>(informative + 1), 0);
    for (int index = 1; index <= informative; ++index) {
        difference_position[index] = result.difference_position[index - 1];
    }
    auto centre = [&](int& beginning, int& ending) {
        if (informative < 1) return;
        auto position_difference = result.position_difference;
        centre_initial_breakpoints(
            beginning, ending, position_difference, difference_position,
            scan_state.sequence_length, informative, settings.circular != 0);
    };
    for (const auto& side : result.sides) {
        if (!side.significant) continue;
        const auto active = choose_active(
            sequences, 8, store_lpv, store_lpv_ub, allocator);
        const int slot = allocator.allocate(active[0], 8, side.probability);
        if (slot < 1) continue;
        auto& event = allocator.event(active[0], slot);
        event.daughter = static_cast<std::int16_t>(active[0]);
        event.minor_parent = static_cast<std::int16_t>(active[1]);
        event.major_parent = static_cast<std::int16_t>(active[2]);
        event.beginning = side.beginning;
        event.ending = side.ending;
        // The event coordinates returned by the current TS prefix are
        // already CentreBP-equivalent; keep the historical initial-path
        // representation here while the FinalTrim side effects are matched
        // against the native fixture separately.
        event.program_flag = 8;
        event.probability = side.probability;
        event.distance_holder = sequences[2];

        // As with the other legacy methods, the second allocation may grow
        // the vector backing the first slot, so snapshot it first.
        auto reverse = event;
        const int companion = allocator.allocate(
            active[0], 8, side.probability);
        if (companion > -1 &&
            allocator.has_strictly_later_slot(companion)) {
            std::swap(reverse.major_parent, reverse.minor_parent);
            std::swap(reverse.beginning, reverse.ending);
            allocator.event(active[0], companion) = reverse;
        }
    }
}

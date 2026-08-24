#include "analysis.hpp"

#include "MathFuncsDll.h"
#include "distance_state.hpp"
#include "tree_state.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct RdpFactorialTables {
    static constexpr int three_way_upper_bound = 97;
    std::vector<double> factorial;
    std::vector<double> three_way;
};

RdpFactorialTables make_factorial_tables() {
    RdpFactorialTables tables;
    tables.factorial.assign(172, 0.0);
    tables.factorial[0] = 1.0;
    for (int value = 1; value <= 170; ++value) {
        tables.factorial[value] = tables.factorial[value - 1] * value;
    }

    constexpr int stride = RdpFactorialTables::three_way_upper_bound + 1;
    std::vector<double> factorial_three(stride, 0.0);
    factorial_three[0] = 1.0;
    for (int value = 1; value < stride; ++value) {
        factorial_three[value] = factorial_three[value - 1] * value;
    }
    tables.three_way.resize(
        static_cast<std::size_t>(stride) * stride * stride);
    for (int first = 0; first < stride; ++first) {
        for (int second = 0; second < stride; ++second) {
            for (int third = 0; third < stride; ++third) {
                const double denominator =
                    factorial_three[second] * factorial_three[third];
                tables.three_way[
                    static_cast<std::size_t>(first) +
                    static_cast<std::size_t>(second) * stride +
                    static_cast<std::size_t>(third) * stride * stride] =
                        factorial_three[first] / denominator;
            }
        }
    }
    return tables;
}

std::vector<double> make_probability_estimate(
    std::vector<double>& factorial) {
    constexpr int upper_bound = 171;
    constexpr int stride = upper_bound + 1;
    constexpr int categories = 51;
    std::vector<double> estimate(
        static_cast<std::size_t>(stride) * stride * categories, 0.0);
    for (int length = 1; length <= upper_bound; ++length) {
        for (int common = 0; common <= upper_bound; ++common) {
            for (int category = 0; category < categories; ++category) {
                estimate[
                    static_cast<std::size_t>(length) +
                    static_cast<std::size_t>(common) * stride +
                    static_cast<std::size_t>(category) * stride * stride] =
                        100.0;
            }
        }
    }
    for (int length = 1; length <= 170; ++length) {
        for (int common = 3; common <= length; ++common) {
            for (int category = 2; category < categories; ++category) {
                estimate[
                    static_cast<std::size_t>(length) +
                    static_cast<std::size_t>(common) * stride +
                    static_cast<std::size_t>(category) * stride * stride] =
                        MathFuncs::MyMathFuncs::ProbCalcP(
                            factorial.data(), length, common,
                            static_cast<double>(category) / 50.0,
                            length * 2);
            }
        }
    }
    return estimate;
}

std::vector<unsigned char> make_fss_rdp(const int upper_bound) {
    const int stride = upper_bound + 1;
    std::vector<unsigned char> table(
        static_cast<std::size_t>(4) * stride * stride * stride, 0);
    const auto cell = [stride](const int plane, const int first,
                               const int second, const int third) {
        return static_cast<std::size_t>(plane) +
            static_cast<std::size_t>(first) * 4 +
            static_cast<std::size_t>(second) * 4 * stride +
            static_cast<std::size_t>(third) * 4 * stride * stride;
    };
    const auto encoded_nucleotide = [](const int encoded, const int site) {
        int value = encoded;
        for (int skipped = 2; skipped > site; --skipped) value /= 5;
        return value % 5;
    };
    for (int first = 0; first < 125; ++first) {
        for (int second = 0; second < 125; ++second) {
            for (int third = 0; third < 125; ++third) {
                int sum = 0;
                for (int site = 0; site < 3; ++site) {
                    const int a = encoded_nucleotide(first, site);
                    const int b = encoded_nucleotide(second, site);
                    const int c = encoded_nucleotide(third, site);
                    unsigned char action = 0;
                    if (a != 0 && b != 0 && c != 0 &&
                        (a != c || a != b)) {
                        if (a == b) action = 1;
                        else if (a == c) action = 2;
                        else if (b == c) action = 3;
                    }
                    table[cell(site, first, second, third)] = action;
                    sum += action;
                }
                table[cell(3, first, second, third)] =
                    static_cast<unsigned char>(sum);
            }
        }
    }
    return table;
}

RdpFactorialTables& shared_factorial_tables() {
    static thread_local RdpFactorialTables tables = make_factorial_tables();
    return tables;
}

std::vector<double>& shared_probability_estimate() {
    static thread_local std::vector<double> estimate =
        make_probability_estimate(shared_factorial_tables().factorial);
    return estimate;
}

std::vector<unsigned char>& shared_fss_rdp() {
    static thread_local std::vector<unsigned char> table = make_fss_rdp(125);
    return table;
}

Dna5ScanPreprocessApi preprocessing_api() {
    return {
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
}

Dna5XoverApi xover_api() {
    return {
        &MathFuncs::MyMathFuncs::FindSubSeqPB3,
        &MathFuncs::MyMathFuncs::XOHomologyP,
        &MathFuncs::MyMathFuncs::FindNextP,
        &MathFuncs::MyMathFuncs::FindFirstCOP,
        &MathFuncs::MyMathFuncs::DefineEventP2,
        &MathFuncs::MyMathFuncs::ProbCalcP2,
        &MathFuncs::MyMathFuncs::ProbCalcP,
        &MathFuncs::MyMathFuncs::FindSubSeqPB4,
        &MathFuncs::MyMathFuncs::CleanXOSNW,
    };
}

RdpInitialAnalysisResult run_rdp_initial_analysis(
    RdpScanState scan_state, const RdpInitialAnalysisOptions& options) {
    if (options.p_value_cutoff <= 0.0 || options.p_value_cutoff > 1.0) {
        throw std::runtime_error("RDP p-value cutoff must be in (0, 1]");
    }
    if (options.window_sites < 2 || options.window_sites > 32767) {
        throw std::runtime_error("RDP window must be between 2 and 32767 sites");
    }

    auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    auto tree_state = build_rdp_upgma_tree_state(
        scan_state.next_no, distance_state);
    const int triplet_count = scan_state.analysis_list_last + 1;
    const int store_lpv_upper_bound = 8;
    const int fss_upper_bound = 125;
    std::vector<double> store_lpv(
        static_cast<std::size_t>(store_lpv_upper_bound + 1) *
            (scan_state.next_no + 1),
        1.0);
    std::vector<unsigned char> redo(triplet_count, 0);
    auto& fss_rdp = shared_fss_rdp();
    auto& factorials = shared_factorial_tables();
    auto& probability_estimate = shared_probability_estimate();
    const int correction_tests = triplet_count;
    const double uncorrected_threshold =
        options.p_value_cutoff / correction_tests;
    const int target = static_cast<int>(std::nearbyint(
        static_cast<double>(scan_state.sequence_length) / 10.0));
    const short full_window = static_cast<short>(options.window_sites);
    const int half_window = static_cast<int>(std::nearbyint(
        static_cast<double>(options.window_sites) / 2.0));
    MathFuncs::MyMathFuncs::AlistRDP4(
        store_lpv_upper_bound, store_lpv.data(),
        scan_state.analysis_list.data(), scan_state.analysis_list_last,
        0, scan_state.analysis_list_last, scan_state.next_no,
        uncorrected_threshold, redo.data(), options.circular ? 1 : 0,
        correction_tests, 0, options.p_value_cutoff, target,
        scan_state.sequence_length, 100, scan_state.next_no,
        distance_state.distance.data(), scan_state.next_no,
        tree_state.tree_distance.data(), fss_upper_bound,
        scan_state.compressed_sequence_ub,
        scan_state.compressed_sequence.data(), scan_state.sequence_data.data(),
        half_window, full_window, fss_rdp.data(), 0, 171, 171,
        probability_estimate.data(),
        RdpFactorialTables::three_way_upper_bound,
        factorials.three_way.data(), factorials.factorial.data());

    const RdpXoverSettings xover_settings{
        100,
        1,
        target,
        options.circular ? 1 : 0,
    };
    const RdpProbabilitySettings probability_settings{
        options.circular ? 1 : 0,
        correction_tests,
        0,
        0,
        171,
        171,
        RdpFactorialTables::three_way_upper_bound,
        options.p_value_cutoff,
    };
    auto events = scan_rdp_redo_triplets(
        scan_state, distance_state, tree_state, redo, fss_rdp, store_lpv,
        store_lpv_upper_bound, fss_upper_bound, half_window, full_window,
        xover_settings, probability_settings, probability_estimate,
        factorials.three_way, factorials.factorial, xover_api());
    return {std::move(scan_state), std::move(events)};
}

}  // namespace

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_text(
    const std::string& fasta_text, const RdpInitialAnalysisOptions& options) {
    return run_rdp_initial_analysis(
        build_rdp_scan_state_from_fasta_text(fasta_text, preprocessing_api()),
        options);
}

RdpInitialAnalysisResult run_rdp_initial_analysis_from_fasta_file(
    const std::string& fasta_path, const RdpInitialAnalysisOptions& options) {
    return run_rdp_initial_analysis(
        build_rdp_scan_state_from_fasta(fasta_path, preprocessing_api()),
        options);
}

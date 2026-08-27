#include "analysis.hpp"

#include "MathFuncsDll.h"
#include "burt_state.hpp"
#include "distance_state.hpp"
#include "legacy_method_state.hpp"
#include "legacy_optional/bootscan.hpp"
#include "legacy_optional/chimaera.hpp"
#include "legacy_optional/geneconv.hpp"
#include "legacy_optional/maxchi.hpp"
#include "legacy_optional/siscan.hpp"
#include "legacy_optional/threeseq.hpp"
#include "legacy_optional_bridge.hpp"
#include "mutation_state.hpp"
#include "rescan_schedule.hpp"
#include "selection_state.hpp"
#include "three_seq_state.hpp"
#include "tree_state.hpp"

#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

int count_events(const RdpRawEventState& events) {
    int count = 0;
    for (const auto& row : events.xover_list) {
        count += static_cast<int>(row.size());
    }
    return count;
}

RdpSlidingWindowProfile make_rdp_sliding_window_profile(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const RdpRawEvent& event,
    const int window_sites, std::vector<unsigned char>& fss_rdp,
    const Dna5XoverApi& api) {
    RdpSlidingWindowProfile profile;
    profile.window_sites = window_sites;
    // DrawPlots deliberately uses the odd, symmetric window width rather
    // than the number of valid observations in each window.
    const int source_half_window = std::max(1, window_sites / 2);
    profile.divisor = source_half_window * 2 + 1;
    profile.alignment_length = scan_state.sequence_length;
    if (event.profile_available == 0) return profile;
    std::array<int, 3> sequences{};
    for (int role = 0; role < 3; ++role) {
        sequences[role] = event.profile_sequences[role];
        if (sequences[role] < 0 || sequences[role] > scan_state.next_no) {
            return profile;
        }
    }

    // XOver uses the compressed PB3/PB4 path for the initial and rescan
    // passes, and the plain FindSubSeqP path in FinalTrim.  Preserve that
    // choice instead of silently reconstructing a different profile.
    const bool use_compress = event.profile_use_compress != 0;
    auto state = build_rdp_first_xover_state(
        scan_state, distance_state, tree_state, 0, 125, fss_rdp,
        source_half_window, static_cast<short>(window_sites), api,
        &sequences, use_compress ? 0 : 1, use_compress);
    if (state.informative_length < source_half_window * 2 ||
        state.homology_length <= 0) {
        return profile;
    }
    if (use_compress) {
        // PB3 writes the compact agreement rows but not XDiffPos.  The
        // source obtains those coordinates lazily through PB4 before
        // DrawPlots is called.
        state.position_map_result = api.find_subsequence_with_positions(
            state.agreement_counts.data(), 125, source_half_window,
            scan_state.compressed_sequence_ub, scan_state.sequence_length,
            scan_state.next_no, sequences[0], sequences[1], sequences[2],
            const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
            state.xover_sequence_ub, state.xover_sequence.data(),
            state.xdiffpos.data(), state.xposdiff.data(), fss_rdp.data());
    }

    const int stride = state.homology_ub + 1;
    profile.positions.reserve(static_cast<std::size_t>(state.homology_length + 1));
    for (auto& row : profile.counts) {
        row.reserve(static_cast<std::size_t>(state.homology_length + 1));
    }
    int lower = 1;
    int upper = 0;
    for (int position = 0; position <= state.homology_length; ++position) {
        int alignment_position = state.xdiffpos[static_cast<std::size_t>(position)];
        // Module30.DrawPlots repairs the sentinel XDiffPos(0) this way,
        // producing the small initial vertical segment in the legacy plot.
        if (position == 0 && alignment_position == 0) {
            alignment_position = state.xdiffpos[1];
        }
        if (alignment_position <= 0) alignment_position = 1;
        profile.positions.push_back(alignment_position);
        for (int pair = 0; pair < 3; ++pair) {
            const int value = state.homology[
                position + pair * stride];
            profile.counts[pair].push_back(value);
            if (position > 0) {
                lower = std::min(lower, value);
                upper = std::max(upper, value);
            }
        }
    }
    profile.minimum = static_cast<double>(lower) /
        static_cast<double>(profile.divisor);
    profile.maximum = static_cast<double>(upper) /
        static_cast<double>(profile.divisor);
    profile.exact = !profile.positions.empty();
    return profile;
}

void report_progress(const RdpInitialAnalysisOptions& options,
                     const int phase, const int round,
                     const int processed_triplets, const int total_triplets,
                     const int event_count) {
    if (options.progress_callback == nullptr) return;
    options.progress_callback(
        phase, round, processed_triplets, total_triplets, event_count,
        options.progress_user);
}

template <typename Function>
void run_initial_rdp_ranges(const int count, Function&& function) {
    if (count <= 0) return;
#if defined(NEXT_RDP_METHOD_MANUAL_THREADS) && !defined(NEXT_RDP_USE_REAL_OPENMP)
    // Emscripten has no libomp in the supported toolchain.  AlistRDP4's
    // per-row work is independent (each row writes only RL[y]), so invoke the
    // same vendored routine over deterministic contiguous ranges on pthreads.
    // Native OpenMP builds call it once and let the literal source pragma
    // select its team, avoiding nested parallel regions.
    const int requested = std::max(1, rdp_method_worker_threads());
    const int workers = std::min(requested, count);
    if (workers <= 1) {
        function(0, count - 1);
        return;
    }
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(workers));
    std::exception_ptr failure;
    std::mutex failure_mutex;
    for (int worker = 0; worker < workers; ++worker) {
        const int first = (count * worker) / workers;
        const int last = (count * (worker + 1)) / workers - 1;
        threads.emplace_back([&, first, last] {
            try {
                if (first <= last) function(first, last);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        });
    }
    for (auto& thread : threads) thread.join();
    if (failure) std::rethrow_exception(failure);
#else
    function(0, count - 1);
#endif
}

void trace_legacy_method_state(const char* label,
                               const RdpRawEventState& events) {
    if (std::getenv("RDP_TRACE_TEMP") == nullptr) return;
    std::array<int, 9> counts{};
    int total = 0;
    for (const auto& row : events.xover_list) {
        for (const auto& event : row) {
            ++total;
            if (event.program_flag < counts.size()) {
                ++counts[event.program_flag];
            }
        }
    }
    std::cerr << "temp-method " << label << " total=" << total;
    for (int program = 0; program < static_cast<int>(counts.size()); ++program) {
        if (counts[program] != 0) std::cerr << " p" << program << '=' << counts[program];
    }
    std::cerr << " rows=";
    for (std::size_t sequence = 0; sequence < events.xover_list.size(); ++sequence) {
        const auto& row = events.xover_list[sequence];
        if (row.empty()) continue;
        std::array<int, 9> row_counts{};
        for (const auto& event : row) {
            if (event.program_flag < row_counts.size()) {
                ++row_counts[event.program_flag];
            }
        }
        std::cerr << sequence << ':' << row.size() << '(';
        for (int program = 0; program < static_cast<int>(row_counts.size()); ++program) {
            if (row_counts[program] != 0) std::cerr << program << '=' << row_counts[program] << ',';
        }
        std::cerr << ')';
    }
    std::cerr << '\n';
    if (std::getenv("RDP_TRACE_ENDPOINTS") != nullptr &&
        std::string(label) == "collection") {
        for (std::size_t row = 0; row < events.xover_list.size(); ++row) {
            for (const auto& event : events.xover_list[row]) {
                if (event.program_flag == 0 &&
                    (event.daughter == 15 || event.daughter == 17 ||
                     event.daughter == 18)) {
                    std::cerr << "temp-endpoint row=" << row << " p0="
                              << event.daughter << ':' << event.major_parent
                              << ':' << event.minor_parent << ':'
                              << event.beginning << '-' << event.ending
                              << " prob=" << event.probability << '\n';
                }
            }
        }
    }
}

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

std::vector<float>& shared_three_seq_table() {
    static thread_local std::vector<float> table =
        make_rdp_three_seq_probability_table(44);
    return table;
}

thread_local const std::vector<float>* configured_three_seq_table = nullptr;
thread_local int configured_three_seq_table_bound = 45;

const std::vector<float>& active_three_seq_table() {
    return configured_three_seq_table != nullptr
        ? *configured_three_seq_table : shared_three_seq_table();
}

std::vector<float> load_three_seq_table_file(
    const std::string& fasta_path, int& table_bound) {
    const auto slash = fasta_path.find_last_of("/\\");
    const std::string path = (slash == std::string::npos
        ? std::string{} : fasta_path.substr(0, slash + 1)) + "3seqTable";
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return {};
    const auto bytes = static_cast<std::size_t>(input.tellg());
    if (bytes < sizeof(std::int32_t) + sizeof(float)) return {};
    input.seekg(0);
    std::int32_t maximum = 0;
    input.read(reinterpret_cast<char*>(&maximum), sizeof(maximum));
    if (maximum < 1) return {};
    std::vector<float> table((bytes - sizeof(maximum)) / sizeof(float));
    input.read(reinterpret_cast<char*>(table.data()),
               static_cast<std::streamsize>(table.size() * sizeof(float)));
    table_bound = maximum + 1;
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
        &MathFuncs::MyMathFuncs::FindSubSeqP,
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

void append_legacy_optional_events(
    RdpFullAnalysisResult& output,
    const RdpScanState& scan_state,
    const RdpInitialAnalysisOptions& options) {
    if (!options.enable_bootscan && !options.enable_siscan) return;
    const auto alignment = next_rdp_legacy_optional_bridge::make_alignment(
        scan_state);
    const std::size_t count = alignment.sequence_count();
    const std::vector<std::uint32_t> origins = [&] {
        std::vector<std::uint32_t> values(count);
        for (std::size_t index = 0; index < count; ++index) {
            values[index] = static_cast<std::uint32_t>(index);
        }
        return values;
    }();
    std::vector<std::uint8_t> disabled(count, 0);
    for (std::size_t index = 0;
         index < count && index < options.disabled_sequences.size(); ++index) {
        disabled[index] = options.disabled_sequences[index] != 0 ? 1 : 0;
    }
    std::vector<std::uint8_t> missing;
    next_rdp_legacy_optional::BootscanWorkspace bootscan_workspace;
    next_rdp_legacy_optional::SiscanWorkspace siscan_workspace;
    const std::uint64_t correction_tests = std::max<std::uint64_t>(
        1, static_cast<std::uint64_t>(output.triplet_count));
    const int list_last = scan_state.analysis_list_last;
    for (int index = 0; index <= list_last; ++index) {
        const std::array<std::uint32_t, 3> triplet{
            static_cast<std::uint32_t>(scan_state.analysis_list[index * 3]),
            static_cast<std::uint32_t>(scan_state.analysis_list[index * 3 + 1]),
            static_cast<std::uint32_t>(scan_state.analysis_list[index * 3 + 2])};
        bool skip = false;
        for (const auto sequence : triplet) {
            if (sequence >= options.disabled_sequences.size()) continue;
            if (options.disabled_sequences[sequence] != 0 ||
                (sequence < options.masked_sequences.size() &&
                 options.masked_sequences[sequence] != 0)) {
                skip = true;
                break;
            }
        }
        if (skip) continue;
        const auto similarities =
            next_rdp_legacy_optional_bridge::pair_similarity(alignment, triplet);
        missing = next_rdp_legacy_optional_bridge::triplet_missing_data(
            alignment, triplet);

        if (options.enable_bootscan) {
            next_rdp_legacy_optional::BootscanDiscoveryOptions discovery;
            discovery.circular = options.circular;
            discovery.bonferroni = options.correction_bonferroni;
            discovery.p_value_cutoff = options.p_value_cutoff;
            discovery.correction_tests = correction_tests;
            discovery.window_sites = std::max(5, options.bootscan_window_sites);
            discovery.step_sites = std::max(1, options.bootscan_step_sites);
            discovery.bootstrap_replicates = std::max(
                10, options.bootscan_bootstrap_replicates);
            discovery.support_cutoff = options.bootscan_support_cutoff;
            discovery.random_seed = options.bootscan_random_seed;
            std::vector<next_rdp_legacy_optional::BootscanDiscoveryCandidate> candidates;
            (void)next_rdp_legacy_optional::bootscan_discover(
                alignment, triplet, missing, similarities, discovery,
                bootscan_workspace, candidates);
            for (const auto& candidate : candidates) {
                RdpFinalEvent event;
                event.event_number = static_cast<int>(output.events.size()) + 1;
                event.program = 5;
                event.winning_role = 0;
                event.probability = candidate.corrected_p_value;
                event.beginning = static_cast<int>(candidate.beginning);
                event.ending = static_cast<int>(candidate.ending);
                event.representative_sequences = {
                    static_cast<int>(triplet[candidate.recombinant_local]),
                    static_cast<int>(triplet[candidate.major_parent_local]),
                    static_cast<int>(triplet[candidate.minor_parent_local])};
                event.profile_sequences = {
                    static_cast<int>(triplet[0]), static_cast<int>(triplet[1]),
                    static_cast<int>(triplet[2])};
                event.profile_sequences_available = true;
                for (int role = 0; role < 3; ++role) {
                    event.sequence_groups[role].push_back(
                        event.representative_sequences[role]);
                }
                event.bootscan_available = true;
                event.bootscan_discovery = candidate;
                output.events.push_back(std::move(event));
            }
        }

        if (options.enable_siscan) {
            next_rdp_legacy_optional::SiscanOptions discovery;
            discovery.circular = options.circular;
            discovery.bonferroni = options.correction_bonferroni;
            discovery.p_value_cutoff = options.p_value_cutoff;
            discovery.correction_tests = correction_tests;
            discovery.window_sites = std::max(5, options.siscan_window_sites);
            discovery.step_sites = std::max(1, options.siscan_step_sites);
            discovery.scan_permutations = std::max(
                10, options.siscan_scan_permutations);
            discovery.p_value_permutations = std::max<std::size_t>(
                discovery.scan_permutations,
                static_cast<std::size_t>(std::max(1, options.siscan_p_value_permutations)));
            discovery.random_seed = options.siscan_random_seed;
            std::vector<next_rdp_legacy_optional::SiscanDiscoveryCandidate> candidates;
            (void)next_rdp_legacy_optional::siscan_discover(
                alignment, triplet, missing, similarities, origins, disabled,
                discovery, siscan_workspace, candidates);
            for (const auto& candidate : candidates) {
                RdpFinalEvent event;
                event.event_number = static_cast<int>(output.events.size()) + 1;
                event.program = 6;
                event.winning_role = 0;
                event.probability = candidate.corrected_p_value;
                event.beginning = static_cast<int>(candidate.beginning);
                event.ending = static_cast<int>(candidate.ending);
                event.representative_sequences = {
                    static_cast<int>(triplet[candidate.recombinant_local]),
                    static_cast<int>(triplet[candidate.major_parent_local]),
                    static_cast<int>(triplet[candidate.minor_parent_local])};
                event.profile_sequences = {
                    static_cast<int>(triplet[0]), static_cast<int>(triplet[1]),
                    static_cast<int>(triplet[2])};
                event.profile_sequences_available = true;
                for (int role = 0; role < 3; ++role) {
                    event.sequence_groups[role].push_back(
                        event.representative_sequences[role]);
                }
                event.siscan_available = true;
                event.siscan_discovery = candidate;
                output.events.push_back(std::move(event));
            }
        }
    }
    output.raw_candidate_count = static_cast<int>(output.events.size());
}

void apply_legacy_optional_rechecks(
    RdpFullAnalysisResult& output,
    const RdpScanState& scan_state,
    const RdpInitialAnalysisOptions& options) {
    if (!options.enable_bootscan_secondary && !options.enable_siscan_secondary) {
        return;
    }
    const auto alignment = next_rdp_legacy_optional_bridge::make_alignment(scan_state);
    const auto origins = [&] {
        std::vector<std::uint32_t> values(alignment.sequence_count());
        for (std::size_t index = 0; index < values.size(); ++index) {
            values[index] = static_cast<std::uint32_t>(index);
        }
        return values;
    }();
    std::vector<std::uint8_t> disabled(alignment.sequence_count(), 0);
    for (std::size_t index = 0;
         index < disabled.size() && index < options.disabled_sequences.size();
         ++index) {
        disabled[index] = options.disabled_sequences[index] != 0 ? 1 : 0;
    }
    const std::uint64_t correction_tests = std::max<std::uint64_t>(
        1, static_cast<std::uint64_t>(output.triplet_count));
    next_rdp_legacy_optional::BootscanWorkspace bootscan_workspace;
    next_rdp_legacy_optional::SiscanWorkspace siscan_workspace;
    for (auto& event : output.events) {
        std::array<std::uint32_t, 3> triplet{};
        bool valid = true;
        for (int role = 0; role < 3; ++role) {
            const int sequence = event.representative_sequences[role];
            if (sequence < 0 || sequence >= static_cast<int>(alignment.sequence_count())) {
                valid = false;
                break;
            }
            triplet[role] = static_cast<std::uint32_t>(sequence);
        }
        if (!valid) continue;
        const auto missing = next_rdp_legacy_optional_bridge::triplet_missing_data(
            alignment, triplet);
        if (options.enable_bootscan_secondary) {
            next_rdp_legacy_optional::BootscanRecheckOptions recheck;
            recheck.circular = options.circular;
            recheck.bonferroni = options.correction_bonferroni;
            recheck.p_value_cutoff = options.p_value_cutoff;
            recheck.correction_tests = correction_tests;
            recheck.window_sites = std::max(5, options.bootscan_window_sites);
            recheck.step_sites = std::max(1, options.bootscan_step_sites);
            recheck.bootstrap_replicates = std::max(10, options.bootscan_bootstrap_replicates);
            recheck.support_cutoff = options.bootscan_support_cutoff;
            recheck.random_seed = options.bootscan_random_seed;
            event.bootscan_recheck = next_rdp_legacy_optional::bootscan_recheck(
                alignment, triplet, missing,
                static_cast<std::size_t>(std::max(1, event.beginning)),
                static_cast<std::size_t>(std::max(1, event.ending)), recheck,
                bootscan_workspace);
        }
        if (options.enable_siscan_secondary) {
            next_rdp_legacy_optional::SiscanOptions recheck;
            recheck.circular = options.circular;
            recheck.bonferroni = options.correction_bonferroni;
            recheck.p_value_cutoff = options.p_value_cutoff;
            recheck.correction_tests = correction_tests;
            recheck.window_sites = std::max(5, options.siscan_window_sites);
            recheck.step_sites = std::max(1, options.siscan_step_sites);
            recheck.scan_permutations = std::max(10, options.siscan_scan_permutations);
            recheck.p_value_permutations = std::max<std::size_t>(
                recheck.scan_permutations,
                static_cast<std::size_t>(std::max(1, options.siscan_p_value_permutations)));
            recheck.random_seed = options.siscan_random_seed;
            event.siscan_recheck = next_rdp_legacy_optional::siscan_recheck(
                alignment, triplet, missing, origins, disabled,
                static_cast<std::size_t>(std::max(1, event.beginning)),
                static_cast<std::size_t>(std::max(1, event.ending)), recheck,
                siscan_workspace);
        }
    }
}

// Module3.MakeAnalysisListQvR builds a deliberately constrained list rather
// than filtering MakeAListP2 after the fact.  References are the enabled
// sequences assigned to a positive group; group zero sequences are queries.
// The source walks references in input order, considers each unordered pair,
// rejects pairs from the same group, and then appends every query in input
// order.  Preserve that order because it is also the order used by the
// source correction/list screens.
void apply_query_reference_analysis_list(
    RdpScanState& scan_state, const RdpInitialAnalysisOptions& options) {
    if (!options.query_reference_mode) return;
    const int sequence_count = scan_state.next_no + 1;
    if (sequence_count < 3) {
        throw std::runtime_error(
            "Query-vs-reference analysis requires at least three sequences");
    }
    if (options.reference_groups.size() <
        static_cast<std::size_t>(sequence_count)) {
        throw std::runtime_error(
            "Query-vs-reference groups do not match the alignment");
    }
    const auto is_excluded = [&](const int sequence) {
        const bool masked = sequence < static_cast<int>(options.masked_sequences.size()) &&
            options.masked_sequences[static_cast<std::size_t>(sequence)] != 0;
        const bool disabled = sequence < static_cast<int>(options.disabled_sequences.size()) &&
            options.disabled_sequences[static_cast<std::size_t>(sequence)] != 0;
        return masked || disabled;
    };
    std::vector<int> references;
    std::vector<int> queries;
    references.reserve(static_cast<std::size_t>(sequence_count));
    queries.reserve(static_cast<std::size_t>(sequence_count));
    for (int sequence = 0; sequence < sequence_count; ++sequence) {
        if (is_excluded(sequence)) continue;
        const auto group = options.reference_groups[static_cast<std::size_t>(sequence)];
        if (group == 0) queries.push_back(sequence);
        else references.push_back(sequence);
    }
    std::vector<short> list;
    const std::size_t reserve_count = references.size() > 1
        ? (references.size() * (references.size() - 1) / 2) * queries.size()
        : 0;
    if (reserve_count >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / 3U) {
        throw std::runtime_error("Query-vs-reference analysis list is too large");
    }
    list.reserve(reserve_count * 3U);
    for (std::size_t first = 0; first < references.size(); ++first) {
        const int reference_a = references[first];
        const auto group_a = options.reference_groups[
            static_cast<std::size_t>(reference_a)];
        for (std::size_t second = first + 1; second < references.size(); ++second) {
            const int reference_b = references[second];
            const auto group_b = options.reference_groups[
                static_cast<std::size_t>(reference_b)];
            if (group_a == group_b) continue;
            for (const int query : queries) {
                list.push_back(static_cast<short>(reference_a));
                list.push_back(static_cast<short>(reference_b));
                list.push_back(static_cast<short>(query));
            }
        }
    }
    if (list.empty()) {
        throw std::runtime_error(
            "Query-vs-reference analysis needs at least one query and a pair of differently grouped references");
    }
    scan_state.analysis_list = std::move(list);
    scan_state.analysis_list_last = static_cast<int>(
        scan_state.analysis_list.size() / 3U) - 1;
}

RdpInitialAnalysisResult run_rdp_initial_analysis(
    RdpScanState scan_state, const RdpInitialAnalysisOptions& options) {
    if (options.p_value_cutoff <= 0.0 || options.p_value_cutoff > 1.0) {
        throw std::runtime_error("RDP p-value cutoff must be in (0, 1]");
    }
    if (options.window_sites < 2 || options.window_sites > 32767) {
        throw std::runtime_error("RDP window must be between 2 and 32767 sites");
    }

    apply_query_reference_analysis_list(scan_state, options);
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
    const auto run_alist_rdp4 = [&](const int first, const int last) {
        (void)MathFuncs::MyMathFuncs::AlistRDP4(
            store_lpv_upper_bound, store_lpv.data(),
            scan_state.analysis_list.data(), scan_state.analysis_list_last,
            first, last, scan_state.next_no,
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
    };
    // AlistRDP4 supplies the shared StoreLPV screening table used by the
    // optional source methods even when RDP itself is not selected.  The
    // actual RDP XOver walk is skipped for an optional-method-only scan.
    run_initial_rdp_ranges(triplet_count, run_alist_rdp4);
    report_progress(options, 0, 1, 0, triplet_count, 0);

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
    RdpRawEventState events;
    if (options.enable_rdp) {
        events = scan_rdp_redo_triplets(
            scan_state, distance_state, tree_state, redo, fss_rdp, store_lpv,
            store_lpv_upper_bound, fss_upper_bound, half_window, full_window,
            xover_settings, probability_settings, probability_estimate,
            factorials.three_way, factorials.factorial, xover_api(), 0, nullptr,
            nullptr, 0, nullptr, false, nullptr, options.progress_callback,
            options.progress_user, 0, 1);
    } else {
        events.current_xover.assign(scan_state.next_no + 1, 0);
        events.xover_list.resize(scan_state.next_no + 1);
        events.triplets_with_events.assign(static_cast<std::size_t>(triplet_count), 0);
    }
    const int selected_method_count =
        static_cast<int>(options.enable_rdp) +
        static_cast<int>(options.enable_geneconv) +
        static_cast<int>(options.enable_maxchi) +
        static_cast<int>(options.enable_chimaera) +
        static_cast<int>(options.enable_three_seq);
    // InnerScan2 schedules the optional methods only after the RDP event
    // selection has entered the full event loop.  Do not inject their
    // recheck records into the initial RDP event state: that state is fed
    // directly to Identify/MakeNJTrees and method records can make the
    // source tree pass invalid (notably MaxChi on Dataset2).  The retained
    // implementation below is kept as a fixture-only reference for the
    // earlier call-order experiments; the live path is the scheduler below.
    if (selected_method_count > static_cast<int>(options.enable_rdp)) {
        const bool has_geneconv_order =
            !options.geneconv_call_order_path.empty() &&
            !options.geneconv_call_count_path.empty();
        const bool has_maxchi_order =
            !options.maxchi_call_order_path.empty() &&
            !options.maxchi_call_count_path.empty();
        const bool has_chimaera_order =
            !options.chimaera_call_order_path.empty() &&
            !options.chimaera_call_count_path.empty();
        // 3SEQ is a post-RDP InnerScan2 method (the source calls TSXover
        // from its shared method loop rather than from the initial AlistRDP4
        // pass), so it is intentionally not added to this initial event
        // state.  Its source-order recheck is run below in the full event
        // loop alongside the other selectable methods.
        RdpLegacyEventAllocator allocator(events, selected_method_count);
        if (options.enable_geneconv) {
            if (has_geneconv_order) {
                std::ifstream count_input(options.geneconv_call_count_path,
                                           std::ios::binary);
                int call_count = 0;
                count_input.read(reinterpret_cast<char*>(&call_count),
                                 sizeof(call_count));
                std::ifstream order_input(options.geneconv_call_order_path,
                                          std::ios::binary);
                for (int call = 0; call < call_count; ++call) {
                    std::array<int, 9> record{};
                    order_input.read(reinterpret_cast<char*>(record.data()),
                                     sizeof(record));
                    if (!order_input) throw std::runtime_error(
                        "truncated GENECONV call-order fixture");
                    std::array<int, 3> input{record[5], record[6], record[7]};
                    run_rdp_geneconv_recheck(
                        scan_state, input, store_lpv, store_lpv_upper_bound,
                        probability_settings, allocator);
                }
            } else {
                const auto screened = screen_rdp_geneconv_candidates(
                    scan_state, store_lpv, store_lpv_upper_bound,
                    correction_tests, options.p_value_cutoff,
                    options.circular ? 1 : 0, 0, target);
                for (const auto& candidate : screened.candidates) {
                    auto input = candidate;
                    run_rdp_geneconv_recheck(
                        scan_state, input,
                        store_lpv, store_lpv_upper_bound, probability_settings,
                        allocator);
                }
            }
        }
        if (options.enable_maxchi) {
            if (has_maxchi_order) {
                std::ifstream count_input(options.maxchi_call_count_path,
                                           std::ios::binary);
                int call_count = 0;
                count_input.read(reinterpret_cast<char*>(&call_count), sizeof(call_count));
                std::ifstream order_input(options.maxchi_call_order_path,
                                          std::ios::binary);
                for (int call = 0; call < call_count; ++call) {
                    std::array<int, 8> record{};
                    order_input.read(reinterpret_cast<char*>(record.data()), sizeof(record));
                    if (!order_input) throw std::runtime_error("truncated MaxChi call-order fixture");
                    std::array<int, 3> input{record[4], record[5], record[6]};
                    run_rdp_maxchi_recheck(
                        scan_state, input, store_lpv, store_lpv_upper_bound,
                        probability_settings, allocator, 0, 0, true);
                }
            } else {
                const auto screened = screen_rdp_maxchi_candidates(
                    scan_state, store_lpv, store_lpv_upper_bound,
                    correction_tests, options.p_value_cutoff,
                    options.circular ? 1 : 0, 0);
                for (const auto& input : screened.candidates) {
                    run_rdp_maxchi_recheck(
                        scan_state, input, store_lpv, store_lpv_upper_bound,
                        probability_settings, allocator, 0, 0, true);
                }
            }
        }
        if (options.enable_chimaera) {
            if (has_chimaera_order) {
                std::ifstream count_input(options.chimaera_call_count_path,
                                           std::ios::binary);
                int call_count = 0;
                count_input.read(reinterpret_cast<char*>(&call_count), sizeof(call_count));
                std::ifstream order_input(options.chimaera_call_order_path,
                                          std::ios::binary);
                for (int call = 0; call < call_count; ++call) {
                    std::array<int, 8> record{};
                    order_input.read(reinterpret_cast<char*>(record.data()), sizeof(record));
                    if (!order_input) throw std::runtime_error("truncated Chimaera call-order fixture");
                    run_rdp_chimaera_recheck(
                        scan_state, {record[3], record[4], record[5]}, store_lpv,
                        store_lpv_upper_bound, probability_settings, allocator,
                        0, 0, true);
                }
            } else {
                const auto screened = screen_rdp_chimaera_candidates(
                    scan_state, store_lpv, store_lpv_upper_bound,
                    correction_tests, options.p_value_cutoff,
                    options.circular ? 1 : 0, 0);
                for (int index = 0; index <= scan_state.analysis_list_last;
                     ++index) {
                    const unsigned char redo_flags = screened.redo[index];
                    const std::array<int, 3> candidate{
                        scan_state.analysis_list[index * 3],
                        scan_state.analysis_list[index * 3 + 1],
                        scan_state.analysis_list[index * 3 + 2]};
                    const std::array<std::array<int, 3>, 3> orientations{{
                        candidate,
                        {candidate[1], candidate[2], candidate[0]},
                        {candidate[2], candidate[0], candidate[1]}}};
                    const std::array<unsigned char, 3> orientation_bits{
                        1, 4, 16};
                    for (int rotation = 0; rotation < 3; ++rotation) {
                        if ((redo_flags & orientation_bits[rotation]) == 0)
                            continue;
                        run_rdp_chimaera_recheck(
                            scan_state, orientations[rotation], store_lpv,
                            store_lpv_upper_bound, probability_settings,
                            allocator, 0, 0, true);
                    }
                }
            }
        }
    }
    report_progress(
        options, 0, 1, triplet_count, triplet_count, count_events(events));
    return {
        std::move(scan_state), std::move(events), std::move(store_lpv),
        store_lpv_upper_bound};
}

void append_rdp_events(RdpRawEventState& destination,
                       const RdpRawEventState& source) {
    if (destination.xover_list.size() < source.xover_list.size()) {
        destination.xover_list.resize(source.xover_list.size());
        destination.current_xover.resize(source.xover_list.size(), 0);
    }
    for (std::size_t row = 0; row < source.xover_list.size(); ++row) {
        destination.xover_list[row].insert(
            destination.xover_list[row].end(),
            source.xover_list[row].begin(), source.xover_list[row].end());
        destination.current_xover[row] = static_cast<std::int16_t>(
            destination.xover_list[row].size());
    }
}

std::vector<int> calculate_initial_sequence_sizes(
    const RdpScanState& scan_state) {
    std::vector<int> sizes(scan_state.next_no + 1, 0);
    const int stride = scan_state.sequence_length + 1;
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        for (int position = 1; position <= scan_state.sequence_length;
             ++position) {
            if (scan_state.sequence_data[position + sequence * stride] != 46) {
                ++sizes[sequence];
            }
        }
    }
    return sizes;
}

RdpFullAnalysisResult run_rdp_full_analysis(
    RdpScanState initial_scan_state,
    const RdpInitialAnalysisOptions& options,
    RdpFullAnalysisTrace* const analysis_trace = nullptr) {
    constexpr std::array<int, 6> comparison{1, 0, 0, 2, 2, 1};
    auto initial = run_rdp_initial_analysis(
        std::move(initial_scan_state), options);
    RdpFullAnalysisResult output;
    output.sequence_count = initial.alignment.next_no + 1;
    output.sequence_length = initial.alignment.sequence_length;
    output.triplet_count = initial.alignment.analysis_list_last + 1;
    for (const auto& row : initial.events.xover_list) {
        output.raw_candidate_count += static_cast<int>(row.size());
    }

    if (!options.enable_rdp) {
        // Optional-method-only scans use the shared initial alignment and
        // method screening tables, but do not seed the RDP cyclic
        // tract-erasure scheduler. Preserve each selected method's source
        // emission order as a directly reviewable event.
        for (const auto& row : initial.events.xover_list) {
            for (const auto& raw : row) {
                if (raw.program_flag == 0) continue;
                RdpFinalEvent event;
                event.event_number = static_cast<int>(output.events.size()) + 1;
                event.program = static_cast<int>(raw.program_flag);
                event.winning_role = 0;
                event.probability = raw.probability;
                event.beginning = raw.beginning;
                event.ending = raw.ending;
                event.representative_sequences = {
                    raw.daughter, raw.major_parent, raw.minor_parent};
                if (raw.profile_available != 0) {
                    for (int role = 0; role < 3; ++role) {
                        event.profile_sequences[role] = raw.profile_sequences[role];
                    }
                    event.profile_sequences_available = true;
                } else {
                    event.profile_sequences = event.representative_sequences;
                }
                event.method_target_role = 0;
                for (int role = 0; role < 3; ++role) {
                    event.sequence_groups[role].push_back(
                        event.representative_sequences[role]);
                }
                if (options.polish_breakpoints_with_burt) {
                    event.burt = run_rdp_burt(
                        initial.alignment, event.representative_sequences,
                        event.beginning, event.ending, options.circular, 20, 0);
                    event.burt_attempted = event.burt.attempted;
                    event.burt_applied = event.burt.trained;
                    if (event.burt.trained) {
                        event.beginning = event.burt.polished_beginning;
                        event.ending = event.burt.polished_ending;
                    }
                }
                output.events.push_back(std::move(event));
            }
        }
        append_legacy_optional_events(output, initial.alignment, options);
        apply_legacy_optional_rechecks(output, initial.alignment, options);
        return output;
    }

    const int permanent_next_no = initial.alignment.next_no;
    const auto permanent_analysis_scan = initial.alignment;
    // RDP derives MinSeqSize from the alignment length (VB CLng rounds the
    // one-percent threshold).  The old port used a fixed 20 here, which
    // allowed short, poorly comparable sequence pairs into the NJ pass.
    const int minimum_sequence_size = static_cast<int>(std::nearbyint(
        static_cast<double>(initial.alignment.sequence_length) / 100.0));
    const auto initially_screened_triplets =
        initial.events.triplets_with_events;
    auto scan_state = initial.alignment;
    auto events = initial.events;
    int shared_xdiffpos0 = 0;
    auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    auto tree_state = build_rdp_upgma_tree_state(
        scan_state.next_no, distance_state);
    std::vector<unsigned char> missing_data(scan_state.sequence_data.size(), 0);
    std::vector<int> trace_sub(scan_state.next_no + 1, 0);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        trace_sub[sequence] = sequence;
    }
    auto actual_sequence_sizes = calculate_initial_sequence_sizes(scan_state);
    std::vector<unsigned char> done_sequence;
    int done_row_upper_bound = -1;

    // The initial AlistRDP4/redo pass has completed at this point. Keep the
    // source triplet total, then let the cyclic loop report its round
    // boundaries without pretending that its internal DNA5 calls are
    // interruptible.
    report_progress(
        options, 0, 1, output.triplet_count, output.triplet_count,
        count_events(events));

    auto& fss_rdp = shared_fss_rdp();
    auto& factorials = shared_factorial_tables();
    auto& probability_estimate = shared_probability_estimate();
    const int correction_tests = output.triplet_count;
    const int target = static_cast<int>(std::nearbyint(
        static_cast<double>(scan_state.sequence_length) / 10.0));
    const short full_window = static_cast<short>(options.window_sites);
    const int half_window = static_cast<int>(std::nearbyint(
        static_cast<double>(options.window_sites) / 2.0));
    const RdpXoverSettings xover_settings{
        100, 1, target, options.circular ? 1 : 0};
    const RdpProbabilitySettings probability_settings{
        options.circular ? 1 : 0, correction_tests, 0, 0,
        171, 171, RdpFactorialTables::three_way_upper_bound,
        options.p_value_cutoff};
    RdpRescanScreenSettings rescan_settings;
    rescan_settings.circular = options.circular ? 1 : 0;
    rescan_settings.correction_tests = correction_tests;
    rescan_settings.probability_cutoff = options.p_value_cutoff;
    rescan_settings.target = target;
    rescan_settings.half_window = half_window;
    rescan_settings.full_window = full_window;

    for (int round_number = 1; round_number <= 1000; ++round_number) {
        report_progress(
            options, round_number == 1 ? 0 : 1, round_number,
            output.triplet_count, output.triplet_count,
            static_cast<int>(output.events.size()));
        auto selection = select_rdp_best_event(
            events, scan_state.next_no, options.p_value_cutoff,
            done_sequence, done_row_upper_bound);
        if (!selection.found) {
            append_legacy_optional_events(output, initial.alignment, options);
            apply_legacy_optional_rechecks(output, initial.alignment, options);
            return output;
        }
        if (selection.trace[0] < 0 ||
            selection.trace[0] >= static_cast<int>(events.xover_list.size()) ||
            selection.trace[1] < 1 ||
            selection.trace[1] > static_cast<int>(
                events.xover_list[selection.trace[0]].size())) {
            throw std::runtime_error("RDP selected an invalid event slot");
        }
        auto& selected_slot = events.xover_list[
            selection.trace[0]][selection.trace[1] - 1];
        if (selected_slot.daughter == selected_slot.minor_parent ||
            selected_slot.daughter == selected_slot.major_parent ||
            selected_slot.minor_parent == selected_slot.major_parent) {
            selected_slot.probability = 1.0;
            done_sequence = std::move(selection.done_sequence);
            done_row_upper_bound = scan_state.next_no;
            continue;
        }
        RdpRawEvent selected = selected_slot;
        if (analysis_trace != nullptr) {
            analysis_trace->selected_slots.push_back({
                selection.trace[0], selection.trace[1]});
            analysis_trace->done_before_selection.push_back(done_sequence);
            analysis_trace->done_row_upper_bounds.push_back(
                done_row_upper_bound);
        }
        auto round = identify_rdp_complete_round(
            scan_state, distance_state, events, selected, missing_data,
            permanent_next_no, minimum_sequence_size);
        if (std::getenv("RDP_TRACE_RLIST") != nullptr) {
            std::cerr << "rlist round=" << round_number << " winner="
                      << round.consensus.winning_role << " actual=";
            for (int role = 0; role < 3; ++role) {
                std::cerr << '[';
                for (int slot = 0;
                     slot <= round.prefix.actual_resolution.candidates.last[role];
                     ++slot) {
                    if (slot) std::cerr << ',';
                    std::cerr << round.prefix.actual_resolution.candidates.list[
                        role + slot * 3];
                }
                std::cerr << ']';
            }
            std::cerr << " final=";
            for (int role = 0; role < 3; ++role) {
                std::cerr << '[';
                for (int slot = 0;
                     slot <= round.final_candidates.candidate_last[role];
                     ++slot) {
                    if (slot) std::cerr << ',';
                    std::cerr << round.final_candidates.candidate_list[
                        role + slot * 3];
                }
                std::cerr << ']';
            }
            std::cerr << '\n';
        }
        if (std::getenv("RDP_TRACE_MUTATION_HEADERS") != nullptr) {
            std::cerr << "mutation-header round=" << round_number
                      << " next=" << scan_state.next_no
                      << " begin=" << selected.beginning
                      << " end=" << selected.ending
                      << " winner=" << round.consensus.winning_role
                      << " last=" << round.final_candidates.candidate_last[0]
                      << ':' << round.final_candidates.candidate_last[1]
                      << ':' << round.final_candidates.candidate_last[2]
                      << " list=";
            for (int role = 0; role < 3; ++role) {
                std::cerr << '[';
                for (int slot = 0;
                     slot <= round.final_candidates.candidate_last[role]; ++slot) {
                    if (slot) std::cerr << ',';
                    std::cerr << round.final_candidates.candidate_list[
                        role + slot * 3];
                }
                std::cerr << ']';
            }
            std::cerr << '\n';
        }
        RdpBurtResult burt_result;
        if (options.polish_breakpoints_with_burt) {
            // RDP calls PolishBP(20, 0, ...): the HMM iteration count is a
            // fixed 20, independent of MinSeqSize.  MinSeqSize is only the
            // sequence-size gate used by the surrounding NJ/rescan logic.
            // Passing minimum_sequence_size here accidentally made long
            // alignments run dozens (or hundreds) of extra HMM cycles.
            burt_result = run_rdp_burt(
                scan_state, round.prefix.sequences, selected.beginning,
                selected.ending, options.circular, 20, 0);
            if (burt_result.trained) {
                selected.beginning = static_cast<std::int16_t>(
                    burt_result.polished_beginning);
                selected.ending = static_cast<std::int16_t>(
                    burt_result.polished_ending);
            }
        }
        if (analysis_trace != nullptr) {
            analysis_trace->complete_round_states.push_back(round);
            analysis_trace->consensus_states.push_back(round.consensus);
            analysis_trace->consensus_candidate_states.push_back(
                round.consensus_candidates);
            analysis_trace->final_candidate_states.push_back(
                round.final_candidates);
        }
        missing_data = round.prefix.missing_data;
        const int winner = round.consensus.winning_role;
        if (std::getenv("RDP_TRACE_WINNER") != nullptr) {
            std::cerr << "winner-state next=" << scan_state.next_no
                      << " winner=" << winner << " consensus="
                      << round.consensus.consensus[0] << ':'
                      << round.consensus.consensus[1] << ':'
                      << round.consensus.consensus[2] << " tbreak="
                      << round.consensus.rounded_inputs.phylpro[0] -
                          round.consensus.rounded_inputs.triplet_score[0] +
                          round.consensus.rounded_inputs.compatibility[0] +
                          round.consensus.rounded_inputs.region_compatibility[0]
                      << ':'
                      << round.consensus.rounded_inputs.phylpro[1] -
                          round.consensus.rounded_inputs.triplet_score[1] +
                          round.consensus.rounded_inputs.compatibility[1] +
                          round.consensus.rounded_inputs.region_compatibility[1]
                      << ':'
                      << round.consensus.rounded_inputs.phylpro[2] -
                          round.consensus.rounded_inputs.triplet_score[2] +
                          round.consensus.rounded_inputs.compatibility[2] +
                          round.consensus.rounded_inputs.region_compatibility[2]
                      << '\n';
        }
        if (std::getenv("RDP_TRACE_CONSENSUS") != nullptr) {
            const auto& ci = round.consensus.rounded_inputs;
            std::cerr << "cons-input next=" << scan_state.next_no
                      << " seqs=" << round.prefix.sequences[0]
                      << ':' << round.prefix.sequences[1] << ':'
                      << round.prefix.sequences[2] << " winner=" << winner
                      << " dmax=" << ci.maximum_distance[0] << ':'
                      << ci.maximum_distance[1] << ':' << ci.maximum_distance[2]
                      << " ss=" << ci.split_distance[0] << ':'
                      << ci.split_distance[1] << ':' << ci.split_distance[2]
                      << " ph=" << ci.phylpro[0] << ':' << ci.phylpro[1]
                      << ':' << ci.phylpro[2] << " ph2="
                      << ci.phylpro_secondary[0] << ':'
                      << ci.phylpro_secondary[1] << ':'
                      << ci.phylpro_secondary[2] << " sub="
                      << ci.subtree_phylpro[0] << ':' << ci.subtree_phylpro[1]
                      << ':' << ci.subtree_phylpro[2] << " list="
                      << ci.list_correlation[0] << ':'
                      << ci.list_correlation[1] << ':'
                      << ci.list_correlation[2] << " bad="
                      << ci.bad_distances[0] << ':' << ci.bad_distances[1]
                      << ':' << ci.bad_distances[2] << " comp="
                      << ci.compatibility[0] << ':' << ci.compatibility[1]
                      << ':' << ci.compatibility[2] << " compS="
                      << ci.region_compatibility[0] << ':'
                      << ci.region_compatibility[1] << ':'
                      << ci.region_compatibility[2] << " post="
                      << ci.post_trim_compatibility[0] << ':'
                      << ci.post_trim_compatibility[1] << ':'
                      << ci.post_trim_compatibility[2] << ','
                      << ci.post_trim_region_compatibility[0] << ':'
                      << ci.post_trim_region_compatibility[1] << ':'
                      << ci.post_trim_region_compatibility[2]
                      << " min=" << round.prefix.minimum_pair[0] << ':'
                      << round.prefix.minimum_pair[1]
                      << " pair=" << round.prefix.sequence_pair[0] << ':'
                      << round.prefix.sequence_pair[1] << ':'
                      << round.prefix.sequence_pair[2] << '\n';
        }

        RdpFinalEvent final_event;
        final_event.event_number = static_cast<int>(output.events.size()) + 1;
        final_event.program = static_cast<int>(selected.program_flag);
        final_event.winning_role = winner;
        final_event.probability = selected.probability;
        final_event.beginning = selected.beginning;
        final_event.ending = selected.ending;
        final_event.consensus = round.consensus.consensus;
        final_event.burt = std::move(burt_result);
        final_event.burt_attempted = final_event.burt.attempted;
        final_event.burt_applied = final_event.burt.trained;
        if (selected.program_flag == 0) {
            // AddjustCXO/MakeCollecteventsC communicate through the legacy
            // XOVERDEFINE ABI, so their bookkeeping copies do not carry the
            // plot provenance fields.  In that case the complete-round
            // prefix is the live triplet that the selected event is about;
            // all post-initial cyclic XOver passes use PB3.
            RdpRawEvent plot_event = selected;
            if (plot_event.profile_available == 0) {
                plot_event.profile_sequences = {
                    static_cast<std::int16_t>(round.prefix.sequences[0]),
                    static_cast<std::int16_t>(round.prefix.sequences[1]),
                    static_cast<std::int16_t>(round.prefix.sequences[2])};
                plot_event.profile_available = 1;
                plot_event.profile_use_compress = 1;
            }
            final_event.rdp_profile = make_rdp_sliding_window_profile(
                scan_state, distance_state, tree_state, plot_event,
                options.window_sites, fss_rdp, xover_api());
            if (plot_event.profile_available != 0) {
                final_event.profile_sequences = {
                    plot_event.profile_sequences[0] < static_cast<int>(trace_sub.size())
                        ? trace_sub[plot_event.profile_sequences[0]]
                        : static_cast<int>(plot_event.profile_sequences[0]),
                    plot_event.profile_sequences[1] < static_cast<int>(trace_sub.size())
                        ? trace_sub[plot_event.profile_sequences[1]]
                        : static_cast<int>(plot_event.profile_sequences[1]),
                    plot_event.profile_sequences[2] < static_cast<int>(trace_sub.size())
                        ? trace_sub[plot_event.profile_sequences[2]]
                        : static_cast<int>(plot_event.profile_sequences[2])};
                final_event.profile_sequences_available = true;
            }
        }
        if (!final_event.profile_sequences_available &&
            selected.profile_available != 0) {
            for (int role = 0; role < 3; ++role) {
                const int working = selected.profile_sequences[role];
                final_event.profile_sequences[role] =
                    working >= 0 && working < static_cast<int>(trace_sub.size())
                        ? trace_sub[working] : working;
            }
            final_event.profile_sequences_available = true;
        }
        for (int role = 0; role < 3; ++role) {
            const int representative = round.prefix.sequences[role];
            final_event.representative_sequences[role] =
                representative < static_cast<int>(trace_sub.size())
                ? trace_sub[representative] : representative;
            for (int slot = 0;
                 slot <= round.final_candidates.candidate_last[role]; ++slot) {
                const int sequence = round.final_candidates.candidate_list[
                    role + slot * 3];
                final_event.sequence_groups[role].push_back(
                    sequence < static_cast<int>(trace_sub.size())
                    ? trace_sub[sequence] : sequence);
            }
        }
        if (!final_event.profile_sequences_available) {
            final_event.profile_sequences = final_event.representative_sequences;
        }
        if (selected.program_flag == 4) {
            for (int role = 0; role < 3; ++role) {
                if (final_event.representative_sequences[role] == selected.daughter) {
                    final_event.method_target_role = role;
                    break;
                }
            }
        }
        output.events.push_back(std::move(final_event));
        report_progress(
            options, round_number == 1 ? 0 : 1, round_number,
            output.triplet_count, output.triplet_count,
            static_cast<int>(output.events.size()));

        auto selected_collection_event = selected;
        selected_collection_event.distance_holder =
            (selected_collection_event.distance_holder + 0.00000001) * -1.0;
        if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
                    std::cerr << "temp-selected next=" << scan_state.next_no
                      << " program="
                      << static_cast<int>(selected_collection_event.program_flag)
                      << " daughter=" << selected_collection_event.daughter
                      << " parents=" << selected_collection_event.major_parent
                      << ':' << selected_collection_event.minor_parent
                      << " interval=" << selected_collection_event.beginning
                      << '-' << selected_collection_event.ending << '\n';
        }
        RdpRawEventState events_after_final_trim = events;
        events_after_final_trim.xover_list[selection.trace[0]][
            selection.trace[1] - 1] = selected_collection_event;
        RdpRawEventState temporary_events;
        temporary_events.current_xover.assign(scan_state.next_no + 1, 0);
        temporary_events.xover_list.resize(scan_state.next_no + 1);
        for (const auto& seed : round.final_candidates.synthetic_event_roles) {
            const int role = seed[0];
            const int sequence = seed[1];
            RdpRawEvent synthetic = selected_collection_event;
            synthetic.daughter = static_cast<std::int16_t>(sequence);
            synthetic.major_parent = static_cast<std::int16_t>(
                round.prefix.sequences[comparison[role]]);
            synthetic.minor_parent = static_cast<std::int16_t>(
                round.prefix.sequences[comparison[role + 3]]);
            synthetic.probability = 0.9;
            synthetic.distance_holder = std::abs(synthetic.distance_holder);
            temporary_events.xover_list[sequence].push_back(synthetic);
            temporary_events.current_xover[sequence] =
                static_cast<std::int16_t>(
                    temporary_events.xover_list[sequence].size());
        }
        std::vector<unsigned char> single_redo(
            static_cast<std::size_t>(scan_state.analysis_list_last + 1), 0);
        if (single_redo.empty()) single_redo.resize(1, 0);
        single_redo[0] = 1;
        auto relaxed_probability_settings = probability_settings;
        relaxed_probability_settings.lowest_probability = std::max<double>({
            probability_settings.lowest_probability,
            static_cast<double>(
                selected_collection_event.probability * 100000.0),
            static_cast<double>(
                probability_settings.lowest_probability * correction_tests)});
        for (int role = 0; role < 3; ++role) {
            for (int slot = 0;
                 slot <= round.final_candidates.candidate_last[role]; ++slot) {
                const int sequence = round.final_candidates.candidate_list[
                    role + slot * 3];
                if (sequence == round.prefix.sequences[role]) continue;
                const std::array<int, 3> triplet{
                    sequence,
                    round.prefix.sequences[comparison[role]],
                    round.prefix.sequences[comparison[role + 3]]};
                temporary_events = scan_rdp_redo_triplets(
                    scan_state, distance_state, tree_state, single_redo,
                    fss_rdp, initial.store_lpv,
                    initial.store_lpv_upper_bound, 125, half_window,
                    full_window, xover_settings,
                    relaxed_probability_settings, probability_estimate,
                    factorials.three_way, factorials.factorial, xover_api(), 1,
                    &temporary_events, &triplet, 1, &missing_data, false,
                    &shared_xdiffpos0);
            }
        }
        trace_legacy_method_state("after-rdp-rescans", temporary_events);
        trace_legacy_method_state("after-rdp-winning-rescan", temporary_events);
        // Module2.FinalTrim has two distinct passes.  Its first pass visits
        // every RList role but calls only the method that produced the
        // selected event.  Its second pass visits the winner list and calls
        // every other method.  Keeping those passes separate matters: the
        // second pass is intentionally run even when RDP.ini disables the
        // method (a long-standing RDP quirk).
        const auto& final_trim_three_seq_table = active_three_seq_table();
        RdpLegacyEventAllocator final_trim_legacy_events(temporary_events, 1);
        // First pass: all roles, selected method only.  The selected RDP
        // method is already represented by the source-faithful RDP rescan
        // above; the optional methods need their direct FinalTrim call.
        for (int final_trim_role = 0; final_trim_role < 3;
             ++final_trim_role) {
            for (int slot = 0;
                 slot <= round.final_candidates.candidate_last[final_trim_role];
                 ++slot) {
                const int sequence = round.final_candidates.candidate_list[
                    final_trim_role + slot * 3];
                if (sequence == round.prefix.sequences[final_trim_role]) {
                    continue;
                }
                std::array<int, 3> triplet{
                    sequence,
                    round.prefix.sequences[comparison[final_trim_role]],
                    round.prefix.sequences[comparison[final_trim_role + 3]]};
                if (selected_collection_event.program_flag == 1) {
                    run_rdp_geneconv_recheck(
                        scan_state, triplet, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        final_trim_legacy_events);
                }
                if (selected_collection_event.program_flag == 3) {
                    run_rdp_maxchi_recheck(
                        scan_state, triplet, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        final_trim_legacy_events, selected.beginning,
                        selected.ending);
                }
                if (selected_collection_event.program_flag == 4) {
                    auto rotated = triplet;
                    for (int rotation = 0; rotation < 3; ++rotation) {
                        run_rdp_chimaera_recheck(
                            scan_state, rotated, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            final_trim_legacy_events, selected.beginning,
                            selected.ending, false,
                            rotation < 2 ? &shared_xdiffpos0 : nullptr);
                        rotated = {rotated[1], rotated[2], rotated[0]};
                    }
                }
                if (selected_collection_event.program_flag == 8 &&
                    !final_trim_three_seq_table.empty()) {
                    auto rotated = triplet;
                    for (int rotation = 0; rotation < 3; ++rotation) {
                        run_rdp_three_seq_recheck(
                            scan_state, rotated, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            final_trim_three_seq_table,
                            configured_three_seq_table != nullptr
                                ? configured_three_seq_table_bound : 45,
                            final_trim_legacy_events);
                        rotated = {rotated[1], rotated[2], rotated[0]};
                    }
                }
            }
        }
        // Second pass: winner list, every method except the selected one.
        for (int slot = 0; slot <= round.final_candidates.candidate_last[winner];
             ++slot) {
            const int sequence = round.final_candidates.candidate_list[
                winner + slot * 3];
            std::array<int, 3> triplet{
                sequence,
                round.prefix.sequences[comparison[winner]],
                round.prefix.sequences[comparison[winner + 3]]};
            if (selected_collection_event.program_flag != 1) {
                run_rdp_geneconv_recheck(
                    scan_state, triplet, initial.store_lpv,
                    initial.store_lpv_upper_bound, probability_settings,
                    final_trim_legacy_events);
            }
            if (selected_collection_event.program_flag != 3) {
                run_rdp_maxchi_recheck(
                    scan_state, triplet, initial.store_lpv,
                    initial.store_lpv_upper_bound, probability_settings,
                    final_trim_legacy_events, selected.beginning,
                    selected.ending);
            }
            if (selected_collection_event.program_flag != 4) {
                auto rotated = triplet;
                for (int rotation = 0; rotation < 3; ++rotation) {
                    run_rdp_chimaera_recheck(
                        scan_state, rotated, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        final_trim_legacy_events, selected.beginning,
                        selected.ending, false,
                        rotation < 2 ? &shared_xdiffpos0 : nullptr);
                    rotated = {rotated[1], rotated[2], rotated[0]};
                }
            }
            if (selected_collection_event.program_flag != 8 &&
                !final_trim_three_seq_table.empty()) {
                auto rotated = triplet;
                for (int rotation = 0; rotation < 3; ++rotation) {
                    run_rdp_three_seq_recheck(
                        scan_state, rotated, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        final_trim_three_seq_table,
                        configured_three_seq_table != nullptr
                            ? configured_three_seq_table_bound : 45,
                        final_trim_legacy_events);
                    rotated = {rotated[1], rotated[2], rotated[0]};
                }
            }
        }
        if (std::getenv("RDP_TRACE_TEMP") != nullptr && winner == 1 &&
            round.final_candidates.candidate_last[winner] >= 1) {
            const auto lpv = [&](int program, int sequence) {
                return initial.store_lpv[static_cast<std::size_t>(sequence) *
                    static_cast<std::size_t>(initial.store_lpv_upper_bound + 1) +
                    static_cast<std::size_t>(program)];
            };
            std::cerr << "temp-lpv winner=" << winner;
            for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
                if (lpv(0, sequence) != 1.0 || lpv(3, sequence) != 1.0 ||
                    lpv(4, sequence) != 1.0 || lpv(8, sequence) != 1.0) {
                    std::cerr << ' ' << sequence << "=" << lpv(0, sequence)
                              << '/' << lpv(3, sequence) << '/' << lpv(4, sequence)
                              << '/' << lpv(8, sequence);
                }
            }
            std::cerr << '\n';
        }
        const int selected_method_count =
            1 + static_cast<int>(options.enable_geneconv) +
            static_cast<int>(options.enable_maxchi) +
            static_cast<int>(options.enable_chimaera) +
            static_cast<int>(options.enable_three_seq);
        RdpLegacyEventAllocator legacy_events(
            temporary_events, selected_method_count);
        const auto& three_seq_table = active_three_seq_table();
        // InnerScan2 does not feed the legacy methods only the final winner
        // list.  It expands the permanent AnalysisList through the selected
        // WinPP RList using MakeAListISP3, then runs each method's Alist* and
        // Xover routine in source order.  Keep the scheduler shared by all
        // methods so random method combinations exercise the same list.
        const int method_bits = 1 |
            (options.enable_geneconv ? 2 : 0) |
            (options.enable_maxchi ? 8 : 0) |
            (options.enable_chimaera ? 16 : 0) |
            (options.enable_three_seq ? 64 : 0);
        // AddjustCXO is the source boundary immediately before InnerScan2.
        // Its DoPairs matrix is not an all-ones matrix: it is derived from
        // the collected PXOList records and then propagated through the
        // selected RList.  Compute that matrix before constructing any
        // method AList so a method combination sees the same pair gates as
        // the VB implementation.  The full adjustment below is still run
        // after method emissions have been appended, as in the existing
        // event-state pipeline.
        const auto method_collection_events = prepare_rdp_collection_event_list(
            scan_state.next_no, winner, round.prefix.sequences, selection.trace,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list,
            round.final_candidates.acceptable_sequences,
            events_after_final_trim);
        auto method_adjusted = adjust_rdp_events_exact(
            scan_state.next_no, winner, options.p_value_cutoff,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, trace_sub,
            method_collection_events, selection.done_sequence,
            scan_state.next_no, selection.slot_upper_bound);
        auto method_pairs = method_adjusted.pairs_to_rescan;
        propagate_rdp_group_pairs(
            scan_state.next_no, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, method_pairs);
        const auto& method_rlist = round.prefix.actual_resolution.candidates;
        const auto method_triplets = [&](const int program) {
            return make_rdp_inner_method_triplets(
                scan_state, method_rlist.last, method_rlist.list, winner,
                trace_sub, actual_sequence_sizes, method_pairs,
                permanent_next_no, minimum_sequence_size, program, method_bits);
        };
        const auto run_method_scan = [&](const int program,
                                         const bool screen_geneconv,
                                         const bool screen_maxchi) {
            auto triplets = method_triplets(program);
            RdpScanState method_scan = scan_state;
            method_scan.analysis_list = flatten_rdp_triplets(triplets);
            method_scan.analysis_list_last =
                static_cast<int>(triplets.size()) - 1;
            if (method_scan.analysis_list_last < 0) return;
            if (screen_geneconv) {
                const auto screened = screen_rdp_geneconv_candidates(
                    method_scan, initial.store_lpv,
                    initial.store_lpv_upper_bound, correction_tests,
                    options.p_value_cutoff, options.circular ? 1 : 0, 0,
                    target);
                for (const auto& candidate : screened.candidates) {
                    auto input = candidate;
                    run_rdp_geneconv_recheck(
                        scan_state, input, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        legacy_events);
                }
            } else if (screen_maxchi) {
                const auto screened = screen_rdp_maxchi_candidates(
                    method_scan, initial.store_lpv,
                    initial.store_lpv_upper_bound, correction_tests,
                    options.p_value_cutoff, options.circular ? 1 : 0, 0);
                for (const auto& candidate : screened.candidates) {
                    run_rdp_maxchi_recheck(
                        scan_state, candidate, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        legacy_events, selected.beginning, selected.ending);
                }
            } else if (program == 4) {
                // CXoverA is entered only for AlistChi redo==1 rows.  The
                // source then gates its three role orientations through
                // ProgBinRead; until that per-redo orientation table is
                // threaded through this API, keep the source rotation order
                // but do not recheck every raw MakeAList row.
                const auto screened = screen_rdp_chimaera_candidates(
                    method_scan, initial.store_lpv,
                    initial.store_lpv_upper_bound, correction_tests,
                    options.p_value_cutoff, options.circular ? 1 : 0, 0);
                for (int index = 0; index <= method_scan.analysis_list_last;
                     ++index) {
                    const unsigned char redo_flags = screened.redo[index];
                    if ((redo_flags & 1U) == 0 &&
                        (redo_flags & 4U) == 0 &&
                        (redo_flags & 16U) == 0) {
                        continue;
                    }
                    const std::array<int, 3> candidate{
                        method_scan.analysis_list[index * 3],
                        method_scan.analysis_list[index * 3 + 1],
                        method_scan.analysis_list[index * 3 + 2]};
                    const std::array<std::array<int, 3>, 3> orientations{{
                        candidate,
                        {candidate[1], candidate[2], candidate[0]},
                        {candidate[2], candidate[0], candidate[1]}}};
                    const std::array<unsigned char, 3> orientation_bits{
                        1, 4, 16};
                    for (int rotation = 0; rotation < 3; ++rotation) {
                        if ((redo_flags & orientation_bits[rotation]) == 0)
                            continue;
                        const auto& rotated = orientations[rotation];
                        run_rdp_chimaera_recheck(
                            scan_state, rotated, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            legacy_events, selected.beginning, selected.ending,
                            false, rotation < 2 ? &shared_xdiffpos0 : nullptr);
                    }
                }
            } else if (program == 6) {
                for (const auto& candidate : triplets) {
                    auto rotated = candidate;
                    for (int rotation = 0; rotation < 3; ++rotation) {
                        run_rdp_three_seq_recheck(
                            scan_state, rotated, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            three_seq_table,
                            configured_three_seq_table != nullptr
                                ? configured_three_seq_table_bound : 45,
                            legacy_events, false, false);
                        rotated = {rotated[1], rotated[2], rotated[0]};
                    }
                }
            }
        };
        if (options.enable_geneconv) run_method_scan(1, true, false);
        if (options.enable_maxchi) run_method_scan(3, false, true);
        if (options.enable_chimaera) run_method_scan(4, false, false);
        if (options.enable_three_seq) run_method_scan(6, false, false);
        /*
        for (int slot = 0;
             slot <= round.final_candidates.candidate_last[winner]; ++slot) {
            const int sequence = round.final_candidates.candidate_list[
                winner + slot * 3];
            if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
                std::cerr << "temp-side-prefix winner=" << winner
                          << " prefix=" << round.prefix.sequences[0] << ':'
                          << round.prefix.sequences[1] << ':'
                          << round.prefix.sequences[2]
                          << " comparison=" << comparison[0] << ','
                          << comparison[1] << ',' << comparison[2] << ','
                          << comparison[3] << ',' << comparison[4] << ','
                          << comparison[5] << " candidate=" << sequence
                          << '\n';
            }
            std::array<int, 3> triplet{
                sequence,
                round.prefix.sequences[comparison[winner]],
                round.prefix.sequences[comparison[winner + 3]]};
            const auto trace_call = [&](const char* method, const int rotation) {
                if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
                    std::cerr << "temp-call " << method << " rotation="
                              << rotation << " triplet=" << triplet[0] << ':'
                              << triplet[1] << ':' << triplet[2] << '\n';
                }
            };
            if (options.enable_geneconv) {
                trace_call("geneconv", 0);
                run_rdp_geneconv_recheck(
                    scan_state, triplet, initial.store_lpv,
                    initial.store_lpv_upper_bound, probability_settings,
                    legacy_events);
                trace_legacy_method_state("after-geneconv", temporary_events);
            }
            if (options.enable_maxchi) {
                trace_call("maxchi", 0);
                run_rdp_maxchi_recheck(
                    scan_state, triplet, initial.store_lpv,
                    initial.store_lpv_upper_bound, probability_settings,
                    legacy_events, selected.beginning, selected.ending);
                trace_legacy_method_state("after-maxchi", temporary_events);
            }
            if (options.enable_chimaera) {
                auto rotated = triplet;
                for (int rotation = 0; rotation < 3; ++rotation) {
                    if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
                        std::cerr << "temp-call chimaera rotation=" << rotation
                                  << " triplet=" << rotated[0] << ':'
                                  << rotated[1] << ':' << rotated[2] << '\n';
                    }
                    // Preserve the native CXoverA map hand-off: its third
                    // orientation does not replace the map used by later XOver.
                    run_rdp_chimaera_recheck(
                        scan_state, rotated, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        legacy_events, selected.beginning, selected.ending, false,
                        rotation < 2 ? &shared_xdiffpos0 : nullptr);
                    trace_legacy_method_state("after-chimaera", temporary_events);
                    rotated = {rotated[1], rotated[2], rotated[0]};
                }
            }
            if (options.enable_three_seq) {
                auto rotated = triplet;
                for (int rotation = 0; rotation < 3; ++rotation) {
                    if (std::getenv("RDP_TRACE_TEMP") != nullptr) {
                        std::cerr << "temp-call three-seq rotation=" << rotation
                                  << " triplet=" << rotated[0] << ':'
                                  << rotated[1] << ':' << rotated[2] << '\n';
                    }
                    run_rdp_three_seq_recheck(
                        scan_state, rotated, initial.store_lpv,
                        initial.store_lpv_upper_bound, probability_settings,
                        three_seq_table, configured_three_seq_table != nullptr
                            ? configured_three_seq_table_bound : 45,
                        legacy_events);
                    trace_legacy_method_state("after-three-seq", temporary_events);
                    rotated = {rotated[1], rotated[2], rotated[0]};
                }
            }
        }
        */
        append_rdp_events(events_after_final_trim, temporary_events);
        trace_legacy_method_state("final-trim-merged", events_after_final_trim);
        // FinalTrim's rescans above use the original selected probability.
        // MakeCollecteventsX also receives that original PXOList record.  The
        // source marks the persistent selected signal as consumed only after
        // collection, before MakeBestXOList/AddjustCXO.
        const std::array<int, 2> trace{
            selection.trace[0], selection.trace[1]};
        // FinalTrim adds a probability-1 copy of the selected PXO record for
        // each accepted member of the winning RList before MakeCollecteventsC
        // is called.  `prepare_rdp_collection_event_list` is the source-order
        // representation of that persistent PXOList mutation; using the raw
        // pre-copy list here silently loses those records and makes every
        // later WinnerPos/StripUnfound decision wrong.
        auto collection_events = prepare_rdp_collection_event_list(
            scan_state.next_no, winner, round.prefix.sequences, trace,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list,
            round.final_candidates.acceptable_sequences,
            events_after_final_trim);
        if (std::getenv("RDP_TRACE_STRIP") != nullptr) {
            const auto collected = make_rdp_parent_collect_events(
                scan_state.sequence_length, scan_state.next_no, winner,
                round.prefix.actual_resolution.region_sizes,
                round.prefix.actual_resolution.event_overlap_mask,
                comparison,
                round.final_candidates.candidate_last,
                round.final_candidates.candidate_list, 0,
                round.prefix.sequences, selection.trace,
                collection_events);
            std::cerr << "strip-candidates round=" << round_number
                      << " winner=" << winner << " source=";
            for (int slot = 0; slot <= round.final_candidates.candidate_last[winner]; ++slot) {
                const auto sequence = round.final_candidates.candidate_list[winner + slot * 3];
                const auto& event = collected.events[slot];
                const auto ok_index = static_cast<std::size_t>(winner) +
                    18U * 3U + static_cast<std::size_t>(sequence) * 57U;
                std::cerr << sequence << ':' << event.probability << ":ok="
                          << round.final_candidates.acceptable_sequences[ok_index] << ',';
            }
            std::cerr << '\n';
        }
        events_after_final_trim.xover_list[selection.trace[0]][
            selection.trace[1] - 1].probability = 1.0;
        if (analysis_trace != nullptr) {
            analysis_trace->collection_events_before_adjustment.push_back(
                collection_events);
            trace_legacy_method_state("collection", collection_events);
        }

        const int selected_count =
            round.final_candidates.candidate_last[winner] + 1;
        std::vector<int> breakpoints(
            static_cast<std::size_t>(selected_count) * 2, 0);
        for (int slot = 0; slot < selected_count; ++slot) {
            breakpoints[slot * 2] = selected.beginning;
            breakpoints[slot * 2 + 1] = selected.ending;
        }
        auto mutation = erase_rdp_recombinant_tracts(
            scan_state.sequence_length, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, selected.beginning,
            selected.ending, breakpoints, scan_state.sequence_data,
            missing_data);
        const int expanded_next_no = scan_state.next_no + selected_count;
        auto expanded_sequences = mutation.sequence_data;
        auto expanded_missing = mutation.missing_data;
        rebuild_rdp_recombinant_tracts(
            scan_state.sequence_length, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, mutation.breakpoints,
            mutation.saved_tracts, expanded_sequences);
        expanded_sequences.resize(
            static_cast<std::size_t>(expanded_next_no + 1) *
                (scan_state.sequence_length + 1), 0);
        expanded_missing.resize(expanded_sequences.size(), 0);
        auto tail_breakpoints = mutation.breakpoints;
        make_rdp_fragment_rows(
            scan_state.sequence_length, expanded_next_no, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, selected.beginning,
            selected.ending, tail_breakpoints, expanded_sequences,
            expanded_missing);
        std::vector<int> compact_candidates(selected_count, 0);
        for (int slot = 0; slot < selected_count; ++slot) {
            compact_candidates[slot] = round.final_candidates.candidate_list[
                winner + slot * 3];
        }
        erase_rdp_original_tracts(
            scan_state.sequence_length, expanded_next_no, winner,
            round.final_candidates.candidate_last, compact_candidates,
            selected.beginning, selected.ending, tail_breakpoints,
            expanded_sequences, expanded_missing);
        const auto expanded_sizes = calculate_rdp_actual_sequence_sizes(
            scan_state.sequence_length, expanded_next_no, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, expanded_sequences);
        if (analysis_trace != nullptr) {
            analysis_trace->expanded_actual_sequence_sizes.push_back(
                expanded_sizes);
        }
        auto expanded_trace_sub = trace_sub;
        expanded_trace_sub.resize(expanded_next_no + 1, 0);
        for (int slot = 0; slot < selected_count; ++slot) {
            const int created = expanded_next_no -
                round.final_candidates.candidate_last[winner] + slot;
            expanded_trace_sub[created] =
                trace_sub[compact_candidates[slot]];
        }
        if (analysis_trace != nullptr) {
            analysis_trace->trace_sub_after_expansion.push_back(
                expanded_trace_sub);
        }

        auto adjusted = adjust_rdp_events_exact(
            scan_state.next_no, winner, options.p_value_cutoff,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, expanded_trace_sub,
            collection_events, selection.done_sequence,
            scan_state.next_no, selection.slot_upper_bound);
        if (analysis_trace != nullptr) {
            analysis_trace->adjusted_events.push_back(adjusted.events);
            analysis_trace->adjusted_pairs_to_rescan.push_back(
                adjusted.pairs_to_rescan);
        }
        auto propagated_pairs = adjusted.pairs_to_rescan;
        propagate_rdp_group_pairs(
            scan_state.next_no, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, propagated_pairs);
        const auto inner_triplets = make_rdp_inner_scan_triplets(
            permanent_analysis_scan, initially_screened_triplets, winner,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list, expanded_trace_sub,
            actual_sequence_sizes, permanent_next_no, minimum_sequence_size,
            propagated_pairs);
        const auto expanded_scan_state = rebuild_rdp_scan_state(
            expanded_next_no, scan_state.sequence_length,
            expanded_sequences, preprocessing_api());
        auto next_actual_sizes = actual_sequence_sizes;
        next_actual_sizes.resize(expanded_next_no + 1, 0);
        for (int slot = 0; slot < selected_count; ++slot) {
            const int source = compact_candidates[slot];
            const int created = expanded_next_no -
                round.final_candidates.candidate_last[winner] + slot;
            next_actual_sizes[source] = expanded_sizes[source];
            next_actual_sizes[created] = expanded_sizes[created];
        }
        const auto expanded_distance = build_rdp_distance_state(
            expanded_scan_state, 1, expanded_scan_state.sequence_length);
        // OuterScan4 re-dimensions SubValid for the post-mutation sequence
        // set and initializes every off-diagonal pair to 100.  It does not
        // reuse FastDistanceCalcZ's raw valid-site matrix here.
        std::vector<float> outer_subvalid(
            static_cast<std::size_t>(expanded_next_no + 1) *
                (expanded_next_no + 1), 0.0F);
        for (int first = 0; first < expanded_next_no; ++first) {
            for (int second = first + 1; second <= expanded_next_no; ++second) {
                outer_subvalid[first + second * (expanded_next_no + 1)] = 100.0F;
                outer_subvalid[second + first * (expanded_next_no + 1)] = 100.0F;
            }
        }
        const auto inner_scan_state = rebuild_rdp_scan_state(
            scan_state.next_no, scan_state.sequence_length,
            mutation.sequence_data, preprocessing_api());
        const auto inner_redo = screen_rdp_rescan_triplets(
            inner_triplets, inner_scan_state, distance_state, tree_state,
            rescan_settings, fss_rdp, probability_estimate,
            factorials.three_way, factorials.factorial);
        auto inner_event_state = inner_scan_state;
        inner_event_state.analysis_list = flatten_rdp_triplets(inner_triplets);
        inner_event_state.analysis_list_last =
            static_cast<int>(inner_triplets.size()) - 1;
        auto inner_events = scan_rdp_redo_triplets(
            inner_event_state, distance_state, tree_state, inner_redo,
            fss_rdp, initial.store_lpv, initial.store_lpv_upper_bound, 125,
            half_window, full_window, xover_settings, probability_settings,
            probability_estimate, factorials.three_way,
            factorials.factorial, xover_api(), 0, nullptr, nullptr, 1,
            &mutation.missing_data, true, &shared_xdiffpos0);

        // Module8.InnerScan2 runs the optional method branches against the
        // same permanent AnalysisList/RList expansion as the RDP branch,
        // after the selected tract has been erased.  The earlier FinalTrim
        // pass above is a separate source pass; omitting this post-mutation
        // pass drops events that only become significant on the edited
        // sequence (for example the MaxChi event at 71--1220 on Dataset0).
        if (selected_method_count > 1) {
            RdpLegacyEventAllocator inner_method_events(
                inner_events, selected_method_count);
            const auto inner_method_triplets = [&](const int program) {
                auto& selected_pairs = propagated_pairs;
                if (std::getenv("RDP_TRACE_INNER_ALLPAIRS") != nullptr &&
                    (round_number == 7 || round_number == 9)) {
                    selected_pairs.assign(
                        static_cast<std::size_t>(scan_state.next_no + 1) *
                            static_cast<std::size_t>(scan_state.next_no + 1), 1);
                }
                return make_rdp_inner_method_triplets(
                    permanent_analysis_scan,
                    round.final_candidates.candidate_last,
                    round.final_candidates.candidate_list, winner,
                    trace_sub, actual_sequence_sizes, selected_pairs,
                    permanent_next_no, minimum_sequence_size, program,
                    method_bits);
            };
            const auto run_inner_method_scan = [&](const int program,
                                                   const bool screen_geneconv,
                                                   const bool screen_maxchi) {
                const auto triplets = inner_method_triplets(program);
                auto method_scan = inner_scan_state;
                method_scan.analysis_list = flatten_rdp_triplets(triplets);
                method_scan.analysis_list_last =
                    static_cast<int>(triplets.size()) - 1;
                if (method_scan.analysis_list_last < 0) return;
                if (screen_geneconv) {
                    const auto screened = screen_rdp_geneconv_candidates(
                        method_scan, initial.store_lpv,
                        initial.store_lpv_upper_bound, correction_tests,
                        options.p_value_cutoff, options.circular ? 1 : 0, 0,
                        target);
                    for (const auto& candidate : screened.candidates) {
                        auto input = candidate;
                        run_rdp_geneconv_recheck(
                            inner_scan_state, input, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            inner_method_events);
                    }
                } else if (screen_maxchi) {
                    const auto screened = screen_rdp_maxchi_candidates(
                        method_scan, initial.store_lpv,
                        initial.store_lpv_upper_bound, correction_tests,
                        options.p_value_cutoff, options.circular ? 1 : 0, 0,
                        round_number, &mutation.missing_data);
                    for (const auto& candidate : screened.candidates) {
                        run_rdp_maxchi_recheck(
                            inner_scan_state, candidate, initial.store_lpv,
                            initial.store_lpv_upper_bound, probability_settings,
                            inner_method_events, 0, 0, false, nullptr, true,
                            &mutation.missing_data);
                    }
                } else if (program == 4) {
                    const auto screened = screen_rdp_chimaera_candidates(
                        method_scan, initial.store_lpv,
                        initial.store_lpv_upper_bound, correction_tests,
                        options.p_value_cutoff, options.circular ? 1 : 0, 0,
                        round_number, &mutation.missing_data);
                    for (int index = 0;
                         index <= method_scan.analysis_list_last; ++index) {
                        const unsigned char redo_flags = screened.redo[index];
                        const std::array<int, 3> candidate{
                            method_scan.analysis_list[index * 3],
                            method_scan.analysis_list[index * 3 + 1],
                            method_scan.analysis_list[index * 3 + 2]};
                        const std::array<std::array<int, 3>, 3> orientations{{
                            candidate,
                            {candidate[1], candidate[2], candidate[0]},
                            {candidate[2], candidate[0], candidate[1]}}};
                        const std::array<unsigned char, 3> orientation_bits{
                            1, 4, 16};
                        for (int rotation = 0; rotation < 3; ++rotation) {
                            if ((redo_flags & orientation_bits[rotation]) == 0)
                                continue;
                            run_rdp_chimaera_recheck(
                                inner_scan_state, orientations[rotation],
                                initial.store_lpv,
                                initial.store_lpv_upper_bound,
                                probability_settings, inner_method_events,
                                0, 0, false,
                                rotation < 2 ? &shared_xdiffpos0 : nullptr,
                                true, &mutation.missing_data);
                        }
                    }
                } else if (program == 6) {
                    for (const auto& candidate : triplets) {
                        auto rotated = candidate;
                        for (int rotation = 0; rotation < 3; ++rotation) {
                            run_rdp_three_seq_recheck(
                                inner_scan_state, rotated, initial.store_lpv,
                                initial.store_lpv_upper_bound,
                                probability_settings, three_seq_table,
                                configured_three_seq_table != nullptr
                                    ? configured_three_seq_table_bound : 45,
                                inner_method_events);
                            rotated = {rotated[1], rotated[2], rotated[0]};
                        }
                    }
                }
            };
            if (options.enable_geneconv) run_inner_method_scan(1, true, false);
            if (options.enable_maxchi) run_inner_method_scan(3, false, true);
            if (options.enable_chimaera) run_inner_method_scan(4, false, false);
            if (options.enable_three_seq) run_inner_method_scan(6, false, false);
            trace_legacy_method_state("inner-methods", inner_events);
        }

        // FindActualEventsVB calls MakeBestXOList, then InnerScan2, then
        // StripUnfound2.  WinnerPos is populated from MakeCollecteventsC;
        // it is not inferred from InnerScan2's emitted events.  Preserve the
        // pre-strip RList for all preceding work and use this source-order
        // copy only when constructing OuterScan4's AList.
        auto outer_candidate_last = round.final_candidates.candidate_last;
        auto outer_candidate_list = round.final_candidates.candidate_list;
        constexpr int source_method_upper_bound = 9;
        const auto collected_for_winner = make_rdp_parent_collect_events(
            scan_state.sequence_length, scan_state.next_no, winner,
            round.prefix.actual_resolution.region_sizes,
            round.prefix.actual_resolution.event_overlap_mask,
            comparison,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list,
            source_method_upper_bound, round.prefix.sequences, selection.trace,
            collection_events);
        const int winner_count =
            round.final_candidates.candidate_last[winner] + 1;
        const int winner_method_count = source_method_upper_bound + 1;
        std::vector<unsigned char> winner_pos(
            static_cast<std::size_t>(std::max(0, winner_count)) *
                static_cast<std::size_t>(winner_method_count), 0);
        for (int slot = 0; slot < winner_count; ++slot) {
            for (int method = 0; method <= source_method_upper_bound; ++method) {
                const auto& event = collected_for_winner.events[
                    static_cast<std::size_t>(slot) +
                    static_cast<std::size_t>(method) *
                        static_cast<std::size_t>(scan_state.next_no + 1)];
                if (event.probability > 0.0) {
                    winner_pos[static_cast<std::size_t>(slot) *
                                   static_cast<std::size_t>(winner_method_count) +
                               static_cast<std::size_t>(method)] = 1;
                }
            }
        }
        int strip_slot = 0;
        while (strip_slot <= outer_candidate_last[winner]) {
            bool found = false;
            for (int method = 0; method <= source_method_upper_bound; ++method) {
                if (winner_pos[static_cast<std::size_t>(strip_slot) *
                                   static_cast<std::size_t>(winner_method_count) +
                               static_cast<std::size_t>(method)] != 0) {
                    found = true;
                    break;
                }
            }
            if (found) {
                ++strip_slot;
                continue;
            }
            const int last_slot = outer_candidate_last[winner];
            if (strip_slot < last_slot) {
                outer_candidate_list[winner + strip_slot * 3] =
                    outer_candidate_list[winner + last_slot * 3];
                for (int method = 0; method <= source_method_upper_bound; ++method) {
                    winner_pos[static_cast<std::size_t>(strip_slot) *
                                   static_cast<std::size_t>(winner_method_count) +
                               static_cast<std::size_t>(method)] =
                        winner_pos[static_cast<std::size_t>(last_slot) *
                                       static_cast<std::size_t>(winner_method_count) +
                                   static_cast<std::size_t>(method)];
                }
            }
            --outer_candidate_last[winner];
        }
        if (std::getenv("RDP_TRACE_STRIP") != nullptr) {
            std::cerr << "strip-output round=" << round_number << " role="
                      << winner << " last=" << outer_candidate_last[winner]
                      << " list=";
            for (int slot = 0; slot <= outer_candidate_last[winner]; ++slot) {
                std::cerr << outer_candidate_list[winner + slot * 3] << ',';
            }
            std::cerr << '\n';
        }
        const auto outer_triplets = make_rdp_outer_scan_triplets(
            permanent_analysis_scan, initially_screened_triplets,
            expanded_next_no,
            scan_state.next_no, winner,
            outer_candidate_last, expanded_trace_sub,
            next_actual_sizes, permanent_next_no, minimum_sequence_size,
            propagated_pairs,
            expanded_next_no, outer_subvalid);
        if (analysis_trace != nullptr) {
            analysis_trace->inner_triplets.push_back(inner_triplets);
            analysis_trace->outer_triplets.push_back(outer_triplets);
        }
        if (std::getenv("RDP_TRACE_WINNERPOS") != nullptr) {
            std::cerr << "winnerpos round=" << round_number << " role="
                      << winner << " candidates=";
            for (int slot = 0;
                 slot <= round.final_candidates.candidate_last[winner];
                 ++slot) {
                const int sequence = round.final_candidates.candidate_list[
                    winner + slot * 3];
                const int count = sequence < static_cast<int>(
                    inner_events.xover_list.size())
                    ? static_cast<int>(inner_events.xover_list[sequence].size()) : 0;
                std::cerr << sequence << ':' << count << ',';
            }
            std::cerr << '\n';
        }
        auto combined_events = adjusted.events;
        append_rdp_events(combined_events, inner_events);
        if (std::getenv("RDP_TRACE_COMBINED_BROAD") != nullptr &&
            round_number <= 3) {
            for (std::size_t row = 0; row < combined_events.xover_list.size(); ++row) {
                for (const auto& event : combined_events.xover_list[row]) {
                    if (event.beginning >= 1400 && event.beginning <= 2300 &&
                        event.ending >= 3500 && event.ending <= 4600) {
                        std::cerr << "combined-broad round=" << round_number
                                  << " row=" << row << " method="
                                  << static_cast<int>(event.program_flag)
                                  << " roles=" << event.daughter << ':'
                                  << event.minor_parent << ':' << event.major_parent
                                  << " interval=" << event.beginning << '-'
                                  << event.ending << " p=" << event.probability
                                  << '\n';
                    }
                }
            }
        }
        combined_events.xover_list.resize(expanded_next_no + 1);
        combined_events.current_xover.resize(expanded_next_no + 1, 0);
        const auto expanded_tree = build_rdp_upgma_tree_state(
            expanded_next_no, expanded_distance);
        const auto outer_redo = screen_rdp_rescan_triplets(
            outer_triplets, expanded_scan_state, expanded_distance,
            expanded_tree, rescan_settings, fss_rdp, probability_estimate,
            factorials.three_way, factorials.factorial);
        if (analysis_trace != nullptr) {
            analysis_trace->inner_redo.push_back(inner_redo);
            analysis_trace->outer_redo.push_back(outer_redo);
        }
        auto outer_event_state = expanded_scan_state;
        outer_event_state.analysis_list = flatten_rdp_triplets(outer_triplets);
        outer_event_state.analysis_list_last =
            static_cast<int>(outer_triplets.size()) - 1;
        std::vector<double> expanded_store_lpv(
            static_cast<std::size_t>(expanded_next_no + 1) *
                (initial.store_lpv_upper_bound + 1), 1.0);
        const auto outer_events = scan_rdp_redo_triplets(
            outer_event_state, expanded_distance, expanded_tree, outer_redo,
            fss_rdp, expanded_store_lpv, initial.store_lpv_upper_bound, 125,
            half_window, full_window, xover_settings, probability_settings,
            probability_estimate, factorials.three_way,
            factorials.factorial, xover_api(), 0, nullptr, nullptr, 1,
            &expanded_missing, true, &shared_xdiffpos0);
        if (std::getenv("RDP_TRACE_OUTER_BROAD") != nullptr &&
            round_number <= 3) {
            for (std::size_t row = 0; row < outer_events.xover_list.size(); ++row) {
                for (const auto& event : outer_events.xover_list[row]) {
                    if (event.beginning >= 1400 && event.beginning <= 2300 &&
                        event.ending >= 3500 && event.ending <= 4600) {
                        std::cerr << "outer-broad round=" << round_number
                                  << " row=" << row << " method="
                                  << static_cast<int>(event.program_flag)
                                  << " roles=" << event.daughter << ':'
                                  << event.minor_parent << ':' << event.major_parent
                                  << " interval=" << event.beginning << '-'
                                  << event.ending << " p=" << event.probability
                                  << '\n';
                    }
                }
            }
        }
        auto dropped = drop_rdp_unused_fragment_events(
            permanent_next_no, scan_state.next_no, expanded_next_no,
            scan_state.sequence_length, minimum_sequence_size,
            expanded_trace_sub,
            next_actual_sizes, expanded_scan_state.sequence_data,
            expanded_missing, combined_events, outer_events);
        if (analysis_trace != nullptr) {
            analysis_trace->events_before_drop.push_back(
                dropped.events_before_drop);
            analysis_trace->post_rescan_events.push_back(dropped.events);
            analysis_trace->retained_actual_sequence_sizes.push_back(
                dropped.actual_sequence_sizes);
            analysis_trace->fragment_reference_counts.push_back(
                dropped.reference_counts);
            analysis_trace->fragment_reference_counts_before_drop.push_back(
                dropped.reference_counts_before_drop);
            analysis_trace->actual_sequence_sizes_before_drop.push_back(
                dropped.actual_sequence_sizes_before_drop);
        }

        const int previous_store_sequence_count = static_cast<int>(
            initial.store_lpv.size() /
            static_cast<std::size_t>(initial.store_lpv_upper_bound + 1));
        if (dropped.next_no + 1 > previous_store_sequence_count) {
            initial.store_lpv.resize(
                static_cast<std::size_t>(dropped.next_no + 1) *
                    (initial.store_lpv_upper_bound + 1),
                0.0);
            // ResortCurrentXOverIII grows StoreLPV on first use of a new
            // fragment row. In an RDP-only run that first program is RDP;
            // VB initializes only that program's newly appended cells to 1,
            // leaving disabled-program cells at ReDim Preserve's zero.
            for (int sequence = previous_store_sequence_count;
                 sequence <= dropped.next_no; ++sequence) {
                initial.store_lpv[static_cast<std::size_t>(sequence) *
                    (initial.store_lpv_upper_bound + 1)] = 1.0;
            }
        }
        scan_state = rebuild_rdp_scan_state(
            dropped.next_no, scan_state.sequence_length,
            dropped.sequence_data,
            preprocessing_api());
        missing_data = std::move(dropped.missing_data);
        events = std::move(dropped.events);
        trace_legacy_method_state("post-round", events);
        trace_sub = std::move(dropped.trace_sub);
        actual_sequence_sizes = std::move(dropped.actual_sequence_sizes);
        done_sequence = std::move(adjusted.done_sequence);
        done_row_upper_bound = adjusted.done_row_upper_bound;
        distance_state = build_rdp_distance_state(
            scan_state, 1, scan_state.sequence_length);
        tree_state = build_rdp_upgma_tree_state(
            scan_state.next_no, distance_state);
    }
    throw std::runtime_error("RDP event loop exceeded 1000 finalized events");
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

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_text(
    const std::string& fasta_text, const RdpInitialAnalysisOptions& options) {
    return run_rdp_full_analysis(
        build_rdp_scan_state_from_fasta_text(fasta_text, preprocessing_api()),
        options, nullptr);
}

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_file(
    const std::string& fasta_path, const RdpInitialAnalysisOptions& options) {
    int table_bound = 45;
    auto table = load_three_seq_table_file(fasta_path, table_bound);
    const auto* prior_table = configured_three_seq_table;
    const int prior_bound = configured_three_seq_table_bound;
    if (!table.empty()) {
        configured_three_seq_table = &table;
        configured_three_seq_table_bound = table_bound;
    }
    auto result = run_rdp_full_analysis(
        build_rdp_scan_state_from_fasta(fasta_path, preprocessing_api()),
        options, nullptr);
    configured_three_seq_table = prior_table;
    configured_three_seq_table_bound = prior_bound;
    return result;
}

RdpFullAnalysisResult run_rdp_full_analysis_from_fasta_file_with_trace(
    const std::string& fasta_path, RdpFullAnalysisTrace& trace,
    const RdpInitialAnalysisOptions& options) {
    trace = {};
    int table_bound = 45;
    auto table = load_three_seq_table_file(fasta_path, table_bound);
    const auto* prior_table = configured_three_seq_table;
    const int prior_bound = configured_three_seq_table_bound;
    if (!table.empty()) {
        configured_three_seq_table = &table;
        configured_three_seq_table_bound = table_bound;
    }
    auto result = run_rdp_full_analysis(
        build_rdp_scan_state_from_fasta(fasta_path, preprocessing_api()),
        options, &trace);
    configured_three_seq_table = prior_table;
    configured_three_seq_table_bound = prior_bound;
    return result;
}

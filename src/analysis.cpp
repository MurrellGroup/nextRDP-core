#include "analysis.hpp"

#include "MathFuncsDll.h"
#include "distance_state.hpp"
#include "legacy_method_state.hpp"
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
#include <stdexcept>
#include <string>
#include <vector>

namespace {

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

    const int permanent_next_no = initial.alignment.next_no;
    const auto permanent_analysis_scan = initial.alignment;
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
        auto selection = select_rdp_best_event(
            events, scan_state.next_no, options.p_value_cutoff,
            done_sequence, done_row_upper_bound);
        if (!selection.found) return output;
        if (selection.trace[0] < 0 ||
            selection.trace[0] >= static_cast<int>(events.xover_list.size()) ||
            selection.trace[1] < 1 ||
            selection.trace[1] > static_cast<int>(
                events.xover_list[selection.trace[0]].size())) {
            throw std::runtime_error("RDP selected an invalid event slot");
        }
        auto& selected_slot = events.xover_list[
            selection.trace[0]][selection.trace[1] - 1];
        if (selected_slot.program_flag != 0 ||
            selected_slot.daughter == selected_slot.minor_parent ||
            selected_slot.daughter == selected_slot.major_parent ||
            selected_slot.minor_parent == selected_slot.major_parent) {
            selected_slot.probability = 1.0;
            done_sequence = std::move(selection.done_sequence);
            done_row_upper_bound = scan_state.next_no;
            continue;
        }
        const RdpRawEvent selected = selected_slot;
        if (analysis_trace != nullptr) {
            analysis_trace->selected_slots.push_back({
                selection.trace[0], selection.trace[1]});
            analysis_trace->done_before_selection.push_back(done_sequence);
            analysis_trace->done_row_upper_bounds.push_back(
                done_row_upper_bound);
        }
        auto round = identify_rdp_complete_round(
            scan_state, distance_state, events, selected, missing_data,
            permanent_next_no, 20);
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
        final_event.winning_role = winner;
        final_event.probability = selected.probability;
        final_event.beginning = selected.beginning;
        final_event.ending = selected.ending;
        final_event.consensus = round.consensus.consensus;
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
        output.events.push_back(std::move(final_event));

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
        RdpLegacyEventAllocator legacy_events(temporary_events, 1);
        const auto& three_seq_table = active_three_seq_table();
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
            trace_call("geneconv", 0);
            run_rdp_geneconv_recheck(
                scan_state, triplet, initial.store_lpv,
                initial.store_lpv_upper_bound, probability_settings,
                legacy_events);
            trace_legacy_method_state("after-geneconv", temporary_events);
            trace_call("maxchi", 0);
            run_rdp_maxchi_recheck(
                scan_state, triplet, initial.store_lpv,
                initial.store_lpv_upper_bound, probability_settings,
                legacy_events, selected.beginning, selected.ending);
            trace_legacy_method_state("after-maxchi", temporary_events);
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
            rotated = triplet;
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
        append_rdp_events(events_after_final_trim, temporary_events);
        // FinalTrim's rescans above use the original selected probability.
        // MakeCollecteventsX also receives that original PXOList record.  The
        // source marks the persistent selected signal as consumed only after
        // collection, before MakeBestXOList/AddjustCXO.
        const std::array<int, 2> trace{
            selection.trace[0], selection.trace[1]};
        auto collection_events = prepare_rdp_collection_event_list(
            scan_state.next_no, winner, round.prefix.sequences, trace,
            round.final_candidates.candidate_last,
            round.final_candidates.candidate_list,
            round.final_candidates.acceptable_sequences,
            events_after_final_trim);
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
            actual_sequence_sizes, permanent_next_no, 20, propagated_pairs);
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
        const auto outer_triplets = make_rdp_outer_scan_triplets(
            permanent_analysis_scan, initially_screened_triplets,
            expanded_next_no,
            scan_state.next_no, winner,
            round.final_candidates.candidate_last, expanded_trace_sub,
            next_actual_sizes, permanent_next_no, 20, propagated_pairs,
            expanded_next_no, expanded_distance.valid_sites);
        if (analysis_trace != nullptr) {
            analysis_trace->inner_triplets.push_back(inner_triplets);
            analysis_trace->outer_triplets.push_back(outer_triplets);
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
        const auto inner_events = scan_rdp_redo_triplets(
            inner_event_state, distance_state, tree_state, inner_redo,
            fss_rdp, initial.store_lpv, initial.store_lpv_upper_bound, 125,
            half_window, full_window, xover_settings, probability_settings,
            probability_estimate, factorials.three_way,
            factorials.factorial, xover_api(), 0, nullptr, nullptr, 1,
            &mutation.missing_data, true, &shared_xdiffpos0);
        auto combined_events = adjusted.events;
        append_rdp_events(combined_events, inner_events);
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
        auto dropped = drop_rdp_unused_fragment_events(
            permanent_next_no, scan_state.next_no, expanded_next_no,
            scan_state.sequence_length, 20, expanded_trace_sub,
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

#include "MathFuncsDll.h"
#include "alist_rdp4_fixture.hpp"
#include "distance_fixture.hpp"
#include "distance_state.hpp"
#include "event_state_fixture.hpp"
#include "find_subseq_pb3_fixture.hpp"
#include "identification_fixture.hpp"
#include "identification_state.hpp"
#include "preprocess_fixture.hpp"
#include "rdp_walk_fixture.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

int self_test() {
    constexpr int next_no = 4;
    std::array<short, next_no + 1> mask{};
    std::array<short, 3 * 10> analysis_list{};

    const int last_triplet = MathFuncs::MyMathFuncs::MakeAListP2(
        1.0F,
        next_no,
        mask.data(),
        2,
        analysis_list.data());

    if (last_triplet != 9) {
        std::cerr << "MakeAListP2 self-test failed: expected last index 9, got "
                  << last_triplet << '\n';
        return 1;
    }

    std::cout << "MakeAListP2: ok (10 triplets)\n";
    return 0;
}

int preprocess_fixture() {
    const PreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    return run_preprocess_fixture(api, std::cout);
}

int distance_fixture() {
    return run_distance_fixture(&MathFuncs::MyMathFuncs::SuperDistP, std::cout);
}

int fasta_preprocess_fixture(
    const std::string& fasta_path, const std::string& fixture_path) {
    const Dna5ScanPreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    const auto actual = build_rdp_scan_state_from_fasta(fasta_path, api);
    const auto fixture = load_alist_rdp4_fixture(fixture_path);
    const auto expected_sequence = alist_rdp4_typed_section<short>(
        fixture, AlistRdp4Section::sequence_data_in);
    const auto expected_analysis = alist_rdp4_typed_section<short>(
        fixture, AlistRdp4Section::analysis_list_in);
    const auto expected_compressed =
        alist_rdp4_typed_section<unsigned char>(
            fixture, AlistRdp4Section::compressed_sequence_in);
    const bool matches = actual.next_no == fixture.header.next_no &&
        actual.sequence_length == fixture.header.sequence_length &&
        actual.analysis_list_last == fixture.header.list_length &&
        actual.compressed_sequence_ub ==
            fixture.header.compressed_sequence_ub &&
        actual.sequence_data == expected_sequence &&
        actual.analysis_list == expected_analysis &&
        actual.compressed_sequence == expected_compressed;
    if (!matches) {
        const auto first_difference = [](const auto& first, const auto& second) {
            const auto common = std::min(first.size(), second.size());
            std::size_t index = 0;
            while (index < common && first[index] == second[index]) ++index;
            return index;
        };
        std::cerr << "FASTA preprocessing parity: FAIL"
                  << " nextNo=" << (actual.next_no == fixture.header.next_no)
                  << " length="
                  << (actual.sequence_length == fixture.header.sequence_length)
                  << " listLast="
                  << (actual.analysis_list_last == fixture.header.list_length)
                  << " compressedUb="
                  << (actual.compressed_sequence_ub ==
                      fixture.header.compressed_sequence_ub)
                  << " sequence="
                  << (actual.sequence_data == expected_sequence)
                  << " sequenceFirst="
                  << first_difference(actual.sequence_data, expected_sequence)
                  << " analysis="
                  << (actual.analysis_list == expected_analysis)
                  << " compressed="
                  << (actual.compressed_sequence == expected_compressed)
                  << " compressedFirst="
                  << first_difference(
                         actual.compressed_sequence, expected_compressed)
                  << '\n';
        return 1;
    }
    std::cout << "FASTA preprocessing parity: PASS ("
              << (actual.next_no + 1) << " sequences, "
              << (actual.analysis_list_last + 1) << " triplets)\n";
    return 0;
}

int fasta_distance_fixture(
    const std::string& fasta_path, const std::string& fixture_path) {
    const Dna5ScanPreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    const auto scan_state = build_rdp_scan_state_from_fasta(fasta_path, api);
    const auto actual = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    const auto fixture = load_alist_rdp4_fixture(fixture_path);
    const auto expected = alist_rdp4_typed_section<float>(
        fixture, AlistRdp4Section::distance_in);
    if (actual.distance != expected) {
        const auto common = std::min(actual.distance.size(), expected.size());
        std::size_t first = 0;
        while (first < common && actual.distance[first] == expected[first]) {
            ++first;
        }
        std::size_t mismatch_count =
            actual.distance.size() > common ? actual.distance.size() - common :
            expected.size() - common;
        std::size_t first_off_diagonal = common;
        const auto matrix_stride =
            static_cast<std::size_t>(scan_state.next_no + 1);
        for (std::size_t index = 0; index < common; ++index) {
            if (actual.distance[index] != expected[index]) {
                ++mismatch_count;
                if (index / matrix_stride != index % matrix_stride &&
                    first_off_diagonal == common) {
                    first_off_diagonal = index;
                }
            }
        }
        std::cerr << "FASTA distance parity: FAIL (actual="
                  << actual.distance.size() << ", expected=" << expected.size()
                  << ", mismatches=" << mismatch_count
                  << ", first difference=" << first
                  << ", first off-diagonal difference="
                  << first_off_diagonal;
        if (first < common) {
            std::cerr << ", actual value=" << actual.distance[first]
                      << ", expected value=" << expected[first];
        }
        std::cerr << ")\n";
        return 1;
    }
    std::cout << "FASTA distance parity: PASS ("
              << (scan_state.next_no + 1) << " sequences)\n";
    return 0;
}

int fasta_tree_distance_fixture(
    const std::string& fasta_path, const std::string& fixture_path) {
    const Dna5ScanPreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    const auto scan_state = build_rdp_scan_state_from_fasta(fasta_path, api);
    const auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    const auto actual =
        build_rdp_upgma_tree_state(scan_state.next_no, distance_state);
    const auto fixture = load_alist_rdp4_fixture(fixture_path);
    const auto expected = alist_rdp4_typed_section<float>(
        fixture, AlistRdp4Section::tree_distance_in);
    if (actual.tree_distance != expected) {
        const auto common =
            std::min(actual.tree_distance.size(), expected.size());
        std::size_t first = 0;
        std::size_t mismatches =
            actual.tree_distance.size() > common ?
            actual.tree_distance.size() - common : expected.size() - common;
        while (first < common && actual.tree_distance[first] == expected[first]) {
            ++first;
        }
        for (std::size_t index = 0; index < common; ++index) {
            if (actual.tree_distance[index] != expected[index]) ++mismatches;
        }
        std::cerr << "FASTA tree-distance parity: FAIL (mismatches="
                  << mismatches << ", first difference=" << first;
        if (first < common) {
            std::cerr << ", actual value=" << actual.tree_distance[first]
                      << ", expected value=" << expected[first];
        }
        std::cerr << ")\n";
        return 1;
    }
    std::cout << "FASTA tree-distance parity: PASS ("
              << (scan_state.next_no + 1) << " sequences)\n";
    return 0;
}

int fasta_alist_rdp4_fixture(
    const std::string& fasta_path, const std::string& fixture_path) {
    const Dna5ScanPreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    auto scan_state = build_rdp_scan_state_from_fasta(fasta_path, api);
    auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    auto tree_state =
        build_rdp_upgma_tree_state(scan_state.next_no, distance_state);
    const auto fixture = load_alist_rdp4_fixture(fixture_path);
    const auto& h = fixture.header;
    if (scan_state.next_no != h.next_no ||
        scan_state.sequence_length != h.sequence_length ||
        scan_state.analysis_list_last != h.list_length ||
        scan_state.compressed_sequence_ub != h.compressed_sequence_ub) {
        std::cerr << "FASTA AlistRDP4 parity: FAIL (state dimensions)\n";
        return 1;
    }

    auto store_lpv = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::store_lpv_in);
    auto redo_list = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::redo_list_in);
    auto fss_rdp = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::fss_rdp_in);
    auto probability_estimate = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::probability_estimate_in);
    auto fact_three = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::fact_three_in);
    auto fact = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::fact_in);
    const auto expected_redo = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::redo_list_out);
    const auto expected_store = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::store_lpv_out);
    const auto expected_result = alist_rdp4_typed_section<int>(
        fixture, AlistRdp4Section::result_out);

    const int result = MathFuncs::MyMathFuncs::AlistRDP4(
        h.store_lpv_ub, store_lpv.data(), scan_state.analysis_list.data(),
        h.list_length, h.start, h.end, h.next_no, h.sub_threshold,
        redo_list.data(), h.circular, h.mc_correction, h.mc_flag,
        h.lowest_probability, h.target_x, h.sequence_length, h.short_output,
        h.distance_ub, distance_state.distance.data(), h.tree_distance_ub,
        tree_state.tree_distance.data(), h.fss_rdp_ub,
        h.compressed_sequence_ub, scan_state.compressed_sequence.data(),
        scan_state.sequence_data.data(), h.xover_window, h.xover_window_x,
        fss_rdp.data(), h.probability_file_flag, h.probability_one_ub,
        h.probability_two_ub, probability_estimate.data(), h.fact_three_ub,
        fact_three.data(), fact.data());
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && redo_list == expected_redo &&
        store_lpv == expected_store;
    if (!matches) {
        std::cerr << "FASTA AlistRDP4 parity: FAIL (result=" << result
                  << ", redo=" << (redo_list == expected_redo)
                  << ", store=" << (store_lpv == expected_store) << ")\n";
        return 1;
    }
    std::cout << "FASTA AlistRDP4 parity: PASS ("
              << (scan_state.analysis_list_last + 1) << " triplets)\n";
    return 0;
}

int fasta_first_xover_fixture(
    const std::string& fasta_path, const std::string& alist_fixture_path,
    const std::string& find_subseq_fixture_path) {
    const Dna5ScanPreprocessApi api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    auto scan_state = build_rdp_scan_state_from_fasta(fasta_path, api);
    const auto alist_fixture = load_alist_rdp4_fixture(alist_fixture_path);
    const auto redo = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::redo_list_out);
    int first_redo = -1;
    for (int triplet = alist_fixture.header.start;
         triplet <= alist_fixture.header.end; ++triplet) {
        if (redo[triplet] == 1) {
            first_redo = triplet;
            break;
        }
    }
    if (first_redo < 0) {
        std::cerr << "FASTA first-XOver parity: FAIL (no redo triplet)\n";
        return 1;
    }

    const std::array<int, 3> sequences{
        scan_state.analysis_list[0 + first_redo * 3],
        scan_state.analysis_list[1 + first_redo * 3],
        scan_state.analysis_list[2 + first_redo * 3],
    };
    const auto fixture =
        load_find_subseq_pb3_fixture(find_subseq_fixture_path);
    const auto& h = fixture.header;
    if (sequences[0] != h.seq1 || sequences[1] != h.seq2 ||
        sequences[2] != h.seq3 || h.next_no != scan_state.next_no ||
        h.sequence_length != scan_state.sequence_length ||
        h.compressed_sequence_ub != scan_state.compressed_sequence_ub) {
        std::cerr << "FASTA first-XOver parity: FAIL (selected triplet/state)\n";
        return 1;
    }
    std::array<int, 3> ah{};
    std::vector<char> xover_sequence(
        static_cast<std::size_t>(h.xover_sequence_ub + 1) * 3, 0);
    auto fss_rdp = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::fss_rdp_in);
    const auto expected_ah = find_subseq_pb3_section<int>(fixture, 101);
    const auto expected_xover = find_subseq_pb3_section<char>(fixture, 102);
    const auto expected_result = find_subseq_pb3_section<int>(fixture, 103);
    const int result = MathFuncs::MyMathFuncs::FindSubSeqPB3(
        ah.data(), h.fss_ub, h.xover_window, h.compressed_sequence_ub,
        h.sequence_length, h.next_no, h.seq1, h.seq2, h.seq3,
        scan_state.compressed_sequence.data(), h.xover_sequence_ub,
        xover_sequence.data(), fss_rdp.data());
    const std::vector<int> actual_ah(ah.begin(), ah.end());
    if (expected_result.size() != 1 || result != expected_result[0] ||
        actual_ah != expected_ah || xover_sequence != expected_xover) {
        std::cerr << "FASTA first-XOver parity: FAIL (result=" << result
                  << ", AH=" << (actual_ah == expected_ah)
                  << ", window=" << (xover_sequence == expected_xover)
                  << ", actual AH=" << actual_ah[0] << ',' << actual_ah[1]
                  << ',' << actual_ah[2] << ", expected AH="
                  << expected_ah[0] << ',' << expected_ah[1] << ','
                  << expected_ah[2] << ")\n";
        return 1;
    }
    std::cout << "FASTA first-XOver parity: PASS (triplet "
              << sequences[0] << ',' << sequences[1] << ',' << sequences[2]
              << ", informative length " << result << ")\n";
    return 0;
}

int fasta_first_xover_walk_fixture(
    const std::string& fasta_path, const std::string& alist_fixture_path,
    const std::string& find_subseq_fixture_path,
    const std::string& homology_fixture_path,
    const std::string& find_next_fixture_path,
    const std::string& define_event_fixture_path,
    const std::string& probability_p2_fixture_path,
    const std::string& probability_p_fixture_path,
    const std::string& pb4_fixture_path,
    const std::string& clean_fixture_path) {
    const Dna5ScanPreprocessApi preprocess_api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    const Dna5XoverApi xover_api{
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
    auto scan_state =
        build_rdp_scan_state_from_fasta(fasta_path, preprocess_api);
    auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    auto tree_state =
        build_rdp_upgma_tree_state(scan_state.next_no, distance_state);
    const auto alist_fixture = load_alist_rdp4_fixture(alist_fixture_path);
    const auto& alist_header = alist_fixture.header;
    auto store_lpv = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::store_lpv_in);
    auto redo = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::redo_list_in);
    auto fss_rdp = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::fss_rdp_in);
    auto probability_estimate = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::probability_estimate_in);
    auto fact_three = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::fact_three_in);
    auto fact = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::fact_in);
    const auto expected_redo = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::redo_list_out);
    const auto expected_store = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::store_lpv_out);
    const auto expected_alist_result = alist_rdp4_typed_section<int>(
        alist_fixture, AlistRdp4Section::result_out);
    const int alist_result = MathFuncs::MyMathFuncs::AlistRDP4(
        alist_header.store_lpv_ub, store_lpv.data(),
        scan_state.analysis_list.data(), alist_header.list_length,
        alist_header.start, alist_header.end, alist_header.next_no,
        alist_header.sub_threshold, redo.data(), alist_header.circular,
        alist_header.mc_correction, alist_header.mc_flag,
        alist_header.lowest_probability, alist_header.target_x,
        alist_header.sequence_length, alist_header.short_output,
        alist_header.distance_ub, distance_state.distance.data(),
        alist_header.tree_distance_ub, tree_state.tree_distance.data(),
        alist_header.fss_rdp_ub, alist_header.compressed_sequence_ub,
        scan_state.compressed_sequence.data(), scan_state.sequence_data.data(),
        alist_header.xover_window, alist_header.xover_window_x,
        fss_rdp.data(), alist_header.probability_file_flag,
        alist_header.probability_one_ub, alist_header.probability_two_ub,
        probability_estimate.data(), alist_header.fact_three_ub,
        fact_three.data(), fact.data());
    if (expected_alist_result.size() != 1 ||
        alist_result != expected_alist_result[0] || redo != expected_redo ||
        store_lpv != expected_store) {
        std::cerr << "FASTA first-XOver walk parity: FAIL AlistRDP4\n";
        return 1;
    }
    int first_redo = -1;
    for (int triplet = alist_fixture.header.start;
         triplet <= alist_fixture.header.end; ++triplet) {
        if (redo[triplet] == 1) {
            first_redo = triplet;
            break;
        }
    }
    const auto pb3_fixture =
        load_find_subseq_pb3_fixture(find_subseq_fixture_path);
    auto state = build_rdp_first_xover_state(
        scan_state, distance_state, tree_state, first_redo,
        pb3_fixture.header.fss_ub, fss_rdp,
        pb3_fixture.header.xover_window,
        alist_fixture.header.xover_window_x, xover_api);
    const auto expected_ah =
        find_subseq_pb3_section<int>(pb3_fixture, 101);
    const auto expected_xover =
        find_subseq_pb3_section<char>(pb3_fixture, 102);
    const auto expected_pb3_result =
        find_subseq_pb3_section<int>(pb3_fixture, 103);

    const char xoh_magic[8] = {'X', 'O', 'H', 'O', 'M', 'P', '\0', '\0'};
    const auto xoh_fixture =
        load_rdp_sectioned_fixture<XOHomologyPCaptureHeader>(
            homology_fixture_path, xoh_magic);
    const auto expected_xoh_sequence =
        rdp_fixture_section<char>(xoh_fixture, 1);
    const auto expected_homology =
        rdp_fixture_section<int>(xoh_fixture, 101);
    const auto expected_xoh_result =
        rdp_fixture_section<int>(xoh_fixture, 102);
    const char find_next_magic[8] = {
        'F', 'I', 'N', 'D', 'N', 'X', 'T', '\0'};
    const auto find_next_fixture =
        load_rdp_sectioned_fixture<FindNextPCaptureHeader>(
            find_next_fixture_path, find_next_magic);
    const auto expected_find_next_homology =
        rdp_fixture_section<int>(find_next_fixture, 1);
    const auto expected_next_position =
        rdp_fixture_section<int>(find_next_fixture, 101);
    const char define_event_magic[8] = {
        'D', 'E', 'F', 'E', 'V', 'P', '2', '\0'};
    const auto define_event_fixture =
        load_rdp_sectioned_fixture<DefineEventP2CaptureHeader>(
            define_event_fixture_path, define_event_magic);
    const auto define_scalars_in =
        rdp_fixture_section<int>(define_event_fixture, 1);
    const auto define_xover_in =
        rdp_fixture_section<char>(define_event_fixture, 2);
    const auto define_homology_in =
        rdp_fixture_section<int>(define_event_fixture, 3);
    const auto expected_define_scalars =
        rdp_fixture_section<int>(define_event_fixture, 101);
    const auto expected_define_result =
        rdp_fixture_section<int>(define_event_fixture, 102);
    const RdpXoverSettings settings{
        define_event_fixture.header.short_output,
        define_event_fixture.header.long_winded,
        define_event_fixture.header.target_x,
        define_event_fixture.header.circular,
    };
    const auto homology_before_event = state.homology;
    const auto xover_before_event = state.xover_sequence;
    define_rdp_first_xover_event(
        state, scan_state, pb3_fixture.header.xover_window, settings,
        xover_api);
    const std::vector<int> actual_define_scalars{
        state.end_flag, state.event_begin, state.event_end,
        state.number_in_common, state.event_length};
    const int first_define_result = state.event_position;
    const bool define_buffers_match =
        state.homology_at_define == define_homology_in &&
        state.xover_sequence_at_define == define_xover_in;
    const bool define_header_matches =
        state.homology_ub == define_event_fixture.header.homology_ub &&
        state.med_homology == define_event_fixture.header.med &&
        state.high_homology == define_event_fixture.header.high &&
        state.low_homology == define_event_fixture.header.low &&
        alist_fixture.header.short_output ==
            define_event_fixture.header.short_output &&
        alist_fixture.header.target_x == define_event_fixture.header.target_x &&
        alist_fixture.header.circular == define_event_fixture.header.circular &&
        state.define_input_position == define_event_fixture.header.xx &&
        pb3_fixture.header.xover_window ==
            define_event_fixture.header.xover_window &&
        scan_state.sequence_length ==
            define_event_fixture.header.sequence_length &&
        state.homology_length == define_event_fixture.header.xover_length &&
        state.sequence_minor ==
            define_event_fixture.header.sequence_daughter &&
        state.sequence_daughter == define_event_fixture.header.sequence_minor;
    const char probability_magic[8] = {
        'P', 'R', 'O', 'B', 'C', 'P', '2', '\0'};
    const auto probability_fixture =
        load_rdp_sectioned_fixture<ProbCalcP2CaptureHeader>(
            probability_p2_fixture_path, probability_magic);
    const auto expected_probability_fact =
        rdp_fixture_section<double>(probability_fixture, 1);
    const auto expected_probability =
        rdp_fixture_section<double>(probability_fixture, 101);
    const char probability_p_magic[8] = {
        'P', 'R', 'O', 'B', 'C', 'P', '\0', '\0'};
    const auto probability_p_fixture =
        load_rdp_sectioned_fixture<ProbCalcPCaptureHeader>(
            probability_p_fixture_path, probability_p_magic);
    const auto expected_probability_fact_p =
        rdp_fixture_section<double>(probability_p_fixture, 1);
    const auto expected_probability_p =
        rdp_fixture_section<double>(probability_p_fixture, 101);
    const RdpProbabilitySettings probability_settings{
        alist_header.circular,
        alist_header.mc_correction,
        alist_header.mc_flag,
        alist_header.probability_file_flag,
        alist_header.probability_one_ub,
        alist_header.probability_two_ub,
        alist_header.fact_three_ub,
        alist_header.lowest_probability,
    };
    calculate_rdp_first_xover_probability(
        state, probability_settings, probability_estimate, fact_three, fact,
        xover_api);
    if (!state.probability_tested) {
        continue_rdp_xover_to_first_probability(
            state, scan_state, pb3_fixture.header.xover_window, settings,
            probability_settings, probability_estimate, fact_three, fact,
            xover_api);
    }
    const auto probability_close = [](const double actual,
                                      const std::vector<double>& expected) {
        if (expected.size() != 1) return false;
        const double scale = std::max(1.0, std::abs(expected[0]));
        return std::abs(actual - expected[0]) <= 1e-14 * scale;
    };
    const bool probability_p2_matches = state.probability_tested &&
        state.used_probability_p2 && fact_three == expected_probability_fact &&
        state.probability_length == probability_fixture.header.xover_length &&
        state.probability_same == probability_fixture.header.number_in_common &&
        state.individual_probability ==
            probability_fixture.header.individual_probability &&
        state.homology_length == probability_fixture.header.informative_length &&
        probability_close(state.event_probability, expected_probability);
    const bool probability_p_matches = state.probability_tested &&
        !state.used_probability_p2 && fact == expected_probability_fact_p &&
        state.probability_length == probability_p_fixture.header.xover_length &&
        state.probability_same ==
            probability_p_fixture.header.number_in_common &&
        state.individual_probability ==
            probability_p_fixture.header.individual_probability &&
        state.homology_length ==
            probability_p_fixture.header.informative_length &&
        probability_close(state.event_probability, expected_probability_p);
    const bool probability_matches =
        probability_p2_matches || probability_p_matches;
    const std::vector<int> actual_ah(
        state.agreement_counts.begin(), state.agreement_counts.end());
    const bool matches = expected_pb3_result.size() == 1 &&
        state.informative_length == expected_pb3_result[0] &&
        actual_ah == expected_ah && xover_before_event == expected_xover &&
        state.initial_high_homology == xoh_fixture.header.inlyer &&
        xover_before_event == expected_xoh_sequence &&
        state.homology_sequence_length ==
            xoh_fixture.header.sequence_length &&
        state.homology_length == xoh_fixture.header.xover_length &&
        expected_xoh_result.size() == 1 &&
        state.homology_start == expected_xoh_result[0] &&
        homology_before_event == expected_homology &&
        homology_before_event == expected_find_next_homology &&
        state.homology_ub == find_next_fixture.header.homology_ub &&
        find_next_fixture.header.start == 1 &&
        state.high_homology == find_next_fixture.header.high &&
        state.med_homology == find_next_fixture.header.med &&
        state.low_homology == find_next_fixture.header.low &&
        state.homology_length == find_next_fixture.header.xover_length &&
        pb3_fixture.header.xover_window ==
            find_next_fixture.header.xover_window &&
        expected_next_position.size() == 1 &&
        state.next_position == expected_next_position[0] &&
        define_scalars_in == std::vector<int>(5, 0) &&
        define_buffers_match && define_header_matches &&
        expected_define_result.size() == 1 &&
        first_define_result == expected_define_result[0] &&
        actual_define_scalars == expected_define_scalars &&
        probability_matches;
    if (!matches) {
        const auto homology_difference = std::mismatch(
            expected_homology.begin(), expected_homology.end(),
            homology_before_event.begin(), homology_before_event.end());
        const auto sequence_difference = std::mismatch(
            expected_xoh_sequence.begin(), expected_xoh_sequence.end(),
            xover_before_event.begin(), xover_before_event.end());
        std::cerr << "FASTA first-XOver walk parity: FAIL"
                  << " PB3="
                  << (expected_pb3_result.size() == 1 &&
                      state.informative_length == expected_pb3_result[0] &&
                      actual_ah == expected_ah &&
                      xover_before_event == expected_xover)
                  << " role="
                  << (state.initial_high_homology ==
                      xoh_fixture.header.inlyer)
                  << " XOH="
                  << (expected_xoh_result.size() == 1 &&
                      state.homology_start == expected_xoh_result[0] &&
                      homology_before_event == expected_homology)
                  << " result=" << state.homology_start << '/'
                  << (expected_xoh_result.empty() ? -1 : expected_xoh_result[0])
                  << " homology-size=" << homology_before_event.size() << '/'
                  << expected_homology.size()
                  << " homology-diff="
                  << (homology_difference.first == expected_homology.end()
                          ? -1
                          : static_cast<long long>(std::distance(
                                expected_homology.begin(),
                                homology_difference.first)))
                  << " xoh-input-size=" << xover_before_event.size() << '/'
                  << expected_xoh_sequence.size()
                  << " xoh-input-diff="
                  << (sequence_difference.first == expected_xoh_sequence.end()
                          ? -1
                          : static_cast<long long>(std::distance(
                                expected_xoh_sequence.begin(),
                                sequence_difference.first)))
                  << " roles=" << state.high_homology << ','
                  << state.med_homology << ',' << state.low_homology << '/'
                  << find_next_fixture.header.high << ','
                  << find_next_fixture.header.med << ','
                  << find_next_fixture.header.low
                  << " next=" << state.next_position << '/'
                  << (expected_next_position.empty()
                          ? -1
                          : expected_next_position[0])
                  << " define=" << first_define_result << '/'
                  << (expected_define_result.empty()
                          ? -1
                          : expected_define_result[0])
                  << " define-scalars="
                  << (actual_define_scalars == expected_define_scalars)
                  << " define-input-scalars="
                  << (define_scalars_in == std::vector<int>(5, 0))
                  << " define-buffers=" << define_buffers_match
                  << " define-header=" << define_header_matches
                  << " passed-daughter/minor=" << state.sequence_minor << ','
                  << state.sequence_daughter << '/'
                  << define_event_fixture.header.sequence_daughter << ','
                  << define_event_fixture.header.sequence_minor
                  << " probability=" << probability_matches
                  << " probability-state=" << state.probability_tested << ','
                  << (state.used_probability_p2 ? "P2" : "P") << ','
                  << state.probability_length << ','
                  << state.probability_same << ','
                  << state.individual_probability << " expected="
                  << (state.used_probability_p2
                          ? probability_fixture.header.xover_length
                          : probability_p_fixture.header.xover_length) << ','
                  << (state.used_probability_p2
                          ? probability_fixture.header.number_in_common
                          : probability_p_fixture.header.number_in_common) << ','
                  << (state.used_probability_p2
                          ? probability_fixture.header.individual_probability
                          : probability_p_fixture.header.individual_probability)
                  << " prefilter=" << state.probability_prefilter_value << '/'
                  << alist_header.lowest_probability
                  << '\n';
        return 1;
    }
    bool reached_significant_event =
        apply_rdp_probability_cutoff(state, probability_settings);
    while (!reached_significant_event) {
        state.probability_tested = false;
        continue_rdp_xover_to_first_probability(
            state, scan_state, pb3_fixture.header.xover_window, settings,
            probability_settings, probability_estimate, fact_three, fact,
            xover_api);
        if (!state.probability_tested) {
            if (!advance_rdp_role_cycle(state, 0)) break;
            scan_rdp_current_roles_to_first_probability(
                state, scan_state, pb3_fixture.header.xover_window, settings,
                probability_settings, probability_estimate, fact_three, fact,
                xover_api);
        }
        if (!state.probability_tested) continue;
        reached_significant_event =
            apply_rdp_probability_cutoff(state, probability_settings);
    }
    if (reached_significant_event) {
        const char pb4_magic[8] = {
            'F', 'S', 'P', 'B', '4', '\0', '\0', '\0'};
        const auto pb4_fixture = load_find_subseq_fixture(
            pb4_fixture_path, pb4_magic, "FindSubSeqPB4");
        const auto expected_pb4_ah =
            find_subseq_pb3_section<int>(pb4_fixture, 101);
        const auto expected_pb4_xover =
            find_subseq_pb3_section<char>(pb4_fixture, 102);
        const auto expected_xdiffpos =
            find_subseq_pb3_section<int>(pb4_fixture, 103);
        const auto expected_xposdiff =
            find_subseq_pb3_section<int>(pb4_fixture, 104);
        const auto expected_pb4_result =
            find_subseq_pb3_section<int>(pb4_fixture, 105);
        const auto expected_pb4_ah_in =
            find_subseq_pb3_section<int>(pb4_fixture, 1);
        const auto expected_pb4_compressed =
            find_subseq_pb3_section<unsigned char>(pb4_fixture, 2);
        const auto expected_pb4_xover_in =
            find_subseq_pb3_section<char>(pb4_fixture, 3);
        const auto expected_xdiffpos_in =
            find_subseq_pb3_section<int>(pb4_fixture, 4);
        const auto expected_xposdiff_in =
            find_subseq_pb3_section<int>(pb4_fixture, 5);
        const auto expected_pb4_fss =
            find_subseq_pb3_section<unsigned char>(pb4_fixture, 6);
        const std::vector<int> pb4_ah_before(
            state.agreement_counts.begin(), state.agreement_counts.end());
        const auto pb4_xover_before = state.xover_sequence;
        const std::vector<int> zero_positions(
            static_cast<std::size_t>(scan_state.sequence_length + 201), 0);
        const bool pb4_inputs_match =
            pb4_fixture.header.fss_ub == alist_header.fss_rdp_ub &&
            pb4_fixture.header.xover_window ==
                pb3_fixture.header.xover_window &&
            pb4_fixture.header.compressed_sequence_ub ==
                scan_state.compressed_sequence_ub &&
            pb4_fixture.header.sequence_length == scan_state.sequence_length &&
            pb4_fixture.header.next_no == scan_state.next_no &&
            pb4_fixture.header.seq1 == state.sequences[0] &&
            pb4_fixture.header.seq2 == state.sequences[1] &&
            pb4_fixture.header.seq3 == state.sequences[2] &&
            pb4_fixture.header.xover_sequence_ub ==
                state.xover_sequence_ub &&
            pb4_ah_before == expected_pb4_ah_in &&
            scan_state.compressed_sequence == expected_pb4_compressed &&
            pb4_xover_before == expected_pb4_xover_in &&
            zero_positions == expected_xdiffpos_in &&
            zero_positions == expected_xposdiff_in &&
            fss_rdp == expected_pb4_fss;
        build_rdp_first_position_maps(
            state, scan_state, alist_header.fss_rdp_ub,
            pb3_fixture.header.xover_window, fss_rdp, xover_api);
        const std::vector<int> pb4_ah_after(
            state.agreement_counts.begin(), state.agreement_counts.end());
        const bool pb4_matches = pb4_inputs_match &&
            expected_pb4_result.size() == 1 &&
            state.position_map_result == expected_pb4_result[0] &&
            pb4_ah_after == expected_pb4_ah &&
            state.xover_sequence == expected_pb4_xover &&
            state.xdiffpos == expected_xdiffpos &&
            state.xposdiff == expected_xposdiff;
        if (!pb4_matches) {
            std::cerr << "FASTA first-XOver PB4 parity: FAIL inputs="
                      << pb4_inputs_match << " result="
                      << state.position_map_result << '/'
                      << (expected_pb4_result.empty()
                              ? -1
                              : expected_pb4_result[0])
                      << '\n';
            return 1;
        }
    }
    while (true) {
        state.probability_tested = false;
        continue_rdp_xover_to_first_probability(
            state, scan_state, pb3_fixture.header.xover_window, settings,
            probability_settings, probability_estimate, fact_three, fact,
            xover_api);
        if (state.probability_tested) {
            apply_rdp_probability_cutoff(state, probability_settings);
            continue;
        }
        if (!advance_rdp_role_cycle(state, 0)) break;
        scan_rdp_current_roles_to_first_probability(
            state, scan_state, pb3_fixture.header.xover_window, settings,
            probability_settings, probability_estimate, fact_three, fact,
            xover_api);
        if (state.probability_tested) {
            apply_rdp_probability_cutoff(state, probability_settings);
        }
    }
    const char clean_magic[8] = {
        'C', 'L', 'N', 'X', 'O', 'S', 'N', 'W'};
    const auto clean_fixture =
        load_rdp_sectioned_fixture<CleanXOSNWCaptureHeader>(
            clean_fixture_path, clean_magic);
    const auto expected_clean_input =
        rdp_fixture_section<char>(clean_fixture, 1);
    const auto expected_clean_output =
        rdp_fixture_section<char>(clean_fixture, 101);
    const auto expected_clean_result =
        rdp_fixture_section<int>(clean_fixture, 102);
    const bool clean_input_matches =
        clean_fixture.header.xover_length ==
            state.homology_length + pb3_fixture.header.xover_window * 2 &&
        clean_fixture.header.xover_window ==
            pb3_fixture.header.xover_window &&
        clean_fixture.header.xover_sequence_ub == state.xover_sequence_ub &&
        state.xover_sequence == expected_clean_input;
    const int clean_result = xover_api.clean_xover_sequence(
        state.homology_length + pb3_fixture.header.xover_window * 2,
        pb3_fixture.header.xover_window, state.xover_sequence_ub,
        state.xover_sequence.data());
    if (!clean_input_matches || expected_clean_result.size() != 1 ||
        clean_result != expected_clean_result[0] ||
        state.xover_sequence != expected_clean_output) {
        std::cerr << "FASTA first-XOver CleanXOSNW parity: FAIL input="
                  << clean_input_matches << " result=" << clean_result << '/'
                  << (expected_clean_result.empty()
                          ? -1
                          : expected_clean_result[0])
                  << '\n';
        return 1;
    }
    std::cout << "FASTA first-XOver walk parity: PASS (FindSubSeqPB3 -> "
                 "XOHomologyP -> role ranking -> FindNextP -> "
                 "DefineEventP2 -> ProbCalcP/P2; significant-in-role="
              << reached_significant_event
              << (reached_significant_event ? " -> FindSubSeqPB4" : "")
              << " -> CleanXOSNW)\n";
    return 0;
}

int fasta_all_redo_events_fixture(
    const std::string& fasta_path, const std::string& alist_fixture_path,
    const std::string& define_event_fixture_path,
    const std::string& make_test_fixture_path,
    const std::string& find_best_fixture_path,
    const std::string& ufdist_fixture_path,
    const std::string& super_dist_p2_fixture_path,
    const std::string& check_matrix_fixture_path,
    const std::string& make_nj_fixture_path,
    const std::string& make_sdmp_fixture_path,
    const std::array<std::string, 3>& fill_rmat_fixture_paths,
    const std::string& calcr_fixture_path,
    const std::string& make_rlist_fixture_path,
    const std::string& find_actual_events_fixture_path,
    const std::string& strip_dup_inv_fixture_path,
    const std::string& rcompat_fixture_path,
    const std::string& phpr_fixture_path,
    const std::string& score_support_fixture_path) {
    const Dna5ScanPreprocessApi preprocess_api{
        &MathFuncs::MyMathFuncs::MakeAListP2,
        &MathFuncs::MyMathFuncs::CountNucs,
        &MathFuncs::MyMathFuncs::RecodeNucs,
        &MathFuncs::MyMathFuncs::DoRecodeP,
        &MathFuncs::MyMathFuncs::MakeCompressSeqP,
    };
    const Dna5XoverApi xover_api{
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
    auto scan_state =
        build_rdp_scan_state_from_fasta(fasta_path, preprocess_api);
    auto distance_state = build_rdp_distance_state(
        scan_state, 1, scan_state.sequence_length);
    auto tree_state =
        build_rdp_upgma_tree_state(scan_state.next_no, distance_state);
    const auto alist_fixture = load_alist_rdp4_fixture(alist_fixture_path);
    const auto& h = alist_fixture.header;
    auto store_lpv = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::store_lpv_in);
    auto redo = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::redo_list_in);
    auto fss_rdp = alist_rdp4_typed_section<unsigned char>(
        alist_fixture, AlistRdp4Section::fss_rdp_in);
    auto probability_estimate = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::probability_estimate_in);
    auto fact_three = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::fact_three_in);
    auto fact = alist_rdp4_typed_section<double>(
        alist_fixture, AlistRdp4Section::fact_in);
    const int redo_count = MathFuncs::MyMathFuncs::AlistRDP4(
        h.store_lpv_ub, store_lpv.data(), scan_state.analysis_list.data(),
        h.list_length, h.start, h.end, h.next_no, h.sub_threshold,
        redo.data(), h.circular, h.mc_correction, h.mc_flag,
        h.lowest_probability, h.target_x, h.sequence_length, h.short_output,
        h.distance_ub, distance_state.distance.data(), h.tree_distance_ub,
        tree_state.tree_distance.data(), h.fss_rdp_ub,
        h.compressed_sequence_ub, scan_state.compressed_sequence.data(),
        scan_state.sequence_data.data(), h.xover_window, h.xover_window_x,
        fss_rdp.data(), h.probability_file_flag, h.probability_one_ub,
        h.probability_two_ub, probability_estimate.data(), h.fact_three_ub,
        fact_three.data(), fact.data());

    const char define_magic[8] = {
        'D', 'E', 'F', 'E', 'V', 'P', '2', '\0'};
    const auto define_fixture =
        load_rdp_sectioned_fixture<DefineEventP2CaptureHeader>(
            define_event_fixture_path, define_magic);
    const RdpXoverSettings xover_settings{
        define_fixture.header.short_output,
        define_fixture.header.long_winded,
        define_fixture.header.target_x,
        define_fixture.header.circular,
    };
    const RdpProbabilitySettings probability_settings{
        h.circular,
        h.mc_correction,
        h.mc_flag,
        h.probability_file_flag,
        h.probability_one_ub,
        h.probability_two_ub,
        h.fact_three_ub,
        h.lowest_probability,
    };
    auto events = scan_rdp_redo_triplets(
        scan_state, distance_state, tree_state, redo, fss_rdp, store_lpv,
        h.store_lpv_ub, h.fss_rdp_ub, h.xover_window, h.xover_window_x,
        xover_settings, probability_settings, probability_estimate,
        fact_three, fact, xover_api);

    const char make_test_magic[8] = {
        'M', 'K', 'T', 'E', 'S', 'T', 'P', 'V'};
    const auto make_test_fixture =
        load_rdp_sectioned_fixture<MakeTestPVsCaptureHeader>(
            make_test_fixture_path, make_test_magic);
    const auto native_current =
        rdp_fixture_section<short>(make_test_fixture, 2);
    const auto native_events =
        rdp_fixture_section<XOVERDEFINE>(make_test_fixture, 3);
    int raw_total = 0;
    int native_total = 0;
    int matching_rows = 0;
    int matching_event_identity = 0;
    int matching_event_probability = 0;
    int compared_events = 0;
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        raw_total += events.current_xover[sequence];
        native_total += native_current[sequence];
        if (events.current_xover[sequence] == native_current[sequence]) {
            ++matching_rows;
        }
        const int common = std::min<int>(
            events.current_xover[sequence], native_current[sequence]);
        for (int slot = 1; slot <= common; ++slot) {
            const auto& actual = events.xover_list[sequence][slot - 1];
            const auto native_index = static_cast<std::size_t>(sequence) +
                static_cast<std::size_t>(slot) *
                    static_cast<std::size_t>(
                        make_test_fixture.header.xover_rows_ub + 1);
            if (native_index >= native_events.size()) {
                throw std::runtime_error(
                    "native XOverList lookup exceeds its bounds");
            }
            const auto& expected = native_events[native_index];
            const bool identity_matches =
                actual.outside_flag == expected.OutsideFlag &&
                actual.misidentify_flag == expected.MissIdentifyFlag &&
                actual.program_flag == expected.ProgramFlag &&
                actual.sbp_flag == expected.SBPFlag &&
                actual.accept == expected.Accept &&
                actual.major_parent == expected.MajorP &&
                actual.minor_parent == expected.MinorP &&
                actual.daughter == expected.Daughter &&
                actual.beginning == expected.Beginning &&
                actual.ending == expected.Ending &&
                actual.length_holder == expected.LHolder &&
                actual.event_number == expected.Eventnumber &&
                actual.permutation_pvalue == expected.PermPVal &&
                actual.begin_parent == expected.BeginP &&
                actual.end_parent == expected.EndP &&
                actual.distance_holder == expected.DHolder;
            if (identity_matches) ++matching_event_identity;
            if (actual.probability == expected.Probability) {
                ++matching_event_probability;
            }
            ++compared_events;
        }
    }
    const bool raw_structure_matches = raw_total == native_total &&
        matching_rows == scan_state.next_no + 1 &&
        matching_event_identity == compared_events &&
        compared_events == raw_total;
    if (!raw_structure_matches) {
        std::cerr << "RDP all-redo raw event structural parity: FAIL\n";
        return 1;
    }

    const int row_count = make_test_fixture.header.xover_rows_ub + 1;
    const int slot_count = make_test_fixture.header.xover_slots_ub + 1;
    std::vector<XOVERDEFINE> generated_xovers(
        static_cast<std::size_t>(row_count) * slot_count);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        for (int slot = 1; slot <= events.current_xover[sequence]; ++slot) {
            const auto& source = events.xover_list[sequence][slot - 1];
            auto& destination = generated_xovers[
                static_cast<std::size_t>(sequence) +
                static_cast<std::size_t>(slot) * row_count];
            destination.OutsideFlag = source.outside_flag;
            destination.MissIdentifyFlag = source.misidentify_flag;
            destination.ProgramFlag = source.program_flag;
            destination.SBPFlag = source.sbp_flag;
            destination.Accept = source.accept;
            destination.MajorP = source.major_parent;
            destination.MinorP = source.minor_parent;
            destination.Daughter = source.daughter;
            destination.Beginning = source.beginning;
            destination.Ending = source.ending;
            destination.LHolder = source.length_holder;
            destination.Eventnumber = source.event_number;
            destination.PermPVal = source.permutation_pvalue;
            destination.BeginP = source.begin_parent;
            destination.EndP = source.end_parent;
            destination.Probability = source.probability;
            destination.DHolder = source.distance_holder;
        }
    }
    auto generated_done =
        rdp_fixture_section<unsigned char>(make_test_fixture, 1);
    auto generated_test_pvs =
        rdp_fixture_section<double>(make_test_fixture, 4);
    const int make_test_result = MathFuncs::MyMathFuncs::MakeTestPVs(
        make_test_fixture.header.done_sequence_ub, generated_done.data(),
        scan_state.next_no, make_test_fixture.header.xover_rows_ub,
        make_test_fixture.header.xover_slots_ub,
        events.current_xover.data(), generated_xovers.data(),
        generated_test_pvs.data());
    const auto native_done_out =
        rdp_fixture_section<unsigned char>(make_test_fixture, 101);
    const auto native_make_test_result =
        rdp_fixture_section<int>(make_test_fixture, 105);
    const bool make_test_structure_matches =
        native_make_test_result.size() == 1 &&
        make_test_result == native_make_test_result[0] &&
        generated_done == native_done_out;

    const char find_best_magic[8] = {
        'F', 'B', 'R', 'S', 'I', 'G', '2', '\0'};
    const auto find_best_fixture =
        load_rdp_sectioned_fixture<FindBestRecSignalP2CaptureHeader>(
            find_best_fixture_path, find_best_magic);
    auto low_probability =
        rdp_fixture_section<double>(find_best_fixture, 1);
    auto trace = rdp_fixture_section<int>(find_best_fixture, 3);
    const int find_best_result = MathFuncs::MyMathFuncs::FindBestRecSignalP2(
        find_best_fixture.header.done_target, scan_state.next_no,
        find_best_fixture.header.probability_rows_ub,
        find_best_fixture.header.probability_columns_ub,
        low_probability.data(),
        reinterpret_cast<char*>(generated_done.data()), trace.data(),
        events.current_xover.data(), generated_test_pvs.data());
    const auto native_trace =
        rdp_fixture_section<int>(find_best_fixture, 103);
    const auto native_find_best_result =
        rdp_fixture_section<int>(find_best_fixture, 106);
    const bool first_selection_matches =
        native_find_best_result.size() == 1 &&
        find_best_result == native_find_best_result[0] &&
        trace == native_trace;
    const char ufdist_magic[8] = {
        'U', 'F', 'D', 'I', 'S', 'T', '\0', '\0'};
    const auto ufdist_fixture =
        load_rdp_sectioned_fixture<UFDistCaptureHeader>(
            ufdist_fixture_path, ufdist_magic);
    const auto& selected = events.xover_list[trace[0]][trace[1] - 1];
    std::vector<int> selected_sequences{
        selected.daughter, selected.minor_parent, selected.major_parent};
    const auto native_valid =
        rdp_fixture_section<float>(ufdist_fixture, 1);
    const auto native_differences =
        rdp_fixture_section<float>(ufdist_fixture, 2);
    const auto native_sequences =
        rdp_fixture_section<int>(ufdist_fixture, 5);
    const auto native_sequence_data =
        rdp_fixture_section<short>(ufdist_fixture, 6);
    bool selected_sequence_data_matches = true;
    const int sequence_stride = scan_state.sequence_length + 1;
    for (const int sequence : selected_sequences) {
        for (int position = 0; position <= scan_state.sequence_length;
             ++position) {
            const auto index = static_cast<std::size_t>(position) +
                static_cast<std::size_t>(sequence) * sequence_stride;
            if (scan_state.sequence_data[index] != native_sequence_data[index]) {
                selected_sequence_data_matches = false;
            }
        }
    }
    const bool ufdist_inputs_match =
        ufdist_fixture.header.sequence_length == scan_state.sequence_length &&
        ufdist_fixture.header.begin == selected.beginning &&
        ufdist_fixture.header.end == selected.ending &&
        ufdist_fixture.header.pair_matrix_ub == scan_state.next_no &&
        ufdist_fixture.header.sequence_data_ub ==
            scan_state.sequence_length &&
        distance_state.valid_sites == native_valid &&
        distance_state.differences == native_differences &&
        selected_sequences == native_sequences &&
        selected_sequence_data_matches;
    std::vector<float> breakpoint_distance(3, 0.0F);
    std::vector<float> remainder_distance(3, 0.0F);
    const int ufdist_result = MathFuncs::MyMathFuncs::UFDist(
        scan_state.sequence_length, selected.beginning, selected.ending,
        scan_state.next_no, distance_state.valid_sites.data(),
        distance_state.differences.data(), breakpoint_distance.data(),
        remainder_distance.data(), selected_sequences.data(),
        scan_state.sequence_length, scan_state.sequence_data.data());
    const auto native_breakpoint =
        rdp_fixture_section<float>(ufdist_fixture, 101);
    const auto native_remainder =
        rdp_fixture_section<float>(ufdist_fixture, 102);
    const auto native_ufdist_result =
        rdp_fixture_section<int>(ufdist_fixture, 103);
    const bool ufdist_matches = ufdist_inputs_match &&
        native_ufdist_result.size() == 1 &&
        ufdist_result == native_ufdist_result[0] &&
        breakpoint_distance == native_breakpoint &&
        remainder_distance == native_remainder;
    if (!ufdist_matches) {
        std::cerr << "UFDist integration mismatch: interval="
                  << selected.beginning << '-' << selected.ending << '/'
                  << ufdist_fixture.header.begin << '-'
                  << ufdist_fixture.header.end << " roles="
                  << selected_sequences[0] << ',' << selected_sequences[1]
                  << ',' << selected_sequences[2] << '/';
        for (std::size_t index = 0; index < native_sequences.size(); ++index) {
            if (index != 0) std::cerr << ',';
            std::cerr << native_sequences[index];
        }
        std::cerr << " valid="
                  << (distance_state.valid_sites == native_valid)
                  << " differences="
                  << (distance_state.differences == native_differences)
                  << " sequence-data="
                  << (scan_state.sequence_data == native_sequence_data)
                  << " selected-sequence-data="
                  << selected_sequence_data_matches
                  << " result=" << ufdist_result << '/'
                  << (native_ufdist_result.empty() ?
                          -1 : native_ufdist_result[0])
                  << " breakpoint="
                  << (breakpoint_distance == native_breakpoint)
                  << " remainder="
                  << (remainder_distance == native_remainder) << '\n';
    }
    const char super_dist_magic[8] = {
        'S', 'U', 'P', 'D', 'I', 'S', 'T', '2'};
    const auto super_dist_fixture =
        load_rdp_sectioned_fixture<SuperDistP2CaptureHeader>(
            super_dist_p2_fixture_path, super_dist_magic);
    const auto native_region_differences =
        rdp_fixture_section<float>(super_dist_fixture, 102);
    const auto native_region_valid =
        rdp_fixture_section<float>(super_dist_fixture, 103);
    const auto native_region_distance =
        rdp_fixture_section<float>(super_dist_fixture, 104);
    auto region_distance_state = build_rdp_distance_state(
        scan_state, selected.beginning, selected.ending, false);
    const bool region_distance_matches =
        super_dist_fixture.header.next_no == scan_state.next_no &&
        region_distance_state.differences == native_region_differences &&
        region_distance_state.valid_sites == native_region_valid &&
        region_distance_state.distance == native_region_distance;

    const char check_matrix_magic[8] = {
        'C', 'H', 'K', 'M', 'A', 'T', 'P', '\0'};
    const auto check_matrix_fixture =
        load_rdp_sectioned_fixture<CheckMatrixPCaptureHeader>(
            check_matrix_fixture_path, check_matrix_magic);
    const auto& check_header = check_matrix_fixture.header;
    auto generated_matrices = finish_rdp_event_distances(
        scan_state.next_no, distance_state, region_distance_state);
    std::vector<unsigned char> minimum_pair{3, 3, 0};
    std::vector<unsigned char> sequence_pair(3, 0);
    const std::array<int, 3> role_outlier{2, 1, 0};
    float minimum_background = 1000000.0F;
    float minimum_region = 1000000.0F;
    int role_pair = 0;
    for (int first = 0; first < 2; ++first) {
        for (int second = first + 1; second < 3; ++second) {
            const auto offset =
                static_cast<std::size_t>(selected_sequences[first]) +
                static_cast<std::size_t>(selected_sequences[second]) *
                    (scan_state.next_no + 1);
            if (generated_matrices.background[offset] < minimum_background) {
                minimum_background = generated_matrices.background[offset];
                minimum_pair[0] = static_cast<unsigned char>(role_pair);
                sequence_pair[0] = static_cast<unsigned char>(first);
                sequence_pair[1] = static_cast<unsigned char>(second);
                sequence_pair[2] =
                    static_cast<unsigned char>(role_outlier[role_pair]);
            }
            if (generated_matrices.event_region[offset] < minimum_region) {
                minimum_region = generated_matrices.event_region[offset];
                minimum_pair[1] = static_cast<unsigned char>(role_pair);
            }
            ++role_pair;
        }
    }
    const int matrix_stride = scan_state.next_no + 1;
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        generated_matrices.background[
            static_cast<std::size_t>(sequence) * matrix_stride + sequence] =
            0.0F;
    }

    // CheckMatrixX's VB prepass, immediately before CheckMatrixP.
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        bool erase = false;
        for (int role = 0; role < 3; ++role) {
            const int selected_sequence = selected_sequences[role];
            if (selected_sequence <= scan_state.next_no &&
                sequence != selected_sequence) {
                const auto offset =
                    static_cast<std::size_t>(selected_sequence) +
                    static_cast<std::size_t>(sequence) * matrix_stride;
                if (distance_state.valid_sites[offset] -
                            region_distance_state.valid_sites[offset] <
                        check_header.minimum_sequence_size ||
                    region_distance_state.valid_sites[offset] <
                        check_header.sco) {
                    erase = true;
                    break;
                }
            }
        }
        if (erase) {
            for (int other = 0; other <= scan_state.next_no; ++other) {
                const auto first = static_cast<std::size_t>(sequence) +
                    static_cast<std::size_t>(other) * matrix_stride;
                const auto second = static_cast<std::size_t>(other) +
                    static_cast<std::size_t>(sequence) * matrix_stride;
                generated_matrices.background[first] = 3.0F;
                generated_matrices.background[second] = 3.0F;
                generated_matrices.event_region[first] = 3.0F;
                generated_matrices.event_region[second] = 3.0F;
            }
        }
    }

    std::vector<int> minimums(matrix_stride, 0);
    std::vector<unsigned char> missing_pair(
        static_cast<std::size_t>(matrix_stride) * matrix_stride, 0);
    std::vector<int> background_total(matrix_stride, 0);
    std::vector<int> region_total(matrix_stride, 0);
    const auto native_minimums_in =
        rdp_fixture_section<int>(check_matrix_fixture, 1);
    const auto native_sequences_in =
        rdp_fixture_section<int>(check_matrix_fixture, 2);
    const auto native_missing_in =
        rdp_fixture_section<unsigned char>(check_matrix_fixture, 3);
    const auto native_check_valid =
        rdp_fixture_section<float>(check_matrix_fixture, 4);
    const auto native_check_sub_valid =
        rdp_fixture_section<float>(check_matrix_fixture, 5);
    const auto native_background_in =
        rdp_fixture_section<float>(check_matrix_fixture, 6);
    const auto native_region_in =
        rdp_fixture_section<float>(check_matrix_fixture, 7);
    const auto native_background_total_in =
        rdp_fixture_section<int>(check_matrix_fixture, 8);
    const auto native_region_total_in =
        rdp_fixture_section<int>(check_matrix_fixture, 9);
    const auto float_vectors_match = [](const std::vector<float>& actual,
                                        const std::vector<float>& expected) {
        if (actual.size() != expected.size()) return false;
        for (std::size_t index = 0; index < actual.size(); ++index) {
            if (std::abs(actual[index] - expected[index]) > 1.0e-6F) {
                return false;
            }
        }
        return true;
    };
    const bool check_matrix_inputs_match =
        check_header.next_no == scan_state.next_no &&
        check_header.missing_pair_ub == scan_state.next_no &&
        check_header.valid_ub == scan_state.next_no &&
        check_header.sub_valid_ub == scan_state.next_no &&
        check_header.matrix_ub == scan_state.next_no &&
        minimums == native_minimums_in &&
        selected_sequences == native_sequences_in &&
        missing_pair == native_missing_in &&
        distance_state.valid_sites == native_check_valid &&
        region_distance_state.valid_sites == native_check_sub_valid &&
        float_vectors_match(
            generated_matrices.background, native_background_in) &&
        float_vectors_match(
            generated_matrices.event_region, native_region_in) &&
        background_total == native_background_total_in &&
        region_total == native_region_total_in;
    const int check_matrix_result = MathFuncs::MyMathFuncs::CheckMatrixP(
        minimums.data(), selected_sequences.data(), scan_state.next_no,
        check_header.sco, check_header.minimum_sequence_size,
        scan_state.next_no, missing_pair.data(), scan_state.next_no,
        distance_state.valid_sites.data(), scan_state.next_no,
        region_distance_state.valid_sites.data(), scan_state.next_no,
        generated_matrices.background.data(),
        generated_matrices.event_region.data(), background_total.data(),
        region_total.data());
    const auto native_minimums_out =
        rdp_fixture_section<int>(check_matrix_fixture, 101);
    const auto native_missing_out =
        rdp_fixture_section<unsigned char>(check_matrix_fixture, 102);
    const auto native_background_out =
        rdp_fixture_section<float>(check_matrix_fixture, 103);
    const auto native_region_out =
        rdp_fixture_section<float>(check_matrix_fixture, 104);
    const auto native_background_total_out =
        rdp_fixture_section<int>(check_matrix_fixture, 105);
    const auto native_region_total_out =
        rdp_fixture_section<int>(check_matrix_fixture, 106);
    const auto native_check_result =
        rdp_fixture_section<int>(check_matrix_fixture, 107);
    const bool check_matrix_matches = check_matrix_inputs_match &&
        native_check_result.size() == 1 &&
        check_matrix_result == native_check_result[0] &&
        minimums == native_minimums_out &&
        missing_pair == native_missing_out &&
        float_vectors_match(
            generated_matrices.background, native_background_out) &&
        float_vectors_match(
            generated_matrices.event_region, native_region_out) &&
        background_total == native_background_total_out &&
        region_total == native_region_total_out;
    if (!check_matrix_matches) {
        const auto report_float_difference = [&](const char* label,
                                                  const auto& actual,
                                                  const auto& expected) {
            std::size_t mismatch_count = 0;
            std::size_t first_mismatch = 0;
            for (std::size_t index = 0;
                 index < std::min(actual.size(), expected.size()); ++index) {
                if (actual[index] != expected[index]) {
                    if (mismatch_count == 0) first_mismatch = index;
                    ++mismatch_count;
                }
            }
            std::cerr << ' ' << label << '=' << mismatch_count;
            if (mismatch_count != 0) {
                std::cerr << '@' << first_mismatch << ':'
                          << actual[first_mismatch] << '/'
                          << expected[first_mismatch];
            }
        };
        std::cerr << "CheckMatrix integration mismatch: inputs="
                  << check_matrix_inputs_match << " minimums-in="
                  << (std::vector<int>(matrix_stride, 0) == native_minimums_in)
                  << " sequences="
                  << (selected_sequences == native_sequences_in)
                  << " valid="
                  << (distance_state.valid_sites == native_check_valid)
                  << " sub-valid="
                  << (region_distance_state.valid_sites ==
                      native_check_sub_valid)
                  << " background-in="
                  << (native_background_in ==
                      finish_rdp_event_distances(
                          scan_state.next_no, distance_state,
                          region_distance_state).background)
                  << " region-in="
                  << (native_region_in == region_distance_state.distance)
                  << " result=" << check_matrix_result << '/'
                  << (native_check_result.empty() ? -1 : native_check_result[0])
                  << ';';
        report_float_difference(
            "background-out", generated_matrices.background,
            native_background_out);
        report_float_difference(
            "region-out", generated_matrices.event_region,
            native_region_out);
        std::cerr << '\n';
    }

    const char make_nj_magic[8] = {
        'M', 'A', 'K', 'E', 'N', 'J', 'P', '2'};
    const auto make_nj_fixture =
        load_rdp_sectioned_fixture<MakeNJTreesP2CaptureHeader>(
            make_nj_fixture_path, make_nj_magic);
    const auto& nj_header = make_nj_fixture.header;
    std::vector<int> outlier{2, 1, 0};
    std::vector<int> redo_list(matrix_stride, 0);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        const auto diagonal = static_cast<std::size_t>(sequence) +
            static_cast<std::size_t>(sequence) * matrix_stride;
        if (generated_matrices.background[diagonal] == 3.0F) {
            redo_list[sequence] = 1;
        }
    }
    int local_last_sequence = -1;
    std::vector<int> trace_sequences(
        static_cast<std::size_t>(2) * matrix_stride, 0);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        if (redo_list[sequence] == 0) {
            ++local_last_sequence;
            trace_sequences[static_cast<std::size_t>(sequence) * 2] =
                local_last_sequence;
            trace_sequences[1 + static_cast<std::size_t>(local_last_sequence) *
                    2] = sequence;
        }
    }
    const int name_length = std::max(
        2, static_cast<int>(std::to_string(local_last_sequence).size()));
    std::vector<float> background_adjusted(
        static_cast<std::size_t>(matrix_stride) * matrix_stride, 0.0F);
    std::vector<float> region_adjusted(
        static_cast<std::size_t>(matrix_stride) * matrix_stride, 0.0F);
    const int local_stride = local_last_sequence + 1;
    std::vector<char> background_holder(
        static_cast<std::size_t>(local_stride) * 80 + 1, 0);
    std::vector<char> region_holder(
        static_cast<std::size_t>(local_stride) * 80 + 1, 0);
    std::vector<float> temporary_background(
        static_cast<std::size_t>(local_stride) * local_stride, 0.0F);
    std::vector<float> temporary_region(
        static_cast<std::size_t>(local_stride) * local_stride, 0.0F);

    const auto native_nj_sequences =
        rdp_fixture_section<int>(make_nj_fixture, 1);
    const auto native_nj_minimum_pair =
        rdp_fixture_section<unsigned char>(make_nj_fixture, 2);
    const auto native_nj_sequence_pair =
        rdp_fixture_section<unsigned char>(make_nj_fixture, 3);
    const auto native_nj_outlier =
        rdp_fixture_section<int>(make_nj_fixture, 4);
    const auto native_nj_trace =
        rdp_fixture_section<int>(make_nj_fixture, 5);
    const auto native_nj_background =
        rdp_fixture_section<float>(make_nj_fixture, 6);
    const auto native_nj_region =
        rdp_fixture_section<float>(make_nj_fixture, 7);
    const auto native_nj_background_adjusted =
        rdp_fixture_section<float>(make_nj_fixture, 8);
    const auto native_nj_region_adjusted =
        rdp_fixture_section<float>(make_nj_fixture, 9);
    const auto native_nj_redo =
        rdp_fixture_section<int>(make_nj_fixture, 10);
    const auto native_nj_background_holder =
        rdp_fixture_section<char>(make_nj_fixture, 11);
    const auto native_nj_region_holder =
        rdp_fixture_section<char>(make_nj_fixture, 12);
    const auto native_nj_temporary_background =
        rdp_fixture_section<float>(make_nj_fixture, 13);
    const auto native_nj_temporary_region =
        rdp_fixture_section<float>(make_nj_fixture, 14);
    const bool make_nj_inputs_match =
        nj_header.resolve_root == 1 &&
        nj_header.nseqs == local_last_sequence &&
        nj_header.next_no == scan_state.next_no &&
        nj_header.name_length == name_length &&
        nj_header.sequence_length == scan_state.sequence_length &&
        nj_header.trace_sequences_ub == 1 &&
        nj_header.first_matrix_ub == scan_state.next_no &&
        nj_header.second_matrix_ub == scan_state.next_no &&
        nj_header.first_adjusted_matrix_ub == scan_state.next_no &&
        nj_header.second_adjusted_matrix_ub == scan_state.next_no &&
        selected_sequences == native_nj_sequences &&
        minimum_pair == native_nj_minimum_pair &&
        sequence_pair == native_nj_sequence_pair &&
        outlier == native_nj_outlier &&
        trace_sequences == native_nj_trace &&
        float_vectors_match(
            generated_matrices.background, native_nj_background) &&
        float_vectors_match(
            generated_matrices.event_region, native_nj_region) &&
        background_adjusted == native_nj_background_adjusted &&
        region_adjusted == native_nj_region_adjusted &&
        redo_list == native_nj_redo &&
        background_holder == native_nj_background_holder &&
        region_holder == native_nj_region_holder &&
        temporary_background == native_nj_temporary_background &&
        temporary_region == native_nj_temporary_region;
    const int make_nj_result = MathFuncs::MyMathFuncs::MakeNJTreesP2(
        1, local_last_sequence, scan_state.next_no,
        selected_sequences.data(), minimum_pair.data(), sequence_pair.data(),
        nj_header.seed, name_length, scan_state.sequence_length, 1,
        outlier.data(), trace_sequences.data(), scan_state.next_no,
        generated_matrices.background.data(), scan_state.next_no,
        generated_matrices.event_region.data(), scan_state.next_no,
        background_adjusted.data(), scan_state.next_no,
        region_adjusted.data(), redo_list.data(), background_holder.data(),
        region_holder.data(), temporary_background.data(),
        temporary_region.data());
    const auto native_nj_result =
        rdp_fixture_section<int>(make_nj_fixture, 115);
    const bool make_nj_matches = make_nj_inputs_match &&
        native_nj_result.size() == 1 && make_nj_result == native_nj_result[0] &&
        selected_sequences ==
            rdp_fixture_section<int>(make_nj_fixture, 101) &&
        minimum_pair ==
            rdp_fixture_section<unsigned char>(make_nj_fixture, 102) &&
        sequence_pair ==
            rdp_fixture_section<unsigned char>(make_nj_fixture, 103) &&
        outlier == rdp_fixture_section<int>(make_nj_fixture, 104) &&
        trace_sequences == rdp_fixture_section<int>(make_nj_fixture, 105) &&
        float_vectors_match(
            generated_matrices.background,
            rdp_fixture_section<float>(make_nj_fixture, 106)) &&
        float_vectors_match(
            generated_matrices.event_region,
            rdp_fixture_section<float>(make_nj_fixture, 107)) &&
        float_vectors_match(
            background_adjusted,
            rdp_fixture_section<float>(make_nj_fixture, 108)) &&
        float_vectors_match(
            region_adjusted,
            rdp_fixture_section<float>(make_nj_fixture, 109)) &&
        redo_list == rdp_fixture_section<int>(make_nj_fixture, 110) &&
        background_holder ==
            rdp_fixture_section<char>(make_nj_fixture, 111) &&
        region_holder == rdp_fixture_section<char>(make_nj_fixture, 112) &&
        float_vectors_match(
            temporary_background,
            rdp_fixture_section<float>(make_nj_fixture, 113)) &&
        float_vectors_match(
            temporary_region,
            rdp_fixture_section<float>(make_nj_fixture, 114));
    if (!make_nj_matches) {
        const auto native_background_holder_out =
            rdp_fixture_section<char>(make_nj_fixture, 111);
        const auto native_region_holder_out =
            rdp_fixture_section<char>(make_nj_fixture, 112);
        const auto first_char_difference = [](const auto& actual,
                                              const auto& expected) {
            const auto common = std::min(actual.size(), expected.size());
            std::size_t index = 0;
            while (index < common && actual[index] == expected[index]) ++index;
            return index;
        };
        const auto holder_text = [](const std::vector<char>& holder) {
            std::string value;
            for (std::size_t index = 0;
                 index < std::min<std::size_t>(holder.size(), 160); ++index) {
                const unsigned char character = holder[index];
                if (character >= 32 && character < 127) {
                    value.push_back(static_cast<char>(character));
                } else {
                    value += '<' + std::to_string(character) + '>';
                }
            }
            return value;
        };
        std::cerr << "MakeNJTreesP2 integration mismatch: inputs="
                  << make_nj_inputs_match << " nseqs="
                  << local_last_sequence << '/' << nj_header.nseqs
                  << " min-pair=" << static_cast<int>(minimum_pair[0])
                  << ',' << static_cast<int>(minimum_pair[1]) << '/';
        for (std::size_t index = 0; index < native_nj_minimum_pair.size();
             ++index) {
            if (index != 0) std::cerr << ',';
            std::cerr << static_cast<int>(native_nj_minimum_pair[index]);
        }
        std::cerr << " seq-pair="
                  << static_cast<int>(sequence_pair[0]) << ','
                  << static_cast<int>(sequence_pair[1]) << ','
                  << static_cast<int>(sequence_pair[2]) << " trace="
                  << (trace_sequences == native_nj_trace) << " redo="
                  << (redo_list == native_nj_redo) << " holders="
                  << (background_holder == native_background_holder_out)
                  << ','
                  << (region_holder == native_region_holder_out)
                  << " holder-first="
                  << first_char_difference(
                      background_holder, native_background_holder_out)
                  << ','
                  << first_char_difference(
                      region_holder, native_region_holder_out)
                  << " region-tree='" << holder_text(region_holder)
                  << "'/'" << holder_text(native_region_holder_out) << '\''
                  << " result=" << make_nj_result << '/'
                  << (native_nj_result.empty() ? -1 : native_nj_result[0])
                  << '\n';
    }

    const char make_sdmp_magic[8] = {
        'M', 'A', 'K', 'E', 'S', 'D', 'M', 'P'};
    const auto make_sdmp_fixture =
        load_rdp_sectioned_fixture<MakeSDMP2CaptureHeader>(
            make_sdmp_fixture_path, make_sdmp_magic);
    const auto breakpoint_flanks = make_rdp_breakpoint_flanks(
        scan_state, selected.beginning, selected.ending,
        {selected_sequences[0], selected_sequences[1],
         selected_sequences[2]});
    std::vector<int> start_positions{
        breakpoint_flanks.positions[0], selected.beginning,
        breakpoint_flanks.positions[2], selected.ending + 1,
        selected.beginning};
    std::vector<int> end_positions{
        selected.beginning - 1, breakpoint_flanks.positions[1],
        selected.ending, breakpoint_flanks.positions[3], selected.ending};
    for (int region = 0; region < 4; ++region) {
        if (start_positions[region] > scan_state.sequence_length) {
            start_positions[region] -= scan_state.sequence_length;
        } else if (start_positions[region] < 1) {
            start_positions[region] += scan_state.sequence_length;
        }
        if (end_positions[region] > scan_state.sequence_length) {
            end_positions[region] -= scan_state.sequence_length;
        } else if (end_positions[region] < 1) {
            end_positions[region] += scan_state.sequence_length;
        }
    }
    std::vector<int> comparison_matrix{1, 0, 0, 2, 2, 1};
    const auto alignment_cells =
        static_cast<std::size_t>(matrix_stride) *
        (scan_state.sequence_length + 1);
    std::vector<unsigned char> missing_data(alignment_cells, 0);
    std::vector<double> summary_matrix(
        static_cast<std::size_t>(9) * matrix_stride, 0.0);
    std::vector<double> regional_distance_matrix(
        static_cast<std::size_t>(45) * matrix_stride, 0.0);
    const auto native_sdmp_sequence_data =
        rdp_fixture_section<short>(make_sdmp_fixture, 6);
    const bool sdmp_sequence_data_matches =
        scan_state.sequence_data.size() >= native_sdmp_sequence_data.size() &&
        std::equal(native_sdmp_sequence_data.begin(),
                   native_sdmp_sequence_data.end(),
                   scan_state.sequence_data.begin());
    const bool make_sdmp_inputs_match =
        make_sdmp_fixture.header.next_no == scan_state.next_no &&
        make_sdmp_fixture.header.sequence_length ==
            scan_state.sequence_length &&
        start_positions ==
            rdp_fixture_section<int>(make_sdmp_fixture, 1) &&
        end_positions == rdp_fixture_section<int>(make_sdmp_fixture, 2) &&
        selected_sequences ==
            rdp_fixture_section<int>(make_sdmp_fixture, 3) &&
        comparison_matrix ==
            rdp_fixture_section<int>(make_sdmp_fixture, 4) &&
        missing_data ==
            rdp_fixture_section<unsigned char>(make_sdmp_fixture, 5) &&
        sdmp_sequence_data_matches &&
        summary_matrix ==
            rdp_fixture_section<double>(make_sdmp_fixture, 7) &&
        regional_distance_matrix ==
            rdp_fixture_section<double>(make_sdmp_fixture, 8);
    const int make_sdmp_result = MathFuncs::MyMathFuncs::MakeSDMP2(
        scan_state.next_no, scan_state.sequence_length,
        start_positions.data(), end_positions.data(),
        selected_sequences.data(), comparison_matrix.data(),
        missing_data.data(), scan_state.sequence_data.data(),
        summary_matrix.data(), regional_distance_matrix.data());
    const auto native_make_sdmp_result =
        rdp_fixture_section<int>(make_sdmp_fixture, 103);
    const bool make_sdmp_matches = make_sdmp_inputs_match &&
        native_make_sdmp_result.size() == 1 &&
        make_sdmp_result == native_make_sdmp_result[0] &&
        summary_matrix ==
            rdp_fixture_section<double>(make_sdmp_fixture, 101) &&
        regional_distance_matrix ==
            rdp_fixture_section<double>(make_sdmp_fixture, 102);
    if (!make_sdmp_matches) {
        int differing_rows = 0;
        int differing_cells = 0;
        int ufdist_to_sdmp_cells = 0;
        std::array<int, 256> recode_counts{};
        for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
            bool row_differs = false;
            for (int position = 0; position <= scan_state.sequence_length;
                 ++position) {
                const auto offset = static_cast<std::size_t>(position) +
                    static_cast<std::size_t>(sequence) *
                        (scan_state.sequence_length + 1);
                if (scan_state.sequence_data[offset] !=
                    native_sdmp_sequence_data[offset]) {
                    row_differs = true;
                    ++differing_cells;
                    const int key =
                        (scan_state.sequence_data[offset] & 15) * 16 +
                        (native_sdmp_sequence_data[offset] & 15);
                    ++recode_counts[key];
                }
                if (native_sequence_data[offset] !=
                    native_sdmp_sequence_data[offset]) {
                    ++ufdist_to_sdmp_cells;
                }
            }
            if (row_differs) ++differing_rows;
        }
        std::cerr << "MakeSDMP2 integration mismatch: inputs="
                  << make_sdmp_inputs_match << " SP=";
        for (const int position : start_positions) {
            std::cerr << position << ',';
        }
        std::cerr << " EP=";
        for (const int position : end_positions) {
            std::cerr << position << ',';
        }
        std::cerr << " missing="
                  << (missing_data == rdp_fixture_section<unsigned char>(
                          make_sdmp_fixture, 5))
                  << " seqnum="
                  << sdmp_sequence_data_matches
                  << " seqnum-diff=" << differing_rows << "rows/"
                  << differing_cells << "cells ufdist-to-sdmp="
                  << ufdist_to_sdmp_cells << " sizes="
                  << scan_state.sequence_data.size() << '/'
                  << native_sdmp_sequence_data.size() << " recodes=";
        for (int key = 0; key < 256; ++key) {
            if (recode_counts[key] != 0) {
                std::cerr << (key / 16) << '>' << (key % 16) << ':'
                          << recode_counts[key] << ',';
            }
        }
        std::cerr
                  << " result=" << make_sdmp_result << '/'
                  << (native_make_sdmp_result.empty() ?
                          -1 : native_make_sdmp_result[0])
                  << '\n';
    }

    bool fill_rmat_matches = make_sdmp_matches;
    std::vector<unsigned char> correlation_positions(2, 0);
    std::array<std::vector<double>, 3> correlation_matrices;
    for (int target = 0; target < 3; ++target) {
        const char fill_rmat_magic[8] = {
            'F', 'I', 'L', 'L', 'R', 'M', 'A', 'T'};
        const auto fill_fixture =
            load_rdp_sectioned_fixture<FillRmatCaptureHeader>(
                fill_rmat_fixture_paths[target], fill_rmat_magic);
        const auto& fill_header = fill_fixture.header;
        auto& correlation_matrix = correlation_matrices[target];
        correlation_matrix.assign(
            static_cast<std::size_t>(3) * 6 * matrix_stride, 0.0);
        const bool fill_inputs_match = fill_header.y == target &&
            fill_header.next_no == scan_state.next_no &&
            fill_header.result_matrix_ub1 == 2 &&
            fill_header.result_matrix_ub2 == 5 &&
            fill_header.distance_matrix_ub1 == 2 &&
            fill_header.distance_matrix_ub2 == 4 &&
            fill_header.distance_matrix_ub3 == scan_state.next_no &&
            correlation_matrix ==
                rdp_fixture_section<double>(fill_fixture, 1) &&
            regional_distance_matrix ==
                rdp_fixture_section<double>(fill_fixture, 2) &&
            correlation_positions ==
                rdp_fixture_section<unsigned char>(fill_fixture, 3);
        const int fill_result = MathFuncs::MyMathFuncs::FillRmat(
            target, scan_state.next_no, 2, 5, 2, 4,
            scan_state.next_no, correlation_matrix.data(),
            regional_distance_matrix.data(), correlation_positions.data());
        const auto native_fill_result =
            rdp_fixture_section<int>(fill_fixture, 104);
        const bool target_matches = fill_inputs_match &&
            native_fill_result.size() == 1 &&
            fill_result == native_fill_result[0] &&
            correlation_matrix ==
                rdp_fixture_section<double>(fill_fixture, 101) &&
            regional_distance_matrix ==
                rdp_fixture_section<double>(fill_fixture, 102) &&
            correlation_positions ==
                rdp_fixture_section<unsigned char>(fill_fixture, 103);
        if (!target_matches) {
            std::cerr << "FillRmat integration mismatch: Y=" << target
                      << " inputs=" << fill_inputs_match << " positions="
                      << static_cast<int>(correlation_positions[0]) << ','
                      << static_cast<int>(correlation_positions[1])
                      << " result=" << fill_result << '/'
                      << (native_fill_result.empty() ?
                              -1 : native_fill_result[0])
                      << '\n';
        }
        fill_rmat_matches = fill_rmat_matches && target_matches;
    }
    const char calcr_magic[8] = {
        'C', 'A', 'L', 'C', 'R', '3', '\0', '\0'};
    const auto calcr_fixture =
        load_rdp_sectioned_fixture<CalCRChainCaptureHeader>(
            calcr_fixture_path, calcr_magic);
    const std::array<int, 3> correlation_sequences{
        selected_sequences[0], selected_sequences[1], selected_sequences[2]};
    const std::array<int, 6> correlation_comparison{1, 0, 0, 2, 2, 1};
    const auto correlation_state = calculate_rdp_correlations(
        scan_state.next_no, correlation_sequences, correlation_comparison,
        correlation_matrices);
    const bool calcr_matches = fill_rmat_matches &&
        calcr_fixture.header.next_no == scan_state.next_no &&
        std::vector<int>(correlation_sequences.begin(),
                         correlation_sequences.end()) ==
            rdp_fixture_section<int>(calcr_fixture, 1) &&
        std::vector<int>(correlation_comparison.begin(),
                         correlation_comparison.end()) ==
            rdp_fixture_section<int>(calcr_fixture, 2) &&
        correlation_matrices[0] ==
            rdp_fixture_section<double>(calcr_fixture, 3) &&
        correlation_matrices[1] ==
            rdp_fixture_section<double>(calcr_fixture, 4) &&
        correlation_matrices[2] ==
            rdp_fixture_section<double>(calcr_fixture, 5) &&
        correlation_state.correlation ==
            rdp_fixture_section<float>(calcr_fixture, 101) &&
        correlation_state.inversion ==
            rdp_fixture_section<float>(calcr_fixture, 102) &&
        correlation_state.tested_correlation ==
            rdp_fixture_section<float>(calcr_fixture, 103) &&
        std::vector<double>(correlation_state.intermediate.begin(),
                            correlation_state.intermediate.end()) ==
            rdp_fixture_section<double>(calcr_fixture, 104) &&
        std::vector<double>(correlation_state.results.begin(),
                            correlation_state.results.end()) ==
            rdp_fixture_section<double>(calcr_fixture, 105);
    if (!calcr_matches) {
        const auto report_calcr_difference = [&](const char* name,
                                                 const auto& actual,
                                                 const auto& expected) {
            std::size_t count = 0;
            std::size_t first = 0;
            for (std::size_t index = 0;
                 index < std::min(actual.size(), expected.size()); ++index) {
                if (actual[index] != expected[index]) {
                    if (count == 0) first = index;
                    ++count;
                }
            }
            std::cerr << ' ' << name << '=' << count;
            if (count != 0) {
                std::cerr << '@' << first << ':' << actual[first] << '/'
                          << expected[first];
            }
        };
        const auto native_calcr_correlation =
            rdp_fixture_section<float>(calcr_fixture, 101);
        const auto native_calcr_inversion =
            rdp_fixture_section<float>(calcr_fixture, 102);
        const auto native_calcr_tested =
            rdp_fixture_section<float>(calcr_fixture, 103);
        std::cerr << "CalCR chain integration mismatch: correlation="
                  << (correlation_state.correlation ==
                      rdp_fixture_section<float>(calcr_fixture, 101))
                  << " inversion="
                  << (correlation_state.inversion ==
                      rdp_fixture_section<float>(calcr_fixture, 102))
                  << " tested="
                  << (correlation_state.tested_correlation ==
                      rdp_fixture_section<float>(calcr_fixture, 103))
                  << " intermediate="
                  << (std::vector<double>(
                          correlation_state.intermediate.begin(),
                          correlation_state.intermediate.end()) ==
                      rdp_fixture_section<double>(calcr_fixture, 104))
                  << ';';
        report_calcr_difference(
            "correlation", correlation_state.correlation,
            native_calcr_correlation);
        report_calcr_difference(
            "inversion", correlation_state.inversion,
            native_calcr_inversion);
        report_calcr_difference(
            "tested", correlation_state.tested_correlation,
            native_calcr_tested);
        for (std::size_t index = 0;
             index < correlation_state.inversion.size(); ++index) {
            if (correlation_state.inversion[index] !=
                native_calcr_inversion[index]) {
                const int sequence = static_cast<int>(index / 9);
                const int remainder = static_cast<int>(index % 9);
                const int region = remainder / 3;
                const int target = remainder % 3;
                std::cerr << " inv-first-context=" << target << ',' << region
                          << ',' << sequence << " tested";
                for (int permutation = 0; permutation < 5; ++permutation) {
                    const auto tested_index = static_cast<std::size_t>(target) +
                        region * 3 + permutation * 9 + sequence * 45;
                    std::cerr << ':'
                              << correlation_state.tested_correlation[
                                     tested_index];
                }
                break;
            }
        }
        std::cerr << '\n';
    }

    const std::array<int, 4> local_starts{
        start_positions[0], start_positions[1], start_positions[2],
        start_positions[3]};
    const std::array<int, 4> local_ends{
        end_positions[0], end_positions[1], end_positions[2],
        end_positions[3]};
    auto correlation_decisions = finalize_rdp_correlations(
        scan_state.next_no, correlation_state, correlation_sequences,
        correlation_comparison, summary_matrix, regional_distance_matrix);
    if (correlation_decisions.warnings[0] != 0 &&
        correlation_decisions.warnings[1] != 0) {
        correlation_decisions.warnings[2] = 0;
    }
    const auto local_distance_panels = make_rdp_local_distance_panels(
        scan_state, local_starts, local_ends, correlation_sequences);
    apply_rdp_distance_warnings(
        scan_state.next_no, correlation_sequences, local_distance_panels,
        correlation_decisions.warnings);
    const auto good_comparisons = make_rdp_good_comparisons(
        scan_state,
        {breakpoint_flanks.positions[0], breakpoint_flanks.positions[1],
         breakpoint_flanks.positions[2], breakpoint_flanks.positions[3]});

    const auto make_small_matrix = [&](const std::vector<float>& matrix) {
        std::vector<float> result(
            static_cast<std::size_t>(3) * matrix_stride, 0.0F);
        for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
            for (int role = 0; role < 3; ++role) {
                result[role + sequence * 3] = matrix[
                    selected_sequences[role] + sequence * matrix_stride];
            }
        }
        return result;
    };
    const auto first_direct_small =
        make_small_matrix(generated_matrices.background);
    const auto second_direct_small =
        make_small_matrix(generated_matrices.event_region);
    auto first_adjusted_small = make_small_matrix(background_adjusted);
    auto second_adjusted_small = make_small_matrix(region_adjusted);

    std::array<unsigned char, 2> final_minimum_pair{};
    std::array<float, 2> final_minimum_distance{
        1000000.0F, 1000000.0F};
    int final_pair = 0;
    for (int first = 0; first < 2; ++first) {
        for (int second = first + 1; second < 3; ++second) {
            const auto offset = static_cast<std::size_t>(first) +
                static_cast<std::size_t>(selected_sequences[second]) * 3;
            if (first_adjusted_small[offset] < final_minimum_distance[0]) {
                final_minimum_distance[0] = first_adjusted_small[offset];
                final_minimum_pair[0] =
                    static_cast<unsigned char>(final_pair);
            }
            if (second_adjusted_small[offset] < final_minimum_distance[1]) {
                final_minimum_distance[1] = second_adjusted_small[offset];
                final_minimum_pair[1] =
                    static_cast<unsigned char>(final_pair);
            }
            ++final_pair;
        }
    }
    const auto role_lists = make_rdp_role_lists(final_minimum_pair);
    auto first_collapsed = first_adjusted_small;
    auto second_collapsed = second_adjusted_small;
    const auto acceptable_correlations = make_rdp_acceptable_correlations(
        scan_state.next_no, correlation_sequences, role_lists.inside,
        first_direct_small, second_direct_small,
        first_adjusted_small, second_adjusted_small, first_collapsed,
        second_collapsed);
    const std::vector<unsigned char> dont_redo(
        static_cast<std::size_t>(3) * matrix_stride, 0);

    const char make_rlist_magic[8] = {
        'M', 'R', 'L', 'I', 'S', 'T', '1', '\0'};
    const auto make_rlist_fixture =
        load_rdp_sectioned_fixture<MakeRListCaptureHeader>(
            make_rlist_fixture_path, make_rlist_magic);
    const bool make_rlist_inputs_match =
        make_rlist_fixture.header.next_no == scan_state.next_no &&
        selected_sequences ==
            rdp_fixture_section<int>(make_rlist_fixture, 1) &&
        good_comparisons ==
            rdp_fixture_section<int>(make_rlist_fixture, 2) &&
        correlation_decisions.correlations.inversion ==
            rdp_fixture_section<float>(make_rlist_fixture, 3) &&
        correlation_decisions.correlations.correlation ==
            rdp_fixture_section<float>(make_rlist_fixture, 4) &&
        dont_redo == rdp_fixture_section<unsigned char>(
            make_rlist_fixture, 5) &&
        acceptable_correlations == rdp_fixture_section<unsigned char>(
            make_rlist_fixture, 6) &&
        std::vector<unsigned char>(
            correlation_decisions.warnings.begin(),
            correlation_decisions.warnings.end()) ==
            rdp_fixture_section<unsigned char>(make_rlist_fixture, 7);
    const auto candidate_lists = make_rdp_candidate_lists(
        scan_state.next_no, good_comparisons, correlation_sequences,
        correlation_decisions, dont_redo, acceptable_correlations);
    const auto probability_vectors_match = [](const auto& actual,
                                                const auto& expected) {
        if (actual.size() != expected.size()) return false;
        for (std::size_t index = 0; index < actual.size(); ++index) {
            if (std::abs(actual[index] - expected[index]) > 1.0e-12) {
                return false;
            }
        }
        return true;
    };
    const bool make_rlist_matches = make_rlist_inputs_match &&
        std::vector<int>(candidate_lists.last.begin(),
                         candidate_lists.last.end()) ==
            rdp_fixture_section<int>(make_rlist_fixture, 101) &&
        candidate_lists.list ==
            rdp_fixture_section<int>(make_rlist_fixture, 102) &&
        candidate_lists.inverse ==
            rdp_fixture_section<int>(make_rlist_fixture, 103) &&
        probability_vectors_match(
            candidate_lists.totals,
            rdp_fixture_section<double>(make_rlist_fixture, 104)) &&
        probability_vectors_match(
            candidate_lists.list_scores,
            rdp_fixture_section<double>(make_rlist_fixture, 105));
    if (!make_rlist_matches) {
        const auto expected_warnings =
            rdp_fixture_section<unsigned char>(make_rlist_fixture, 7);
        const auto expected_last =
            rdp_fixture_section<int>(make_rlist_fixture, 101);
        const auto count_differences = [](const auto& actual,
                                          const auto& expected) {
            std::size_t count = actual.size() == expected.size() ? 0 : 1;
            for (std::size_t index = 0;
                 index < std::min(actual.size(), expected.size()); ++index) {
                if (actual[index] != expected[index]) ++count;
            }
            return count;
        };
        std::cerr << "MakeRList integration mismatch: inputs="
                  << make_rlist_inputs_match << " selected="
                  << (selected_sequences ==
                      rdp_fixture_section<int>(make_rlist_fixture, 1))
                  << " good="
                  << count_differences(
                      good_comparisons,
                      rdp_fixture_section<int>(make_rlist_fixture, 2))
                  << " rcorr="
                  << count_differences(
                      correlation_decisions.correlations.correlation,
                      rdp_fixture_section<float>(make_rlist_fixture, 4))
                  << " rinv="
                  << count_differences(
                      correlation_decisions.correlations.inversion,
                      rdp_fixture_section<float>(make_rlist_fixture, 3))
                  << " acceptable="
                  << count_differences(
                      acceptable_correlations,
                      rdp_fixture_section<unsigned char>(
                          make_rlist_fixture, 6))
                  << " warnings=";
        for (const auto warning : correlation_decisions.warnings) {
            std::cerr << static_cast<int>(warning);
        }
        std::cerr << '/';
        for (const auto warning : expected_warnings) {
            std::cerr << static_cast<int>(warning);
        }
        std::cerr << " last=";
        for (const auto value : candidate_lists.last) std::cerr << value << ',';
        std::cerr << '/';
        for (const auto value : expected_last) std::cerr << value << ',';
        std::cerr << " list="
                  << count_differences(
                      candidate_lists.list,
                      rdp_fixture_section<int>(make_rlist_fixture, 102))
                  << " inverse="
                  << count_differences(
                      candidate_lists.inverse,
                      rdp_fixture_section<int>(make_rlist_fixture, 103))
                  << " totals="
                  << count_differences(
                      candidate_lists.totals,
                      rdp_fixture_section<double>(make_rlist_fixture, 104))
                  << " scores="
                  << count_differences(
                      candidate_lists.list_scores,
                      rdp_fixture_section<double>(make_rlist_fixture, 105));
        const auto expected_totals =
            rdp_fixture_section<double>(make_rlist_fixture, 104);
        for (std::size_t index = 0;
             index < candidate_lists.totals.size(); ++index) {
            if (candidate_lists.totals[index] != expected_totals[index]) {
                std::cerr << " first-total=" << index << ':'
                          << candidate_lists.totals[index] << '/'
                          << expected_totals[index];
                break;
            }
        }
        std::cerr << '\n';
    }

    const char find_actual_magic[8] = {
        'F', 'A', 'E', 'V', 'E', 'N', 'T', '1'};
    const auto find_actual_fixture =
        load_rdp_sectioned_fixture<FindActualEventsCaptureHeader>(
            find_actual_events_fixture_path, find_actual_magic);
    std::array<int, 6> actual_starts{};
    std::array<int, 6> actual_ends{};
    for (int index = 0; index < 5; ++index) {
        actual_starts[index] = start_positions[index];
        actual_ends[index] = end_positions[index];
    }
    const auto actual_resolution = resolve_rdp_actual_events(
        scan_state.sequence_length, scan_state.next_no,
        correlation_sequences, correlation_comparison,
        actual_starts, actual_ends, correlation_decisions, candidate_lists,
        dont_redo, events, scan_state.next_no);
    const auto as_vector = [](const auto& values) {
        return std::vector<typename std::decay_t<decltype(values)>::value_type>(
            values.begin(), values.end());
    };
    bool find_actual_matches = make_rlist_matches &&
        find_actual_fixture.header.sequence_length ==
            scan_state.sequence_length &&
        find_actual_fixture.header.next_no == scan_state.next_no;
    for (int role = 0; role < 3; ++role) {
        const int base = role * 1000;
        const auto& call = actual_resolution.calls[role];
        const auto expected_inversion =
            rdp_fixture_section<unsigned char>(
                find_actual_fixture, base + 18);
        const auto expected_membership =
            rdp_fixture_section<unsigned char>(
                find_actual_fixture, base + 24);
        const bool inputs_match =
            selected_sequences ==
                rdp_fixture_section<int>(find_actual_fixture, base + 1) &&
            as_vector(correlation_comparison) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 2) &&
            as_vector(call.region_sizes_before) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 3) &&
            dont_redo == rdp_fixture_section<unsigned char>(
                find_actual_fixture, base + 4) &&
            call.breakpoint_matches_before ==
                rdp_fixture_section<int>(find_actual_fixture, base + 5) &&
            call.best_matches_before ==
                rdp_fixture_section<float>(find_actual_fixture, base + 6) &&
            probability_vectors_match(call.acceptable_sequences_before,
                rdp_fixture_section<double>(find_actual_fixture, base + 7)) &&
            call.found_before ==
                rdp_fixture_section<int>(find_actual_fixture, base + 8) &&
            as_vector(actual_starts) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 9) &&
            as_vector(actual_ends) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 10) &&
            actual_resolution.correlations.correlations.correlation ==
                rdp_fixture_section<float>(find_actual_fixture, base + 11) &&
            actual_resolution.event_overlap_mask ==
                rdp_fixture_section<int>(find_actual_fixture, base + 12) &&
            actual_resolution.beginning_overlap_mask ==
                rdp_fixture_section<int>(find_actual_fixture, base + 13) &&
            actual_resolution.ending_overlap_mask ==
                rdp_fixture_section<int>(find_actual_fixture, base + 14) &&
            as_vector(call.candidate_scratch_before) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 15) &&
            as_vector(call.candidate_last_before) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 16) &&
            call.candidate_list_before ==
                rdp_fixture_section<int>(find_actual_fixture, base + 17) &&
            call.inversion_state_before == expected_inversion &&
            as_vector(call.trace_sequences_before) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 19) &&
            as_vector(call.match_before) ==
                rdp_fixture_section<double>(find_actual_fixture, base + 20) &&
            std::vector<short>(events.current_xover.begin(),
                               events.current_xover.end()) ==
                rdp_fixture_section<short>(find_actual_fixture, base + 21) &&
            as_vector(call.sequence_scratch_before) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 22) &&
            as_vector(call.tried_permutations_before) ==
                rdp_fixture_section<unsigned char>(
                    find_actual_fixture, base + 23) &&
            call.role_membership_before == expected_membership;
        const bool outputs_match =
            std::vector<int>{call.result} ==
                rdp_fixture_section<int>(find_actual_fixture, base + 101) &&
            as_vector(call.region_sizes_after) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 102) &&
            call.breakpoint_matches_after ==
                rdp_fixture_section<int>(find_actual_fixture, base + 103) &&
            call.best_matches_after ==
                rdp_fixture_section<float>(find_actual_fixture, base + 104) &&
            probability_vectors_match(call.acceptable_sequences_after,
                rdp_fixture_section<double>(
                    find_actual_fixture, base + 105)) &&
            call.found_after ==
                rdp_fixture_section<int>(find_actual_fixture, base + 106) &&
            as_vector(call.candidate_scratch_after) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 107) &&
            as_vector(call.trace_sequences_after) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 108) &&
            as_vector(call.match_after) ==
                rdp_fixture_section<double>(find_actual_fixture, base + 109) &&
            as_vector(call.sequence_scratch_after) ==
                rdp_fixture_section<int>(find_actual_fixture, base + 110) &&
            as_vector(call.tried_permutations_after) ==
                rdp_fixture_section<unsigned char>(
                    find_actual_fixture, base + 111);
        if (!inputs_match || !outputs_match) {
            const auto differences = [](const auto& actual,
                                        const auto& expected) {
                std::size_t result = actual.size() == expected.size() ? 0 : 1;
                for (std::size_t index = 0;
                     index < std::min(actual.size(), expected.size()); ++index) {
                    if (actual[index] != expected[index]) ++result;
                }
                return result;
            };
            const auto maximum_difference = [](const auto& actual,
                                               const auto& expected) {
                double result = 0.0;
                for (std::size_t index = 0;
                     index < std::min(actual.size(), expected.size()); ++index) {
                    const double difference = static_cast<double>(std::abs(
                        static_cast<double>(actual[index]) -
                        static_cast<double>(expected[index])));
                    result = std::max(result, difference);
                }
                return result;
            };
            const auto expected_ok_in = rdp_fixture_section<double>(
                find_actual_fixture, base + 7);
            const auto expected_ok_out = rdp_fixture_section<double>(
                find_actual_fixture, base + 105);
            std::cerr << "FindActualEvents role " << role
                      << " mismatch: inputs=" << inputs_match
                      << " outputs=" << outputs_match
                      << " rnum=" << differences(
                             as_vector(call.candidate_last_before),
                             rdp_fixture_section<int>(
                                 find_actual_fixture, base + 16))
                      << " rlist=" << differences(
                             call.candidate_list_before,
                             rdp_fixture_section<int>(
                                 find_actual_fixture, base + 17))
                      << " invs=" << differences(
                             call.inversion_state_before, expected_inversion)
                      << " ok-in=" << differences(
                             call.acceptable_sequences_before,
                             expected_ok_in)
                      << '/' << maximum_difference(
                             call.acceptable_sequences_before,
                             expected_ok_in)
                      << " found=" << differences(
                             call.found_after,
                             rdp_fixture_section<int>(
                                 find_actual_fixture, base + 106))
                      << " ok-out=" << differences(
                             call.acceptable_sequences_after,
                             expected_ok_out)
                      << '/' << maximum_difference(
                             call.acceptable_sequences_after,
                             expected_ok_out)
                      << " bp=" << differences(
                             call.breakpoint_matches_after,
                             rdp_fixture_section<int>(
                                 find_actual_fixture, base + 103))
                      << " best=" << differences(
                             call.best_matches_after,
                             rdp_fixture_section<float>(
                                 find_actual_fixture, base + 104))
                      << '\n';
        }
        find_actual_matches =
            find_actual_matches && inputs_match && outputs_match;
    }
    const char strip_dup_magic[8] = {
        'S', 'T', 'R', 'I', 'P', 'D', 'I', '1'};
    const auto strip_dup_fixture =
        load_rdp_sectioned_fixture<StripDupInvCaptureHeader>(
            strip_dup_inv_fixture_path, strip_dup_magic);
    const bool strip_dup_matches = find_actual_matches &&
        strip_dup_fixture.header.next_no == scan_state.next_no &&
        as_vector(actual_resolution.candidate_last_before_strip) ==
            rdp_fixture_section<int>(strip_dup_fixture, 1) &&
        actual_resolution.candidate_list_before_strip ==
            rdp_fixture_section<int>(strip_dup_fixture, 2) &&
        actual_resolution.candidate_inverse_before_strip ==
            rdp_fixture_section<int>(strip_dup_fixture, 3) &&
        as_vector(actual_resolution.candidates.last) ==
            rdp_fixture_section<int>(strip_dup_fixture, 101) &&
        actual_resolution.candidates.list ==
            rdp_fixture_section<int>(strip_dup_fixture, 102) &&
        actual_resolution.candidates.inverse ==
            rdp_fixture_section<int>(strip_dup_fixture, 103) &&
        as_vector(actual_resolution.inversion_penalty) ==
            rdp_fixture_section<int>(strip_dup_fixture, 104);
    if (!strip_dup_matches) {
        std::cerr << "StripDupInv mismatch: before-last="
                  << (as_vector(actual_resolution.candidate_last_before_strip) ==
                      rdp_fixture_section<int>(strip_dup_fixture, 1))
                  << " before-list="
                  << (actual_resolution.candidate_list_before_strip ==
                      rdp_fixture_section<int>(strip_dup_fixture, 2))
                  << " after-last="
                  << (as_vector(actual_resolution.candidates.last) ==
                      rdp_fixture_section<int>(strip_dup_fixture, 101))
                  << " after-list="
                  << (actual_resolution.candidates.list ==
                      rdp_fixture_section<int>(strip_dup_fixture, 102))
                  << " penalty="
                  << (as_vector(actual_resolution.inversion_penalty) ==
                      rdp_fixture_section<int>(strip_dup_fixture, 104))
                  << '\n';
    }
    const char rcompat_magic[8] = {
        'R', 'C', 'O', 'M', 'P', 'A', 'T', '1'};
    const auto rcompat_fixture =
        load_rdp_sectioned_fixture<RCompatCaptureHeader>(
            rcompat_fixture_path, rcompat_magic);
    const auto initial_tree_compatibility = evaluate_rdp_tree_compatibility(
        scan_state.next_no, correlation_sequences, correlation_comparison,
        actual_resolution.inversion_penalty,
        actual_resolution.candidates.last,
        actual_resolution.candidates.list, good_comparisons,
        background_adjusted, region_adjusted);
    const auto tied_compatibility = [](const std::array<int, 3>& values) {
        return values[0] == values[1] && values[0] == values[2];
    };
    std::vector<std::vector<float>> collapsed_matrices;
    bool repeated_background_primary = false;
    bool repeated_region_primary = false;
    for (unsigned int call_index = 6;
         call_index < rcompat_fixture.header.calls; call_index += 3) {
        const int base = static_cast<int>(call_index) * 1000;
        const auto call_last =
            rdp_fixture_section<int>(rcompat_fixture, base + 7);
        const auto call_list =
            rdp_fixture_section<int>(rcompat_fixture, base + 11);
        const auto matrix =
            rdp_fixture_section<float>(rcompat_fixture, base + 13);
        if (call_last != as_vector(actual_resolution.candidates.last) ||
            call_list != actual_resolution.candidates.list) {
            continue;
        }
        if (matrix == background_adjusted) {
            repeated_background_primary = true;
            continue;
        }
        if (matrix == region_adjusted) {
            repeated_region_primary = true;
            continue;
        }
        if (std::find(collapsed_matrices.begin(), collapsed_matrices.end(),
                      matrix) == collapsed_matrices.end()) {
            collapsed_matrices.push_back(matrix);
        }
    }
    std::vector<float> collapsed_background;
    std::vector<float> collapsed_region;
    auto collapsed = collapsed_matrices.begin();
    if (tied_compatibility(
            initial_tree_compatibility.background_compatibility)) {
        if (repeated_background_primary) {
            collapsed_background = background_adjusted;
        } else if (collapsed != collapsed_matrices.end()) {
            collapsed_background = *collapsed++;
        }
    }
    if (tied_compatibility(initial_tree_compatibility.region_compatibility)) {
        if (repeated_region_primary) {
            collapsed_region = region_adjusted;
        } else if (collapsed != collapsed_matrices.end()) {
            collapsed_region = *collapsed++;
        }
    }
    const auto tree_compatibility_flow = run_rdp_tree_compatibility_flow(
        scan_state.sequence_length, scan_state.next_no, selected.beginning,
        selected.ending, correlation_sequences, correlation_comparison,
        actual_resolution.inversion_penalty,
        actual_resolution.candidates.last,
        actual_resolution.candidates.list, good_comparisons,
        background_adjusted, region_adjusted, collapsed_background,
        collapsed_region, events);
    bool rcompat_matches = strip_dup_matches &&
        rcompat_fixture.header.next_no == scan_state.next_no &&
        rcompat_fixture.header.calls >= 6;
    bool rcompat_flow_matches =
        tree_compatibility_flow.calls.size() <= rcompat_fixture.header.calls &&
        (tree_compatibility_flow.calls.size() == rcompat_fixture.header.calls ||
         tree_compatibility_flow.calls.size() + 6 ==
             rcompat_fixture.header.calls);
    for (unsigned int call_index = 0;
         call_index < rcompat_fixture.header.calls; ++call_index) {
        const int base = call_index * 1000;
        const auto metadata =
            rdp_fixture_section<int>(rcompat_fixture, base + 1);
        const auto call_sequences =
            rdp_fixture_section<int>(rcompat_fixture, base + 2);
        const auto call_comparison =
            rdp_fixture_section<int>(rcompat_fixture, base + 3);
        const auto compatibility_before =
            rdp_fixture_section<int>(rcompat_fixture, base + 4);
        const auto reverse_before =
            rdp_fixture_section<int>(rcompat_fixture, base + 5);
        const auto call_penalty =
            rdp_fixture_section<int>(rcompat_fixture, base + 6);
        const auto call_last =
            rdp_fixture_section<int>(rcompat_fixture, base + 7);
        const auto nonrecombinant_before =
            rdp_fixture_section<int>(rcompat_fixture, base + 8);
        const auto call_good =
            rdp_fixture_section<int>(rcompat_fixture, base + 9);
        const auto expected_done =
            rdp_fixture_section<int>(rcompat_fixture, base + 10);
        const auto call_list =
            rdp_fixture_section<int>(rcompat_fixture, base + 11);
        auto nonrecombinant_list =
            rdp_fixture_section<int>(rcompat_fixture, base + 12);
        const auto matrix =
            rdp_fixture_section<float>(rcompat_fixture, base + 13);
        const auto list_distances =
            rdp_fixture_section<double>(rcompat_fixture, base + 14);
        const bool dimensions_match = metadata.size() == 2 &&
            call_sequences.size() == 3 && call_comparison.size() == 6 &&
            compatibility_before.size() == 3 && reverse_before.size() == 3 &&
            call_penalty.size() == 3 && call_last.size() == 3 &&
            nonrecombinant_before.size() == 3 && list_distances.size() == 3;
        RdpTreeCompatibilityCallState call;
        if (dimensions_match) {
            std::array<int, 3> compatibility{
                compatibility_before[0], compatibility_before[1],
                compatibility_before[2]};
            std::array<int, 3> reverse{
                reverse_before[0], reverse_before[1], reverse_before[2]};
            std::array<int, 3> nonrecombinant_last{
                nonrecombinant_before[0], nonrecombinant_before[1],
                nonrecombinant_before[2]};
            call = make_rdp_tree_compatibility_call(
                scan_state.next_no,
                {call_sequences[0], call_sequences[1], call_sequences[2]},
                {call_comparison[0], call_comparison[1], call_comparison[2],
                 call_comparison[3], call_comparison[4], call_comparison[5]},
                metadata[1],
                {call_penalty[0], call_penalty[1], call_penalty[2]},
                {call_last[0], call_last[1], call_last[2]}, call_list,
                call_good, matrix,
                {list_distances[0], list_distances[1], list_distances[2]},
                compatibility, reverse, nonrecombinant_last,
                nonrecombinant_list);
        }
        bool integrated_matches = true;
        if (call_index < tree_compatibility_flow.calls.size() &&
            dimensions_match) {
            const auto& integrated =
                tree_compatibility_flow.calls[call_index];
            integrated_matches =
                integrated.role == metadata[1] &&
                as_vector(integrated.compatibility_before) ==
                    compatibility_before &&
                as_vector(integrated.reverse_compatibility_before) ==
                    reverse_before &&
                as_vector(integrated.nonrecombinant_last_before) ==
                    nonrecombinant_before &&
                integrated.done_before == expected_done &&
                integrated.nonrecombinant_list_before ==
                    rdp_fixture_section<int>(rcompat_fixture, base + 12) &&
                as_vector(integrated.list_distances) == list_distances &&
                integrated.compatibility_after == call.compatibility_after &&
                integrated.reverse_compatibility_after ==
                    call.reverse_compatibility_after &&
                integrated.nonrecombinant_last_after ==
                    call.nonrecombinant_last_after;
            rcompat_flow_matches =
                rcompat_flow_matches && integrated_matches;
        }
        const bool call_matches =
            dimensions_match && integrated_matches &&
            as_vector(call.compatibility_before) ==
                compatibility_before &&
            as_vector(call.reverse_compatibility_before) ==
                reverse_before &&
            as_vector(call.nonrecombinant_last_before) ==
                nonrecombinant_before &&
            call.done_before == expected_done &&
            call.nonrecombinant_list_before ==
                rdp_fixture_section<int>(rcompat_fixture, base + 12) &&
            as_vector(call.list_distances) ==
                list_distances &&
            as_vector(call.compatibility_after) ==
                rdp_fixture_section<int>(rcompat_fixture, base + 101) &&
            as_vector(call.reverse_compatibility_after) ==
                rdp_fixture_section<int>(rcompat_fixture, base + 102) &&
            as_vector(call.nonrecombinant_last_after) ==
                rdp_fixture_section<int>(rcompat_fixture, base + 103);
        if (!call_matches) {
            std::cerr << "MakeRCompat call " << call_index
                      << " mismatch: dimensions=" << dimensions_match
                      << " integrated=" << integrated_matches
                      << " ldist="
                      << (as_vector(call.list_distances) ==
                          rdp_fixture_section<double>(
                              rcompat_fixture, base + 14))
                      << " before="
                      << (as_vector(call.compatibility_before) ==
                          rdp_fixture_section<int>(
                              rcompat_fixture, base + 4))
                      << " after="
                      << (as_vector(call.compatibility_after) ==
                          rdp_fixture_section<int>(
                              rcompat_fixture, base + 101))
                      << " reverse="
                      << (as_vector(call.reverse_compatibility_after) ==
                          rdp_fixture_section<int>(
                              rcompat_fixture, base + 102))
                      << " nrnum="
                      << (as_vector(call.nonrecombinant_last_after) ==
                          rdp_fixture_section<int>(
                              rcompat_fixture, base + 103))
                      << '\n';
        }
        rcompat_matches = rcompat_matches && call_matches;
    }
    if (!rcompat_flow_matches) {
        const auto diagnostic_matrix = rdp_fixture_section<float>(
            rcompat_fixture, 6013);
        const auto diagnostic_last = rdp_fixture_section<int>(
            rcompat_fixture, 6007);
        const auto diagnostic_list = rdp_fixture_section<int>(
            rcompat_fixture, 6011);
        std::cerr << "RDP compatibility flow mismatch: generated "
                  << tree_compatibility_flow.calls.size() << " of "
                  << rcompat_fixture.header.calls << " captured calls; F="
                  << tree_compatibility_flow.background[0] << ','
                  << tree_compatibility_flow.background[1] << ','
                  << tree_compatibility_flow.background[2] << " S="
                  << tree_compatibility_flow.region[0] << ','
                  << tree_compatibility_flow.region[1] << ','
                  << tree_compatibility_flow.region[2] << " FC="
                  << tree_compatibility_flow.background_secondary[0] << ','
                  << tree_compatibility_flow.background_secondary[1] << ','
                  << tree_compatibility_flow.background_secondary[2] << " SC="
                  << tree_compatibility_flow.region_secondary[0] << ','
                  << tree_compatibility_flow.region_secondary[1] << ','
                  << tree_compatibility_flow.region_secondary[2]
                  << "; call6 matrix F/FC/S/SC="
                  << (diagnostic_matrix == background_adjusted) << '/'
                  << (diagnostic_matrix == temporary_background) << '/'
                  << (diagnostic_matrix == region_adjusted) << '/'
                  << (diagnostic_matrix == temporary_region)
                  << " lists primary/sets="
                  << (diagnostic_last ==
                      as_vector(actual_resolution.candidates.last) &&
                      diagnostic_list == actual_resolution.candidates.list)
                  << '/'
                  << (diagnostic_last ==
                      as_vector(tree_compatibility_flow.event_sets.candidate_last) &&
                      diagnostic_list ==
                          tree_compatibility_flow.event_sets.candidate_list)
                  << '\n';
    }
    const char phpr_magic[8] = {
        'P', 'H', 'P', 'R', 'S', 'C', 'O', '1'};
    const auto phpr_fixture =
        load_rdp_sectioned_fixture<PhPrScoreCaptureHeader>(
            phpr_fixture_path, phpr_magic);
    bool phpr_matches =
        phpr_fixture.header.next_no == scan_state.next_no &&
        phpr_fixture.header.calls >= 2 &&
        phpr_fixture.header.calls <= 3;
    for (unsigned int call_index = 0;
         call_index < phpr_fixture.header.calls; ++call_index) {
        const int base = static_cast<int>(call_index) * 1000;
        const auto offset_values =
            rdp_fixture_section<double>(phpr_fixture, base + 2);
        const auto call_sequences =
            rdp_fixture_section<int>(phpr_fixture, base + 3);
        const auto done_this =
            rdp_fixture_section<int>(phpr_fixture, base + 4);
        const auto expected_trace =
            rdp_fixture_section<int>(phpr_fixture, base + 5);
        const auto first_matrix =
            rdp_fixture_section<float>(phpr_fixture, base + 6);
        const auto second_matrix =
            rdp_fixture_section<float>(phpr_fixture, base + 7);
        const auto expected_scores =
            rdp_fixture_section<double>(phpr_fixture, base + 8);
        const auto expected_sub_scores =
            rdp_fixture_section<double>(phpr_fixture, base + 9);
        const auto expected_sub_distances =
            rdp_fixture_section<double>(phpr_fixture, base + 10);
        const bool dimensions_match = offset_values.size() == 1 &&
            call_sequences.size() == 3;
        RdpPhylProScoreState actual;
        if (dimensions_match) {
            actual = make_rdp_phylpro_scores(
                scan_state.next_no, offset_values[0], done_this,
                {call_sequences[0], call_sequences[1], call_sequences[2]},
                first_matrix, second_matrix);
        }
        const bool call_matches = dimensions_match &&
            actual.trace_involved == expected_trace &&
            probability_vectors_match(
                as_vector(actual.scores), expected_scores) &&
            probability_vectors_match(
                as_vector(actual.sub_scores), expected_sub_scores) &&
            probability_vectors_match(
                as_vector(actual.sub_distance_scores),
                expected_sub_distances);
        if (!call_matches) {
            std::cerr << "MakePhPrScore call " << call_index
                      << " mismatch: dimensions=" << dimensions_match
                      << " trace="
                      << (dimensions_match &&
                          actual.trace_involved == expected_trace)
                      << " score="
                      << (dimensions_match &&
                          probability_vectors_match(
                              as_vector(actual.scores), expected_scores))
                      << " sub="
                      << (dimensions_match &&
                          probability_vectors_match(
                              as_vector(actual.sub_scores),
                              expected_sub_scores))
                      << " distance="
                      << (dimensions_match &&
                          probability_vectors_match(
                              as_vector(actual.sub_distance_scores),
                              expected_sub_distances))
                      << " values=";
            if (dimensions_match) {
                for (int role = 0; role < 3; ++role) {
                    std::cerr << actual.scores[role] << '/'
                              << expected_scores[role] << ',';
                }
                std::cerr << " sub=";
                for (int role = 0; role < 3; ++role) {
                    std::cerr << actual.sub_scores[role] << '/'
                              << expected_sub_scores[role] << ',';
                }
            }
            std::cerr << '\n';
        }
        phpr_matches = phpr_matches && call_matches;
    }
    const char score_support_magic[8] = {
        'S', 'C', 'O', 'R', 'E', 'S', 'P', '1'};
    const auto score_support_fixture =
        load_rdp_sectioned_fixture<ScoreSupportCaptureHeader>(
            score_support_fixture_path, score_support_magic);
    bool score_support_matches =
        score_support_fixture.header.next_no == scan_state.next_no &&
        score_support_fixture.header.done_calls == 2 &&
        score_support_fixture.header.group_calls == 3 &&
        score_support_fixture.header.score_calls == 3;
    for (unsigned int call_index = 0; call_index < 2; ++call_index) {
        const int base = static_cast<int>(call_index) * 1000;
        const auto call_sequences =
            rdp_fixture_section<int>(score_support_fixture, base + 2);
        const auto done_before =
            rdp_fixture_section<int>(score_support_fixture, base + 3);
        const auto raw_background =
            rdp_fixture_section<float>(score_support_fixture, base + 4);
        const auto ancestor_background =
            rdp_fixture_section<float>(score_support_fixture, base + 5);
        const auto ancestor_region =
            rdp_fixture_section<float>(score_support_fixture, base + 6);
        const auto expected_result =
            rdp_fixture_section<int>(score_support_fixture, base + 7);
        const auto expected_done =
            rdp_fixture_section<int>(score_support_fixture, base + 8);
        const bool dimensions_match = call_sequences.size() == 3 &&
            expected_result.size() == 1;
        std::vector<int> actual;
        if (dimensions_match) {
            actual = make_rdp_score_filter(
                scan_state.next_no,
                {call_sequences[0], call_sequences[1], call_sequences[2]},
                raw_background, ancestor_background, ancestor_region);
        }
        const bool call_matches = dimensions_match &&
            std::all_of(done_before.begin(), done_before.end(),
                        [](const int value) { return value == 0; }) &&
            expected_result[0] == 1 && actual == expected_done;
        if (!call_matches) {
            std::cerr << "MakeDoneThis3 call " << call_index
                      << " mismatch: dimensions=" << dimensions_match
                      << " output=" << (actual == expected_done) << '\n';
        }
        score_support_matches = score_support_matches && call_matches;
    }
    std::array<RdpTripletGroupState, 3> triplet_groups;
    for (unsigned int call_index = 0; call_index < 3; ++call_index) {
        const int base = 10000 + static_cast<int>(call_index) * 1000;
        const auto raw_header =
            rdp_fixture_section<unsigned int>(score_support_fixture, base + 1);
        const auto call_sequences =
            rdp_fixture_section<int>(score_support_fixture, base + 2);
        const auto call_comparison =
            rdp_fixture_section<int>(score_support_fixture, base + 3);
        const auto counts_before =
            rdp_fixture_section<int>(score_support_fixture, base + 4);
        const auto done_before =
            rdp_fixture_section<int>(score_support_fixture, base + 5);
        const auto groups_before =
            rdp_fixture_section<int>(score_support_fixture, base + 6);
        const auto minimum_before =
            rdp_fixture_section<double>(score_support_fixture, base + 7);
        const auto matrix =
            rdp_fixture_section<float>(score_support_fixture, base + 8);
        const auto expected_result =
            rdp_fixture_section<int>(score_support_fixture, base + 9);
        const auto expected_counts =
            rdp_fixture_section<int>(score_support_fixture, base + 10);
        const auto expected_done =
            rdp_fixture_section<int>(score_support_fixture, base + 11);
        const auto expected_groups =
            rdp_fixture_section<int>(score_support_fixture, base + 12);
        const auto expected_minimum =
            rdp_fixture_section<double>(score_support_fixture, base + 13);
        const int role = raw_header.size() == 4
            ? static_cast<int>(raw_header[2]) : -1;
        const bool dimensions_match = role >= 0 && role < 3 &&
            call_sequences.size() == 3 && call_comparison.size() == 6 &&
            expected_result.size() == 1;
        if (dimensions_match) {
            triplet_groups[role] = make_rdp_triplet_groups(
                role, scan_state.next_no,
                {call_sequences[0], call_sequences[1], call_sequences[2]},
                {call_comparison[0], call_comparison[1], call_comparison[2],
                 call_comparison[3], call_comparison[4], call_comparison[5]},
                matrix,
                {minimum_before[0], minimum_before[1], minimum_before[2]});
        }
        const auto minimum_actual = dimensions_match
            ? as_vector(triplet_groups[role].minimum_distances)
            : std::vector<double>{};
        const bool call_matches = dimensions_match &&
            std::all_of(counts_before.begin(), counts_before.end(),
                        [](const int value) { return value == 0; }) &&
            std::all_of(done_before.begin(), done_before.end(),
                        [](const int value) { return value == 0; }) &&
            std::all_of(groups_before.begin(), groups_before.end(),
                        [](const int value) { return value == 0; }) &&
            triplet_groups[role].result == expected_result[0] &&
            triplet_groups[role].counts == expected_counts &&
            triplet_groups[role].done == expected_done &&
            triplet_groups[role].groups == expected_groups &&
            probability_vectors_match(minimum_actual, expected_minimum);
        if (!call_matches) {
            std::cerr << "MakeTrpGroups2 call " << call_index
                      << " mismatch: role=" << role
                      << " counts=" << (dimensions_match &&
                          triplet_groups[role].counts == expected_counts)
                      << " done=" << (dimensions_match &&
                          triplet_groups[role].done == expected_done)
                      << " groups=" << (dimensions_match &&
                          triplet_groups[role].groups == expected_groups)
                      << " minimum=" << (dimensions_match &&
                          probability_vectors_match(
                              minimum_actual, expected_minimum))
                      << '\n';
        }
        score_support_matches = score_support_matches && call_matches;
    }
    for (unsigned int call_index = 0; call_index < 3; ++call_index) {
        const int base = 20000 + static_cast<int>(call_index) * 1000;
        const auto raw_header =
            rdp_fixture_section<unsigned int>(score_support_fixture, base + 1);
        const auto call_sequences =
            rdp_fixture_section<int>(score_support_fixture, base + 2);
        auto scores =
            rdp_fixture_section<double>(score_support_fixture, base + 3);
        const auto counts =
            rdp_fixture_section<int>(score_support_fixture, base + 4);
        const auto groups =
            rdp_fixture_section<int>(score_support_fixture, base + 5);
        const auto first_matrix =
            rdp_fixture_section<float>(score_support_fixture, base + 6);
        const auto second_matrix =
            rdp_fixture_section<float>(score_support_fixture, base + 7);
        const auto expected_result =
            rdp_fixture_section<int>(score_support_fixture, base + 8);
        const auto expected_scores =
            rdp_fixture_section<double>(score_support_fixture, base + 9);
        const int role = raw_header.size() == 4
            ? static_cast<int>(raw_header[2]) : -1;
        const bool dimensions_match = role >= 0 && role < 3 &&
            call_sequences.size() == 3 && scores.size() == 4 &&
            expected_result.size() == 1;
        if (dimensions_match) {
            RdpTripletGroupState inputs = triplet_groups[role];
            inputs.counts = counts;
            inputs.groups = groups;
            scores[role] = make_rdp_triplet_tree_score(
                role, scan_state.next_no,
                {call_sequences[0], call_sequences[1], call_sequences[2]},
                first_matrix, second_matrix, inputs);
        }
        const bool call_matches = dimensions_match &&
            expected_result[0] == 1 &&
            probability_vectors_match(scores, expected_scores);
        if (!call_matches) {
            std::cerr << "MakeTrpScore2 call " << call_index
                      << " mismatch: role=" << role
                      << " output="
                      << (dimensions_match &&
                          probability_vectors_match(scores, expected_scores))
                      << '\n';
        }
        score_support_matches = score_support_matches && call_matches;
    }
    std::cout << "RDP all-redo raw event scan: " << events.scanned_triplets
              << " triplets (Alist return " << redo_count << "), "
              << events.significant_candidates << " significant intervals, "
              << raw_total << " stored candidates; native MakeTestPVs boundary "
              << native_total << " candidates, " << matching_rows << '/'
              << (scan_state.next_no + 1) << " row counts equal, "
              << matching_event_identity << '/' << compared_events
              << " event identities exact, " << matching_event_probability
              << '/' << compared_events << " probabilities exact; MakeTestPVs "
              << (make_test_structure_matches ? "PASS" : "FAIL")
              << ", first selection "
              << (first_selection_matches ? "PASS" : "FAIL") << " ("
              << trace[0] << ',' << trace[1] << "), UFDist "
              << (ufdist_matches ? "PASS" : "FAIL")
              << ", region distance "
              << (region_distance_matches ? "PASS" : "FAIL")
              << ", CheckMatrix "
              << (check_matrix_matches ? "PASS" : "FAIL")
              << ", first NJ tree "
              << (make_nj_matches ? "PASS" : "FAIL")
              << ", MakeSDMP2 "
              << (make_sdmp_matches ? "PASS" : "FAIL")
              << ", FillRmat "
              << (fill_rmat_matches ? "PASS" : "FAIL")
              << ", CalCR " << (calcr_matches ? "PASS" : "FAIL")
              << ", MakeRList "
              << (make_rlist_matches ? "PASS" : "FAIL")
              << ", FindActualEvents "
              << (find_actual_matches ? "PASS" : "FAIL")
              << ", StripDupInv "
              << (strip_dup_matches ? "PASS" : "FAIL")
              << ", MakeRCompat "
              << (rcompat_matches ? "PASS" : "FAIL")
              << ", compatibility flow "
              << (rcompat_flow_matches ? "PASS" : "FAIL")
              << ", MakePhPrScore "
              << (phpr_matches ? "PASS" : "FAIL")
              << ", score support "
              << (score_support_matches ? "PASS" : "FAIL") << "\n";

    return make_test_structure_matches && first_selection_matches &&
            ufdist_matches && region_distance_matches &&
            check_matrix_matches && make_nj_matches ?
        (make_sdmp_matches && fill_rmat_matches && calcr_matches &&
         make_rlist_matches && find_actual_matches && strip_dup_matches &&
             rcompat_matches && rcompat_flow_matches && phpr_matches &&
             score_support_matches
             ? 0 : 1) : 1;
}

int alist_rdp4_fixture(const std::string& path) {
    return run_alist_rdp4_fixture(
        &MathFuncs::MyMathFuncs::AlistRDP4, path, std::cout, std::cerr);
}

int find_subseq_pb3_fixture(const std::string& path) {
    return run_find_subseq_pb3_fixture(
        &MathFuncs::MyMathFuncs::FindSubSeqPB3, path, std::cout, std::cerr);
}

int find_subseq_pb4_fixture(const std::string& path) {
    return run_find_subseq_pb4_fixture(
        &MathFuncs::MyMathFuncs::FindSubSeqPB4, path, std::cout, std::cerr);
}

int xohomology_p_fixture(const std::string& path) {
    return run_xohomology_p_fixture(
        &MathFuncs::MyMathFuncs::XOHomologyP, path, std::cout, std::cerr);
}

int find_next_p_fixture(const std::string& path) {
    return run_find_next_p_fixture(
        &MathFuncs::MyMathFuncs::FindNextP, path, std::cout, std::cerr);
}

int define_event_p2_fixture(const std::string& path) {
    return run_define_event_p2_fixture(
        &MathFuncs::MyMathFuncs::DefineEventP2, path, std::cout, std::cerr);
}

int prob_calc_p2_fixture(const std::string& path) {
    return run_prob_calc_p2_fixture(
        &MathFuncs::MyMathFuncs::ProbCalcP2, path, std::cout, std::cerr);
}

int find_first_co_p_fixture(const std::string& path) {
    return run_find_first_co_p_fixture(
        &MathFuncs::MyMathFuncs::FindFirstCOP, path, std::cout, std::cerr);
}

int prob_calc_p_fixture(const std::string& path) {
    return run_prob_calc_p_fixture(
        &MathFuncs::MyMathFuncs::ProbCalcP, path, std::cout, std::cerr);
}

int clean_xosnw_fixture(const std::string& path) {
    return run_clean_xosnw_fixture(
        &MathFuncs::MyMathFuncs::CleanXOSNW, path, std::cout, std::cerr);
}

int make_test_pvs_fixture(const std::string& path) {
    return run_make_test_pvs_fixture<XOVERDEFINE>(
        &MathFuncs::MyMathFuncs::MakeTestPVs, path, std::cout, std::cerr);
}

int find_best_rec_signal_p2_fixture(const std::string& path) {
    return run_find_best_rec_signal_p2_fixture(
        &MathFuncs::MyMathFuncs::FindBestRecSignalP2, path, std::cout,
        std::cerr);
}

int ufdist_fixture(const std::string& path) {
    return run_ufdist_fixture(
        &MathFuncs::MyMathFuncs::UFDist, path, std::cout, std::cerr);
}

int super_dist_p2_fixture(const std::string& path) {
    return run_super_dist_p2_fixture(
        &MathFuncs::MyMathFuncs::SuperDistP2, path, std::cout, std::cerr);
}

int super_dist_p_fixture(const std::string& path) {
    return run_super_dist_p_fixture(
        &MathFuncs::MyMathFuncs::SuperDistP, path, std::cout, std::cerr);
}

int check_matrix_p_fixture(const std::string& path) {
    return run_check_matrix_p_fixture(
        &MathFuncs::MyMathFuncs::CheckMatrixP, path, std::cout, std::cerr);
}

int make_nj_trees_p2_fixture(const std::string& path) {
    return run_make_nj_trees_p2_fixture(
        &MathFuncs::MyMathFuncs::MakeNJTreesP2, path, std::cout, std::cerr);
}

int mark_outsides_fixture(const std::string& path) {
    return run_mark_outsides_fixture<XOVERDEFINE>(
        &MathFuncs::MyMathFuncs::MarkOutsides, path, std::cout, std::cerr);
}

int make_sdmp2_fixture(const std::string& path) {
    return run_make_sdmp2_fixture(
        &MathFuncs::MyMathFuncs::MakeSDMP2, path, std::cout, std::cerr);
}

int fill_rmat_fixture(const std::string& path) {
    return run_fill_rmat_fixture(
        &MathFuncs::MyMathFuncs::FillRmat, path, std::cout, std::cerr);
}

int oracle_fixture_chain(const std::string& directory) {
    const auto fixture = [&directory](std::string_view name) {
        return directory + "/" + std::string(name);
    };
    if (super_dist_p_fixture(fixture("super-dist-p-v1.bin")) != 0) return 1;
    if (alist_rdp4_fixture(fixture("alist-rdp4-v1.bin")) != 0) return 1;
    if (find_subseq_pb3_fixture(fixture("find-subseq-pb3-v1.bin")) != 0) {
        return 1;
    }
    if (xohomology_p_fixture(fixture("xohomology-p-v1.bin")) != 0) return 1;
    if (find_next_p_fixture(fixture("find-next-p-v1.bin")) != 0) return 1;
    if (define_event_p2_fixture(fixture("define-event-p2-v1.bin")) != 0) {
        return 1;
    }
    if (prob_calc_p2_fixture(fixture("prob-calc-p2-v1.bin")) != 0) return 1;
    if (find_subseq_pb4_fixture(fixture("find-subseq-pb4-v1.bin")) != 0) {
        return 1;
    }
    if (find_first_co_p_fixture(fixture("find-first-co-p-v1.bin")) != 0) {
        return 1;
    }
    if (prob_calc_p_fixture(fixture("prob-calc-p-v1.bin")) != 0) return 1;
    if (clean_xosnw_fixture(fixture("clean-xosnw-v1.bin")) != 0) return 1;
    if (make_test_pvs_fixture(fixture("make-test-pvs-v1.bin")) != 0) return 1;
    if (find_best_rec_signal_p2_fixture(
            fixture("find-best-rec-signal-p2-v1.bin")) != 0) {
        return 1;
    }
    if (ufdist_fixture(fixture("ufdist-v1.bin")) != 0) return 1;
    if (super_dist_p2_fixture(fixture("super-dist-p2-v1.bin")) != 0) return 1;
    if (check_matrix_p_fixture(fixture("check-matrix-p-v1.bin")) != 0) return 1;
    if (make_nj_trees_p2_fixture(fixture("make-nj-trees-p2-v1.bin")) != 0) {
        return 1;
    }
    if (mark_outsides_fixture(fixture("mark-outsides-v1.bin")) != 0) return 1;
    if (make_sdmp2_fixture(fixture("make-sdmp2-v1.bin")) != 0) return 1;
    if (fill_rmat_fixture(fixture("fill-rmat-y0-v1.bin")) != 0) return 1;
    if (fill_rmat_fixture(fixture("fill-rmat-y1-v1.bin")) != 0) return 1;
    if (fill_rmat_fixture(fixture("fill-rmat-y2-v1.bin")) != 0) return 1;
    std::cout << "oracle fixture chain: PASS (22 live calls in RDP order)\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "self-test") {
        return self_test();
    }
    if (argc == 2 && std::string_view(argv[1]) == "preprocess-fixture") {
        return preprocess_fixture();
    }
    if (argc == 2 && std::string_view(argv[1]) == "distance-fixture") {
        return distance_fixture();
    }
    if (argc == 4 &&
        std::string_view(argv[1]) == "fasta-preprocess-fixture") {
        return fasta_preprocess_fixture(argv[2], argv[3]);
    }
    if (argc == 4 &&
        std::string_view(argv[1]) == "fasta-distance-fixture") {
        return fasta_distance_fixture(argv[2], argv[3]);
    }
    if (argc == 4 &&
        std::string_view(argv[1]) == "fasta-tree-distance-fixture") {
        return fasta_tree_distance_fixture(argv[2], argv[3]);
    }
    if (argc == 4 &&
        std::string_view(argv[1]) == "fasta-alist-rdp4-fixture") {
        return fasta_alist_rdp4_fixture(argv[2], argv[3]);
    }
    if (argc == 5 &&
        std::string_view(argv[1]) == "fasta-first-xover-fixture") {
        return fasta_first_xover_fixture(argv[2], argv[3], argv[4]);
    }
    if (argc == 12 &&
        std::string_view(argv[1]) == "fasta-first-xover-walk-fixture") {
        return fasta_first_xover_walk_fixture(
            argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8],
            argv[9], argv[10], argv[11]);
    }
    if (argc == 22 &&
        std::string_view(argv[1]) == "fasta-all-redo-events-fixture") {
        return fasta_all_redo_events_fixture(
            argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8],
            argv[9], argv[10], argv[11], {argv[12], argv[13], argv[14]},
            argv[15], argv[16], argv[17], argv[18], argv[19], argv[20],
            argv[21]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "alist-rdp4-fixture") {
        return alist_rdp4_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "find-subseq-pb3-fixture") {
        return find_subseq_pb3_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "find-subseq-pb4-fixture") {
        return find_subseq_pb4_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "xohomology-p-fixture") {
        return xohomology_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "find-next-p-fixture") {
        return find_next_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "define-event-p2-fixture") {
        return define_event_p2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "prob-calc-p2-fixture") {
        return prob_calc_p2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "find-first-co-p-fixture") {
        return find_first_co_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "prob-calc-p-fixture") {
        return prob_calc_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "clean-xosnw-fixture") {
        return clean_xosnw_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "make-test-pvs-fixture") {
        return make_test_pvs_fixture(argv[2]);
    }
    if (argc == 3 &&
        std::string_view(argv[1]) == "find-best-rec-signal-p2-fixture") {
        return find_best_rec_signal_p2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "oracle-fixture-chain") {
        return oracle_fixture_chain(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "ufdist-fixture") {
        return ufdist_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "super-dist-p2-fixture") {
        return super_dist_p2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "super-dist-p-fixture") {
        return super_dist_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "check-matrix-p-fixture") {
        return check_matrix_p_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "make-nj-trees-p2-fixture") {
        return make_nj_trees_p2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "mark-outsides-fixture") {
        return mark_outsides_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "make-sdmp2-fixture") {
        return make_sdmp2_fixture(argv[2]);
    }
    if (argc == 3 && std::string_view(argv[1]) == "fill-rmat-fixture") {
        return fill_rmat_fixture(argv[2]);
    }

    std::cerr
        << "usage: rdp-core <self-test|distance-fixture|preprocess-fixture>\n"
        << "       rdp-core fasta-preprocess-fixture <alignment.fasta> <alist-capture.bin>\n"
        << "       rdp-core fasta-distance-fixture <alignment.fasta> <alist-capture.bin>\n"
        << "       rdp-core fasta-tree-distance-fixture <alignment.fasta> <alist-capture.bin>\n"
        << "       rdp-core fasta-alist-rdp4-fixture <alignment.fasta> <alist-capture.bin>\n"
        << "       rdp-core fasta-first-xover-fixture <alignment.fasta> <alist-capture.bin> <find-subseq-pb3-capture.bin>\n"
        << "       rdp-core fasta-first-xover-walk-fixture <alignment.fasta> <alist-capture.bin> <find-subseq-pb3-capture.bin> <xohomology-capture.bin> <find-next-capture.bin> <define-event-capture.bin> <prob-calc-p2-capture.bin> <prob-calc-p-capture.bin> <find-subseq-pb4-capture.bin> <clean-xosnw-capture.bin>\n"
        << "       rdp-core alist-rdp4-fixture <capture.bin>\n"
        << "       rdp-core find-subseq-pb3-fixture <capture.bin>\n"
        << "       rdp-core find-subseq-pb4-fixture <capture.bin>\n"
        << "       rdp-core xohomology-p-fixture <capture.bin>\n"
        << "       rdp-core find-next-p-fixture <capture.bin>\n"
        << "       rdp-core define-event-p2-fixture <capture.bin>\n"
        << "       rdp-core prob-calc-p2-fixture <capture.bin>\n"
        << "       rdp-core find-first-co-p-fixture <capture.bin>\n"
        << "       rdp-core prob-calc-p-fixture <capture.bin>\n"
        << "       rdp-core clean-xosnw-fixture <capture.bin>\n"
        << "       rdp-core make-test-pvs-fixture <capture.bin>\n"
        << "       rdp-core find-best-rec-signal-p2-fixture <capture.bin>\n"
        << "       rdp-core ufdist-fixture <capture.bin>\n"
        << "       rdp-core super-dist-p2-fixture <capture.bin>\n"
        << "       rdp-core super-dist-p-fixture <capture.bin>\n"
        << "       rdp-core check-matrix-p-fixture <capture.bin>\n"
        << "       rdp-core make-nj-trees-p2-fixture <capture.bin>\n"
        << "       rdp-core mark-outsides-fixture <capture.bin>\n"
        << "       rdp-core make-sdmp2-fixture <capture.bin>\n"
        << "       rdp-core fill-rmat-fixture <capture.bin>\n"
        << "       rdp-core oracle-fixture-chain <fixture-directory>\n";
    return 2;
}

#include "MathFuncsDll.h"
#include "alist_rdp4_fixture.hpp"
#include "distance_fixture.hpp"
#include "distance_state.hpp"
#include "event_state_fixture.hpp"
#include "find_subseq_pb3_fixture.hpp"
#include "identification_fixture.hpp"
#include "preprocess_fixture.hpp"
#include "rdp_walk_fixture.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <string_view>
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
    const std::string& find_best_fixture_path) {
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
              << trace[0] << ',' << trace[1] << ")\n";
    return make_test_structure_matches && first_selection_matches ? 0 : 1;
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
    if (argc == 7 &&
        std::string_view(argv[1]) == "fasta-all-redo-events-fixture") {
        return fasta_all_redo_events_fixture(
            argv[2], argv[3], argv[4], argv[5], argv[6]);
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

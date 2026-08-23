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

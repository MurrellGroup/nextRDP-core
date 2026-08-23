#include "MathFuncsDll.h"
#include "alist_rdp4_fixture.hpp"
#include "distance_fixture.hpp"
#include "event_state_fixture.hpp"
#include "find_subseq_pb3_fixture.hpp"
#include "identification_fixture.hpp"
#include "preprocess_fixture.hpp"
#include "rdp_walk_fixture.hpp"

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
    std::cout << "oracle fixture chain: PASS (21 live calls in RDP order)\n";
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
        << "       rdp-core check-matrix-p-fixture <capture.bin>\n"
        << "       rdp-core make-nj-trees-p2-fixture <capture.bin>\n"
        << "       rdp-core mark-outsides-fixture <capture.bin>\n"
        << "       rdp-core make-sdmp2-fixture <capture.bin>\n"
        << "       rdp-core fill-rmat-fixture <capture.bin>\n"
        << "       rdp-core oracle-fixture-chain <fixture-directory>\n";
    return 2;
}

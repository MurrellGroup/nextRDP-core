#pragma once

#include "rdp_walk_fixture.hpp"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

#pragma pack(push, 1)
struct UFDistCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t sequence_length;
    std::int32_t begin;
    std::int32_t end;
    std::int32_t pair_matrix_ub;
    std::int32_t sequence_data_ub;
};

struct SuperDistP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t x;
    std::int32_t next_no;
    std::int32_t ub14;
    std::int32_t ub04;
    std::int32_t ub13;
    std::int32_t ub03;
    std::int32_t ub12;
    std::int32_t ub02;
    std::int32_t ub11;
};

struct CheckMatrixPCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::int32_t sco;
    std::int32_t minimum_sequence_size;
    std::int32_t missing_pair_ub;
    std::int32_t valid_ub;
    std::int32_t sub_valid_ub;
    std::int32_t matrix_ub;
};

struct MakeNJTreesP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t resolve_root;
    std::int32_t nseqs;
    std::int32_t next_no;
    std::int32_t seed;
    std::int32_t name_length;
    std::int32_t sequence_length;
    std::int32_t trace_sequences_ub;
    std::int32_t first_matrix_ub;
    std::int32_t second_matrix_ub;
    std::int32_t first_adjusted_matrix_ub;
    std::int32_t second_adjusted_matrix_ub;
};

struct MarkOutsidesCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t done_sequence_ub;
    std::int32_t next_no;
    std::int32_t xover_rows_ub;
    std::int32_t maximum_current_xover;
    std::int32_t xover_struct_bytes;
};

struct MakeSDMP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::int32_t sequence_length;
};

struct FillRmatCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t y;
    std::int32_t next_no;
    std::int32_t result_matrix_ub1;
    std::int32_t result_matrix_ub2;
    std::int32_t distance_matrix_ub1;
    std::int32_t distance_matrix_ub2;
    std::int32_t distance_matrix_ub3;
};

struct CalCRChainCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
};

struct MakeRListCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
};

struct FindActualEventsCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t sequence_length;
    std::int32_t next_no;
    std::int32_t event_list_ub;
};

struct StripDupInvCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
};

struct RCompatCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::uint32_t calls;
};

struct PhPrScoreCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::uint32_t calls;
};

struct ScoreSupportCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::uint32_t done_calls;
    std::uint32_t group_calls;
    std::uint32_t score_calls;
};

struct CheckPatternCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::int32_t sequence_length;
};

struct CMaxD2P3CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::int32_t sequence_length;
};

struct ConsensusCsvCaptureHeader {
    char magic[8];
    std::uint32_t version;
};

struct CollectEventsCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t next_no;
    std::int32_t sequence_length;
};
#pragma pack(pop)

static_assert(sizeof(UFDistCaptureHeader) == 32);
static_assert(sizeof(SuperDistP2CaptureHeader) == 48);
static_assert(sizeof(CheckMatrixPCaptureHeader) == 40);
static_assert(sizeof(MakeNJTreesP2CaptureHeader) == 56);
static_assert(sizeof(MarkOutsidesCaptureHeader) == 32);
static_assert(sizeof(MakeSDMP2CaptureHeader) == 20);
static_assert(sizeof(FillRmatCaptureHeader) == 40);
static_assert(sizeof(CalCRChainCaptureHeader) == 16);
static_assert(sizeof(MakeRListCaptureHeader) == 16);
static_assert(sizeof(FindActualEventsCaptureHeader) == 24);
static_assert(sizeof(StripDupInvCaptureHeader) == 16);
static_assert(sizeof(RCompatCaptureHeader) == 20);
static_assert(sizeof(PhPrScoreCaptureHeader) == 20);
static_assert(sizeof(ScoreSupportCaptureHeader) == 28);
static_assert(sizeof(CheckPatternCaptureHeader) == 20);
static_assert(sizeof(ConsensusCsvCaptureHeader) == 12);
static_assert(sizeof(CollectEventsCaptureHeader) == 20);

template <typename UFDistFn>
inline int run_ufdist_fixture(
    UFDistFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'U', 'F', 'D', 'I', 'S', 'T', '\0', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<UFDistCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto valid = rdp_fixture_section<float>(fixture, 1);
    auto differences = rdp_fixture_section<float>(fixture, 2);
    auto breakpoint_distance = rdp_fixture_section<float>(fixture, 3);
    auto remainder_distance = rdp_fixture_section<float>(fixture, 4);
    auto sequences = rdp_fixture_section<int>(fixture, 5);
    auto sequence_data = rdp_fixture_section<short>(fixture, 6);
    const auto expected_breakpoint = rdp_fixture_section<float>(fixture, 101);
    const auto expected_remainder = rdp_fixture_section<float>(fixture, 102);
    const auto expected_result = rdp_fixture_section<int>(fixture, 103);

    const int result = function(
        h.sequence_length, h.begin, h.end, h.pair_matrix_ub, valid.data(),
        differences.data(), breakpoint_distance.data(), remainder_distance.data(),
        sequences.data(), h.sequence_data_ub, sequence_data.data());
    const bool breakpoint_matches =
        breakpoint_distance.size() == expected_breakpoint.size() &&
        std::memcmp(breakpoint_distance.data(), expected_breakpoint.data(),
                    breakpoint_distance.size() * sizeof(float)) == 0;
    const bool remainder_matches =
        remainder_distance.size() == expected_remainder.size() &&
        std::memcmp(remainder_distance.data(), expected_remainder.data(),
                    remainder_distance.size() * sizeof(float)) == 0;
    if (expected_result.size() != 1 || result != expected_result[0] ||
        !breakpoint_matches || !remainder_matches) {
        error << "UFDist parity: FAIL\n";
        return 1;
    }
    output << "UFDist parity: PASS (event interval " << h.begin << '-' << h.end
           << ")\n";
    return 0;
}

template <typename SuperDistFn>
inline int run_super_dist_fixture(
    SuperDistFn function, const std::string& path, const char (&magic)[8],
    const char* routine_name, std::ostream& output, std::ostream& error) {
    const auto fixture =
        load_rdp_sectioned_fixture<SuperDistP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto average = rdp_fixture_section<double>(fixture, 1);
    auto pair_diff = rdp_fixture_section<float>(fixture, 2);
    auto pair_valid = rdp_fixture_section<float>(fixture, 3);
    auto distance = rdp_fixture_section<float>(fixture, 4);
    auto redo = rdp_fixture_section<short>(fixture, 5);
    auto category_count = rdp_fixture_section<int>(fixture, 6);
    auto iseq14 = rdp_fixture_section<short>(fixture, 7);
    auto iseq04 = rdp_fixture_section<short>(fixture, 8);
    auto iseq13 = rdp_fixture_section<short>(fixture, 9);
    auto iseq03 = rdp_fixture_section<short>(fixture, 10);
    auto iseq12 = rdp_fixture_section<short>(fixture, 11);
    auto iseq02 = rdp_fixture_section<short>(fixture, 12);
    auto iseq11 = rdp_fixture_section<short>(fixture, 13);
    auto valid14 = rdp_fixture_section<char>(fixture, 14);
    auto diff14 = rdp_fixture_section<char>(fixture, 15);
    auto valid13 = rdp_fixture_section<char>(fixture, 16);
    auto diff13 = rdp_fixture_section<char>(fixture, 17);
    auto valid12 = rdp_fixture_section<char>(fixture, 18);
    auto diff12 = rdp_fixture_section<char>(fixture, 19);
    auto valid11 = rdp_fixture_section<char>(fixture, 20);
    auto diff11 = rdp_fixture_section<char>(fixture, 21);
    auto diff04 = rdp_fixture_section<char>(fixture, 22);
    auto diff03 = rdp_fixture_section<char>(fixture, 23);
    auto diff02 = rdp_fixture_section<char>(fixture, 24);
    const auto expected_average = rdp_fixture_section<double>(fixture, 101);
    const auto expected_pair_diff = rdp_fixture_section<float>(fixture, 102);
    const auto expected_pair_valid = rdp_fixture_section<float>(fixture, 103);
    const auto expected_distance = rdp_fixture_section<float>(fixture, 104);
    const auto expected_result = rdp_fixture_section<double>(fixture, 105);
    if (average.size() != 1) {
        throw std::runtime_error(
            std::string(routine_name) + " average section must be scalar");
    }
    const double result = function(
        h.x, h.next_no, h.ub14, h.ub04, h.ub13, h.ub03, h.ub12, h.ub02,
        h.ub11, average.data(), pair_diff.data(), pair_valid.data(),
        distance.data(), redo.data(), category_count.data(), iseq14.data(),
        iseq04.data(), iseq13.data(), iseq03.data(), iseq12.data(),
        iseq02.data(), iseq11.data(), valid14.data(), diff14.data(),
        valid13.data(), diff13.data(), valid12.data(), diff12.data(),
        valid11.data(), diff11.data(), diff04.data(), diff03.data(),
        diff02.data());
    const auto bytes_equal = [](const auto& actual, const auto& expected) {
        using Value = typename std::decay_t<decltype(actual)>::value_type;
        return actual.size() == expected.size() &&
            std::memcmp(actual.data(), expected.data(),
                        actual.size() * sizeof(Value)) == 0;
    };
    const bool matches = expected_result.size() == 1 &&
        std::memcmp(&result, expected_result.data(), sizeof(result)) == 0 &&
        bytes_equal(average, expected_average) &&
        bytes_equal(pair_diff, expected_pair_diff) &&
        bytes_equal(pair_valid, expected_pair_valid) &&
        bytes_equal(distance, expected_distance);
    if (!matches) {
        error << routine_name << " parity: FAIL\n";
        return 1;
    }
    output << routine_name << " parity: PASS (" << (h.next_no + 1)
           << " sequences)\n";
    return 0;
}

template <typename SuperDistPFn>
inline int run_super_dist_p_fixture(
    SuperDistPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'S', 'U', 'P', 'D', 'I', 'S', 'T', '1'};
    return run_super_dist_fixture(
        function, path, magic, "SuperDistP", output, error);
}

template <typename SuperDistP2Fn>
inline int run_super_dist_p2_fixture(
    SuperDistP2Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'S', 'U', 'P', 'D', 'I', 'S', 'T', '2'};
    return run_super_dist_fixture(
        function, path, magic, "SuperDistP2", output, error);
}

template <typename CheckMatrixPFn>
inline int run_check_matrix_p_fixture(
    CheckMatrixPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'C', 'H', 'K', 'M', 'A', 'T', 'P', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<CheckMatrixPCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto minimums = rdp_fixture_section<int>(fixture, 1);
    auto sequences = rdp_fixture_section<int>(fixture, 2);
    auto missing_pair = rdp_fixture_section<unsigned char>(fixture, 3);
    auto valid = rdp_fixture_section<float>(fixture, 4);
    auto sub_valid = rdp_fixture_section<float>(fixture, 5);
    auto first_matrix = rdp_fixture_section<float>(fixture, 6);
    auto second_matrix = rdp_fixture_section<float>(fixture, 7);
    auto first_total = rdp_fixture_section<int>(fixture, 8);
    auto second_total = rdp_fixture_section<int>(fixture, 9);
    const auto expected_minimums = rdp_fixture_section<int>(fixture, 101);
    const auto expected_missing =
        rdp_fixture_section<unsigned char>(fixture, 102);
    const auto expected_first_matrix = rdp_fixture_section<float>(fixture, 103);
    const auto expected_second_matrix = rdp_fixture_section<float>(fixture, 104);
    const auto expected_first_total = rdp_fixture_section<int>(fixture, 105);
    const auto expected_second_total = rdp_fixture_section<int>(fixture, 106);
    const auto expected_result = rdp_fixture_section<int>(fixture, 107);
    const int result = function(
        minimums.data(), sequences.data(), h.next_no, h.sco,
        h.minimum_sequence_size, h.missing_pair_ub, missing_pair.data(),
        h.valid_ub, valid.data(), h.sub_valid_ub, sub_valid.data(), h.matrix_ub,
        first_matrix.data(), second_matrix.data(), first_total.data(),
        second_total.data());
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && minimums == expected_minimums &&
        missing_pair == expected_missing && first_matrix == expected_first_matrix &&
        second_matrix == expected_second_matrix && first_total == expected_first_total &&
        second_total == expected_second_total;
    if (!matches) {
        error << "CheckMatrixP parity: FAIL\n";
        return 1;
    }
    output << "CheckMatrixP parity: PASS (" << (h.next_no + 1)
           << " sequences)\n";
    return 0;
}

template <typename MakeNJTreesP2Fn>
inline int run_make_nj_trees_p2_fixture(
    MakeNJTreesP2Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'M', 'A', 'K', 'E', 'N', 'J', 'P', '2'};
    const auto fixture =
        load_rdp_sectioned_fixture<MakeNJTreesP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto sequences = rdp_fixture_section<int>(fixture, 1);
    auto min_pair = rdp_fixture_section<unsigned char>(fixture, 2);
    auto sequence_pair = rdp_fixture_section<unsigned char>(fixture, 3);
    auto outlier = rdp_fixture_section<int>(fixture, 4);
    auto trace_sequences = rdp_fixture_section<int>(fixture, 5);
    auto first_matrix = rdp_fixture_section<float>(fixture, 6);
    auto second_matrix = rdp_fixture_section<float>(fixture, 7);
    auto first_adjusted = rdp_fixture_section<float>(fixture, 8);
    auto second_adjusted = rdp_fixture_section<float>(fixture, 9);
    auto redo_list = rdp_fixture_section<int>(fixture, 10);
    auto first_holder = rdp_fixture_section<char>(fixture, 11);
    auto second_holder = rdp_fixture_section<char>(fixture, 12);
    auto temp_first = rdp_fixture_section<float>(fixture, 13);
    auto temp_second = rdp_fixture_section<float>(fixture, 14);
    const int result = function(
        h.resolve_root, h.nseqs, h.next_no, sequences.data(), min_pair.data(),
        sequence_pair.data(), h.seed, h.name_length, h.sequence_length,
        h.trace_sequences_ub, outlier.data(), trace_sequences.data(),
        h.first_matrix_ub, first_matrix.data(), h.second_matrix_ub,
        second_matrix.data(), h.first_adjusted_matrix_ub, first_adjusted.data(),
        h.second_adjusted_matrix_ub, second_adjusted.data(), redo_list.data(),
        first_holder.data(), second_holder.data(), temp_first.data(),
        temp_second.data());
    const auto bytes_equal = [](const auto& actual, const auto& expected) {
        using Value = typename std::decay_t<decltype(actual)>::value_type;
        return actual.size() == expected.size() &&
            std::memcmp(actual.data(), expected.data(),
                        actual.size() * sizeof(Value)) == 0;
    };
    bool matches = true;
    const auto check = [&](const char* name, const auto& actual,
                           const auto& expected) {
        if (!bytes_equal(actual, expected)) {
            std::size_t first = 0;
            const auto common = std::min(actual.size(), expected.size());
            while (first < common &&
                   std::memcmp(&actual[first], &expected[first],
                               sizeof(actual[first])) == 0) {
                ++first;
            }
            error << "MakeNJTreesP2 mismatch: " << name
                  << " (first element " << first;
            if (first < common) {
                error << ", actual bytes";
                const auto* actual_bytes = reinterpret_cast<const unsigned char*>(
                    &actual[first]);
                const auto* expected_bytes = reinterpret_cast<const unsigned char*>(
                    &expected[first]);
                for (std::size_t byte = 0; byte < sizeof(actual[first]); ++byte) {
                    error << ' ' << static_cast<unsigned int>(actual_bytes[byte]);
                }
                error << ", expected bytes";
                for (std::size_t byte = 0; byte < sizeof(expected[first]); ++byte) {
                    error << ' ' << static_cast<unsigned int>(expected_bytes[byte]);
                }
            }
            error << ")\n";
            matches = false;
        }
    };
    if (result != rdp_fixture_section<int>(fixture, 115).at(0)) {
        error << "MakeNJTreesP2 mismatch: result\n";
        matches = false;
    }
    check("sequences", sequences, rdp_fixture_section<int>(fixture, 101));
    check("min_pair", min_pair,
          rdp_fixture_section<unsigned char>(fixture, 102));
    check("sequence_pair", sequence_pair,
          rdp_fixture_section<unsigned char>(fixture, 103));
    check("outlier", outlier, rdp_fixture_section<int>(fixture, 104));
    check("trace_sequences", trace_sequences,
          rdp_fixture_section<int>(fixture, 105));
    check("first_matrix", first_matrix,
          rdp_fixture_section<float>(fixture, 106));
    check("second_matrix", second_matrix,
          rdp_fixture_section<float>(fixture, 107));
    check("first_adjusted", first_adjusted,
          rdp_fixture_section<float>(fixture, 108));
    check("second_adjusted", second_adjusted,
          rdp_fixture_section<float>(fixture, 109));
    check("redo_list", redo_list, rdp_fixture_section<int>(fixture, 110));
    check("first_holder", first_holder,
          rdp_fixture_section<char>(fixture, 111));
    check("second_holder", second_holder,
          rdp_fixture_section<char>(fixture, 112));
    check("temp_first", temp_first,
          rdp_fixture_section<float>(fixture, 113));
    check("temp_second", temp_second,
          rdp_fixture_section<float>(fixture, 114));
    if (!matches) {
        error << "MakeNJTreesP2 parity: FAIL\n";
        return 1;
    }
    output << "MakeNJTreesP2 parity: PASS (" << (h.nseqs + 1)
           << " local sequences)\n";
    return 0;
}

template <typename XOverType, typename MarkOutsidesFn>
inline int run_mark_outsides_fixture(
    MarkOutsidesFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'M', 'A', 'R', 'K', 'O', 'U', 'T', 'S'};
    const auto fixture =
        load_rdp_sectioned_fixture<MarkOutsidesCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    if (h.xover_struct_bytes != static_cast<int>(sizeof(XOverType))) {
        throw std::runtime_error("XOVERDEFINE layout differs from Windows oracle");
    }
    auto done_sequence = rdp_fixture_section<unsigned char>(fixture, 1);
    auto current_xover = rdp_fixture_section<short>(fixture, 2);
    auto xover_list = rdp_fixture_section<XOverType>(fixture, 3);
    const int result = function(
        h.done_sequence_ub, done_sequence.data(), h.next_no,
        h.xover_rows_ub, current_xover.data(), xover_list.data());
    const auto expected_done =
        rdp_fixture_section<unsigned char>(fixture, 101);
    const auto expected_current = rdp_fixture_section<short>(fixture, 102);
    const auto expected_xovers = rdp_fixture_section<XOverType>(fixture, 103);
    const auto expected_result = rdp_fixture_section<int>(fixture, 104);
    const bool xovers_match = xover_list.size() == expected_xovers.size() &&
        std::memcmp(xover_list.data(), expected_xovers.data(),
                    xover_list.size() * sizeof(XOverType)) == 0;
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && done_sequence == expected_done &&
        current_xover == expected_current && xovers_match;
    if (!matches) {
        error << "MarkOutsides parity: FAIL\n";
        return 1;
    }
    output << "MarkOutsides parity: PASS (" << (h.next_no + 1)
           << " sequence rows, " << h.maximum_current_xover
           << " event slots)\n";
    return 0;
}

template <typename MakeSDMP2Fn>
inline int run_make_sdmp2_fixture(
    MakeSDMP2Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'M', 'A', 'K', 'E', 'S', 'D', 'M', 'P'};
    const auto fixture =
        load_rdp_sectioned_fixture<MakeSDMP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto start_positions = rdp_fixture_section<int>(fixture, 1);
    auto end_positions = rdp_fixture_section<int>(fixture, 2);
    auto sequences = rdp_fixture_section<int>(fixture, 3);
    auto comparison_matrix = rdp_fixture_section<int>(fixture, 4);
    auto missing_data = rdp_fixture_section<unsigned char>(fixture, 5);
    auto sequence_data = rdp_fixture_section<short>(fixture, 6);
    auto summary_matrix = rdp_fixture_section<double>(fixture, 7);
    auto distance_matrix = rdp_fixture_section<double>(fixture, 8);
    const int result = function(
        h.next_no, h.sequence_length, start_positions.data(),
        end_positions.data(), sequences.data(), comparison_matrix.data(),
        missing_data.data(), sequence_data.data(), summary_matrix.data(),
        distance_matrix.data());
    const auto expected_summary =
        rdp_fixture_section<double>(fixture, 101);
    const auto expected_distance =
        rdp_fixture_section<double>(fixture, 102);
    const auto expected_result = rdp_fixture_section<int>(fixture, 103);
    const bool summary_matches = summary_matrix.size() == expected_summary.size() &&
        std::memcmp(summary_matrix.data(), expected_summary.data(),
                    summary_matrix.size() * sizeof(double)) == 0;
    const bool distance_matches =
        distance_matrix.size() == expected_distance.size() &&
        std::memcmp(distance_matrix.data(), expected_distance.data(),
                    distance_matrix.size() * sizeof(double)) == 0;
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && summary_matches && distance_matches;
    if (!matches) {
        error << "MakeSDMP2 parity: FAIL\n";
        return 1;
    }
    output << "MakeSDMP2 parity: PASS (" << (h.next_no + 1)
           << " sequences, length " << h.sequence_length << ")\n";
    return 0;
}

template <typename FillRmatFn>
inline int run_fill_rmat_fixture(
    FillRmatFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'F', 'I', 'L', 'L', 'R', 'M', 'A', 'T'};
    const auto fixture =
        load_rdp_sectioned_fixture<FillRmatCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto result_matrix = rdp_fixture_section<double>(fixture, 1);
    auto distance_matrix = rdp_fixture_section<double>(fixture, 2);
    auto positions = rdp_fixture_section<unsigned char>(fixture, 3);
    const int result = function(
        h.y, h.next_no, h.result_matrix_ub1, h.result_matrix_ub2,
        h.distance_matrix_ub1, h.distance_matrix_ub2,
        h.distance_matrix_ub3, result_matrix.data(), distance_matrix.data(),
        positions.data());
    const auto expected_matrix = rdp_fixture_section<double>(fixture, 101);
    const auto expected_distance = rdp_fixture_section<double>(fixture, 102);
    const auto expected_positions =
        rdp_fixture_section<unsigned char>(fixture, 103);
    const auto expected_result = rdp_fixture_section<int>(fixture, 104);
    const bool matrix_matches = result_matrix.size() == expected_matrix.size() &&
        std::memcmp(result_matrix.data(), expected_matrix.data(),
                    result_matrix.size() * sizeof(double)) == 0;
    const bool distance_matches =
        distance_matrix.size() == expected_distance.size() &&
        std::memcmp(distance_matrix.data(), expected_distance.data(),
                    distance_matrix.size() * sizeof(double)) == 0;
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && matrix_matches && distance_matches &&
        positions == expected_positions;
    if (!matches) {
        error << "FillRmat parity: FAIL (Y=" << h.y << ")\n";
        return 1;
    }
    output << "FillRmat parity: PASS (Y=" << h.y << ")\n";
    return 0;
}

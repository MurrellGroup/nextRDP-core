#pragma once

#include "rdp_walk_fixture.hpp"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#pragma pack(push, 1)
struct MakeTestPVsCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t done_sequence_ub;
    std::int32_t next_no;
    std::int32_t xover_rows_ub;
    std::int32_t xover_slots_ub;
    std::int32_t xover_struct_bytes;
};

struct FindBestRecSignalP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    char done_target;
    std::int32_t next_no;
    std::int32_t probability_rows_ub;
    std::int32_t probability_columns_ub;
};
#pragma pack(pop)

static_assert(sizeof(MakeTestPVsCaptureHeader) == 32);
static_assert(sizeof(FindBestRecSignalP2CaptureHeader) == 25);

template <typename XOverType, typename MakeTestPVsFn>
inline int run_make_test_pvs_fixture(
    MakeTestPVsFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'M', 'K', 'T', 'E', 'S', 'T', 'P', 'V'};
    const auto fixture =
        load_rdp_sectioned_fixture<MakeTestPVsCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    if (h.xover_struct_bytes != static_cast<int>(sizeof(XOverType))) {
        throw std::runtime_error("XOVERDEFINE layout differs from Windows oracle");
    }
    auto done_sequence = rdp_fixture_section<unsigned char>(fixture, 1);
    auto current_xover = rdp_fixture_section<short>(fixture, 2);
    auto xover_list = rdp_fixture_section<XOverType>(fixture, 3);
    auto test_pvs = rdp_fixture_section<double>(fixture, 4);
    const auto expected_done = rdp_fixture_section<unsigned char>(fixture, 101);
    const auto expected_current = rdp_fixture_section<short>(fixture, 102);
    const auto expected_xovers = rdp_fixture_section<XOverType>(fixture, 103);
    const auto expected_test_pvs = rdp_fixture_section<double>(fixture, 104);
    const auto expected_result = rdp_fixture_section<int>(fixture, 105);

    const int result = function(
        h.done_sequence_ub, done_sequence.data(), h.next_no, h.xover_rows_ub,
        h.xover_slots_ub, current_xover.data(), xover_list.data(),
        test_pvs.data());
    const bool xovers_match = xover_list.size() == expected_xovers.size() &&
        std::memcmp(xover_list.data(), expected_xovers.data(),
                    xover_list.size() * sizeof(XOverType)) == 0;
    const bool test_pvs_match = test_pvs.size() == expected_test_pvs.size() &&
        std::memcmp(test_pvs.data(), expected_test_pvs.data(),
                    test_pvs.size() * sizeof(double)) == 0;
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && done_sequence == expected_done &&
        current_xover == expected_current && xovers_match && test_pvs_match;
    if (!matches) {
        error << "MakeTestPVs parity: FAIL\n";
        return 1;
    }
    output << "MakeTestPVs parity: PASS (" << (h.next_no + 1)
           << " sequence rows)\n";
    return 0;
}

template <typename FindBestRecSignalP2Fn>
inline int run_find_best_rec_signal_p2_fixture(
    FindBestRecSignalP2Fn function, const std::string& path,
    std::ostream& output, std::ostream& error) {
    const char magic[8] = {'F', 'B', 'R', 'S', 'I', 'G', '2', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<FindBestRecSignalP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto low_p = rdp_fixture_section<double>(fixture, 1);
    auto done_sequence = rdp_fixture_section<char>(fixture, 2);
    auto trace = rdp_fixture_section<int>(fixture, 3);
    auto current_xover = rdp_fixture_section<short>(fixture, 4);
    auto test_pvs = rdp_fixture_section<double>(fixture, 5);
    const auto expected_low_p = rdp_fixture_section<double>(fixture, 101);
    const auto expected_done = rdp_fixture_section<char>(fixture, 102);
    const auto expected_trace = rdp_fixture_section<int>(fixture, 103);
    const auto expected_current = rdp_fixture_section<short>(fixture, 104);
    const auto expected_test_pvs = rdp_fixture_section<double>(fixture, 105);
    const auto expected_result = rdp_fixture_section<int>(fixture, 106);
    if (low_p.size() != 1 || trace.size() != 2) {
        throw std::runtime_error("invalid FindBestRecSignalP2 scalar sections");
    }
    const int result = function(
        h.done_target, h.next_no, h.probability_rows_ub,
        h.probability_columns_ub, low_p.data(), done_sequence.data(),
        trace.data(), current_xover.data(), test_pvs.data());
    const bool low_p_matches = expected_low_p.size() == 1 &&
        std::memcmp(low_p.data(), expected_low_p.data(), sizeof(double)) == 0;
    const bool test_pvs_match = test_pvs.size() == expected_test_pvs.size() &&
        std::memcmp(test_pvs.data(), expected_test_pvs.data(),
                    test_pvs.size() * sizeof(double)) == 0;
    const bool matches = expected_result.size() == 1 &&
        result == expected_result[0] && low_p_matches &&
        done_sequence == expected_done && trace == expected_trace &&
        current_xover == expected_current && test_pvs_match;
    if (!matches) {
        error << "FindBestRecSignalP2 parity: FAIL\n";
        return 1;
    }
    output.precision(17);
    output << "FindBestRecSignalP2 parity: PASS (selected " << trace[0] << ','
           << trace[1] << ", p=" << low_p[0] << ")\n";
    return 0;
}

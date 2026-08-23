#pragma once

#include "alist_rdp4_fixture.hpp"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma pack(push, 1)
struct XOHomologyPCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int16_t inlyer;
    std::int32_t sequence_length;
    std::int32_t xover_length;
    std::int16_t xover_window;
};

struct FindNextPCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t homology_ub;
    std::int32_t start;
    std::int32_t high;
    std::int32_t med;
    std::int32_t low;
    std::int32_t xover_length;
    std::int32_t xover_window;
};

struct DefineEventP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t homology_ub;
    std::int32_t short_output;
    std::int32_t long_winded;
    std::int32_t med;
    std::int32_t high;
    std::int32_t low;
    std::int32_t target_x;
    std::int32_t circular;
    std::int32_t xx;
    std::int32_t xover_window;
    std::int32_t sequence_length;
    std::int32_t xover_length;
    std::int32_t sequence_daughter;
    std::int32_t sequence_minor;
};

struct ProbCalcP2CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t fact_three_ub;
    std::int32_t xover_length;
    std::int32_t number_in_common;
    double individual_probability;
    std::int32_t informative_length;
};

struct FindFirstCOPCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t x;
    std::int32_t med;
    std::int32_t high;
    std::int32_t xover_length;
    std::int32_t homology_ub;
};

struct ProbCalcPCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t xover_length;
    std::int32_t number_in_common;
    double individual_probability;
    std::int32_t informative_length;
};

struct CleanXOSNWCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t xover_length;
    std::int32_t xover_window;
    std::int32_t xover_sequence_ub;
};
#pragma pack(pop)

static_assert(sizeof(XOHomologyPCaptureHeader) == 24);
static_assert(sizeof(FindNextPCaptureHeader) == 40);
static_assert(sizeof(DefineEventP2CaptureHeader) == 68);
static_assert(sizeof(ProbCalcP2CaptureHeader) == 36);
static_assert(sizeof(FindFirstCOPCaptureHeader) == 32);
static_assert(sizeof(ProbCalcPCaptureHeader) == 32);
static_assert(sizeof(CleanXOSNWCaptureHeader) == 24);

template <typename Header>
struct RdpSectionedFixture {
    Header header{};
    std::unordered_map<std::uint32_t, std::vector<std::uint8_t>> sections;
};

template <typename Header>
inline RdpSectionedFixture<Header> load_rdp_sectioned_fixture(
    const std::string& path, const char (&magic)[8]) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open fixture: " + path);
    RdpSectionedFixture<Header> fixture;
    input.read(reinterpret_cast<char*>(&fixture.header), sizeof(Header));
    if (!input || std::memcmp(fixture.header.magic, magic, 8) != 0 ||
        fixture.header.version != 1) {
        throw std::runtime_error("invalid routine fixture header: " + path);
    }
    while (true) {
        AlistRdp4SectionHeader section{};
        input.read(reinterpret_cast<char*>(&section), sizeof(section));
        if (!input) throw std::runtime_error("truncated routine fixture");
        if (section.id == static_cast<std::uint32_t>(AlistRdp4Section::end_marker)) {
            if (section.bytes != 0) throw std::runtime_error("invalid end marker");
            break;
        }
        auto& bytes = fixture.sections[section.id];
        if (!bytes.empty()) throw std::runtime_error("duplicate fixture section");
        bytes.resize(section.bytes);
        input.read(reinterpret_cast<char*>(bytes.data()), section.bytes);
        if (!input) throw std::runtime_error("truncated fixture section");
    }
    return fixture;
}

template <typename T, typename Header>
inline std::vector<T> rdp_fixture_section(
    const RdpSectionedFixture<Header>& fixture, std::uint32_t id) {
    const auto found = fixture.sections.find(id);
    if (found == fixture.sections.end() || found->second.size() % sizeof(T) != 0) {
        throw std::runtime_error("missing or mis-sized section " +
                                 std::to_string(id));
    }
    std::vector<T> values(found->second.size() / sizeof(T));
    std::memcpy(values.data(), found->second.data(), found->second.size());
    return values;
}

template <typename XOHomologyPFn>
inline int run_xohomology_p_fixture(
    XOHomologyPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'X', 'O', 'H', 'O', 'M', 'P', '\0', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<XOHomologyPCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto xover_sequence = rdp_fixture_section<char>(fixture, 1);
    auto homology = rdp_fixture_section<int>(fixture, 2);
    const auto expected_homology = rdp_fixture_section<int>(fixture, 101);
    const auto expected_result = rdp_fixture_section<int>(fixture, 102);
    const int result = function(
        h.inlyer, h.sequence_length, h.xover_length, h.xover_window,
        xover_sequence.data(), homology.data());
    if (expected_result.size() != 1 || result != expected_result[0] ||
        homology != expected_homology) {
        error << "XOHomologyP parity: FAIL\n";
        return 1;
    }
    output << "XOHomologyP parity: PASS (start position " << result << ")\n";
    return 0;
}

template <typename FindNextPFn>
inline int run_find_next_p_fixture(
    FindNextPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'F', 'I', 'N', 'D', 'N', 'X', 'T', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<FindNextPCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto homology = rdp_fixture_section<int>(fixture, 1);
    const auto expected_result = rdp_fixture_section<int>(fixture, 101);
    const int result = function(
        h.homology_ub, h.start, h.high, h.med, h.low, h.xover_length,
        h.xover_window, homology.data());
    if (expected_result.size() != 1 || result != expected_result[0]) {
        error << "FindNextP parity: FAIL\n";
        return 1;
    }
    output << "FindNextP parity: PASS (next position " << result << ")\n";
    return 0;
}

template <typename DefineEventP2Fn>
inline int run_define_event_p2_fixture(
    DefineEventP2Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'D', 'E', 'F', 'E', 'V', 'P', '2', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<DefineEventP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto scalar_state = rdp_fixture_section<int>(fixture, 1);
    auto xover_sequence = rdp_fixture_section<char>(fixture, 2);
    auto homology = rdp_fixture_section<int>(fixture, 3);
    const auto expected_scalars = rdp_fixture_section<int>(fixture, 101);
    const auto expected_result = rdp_fixture_section<int>(fixture, 102);
    if (scalar_state.size() != 5) {
        throw std::runtime_error("DefineEventP2 scalar state must have 5 values");
    }
    const int result = function(
        h.homology_ub, h.short_output, h.long_winded, h.med, h.high, h.low,
        h.target_x, h.circular, h.xx, h.xover_window, h.sequence_length,
        h.xover_length, h.sequence_daughter, h.sequence_minor,
        &scalar_state[0], &scalar_state[1], &scalar_state[2], &scalar_state[3],
        &scalar_state[4], xover_sequence.data(), homology.data());
    if (expected_result.size() != 1 || result != expected_result[0] ||
        scalar_state != expected_scalars) {
        error << "DefineEventP2 parity: FAIL\n";
        return 1;
    }
    output << "DefineEventP2 parity: PASS (begin " << scalar_state[1]
           << ", end " << scalar_state[2] << ", next " << result << ")\n";
    return 0;
}

template <typename ProbCalcP2Fn>
inline int run_prob_calc_p2_fixture(
    ProbCalcP2Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'P', 'R', 'O', 'B', 'C', 'P', '2', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<ProbCalcP2CaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto fact_three = rdp_fixture_section<double>(fixture, 1);
    const auto expected_result = rdp_fixture_section<double>(fixture, 101);
    const double result = function(
        fact_three.data(), h.fact_three_ub, h.xover_length,
        h.number_in_common, h.individual_probability, h.informative_length);
    const bool bit_exact = expected_result.size() == 1 &&
                           std::memcmp(&result, expected_result.data(),
                                       sizeof(result)) == 0;
    const double scale = expected_result.empty()
                             ? 1.0
                             : std::max(1.0, std::abs(expected_result[0]));
    const bool numerically_equivalent = expected_result.size() == 1 &&
        std::abs(result - expected_result[0]) <= 1e-14 * scale;
    if (!numerically_equivalent) {
        error.precision(17);
        error << "ProbCalcP2 parity: FAIL (result=" << result << ", expected="
              << (expected_result.empty() ? -1.0 : expected_result[0]) << ")\n";
        return 1;
    }
    output.precision(17);
    output << "ProbCalcP2 parity: PASS (p=" << result;
    if (!bit_exact) output << ", accepted platform-math rounding";
    output << ")\n";
    return 0;
}

template <typename FindFirstCOPFn>
inline int run_find_first_co_p_fixture(
    FindFirstCOPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'F', 'F', 'I', 'R', 'S', 'T', '\0', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<FindFirstCOPCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto homology = rdp_fixture_section<int>(fixture, 1);
    const auto expected_result = rdp_fixture_section<int>(fixture, 101);
    const int result = function(
        h.x, h.med, h.high, h.xover_length, h.homology_ub, homology.data());
    if (expected_result.size() != 1 || result != expected_result[0]) {
        error << "FindFirstCOP parity: FAIL\n";
        return 1;
    }
    output << "FindFirstCOP parity: PASS (position " << result << ")\n";
    return 0;
}

template <typename ProbCalcPFn>
inline int run_prob_calc_p_fixture(
    ProbCalcPFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'P', 'R', 'O', 'B', 'C', 'P', '\0', '\0'};
    const auto fixture =
        load_rdp_sectioned_fixture<ProbCalcPCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto fact = rdp_fixture_section<double>(fixture, 1);
    const auto expected_result = rdp_fixture_section<double>(fixture, 101);
    const double result = function(
        fact.data(), h.xover_length, h.number_in_common,
        h.individual_probability, h.informative_length);
    const bool bit_exact = expected_result.size() == 1 &&
        std::memcmp(&result, expected_result.data(), sizeof(result)) == 0;
    const double scale = expected_result.empty()
                             ? 1.0
                             : std::max(1.0, std::abs(expected_result[0]));
    const bool numerically_equivalent = expected_result.size() == 1 &&
        std::abs(result - expected_result[0]) <= 1e-14 * scale;
    if (!numerically_equivalent) {
        error.precision(17);
        error << "ProbCalcP parity: FAIL (result=" << result << ", expected="
              << (expected_result.empty() ? -1.0 : expected_result[0]) << ")\n";
        return 1;
    }
    output.precision(17);
    output << "ProbCalcP parity: PASS (p=" << result;
    if (!bit_exact) output << ", accepted platform-math rounding";
    output << ")\n";
    return 0;
}

template <typename CleanXOSNWFn>
inline int run_clean_xosnw_fixture(
    CleanXOSNWFn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const char magic[8] = {'C', 'L', 'N', 'X', 'O', 'S', 'N', 'W'};
    const auto fixture =
        load_rdp_sectioned_fixture<CleanXOSNWCaptureHeader>(path, magic);
    const auto& h = fixture.header;
    auto xover_sequence = rdp_fixture_section<char>(fixture, 1);
    const auto expected_xover = rdp_fixture_section<char>(fixture, 101);
    const auto expected_result = rdp_fixture_section<int>(fixture, 102);
    const int result = function(
        h.xover_length, h.xover_window, h.xover_sequence_ub,
        xover_sequence.data());
    if (expected_result.size() != 1 || result != expected_result[0] ||
        xover_sequence != expected_xover) {
        error << "CleanXOSNW parity: FAIL\n";
        return 1;
    }
    output << "CleanXOSNW parity: PASS\n";
    return 0;
}

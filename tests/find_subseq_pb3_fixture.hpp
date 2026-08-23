#pragma once

#include "alist_rdp4_fixture.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma pack(push, 1)
struct FindSubSeqPB3CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t fss_ub;
    std::int32_t xover_window;
    std::int32_t compressed_sequence_ub;
    std::int32_t sequence_length;
    std::int32_t next_no;
    std::int32_t seq1;
    std::int32_t seq2;
    std::int32_t seq3;
    std::int32_t xover_sequence_ub;
};
#pragma pack(pop)

static_assert(sizeof(FindSubSeqPB3CaptureHeader) == 48);

struct FindSubSeqPB3Fixture {
    FindSubSeqPB3CaptureHeader header{};
    std::unordered_map<std::uint32_t, std::vector<std::uint8_t>> sections;
};

inline FindSubSeqPB3Fixture load_find_subseq_pb3_fixture(
    const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open FindSubSeqPB3 fixture: " + path);
    }
    FindSubSeqPB3Fixture fixture;
    input.read(reinterpret_cast<char*>(&fixture.header), sizeof(fixture.header));
    const char expected_magic[8] = {'F', 'S', 'P', 'B', '3', '\0', '\0', '\0'};
    if (!input || std::memcmp(fixture.header.magic, expected_magic, 8) != 0 ||
        fixture.header.version != 1) {
        throw std::runtime_error("invalid FindSubSeqPB3 fixture header");
    }
    while (true) {
        AlistRdp4SectionHeader section{};
        input.read(reinterpret_cast<char*>(&section), sizeof(section));
        if (!input) throw std::runtime_error("truncated FindSubSeqPB3 fixture");
        if (section.id == static_cast<std::uint32_t>(AlistRdp4Section::end_marker)) {
            if (section.bytes != 0) {
                throw std::runtime_error("invalid FindSubSeqPB3 end marker");
            }
            break;
        }
        auto& bytes = fixture.sections[section.id];
        if (!bytes.empty()) {
            throw std::runtime_error("duplicate FindSubSeqPB3 section");
        }
        bytes.resize(section.bytes);
        input.read(reinterpret_cast<char*>(bytes.data()), section.bytes);
        if (!input) {
            throw std::runtime_error("truncated FindSubSeqPB3 section " +
                                     std::to_string(section.id));
        }
    }
    return fixture;
}

template <typename T>
inline std::vector<T> find_subseq_pb3_section(
    const FindSubSeqPB3Fixture& fixture, std::uint32_t id) {
    const auto found = fixture.sections.find(id);
    if (found == fixture.sections.end() || found->second.size() % sizeof(T) != 0) {
        throw std::runtime_error("missing or mis-sized FindSubSeqPB3 section " +
                                 std::to_string(id));
    }
    std::vector<T> values(found->second.size() / sizeof(T));
    std::memcpy(values.data(), found->second.data(), found->second.size());
    return values;
}

template <typename FindSubSeqPB3Fn>
inline int run_find_subseq_pb3_fixture(
    FindSubSeqPB3Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const auto fixture = load_find_subseq_pb3_fixture(path);
    const auto& h = fixture.header;
    auto ah = find_subseq_pb3_section<int>(fixture, 1);
    auto compressed_sequence =
        find_subseq_pb3_section<unsigned char>(fixture, 2);
    auto xover_sequence = find_subseq_pb3_section<char>(fixture, 3);
    auto fss_rdp = find_subseq_pb3_section<unsigned char>(fixture, 4);
    const auto expected_ah = find_subseq_pb3_section<int>(fixture, 101);
    const auto expected_xover = find_subseq_pb3_section<char>(fixture, 102);
    const auto expected_result = find_subseq_pb3_section<int>(fixture, 103);

    const int result = function(
        ah.data(), h.fss_ub, h.xover_window, h.compressed_sequence_ub,
        h.sequence_length, h.next_no, h.seq1, h.seq2, h.seq3,
        compressed_sequence.data(), h.xover_sequence_ub,
        xover_sequence.data(), fss_rdp.data());

    const bool matches = expected_result.size() == 1 &&
                         result == expected_result[0] && ah == expected_ah &&
                         xover_sequence == expected_xover;
    if (!matches) {
        error << "FindSubSeqPB3 parity: FAIL (result=" << result
              << ", expected="
              << (expected_result.empty() ? -1 : expected_result[0]) << ")\n";
        return 1;
    }
    output << "FindSubSeqPB3 parity: PASS (triplet " << h.seq1 << ',' << h.seq2
           << ',' << h.seq3 << ", informative length " << result << ")\n";
    return 0;
}

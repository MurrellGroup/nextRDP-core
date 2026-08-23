#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#pragma pack(push, 1)
struct AlistRdp4CaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::int32_t store_lpv_ub;
    std::int32_t list_length;
    std::int32_t start;
    std::int32_t end;
    std::int32_t next_no;
    double sub_threshold;
    std::int32_t circular;
    std::int32_t mc_correction;
    std::int32_t mc_flag;
    double lowest_probability;
    std::int32_t target_x;
    std::int32_t sequence_length;
    std::int32_t short_output;
    std::int32_t distance_ub;
    std::int32_t tree_distance_ub;
    std::int32_t fss_rdp_ub;
    std::int32_t compressed_sequence_ub;
    std::int32_t xover_window;
    std::int16_t xover_window_x;
    std::int32_t probability_file_flag;
    std::int32_t probability_one_ub;
    std::int32_t probability_two_ub;
    std::int32_t fact_three_ub;
};

struct AlistRdp4SectionHeader {
    std::uint32_t id;
    std::uint32_t bytes;
};
#pragma pack(pop)

static_assert(sizeof(AlistRdp4CaptureHeader) == 110);
static_assert(sizeof(AlistRdp4SectionHeader) == 8);

enum class AlistRdp4Section : std::uint32_t {
    store_lpv_in = 1,
    analysis_list_in = 2,
    redo_list_in = 3,
    distance_in = 4,
    tree_distance_in = 5,
    fss_rdp_in = 6,
    compressed_sequence_in = 7,
    sequence_data_in = 8,
    probability_estimate_in = 9,
    fact_three_in = 10,
    fact_in = 11,
    redo_list_out = 101,
    store_lpv_out = 102,
    result_out = 103,
    end_marker = 0xffffffffU,
};

struct AlistRdp4Fixture {
    AlistRdp4CaptureHeader header{};
    std::unordered_map<std::uint32_t, std::vector<std::uint8_t>> sections;
};

inline AlistRdp4Fixture load_alist_rdp4_fixture(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open AlistRDP4 fixture: " + path);
    }

    AlistRdp4Fixture fixture;
    input.read(reinterpret_cast<char*>(&fixture.header), sizeof(fixture.header));
    const char expected_magic[8] = {'A', 'L', 'R', 'D', 'P', '4', '\0', '\0'};
    if (!input || std::memcmp(fixture.header.magic, expected_magic, 8) != 0 ||
        fixture.header.version != 1) {
        throw std::runtime_error("invalid AlistRDP4 fixture header");
    }

    while (true) {
        AlistRdp4SectionHeader section{};
        input.read(reinterpret_cast<char*>(&section), sizeof(section));
        if (!input) {
            throw std::runtime_error("truncated AlistRDP4 section header");
        }
        if (section.id == static_cast<std::uint32_t>(AlistRdp4Section::end_marker)) {
            if (section.bytes != 0) {
                throw std::runtime_error("invalid AlistRDP4 end marker");
            }
            break;
        }
        auto& bytes = fixture.sections[section.id];
        if (!bytes.empty()) {
            throw std::runtime_error("duplicate AlistRDP4 section");
        }
        bytes.resize(section.bytes);
        input.read(reinterpret_cast<char*>(bytes.data()), section.bytes);
        if (!input) {
            throw std::runtime_error(
                "truncated AlistRDP4 section body " + std::to_string(section.id) +
                " (declared " + std::to_string(section.bytes) + " bytes)");
        }
    }
    if (input.peek() != std::ifstream::traits_type::eof()) {
        throw std::runtime_error("trailing bytes after AlistRDP4 fixture");
    }
    return fixture;
}

template <typename T>
inline std::vector<T> alist_rdp4_typed_section(
    const AlistRdp4Fixture& fixture, AlistRdp4Section section) {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto id = static_cast<std::uint32_t>(section);
    const auto found = fixture.sections.find(id);
    if (found == fixture.sections.end() || found->second.size() % sizeof(T) != 0) {
        throw std::runtime_error("missing or mis-sized AlistRDP4 section " +
                                 std::to_string(id));
    }
    std::vector<T> values(found->second.size() / sizeof(T));
    std::memcpy(values.data(), found->second.data(), found->second.size());
    return values;
}

template <typename AlistRdp4Fn>
inline int run_alist_rdp4_fixture(
    AlistRdp4Fn function, const std::string& path, std::ostream& output,
    std::ostream& error) {
    const auto fixture = load_alist_rdp4_fixture(path);
    const auto& h = fixture.header;
    auto store_lpv = alist_rdp4_typed_section<double>(
        fixture, AlistRdp4Section::store_lpv_in);
    auto analysis_list = alist_rdp4_typed_section<short>(
        fixture, AlistRdp4Section::analysis_list_in);
    auto redo_list = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::redo_list_in);
    const auto redo_list_before = redo_list;
    auto distance = alist_rdp4_typed_section<float>(
        fixture, AlistRdp4Section::distance_in);
    auto tree_distance = alist_rdp4_typed_section<float>(
        fixture, AlistRdp4Section::tree_distance_in);
    auto fss_rdp = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::fss_rdp_in);
    auto compressed_sequence = alist_rdp4_typed_section<unsigned char>(
        fixture, AlistRdp4Section::compressed_sequence_in);
    auto sequence_data = alist_rdp4_typed_section<short>(
        fixture, AlistRdp4Section::sequence_data_in);
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

    const int result = function(
        h.store_lpv_ub, store_lpv.data(), analysis_list.data(), h.list_length,
        h.start, h.end, h.next_no, h.sub_threshold, redo_list.data(),
        h.circular, h.mc_correction, h.mc_flag, h.lowest_probability,
        h.target_x, h.sequence_length, h.short_output, h.distance_ub,
        distance.data(), h.tree_distance_ub, tree_distance.data(),
        h.fss_rdp_ub, h.compressed_sequence_ub, compressed_sequence.data(),
        sequence_data.data(), h.xover_window, h.xover_window_x,
        fss_rdp.data(), h.probability_file_flag, h.probability_one_ub,
        h.probability_two_ub, probability_estimate.data(), h.fact_three_ub,
        fact_three.data(), fact.data());

    bool matches = expected_result.size() == 1 && result == expected_result[0];
    matches = matches && redo_list == expected_redo;
    matches = matches && store_lpv.size() == expected_store.size() &&
              std::memcmp(store_lpv.data(), expected_store.data(),
                          store_lpv.size() * sizeof(double)) == 0;

    std::size_t changed = 0;
    for (std::size_t i = 0; i < redo_list.size(); ++i) {
        if (redo_list[i] != redo_list_before[i]) {
            ++changed;
        }
    }

    if (!matches) {
        error << "AlistRDP4 parity: FAIL (result=" << result
              << ", expected="
              << (expected_result.empty() ? -1 : expected_result[0]) << ")\n";
        return 1;
    }
    output << "AlistRDP4 parity: PASS (" << (h.end - h.start + 1)
           << " triplets, " << changed << " redo states changed)\n";
    return 0;
}

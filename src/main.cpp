#include "MathFuncsDll.h"
#include "alist_rdp4_fixture.hpp"
#include "distance_fixture.hpp"
#include "distance_state.hpp"
#include "event_state_fixture.hpp"
#include "find_subseq_pb3_fixture.hpp"
#include "identification_fixture.hpp"
#include "identification_state.hpp"
#include "mutation_state.hpp"
#include "preprocess_fixture.hpp"
#include "rdp_walk_fixture.hpp"
#include "rescan_schedule.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"
#include "xover_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

#pragma pack(push, 1)
struct ModSeqNumYCaptureHeader {
    char magic[8];
    std::uint32_t version;
    std::uint32_t calls;
    std::uint32_t reserved;
};
using MutationTailCaptureHeader = ModSeqNumYCaptureHeader;
#pragma pack(pop)
static_assert(sizeof(ModSeqNumYCaptureHeader) == 20);

struct AlistRdp3TraceCall {
    int next_no = -1;
    std::vector<std::array<int, 3>> triplets;
    std::vector<unsigned char> redo;
};

struct NativeXoverDefine {
    std::uint8_t outside_flag;
    std::uint8_t misidentify_flag;
    std::uint8_t program_flag;
    std::uint8_t sbp_flag;
    std::uint8_t accept;
    std::int16_t major_parent;
    std::int16_t minor_parent;
    std::int16_t daughter;
    std::int32_t beginning;
    std::int32_t ending;
    std::int32_t length_holder;
    std::int32_t event_number;
    float permutation_pvalue;
    std::int32_t begin_parent;
    std::int32_t end_parent;
    double probability;
    double distance_holder;
};
static_assert(sizeof(NativeXoverDefine) == 56);
static_assert(offsetof(NativeXoverDefine, major_parent) == 6);
static_assert(offsetof(NativeXoverDefine, probability) == 40);

RdpRawEvent to_raw_event(const NativeXoverDefine& source) {
    RdpRawEvent event;
    event.outside_flag = source.outside_flag;
    event.misidentify_flag = source.misidentify_flag;
    event.program_flag = source.program_flag;
    event.sbp_flag = source.sbp_flag;
    event.accept = source.accept;
    event.major_parent = source.major_parent;
    event.minor_parent = source.minor_parent;
    event.daughter = source.daughter;
    event.beginning = source.beginning;
    event.ending = source.ending;
    event.length_holder = source.length_holder;
    event.event_number = source.event_number;
    event.permutation_pvalue = source.permutation_pvalue;
    event.begin_parent = source.begin_parent;
    event.end_parent = source.end_parent;
    event.probability = source.probability;
    event.distance_holder = source.distance_holder;
    return event;
}

bool raw_event_equivalent(const RdpRawEvent& first,
                          const RdpRawEvent& second) {
    const auto close = [](const double actual, const double expected) {
        return std::abs(actual - expected) <= 1e-12 * std::max(
            {1.0, std::abs(actual), std::abs(expected)});
    };
    return first.outside_flag == second.outside_flag &&
        first.misidentify_flag == second.misidentify_flag &&
        first.program_flag == second.program_flag &&
        first.sbp_flag == second.sbp_flag && first.accept == second.accept &&
        first.major_parent == second.major_parent &&
        first.minor_parent == second.minor_parent &&
        first.daughter == second.daughter &&
        first.beginning == second.beginning && first.ending == second.ending &&
        first.length_holder == second.length_holder &&
        first.event_number == second.event_number &&
        std::memcmp(&first.permutation_pvalue, &second.permutation_pvalue,
                    sizeof(float)) == 0 &&
        first.begin_parent == second.begin_parent &&
        first.end_parent == second.end_parent &&
        close(first.probability, second.probability) &&
        close(first.distance_holder, second.distance_holder);
}

std::vector<AlistRdp3TraceCall> load_alist_rdp3_trace(
    const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("cannot open AlistRDP3 trace: " + path);
    const auto length = input.tellg();
    input.seekg(0);
    std::vector<unsigned char> bytes(static_cast<std::size_t>(length));
    input.read(reinterpret_cast<char*>(bytes.data()), length);
    std::size_t offset = 0;
    const auto read_int = [&]() {
        if (offset + sizeof(int) > bytes.size()) {
            throw std::runtime_error("truncated AlistRDP3 trace");
        }
        int value = 0;
        std::memcpy(&value, bytes.data() + offset, sizeof(value));
        offset += sizeof(value);
        return value;
    };
    std::vector<AlistRdp3TraceCall> calls;
    while (offset < bytes.size()) {
        const int magic = read_int();
        const int invocation = read_int();
        const int list_length = read_int();
        const int start = read_int();
        const int end = read_int();
        const int next_no = read_int();
        const int result = read_int();
        const int sequence_length = read_int();
        (void)invocation;
        (void)result;
        (void)sequence_length;
        if (magic != 0x41523343 || start < 0 || end < start ||
            list_length < end) {
            throw std::runtime_error("AlistRDP3 trace header differs");
        }
        AlistRdp3TraceCall call;
        call.next_no = next_no;
        for (int item = start; item <= end; ++item) {
            std::array<int, 3> triplet{};
            for (int role = 0; role < 3; ++role) {
                if (offset + sizeof(short) > bytes.size()) {
                    throw std::runtime_error("truncated AlistRDP3 triplet");
                }
                short value = 0;
                std::memcpy(&value, bytes.data() + offset, sizeof(value));
                offset += sizeof(value);
                triplet[role] = value;
            }
            if (offset >= bytes.size()) {
                throw std::runtime_error("truncated AlistRDP3 redo list");
            }
            call.triplets.push_back(triplet);
            call.redo.push_back(bytes[offset++]);
        }
        calls.push_back(std::move(call));
    }
    return calls;
}

std::vector<unsigned char> load_first_addjust_pairs(
    const std::string& path, int expected_next_no) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open AddjustCXO pair trace");
    std::array<int, 5> header{};
    input.read(reinterpret_cast<char*>(header.data()), sizeof(header));
    if (header[0] != 0x41445052 || header[2] != expected_next_no) {
        throw std::runtime_error("AddjustCXO pair trace header differs");
    }
    std::vector<unsigned char> pairs(
        static_cast<std::size_t>(expected_next_no + 1) *
        (expected_next_no + 1));
    input.read(reinterpret_cast<char*>(pairs.data()), pairs.size());
    if (!input) throw std::runtime_error("truncated AddjustCXO pair trace");
    return pairs;
}

RdpRawEventState load_first_addjust_events(const std::string& pairs_path) {
    std::string path = pairs_path;
    const std::string suffix = "addjust-dopairs.bin";
    const auto suffix_position = path.rfind(suffix);
    if (suffix_position == std::string::npos) {
        throw std::runtime_error("cannot derive AddjustCXO event trace path");
    }
    path.replace(suffix_position, suffix.size(), "addjust-cxo-temp.bin");
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open AddjustCXO event trace");
    std::array<int, 12> header{};
    input.read(reinterpret_cast<char*>(header.data()), sizeof(header));
    if (!input || header[0] != 0x41435854 || header[1] != 1 ||
        header[11] != static_cast<int>(sizeof(NativeXoverDefine))) {
        throw std::runtime_error("AddjustCXO event trace header differs");
    }
    const int next_no = header[2];
    const int trace_ub = header[5];
    const int row_ub = header[6];
    const int slot_ub = header[7];
    std::vector<int> trace_sub(trace_ub + 1);
    input.read(reinterpret_cast<char*>(trace_sub.data()),
               static_cast<std::streamsize>(trace_sub.size() * sizeof(int)));
    std::vector<std::int16_t> current(next_no + 1);
    input.read(reinterpret_cast<char*>(current.data()),
               static_cast<std::streamsize>(current.size() * sizeof(short)));
    std::vector<NativeXoverDefine> matrix(
        static_cast<std::size_t>(row_ub + 1) * (slot_ub + 1));
    input.read(reinterpret_cast<char*>(matrix.data()),
               static_cast<std::streamsize>(
                   matrix.size() * sizeof(NativeXoverDefine)));
    if (!input) throw std::runtime_error("truncated AddjustCXO event trace");
    RdpRawEventState state;
    state.current_xover = current;
    state.xover_list.resize(next_no + 1);
    for (int row = 0; row <= next_no; ++row) {
        for (int slot = 1; slot <= current[row]; ++slot) {
            state.xover_list[row].push_back(to_raw_event(
                matrix[row + slot * (row_ub + 1)]));
        }
    }
    return state;
}

RdpRawEventState load_findbetter_events(
    const std::string& pairs_path, const int wanted_invocation) {
    std::string path = pairs_path;
    const std::string suffix = "addjust-dopairs.bin";
    const auto suffix_position = path.rfind(suffix);
    if (suffix_position == std::string::npos) {
        throw std::runtime_error("cannot derive MakeTestPVs trace path");
    }
    path.replace(suffix_position, suffix.size(), "findbetter-pxolist.bin");
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open MakeTestPVs trace");
    while (input) {
        std::array<int, 8> header{};
        input.read(reinterpret_cast<char*>(header.data()), sizeof(header));
        if (!input) break;
        if (header[0] != 0x4658504c ||
            header[7] != static_cast<int>(sizeof(NativeXoverDefine))) {
            throw std::runtime_error("MakeTestPVs trace header differs");
        }
        const int invocation = header[1];
        const int done_row_ub = header[2];
        const int next_no = header[3];
        const int xover_row_ub = header[4];
        const int xover_slot_ub = header[5];
        std::vector<std::int16_t> current(next_no + 1);
        input.read(reinterpret_cast<char*>(current.data()),
                   static_cast<std::streamsize>(
                       current.size() * sizeof(std::int16_t)));
        input.seekg(
            static_cast<std::streamoff>(done_row_ub + 1) *
                (xover_slot_ub + 1),
            std::ios::cur);
        std::vector<NativeXoverDefine> matrix(
            static_cast<std::size_t>(xover_row_ub + 1) *
            (xover_slot_ub + 1));
        input.read(reinterpret_cast<char*>(matrix.data()),
                   static_cast<std::streamsize>(
                       matrix.size() * sizeof(NativeXoverDefine)));
        input.seekg(
            static_cast<std::streamoff>(next_no + 1) *
                (xover_slot_ub + 1) * sizeof(double),
            std::ios::cur);
        if (!input) throw std::runtime_error("truncated MakeTestPVs trace");
        if (invocation != wanted_invocation) continue;
        RdpRawEventState state;
        state.current_xover = current;
        state.xover_list.resize(next_no + 1);
        for (int row = 0; row <= next_no; ++row) {
            for (int slot = 1; slot <= current[row]; ++slot) {
                state.xover_list[row].push_back(to_raw_event(
                    matrix[row + slot * (xover_row_ub + 1)]));
            }
        }
        return state;
    }
    throw std::runtime_error("requested MakeTestPVs invocation is absent");
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

std::pair<bool, std::size_t> compare_rdp_event_states(
    const RdpRawEventState& actual, const RdpRawEventState& expected) {
    bool matches = actual.current_xover == expected.current_xover &&
        actual.xover_list.size() == expected.xover_list.size();
    std::size_t mismatches = 0;
    const auto rows = std::min(
        actual.xover_list.size(), expected.xover_list.size());
    for (std::size_t row = 0; row < rows; ++row) {
        const auto& actual_row = actual.xover_list[row];
        const auto& expected_row = expected.xover_list[row];
        if (actual_row.size() != expected_row.size()) {
            matches = false;
            mismatches += actual_row.size() > expected_row.size()
                ? actual_row.size() - expected_row.size()
                : expected_row.size() - actual_row.size();
        }
        const auto slots = std::min(actual_row.size(), expected_row.size());
        for (std::size_t slot = 0; slot < slots; ++slot) {
            if (!raw_event_equivalent(actual_row[slot], expected_row[slot])) {
                matches = false;
                ++mismatches;
            }
        }
    }
    return {matches, mismatches};
}

std::vector<short> flatten_rdp_triplets(
    const std::vector<std::array<int, 3>>& triplets) {
    std::vector<short> output;
    output.reserve(triplets.size() * 3);
    for (const auto& triplet : triplets) {
        for (const int sequence : triplet) {
            output.push_back(static_cast<short>(sequence));
        }
    }
    return output;
}

std::vector<unsigned char> screen_rdp_rescan_triplets(
    const std::vector<std::array<int, 3>>& triplets,
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const AlistRdp4CaptureHeader& settings,
    std::vector<unsigned char>& fss_rdp,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact) {
    if (triplets.empty()) return {};
    auto analysis_list = flatten_rdp_triplets(triplets);
    std::vector<unsigned char> redo(triplets.size(), 0);
    const double uncorrected_threshold = settings.mc_flag == 0
        ? settings.lowest_probability / settings.mc_correction
        : settings.lowest_probability;
    MathFuncs::MyMathFuncs::AlistRDP3(
        analysis_list.data(), static_cast<int>(triplets.size()) - 1, 0,
        static_cast<int>(triplets.size()) - 1, scan_state.next_no,
        uncorrected_threshold, redo.data(), settings.circular,
        settings.mc_correction, settings.mc_flag, settings.lowest_probability,
        settings.target_x, scan_state.sequence_length, settings.short_output,
        scan_state.next_no, const_cast<float*>(distance_state.distance.data()),
        scan_state.next_no,
        const_cast<float*>(tree_state.tree_distance.data()),
        settings.fss_rdp_ub, scan_state.compressed_sequence_ub,
        const_cast<unsigned char*>(scan_state.compressed_sequence.data()),
        const_cast<short*>(scan_state.sequence_data.data()),
        settings.xover_window, settings.xover_window_x, fss_rdp.data(),
        settings.probability_file_flag, settings.probability_one_ub,
        settings.probability_two_ub, probability_estimate.data(),
        settings.fact_three_ub, fact_three.data(), fact.data());
    return redo;
}

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
    const std::string& score_support_fixture_path,
    const std::string& check_pattern_fixture_path,
    const std::string& cmaxd_fixture_path,
    const std::string& consensus_fixture_path,
    const std::string& collect_events_fixture_path,
    const std::string& modseqnumy_fixture_path,
    const std::string& mutation_tail_fixture_path,
    const std::string& alist_rdp3_trace_path,
    const std::string& addjust_pairs_trace_path) {
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
            const bool integrated_matrix_matches = call_index < 3
                ? matrix == background_adjusted
                : (call_index < 6 ? matrix == region_adjusted : true);
            integrated_matches =
                integrated.role == metadata[1] &&
                integrated_matrix_matches &&
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
        std::cerr << "RDP compatibility flow mismatch: generated "
                  << tree_compatibility_flow.calls.size() << " of "
                  << rcompat_fixture.header.calls << " captured calls\n";
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
    std::array<RdpPhylProScoreState, 3> phpr_states;
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
            phpr_states[call_index] = actual;
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
    std::array<std::vector<int>, 2> score_filters;
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
            score_filters[call_index] = actual;
        }
        const bool call_matches = dimensions_match &&
            (call_index != 0 ||
             (raw_background == first_direct_small &&
              ancestor_background == first_direct_small &&
              ancestor_region == second_direct_small)) &&
            std::all_of(done_before.begin(), done_before.end(),
                        [](const int value) { return value == 0; }) &&
            expected_result[0] == 1 && actual == expected_done;
        if (!call_matches) {
            std::cerr << "MakeDoneThis3 call " << call_index
                      << " mismatch: dimensions=" << dimensions_match
                      << " output=" << (actual == expected_done)
                      << " matrices="
                      << (raw_background == first_direct_small)
                      << (ancestor_background == first_direct_small)
                      << (ancestor_region == second_direct_small) << '\n';
        }
        score_support_matches = score_support_matches && call_matches;
    }
    std::array<RdpTripletGroupState, 3> triplet_groups;
    std::array<double, 3> triplet_scores{};
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
            matrix == first_adjusted_small &&
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
            triplet_scores[role] = scores[role];
        }
        const bool call_matches = dimensions_match &&
            first_matrix == first_adjusted_small &&
            second_matrix == second_adjusted_small &&
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
    const char check_pattern_magic[8] = {
        'C', 'H', 'K', 'P', 'A', 'T', '1', '\0'};
    const auto check_pattern_fixture =
        load_rdp_sectioned_fixture<CheckPatternCaptureHeader>(
            check_pattern_fixture_path, check_pattern_magic);
    const std::array<int, 3> pattern_starts{
        actual_starts[0], actual_starts[2], actual_starts[4]};
    const std::array<int, 3> pattern_ends{
        actual_ends[1], actual_ends[3], actual_ends[4]};
    const auto expected_pattern_before =
        rdp_fixture_section<double>(check_pattern_fixture, 5);
    const auto expected_done_before =
        rdp_fixture_section<unsigned char>(check_pattern_fixture, 6);
    const auto expected_pattern =
        rdp_fixture_section<double>(check_pattern_fixture, 102);
    const auto expected_pattern_done =
        rdp_fixture_section<unsigned char>(check_pattern_fixture, 103);
    const auto expected_pattern_sequences =
        rdp_fixture_section<short>(check_pattern_fixture, 4);
    const bool pattern_sequences_match =
        scan_state.sequence_data.size() >= expected_pattern_sequences.size() &&
        std::equal(expected_pattern_sequences.begin(),
                   expected_pattern_sequences.end(),
                   scan_state.sequence_data.begin());
    const auto pattern_state = check_rdp_sequence_patterns(
        scan_state.sequence_length, scan_state.next_no,
        correlation_sequences, pattern_starts, pattern_ends,
        correlation_comparison, scan_state.sequence_data,
        actual_resolution.acceptable_sequences);
    const bool check_pattern_matches =
        check_pattern_fixture.header.next_no == scan_state.next_no &&
        check_pattern_fixture.header.sequence_length ==
            scan_state.sequence_length &&
        as_vector(correlation_sequences) ==
            rdp_fixture_section<int>(check_pattern_fixture, 1) &&
        as_vector(pattern_starts) ==
            rdp_fixture_section<int>(check_pattern_fixture, 2) &&
        as_vector(pattern_ends) ==
            rdp_fixture_section<int>(check_pattern_fixture, 3) &&
        pattern_sequences_match &&
        std::all_of(expected_pattern_before.begin(),
                    expected_pattern_before.end(),
                    [](const double value) { return value == 0.0; }) &&
        std::all_of(expected_done_before.begin(), expected_done_before.end(),
                    [](const unsigned char value) { return value == 0; }) &&
        rdp_fixture_section<int>(check_pattern_fixture, 101) ==
            std::vector<int>{1} &&
        pattern_state.pattern == expected_pattern &&
        pattern_state.done == expected_pattern_done;
    if (!check_pattern_matches) {
        std::cerr << "CheckPatternX mismatch: inputs="
                  << pattern_sequences_match
                  << " pattern="
                  << (pattern_state.pattern == expected_pattern)
                  << " done="
                  << (pattern_state.done == expected_pattern_done) << '\n';
    }
    bool final_trim_prefix_matches = true;
    bool final_trim_runs = false;
    RdpFinalTrimState downstream_candidates;
    downstream_candidates.candidate_last = actual_resolution.candidates.last;
    downstream_candidates.candidate_list = actual_resolution.candidates.list;
    downstream_candidates.acceptable_sequences =
        pattern_state.acceptable_sequences;
    std::array<int, 3> expected_trim_last{};
    std::vector<int> expected_trim_list;
    for (unsigned int call_index = 6;
         call_index < rcompat_fixture.header.calls; ++call_index) {
        const int base = static_cast<int>(call_index) * 1000;
        const auto call_last =
            rdp_fixture_section<int>(rcompat_fixture, base + 7);
        const auto call_list =
            rdp_fixture_section<int>(rcompat_fixture, base + 11);
        if (call_last != as_vector(actual_resolution.candidates.last) ||
            call_list != actual_resolution.candidates.list) {
            if (call_last.size() == 3) {
                expected_trim_last = {
                    call_last[0], call_last[1], call_last[2]};
                expected_trim_list = call_list;
                final_trim_runs = true;
            }
            break;
        }
    }
    if (final_trim_runs) {
        auto trim_prefix = run_rdp_final_trim_candidate_maintenance(
            scan_state.next_no, correlation_sequences,
            correlation_comparison, final_minimum_pair, role_lists.inside,
            correlation_decisions.warnings, actual_resolution.unfound,
            actual_resolution.correlations.correlations.correlation,
            actual_resolution.correlations.correlations.inversion,
            local_distance_panels, first_adjusted_small,
            second_adjusted_small, first_collapsed, second_collapsed,
            actual_resolution.candidates.last,
            actual_resolution.candidates.list,
            pattern_state.acceptable_sequences);
        const auto trim_maintenance = trim_prefix;
        trim_prefix.acceptable_sequences = calculate_rdp_match_evidence(
            scan_state.sequence_length, scan_state.next_no,
            selected.beginning, selected.ending, correlation_sequences,
            correlation_comparison, scan_state.sequence_data,
            trim_prefix.acceptable_sequences, true);
        trim_prefix = make_rdp_consensus_candidates(
            scan_state.next_no, correlation_sequences,
            correlation_comparison, correlation_decisions.warnings,
            actual_resolution.correlations.correlations.correlation,
            actual_resolution.correlations.correlations.inversion,
            generated_matrices.background, generated_matrices.event_region,
            background_adjusted, region_adjusted, first_direct_small,
            second_direct_small, first_adjusted_small, second_adjusted_small,
            first_collapsed, second_collapsed, std::move(trim_prefix), true);
        downstream_candidates = trim_prefix;
        final_trim_prefix_matches =
            trim_prefix.candidate_last == expected_trim_last;
        for (int role = 0; role < 3 && final_trim_prefix_matches; ++role) {
            for (int slot = 0; slot <= expected_trim_last[role]; ++slot) {
                if (trim_prefix.candidate_list[role + slot * 3] !=
                    expected_trim_list[role + slot * 3]) {
                    final_trim_prefix_matches = false;
                    break;
                }
            }
        }
        if (!final_trim_prefix_matches) {
            std::cerr << "FinalTrim/ConsensusOK checkpoint: event="
                      << selected.beginning << '-' << selected.ending
                      << " maintenance=";
            for (const int value : trim_maintenance.candidate_last) {
                std::cerr << value << ',';
            }
            std::cerr << ':';
            for (int role = 0; role < 3; ++role) {
                std::cerr << '[';
                for (int slot = 0;
                     slot <= trim_maintenance.candidate_last[role]; ++slot) {
                    std::cerr << trim_maintenance.candidate_list[
                        role + slot * 3] << ',';
                }
                std::cerr << ']';
            }
            std::cerr << " last=";
            for (const int value : trim_prefix.candidate_last) {
                std::cerr << value << ',';
            }
            std::cerr << '/';
            for (const int value : expected_trim_last) {
                std::cerr << value << ',';
            }
            std::cerr << " lists=";
            for (int role = 0; role < 3; ++role) {
                std::cerr << '[';
                for (int slot = 0; slot <= trim_prefix.candidate_last[role];
                     ++slot) {
                    std::cerr << trim_prefix.candidate_list[
                        role + slot * 3] << ',';
                }
                std::cerr << "]/[";
                for (int slot = 0; slot <= expected_trim_last[role]; ++slot) {
                    std::cerr << expected_trim_list[role + slot * 3] << ',';
                }
                std::cerr << "] ";
            }
            for (int role = 0; role < 3; ++role) {
                for (int slot = 0; slot <= expected_trim_last[role]; ++slot) {
                    const int sequence = expected_trim_list[role + slot * 3];
                    bool present = false;
                    for (int actual_slot = 0;
                         actual_slot <= trim_prefix.candidate_last[role];
                         ++actual_slot) {
                        present = present || trim_prefix.candidate_list[
                            role + actual_slot * 3] == sequence;
                    }
                    if (present) continue;
                    std::cerr << " missing(" << role << ',' << sequence
                              << ") ok=";
                    for (int category = 0; category <= 18; ++category) {
                        std::cerr << trim_prefix.acceptable_sequences[
                            role + category * 3 + sequence * 57] << ',';
                    }
                    std::cerr << " F/FA/FC="
                              << first_direct_small[role + sequence * 3]
                              << '/'
                              << first_adjusted_small[role + sequence * 3]
                              << '/' << first_collapsed[role + sequence * 3]
                              << " S/SA/SC="
                              << second_direct_small[role + sequence * 3]
                              << '/'
                              << second_adjusted_small[role + sequence * 3]
                              << '/' << second_collapsed[role + sequence * 3]
                              << " rc=";
                    for (int region = 0; region < 3; ++region) {
                        std::cerr << actual_resolution.correlations.correlations
                            .correlation[role + region * 3 + sequence * 9]
                                  << ',';
                    }
                    const auto full_index = [matrix_stride](
                        const int first, const int second) {
                        return static_cast<std::size_t>(first) +
                            static_cast<std::size_t>(second) * matrix_stride;
                    };
                    std::cerr << " stragglers=";
                    for (int actual_slot = 0;
                         actual_slot <= trim_prefix.candidate_last[role];
                         ++actual_slot) {
                        const int candidate = trim_prefix.candidate_list[
                            role + actual_slot * 3];
                        std::cerr << candidate << '{'
                                  << second_adjusted_small[
                                         role + candidate * 3]
                                  << ">=" << region_adjusted[
                                         full_index(candidate, sequence)]
                                  << ',' << first_adjusted_small[
                                         role + candidate * 3]
                                  << ">=" << background_adjusted[
                                         full_index(candidate, sequence)]
                                  << ',' << second_direct_small[
                                         role + candidate * 3]
                                  << ">=" << generated_matrices.event_region[
                                         full_index(candidate, sequence)]
                                  << ',' << first_direct_small[
                                         role + candidate * 3]
                                  << ">=" << generated_matrices.background[
                                         full_index(candidate, sequence)]
                                  << "},";
                    }
                }
            }
            std::cerr << '\n';
        }
    }
    const char cmaxd_magic[8] = {
        'C', 'M', 'A', 'X', 'D', '3', 'V', '1'};
    const auto cmaxd_fixture =
        load_rdp_sectioned_fixture<CMaxD2P3CaptureHeader>(
            cmaxd_fixture_path, cmaxd_magic);
    const auto cmaxd_metadata =
        rdp_fixture_section<unsigned int>(cmaxd_fixture, 1);
    const auto cmaxd_expected_sequence =
        rdp_fixture_section<short>(cmaxd_fixture, 2);
    const std::vector<unsigned char> maximum_distance_mask(
        static_cast<std::size_t>(scan_state.next_no + 1), 0);
    const auto maximum_distance = calculate_rdp_maximum_distances(
        scan_state.sequence_length, scan_state.next_no,
        correlation_sequences, selected.beginning, selected.ending,
        scan_state.sequence_data, maximum_distance_mask);
    const bool cmaxd_sequence_matches =
        scan_state.sequence_data.size() >= cmaxd_expected_sequence.size() &&
        std::equal(cmaxd_expected_sequence.begin(),
                   cmaxd_expected_sequence.end(),
                   scan_state.sequence_data.begin());
    const bool cmaxd_matches =
        cmaxd_fixture.header.next_no == scan_state.next_no &&
        cmaxd_fixture.header.sequence_length == scan_state.sequence_length &&
        cmaxd_metadata.size() == 10 &&
        cmaxd_metadata[2] ==
            static_cast<unsigned int>(maximum_distance.included_last) &&
        std::equal(correlation_sequences.begin(), correlation_sequences.end(),
                   cmaxd_metadata.begin() + 3) &&
        cmaxd_metadata[6] == static_cast<unsigned int>(selected.beginning) &&
        cmaxd_metadata[7] == static_cast<unsigned int>(selected.ending) &&
        cmaxd_sequence_matches &&
        maximum_distance.informative_to_position ==
            rdp_fixture_section<int>(cmaxd_fixture, 3) &&
        maximum_distance.position_to_informative ==
            rdp_fixture_section<int>(cmaxd_fixture, 4) &&
        maximum_distance.nucleotide_map ==
            rdp_fixture_section<unsigned char>(cmaxd_fixture, 5) &&
        maximum_distance.included_sequences ==
            rdp_fixture_section<int>(cmaxd_fixture, 6) &&
        maximum_distance.included_mask ==
            rdp_fixture_section<unsigned char>(cmaxd_fixture, 7) &&
        maximum_distance.representative_mask ==
            rdp_fixture_section<unsigned char>(cmaxd_fixture, 8) &&
        maximum_distance.split_scores ==
            rdp_fixture_section<float>(cmaxd_fixture, 9) &&
        std::vector<int>{maximum_distance.result} ==
            rdp_fixture_section<int>(cmaxd_fixture, 10) &&
        as_vector(maximum_distance.distance_totals) ==
            rdp_fixture_section<float>(cmaxd_fixture, 11) &&
        as_vector(maximum_distance.distance_counts) ==
            rdp_fixture_section<int>(cmaxd_fixture, 12);
    if (!cmaxd_matches) {
        std::cerr << "CalcMaxD/CMaxD2P3 mismatch: sequence="
                  << cmaxd_sequence_matches << " map="
                  << (maximum_distance.informative_to_position ==
                      rdp_fixture_section<int>(cmaxd_fixture, 3))
                  << '/' << (maximum_distance.position_to_informative ==
                      rdp_fixture_section<int>(cmaxd_fixture, 4))
                  << " included="
                  << (maximum_distance.included_sequences ==
                      rdp_fixture_section<int>(cmaxd_fixture, 6))
                  << " totals=";
        const auto expected_totals =
            rdp_fixture_section<float>(cmaxd_fixture, 11);
        for (int role = 0; role < 3; ++role) {
            std::cerr << maximum_distance.distance_totals[role] << '/'
                      << expected_totals[role] << ',';
        }
        std::cerr << '\n';
    }
    const char consensus_magic[8] = {
        'C', 'O', 'N', 'S', 'C', 'V', '1', '\0'};
    const auto consensus_fixture =
        load_rdp_sectioned_fixture<ConsensusCsvCaptureHeader>(
            consensus_fixture_path, consensus_magic);
    const auto consensus_values =
        rdp_fixture_section<double>(consensus_fixture, 1);
    if (consensus_values.size() != 120) {
        throw std::runtime_error("consensus fixture must contain 120 values");
    }
    const auto csv_value = [&](const int metric, const int role) {
        return consensus_values[metric * 3 + role];
    };
    RdpConsensusInputs consensus_inputs;
    consensus_inputs.next_no = scan_state.next_no;
    consensus_inputs.permanent_next_no = scan_state.next_no;
    consensus_inputs.comparison_matrix = correlation_comparison;
    for (int role = 0; role < 3; ++role) {
        consensus_inputs.list_correlation[role] = csv_value(0, role);
        consensus_inputs.simple_distance_strength[role] = csv_value(1, role);
        consensus_inputs.simple_distance_score[role] =
            static_cast<int>(csv_value(2, role));
        consensus_inputs.phylpro[role] = csv_value(3, role);
        consensus_inputs.phylpro_secondary[role] = csv_value(4, role);
        consensus_inputs.phylpro_collapsed[role] = csv_value(5, role);
        consensus_inputs.subtree_score[role] = csv_value(6, role);
        consensus_inputs.split_distance[role] = csv_value(7, role);
        consensus_inputs.outlier_index[role] =
            static_cast<int>(csv_value(8, role));
        consensus_inputs.subtree_phylpro[role] = csv_value(9, role);
        consensus_inputs.subtree_score_secondary[role] = csv_value(10, role);
        consensus_inputs.subtree_phylpro_secondary[role] = csv_value(11, role);
        consensus_inputs.compatibility[role] =
            static_cast<int>(csv_value(14, role));
        consensus_inputs.compatibility_secondary[role] =
            static_cast<int>(csv_value(15, role));
        consensus_inputs.compatibility_tertiary[role] =
            static_cast<int>(csv_value(16, role));
        consensus_inputs.compatibility_quaternary[role] =
            static_cast<int>(csv_value(17, role));
        consensus_inputs.region_compatibility[role] =
            static_cast<int>(csv_value(18, role));
        consensus_inputs.region_compatibility_secondary[role] =
            static_cast<int>(csv_value(19, role));
        consensus_inputs.region_compatibility_tertiary[role] =
            static_cast<int>(csv_value(20, role));
        consensus_inputs.region_compatibility_quaternary[role] =
            static_cast<int>(csv_value(21, role));
        consensus_inputs.post_trim_compatibility[role] =
            static_cast<int>(csv_value(24, role));
        consensus_inputs.post_trim_region_compatibility[role] =
            static_cast<int>(csv_value(25, role));
        consensus_inputs.triplet_score[role] = csv_value(26, role);
        consensus_inputs.bad_distances[role] = csv_value(27, role);
        consensus_inputs.outside_list[role] =
            static_cast<int>(csv_value(28, role));
        consensus_inputs.list_correlation_secondary[role] =
            csv_value(29, role);
        consensus_inputs.list_correlation_tertiary[role] =
            csv_value(30, role);
        consensus_inputs.outlier_check[role] =
            static_cast<int>(csv_value(34, role));
        consensus_inputs.ranks[role] = {
            static_cast<int>(csv_value(37, role)),
            static_cast<int>(csv_value(38, role))};
        consensus_inputs.maximum_distance[role] =
            static_cast<float>(csv_value(39, role));
    }
    const auto split_distances = calculate_rdp_split_distances(
        scan_state.next_no, correlation_sequences, role_lists.inside,
        first_adjusted_small, generated_matrices.background,
        generated_matrices.event_region, score_filters[0]);
    const auto simple_distances = calculate_rdp_simple_distances(
        scan_state.next_no, correlation_sequences, role_lists.inside,
        actual_resolution.candidates.last,
        actual_resolution.candidates.list, generated_matrices.background,
        generated_matrices.event_region);
    std::array<int, 3> outlier_checks{};
    if (final_minimum_pair[0] != final_minimum_pair[1]) {
        outlier_checks = calculate_rdp_outlier_checks(
            scan_state.next_no, correlation_sequences, role_lists.inside,
            first_adjusted_small, second_adjusted_small);
    }
    const auto bad_distances = calculate_rdp_bad_distances(
        scan_state.next_no, correlation_sequences, correlation_comparison,
        actual_resolution.candidates.last,
        actual_resolution.candidates.list, actual_resolution.unfound,
        actual_resolution.correlations.correlations.correlation,
        first_adjusted_small, local_distance_panels);
    RdpListCorrelationState list_correlations;
    if (final_minimum_pair[0] != final_minimum_pair[1]) {
        list_correlations = calculate_rdp_list_correlations(
            scan_state.next_no, correlation_sequences, role_lists.inside,
            correlation_decisions.warnings,
            candidate_lists.last, candidate_lists.list,
            correlation_decisions.correlations.inversion,
            correlation_decisions.correlations.tested_correlation,
            first_adjusted_small, second_adjusted_small);
    }
    std::array<int, 3> post_trim_background{};
    std::array<int, 3> post_trim_region{};
    if (final_trim_runs) {
        const auto run_post_trim_compatibility = [&](
            const std::vector<float>& matrix,
            const std::array<double, 3>& list_distances,
            std::array<int, 3>& compatibility) {
            std::array<int, 3> reverse{};
            for (int role = 0; role < 3; ++role) {
                std::array<int, 3> nonrecombinant_last{};
                std::vector<int> nonrecombinant_list(
                    static_cast<std::size_t>(3) *
                        (scan_state.next_no + 1), 0);
                make_rdp_tree_compatibility_call(
                    scan_state.next_no, correlation_sequences,
                    correlation_comparison, role,
                    actual_resolution.inversion_penalty,
                    downstream_candidates.candidate_last,
                    downstream_candidates.candidate_list, good_comparisons,
                    matrix, list_distances, compatibility, reverse,
                    nonrecombinant_last, nonrecombinant_list);
            }
        };
        run_post_trim_compatibility(
            background_adjusted,
            tree_compatibility_flow.calls[0].list_distances,
            post_trim_background);
        run_post_trim_compatibility(
            region_adjusted,
            tree_compatibility_flow.calls[3].list_distances,
            post_trim_region);
    }
    const auto metric_close = [&](const double actual, const int metric,
                                  const int role) {
        return std::abs(actual - csv_value(metric, role)) < 0.00051;
    };
    bool generated_consensus_inputs_match = true;
    for (int role = 0; role < 3; ++role) {
        generated_consensus_inputs_match = generated_consensus_inputs_match &&
            metric_close(split_distances.distances[role], 7, role) &&
            split_distances.outlier_index[role] ==
                static_cast<int>(csv_value(8, role)) &&
            metric_close(simple_distances.strengths[role], 1, role) &&
            simple_distances.scores[role] ==
                static_cast<int>(csv_value(2, role)) &&
            outlier_checks[role] == static_cast<int>(csv_value(34, role)) &&
            metric_close(bad_distances[role], 27, role) &&
            metric_close(list_correlations.mismatches[role], 0, role) &&
            metric_close(list_correlations.expected_strength[role], 29, role) &&
            metric_close(list_correlations.absent_strength[role], 30, role) &&
            post_trim_background[role] ==
                static_cast<int>(csv_value(24, role)) &&
            post_trim_region[role] ==
                static_cast<int>(csv_value(25, role));
    }
    if (!generated_consensus_inputs_match) {
        std::cerr << "Generated consensus support mismatch:";
        for (int role = 0; role < 3; ++role) {
            std::cerr << " role" << role << " SSD="
                      << split_distances.distances[role] << '/'
                      << csv_value(7, role) << " OUI="
                      << split_distances.outlier_index[role] << '/'
                      << csv_value(8, role) << " Sim="
                      << simple_distances.strengths[role] << '/'
                      << csv_value(1, role) << ':'
                      << simple_distances.scores[role] << '/'
                      << csv_value(2, role) << " OU="
                      << outlier_checks[role] << '/'
                      << csv_value(34, role) << " Bad="
                      << bad_distances[role] << '/'
                      << csv_value(27, role) << " LC="
                      << list_correlations.mismatches[role] << '/'
                      << csv_value(0, role) << ':'
                      << list_correlations.expected_strength[role] << '/'
                      << csv_value(29, role) << ':'
                      << list_correlations.absent_strength[role] << '/'
                      << csv_value(30, role) << " post="
                      << post_trim_background[role] << '/'
                      << csv_value(24, role) << ':'
                      << post_trim_region[role] << '/'
                      << csv_value(25, role);
        }
        std::cerr << '\n';
    }
    auto connected_consensus_inputs = consensus_inputs;
    connected_consensus_inputs.list_correlation =
        list_correlations.mismatches;
    connected_consensus_inputs.list_correlation_secondary =
        list_correlations.expected_strength;
    connected_consensus_inputs.list_correlation_tertiary =
        list_correlations.absent_strength;
    connected_consensus_inputs.post_trim_compatibility =
        post_trim_background;
    connected_consensus_inputs.post_trim_region_compatibility =
        post_trim_region;
    connected_consensus_inputs.simple_distance_strength =
        simple_distances.strengths;
    connected_consensus_inputs.simple_distance_score = simple_distances.scores;
    connected_consensus_inputs.phylpro = phpr_states[0].scores;
    connected_consensus_inputs.phylpro_secondary = phpr_states[1].scores;
    if (phpr_fixture.header.calls >= 3) {
        connected_consensus_inputs.phylpro_collapsed = phpr_states[2].scores;
    }
    connected_consensus_inputs.subtree_score =
        phpr_states[0].sub_distance_scores;
    connected_consensus_inputs.subtree_phylpro = phpr_states[0].sub_scores;
    connected_consensus_inputs.subtree_score_secondary =
        phpr_states[1].sub_distance_scores;
    connected_consensus_inputs.subtree_phylpro_secondary =
        phpr_states[1].sub_scores;
    connected_consensus_inputs.split_distance = split_distances.distances;
    connected_consensus_inputs.outlier_index = split_distances.outlier_index;
    connected_consensus_inputs.compatibility = tree_compatibility_flow.background;
    connected_consensus_inputs.compatibility_secondary =
        tree_compatibility_flow.background_secondary;
    connected_consensus_inputs.compatibility_tertiary =
        tree_compatibility_flow.background_sets;
    connected_consensus_inputs.compatibility_quaternary =
        tree_compatibility_flow.background_secondary_sets;
    connected_consensus_inputs.region_compatibility =
        tree_compatibility_flow.region;
    connected_consensus_inputs.region_compatibility_secondary =
        tree_compatibility_flow.region_secondary;
    connected_consensus_inputs.region_compatibility_tertiary =
        tree_compatibility_flow.region_sets;
    connected_consensus_inputs.region_compatibility_quaternary =
        tree_compatibility_flow.region_secondary_sets;
    connected_consensus_inputs.triplet_score = triplet_scores;
    connected_consensus_inputs.bad_distances = bad_distances;
    connected_consensus_inputs.outlier_check = outlier_checks;
    connected_consensus_inputs.maximum_distance =
        maximum_distance.maximum_distances;
    connected_consensus_inputs.ranks = simple_distances.ranks;
    bool connected_metrics_match = generated_consensus_inputs_match;
    for (int role = 0; role < 3; ++role) {
        connected_metrics_match = connected_metrics_match &&
            metric_close(connected_consensus_inputs.phylpro[role], 3, role) &&
            metric_close(connected_consensus_inputs.phylpro_secondary[role],
                         4, role) &&
            metric_close(connected_consensus_inputs.phylpro_collapsed[role],
                         5, role) &&
            metric_close(connected_consensus_inputs.subtree_score[role],
                         6, role) &&
            metric_close(connected_consensus_inputs.subtree_phylpro[role],
                         9, role) &&
            metric_close(
                connected_consensus_inputs.subtree_score_secondary[role],
                10, role) &&
            metric_close(
                connected_consensus_inputs.subtree_phylpro_secondary[role],
                11, role) &&
            connected_consensus_inputs.compatibility[role] ==
                static_cast<int>(csv_value(14, role)) &&
            connected_consensus_inputs.compatibility_secondary[role] ==
                static_cast<int>(csv_value(15, role)) &&
            connected_consensus_inputs.compatibility_tertiary[role] ==
                static_cast<int>(csv_value(16, role)) &&
            connected_consensus_inputs.compatibility_quaternary[role] ==
                static_cast<int>(csv_value(17, role)) &&
            connected_consensus_inputs.region_compatibility[role] ==
                static_cast<int>(csv_value(18, role)) &&
            connected_consensus_inputs.region_compatibility_secondary[role] ==
                static_cast<int>(csv_value(19, role)) &&
            connected_consensus_inputs.region_compatibility_tertiary[role] ==
                static_cast<int>(csv_value(20, role)) &&
            connected_consensus_inputs.region_compatibility_quaternary[role] ==
                static_cast<int>(csv_value(21, role)) &&
            metric_close(connected_consensus_inputs.triplet_score[role],
                         26, role) &&
            connected_consensus_inputs.ranks[role][0] ==
                static_cast<int>(csv_value(37, role)) &&
            connected_consensus_inputs.ranks[role][1] ==
                static_cast<int>(csv_value(38, role));
    }
    const auto consensus_state = make_rdp_consensus(consensus_inputs);
    const auto connected_consensus_state =
        make_rdp_consensus(connected_consensus_inputs);
    bool consensus_scores_match = true;
    for (int role = 0; role < 3; ++role) {
        // NN_inputs formats source statistics to four decimal places, so this
        // gate admits only the resulting last-decimal dMax contribution loss.
        consensus_scores_match = consensus_scores_match &&
            std::abs(consensus_state.consensus[role] -
                     csv_value(33, role)) < 0.01;
    }
    const char collect_magic[8] = {
        'C', 'O', 'L', 'L', 'E', 'C', 'T', '1'};
    const auto collect_fixture =
        load_rdp_sectioned_fixture<CollectEventsCaptureHeader>(
            collect_events_fixture_path, collect_magic);
    const auto collect_first =
        rdp_fixture_section<unsigned int>(collect_fixture, 1);
    const auto collect_second =
        rdp_fixture_section<unsigned int>(collect_fixture, 1001);
    if (collect_first.size() != 17 || collect_second.size() != 17) {
        throw std::runtime_error("collect-events fixture header differs");
    }
    const int first_parent = static_cast<int>(collect_first[4]);
    const int second_parent = static_cast<int>(collect_second[4]);
    int expected_winner = -1;
    for (int role = 0; role < 3; ++role) {
        if (role != first_parent && role != second_parent) expected_winner = role;
    }
    const auto expected_final_list =
        rdp_fixture_section<int>(collect_fixture, 5);
    const auto expected_final_last =
        rdp_fixture_section<int>(collect_fixture, 6);
    const auto expected_collect_region_rows =
        rdp_fixture_section<float>(collect_fixture, 7);
    if (expected_collect_region_rows != second_direct_small) {
        std::size_t differences = 0;
        for (std::size_t index = 0;
             index < std::min(expected_collect_region_rows.size(),
                              second_direct_small.size()); ++index) {
            if (expected_collect_region_rows[index] !=
                second_direct_small[index]) ++differences;
        }
        std::cerr << "Second FinalTrim SMatSmall input mismatch: "
                  << differences << " cells\n";
    }
    const auto final_candidates = apply_rdp_strict_group_constraints(
        scan_state.next_no, correlation_sequences, correlation_comparison,
        generated_matrices.background, generated_matrices.event_region,
        first_direct_small, second_direct_small,
        first_adjusted_small, second_adjusted_small,
        downstream_candidates);
    bool second_final_trim_matches = expected_final_last.size() == 3;
    for (int role = 0; role < 3 && second_final_trim_matches; ++role) {
        second_final_trim_matches = expected_final_last[role] ==
            final_candidates.candidate_last[role];
        for (int slot = 0;
             slot <= expected_final_last[role] && second_final_trim_matches;
             ++slot) {
            second_final_trim_matches = expected_final_list[role + slot * 3] ==
                final_candidates.candidate_list[role + slot * 3];
        }
    }
    if (!second_final_trim_matches) {
        std::cerr << "Second FinalTrim lists mismatch:";
        for (int role = 0; role < 3; ++role) {
            std::cerr << " role" << role << '[';
            for (int slot = 0;
                 slot <= final_candidates.candidate_last[role]; ++slot) {
                std::cerr << final_candidates.candidate_list[
                    role + slot * 3] << ',';
            }
            std::cerr << "]/[";
            for (int slot = 0; slot <= expected_final_last[role]; ++slot) {
                std::cerr << expected_final_list[role + slot * 3] << ',';
            }
            std::cerr << ']';
            for (int expected_slot = 0;
                 expected_slot <= expected_final_last[role]; ++expected_slot) {
                const int sequence =
                    expected_final_list[role + expected_slot * 3];
                bool present = false;
                for (int actual_slot = 0;
                     actual_slot <= final_candidates.candidate_last[role];
                     ++actual_slot) {
                    present = present || final_candidates.candidate_list[
                        role + actual_slot * 3] == sequence;
                }
                if (!present) {
                    const int p0 = correlation_comparison[role];
                    const int p1 = correlation_comparison[role + 3];
                    std::cerr << " missing" << sequence << " rows="
                              << second_adjusted_small[role + sequence * 3]
                              << '<' << second_adjusted_small[p0 + sequence * 3]
                              << ',' << second_adjusted_small[p1 + sequence * 3]
                              << ';' << second_direct_small[role + sequence * 3]
                              << '<' << second_direct_small[p0 + sequence * 3]
                              << ',' << second_direct_small[p1 + sequence * 3]
                              << ';' << first_adjusted_small[role + sequence * 3]
                              << '<' << first_adjusted_small[p0 + sequence * 3]
                              << ',' << first_adjusted_small[p1 + sequence * 3]
                              << ';' << first_direct_small[role + sequence * 3]
                              << '<' << first_direct_small[p0 + sequence * 3]
                              << ',' << first_direct_small[p1 + sequence * 3];
                }
            }
        }
        std::cerr << '\n';
    }
    const auto relevant_sequences = make_rdp_relevant_sequences(
        scan_state.next_no, final_candidates.candidate_last,
        final_candidates.candidate_list);
    const std::array<int, 2> collection_trace{trace[0], trace[1]};
    auto selected_collection_event =
        events.xover_list[collection_trace[0]][collection_trace[1] - 1];
    selected_collection_event.distance_holder =
        (selected_collection_event.distance_holder + 0.00000001) * -1.0;
    const auto add_synthetic_events = [&] (
        const RdpFinalTrimState& pass, RdpRawEventState& target_events) {
      for (const auto& seed : pass.synthetic_event_roles) {
        const int role = seed[0];
        const int sequence = seed[1];
        RdpRawEvent synthetic = selected_collection_event;
        synthetic.daughter = static_cast<std::int16_t>(sequence);
        synthetic.major_parent = static_cast<std::int16_t>(
            correlation_sequences[correlation_comparison[role]]);
        synthetic.minor_parent = static_cast<std::int16_t>(
            correlation_sequences[correlation_comparison[role + 3]]);
        synthetic.probability = 0.9;
        synthetic.distance_holder = std::abs(synthetic.distance_holder);
        target_events.xover_list[sequence].push_back(synthetic);
        target_events.current_xover[sequence] = static_cast<std::int16_t>(
            target_events.xover_list[sequence].size());
      }
    };
    std::vector<unsigned char> single_redo(
        static_cast<std::size_t>(scan_state.analysis_list_last + 1), 0);
    single_redo[0] = 1;
    auto rescan_probability_settings = probability_settings;
    rescan_probability_settings.lowest_probability = std::max<double>({
        probability_settings.lowest_probability,
        static_cast<double>(
            selected_collection_event.probability * 100000.0),
        static_cast<double>(probability_settings.lowest_probability *
            probability_settings.mc_correction)});
    const auto run_rescan_pass = [&] (
        const RdpFinalTrimState& pass, RdpRawEventState& target_events) {
      for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= pass.candidate_last[role]; ++slot) {
            const int sequence =
                pass.candidate_list[role + slot * 3];
            if (sequence == correlation_sequences[role]) continue;
            const std::array<int, 3> rescan_triplet{
                sequence,
                correlation_sequences[correlation_comparison[role]],
                correlation_sequences[correlation_comparison[role + 3]]};
            target_events = scan_rdp_redo_triplets(
                scan_state, distance_state, tree_state, single_redo, fss_rdp,
                store_lpv, h.store_lpv_ub, h.fss_rdp_ub, h.xover_window,
                h.xover_window_x, xover_settings, rescan_probability_settings,
                probability_estimate, fact_three, fact, xover_api, 1,
                &target_events, &rescan_triplet);
        }
      }
    };
    // The source FinalTrim invocation clears its temporary XOverList,
    // populates it, and then copies every active slot into persistent PXOList.
    RdpRawEventState events_after_rescan = events;
    events_after_rescan.xover_list[collection_trace[0]][
        collection_trace[1] - 1] = selected_collection_event;
    RdpRawEventState rescan_events;
    const auto run_final_trim_side_effects = [&] (
        const RdpFinalTrimState& pass) {
        RdpRawEventState temporary;
        temporary.current_xover.assign(scan_state.next_no + 1, 0);
        temporary.xover_list.resize(scan_state.next_no + 1);
        add_synthetic_events(pass, temporary);
        run_rescan_pass(pass, temporary);
        append_rdp_events(events_after_rescan, temporary);
        rescan_events = std::move(temporary);
    };
    // FinalTrim temporarily relaxes LowestProb for each RDP rescan tail.
    run_final_trim_side_effects(final_candidates);
    const auto collection_event_list = prepare_rdp_collection_event_list(
        scan_state.next_no, expected_winner, correlation_sequences,
        collection_trace, final_candidates.candidate_last,
        final_candidates.candidate_list,
        final_candidates.acceptable_sequences, events_after_rescan);
    const auto probability_matches = [](const double actual,
                                        const double expected) {
        if (actual == expected) return true;
        const double scale = std::max({std::abs(actual), std::abs(expected),
                                       std::numeric_limits<double>::min()});
        return std::abs(actual - expected) <= scale * 1.0e-12;
    };
    const auto event_matches = [&](const RdpRawEvent& actual,
                                   const XOVERDEFINE& expected) {
        return actual.outside_flag == expected.OutsideFlag &&
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
            probability_matches(actual.permutation_pvalue,
                                expected.PermPVal) &&
            actual.begin_parent == expected.BeginP &&
            actual.end_parent == expected.EndP &&
            probability_matches(actual.probability, expected.Probability) &&
            actual.distance_holder == expected.DHolder;
    };
    const auto differing_event_field = [&](const RdpRawEvent& actual,
                                            const XOVERDEFINE& expected) {
#define RDP_EVENT_FIELD(actual_field, expected_field) \
        if (actual.actual_field != expected.expected_field) return #actual_field
        RDP_EVENT_FIELD(outside_flag, OutsideFlag);
        RDP_EVENT_FIELD(misidentify_flag, MissIdentifyFlag);
        RDP_EVENT_FIELD(program_flag, ProgramFlag);
        RDP_EVENT_FIELD(sbp_flag, SBPFlag);
        RDP_EVENT_FIELD(accept, Accept);
        RDP_EVENT_FIELD(major_parent, MajorP);
        RDP_EVENT_FIELD(minor_parent, MinorP);
        RDP_EVENT_FIELD(daughter, Daughter);
        RDP_EVENT_FIELD(beginning, Beginning);
        RDP_EVENT_FIELD(ending, Ending);
        RDP_EVENT_FIELD(length_holder, LHolder);
        RDP_EVENT_FIELD(event_number, Eventnumber);
        if (!probability_matches(actual.permutation_pvalue,
                                 expected.PermPVal)) {
            return "permutation_pvalue";
        }
        RDP_EVENT_FIELD(begin_parent, BeginP);
        RDP_EVENT_FIELD(end_parent, EndP);
        if (!probability_matches(actual.probability, expected.Probability)) {
            return "probability";
        }
        RDP_EVENT_FIELD(distance_holder, DHolder);
#undef RDP_EVENT_FIELD
        return "none";
    };
    bool make_relevant_matches =
        relevant_sequences.size() ==
            static_cast<std::size_t>(scan_state.next_no + 1);
    for (int role = 0; role < 3; ++role) {
        for (int slot = 0; slot <= final_candidates.candidate_last[role];
             ++slot) {
            make_relevant_matches = make_relevant_matches &&
                relevant_sequences[final_candidates.candidate_list[
                    role + slot * 3]] == 1;
        }
    }
    const auto expected_rdp_locations =
        rdp_fixture_section<int>(collect_fixture, 11);
    const auto expected_rdp_events =
        rdp_fixture_section<XOVERDEFINE>(collect_fixture, 12);
    bool rdp_working_list_matches =
        expected_rdp_locations.size() == expected_rdp_events.size() * 2;
    bool reported_working_list_difference = false;
    std::size_t working_event_index = 0;
    for (int row = 0; row <= scan_state.next_no; ++row) {
        int rdp_ordinal = 0;
        for (std::size_t slot = 0;
             slot < collection_event_list.xover_list[row].size(); ++slot) {
            const auto& event = collection_event_list.xover_list[row][slot];
            if (event.program_flag != 0) continue;
            const bool event_ok =
                working_event_index < expected_rdp_events.size() &&
                expected_rdp_locations[working_event_index * 2] == row &&
                event_matches(event, expected_rdp_events[working_event_index]);
            if (!event_ok && !reported_working_list_difference) {
                reported_working_list_difference = true;
                rdp_working_list_matches = false;
                if (working_event_index >= expected_rdp_events.size()) {
                    std::cerr << "RDP PXOList has extra event at row " << row
                              << " RDP ordinal " << rdp_ordinal << '\n';
                } else {
                const auto& expected = expected_rdp_events[working_event_index];
                std::cerr << "RDP PXOList first mismatch at row " << row
                          << " RDP ordinal " << rdp_ordinal
                          << " (physical slot " << slot + 1
                          << ") expected-location="
                          << expected_rdp_locations[working_event_index * 2]
                          << ',' << expected_rdp_locations[
                              working_event_index * 2 + 1]
                          << " field="
                          << differing_event_field(event, expected)
                          << " event=" << event.daughter << ','
                          << event.major_parent << ',' << event.minor_parent
                          << ',' << event.beginning << ',' << event.ending
                          << ',' << std::setprecision(17)
                          << event.probability << ','
                          << event.distance_holder << " expected="
                          << expected.Daughter << ',' << expected.MajorP << ','
                          << expected.MinorP << ',' << expected.Beginning << ','
                          << expected.Ending << ',' << std::setprecision(17)
                          << expected.Probability
                          << ',' << expected.DHolder << '\n';
                }
            }
            ++working_event_index;
            ++rdp_ordinal;
        }
        if (working_event_index < expected_rdp_events.size() &&
            expected_rdp_locations[working_event_index * 2] == row) {
            if (!reported_working_list_difference) {
                reported_working_list_difference = true;
                const auto& expected = expected_rdp_events[working_event_index];
                std::cerr << "RDP PXOList missing event at row " << row
                          << " RDP ordinal " << rdp_ordinal
                          << " expected-location="
                          << expected_rdp_locations[working_event_index * 2]
                          << ',' << expected_rdp_locations[
                              working_event_index * 2 + 1]
                          << " expected=" << expected.Daughter << ','
                          << expected.MajorP << ',' << expected.MinorP << ','
                          << expected.Beginning << ',' << expected.Ending << ','
                          << expected.Probability << ',' << expected.DHolder
                          << '\n';
            }
            rdp_working_list_matches = false;
            while (working_event_index < expected_rdp_events.size() &&
                   expected_rdp_locations[working_event_index * 2] == row) {
                ++working_event_index;
            }
        }
    }
    rdp_working_list_matches = rdp_working_list_matches &&
        working_event_index == expected_rdp_events.size();
    bool collect_events_matches =
        make_relevant_matches && rdp_working_list_matches;
    const auto check_collect_call = [&](const int call_index,
                                        const auto& raw_header) {
        const int base = call_index * 1000;
        const int role = static_cast<int>(raw_header[4]);
        const int add_num = static_cast<int>(raw_header[7]);
        const auto expected_region_sizes =
            rdp_fixture_section<int>(collect_fixture, base + 2);
        const auto expected_overlap =
            rdp_fixture_section<int>(collect_fixture, base + 3);
        const auto expected_current =
            rdp_fixture_section<short>(collect_fixture, base + 8);
        const auto expected_rdp_current =
            rdp_fixture_section<short>(collect_fixture, base + 9);
        const auto expected_output =
            rdp_fixture_section<XOVERDEFINE>(collect_fixture, base + 10);
        const auto actual = make_rdp_parent_collect_events(
            scan_state.sequence_length, scan_state.next_no, role,
            actual_resolution.region_sizes,
            actual_resolution.event_overlap_mask, correlation_comparison,
            final_candidates.candidate_last,
            final_candidates.candidate_list, add_num,
            correlation_sequences, collection_trace, collection_event_list);
        bool output_matches = actual.result ==
                static_cast<int>(raw_header[11]) &&
            expected_output.size() ==
                static_cast<std::size_t>(scan_state.next_no + 1);
        std::size_t first_output_difference = expected_output.size();
        for (std::size_t slot = 0;
             slot < expected_output.size() && output_matches; ++slot) {
            output_matches = event_matches(actual.events[slot],
                                           expected_output[slot]);
            if (!output_matches) first_output_difference = slot;
        }
        bool active_lists_match = true;
        const auto expected_lists =
            rdp_fixture_section<int>(collect_fixture, base + 5);
        const auto expected_last =
            rdp_fixture_section<int>(collect_fixture, base + 6);
        for (int candidate_role = 0;
             candidate_role < 3 && active_lists_match; ++candidate_role) {
            active_lists_match = expected_last[candidate_role] ==
                final_candidates.candidate_last[candidate_role];
            for (int slot = 0;
                 slot <= expected_last[candidate_role] && active_lists_match;
                 ++slot) {
                active_lists_match = expected_lists[
                    candidate_role + slot * 3] ==
                    final_candidates.candidate_list[
                        candidate_role + slot * 3];
            }
        }
        std::vector<short> actual_rdp_current(scan_state.next_no + 1, 0);
        for (int row = 0; row <= scan_state.next_no; ++row) {
            for (const auto& event : collection_event_list.xover_list[row]) {
                if (event.program_flag == 0) ++actual_rdp_current[row];
            }
        }
        const bool call_matches =
            raw_header[2] == static_cast<unsigned int>(scan_state.next_no) &&
            raw_header[3] ==
                static_cast<unsigned int>(scan_state.sequence_length) &&
            raw_header[12] ==
                static_cast<unsigned int>(correlation_sequences[0]) &&
            raw_header[13] ==
                static_cast<unsigned int>(correlation_sequences[1]) &&
            raw_header[14] ==
                static_cast<unsigned int>(correlation_sequences[2]) &&
            raw_header[15] == static_cast<unsigned int>(collection_trace[0]) &&
            raw_header[16] == static_cast<unsigned int>(collection_trace[1]) &&
            expected_region_sizes[0] == actual_resolution.region_sizes[0] &&
            expected_overlap == actual_resolution.event_overlap_mask &&
            rdp_fixture_section<int>(collect_fixture, base + 4) ==
                as_vector(correlation_comparison) &&
            active_lists_match &&
            expected_rdp_current == actual_rdp_current &&
            output_matches;
        if (!call_matches) {
            std::cerr << "MakeCollecteventsC call " << call_index
                      << " mismatch: role=" << role
                      << " region="
                      << (expected_region_sizes[0] ==
                          actual_resolution.region_sizes[0])
                      << " overlap="
                      << (expected_overlap ==
                          actual_resolution.event_overlap_mask)
                      << " lists="
                      << active_lists_match
                      << " current="
                      << (expected_rdp_current == actual_rdp_current)
                      << " trace=" << collection_trace[0] << ','
                      << collection_trace[1] << '/' << raw_header[15] << ','
                      << raw_header[16] << " output=" << output_matches;
            if (expected_rdp_current != actual_rdp_current) {
                std::cerr << " current-diffs=";
                for (int row = 0; row <= scan_state.next_no; ++row) {
                    if (expected_rdp_current[row] != actual_rdp_current[row]) {
                        std::cerr << row << ':'
                                  << actual_rdp_current[row]
                                  << '/' << expected_rdp_current[row] << ',';
                    }
                }
            }
            if (!output_matches &&
                first_output_difference < expected_output.size()) {
                const auto& observed = actual.events[first_output_difference];
                const auto& expected = expected_output[first_output_difference];
                std::cerr << " first-output=" << first_output_difference
                          << " d=" << observed.daughter << '/'
                          << expected.Daughter << " ma="
                          << observed.major_parent << '/' << expected.MajorP
                          << " mi=" << observed.minor_parent << '/'
                          << expected.MinorP << " p="
                          << observed.probability << '/'
                          << expected.Probability << " pp="
                          << observed.permutation_pvalue << '/'
                          << expected.PermPVal << " dh="
                          << observed.distance_holder << '/'
                          << expected.DHolder << " flags="
                          << static_cast<int>(observed.outside_flag) << ','
                          << static_cast<int>(observed.misidentify_flag) << ','
                          << static_cast<int>(observed.sbp_flag) << '/'
                          << static_cast<int>(expected.OutsideFlag) << ','
                          << static_cast<int>(expected.MissIdentifyFlag) << ','
                          << static_cast<int>(expected.SBPFlag) << " bp="
                          << observed.beginning << ',' << observed.ending << '/'
                          << expected.Beginning << ',' << expected.Ending;
            }
            std::cerr << '\n';
        }
        return call_matches;
    };
    const bool collect_first_matches = check_collect_call(0, collect_first);
    const bool collect_second_matches = check_collect_call(1, collect_second);
    collect_events_matches = collect_events_matches &&
        collect_first_matches && collect_second_matches;

    const char modseq_magic[8] = {
        'M', 'S', 'E', 'Q', 'Y', 'V', '1', '\0'};
    const auto modseq_fixture =
        load_rdp_sectioned_fixture<ModSeqNumYCaptureHeader>(
            modseqnumy_fixture_path, modseq_magic);
    const auto modseq_header = rdp_fixture_section<int>(modseq_fixture, 1);
    const auto modseq_inputs = rdp_fixture_section<int>(modseq_fixture, 2);
    const auto modseq_result = rdp_fixture_section<int>(modseq_fixture, 3);
    const auto modseq_outputs = rdp_fixture_section<int>(modseq_fixture, 4);
    const auto expected_mutated_sequences =
        rdp_fixture_section<short>(modseq_fixture, 5);
    const auto expected_saved_tracts =
        rdp_fixture_section<short>(modseq_fixture, 6);
    const auto expected_mutation_missing =
        rdp_fixture_section<unsigned char>(modseq_fixture, 7);
    bool mutation_metadata_matches = modseq_fixture.header.calls > 0 &&
        modseq_header.size() == 10 && modseq_result == std::vector<int>{1} &&
        modseq_header[2] == selected.beginning &&
        modseq_header[3] == selected.ending &&
        modseq_header[4] == scan_state.sequence_length &&
        modseq_header[5] == expected_winner;
    std::vector<int> generated_breakpoints(
        static_cast<std::size_t>(2 *
            (final_candidates.candidate_last[expected_winner] + 1)), 0);
    for (int slot = 0;
         slot <= final_candidates.candidate_last[expected_winner]; ++slot) {
        // In the RDP-only MakeBreaks path every retained winner has a
        // legitimate ProgramFlag-zero record. The literal assignment uses
        // the selected BPos/EPos, not that record's own interval.
        generated_breakpoints[slot * 2] = selected.beginning;
        generated_breakpoints[slot * 2 + 1] = selected.ending;
    }
    for (std::size_t offset = 0; offset + 4 < modseq_inputs.size();
         offset += 5) {
        const int role = modseq_inputs[offset];
        const int slot = modseq_inputs[offset + 1];
        mutation_metadata_matches = mutation_metadata_matches &&
            role >= 0 && role < 3 && slot >= 0 &&
            slot <= final_candidates.candidate_last[role] &&
            modseq_inputs[offset + 2] ==
                final_candidates.candidate_list[role + slot * 3];
        if (role == expected_winner) {
            mutation_metadata_matches = mutation_metadata_matches &&
                generated_breakpoints[slot * 2] ==
                    modseq_inputs[offset + 3] &&
                generated_breakpoints[slot * 2 + 1] ==
                    modseq_inputs[offset + 4];
        }
    }
    std::vector<unsigned char> initial_missing(scan_state.sequence_data.size(), 0);
    const auto mutation_state = erase_rdp_recombinant_tracts(
        scan_state.sequence_length, expected_winner,
        final_candidates.candidate_last, final_candidates.candidate_list,
        selected.beginning, selected.ending, generated_breakpoints,
        scan_state.sequence_data, initial_missing);
    const int selected_count =
        final_candidates.candidate_last[expected_winner] + 1;
    std::vector<short> actual_mutated_sequences;
    actual_mutated_sequences.reserve(
        static_cast<std::size_t>(selected_count) *
            (scan_state.sequence_length + 1));
    std::vector<unsigned char> actual_mutation_missing;
    actual_mutation_missing.reserve(
        static_cast<std::size_t>(selected_count) *
            (scan_state.sequence_length + 1));
    for (int slot = 0; slot < selected_count; ++slot) {
        const int sequence =
            final_candidates.candidate_list[expected_winner + slot * 3];
        const auto begin = mutation_state.sequence_data.begin() +
            static_cast<std::ptrdiff_t>(sequence) *
                (scan_state.sequence_length + 1);
        actual_mutated_sequences.insert(
            actual_mutated_sequences.end(), begin,
            begin + scan_state.sequence_length + 1);
        const auto missing_begin = mutation_state.missing_data.begin() +
            static_cast<std::ptrdiff_t>(sequence) *
                (scan_state.sequence_length + 1);
        actual_mutation_missing.insert(
            actual_mutation_missing.end(), missing_begin,
            missing_begin + scan_state.sequence_length + 1);
        mutation_metadata_matches = mutation_metadata_matches &&
            modseq_outputs[slot * 5] == expected_winner &&
            modseq_outputs[slot * 5 + 1] == slot &&
            modseq_outputs[slot * 5 + 2] == sequence &&
            modseq_outputs[slot * 5 + 3] ==
                mutation_state.breakpoints[slot * 2] &&
            modseq_outputs[slot * 5 + 4] ==
                mutation_state.breakpoints[slot * 2 + 1];
    }
    const bool modseqnumy_matches = mutation_metadata_matches &&
        actual_mutated_sequences == expected_mutated_sequences &&
        mutation_state.saved_tracts == expected_saved_tracts &&
        actual_mutation_missing == expected_mutation_missing;
    if (!modseqnumy_matches) {
        std::cerr << "ModSeqNumY mismatch: metadata="
                  << mutation_metadata_matches << " sequence="
                  << (actual_mutated_sequences == expected_mutated_sequences)
                  << " saved="
                  << (mutation_state.saved_tracts == expected_saved_tracts)
                  << " missing="
                  << (actual_mutation_missing == expected_mutation_missing)
                  << '\n';
    }
    const char mutation_tail_magic[8] = {
        'M', 'T', 'A', 'I', 'L', 'V', '1', '\0'};
    const auto mutation_tail_fixture =
        load_rdp_sectioned_fixture<MutationTailCaptureHeader>(
            mutation_tail_fixture_path, mutation_tail_magic);
    const auto msn_header =
        rdp_fixture_section<int>(mutation_tail_fixture, 1);
    const auto msn_records =
        rdp_fixture_section<int>(mutation_tail_fixture, 2);
    const auto expected_fragment_sequences =
        rdp_fixture_section<short>(mutation_tail_fixture, 3);
    const auto expected_fragment_missing =
        rdp_fixture_section<unsigned char>(mutation_tail_fixture, 4);
    const auto modz_header =
        rdp_fixture_section<int>(mutation_tail_fixture, 5);
    const auto modz_records =
        rdp_fixture_section<int>(mutation_tail_fixture, 6);
    const auto expected_original_sequences =
        rdp_fixture_section<short>(mutation_tail_fixture, 7);
    const auto expected_original_missing =
        rdp_fixture_section<unsigned char>(mutation_tail_fixture, 8);
    const auto size_header =
        rdp_fixture_section<int>(mutation_tail_fixture, 9);
    const auto expected_sizes =
        rdp_fixture_section<int>(mutation_tail_fixture, 10);
    const int expanded_next_no = scan_state.next_no + selected_count;
    auto expanded_sequences = mutation_state.sequence_data;
    auto expanded_missing = mutation_state.missing_data;
    rebuild_rdp_recombinant_tracts(
        scan_state.sequence_length, expected_winner,
        final_candidates.candidate_last, final_candidates.candidate_list,
        mutation_state.breakpoints, mutation_state.saved_tracts,
        expanded_sequences);
    expanded_sequences.resize(
        static_cast<std::size_t>(expanded_next_no + 1) *
            (scan_state.sequence_length + 1), 0);
    expanded_missing.resize(expanded_sequences.size(), 0);
    auto tail_breakpoints = mutation_state.breakpoints;
    make_rdp_fragment_rows(
        scan_state.sequence_length, expanded_next_no, expected_winner,
        final_candidates.candidate_last, final_candidates.candidate_list,
        selected.beginning, selected.ending, tail_breakpoints,
        expanded_sequences, expanded_missing);
    std::vector<short> actual_fragment_sequences;
    std::vector<unsigned char> actual_fragment_missing;
    bool mutation_tail_metadata_matches =
        mutation_tail_fixture.header.calls > 0 && msn_header.size() == 11 &&
        msn_header[2] == expanded_next_no &&
        msn_header[3] == scan_state.sequence_length &&
        msn_header[4] == selected.beginning &&
        msn_header[5] == selected.ending &&
        msn_header[6] == expected_winner && msn_header[10] == 1;
    for (int slot = 0; slot < selected_count; ++slot) {
        const int source =
            final_candidates.candidate_list[expected_winner + slot * 3];
        const int created = expanded_next_no -
            final_candidates.candidate_last[expected_winner] + slot;
        mutation_tail_metadata_matches = mutation_tail_metadata_matches &&
            msn_records[slot * 5] == slot &&
            msn_records[slot * 5 + 1] == source &&
            msn_records[slot * 5 + 2] == created &&
            msn_records[slot * 5 + 3] == tail_breakpoints[slot * 2] &&
            msn_records[slot * 5 + 4] == tail_breakpoints[slot * 2 + 1];
        const auto sequence_begin = expanded_sequences.begin() +
            static_cast<std::ptrdiff_t>(created) *
                (scan_state.sequence_length + 1);
        actual_fragment_sequences.insert(
            actual_fragment_sequences.end(), sequence_begin,
            sequence_begin + scan_state.sequence_length + 1);
        const auto missing_begin = expanded_missing.begin() +
            static_cast<std::ptrdiff_t>(created) *
                (scan_state.sequence_length + 1);
        actual_fragment_missing.insert(
            actual_fragment_missing.end(), missing_begin,
            missing_begin + scan_state.sequence_length + 1);
    }
    std::vector<int> compact_candidates(selected_count, 0);
    for (int slot = 0; slot < selected_count; ++slot) {
        compact_candidates[slot] =
            final_candidates.candidate_list[expected_winner + slot * 3];
    }
    erase_rdp_original_tracts(
        scan_state.sequence_length, expanded_next_no, expected_winner,
        final_candidates.candidate_last, compact_candidates,
        selected.beginning, selected.ending, tail_breakpoints,
        expanded_sequences, expanded_missing);
    std::vector<short> actual_original_sequences;
    std::vector<unsigned char> actual_original_missing;
    for (int slot = 0; slot < selected_count; ++slot) {
        const int source = compact_candidates[slot];
        mutation_tail_metadata_matches = mutation_tail_metadata_matches &&
            modz_records[slot * 4] == slot &&
            modz_records[slot * 4 + 1] == source &&
            modz_records[slot * 4 + 2] == tail_breakpoints[slot * 2] &&
            modz_records[slot * 4 + 3] == tail_breakpoints[slot * 2 + 1];
        const auto sequence_begin = expanded_sequences.begin() +
            static_cast<std::ptrdiff_t>(source) *
                (scan_state.sequence_length + 1);
        actual_original_sequences.insert(
            actual_original_sequences.end(), sequence_begin,
            sequence_begin + scan_state.sequence_length + 1);
        const auto missing_begin = expanded_missing.begin() +
            static_cast<std::ptrdiff_t>(source) *
                (scan_state.sequence_length + 1);
        actual_original_missing.insert(
            actual_original_missing.end(), missing_begin,
            missing_begin + scan_state.sequence_length + 1);
    }
    const auto actual_sizes = calculate_rdp_actual_sequence_sizes(
        scan_state.sequence_length, expanded_next_no, expected_winner,
        final_candidates.candidate_last, final_candidates.candidate_list,
        expanded_sequences);
    mutation_tail_metadata_matches = mutation_tail_metadata_matches &&
        modz_header.size() == 11 && modz_header[2] == expanded_next_no &&
        modz_header[6] == expected_winner && modz_header[10] == 1 &&
        size_header.size() == 9 && size_header[2] == expanded_next_no &&
        size_header[4] == expected_winner && size_header[8] == 1;
    for (int slot = 0; slot < selected_count; ++slot) {
        const int source = compact_candidates[slot];
        const int created = expanded_next_no -
            final_candidates.candidate_last[expected_winner] + slot;
        mutation_tail_metadata_matches = mutation_tail_metadata_matches &&
            expected_sizes[slot * 5] == slot &&
            expected_sizes[slot * 5 + 1] == source &&
            expected_sizes[slot * 5 + 2] == created &&
            expected_sizes[slot * 5 + 3] == actual_sizes[source] &&
            expected_sizes[slot * 5 + 4] == actual_sizes[created];
    }
    const bool mutation_tail_matches = mutation_tail_metadata_matches &&
        actual_fragment_sequences == expected_fragment_sequences &&
        actual_fragment_missing == expected_fragment_missing &&
        actual_original_sequences == expected_original_sequences &&
        actual_original_missing == expected_original_missing;
    if (!mutation_tail_matches) {
        std::cerr << "Mutation tail mismatch: metadata="
                  << mutation_tail_metadata_matches << " new-sequence="
                  << (actual_fragment_sequences == expected_fragment_sequences)
                  << " new-missing="
                  << (actual_fragment_missing == expected_fragment_missing)
                  << " original-sequence="
                  << (actual_original_sequences == expected_original_sequences)
                  << " original-missing="
                  << (actual_original_missing == expected_original_missing)
                  << '\n';
    }
    std::vector<int> trace_sub(expanded_next_no + 1, 0);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        trace_sub[sequence] = sequence;
    }
    for (int slot = 0; slot < selected_count; ++slot) {
        const int created = expanded_next_no -
            final_candidates.candidate_last[expected_winner] + slot;
        trace_sub[created] = compact_candidates[slot];
    }
    std::vector<int> initial_actual_sizes(scan_state.next_no + 1, 0);
    for (int sequence = 0; sequence <= scan_state.next_no; ++sequence) {
        for (int position = 1; position <= scan_state.sequence_length;
             ++position) {
            if (scan_state.sequence_data[
                    position + sequence * (scan_state.sequence_length + 1)] !=
                46) {
                ++initial_actual_sizes[sequence];
            }
        }
    }
    const auto redistributed = redistribute_rdp_events(
        scan_state.next_no, expected_winner, h.lowest_probability,
        final_candidates.candidate_last, final_candidates.candidate_list,
        trace_sub, collection_event_list);
    const auto native_addjust_pairs = load_first_addjust_pairs(
        addjust_pairs_trace_path, scan_state.next_no);
    const auto native_redistributed =
        load_first_addjust_events(addjust_pairs_trace_path);
    const bool redistribution_pairs_match =
        redistributed.pairs_to_rescan == native_addjust_pairs;
    bool redistribution_events_match =
        redistributed.events.current_xover ==
            native_redistributed.current_xover &&
        redistributed.events.xover_list.size() ==
            native_redistributed.xover_list.size();
    std::size_t redistribution_event_mismatches = 0;
    for (std::size_t row = 0;
         row < redistributed.events.xover_list.size() &&
         row < native_redistributed.xover_list.size(); ++row) {
        const auto& actual_row = redistributed.events.xover_list[row];
        const auto& expected_row = native_redistributed.xover_list[row];
        if (actual_row.size() != expected_row.size()) {
            redistribution_events_match = false;
            redistribution_event_mismatches +=
                actual_row.size() > expected_row.size()
                ? actual_row.size() - expected_row.size()
                : expected_row.size() - actual_row.size();
        }
        for (std::size_t slot = 0;
             slot < actual_row.size() && slot < expected_row.size(); ++slot) {
            if (!raw_event_equivalent(actual_row[slot], expected_row[slot])) {
                redistribution_events_match = false;
                ++redistribution_event_mismatches;
            }
        }
    }
    if (!redistribution_events_match) {
        const auto count_events = [](const RdpRawEventState& state) {
            std::size_t count = 0;
            for (const auto& row : state.xover_list) count += row.size();
            return count;
        };
        std::cerr << "AddjustCXO event mismatch: "
                  << redistribution_event_mismatches << " records, counts="
                  << (redistributed.events.current_xover ==
                      native_redistributed.current_xover)
                  << " input totals=" << count_events(events) << '+'
                  << (count_events(events_after_rescan) - count_events(events))
                  << '=' << count_events(collection_event_list)
                  << " scan significant=" << rescan_events.significant_candidates
                  << " threshold="
                  << rescan_probability_settings.lowest_probability
                  << " candidates=" << final_candidates.candidate_last[0] + 1
                  << ',' << final_candidates.candidate_last[1] + 1
                  << ',' << final_candidates.candidate_last[2] + 1
                  << " first=" << downstream_candidates.candidate_last[0] + 1
                  << ',' << downstream_candidates.candidate_last[1] + 1
                  << ',' << downstream_candidates.candidate_last[2] + 1
                  << " winner=" << expected_winner
                  << " triplet=" << correlation_sequences[0] << ','
                  << correlation_sequences[1] << ','
                  << correlation_sequences[2]
                  << '\n';
        for (std::size_t row = 0;
             row < redistributed.events.xover_list.size(); ++row) {
            if (redistributed.events.xover_list[row].size() !=
                native_redistributed.xover_list[row].size()) {
                std::cerr << "  input row " << row
                          << " initial=" << events.xover_list[row].size()
                          << " rescanned="
                          << (events_after_rescan.xover_list[row].size() -
                              events.xover_list[row].size())
                          << " collected="
                          << collection_event_list.xover_list[row].size()
                          << '\n';
            }
        }
    }
    auto propagated_pairs = native_addjust_pairs;
    propagate_rdp_group_pairs(
        scan_state.next_no, expected_winner,
        final_candidates.candidate_last, final_candidates.candidate_list,
        propagated_pairs);
    const auto inner_triplets = make_rdp_inner_scan_triplets(
        scan_state, events.triplets_with_events, expected_winner,
        final_candidates.candidate_last,
        final_candidates.candidate_list, trace_sub, initial_actual_sizes,
        scan_state.next_no, 20, propagated_pairs);
    const auto expanded_scan_state = rebuild_rdp_scan_state(
        expanded_next_no, scan_state.sequence_length, expanded_sequences,
        preprocess_api);
    auto expanded_actual_sizes = initial_actual_sizes;
    expanded_actual_sizes.resize(expanded_next_no + 1, 0);
    for (int slot = 0; slot < selected_count; ++slot) {
        const int source = compact_candidates[slot];
        const int created = expanded_next_no -
            final_candidates.candidate_last[expected_winner] + slot;
        expanded_actual_sizes[source] = actual_sizes[source];
        expanded_actual_sizes[created] = actual_sizes[created];
    }
    const auto expanded_distance_state = build_rdp_distance_state(
        expanded_scan_state, 1, expanded_scan_state.sequence_length);
    const auto outer_triplets = make_rdp_outer_scan_triplets(
        scan_state, events.triplets_with_events, expanded_next_no,
        scan_state.next_no,
        expected_winner, final_candidates.candidate_last, trace_sub,
        expanded_actual_sizes, scan_state.next_no, 20, propagated_pairs,
        expanded_next_no, expanded_distance_state.valid_sites);
    const auto alist_rdp3_calls =
        load_alist_rdp3_trace(alist_rdp3_trace_path);
    const bool inner_schedule_matches = !alist_rdp3_calls.empty() &&
        alist_rdp3_calls[0].next_no == scan_state.next_no &&
        alist_rdp3_calls[0].triplets == inner_triplets;
    const bool outer_schedule_matches = alist_rdp3_calls.size() > 1 &&
        alist_rdp3_calls[1].next_no == expanded_next_no &&
        alist_rdp3_calls[1].triplets == outer_triplets;
    const bool rescan_schedule_matches =
        inner_schedule_matches && outer_schedule_matches;

    auto rescan_fss_rdp = fss_rdp;
    auto rescan_probability_estimate = probability_estimate;
    auto rescan_fact_three = fact_three;
    auto rescan_fact = fact;
    const auto inner_scan_state = rebuild_rdp_scan_state(
        scan_state.next_no, scan_state.sequence_length,
        mutation_state.sequence_data, preprocess_api);
    const auto inner_redo = screen_rdp_rescan_triplets(
        inner_triplets, inner_scan_state, distance_state, tree_state, h,
        rescan_fss_rdp, rescan_probability_estimate, rescan_fact_three,
        rescan_fact);
    auto inner_event_scan_state = inner_scan_state;
    inner_event_scan_state.analysis_list = flatten_rdp_triplets(inner_triplets);
    inner_event_scan_state.analysis_list_last =
        static_cast<int>(inner_triplets.size()) - 1;
    const auto inner_rescan_events = scan_rdp_redo_triplets(
        inner_event_scan_state, distance_state, tree_state, inner_redo,
        rescan_fss_rdp, store_lpv, h.store_lpv_ub, h.fss_rdp_ub,
        h.xover_window, h.xover_window_x, xover_settings,
        probability_settings, rescan_probability_estimate,
        rescan_fact_three, rescan_fact, xover_api, 0, nullptr, nullptr, 1,
        &mutation_state.missing_data);
    auto combined_rescan_events = redistributed.events;
    append_rdp_events(combined_rescan_events, inner_rescan_events);
    combined_rescan_events.xover_list.resize(expanded_next_no + 1);
    combined_rescan_events.current_xover.resize(expanded_next_no + 1, 0);
    const auto expanded_tree_state = build_rdp_upgma_tree_state(
        expanded_next_no, expanded_distance_state);
    const auto outer_redo = screen_rdp_rescan_triplets(
        outer_triplets, expanded_scan_state, expanded_distance_state,
        expanded_tree_state, h, rescan_fss_rdp,
        rescan_probability_estimate, rescan_fact_three, rescan_fact);
    auto outer_event_scan_state = expanded_scan_state;
    outer_event_scan_state.analysis_list = flatten_rdp_triplets(outer_triplets);
    outer_event_scan_state.analysis_list_last =
        static_cast<int>(outer_triplets.size()) - 1;
    std::vector<double> expanded_store_lpv(
        static_cast<std::size_t>(expanded_next_no + 1) *
            (h.store_lpv_ub + 1),
        1.0);
    const auto outer_rescan_events = scan_rdp_redo_triplets(
        outer_event_scan_state, expanded_distance_state, expanded_tree_state,
        outer_redo, rescan_fss_rdp, expanded_store_lpv, h.store_lpv_ub,
        h.fss_rdp_ub, h.xover_window, h.xover_window_x, xover_settings,
        probability_settings, rescan_probability_estimate,
        rescan_fact_three, rescan_fact, xover_api, 0, nullptr, nullptr, 1,
        &expanded_missing);
    const auto dropped_rescan_state = drop_rdp_unused_fragment_events(
        scan_state.next_no, expanded_next_no, 20, trace_sub,
        expanded_actual_sizes, combined_rescan_events, outer_rescan_events);
    const auto& post_rescan_events = dropped_rescan_state.events;
    const auto native_post_rescan =
        load_findbetter_events(addjust_pairs_trace_path, 2);
    const auto [post_rescan_events_match, post_rescan_mismatches] =
        compare_rdp_event_states(post_rescan_events, native_post_rescan);
    if (!post_rescan_events_match) {
        std::cerr << "Post-rescan event mismatch: "
                  << post_rescan_mismatches << " records, counts="
                  << (post_rescan_events.current_xover ==
                      native_post_rescan.current_xover)
                  << " totals=";
        const auto total_events = [](const RdpRawEventState& state) {
            std::size_t total = 0;
            for (const auto& row : state.xover_list) total += row.size();
            return total;
        };
        std::cerr << total_events(post_rescan_events) << '/'
                  << total_events(native_post_rescan)
                  << " retained=" << total_events(redistributed.events)
                  << " inner=" << total_events(inner_rescan_events)
                  << " outer=" << total_events(outer_rescan_events)
                  << " next=" << expanded_next_no << "->"
                  << dropped_rescan_state.next_no << " refs=";
        for (int sequence = scan_state.next_no + 1;
             sequence <= dropped_rescan_state.next_no; ++sequence) {
            std::cerr << sequence << ':'
                      << dropped_rescan_state.reference_counts[sequence]
                      << ',';
        }
        std::cerr << '\n';
        for (std::size_t row = 0;
             row < post_rescan_events.xover_list.size() &&
             row < native_post_rescan.xover_list.size(); ++row) {
            const auto actual_count = post_rescan_events.xover_list[row].size();
            const auto expected_count = native_post_rescan.xover_list[row].size();
            if (actual_count != expected_count) {
                std::cerr << "  counts row " << row << " actual="
                          << actual_count << " native=" << expected_count
                          << " retained="
                          << (row < redistributed.events.xover_list.size()
                                  ? redistributed.events.xover_list[row].size()
                                  : 0)
                          << " inner="
                          << (row < inner_rescan_events.xover_list.size()
                                  ? inner_rescan_events.xover_list[row].size()
                                  : 0)
                          << " outer="
                          << (row < outer_rescan_events.xover_list.size()
                                  ? outer_rescan_events.xover_list[row].size()
                                  : 0)
                          << '\n';
                const std::size_t retained_count =
                    row < redistributed.events.xover_list.size()
                    ? redistributed.events.xover_list[row].size() : 0;
                std::cerr << "    actual-new";
                for (std::size_t slot = retained_count;
                     slot < post_rescan_events.xover_list[row].size(); ++slot) {
                    const auto& event = post_rescan_events.xover_list[row][slot];
                    std::cerr << " [" << event.major_parent << ','
                              << event.minor_parent << ':' << event.beginning
                              << '-' << event.ending << ':'
                              << event.probability << ']';
                }
                std::cerr << "\n    native-new";
                for (std::size_t slot = retained_count;
                     slot < native_post_rescan.xover_list[row].size(); ++slot) {
                    const auto& event = native_post_rescan.xover_list[row][slot];
                    std::cerr << " [" << event.major_parent << ','
                              << event.minor_parent << ':' << event.beginning
                              << '-' << event.ending << ':'
                              << event.probability << ']';
                }
                std::cerr << '\n';
            }
        }
        int reports = 0;
        for (std::size_t row = 0;
             row < post_rescan_events.xover_list.size() &&
             row < native_post_rescan.xover_list.size(); ++row) {
            const auto& actual_row = post_rescan_events.xover_list[row];
            const auto& expected_row = native_post_rescan.xover_list[row];
            if (actual_row.size() != expected_row.size() && reports < 8) {
                std::cerr << "  row " << row << " count "
                          << actual_row.size() << '/' << expected_row.size()
                          << '\n';
                ++reports;
            }
            for (std::size_t slot = 0;
                 slot < actual_row.size() && slot < expected_row.size() &&
                 reports < 8; ++slot) {
                if (raw_event_equivalent(actual_row[slot],
                                         expected_row[slot])) {
                    continue;
                }
                const auto& actual = actual_row[slot];
                const auto& expected = expected_row[slot];
                std::cerr << "  row " << row << " slot " << slot + 1
                          << " roles " << actual.daughter << ','
                          << actual.major_parent << ',' << actual.minor_parent
                          << '/' << expected.daughter << ','
                          << expected.major_parent << ','
                          << expected.minor_parent << " bp "
                          << actual.beginning << '-' << actual.ending << '/'
                          << expected.beginning << '-' << expected.ending
                          << " event " << actual.event_number << '/'
                          << expected.event_number << " p "
                          << actual.probability << '/'
                          << expected.probability << " flags(out,mis,prog,sbp,acc) "
                          << static_cast<int>(actual.outside_flag) << ','
                          << static_cast<int>(actual.misidentify_flag) << ','
                          << static_cast<int>(actual.program_flag) << ','
                          << static_cast<int>(actual.sbp_flag) << ','
                          << static_cast<int>(actual.accept) << '/'
                          << static_cast<int>(expected.outside_flag) << ','
                          << static_cast<int>(expected.misidentify_flag) << ','
                          << static_cast<int>(expected.program_flag) << ','
                          << static_cast<int>(expected.sbp_flag) << ','
                          << static_cast<int>(expected.accept)
                          << " holders(len,perm,bp,ep,dist) "
                          << actual.length_holder << ','
                          << actual.permutation_pvalue << ','
                          << actual.begin_parent << ',' << actual.end_parent
                          << ',' << actual.distance_holder << '/'
                          << expected.length_holder << ','
                          << expected.permutation_pvalue << ','
                          << expected.begin_parent << ','
                          << expected.end_parent << ','
                          << expected.distance_holder << '\n';
                ++reports;
            }
        }
    }
    const bool inner_screen_matches = !alist_rdp3_calls.empty() &&
        inner_redo == alist_rdp3_calls[0].redo;
    const bool outer_screen_matches = alist_rdp3_calls.size() > 1 &&
        outer_redo == alist_rdp3_calls[1].redo;
    const bool rescan_screen_matches =
        inner_screen_matches && outer_screen_matches;
    if (!rescan_schedule_matches) {
        std::cerr << "RDP rescan schedule mismatch: pairs="
                  << redistribution_pairs_match << " inner="
                  << inner_triplets.size() << '/'
                  << (alist_rdp3_calls.empty() ? 0 :
                      alist_rdp3_calls[0].triplets.size())
                  << " outer=" << outer_triplets.size() << '/'
                  << (alist_rdp3_calls.size() < 2 ? 0 :
                      alist_rdp3_calls[1].triplets.size());
        const auto report_first = [&](const char* label, const auto& actual,
                                      const auto& expected) {
            std::size_t first = 0;
            while (first < actual.size() && first < expected.size() &&
                   actual[first] == expected[first]) {
                ++first;
            }
            std::cerr << ' ' << label << "-first=" << first;
            if (first < actual.size()) {
                std::cerr << " actual=" << actual[first][0] << ','
                          << actual[first][1] << ',' << actual[first][2];
            }
            if (first < expected.size()) {
                std::cerr << " expected=" << expected[first][0] << ','
                          << expected[first][1] << ',' << expected[first][2];
            }
        };
        if (!alist_rdp3_calls.empty()) {
            report_first("inner", inner_triplets,
                         alist_rdp3_calls[0].triplets);
        }
        if (alist_rdp3_calls.size() > 1) {
            report_first("outer", outer_triplets,
                         alist_rdp3_calls[1].triplets);
        }
        std::cerr << '\n';
    }
    if (!rescan_screen_matches) {
        const auto mismatch_count = [](const auto& actual,
                                       const auto& expected) {
            std::size_t count = actual.size() > expected.size()
                ? actual.size() - expected.size()
                : expected.size() - actual.size();
            for (std::size_t index = 0;
                 index < actual.size() && index < expected.size(); ++index) {
                count += actual[index] != expected[index];
            }
            return count;
        };
        std::cerr << "AlistRDP3 rescan mismatch: inner="
                  << mismatch_count(
                         inner_redo, alist_rdp3_calls.empty()
                             ? std::vector<unsigned char>{}
                             : alist_rdp3_calls[0].redo)
                  << " outer="
                  << mismatch_count(
                         outer_redo, alist_rdp3_calls.size() < 2
                             ? std::vector<unsigned char>{}
                             : alist_rdp3_calls[1].redo)
                  << '\n';
    }
    const bool consensus_matches = consensus_scores_match &&
        consensus_state.winning_role == expected_winner &&
        connected_metrics_match &&
        connected_consensus_state.winning_role == expected_winner;
    if (!consensus_matches) {
        std::cerr << "MakeConsensusC mismatch: scores="
                  << consensus_scores_match << " winner="
                  << consensus_state.winning_role << '/' << expected_winner
                  << " connected=" << connected_metrics_match << ':'
                  << connected_consensus_state.winning_role
                  << " second-trim=" << second_final_trim_matches
                  << " raw=";
        for (int role = 0; role < 3; ++role) {
            std::cerr << consensus_state.consensus[role] << '/'
                      << csv_value(33, role) << ',';
        }
        std::cerr << '\n';
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
              << (score_support_matches ? "PASS" : "FAIL")
              << ", CheckPatternX "
              << (check_pattern_matches ? "PASS" : "FAIL")
              << ", FinalTrim prefix "
              << (final_trim_prefix_matches ? "PASS" : "FAIL")
              << ", CalcMaxD " << (cmaxd_matches ? "PASS" : "FAIL")
              << ", MakeConsensusC "
              << (consensus_matches ? "PASS" : "FAIL")
              << ", second FinalTrim "
              << (second_final_trim_matches ? "PASS" : "FAIL")
              << ", MakeRelevant/MakeCollecteventsC "
              << (collect_events_matches ? "PASS" : "FAIL")
              << ", ModSeqNumY "
              << (modseqnumy_matches ? "PASS" : "FAIL")
              << ", ModSN/ModSeqNumZ "
              << (mutation_tail_matches ? "PASS" : "FAIL")
              << ", AddjustCXO events "
              << (redistribution_events_match ? "PASS" : "FAIL")
              << ", rescan schedules "
              << (rescan_schedule_matches ? "PASS" : "FAIL")
              << ", AlistRDP3 rescans "
              << (rescan_screen_matches ? "PASS" : "FAIL")
              << ", post-rescan events "
              << (post_rescan_events_match ? "PASS" : "FAIL")
              << "\n";

    return make_test_structure_matches && first_selection_matches &&
            ufdist_matches && region_distance_matches &&
            check_matrix_matches && make_nj_matches ?
        (make_sdmp_matches && fill_rmat_matches && calcr_matches &&
         make_rlist_matches && find_actual_matches && strip_dup_matches &&
             rcompat_matches && rcompat_flow_matches && phpr_matches &&
             score_support_matches && check_pattern_matches &&
             final_trim_prefix_matches && cmaxd_matches && consensus_matches &&
             second_final_trim_matches && collect_events_matches
             && modseqnumy_matches
             && mutation_tail_matches
             && redistribution_events_match
             && rescan_schedule_matches
             && rescan_screen_matches
             && post_rescan_events_match
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
    if (argc == 30 &&
        std::string_view(argv[1]) == "fasta-all-redo-events-fixture") {
        return fasta_all_redo_events_fixture(
            argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8],
            argv[9], argv[10], argv[11], {argv[12], argv[13], argv[14]},
            argv[15], argv[16], argv[17], argv[18], argv[19], argv[20],
            argv[21], argv[22], argv[23], argv[24], argv[25], argv[26],
            argv[27], argv[28], argv[29]);
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

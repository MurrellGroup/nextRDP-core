#pragma once

#include <array>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#if defined(_WIN32)
#define RDP_CALL __stdcall
#else
#define RDP_CALL
#endif

struct PreprocessApi {
    int(RDP_CALL* make_a_list_p2)(float, int, short*, int, short*);
    int(RDP_CALL* count_nucs)(int, int, int, short*, unsigned char*, int, int*);
    int(RDP_CALL* recode_nucs)(int, int, int, int*, int, unsigned char*);
    int(RDP_CALL* do_recode_p)(int, int, int, short*, int, unsigned char*, unsigned char*, int, unsigned char*);
    int(RDP_CALL* make_compress_seq_p)(int, int, unsigned char*, int, unsigned char*, int, int, unsigned char*);
};

template <typename T>
void write_json_array(std::ostream& output, const std::vector<T>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << static_cast<long long>(values[index]);
    }
    output << ']';
}

inline int run_preprocess_fixture(const PreprocessApi& api, std::ostream& output) {
    constexpr int next_no = 4;
    constexpr int sequence_length = 12;
    constexpr int sequence_stride = sequence_length + 1;
    constexpr int nucleotide_categories = 4;

    const std::array<std::string, next_no + 1> sequences = {
        "ACGTACGTACGT",
        "ACGTTCGTACGA",
        "GCGTACATACGT",
        "ACCTACGTTCGT",
        "TCGTACGTACGA",
    };

    std::vector<short> sequence_numbers(sequence_stride * (next_no + 1), 0);
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        for (int site = 1; site <= sequence_length; ++site) {
            sequence_numbers[site + sequence * sequence_stride] =
                static_cast<short>(sequences[sequence][site - 1]);
        }
    }

    std::vector<unsigned char> nucleotide_map(256, 0);
    nucleotide_map[static_cast<unsigned char>('A')] = 1;
    nucleotide_map[static_cast<unsigned char>('C')] = 2;
    nucleotide_map[static_cast<unsigned char>('G')] = 3;
    nucleotide_map[static_cast<unsigned char>('T')] = 4;

    std::vector<short> mask(next_no + 1, 0);
    std::vector<short> analysis_list(3 * 10, -1);
    const int last_triplet = api.make_a_list_p2(
        1.0F, next_no, mask.data(), 2, analysis_list.data());

    std::vector<int> nucleotide_counts(
        sequence_stride * (nucleotide_categories + 1), 0);
    const int count_result = api.count_nucs(
        next_no,
        sequence_length,
        sequence_length,
        sequence_numbers.data(),
        nucleotide_map.data(),
        sequence_length,
        nucleotide_counts.data());
    const std::vector<int> counted_nucleotides = nucleotide_counts;

    std::vector<unsigned char> replacements(
        sequence_stride * (nucleotide_categories + 1), 0);
    const int recode_result = api.recode_nucs(
        next_no,
        sequence_length,
        sequence_length,
        nucleotide_counts.data(),
        sequence_length,
        replacements.data());

    std::vector<unsigned char> recoded(sequence_stride * (next_no + 1), 0);
    const int do_recode_result = api.do_recode_p(
        next_no,
        sequence_length,
        sequence_length,
        sequence_numbers.data(),
        sequence_length,
        recoded.data(),
        nucleotide_map.data(),
        sequence_length,
        replacements.data());

    constexpr int compressor_bound = 4;
    std::vector<unsigned char> compressor(125, 0);
    for (int third = 0; third <= compressor_bound; ++third) {
        for (int second = 0; second <= compressor_bound; ++second) {
            for (int first = 0; first <= compressor_bound; ++first) {
                const int offset = first + second * 5 + third * 25;
                compressor[offset] = static_cast<unsigned char>(offset);
            }
        }
    }

    constexpr int compressed_bound = 4;
    std::vector<unsigned char> compressed(
        (compressed_bound + 1) * (next_no + 1), 0);
    const int compress_result = api.make_compress_seq_p(
        next_no,
        sequence_length,
        recoded.data(),
        compressed_bound,
        compressed.data(),
        compressor_bound,
        compressor_bound,
        compressor.data());

    output << '{';
    output << "\"lastTriplet\":" << last_triplet;
    output << ",\"returnCodes\":[" << count_result << ',' << recode_result << ','
           << do_recode_result << ',' << compress_result << ']';
    output << ",\"analysisList\":";
    write_json_array(output, analysis_list);
    output << ",\"countedNucleotides\":";
    write_json_array(output, counted_nucleotides);
    output << ",\"mutatedNucleotideCounts\":";
    write_json_array(output, nucleotide_counts);
    output << ",\"replacements\":";
    write_json_array(output, replacements);
    output << ",\"recoded\":";
    write_json_array(output, recoded);
    output << ",\"compressed\":";
    write_json_array(output, compressed);
    output << "}\n";
    return 0;
}


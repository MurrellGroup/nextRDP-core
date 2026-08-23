#include "scan_state.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> read_fasta_sequences(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open FASTA: " + path);
    std::vector<std::string> sequences;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line.front() == '>') {
            sequences.emplace_back();
            continue;
        }
        if (sequences.empty()) {
            if (line.find_first_not_of(" \t") != std::string::npos) {
                throw std::runtime_error("FASTA sequence precedes first header");
            }
            continue;
        }
        for (const unsigned char character : line) {
            if (character != ' ' && character != '\t') {
                sequences.back().push_back(
                    static_cast<char>(std::toupper(character)));
            }
        }
    }
    if (sequences.empty() || sequences.front().empty()) {
        throw std::runtime_error("FASTA contains no sequence data");
    }
    const auto length = sequences.front().size();
    for (const auto& sequence : sequences) {
        if (sequence.size() != length) {
            throw std::runtime_error("RDP requires an aligned equal-length FASTA");
        }
    }
    if (sequences.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("too many FASTA sequences");
    }
    return sequences;
}

short rdp_sequence_character(const char character) {
    switch (character) {
        case 'A': return 66;
        case 'C': return 68;
        case 'G': return 72;
        case 'T':
        case 'U': return 85;
        case '-':
        case '.': return 46;
        default: return 46;
    }
}

std::vector<unsigned char> make_compressor() {
    constexpr int upper_bound = 125;
    constexpr int stride = upper_bound + 1;
    std::vector<unsigned char> compressor(
        static_cast<std::size_t>(stride) * stride * stride, 0);
    int position = -1;
    for (int first = 0; first <= 4; ++first) {
        for (int second = 0; second <= 4; ++second) {
            for (int third = 0; third <= 4; ++third) {
                ++position;
                compressor[first + second * stride + third * stride * stride] =
                    static_cast<unsigned char>(position);
            }
        }
    }
    return compressor;
}

}  // namespace

RdpScanState build_rdp_scan_state_from_fasta(
    const std::string& fasta_path, const Dna5ScanPreprocessApi& api) {
    const auto fasta = read_fasta_sequences(fasta_path);
    RdpScanState state;
    state.next_no = static_cast<int>(fasta.size()) - 1;
    state.sequence_length = static_cast<int>(fasta.front().size());
    const int sequence_stride = state.sequence_length + 1;
    state.sequence_data.assign(
        static_cast<std::size_t>(sequence_stride) * (fasta.size() + 1), 0);
    for (int sequence = 0; sequence <= state.next_no; ++sequence) {
        for (int site = 1; site <= state.sequence_length; ++site) {
            state.sequence_data[site + sequence * sequence_stride] =
                rdp_sequence_character(fasta[sequence][site - 1]);
        }
    }

    std::array<unsigned char, 256> nucleotide_map{};
    nucleotide_map[66] = 1;
    nucleotide_map[68] = 2;
    nucleotide_map[72] = 3;
    nucleotide_map[85] = 4;
    nucleotide_map[46] = 0;
    const int column_stride = state.sequence_length + 1;
    std::vector<int> counts(static_cast<std::size_t>(column_stride) * 5, 0);
    if (api.count_nucleotides(
            state.next_no, state.sequence_length, state.sequence_length,
            state.sequence_data.data(), nucleotide_map.data(),
            state.sequence_length, counts.data()) != 1) {
        throw std::runtime_error("CountNucs failed");
    }
    std::vector<unsigned char> replacements(
        static_cast<std::size_t>(column_stride) * 5, 0);
    if (api.rank_nucleotides(
            state.next_no, state.sequence_length, state.sequence_length,
            counts.data(), state.sequence_length, replacements.data()) != 1) {
        throw std::runtime_error("RecodeNucs failed");
    }

    const int recoded_ub = state.sequence_length + 3;
    std::vector<unsigned char> recoded(
        static_cast<std::size_t>(recoded_ub + 1) * fasta.size(), 0);
    if (api.apply_recoding(
            state.next_no, state.sequence_length, state.sequence_length,
            state.sequence_data.data(), recoded_ub, recoded.data(),
            nucleotide_map.data(), state.sequence_length,
            replacements.data()) != 1) {
        throw std::runtime_error("DoRecodeP failed");
    }

    state.compressed_sequence_ub = static_cast<int>(std::nearbyint(
        static_cast<double>(recoded_ub) / 3.0 + 1.0));
    state.compressed_sequence.assign(
        static_cast<std::size_t>(state.compressed_sequence_ub + 1) *
            fasta.size(),
        0);
    auto compressor = make_compressor();
    if (api.compress_sequences(
            state.next_no, recoded_ub, recoded.data(),
            state.compressed_sequence_ub, state.compressed_sequence.data(),
            125, 125, compressor.data()) != 1) {
        throw std::runtime_error("MakeCompressSeqP failed");
    }

    const auto sequence_count = static_cast<std::uint64_t>(fasta.size());
    const auto triplet_count =
        sequence_count * (sequence_count - 1) * (sequence_count - 2) / 6;
    if (triplet_count >
        static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("too many sequence triplets");
    }
    std::vector<short> mask(fasta.size(), 0);
    state.analysis_list.assign(
        static_cast<std::size_t>(triplet_count) * 3, 0);
    state.analysis_list_last = api.make_analysis_list(
        1.0F, state.next_no, mask.data(), 2, state.analysis_list.data());
    if (state.analysis_list_last + 1 != static_cast<int>(triplet_count)) {
        throw std::runtime_error("MakeAListP2 returned an incomplete list");
    }
    return state;
}

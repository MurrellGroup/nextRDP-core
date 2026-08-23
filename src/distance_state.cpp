#include "distance_state.hpp"

#include "MathFuncsDll.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {

struct CategorySpec {
    int base;
    int sites_per_word;
    int table_upper_bound;
    bool gaps_are_missing;
};

constexpr CategorySpec category_14{5, 4, 625, true};
constexpr CategorySpec category_13{4, 5, 1024, true};
constexpr CategorySpec category_04{4, 5, 1024, false};
constexpr CategorySpec category_12{3, 6, 729, true};
constexpr CategorySpec category_03{3, 6, 729, false};
constexpr CategorySpec category_11{2, 10, 1024, true};
constexpr CategorySpec category_02{2, 10, 1024, false};

int power(int base, int exponent) {
    int value = 1;
    while (exponent-- > 0) value *= base;
    return value;
}

struct CompressionTable {
    CategorySpec spec{};
    std::vector<short> sequence_compressor;
    std::vector<char> valid;
    std::vector<char> differences;
};

CompressionTable make_compression_table(const CategorySpec spec) {
    CompressionTable table;
    table.spec = spec;
    const int states = power(spec.base, spec.sites_per_word);
    const int stride = spec.table_upper_bound + 1;
    table.sequence_compressor.assign(states, 0);
    table.valid.assign(static_cast<std::size_t>(stride) * stride, 0);
    table.differences.assign(static_cast<std::size_t>(stride) * stride, 0);

    // This is MakeSeqCompressors from Module4.bas, retaining its two distinct
    // orders: VB column-major subscript order and its lexicographic loop order.
    for (int position = 0; position < states; ++position) {
        int remainder = position;
        int vb_subscript = 0;
        std::array<int, 10> digits{};
        for (int site = spec.sites_per_word - 1; site >= 0; --site) {
            digits[site] = remainder % spec.base;
            remainder /= spec.base;
        }
        for (int site = 0; site < spec.sites_per_word; ++site) {
            vb_subscript += digits[site] * power(spec.base, site);
        }
        table.sequence_compressor[vb_subscript] =
            static_cast<short>(position);
    }

    for (int first = 0; first < states; ++first) {
        for (int second = 0; second < states; ++second) {
            int first_remainder = first;
            int second_remainder = second;
            int valid = 0;
            int differences = 0;
            for (int site = spec.sites_per_word - 1; site >= 0; --site) {
                const int first_digit = first_remainder % spec.base;
                const int second_digit = second_remainder % spec.base;
                first_remainder /= spec.base;
                second_remainder /= spec.base;
                const bool comparable = !spec.gaps_are_missing ||
                    (first_digit != 0 && second_digit != 0);
                if (comparable) {
                    ++valid;
                    if (first_digit != second_digit) ++differences;
                }
            }
            const auto offset = static_cast<std::size_t>(first) +
                static_cast<std::size_t>(second) * stride;
            table.valid[offset] = static_cast<char>(valid);
            table.differences[offset] = static_cast<char>(differences);
        }
    }
    return table;
}

struct CompressionTables {
    CompressionTable c14 = make_compression_table(category_14);
    CompressionTable c13 = make_compression_table(category_13);
    CompressionTable c04 = make_compression_table(category_04);
    CompressionTable c12 = make_compression_table(category_12);
    CompressionTable c03 = make_compression_table(category_03);
    CompressionTable c11 = make_compression_table(category_11);
    CompressionTable c02 = make_compression_table(category_02);
};

CompressionTables& compression_tables() {
    static CompressionTables tables;
    return tables;
}

// Exact scalar body of legacy DNA.dll!MakeSeqCatCount2, expressed with
// descriptive names but retaining its indexing, mutation, and wrap behavior.
void make_sequence_category_counts(
    int next_no, int sequence_length, int sequence_ub,
    int start_position, int end_position, const short* sequence_data,
    std::array<int, 10>& counts, std::vector<unsigned char>& flp,
    std::vector<unsigned char>& missing_level,
    std::vector<unsigned char>& nucleotide_level) {
    std::array<int, 5> present{};
    std::array<unsigned char, 256> nucleotide_map{};
    std::array<unsigned char, 5> nucleotide_code{};
    std::array<unsigned char, 256> recode{};
    nucleotide_map[66] = 1;
    nucleotide_map[68] = 2;
    nucleotide_map[72] = 3;
    nucleotide_map[85] = 4;
    nucleotide_code[1] = 66;
    nucleotide_code[2] = 68;
    nucleotide_code[3] = 72;
    nucleotide_code[4] = 85;
    counts.fill(0);

    const int sequence_stride = sequence_ub + 1;
    const int flp_stride = next_no + 1;
    const auto visit = [&](const int position) {
        present.fill(0);
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            present[nucleotide_map[sequence_data[
                position + sequence * sequence_stride]]] = 1;
        }
        const int nucleotide_types =
            present[1] + present[2] + present[3] + present[4];
        int next_recode = present[0];
        ++counts[present[0] + nucleotide_types * 2];
        for (int nucleotide = 1; nucleotide <= 4; ++nucleotide) {
            if (present[nucleotide] == 1) {
                recode[nucleotide_code[nucleotide]] =
                    static_cast<unsigned char>(next_recode++);
            }
        }
        missing_level[position] = static_cast<unsigned char>(present[0]);
        nucleotide_level[position] =
            static_cast<unsigned char>(nucleotide_types);
        for (int sequence = 0; sequence <= next_no; ++sequence) {
            flp[sequence + position * flp_stride] = recode[sequence_data[
                position + sequence * sequence_stride]];
        }
    };

    if (start_position <= end_position) {
        for (int position = start_position; position <= end_position;
             ++position) {
            visit(position);
        }
    } else {
        for (int position = start_position; position <= sequence_length;
             ++position) {
            visit(position);
        }
        for (int position = 1; position <= end_position; ++position) {
            visit(position);
        }
    }
}

struct PackedCategory {
    CategorySpec spec{};
    int site_count = 0;
    int word_ub = 0;
    std::vector<unsigned char> sites;
    std::vector<short> words;
};

PackedCategory make_category(
    const CategorySpec spec, int site_count, int next_no) {
    PackedCategory category;
    category.spec = spec;
    category.site_count = site_count;
    category.sites.assign(
        static_cast<std::size_t>(site_count + 1) * (next_no + 2), 0);
    category.word_ub = static_cast<int>(std::nearbyint(
        static_cast<double>(site_count) / spec.sites_per_word)) + 1;
    category.words.assign(
        static_cast<std::size_t>(category.word_ub + 1) * (next_no + 1), 0);
    return category;
}

void fill_category_site(
    PackedCategory& category, int category_position, int source_position,
    int next_no, const std::vector<unsigned char>& flp) {
    const int category_stride = category.site_count + 1;
    const int flp_stride = next_no + 1;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        category.sites[category_position + sequence * category_stride] =
            flp[sequence + source_position * flp_stride];
    }
}

void pack_category(
    PackedCategory& category, int next_no,
    const CompressionTable& table) {
    const int site_stride = category.site_count + 1;
    const int word_stride = category.word_ub + 1;
    for (int sequence = 0; sequence <= next_no; ++sequence) {
        int word_position = 0;
        int site_position = 1;
        for (; site_position <= category.site_count - category.spec.sites_per_word;
             site_position += category.spec.sites_per_word) {
            ++word_position;
            int compressor_index = 0;
            int multiplier = 1;
            for (int site = 0; site < category.spec.sites_per_word; ++site) {
                compressor_index += category.sites[
                    site_position + site + sequence * site_stride] * multiplier;
                multiplier *= category.spec.base;
            }
            category.words[word_position + sequence * word_stride] =
                table.sequence_compressor[compressor_index];
        }

        ++word_position;
        int compressor_index = 0;
        int multiplier = 1;
        for (int site = 0; site < category.spec.sites_per_word; ++site) {
            if (site_position + site <= category.site_count) {
                compressor_index += category.sites[
                    site_position + site + sequence * site_stride] * multiplier;
            }
            multiplier *= category.spec.base;
        }
        category.words[word_position + sequence * word_stride] =
            table.sequence_compressor[compressor_index];
    }
}

}  // namespace

RdpDistanceState build_rdp_distance_state(
    const RdpScanState& scan_state, const int start_position,
    const int end_position) {
    if (scan_state.next_no < 0 || scan_state.sequence_length < 1) {
        throw std::runtime_error("distance calculation requires sequences");
    }
    if (start_position < 1 || start_position > scan_state.sequence_length ||
        end_position < 1 || end_position > scan_state.sequence_length) {
        throw std::runtime_error("distance interval is outside the alignment");
    }

    const int next_no = scan_state.next_no;
    const int sequence_count = next_no + 1;
    std::array<int, 10> counts{};
    std::vector<unsigned char> flp(
        static_cast<std::size_t>(scan_state.sequence_length + 1) *
            sequence_count,
        0);
    std::vector<unsigned char> missing_level(scan_state.sequence_length + 1, 0);
    std::vector<unsigned char> nucleotide_level(
        scan_state.sequence_length + 1, 0);
    make_sequence_category_counts(
        next_no, scan_state.sequence_length, scan_state.sequence_length,
        start_position, end_position, scan_state.sequence_data.data(), counts,
        flp, missing_level, nucleotide_level);

    PackedCategory c11 = make_category(category_11, counts[3], next_no);
    PackedCategory c02 = make_category(category_02, counts[4], next_no);
    PackedCategory c12 = make_category(category_12, counts[5], next_no);
    PackedCategory c03 = make_category(category_03, counts[6], next_no);
    PackedCategory c13 = make_category(category_13, counts[7], next_no);
    PackedCategory c04 = make_category(category_04, counts[8], next_no);
    PackedCategory c14 = make_category(category_14, counts[9], next_no);

    std::array<int, 10> category_positions{};
    int position = start_position;
    while (true) {
        const int category_index =
            missing_level[position] + nucleotide_level[position] * 2;
        if (category_index >= 3 && category_index <= 9) {
            ++category_positions[category_index];
            PackedCategory* target = nullptr;
            switch (category_index) {
                case 3: target = &c11; break;
                case 4: target = &c02; break;
                case 5: target = &c12; break;
                case 6: target = &c03; break;
                case 7: target = &c13; break;
                case 8: target = &c04; break;
                case 9: target = &c14; break;
            }
            fill_category_site(
                *target, category_positions[category_index], position,
                next_no, flp);
        }
        if (position == end_position) break;
        ++position;
        if (position > scan_state.sequence_length) position = 1;
    }

    auto& tables = compression_tables();
    pack_category(c14, next_no, tables.c14);
    pack_category(c13, next_no, tables.c13);
    pack_category(c04, next_no, tables.c04);
    pack_category(c12, next_no, tables.c12);
    pack_category(c03, next_no, tables.c03);
    pack_category(c11, next_no, tables.c11);
    pack_category(c02, next_no, tables.c02);

    RdpDistanceState result;
    const auto matrix_size =
        static_cast<std::size_t>(sequence_count) * sequence_count;
    result.differences.assign(matrix_size, 0.0F);
    result.valid_sites.assign(matrix_size, 0.0F);
    result.distance.assign(matrix_size, 0.0F);
    result.redo_distance.assign(sequence_count, 1);
    // ShowProg=0 in the observed initial path, so FastDistanceCalcZ calls
    // SuperDistP once for each X instead of taking the SuperDistP2 branch.
    for (int sequence = 0; sequence < next_no; ++sequence) {
        const double upper = MathFuncs::MyMathFuncs::SuperDistP(
            sequence, next_no, c14.word_ub, c04.word_ub, c13.word_ub,
            c03.word_ub, c12.word_ub, c02.word_ub, c11.word_ub,
            &result.average_distance_accumulator, result.differences.data(),
            result.valid_sites.data(), result.distance.data(),
            result.redo_distance.data(), counts.data(), c14.words.data(),
            c04.words.data(), c13.words.data(), c03.words.data(),
            c12.words.data(), c02.words.data(), c11.words.data(),
            tables.c14.valid.data(), tables.c14.differences.data(),
            tables.c13.valid.data(), tables.c13.differences.data(),
            tables.c12.valid.data(), tables.c12.differences.data(),
            tables.c11.valid.data(), tables.c11.differences.data(),
            tables.c04.differences.data(), tables.c03.differences.data(),
            tables.c02.differences.data());
        if (result.upper_distance < upper) result.upper_distance = upper;
    }

    // CalcDistances performs this immediately after FastDistanceCalcZ. It is
    // intentionally not folded into the distance kernel because RDP does not.
    int cutoff = scan_state.sequence_length / 100;
    if (cutoff < 30) cutoff = 30;
    if (cutoff > 50) cutoff = 50;
    for (int first = 0; first <= next_no; ++first) {
        result.distance[first + first * sequence_count] = 1.0F;
        for (int second = first + 1; second <= next_no; ++second) {
            const auto first_offset = first + second * sequence_count;
            const auto second_offset = second + first * sequence_count;
            if (result.valid_sites[first_offset] < cutoff) {
                result.valid_sites[first_offset] = 0.0F;
                result.valid_sites[second_offset] = 0.0F;
                result.differences[first_offset] = 0.0F;
                result.differences[second_offset] = 0.0F;
                result.distance[first_offset] = 0.0F;
                result.distance[second_offset] = 0.0F;
            }
        }
    }
    return result;
}

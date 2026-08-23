#pragma once

#include <cstdint>
#include <cstring>
#include <ostream>
#include <vector>

#if defined(_WIN32)
#define RDP_DISTANCE_CALL __stdcall
#else
#define RDP_DISTANCE_CALL
#endif

using SuperDistPFunction = double(RDP_DISTANCE_CALL*)(
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    int,
    double*,
    float*,
    float*,
    float*,
    short*,
    int*,
    short*,
    short*,
    short*,
    short*,
    short*,
    short*,
    short*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*,
    char*);

template <typename To, typename From>
To copy_bits(const From& value) {
    static_assert(sizeof(To) == sizeof(From));
    To result{};
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

inline void write_float_bits(std::ostream& output, const std::vector<float>& values) {
    output << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << copy_bits<std::uint32_t>(values[index]);
    }
    output << ']';
}

inline int run_distance_fixture(SuperDistPFunction super_dist_p, std::ostream& output) {
    constexpr int next_no = 4;
    constexpr int upper_bound = 2;
    constexpr int stride = upper_bound + 1;
    constexpr int matrix_stride = next_no + 1;

    std::vector<short> category14(stride * matrix_stride, 0);
    std::vector<short> category04(stride * matrix_stride, 0);
    std::vector<short> category13(stride * matrix_stride, 0);
    std::vector<short> category03(stride * matrix_stride, 0);
    std::vector<short> category12(stride * matrix_stride, 0);
    std::vector<short> category02(stride * matrix_stride, 0);
    std::vector<short> category11(stride * matrix_stride, 0);

    for (int sequence = 0; sequence <= next_no; ++sequence) {
        for (int site = 1; site <= upper_bound; ++site) {
            const int offset = site + sequence * stride;
            category14[offset] = static_cast<short>((sequence + site) % 4);
            category04[offset] = static_cast<short>((sequence + site * 2) % 4);
            category13[offset] = static_cast<short>((sequence * 2 + site) % 4);
            category03[offset] = static_cast<short>((sequence * 3 + site) % 4);
            category12[offset] = static_cast<short>((sequence + site * 3) % 4);
            category02[offset] = static_cast<short>((sequence * 2 + site * 3) % 4);
            category11[offset] = static_cast<short>((sequence * 3 + site * 2) % 4);
        }
    }

    const auto make_lookup = [](int categories, bool valid) {
        std::vector<char> lookup(categories * categories, 0);
        for (int second = 0; second < 4; ++second) {
            for (int first = 0; first < 4; ++first) {
                lookup[first + second * categories] = static_cast<char>(
                    valid ? 1 : (first == second ? 0 : 1));
            }
        }
        return lookup;
    };

    auto valid14 = make_lookup(626, true);
    auto diff14 = make_lookup(626, false);
    auto valid13 = make_lookup(1025, true);
    auto diff13 = make_lookup(1025, false);
    auto valid12 = make_lookup(730, true);
    auto diff12 = make_lookup(730, false);
    auto valid11 = make_lookup(1025, true);
    auto diff11 = make_lookup(1025, false);
    auto diff04 = make_lookup(1025, false);
    auto diff03 = make_lookup(730, false);
    auto diff02 = make_lookup(1025, false);

    std::vector<short> redo_distance(matrix_stride, 1);
    std::vector<int> sequence_category_count(9, 0);
    sequence_category_count[2] = upper_bound;
    sequence_category_count[4] = upper_bound;
    sequence_category_count[6] = upper_bound;
    sequence_category_count[8] = upper_bound;

    std::vector<float> pair_differences(matrix_stride * matrix_stride, 0.0F);
    std::vector<float> pair_valid(matrix_stride * matrix_stride, 0.0F);
    std::vector<float> pair_distance(matrix_stride * matrix_stride, 0.0F);
    std::vector<std::uint64_t> returned_upper_bits;
    double average_distance_sum = 0.0;

    for (int sequence = 0; sequence < next_no; ++sequence) {
        const double upper = super_dist_p(
            sequence,
            next_no,
            upper_bound,
            upper_bound,
            upper_bound,
            upper_bound,
            upper_bound,
            upper_bound,
            upper_bound,
            &average_distance_sum,
            pair_differences.data(),
            pair_valid.data(),
            pair_distance.data(),
            redo_distance.data(),
            sequence_category_count.data(),
            category14.data(),
            category04.data(),
            category13.data(),
            category03.data(),
            category12.data(),
            category02.data(),
            category11.data(),
            valid14.data(),
            diff14.data(),
            valid13.data(),
            diff13.data(),
            valid12.data(),
            diff12.data(),
            valid11.data(),
            diff11.data(),
            diff04.data(),
            diff03.data(),
            diff02.data());
        returned_upper_bits.push_back(copy_bits<std::uint64_t>(upper));
    }

    output << '{';
    output << "\"averageDistanceBits\":"
           << copy_bits<std::uint64_t>(average_distance_sum);
    output << ",\"returnedUpperBits\":[";
    for (std::size_t index = 0; index < returned_upper_bits.size(); ++index) {
        if (index != 0) {
            output << ',';
        }
        output << returned_upper_bits[index];
    }
    output << ']';
    output << ",\"pairDifferenceBits\":";
    write_float_bits(output, pair_differences);
    output << ",\"pairValidBits\":";
    write_float_bits(output, pair_valid);
    output << ",\"pairDistanceBits\":";
    write_float_bits(output, pair_distance);
    output << "}\n";
    return 0;
}

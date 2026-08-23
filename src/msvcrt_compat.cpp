#include <cstdint>

#if !defined(_WIN32)
namespace {
std::uint32_t msvcrt_random_state = 1;
}

extern "C" void srand(unsigned int seed) {
    msvcrt_random_state = seed;
}

extern "C" int rand() {
    msvcrt_random_state =
        msvcrt_random_state * UINT32_C(214013) + UINT32_C(2531011);
    return static_cast<int>((msvcrt_random_state >> 16U) & UINT32_C(0x7fff));
}
#endif

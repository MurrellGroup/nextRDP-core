#pragma once

inline int __popcnt(unsigned int value) {
    return __builtin_popcount(value);
}


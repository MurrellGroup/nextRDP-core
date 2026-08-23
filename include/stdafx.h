#pragma once

#include <cstring>

#include "rdp/platform_compat.hpp"

// The pinned Windows oracle links the MSVCRT pseudo-random generator.  RDP's
// tree code also uses RAND_MAX in its scaling, so both parts of that ABI are
// required when the unchanged source is compiled for another platform.
#if !defined(_WIN32)
#undef RAND_MAX
#define RAND_MAX 32767
#endif

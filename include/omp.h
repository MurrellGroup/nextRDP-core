#pragma once

// The source DLL uses OpenMP for the independent per-triplet method passes.
// Keep the tiny single-thread fallback for Emscripten builds that do not opt
// into pthreads, but let native and pthread-enabled builds include the real
// runtime header.  `include_next` is intentional: this compatibility header
// is first on the target's include path, while GCC/Clang provide the actual
// omp.h later in their system paths.
#if defined(NEXT_RDP_USE_REAL_OPENMP)
#  include_next <omp.h>
#else
inline int omp_get_num_procs() { return 1; }
inline int omp_get_thread_num() { return 0; }
inline void omp_set_num_threads(int) {}
#endif

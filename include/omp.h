#pragma once

inline int omp_get_num_procs() { return 1; }
inline int omp_get_thread_num() { return 0; }
inline void omp_set_num_threads(int) {}

#pragma once

#include "scan_state.hpp"

#include <array>
#include <vector>

// The BURT implementation is kept separate from the event allocator because
// PolishBP/BenHMM operate on one three-sequence event at a time.  Coordinates
// are the one-based coordinates used by the VB source.
struct RdpBurtInterval {
    std::array<int, 5> values{}; // 99% left/right, HMM point, 95% left/right
};

struct RdpBurtResult {
    std::array<int, 10> confidence{};
    std::vector<RdpBurtInterval> intervals;
    int input_beginning = 0;
    int input_ending = 0;
    int polished_beginning = 0;
    int polished_ending = 0;
    int information_rich_sites = 0;
    double best_log_likelihood = -1000000.0;
    bool attempted = false;
    bool trained = false;
    bool available = false;
    bool single_transition_assignment = false;
    bool insufficient_inside_or_outside_reverted = false;
};

// Source Module4.BenHMM followed by Module4.PolishBP(HMMCycles,0,...).
// This is intentionally a direct source-order routine: it uses the exact
// DNA5 HMM helper bodies linked into next_rdp_core rather than a browser-side
// statistical approximation.
RdpBurtResult run_rdp_burt(
    const RdpScanState& scan_state,
    const std::array<int, 3>& representatives,
    int beginning,
    int ending,
    bool circular = true,
    int hmm_cycles = 20,
    int repos_flag = 0);


#pragma once

#include "scan_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

struct GeneconvEmissionTrace {
    std::array<int, 3> input{};
    std::array<int, 3> counts{};
    std::array<double, 3> store_lpv{};
    std::array<int, 3> active{};
    RdpRawEvent event{};
};

struct MaxchiPeakTrace {
    int repetition = 0;
    int maximum_position = 0;
    int comparison = 0;
    int test_window = 0;
    int initial_left = 0;
    int initial_right = 0;
    int best_window = 0;
    int top_left = 0;
    int top_right = 0;
    int top_left_position = 0;
    int top_right_position = 0;
    double initial_maximum = 0.0;
    double grown_maximum = 0.0;
    double probability = 1.0;
    std::uint64_t difference_position_hash = 0;
    std::uint64_t score_hash = 0;
    std::array<std::uint64_t, 3> score_plane_hash{};
    std::uint64_t chi_table_hash = 0;
    std::uint64_t chi_map_hash = 0;
    bool accepted = false;
};

// State shared by the methods called from Module2.FinalTrim.  In VB6 the
// methods all append through UpdateXOList3 into one rectangular XOverList;
// its second-dimension upper bound is observable (TSXOver has a strict-bound
// copy check), so it cannot be represented by independent growable rows.
class RdpLegacyEventAllocator {
public:
    RdpLegacyEventAllocator(RdpRawEventState& events, int selected_program_count,
                            int max_events = 32767);

    int allocate(int active_sequence, int program, double probability);
    bool has_strictly_later_slot(int slot) const;
    RdpRawEvent& event(int sequence, int slot);
    int count(int sequence) const;
    int upper_bound() const { return upper_bound_; }

private:
    RdpRawEventState& events_;
    int selected_program_count_;
    static constexpr int method_count_ = 9;
    int max_events_;
    int upper_bound_ = 10;
    std::vector<int> max_xop_;
};

// Literal computational path of Module30.GCXoverD(1), with the plotting,
// UI, redo-list bookkeeping and missing-data split branch omitted because
// Module2.FinalTrim sets DontWorryAboutSplitsFlag=1 before this call.
void run_rdp_geneconv_recheck(
    const RdpScanState& scan_state, std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, bool long_winded = true,
    std::vector<GeneconvEmissionTrace>* trace = nullptr);

// Module5.MCXoverF(1), the MaxChi recheck used by FinalTrim.  This is the
// enumerating VB path (not DNA5.FastRecCheckMC2, which deliberately returns
// after its first significant peak).
void run_rdp_maxchi_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, int event_beginning, int event_ending,
    bool initial_scan = false,
    std::vector<MaxchiPeakTrace>* trace = nullptr);

// Module30.CXoverA(2), the Chimaera recheck used in each of FinalTrim's three
// source-order role rotations.
void run_rdp_chimaera_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, int event_beginning, int event_ending);

// Literal allocation/order portion of Module3.TSXOver(1).  The two tested
// excursions are emitted in source order and each accepted event requests
// the source's reverse companion slot.
void run_rdp_three_seq_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    const std::vector<float>& probability_table, int table_bound,
    RdpLegacyEventAllocator& allocator);

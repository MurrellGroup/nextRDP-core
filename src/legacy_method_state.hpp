#pragma once

#include "scan_state.hpp"
#include "xover_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

// Source-order screening result from DNA5.AlistGC2.  The compatibility
// runner keeps this separate from event emission: the VB scanner first builds
// an AList, calls AlistGC2, then invokes GCXoverD only for RL==1 entries.
struct RdpMethodScreenResult {
    std::vector<std::array<int, 3>> candidates;
    std::vector<unsigned char> redo;
};

// Literal DNA5.MakeAListISP3 scheduling for one InnerScan2 method pass.
// The method kernels do not scan the permanent AnalysisList directly once an
// event has been selected: they expand the selected WinPP RList through the
// source's TraceSub/ActualSeqSize/DoPairs gates first.  Keeping this routine
// separate makes that source-order boundary testable for every method.
std::vector<std::array<int, 3>> make_rdp_inner_method_triplets(
    const RdpScanState& scan_state, const std::array<int, 3>& rnum,
    const std::vector<int>& rlist, int win_pp,
    const std::vector<int>& trace_sub,
    const std::vector<int>& actual_sequence_sizes,
    const std::vector<unsigned char>& do_pairs,
    int permanent_next_no, int min_sequence_size, int method_program,
    int selected_program_bits, float probability_step = 1.1F);

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
    std::uint64_t position_difference_hash = 0;
    std::uint64_t score_hash = 0;
    std::array<std::uint64_t, 3> score_plane_hash{};
    std::uint64_t chi_table_hash = 0;
    std::uint64_t chi_map_hash = 0;
    bool accepted = false;
    int raw_beginning = 0;
    int raw_ending = 0;
    int first_centered_beginning = 0;
    int first_centered_ending = 0;
    int polished_beginning = 0;
    int polished_ending = 0;
    int polish_begin_index = 0;
    int polish_left_stop_index = 0;
    int polish_begin_candidate_index = 0;
    int polish_begin_candidate = 0;
    bool polish_begin_accepted = false;
    int opt_initial_score = 0;
    int opt_first_score = 0;
    int opt_first_position = 0;
    bool opt_first_improved = false;
    int centered_beginning = 0;
    int centered_ending = 0;
    int destroy_left = 0;
    int destroy_right = 0;
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

// Literal DNA5.AlistGC2 screening over the supplied analysis list.  This is
// deliberately only a scanner: source scheduling (RList/Worthwhilescan and
// MakeAListISP3) is layered above it by the analysis loop.
RdpMethodScreenResult screen_rdp_geneconv_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    int store_lpv_ub, int correction_tests, double lowest_probability,
    int circular, int mc_flag = 0, int target = 0);

RdpMethodScreenResult screen_rdp_maxchi_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    int store_lpv_ub, int correction_tests, double lowest_probability,
    int circular, int mc_flag = 0, int event_number = 0,
    const std::vector<unsigned char>* missing_data = nullptr);

RdpMethodScreenResult screen_rdp_chimaera_candidates(
    const RdpScanState& scan_state, const std::vector<double>& store_lpv,
    int store_lpv_ub, int correction_tests, double lowest_probability,
    int circular, int mc_flag = 0, int event_number = 0,
    const std::vector<unsigned char>* missing_data = nullptr);

// Module5.MCXoverF(1), the MaxChi recheck used by FinalTrim.  This is the
// enumerating VB path (not DNA5.FastRecCheckMC2, which deliberately returns
// after its first significant peak).
void run_rdp_maxchi_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, int event_beginning, int event_ending,
    bool initial_scan = false,
    std::vector<MaxchiPeakTrace>* trace = nullptr,
    bool inner_scan = false,
    const std::vector<unsigned char>* missing_data = nullptr);

// Module30.CXoverA(2), the Chimaera recheck used in each of FinalTrim's three
// source-order role rotations.
void run_rdp_chimaera_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    RdpLegacyEventAllocator& allocator, int event_beginning, int event_ending,
    bool initial_scan = false, int* shared_xdiffpos0 = nullptr,
    bool inner_scan = false,
    const std::vector<unsigned char>* missing_data = nullptr);

// Literal allocation/order portion of Module3.TSXOver(1).  The two tested
// excursions are emitted in source order and each accepted event requests
// the source's reverse companion slot.
void run_rdp_three_seq_recheck(
    const RdpScanState& scan_state, const std::array<int, 3>& sequences,
    const std::vector<double>& store_lpv, int store_lpv_ub,
    const RdpProbabilitySettings& settings,
    const std::vector<float>& probability_table, int table_bound,
    RdpLegacyEventAllocator& allocator);

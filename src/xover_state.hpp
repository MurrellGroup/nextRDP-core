#pragma once

#include "distance_state.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"

#include <array>
#include <cstdint>
#include <vector>

#if defined(_WIN32)
#define RDP_XOVER_CALL __stdcall
#else
#define RDP_XOVER_CALL
#endif

struct Dna5XoverApi {
    int(RDP_XOVER_CALL* find_subsequence_plain)(
        int*, short, short, int, int, int, int, int, int, short, short*,
        short*, char*, short*, int*, int*, short*);
    int(RDP_XOVER_CALL* find_subsequence)(
        int*, int, int, int, int, int, int, int, int, unsigned char*, int,
        char*, unsigned char*);
    int(RDP_XOVER_CALL* calculate_homology)(
        short, int, int, short, char*, int*);
    int(RDP_XOVER_CALL* find_next)(
        int, int, int, int, int, int, int, int*);
    int(RDP_XOVER_CALL* find_first)(int, int, int, int, int, int*);
    int(RDP_XOVER_CALL* define_event)(
        int, int, int, int, int, int, int, int, int, int, int, int, int,
        int, int*, int*, int*, int*, int*, char*, int*);
    double(RDP_XOVER_CALL* probability_p2)(
        double*, int, int, int, double, int);
    double(RDP_XOVER_CALL* probability_p)(double*, int, int, double, int);
    int(RDP_XOVER_CALL* find_subsequence_with_positions)(
        int*, int, int, int, int, int, int, int, int, unsigned char*, int,
        char*, int*, int*, unsigned char*);
    int(RDP_XOVER_CALL* clean_xover_sequence)(int, int, int, char*);
};

struct RdpXoverSettings {
    int short_output = 0;
    int long_winded = 0;
    int target = 0;
    int circular = 0;
};

struct RdpProbabilitySettings {
    int circular = 0;
    int mc_correction = 1;
    int mc_flag = 0;
    int probability_file_flag = 0;
    int probability_one_ub = 0;
    int probability_two_ub = 0;
    int fact_three_ub = 0;
    double lowest_probability = 0.0;
};

// Optional observer used by the source-faithful scan wrapper.  This is kept
// in the xover header (rather than depending on analysis.hpp) so the low-level
// triplet walker can report real work while a vendored DNA5 call is running.
using RdpXoverProgressCallback = void (*) (
    int phase, int round, int processed_triplets, int total_triplets,
    int event_count, void* user);

struct RdpFirstXoverState {
    std::array<int, 3> sequences{};
    std::array<int, 3> agreement_counts{};
    int informative_length = 0;
    int homology_length = 0;
    int homology_sequence_length = 0;
    int initial_high_homology = 0;
    int high_homology = 0;
    int med_homology = 0;
    int low_homology = 0;
    int homology_start = 0;
    int next_position = -1;
    int define_input_position = -1;
    int event_position = -1;
    int old_find_position = -1;
    int find_cycle = 0;
    int end_flag = 0;
    int event_begin = 0;
    int event_end = 0;
    int number_in_common = 0;
    int event_length = 0;
    bool used_find_first = false;
    bool probability_tested = false;
    bool used_probability_p2 = false;
    int probability_length = 0;
    int probability_same = 0;
    int probability_different = 0;
    double probability_scale = 1.0;
    double individual_probability = 0.0;
    double probability_prefilter_value = 0.0;
    double event_probability = 0.0;
    double adjusted_event_probability = 0.0;
    bool significant_event = false;
    int position_map_result = 0;
    int active_sequence = -1;
    int active_major_parent = -1;
    int active_minor_parent = -1;
    int sequence_daughter = -1;
    int sequence_minor = -1;
    int xover_sequence_ub = 0;
    int homology_ub = 0;
    std::array<double, 3> average_homology{};
    std::vector<char> xover_sequence;
    std::vector<int> homology;
    std::vector<char> xover_sequence_at_define;
    std::vector<int> homology_at_define;
    std::vector<int> xdiffpos;
    std::vector<int> xposdiff;
};

struct RdpRawEvent {
    std::uint8_t outside_flag = 0;
    std::uint8_t misidentify_flag = 0;
    std::uint8_t program_flag = 0;
    std::uint8_t sbp_flag = 0;
    std::uint8_t accept = 0;
    std::int16_t major_parent = 0;
    std::int16_t minor_parent = 0;
    std::int16_t daughter = 0;
    std::int32_t beginning = 0;
    std::int32_t ending = 0;
    std::int32_t length_holder = 0;
    std::int32_t event_number = 0;
    float permutation_pvalue = 0.0F;
    std::int32_t begin_parent = 0;
    std::int32_t end_parent = 0;
    double probability = 0.0;
    double distance_holder = 0.0;
    // The source keeps the triplet and whether XOver used the compressed
    // FindSubSeqPB3 path implicit in the live buffers.  Retain that small
    // bit of provenance with each RDP record so the review plot can replay
    // the same XOverHomologyP profile after the scan has moved on.
    std::array<std::int16_t, 3> profile_sequences{};
    std::uint8_t profile_available = 0;
    std::uint8_t profile_use_compress = 0;
};

struct RdpRawEventState {
    std::vector<std::int16_t> current_xover;
    std::vector<std::vector<RdpRawEvent>> xover_list;
    int scanned_triplets = 0;
    int significant_candidates = 0;
    std::vector<unsigned char> triplets_with_events;
};

RdpFirstXoverState build_rdp_first_xover_state(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, int triplet_index, int fss_ub,
    std::vector<unsigned char>& fss_rdp, int xover_window,
    short xover_window_x, const Dna5XoverApi& api,
    const std::array<int, 3>* explicit_sequences = nullptr,
    int sequence_event_number = 0,
    bool use_compress = false);

void define_rdp_first_xover_event(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    int xover_window, const RdpXoverSettings& settings,
    const Dna5XoverApi& api);

void calculate_rdp_first_xover_probability(
    RdpFirstXoverState& state, const RdpProbabilitySettings& settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api);

void continue_rdp_xover_to_first_probability(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    int xover_window, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api);

bool apply_rdp_probability_cutoff(
    RdpFirstXoverState& state, const RdpProbabilitySettings& settings);

void build_rdp_first_position_maps(
    RdpFirstXoverState& state, const RdpScanState& scan_state, int fss_ub,
    int xover_window, std::vector<unsigned char>& fss_rdp,
    const Dna5XoverApi& api);

bool advance_rdp_role_cycle(RdpFirstXoverState& state, int do_all);

void scan_rdp_current_roles_to_first_probability(
    RdpFirstXoverState& state, const RdpScanState& scan_state,
    int xover_window, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api);

RdpRawEventState scan_rdp_redo_triplets(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, const std::vector<unsigned char>& redo,
    std::vector<unsigned char>& fss_rdp, const std::vector<double>& store_lpv,
    int store_lpv_ub, int fss_ub, int xover_window,
    short xover_window_x, const RdpXoverSettings& xover_settings,
    const RdpProbabilitySettings& probability_settings,
    std::vector<double>& probability_estimate,
    std::vector<double>& fact_three, std::vector<double>& fact,
    const Dna5XoverApi& api, int do_all = 0,
    const RdpRawEventState* initial_events = nullptr,
    const std::array<int, 3>* explicit_sequences = nullptr,
    int sequence_event_number = 0,
    const std::vector<unsigned char>* missing_data = nullptr,
    bool use_compress = false,
    int* shared_xdiffpos0 = nullptr,
    RdpXoverProgressCallback progress_callback = nullptr,
    void* progress_user = nullptr,
    int progress_phase = 0,
    int progress_round = 1);

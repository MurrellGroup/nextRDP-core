#pragma once

#include "distance_state.hpp"
#include "scan_state.hpp"
#include "tree_state.hpp"

#include <array>
#include <vector>

#if defined(_WIN32)
#define RDP_XOVER_CALL __stdcall
#else
#define RDP_XOVER_CALL
#endif

struct Dna5XoverApi {
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

RdpFirstXoverState build_rdp_first_xover_state(
    const RdpScanState& scan_state, const RdpDistanceState& distance_state,
    const RdpTreeState& tree_state, int triplet_index, int fss_ub,
    std::vector<unsigned char>& fss_rdp, int xover_window,
    short xover_window_x, const Dna5XoverApi& api);

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

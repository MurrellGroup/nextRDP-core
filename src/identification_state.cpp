#include "identification_state.hpp"

#include <stdexcept>

RdpBreakpointFlanks make_rdp_breakpoint_flanks(
    const RdpScanState& scan_state, const int beginning, const int ending,
    const std::array<int, 3>& sequences, const int variable_site_target,
    const int total_site_target) {
    if (beginning < 1 || beginning > scan_state.sequence_length ||
        ending < 1 || ending > scan_state.sequence_length) {
        throw std::runtime_error("MakeBPosLR interval is outside alignment");
    }
    for (const int sequence : sequences) {
        if (sequence < 0 || sequence > scan_state.next_no) {
            throw std::runtime_error("MakeBPosLR sequence is outside alignment");
        }
    }

    RdpBreakpointFlanks result;
    const int length = scan_state.sequence_length;
    const int stride = length + 1;
    const std::array<int, 4> increments{-1, 1, -1, 1};
    for (int flank = 0; flank < 4; ++flank) {
        int position;
        int stop;
        if (flank == 0) {
            position = beginning - 1;
            stop = ending + 1;
        } else if (flank == 1) {
            position = beginning;
            stop = ending;
        } else if (flank == 2) {
            position = ending;
            stop = beginning;
        } else {
            position = ending + 1;
            stop = beginning - 1;
        }
        if (position > length) position -= length;
        if (stop > length) stop -= length;

        int variable_count = 0;
        int total_count = 0;
        while (position != stop) {
            const short first = scan_state.sequence_data[
                position + sequences[0] * stride];
            const short second = scan_state.sequence_data[
                position + sequences[1] * stride];
            const short third = scan_state.sequence_data[
                position + sequences[2] * stride];
            if (first != 46 && second != 46 && third != 46) {
                ++total_count;
                if ((first != second || first != third) &&
                    (first == second || first == third || second == third)) {
                    ++variable_count;
                    if (variable_count == variable_site_target &&
                        total_count > total_site_target) {
                        break;
                    }
                }
            }
            position += increments[flank];
            if (position < 1) {
                position = length;
            } else if (position > length) {
                position = 1;
            }
        }
        result.positions[flank] = position;
        result.informative_counts[flank] =
            variable_count == 0 ? 1.0 : variable_count;
    }
    return result;
}


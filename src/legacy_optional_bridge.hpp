#pragma once

#include "legacy_optional/alignment.hpp"
#include "scan_state.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace next_rdp_legacy_optional_bridge {

// Adapt the source-faithful scan state (ASCII DNA5 symbols, column-major
// sentinel layout) to the compact alignment view consumed by the optional
// BootScan/SISCAN kernels.  This is a representation adapter only; it does
// not alter coordinates or missing-state semantics.
next_rdp_legacy_optional::Alignment make_alignment(
    const RdpScanState& scan_state);

next_rdp_legacy_optional::Alignment make_alignment(
    const std::vector<std::string>& sequences);

std::array<double, 3> pair_similarity(
    const next_rdp_legacy_optional::Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet);

std::vector<std::uint8_t> triplet_missing_data(
    const next_rdp_legacy_optional::Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet);

}  // namespace next_rdp_legacy_optional_bridge


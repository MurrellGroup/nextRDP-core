#include "legacy_optional_bridge.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace next_rdp_legacy_optional_bridge {
namespace {

std::uint8_t state_for_code(short code) {
  switch (code) {
    case 66: return 1;  // A
    case 68: return 2;  // C
    case 72: return 3;  // G
    case 85: return 4;  // T/U
    default: return 0;
  }
}

std::uint8_t state_for_character(char character) {
  switch (static_cast<unsigned char>(std::toupper(
      static_cast<unsigned char>(character)))) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 3;
    case 'T':
    case 'U': return 4;
    default: return 0;
  }
}

void finish_alignment(next_rdp_legacy_optional::Alignment& alignment) {
  const std::size_t count = alignment.names.size();
  if (count == 0 || alignment.length == 0) {
    throw std::runtime_error("optional-method alignment is empty");
  }
  alignment.pair_similarity.assign(count * count, 0.0F);
  for (std::size_t first = 0; first < count; ++first) {
    alignment.pair_similarity[first * count + first] = 1.0F;
    for (std::size_t second = first + 1; second < count; ++second) {
      std::size_t valid = 0;
      std::size_t matches = 0;
      for (std::size_t position = 0; position < alignment.length; ++position) {
        const auto left = alignment.at(first, position);
        const auto right = alignment.at(second, position);
        if (left == 0 || right == 0) continue;
        ++valid;
        if (left == right) ++matches;
      }
      const float similarity = valid == 0
          ? 0.0F
          : static_cast<float>(matches) / static_cast<float>(valid);
      alignment.pair_similarity[first * count + second] = similarity;
      alignment.pair_similarity[second * count + first] = similarity;
    }
  }
}

}  // namespace

next_rdp_legacy_optional::Alignment make_alignment(
    const RdpScanState& scan_state) {
  if (scan_state.next_no < 0 || scan_state.sequence_length < 1) {
    throw std::runtime_error("optional-method scan state is empty");
  }
  const std::size_t count = static_cast<std::size_t>(scan_state.next_no + 1);
  const std::size_t length = static_cast<std::size_t>(scan_state.sequence_length);
  const std::size_t stride = length + 1;
  if (scan_state.sequence_data.size() < count * stride) {
    throw std::runtime_error("optional-method scan state dimensions differ");
  }
  next_rdp_legacy_optional::Alignment alignment;
  alignment.length = length;
  alignment.names.reserve(count);
  alignment.states.assign(count * length, 0);
  for (std::size_t sequence = 0; sequence < count; ++sequence) {
    alignment.names.push_back("Sequence " + std::to_string(sequence + 1));
    for (std::size_t position = 0; position < length; ++position) {
      alignment.states[sequence * length + position] = state_for_code(
          scan_state.sequence_data[sequence * stride + position + 1]);
    }
  }
  finish_alignment(alignment);
  return alignment;
}

next_rdp_legacy_optional::Alignment make_alignment(
    const std::vector<std::string>& sequences) {
  if (sequences.empty() || sequences.front().empty()) {
    throw std::runtime_error("optional-method alignment is empty");
  }
  const std::size_t length = sequences.front().size();
  next_rdp_legacy_optional::Alignment alignment;
  alignment.length = length;
  alignment.names.reserve(sequences.size());
  alignment.states.assign(sequences.size() * length, 0);
  for (std::size_t sequence = 0; sequence < sequences.size(); ++sequence) {
    if (sequences[sequence].size() != length) {
      throw std::runtime_error("optional-method alignment lengths differ");
    }
    alignment.names.push_back("Sequence " + std::to_string(sequence + 1));
    for (std::size_t position = 0; position < length; ++position) {
      alignment.states[sequence * length + position] = state_for_character(
          sequences[sequence][position]);
    }
  }
  finish_alignment(alignment);
  return alignment;
}

std::array<double, 3> pair_similarity(
    const next_rdp_legacy_optional::Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet) {
  if (triplet[0] >= alignment.sequence_count() ||
      triplet[1] >= alignment.sequence_count() ||
      triplet[2] >= alignment.sequence_count()) {
    return {};
  }
  const auto count = alignment.sequence_count();
  return {
      alignment.pair_similarity[triplet[0] * count + triplet[1]],
      alignment.pair_similarity[triplet[0] * count + triplet[2]],
      alignment.pair_similarity[triplet[1] * count + triplet[2]],
  };
}

std::vector<std::uint8_t> triplet_missing_data(
    const next_rdp_legacy_optional::Alignment& alignment,
    const std::array<std::uint32_t, 3>& triplet) {
  std::vector<std::uint8_t> missing(alignment.length, 0);
  // Literal alignment gaps are already skipped by both optional kernels. The
  // missing-data vector is reserved for cyclic erasures; initial scans have
  // no erased working tract, so keep it all-zero here.
  (void)triplet;
  return missing;
}

}  // namespace next_rdp_legacy_optional_bridge


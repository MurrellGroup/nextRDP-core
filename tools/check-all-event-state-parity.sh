#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$workspace_dir/software/cmake/usr/bin:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null

printf '| Dataset | Redo triplets | Stored events | Row counts | Event identities | MakeTestPVs | First selection |\n'
printf '|---|---:|---:|---:|---:|---:|---:|\n'
for number in {0..9}; do
  dataset="Dataset$number"
  if [[ $number == 0 ]]; then
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
  else
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/$dataset"
  fi
  output=$(node "$project_dir/build/wasm/rdp-core.js" \
    fasta-all-redo-events-fixture \
    "$fixture_dir/$dataset.fas" \
    "$fixture_dir/alist-rdp4-v1.bin" \
    "$fixture_dir/define-event-p2-v1.bin" \
    "$fixture_dir/make-test-pvs-v1.bin" \
    "$fixture_dir/find-best-rec-signal-p2-v1.bin")
  if [[ $output =~ scan:\ ([0-9]+)\ triplets.*\ ([0-9]+)\ stored\ candidates.*\ ([0-9]+/[0-9]+)\ row\ counts\ equal,\ ([0-9]+/[0-9]+)\ event\ identities\ exact.*MakeTestPVs\ (PASS|FAIL),\ first\ selection\ (PASS|FAIL) ]]; then
    printf '| %s | %s | %s | %s | %s | %s | %s |\n' \
      "$dataset" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" \
      "${BASH_REMATCH[3]}" "${BASH_REMATCH[4]}" \
      "${BASH_REMATCH[5]}" "${BASH_REMATCH[6]}"
  else
    printf '%s\n' "$output" >&2
    exit 1
  fi
done

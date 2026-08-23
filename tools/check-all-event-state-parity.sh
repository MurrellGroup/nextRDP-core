#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$workspace_dir/software/cmake/usr/bin:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null
if [[ ! -f "$project_dir/sandbox/calcr-capture/Dataset9-calcr3-v1.bin" ]]; then
  "$project_dir/tools/capture-calcr-chain.sh"
fi

if [[ ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-make-rlist-v1.bin" ]]; then
  "$project_dir/tools/capture-make-rlist-chain.sh"
fi
if [[ ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-find-actual-events-v1.bin" ]]; then
  "$project_dir/tools/capture-find-actual-events.sh"
fi

printf '| Dataset | Redo triplets | Stored events | Row counts | Event identities | MakeTestPVs | First selection | UFDist | Region distance | CheckMatrix | First NJ tree | MakeSDMP2 | FillRmat | CalCR | MakeRList | FindActualEvents | StripDupInv | MakeRCompat |\n'
printf '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n'
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
    "$fixture_dir/find-best-rec-signal-p2-v1.bin" \
    "$fixture_dir/ufdist-v1.bin" \
    "$fixture_dir/super-dist-p2-v1.bin" \
    "$fixture_dir/check-matrix-p-v1.bin" \
    "$fixture_dir/make-nj-trees-p2-v1.bin" \
    "$fixture_dir/make-sdmp2-v1.bin" \
    "$fixture_dir/fill-rmat-y0-v1.bin" \
    "$fixture_dir/fill-rmat-y1-v1.bin" \
    "$fixture_dir/fill-rmat-y2-v1.bin" \
    "$project_dir/sandbox/calcr-capture/$dataset-calcr3-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-make-rlist-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-find-actual-events-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-strip-dup-inv-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-rcompat-v1.bin")
  if [[ $output =~ scan:\ ([0-9]+)\ triplets.*\ ([0-9]+)\ stored\ candidates.*\ ([0-9]+/[0-9]+)\ row\ counts\ equal,\ ([0-9]+/[0-9]+)\ event\ identities\ exact.*MakeTestPVs\ (PASS|FAIL),\ first\ selection\ (PASS|FAIL).*UFDist\ (PASS|FAIL),\ region\ distance\ (PASS|FAIL),\ CheckMatrix\ (PASS|FAIL),\ first\ NJ\ tree\ (PASS|FAIL),\ MakeSDMP2\ (PASS|FAIL),\ FillRmat\ (PASS|FAIL),\ CalCR\ (PASS|FAIL),\ MakeRList\ (PASS|FAIL),\ FindActualEvents\ (PASS|FAIL),\ StripDupInv\ (PASS|FAIL),\ MakeRCompat\ (PASS|FAIL) ]]; then
    printf '| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' \
      "$dataset" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" \
      "${BASH_REMATCH[3]}" "${BASH_REMATCH[4]}" \
      "${BASH_REMATCH[5]}" "${BASH_REMATCH[6]}" \
      "${BASH_REMATCH[7]}" "${BASH_REMATCH[8]}" \
      "${BASH_REMATCH[9]}" "${BASH_REMATCH[10]}" \
      "${BASH_REMATCH[11]}" "${BASH_REMATCH[12]}" \
      "${BASH_REMATCH[13]}" "${BASH_REMATCH[14]}" \
      "${BASH_REMATCH[15]}" "${BASH_REMATCH[16]}" \
      "${BASH_REMATCH[17]}"
  else
    printf '%s\n' "$output" >&2
    exit 1
  fi
done

#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/mutation-capture"
run_dir="$capture_dir/run"
trace_path="$workspace_dir/sandbox/native-trace/modseq-d6.bin"
tail_trace_path="$workspace_dir/sandbox/native-trace/mutation-tail.bin"
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
proxy_dir="$workspace_dir/sandbox/native-trace/dna-proxy"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"
node_bin="$workspace_dir/software/emsdk/node/24.19.0_64bit/bin/node"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$proxy_dir/DNA.mutation-all.dll" "$run_dir/DNA.dll"
cp "$workspace_dir/software/rdp5/prefix/drive_c/windows/syswow64/dna.dll" \
  "$run_dir/dna_original.dll"
: > "$trace_path"
: > "$tail_trace_path"

dataset_numbers=${DATASET_NUMBERS:-"0 1 2 3 4 5 6 7 8 9"}
for number in $dataset_numbers; do
  dataset="Dataset$number"
  before=$(stat -c %s "$trace_path")
  tail_before=$(stat -c %s "$tail_trace_path")
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.modseq.out" 2> "$dataset.modseq.err"
  )
  after=$(stat -c %s "$trace_path")
  tail_after=$(stat -c %s "$tail_trace_path")
  if [[ $after -le $before ]]; then
    printf 'no ModSeqNumY trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  if [[ $tail_after -le $tail_before ]]; then
    printf 'no mutation-tail trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  dd if="$trace_path" of="$capture_dir/$dataset-modseqnumy-trace.bin" \
    bs=1 skip="$before" count="$((after-before))" status=none
  "$node_bin" "$project_dir/tools/convert-modseqnumy-trace.mjs" \
    "$capture_dir/$dataset-modseqnumy-trace.bin" \
    "$capture_dir/$dataset-modseqnumy-v1.bin"
  dd if="$tail_trace_path" \
    of="$capture_dir/$dataset-mutation-tail-trace.bin" \
    bs=1 skip="$tail_before" count="$((tail_after-tail_before))" status=none
  "$node_bin" "$project_dir/tools/convert-mutation-tail-trace.mjs" \
    "$capture_dir/$dataset-mutation-tail-trace.bin" \
    "$capture_dir/$dataset-mutation-tail-v1.bin"
done

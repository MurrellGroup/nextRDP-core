#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/make-rlist-capture"
run_dir="$capture_dir/run"
trace_path="$workspace_dir/sandbox/native-trace/make-rlist-trace.bin"
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
proxy_dir="$workspace_dir/sandbox/native-trace/addjust-d7-run"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"
node_bin="$workspace_dir/software/emsdk/node/24.19.0_64bit/bin/node"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$proxy_dir/DNA.dll" "$run_dir/DNA.dll"
cp "$proxy_dir/DNA_original.dll" "$run_dir/DNA_original.dll"

for number in {0..9}; do
  dataset="Dataset$number"
  before=$(stat -c %s "$trace_path")
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.trace.out" 2> "$dataset.trace.err"
  )
  after=$(stat -c %s "$trace_path")
  if [[ $after -le $before ]]; then
    printf 'no MakeRList trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  dd if="$trace_path" of="$capture_dir/$dataset-make-rlist-trace.bin" \
    bs=1 skip="$before" count="$((after-before))" status=none
  "$node_bin" "$project_dir/tools/convert-make-rlist-trace.mjs" \
    "$capture_dir/$dataset-make-rlist-trace.bin" \
    "$capture_dir/$dataset-make-rlist-v1.bin"
  printf '%s captured\n' "$dataset"
done

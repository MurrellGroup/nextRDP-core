#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/make-rlist-capture"
run_dir="$capture_dir/run"
trace_path="$workspace_dir/sandbox/native-trace/find-actual-events-trace.bin"
strip_trace_path="$workspace_dir/sandbox/native-trace/strip-dup-inv-trace.bin"
rcompat_trace_path="$workspace_dir/sandbox/native-trace/rcompat-trace-current-d0.bin"
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
proxy_dir="$workspace_dir/sandbox/native-trace/dna-proxy"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"
node_bin="$workspace_dir/software/emsdk/node/24.19.0_64bit/bin/node"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$proxy_dir/DNA-findactual.dll" "$run_dir/DNA.dll"
cp "$workspace_dir/software/rdp5/prefix/drive_c/windows/syswow64/dna.dll" \
  "$run_dir/DNA_original.dll"
touch "$trace_path"
touch "$strip_trace_path"
touch "$rcompat_trace_path"

for number in {0..9}; do
  dataset="Dataset$number"
  before=$(stat -c %s "$trace_path")
  strip_before=$(stat -c %s "$strip_trace_path")
  rcompat_before=$(stat -c %s "$rcompat_trace_path")
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.find-actual.out" 2> "$dataset.find-actual.err"
  )
  after=$(stat -c %s "$trace_path")
  strip_after=$(stat -c %s "$strip_trace_path")
  rcompat_after=$(stat -c %s "$rcompat_trace_path")
  if [[ $after -le $before ]]; then
    printf 'no FindActualEvents trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  if [[ $strip_after -le $strip_before ]]; then
    printf 'no StripDupInv trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  if [[ $rcompat_after -le $rcompat_before ]]; then
    printf 'no MakeRCompat trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  dd if="$trace_path" \
    of="$capture_dir/$dataset-find-actual-events-trace.bin" \
    bs=1 skip="$before" count="$((after-before))" status=none
  "$node_bin" "$project_dir/tools/convert-find-actual-events-trace.mjs" \
    "$capture_dir/$dataset-find-actual-events-trace.bin" \
    "$capture_dir/$dataset-find-actual-events-v1.bin"
  dd if="$strip_trace_path" \
    of="$capture_dir/$dataset-strip-dup-inv-trace.bin" \
    bs=1 skip="$strip_before" count="$((strip_after-strip_before))" \
    status=none
  "$node_bin" "$project_dir/tools/convert-strip-dup-inv-trace.mjs" \
    "$capture_dir/$dataset-strip-dup-inv-trace.bin" \
    "$capture_dir/$dataset-strip-dup-inv-v1.bin"
  dd if="$rcompat_trace_path" \
    of="$capture_dir/$dataset-rcompat-trace.bin" \
    bs=1 skip="$rcompat_before" count="$((rcompat_after-rcompat_before))" \
    status=none
  "$node_bin" "$project_dir/tools/convert-rcompat-trace.mjs" \
    "$capture_dir/$dataset-rcompat-trace.bin" \
    "$capture_dir/$dataset-rcompat-v1.bin"
  printf '%s captured\n' "$dataset"
done

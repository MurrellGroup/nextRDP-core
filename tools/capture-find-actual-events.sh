#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/make-rlist-capture"
run_dir="$capture_dir/run"
trace_path="$workspace_dir/sandbox/native-trace/find-actual-events-trace.bin"
strip_trace_path="$workspace_dir/sandbox/native-trace/strip-dup-inv-trace.bin"
rcompat_trace_path="$workspace_dir/sandbox/native-trace/rcompat-trace-current-d0.bin"
collect_trace_path="$workspace_dir/sandbox/native-trace/collectevents-boundary.bin"
phpr_trace_path="$workspace_dir/sandbox/native-trace/make-phpr-trace.bin"
done_trace_path="$workspace_dir/sandbox/native-trace/make-done-this3-trace.bin"
trp_group_trace_path="$workspace_dir/sandbox/native-trace/make-trp-groups2-trace.bin"
trp_score_trace_path="$workspace_dir/sandbox/native-trace/make-trp-score2-trace.bin"
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
cp "$runtime_dir/DNA5.dll" "$run_dir/DNA5_original.dll"
cp "$workspace_dir/sandbox/native-trace/dna5-proxy/DNA5.collect-boundary.dll" \
  "$run_dir/DNA5.dll"
touch "$trace_path"
touch "$strip_trace_path"
touch "$rcompat_trace_path"
touch "$collect_trace_path"
touch "$phpr_trace_path"
touch "$done_trace_path" "$trp_group_trace_path" "$trp_score_trace_path"

dataset_numbers=${DATASET_NUMBERS:-"0 1 2 3 4 5 6 7 8 9"}
for number in $dataset_numbers; do
  dataset="Dataset$number"
  before=$(stat -c %s "$trace_path")
  strip_before=$(stat -c %s "$strip_trace_path")
  rcompat_before=$(stat -c %s "$rcompat_trace_path")
  collect_before=$(stat -c %s "$collect_trace_path")
  phpr_before=$(stat -c %s "$phpr_trace_path")
  done_before=$(stat -c %s "$done_trace_path")
  trp_group_before=$(stat -c %s "$trp_group_trace_path")
  trp_score_before=$(stat -c %s "$trp_score_trace_path")
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.find-actual.out" 2> "$dataset.find-actual.err"
  )
  after=$(stat -c %s "$trace_path")
  strip_after=$(stat -c %s "$strip_trace_path")
  rcompat_after=$(stat -c %s "$rcompat_trace_path")
  collect_after=$(stat -c %s "$collect_trace_path")
  phpr_after=$(stat -c %s "$phpr_trace_path")
  done_after=$(stat -c %s "$done_trace_path")
  trp_group_after=$(stat -c %s "$trp_group_trace_path")
  trp_score_after=$(stat -c %s "$trp_score_trace_path")
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
  if [[ $collect_after -le $collect_before ]]; then
    printf 'no MakeCollecteventsC boundary captured for %s\n' "$dataset" >&2
    exit 1
  fi
  if [[ $phpr_after -le $phpr_before ]]; then
    printf 'no MakePhPrScore trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  if [[ $done_after -le $done_before || $trp_group_after -le $trp_group_before ||
        $trp_score_after -le $trp_score_before ]]; then
    printf 'incomplete score support traces for %s\n' "$dataset" >&2
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
  dd if="$collect_trace_path" \
    of="$capture_dir/$dataset-collectevents-boundary.bin" \
    bs=1 skip="$collect_before" count="$((collect_after-collect_before))" \
    status=none
  "$node_bin" "$project_dir/tools/convert-collectevents-boundary.mjs" \
    "$capture_dir/$dataset-collectevents-boundary.bin" \
    "$capture_dir/$dataset-collectevents-v1.bin"
  dd if="$phpr_trace_path" \
    of="$capture_dir/$dataset-phpr-trace.bin" \
    bs=1 skip="$phpr_before" count="$((phpr_after-phpr_before))" \
    status=none
  "$node_bin" "$project_dir/tools/convert-phpr-trace.mjs" \
    "$capture_dir/$dataset-phpr-trace.bin" \
    "$capture_dir/$dataset-phpr-v1.bin"
  dd if="$done_trace_path" of="$capture_dir/$dataset-done-this3-trace.bin" \
    bs=1 skip="$done_before" count="$((done_after-done_before))" status=none
  dd if="$trp_group_trace_path" of="$capture_dir/$dataset-trp-groups2-trace.bin" \
    bs=1 skip="$trp_group_before" count="$((trp_group_after-trp_group_before))" status=none
  dd if="$trp_score_trace_path" of="$capture_dir/$dataset-trp-score2-trace.bin" \
    bs=1 skip="$trp_score_before" count="$((trp_score_after-trp_score_before))" status=none
  "$node_bin" "$project_dir/tools/convert-score-support-traces.mjs" \
    "$capture_dir/$dataset-done-this3-trace.bin" \
    "$capture_dir/$dataset-trp-groups2-trace.bin" \
    "$capture_dir/$dataset-trp-score2-trace.bin" \
    "$capture_dir/$dataset-score-support-v1.bin"
  printf '%s captured\n' "$dataset"
done

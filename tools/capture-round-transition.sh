#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/round-transition-capture"
run_dir="$capture_dir/run"
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
trace_dir="$workspace_dir/sandbox/native-trace"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$workspace_dir/software/rdp5/prefix/drive_c/windows/syswow64/dna.dll" \
  "$run_dir/DNA.dll"
cp "$runtime_dir/DNA5.dll" "$run_dir/DNA5_original.dll"
cp "$trace_dir/dna5-proxy/DNA5.round-trace.dll" "$run_dir/DNA5.dll"

trace_names=(
  addjust-input-records.bin
  addjust-cxo.bin
  addjust-cxo-temp.bin
  addjust-rlist.bin
  addjust-dopairs.bin
  alist-rdp3-calls.bin
  findbetter-pxolist.bin
  update-done-pvco.bin
)
for trace_name in "${trace_names[@]}"; do
  touch "$trace_dir/$trace_name"
done

dataset_numbers=${DATASET_NUMBERS:-"0 1 2 3 4 5 6 7 8 9"}
for number in $dataset_numbers; do
  dataset="Dataset$number"
  declare -A before
  for trace_name in "${trace_names[@]}"; do
    before[$trace_name]=$(stat -c %s "$trace_dir/$trace_name")
  done
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.round-transition.out" \
      2> "$dataset.round-transition.err"
  )
  for trace_name in "${trace_names[@]}"; do
    after=$(stat -c %s "$trace_dir/$trace_name")
    if [[ $after -le ${before[$trace_name]} ]]; then
      printf 'no %s trace captured for %s\n' "$trace_name" "$dataset" >&2
      exit 1
    fi
    dd if="$trace_dir/$trace_name" \
      of="$capture_dir/$dataset-$trace_name" bs=1 \
      skip="${before[$trace_name]}" \
      count="$((after-before[$trace_name]))" status=none
  done
  printf '%s round transitions captured\n' "$dataset"
done

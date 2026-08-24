#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/findbest-series-capture"
run_dir="$capture_dir/run"
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
trace_dir="$workspace_dir/sandbox/native-trace"
proxy_dir="$trace_dir/dna5-proxy"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"
trace_path="$trace_dir/findbest-series.bin"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$workspace_dir/software/rdp5/prefix/drive_c/windows/syswow64/dna.dll" \
  "$run_dir/DNA.dll"
cp "$runtime_dir/DNA5.dll" "$run_dir/DNA5_original.dll"
cp "$proxy_dir/DNA5.findbest-series.dll" "$run_dir/DNA5.dll"
touch "$trace_path"

dataset_numbers=${DATASET_NUMBERS:-"0 1 2 3 4 5 6 7 8 9"}
for number in $dataset_numbers; do
  dataset="Dataset$number"
  before=$(stat -c %s "$trace_path")
  (
    cd "$run_dir"
    PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" \
      timeout 45s xvfb-run -a wine RDP5CL.exe "-f$dataset.fas" -ds \
      > "$dataset.findbest-series.out" \
      2> "$dataset.findbest-series.err"
  )
  after=$(stat -c %s "$trace_path")
  if [[ $after -le $before ]]; then
    printf 'no FindBestRecSignalP2 trace captured for %s\n' "$dataset" >&2
    exit 1
  fi
  dd if="$trace_path" \
    of="$capture_dir/$dataset-findbest-series.bin" bs=1 \
    skip="$before" count="$((after-before))" status=none
  printf '%s selection series captured\n' "$dataset"
done

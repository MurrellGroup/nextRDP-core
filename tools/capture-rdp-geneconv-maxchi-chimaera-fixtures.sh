#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
runtime=${1:-$project_dir/sandbox/method-combos/rdp_gc_mc_ch_all_runtime}
fixture_root=${2:-$project_dir/sandbox/method-combos/rdp_gc_mc_ch_all}
first_dataset=${3:-0}
last_dataset=${4:-9}

export PATH="$workspace_dir/software/wine-11.13/bin:$PATH"
export WINEPREFIX="$workspace_dir/sandbox/rdp5-trace-prefix"
export WINEDEBUG=-all

mkdir -p "$fixture_root"
for number in $(seq "$first_dataset" "$last_dataset"); do
  name="Dataset$number"
  (
    cd "$runtime"
    timeout 120s xvfb-run -a wine RDP5CL.exe "-f$name.fas" -ds \
      > "$name.chimaera-capture.out" \
      2> "$name.chimaera-capture.err"
  )
  fixture="$fixture_root/$name"
  mkdir -p "$fixture"
  for file in \
      "$name.fas" \
      alist-rdp4-v1.bin \
      define-event-p2-v1.bin \
      make-test-pvs-v1.bin \
      geneconv-call-order.bin \
      geneconv-count-at-first-make-test.bin \
      maxchi-call-order.bin \
      maxchi-count-at-first-make-test.bin \
      chimaera-call-order.bin \
      chimaera-count-at-first-make-test.bin; do
    test -s "$runtime/$file"
    cp "$runtime/$file" "$fixture/$file"
  done
  printf '%s captured\n' "$name"
done

#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 4 )); then
  echo "usage: $0 METHOD_MASK DATASET [RUNTIME] [FIXTURE_ROOT]" >&2
  exit 2
fi

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
mask=$1
name=$2
runtime=${3:-$project_dir/sandbox/method-combos/rdp_gc_mc_ch_all_runtime}
fixture_root=${4:-$project_dir/sandbox/method-combos/random_samples}

geneconv=$((mask & 1 ? 1 : 0))
maxchi=$((mask & 2 ? 1 : 0))
chimaera=$((mask & 4 ? 1 : 0))
perl -pi -e \
  "if (\$. == 13) { \$_ = \"1,$geneconv,0,$maxchi,$chimaera,0\\r\\n\" }" \
  "$runtime/RDP.ini"

export PATH="$workspace_dir/software/wine-11.13/bin:$PATH"
export WINEPREFIX="$workspace_dir/sandbox/rdp5-trace-prefix"
export WINEDEBUG=-all
(
  cd "$runtime"
  timeout 120s xvfb-run -a wine RDP5CL.exe "-f$name.fas" -ds \
    > "$name-mask$mask-capture.out" \
    2> "$name-mask$mask-capture.err"
)

fixture="$fixture_root/mask$mask-$name"
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
printf '%s mask %d captured in %s\n' "$name" "$mask" "$fixture"

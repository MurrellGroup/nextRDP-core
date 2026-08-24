#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
fixture_root=${1:-$project_dir/sandbox/method-combos/rdp_gc_mc_all}

export PATH="$workspace_dir/software/cmake/usr/bin:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
export EMSDK_QUIET=1
source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null

cmake -S "$project_dir" -B "$project_dir/build/wasm" >/dev/null
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null

failures=0
printf '%-9s %s\n' Dataset Result
for number in $(seq 0 9); do
  name="Dataset$number"
  fixture="$fixture_root/$name"
  set +e
  output=$(
    "$workspace_dir/software/emsdk/node/24.19.0_64bit/bin/node" \
      "$project_dir/build/wasm/rdp-core.js" \
      fasta-geneconv-maxchi-events-fixture \
      "$fixture/$name.fas" \
      "$fixture/alist-rdp4-v1.bin" \
      "$fixture/define-event-p2-v1.bin" \
      "$fixture/make-test-pvs-v1.bin" \
      "$fixture/maxchi-call-order.bin" \
      "$fixture/maxchi-count-at-first-make-test.bin" 2>&1)
  status=$?
  set -e
  line=$(printf '%s\n' "$output" | tail -n 1)
  printf '%-9s %s\n' "$name" "$line"
  if ((status != 0)); then
    failures=$((failures + 1))
    printf '%s\n' "$output" > "$fixture/wasm-parity.log"
  fi
done

if ((failures != 0)); then
  printf '%d dataset(s) differ; detailed logs are in the fixture directories.\n' \
    "$failures" >&2
  exit 1
fi

#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
fixture_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
cmake_dir="$workspace_dir/software/cmake/usr/bin"

"$project_dir/tools/check-preprocess-parity.sh"
"$project_dir/tools/capture-alist-rdp4.sh"

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$cmake_dir:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null

node "$project_dir/build/wasm/rdp-core.js" \
  fasta-preprocess-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  fasta-distance-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  fasta-tree-distance-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  fasta-alist-rdp4-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  fasta-first-xover-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin" \
  "$fixture_dir/find-subseq-pb3-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  fasta-first-xover-walk-fixture "$fixture_dir/Dataset0.fas" \
  "$fixture_dir/alist-rdp4-v1.bin" \
  "$fixture_dir/find-subseq-pb3-v1.bin" \
  "$fixture_dir/xohomology-p-v1.bin" \
  "$fixture_dir/find-next-p-v1.bin" \
  "$fixture_dir/define-event-p2-v1.bin" \
  "$fixture_dir/prob-calc-p2-v1.bin" \
  "$fixture_dir/prob-calc-p-v1.bin" \
  "$fixture_dir/find-subseq-pb4-v1.bin"
node "$project_dir/build/wasm/rdp-core.js" \
  oracle-fixture-chain "$fixture_dir"

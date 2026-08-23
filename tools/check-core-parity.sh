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
  oracle-fixture-chain "$fixture_dir"

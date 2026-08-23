#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$workspace_dir/software/cmake/usr/bin:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null

printf '| Dataset | FASTA state | Distance | Tree distance | AlistRDP4 |\n'
printf '|---|---:|---:|---:|---:|\n'
for number in {0..9}; do
  dataset="Dataset$number"
  "$project_dir/tools/capture-alist-rdp4.sh" "$dataset" >/dev/null
  if [[ $dataset == Dataset0 ]]; then
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
  else
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/$dataset"
  fi
  fasta="$fixture_dir/$dataset.fas"
  fixture="$fixture_dir/alist-rdp4-v1.bin"
  node "$project_dir/build/wasm/rdp-core.js" \
    fasta-preprocess-fixture "$fasta" "$fixture" >/dev/null
  node "$project_dir/build/wasm/rdp-core.js" \
    fasta-distance-fixture "$fasta" "$fixture" >/dev/null
  node "$project_dir/build/wasm/rdp-core.js" \
    fasta-tree-distance-fixture "$fasta" "$fixture" >/dev/null
  node "$project_dir/build/wasm/rdp-core.js" \
    fasta-alist-rdp4-fixture "$fasta" "$fixture" >/dev/null
  printf '| %s | PASS | PASS | PASS | PASS |\n' "$dataset"
done

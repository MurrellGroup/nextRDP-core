#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
oracle_dll="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt/DNA5.dll"
mingw_bin="$workspace_dir/software/mingw/root/usr/bin"
wine_dir="$workspace_dir/software/wine-11.13"
cmake_dir="$workspace_dir/software/cmake/usr/bin"
run_dir="$project_dir/sandbox/preprocess-parity"

mkdir -p "$run_dir"

"$mingw_bin/i686-w64-mingw32-g++-posix" \
    -std=c++20 -O2 -static \
    -I"$project_dir/tests" \
    "$project_dir/tests/oracle_preprocess.cpp" \
    -o "$run_dir/oracle-preprocess.exe"
"$mingw_bin/i686-w64-mingw32-g++-posix" \
    -std=c++20 -O2 -static \
    -I"$project_dir/tests" \
    "$project_dir/tests/oracle_distance.cpp" \
    -o "$run_dir/oracle-distance.exe"

export PATH="$wine_dir/bin:$cmake_dir:$PATH"
export LD_LIBRARY_PATH="$wine_dir/lib:$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
export WINEDLLPATH="$wine_dir/lib/wine"
export WINEPREFIX="$workspace_dir/sandbox/rdp5-trace-prefix"

oracle_windows_path=$(winepath -w "$oracle_dll")
wine "$run_dir/oracle-preprocess.exe" "$oracle_windows_path" > "$run_dir/oracle.json" 2> "$run_dir/oracle.stderr"
tr -d '\r' < "$run_dir/oracle.json" > "$run_dir/oracle.normalized.json"
wine "$run_dir/oracle-distance.exe" "$oracle_windows_path" > "$run_dir/oracle-distance.json" 2> "$run_dir/oracle-distance.stderr"
tr -d '\r' < "$run_dir/oracle-distance.json" > "$run_dir/oracle-distance.normalized.json"

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$cmake_dir:$PATH"
emcmake cmake -S "$project_dir" -B "$project_dir/build/wasm" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null
node "$project_dir/build/wasm/rdp-core.js" preprocess-fixture > "$run_dir/wasm.json"
node "$project_dir/build/wasm/rdp-core.js" distance-fixture > "$run_dir/wasm-distance.json"

cmp "$run_dir/oracle.normalized.json" "$run_dir/wasm.json"
cmp "$run_dir/oracle-distance.normalized.json" "$run_dir/wasm-distance.json"
echo "distance + preprocess parity: PASS"

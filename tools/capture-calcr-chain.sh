#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
capture_dir="$project_dir/sandbox/calcr-capture"
mingw_bin="$workspace_dir/software/mingw/root/usr/bin"
wine_bin="$workspace_dir/software/wine-11.13/bin"
wine_prefix="$workspace_dir/software/rdp5/prefix"

mkdir -p "$capture_dir"
"$mingw_bin/i686-w64-mingw32-gcc-posix" -O2 -static \
  -o "$capture_dir/capture-calcr-chain.exe" \
  "$project_dir/tools/capture-calcr-chain.c"

to_windows_path() {
  local path=$1
  printf 'Z:%s' "${path//\//\\}"
}

for number in {0..9}; do
  dataset="Dataset$number"
  fixture_dir="$project_dir/sandbox/alist-rdp4-capture/$dataset"
  if [[ $number == 0 ]]; then
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
  fi
  PATH="$wine_bin:$PATH" WINEPREFIX="$wine_prefix" wine \
    "$capture_dir/capture-calcr-chain.exe" \
    "$(to_windows_path "$fixture_dir/make-sdmp2-v1.bin")" \
    "$(to_windows_path "$fixture_dir/fill-rmat-y0-v1.bin")" \
    "$(to_windows_path "$fixture_dir/fill-rmat-y1-v1.bin")" \
    "$(to_windows_path "$fixture_dir/fill-rmat-y2-v1.bin")" \
    "$(to_windows_path "$capture_dir/$dataset-calcr3-v1.bin")"
done


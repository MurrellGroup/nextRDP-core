#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
runtime_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
capture_dir="$project_dir/sandbox/calcmatch-capture"
run_dir="$capture_dir/runtime"
source_def="$workspace_dir/sandbox/source-build/dna-old/threshold.DEF"
mingw_bin="$workspace_dir/software/mingw/root/usr/bin"
wine_dir="$workspace_dir/software/wine-11.13"
wine_prefix="$workspace_dir/software/rdp5/prefix"

mkdir -p "$run_dir"
cp -a "$runtime_dir/." "$run_dir/"
cp "$wine_prefix/drive_c/windows/syswow64/dna.dll" \
  "$run_dir/dna_original.dll"

awk '
  BEGIN { print "LIBRARY dna"; print "EXPORTS" }
  /^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*[[:space:]]+@[0-9]+/ {
    name=$1
    ordinal=$2
    if (name == "MakeVarMap2")
      print "    MakeVarMap2=MakeVarMap2Capture@36 " ordinal
    else if (name == "MakeCntHit2")
      print "    MakeCntHit2=MakeCntHit2Capture@44 " ordinal
    else
      print "    " name "=dna_original." name " " ordinal
  }
' "$source_def" > "$capture_dir/dna.capture.def"

"$mingw_bin/i686-w64-mingw32-gcc-posix" \
  -std=c11 -O2 -shared -static \
  "$project_dir/tools/dna_calcmatch_capture.c" \
  "$capture_dir/dna.capture.def" \
  -Wl,--out-implib,"$capture_dir/libdna-capture.a" \
  -o "$run_dir/dna.dll"

rm -f "$run_dir/make-var-map2-trace.bin" \
  "$run_dir/make-cnt-hit2-trace.bin"

export PATH="$wine_dir/bin:$PATH"
export LD_LIBRARY_PATH="$wine_dir/lib:${LD_LIBRARY_PATH:-}"
export WINEDLLPATH="$wine_dir/lib/wine"
export WINEPREFIX="$wine_prefix"

(
  cd "$run_dir"
  timeout 180s xvfb-run -a wine RDP5CL.exe -fDataset0.fas -ds \
    > capture.stdout 2> capture.stderr
)

test -s "$run_dir/make-var-map2-trace.bin"
test -s "$run_dir/make-cnt-hit2-trace.bin"
printf 'captured %s MakeVarMap2 calls and %s MakeCntHit2 calls\n' \
  "$(od -An -tu4 -j4 -N4 "$run_dir/make-var-map2-trace.bin" | tr -d ' ')" \
  "$(od -An -tu4 -j4 -N4 "$run_dir/make-cnt-hit2-trace.bin" | tr -d ' ')"

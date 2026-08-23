#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)
oracle_dir="$workspace_dir/sandbox/source-build/rdp-dll-smoke/rebuilt"
mingw_bin="$workspace_dir/software/mingw/root/usr/bin"
wine_dir="$workspace_dir/software/wine-11.13"
run_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
proxy_def="$project_dir/sandbox/alist-rdp4-capture/DNA5.capture.def"

mkdir -p "$run_dir"

# The pinned DLL's existing definition file has its complete export list and
# ordinals. Generate an all-forwarding proxy, replacing only AlistRDP4.
awk '
  BEGIN { print "LIBRARY DNA5"; print "EXPORTS" }
  /^[[:space:]]+[A-Za-z_]/ {
    name=$1
    sub(/=.*/, "", name)
    ordinal=$NF
    if (name == "AlistRDP4")
      print "    AlistRDP4=AlistRDP4Capture@144 " ordinal
    else if (name == "FindSubSeqPB3")
      print "    FindSubSeqPB3=FindSubSeqPB3Capture@52 " ordinal
    else if (name == "XOHomologyP")
      print "    XOHomologyP=XOHomologyPCapture@24 " ordinal
    else if (name == "FindNextP")
      print "    FindNextP=FindNextPCapture@32 " ordinal
    else if (name == "DefineEventP2")
      print "    DefineEventP2=DefineEventP2Capture@84 " ordinal
    else if (name == "ProbCalcP2")
      print "    ProbCalcP2=ProbCalcP2Capture@28 " ordinal
    else
      print "    " name "=DNA5_original." name " " ordinal
  }
' "$workspace_dir/sandbox/native-trace/dna5-proxy/DNA5.collect-input.def" > "$proxy_def"

"$mingw_bin/i686-w64-mingw32-gcc-posix" \
    -std=c11 -O2 -shared -static \
    "$project_dir/tools/alist_rdp4_capture.c" "$proxy_def" \
    -Wl,--out-implib,"$project_dir/sandbox/alist-rdp4-capture/libDNA5-capture.a" \
    -o "$project_dir/sandbox/alist-rdp4-capture/DNA5.dll"

# This runtime is disposable. Copy explicitly into this project's sandbox.
cp -a "$oracle_dir/." "$run_dir/"
cp "$run_dir/DNA5.dll" "$run_dir/DNA5_original.dll"
cp "$project_dir/sandbox/alist-rdp4-capture/DNA5.dll" "$run_dir/DNA5.dll"
rm -f "$run_dir/alist-rdp4-v1.bin"
rm -f "$run_dir/find-subseq-pb3-v1.bin"
rm -f "$run_dir/xohomology-p-v1.bin" "$run_dir/find-next-p-v1.bin"
rm -f "$run_dir/define-event-p2-v1.bin"
rm -f "$run_dir/prob-calc-p2-v1.bin"

export PATH="$wine_dir/bin:$PATH"
export LD_LIBRARY_PATH="$wine_dir/lib:${LD_LIBRARY_PATH:-}"
export WINEDLLPATH="$wine_dir/lib/wine"
export WINEPREFIX="$workspace_dir/sandbox/rdp5-trace-prefix"

(
  cd "$run_dir"
  xvfb-run -a wine RDP5CL.exe -fDataset0.fas -ds > capture.stdout 2> capture.stderr
)

test -s "$run_dir/alist-rdp4-v1.bin"
test -s "$run_dir/find-subseq-pb3-v1.bin"
test -s "$run_dir/xohomology-p-v1.bin"
test -s "$run_dir/find-next-p-v1.bin"
test -s "$run_dir/define-event-p2-v1.bin"
test -s "$run_dir/prob-calc-p2-v1.bin"
echo "captured $run_dir/alist-rdp4-v1.bin ($(stat -c %s "$run_dir/alist-rdp4-v1.bin") bytes)"
echo "captured $run_dir/find-subseq-pb3-v1.bin ($(stat -c %s "$run_dir/find-subseq-pb3-v1.bin") bytes)"
echo "captured $run_dir/xohomology-p-v1.bin ($(stat -c %s "$run_dir/xohomology-p-v1.bin") bytes)"
echo "captured $run_dir/find-next-p-v1.bin ($(stat -c %s "$run_dir/find-next-p-v1.bin") bytes)"
echo "captured $run_dir/define-event-p2-v1.bin ($(stat -c %s "$run_dir/define-event-p2-v1.bin") bytes)"
echo "captured $run_dir/prob-calc-p2-v1.bin ($(stat -c %s "$run_dir/prob-calc-p2-v1.bin") bytes)"

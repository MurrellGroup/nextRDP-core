#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "$0")/.." && pwd)
workspace_dir=$(cd "$project_dir/.." && pwd)

source "$workspace_dir/software/emsdk/emsdk_env.sh" >/dev/null
export PATH="$workspace_dir/software/cmake/usr/bin:$PATH"
export LD_LIBRARY_PATH="$workspace_dir/software/cmake/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
cmake --build "$project_dir/build/wasm" --parallel 1 >/dev/null
if [[ ! -f "$project_dir/sandbox/calcr-capture/Dataset9-calcr3-v1.bin" ]]; then
  "$project_dir/tools/capture-calcr-chain.sh"
fi

if [[ ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-make-rlist-v1.bin" ]]; then
  "$project_dir/tools/capture-make-rlist-chain.sh"
fi
if [[ ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-find-actual-events-v1.bin" ||
      ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-phpr-v1.bin" ||
      ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-score-support-v1.bin" ||
      ! -f "$project_dir/sandbox/make-rlist-capture/Dataset9-check-pattern-v1.bin" ]]; then
  "$project_dir/tools/capture-find-actual-events.sh"
fi
if [[ ! -f "$project_dir/sandbox/findbest-series-capture/Dataset9-findbest-series.bin" ]]; then
  "$project_dir/tools/capture-findbest-series.sh"
fi

printf '| Dataset | Redo triplets | Stored events | Row counts | Event identities | MakeTestPVs | First selection | UFDist | Region distance | CheckMatrix | First NJ tree | MakeSDMP2 | FillRmat | CalCR | MakeRList | FindActualEvents | StripDupInv | MakeRCompat | Compatibility flow | MakePhPrScore | Score support | CheckPatternX | FinalTrim prefix | CalcMaxD | MakeConsensusC | Second FinalTrim | MakeRelevant/Collect | ModSeqNumY | ModSN/ModSeqNumZ | AddjustCXO events | Rescan schedules | AlistRDP3 rescans | Post-rescan events | Second MakeTestPVs | Second selection | Second round prefix | Reusable complete round |\n'
printf '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n'
dataset_numbers=${DATASET_NUMBERS:-"0 1 2 3 4 5 6 7 8 9"}
overall_status=0
for number in $dataset_numbers; do
  dataset="Dataset$number"
  if [[ $number == 0 ]]; then
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/runtime"
  else
    fixture_dir="$project_dir/sandbox/alist-rdp4-capture/$dataset"
  fi
  set +e
  output=$(node "$project_dir/build/wasm/rdp-core.js" \
    fasta-all-redo-events-fixture \
    "$fixture_dir/$dataset.fas" \
    "$fixture_dir/alist-rdp4-v1.bin" \
    "$fixture_dir/define-event-p2-v1.bin" \
    "$fixture_dir/make-test-pvs-v1.bin" \
    "$fixture_dir/find-best-rec-signal-p2-v1.bin" \
    "$fixture_dir/ufdist-v1.bin" \
    "$fixture_dir/super-dist-p2-v1.bin" \
    "$fixture_dir/check-matrix-p-v1.bin" \
    "$fixture_dir/make-nj-trees-p2-v1.bin" \
    "$fixture_dir/make-sdmp2-v1.bin" \
    "$fixture_dir/fill-rmat-y0-v1.bin" \
    "$fixture_dir/fill-rmat-y1-v1.bin" \
    "$fixture_dir/fill-rmat-y2-v1.bin" \
    "$project_dir/sandbox/calcr-capture/$dataset-calcr3-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-make-rlist-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-find-actual-events-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-strip-dup-inv-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-rcompat-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-phpr-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-score-support-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-check-pattern-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-cmaxd2p3-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-consensus-v1.bin" \
    "$project_dir/sandbox/make-rlist-capture/$dataset-collectevents-v1.bin" \
    "$project_dir/sandbox/mutation-capture/$dataset-modseqnumy-v1.bin" \
    "$project_dir/sandbox/mutation-capture/$dataset-mutation-tail-v1.bin" \
    "$project_dir/sandbox/round-transition-capture/$dataset-alist-rdp3-calls.bin" \
    "$project_dir/sandbox/round-transition-capture/$dataset-addjust-dopairs.bin" \
    "$project_dir/sandbox/findbest-series-capture/$dataset-findbest-series.bin")
  command_status=$?
  set -e
  if [[ $command_status -ne 0 ]]; then overall_status=1; fi
  if [[ $output =~ scan:\ ([0-9]+)\ triplets.*\ ([0-9]+)\ stored\ candidates.*\ ([0-9]+/[0-9]+)\ row\ counts\ equal,\ ([0-9]+/[0-9]+)\ event\ identities\ exact.*MakeTestPVs\ (PASS|FAIL),\ first\ selection\ (PASS|FAIL).*UFDist\ (PASS|FAIL),\ region\ distance\ (PASS|FAIL),\ CheckMatrix\ (PASS|FAIL),\ first\ NJ\ tree\ (PASS|FAIL),\ MakeSDMP2\ (PASS|FAIL),\ FillRmat\ (PASS|FAIL),\ CalCR\ (PASS|FAIL),\ MakeRList\ (PASS|FAIL),\ FindActualEvents\ (PASS|FAIL),\ StripDupInv\ (PASS|FAIL),\ MakeRCompat\ (PASS|FAIL),\ compatibility\ flow\ (PASS|FAIL),\ MakePhPrScore\ (PASS|FAIL),\ score\ support\ (PASS|FAIL),\ CheckPatternX\ (PASS|FAIL),\ FinalTrim\ prefix\ (PASS|FAIL),\ CalcMaxD\ (PASS|FAIL),\ MakeConsensusC\ (PASS|FAIL),\ second\ FinalTrim\ (PASS|FAIL),\ MakeRelevant/MakeCollecteventsC\ (PASS|FAIL),\ ModSeqNumY\ (PASS|FAIL),\ ModSN/ModSeqNumZ\ (PASS|FAIL),\ AddjustCXO\ events\ (PASS|FAIL),\ rescan\ schedules\ (PASS|FAIL),\ AlistRDP3\ rescans\ (PASS|FAIL),\ post-rescan\ events\ (PASS|FAIL),\ second\ MakeTestPVs\ (PASS|FAIL),\ second\ selection\ (PASS|FAIL),\ second\ round\ prefix\ (PASS|FAIL),\ reusable\ complete\ round\ (PASS|FAIL) ]]; then
    printf '| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' \
      "$dataset" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" \
      "${BASH_REMATCH[3]}" "${BASH_REMATCH[4]}" \
      "${BASH_REMATCH[5]}" "${BASH_REMATCH[6]}" \
      "${BASH_REMATCH[7]}" "${BASH_REMATCH[8]}" \
      "${BASH_REMATCH[9]}" "${BASH_REMATCH[10]}" \
      "${BASH_REMATCH[11]}" "${BASH_REMATCH[12]}" \
      "${BASH_REMATCH[13]}" "${BASH_REMATCH[14]}" \
      "${BASH_REMATCH[15]}" "${BASH_REMATCH[16]}" \
      "${BASH_REMATCH[17]}" "${BASH_REMATCH[18]}" \
      "${BASH_REMATCH[19]}" "${BASH_REMATCH[20]}" \
      "${BASH_REMATCH[21]}" "${BASH_REMATCH[22]}" \
      "${BASH_REMATCH[23]}" "${BASH_REMATCH[24]}" \
      "${BASH_REMATCH[25]}" "${BASH_REMATCH[26]}" \
      "${BASH_REMATCH[27]}" "${BASH_REMATCH[28]}" \
      "${BASH_REMATCH[29]}" "${BASH_REMATCH[30]}" \
      "${BASH_REMATCH[31]}" "${BASH_REMATCH[32]}" \
      "${BASH_REMATCH[33]}" "${BASH_REMATCH[34]}" \
      "${BASH_REMATCH[35]}" "${BASH_REMATCH[36]}"
  else
    printf '%s\n' "$output" >&2
    exit 1
  fi
done
exit "$overall_status"

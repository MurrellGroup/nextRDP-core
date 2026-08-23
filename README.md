# nextRDP core parity port

This directory deliberately contains no UI. It has three layers:

1. `vendor/dna5`: the exact supplied C++ computational source.
2. `src`: a WASM-compatible command-line orchestration port of the VB execution path.
3. `tests`: per-routine parity runners that compare the local Windows oracle with the WASM build.

The initial execution order was recovered from the local RDP-only run by observing the first resolution of DNA5 exports. It is recorded in `reference/rdp-only-dna5-order.txt`. The list is a lower bound: VB routines between these calls are added to the manifest as their boundaries are instrumented.

Build the WASM command-line target:

```sh
source ../software/emsdk/emsdk_env.sh
emcmake cmake -S . -B build/wasm
cmake --build build/wasm --parallel 1
node build/wasm/rdp-core.js self-test
```

Generated oracle traces and comparison output belong in `sandbox/`, which is ignored by git.

Run every implemented oracle gate in execution order:

```sh
./tools/check-core-parity.sh
```

See `PARITY.md` for the current routine-by-routine ledger.

The live captured routine boundaries can also be replayed through one WASM
command, in their observed RDP execution order:

```sh
node build/wasm/rdp-core.js oracle-fixture-chain sandbox/alist-rdp4-capture/runtime
```

This currently checks 25 distinct routines (22 live calls because `FillRmat`
is exercised for all three `Y` branches). The pinned DLL is a 32-bit x87 build;
the WASM compilation uses Clang's extended evaluation mode to preserve that
source-level floating-point ABI without editing the supplied routine bodies.
The compatibility layer also supplies the Microsoft CRT 15-bit `rand()` stream.

The first shared-state orchestration segment now starts from FASTA and calls the
same preprocessing routines in order. Its live oracle gate is:

```sh
node build/wasm/rdp-core.js fasta-preprocess-fixture \
  sandbox/alist-rdp4-capture/runtime/Dataset0.fas \
  sandbox/alist-rdp4-capture/runtime/alist-rdp4-v1.bin
```

The next shared-state segment reconstructs RDP's compressed distance pipeline
from the same FASTA, calls `SuperDistP` row-by-row as observed, and compares the
complete matrix passed into `AlistRDP4`:

```sh
node build/wasm/rdp-core.js fasta-distance-fixture \
  sandbox/alist-rdp4-capture/runtime/Dataset0.fas \
  sandbox/alist-rdp4-capture/runtime/alist-rdp4-v1.bin
```

The following source-faithful UPGMA block is gated separately:

```sh
node build/wasm/rdp-core.js fasta-tree-distance-fixture \
  sandbox/alist-rdp4-capture/runtime/Dataset0.fas \
  sandbox/alist-rdp4-capture/runtime/alist-rdp4-v1.bin
```

The generated states are then wired into the actual `AlistRDP4` call and the
first selected `XOver` path through `FindSubSeqPB3`, `XOHomologyP`, the VB role
ranking block, `FindNextP`, circular-start handling, `DefineEventP2`, the
`ProbEstimate` prefilter, `ProbCalcP`/`ProbCalcP2`, the source role cycles, and
the first significant event's `FindSubSeqPB4` position-map reconstruction. It
then exhausts the remaining interval/role scan and reaches the exact
end-of-triplet `CleanXOSNW` state. The redo list used by this walk is produced
by the generated `AlistRDP4` call in the same process. Run the all-dataset
matrix for this connected prefix with:

```sh
./tools/check-all-pre-scan-parity.sh
```

Run the all-redo event-state, first-selection, regional-distance,
`CheckMatrixP`, first-NJ-tree, `MakeSDMP2`, and `FillRmat` parity table with:

```sh
./tools/check-all-event-state-parity.sh
```

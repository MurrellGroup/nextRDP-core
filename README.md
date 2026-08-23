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

This currently checks 25 distinct routines (21 live calls because `FillRmat`
is exercised for all three `Y` branches). The pinned DLL is a 32-bit x87 build;
the WASM compilation uses Clang's extended evaluation mode to preserve that
source-level floating-point ABI without editing the supplied routine bodies.
The compatibility layer also supplies the Microsoft CRT 15-bit `rand()` stream.

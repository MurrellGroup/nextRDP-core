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

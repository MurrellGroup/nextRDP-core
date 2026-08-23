# Routine parity ledger

Oracle means the pinned locally compiled `DNA5.dll` invoked by the pinned
`RDP5CL.exe`. All discrete state comparisons are byte-for-byte. Probability
outputs use the user's accepted platform-math tolerance when Windows and WASM
`pow` differ only in final rounding.

| Execution order | Routine | Oracle boundary | Current result |
|---:|---|---|---|
| 1 | `SuperDistP` | synthetic input, Windows DLL vs WASM | byte-identical |
| 2 | `MakeAListP2` | synthetic input, Windows DLL vs WASM | byte-identical |
| 3 | `CountNucs` | chained synthetic state | byte-identical |
| 4 | `RecodeNucs` | chained synthetic state | byte-identical |
| 5 | `DoRecodeP` | chained synthetic state | byte-identical |
| 6 | `MakeCompressSeqP` | chained synthetic state | byte-identical |
| 7 | `AlistRDP4` | live Dataset0 call from `RDP5CL.exe` | byte-identical result, `RL`, and `StoreLPV` |
| 8 | `FindSubSeqPB3` | live Dataset0 call from `RDP5CL.exe` | byte-identical result, `AH`, and window state |
| 9 | `XOHomologyP` | live Dataset0 call from `RDP5CL.exe` | byte-identical result and homology state |
| 10 | `FindNextP` | live Dataset0 call from `RDP5CL.exe` | byte-identical result |
| 11 | `DefineEventP2` | live Dataset0 call from `RDP5CL.exe` | byte-identical result and all scalar output state |
| 12 | `ProbCalcP2` | live Dataset0 call from `RDP5CL.exe` | numerically equivalent; one final-bit `pow` rounding difference |

The live `AlistRDP4` fixture covers 2,300 triplets and changes 709 redo states.
Because the exact C++ implementation calls `FastRecCheckPB` internally, that
passing outer boundary also validates its complete internal chain collectively;
the later rows isolate the same routines at their separately exported calls.

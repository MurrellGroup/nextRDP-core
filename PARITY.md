# Routine parity ledger

Oracle means the pinned locally compiled `DNA5.dll` invoked by the pinned
`RDP5CL.exe`. All discrete state comparisons are byte-for-byte. Probability
outputs use the user's accepted platform-math tolerance when Windows and WASM
`pow` differ only in final rounding.

| Execution order | Routine | Oracle boundary | Current result |
|---:|---|---|---|
| 1 | `SuperDistP` | first live Dataset0 row call from `RDP5CL.exe` | byte-identical result and distance state |
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
| 13 | `FindSubSeqPB4` | live Dataset0 call from `RDP5CL.exe` | byte-identical result and position maps |
| 14 | `FindFirstCOP` | live Dataset0 call from `RDP5CL.exe` | byte-identical result |
| 15 | `ProbCalcP` | live Dataset0 call from `RDP5CL.exe` | numeric result gate |
| 16 | `CleanXOSNW` | live Dataset0 call from `RDP5CL.exe` | byte-identical cleared state |
| 17 | `MakeTestPVs` | live Dataset0 call from `RDP5CL.exe` | byte-identical event-test state |
| 18 | `FindBestRecSignalP2` | live Dataset0 call from `RDP5CL.exe` | byte-identical selected event and state |
| 19 | `UFDist` | live Dataset0 call from `RDP5CL.exe` | byte-identical interval distances |
| 20 | `SuperDistP2` | live Dataset0 call from `RDP5CL.exe` | byte-identical distance matrices |
| 21 | `CheckMatrixP` | live Dataset0 call from `RDP5CL.exe` | byte-identical matrix filtering state |
| 22 | `MakeNJTreesP2` | live Dataset0 call from `RDP5CL.exe` | byte-identical tree-building state |
| 23 | `MarkOutsides` | live Dataset0 call from `RDP5CL.exe` | byte-identical event marking state |
| 24 | `MakeSDMP2` | live Dataset0 call from `RDP5CL.exe` | byte-identical regional distance state |
| 25 | `FillRmat` | all three live Dataset0 calls from `RDP5CL.exe` | byte-identical correlation-matrix state |

The live `AlistRDP4` fixture covers 2,300 triplets and changes 709 redo states.
Because the exact C++ implementation calls `FastRecCheckPB` internally, that
passing outer boundary also validates its complete internal chain collectively;
the later rows isolate the same routines at their separately exported calls.

The CLI now constructs `Seqnum`, the complete triplet list, and `CompressSeq`
directly from the supplied Dataset0 FASTA. Those three arrays are byte-identical
to the live state passed by RDP into `AlistRDP4`; they are no longer supplied to
that preprocessing gate by a fixture.

The CLI also ports the scalar VB/DNA distance setup preceding the first
`SuperDistP`: site categorization, the seven `NucXX` categories, RDP's unusual
fixed-width sequence compressors, the sequential `SuperDistP` row loop, and
the caller's diagonal/low-validity normalization. Starting from FASTA, the
resulting complete 25-by-25 distance matrix is byte-identical to the live
matrix passed into `AlistRDP4`. The legacy helper bodies are implemented but
remain collectively gated at this shared-state boundary until their individual
`DNA.dll` calls are captured.

The following UPGMA block now starts from that generated distance matrix and
ports `MakeDistanceBakB`, `MakeDistMapX`, `ShortestDistB`,
`AddSeqToUPGMA`, `UpdateDistMapX`, and `TreeDist2` in their VB call order.
The resulting complete `TreeDistance` matrix is byte-identical to the live
matrix passed into `AlistRDP4`.

With those generated matrices substituted for captured state, `AlistRDP4`
still produces the byte-identical redo list and `StoreLPV`. The first redo
triplet selected from that list drives generated `CompressSeq` through
`FindSubSeqPB3`, `XOHomologyP`, the intervening VB average-homology role
ranking and distance tie breaks, and `FindNextP`. The connected prefix is
byte-identical at every captured state boundary on all ten supplied datasets;
run `tools/check-all-pre-scan-parity.sh` for that table.

The orchestration preserves two otherwise easy-to-collapse state values:
RDP passes the FASTA length to `FindSubSeqPB3` but the FASTA length plus one to
`XOHomologyP`, and its initial tree-selected `HighHomol` can differ from the
average-selected `HighHomol` later passed to `FindNextP`.

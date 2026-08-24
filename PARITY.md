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
| 26 | `CalCR` | all three first-event calls, fresh local-RDP capture on all ten datasets | byte-identical `RCorr`, `RInv`, and `tRCorr` state |
| 27 | `MakeProperRCorr` caller tail | first-event `MakeRList` input, fresh local-RDP capture on all ten datasets | byte-identical filtered correlations, inversion categories, and warning mask |
| 28 | `QuickDist6` / `MakeDMatS` | first-event `MakeRList` warning input on all ten datasets | byte-identical downstream warning state; literal inclusive-window source port |
| 29 | `MakeGoodC` | first-event `MakeRList` input on all ten datasets | byte-identical overlap eligibility |
| 30 | `MakeINList` / `MakeACOR` | first-event `MakeRList` input on all ten datasets | byte-identical role map and acceptable-correlation mask |
| 31 | `MakeRList` | first live call from each of ten fresh local-RDP runs | exact row bounds, candidate membership, list order, and inversion flags; accepted platform-math tolerance for derived probability cells |
| 32 | `FindActualEvents` / `StripUnfound` | all three first-event role calls from fresh local-RDP runs on all ten datasets | exact candidate inputs, found-slot sets, overlap scores, breakpoint matches, swap-last pruning, and structural `OKSeq` state; accepted p-value tolerance |
| 33 | `StripDupInv` | first live call from each of ten fresh local-RDP runs | exact pre/post row bounds, list order, inverse removal, and inversion penalties |
| 34 | `MakeLDist` / `MakeRCompat` | every first-event call (6–24 calls per run) from fresh local-RDP runs on all ten datasets | byte-identical list distances, compatibility/reverse-compatibility scores, and non-recombinant list bounds across primary, tie-break, alternate-list, and retrim families |
| 35 | `MakeDoneThis3` / `MakePhPrScore` | all first-event score calls from fresh local-RDP runs on all ten datasets | exact filter/trace state and numerically equivalent correlation/subscores |
| 36 | `MakeTrpGroups2` / `MakeTrpScore2` | all three role calls from fresh local-RDP runs on all ten datasets | byte-identical group membership/counts and numerically equivalent triplet tree scores |
| 37 | `FindSets` and pre-retrim compatibility orchestration | complete first-event path from fresh local-RDP runs on all ten datasets | exact alternate candidate sets and exact conditional 6-, 9-, or 24-call `FAMat`/`FCMat`/`SAMat`/`SCMat` family order |
| 38 | `CheckPatternX` / `CheckPattern` score panel | first live call from each of ten fresh local-RDP runs | byte-identical per-sequence pattern counts and completion mask; integrated VB proportion-panel update |
| 39 | `FinalTrim` candidate maintenance | first retrim path from fresh local-RDP runs on all ten datasets | exact active candidate bounds and membership on every dataset that retrims |
| 40 | `MakeVarMap2` / `MakeCntHit2` / `CalcMatchY` | native DLL trace plus connected first-event path | byte-identical hit products and breakpoint smoothing state; both native grouping-threshold branches ported |
| 41 | `ConsensusOK` | post-retrim candidate state from fresh local-RDP runs on all ten datasets | exact active candidate bounds, membership, and list order |
| 42 | second `FinalTrim` strict-group block | first-event path from fresh local-RDP runs on all ten datasets | exact active candidate bounds, list order, synthetic RDP group records, and RDP rescan tail |
| 43 | `MakeRelevant` / `MakeCollecteventsC` | both parent calls for the first event from fresh local-RDP runs on all ten datasets | exact relevant membership, nonblank RDP working records, per-row counts, and collected parent-event outputs; accepted platform-math tolerance for p-values |
| 44 | `FillFSSGC` / `FindSubSeqGCAP7` | every initial GENECONV triplet through the first `MakeTestPVs` boundary | exact compressed informative-site stream, difference counts, and coordinate maps on all ten datasets |
| 45 | `GetFragsP` / `GetMaxFragScoreP` / `CalcKMaxP` / `GCCalcPValP2` | targeted native peak traces plus the complete initial GENECONV scan | exact peak order and accepted-fragment topology, preserving the independent `GCDimSize=2999` fragment stride |
| 46 | `GCGetHiPValP` / `DelPValsP` / `MakeDeleteArrayP` / `GCXoverD` caller loop | complete first `MakeTestPVs` state with only RDP and GENECONV enabled on all ten datasets | exact counts, rows, event identities, ordering, and p-values except one accepted platform-math difference on Dataset8 |

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
ranking and distance tie breaks, `FindNextP`, `FindFirstCOP` when the circular
start branch applies, and `DefineEventP2`. Rejected intervals advance through
the source `NextPosX` loop and the probability-table prefilter until the first
literal `ProbCalcP` or `ProbCalcP2` branch is reached. The scan then preserves
the source's two role transitions and its `OldX` state across cycles until the
first corrected-significant interval is found. Generated `AH`, `XoverSeqNumW`,
`XDiffPos`, and `XPosDiff` at the resulting `FindSubSeqPB4` call are
byte-identical on all ten datasets (with the accepted platform-math tolerance
for final p-values). The port then exhausts the remaining two role scans and
matches the complete `XoverSeqNumW` input and cleared output at the final
`CleanXOSNW` call on all ten datasets; run
`tools/check-all-pre-scan-parity.sh` for that table.

The same source-order walk now runs every triplet selected by the generated
`AlistRDP4` redo list. It ports the late `FindSubSeqPB4` map construction,
compressed-coordinate conversion, the `XPDDone` branch of `CentreBP`, RDP's
`CurrentXOver`/`StoreLPV` event-owner choice, and `UpdateXOList3`'s ordinary
append path. Across all ten datasets, all 16,431 generated active records have
the exact native daughter, parents, centered breakpoints, program/SBP flags,
and zero-initialized metadata at the first `MakeTestPVs` boundary. Every
per-sequence `CurrentXOver` count also matches. The generated arrays then pass
through the exact `MakeTestPVs` and `FindBestRecSignalP2` routines and select
the same first event on all ten datasets. That selection is now exposed as a
fixture-independent library routine, including the source's fallback from
topology-changing (`DoneTarget=0`) to non-changing (`DoneTarget=1`) signals.
Run
`tools/check-all-event-state-parity.sh` for the current table. Probability
cells retain the accepted Windows/WASM final-rounding differences and are not
part of the structural identity gate.

The initial mixed-method path now follows the RDP scan with the literal
compressed `GCXoverD` path.  Its generated event state matches the locally
compiled oracle on all ten datasets: totals are 1,400, 663, 987, 1,144, 775,
3,366, 6,008, 3,799, 3,797, and 4,556 respectively, with every row count and
all 26,495 structural event identities exact.  Dataset8 has one accepted
last-bit probability difference; every other probability cell is exact.

The selected event now continues through the first phylogenetic test. The
generated full-alignment state and selected-region state match the native
`UFDist` and `SuperDistP2` boundaries on all ten datasets. The port then
executes the installed `DNA.dll!FinishDists` evaluation order, the VB
`CheckMatrixX` masking prepass, `CheckMatrixP`, sequence removal and trace-list
construction, and `MakeNJTreesP2`. Both generated NJ tree strings and all
discrete outputs match on all ten datasets. Tiny platform `log` differences in
otherwise equivalent stored Jukes-Cantor cells are compared at `1e-6`; the
delayed-rounding implementation is required because rounding the apparent
source `float` temporaries changes the native tied-tree merge order.

From the same generated event and alignment, the port now also executes the
literal circular `MakeBPosLR` flank walk, constructs the five VB correlation
regions and fixed triplet comparison matrix, and reaches `MakeSDMP2` with the
exact native input state. Its complete `SDM` and `DistMat` outputs and the
three sequential `FillRmat` branches are byte-identical on all ten datasets.
The resulting matrices then pass through a direct port of all three `CalCR`
calls. Its complete `RCorr`, `RInv`, and five-permutation `tRCorr` state is
byte-identical to a reproducible installed-`DNA.dll` oracle capture on all ten
datasets.

The connected path now continues through the literal VB tail of
`MakeProperRCorr`, including SDM suppression, both warning passes, six-decimal
tested-correlation rounding, and inversion-category remapping. It then ports
the old-DNA `QuickDist6`, `MakeDMatS`, and `MakeGoodC` loops, constructs the
compact three-row matrix panels used by `MakeINList`/`MakeACOR`, and executes a
direct source port of `MakeRList`. Fresh first-call traces from the locally
compiled RDP match all structural inputs and outputs on all ten datasets.
The same connected path now executes Module3's active `FindActualEventsVB`
wrapper, the compiled `FindActualEvents` body, `StripUnfound`, and
`StripDupInv`, preserving their swap-last ordering and two confirmed legacy
indexing bugs. It then runs `MakeLDist` and directly checks every observed
`MakeRCompat` call, including the optional FCMat/SCMat, alternate-list, and
post-retrim families (6–24 calls per dataset). The connected score path also
ports `MakeDoneThis3`, `MakePhPrScore`, `MakeTrpGroups2`, and `MakeTrpScore2`.
Every routine boundary matches fresh oracle runs on all ten datasets; carried
probability and correlation-score cells use the approved `1e-12` tolerance.

The orchestration preserves two otherwise easy-to-collapse state values:
RDP passes the FASTA length to `FindSubSeqPB3` but the FASTA length plus one to
`XOHomologyP`, and its initial tree-selected `HighHomol` can differ from the
average-selected `HighHomol` later passed to `FindNextP`.

The connected path now also enters `FinalTrim` when the native compatibility
scores request a retrim. It ports the active candidate-maintenance loop,
`MakeVarMap2`, `MakeCntHit2`, `CalcMatchY`, and `ConsensusOK`, including the
configured conservative grouping thresholds. The resulting active `RNum` and
`RList` ranges match the post-retrim `MakeRCompat` boundary on all four demo
datasets that take this path; the other six correctly skip it. As in VB, cells
beyond each `RNum` are stale backing-array storage and are not observable list
members.

Immediately after that retrim path, the clean port now reproduces the active
`CalcMaxD` prelude and its `CMaxD2P3` boundary. It constructs the same
representative and included-sequence masks, informative-site nucleotide maps,
pairwise V-score matrix, and accumulated maximum-distance values. Every input
and output cell at that boundary matches fresh locally compiled RDP traces on
all ten datasets.

The connected first-event path now also ports `MakeSSDistB`, `SimpleDist`,
`MakeOUCheck`, `MakeEList`, `MakeListCorr`, `GetBadDists`, the two post-retrim
`MakeRCompat` families, and the active decision-tree path of
`MakeConsensusC`. Every generated input statistic agrees with the native
`NN_inputs` boundary on all ten datasets at its printed precision, and the
raw consensus scores and selected recombinant role agree on all ten. The CSV
oracle is used only as a boundary assertion; no value from it drives the
connected consensus calculation.

The first-event path now continues through the second, strict-group
`FinalTrim`, its synthetic group-member records and relaxed primary-RDP
rescans, `MakeRelevant`, and both parent-role `MakeCollecteventsC` calls.
Every discrete RDP record and collected event agrees across all ten datasets.
The exact `AddjustCXO` redistribution boundary is now also a reusable library
routine rather than fixture-runner-only code; its retained event matrix,
done-state matrix, and pair-rescan schedule remain covered by the same ten
dataset gate.
The immediately following inner- and outer-pair rescan screening is now also
fixture-independent library code. It preserves the source `AlistRDP3` call,
including the flattened triplet order and correction threshold, and remains
exact through both post-rescan event sets, the second `MakeTestPVs`/selection,
and the second resolved-round prefix on all ten datasets.
Native zero-filled slots created by disabled-method bookkeeping are excluded
from the RDP-event fixture because they contain no event and cannot overlap a
real interval. The enclosing post-event sequence-mutation and cyclic rescan
loop remains the next parity boundary.

The initial method stream now also enters MaxChi using the native dispatcher’s
captured triplet order. Runtime traces establish that this path uses the fixed
70-site window and `FindAllFlag=0`: the first `GrowMChiWinP2` call gets
`TWin=10`, not the `TWin=36` implied by the supplied VB caller’s assignment of
`FindAllFlagX=1`. Reproducing the observed flag, the source `CalcChiVals4P3`
lookup path, compressed informative-site map, and initial `CentreBP` branch
reduces Dataset0 MaxChi output from 1,286 port records to 636 versus 634 in
the oracle. A direct trace of the first remaining divergent growth call then
showed that circular `MCXoverF` does not use that lookup path for its initial
scan: it calls formula-based `CalcChiValsP`, preserving a full-precision
double peak value. The port's rounded float-table peak compared equal to a
later window and reset `GrowMChiWinP2`'s failure count, while the native peak
compared strictly greater. Restoring the source's circular/formula branch
removes that event and advances the first count divergence from call 41 to
call 341. Dataset0 is now 2,033 port records versus 2,034 in the oracle, with
24 of 25 event-bearing sequence rows exact. The remaining call-341 miss came
from collapsing MCXoverF's explicit `MaxX = 1` opposite-boundary branch into a
generic circular-position helper. That changed the `DestroyPeaks` interval
from `(1203,1)` to `(1,1492)` and erased an adjacent peak. Reinstating the
literal source branch restores Dataset0 to 2,034/2,034 records, 25/25 row
counts, and 2,034/2,034 probabilities. Breakpoint/window structural identity
then reaches 2,034/2,034 only after the active program-3 branch of
`FindDaughter` and MCXoverF's second `CentreBP` call are preserved. Despite
its name, this branch does not change the selected roles: it walks the
informative-site map according to `MCMaxY` and `LoHiFlag`, conditionally
polishes both endpoints, and recentres them. Dataset0's complete initial
RDP+GENECONV+MAXCHI stream is now exact in total, row counts, structural
identity, and probability.

Fresh mixed-method fixtures now cover all ten demos and can be checked with
`tools/check-all-rdp-geneconv-maxchi-parity.sh`. Dataset0, Dataset2, and
Dataset7 are completely exact. Current signal totals are: Dataset0
2,034/2,034; Dataset1 1,061/1,062; Dataset2 1,457/1,457; Dataset3
1,741/1,744; Dataset4 1,103/1,104; Dataset5 4,310/4,313; Dataset6
8,045/8,045; Dataset7 4,672/4,672; Dataset8 5,318/5,326; and Dataset9
5,763/5,764. Dataset6 has equal totals but 8,010 of 8,044 comparable
structural identities; the other non-passing datasets still have one to eight
missing signals plus later coordinate differences.

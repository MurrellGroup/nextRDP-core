# RDP bugs and quirks

This file records behaviour that the WASM port must preserve while parity is
the goal.  It is intentionally separate from the eventual “clean” semantics
we may want to expose to users.

- `InnerScan2` is a post-mutation pass.  After `FindActualEvents` erases the
  selected tract and adjusts the event state, RDP expands the selected RList
  again and runs each enabled legacy method in source order.  Running an
  optional method only during the initial scan misses events that appear after
  the mutation.
- The legacy method scheduler uses source program bits and `MakeAListISP3`
  filtering.  It is not equivalent to iterating the permanent analysis list
  directly; `Worthwhilescan`, `ProgBinRead`, `TraceSub`, `DoPairs`, and the
  `ProbDo` sampling state all affect which triplets are called.
- The optional-method role allocator has a strict-count branch before its
  StoreLPV fallback.  Equal sequence counts therefore do not select a role in
  that branch; the fallback then decides the orientation.  This differs from
  the superficially similar inclusive comparisons in some VB helper routines,
  and must not be normalized without checking the complete caller state.
- Inner MaxChi and Chimaera scans rebuild `MDMap`/`BanWin` with `MakeBanWinP`
  in informative-site coordinates after a mutation.  Passing raw alignment
  coordinates, or using the unmasked `CalcChiVals3` path, changes the source
  `CalcChiVals5`/`GrowMChiWin2` window decisions.  The port preserves this
  source ordering even when a demo alignment has no missing sites.
- The source can enter round-prefix identification with every sequence marked
  for redo.  In that state it has no NJ panel; passing a negative local upper
  bound into `MakeNJTreesP2` is an underflow/crash in the port.  The compatible
  path skips the NJ call for an empty panel and keeps a neutral one-cell state
  for downstream code.
- 3Seq `GetTSPVal` always uses the configured lookup table.  At or beyond the
  table boundary it rescales `(nM,nN,nK)`, truncates with VB `Int`, then raises
  the looked-up value by the rescale factor; it does not switch to a direct
  Siegmund approximation when a table lookup is available.
- 3Seq has two source call modes: the first analysis pass calls `TSXOver(0)`
  (`FindSubSeqTS`, one selected excursion), while post-event rechecks call
  `TSXOver(1)` (`FindSubSeqTS2`, allowing the reverse excursion).  Treating
  every call as the latter changes event order even when the walk itself is
  identical.
- `StripUnfound2` in the DNA5 path uses the legacy informative-site bounds and
  can leave the final informative site out of the strip loop.  Do not “fix” the
  bound while reproducing RDP; record it here instead.

- `PolishBP` invokes `BenHMM` with a fixed `HMMCycles = 20`.  This is
  independent of `MinSeqSize`, which is a separate sequence-size gate used by
  the surrounding event/tree logic.  Passing `MinSeqSize` as the HMM cycle
  count is a port defect, not an RDP optimization: on a 9,594-site alignment
  it changed 20 cycles to 96 and materially increased BURT runtime.  The
  compatibility port keeps the source's fixed 20-cycle call.

- The DNA5 DLL's `AlistGC2`, `AlistMC3`, and `AlistChi` routines use OpenMP
  worksharing over triplet rows.  Their per-row scratch buffers and `RL[y]`
  outputs are independent, so deterministic range partitioning reproduces the
  serial results while restoring the intended parallelism in pthread WASM.
  `BenHMM` is different: the active source path explicitly calls
  `DoHMMCyclesSerial`, so BURT must remain serial even when method screens are
  threaded.  A cleanup should expose these two policies instead of treating
  every source `#pragma omp` as interchangeable.

- `AlistRDP4` has the same independent-row OpenMP structure.  In the
  Emscripten pthread build it is invoked over deterministic contiguous row
  ranges; native OpenMP builds invoke the vendored routine once so its own
  source workshare is retained.  The resulting `Redo` bytes and downstream
  event order are unchanged.

- RDP's optional methods can be selected without the RDP method.  The source
  still builds the shared StoreLPV screening table, but it does not run the
  RDP XOver walk or the cyclic tract-erasure scheduler in that mode.  The
  compatibility path therefore emits the selected GENECONV/MaxChi/CHIMAERA/
  3SEQ records directly in their source method order; enabling RDP must not
  be inferred merely because an optional method was selected.

- The RDP `DrawPlots`/`XOver` profile is an informative-site homology plot,
  not a percent-identity scan over every nucleotide.  `FindSubSeqP` or PB3
  retains only the three-sequence informative sites, `XOHomologyP` writes
  rolling integer agreement counts, and `DrawPlots` divides them by the odd
  width `2 * Int(XOverWindowX / 2) + 1` even when gaps or invalid bases mean
  fewer sites were compared.  The x-axis is `XDiffPos` in the full alignment;
  `XDiffPos(0)` is replaced with `XDiffPos(1)` when it is still the zero
  sentinel.  Replacing this with raw alignment identities makes the plot look
  smooth/high and is visibly unlike the original RDP plot.

- BootScan's automated `BSXoverR` screen counts only strict, unique closest-pair
  votes. A window tied for the closest distance is not a vote for either pair;
  treating ties as a shared vote inflates support and can create a false tract.
  The review curve must use the same strict rule as discovery.

- SISCAN's fast `QuickCheckB` window pass deliberately uses the source's
  one-site-short window bound before the full `ShrinkRegionC` pass. It also
  consumes a single seeded flat `MakeVRand` prefix across windows and triplets.
  Recomputing a fresh random stream per window or changing the bound to an
  inclusive window silently changes both the selected sister pair and its
  permutation P values.

- The browser's initial BootScan/SISCAN implementation emits their direct
  source-kernel discoveries after the shared initial triplet screen. They do
  not seed the RDP cyclic tract-erasure scheduler yet; enabling RDP therefore
  preserves the RDP event order and appends optional-lane records separately.
  This is an integration boundary, not permission to infer RDP events from an
  optional method.

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
  `DoHMMCyclesSerial`. Its 21 seeded restarts are computationally independent,
  but two pieces of state make naive fan-out change results: one MSVCRT
  `rand()` stream seeds them in source order, and `PathMax` is not reset between
  restarts, so a likelihood encountered at any iteration can stop the next
  restart early. The self-transition initialization also evaluates `1 - iVal`
  at `float` precision before converting to `double` for `log`. The pthread
  port therefore captures the seed stream first, evaluates every restart's
  iteration states independently, and replays the carried-`PathMax`/first-best
  selection in source order. Exact serial/parallel equality was checked for
  every BURT call in the full cyclic runs on all ten supplied datasets.

- `AlistRDP4` has the same independent-row OpenMP structure.  In the
  Emscripten pthread build it is invoked over deterministic contiguous row
  ranges; native OpenMP builds invoke the vendored routine once so its own
  source workshare is retained.  The resulting `Redo` bytes and downstream
  event order are unchanged.

- The source probability-table setup calls `ProbCalcP` independently for
  every `(length, lower-bound, probability)` cell. All lower bounds for a
  fixed length/probability recalculate the same factorial and `pow` terms.
  The compatibility port caches those long-double terms and parallelizes the
  independent probability columns, but deliberately repeats `ProbCalcP`'s
  ascending per-bound addition with a double rounding after every term.
  Replacing it with an ordinary reverse cumulative sum is faster but changes
  floating-point grouping. The retained implementation was checked against
  every source table cell and keeps the full cyclic result hash unchanged.

- RDP's optional methods can be selected without the RDP method.  `DoRDP` is
  still called and its common cyclic selection/tract bookkeeping is entered;
  only the RDP-specific `XOver` walk is guarded by `DoScans(0, 0)`.  During
  selection, records for disabled program bits are assigned probability 1 and
  skipped before the next candidate is considered.  A port must therefore
  keep the scheduler active for optional-only runs without exposing disabled
  RDP records as results (the old direct-emission shortcut was not source
  faithful).

- The RDP `DrawPlots`/`XOver` profile is an informative-site homology plot,
  not a percent-identity scan over every nucleotide.  `FindSubSeqP` or PB3
  retains only the three-sequence informative sites, `XOHomologyP` writes
  rolling integer agreement counts, and `DrawPlots` divides them by the odd
  width `2 * Int(XOverWindowX / 2) + 1` even when gaps or invalid bases mean
  fewer sites were compared.  The x-axis is `XDiffPos` in the full alignment;
  `XDiffPos(0)` is replaced with `XDiffPos(1)` when it is still the zero
  sentinel.  The plotting loop nevertheless starts at index 1, so that
  repaired index-0 sentinel is not itself a displayed point.  Replacing this
  with raw alignment identities makes the plot look smooth/high and is
  visibly unlike the original RDP plot.

- BootScan's automated `BSXoverR` screen counts only strict, unique closest-pair
  votes. A window tied for the closest distance is not a vote for either pair;
  treating ties as a shared vote inflates support and can create a false tract.
  The review curve must use the same strict rule as discovery. Fixed-region
  rechecks also repeat the same seeded pair/window distance profiles used by
  discovery; those profiles may be reused by sequence-pair key, but the cache
  must be invalidated after tract erasure. Reusing that exact profile path
  removed the dominant secondary-check cost without changing any serialized
  BootScan result at one, four, or eight requested workers.

- SISCAN's fast `QuickCheckB` window pass deliberately uses the source's
  one-site-short window bound before the full `ShrinkRegionC` pass. It also
  consumes a single seeded flat `MakeVRand` prefix across windows and triplets.
  Recomputing a fresh random stream per window or changing the bound to an
  inclusive window silently changes both the selected sister pair and its
  permutation P values. The vertical remap depends only on the source category
  family and the 1..12 random value, so a lookup table is equivalent to the
  source modulo branches only while the flat prefix cursor remains in exactly
  the same category/occurrence/permutation order.

- The browser's initial BootScan/SISCAN implementation emits their direct
  source-kernel discoveries after the shared initial triplet screen. They do
  not seed the RDP cyclic tract-erasure scheduler yet; enabling RDP therefore
  preserves the RDP event order and appends optional-lane records separately.
  This is an integration boundary, not permission to infer RDP events from an
  optional method.

- PHYLPRO's intended profile is indexed by polymorphic alignment columns, not
  by every nucleotide coordinate.  Each target keeps rolling mismatch counts
  against the non-disabled context rows on the two sides of a moving partition;
  Pearson is then computed after excluding the target's own context row unless
  the explicit self option is enabled.  The source uses VB CInt half-window
  rounding, circular wrap-around, and source-default correlation `1` for
  degenerate rows.  A browser plot may downsample the returned points, but it
  must retain both breakpoint-nearest points and each role's global minimum.

- Final event records retain representative sequences in discovery order while
  `winningRole` identifies the source's selected recombinant row.  Event
  alignments, trees, and PHYLPRO must rotate that prefix to recombinant/major/
  minor before assigning labels; using the stored prefix directly silently
  relabels evidence whenever the winning role is not zero.

- The old review page has two deliberately different triplet contracts.  The
  event alignment/tree/PHYLPRO panels use the reconciled recombinant/major/
  minor roles, but the signal plot is keyed to the original canonical
  discovery triplet and its pair slots.  Reusing the winner-rotated event
  roles for the plot swaps coloured traces (especially when the winner is
  slot 1 or 2), even though the underlying pairwise counts are unchanged.

- The browser tree endpoint is a source-shaped display reconstruction.  Its
  neighbor-joining edges are useful for review, but they are not a persisted
  RDP5 tree and must not be interpreted as source tree-parity evidence.

- `MakeAnalysisListQvR` is not an ordinary exploratory triplet filter. It
  preserves input order, emits every group-0 query for each unordered pair of
  enabled references from different positive groups, and omits same-group
  reference pairs entirely. The desktop correction bookkeeping still derives
  its opportunity count from reference-group pairs × queries, which can differ
  from the number of emitted triplets when groups have uneven sizes.

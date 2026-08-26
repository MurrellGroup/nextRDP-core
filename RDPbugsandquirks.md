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
- The source can enter round-prefix identification with every sequence marked
  for redo.  In that state it has no NJ panel; passing a negative local upper
  bound into `MakeNJTreesP2` is an underflow/crash in the port.  The compatible
  path skips the NJ call for an empty panel and keeps a neutral one-cell state
  for downstream code.
- `StripUnfound2` in the DNA5 path uses the legacy informative-site bounds and
  can leave the final informative site out of the strip loop.  Do not “fix” the
  bound while reproducing RDP; record it here instead.

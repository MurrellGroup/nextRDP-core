# Core parity rules

- This project is a clean-room orchestration port around the supplied RDP source snapshot.
- The compatibility target is the pinned local `RDP5CL.exe` plus locally rebuilt `DNA5.dll` recorded in `ORACLE.md`.
- Port routines in the order the RDP-only command-line execution reaches them.
- A routine is not considered ported until its observable outputs match a captured oracle fixture.
- Do not add dataset-specific conditions, lookup tables, role swaps, or breakpoint corrections.
- Preserve native RDP bugs and quirks. Record confirmed quirks in `RDPbugsandquirks.md`.
- For every confirmed source/runtime behavior needed for parity, record the mechanism and the
  validating dataset/oracle evidence in that ledger; do not hide a compatibility quirk in a
  dataset-specific workaround.
- Keep generated traces, fixtures under construction, and runtime output in this project's `sandbox/` directory.
- Keep the command-line core authoritative; web integration may consume it when the
  user explicitly requests interface functionality, but must not introduce separate
  algorithmic behavior.

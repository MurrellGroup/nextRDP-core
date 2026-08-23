# Core parity rules

- This project is a clean-room orchestration port around the supplied RDP source snapshot.
- The compatibility target is the pinned local `RDP5CL.exe` plus locally rebuilt `DNA5.dll` recorded in `ORACLE.md`.
- Port routines in the order the RDP-only command-line execution reaches them.
- A routine is not considered ported until its observable outputs match a captured oracle fixture.
- Do not add dataset-specific conditions, lookup tables, role swaps, or breakpoint corrections.
- Preserve native RDP bugs and quirks. Record confirmed quirks in `../nextRDP-wasm/RDPbugsandquirks.md`.
- Keep generated traces, fixtures under construction, and runtime output in this project's `sandbox/` directory.
- Do not add UI or web integration until the command-line core reaches parity.


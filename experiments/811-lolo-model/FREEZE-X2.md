# FREEZE (coordinator, 2026-10-22T08:13Z by the clock): do not touch work/banks or hook811.py until X2 is submitted or retired

GO-X2 written at c3f24c91. X2 = C0 (R0 + 811 via dumphook811c0.sh) builds after X1, ~11:30Z → ~16:30Z on 10-22. Its hook reads
`work/banks` (manifest 4cc6d82a, committed by the coordinator at 77c50953) and `hook811.py` (c7cff26d). Any change to either
breaks the build's pins. Lift this note only after the coordinator records X2's submission or retirement.

## Extended 2026-10-22T17:50Z by the clock: X2 SHIPS — freeze until the final closes (Nov 5)
The coordinator shipped X2 = C0 as the replacement final entry (6a decision 63803087; package READY e39d2ba6; the dump used == the hook
on R0's stack dump by content). 811's banks (work/banks, manifest 4cc6d82a), hook811.py (c7cff26d), dumphook811c0.sh (297c3832) and
811's records (RESULTS, PREREG-*, wear811.jsonl, run-comp811kit-*, run-comp811-LZ, defbanks/, kamac/) stay unchanged until the final
closes on 2026-11-05. The old bank set work/banks.pre-def-37e34335 is kept too.

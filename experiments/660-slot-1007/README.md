# 660 — the 10-07 slot chain (coordinator, 2026-10-02)

- `chain1007.sh` waits for 10-06 20:00Z.
  - It checks the sha of the frozen rule (`experiments/653-rule1007-review/rule1007_v4_653.py`).
  - It reads the captured entries of 446bj2zcp (10-01), 495 (10-03), 496a (10-05) and 446bj2zk16 / 446bj2zk20 (10-06).
  - It runs the rule WITHOUT `--arc-hook` and arms `arm1007.sh` with the first READY package of each slot's list.
- `arm1007.sh` is arm1004.sh with the dates moved to 10-07: submit after the 00:01Z reset, cap 12:00Z, 536 retry block kept.
- **R6:** if the coordinator's read of 524dr (10-06 s1) is LAG or hits a falsifier, write `work/HOLD-1007-s1` before 20:00Z
  10-06. Slot 1 is then not armed, and slot 2's pick goes alone.
- `NOGO-1007-s<n>` in `work/` stops a slot before its submit, as in every arm script.
- Names were claimed from docs/MODEL-NAMES.md's queue (446bj2zcph at 10-02 20:01Z; 446bj2zk20 and 459dr earlier). `-d` is empty.
- `DRY=1 bash chain1007.sh` prints the plan without waiting or arming.
  - 10-02 20:03Z dry run: s1 446bj2zcph, s2 459dr (k16 not yet read). Expected on 10-06: s1 CPH; s2 K20 if k16 ≥ MID, else 459dr.
- Armed 10-02 20:03Z (pid in `work/chain1007.pid`). After a reboot, re-arm with the same command:
  `nohup bash experiments/660-slot-1007/chain1007.sh > experiments/660-slot-1007/work/chain1007.nohup 2>&1 &`

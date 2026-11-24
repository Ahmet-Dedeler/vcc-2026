# 659 — the 10-08 slot chain (DRAFT for the coordinator; UNARMED)

[7772dc], assigned 10-04 07:3xZ. Modelled on 660's chain1007 / arm1007. The coordinator reviews and arms it.

## Files
- **`chain1008.sh`** waits for 10-07 20:00Z (a late wake after 10-08 06:00Z arms nothing). Then:
  1. It checks the frozen rule's sha (653 v4, `e00af6c7…`, the same pin as chain1007).
  2. It reads the captured entries with `slot_entry_id.py` over `experiments/*/work/submit100[1-7]s[12].log`:
     - cp 446bj2zcp, h 496a, k16, k20 (scored 10-04 or later), c 495, cph 446bj2zcph.
  3. It derives `--scored=` from 660's `submit1007s1/s2.log`.
     - A slot whose LAST armed run reached `SUBMITTED ok` or `CAPTURED own entry` scored its tag.
     - A submitted but uncaptured entry counts. That is stricter than ARMED+CAPTURED, so it never re-arms a scored package.
  4. It runs the rule WITHOUT `--arc-hook` and WITH `--scored=` (empty if nothing scored), then builds each slot's
     candidate list:
     - **Default:** the rule's `10-08:` list (rule tag → build tag: CPH 446bj2zcph, K20 446bj2zk20, 459dr, ALLSW allsw658,
       CPHK20 446bj2zcphk20), then the queue tail from `tail1008.tsv` in row order (630arm → 524dr2 → allsw658), deduplicated.
     - **Opt-in `work/QUEUEFIRST-1008`:** the rule's DE follow-ups (CPHK20, CPMZ, K20, MZ) first, then the tail, then the rule's
       own fallbacks (459dr, ALLSW, S1).
  5. It arms `arm1008.sh` with the first READY candidate of each slot. Slot 2 never takes slot 1's pick.
- **READY:**
  - **Rule-list arm:** package present; `[tag] PASS`, `[tag] CLIFF OK` and `novelty rc=0` in its gate log; a claimed name
    (`docs/model-names-claimed.tsv`); the novelty check passes now.
  - **Tail arm:**
    - its driver log's last line is `<ts> READY <tag> <sha>`;
    - the driver's `package <path> sha256 <sha>` line names the package, and shasum == sha;
    - the gate log has `[gate_tag] PASS`, `[gate_tag] CLIFF OK`, `novelty rc=0` and `package sha256: <sha>`;
    - a claimed name exists, and the novelty check passes.
- **`tail1008.tsv`** (tag, rule_tag, driver_log, gate_log, gate_tag) is read at 20:00Z 10-07, so the rows can change until
  then without touching the live script.
  - 630arm: [d4cd5f]'s driver.
  - 524dr2 and allsw658: 658's drivers. 524dr2's gate tag is `524dr`, because gate524.py writes `[524dr]`.
- **`arm1008.sh`** is arm1007.sh with every date moved to 10-08 (reset 00:01Z, cap 12:00Z; the 536b retry block kept).
  - The one functional change: an optional 7th argument, the gate log's bracket tag (default = the arm's tag), for 524dr2.
- **Files that stop a slot:**
  - `work/HOLD-1008-s1` or `-s2` drops that slot's list (R6 style). If slot 1 is held, slot 2's pick goes as slot 1.
  - `work/NOGO-1008-s<n>` stops a slot inside arm1008.sh.
- **Names** come from the MODEL-NAMES queue, claimed 10-04 09:04Z for 630arm, 524dr2 and allsw658; 459dr, K20 and CPH had
  theirs already. `-d` stays empty.

## The order question (for the coordinator)
- The frozen rule's 10-08 list ends with its own fallbacks: 459dr, then ALLSW.
- The likely state on 10-07: K20 is scored 10-04 (if its resume lands), so 10-07 s2 = 459dr.
- The default order then gives 10-08 = **allsw658, then 630arm**, and 524dr2 waits for 10-09.
- **With `QUEUEFIRST-1008` it is 630arm, then 524dr2**, the queue's order. That reads 630 first (DNMAG waits on 630's entry)
  and the October rehearsal early, and ALLSW still lands before Oct 22.
- Create the file before 20:00Z 10-07 if you want the queue order.

## Tests and dry runs
- `uv run python experiments/659-chain1008/test659.py`: **13/13**. It uses synthetic 10-07 logs, tail drivers, gate logs and
  packages, with novelty stubbed (DRY only). Covered:
  - scored derivation, candidate order and dedup;
  - HOLD s1 and both;
  - sha mismatch, missing gate sha, a NOT READY driver, a TBD row, the gate tag;
  - QUEUEFIRST;
  - `bash -n` on both scripts.
- **arm1008** with a wrong gate tag stops at the gate check (`NOT SUBMITTING: config diff did not PASS`). The test's slot log
  was deleted.
- **Live DRY on 10-04 09:0xZ** (`DRY=1 bash chain1008.sh` → `work/dry-chain1008.log`):
  - inputs cp K1Q0…, c 9Hxef…, everything else none; `--scored=` empty;
  - rule `10-08: 459dr ALLSW`; would arm s1 459dr, s2 none, because no tail arm is built yet.
- DRY mode writes only `work/dry-chain1008.log`.

## Arm (coordinator)
```bash
nohup bash experiments/659-chain1008/chain1008.sh > experiments/659-chain1008/work/chain1008.nohup 2>&1 &
```

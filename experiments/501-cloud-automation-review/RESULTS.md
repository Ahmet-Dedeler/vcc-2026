# 501 — adversarial review of the 2026-09-27 slot and build automation

Reviewed on current `master` (`14695cef`), from the chain sources only. Reproductions are under `repro/` and do not read gitignored data. Nothing here was applied to the armed scripts.

| # | sev | where | failure |
|---|---|---|---|
| 1 | high | `arm1002.sh:20`, `arm1003.sh:20`, `arm1004.sh:20` | Slot-1 noon cap never fires. A late wake submits that day's package on a later day. |
| 2 | high | `chain1004.sh:15-19` | `eid()` returns on the first ARMED hit. A failed earlier arm hides the real entry, or a later CAPTURED line is attributed to the wrong tag. `rule1004` then treats 495/479/453 as not submitted. |
| 3 | high | queue deadlines, see §3 | Five caps are at or before the predecessor they wait on. The successor logs `END` and exits while that predecessor is still allowed to run. Nobody restarts a dead committer. |
| 4 | high | `relaunch479b.sh:9-10`, `chain495.sh:19-20` | The relauncher treats "495 is not running" or a stale `patched:` line as a pass and starts 479b. 495's start gate matches `chain479b.sh`. Both then wait for each other. |
| 5 | high | `commit_after_chain.sh:45-47` | The committer exits 0 on the first FILE that already equals HEAD, before `git commit`. Every other dirty FILE is left uncommitted. The next chain waits for HEAD until its own cap. |
| 6 | high | `patch_builder479.py` `sharpen_q_by_tercile` | On a tie-free panel the split stays 100/100/100. A cap block that contains the upper cut is demoted to MID, and a stercile tie at the lower cut or a collapsed cut empties LOW or sends every row to the mildest q. |

## 1. Slot-1 noon cap is dead code

`arm1002.sh:20-23` (same shape in `arm1003.sh:20-23`, `arm1004.sh:20-23`):

```bash
until [ "$(date -u +%Y%m%d%H%M)" -ge 202610020001 ]; do
  [ -f "$W/NOGO-1002-s1" ] && stop "..."
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610021200 ] && stop "clock passed 12:00Z ..."
  sleep 20; done
```

The body runs only while the clock is **before** 10-02 00:01Z. The cap is 10-02 12:00Z, nine hours later, so the cap test is false on every iteration that actually runs. When the clock is already past 00:01 — a wake at 15:00, or on 10-03 — the `until` is already true, the body is skipped, and the script submits. Slot 2's cap is a different loop (it waits for a log line) and does fire.

CONFIRMED. `bash experiments/501-cloud-automation-review/repro/slot_cap.sh`:

```
now=202610012359 ... cap_would_stop=0
now=202610020001 ... cap_checked=no  WOULD_SUBMIT=yes
now=202610021500 ... cap_checked=no  WOULD_SUBMIT=yes
now=202610030900 ... cap_checked=no  WOULD_SUBMIT=yes
now=202610051200 slot=202610040001 ... WOULD_SUBMIT=yes
slot2: CAPTURED absent and now=202610021500 >= noon -> STOP (cap fires)
```

A machine that sleeps across noon, or a re-arm after a reboot on the next morning, spends the new day's slot on the previous day's package and can collide with that day's own waiter (one submission in flight).

Fix: `patches/f1-slot-cap.patch`. Check the cap on every wake, and break out only when the clock is inside `[00:01, 12:00)`.

## 2. `eid()` attributes the wrong entry, or none

`chain1004.sh:15-19` walks the 09-28, 10-02 and 10-03 slot logs. On the first file whose ARMED line names the tag, it prints that file's last `CAPTURED own entry` and returns, even when the capture belongs to a different tag or is absent.

That id is the only input `rule1004.py` has. An empty result becomes `none` (`${I495:-none}`), and `C`/`P`/`O` come back `UNKNOWN`. `UNKNOWN` is not `FAILS`, so s1 still takes `498<P>` (`rule1004.py:93-94`). A 495 that actually failed the board is then submitted as the combined arm. The other direction: a re-arm that appends to `submit1002s1.log` makes `eid 479p2g` return 495's id, and the pair pick is scored off the wrong entry.

The arm scripts' own extractor is fine. On `submit0927s1.log` the primary regex `entry [A-Za-z0-9]{20}` matches `submitted — entry YvVru3DAF89Jiyr34olR`. The `entry:` colon form does not match; it does not have to, because the success line does.

CONFIRMED. `bash experiments/501-cloud-automation-review/repro/eid.sh`:

```
eid 495 -> []                          # 10-02 ARMED, no capture; 10-03 capture ignored
eid 479p2g -> [aaaaaaaaaaBBBBBBBBBB]   # that id is 495's, appended later in the same log
A fixed eid 495 -> [ccccccccccDDDDDDDDDD]
B fixed eid 479p2g -> [eeeeeeeeeeFFFFFFFFFF]
```

Fix: `patches/f2-eid.patch`. Keep the capture that follows this tag's ARMED line, ignore a later arm's capture, and continue to the next log when this arm was not captured.

## 3. Deadline inversions drop the rest of the queue

These numbers are wait caps. The happy path does not sleep until them. A successor that gives up while its predecessor is still inside its own cap never retries, and a committer that exits is not re-armed. Later chains then wait for `builder == HEAD` until their own cap and log `END`.

CONFIRMED. `python3 experiments/501-cloud-automation-review/repro/deadlines.py`:

```
INVERSION: commit479 202609271600  vs  chain479 202609271600
         committer dies at chain479's release cap, while the three KAs run after that cap
INVERSION: chain476 202609281200  vs  chain479b 202609281800
         476's awk waits until chain479b.sh is gone; 479b may keep waiting until its own later cap, then still has to build
INVERSION: commit495 202609281200  vs  chain479b 202609281800
         479b will not build until the 495 builder is committed; the committer can die while 495 is still building
INVERSION: chain497 202609281800  vs  commit476 202609281800
         497 requires 476's files at HEAD; both give up at the same minute, and 476's rehearsals run after its release
INVERSION: chain496p 202609291200  vs  chain498 202609291200
         496p waits for chain498 to exit, and both give up at the same minute — 498 is allowed to start at that minute
```

The same shape, one step earlier: `chain495.sh:20` matches `chain479.sh` (`chain479b?\.sh`) and gives up at 06:00Z 09-28 (`chain495.sh:21,24`). `commit479.spec` has already exited at 16:00Z 09-27 if the KAs were still running. From an arm time of 09-27 ~07Z, both of those caps are inside the unattended window. 10-02's primary arms (479) and everything downstream fall through to whatever package already exists.

`chain479.sh`'s own 16:00Z 09-27 start cap is left as it is. The replacement caps are still caps, in dependency order:

| waiter | was | cap |
|---|---|---|
| commit479 | 09-27 16:00 | 09-28 10:00 |
| chain495 (both waits) | 09-28 06:00 | 09-28 12:00 |
| commit495 | 09-28 12:00 | 09-29 06:00 |
| chain479b, relaunch479b | 09-28 18:00 | 09-29 08:00 |
| chain476 (both waits) | 09-28 12:00 | 09-30 02:00 |
| commit476 | 09-28 18:00 | 09-30 14:00 |
| chain497 | 09-28 18:00 | 09-30 18:00 |
| commit497 | 09-29 18:00 | 10-01 10:00 |
| chain498 | 09-29 12:00 | 10-01 12:00 |
| chain496p | 09-29 12:00 | 10-02 02:00 |
| commit496 | 09-30 12:00 | 10-02 10:00 |
| chain497e | 09-29 18:00 | 10-02 14:00 |

479b's packages are read on 10-02, 497's on 10-03 as a fallback, 498's on 10-04. Each of those slots is still after the new cap, and the build itself is not capped once the wait ends. If 495's own KAs and dumps run past 06:00Z 09-29, commit495 still dies and is not restarted; that fuse is 18 h after 495's latest start, not an inversion with 479b.

Fix: `patches/f3-deadlines.patch` (dry-run applies). Stop-line text matches the new clock. `chain495.sh`'s stop line names chain479, because that is the process the regex actually waits on at the start of the queue. `relaunch479b.sh` is also rewritten by §4; apply that hunk instead of this patch's hunk for that one file.

## 4. Relauncher starts 479b before 495 has left its gate

`chain495.sh:19-20` waits while any bash argv matches `chain479b?\.sh`. That matches `chain479b.sh` and does not match `chain495.sh`. `chain479b.sh:14-15` waits while any bash argv matches `chain495\.sh`. It does not match itself. `relaunch479b.sh` matches neither gate.

`relaunch479b.sh:8-10` reads the pid once, then loops until the log contains `patched:` **or** that pid is not running `chain495`. The log is appended (`>>`). Two states that are true at arm time both end the loop:

- 495 has not been started, so `ps -p` fails and the `|| ! ps` term is true.
- A previous run left `patched:` in `chain495.log`.

Either one launches `chain479b` while 495 is still in its start `until`. The two waits then never clear.

Once 495 has logged `patched:` and left that `until`, starting 479b beside it is the intended handoff: 479b keeps waiting until 495's pid exits and no `run_ka_4` / `run_dump4` / `run_arm_generic.sh` is visible, so the two builds do not run together. The bug is the false start.

CONFIRMED. `bash experiments/501-cloud-automation-review/repro/relauncher_deadlock.sh`:

```
495-gate   vs 479b-proc    MATCH
479b-gate  vs 495-proc     MATCH
495-gate   vs 495-proc     no
479b-gate  vs 479b-proc    no
495-gate   vs relauncher   no
495 not started, no patched : START_479b   (BUG: launches 479b)
495 not started, STALE patched in the appended log : START_479b   (BUG)
after 3s, 495 pid 3014 still RUNNING
after 3s, 479b pid 3015 still RUNNING
```

The deadline exit of the current relauncher writes no `END` line (`relaunch479b.sh:11`).

Fix: `patches/f4-relauncher.patch`. 495 deletes `work/chain495.past-gate` on start and in `back()`, and writes its own pid there immediately after `patched:`. The relauncher starts 479b only after it has seen that pid alive and the sentinel equals it, or after that same pid has exited (a 495 failure must not block the 10-02 arms; they only need `builder == HEAD`). A missing pid and a stale log line do not count. The failure path writes `END`, and the cap is 08:00Z 09-29, the same cap §3 gives `chain479b`. The chain495 hunks do not touch the deadline lines, so they still apply after §3. The relauncher hunk replaces §3's hunk for that file.

The replacement predicate, stubbed without `ps` (`repro/relauncher_predicate.sh`):

```
not started (no pid, stale gate ignored) : WAIT
running, gate not written yet            : WAIT
running, past-gate is this pid           : START
seen, then exited before patch (back)    : START
stale gate, a new pid is what we saw     : WAIT
never saw it, gate equals some old pid   : WAIT
```

## 5. One unchanged file aborts the whole commit

`scripts/commit_after_chain.sh:45-47`:

```bash
for p in "${FILES[@]}"; do
  [ "$(sha of $p)" != "$(sha of HEAD:$p)" ] || { say "$p equals HEAD -- nothing to commit"; exit 0; }
done
```

The `exit 0` is inside the loop. The first path that already matches HEAD ends the script before `git commit --only` and before the EXPECT check. Exit status 0, so a supervisor that restarts only on failure leaves it there.

On a first pass where every listed file differs, the loop finishes and the commit runs. `commit476.spec` lists five files and `chain476.sh:49-57` refuses, via `back()`, if the pin anchor is already gone, so that particular re-run does not produce a mixed set. The mixed set is any later state in which one listed path is already HEAD and another is not: a committer re-armed after one of the files was committed on its own, or a two-file spec (`commit497.spec`, `commit496.spec`) whose second path was edited after the first was committed. The script then reports success and leaves the dirty file uncommitted. The next chain waits for that path to equal HEAD until its cap (§3) and stops. A one-file spec (`commit479`, `commit495`) is the case the `exit 0` was written for.

The RESTORED refusal (`commit_after_chain.sh:38`) is before this loop and is not the bug. It still refuses a log that contains `RESTORED`.

CONFIRMED. `bash experiments/501-cloud-automation-review/repro/committer_head.sh` (temp repo, a.txt equals HEAD, b.txt differs):

```
EXIT 0 at a.txt (equals HEAD) — commit not reached
HEAD after the loop: f698eb2 init
worktree b.txt=TWO  (uncommitted if the loop exited early)
status:  M b.txt
```

Fix: `patches/f5-committer-head.patch`. Count differing files, and exit 0 only when every listed file equals HEAD. The EXPECT check and `git commit --only` then run for the files that differ.

## 6. The tie rule demotes a cap block; stercile ties miss LOW

`sharpen_q_by_tercile` is not in the HEAD builder. `chain479` inserts it from `experiments/479-pertarget-sharpen/patch_builder479.py`. `np.digitize` defaults to `right=False`: a value equal to the lower edge goes to the next bin (MID), a value equal to the upper edge goes to HIGH, and duplicate edges `[c, c]` send both `x == c` and `x > c` to HIGH. Flipping `right` moves single boundary points and can change the priced 100/100/100 split, so the fix does not.

The tercile tie loop (`patch_builder479.py` inside the inserted function, the `for e in edges` block) marks every block of ≥2 rows within 1e-9 relative of **either** cut as MID, then the floor guard overwrites rows at `min(r)` back to LOW. The only assignment that loop changes, once the floor has run, is an upper-edge block. That is the cap. 493h says not to infer a tie from `r`, because a cap block (`w` at the cap) exists, and 540's reading is that capped targets sit in HIGH.

CONFIRMED. `python3 experiments/501-cloud-automation-review/repro/tercile.py`:

```
CONTROL distinct ratios          counts [100, 100, 100]  tie_to_mid=0
TERCILE cap of 49 (cut below the cap)   cap rows HIGH, counts [100, 100, 100]  tie_to_mid=0
TERCILE cap of 150 (cut inside the cap) counts [100, 200, 0]  tie_to_mid=150  cap rows all MID
STERCILE 140 at the lower edge   counts [0, 200, 100]   digitize(lower)=1
STERCILE both cuts collapsed     counts [0, 0, 300]     digitize(tied value)=2  (mildest q)
STERCILE 140 exact zeros         counts [0, 200, 100]   LOW empty
TERCILE 10 zero-carrier rows     bins all 1 (median fill)
TERCILE one 1e6 tail             counts [100, 100, 100]
```

A 49-row cap does not contain the 2/3 cut and is already HIGH. A single huge tail does not move the quantiles. The floor collapse (220 rows at `min(r)`, both cuts equal to that value) is already LOW/HIGH 220/80, because the floor guard overwrites the tie. Those three are not the bug.

`stercile` has no floor and no tie rule. A tie of size ≥2 on the lower edge empties LOW: the weakest rows get q 2.0 instead of 2.5. When more than two thirds of the rows are identical, both quantiles collapse and every row gets q 1.75, the mildest, which is the opposite of sharpening the weak rows. S_P2 is the October default candidate, and 493 measured it tie-free; the failure shows up on a panel that is not.

Non-finite ratios (zero carrier norm) are filled with the median and then digitized, so they land in MID and they move the cuts. The A/B/C known answer has `nonfinite == 0` (`test_pt479.py`).

Fix: `patches/f6-tercile.patch`, applied to `patch_builder479.py` (the chain patches the builder from that file).

- Tercile: run the tie-to-MID loop only for an edge at the floor. Upper-edge rows stay at digitize's HIGH. Floor still wins. Control and the 49-cap stay `[100, 100, 100]`, `tie_to_mid` 0. The 150-cap becomes `[100, 50, 150]`.
- Stercile: a lower-edge block of size ≥2 goes LOW (140 zeros → `[140, 60, 100]`). Equal edges send the tied mass LOW and everything strictly above HIGH (`[220, 0, 80]` instead of `[0, 0, 300]`). A single boundary point stays where digitize put it, so the distinct-strength control stays `[100, 100, 100]`.
- Quantiles are computed on finite rows only. Non-finite rows are binned MID and do not enter the cuts. With `nonfinite == 0` the edges match the old ones. The synthetic 10-nan panel moves from `[100, 100, 100]` to `[97, 106, 97]` because those 10 no longer take tercile slots; the unit test forbids that panel.

`test_pt479.py` still requires 100/100/100, `tie_to_mid == 0`, `nonfinite == 0` on 446bj2's LODO dump. That dump is gitignored here, so this review did not re-run it. 493h records that dump as tie-free. If a lower-edge block of size ≥2 is actually present, the stercile assertion fails and this patch must not be applied.

## Checked, no finding

- **-m survives the hop.** `chain1002.sh` / `chain1003.sh` / `chain1004.sh` pass `"$(name $tag)"` as the arm's `$5`, and the arm runs `vcc submit -m "$NAME"` and records `--note "$NOTE: $NAME"`. `repro/names.sh` pulls those strings out of `name()` and execs that hop. Apostrophes, parentheses, the em dash and `!!` each arrive as one argv. Non-interactive bash has `histexpand` off, so `!!` is literal:

```
NAME=taking action is the best action—anytime. dont let it stay in your mind, go grind!!
ARGV <vcc> <submit> <pkg.vcc> <-m> <taking action is the best action—anytime. dont let it stay in your mind, go grind!!> <--wait>
histexpand      off
```

- **s2 is never s1.** `chain1002.sh:47` and `chain1003.sh:46` take the first two distinct ready tags. `chain1004.sh` skips s1's tag inside s2's list and promotes s2 only when s1's list produced nothing. `repro/select.sh`: both ready → `479p2g`/`479sp2`; only fallbacks → `495`/`471c`; one ready → s2 `<none>`; `498sp2` on both of 1004's lists → `498sp2`/`460dr`. No case sets them equal.

- **An already-scored package is not re-submitted.** `ready()` and the arm (`arm1002.sh:18,36`, same in 1003 and 1004) both require `check_submission_novelty.py check` to exit 0. On a ledger hit that function prints `ALREADY BEEN SCORED` and returns 1 (`scripts/check_submission_novelty.py:311-326`), for the `.vcc` digest, the prediction digest, or the content digest. The arm stops before `vcc submit`.

- **Rho band.** Production re-levelled arms use `abs(x / 0.9793 - 1) <= 0.005` on all three contexts (`chain1002.sh:38`, `chain1003.sh:35`, `chain1004.sh:56`). `0.9745 / 0.9793 - 1` is 0.00490, inside the band. `497dr` skips it (`chain1003.sh:30`); 1004 checks only `498p2g` and `498sp2`. The fallbacks 458e/460/461/469 are not on that band. Their gate is the package, the PASS/CLIFF lines and novelty.

- **REQ lines match the logs.** `KA-$T PASS` with `T` in `kas kat kast` (`chain479.sh:40`) matches `KA-kas PASS` and does not match `KA-kast PASS`. `END builder ${NEW:0:16}` matches `END builder d09a986c`. Applying `patch_builder479.py` to a copy of the builder hashes to `d09a986c3eb242c0`, which is `commit479.spec`'s EXPECT. `patched: ${OLD:0:16} -> ${NEW:0:16}`, `KA-1 PASS`, `dump KA PASS`, `(a) PASS` … `(d) PASS`, `END all PASS`, `applied: def_day` and `KA-H1 PASS` are the lines those chains log. `back()` logs `RESTORED` and `END`; the committer refuses a log containing `RESTORED` (`commit_after_chain.sh:38`). A step-5 build failure in `chain495.sh:88-95` logs `END` and leaves the builder patched on purpose.

- **Queue argv.** `repro/gate_matrix.py`: no chain's gate matches its own argv. `commit_after_chain.sh`, `chain1002/1003/1004`, `arm1002/1003/1004` and `relaunch479b.sh` match none of the gates. `run_arm_generic.sh`, `run_ka_4*`, `def_day.sh` and `run_dump4*` match while a build is in progress, so a second chain does not leave its `until` during them. The one cross-match that deadlocks is §4 (`chain495` vs `chain479b`). `chain496p.sh` does not match `chain496\.sh`.

- **Committer parse.** The spec reader splits on the first `=`, appends `REQ` and `FILE` with arrays, and has no `declare -A`. Those are bash 3.2 features. The EXPECT value `path=sha16` survives the first-equals split. These repros ran on this VM's bash 5.2.21; there is no bash 3.2 binary here. Nothing in the spec grammar depends on a 4.x feature. The `index.lock` loop retries the commit up to 20 times. `seq` is the same command `scripts/def_day.sh` already runs on the armed machine.

- **Entry-id regex on `vcc`'s own line.** `arm1002.sh`'s `grep -o 'entry [A-Za-z0-9]\{20\}'` matches `submitted — entry YvVru3DAF89Jiyr34olR` in `submit0927s1.log`. The broken lookup is `eid()` (§2), not that regex.

- **Digitize plus one huge tail, and a 49-row cap,** do not move the tercile counts (§6). The floor guard already keeps a collapsed floor mass in LOW.

## Coordinator's disposition (2026-09-27 07:56Z)
All six accepted; the first run (bc-89f0f1e2) errored before pushing, and this report is the relaunch (bc-5b11face).
1. **Slot-1 noon cap: fixed** (daa69b7d). The 12:00Z check is repeated after the wait in arm0929/0930/1001–1005.
   arm0928 is running and was left alone: caffeinate is on and the Mac is on AC.
2. **eid(): fixed before this report landed** (82bc4e39). `scripts/slot_entry_id.py` takes the capture after the tag's LAST
   ARMED line; chain1004/1005 use it.
3. **Deadline inversions: fixed** (bc6c23ea). Every waiting chain and committer was re-armed with later caps, in dependency
   order (the table above, plus chain510 at 18:00Z 10-02 and commit511).
4. **Relauncher: moot for this run.** It fired at 06:56:52Z, after chain495 logged `patched:` (a fresh log, so no stale
   line), and chain479b is now waiting correctly behind chain495. The sentinel design is noted for any future relauncher.
5. **Committer early exit: fixed before this report landed** (82bc4e39). It exits 0 only when every file equals HEAD, and
   it re-reads the pid file each minute.
6. **Tie rules: queued as 511.** The patch is tolerance-based: the review's exact-tie version split near-ties of ~1e-14
   in our tests. It is unchanged on A/B/C and efr.

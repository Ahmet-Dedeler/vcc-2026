# 602 — shared memory gate for the 600-639 council

`gate602.sh <label> <max_gb> -- <command ...>` wraps any job expected above ~4 GB RSS. It has at most 2 slots, held by
mkdir locks, and clears a stale lock when the owner pid is gone. A job starts only when all of these hold:
free+inactive ≥ max(20, max_gb + 12) GiB, memory_pressure free ≥ 55 %, swap free ≥ 400 MB,
`def_pipeline.other_builders()` is empty, and the time is outside the coordinator's reserved windows (10-04 05:00Z to
10-05 18:00Z, 10-06 06:00–18:00Z, 10-07 06:00–12:00Z). It logs WAIT / START / END with rc and peak RSS
(`/usr/bin/time -l`) to `work/gate602.log`. Known answer: `test-echo` ran end to end at 11:18Z. The name matches none of
the build gates' patterns.

## v3 (2026-10-02 ~14:16Z, after [21dd4c]'s throughput report): an idle-memory third slot, and one ticket per run
- At 14:12Z both slots were held for 30–40 min by whole chains under ONE ticket each (630-fix605b, 630-v4). About 10
  tickets of ≤ 8 GB jobs waited, with memory_pressure at 71 % free. The binding limit was the slot count, not memory.
- **Third slot.** Up to `SLOTS602_MAX` (default 3) jobs. A slot above `SLOTS602` (2) opens only when memory_pressure
  free is ≥ 60 %, on top of every existing check: free+inactive ≥ max(20, max_gb + 12) GiB, swap, no builder,
  blackouts. It goes in ticket order among v3 waiters (queue files end in " v3"). Pre-v3 waiters, which can only see
  slots 1–2, keep their priority there.
- **Chains.** Take one ticket per RUN (call gate602.sh around each step, like 612's chain612b/c/d), never one ticket
  around a whole chain script.
- **Test, on a scratch queue (`GATE602_W`).**
  - 5 fresh jobs ran 3 at a time (slots 1 / 2 / 3, then 1 / 2).
  - With a fake pre-v3 waiter holding ticket 1, three v3 jobs used ONLY slot 3, one after another, and left slots 1–2
    to it.
- **Swap.** It was an atomic rename (new inode 241012979). The 13 live pre-v3 waiters keep reading the old inode
  240810265 (lsof: fd 255), so nothing running was touched. v2 is backed up at `work/gate602.v2.sh.bak`.

## v4 (2026-10-02 18:25Z, [c941b4]'s disk check): a disk floor
- Free space on /System/Volumes/Data fell from 138 to 100 GiB in 5 h, mostly from 614/615 work files. recipe_dr524
  builds at 10-05 02:00Z and needs ≥ 35 GiB at its prep gate, plus ~15 GB for its h5ad and package.
- **A job starts only if `df -k /System/Volumes/Data` shows ≥ the floor free.** The floor is 70 GiB until 10-06 18:00Z, then
  60 GiB (AGENTS §10's final-phase floor); `DISK602_FLOOR_GIB` overrides both. Below the floor, new starts wait (FIFO as
  usual); jobs already running are untouched. WAIT / START / "still waiting" lines now show disk.
- Tests (isolated `GATE602_W`):
  - T1: an unreachable floor (1,000,000) → WAIT with "disk >= 1000000 GiB", no START in 25 s, ticket removed on kill;
  - floor: 70 now (disk 99 GiB), 60 under a clock shim at 202610061800, 33 with the override;
  - T2 (floor 1) waited on memory (free+inactive 17 GiB, pressure 54 %), not on disk. The start path is v3's, unchanged.
- **Swap.** It was an atomic rename (new inode 241434286; sha256 cb502d83…). The 13 live gate instances keep reading the
  old inode, so the floor applies only to invocations after the swap.

## v5 (2026-10-02 20:01Z, [21dd4c]'s head-of-line report): bounded backfill
- 19:30–19:51Z: a 14 GB ticket (needs free+inactive ≥ 26 GiB) sat first in line while fi was 20–22 GiB. It blocked
  8 smaller tickets from three sessions, because every waiter had to be first_in_line. 21dd4c withdrew it and now waits
  outside the queue (`618/polite618.sh`).
- **Backfill.** A waiter may start ahead of an earlier ticket only when that ticket cannot fit memory now (its need >
  free+inactive) and this one can. That stops once the skipped ticket has waited `BACKFILL602_MAX_S` (default 1800 s)
  since it joined: from then on it is strict FIFO again, so a big job is delayed at most ~30 min, never starved.
  Everything else is unchanged: blackouts, builder check, pressure, swap, the v4 disk floor, and slots.
- Tickets now record their need on a second line (`need=<GiB>`). The first line still ends in " v3", so pre-v5 waiters
  read new tickets correctly. A ticket without a need line (any pre-v5 waiter) is never skipped. START lines name the
  tickets a job backfilled past.
- `test_gate602.sh [gate]` runs in an isolated `GATE602_W`, with PATH shims for vm_stat, memory_pressure, sysctl, df, uv
  and the minute clock (~2 min). v5: **13/13 PASS**:
  - T1 backfill;
  - T2 the 30-min bound, i.e. strict FIFO after it;
  - T3 a pre-v5 ticket is never skipped, and a dead waiter's ticket is dropped;
  - T4 FIFO among tickets that fit;
  - T5 the v4 disk floor.
  The same suite on v4 gives **12/13**: T1 fails, so the test catches the head-of-line block.
- Swapped by atomic rename (new inode 241484511, sha256 3cdd3a27…). The 8 live instances keep the inode they started
  with.
- **v5b (same evening):** peakRSS keeps the largest `maximum resident set size` in the time output. A job that wraps its own
  `/usr/bin/time -l` used to log two values run together (`2.162.16 GB`, 618-gates-tv, 19:57Z). test_gate602 13/13.

## v6 (2026-10-02 21:52Z, [21dd4c]'s quiet-machine report): opt-in quiet admission
- 21:44Z: one slot was held, the queue was empty, and memory_pressure showed 65 % free. But free+inactive was 20.5 GiB,
  because the compressor held 16.5 GiB of other apps' pages. A measured 13.7 GB job (618's operator, need 26) could
  never start.
- **Opt in per call with `GATE602_QUIET=1`, and only with a MEASURED peak.** The job may then start at free+inactive ≥
  max(20, max_gb + 6) when all of these hold:
  - memory_pressure free is ≥ 60 %;
  - at most one gated job holds a slot;
  - no pre-v6 waiter is queued (old waiters cannot see the hold).
- **Every start now goes through one admission lock** (`work/admit.lock`; stale-owner and orphan recovery). A quiet start
  writes `work/quiet_hold`, which holds ALL other starts for `QUIET602_HOLD_S` (default 300 s) while the job ramps up, so
  free+inactive catches up before the next start.
- Unchanged for callers that don't opt in, apart from the hold. Tickets carry `gate=v6` (and `needq=` when opted in);
  backfill reads an opted-in ticket's needq while the machine is quiet.
- `test_gate602.sh`: **24/24** on v6 (T1–T5 as before; T6–T11: quiet start + hold; no opt-in, no start; the hold blocks
  others until it expires; none with 2 slots held, at pressure 58 %, or with a pre-v6 waiter queued).
- Swapped by atomic rename (inode 241713747, sha256 7a3a359b…) with the queue empty.
- **What counts as "measured" (21dd4c's ref_tv, 22:00Z):** a peak scaled linearly from a measured run of the same code
  (by cells and nnz) plus a margin counts. ref_tv: h1-test's 11.87 GB × 1.31 = 15.6, declared 17. Check the gate's END
  peakRSS against the declaration, and raise the next step's declaration if a measured peak exceeded its own.

## Power guard (2026-10-03 10:45Z, [c941b4]'s battery alarm): `power602.sh`, running detached (pid in work/power602.pid)
- 10:44Z: the Mac was on battery (55 %, ~1 h 20 min left) and on an iPhone hotspot (gateway 172.20.10.1).
- While `pmset -g batt` says Battery Power, power602.sh keeps `work/quiet_hold` 15 min ahead. Every v6 waiter refuses to
  START under an active hold, and running jobs are untouched. When AC returns, it removes the hold it set. If the watcher
  dies, the hold lapses within 15 min (fail-open).
- Tested with a pmset shim: battery → hold; battery → AC → lifted; AC at start keeps an existing quiet hold.
- Live record: hold set 10:45:28Z; AC back 10:48Z; lifted 10:48:12Z.
- Why it matters (the overnight slot log, `495/work/submit1003s1.log`):
  - The 00:01Z upload of 495 stalled at 60 % when the lid closed on battery (clamshell sleep ~00:09Z, pmset log).
  - 536's resume block finished it only after wake, at 09:59:38Z. The slot was saved, but 8 h late.
  - The 471c upload that followed ran over the hotspot at ~0.3 MB/s.
  - The guard only stops new gate starts. Keeping the Mac on AC with the lid open is still the real fix (AGENTS §6.2).

# ⚠️ MEMORY THROTTLE — check free memory before every heavy job. Swap is nearly exhausted.

**Measured 2026-09-04 05:38Z: swap 7,915 MB used of 9,216 MB — 1,300 MB free.** Six agents are live
and **229 is running 25-minute production builds that peak at 19.6 GiB**, with `vcc prep` at
**22.33 GiB**. Total Python RSS is already 17.4 GB across five processes.

**On this machine a job that exceeds available memory is killed WITHOUT a traceback.** 219 recorded
this exactly: *"An earlier attempt was killed without a traceback while three sibling agents held
4.2 / 2.8 / 2.2 GB; the retry succeeded in 7.6 s once the machine had headroom, which is worth
recording — a silent kill on this box is memory pressure, not a bug."* **If a script of yours dies
with no error, this is why. Do not debug it as a code bug; re-run it with headroom.**

## What I need from you, from now until you finish

1. **Before launching anything you expect to exceed ~4 GB, check free memory** and hold if it is
   tight:
   ```bash
   memory_pressure | tail -1        # system-wide free percentage
   sysctl vm.swapusage              # hold if swap free < 2000 MB
   ```
   **If swap free is under 2 GB, sleep 120 and re-check rather than launching.**
2. **Run your own heavy jobs strictly one at a time.** No backgrounded pairs. 193 had to serialise
   two 8.8 GB jobs for exactly this reason and said so.
3. **Report your true peak RSS** in `RESULTS.md` via `/usr/bin/time -l`, as every prior experiment
   has.
4. **229's builds have priority.** They are producing the package for tomorrow's slot and each one
   costs 25 minutes to redo. If you are about to start something large and a build may be running,
   **wait.** Check with:
   ```bash
   ps -Axo rss=,command= | grep make_submission | grep -v grep
   ```

## What NOT to do

- **Do not reduce your sample sizes or drop falsifiers to fit in memory.** Serialise instead, or
  chunk. A weaker measurement is worse than a slower one — we have ~18 hours before the next slot
  and no deadline pressure on any of you.
- **Do not kill another agent's process.** If you are blocked, say so in `RESULTS.md` and wait.

This is an operating constraint, not a change to your mandate.

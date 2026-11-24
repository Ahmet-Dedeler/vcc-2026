"""300 -- the run scheduler.  The arm ROUNDS are frozen in PREREGISTRATION.md section 7; this walks them in order,
levelling the three cancer folds after every round, with AT MOST TWO fold processes in flight and never two
processes on the same fold (they would race on the same sweep json).  Every arm is flushed by `run300.py` as soon
as it is measured, so stopping at the deadline loses only the arms in flight.

    PYTHONDONTWRITEBYTECODE=1 uv run python sched300.py <deadline HH:MM UTC>
"""
from __future__ import annotations

import datetime as dt
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps")
ROUNDS = [
    ["cd4svd5", "STACK", "tilt:d1.0"],      # KA-3 / KA-4 / R2 / KA-W1 + the primary arm
    ["tilt:d-1.0"],                          # KA-W2 + the primary arm's mirror control
    ["tilt:d2.5"],                           # the A ladder (P1)
    ["tilt:d1.5"],
    ["tilt:d0.5"],
    ["tilt:d2.0"],
    ["prune:d0.0"],                          # B: the prune alone (KA-B)
    ["spread:d1.0"],                         # C (P3)
    ["prune:d1.0"],                          # B at the primary dose (P2)
    ["spread:d2.0"],
    ["rand:d1.0:s300001"],                   # the random-weight control
    ["tilt:d1.8"],                           # KA-W3
    ["spread:d-1.0"],                        # C's mirror
    ["rand:d1.0:s300002"],
    ["prune:d2.5"],
    ["tilt:d-2.5"],
    ["tilt:d0.0"],                           # KA-2a / KA-2e / KA-2f
    ["tilt:d-1.5"],
    ["prune:d1.5"],
    ["spread:d-2.0"],
]
MAXPROC = 2


def free_pct() -> float:
    try:
        out = subprocess.run(["memory_pressure"], capture_output=True, text=True, check=True).stdout
        for ln in out.splitlines():
            if "free percentage" in ln:
                return float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    except Exception:
        pass
    return 100.0


def main(deadline: str) -> None:
    hh, mm = deadline.split(":")
    now = dt.datetime.now(dt.timezone.utc)
    dl = now.replace(hour=int(hh), minute=int(mm), second=0, microsecond=0)
    if dl < now:
        dl += dt.timedelta(days=1)
    log = WORK / "sched300.log"

    def W(s: str) -> None:
        line = f"{dt.datetime.now(dt.timezone.utc):%H:%M:%S}Z {s}"
        print(line, flush=True)
        with log.open("a") as fh:
            fh.write(line + "\n")

    W(f"scheduler up, deadline {dl:%Y-%m-%dT%H:%M:%S}Z, {len(ROUNDS)} rounds x {len(FOLDS)} folds, max {MAXPROC} in flight")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    # one global queue in the frozen round order; a slot takes the FIRST job whose fold is not already running, so
    # the folds stay levelled without the idle half-wave a strict round barrier costs
    queue = [(ri, f, arms) for ri, arms in enumerate(ROUNDS) for f in FOLDS]
    running: dict[str, tuple] = {}
    while queue or running:
        for f in list(running):
            p, t0, fh, ri = running[f]
            if p.poll() is not None:
                fh.close()
                W(f"  round {ri} {f} exit {p.returncode} after {(time.time()-t0)/60:.1f} min")
                del running[f]
        if dt.datetime.now(dt.timezone.utc) >= dl:
            if not running:
                W("DEADLINE reached, nothing in flight; stopping")
                return
            time.sleep(10)
            continue
        while len(running) < MAXPROC:
            nxt = next((j for j in queue if j[1] not in running), None)
            if nxt is None:
                break
            fp = free_pct()
            if fp < 40.0:
                W(f"  memory gate: {fp:.0f} % free < 40 %, holding")
                break
            queue.remove(nxt)
            ri, f, arms = nxt
            fh = (WORK / f"run300-{f}.nohup").open("a")
            cmd = ["uv", "run", "python", str(HERE / "run300.py"), f, *arms]
            pr = subprocess.Popen(cmd, cwd=str(HERE), env=env, stdout=fh, stderr=subprocess.STDOUT)
            running[f] = (pr, time.time(), fh, ri)
            W(f"round {ri} [{'|'.join(arms)}] -> {f} pid {pr.pid} ({fp:.0f} % free, {len(queue)} queued)")
        time.sleep(10)
    W("all rounds complete")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "23:00")

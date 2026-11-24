"""534 Part C -- the driver: every heavy step strictly serial behind gate534c (MANDATE-C's gate and windows), resumable.

  abc, real, efp:  full534c <panel> (energies + levels, production known answers first)
                   -> comp534c <panel> <fold> --level (one process per fold)
  efr (fold route): comp534c efr <fold> (unlevelled first pass, x3) -> foldlevel534c efr (light: fold-route factors,
                   KA vs 438's efr factors) -> comp534c efr <fold> --level (second pass, KA-RESUME)
Order: abc, real, efp, efr. State: work/partC/state534c.json. A closed window writes work/partC/RESUME.md and exits 3.

    uv run python drive534c.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import designs534c as D  # noqa: E402
import gate534c as G  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REL = HERE.relative_to(ROOT)
STATE = HERE / "state534c.json"
PY = ["uv", "run", "python"]


def state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def done(step, rec):
    s = state()
    s[step] = dict(rec, done=G.utc())
    STATE.write_text(json.dumps(s, indent=1))


def say(*a):
    print(f"[drive534c {G.utc()}]", *a, flush=True)


def heavy(step, cmd, log):
    if step in state():
        return say(f"skip {step} (done)")
    say(f"{step}: waiting for the gate")
    r = G.heavy(cmd, step, log)
    say(f"{step}: rc {r['rc']}, {r['wall_s']:.0f} s, peak RSS {r['peak_rss_gib']:.2f} GiB")
    if r["rc"]:
        raise SystemExit(f"{step} failed (rc {r['rc']}); see {log}")
    done(step, r)


def light(step, cmd, log):
    if step in state():
        return say(f"skip {step} (done)")
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    say(f"{step}: rc {rc} ({time.time() - t0:.0f} s)")
    if rc:
        raise SystemExit(f"{step} failed (rc {rc}); see {log}")
    done(step, dict(rc=rc, seconds=time.time() - t0, log=str(log.relative_to(ROOT))))


def main():
    pk = "--pk" in sys.argv                     # secondary: the same runs plus B1pk..B3pk (toggled base, production key)
    x, fl = ("pk:", ["--pk"]) if pk else ("", [])
    s = "pk" if pk else ""
    try:
        for p in D.PANELS:
            if p in D.FOLD_ROUTE:
                for f in D.CANCER:
                    heavy(f"{x}{p}:comp534c:{f}:pass1", [*PY, str(REL / "comp534c.py"), p, f, *fl],
                          HERE / f"logs/comp534c{s}-{p}-{f}-pass1.log")
                light(f"{x}{p}:foldlevel", [*PY, str(REL / "foldlevel534c.py"), p, *fl], HERE / f"logs/foldlevel534c{s}-{p}.log")
            else:
                heavy(f"{x}{p}:full534c", [*PY, str(REL / "full534c.py"), p, *fl], HERE / f"logs/full534c{s}-{p}.log")
            lv = str(REL / p / f"level534c{s}.json")
            for f in D.CANCER:
                heavy(f"{x}{p}:comp534c:{f}", [*PY, str(REL / "comp534c.py"), p, f, "--level", lv, *fl],
                      HERE / f"logs/comp534c{s}-{p}-{f}.log")
        say("all panels priced" + (" (pk)" if pk else ""))
    except G.WindowClosed as e:
        s = state()
        (HERE / "RESUME.md").write_text("\n".join(
            ["# 534 Part C -- RESUME", "", f"Written {G.utc()} by drive534c.py: a MANDATE-C window closed ({e}).", "",
             "Done steps (work/partC/state534c.json):", ""] + [f"- {k} ({v.get('done')})" for k, v in s.items()]
            + ["", "Next: after the window, re-run `uv run python experiments/534-generalization-audit/work/partC/drive534c.py`"
               " (done steps are skipped), then `aggregate534c.py`, then write RESULTS-C.md.", ""]))
        say(f"WINDOW CLOSED ({e}); RESUME.md written")
        raise SystemExit(3)


if __name__ == "__main__":
    main()

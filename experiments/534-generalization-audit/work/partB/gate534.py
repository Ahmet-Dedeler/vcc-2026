"""534 Part B -- the memory/builder/chain gate (PREREGISTRATION-B section 9, MANDATE-B "Rules") and the heavy runner.

A heavy process (an operator fit or a scoring run) starts only when ALL hold:
  * `memory_pressure` reports >= 40 % free;
  * `vm_stat` free + inactive >= 33 GiB;
  * `def_pipeline.other_builders()` is empty (run exactly as the mandate words it, in a fresh interpreter);
  * no bash/sh/zsh process is RUNNING (its script argument, not its mandate text) run_arm_generic.sh, run_ka_*,
    def_day.sh, chain4NN* / chain5NN*.sh, relevel_sp2.sh, run_dump4*, build527.sh or fix526.sh (AGENTS 6.2: match the
    interpreter plus the script path; a worker whose prompt MENTIONS a script is not running it);
  * the coordinator's window: no heavy START after 2026-09-29 23:45Z or during 2026-09-30 00:30-06:00Z.
Closed -> poll every 60 s; each wait is capped at 2 h, after which it is logged and re-armed (re-checked). A window
closure is terminal: the caller writes RESUME.md.

`heavy(cmd, tag)` runs cmd under /usr/bin/time -l (cwd = repo root), strictly serially, appending its start and end
(UTC, command, rc, wall, peak RSS) to work/partB/heavy.log.
"""
from __future__ import annotations

import datetime as dt
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
HEAVY_LOG = HERE / "heavy.log"
GATE_LOG = HERE / "gate.log"
MIN_FREE_PCT, MIN_FREE_GIB = 40.0, 33.0
CHAIN = re.compile(r"(^|/)(run_arm_generic\.sh|run_ka_[^/]*|def_day\.sh|chain[45][0-9][0-9][^/]*\.sh|relevel_sp2\.sh|"
                   r"run_dump4[^/]*|build527\.sh|fix526\.sh)$")
SHELLS = {"bash", "sh", "zsh"}
NO_START_AFTER = dt.datetime(2026, 9, 29, 23, 45, tzinfo=dt.timezone.utc)
BLACKOUT = (dt.datetime(2026, 9, 30, 0, 30, tzinfo=dt.timezone.utc), dt.datetime(2026, 9, 30, 6, 0, tzinfo=dt.timezone.utc))


class WindowClosed(Exception):
    pass


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def window_ok(now=None) -> tuple[bool, str]:
    now = now or dt.datetime.now(dt.timezone.utc)
    if BLACKOUT[0] <= now < BLACKOUT[1]:
        return False, "coordinator blackout 2026-09-30 00:30-06:00Z"
    if NO_START_AFTER <= now < BLACKOUT[0]:
        return False, "no heavy start after 2026-09-29 23:45Z"
    return True, "window open"


def memory() -> tuple[float, float]:
    mp = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    pct = 0.0
    for ln in mp.splitlines():
        if "free percentage" in ln:
            pct = float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = 16384
    m = re.search(r"page size of (\d+) bytes", vm)
    if m:
        page = int(m.group(1))

    def pages(tag):
        for ln in vm.splitlines():
            if ln.startswith(tag):
                return int(ln.split(":")[1].strip().rstrip("."))
        return 0
    return pct, (pages("Pages free") + pages("Pages inactive")) * page / 2**30


def builders() -> list[str]:
    r = subprocess.run(["uv", "run", "python", "-c",
                        "import sys; sys.path.insert(0,'scripts'); import def_pipeline as d; print(repr(d.other_builders()))"],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        return [f"other_builders() failed rc {r.returncode}: {r.stderr.strip()[-200:]}"]   # fail closed
    import ast
    return list(ast.literal_eval(r.stdout.strip().splitlines()[-1]))


def chains() -> list[str]:
    out = subprocess.run(["ps", "-axo", "pid=,comm=,args="], capture_output=True, text=True).stdout
    hits = []
    for ln in out.splitlines():
        parts = ln.split(None, 2)
        if len(parts) < 3:
            continue
        pid, comm, args = parts
        if os.path.basename(comm) not in SHELLS:
            continue
        argv = args.split()
        # the script a shell RUNS: its first non-option argument (bash [-opts] script ...); `-c CMD`: the command's
        # first word. A prompt text that merely mentions a name is never argv[1] of a shell running it.
        script = None
        i = 1
        while i < len(argv):
            a = argv[i]
            if a == "-c":
                script = argv[i + 1] if i + 1 < len(argv) else None
                break
            if a.startswith("-"):
                i += 1
                continue
            script = a
            break
        if script and CHAIN.search(script):
            hits.append(f"{pid} {args[:160]}")
    return hits


OVERRIDE = HERE / "GATE-OVERRIDE.json"


def floor_gib() -> tuple[float, str]:
    """The free+inactive floor: the preregistered 33 GiB, unless the OWNER of 534 has written GATE-OVERRIDE.json
    ({"min_free_gib": x, "by": ..., "reason": ...}; 10 <= x <= 33). This worker never writes that file: it is the owner's
    amendment of PREREGISTRATION-B section 9, recorded as a deviation in RESULTS-B.md if it is ever used."""
    if OVERRIDE.exists():
        import json
        o = json.loads(OVERRIDE.read_text())
        x = float(o["min_free_gib"])
        assert 10.0 <= x <= MIN_FREE_GIB and o.get("by") and o.get("reason"), f"malformed {OVERRIDE}: {o}"
        return x, f"owner override {x:g} GiB by {o['by']}"
    return MIN_FREE_GIB, "preregistered 33 GiB"


def gate() -> tuple[bool, str]:
    ok_w, wmsg = window_ok()
    pct, fi = memory()
    b = builders()
    c = chains()
    floor, fsrc = floor_gib()
    ok = pct >= MIN_FREE_PCT and fi >= floor and not b and not c
    msg = f"memory_pressure free {pct:.0f}%, free+inactive {fi:.1f} GiB (floor {floor:g}: {fsrc}), builders {len(b)}, chains {len(c)}"
    if b:
        msg += " | builders: " + " ;; ".join(x[:140] for x in b)
    if c:
        msg += " | chains: " + " ;; ".join(c)
    return ok and ok_w, msg + f" | {wmsg}"


def wait_gate(tag: str) -> str:
    t0 = time.time()
    n = 0
    while True:
        ok_w, wmsg = window_ok()
        if not ok_w:
            with open(GATE_LOG, "a") as fh:
                fh.write(f"{utc()} {tag}: WINDOW CLOSED ({wmsg})\n")
            raise WindowClosed(wmsg)
        ok, msg = gate()
        if ok:
            with open(GATE_LOG, "a") as fh:
                fh.write(f"{utc()} {tag}: OPEN after {time.time() - t0:.0f}s ({msg})\n")
            return msg
        if n % 10 == 0:
            with open(GATE_LOG, "a") as fh:
                fh.write(f"{utc()} {tag}: closed ({msg})\n")
        n += 1
        if time.time() - t0 > 7200:                    # a 2 h wait: log it and re-arm (re-check from scratch)
            with open(GATE_LOG, "a") as fh:
                fh.write(f"{utc()} {tag}: gate closed for 2 h; re-arming the wait\n")
            t0 = time.time()
        time.sleep(60)


def heavy(cmd: list[str], tag: str, log: Path) -> dict:
    """Wait for the gate, run cmd under /usr/bin/time -l (cwd ROOT), log start/end/peak RSS to heavy.log."""
    gmsg = wait_gate(tag)
    start = utc()
    with open(HEAVY_LOG, "a") as fh:
        fh.write(f"{start} START {tag} :: {' '.join(cmd)} :: gate {gmsg}\n")
    t0 = time.time()
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        rc = subprocess.run(["/usr/bin/time", "-l", *cmd], cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT,
                            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}).returncode
    wall = time.time() - t0
    txt = log.read_text()
    rss = [ln for ln in txt.splitlines() if "maximum resident set size" in ln]
    peak = int(rss[-1].split()[0]) / 2**30 if rss else float("nan")
    with open(HEAVY_LOG, "a") as fh:
        fh.write(f"{utc()} END   {tag} :: rc {rc} :: wall {wall:.0f}s :: peak RSS {peak:.2f} GiB (ru_maxrss) :: log {log.relative_to(ROOT)}\n")
    return dict(tag=tag, rc=rc, wall_s=wall, peak_rss_gib=peak, start=start, log=str(log.relative_to(ROOT)))


if __name__ == "__main__":
    ok, msg = gate()
    print(("OPEN " if ok else "CLOSED ") + msg)
    sys.exit(0 if ok else 1)

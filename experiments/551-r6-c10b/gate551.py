"""551 (C10b) -- gate534.py with 551's rules (PREREGISTRATION.md section 9, agreed with the coordinator 2026-09-30):
floor 27 GiB free+inactive (534's owner override, 33 GiB is unreachable with the desktop apps open), memory_pressure
>= 40 %, no builder, no build script running (the current family: run_arm_generic.sh, def_day.sh, run_ka_*, build5NN.sh,
rebuild5NN.sh, relevel*.sh, run_dump*, fix5NN.sh, chain4NN/5NN*.sh builds; NOT the long-lived slot
waiters chain1NNN.sh / arm1NNN.sh, nor recipe_dr524.sh while it waits: once it builds, its def_day.sh
child and the python builders close the gate), data volume >= 25 GiB free, and no heavy START during 2026-10-05 00:00-18:00Z (the
coordinator's crowded day) or after 2026-10-06 12:00Z (the ledger's limit). Amended 2026-09-30 ~21:55Z after the
coordinator's message: its windows are 10-04 12:00Z - 10-05 18:00Z, and C10b's access starts heavy steps only in
10-04 05:00-16:00Z (PHASE 'access'). gate534's docstring follows.

534 Part B -- the memory/builder/chain gate (PREREGISTRATION-B section 9, MANDATE-B "Rules") and the heavy runner.

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
ROOT = HERE.parents[1]
HEAVY_LOG = HERE / "work" / "heavy.log"
GATE_LOG = HERE / "work" / "gate.log"
MIN_FREE_PCT, MIN_FREE_GIB, MIN_DISK_GIB = 40.0, 27.0, 25.0
CHAIN = re.compile(r"(^|/)(run_arm_generic\.sh|run_ka_[^/]*|def_day\.sh|chain[45][0-9][0-9][^/]*\.sh|relevel[^/]*\.sh|"
                   r"run_dump[^/]*|build5[0-9][0-9][^/]*\.sh|rebuild5[0-9][0-9][^/]*\.sh|fix5[0-9][0-9][^/]*\.sh|"
                   r"prep5[0-9][0-9][^/]*\.sh)$")
# The coordinator's approval (message 2026-09-30 ~22:30Z): <= 8 GB jobs from 10-01 06:00Z at free+inactive >= 20 GiB AND
# memory_pressure >= 55 % AND no builder AND nothing of its building; C10b's 13 GiB operator fits on 10-04 at >= 22 GiB AND
# memory_pressure >= 55 % AND swap free >= 400 MiB. (free+inactive does not see the compressor, which holds ~12 GB of
# desktop-app pages on this Mac, while memory_pressure reads ~60 % free: the pressure floor is the real guard.)
RELAX_FROM = dt.datetime(2026, 10, 1, 6, 0, tzinfo=dt.timezone.utc)
SHELLS = {"bash", "sh", "zsh"}
NO_START_AFTER = dt.datetime(2026, 10, 6, 12, 0, tzinfo=dt.timezone.utc)
# the coordinator's build windows (message 2026-09-30 ~21:50Z): 10-04 12:00Z - 10-05 18:00Z are its own
BLACKOUT = (dt.datetime(2026, 10, 4, 12, 0, tzinfo=dt.timezone.utc), dt.datetime(2026, 10, 5, 18, 0, tzinfo=dt.timezone.utc))
# C10b's agreed window: 10-04 05:00Z - 12:00Z, slipping with the coordinator's read but ending by 16:00Z at the latest
ACCESS = (dt.datetime(2026, 10, 4, 5, 0, tzinfo=dt.timezone.utc), dt.datetime(2026, 10, 4, 16, 0, tzinfo=dt.timezone.utc))
PHASE = {"name": None}                          # drive551 sets "access" for the C10b run


class WindowClosed(Exception):
    pass


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def window_ok(now=None) -> tuple[bool, str]:
    now = now or dt.datetime.now(dt.timezone.utc)
    if PHASE["name"] == "access":
        if not (ACCESS[0] <= now < ACCESS[1]):
            return False, "C10b's window is 2026-10-04 05:00-16:00Z (the coordinator's; 05:00-12:00Z unless its read slips)"
        return True, "C10b window open"
    if BLACKOUT[0] <= now < BLACKOUT[1]:
        return False, "coordinator build window 2026-10-04 12:00Z - 10-05 18:00Z"
    if now >= NO_START_AFTER:
        return False, "no heavy start after 2026-10-06 12:00Z (ledger limit)"
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


def thresholds(now=None) -> tuple[float, float, float, str]:
    """(free+inactive floor GiB, memory_pressure floor %, swap-free floor MiB, source) for the current phase and time."""
    now = now or dt.datetime.now(dt.timezone.utc)
    if PHASE["name"] == "access":
        return 22.0, 55.0, 400.0, "C10b access: 22 GiB + pressure 55 % + swap 400 MiB (coordinator, 2026-09-30)"
    if now >= RELAX_FROM:
        return 20.0, 55.0, 0.0, "<= 8 GB job: 20 GiB + pressure 55 % (coordinator, from 10-01 06:00Z)"
    return MIN_FREE_GIB, MIN_FREE_PCT, 0.0, "551's 27 GiB + pressure 40 % (before 10-01 06:00Z)"


def swap_free_mib() -> float:
    out = subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout
    if re.search(r"total = 0\.00M", out):
        return float("inf")                          # 670 (10-16): no swap file yet = nothing swapped, not exhausted
    m = re.search(r"free = ([0-9.]+)([MG])", out)
    if not m:
        return 0.0                                   # fail closed
    return float(m.group(1)) * (1024.0 if m.group(2) == "G" else 1.0)


def vcc_preps() -> int:
    out = subprocess.run(["ps", "-axo", "args="], capture_output=True, text=True).stdout
    return sum(1 for ln in out.splitlines() if re.search(r"(^|/| )vcc prep( |$)", ln) and "claude" not in ln)


def floor_gib() -> tuple[float, str]:
    """The free+inactive floor: 551's 27 GiB, unless the coordinator asks for another value (GATE-OVERRIDE.json
    ({"min_free_gib": x, "by": ..., "reason": ...}; 10 <= x <= 33). This worker never writes that file: it is the owner's
    amendment of PREREGISTRATION-B section 9, recorded as a deviation in RESULTS-B.md if it is ever used."""
    if OVERRIDE.exists():
        import json
        o = json.loads(OVERRIDE.read_text())
        x = float(o["min_free_gib"])
        assert 10.0 <= x <= MIN_FREE_GIB and o.get("by") and o.get("reason"), f"malformed {OVERRIDE}: {o}"
        return x, f"owner override {x:g} GiB by {o['by']}"
    f, _, _, src = thresholds()
    return f, src


def gate() -> tuple[bool, str]:
    ok_w, wmsg = window_ok()
    pct, fi = memory()
    b = builders()
    c = chains()
    floor, fsrc = floor_gib()
    _, min_pct, min_swap, _ = thresholds()
    import shutil
    disk = shutil.disk_usage(ROOT).free / 2**30
    sw = swap_free_mib()
    vp = vcc_preps()
    ok = (pct >= min_pct and fi >= floor and not b and not c and disk >= MIN_DISK_GIB and sw >= min_swap and vp == 0)
    msg = (f"memory_pressure free {pct:.0f}% (floor {min_pct:g}), free+inactive {fi:.1f} GiB (floor {floor:g}: {fsrc}), "
           f"swap free {sw:.0f} MiB (floor {min_swap:g}), builders {len(b)}, chains {len(c)}, vcc preps {vp}, disk {disk:.0f} GiB")
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

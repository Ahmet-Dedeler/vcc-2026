"""534 Part B -- the driver: every step in the mandated order, heavy steps strictly serial behind gate534, resumable.

  ka        efp known answers (preregistration section 5), before any draw:
              full534 efp --ka  ->  comp534 efp <3 folds, one process> --level --ka  ->  ka_price534 (light)
  draw      draw534 (light): R6's four panels + bundles; then each panel's spec.json
  panel P   (P in S1, S2, S3, S4, in that order)
              cd4_534 (light) -> the four operator fits (heavy; build_union_whitener.py, 238's flags, --max-rss-gib 14):
              held-out K562 / HCT116 / HEK293T (the other two donors) and the full 3-donor operator -> axes534.check_ops
              -> full534 P (heavy) -> comp534 P <3 folds, one process> --level (heavy) -> operators deleted (sha256 +
              command kept)
  abc       supplementary: the exact designs on A/B/C's own panel (after R6)
  all       ka, draw, S1, S2, S3, S4, abc   (r6: S1..S4 only)

State: work/partB/state534.json (a step is skipped once recorded done). A closed window (no heavy start after
2026-09-29 23:45Z or in 2026-09-30 00:30-06:00Z) writes work/partB/RESUME.md and exits 3.

    uv run python drive534.py all
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import gate534 as G  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REL = HERE.relative_to(ROOT)
STATE = HERE / "state534.json"
PANELS = ("S1", "S2", "S3", "S4")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FITS = {"held_k562_gwps": ("orion_hct116", "orion_hek293t"), "held_orion_hct116": ("k562_gwps", "orion_hek293t"),
        "held_orion_hek293t": ("k562_gwps", "orion_hct116"), "full3": ("k562_gwps", "orion_hct116", "orion_hek293t")}
PY = ["uv", "run", "python"]


def state() -> dict:
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def done(step: str, rec: dict) -> None:
    s = state()
    s[step] = dict(rec, done=G.utc())
    STATE.write_text(json.dumps(s, indent=1))


def say(*a):
    print(f"[drive534 {G.utc()}]", *a, flush=True)


def light(step: str, cmd: list[str], log: Path) -> None:
    if step in state():
        return say(f"skip {step} (done)")
    t0 = time.time()
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
    say(f"{step}: rc {rc} ({time.time() - t0:.0f} s) -> {log.relative_to(ROOT)}")
    if rc:
        raise SystemExit(f"{step} failed (rc {rc}); see {log}")
    done(step, dict(rc=rc, seconds=time.time() - t0, log=str(log.relative_to(ROOT))))


def heavy(step: str, cmd: list[str], log: Path) -> dict:
    if step in state():
        say(f"skip {step} (done)")
        return state()[step]
    say(f"{step}: waiting for the gate")
    r = G.heavy(cmd, step, log)
    say(f"{step}: rc {r['rc']}, {r['wall_s']:.0f} s, peak RSS {r['peak_rss_gib']:.2f} GiB")
    if r["rc"]:
        raise SystemExit(f"{step} failed (rc {r['rc']}); see {log}")
    done(step, r)
    return r


def disk_free_gib() -> float:
    return shutil.disk_usage(ROOT).free / 2**30


def ka():
    lv = str(REL / "efp/level534.json")
    heavy("ka:full534:efp", [*PY, str(REL / "full534.py"), "efp", "--ka"], HERE / "logs/full534-efp.log")
    heavy("ka:comp534:efp", [*PY, str(REL / "comp534.py"), "efp", ",".join(CANCER), "--level", lv, "--ka"],
          HERE / "logs/comp534-efp.log")
    light("ka:price", [*PY, str(REL / "ka_price534.py")], HERE / "logs/ka_price534.log")


def abc():
    """Supplementary (not preregistered): the exact milestone designs on A/B/C's own panel (438's abc inputs, rules abc,
    187's claims), so retention can also be read against an exact-design A/B/C value. Its KA: 438's abc rows and 493's."""
    lv = str(REL / "abc/level534.json")
    heavy("abc:full534", [*PY, str(REL / "full534.py"), "abc", "--ka"], HERE / "logs/full534-abc.log")
    heavy("abc:comp534", [*PY, str(REL / "comp534.py"), "abc", ",".join(CANCER), "--level", lv, "--ka"], HERE / "logs/comp534-abc.log")


def draw():
    light("draw", [*PY, str(REL / "draw534.py")], HERE / "logs/draw534.log")
    for p in PANELS:
        spec = dict(bundle=str(HERE / p / "bundle"), rules="def", cd4_root=str(HERE / p / "cd4/cd4_compact"),
                    full_op=str(HERE / p / "ops" / f"{p}__full3.npz"), folds=list(CANCER), synthetic=True, meta=None,
                    fold_ops={h: str(HERE / p / "ops" / f"{p}__held_{h}.npz") for h in CANCER})
        sp = HERE / p / "spec.json"
        if not sp.exists():
            sp.write_text(json.dumps(spec, indent=1))


def panel(p: str):
    light(f"{p}:cd4", [*PY, str(REL / "cd4_534.py"), p], HERE / f"logs/cd4_534-{p}.log")
    ops = HERE / p / "ops"
    ops.mkdir(parents=True, exist_ok=True)
    paths = {}
    for tag, donors in FITS.items():
        out = ops / f"{p}__{tag}.npz"
        paths[tag] = out
        step = f"{p}:fit:{tag}"
        if step in state():
            continue
        free = disk_free_gib()
        assert free > 10.0, f"disk free {free:.1f} GiB: a fit needs ~1.3 GB + a {'4.2' if tag == 'full3' else '3.1'} GB scratch"
        cmd = [*PY, "scripts/build_union_whitener.py", "--bundle", str(HERE / p / "bundle"), "--donors", ",".join(donors),
               "--basis-prune", "essential", "--out", str(out), "--scratch", str(ops / ".scratch"), "--max-rss-gib", "14"]
        heavy(step, cmd, HERE / f"logs/fit534-{p}-{tag}.log")
        s = state()
        s[step]["op_sha256"] = json.loads(out.with_suffix(".json").read_text())["sha256"]
        s[step]["command"] = " ".join(cmd)
        STATE.write_text(json.dumps(s, indent=1))
    if f"{p}:axes" not in state():
        import axes534
        rec = axes534.check_ops(p, paths)
        (HERE / p / "axes534.json").write_text(json.dumps(rec, indent=1))
        done(f"{p}:axes", rec)
        say(f"{p}: axes OK " + ", ".join(f"{t} {r['axis_sha256'][:12]}" for t, r in rec.items()))
    lv = str(REL / p / "level534.json")
    heavy(f"{p}:full534", [*PY, str(REL / "full534.py"), p], HERE / f"logs/full534-{p}.log")
    heavy(f"{p}:comp534", [*PY, str(REL / "comp534.py"), p, ",".join(CANCER), "--level", lv], HERE / f"logs/comp534-{p}.log")
    # verify, then free the disk: the four operators (~5.1 GB) are deterministic re-fits (seed 0), sha256 + command kept
    for f in CANCER:
        c = json.loads((HERE / p / f"comp534-{f}.json").read_text())
        assert "peak_rss_gib" in c and c["op_axis_checked"], (p, f)
    if f"{p}:ops_deleted" not in state():
        rec = {}
        for tag, path in paths.items():
            rec[tag] = dict(sha256=json.loads(path.with_suffix(".json").read_text())["sha256"], bytes=path.stat().st_size)
            path.unlink()
        done(f"{p}:ops_deleted", rec)
        say(f"{p}: priced; operators deleted to free {sum(r['bytes'] for r in rec.values()) / 2**30:.1f} GiB (jsons kept)")


def resume_note(msg: str):
    s = state()
    lines = ["# 534 Part B -- RESUME", "", f"Written {G.utc()} by drive534.py: {msg}.", "",
             "Done steps (work/partB/state534.json):", ""] + [f"- {k} ({v.get('done')})" for k, v in s.items()] + [
             "", "Next: re-run `uv run python experiments/534-generalization-audit/work/partB/drive534.py all` after "
             "2026-09-30 06:00Z. Done steps are skipped; the next heavy step waits for the gate.", ""]
    (HERE / "RESUME.md").write_text("\n".join(lines))


def main(what: str):
    try:
        if what in ("ka", "all"):
            ka()
        if what in ("draw", "all"):
            draw()
        for p in PANELS:
            if what in (p, "all", "r6"):
                panel(p)
        if what in ("abc", "all"):
            abc()
        say(f"{what}: complete")
    except G.WindowClosed as e:
        resume_note(f"the coordinator's window closed ({e})")
        say(f"WINDOW CLOSED ({e}); RESUME.md written")
        raise SystemExit(3)


if __name__ == "__main__":
    main(sys.argv[1])

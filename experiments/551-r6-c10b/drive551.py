"""551 (C10b) -- the driver: heavy steps strictly serial behind gate551, resumable (work/state551.json; a done step is skipped).

  ka      efp known answers (PREREGISTRATION section 4), before any R6 access and before F is known:
            full551 efp --ka  ->  comp551 efp <3 folds, one process> --level --ka
  abc     supplementary (not an R6 access; A/B/C's own panel is tuned-on): F's designs on abc, 534's abc KA rows
            full551 abc --ka --config F  ->  comp551 abc <3 folds> --level --ka --config F
  access  THE C10b ACCESS (needs work/c10b-config.json, recorded in the ledger first). Per R6 panel, S1..S4 in order:
            spec.json (534's bundle and CD4 banks, read-only; operators under 551/work) -> the four operator re-fits with
            534's exact commands and --expect-sha256 = C10a's sha256 -> axis check vs 534's axes534.json -> full551 P
            --config F -> comp551 P <3 folds> --level --repro --config F (KA-REPRO before any new number) -> operators
            deleted (sha256 + bytes kept). Then aggregate551 (light).

    uv run python drive551.py ka | abc | access
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard551  # noqa: E402,F401
import json  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import gate551 as G  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REL = HERE.relative_to(ROOT)
WORK = HERE / "work"
STATE = WORK / "state551.json"
W534 = ROOT / "experiments/534-generalization-audit/work/partB"
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
    print(f"[drive551 {G.utc()}]", *a, flush=True)


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


def ka():
    for bank in ("k562_gwps", "orion_hct116", "orion_hek293t", "perturbmatch"):   # one per process (062's 9 GiB guard)
        heavy(f"ka:bankstats:{bank}", [*PY, str(REL / "bankstats551.py"), bank], WORK / f"logs/bankstats551-{bank}.log")
    lv = str(REL / "work/efp/level551.json")
    heavy("ka:full551:efp", [*PY, str(REL / "full551.py"), "efp", "--ka"], WORK / "logs/full551-efp.log")
    heavy("ka:comp551:efp", [*PY, str(REL / "comp551.py"), "efp", ",".join(CANCER), "--level", lv, "--ka"],
          WORK / "logs/comp551-efp.log")


def config() -> str:
    c = json.loads((WORK / "c10b-config.json").read_text())
    return c["config"]


def abc():
    F = config()
    lv = str(REL / "work/abc/level551.json")
    heavy("abc:full551", [*PY, str(REL / "full551.py"), "abc", "--ka", "--config", F], WORK / "logs/full551-abc.log")
    heavy("abc:comp551", [*PY, str(REL / "comp551.py"), "abc", ",".join(CANCER), "--level", lv, "--ka", "--config", F],
          WORK / "logs/comp551-abc.log")


def panel(p: str, F: str):
    pw = WORK / p
    ops = pw / "ops"
    ops.mkdir(parents=True, exist_ok=True)
    spec = dict(bundle=str(W534 / p / "bundle"), rules="def", cd4_root=str(W534 / p / "cd4/cd4_compact"),
                full_op=str(ops / f"{p}__full3.npz"), folds=list(CANCER), synthetic=True, meta=None,
                fold_ops={h: str(ops / f"{p}__held_{h}.npz") for h in CANCER})
    (pw / "spec.json").write_text(json.dumps(spec, indent=1))
    s534 = json.loads((W534 / "state534.json").read_text())
    axes534 = json.loads((W534 / p / "axes534.json").read_text())
    paths = {}
    for tag, donors in FITS.items():
        out = ops / f"{p}__{tag}.npz"
        paths[tag] = out
        step = f"{p}:fit:{tag}"
        if step in state():
            continue
        want = s534[f"{p}:fit:{tag}"]["op_sha256"]
        free = shutil.disk_usage(ROOT).free / 2**30
        assert free > 10.0, f"disk free {free:.1f} GiB"
        cmd = [*PY, "scripts/build_union_whitener.py", "--bundle", str(W534 / p / "bundle"), "--donors", ",".join(donors),
               "--basis-prune", "essential", "--out", str(out), "--scratch", str(ops / ".scratch"), "--max-rss-gib", "14",
               "--expect-sha256", want]
        c534 = s534[f"{p}:fit:{tag}"]["command"].split()
        assert c534[c534.index("--donors") + 1] == ",".join(donors), (p, tag, "not 534's donors")
        heavy(step, cmd, WORK / f"logs/fit551-{p}-{tag}.log")
        got = json.loads(out.with_suffix(".json").read_text())
        assert got["sha256"] == want, f"{p} {tag}: re-fitted operator {got['sha256'][:12]} != C10a's {want[:12]}"
        assert got["axis_sha256"] == axes534[tag]["axis_sha256"], f"{p} {tag}: axis differs from C10a's"
        s = state()
        s[step]["op_sha256"] = got["sha256"]
        s[step]["command"] = " ".join(cmd)
        STATE.write_text(json.dumps(s, indent=1))
        say(f"{p} {tag}: operator re-fitted, sha256 {want[:12]} == C10a's, axis {got['axis_sha256'][:12]} == C10a's")
    lv = str(REL / f"work/{p}/level551.json")
    heavy(f"{p}:full551", [*PY, str(REL / "full551.py"), p, "--config", F], WORK / f"logs/full551-{p}.log")
    heavy(f"{p}:comp551", [*PY, str(REL / "comp551.py"), p, ",".join(CANCER), "--level", lv, "--repro", "--config", F],
          WORK / f"logs/comp551-{p}.log")
    for f in CANCER:
        c = json.loads((pw / f"comp551-{f}.json").read_text())
        assert "peak_rss_gib" in c and c["op_axis_checked"] and c["ka"]["repro"]["worst"] <= 1e-12, (p, f)
    if f"{p}:ops_deleted" not in state():
        rec = {}
        for tag, path in paths.items():        # resumable (555): an operator an interrupted run already deleted is recorded
            size = path.stat().st_size if path.exists() else None
            rec[tag] = dict(sha256=json.loads(path.with_suffix(".json").read_text())["sha256"], bytes=size)
            if path.exists():
                path.unlink()
        done(f"{p}:ops_deleted", rec)
        say(f"{p}: priced; operators deleted ({sum(r['bytes'] or 0 for r in rec.values()) / 2**30:.1f} GiB freed; jsons kept)")


def access():
    G.PHASE["name"] = "access"                  # heavy starts only inside the coordinator's C10b window
    F = config()
    import designs551 as D
    if D.CONFIGS[F][1] is None:
        say(f"F = {F} is a C10a design: C10b does not access R6 (PREREGISTRATION section 2). Nothing done.")
        return
    for p in PANELS:
        panel(p, F)
    light("aggregate551", [*PY, str(REL / "aggregate551.py")], WORK / "logs/aggregate551.log")


def main(what: str):
    WORK.mkdir(parents=True, exist_ok=True)
    try:
        {"ka": ka, "abc": abc, "access": access}[what]()
        say(f"{what}: complete")
    except G.WindowClosed as e:
        say(f"WINDOW CLOSED ({e}); re-run `uv run python {REL}/drive551.py {what}` when it reopens (done steps are skipped)")
        raise SystemExit(3)


if __name__ == "__main__":
    main(sys.argv[1])

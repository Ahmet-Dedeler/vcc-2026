"""233 Phase 0.1 / Phase 1.1 — re-fit the whitening operator for every (panel, fold).

AGENTS.md §8.3 / DEF-RUNBOOK §0.1: the operator's gene axis is `union(donor coverage) - panel`.
It is a FUNCTION OF THE PANEL and it fails SILENTLY on a new one -- every check the production
builder runs passes (n_cols, eigenvector shape, file sha256) while 98.8% of coordinates are
paired with the wrong gene.  212 priced that trap at -0.0429 overall on `pds` alone.  So a
retention measurement that reused today's operator on a new panel would measure the trap, not
the panel.

Ten fits, all through the production script `scripts/build_union_whitener.py` with identical
flags (`--basis-prune essential`, cap 12,000, seed 0) so nothing but `--bundle` (the panel) and
`--donors` (leave-one-donor-out) varies:

  KA  abc x {k562_gwps, orion_hct116, orion_hek293t}   -- the production 3-donor fit; MUST
      reproduce `e2f65ecf...` (DEF-RUNBOOK §4's known-answer row).
  9   {abc, real, pess} x each 2-donor leave-one-out basis, for the three PLOLO cancer folds.

Each fit's sidecar json is then given the field production's own metadata LACKS and whose
absence is the trap: `axis_sha256`, the sha256 of the operator's gene-axis column vector.
`panels233.load_frozen_operator_233` refuses to load an operator whose axis hash does not
match the panel it is asked for -- and `trap233.py` shows that assertion firing while every
production-side check still passes.

    uv run python whiten233.py --plan                 # dry-run every fit, no eigensolve
    uv run python whiten233.py --jobs ka abc:k562_gwps ...
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import panels233 as P  # noqa: E402

BUILDER = ROOT / "scripts/build_union_whitener.py"
KA_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
OPS = P.WORK / "ops"
OPS.mkdir(parents=True, exist_ok=True)


def jobs() -> dict[str, dict]:
    out = {"ka": dict(panel="abc", donors=("k562_gwps", "orion_hct116", "orion_hek293t"),
                      out=OPS / "abc__production3.npz", expect=KA_SHA)}
    for panel in ("abc", "real", "pess"):
        for held, srcs in P.FOLD_SOURCES.items():
            out[f"{panel}:{held}"] = dict(panel=panel, donors=srcs,
                                          out=P.op_path(panel, held), expect=None)
    return out


def memory_ok(min_free_pct: float = 40.0, min_free_gib: float = 20.0,
              min_swap_free_mb: float = 400.0) -> tuple[bool, str]:
    """Hold if the machine cannot absorb an ~9 GiB anonymous working set.

    Three other agents are live and 229 runs 20 GiB builds; 219 recorded that a job which
    exceeds memory here is killed with NO traceback, so this is checked before every fit.

    PRIMARY SIGNAL: `memory_pressure`'s system-wide free percentage plus free+inactive pages.
    `vm.swapusage`'s `free` field is the unused TAIL OF A DYNAMICALLY SIZED SWAP FILE, not a
    pressure reading -- macOS grows the file on demand while ~145 GiB of disk is free.  It sat
    at 0.8-1.5 GiB for this whole session while `memory_pressure` reported 86% free, and six
    9 GiB fits completed in that state with no eviction.  Gating on it alone would have
    stalled indefinitely without measuring anything.  It is still checked, at a floor.
    """
    mp = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    free_pct = 100.0
    for ln in mp.splitlines():
        if "free percentage" in ln:
            free_pct = float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    sw = subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout
    free_sw = float(sw.split("free =")[1].split("M")[0].strip())
    vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = 16384
    def pages(tag):
        for ln in vm.splitlines():
            if ln.startswith(tag):
                return int(ln.split(":")[1].strip().rstrip("."))
        return 0
    free_gib = (pages("Pages free") + pages("Pages inactive")
                + pages("Pages speculative") + pages("Pages purgeable")) * page / 2**30
    msg = (f"memory_pressure free {free_pct:.0f}%, free+inactive {free_gib:.1f} GiB, "
           f"swap_free {free_sw:.0f} MiB")
    return (free_pct >= min_free_pct and free_gib >= min_free_gib
            and free_sw >= min_swap_free_mb), msg


def fit(name: str, spec: dict, dry: bool, max_rss: float) -> dict:
    cmd = ["uv", "run", "python", str(BUILDER),
           "--bundle", str(P.BUNDLES[spec["panel"]]),
           "--donors", ",".join(spec["donors"]),
           "--basis-prune", "essential",
           "--out", str(spec["out"]),
           "--scratch", str(OPS / ".scratch"),
           "--max-rss-gib", str(max_rss)]
    if spec["expect"]:
        cmd += ["--expect-sha256", spec["expect"]]
    if dry:
        cmd += ["--dry-run"]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    tag = f"{name.replace(':', '__')}{'-dry' if dry else ''}"
    (P.WORK / f"whiten-{tag}.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    rec = dict(job=name, panel=spec["panel"], donors=list(spec["donors"]),
               rc=r.returncode, seconds=round(time.time() - t0, 1),
               out=str(spec["out"].relative_to(ROOT)), dry=dry,
               stdout=r.stdout.strip().splitlines()[-8:],
               stderr=r.stderr.strip().splitlines()[-4:])
    if spec["expect"]:
        rec["known_answer"] = "MATCH" if "MATCH" in r.stdout else (
            "MISMATCH" if "MISMATCH" in r.stdout else "n/a (dry)")
    if not dry and r.returncode == 0:
        rec.update(stamp(spec))
    return rec


def stamp(spec: dict) -> dict:
    """Add the axis hash production's metadata lacks, and the panel it was fitted for."""
    from plolo import E062  # noqa: F401  (only for Bank.load)
    import plolo
    banks = [plolo.E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in spec["donors"]]
    sha, n = P.axis_sha256(set(P.read_panel(spec["panel"])), banks)
    jpath = spec["out"].with_suffix(".json")
    meta = json.loads(jpath.read_text())
    assert int(meta["n_cols"]) == n, f"n_cols {meta['n_cols']} != recomputed axis {n}"
    meta["axis_sha256"] = sha
    meta["panel"] = spec["panel"]
    meta["panel_csv_sha256"] = P.hashlib.sha256(
        (P.BUNDLES[spec["panel"]] / "pert_counts.csv").read_bytes()).hexdigest()
    jpath.write_text(json.dumps(meta, indent=2))
    return dict(n_cols=n, axis_sha256=sha[:16], operator_sha256=meta.get("sha256", "")[:16])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", action="store_true", help="dry-run every fit")
    ap.add_argument("--jobs", nargs="*", default=None)
    ap.add_argument("--max-rss-gib", type=float, default=13.0)
    ap.add_argument("--no-gate", action="store_true")
    args = ap.parse_args()

    all_jobs = jobs()
    names = args.jobs if args.jobs else list(all_jobs)
    out_path = P.WORK / ("whiten233-plan.json" if args.plan else "whiten233.json")
    rows = json.loads(out_path.read_text()) if out_path.exists() else {}

    for name in names:
        spec = all_jobs[name]
        if not args.plan and spec["out"].exists() and name in rows and rows[name]["rc"] == 0:
            print(f"[{name}] already fitted -> {spec['out'].name}")
            continue
        if not args.plan and not args.no_gate:
            for attempt in range(120):
                ok, msg = memory_ok()
                if ok:
                    break
                print(f"[{name}] HOLD ({msg}) — 229 runs 20 GiB builds; waiting 60 s",
                      flush=True)
                time.sleep(60)
            else:
                raise SystemExit(f"[{name}] memory never freed: {msg}")
            print(f"[{name}] memory OK ({msg}); fitting {spec['donors']} on panel "
                  f"{spec['panel']}", flush=True)
        rec = fit(name, spec, args.plan, args.max_rss_gib)
        rows[name] = rec
        out_path.write_text(json.dumps(rows, indent=2))
        head = [l for l in rec["stdout"] if "union axis" in l or "basis rows" in l
                or "sha256" in l or "MATCH" in l]
        print(f"[{name}] rc {rec['rc']} {rec['seconds']}s  " + " | ".join(head), flush=True)
        if rec["rc"] != 0:
            for l in rec["stderr"]:
                print("   ERR", l[:200])
            raise SystemExit(f"fit failed: {name}")
    print(f"-> {out_path.name}")


if __name__ == "__main__":
    main()

"""238 Phase 1c — re-fit the whitening operator for every (essential-free panel, fold).

AGENTS.md §8.3 / `DEF-RUNBOOK` §0.1: the operator's gene axis is `union(donor coverage) -
panel`.  It is a FUNCTION OF THE PANEL and it fails SILENTLY on a new one -- every check the
production builder runs passes (`n_cols` is 17,942 on every panel on disk) while 98.8% of the
coordinates are paired with the wrong gene.  233 §2 priced that at **-0.021 overall on `pds`
alone**.  This is falsifier 1 of the mandate: *"operator re-fit per panel with the axis
asserted -- it fails silently otherwise."*

Seven fits, all through the production script `scripts/build_union_whitener.py` with flags
identical to 233's, so nothing but `--bundle` (the panel) and `--donors` (leave-one-donor-out)
varies:

  ka   abc x all three donors           -- MUST reproduce `e2f65ecf...` (DEF-RUNBOOK §4)
  6    {efp, efr} x each 2-donor leave-one-out basis, for the three PLOLO cancer folds

Each fit's sidecar json is then given the field production's own metadata lacks and whose
absence is the trap: `axis_sha256`.  `panels238.load_frozen_operator_238` refuses to load an
operator whose axis hash does not match the panel it is asked for.

MEMORY.  ~8 GiB anonymous at the eigensolve, 13-14 GiB process high-water with the 4.2 GB
basis memmap resident.  Other agents are live, so every fit is gated on
`memory_pressure >= 40% free AND free+inactive >= 20 GiB AND swap free >= 400 MiB` and run
STRICTLY SERIALLY.  A job that exceeds memory on this box is killed without a traceback.

    uv run python whiten238.py --plan          # dry-run every fit, no eigensolve
    uv run python whiten238.py                 # fit them, serially, gated
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import panels238 as P  # noqa: E402

BUILDER = ROOT / "scripts/build_union_whitener.py"
KA_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_AXIS_SHA = "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
OPS = P.WORK / "ops"
OPS.mkdir(parents=True, exist_ok=True)


def jobs() -> dict[str, dict]:
    out = {"ka": dict(panel="abc", donors=("k562_gwps", "orion_hct116", "orion_hek293t"),
                      out=OPS / "abc__production3.npz", expect=KA_SHA)}
    for panel in P.OWN:
        for held, srcs in P.FOLD_SOURCES.items():
            out[f"{panel}:{held}"] = dict(panel=panel, donors=srcs,
                                          out=P.op_path(panel, held), expect=None)
    return out


def memory_ok(min_free_pct: float = 40.0, min_free_gib: float = 20.0,
              min_swap_free_mb: float = 400.0) -> tuple[bool, str]:
    """The mandate's gate.  `vm.swapusage free` is the unused tail of a DYNAMICALLY SIZED swap
    file, not a pressure reading (233 measured it at 0.8-1.7 GiB while `memory_pressure`
    reported 86% free), so it is checked at a floor and the pressure reading is primary."""
    mp = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    free_pct = 100.0
    for ln in mp.splitlines():
        if "free percentage" in ln:
            free_pct = float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    sw = subprocess.run(["sysctl", "vm.swapusage"], capture_output=True, text=True).stdout
    free_sw = float(sw.split("free =")[1].split("M")[0].strip())
    vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout

    def pages(tag: str) -> int:
        for ln in vm.splitlines():
            if ln.startswith(tag):
                return int(ln.split(":")[1].strip().rstrip("."))
        return 0

    free_gib = (pages("Pages free") + pages("Pages inactive") + pages("Pages speculative")
                + pages("Pages purgeable")) * 16384 / 2**30
    msg = (f"memory_pressure free {free_pct:.0f}%, free+inactive {free_gib:.1f} GiB, "
           f"swap_free {free_sw:.0f} MiB")
    return (free_pct >= min_free_pct and free_gib >= min_free_gib
            and free_sw >= min_swap_free_mb), msg


def fit(name: str, spec: dict, dry: bool, max_rss: float) -> dict:
    cmd = ["/usr/bin/time", "-l", "uv", "run", "python", str(BUILDER),
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
    peak = None
    for ln in r.stderr.splitlines():
        if "maximum resident set size" in ln:
            peak = round(int(ln.strip().split()[0]) / 2**30, 2)
    rec = dict(job=name, panel=spec["panel"], donors=list(spec["donors"]), rc=r.returncode,
               seconds=round(time.time() - t0, 1), peak_rss_gib=peak,
               out=str(spec["out"]).replace(str(ROOT) + "/", ""), dry=dry,
               stdout=r.stdout.strip().splitlines()[-10:],
               stderr=[l for l in r.stderr.strip().splitlines()[-40:]
                       if "resident" in l or "Error" in l or "error" in l][-4:])
    if spec["expect"]:
        rec["known_answer"] = ("MATCH" if "MATCH" in r.stdout else
                               "MISMATCH" if "MISMATCH" in r.stdout else "n/a (dry)")
    if not dry and r.returncode == 0:
        rec.update(stamp(spec))
    return rec


def stamp(spec: dict) -> dict:
    """Add the axis hash production's metadata lacks, and the panel it was fitted for."""
    import plolo
    from plolo import E062
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in spec["donors"]]
    sha, n = P.axis_sha256(set(P.read_panel(spec["panel"])), banks)
    jpath = spec["out"].with_suffix(".json")
    meta = json.loads(jpath.read_text())
    assert int(meta["n_cols"]) == n, f"n_cols {meta['n_cols']} != recomputed axis {n}"
    if spec["panel"] in P.OWN and sha == ABC_AXIS_SHA:
        raise SystemExit(f"{spec['out']}: axis hash equals A/B/C's — the bundle is wrong")
    meta["axis_sha256"] = sha
    meta["panel"] = spec["panel"]
    meta["panel_csv_sha256"] = hashlib.sha256(
        (P.BUNDLES[spec["panel"]] / "pert_counts.csv").read_bytes()).hexdigest()
    jpath.write_text(json.dumps(meta, indent=2))
    return dict(n_cols=n, axis_sha256=sha[:16], operator_sha256=meta.get("sha256", "")[:16])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--jobs", nargs="*", default=None)
    ap.add_argument("--max-rss-gib", type=float, default=14.0)
    args = ap.parse_args()

    all_jobs = jobs()
    names = args.jobs if args.jobs else list(all_jobs)
    out_path = P.WORK / ("whiten238-plan.json" if args.plan else "whiten238.json")
    rows = json.loads(out_path.read_text()) if out_path.exists() else {}

    for name in names:
        spec = all_jobs[name]
        if not args.plan and spec["out"].exists() and name in rows and rows[name]["rc"] == 0:
            print(f"[{name}] already fitted -> {spec['out'].name}", flush=True)
            continue
        if not args.plan:
            for _ in range(180):
                ok, msg = memory_ok()
                if ok:
                    break
                print(f"[{name}] HOLD ({msg}) — waiting 60 s", flush=True)
                time.sleep(60)
            else:
                raise SystemExit(f"[{name}] memory never freed: {msg}")
            print(f"[{name}] memory OK ({msg}); fitting {spec['donors']} on panel "
                  f"{spec['panel']}", flush=True)
        rec = fit(name, spec, args.plan, args.max_rss_gib)
        rows[name] = rec
        out_path.write_text(json.dumps(rows, indent=2))
        head = [l for l in rec["stdout"] if "union axis" in l or "basis rows" in l
                or "axis sha256" in l or "MATCH" in l]
        print(f"[{name}] rc {rec['rc']} {rec['seconds']}s peak {rec['peak_rss_gib']} GiB  "
              + " | ".join(head), flush=True)
        if rec["rc"] != 0:
            for l in rec["stderr"]:
                print("   ERR", l[:200])
            raise SystemExit(f"fit failed: {name}")
    print(f"-> {out_path.name}")


if __name__ == "__main__":
    main()

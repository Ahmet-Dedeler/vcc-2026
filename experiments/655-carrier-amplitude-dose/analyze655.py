#!/usr/bin/env python3
"""655 — the dose curve under zero repair: x1 (654's base) / x1.25 (654's zr) / x1.5 / x2.0 (655's zr runs), seeds 0-2.
Uses 654's analyze654 functions verbatim (members at 521's regimes with the realised order, overall-eq by lib521; 631's
regions on lib521.nmae_rows' conventions), so 654 and 655 are read by identical code.
usage: uv run python experiments/655-carrier-amplitude-dose/analyze655.py"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.argv = [sys.argv[0], "--zr"]                       # 654's module reads --zr at import (the production-faithful files)
spec = importlib.util.spec_from_file_location("a654", ROOT / "experiments/654-carrier-amplitude/analyze654.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
SEEDS = (0, 1, 2)
DOSES = {"x125": 1.25, "x150": 1.5, "x200": 2.0}


def stats_path(design, s):
    if design in ("base", "x125"):
        return A.path(design, s)                        # 654's production-faithful files (base_zr == base, A1's KA)
    hits = sorted((HERE / "work/h1_cpm").glob(f"stats-{design}_s{s}_bps*zr-5*_dnk250zs.npz"))
    return hits[0] if len(hits) == 1 else None


def load(design, s):
    p = stats_path(design, s)
    if p is None or not Path(p).exists():
        return None
    z = np.load(p, allow_pickle=True)
    out = {k: np.asarray(z[k], float) for k in ("q", "p", "lfc")}
    out["sums"] = np.asarray(z["sums"])
    out["zero_total"] = int((out["lfc"] < -30).sum())
    return out


def main() -> int:
    lines, res = [], {"deltas": {}, "regions": {}, "calls": {}, "zero_total": {}}

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    runs = {(d, s): load(d, s) for d in ("base", *DOSES) for s in SEEDS}
    runs = {k: v for k, v in runs.items() if v is not None}
    say(f"runs present: {sorted(runs)}")
    for (d, s), r in sorted(runs.items()):
        res["calls"][f"{d}|{s}"] = float(((r["q"] < 0.05) & ~A.OWN).sum(1).mean())
        res["zero_total"][f"{d}|{s}"] = r["zero_total"]
    say("calls/target and zero-total cells (own targets included; zero repair never touches the own target):")
    for d in ("base", *DOSES):
        say(f"  {d}: " + ", ".join(f"s{s} {res['calls'].get(f'{d}|{s}', float('nan')):.1f} ({res['zero_total'].get(f'{d}|{s}', '-')})" for s in SEEDS))
    say("\nmembers vs x1 (paired per seed), overall-eq (192 INV); 611's call-set rule noted in RESULTS:")
    for d in DOSES:
        for reg in ("1.0", "0.5", "0.35", "cmp"):
            rows = []
            for s in SEEDS:
                if (d, s) in runs and ("base", s) in runs:
                    m1, m0 = A.members(runs[(d, s)], reg), A.members(runs[("base", s)], reg)
                    dd = {k: m1[k] - m0[k] for k in m1}
                    dd["overall"] = A.L.overall_delta(m1, m0)
                    rows.append(dd)
                    res["deltas"][f"{d}|{reg}|{s}"] = dd
            if rows:
                mean = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
                res["deltas"][f"{d}|{reg}|mean"] = mean | {"n": len(rows), "seeds": [r["overall"] for r in rows]}
                say(f"  {d} f {reg:4s}: overall {mean['overall']:+.5f} seeds {[round(r['overall'], 5) for r in rows]}  "
                    f"fid {mean['fid']:+.5f} jac {mean['jac']:+.5f} reach {mean['reach']:+.5f} nmae {mean['nmae']:+.5f}")
    say("\nhead / block / uncalled x direction at f 0.35: nmae share and mean delivered |lfc| (mean over seeds)")
    for d in ("base", *DOSES):
        tabs = [A.regions_table(runs[(d, s)], runs[("base", s)], 0.35) for s in SEEDS if (d, s) in runs and ("base", s) in runs]
        if not tabs:
            continue
        for key in ("head|DOWN", "head|UP", "block|DOWN", "block|UP", "uncalled|DOWN", "uncalled|UP"):
            sh = float(np.mean([t[key]["share"] for t in tabs]))
            mg = float(np.nanmean([t[key]["mag"] for t in tabs]))
            res["regions"][f"{d}|{key}"] = dict(share=sh, mag=mg)
        say(f"  {d:5s}: " + "  ".join(f"{k} {res['regions'][f'{d}|{k}']['share']:.4f}/{res['regions'][f'{d}|{k}']['mag']:.3f}"
                                      for k in ("head|DOWN", "block|DOWN", "uncalled|DOWN", "uncalled|UP")))
    say("\npreregistered reading (655 bins: s* = argmax over {1.25, 1.5, 2.0} of the f 0.35 mean; read at f 0.35 AND 0.5)")
    g = lambda d, reg: res["deltas"].get(f"{d}|{reg}|mean", {}).get("overall", float("nan"))  # noqa: E731
    seeds_same = lambda d: all(np.sign(res["deltas"].get(f"{d}|{reg}|{s}", {}).get("overall", np.nan)) == np.sign(g(d, reg))  # noqa: E731
                               for reg in ("0.35", "0.5") for s in SEEDS)
    cand = [d for d in DOSES if np.isfinite(g(d, "0.35"))]
    if cand:
        best = max(cand, key=lambda d: g(d, "0.35"))
        lo = min(g(best, "0.35"), g(best, "0.5"))
        same = seeds_same(best)
        b = ("CANDIDATE (>= +0.0015 at f 0.35 and 0.5, seeds agree)" if same and lo >= 0.0015 else
             "WEAK CANDIDATE (+0.0005 ... +0.0015)" if same and lo >= 0.0005 else "SMALL / CLOSE (< +0.0005 or seeds disagree)")
        d1 = (g("x150", "0.35") > g("x125", "0.35")) and (g("x200", "0.35") < g("x150", "0.35"))
        say(f"  s* = {best} (x{DOSES[best]}): f 0.35 {g(best, '0.35'):+.5f}, f 0.5 {g(best, '0.5'):+.5f}, seeds agree {same} -> {b}")
        say(f"  D1 (interior peak above 1.25: x1.5 > x1.25 and x2.0 < x1.5 at f 0.35): {'HOLDS' if d1 else 'FAILS'} "
            f"({g('x125', '0.35'):+.5f} / {g('x150', '0.35'):+.5f} / {g('x200', '0.35'):+.5f})")
        res["verdict"] = dict(s_star=best, bin=b, D1=bool(d1))
    (HERE / "work/analyze655.json").write_text(json.dumps(res, indent=1, default=float))
    (HERE / "work/analyze655.out").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

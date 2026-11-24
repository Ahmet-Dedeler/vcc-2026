#!/usr/bin/env python3
"""599: price the H1 real-emitter runs. DE members per run with its OWN realized order (scorer's rule) by 521's exact member
functions at its four regimes (f 1 / 0.5 / 0.35 / cmp, table521.regime_pool), as 612 does; overall-equivalent of the four
DE members at 192's board-fitted INV (lib521.overall_delta). Bulk members (pds, mse) are the runs' own raw values (each
run's pseudobulk vs h1-test); their overall-equivalent uses the BRIEF's mid anchors (pds 1/(0.9525 - 0.5); mse
-1/(0.989 - 0.0365)), approximate. Pairs: tc - t (the read), t - s175 and tc - s175 (each rule against the shipped scalar
render). Post-zs (the October DE design) is primary; pre-zs is reported beside it. The ballast K each run used is checked
against prep599's registered K (it is deterministic).
usage: uv run python experiments/599-tercilec-h1-de/score599.py"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import table521 as TB  # noqa: E402

TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
INV_BULK = {"pds": 1.0 / (0.9525 - 0.5), "mse": -1.0 / (0.989 - 0.0365)}
SEEDS = (0, 1, 2, 3, 4)   # 3 and 4: the preregistered extension (t and tc only), when present
S = L.load_state(0)
own, rl = S["own"], S["rl"]
adj = L.adjud(S)
POOLS = {reg: TB.regime_pool(S, reg)[0] for reg in TB.REGIMES}
PREP = dict(np.load(HERE / "work/prep599.npz", allow_pickle=False))
if (HERE / "work/prepx599.npz").exists():
    PREP.update(dict(np.load(HERE / "work/prepx599.npz", allow_pickle=False)))


def load(mode, s, post):
    t = TAG.format(s=s) + ("_dnk250zs" if post else "")
    p = HERE / f"work/h1_{mode}/stats-{t}.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    m = json.loads((HERE / f"work/h1_{mode}/metrics-{t}.json").read_text())["raw"]
    rec = json.loads((HERE / f"work/h1_{mode}/run599_s{s}.json").read_text())
    K = np.array([rec["per_target"][str(i)]["K"] for i in range(len(rec["per_target"]))])
    return dict(q=np.asarray(z["q"], float), p=np.asarray(z["p"], float), lfc=np.asarray(z["lfc"], float),
                pds=float(m["pds"]), mse=float(m["mse"]), K=K)


def members(run, reg):
    pool = POOLS[reg]
    call = (run["q"] < 0.05) & ~own
    St = dict(S); St["oq"], St["op"], St["ol"] = run["q"], run["p"], run["lfc"]
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], adj[i], run["lfc"][i], rl[i]) for i in range(len(S["targets"]))])
    return dict(fid=float(np.nanmean(L.fid_rows(call, run["lfc"], pool, adj, rl))), jac=float(np.mean(L.jac_rows(call, pool))),
                reach=float(np.nanmean(rr)), nmae=float(np.nanmean(L.nmae_rows(run["lfc"], pool, rl))))


out = {"K_matches_prep": {}, "pairs": {}}
for mode in ("s175", "t", "tc"):
    for s in SEEDS:
        run = load(mode, s, True)
        if run is not None:
            want = PREP.get(f"K_{mode}_s{s}")
            out["K_matches_prep"][f"{mode}_s{s}"] = bool(want is not None and np.array_equal(run["K"], want))
print(f"ballast K == prep599's registered K: {out['K_matches_prep']}")
for post in (True, False):
    for arm, ref in (("tc", "t"), ("t", "s175"), ("tc", "s175")):
        rows = []
        for s in SEEDS:
            a, b = load(arm, s, post), load(ref, s, post)
            if a is None or b is None:
                continue
            rows.append(dict(s=s, de={reg: (members(a, reg), members(b, reg)) for reg in TB.REGIMES},
                             pds=a["pds"] - b["pds"], mse=a["mse"] - b["mse"], dK=float(a["K"].mean() - b["K"].mean())))
        if not rows:
            continue
        key = f"{arm}-{ref}{'' if post else ' (pre-zs)'}"
        bulk_ov = [(r["pds"] * INV_BULK["pds"] + r["mse"] * INV_BULK["mse"]) / 6.0 for r in rows]
        res = {"seeds": [r["s"] for r in rows], "dK_mean": float(np.mean([r["dK"] for r in rows])),
               "pds": [r["pds"] for r in rows], "mse": [r["mse"] for r in rows],
               "bulk_overall_eq_approx": {"mean": float(np.mean(bulk_ov)), "per_seed": bulk_ov}, "regimes": {}}
        print(f"\n{key} ({len(rows)} seeds): mean dK {res['dK_mean']:+.3f}; pds {np.mean(res['pds']):+.5f} "
              f"{[round(x, 5) for x in res['pds']]} mse {np.mean(res['mse']):+.5f} {[round(x, 5) for x in res['mse']]} "
              f"| bulk overall-eq ~{np.mean(bulk_ov):+.5f}")
        for reg in TB.REGIMES:
            d = {m: float(np.mean([r["de"][reg][0][m] - r["de"][reg][1][m] for r in rows])) for m in ("fid", "jac", "reach", "nmae")}
            ov = [L.overall_delta(r["de"][reg][0], r["de"][reg][1]) for r in rows]
            per = {m: [float(r["de"][reg][0][m] - r["de"][reg][1][m]) for r in rows] for m in ("nmae", "fid", "jac", "reach")}
            res["regimes"][reg] = dict(d, de_overall=float(np.mean(ov)), de_overall_sd=float(np.std(ov)) if len(ov) > 1 else None,
                                       de_overall_per_seed=[float(x) for x in ov], per_seed=per)
            print(f"  {reg:4s} fid {d['fid']:+.5f} jac {d['jac']:+.5f} reach {d['reach']:+.5f} nmae {d['nmae']:+.5f} "
                  f"| DE overall {np.mean(ov):+.5f} (per seed {[round(x, 5) for x in ov]})")
        out["pairs"][key] = res
(HERE / "work/score599.json").write_text(json.dumps(out, indent=1) + "\n")

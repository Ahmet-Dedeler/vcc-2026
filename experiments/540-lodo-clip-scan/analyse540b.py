"""540b pricer — the cap extension (reads abc540b.json / arms540b-<fold>.json).

540 pricer — each LODO clip vs 446bj2 (clip 0.5, 3), board-calibrated on the board's own key step.

analyse492g.py's arithmetic (446bj2-anchored 187 cliff model; X from 446bj2's board raw mse; base and arm each re-levelled
at their EV-optimal factors on their own realised energy maps), with the clip arms' energies and maps from abc540b.json and
the PLOLO direction members from arms540-<fold>.json (realised 16, vs V1_0.5_3).
Calibration: this instrument reproduces the board's key step 446b -> 446bj2 (K3 clip 0.3 -> LODO clip 0.5) at k_X 1.31
(analyse492g's retrodiction), and its literal pds x 0.34 reproduced that step's pds (+0.00024 vs +0.00021) where the
board-emulated pds did not (7x too high). A clip change is the same class, so: pds = 0.34 x literal PLOLO dpds; mse at
k_X 1.31 (1.0 as the conservative end); minus the ballast's nmae cost (K model and displacement model, ballast492g's
slopes), -0.2808 overall per unit raw nmae.

    uv run python experiments/540-lodo-clip-scan/analyse540b.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/418-bulk-composition"))
sys.path.insert(0, str(ROOT / "experiments/492-per-target-sharpening"))
import analyse418 as A8  # noqa: E402
from analyse492g import CANCER, CTX, RAW_446BJ2, level_opt  # noqa: E402

FOLDS = (*CANCER, "kolf2")
BASE = "V1_0.5_3"
K_LIT, K_KEY_FIT = 0.34, 1.31
KS = (0.61, 1.0, K_KEY_FIT, A8.K_KEY)
NMAE_OV = -0.2808


def main():
    cl = A8.Cliff()
    g = json.loads((ROOT / "experiments/492-per-target-sharpening/work/analyse492g.json").read_text())
    bal = json.loads((ROOT / "experiments/492-per-target-sharpening/work/ballast492g.json").read_text())
    slope_K, slope_D = bal["fits"]["dK_mean"]["slope"], bal["fits"]["D_from_design"]["slope"]
    AB0 = g["AB_446bj2"]
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    assert abs(X0 - g["X_446bj2"]) < 1e-12
    abc = json.loads((WORK / "abc540b.json").read_text())
    fs = [0.92, 0.95, 0.98, 1.0, 1.02, 1.05, 1.08]

    def ratio_map(tag, key):
        mp = abc["maps"][tag]
        return {c: [mp[f"{f:g}"][c][key] / mp["1"][c][key] for f in fs] for c in CTX}

    f0, raw0 = level_opt(cl, AB0, X0, fs, ratio_map(BASE, "A_own"), ratio_map(BASE, "B_own"))
    R = {f: json.loads((WORK / f"arms540b-{f}.json").read_text()) for f in FOLDS if (WORK / f"arms540b-{f}.json").exists()}
    print(f"446bj2 anchor X {X0:.6f}; base re-level {f0}; EV raw {raw0:.5f}; folds {list(R)}")
    out = dict(base_level=f0, base_raw=raw0, arms={})
    for tg, a in abc["arms"].items():
        if tg == BASE or not R:
            continue
        ABa = {c: dict(A=AB0[c]["A"] * a[c]["A_ratio"], B=AB0[c]["B"] * a[c]["B_ratio"]) for c in CTX}
        amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
        ram, rbm = ratio_map(tg, "A_own"), ratio_map(tg, "B_own")
        rec = dict(lo=a["lo"], hi=a["hi"], folds={})
        for f, Rf in R.items():
            b, x = Rf["arms"][BASE]["realised"]["16"], Rf["arms"][tg]["realised"]["16"]
            r = x["cos_agg"] / b["cos_agg"] - 1
            row = dict(r=r, dpds_bd=x["pds_lam0.5"] - b["pds_lam0.5"], dpds_lit=x["pds"] - b["pds"])
            for k in KS:
                fa, rawa = level_opt(cl, ABa, X0 * (1 + k * r) * amp, fs, ram, rbm)
                row[f"dmse_k{k:.3f}"] = rawa - raw0
            rec["folds"][f] = row
        canc = [f for f in CANCER if f in rec["folds"]]
        mean = lambda key: float(np.mean([rec["folds"][f][key] for f in canc]))
        rec["cancer"] = {key: mean(key) for key in ("r", "dpds_bd", "dpds_lit", *[f"dmse_k{k:.3f}" for k in KS])}
        rec["cancer_pos_r"] = int(sum(rec["folds"][f]["r"] > 0 for f in canc))
        nK, nD = slope_K * a["dK_mean"], slope_D * a["D_from_design"]
        dp = K_LIT * rec["cancer"]["dpds_lit"]
        cal = {}
        for kk in (K_KEY_FIT, 1.0):
            dm = rec["cancer"][f"dmse_k{kk:.3f}"]
            gross = (A8.PDS_SC * dp + A8.MSE_SC * dm) / 6
            cal[f"k_X {kk:g}"] = dict(dpds_raw=dp, dmse_raw=dm, gross=gross, net_K=gross + NMAE_OV * nK, net_D=gross + NMAE_OV * nD)
        rec["calibrated"] = dict(nmae_raw_K=nK, nmae_raw_D=nD, **cal)
        out["arms"][tg] = rec
        c1, c2 = cal[f"k_X {K_KEY_FIT:g}"], cal["k_X 1"]
        kol = rec["folds"].get("kolf2", {})
        print(f"{tg:12s} r {100 * rec['cancer']['r']:+.2f} % ({rec['cancer_pos_r']}/{len(canc)}) kolf2 {100 * kol.get('r', np.nan):+.2f} % | "
              f"dpds lit {rec['cancer']['dpds_lit']:+.5f} bd {rec['cancer']['dpds_bd']:+.5f} | E {[round(a[c]['A_ratio'], 4) for c in CTX]} | "
              f"CAL pds {dp:+.5f} mse {c1['dmse_raw']:+.5f} ({c2['dmse_raw']:+.5f} at k1) nmae K {nK:+.5f} D {nD:+.5f} -> overall "
              f"{c1['net_K']:+.5f} / {c1['net_D']:+.5f} (k1: {c2['net_K']:+.5f} / {c2['net_D']:+.5f})", flush=True)
    (WORK / "analyse540b.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

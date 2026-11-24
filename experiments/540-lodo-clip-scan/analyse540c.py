"""540c pricer — do 479 (P2g / S_P2) and 495 (cap 6) stack? PLOLO direction members from arms540c-<fold>.json (realised 16),
energies and level maps from abc540c.json, 446bj2-anchored cliff model (analyse492g's arithmetic).

Each combination is priced as two steps, each at its own board-calibrated class:
  step 1  446bj2 -> cap 6            (key class: k_X 1.31; 540's price, from analyse540b.json)
  step 2  cap 6 -> cap 6 + P2g/S_P2  (sharpening class: k_X 0.61 / 0.12, the per-target arms' calibration)
with step 2's cliff anchored at the cap-6 state: energies x abc540c's cap-6 ratios, X = X_446bj2 (1 + 1.31 r_cap6) amp_cap6.
The same step-2 arithmetic on the cap-3 base reproduces 492g's P2g / S_P2 gross prices (the known answer of this pricer).
pds = 0.34 x literal PLOLO dpds; ballast nmae by the K model (ballast492g's MDE-era slope) on abc540c's K means.

    uv run python experiments/540-lodo-clip-scan/analyse540c.py
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
K_LIT, K_CAP = 0.34, 1.31
NMAE_OV = -0.2808


def main():
    cl = A8.Cliff()
    g = json.loads((ROOT / "experiments/492-per-target-sharpening/work/analyse492g.json").read_text())
    slope_K = json.loads((ROOT / "experiments/492-per-target-sharpening/work/ballast492g.json").read_text())["fits"]["dK_mean"]["slope"]
    abc = json.loads((WORK / "abc540c.json").read_text())
    AB0 = g["AB_446bj2"]
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    fs = [float(f) for f in abc["level_maps"]["factors"]]
    R = {f: json.loads((WORK / f"arms540c-{f}.json").read_text()) for f in FOLDS if (WORK / f"arms540c-{f}.json").exists()}
    print(f"folds {list(R)}")
    ea = {k: {c: abc["arms"][k][c]["energy_ratio_realised"] for c in CTX} for k in abc["arms"]}
    eb = {k: {c: abc["arms"][k][c]["B_ratio_realised"] for c in CTX} for k in abc["arms"]}
    ea["cap3"] = {c: 1.0 for c in CTX}
    eb["cap3"] = {c: 1.0 for c in CTX}
    Kmean = {k: abc["arms"][k]["K_mean"] for k in abc["arms"]}
    Kmean["cap3"] = abc["base_K_mean"]
    b540 = json.loads((WORK / "analyse540b.json").read_text())["arms"]["V1_0.5_6"]
    r_cap6 = b540["cancer"]["r"]

    m492 = json.loads((ROOT / "experiments/492-per-target-sharpening/work/abc492g.json").read_text())["level_maps"]["maps"]

    def maps_of(tag):
        # cap-3 states use 492g's own maps (so the cap-3 increments reproduce 492g's prices), cap-6 states abc540c's
        mp = {"cap3": m492["base"], "cap3 + P2g": m492["P2g"], "cap3 + S_P2": m492["S_P2"]}.get(tag) or abc["level_maps"]["maps"][tag]
        ra = {c: [mp[f"{f:g}"][c]["realised"]["A_own"] / mp["1"][c]["realised"]["A_own"] for f in fs] for c in CTX}
        rb = {c: [mp[f"{f:g}"][c]["realised"]["B_own"] / mp["1"][c]["realised"]["B_own"] for f in fs] for c in CTX}
        return ra, rb

    out = {}
    for base, amp_r in (("cap3", 0.0), ("cap6", r_cap6)):
        ABb = {c: dict(A=AB0[c]["A"] * ea[base][c], B=AB0[c]["B"] * eb[base][c]) for c in CTX}
        ampb = np.sqrt(np.mean([ABb[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
        Xb = X0 * (1 + K_CAP * amp_r) * ampb
        rab, rbb = maps_of(base)
        _, rawb = level_opt(cl, ABb, Xb, fs, rab, rbb)
        for arm, maptag in (("P2g", "cap6 + P2g"), ("S_P2", "cap6 + S_P2")):
            key = f"{base} + {arm}"
            base_name = "V1_0.5_3" if base == "cap3" else "V1_0.5_6"
            arm_name = f"{base_name} + {arm}"
            rs, dps = [], []
            for f in CANCER:
                if f not in R:
                    continue
                b, a = R[f]["arms"][base_name]["realised"]["16"], R[f]["arms"][arm_name]["realised"]["16"]
                rs.append(a["cos_agg"] / b["cos_agg"] - 1)
                dps.append(a["pds"] - b["pds"])
            r, dpl = float(np.mean(rs)), float(np.mean(dps))
            ABa = {c: dict(A=AB0[c]["A"] * ea[key][c], B=AB0[c]["B"] * eb[key][c]) for c in CTX}
            ampa = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([ABb[c]["A"] for c in CTX]))
            raa, rba = maps_of(key)
            res = dict(r=r, dpds_lit=dpl, folds_pos=int(sum(x > 0 for x in rs)))
            for kx in (0.61, 0.12):
                _, rawa = level_opt(cl, ABa, Xb * (1 + kx * r) * ampa, fs, raa, rba)
                dm = rawa - rawb
                nm = slope_K * (Kmean[key] - Kmean[base])
                gross = (A8.PDS_SC * K_LIT * dpl + A8.MSE_SC * dm) / 6
                res[f"k_X {kx}"] = dict(dmse=dm, gross=gross, nmae_K=nm, net=gross + NMAE_OV * nm)
            if base == "cap3":          # known answer: 492g's gross (literal pds, k_X 0.61) for the same arm
                ref = g["arms"][arm]["calibrated"]["k_X 0.61"]["gross_lit"]
                res["ka_gross_vs_492g"] = res["k_X 0.61"]["gross"] - ref
                print(f"KA {key}: gross {res['k_X 0.61']['gross']:+.6f} vs 492g {ref:+.6f} (|d| {abs(res['ka_gross_vs_492g']):.1e})")
            kol = R.get("kolf2")
            if kol:
                res["kolf2_r"] = kol["arms"][arm_name]["realised"]["16"]["cos_agg"] / kol["arms"][base_name]["realised"]["16"]["cos_agg"] - 1
            out[key] = res
            print(f"{key:14s} (increment over {base}): r {100 * r:+.2f} % ({res['folds_pos']}/3) kolf2 {100 * res.get('kolf2_r', np.nan):+.2f} % "
                  f"dpds lit {dpl:+.5f} | k_X .61: dmse {res['k_X 0.61']['dmse']:+.5f} net {res['k_X 0.61']['net']:+.5f} | "
                  f"k_X .12: net {res['k_X 0.12']['net']:+.5f}", flush=True)
    c6 = b540["calibrated"]
    for arm in ("P2g", "S_P2"):
        tot = {kk: min(c6[f"k_X {kk}"]["net_K"], c6[f"k_X {kk}"]["net_D"]) + out[f"cap6 + {arm}"]["k_X 0.61"]["net"] for kk in ("1.31", "1")}
        out[f"total 446bj2 -> cap6 + {arm}"] = tot
        print(f"TOTAL 446bj2 -> cap 6 + {arm}: cap step (worse ballast model) + increment (k_X .61) = {tot['1.31']:+.5f} (key k_X 1.31) / "
              f"{tot['1']:+.5f} (key k 1); vs the parts alone: cap 6 {min(c6['k_X 1.31']['net_K'], c6['k_X 1.31']['net_D']):+.5f}, "
              f"{arm} on cap 3 {out[f'cap3 + {arm}']['k_X 0.61']['net']:+.5f}")
    (WORK / "analyse540c.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

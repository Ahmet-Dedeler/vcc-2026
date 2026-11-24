"""518 pricer — the strength-key arms vs 446bj2, board-calibrated exactly as 540 priced the cap (analyse540's arithmetic).

Inputs: the arms' exact PLOLO rows (skey518.json: r = cos_agg / cos_agg(446bj2's V1_0.5_3) - 1, literal d pds) and their A/B/C
energies, level maps and ballast (abc518.json). Arithmetic (analyse540 / analyse492g): 446bj2-anchored cliff model (X from the
board raw 0.6264), base and arm each re-levelled at their EV-optimal factors on their own realised energy maps;
pds = 0.34 x literal PLOLO d pds (the key class's measured pds transfer); mse at k_X 1.31 (the key class's fitted transfer,
446b -> 446bj2) and 1.0; minus the ballast's nmae cost by the K-count and displacement models (ballast492g's slopes),
-0.2808 overall per unit raw nmae.
KNOWN ANSWER: the cap-6 arm reproduces analyse540.json's recorded V1_0.5_6 dmse at k_X 1.31 (from arms540's r) when fed
arms540's r; here r comes from skey518's own exact rows, whose cap-6 r equals arms540's (asserted, 1e-12).

    uv run python experiments/518-bulk-mse-gap/analyse518.py
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
from analyse492g import CTX, RAW_446BJ2, level_opt  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
K_LIT, K_KEY_FIT = 0.34, 1.31
NMAE_OV = -0.2808
ARMS = {"cap6": ("cap6", None), "meta th1.25 on cap3": ("meta|theta1.25|on_cap3", None),
        "meta th1 on cap3": ("meta|theta1|on_cap3", None), "src th1 on cap3": ("src|theta1|on_cap3", None),
        "meta th0.5 on cap6": ("meta|theta0.5|on_cap6", "cap6"), "meta th0.75 on cap6": ("meta|theta0.75|on_cap6", "cap6")}


def main():
    cl = A8.Cliff()
    g = json.loads((ROOT / "experiments/492-per-target-sharpening/work/analyse492g.json").read_text())
    bal = json.loads((ROOT / "experiments/492-per-target-sharpening/work/ballast492g.json").read_text())
    slope_K, slope_D = bal["fits"]["dK_mean"]["slope"], bal["fits"]["D_from_design"]["slope"]
    AB0 = g["AB_446bj2"]
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    abc = json.loads((WORK / "abc518.json").read_text())
    sk = json.loads((WORK / "skey518.json").read_text())
    a540 = json.loads((ROOT / "experiments/540-lodo-clip-scan/work/analyse540.json").read_text())["arms"]["V1_0.5_6"]
    fs = [0.92, 0.95, 0.98, 1.0, 1.02, 1.05, 1.08]

    def ratio_map(tag, key):
        mp = abc["maps"][tag]
        return {c: [mp[f"{f:g}"][c][key] / mp["1"][c][key] for f in fs] for c in CTX}
    f0, raw0 = level_opt(cl, AB0, X0, fs, ratio_map("cap3", "A_own"), ratio_map("cap3", "B_own"))
    out = dict(base_level=f0, base_raw=raw0, arms={})
    # KA: skey518's cap-6 r per fold == arms540's recorded r (via analyse540.json)
    for f in CANCER:
        assert abs(sk[f]["cap6"]["r"] - a540["folds"][f]["r"]) < 1e-12, f"KA cap6 r {f}"
    print(f"KA: skey518 cap-6 r == arms540 on the cancer folds; base re-level {f0}, EV raw {raw0:.5f}")
    for tag, (skname, on) in ARMS.items():
        a = abc["arms"][tag]
        ABa = {c: dict(A=AB0[c]["A"] * a[c]["A_ratio"], B=AB0[c]["B"] * a[c]["B_ratio"]) for c in CTX}
        amp = np.sqrt(np.mean([ABa[c]["A"] for c in CTX]) / np.mean([AB0[c]["A"] for c in CTX]))
        ram, rbm = ratio_map(tag, "A_own"), ratio_map(tag, "B_own")
        rows = {}
        for f in (*CANCER, "kolf2"):
            if tag == "cap6":
                r, dp = sk[f]["cap6"]["r"], sk[f]["cap6"]["dpds"]
            else:
                arm = sk[f]["arms"][skname]
                r = arm["r_vs_446bj2"]
                dp = arm["dpds"] + (sk[f]["cap6"]["dpds"] if on == "cap6" else 0.0)
            row = dict(r=r, dpds_lit=dp)
            for k in (1.0, K_KEY_FIT):
                row[f"dmse_k{k:g}"] = level_opt(cl, ABa, X0 * (1 + k * r) * amp, fs, ram, rbm)[1] - raw0
            rows[f] = row
        mean = lambda key: float(np.mean([rows[f][key] for f in CANCER]))
        dp = K_LIT * mean("dpds_lit")
        nK, nD = slope_K * a["dK_mean"], slope_D * a["D_from_design"]
        cal = {}
        for k in (K_KEY_FIT, 1.0):
            dm = mean(f"dmse_k{k:g}")
            gross = (A8.PDS_SC * dp + A8.MSE_SC * dm) / 6
            cal[f"k_X {k:g}"] = dict(dmse_raw=dm, gross=gross, net_K=gross + NMAE_OV * nK, net_D=gross + NMAE_OV * nD)
        kol = rows["kolf2"]
        kol_gross = (A8.PDS_SC * K_LIT * kol["dpds_lit"] + A8.MSE_SC * kol["dmse_k1"]) / 6
        out["arms"][tag] = dict(folds=rows, r_cancer=mean("r"), n_pos=int(sum(rows[f]["r"] > 0 for f in CANCER)), dpds_raw=dp,
                                nmae_raw_K=nK, nmae_raw_D=nD, calibrated=cal, kolf2_gross_k1=kol_gross,
                                energy={c: a[c]["A_ratio"] for c in CTX}, dK=a["dK_mean"])
        c1, c2 = cal["k_X 1.31"], cal["k_X 1"]
        print(f"{tag:22s} r {100 * mean('r'):+.2f} % ({out['arms'][tag]['n_pos']}/3), kolf2 {100 * kol['r']:+.2f} % | pds {dp:+.5f} | "
              f"mse k1.31 {c1['dmse_raw']:+.5f} (k1 {c2['dmse_raw']:+.5f}) | nmae K {nK:+.5f} D {nD:+.5f} | OVERALL vs 446bj2 "
              f"k1.31 {c1['net_K']:+.5f}/{c1['net_D']:+.5f}  k1 {c2['net_K']:+.5f}/{c2['net_D']:+.5f} | kolf2 gross k1 {kol_gross:+.5f}")
    base6 = out["arms"]["cap6"]["calibrated"]
    for tag in ("meta th0.5 on cap6", "meta th0.75 on cap6"):
        inc = {k: {kk: out["arms"][tag]["calibrated"][k][kk] - base6[k][kk] for kk in ("net_K", "net_D")} for k in base6}
        out["arms"][tag]["increment_over_cap6"] = inc
        print(f"   increment of {tag} over cap 6: k1.31 {inc['k_X 1.31']['net_K']:+.5f}/{inc['k_X 1.31']['net_D']:+.5f}, "
              f"k1 {inc['k_X 1']['net_K']:+.5f}/{inc['k_X 1']['net_D']:+.5f}")
    (WORK / "analyse518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

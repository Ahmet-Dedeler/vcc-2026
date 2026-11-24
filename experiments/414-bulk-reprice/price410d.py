"""414: 410d as preregistered (228's gain_reg t 0.25 on the final bulk row, Frobenius-preserved in EFFECT space) vs 410c,
through 187's per-context cliff model, with the render-space energy ratio read on
  (i) the three cancer folds (fold-specific gains, the PLOLO arm `Kinf+gf@0.25`), and
  (ii) kolf2 — whose multiplier IS the shipped 410d artifact (max |d| = 0 checked here) on a production-operator,
       three-production-donor bulk: the closest analogue of the real build.
Also the amplitude-matched variant (bulk rescaled per context so its render-space energy equals 410c's).

    uv run python price410d.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from analyse414 import CTX, MSE_SC, PDS_SC, joint187, fold_rows, mse_c  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def main():
    z = np.load(ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz")
    g = np.load(ROOT / "experiments/228-theta-gain-plolo/work/gain-kolf2.npz", allow_pickle=True)
    live, gr = np.asarray(g["live"], bool), np.asarray(g["gain_reg"], float)
    m = 1.0 + 0.25 * (gr / float(np.median(gr[live])) - 1.0)
    d_art = float(np.abs(np.asarray(z["m"]) - m).max())
    print(f"410d artifact m == kolf2-fold multiplier at t 0.25: max |d| {d_art:.1e}; m q0/1/50/99/100 "
          f"{np.round(np.quantile(z['m'], [0, .01, .5, .99, 1]), 3).tolist()}; genes with m < 0: {(z['m'] < 0).sum()}")
    an = json.loads((WORK / "analyse414.json").read_text())
    dec = json.loads((WORK / "decomp414.json").read_text())
    mdl = joint187.Model(joint187.F_FIT)
    claim = {c: mdl.t[c]["claim"] for c in CTX}
    Kc = {c: mdl.t[c]["Kc"] for c in CTX}
    sumD = joint187.sumD0
    Xb, k = an["Xb"], an["k_key"]
    T = an["plolo_cancer"]
    A410b = {c: dec["410b"][c] for c in CTX}
    # 410c bulk, as in analyse414 (cancer-fold energy ratio Kinf/K3 on 410b's measured A_c, B_c)
    AB_c = {c: dict(A=A410b[c]["A"] * T["Kinf"]["Ar"][c] / T["K3"]["Ar"][c],
                    B=A410b[c]["B"] * T["Kinf"]["Br"][c] / T["K3"]["Br"][c]) for c in CTX}
    X_c = an["amplitude"]["X410c"]
    raw_c = float(np.mean([mse_c(AB_c[c]["A"], AB_c[c]["B"], Kc[c], claim[c], X_c, sumD) for c in CTX]))
    out = dict(artifact_match=d_art, raw410c=raw_c)
    ARMS = sys.argv[1:] or ["Kinf+gf@0.25"]
    for arm, label, folds in [(a, lab, fl) for a in ARMS for lab, fl in
                              (("cancer folds", CANCER), ("kolf2 (= the shipped 410d multiplier for the gain arm)", ("kolf2",)))]:
        rows = [fold_rows(f) for f in folds]
        ar = {c: float(np.mean([r[arm]["AB"][c]["A"] / r["Kinf"]["AB"][c]["A"] for r in rows])) for c in CTX}
        br = {c: float(np.mean([(r[arm]["AB"][c]["B"] if np.isfinite(r[arm]["AB"][c]["B"]) else r[arm]["AB"][c]["A"])
                                / r["Kinf"]["AB"][c]["B"] for r in rows])) for c in CTX}
        cr = float(np.mean([r[arm]["cos_agg"] / r["Kinf"]["cos_agg"] for r in rows]))
        dp = float(np.mean([r[arm]["pds"] - r["Kinf"]["pds"] for r in rows]))
        dpb = float(np.mean([r[arm]["pds_l05"] - r["Kinf"]["pds_l05"] for r in rows]))
        label = f"{arm} | {label}"
        AB_d = {c: dict(A=AB_c[c]["A"] * ar[c], B=AB_c[c]["B"] * br[c]) for c in CTX}
        # X: the calibrated relative cosine increment (k_key) times the amplitude factor
        rel_c = T["Kinf"]["rel"]
        X_d = X_c * (1 + k * ((1 + rel_c) * cr - 1)) / (1 + k * rel_c) * np.sqrt(np.mean([ar[c] for c in CTX]))
        per_d = {c: mse_c(AB_d[c]["A"], AB_d[c]["B"], Kc[c], claim[c], X_d, sumD) for c in CTX}
        raw_d = float(np.mean(list(per_d.values())))
        # amplitude-matched: per-context scale s_c = 1/sqrt(ar_c) -> render-space energy equal to 410c's
        s = {c: 1 / np.sqrt(ar[c]) for c in CTX}
        X_m = X_c * (1 + k * ((1 + rel_c) * cr - 1)) / (1 + k * rel_c)
        raw_m = float(np.mean([mse_c(AB_d[c]["A"] * s[c] ** 2, AB_d[c]["B"] * s[c] ** 2, Kc[c], claim[c], X_m, sumD)
                               for c in CTX]))
        res = dict(A_ratio=ar, B_ratio=br, cos_ratio=cr, dpds=dp, dpds_board=dpb, raw410d=raw_d, d_raw=raw_d - raw_c,
                   per_ctx_rho={c: AB_d[c]["B"] / claim[c] for c in CTX},
                   overall=(PDS_SC * dpb + MSE_SC * (raw_d - raw_c)) / 6,
                   matched_d_raw=raw_m - raw_c, matched_overall=(PDS_SC * dpb + MSE_SC * (raw_m - raw_c)) / 6)
        out[label] = res
        print(f"\n[{label}] vs 410c: render-space energy ratio A/B/C {[round(ar[c], 4) for c in CTX]}; cos ratio {cr:.4f}; "
              f"dpds literal {dp:+.5f}, board-emulated {dpb:+.5f}")
        print(f"   rho_c after arm {[round(AB_d[c]['B'] / claim[c], 4) for c in CTX]} (410c {[round(AB_c[c]['B'] / claim[c], 4) for c in CTX]})")
        print(f"   predicted raw mse 410c {raw_c:.5f} -> arm {raw_d:.5f} (d {raw_d - raw_c:+.5f}); overall {res['overall']:+.5f}")
        print(f"   AMPLITUDE-MATCHED arm (bulk x {[round(s[c], 4) for c in CTX]}): d raw {raw_m - raw_c:+.5f}; overall "
              f"{res['matched_overall']:+.5f}")
    (WORK / "price410d.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

"""518 — the board's own mse identity for 446bj2, term by term, and what the rivals' raw mse implies in the same units.

187's per-context identity (096 (a), verified against cell_eval2 to 8.9e-16):
    MSE_c - 1 = [A_c + Kc_c - min(claim_c, B_c) - 2 X] / sumD
A_c = sum_p |h_p|^2 / G (our realised rows), B_c = the centred (across-target) part, so A_c - B_c = n |mean_p h_p|^2 / G is
our MEAN-ROW energy; claim_c = sum_p min(C^_p, C_p)/G (NTC-proxy); Kc_c = the observed control's jackknife carried by every
row; X = sum_p <h_p, delta*_p>/G (one X, solved from the board raw, KA1 of ka518); sumD = sum_p |delta*_p|^2/G.
Below the cliff (B < claim, where the re-level puts every context at rho* 0.9793 of the NTC-proxy claim) the credit is B, so
    MSE_c = 1 + [(A_c - B_c) + Kc_c - 2 X] / sumD,
i.e. the spread energy is FREE up to the cliff and only the mean row costs; at the cliff the level is pinned by the truth's
own jackknife, and the member is 1.03 - 2X/sumD: at pinned energy, raw mse moves by -(2X/sumD) per unit relative change in
X = the aggregate cosine (every coordinate of every target, energy-weighted). Rivals are placed on the same identity with
OUR mean-row/Kc terms (an assumption, stated) to read the aggregate-cosine gap their raw implies.
Also: the mean-response share of the truth from the metric's own baseline (every row = the pooled perturbed profile, so
B = 0 and no credit): b_c = 1 - n|mean delta*|^2 / sum|delta*|^2 up to the Jensen gap, b in 0.986-0.992 (brief §8).

    uv run python experiments/518-bulk-mse-gap/board518.py
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
from analyse492g import CTX, RAW_446BJ2  # noqa: E402

# raw (pds_cosine, expr_mse_unbiased_capped_norm) from the board (mandate, 2026-09-28)
RIVALS = {"446bj2 (us)": (0.8674284652545521, 0.6264081539046438), "cqawesome": (0.8747, 0.5704),
          "Testh": (0.8587, 0.5723), "Illumina AI": (0.8823, 0.4899)}
B_RANGE = (0.986, 0.992)       # the metric's mse baseline anchors over A/B/C (docs/reference brief §8)
JK_BIAS = 0.0032               # a control-emitting submission reads ~1.0032 (brief §2): the jackknife is 0.32 % high


def main():
    cl = A8.Cliff()
    g = json.loads((ROOT / "experiments/492-per-target-sharpening/work/analyse492g.json").read_text())
    AB0 = g["AB_446bj2"]
    X0 = cl.solve_X(AB0, RAW_446BJ2)
    sD = cl.sumD
    rows = {}
    for c in CTX:
        A, B, Kc, claim = AB0[c]["A"], AB0[c]["B"], cl.Kc[c], cl.claim[c]
        raw = A8.mse_c(A, B, Kc, claim, X0, sD)
        rows[c] = dict(A=A, B=B, claim=claim, rho=B / claim, Kc=Kc, raw=raw,
                       A_over_sumD=A / sD, meanrow_over_sumD=(A - B) / sD, Kc_over_sumD=Kc / sD,
                       credit_over_sumD=min(claim, B) / sD, twoX_over_sumD=2 * X0 / sD,
                       cos_agg=X0 / np.sqrt(A * sD), energy_opt_scale=X0 / A,
                       check=1 + (A - B) / sD + Kc / sD - 2 * X0 / sD if B < claim else None)
    const = float(np.mean([1 + (AB0[c]["A"] - AB0[c]["B"] + cl.Kc[c]) / sD for c in CTX]))
    twoX = 2 * X0 / sD
    riv = {}
    for name, (pds, mse) in RIVALS.items():
        tx = const - mse                        # 2X/sumD they need at our constants
        riv[name] = dict(pds=pds, mse=mse, twoX_over_sumD=tx, X_ratio_vs_us=tx / twoX,
                         dmse_vs_us=mse - RAW_446BJ2, doverall_mse_only=A8.MSE_SC * (mse - RAW_446BJ2) / 6)
    mean_share = dict(lo=1 - B_RANGE[1], hi=1 - B_RANGE[0], lo_jk=1 + JK_BIAS - B_RANGE[1], hi_jk=1 + JK_BIAS - B_RANGE[0])
    ours_mean = {c: rows[c]["meanrow_over_sumD"] for c in CTX}
    # the mean-component pot: |mean h - mean delta*|^2 / sumD for alignment cos in {1, 0.5, 0, -1}
    pot = {}
    for cosm in (1.0, 0.5, 0.0, -1.0):
        vals = []
        for c in CTX:
            a, t = ours_mean[c], 0.5 * (mean_share["lo_jk"] + mean_share["hi_jk"])
            vals.append(a + t - 2 * cosm * np.sqrt(a * t))
        pot[f"cos_{cosm:+g}"] = float(np.mean(vals))
    out = dict(X=X0, sumD=sD, per_context=rows, const_1_plus_meanrow_Kc=const, rivals=riv,
               truth_mean_response_share=mean_share, mean_component_pot_raw=pot,
               d_raw_per_1pct_X=-twoX * 0.01, d_overall_per_1pct_X=A8.MSE_SC * (-twoX * 0.01) / 6)
    print(f"X {X0:.6f} sumD {sD:.6f}  2X/sumD {twoX:.4f}  const (1+(A-B+Kc)/sumD, mean) {const:.4f}")
    print(f"{'ctx':4}{'A/sD':>8}{'meanrow':>9}{'Kc/sD':>8}{'credit':>8}{'rho':>7}{'2X/sD':>8}{'cos_agg':>9}{'X/A':>7}{'raw':>9}")
    for c in CTX:
        r = rows[c]
        print(f"{c:4}{r['A_over_sumD']:8.4f}{r['meanrow_over_sumD']:9.4f}{r['Kc_over_sumD']:8.4f}{r['credit_over_sumD']:8.4f}"
              f"{r['rho']:7.4f}{r['twoX_over_sumD']:8.4f}{r['cos_agg']:9.4f}{r['energy_opt_scale']:7.3f}{r['raw']:9.5f}")
    print(f"at the cliff: d raw per +1 % X = {out['d_raw_per_1pct_X']:+.5f} ({out['d_overall_per_1pct_X']:+.5f} overall)")
    for n, v in riv.items():
        print(f"{n:14s} pds {v['pds']:.4f} mse {v['mse']:.4f}: 2X/sumD {v['twoX_over_sumD']:.4f} = {v['X_ratio_vs_us']:.3f}x ours; "
              f"mse-only overall vs us {v['doverall_mse_only']:+.4f}")
    print(f"truth mean-response share (1 - b): {mean_share['lo']:.3f}-{mean_share['hi']:.3f} "
          f"({mean_share['lo_jk']:.4f}-{mean_share['hi_jk']:.4f} with the jackknife's +0.32 %); our mean row / sumD "
          + " ".join(f"{c} {ours_mean[c]:.4f}" for c in CTX))
    print("mean-component pot (raw mse at an ORACLE mean) by alignment of our mean with the truth mean: "
          + " ".join(f"{k} {v:.4f}" for k, v in pot.items()))
    (WORK / "board518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

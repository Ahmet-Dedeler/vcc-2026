"""518 — where the aggregate cosine is lost on the PLOLO folds: per-target magnitude allocation vs direction, from 420's saved
per-target arrays (no truth beyond PLOLO's; KA3 of ka518 rebuilds comp420's cos_agg from them exactly).

Arm: 446bj2's bulk on the folds = `stack|V1_5_3+sh1.75` levelled at A 1.101 / B 1.155 / C 1.183 (420's level run; 446bj2
shipped 1.0962 / 1.1461 / 1.1681), realised by the ballast at KMAX 16. Per target p and frame c (the challenge context whose
control profile and gate render the row; the truth is the held line's row, shared by the three frames):
    x_pc = <h_pc, t_p>,  a_pc = |h_pc|^2,  d_p = |t_p|^2   over the fold's scored (non-panel) genes.
At the board's cliff the spread energy is pinned, so the mse reads X per unit energy = cos_agg = sum x / sqrt(sum a sum d)
(board518: d raw = -(2X/sumD) r, 2X/sumD = 0.406).
    s_p = sqrt(a_p) our magnitude,  m_p = x_p / sqrt(a_p) the truth's projection on our unit direction.
MAGNITUDE ORACLE (nonnegative per-target scale, energy re-pinned): cos* = sqrt(sum m+^2 / sum d); r_mag = cos*/cos - 1.
    Split: r_neg = zero only the targets with m < 0; the rest is reallocation among m > 0.
    Additive attribution of the oracle's X gain at pinned energy (lambda = sqrt(sum a / sum m+^2)):
        X* - X = sum_p [ (s_p - lambda m+_p)^2 / (2 lambda) + s_p max(-m_p, 0) ]  (each term >= 0),
    reported by decile of OUR energy (a), of the truth's energy (d, noise-inflated on the folds), of the LODO key, and by
    over- (s > lambda m) vs under-scaled.
Noise: the fold truth is a donor pseudobulk; E[x^2/a] carries one dimension of truth noise per target. Its ceiling here
(all truth energy treated as noise spread over the scored genes) is reported as `noise_bias_ceiling` next to the oracle.

    uv run python experiments/518-bulk-mse-gap/folds518.py
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
W420 = ROOT / "experiments/420-bulk-stack/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
ARM = "stack|V1_5_3+sh1.75:A=1.101,B=1.155,C=1.183@R16"
TWO_X = 0.40587             # board518: 2X/sumD at 446bj2
NSCORED = {"k562_gwps": None}


def deciles(v, n=10):
    q = np.quantile(v, np.linspace(0, 1, n + 1))
    idx = np.clip(np.searchsorted(q, v, side="right") - 1, 0, n - 1)
    return idx


def frame_decomp(x, a, d):
    s = np.sqrt(a)
    m = x / np.maximum(s, 1e-300)
    mp = np.maximum(m, 0.0)
    D = d.sum()
    cos = x.sum() / np.sqrt(a.sum() * D)
    cos_star = np.sqrt((mp ** 2).sum() / D)
    keep = m > 0
    cos_neg = x[keep].sum() / np.sqrt(a[keep].sum() * D)        # zero the m<0 targets, re-pin energy
    lam = np.sqrt(a.sum() / (mp ** 2).sum())
    Xnow, Xstar = (s * m).sum(), np.sqrt(a.sum() * (mp ** 2).sum())
    loss = (s - lam * mp) ** 2 / (2 * lam) + s * np.maximum(-m, 0.0)     # + the X an anti-aligned target costs
    assert abs(loss.sum() - (Xstar - Xnow)) < 1e-9 * abs(Xstar), "attribution identity"
    return dict(cos=cos, cos_star=cos_star, r_mag=cos_star / cos - 1, r_neg=cos_neg / cos - 1, frac_neg=float((m < 0).mean()),
                x_neg_share=float(-x[m < 0].sum() / x[m > 0].sum()), loss=loss, s=s, m=m, lam=lam, over=s > lam * mp)


def main():
    out = dict(arm=ARM, folds={})
    for f in FOLDS:
        z = np.load(W420 / f"arrays420L-{f}.npz")
        x, a, d = z[f"{ARM}__x"], z[f"{ARM}__a"], z["d"]
        n = d.size
        comp = json.loads((W420 / f"comp420-{f}.json").read_text())
        w = np.maximum(np.load(W418 / f"arrays418-{f}.npz")["V1__PP"], 0.5)          # the fold's V1_5_3 key
        rec = dict(n=int(n), frames={})
        tot_loss = np.zeros(n)
        for j, c in enumerate("ABC"):
            fd = frame_decomp(x[j], a[j], d)
            tot_loss += fd["loss"] / fd["loss"].sum()
            rec["frames"][c] = {k: float(v) for k, v in fd.items() if np.isscalar(v)}
            if j == 0:
                F0 = fd
        # frame-mean summaries
        for k in ("cos", "cos_star", "r_mag", "r_neg", "frac_neg", "x_neg_share"):
            rec[k] = float(np.mean([rec["frames"][c][k] for c in "ABC"]))
        # noise ceiling: all truth energy as isotropic noise over the scored genes -> bias sum_p d_p/G in sum m^2
        G = 17942 - 0          # scored columns are ~17.9k non-panel genes on every fold (fold audit); ceiling only
        mp2 = np.mean([(np.maximum(x[j] / np.sqrt(a[j]), 0) ** 2).sum() for j in range(3)])
        rec["noise_bias_ceiling_frac_of_sum_m2"] = float(d.sum() / G / mp2)
        # attribution tables (frame-mean of each frame's normalised loss shares)
        share = tot_loss / 3.0
        a_mean, x_mean = a.mean(0), x.mean(0)
        m_mean = x_mean / np.sqrt(a_mean)
        cosp = x_mean / np.sqrt(a_mean * d)
        tabs = {}
        for name, key in (("our_energy", a_mean), ("truth_energy", d), ("truth_proj_m", m_mean), ("key_w", w + 1e-9 * a_mean),
                          ("cos_p", cosp)):
            dec = deciles(key)
            tabs[name] = [dict(decile=int(k), n=int((dec == k).sum()), loss_share=float(share[dec == k].sum()),
                               x_share=float(x_mean[dec == k].sum() / x_mean.sum()), a_share=float(a_mean[dec == k].sum() / a_mean.sum()),
                               d_share=float(d[dec == k].sum() / d.sum()),
                               over_frac=float(np.mean(F0["over"][dec == k])),
                               ratio_med=float(np.median((F0["lam"] * np.maximum(F0["m"], 0) / F0["s"])[dec == k])))
                          for k in range(10)]
        rec["tables"] = tabs
        rec["over_loss_share"] = float(share[F0["over"]].sum())
        # the oracle's ideal relative scale (lambda m+ / s), frame A: how far from 1, and its relation to the key
        ratio = F0["lam"] * np.maximum(F0["m"], 0) / F0["s"]
        rec["ideal_ratio_q"] = np.quantile(ratio, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()
        rho = lambda u, v: float(np.corrcoef(np.argsort(np.argsort(u)), np.argsort(np.argsort(v)))[0, 1])
        rec["spearman_ideal_vs"] = dict(key_w=rho(ratio, w), our_energy=rho(ratio, a_mean), truth_energy=rho(ratio, d),
                                        m=rho(ratio, m_mean))
        rec["pds_lit"] = comp["arms"]["stack|V1_5_3+sh1.75"]["realised"]["16"]["pds"]
        out["folds"][f] = rec
        print(f"{f:14s} n {n} cos {rec['cos']:.5f} -> oracle {rec['cos_star']:.5f}: r_mag {100 * rec['r_mag']:+.1f} % "
              f"(zero m<0 only {100 * rec['r_neg']:+.1f} %; m<0 on {100 * rec['frac_neg']:.0f} % of targets) | "
              f"over-scaled share of the loss {100 * rec['over_loss_share']:.0f} % | noise ceiling {100 * rec['noise_bias_ceiling_frac_of_sum_m2']:.1f} % "
              f"| board d raw at k 1: {-TWO_X * rec['r_mag']:+.4f}")
        print("   ideal/current scale quantiles (10/25/50/75/90):", " ".join(f"{v:.2f}" for v in rec["ideal_ratio_q"]),
              "| spearman(ideal, key)", f"{rec['spearman_ideal_vs']['key_w']:+.2f}", "(ideal, our energy)",
              f"{rec['spearman_ideal_vs']['our_energy']:+.2f}")
        for name in ("our_energy", "key_w", "cos_p"):
            t = tabs[name]
            print(f"   by {name:12s} decile: loss " + " ".join(f"{100 * r['loss_share']:4.0f}" for r in t)
                  + " | X " + " ".join(f"{100 * r['x_share']:4.0f}" for r in t) + " | a " + " ".join(f"{100 * r['a_share']:3.0f}" for r in t)
                  + " | median ideal/cur " + " ".join(f"{r['ratio_med']:.2f}" for r in t))
    canc = [out["folds"][f] for f in CANCER]
    out["cancer_mean"] = {k: float(np.mean([r[k] for r in canc])) for k in ("cos", "cos_star", "r_mag", "r_neg", "frac_neg")}
    cm = out["cancer_mean"]
    print(f"CANCER MEAN r_mag {100 * cm['r_mag']:+.1f} % (zero-negatives part {100 * cm['r_neg']:+.1f} %); board d raw "
          f"k 1 {-TWO_X * cm['r_mag']:+.4f}, k 1.31 {-TWO_X * 1.31 * cm['r_mag']:+.4f}; kolf2 r_mag {100 * out['folds']['kolf2']['r_mag']:+.1f} %")
    (WORK / "folds518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

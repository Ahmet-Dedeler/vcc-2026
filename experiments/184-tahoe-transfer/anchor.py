"""184 Phase 2a -- the anchor rows: reliability, transfer level, decoy level, the consensus law.

These are the Tahoe analogues of every number 159 measured on our CRISPRi donors.  They are the
calibration for Phase 2b: if Tahoe's cross-line structure does not look like our donors', the
transfer operator learned on it cannot be expected to carry over.

Also carries the internal known-answer rows (rule 10):
  * the ORACLE prediction (the held line's own truth) must score cosine 1 and sign accuracy 1;
  * a zero prediction must score sign accuracy 0.5;
  * the different-drug cross-line cosine is the decoy level and must sit near 0 after centring
    (if it does not, the surface carries a blob and nothing else here is readable -- 159 sec.2).

Writes work/anchor.json only.
"""
from __future__ import annotations

import json
import os

import numpy as np
import polars as pl

from common import HERE, SEED, WORK, load, rowcos, sign_acc, unit

REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
TAHOE = os.path.join(REPO, "data/external/tahoe")


def plate_means(genes, lines):
    """mean delta over the drugs on each (line, plate) -- the control-noise blob per plate."""
    return np.load(os.path.join(WORK, "platemean.npy"))


def main():
    rng = np.random.default_rng(SEED)
    R, D, genes, lines, treats, basal, group = load("plate")
    T, L, G = R.shape
    out = {"n_treatments": T, "n_lines": L, "n_genes": G}

    # ---------- 0. known-answer rows on the endpoint code itself
    sub = R[:200, 0, :]
    out["known"] = {
        "oracle_cos": float(np.nanmean(rowcos(sub, sub))),
        "oracle_sign_top200": float(np.nanmean(sign_acc(sub, sub, 200))),
        "zero_sign_top200": float(np.nanmean(sign_acc(np.zeros_like(sub), sub, 200))),
        "flip_sign_top200": float(np.nanmean(sign_acc(-sub, sub, 200))),
    }
    print("known-answer:", json.dumps(out["known"]))
    assert abs(out["known"]["oracle_cos"] - 1) < 1e-5
    assert abs(out["known"]["oracle_sign_top200"] - 1) < 1e-9
    assert abs(out["known"]["zero_sign_top200"] - 0.5) < 1e-9

    # ---------- 1. same-line replicate reliability (the attenuation)
    z = np.load(os.path.join(WORK, "reps.npz"))
    PM = plate_means(genes, lines)
    ti, lj, RV, RP = z["t_idx"], z["l_idx"], z["delta"], z["plate"]
    cos_rep, cos_rep_uncentred = [], []
    for i in range(len(ti)):
        ps = [j for j in range(RP.shape[1]) if RP[i, j] > 0]
        if len(ps) < 2:
            continue
        a, b = ps[0], ps[1]
        va, vb = RV[i, a].astype(np.float64), RV[i, b].astype(np.float64)
        cos_rep_uncentred.append(float(unit(va) @ unit(vb)))
        ca = va - PM[RP[i, a], lj[i]]
        cb = vb - PM[RP[i, b], lj[i]]
        cos_rep.append(float(unit(ca) @ unit(cb)))
    out["replicate"] = {
        "n_pairs": len(cos_rep),
        "same_line_same_drug_cos_centred": float(np.mean(cos_rep)),
        "median": float(np.median(cos_rep)),
        "same_line_same_drug_cos_uncentred": float(np.mean(cos_rep_uncentred)),
    }
    rho = out["replicate"]["same_line_same_drug_cos_centred"]
    print(f"replicate reliability (same line, same drug, two plates): centred {rho:+.4f}  "
          f"uncentred {out['replicate']['same_line_same_drug_cos_uncentred']:+.4f}  "
          f"n={len(cos_rep)}")

    # ---------- 2. cross-line same-drug vs different-drug, on the full grid
    ok = np.isfinite(R[:, :, 0])
    idx = rng.permutation(T)[:400]
    same, diff, within = [], [], []
    pairs = [(a, b) for a in range(L) for b in range(a + 1, L)]
    rng.shuffle(pairs)
    for a, b in pairs[:200]:
        m = ok[idx, a] & ok[idx, b]
        A, B = unit(R[idx[m], a, :]), unit(R[idx[m], b, :])
        S = A @ B.T
        n = S.shape[0]
        same.append(float(np.diag(S).mean()))
        diff.append(float(S[~np.eye(n, dtype=bool)].mean()))
    for a in range(L):
        m = ok[idx, a]
        A = unit(R[idx[m], a, :])
        S = A @ A.T
        n = S.shape[0]
        within.append(float(S[~np.eye(n, dtype=bool)].mean()))
    out["crossline"] = {
        "n_line_pairs": 200, "n_drugs": 400,
        "same_drug": float(np.mean(same)), "diff_drug": float(np.mean(diff)),
        "excess": float(np.mean(same) - np.mean(diff)),
        "ratio": float(np.mean(same) / np.mean(diff)) if np.mean(diff) else float("nan"),
        "within_line_diff_drug": float(np.mean(within)),
        "disattenuated_same_drug": float(np.mean(same) / rho),
    }
    print(f"cross-line: same {out['crossline']['same_drug']:+.4f}  "
          f"diff {out['crossline']['diff_drug']:+.4f}  "
          f"within-line diff-drug {out['crossline']['within_line_diff_drug']:+.4f}  "
          f"disattenuated same {out['crossline']['disattenuated_same_drug']:+.4f}")

    # ---------- 3. consensus-of-k law, 1/cos^2(k) = A + B/k
    ks = [1, 2, 3, 5, 8, 12, 20, 30, 49]
    curve = {}
    for k in ks:
        vals = []
        for held in range(L):
            others = [l for l in range(L) if l != held]
            for _ in range(4 if k < 49 else 1):
                sel = rng.choice(others, size=k, replace=False)
                m = ok[:, held] & ok[:, sel].all(1)
                Cm = R[:, sel, :][m].mean(1)
                vals.append(float(np.nanmean(rowcos(Cm, R[m, held, :]))))
        curve[k] = float(np.mean(vals))
        print(f"  consensus k={k:2d}: cos {curve[k]:.5f}")
    kk = np.asarray(ks, dtype=float)
    y = np.asarray([1.0 / curve[k] ** 2 for k in ks])
    Amat = np.stack([np.ones_like(kk), 1.0 / kk], 1)
    coef, *_ = np.linalg.lstsq(Amat, y, rcond=None)
    pred = Amat @ coef
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    out["consensus"] = {
        "curve": curve, "A": float(coef[0]), "B": float(coef[1]), "R2": float(r2),
        "cos_inf": float(1 / np.sqrt(coef[0])),
        "cos_inf_disattenuated": float(1 / np.sqrt(coef[0]) / np.sqrt(rho)),
        "lines_for_half_ceiling": float(coef[1] / (3 * coef[0])),
    }
    print(f"1/cos^2 = {coef[0]:.2f} + {coef[1]:.2f}/k   R2 {r2:.5f}   "
          f"cos_inf {out['consensus']['cos_inf']:.4f}  "
          f"(disattenuated {out['consensus']['cos_inf_disattenuated']:.4f})  "
          f"lines for half the ceiling {out['consensus']['lines_for_half_ceiling']:.1f}")

    # kappa^2: the fraction of a drug's drug-specific energy that is line-independent
    tau = out["crossline"]["disattenuated_same_drug"]
    out["tau"] = {"from_pair_disattenuation": tau,
                  "from_consensus": float(out["consensus"]["cos_inf_disattenuated"] ** 2)}
    print(f"tau (shared fraction of drug-specific response energy): "
          f"pair route {tau:.3f}   consensus route {out['tau']['from_consensus']:.3f}")

    with open(os.path.join(WORK, "anchor.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

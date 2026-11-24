#!/usr/bin/env python3
"""627 — Is identity erosion UNIVERSAL across perturbation modality? Tahoe-100M: 50 cancer lines x 1,137 drug treatments.

Claim tested (from CRISPRi, 622/639): under perturbation, genes a line expresses above other lines go DOWN.
Data: 184's tensor (tahoebio/tahoe-de-rhaister): delta[d, L, g] = mean(treated) - mean(plate-matched DMSO), library-
normalised, on 2,000 genes; basal_plate[L, p, g] = the line's DMSO control as a library fraction per plate; plate_of[d, L] =
bitmask of the plates (d, L) was measured on.
Shared-noise trap: delta subtracts the plate's DMSO mean, and an identity axis built from the SAME controls would carry the
same sampling noise with the opposite sign, i.e. a spurious "identity genes go down". So the identity axis of every
(d, L) is built ONLY from plates NOT in plate_of[d, L] (leave-plate-out), for all 50 lines alike.
  I_L(d) = log1p(1e6 * basal_L, other plates) - mean over the other 49 lines (same plates), genes expressed in L.
Tests:
  A (line-specific, drug signature removed): r[d, L] = delta[d, L] - mean_L' delta[d, L'];
     stat = Spearman(r[d, L], -I_L(d)) over genes; and the sign agreement of r with -I on the top-100 |r| genes.
     Null: the same with a RANDOM OTHER line's identity (I_L'(d)).
  B (generic, per line and plate): m[L, p] = mean over treatments measured only on plate p of delta[d, L];
     stat = Spearman(m[L, p], -I_L(not p)); null: another line's identity.
  Strata: response strength |delta[d, L]| quintiles (cytotoxicity confound) and dose (the treatment's concentration).
Analysis only; < 4 GB.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W184 = ROOT / "experiments/184-tahoe-transfer/work"
RNG = np.random.default_rng(6270)


def spearman_rows(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Row-wise Spearman between A and B (same shape), ignoring NaN pairs per row."""
    out = np.full(A.shape[0], np.nan)
    for i in range(A.shape[0]):
        ok = np.isfinite(A[i]) & np.isfinite(B[i])
        if ok.sum() > 50:
            a = rankdata(A[i, ok]); b = rankdata(B[i, ok])
            out[i] = np.corrcoef(a, b)[0, 1]
    return out


def main() -> None:
    z = np.load(W184 / "tensor.npz", allow_pickle=True)
    delta = z["delta"].astype(np.float64)                 # T x L x G
    basal = z["basal_plate"].astype(np.float64)           # L x P x G (library fractions)
    cells = z["ctrl_cells"]                               # L x P
    treat = z["treatments"].astype(str)
    plate_of = np.load(W184 / "plate_of.npy", allow_pickle=True).astype(np.int64)   # T x L bitmask
    T, L, G = delta.shape
    P = basal.shape[1]
    valid = cells > 0
    # identity per plate-exclusion mask (cache)
    cacheI: dict[int, np.ndarray] = {}

    def identity(mask: int) -> np.ndarray:
        if mask not in cacheI:
            keep = np.array([not (mask >> p) & 1 for p in range(P)])
            num = np.zeros((L, G)); den = np.zeros((L, 1))
            for p in range(P):
                if not keep[p]:
                    continue
                w = valid[:, p][:, None].astype(float)
                num += np.nan_to_num(basal[:, p, :]) * w; den += w
            frac = num / np.maximum(den, 1)
            lg = np.log1p(1e6 * frac)
            I = np.empty_like(lg)
            for l in range(L):
                I[l] = lg[l] - np.delete(lg, l, axis=0).mean(0)
            expressed = 1e6 * frac > 5
            cacheI[mask] = np.where(expressed, I, np.nan)
        return cacheI[mask]

    dose = np.array([float(m.group(1)) if (m := re.search(r",\s*([0-9.]+)\s*,\s*'uM'", t)) else np.nan for t in treat])
    # ---- A: line-specific deviation vs -I(other plates)
    stats, nulls, signagree, strength, dosev = [], [], [], [], []
    for d in range(T):
        dd = delta[d]                                      # L x G
        ok_lines = np.isfinite(dd).all(1) & (np.abs(dd).sum(1) > 0)
        if ok_lines.sum() < 20:
            continue
        r = dd - np.nanmean(dd[ok_lines], 0)
        for l in np.flatnonzero(ok_lines):
            I = identity(int(plate_of[d, l]))
            other = RNG.choice([x for x in range(L) if x != l])
            s = spearman_rows(r[l][None, :], -I[l][None, :])[0]
            n = spearman_rows(r[l][None, :], -I[other][None, :])[0]
            ok = np.isfinite(r[l]) & np.isfinite(I[l])
            top = np.argsort(-np.abs(np.where(ok, r[l], 0)))[:100]
            sa = float(np.mean(np.sign(r[l][top]) == np.sign(-I[l][top])))
            stats.append(s); nulls.append(n); signagree.append(sa)
            strength.append(float(np.linalg.norm(dd[l]))); dosev.append(dose[d])
    stats, nulls, signagree = np.array(stats), np.array(nulls), np.array(signagree)
    strength, dosev = np.array(strength), np.array(dosev)
    A = dict(n_pairs=int(np.isfinite(stats).sum()), mean_rho=float(np.nanmean(stats)), frac_rho_pos=float(np.nanmean(stats > 0)),
             mean_rho_null=float(np.nanmean(nulls)), frac_null_pos=float(np.nanmean(nulls > 0)),
             top100_sign_agree=float(np.nanmean(signagree)))
    q = np.nanquantile(strength, [0.2, 0.4, 0.6, 0.8])
    A["by_strength_quintile"] = []
    for k in range(5):
        sel = np.digitize(strength, q) == k
        A["by_strength_quintile"].append(dict(q=k + 1, n=int(sel.sum()), mean_rho=float(np.nanmean(stats[sel])),
                                               null=float(np.nanmean(nulls[sel])), top100=float(np.nanmean(signagree[sel]))))
    A["by_dose"] = {}
    for dv in sorted(set(dosev[np.isfinite(dosev)])):
        sel = dosev == dv
        A["by_dose"][str(dv)] = dict(n=int(sel.sum()), mean_rho=float(np.nanmean(stats[sel])), top100=float(np.nanmean(signagree[sel])))
    print(f"[A] line-specific deviation vs -I (other plates): pairs {A['n_pairs']}, mean Spearman {A['mean_rho']:+.4f} "
          f"(positive in {A['frac_rho_pos']:.3f}); NULL (other line's I) {A['mean_rho_null']:+.4f} ({A['frac_null_pos']:.3f}); "
          f"top-100 |r| sign agreement with -I {A['top100_sign_agree']:.4f}", flush=True)
    for row in A["by_strength_quintile"]:
        print(f"    strength Q{row['q']}: rho {row['mean_rho']:+.4f} (null {row['null']:+.4f}), top-100 sign {row['top100']:.4f}, n {row['n']}")
    for dv, row in A["by_dose"].items():
        print(f"    dose {dv} uM: rho {row['mean_rho']:+.4f}, top-100 sign {row['top100']:.4f}, n {row['n']}")
    # ---- B: generic per (line, plate) mean response vs -I(not that plate)
    Bs, Bn = [], []
    for p in range(P):
        onlyp = np.array([(plate_of[d] == (1 << p)).all() for d in range(T)])
        if onlyp.sum() < 10:
            continue
        m = np.nanmean(delta[onlyp], 0)                    # L x G
        I = identity(1 << p)
        for l in range(L):
            if not valid[l, p]:
                continue
            other = RNG.choice([x for x in range(L) if x != l])
            Bs.append(spearman_rows(m[l][None, :], -I[l][None, :])[0]); Bn.append(spearman_rows(m[l][None, :], -I[other][None, :])[0])
    Bs, Bn = np.array(Bs), np.array(Bn)
    Bres = dict(n=int(np.isfinite(Bs).sum()), mean_rho=float(np.nanmean(Bs)), frac_pos=float(np.nanmean(Bs > 0)),
                null=float(np.nanmean(Bn)), null_frac_pos=float(np.nanmean(Bn > 0)))
    print(f"[B] generic (line x plate mean response) vs -I(not that plate): n {Bres['n']}, mean Spearman {Bres['mean_rho']:+.4f} "
          f"(positive {Bres['frac_pos']:.3f}); NULL {Bres['null']:+.4f} ({Bres['null_frac_pos']:.3f})", flush=True)
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/tahoe627.json").write_text(json.dumps({"A": A, "B": Bres}, indent=1))


if __name__ == "__main__":
    main()

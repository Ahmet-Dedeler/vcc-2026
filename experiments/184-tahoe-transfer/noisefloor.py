"""184 Phase 2e -- is the CRISPRi null a NOISE FLOOR or a real absence?

Tahoe's pseudobulks are 3-5x more reliable than our donor tables (same-line replicate cosine
0.281 here against a residual split-half reliability of 0.05-0.12 for the Orion donors at 150
cells, 159 sec.1.3).  So "Tahoe shows the law, CRISPRi does not" has two readings, and they lead
to opposite decisions:

  A  chemical-vs-genetic: drug responses carry a recipient-conditional per-gene amplitude that
     knockdown responses do not.  Then no amount of extra CRISPRi depth helps.
  B  depth: the same law is there in CRISPRi and is buried under measurement noise.  Then deeper
     references (the challenge's own 400-cell references, or pooled donors) would reveal it.

Discriminated by degrading Tahoe to CRISPRi-like reliability and re-measuring.  With observed
energy V split as signal S = rho*V and noise (1-rho)*V, adding independent noise of energy
A = S*(1/rho' - 1/rho) takes the reliability from rho to rho'.  The added noise keeps each gene's
own residual sd, so only the SNR changes, not the structure of the problem.

Writes work/noisefloor.json only.
"""
from __future__ import annotations

import json
import os

import numpy as np

from common import SEED, WORK, load, unit
from phase2 import cos_rows, fit_gain, gain_features, orth, sign_acc_idx, top_index

RHO0 = 0.2808          # measured same-line, same-drug, cross-plate cosine (anchor.py)


def measure_rho(scale, rng):
    """Check the achieved reliability by adding two independent noise draws to the replicate
    store's two plate measurements of the same (line, drug)."""
    z = np.load(os.path.join(WORK, "reps.npz"))
    PM = np.load(os.path.join(WORK, "platemean.npy"))
    ti, lj, RV, RP = z["t_idx"], z["l_idx"], z["delta"], z["plate"]
    out = []
    for i in range(0, len(ti), 4):                      # every 4th pair, enough for 3 decimals
        ps = [j for j in range(RP.shape[1]) if RP[i, j] > 0]
        if len(ps) < 2:
            continue
        a, b = ps[0], ps[1]
        va = RV[i, a] - PM[RP[i, a], lj[i]]
        vb = RV[i, b] - PM[RP[i, b], lj[i]]
        s = scale * np.sqrt(0.5 * (va.var() + vb.var()))
        va = va + rng.normal(0, s, va.shape)
        vb = vb + rng.normal(0, s, vb.shape)
        out.append(float(unit(va.astype(np.float64)) @ unit(vb.astype(np.float64))))
    return float(np.mean(out))


def main():
    rng = np.random.default_rng(SEED)
    R, D, genes, lines, treats, basal, mask = load("plate")
    T, L, G = R.shape
    ok = np.isfinite(R[:, :, 0])
    R = np.nan_to_num(R, nan=0.0)
    logb = np.log10(basal + 1e-8)

    out = {"rho0": RHO0, "levels": {}}
    for rho1 in [RHO0, 0.12, 0.08, 0.05]:
        # energy multiplier for the added noise, relative to the observed residual energy
        amp = 0.0 if rho1 >= RHO0 else np.sqrt(RHO0 * (1.0 / rho1 - 1.0 / RHO0))
        ach = measure_rho(amp, np.random.default_rng(SEED + 1)) if amp else RHO0
        RN = R.copy()
        if amp:
            sd = R.std(0, keepdims=True)                # per (line, gene) residual sd
            RN = RN + rng.normal(0, 1, RN.shape).astype(np.float32) * (amp * sd)
        lvl = {"rho_target": rho1, "noise_amp": float(amp), "rho_achieved": ach, "k": {}}
        for k in [5, 49]:
            acc = {a: [] for a in ["gain-mean", "gain-basal", "gain-basal-clip",
                                   "gain-basal/perm-line", "gain-ORACLE"]}
            sg = {a: [] for a in acc}
            sgM, cosM = [], []
            for held in range(L):
                others = [l for l in range(L) if l != held]
                S = np.asarray(others if k >= L - 1 else rng.choice(others, k, replace=False))
                truth = RN[:, held, :]
                idx = top_index(truth, 200)
                M = RN[:, S, :].mean(1)
                mu, sd_ = logb[S].mean(0), logb[S].std(0) + 1e-6
                Z = (logb - mu) / sd_
                gh = np.asarray([fit_gain(RN[:, l, :], RN[:, np.delete(S, i), :].mean(1))
                                 for i, l in enumerate(S)])
                gbar = gh.mean(0)
                X = np.concatenate([gain_features(gbar, Z[l]) for l in S])
                y = np.concatenate(list(gh))
                cf = np.linalg.lstsq(X[:, :2], y, rcond=None)[0]
                cx = np.linalg.lstsq(X, y, rcond=None)[0]
                gb = gain_features(gbar, Z[held]) @ cx
                zn = Z[rng.choice([l for l in range(L) if l != held])]
                preds = {"gain-mean": M * (gain_features(gbar, Z[held])[:, :2] @ cf),
                         "gain-basal": M * gb,
                         "gain-basal-clip": M * np.maximum(gb, 0),
                         "gain-basal/perm-line": M * (gain_features(gbar, zn) @ cx),
                         "gain-ORACLE": M * fit_gain(truth, M)}
                sgM.append(float(np.nanmean(sign_acc_idx(M, truth, idx))))
                cosM.append(float(np.nanmean(cos_rows(M, truth))))
                for a, P in preds.items():
                    acc[a].append(float(np.nanmean(cos_rows(orth(P, M), orth(truth, M)))))
                    sg[a].append(float(np.nanmean(sign_acc_idx(P, truth, idx))))
            lvl["k"][str(k)] = {"M_sign": float(np.mean(sgM)), "M_cos": float(np.mean(cosM)),
                                "cos_orth": {a: float(np.mean(v)) for a, v in acc.items()},
                                "sign": {a: float(np.mean(v)) for a, v in sg.items()}}
            c = lvl["k"][str(k)]["cos_orth"]
            print(f"rho {rho1:.3f} (achieved {ach:.3f}, amp {amp:.2f})  k={k:2d}  "
                  f"M cos {lvl['k'][str(k)]['M_cos']:.3f}  "
                  f"basal-clip {c['gain-basal-clip']:+.4f}  "
                  f"ctx-free {c['gain-mean']:+.4f}  null {c['gain-basal/perm-line']:+.4f}  "
                  f"ORACLE {c['gain-ORACLE']:+.4f}")
        out["levels"][str(rho1)] = lvl
        del RN

    with open(os.path.join(WORK, "noisefloor.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

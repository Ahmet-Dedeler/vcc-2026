#!/usr/bin/env python3
"""622k — independent review of [4d0207]'s 604 artifact (ident604_flex_lam0.5_ABC.npz): recompute every context's
identity-augmented P with 622's own code path and compare gene by gene.
(1) the pinned A/B/C control CPM (604's flex_refs604.npz) vs 622's abc_cpm622.npz;
(2) I_c with references = the other two contexts + H1 Flex (604's rule: never the context's own call);
(3) S = z(P_C4 - .5) + 0.5 z(I) (missing part = 0; each z over the genes where it is defined);
(4) quantile map onto CLEAN4's P: covered genes by rank, identity-only genes by monotone interpolation;
compare with 604's pdown_<c>: sign agreement (P > .5), Spearman, max |dP| on covered genes. Analysis only, < 1 GB.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
A604 = np.load(ROOT / "experiments/604-identity-prior/work/ident604_flex_lam0.5_ABC.npz", allow_pickle=True)
R604 = np.load(ROOT / "experiments/604-identity-prior/work/flex_refs604.npz", allow_pickle=True)
mine = np.load(HERE / "work/abc_cpm622.npz", allow_pickle=True)
genes = A604["genes"].astype(str)
assert list(mine["genes"].astype(str)) == list(genes) == list(R604["genes"].astype(str))
out = {}
for c in "ABC":
    d = np.nanmax(np.abs(R604[f"cpm_{c}"] - mine[f"ctx{c}"]))
    out[f"cpm_{c}_maxabs_vs_622"] = float(d)
    print(f"(1) cpm_{c}: max |604 - 622| = {d:.3e}")
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
assert list(c4["genes"].astype(str)) == list(genes)
P4 = c4["P_down"].astype(float)
L = {c: np.log1p(R604[f"cpm_{c}"]) for c in "ABC"}
L["h1"] = np.log1p(R604["cpm_h1"])
for c in "ABC":
    refs = [o for o in "ABC" if o != c] + ["h1"]
    ok = np.isfinite(L[c]) & np.all([np.isfinite(L[r]) for r in refs], 0)
    I = np.where(ok, L[c] - np.mean([L[r] for r in refs], 0), np.nan)
    zc = np.where(np.isfinite(P4), (P4 - 0.5 - np.nanmean(P4 - 0.5)) / np.nanstd(P4 - 0.5), np.nan)
    zi = np.where(np.isfinite(I), (I - np.nanmean(I)) / np.nanstd(I), np.nan)
    S = np.nan_to_num(zc) + 0.5 * np.nan_to_num(zi)
    S = np.where(np.isfinite(zc) | np.isfinite(zi), S, np.nan)
    cov = np.isfinite(P4) & np.isfinite(S)
    Ps = np.sort(P4[cov]); order = np.argsort(np.argsort(S[cov], kind="stable"), kind="stable")
    P = np.full(len(genes), np.nan); P[cov] = Ps[order]
    only = ~np.isfinite(P4) & np.isfinite(S)
    Ss = np.sort(S[cov]); P[only] = np.interp(S[only], Ss, Ps)
    Q = A604[f"pdown_{c}"].astype(float)
    both = np.isfinite(P) & np.isfinite(Q)
    sign_agree = float(((P[both] > 0.5) == (Q[both] > 0.5)).mean())
    rho = float(spearmanr(P[both], Q[both]).correlation)
    dmax = float(np.max(np.abs(P[both & cov] - Q[both & cov])))
    out[c] = dict(n_both=int(both.sum()), n_mine=int(np.isfinite(P).sum()), n_604=int(np.isfinite(Q).sum()),
                  sign_agree=sign_agree, spearman=rho, max_abs_dP_covered=dmax,
                  frac_flipped_vs_C4_mine=float(((P[cov] > 0.5) != (P4[cov] > 0.5)).mean()),
                  frac_flipped_vs_C4_604=float(((Q[cov & np.isfinite(Q)] > 0.5) != (P4[cov & np.isfinite(Q)] > 0.5)).mean()))
    print(f"(2-4) {c}: genes mine {out[c]['n_mine']} / 604 {out[c]['n_604']} / both {out[c]['n_both']}; sign agreement {sign_agree:.5f}; "
          f"Spearman {rho:.6f}; max |dP| covered {dmax:.2e}; flipped vs CLEAN4 mine {out[c]['frac_flipped_vs_C4_mine']:.3f} / 604 {out[c]['frac_flipped_vs_C4_604']:.3f}")
(HERE / "work/review604.json").write_text(json.dumps(out, indent=1))

"""Confound and build-detail checks on `donor_cos` before it is recommended.

(i)  Is it a disguised MAGNITUDE key?  096 killed `||dh_p||` and `kg_density` for exactly that.
(ii) Does it correlate with our OWN predicted magnitude, which would make it self-referential?
(iii) BUILD DETAIL: the donors cover different gene sets (k562_gwps 7,679 of 18,533; the two
     Orion screens 18,106).  A cosine taken on the union charges the sparse donor for its zeros,
     which is a constant per PAIR but not across pairs -- so on the challenge, where a target may
     have 3 pairs or 1, averaging raw pairwise cosines mixes scales and can inject a spurious
     dependence on donor COUNT.  Measured here both ways.
"""
from __future__ import annotations
import sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
from vcc2026.donor_transfer import Donor

D = json.load(open(f"{ROOT}/experiments/097-direction/work/perpert.json"))
print("### (i)/(ii) confound checks on the harnesses")
print(f"  {'harness':<10}{'rho(dcos, true|d_p|)':>22}{'rho(dcos, our|dh_p|)':>22}"
      f"{'rho(dcos, donor096)':>21}")
for nm in ("h1-test", "orion"):
    d = D[nm]
    dc = np.array(d["scores"]["donor_cos"]); sn = np.array(d["scores"]["self_norm"])
    d96 = np.array(d["scores"]["donor096"]); a = np.array(d["a"]); x = np.array(d["x"])
    cos_p = np.array(d["cos_p"])
    nd = np.sqrt(np.maximum(a, 1e-30)) * 0 + np.array(d["a"])  # placeholder
    # true |d_p| recovered from cos_p and x: |d_p| = x /(cos_p*|dh_p|)
    true_nd = np.where(np.abs(cos_p) > 1e-9, x / np.maximum(cos_p * np.sqrt(a), 1e-30), np.nan)
    m = np.isfinite(true_nd)
    print(f"  {nm:<10}{spearmanr(dc[m], true_nd[m]).statistic:>+22.3f}"
          f"{spearmanr(dc, sn).statistic:>+22.3f}{spearmanr(dc, d96).statistic:>+21.3f}")

print("\n### (iii) the challenge key, union-axis vs pair-shared-axis")
PROD = ("k562_gwps", "orion_hct116", "orion_hek293t")
targets = sorted(set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str)))
genes = np.asarray(pd.read_csv(f"{ROOT}/data/challenge/gene_names.csv").iloc[:, 0].astype(str))
Ms, hv, cov = {}, {}, {}
for nm in PROD:
    dn = Donor.load(f"{ROOT}/data/donors/{nm}")
    m, s = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cols = np.array([gi.get(g, -1) for g in genes]); ok = cols >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    M = np.zeros((len(targets), genes.size), np.float32); h = np.zeros(len(targets), bool)
    for r, t in enumerate(targets):
        if t in pi:
            M[r, ok] = s[pi[t]][cols[ok]]; h[r] = True
    Ms[nm] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
    hv[nm] = h; cov[nm] = ok
    del dn, s
ndon = sum(hv[n].astype(int) for n in PROD)
for mode in ("union", "pair-shared"):
    tot, cnt = np.zeros(len(targets)), np.zeros(len(targets))
    per_pair = {}
    for i in range(3):
        for j in range(i + 1, 3):
            a, b = Ms[PROD[i]], Ms[PROD[j]]
            if mode == "pair-shared":
                sh = cov[PROD[i]] & cov[PROD[j]]
                a, b = a[:, sh], b[:, sh]
            na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
            c = np.einsum("ij,ij->i", a, b) / np.maximum(na * nb, 1e-12)
            msk = hv[PROD[i]] & hv[PROD[j]]
            per_pair[f"{PROD[i][:6]}x{PROD[j][:6]}"] = (float(np.median(c[msk])), int(msk.sum()),
                                                        int(a.shape[1]))
            tot += np.where(msk, c, 0.0); cnt += msk
    k = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    ok = np.isfinite(k)
    print(f"\n  mode={mode}")
    for p, (med, n, g) in per_pair.items():
        print(f"    pair {p:<16} median cos {med:+.4f}  on {n:>3} targets, {g:>6} genes")
    print(f"    combined: median {np.nanmedian(k):+.4f}  "
          f"rho(key, n_donors) = {spearmanr(k[ok], ndon[ok]).statistic:+.3f}  "
          f"(a large value here means the key is partly counting DONORS, not agreement)")
    np.save(f"{ROOT}/experiments/097-direction/work/challenge_key_{mode}.npy", k)
kk = [np.load(f"{ROOT}/experiments/097-direction/work/challenge_key_{m}.npy")
      for m in ("union", "pair-shared")]
ok = np.isfinite(kk[0]) & np.isfinite(kk[1])
print(f"\n  rank correlation between the two build variants: "
      f"{spearmanr(kk[0][ok], kk[1][ok]).statistic:+.4f}")

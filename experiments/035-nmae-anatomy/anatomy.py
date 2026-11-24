"""035 — where does our `nmae` error actually come from?

`nmae` is our largest genuinely accessible gap (ours +0.065 scaled, Vivai +0.131, and `mse`
above us is closed at our cosine). It has never been decomposed, so it has been attacked only
with global knobs — alpha, the CPM gate, gamma. This asks which *genes* carry the error mass.

    nmae_p = sum_{g in DE} |pred_g - true_g| / sum_{g in DE} |true_g|

Every gene contributes a ratio to that sum, and there are only three ways to contribute:

  silent      pred_g = 0 (gated out, or no donor covers it) -> contributes exactly |true_g|,
              i.e. a ratio of 1. Cannot be worse than doing nothing, cannot be better.
  wrong sign  contributes |pred| + |true| > |true| -> a ratio ABOVE 1. This is the only way to
              be worse than silence, and the only category worth attacking directly.
  right sign  contributes ||pred| - |true|| -> below 1 whenever |pred| < 2|true|.

Knowing the split says which lever could matter. If the mass is in `silent`, coverage is the
problem. If it is in `wrong sign`, we are paying for confident errors and should shrink those
genes. If it is in `right sign` with |pred| far from |true|, it is a calibration problem.

Measured on the somatic HCT116 reference (132 real VCC-panel targets), against the shipped v8
effect matrix.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, TS          # noqa: E402
from vcc2026.multi_donor import MultiDonorModel        # noqa: E402
from vcc2026.whiten import Whitener                    # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
REAL = "data/harness/orion-hct116.h5ad"

a = ad.read_h5ad(REAL)
X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
lab = a.obs[PERT_COL].astype(str).to_numpy()
genes = np.asarray(a.var_names)
Pc = np.asarray(X[lab == CONTROL].sum(0)).ravel().astype(np.float64)
b_ctrl = np.log1p(TS * Pc / Pc.sum())
perts = sorted(set(lab) - {CONTROL})

donors = {n: Donor.load(Path(p)) for n, p in
          (("k562_gwps", "data/donors/k562_gwps"),
           ("orion_hek293t", "data/donors/orion_hek293t"))}   # HCT116 is the target
model = MultiDonorModel(donors, genes, gene_norm=False, eb=False)
shared = set(model.target_genes[model.tcols[model.names[0]]])
for n in model.names[1:]:
    shared &= set(model.target_genes[model.tcols[n]])
# The cached operator was fitted with the panel targets removed, as the runner does.
import pandas as pd  # noqa: E402
shared -= set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
ti = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
W = Whitener.load("experiments/031-somatic/out/whitener.npz")

w = {n: 1.0 for n in model.names}
E = np.vstack([model.effect(p, 1.0796, 0.0, weights=w) for p in perts]).astype(np.float64)
pre = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
E[:, wcols] = W.apply(E[:, wcols], 0.85, 0.0, renorm=False)
E[:, wcols] *= pre / (float(np.median(np.linalg.norm(E[:, wcols], axis=1))) or 1.0)
e = np.expm1(b_ctrl)
gate = (1e6 * e / e.sum()) > 5.0
E = np.where(gate & (b_ctrl > 0), E, 0.0)

K = 300
acc = {k: 0.0 for k in ("silent", "wrong", "right")}
cnt = {k: 0 for k in ("silent", "wrong", "right")}
den = 0.0
wrong_if_silent = 0.0
for i, p in enumerate(perts):
    Pp = np.asarray(X[lab == p].sum(0)).ravel().astype(np.float64)
    true = np.log1p(TS * Pp / Pp.sum()) - b_ctrl
    idx = np.argpartition(-np.abs(true), K - 1)[:K]
    t, pr = true[idx], E[i][idx]
    den += np.abs(t).sum()
    silent = pr == 0
    wrong = (~silent) & (np.sign(pr) != np.sign(t))
    right = (~silent) & (~wrong)
    for k, m in (("silent", silent), ("wrong", wrong), ("right", right)):
        acc[k] += np.abs(pr[m] - t[m]).sum()
        cnt[k] += int(m.sum())
    wrong_if_silent += np.abs(t[wrong]).sum()

tot = sum(acc.values())
print(f"reference: {len(perts)} perturbations, top-{K} genes by |true effect|\n")
print(f"{'category':12} {'genes':>8} {'share of genes':>15} {'share of error':>15} "
      f"{'nmae contribution':>18}")
for k in ("silent", "wrong", "right"):
    n = cnt[k]
    print(f"{k:12} {n:8,} {100*n/sum(cnt.values()):14.1f}% {100*acc[k]/tot:14.1f}% "
          f"{acc[k]/den:18.4f}")
print(f"\ntotal nmae = {tot/den:.4f}")
print(f"  a fully silent prediction would score exactly 1.0000")
print(f"\ncounterfactual: silencing every wrong-sign gene would move nmae to "
      f"{(acc['silent'] + acc['right'] + wrong_if_silent) / den:.4f}")
Path("experiments/035-nmae-anatomy/results.json").write_text(
    json.dumps({"acc": acc, "cnt": cnt, "den": den, "nmae": tot / den}, indent=2))

# ---------------------------------------------------------------- sign accuracy by agreement
#
# The decomposition says wrong-sign genes carry 53% of the error and that silencing them would
# take `nmae` to 0.909 -- past Vivai's 0.921. We cannot identify them directly, but if
# cross-donor sign agreement separates right from wrong at all, shrinking the disagreeing genes
# is a targeted lever rather than a global one.
print("\n" + "=" * 76)
print("Sign accuracy on the reference's DE genes, split by cross-donor agreement\n")
A = np.zeros((len(perts), len(genes)))
B = np.zeros((len(perts), len(genes)))
for i, p in enumerate(perts):
    for M, n in ((A, model.names[0]), (B, model.names[1])):
        if p in model.pert_index[n]:
            M[i, model.tcols[n]] = model.s[n][model.pert_index[n][p]][model.dcols[n]]

tot_ag = tot_dis = hit_ag = hit_dis = 0
err_ag = err_dis = 0.0
for i, p in enumerate(perts):
    Pp = np.asarray(X[lab == p].sum(0)).ravel().astype(np.float64)
    true = np.log1p(TS * Pp / Pp.sum()) - b_ctrl
    idx = np.argpartition(-np.abs(true), K - 1)[:K]
    t, pr = true[idx], E[i][idx]
    live = pr != 0
    agree = (np.sign(A[i][idx]) == np.sign(B[i][idx])) & (A[i][idx] != 0) & live
    dis = live & ~agree
    ok = np.sign(pr) == np.sign(t)
    tot_ag += int(agree.sum()); hit_ag += int((agree & ok).sum())
    tot_dis += int(dis.sum()); hit_dis += int((dis & ok).sum())
    err_ag += np.abs(pr[agree] - t[agree]).sum()
    err_dis += np.abs(pr[dis] - t[dis]).sum()

print(f"{'group':22} {'genes':>9} {'sign accuracy':>15} {'share of nmae':>15}")
print(f"{'donors agree':22} {tot_ag:9,} {100*hit_ag/max(tot_ag,1):14.1f}% "
      f"{err_ag/den:15.4f}")
print(f"{'donors disagree':22} {tot_dis:9,} {100*hit_dis/max(tot_dis,1):14.1f}% "
      f"{err_dis/den:15.4f}")
print(f"\nseparation: {100*hit_ag/max(tot_ag,1) - 100*hit_dis/max(tot_dis,1):+.1f} "
      f"percentage points")
print("A lever exists only if that separation is large enough to pay for shrinking a group "
      "that is still\nmostly right-signed -- shrinking a 50%-accurate group is free, "
      "shrinking a 55%-accurate one is not.")

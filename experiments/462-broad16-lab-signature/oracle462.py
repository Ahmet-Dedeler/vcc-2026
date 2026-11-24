"""462d — how much of 441's "recipient programme" oracle is the recipient's control-sample noise?

441's oracle predicts a held-out target's gate signs on h1-test from the leave-one-target-out (LOTO) mean of the other
84 targets' reference log2 fold changes. Every one of those fold changes is taken against the SAME 18,400 control cells
as the held-out target's truth. So the oracle also "predicts" the control sample's per-gene deviation, which is a common
offset in every target's truth and which no external source can know. This splits the controls into two random halves:
  shared   oracle from lfc(others vs controls A)  against truth lfc(target vs controls A)   (441's construction, half-size)
  disjoint oracle from lfc(others vs controls B)  against truth lfc(target vs controls A)   (no shared noise)
  441      oracle from lfc(others vs all)         against truth lfc(target vs all)          (as 441 measured)
Gate: the scorer's reference table (q < 0.05, finite lfc, own gene out, targets with >= 10 gate genes; 441's 58).
Fold change = the scorer's: log2 of arithmetic means of per-cell CPM (1e6), epsilon 1e-9, on the table's kept genes.
KA: lfc(target vs all controls) must equal the reference table's lfc.
Also CLEAN4 (391's pinned per-gene prior), which does not depend on the recipient's controls, as a fixed comparator.
usage: uv run python experiments/462-broad16-lab-signature/oracle462.py   -> work/oracle462.json
"""
import csv
import glob
import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
EPS = 1e-9
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str))
q, lfc_ref = z["q"], z["lfc"]
kg = genes[keep]
with h5py.File(ROOT / "data/harness/h1-test.h5ad", "r") as f:
    X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=tuple(f["X"].attrs["shape"]))
    tg = f["obs/target_gene"]; tv = tg["categories"][:].astype(str)[tg["codes"][:]]
    vn = f["var/_index"]; hg = (vn[:] if not isinstance(vn, h5py.Group) else vn["values"][:]).astype(str)
assert np.array_equal(hg, genes), "h1-test gene axis != reference genes"
lib = np.asarray(X.sum(1)).ravel()
C = sp.diags(1e6 / lib) @ X                              # per-cell CPM
C = C[:, np.flatnonzero(keep)].tocsr()
ctrl = np.flatnonzero(tv == "non-targeting")
rng = np.random.default_rng(462)
perm = rng.permutation(ctrl); A, B = perm[: ctrl.size // 2], perm[ctrl.size // 2:]
mean = lambda idx: np.asarray(C[idx].mean(0)).ravel()
mF, mA, mB = mean(ctrl), mean(A), mean(B)
M = np.stack([mean(np.flatnonzero(tv == t)) for t in targets])
lg = lambda m, c: np.log2((m + EPS) / (c + EPS))
LF, LA, LB = lg(M, mF), lg(M, mA), lg(M, mB)
fin = np.isfinite(lfc_ref) & np.isfinite(LF)
ka = float(np.nanmax(np.abs(LF[fin] - lfc_ref[fin])))
dif = np.abs(LF[fin] - lfc_ref[fin]); mc = np.broadcast_to(mF, LF.shape)[fin]
print(f"KA: max |lfc(all controls) - reference lfc| = {ka:.3e}; median {np.median(dif):.1e}; 99.9th pct {np.quantile(dif, .999):.1e}; "
      f"at the max, control mean CPM {mc[np.argmax(dif)]:.3g}", flush=True)
assert ka < 1e-4 and np.median(dif) < 1e-5, "KA failed"   # the scorer accumulates float32 CPM: ~1e-6 in log2
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc_ref) & (lfc_ref != 0) & ~own
use = np.flatnonzero(gate.sum(1) >= 10)
print(f"gated targets {use.size} (441: 58); gate pairs {int(gate[use].sum())}", flush=True)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
pos = {g: i for i, g in enumerate(offi)}
c4s = np.asarray(c4["psign"], float)
c4k = np.array([c4s[pos[g]] if g in pos else np.nan for g in kg])
c4cov = np.array([np.isfinite(np.asarray(c4["P_down"], float)[pos[g]]) if g in pos else False for g in kg])


def acc(pred, truth, weight_by):
    """gate sign accuracy over the gated targets: pooled unweighted and |l|-weighted (l = the reference lfc)."""
    hit, wsum, n, whit = 0.0, 0.0, 0, 0.0
    for t in use:
        g = gate[t] & np.isfinite(pred[t]) & (pred[t] != 0)
        h = np.sign(pred[t][g]) == np.sign(truth[t][g])
        w = np.abs(weight_by[t][g])
        hit += h.sum(); n += g.sum(); whit += (h * w).sum(); wsum += w.sum()
    return round(hit / n, 4), round(whit / wsum, 4)


def loto(L):
    tot = L.sum(0)
    return (tot[None, :] - L) / (L.shape[0] - 1)


res = {"ka_max_abs_lfc": ka, "gated_targets": int(use.size), "gate_pairs": int(gate[use].sum())}
res["oracle_441_all_controls"] = acc(loto(LF), LF, lfc_ref)
res["oracle_shared_half_A"] = acc(loto(LA), LA, lfc_ref)
res["oracle_disjoint_B_to_A"] = acc(loto(LB), LA, lfc_ref)
res["oracle_disjoint_A_to_B"] = acc(loto(LA), LB, lfc_ref)
res["oracle_shared_half_B"] = acc(loto(LB), LB, lfc_ref)
c4pred = np.tile(np.where(c4cov, c4k, np.nan), (len(targets), 1))
res["clean4_vs_truth_all"] = acc(c4pred, LF, lfc_ref)
res["clean4_vs_truth_A"] = acc(c4pred, LA, lfc_ref)
res["clean4_coverage_of_gate"] = float(np.mean([c4cov[gate[t]].mean() for t in use]))
d = lambda k: (res[k][0] + res[f"{k.replace('_B_to_A', '_A_to_B')}"][0]) / 2
res["summary"] = dict(
    oracle_shared_half=round((res["oracle_shared_half_A"][0] + res["oracle_shared_half_B"][0]) / 2, 4),
    oracle_disjoint=round((res["oracle_disjoint_B_to_A"][0] + res["oracle_disjoint_A_to_B"][0]) / 2, 4),
    note="shared - disjoint = the part of the LOTO oracle that is the recipient control sample's noise (half controls)")
print(json.dumps(res, indent=1))
(HERE / "work/oracle462.json").write_text(json.dumps(res, indent=1))

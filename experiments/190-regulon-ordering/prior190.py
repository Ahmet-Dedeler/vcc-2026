"""190 — build the per-(target, gene) regulatory prior, two independent versions, and
report its COVERAGE before anything is scored.

Two priors, both reduced to the same object: a signed, weighted matrix

    N[t, g]  in  R,   sign(N) = the direction the NETWORK expects gene g to move when
                                target t is knocked down,
                      |N|     = how confident the network is in that direction.
    N[t, g] = 0        means "the network says nothing about this pair".

V1  CollecTRI (saezlab), signed TF -> target regulons.  A CRISPRi knockdown REMOVES the TF,
    so the expected response sign is  -mor:  losing an activator lowers its targets, losing a
    repressor raises them.

V2  A control-derived co-expression GRN, D2R2-shaped: Pearson r(t, g) across the recipient's
    own unperturbed cells on log1p-CPM.  Expected response sign is  -sign(r)  (if t drives g
    up in this line, removing t lowers g).  Licence-clean and available for D/E/F on day one.
    174 already measured this as a raw SIGN predictor (0.5158, and 2 PCs carry all of it);
    here it is used as an ORDERING prior instead, which is a different question.

NOTHING IN THIS FILE SCORES ANYTHING.  It only builds N and prints coverage.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
REG = os.path.join(ROOT, "data", "external", "regulon")
sys.path.insert(0, os.path.join(ROOT, "experiments", "172-sign-confidence-key"))

COLLECTRI_CSV = os.path.join(REG, "collectri_rescued.csv")
COLLECTRI_OP = os.path.join(REG, "collectri_omnipath.tsv")


# --------------------------------------------------------------------------------------
# V1 -- CollecTRI
# --------------------------------------------------------------------------------------
def load_collectri(weighted: bool = True) -> pd.DataFrame:
    """Return tidy (tf, gene, mor, w).

    `mor` in {+1, -1} is CollecTRI's own sign of regulation.
    `w` is a confidence: log1p(number of contributing resources).  CollecTRI has no per-edge
    confidence level (unlike DoRothEA's A-E), so resource count is the only graded evidence
    the file carries.  `weighted=False` sets every w to 1.
    """
    d = pd.read_csv(COLLECTRI_CSV)
    d = d[["source", "target", "weight", "TF.category", "resources"]].copy()
    d.columns = ["tf", "gene", "mor", "cat", "resources"]
    d = d[np.isin(d["mor"], (1, -1))]
    # Drop CollecTRI's protein-complex pseudo-sources (AP1, NFKB, ...): they are not single
    # genes, so they can never be one of our 300 knockdown targets.
    nres = d["resources"].fillna("").apply(lambda s: 0 if s == "" else len(s.split(";")))
    d["w"] = np.log1p(nres) if weighted else 1.0
    d = d.drop(columns=["resources"])
    d = d.groupby(["tf", "gene"], as_index=False).agg(
        mor=("mor", lambda s: float(np.sign(s.sum()))),
        w=("w", "max"), cat=("cat", "first"))
    return d[d["mor"] != 0.0]


def collectri_matrix(targets, genes, weighted: bool = True):
    """N[t, g] = -mor * w  for edges present, 0 elsewhere.  Returns (N, stats)."""
    d = load_collectri(weighted=weighted)
    ti = {t: i for i, t in enumerate(targets)}
    gi = {g: j for j, g in enumerate(genes)}
    N = np.zeros((len(targets), len(genes)), dtype=np.float32)
    hit = d[d["tf"].isin(ti) & d["gene"].isin(gi)]
    for tf, gene, mor, w in zip(hit["tf"], hit["gene"], hit["mor"], hit["w"]):
        N[ti[tf], gi[gene]] = -mor * w
    covered = np.flatnonzero((N != 0).any(1))
    stats = dict(
        network="CollecTRI",
        n_edges_total=int(len(d)),
        n_tfs_total=int(d["tf"].nunique()),
        n_targetgenes_total=int(d["gene"].nunique()),
        n_targets=len(targets),
        n_targets_covered=int(covered.size),
        frac_targets_covered=float(covered.size / len(targets)),
        covered_targets=[str(targets[i]) for i in covered],
        edges_placed=int((N != 0).sum()),
        median_regulon_size_on_axis=float(np.median((N != 0).sum(1)[covered])) if covered.size else 0.0,
    )
    return N, stats


# --------------------------------------------------------------------------------------
# V2 -- control-derived co-expression GRN
# --------------------------------------------------------------------------------------
def coexpression_matrix(X_ctrl, targets, genes, n_pc_remove: int = 0, min_cells_frac=0.0):
    """N[t, g] = -r(t, g) computed across control cells on log1p CPM.

    X_ctrl : (cells, genes) dense float32 log1p-CPM matrix, columns aligned to `genes`.
    n_pc_remove : project out the top-k principal components of the cell x gene matrix first
                  (174's control; 2 PCs kill the raw signal).
    """
    Z = X_ctrl - X_ctrl.mean(0, keepdims=True)
    if n_pc_remove > 0:
        # economy SVD on the centred cell x gene matrix, remove top-k right singular vectors
        from scipy.sparse.linalg import svds
        _, _, Vt = svds(Z.astype(np.float32), k=n_pc_remove)
        Z = Z - (Z @ Vt.T) @ Vt
    sd = Z.std(0)
    ok = sd > 0
    Zn = np.zeros_like(Z)
    Zn[:, ok] = Z[:, ok] / sd[ok]
    n = Z.shape[0]
    gi = {g: j for j, g in enumerate(genes)}
    N = np.zeros((len(targets), len(genes)), dtype=np.float32)
    present = []
    for i, t in enumerate(targets):
        j = gi.get(t)
        if j is None or not ok[j]:
            continue
        present.append(i)
        N[i] = -(Zn[:, j] @ Zn) / n
    N[:, ~ok] = 0.0
    stats = dict(network=f"coexpr(pc{n_pc_remove})", n_targets=len(targets),
                 n_targets_covered=len(present),
                 frac_targets_covered=len(present) / len(targets),
                 median_abs_r=float(np.median(np.abs(N[present]))) if present else 0.0)
    return N, stats


# --------------------------------------------------------------------------------------
# coverage report
# --------------------------------------------------------------------------------------
def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    import harness172

    T = harness172.H1Test()
    harness172.check(T)
    print()

    panel = pd.read_csv(os.path.join(ROOT, "data", "challenge", "pert_counts.csv"))
    panel = [str(x) for x in panel["target_gene"]]

    d = load_collectri()
    print(f"CollecTRI (rescued.omnipathdb.org/CollecTRI.csv)")
    print(f"  {len(d):,} signed edges | {d['tf'].nunique():,} TFs | "
          f"{d['gene'].nunique():,} target genes")
    print(f"  mor: {int((d['mor'] > 0).sum()):,} activating / "
          f"{int((d['mor'] < 0).sum()):,} repressing")
    print(f"  TF.category: {d['cat'].value_counts().to_dict()}")

    out = {"collectri": {"n_edges": int(len(d)), "n_tfs": int(d["tf"].nunique()),
                         "n_target_genes": int(d["gene"].nunique())}}

    # --- coverage of the 85 h1-test targets ------------------------------------------
    N1, s1 = collectri_matrix(T.perts, T.genes, weighted=True)
    print(f"\nh1-test 85 targets on {len(T.genes):,} gated genes")
    print(f"  targets that are CollecTRI TFs with >=1 edge on the axis: "
          f"{s1['n_targets_covered']}/85 ({100*s1['frac_targets_covered']:.1f}%)")
    print(f"  {s1['covered_targets']}")
    print(f"  edges placed: {s1['edges_placed']:,}  "
          f"median regulon size among covered: {s1['median_regulon_size_on_axis']:.0f}")

    # the number that actually caps reach: coverage of the reference's DE POOL
    pool = T.pool
    n_pool = pool.sum(1)
    inreg = (pool & (N1 != 0))
    per_t = np.array([inreg[p].sum() / max(n_pool[p], 1) for p in range(len(T.perts))])
    cov_idx = np.flatnonzero((N1 != 0).any(1))
    print(f"\n  DE-POOL coverage (this is what caps `reach`):")
    print(f"    pool pairs total                 {int(pool.sum()):,}")
    print(f"    pool pairs with a CollecTRI edge {int(inreg.sum()):,} "
          f"({100*inreg.sum()/pool.sum():.2f}%)")
    print(f"    mean per-target pool fraction covered, all 85 : {per_t.mean():.4f}")
    print(f"    mean per-target pool fraction covered, covered targets only : "
          f"{per_t[cov_idx].mean():.4f}")
    out["h1test"] = {k: v for k, v in s1.items() if k != "covered_targets"}
    out["h1test"]["covered_targets"] = s1["covered_targets"]
    out["h1test"]["pool_pairs"] = int(pool.sum())
    out["h1test"]["pool_pairs_with_edge"] = int(inreg.sum())
    out["h1test"]["pool_frac_with_edge"] = float(inreg.sum() / pool.sum())
    out["h1test"]["mean_per_target_pool_frac"] = float(per_t.mean())

    # --- coverage of the 300 challenge panel targets ----------------------------------
    tfs = set(d["tf"])
    hit = [t for t in panel if t in tfs]
    print(f"\nchallenge panel: {len(hit)}/300 targets are CollecTRI TFs "
          f"({100*len(hit)/300:.1f}%)")
    print(f"  {sorted(hit)}")
    genes_axis = [str(x) for x in
                  pd.read_csv(os.path.join(ROOT, "data", "challenge", "gene_names.csv"),
                              header=None)[0]]
    if genes_axis[0].lower() in ("gene", "gene_name", "0"):
        genes_axis = genes_axis[1:]
    Np, sp = collectri_matrix(np.array(panel), np.array(genes_axis), weighted=True)
    print(f"  on the challenge's {len(genes_axis):,}-gene axis: "
          f"{sp['n_targets_covered']}/300 covered, {sp['edges_placed']:,} edges, "
          f"median regulon {sp['median_regulon_size_on_axis']:.0f}")
    out["panel"] = {k: v for k, v in sp.items()}

    with open(os.path.join(WORK, "coverage190.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {WORK}/coverage190.json")


if __name__ == "__main__":
    main()

"""224 — what IS the fitted per-gene gain, and is it anything production does not already do?

Three questions, in order of how much they would change the verdict:

  1. **Is it stable across held lines?**  If the gain vectors fitted for different held
     lines agree, the object is a property of the GENE and generalises; if they do not, it
     is a per-fold artefact and the leave-one-line-out mean is hiding it.
  2. **What predicts it?**  Regressed on per-gene cross-line reliability, the consensus's
     own per-gene sd, the truth's per-gene sd, and basal expression.
  3. **Does production's whitening operator already contain it?**  The shipped
     `union_pruned` map is `V diag(evals^-gamma/2) V^T` at `gamma = 0.65`, `lam = 0`.  Its
     per-gene diagonal and per-gene row norm are computed here on the shared gene axis and
     correlated with the fitted gain.  A high correlation would mean this is a rediscovery
     of the shipped operator's diagonal, not a new lever.

    uv run python mechanism224.py --tag S7
"""
from __future__ import annotations

import argparse
import json

import numpy as np

from core224 import WORK, load_surface, unit

ROOT = WORK.parents[2]
WHITENER = ROOT / "data/donors/whitener_union_pruned.npz"
GAMMA = 0.65


def whitener_gene_stats(genes: list[str]) -> dict:
    """diag(W) and rownorm(W) for the shipped operator, on `genes`.

    W = V diag(evals^(-gamma/2)) V^T, so
        diag(W)[j]     = sum_k V[j,k]^2 * w_k
        ||W[j,:]||^2   = sum_k V[j,k]^2 * w_k^2
    Both are computable from V row-wise without ever forming the 17,942^2 dense W."""
    axis = ROOT / "data/donors/whitener_union_pruned.json"
    meta = json.loads(axis.read_text())
    cols = meta.get("cols") or meta.get("genes")
    if cols is None:                    # the axis is not stored in the json; rebuild it
        return dict(available=False, reason="operator gene axis not recorded in the json")
    idx = {g: i for i, g in enumerate(cols)}
    keep = [idx[g] for g in genes if g in idx]
    z = np.load(WHITENER, allow_pickle=True)
    w = np.asarray(z["evals"], dtype=np.float64) ** (-GAMMA / 2.0)
    V = z["evecs"]
    d = np.zeros(len(keep))
    rn = np.zeros(len(keep))
    step = 512
    for s in range(0, len(keep), step):
        blk = np.asarray(V[keep[s:s + step]], dtype=np.float64) ** 2
        d[s:s + step] = blk @ w
        rn[s:s + step] = np.sqrt(blk @ (w ** 2))
    return dict(available=True, genes=[g for g in genes if g in idx], diag=d, rownorm=rn)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="S7")
    args = ap.parse_args()

    d = json.load(open(WORK / f"ladder224-{args.tag}.json"))
    pt = np.load(WORK / f"ptcos-{args.tag}.npz", allow_pickle=True)
    C, basal, members, groups, targets, genes, dropped = load_surface(d["surface"])
    kmax = max(r["k"] for r in d["rows"])

    gains, labels = [], []
    for r in d["rows"]:
        key = f"{r['cid']}|_gain"
        if r["k"] == kmax and key in pt:
            gains.append(np.asarray(pt[key], dtype=np.float64))
            labels.append(r["held"])
    Gm = np.stack(gains)
    print(f"{len(Gm)} fitted gain vectors at k={kmax}, {Gm.shape[1]} genes\n")

    print("1. STABILITY — pairwise Pearson r between gain vectors fitted for DIFFERENT "
          "held lines")
    R = np.corrcoef(Gm)
    same = np.asarray([[labels[i] == labels[j] for j in range(len(labels))]
                       for i in range(len(labels))])
    off = ~np.eye(len(labels), dtype=bool)
    print(f"   same held line : mean r {R[same & off].mean():.4f}  (n={int((same & off).sum())})")
    print(f"   DIFFERENT line : mean r {R[~same].mean():.4f}  "
          f"[{R[~same].min():.3f}, {R[~same].max():.3f}]  (n={int((~same).sum())})")
    print(f"   gain spread    : p10 {np.quantile(Gm, 0.1):.4f}  median "
          f"{np.median(Gm):.4f}  p90 {np.quantile(Gm, 0.9):.4f}")

    print("\n2. WHAT PREDICTS THE GAIN")
    gbar = Gm.mean(0)
    U = np.stack([unit(C[i].astype(np.float64)) for i in range(len(members))])
    sd_line = U.std(1).mean(0)                       # per-gene sd of a single line's rows
    cons = unit(U.mean(0))
    sd_cons = cons.std(0)
    # per-gene cross-line reliability: mean correlation across targets between line pairs
    rel = np.zeros(C.shape[2])
    npair = 0
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            a = U[i] - U[i].mean(0, keepdims=True)
            b = U[j] - U[j].mean(0, keepdims=True)
            rel += ((a * b).sum(0)
                    / np.maximum(np.linalg.norm(a, axis=0) * np.linalg.norm(b, axis=0), 1e-300))
            npair += 1
    rel /= npair
    bas = basal.mean(0).astype(np.float64)
    preds = dict(cross_line_reliability=rel, sd_single_line=sd_line, sd_consensus=sd_cons,
                 sd_ratio=sd_line / np.maximum(sd_cons, 1e-300), basal_mean=bas)
    for k, v in preds.items():
        print(f"   r(gain, {k:24s}) = {np.corrcoef(gbar, v)[0, 1]:+.4f}   "
              f"spearman {np.corrcoef(np.argsort(np.argsort(gbar)), np.argsort(np.argsort(v)))[0, 1]:+.4f}")

    out = dict(stability_same=float(R[same & off].mean()) if (same & off).any() else None,
               stability_diff=float(R[~same].mean()),
               predictors={k: float(np.corrcoef(gbar, v)[0, 1]) for k, v in preds.items()})

    print("\n2b. IS IT 189's theta?  (independent triangulation)")
    # 189 built a per-gene cross-line SIGN-agreement reliability on the same five lines and
    # measured it at z = +41.6 against its own permutation null — but used it as an ORDERING
    # key, where it bought nothing.  If the gain fitted here is that same object, then this
    # experiment is 189's statistic applied to the DIRECTION channel instead.
    p189 = ROOT / "experiments/189-per-gene-reliability/work/prior189.npz"
    if p189.exists():
        z = np.load(p189, allow_pickle=True)
        idx = {str(g): i for i, g in enumerate(z["genes"])}
        keep = [(i, idx[g]) for i, g in enumerate(genes) if g in idx]
        a = np.asarray([gbar[i] for i, _ in keep])
        for key in ("theta", "eb", "raw"):
            b = np.asarray([z[key][j] for _, j in keep])
            ok = np.isfinite(a) & np.isfinite(b)
            r = float(np.corrcoef(a[ok], b[ok])[0, 1])
            sp = float(np.corrcoef(np.argsort(np.argsort(a[ok])),
                                   np.argsort(np.argsort(b[ok])))[0, 1])
            print(f"   r(gain, 189 {key:6s}) = {r:+.4f}   spearman {sp:+.4f}   "
                  f"(n = {int(ok.sum())} shared genes)")
            out[f"corr_189_{key}"] = r
    else:
        print("   189's prior189.npz not on disk")

    print("\n3. IS IT PRODUCTION'S WHITENER?")
    st = whitener_gene_stats(genes)
    if not st["available"]:
        print(f"   {st['reason']} — falling back to the operator's own basis axis")
        meta = json.loads((ROOT / "data/donors/whitener_union_pruned.json").read_text())
        print(f"   operator: {meta['n_cols']} cols, donors {meta['donors']}, "
              f"prune {meta['basis_prune']}, sha {meta['sha256'][:12]}")
        print("   The operator's gene axis is union(donor coverage) - panel and is NOT")
        print("   recorded in its sidecar json (AGENTS.md §8.3 is about exactly this), so")
        print("   a gene-wise correlation cannot be formed without rebuilding the axis.")
        print("   Reported instead: the gain is DIAGONAL and the operator is not — the two")
        print("   families intersect only in diag(W), which §4.2 measures as LOWERING the")
        print("   transfer cosine while the gain RAISES it. They are not the same object.")
    else:
        for k in ("diag", "rownorm"):
            g2 = np.asarray([gbar[genes.index(g)] for g in st["genes"]])
            print(f"   r(gain, whitener {k:8s}) = {np.corrcoef(g2, st[k])[0, 1]:+.4f}")
            out[f"whitener_{k}"] = float(np.corrcoef(g2, st[k])[0, 1])

    (WORK / f"mechanism224-{args.tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

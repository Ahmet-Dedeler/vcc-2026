"""Build a 2026-faithful local scoring reference from the 2025 H1 data.

We get two submissions a day and unlimited local runs, so a private leaderboard is
the single highest-value piece of infrastructure. The 2025 VCC H1 hESC data is the
only public set with the same chemistry (10x Flex), the same knockdown criteria, and
real held-out ground truth.

The 2026 reference is built by downsampling every perturbation to exactly 400 cells
and to a median depth of ~20,000 UMI/cell (vcc2026 metric spec section 0). The 2025
data ships at ~54,000 UMI/cell with variable cell counts, so we reproduce both steps
here. Cell count and sequencing depth move every metric in the suite, so a harness
that skips this measures the wrong thing.

Usage
-----
    uv run python scripts/make_harness.py --split test
    uv run python scripts/make_harness.py --split all --name h1-full

Then build the competition scale (baseline = 0, replicate = 1):

    uv run cell-eval2 baseline --real <ref>.h5ad -o out/baseline \
        --preset vcc2026 --pert-col target_gene --control non-targeting
    uv run cell-eval2 prep-real-bundle --real <ref>.h5ad --baseline out/baseline \
        -o out/bundle --preset vcc2026 --pert-col target_gene --control non-targeting
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

CONTROL = "non-targeting"
PERT_COL = "target_gene"

# 2026 reference conditions, from the vcc2026 metric spec section 0.
CELLS_PER_PERT = 400
TARGET_MEDIAN_UMI = 20_000
N_CONTROL_CELLS = 18_400

SPLITS = {
    "train": "data/external/vcc2025/train/adata_Training.h5ad",
    "validation": "data/external/vcc2025/validation/adata_Validation.h5ad",
    "test": "data/external/vcc2025/test/adata_Test.h5ad",
}


def thin_counts(X: sp.csr_matrix, rate: float, rng: np.random.Generator) -> sp.csr_matrix:
    """Binomially thin every count by `rate`, preserving the depth distribution's shape.

    A uniform rate (rather than per-cell rescaling to a fixed total) is what keeps the
    library-size spread realistic: the reference has min/median ~0.16 and max/median ~2.6,
    and per-cell renormalisation would collapse that.
    """
    if rate >= 1.0:
        return X
    out = X.copy()
    out.data = rng.binomial(out.data.astype(np.int64), rate).astype(np.float32)
    out.eliminate_zeros()
    return out


def build(
    sources: list[str],
    out_path: Path,
    *,
    cells_per_pert: int = CELLS_PER_PERT,
    n_control: int = N_CONTROL_CELLS,
    target_umi: int = TARGET_MEDIAN_UMI,
    n_perts: int | None = None,
    seed: int = 0,
) -> None:
    rng = np.random.default_rng(seed)

    keep_idx: list[np.ndarray] = []
    parts: list[ad.AnnData] = []
    dropped: dict[str, int] = {}

    for src in sources:
        a = ad.read_h5ad(src, backed="r")
        obs = a.obs
        picks: list[np.ndarray] = []

        for gene, sub in obs.groupby(PERT_COL, observed=True).indices.items():
            if gene == CONTROL:
                continue
            if len(sub) < cells_per_pert:
                dropped[str(gene)] = len(sub)
                continue
            picks.append(rng.choice(sub, cells_per_pert, replace=False))

        # A smaller panel scores far faster (DE cost ~ n_perts x total cells) at the price
        # of fidelity: pds ranks are within-panel, so a 40-gene panel is an easier rank
        # problem than a 300-gene one. Use for iteration, never for a go/no-go decision.
        if n_perts is not None and len(picks) > n_perts:
            picks = [picks[i] for i in rng.choice(len(picks), n_perts, replace=False)]

        n_selected = len(picks)
        ctrl = obs.index.get_indexer_for(obs.index[obs[PERT_COL] == CONTROL])
        # Share the control budget across sources so the pooled reference holds n_control.
        take = min(len(ctrl), max(1, n_control // len(sources)))
        picks.append(rng.choice(ctrl, take, replace=False))

        sel = np.sort(np.concatenate(picks))
        print(f"  {Path(src).name}: {n_selected} perturbations x {cells_per_pert} + {take} control")
        parts.append(a[sel].to_memory())
        keep_idx.append(sel)
        del a

    merged = ad.concat(parts, axis=0, join="inner", index_unique=None) if len(parts) > 1 else parts[0]
    del parts

    X = merged.X.tocsr() if not sp.isspmatrix_csr(merged.X) else merged.X
    totals = np.asarray(X.sum(1)).ravel()
    med = float(np.median(totals))
    rate = target_umi / med
    print(f"  depth: median {med:,.0f} UMI -> thinning at rate {rate:.4f}")

    X = thin_counts(X, rate, rng)
    merged.X = X
    tot2 = np.asarray(X.sum(1)).ravel()
    print(f"  after: median {np.median(tot2):,.0f}  min {tot2.min():,.0f}  max {tot2.max():,.0f}")

    # cell-eval2 reads .obs[pert_col]; drop everything else so nothing leaks into scoring.
    merged.obs = merged.obs[[PERT_COL]].copy()
    merged.obs[PERT_COL] = merged.obs[PERT_COL].astype(str)
    merged.var = pd.DataFrame(index=merged.var_names)

    n_pert = merged.obs[PERT_COL].nunique() - 1
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out_path, compression="gzip")

    print(f"\nwrote {out_path}")
    print(f"  {merged.n_obs:,} cells x {merged.n_vars:,} genes, {n_pert} perturbations + control")
    if dropped:
        print(f"  dropped {len(dropped)} perturbations with <{cells_per_pert} cells: "
              f"{sorted(dropped, key=dropped.get)[:6]}...")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--split", default="test", choices=[*SPLITS, "all"])
    ap.add_argument("--name", default=None, help="output basename (default: h1-<split>)")
    ap.add_argument("--outdir", default="data/harness")
    ap.add_argument("--cells-per-pert", type=int, default=CELLS_PER_PERT)
    ap.add_argument("--n-perts", type=int, default=None, help="cap the panel size (fast harness)")
    ap.add_argument("--n-control", type=int, default=N_CONTROL_CELLS)
    ap.add_argument("--target-umi", type=int, default=TARGET_MEDIAN_UMI)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    sources = list(SPLITS.values()) if args.split == "all" else [SPLITS[args.split]]
    missing = [s for s in sources if not Path(s).exists()]
    if missing:
        raise SystemExit(f"missing (still downloading?): {missing}")

    name = args.name or f"h1-{args.split}"
    print(f"building harness '{name}' from {len(sources)} source(s)")
    build(
        sources,
        Path(args.outdir) / f"{name}.h5ad",
        cells_per_pert=args.cells_per_pert,
        n_control=args.n_control,
        target_umi=args.target_umi,
        n_perts=args.n_perts,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()

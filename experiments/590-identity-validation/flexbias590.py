"""590 (b) — a per-gene "Arc Flex minus 3'" log2 bias, so 3' profiles can be scored as if Arc had measured them in Flex.

Each profile's deviation from its own line's CCLE profile cancels the line identity:
    dev_q[g] = q[g] − ccle_true(q)[g]        (log2 CPM+1 vs log2 TPM+1 on the official ∩ CCLE axis)
The Flex side is A/B/C (Arc, Flex, 18,400 cells; true lines JURKAT / HELA / CAL33 per 041 + 356); the 3' side is the
MIX-seq control profiles whose consensus call was correct. bias = mean_Flex dev − mean_3' dev, per gene, set to 0 where
the gene is not expressed in CCLE on average (mean log2 TPM ≤ 1) or missing on either side.

Two-fold cross-fitting over MIX-seq lines (by DepMap ID parity of the sorted list): fold k's bias uses only the other
fold's 3' lines, and `match590.py --flexbias work/flexbias590_fold<k>.npz` is read only on fold k's MIX-seq rows (plus
every non-MIX-seq 3' row, which neither fold fits on). Writes work/flexbias590_fold{0,1}.npz.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE))
import def_identity as DI  # noqa: E402
from match590 import counts_query  # noqa: E402

import anndata as ad  # noqa: E402


def main() -> None:
    Xfull, ccle_genes, samples, clinical = DI.load_cache()
    official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    g2c = {g: i for i, g in enumerate(ccle_genes)}
    genes = np.asarray([g for g in official if g in g2c])
    X = Xfull[:, np.array([g2c[g] for g in genes])].astype(np.float64)
    del Xfull
    srow = {s: i for i, s in enumerate(samples)}
    expressed = X.mean(axis=0) > 1.0

    flex_dev = []
    for ctx, true in DI.CALL_041.items():
        h = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad", backed="r")
        pos = {g: i for i, g in enumerate(h.var_names.astype(str))}
        q = DI.log2_cpm_from_counts(h.X)[[pos[g] for g in genes]]
        h.file.close()
        flex_dev.append(q - X[srow[true]])
    flex_dev = np.vstack(flex_dev)

    match = pd.read_csv(HERE / "work/match590.csv")
    ok = match[(match.dataset == "mcfarland_ctrl") & match.correct]
    z = np.load(HERE / "work/profiles/mcfarland_ctrl.npz")
    meta = json.loads(str(z["meta"]))
    rows = {r["query"]: (r, c) for r, c in zip(meta, z["counts"])}
    # folds over EVERY MIX-seq line (ctrl or all-conditions profile, right or wrong call), so no row is dropped;
    # the bias of fold k is fitted only on the other fold's correctly-called control profiles
    all_ids = sorted({q.split(":", 1)[1] for q in match["query"] if q.startswith("mcfarland")})
    folds = {0: all_ids[0::2], 1: all_ids[1::2]}
    ok_ids = {q.split(":", 1)[1] for q in ok["query"]}
    for k in (0, 1):
        fit_on = [f"mcfarland_ctrl:{i}" for i in folds[1 - k] if i in ok_ids]
        dev = []
        for qn in fit_on:
            r, c = rows[qn]
            dev.append(counts_query(c, official, genes) - X[srow[r["true"]]])
        dev = np.vstack(dev)
        b = np.nanmean(flex_dev, axis=0) - np.nanmean(dev, axis=0)
        b = np.where(expressed & np.isfinite(b), b, 0.0)
        np.savez(HERE / f"work/flexbias590_fold{k}.npz", genes=genes.astype("U"), bias=b,
                 fitted_on=np.asarray(fit_on + ["A", "B", "C"], dtype="U"),
                 apply_to=np.asarray([f"mcfarland:{i}" for i in folds[k]], dtype="U"))
        print(f"fold {k}: fitted on {len(fit_on)} MIX-seq lines + A/B/C; applies to {len(folds[k])}; "
              f"|bias| median {np.median(np.abs(b[expressed])):.3f}, sd {np.std(b[expressed]):.3f} log2, "
              f"genes with bias {int((b != 0).sum())}", flush=True)
    # how much of the Flex deviation is shared across A/B/C (a chemistry/lab signal) vs line-specific
    cc = np.corrcoef(np.where(expressed, flex_dev, 0))
    print("A/B/C Flex-deviation correlations (expressed genes):", np.round(cc[np.triu_indices(3, 1)], 3).tolist())


if __name__ == "__main__":
    main()

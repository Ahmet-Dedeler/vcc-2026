"""Explain the known-answer failure on MLLT1, without scoring any new arm."""
import time

import anndata as ad
import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, norm

import evaluate323 as E
import run323 as R
from cell_eval2.de_compute import compute_de, _to_linear
from vcc2026 import rc_emit


def main():
    start = time.monotonic()
    all_data = ad.read_h5ad(R.REAL)
    a = all_data[all_data.obs.target_gene.isin(["MLLT1", E.CTRL])].copy()
    del all_data
    official = compute_de(
        a, backend="scanpy", groupby="target_gene", reference=E.CTRL,
        mean_calc="arithmetic", epsilon=1e-9, input_type="counts",
        target_sum=1e6, clip_value=None, fdr_scope="per_pert",
        filter_gene_min_cpm_cell=5., threads=1)
    official.write_parquet(R.WORK / "official-de-MLLT1.parquet")
    recorded = pl.read_parquet(R.WORK / "reference_de.parquet").filter(
        pl.col("target") == "MLLT1")
    joined = official.join(recorded, on=["target", "feature"], suffix="_fast")
    changed = joined.filter((pl.col("p_adj") < .05) != (pl.col("p_adj_fast") < .05))
    print("changed calls:", changed, flush=True)
    linear = _to_linear(a, "counts", 1e6)
    ctrl = linear.obs.target_gene.to_numpy() == E.CTRL
    means = np.asarray(linear.X[ctrl].mean(0)).ravel()
    keep = means > 5.
    cs = np.log1p(linear.X[ctrl][:, keep].toarray()).T.copy()
    cs.sort(axis=1)
    pert = np.log1p(linear.X[~ctrl][:, keep].toarray()).T.copy()
    z = rc_emit.wilcoxon_z(cs, rc_emit.ctrl_tiesums(cs), pert)
    p = 2 * norm.sf(np.abs(z))
    q = false_discovery_control(p, method="bh")
    exact = pl.DataFrame({"feature": np.array(a.var_names)[keep],
                          "p_exact": p, "q_exact": q})
    compare = official.join(exact, on="feature")
    report = {
        "old_fast_changed_calls": changed.to_dicts(),
        "old_p_max_abs": float((joined["p_value"]-joined["p_value_fast"]).abs().max()),
        "scanpy_normalized_fast_p_max_abs": float(
            (compare["p_value"]-compare["p_exact"]).abs().max()),
        "scanpy_normalized_fast_q_max_abs": float(
            (compare["p_adj"]-compare["q_exact"]).abs().max()),
        "scanpy_normalized_fast_changed_calls": compare.filter(
            (pl.col("p_adj") < .05) != (pl.col("q_exact") < .05)).height,
        "seconds": time.monotonic()-start,
    }
    R.write_json(R.WORK / "de-parity-diagnosis.json", report)
    print(report, flush=True)


if __name__ == "__main__":
    main()

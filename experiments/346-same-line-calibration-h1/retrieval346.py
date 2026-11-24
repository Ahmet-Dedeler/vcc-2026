#!/usr/bin/env python3
"""346 step 8 — does 341's retrieval endpoint predict the delivered PDS member?

341 was killed by a retrieval guard: its same-line ridge raised weighted sign accuracy but
lowered cosine retrieval against the recipient's own response. This computes the SAME kind of
endpoint at the design stage here — cosine rank retrieval of each arm's designed effect
against the reference's real bulk_lognorm delta, with the panel target genes excluded exactly
as PDS's `exclusion_scope="panel"` does, and the same midrank / (n-1) rank denominator PDS
uses — and lines it up against the delivered `pds` member measured on actual emitted counts.
If they track, 341's guard was a faithful early-warning signal for a priced member.
Design-stage arithmetic plus already-scored per-target values; no emission. Writes to work/.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(BASE), str(ROOT / "src"), str(ROOT / "experiments/334-text-conditioned-transfer")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import run323 as R           # noqa: E402
import evaluate323 as E      # noqa: E402
import core334 as P          # noqa: E402
from cell_eval2.prep import bulk_lognorm_means   # noqa: E402

ARMS = ("cap15", "diag_cross_b1", "same_ridge_b05", "same_ridge_b1", "shuffled_b1")
EXTRA = ("source_mean_b1", "crossline_full_b1", "same_diag_b1")


def main() -> None:
    with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
        targets = np.asarray(z["targets"]).astype(str)
        genes = np.asarray(z["genes"]).astype(str)
        sums = z["sums"]
    logmeans = bulk_lognorm_means(sums.astype(float), E.TS)
    truth = logmeans[:-1] - logmeans[-1]           # real bulk_lognorm delta, 85 x 18,080
    panel = set(targets)
    ranked = np.array([g not in panel for g in genes])   # PDS exclusion_scope="panel"
    with np.load(WORK / "effects346.npz", allow_pickle=False) as z:
        egenes = np.asarray(z["genes"]).astype(str)
        covered = set(egenes[np.asarray(z["covered_columns"], dtype=np.int64)])
        arms = {a: np.asarray(z[a], dtype=np.float64) for a in (*ARMS, *EXTRA) if a in z.files}
    on_covered = np.array([g in covered for g in genes]) & ranked

    delivered = {}
    for a in ARMS:
        table = pl.read_csv(WORK / f"per-target-{a}.csv").filter(
            pl.col("metric") == E.MEMBERS["pds"])
        mapping = dict(zip(table["perturbation"], table["value"]))
        delivered[a] = np.array([mapping[str(t)] for t in targets], dtype=float)

    out = {"n_targets": int(len(targets)),
           "ranked_genes_after_panel_exclusion": int(ranked.sum()),
           "ranked_genes_on_covered_axis": int(on_covered.sum()), "arms": {}}
    for a, effect in arms.items():
        mapped = R.map_axis(effect, egenes, genes)
        full = P.endpoints(mapped[:, ranked], truth[:, ranked])
        cov = P.endpoints(mapped[:, on_covered], truth[:, on_covered])
        row = {"design_retrieval": float(full["retrieval"].mean()),
               "design_own_cosine": float(full["cosine"].mean()),
               "design_decoy_cosine": float(full["decoy_cosine"].mean()),
               "design_retrieval_covered_axis": float(cov["retrieval"].mean()),
               "design_own_cosine_covered_axis": float(cov["cosine"].mean())}
        if a in delivered:
            row["delivered_pds"] = float(delivered[a].mean())
            per_target = full["retrieval"]
            row["design_retrieval_vs_delivered_pds_pearson"] = float(
                np.corrcoef(per_target, delivered[a])[0, 1])
        out["arms"][a] = row
    scored = [a for a in ARMS if a in out["arms"]]
    dr = np.array([out["arms"][a]["design_retrieval"] - out["arms"]["cap15"]["design_retrieval"]
                   for a in scored])
    dp = np.array([out["arms"][a]["delivered_pds"] - out["arms"]["cap15"]["delivered_pds"]
                   for a in scored])
    out["across_arm_agreement"] = {
        "arms": scored,
        "design_retrieval_delta_vs_cap15": dr.tolist(),
        "delivered_pds_delta_vs_cap15": dp.tolist(),
        "pearson": float(np.corrcoef(dr, dp)[0, 1]),
        "slope_delivered_per_design": float(np.polyfit(dr, dp, 1)[0]),
        "signs_agree": int(sum(np.sign(dr) == np.sign(dp))),
    }
    out["definition"] = ("design_retrieval: fraction of the 84 decoy reference deltas that the "
                         "designed effect is farther from than from its own, cosine, midrank "
                         "ties, denominator n-1 — the same rank score PDS reports, evaluated on "
                         "the designed effect instead of the emitted pseudobulk")
    (WORK / "retrieval.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

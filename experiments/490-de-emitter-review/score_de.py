#!/usr/bin/env python3
"""490 — the four DE members on the H1 surface, through cell-eval2's own `prepare_de` + `dispatch_de_metrics`, exactly as
experiments/323-current-family-validation/evaluate323.score does (DE part only). Pure table scoring: no cell-eval2
CLI, no reference bundle, small memory."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
from cell_eval2.config import EvalConfig
from cell_eval2.de import prepare_de
from cell_eval2.run import dispatch_de_metrics

ROOT = Path(__file__).resolve().parents[2]
REAL_DE = ROOT / "experiments/323-current-family-validation/work/reference_de.parquet"
MEMBERS = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
           "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
CTRL = "non-targeting"
_REAL = None


def de_frame(targets, genes, p, q, lfc):
    """evaluate323.de_frame verbatim."""
    return pl.DataFrame({
        "target": np.repeat(targets, len(genes)),
        "feature": np.tile(genes, len(targets)),
        "p_value": np.asarray(p).ravel(),
        "p_adj": np.asarray(q).ravel(),
        "log2_fold_change": np.asarray(lfc).ravel(),
    })


def score(targets, genes_gated, p, q, lfc):
    """Returns (raw means per member, per-target dict member -> {target: value})."""
    global _REAL
    if _REAL is None:
        _REAL = pl.read_parquet(REAL_DE)
    pred = de_frame(np.asarray(targets).astype(str), np.asarray(genes_gated).astype(str), p, q, lfc)
    prepared = prepare_de(pred, _REAL, control=CTRL, sort_by="abs_log2_fold_change",
                          p_adj_threshold=.05, nan_lfc_policy="mask")
    cfg = EvalConfig.from_preset("vcc2026")
    rows = dispatch_de_metrics(list(MEMBERS.values()), prepared, cfg)
    raw, per = {}, {}
    for short, full in MEMBERS.items():
        vals = {row["perturbation"]: row["value"] for row in rows if row["metric"] == full}
        per[short] = vals
        raw[short] = float(np.nanmean(list(vals.values())))
    return raw, per

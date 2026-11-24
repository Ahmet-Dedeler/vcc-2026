#!/usr/bin/env python3
"""Eligibility from existing metadata; target-coordinate mixing without truth."""
import argparse
import json

import numpy as np
import pandas as pd

import run323 as R


def inventory():
    R.verify_inputs()
    depth = pd.read_csv(R.ROOT / "experiments/278-depth-inventory/work/depth278.csv")
    bench = pd.read_csv(R.ROOT / "experiments/285-reliability-benchmark/work/bench285.csv")
    sources = ("k562_gwps", "orion_hct116", "orion_hek293t",
               "cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
    rows = []
    for source in sources:
        d = depth[depth.source == source]
        if d.target.duplicated().any():
            raise ValueError(f"duplicate target metadata: {source}")
        cells = pd.to_numeric(d.cells, errors="coerce")
        rows.append({
            "source": source, "panel_targets": len(d),
            "cells_median": float(cells.median()),
            "at_least_200_cells": int((cells >= 200).sum()),
            "at_least_400_cells": int((cells >= 400).sum()),
            "metadata_source": "278/work/depth278.csv; inherited, not new cells",
        })
    pd.DataFrame(rows).to_csv(R.WORK / "coverage.csv", index=False)
    bench.to_csv(R.WORK / "inherited-reliability.csv", index=False)
    R.write_json(R.WORK / "eligibility.json", {
        "basis": "existing metadata and previously reported measurements",
        "panel_coverage": rows,
        "surfaces": [
            {"recipient": "H1 hESC", "targets": 85, "reference_cells_per_target": 400,
             "use": "current-family union-OFF emitted-count regression panel",
             "exclude": ["h1_2025", "gse295214_h1"],
             "limitations": ["one non-cancer recipient", "different target panel",
                            "cannot test H1 union", "must refit operator for H1 exclusions"]},
            {"recipient": "Orion HCT116 / HEK293T", "targets": "290 / 294",
             "use": "relative direction, not calibrated all-six board selection",
             "exclude": ["whole held line from donors and whitening fit"],
             "limitations": ["few targets reach 400 cells", "documented batch/depth shifts",
                            "remove held line from every inner training example"]},
            {"recipient": "K562-GWPS", "targets": 269,
             "use": "relative direction with held-K562 source exclusion",
             "exclude": ["all K562 screens", "K562-derived PerturbMatch"],
             "limitations": ["sparse reference DE on challenge panel", "low cell budget"]},
            {"recipient": "CD4", "targets": "283-287 by activation state",
             "use": "saved donor/guide reliability statistics; no new cell extraction",
             "exclude": ["whole CD4 source for an independent-recipient test"],
             "limitations": ["saved pseudobulk is not raw cells for six-metric scoring",
                            "conditioning on another donor is not another tissue"]},
            {"recipient": "Nadig Jurkat / HepG2", "targets_on_challenge_panel": 0,
             "use": "non-panel exploratory whole-line/construct development only",
             "limitations": ["essential-target regime mismatch",
                            "no direct 300-panel retrodiction"]},
            {"recipient": "Broad-16 working 12 lines", "targets_on_challenge_panel": 0,
             "use": "exploratory recipient-conditioned response development",
             "exclude": ["JVE127", "MG63", "SLR23", "UMRC3"],
             "limitations": ["KO rather than CRISPRi", "essential-target regime",
                            "cell depth varies by line", "study and line not fully crossed"]},
        ],
        "reserve_access": "NONE; spent/reserved source files not opened",
    })


def mixing():
    """Quantify excluded-input contribution, NOT its predictive benefit."""
    R.verify_inputs()
    b = R.builder()
    genes, targets, _, old_cols, _ = R.axes()
    old_info = json.loads((R.HERE / "inputs.json").read_text())
    expected = old_info["inputs_sha256"]["data/donors/whitener_union_pruned.npz"]
    if R.sha(R.OLD_OPERATOR) != expected:
        raise ValueError("old operator drift")
    op = b.Whitener.load(R.OLD_OPERATOR)
    with np.load(R.WORK / "stage.npz", allow_pickle=False) as z:
        pre = z["uniform"]
    h1_in_old = np.isin(genes[old_cols], targets)
    x = pre[:, old_cols]
    only_targets = np.where(h1_in_old[None, :], x, 0.)
    full = op.apply(x, .65, 0., renorm=False)
    contribution = op.apply(only_targets, .65, 0., renorm=False)
    corrected = op.apply(x-only_targets, .65, 0., renorm=False)
    linearity = float(np.max(np.abs(full-corrected-contribution)))
    if linearity > 1e-10:
        raise AssertionError(f"linearity check failed: {linearity}")
    scored = ~h1_in_old
    numerator = np.linalg.norm(contribution[:, scored], axis=1)
    denominator = np.linalg.norm(full[:, scored], axis=1)
    ratios = numerator / np.maximum(denominator, 1e-300)
    R.write_json(R.WORK / "target-coordinate-mixing.json", {
        "population": "85 H1 target predictions from existing production donors",
        "stage": "historical dense whitening, BEFORE global preserve and emission",
        "statistic": "norm of additive excluded-input contribution / norm of total output",
        "uses_held_response_truth": False,
        "old_operator_sha256": expected, "h1_target_columns_in_old_operator": int(h1_in_old.sum()),
        "linearity_max_abs": linearity,
        "frobenius_ratio": float(np.linalg.norm(contribution[:, scored]) /
                                 np.linalg.norm(full[:, scored])),
        "row_ratio_quantiles": np.percentile(ratios, [0, 10, 50, 90, 100]).tolist(),
        "interpretation": "nonzero target-coordinate mixing; NOT a score gain or leakage of held H1 responses",
    })


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["inventory", "mixing"])
    globals()[ap.parse_args().mode]()

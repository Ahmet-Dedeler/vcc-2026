#!/usr/bin/env python3
"""346 step 1 — recount the target/gene coverage of the same-line calibration regime.

Metadata only: reads perturbation and gene name lists plus the challenge panel and
323's evaluation-target list. No response value is read here, and nothing outside
this experiment's work/ directory is written.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"

# The four donor banks named in the registration (K562-GWPS, both Orion banks,
# Nadig HepG2) plus the two extra CRISPRi banks 334/341 used, for context only.
DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
    "nadig_hepg2": "experiments/166-donor-acquisition/donors/nadig_hepg2",
}
EXTRA = {
    "rpe1": "data/donors/rpe1",
    "nadig_jurkat": "data/donors/nadig_jurkat",
}
PAN = "experiments/115-gse295214-h1/work/compact-donor"   # source id gse295214_h1


def names(path: Path) -> list[str]:
    return path.read_text().splitlines()


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    panel = sorted({r["target_gene"] for r in
                    csv.DictReader((ROOT / "data/challenge/pert_counts.csv").open())})
    base = ROOT / "experiments/323-current-family-validation"
    with np.load(base / "work/reference.npz", allow_pickle=False, mmap_mode="r") as z:
        evaluation = [str(t) for t in z["targets"]]
        scored_genes = [str(g) for g in z["genes"]]          # 18,080 H1 test axis
        keep = np.asarray(z["keep"], dtype=bool)             # 10,779-gene DE universe
    frozen = json.loads((base / "inputs.json").read_text())["h1_targets"]
    assert evaluation == list(frozen), "323 evaluation-target drift"
    canonical = [r["gene_name"] for r in
                 csv.DictReader((ROOT / "data/challenge/gene_names.csv").open())]
    with np.load(base / "work/effects.npz", allow_pickle=False, mmap_mode="r") as z:
        design_genes = [str(g) for g in z["genes"]]          # 18,533 design axis
        assert [str(t) for t in z["targets"]] == evaluation, "effects target order drift"
    assert design_genes == canonical, "effects axis is not the canonical challenge axis"

    pan_perts = names(ROOT / PAN / "perts.txt")
    pan_genes = names(ROOT / PAN / "genes.txt")
    assert pan_genes == canonical, "Pan bank is not on the canonical 18,533-gene axis"
    assert len(set(pan_perts)) == len(pan_perts), "duplicate Pan perturbation row"

    donor_perts = {k: set(names(ROOT / v / "perts.txt")) for k, v in DONORS.items()}
    donor_genes = {k: set(names(ROOT / v / "genes.txt")) for k, v in DONORS.items()}
    extra_perts = {k: set(names(ROOT / v / "perts.txt")) for k, v in EXTRA.items()}
    extra_genes = {k: set(names(ROOT / v / "genes.txt")) for k, v in EXTRA.items()}

    excluded = set(panel) | set(evaluation)
    shared_perts = set.intersection(*donor_perts.values()) & set(pan_perts)
    training = sorted(shared_perts - excluded)

    shared_genes = set.intersection(*donor_genes.values()) & set(canonical)
    axis_no_exclusion = sorted(shared_genes)
    axis = sorted(shared_genes - set(training) - excluded)

    # The measured mask of the Pan bank must cover the whole fitted axis.
    mask = np.load(ROOT / PAN / "measured_mask.npy", mmap_mode="r")
    gi = {g: i for i, g in enumerate(canonical)}
    axis_ix = np.array([gi[g] for g in axis])
    pi = {p: i for i, p in enumerate(pan_perts)}
    rows = np.array([pi[t] for t in training])
    sub = np.asarray(mask[np.ix_(rows, axis_ix)])
    measured_fraction = float(sub.mean())
    fully_measured_rows = int(sub.all(1).sum())

    # Test-side availability: donor rows for the 85 evaluation targets.
    test_in_all_donors = [t for t in evaluation
                          if all(t in donor_perts[k] for k in DONORS)]
    test_in_pan = [t for t in evaluation if t in set(pan_perts)]
    per_donor = {k: sum(t in donor_perts[k] for t in evaluation) for k in DONORS}

    # How much of the scored surface the covered axis actually reaches.
    universe = [g for g, k in zip(scored_genes, keep) if k]
    axis_set = set(axis)
    universe_covered = sum(g in axis_set for g in universe)
    universe_on_design_axis = sum(g in set(canonical) for g in universe)

    out = {
        "registered_donors": list(DONORS),
        "pan_source_id": "gse295214_h1",
        "pan_rows": len(pan_perts),
        "n_panel_targets": len(panel),
        "n_evaluation_targets": len(evaluation),
        "targets": {
            "shared_four_donors_and_pan": len(shared_perts),
            "training_after_exclusions": len(training),
            "next_method_screen_claim": 1178,
            "with_rpe1_and_jurkat_too": len(sorted(
                (shared_perts & set.intersection(*extra_perts.values())) - excluded)),
        },
        "genes": {
            "four_donor_common_on_canonical_axis": len(axis_no_exclusion),
            "after_excluding_training_and_evaluation_targets": len(axis),
            "next_method_screen_claim": 5567,
            "with_rpe1_and_jurkat_too": len(sorted(
                (shared_genes & set.intersection(*extra_genes.values()))
                - set(training) - excluded)),
        },
        "pan_measured_mask": {
            "fraction_of_fitted_coordinates_measured": measured_fraction,
            "training_rows_fully_measured": fully_measured_rows,
            "training_rows": len(training),
        },
        "test_side": {
            "evaluation_targets_present_in_every_registered_donor": len(test_in_all_donors),
            "evaluation_targets_with_a_direct_pan_row": len(test_in_pan),
            "next_method_screen_direct_pan_claim": 30,
            "per_donor_evaluation_rows": per_donor,
            "evaluation_targets_missing_from_some_donor": sorted(
                set(evaluation) - set(test_in_all_donors)),
        },
        "scored_surface": {
            "de_universe_genes": int(keep.sum()),
            "de_universe_genes_inside_covered_axis": universe_covered,
            "de_universe_genes_on_the_design_axis": universe_on_design_axis,
            "design_axis_genes": len(canonical),
            "scored_axis_genes": len(scored_genes),
        },
    }
    (WORK / "coverage.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    # These are the FOUR-donor registered counts, kept as the metadata known answer.
    # The axes actually fitted (three production donors) are written by fit346.py.
    (WORK / "fourdonor-training-targets.txt").write_text("\n".join(training) + "\n")
    (WORK / "fourdonor-covered-genes.txt").write_text("\n".join(axis) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

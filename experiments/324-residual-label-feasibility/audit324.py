"""Metadata-only crossed-design inventory. No expression-data dependencies."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BROWSE = ROOT / "experiments/270-data-reconnaissance/work/gpsadb_metadata_browse.csv"
PANEL = ROOT / "data/challenge/pert_counts.csv"
INVENTORY = ROOT / "experiments/285-reliability-benchmark/work/inv285.json"
UNITS = ROOT / "experiments/268-broad16-context-residuals/work/units268.md"
SOURCES = (BROWSE, PANEL, INVENTORY, UNITS, HERE / "PREREGISTRATION.md",
           HERE / "audit324.py", HERE / "test324.py")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def freeze() -> None:
    out = HERE / "inputs.json"
    if out.exists():
        raise FileExistsError("Do not overwrite a frozen manifest")
    dump(out, {"inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in SOURCES}})


def verify() -> None:
    frozen = json.loads((HERE / "inputs.json").read_text())
    expected_paths = {str(p.relative_to(ROOT)) for p in SOURCES}
    if set(frozen["inputs_sha256"]) != expected_paths:
        raise ValueError("Input allowlist mismatch")
    for path, expected in frozen["inputs_sha256"].items():
        if sha(ROOT / path) != expected:
            raise ValueError(f"Frozen input changed: {path}")
    tracked = [str(p.relative_to(ROOT)) for p in SOURCES
               if p.parent == HERE] + [str((HERE / "inputs.json").relative_to(ROOT))]
    for path in tracked:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    changes = subprocess.check_output(["git", "status", "--porcelain", "--", *tracked],
                                      cwd=ROOT, text=True)
    if changes:
        raise ValueError("Commit the frozen protocol and code before inventory")


def context_key(row: dict, mode: str) -> str | None:
    name = row["CellLineName"].strip()
    accession = row["Cellosaurus"].strip()
    if mode == "raw":
        return name or None
    valid_id = accession.startswith("CVCL_") and ";" not in accession and "," not in accession
    if valid_id:
        return accession
    if mode == "id_only":
        return None
    return f"NAME:{name}" if name else None


def design(rows: list[dict], mode: str, same_method: bool):
    """Unique edges, same-context repeated accessions, and 2x2 rectangles."""
    edges = defaultdict(set)
    for row in rows:
        ctx = context_key(row, mode)
        study = row["GSE"].strip()
        target = row["pbgene"].strip()
        method = row["method"].strip().lower() if same_method else "ALL_METHODS"
        if ctx and study and target and (method or not same_method):
            edges[(target, method, ctx, study)].add(row["Dataset"])
    by_context = defaultdict(dict)
    by_target = defaultdict(dict)
    for (target, method, ctx, study), datasets in edges.items():
        by_context[(target, method, ctx)][study] = datasets
        by_target[(target, method)].setdefault(study, set()).add(ctx)

    pairs = []
    for (target, method, ctx), studies in sorted(by_context.items()):
        for sa, sb in itertools.combinations(sorted(studies), 2):
            pairs.append({"target": target, "method": method, "context": ctx,
                          "accession_a": sa, "accession_b": sb,
                          "datasets_a": "|".join(sorted(studies[sa])),
                          "datasets_b": "|".join(sorted(studies[sb]))})
    rectangles = []
    for (target, method), studies in sorted(by_target.items()):
        for sa, sb in itertools.combinations(sorted(studies), 2):
            for ca, cb in itertools.combinations(sorted(studies[sa] & studies[sb]), 2):
                rectangles.append({
                    "target": target, "method": method, "context_a": ca, "context_b": cb,
                    "accession_a": sa, "accession_b": sb,
                    "datasets_aa": "|".join(sorted(edges[(target, method, ca, sa)])),
                    "datasets_ab": "|".join(sorted(edges[(target, method, ca, sb)])),
                    "datasets_ba": "|".join(sorted(edges[(target, method, cb, sa)])),
                    "datasets_bb": "|".join(sorted(edges[(target, method, cb, sb)])),
                })
    summary = {
        "unique_target_method_context_accession_edges": len(edges),
        "same_context_accession_pairs": len(pairs),
        "targets_with_same_context_replication": len({p["target"] for p in pairs}),
        "replicated_target_method_contexts": len({(p["target"], p["method"], p["context"])
                                                 for p in pairs}),
        "crossed_rectangles": len(rectangles),
        "targets_with_crossed_rectangles": len({r["target"] for r in rectangles}),
        "contexts_with_crossed_rectangles": len({r[k] for r in rectangles
                                                for k in ("context_a", "context_b")}),
    }
    return summary, pairs, rectangles


def write_csv(path: Path, records: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "run"))
    args = parser.parse_args()
    if args.mode == "freeze":
        freeze()
        return
    verify()
    with BROWSE.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    with PANEL.open(newline="") as handle:
        panel = {r["target_gene"] for r in csv.DictReader(handle)}
    panel_rows = [r for r in rows if r["pbgene"] in panel]
    studies = defaultdict(set)
    contexts = defaultdict(set)
    for row in panel_rows:
        studies[row["pbgene"]].add(row["GSE"])
        contexts[row["pbgene"]].add(row["CellLineName"])
    known = {
        "dataset_groups": (len(rows), 7665),
        "accessions": (len({r["GSE"] for r in rows}), 4623),
        "panel_targets": (len({r["pbgene"] for r in panel_rows}), 96),
        "panel_groups": (len(panel_rows), 329),
        "panel_multi_accession_targets": (sum(len(v) >= 2 for v in studies.values()), 50),
        "unique_dataset_ids": (len({r["Dataset"] for r in rows}), len(rows)),
        "panel_size": (len(panel), 300),
    }
    if any(got != expected for got, expected in known.values()):
        dump(WORK / "known-answer-failure.json", known)
        raise ValueError(f"Known-answer mismatch: {known}")
    result = {"known_answer": known, "classification": "metadata feasibility only",
              "input_manifest_sha256": sha(HERE / "inputs.json"),
              "historical_panel_multi_accession_and_multiline_targets":
                  sum(len(s) >= 2 and len(contexts[t]) >= 2 for t, s in studies.items()),
              "metadata_missing_fields": [field for field in
                  ("publication_id", "sample_ids", "control_sample_ids", "construct_id",
                   "guide_sequence", "treatment_duration") if field not in rows[0]],
              "panel_method_groups": dict(Counter(r["method"] for r in panel_rows)),
              "designs": {}}
    WORK.mkdir(exist_ok=True)
    for population, subset in (("all_metadata", rows), ("panel", panel_rows)):
        for mode, same_method in (("raw", False), ("id_fallback", False),
                                  ("id_fallback", True), ("id_only", True)):
            key = f"{population}__{mode}__{'same_method' if same_method else 'mixed_methods'}"
            stats, pairs, rectangles = design(subset, mode, same_method)
            result["designs"][key] = stats
            if population == "panel":
                write_csv(WORK / f"{key}__pairs.csv", pairs,
                          ["target", "method", "context", "accession_a", "accession_b",
                           "datasets_a", "datasets_b"])
                write_csv(WORK / f"{key}__rectangles.csv", rectangles,
                          ["target", "method", "context_a", "context_b", "accession_a",
                           "accession_b", "datasets_aa", "datasets_ab",
                           "datasets_ba", "datasets_bb"])

    inv = json.loads(INVENTORY.read_text())
    broad_lines = sorted(inv["broad16_halves"])
    if set(broad_lines) & {"JVE127", "MG63", "SLR23", "UMRC3"}:
        raise ValueError("Reserved lines in working-line inventory")
    result["single_cell_metadata"] = {
        "k562_two_experiment_target_overlap": inv["k562_study_overlap"]["shared"],
        "k562_two_experiment_panel_overlap": len(inv["k562_study_overlap"]["shared_on_panel"]),
        "k562_two_experiment_publications": 1,
        "broad_working_lines": broad_lines,
        "broad_publications": 1,
        "nadig_jurkat_panel_overlap": inv["nadig_jurkat"]["n_panel"],
        "nadig_hepg2_panel_overlap": inv["nadig_hepg2"]["n_panel"],
        "cd4_donors": inv["cd4_bank280"]["donor_n_unique"],
        "cd4_states": inv["cd4_bank280"]["cond_n_unique"],
        "cd4_publications": 1,
        "source_note": "publication groupings inherited from 268/285 RESULTS, not new source attestation",
    }
    result["verdict"] = "PROVENANCE_UNRESOLVED_NO_RESPONSE_LABEL_PILOT"
    result["six_metrics"] = {m: "UNMEASURED_NO_PREDICTION_CHANGED"
                             for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}
    result["reserve_access"] = False
    dump(WORK / "results.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

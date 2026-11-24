"""Frozen pre-outcome depth sensitivity, separate from the running main code."""

import csv
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.stats import spearmanr
import core332 as C
from run332 import ROOT, HERE, WORK, sha, write

META = ROOT / "data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv"
META_SHA = "1cf07cb227aa18821a6c08942bbcaaeead337a568e49081fa8c644cb27f7b6bf"
LOW = ("C4I", "SNU8", "KMRC20", "U343", "T3M5", "HCT15")
HIGH = ("ONE58", "RPMI7951", "KYSE140", "HKGZCC", "SKGII", "KYSE450")
REQUIRED = ("source_mean", "zero", "no_recipient", "wrong_recipient",
            "shuffled_labels", "absolute_boost")


def main():
    if (WORK / "depth.json").exists():
        raise FileExistsError("depth analysis already terminal")
    files = [HERE / "depth332.py", HERE / "DEPTH-ADDENDUM.md"]
    for p in files:
        subprocess.run(["git", "ls-files", "--error-unmatch", str(p)], check=True,
                       stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *map(str, files)])
    assert sha(META) == META_SHA
    result = json.loads((WORK / "results.json").read_text())
    assert sha(WORK / "per-target.csv") == result["per_target_sha256"]
    with META.open() as f:
        depth = {r["cell_line"]: float(r["total_umis_per_ko"])
                 for r in csv.DictReader(f) if r["cell_line"] in C.LINES}
    assert tuple(sorted(depth, key=depth.get)[:6]) == LOW
    assert tuple(sorted(depth, key=depth.get)[6:]) == HIGH
    with (WORK / "per-target.csv").open() as f:
        rows = list(csv.DictReader(f))
    by_key = {(r["line"], r["target"], r["arm"]): r for r in rows}
    targets = sorted({r["target"] for r in rows})
    def contrast(lines, arm, metric):
        sign = 1 if metric in ("absolute_error", "squared_error") else -1
        return sign * np.array([[float(by_key[(l, t, arm)][metric])
                                 - float(by_key[(l, t, "tabicl_median")][metric])
                                 for t in targets] for l in lines])
    output, gates = {}, {}
    for group, lines in (("lower_depth", LOW), ("higher_depth", HIGH)):
        output[group] = {
            "lines": lines,
            "versus_source_mean": {m: C.paired_interval(contrast(lines, "source_mean", m))
                                   for m in ("absolute_error", "squared_error", "cosine",
                                             "retrieval", "orthogonal_cosine")},
            "absolute_error_comparisons": {a: C.paired_interval(
                contrast(lines, a, "absolute_error")) for a in REQUIRED}}
        for a in REQUIRED:
            gates[f"{group}_beats_{a}_point"] = (
                output[group]["absolute_error_comparisons"][a]["mean"] >= 0)
    delta = contrast(C.LINES, "no_recipient", "absolute_error").mean(axis=1)
    correlation = spearmanr(delta, np.log([depth[l] for l in C.LINES]))
    write(WORK / "depth.json", {
        "strata": output, "gates": gates, "depth_gate_pass": all(gates.values()),
        "combined_status": "EARNS_SEPARATE_EFFECT_STAGE" if (
            result["status"] == "EARNS_SEPARATE_EFFECT_STAGE" and all(gates.values()))
            else "DOES_NOT_EARN_FURTHER_DEVELOPMENT",
        "source_results_sha256": sha(WORK / "results.json"),
        "depth_metadata_sha256": META_SHA,
        "registration_sha256": sha(HERE / "DEPTH-ADDENDUM.md"),
        "script_sha256": sha(Path(__file__)),
        "depth_correlation_descriptive": {"spearman_rho": float(correlation.statistic),
                                         "pvalue": float(correlation.pvalue), "n_lines": 12}})


if __name__ == "__main__":
    main()

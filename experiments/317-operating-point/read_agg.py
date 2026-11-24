"""317 — read the six RAW members out of a cell-eval2 `agg_results.csv`."""
import csv, sys
from pathlib import Path

MEMBERS = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm",
           "nmae": "de_wilcoxon_lfc_nmae",
           "fid": "de_wilcoxon_direction_fidelity_yield_raw",
           "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}


def read(path):
    rows = list(csv.DictReader(Path(path).open()))
    mean = next(r for r in rows if r["statistic"] == "mean")
    return {m: float(mean[c]) for m, c in MEMBERS.items()}


if __name__ == "__main__":
    for p in sys.argv[1:]:
        v = read(p)
        print(Path(p).parts[-3], " ".join(f"{m}={v[m]:.10f}" for m in MEMBERS))

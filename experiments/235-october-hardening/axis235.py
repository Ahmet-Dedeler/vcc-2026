"""235 step 1 — measure the whitening operator's AXIS on every panel we hold.

The operator's column axis is `sorted_by_official_index(union(donor coverage) - panel)`
(`scripts/build_union_whitener.py`, and `production_operator_columns` in the builder).  It is
therefore a FUNCTION OF THE PANEL, and `DEF-RUNBOOK` section 0.1 records that on a new 300-target
panel every check the builder runs still passes while 98.8% of the coordinates pair with the wrong
gene.

This script measures, for each panel on disk:

  * `n_cols`  -- the ONLY quantity the builder currently checks (`make_submission_pm_cd4.py:1248`,
    `build_union_whitener.py:363`).  If it is 17,942 on every panel, the existing check is blind
    by construction and this experiment's guard is necessary, not decorative.
  * `axis_sha256` -- sha256 over the ORDERED column gene symbols, newline-joined UTF-8.  This is
    the quantity the guard compares.
  * the mis-pairing rate: the fraction of positions `i` where the A/B/C axis and the other
    panel's axis carry a DIFFERENT gene.  That is the number 199 reported as 98.8%.

Reads only `genes.txt` from each donor plus each bundle's two CSVs -- no donor matrix is loaded,
so this costs about a second and a few hundred MiB.  The authoritative check runs through the
production `load_production_operator` in `fire235.py`.

    uv run python experiments/235-october-hardening/axis235.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work" / "axis235.json"

DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BUNDLES = {
    "abc": ROOT / "data/challenge",
    "def199_pessimistic": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "def199_realistic": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "def199_badaxis": ROOT / "experiments/199-def-runbook/work/def-bundle-badaxis",
}


def axis_sha256(symbols) -> str:
    """The canonical operator-axis digest.  MUST match both production implementations."""
    return hashlib.sha256("\n".join(map(str, symbols)).encode("utf-8")).hexdigest()


def donor_symbols() -> dict[str, list[str]]:
    return {
        name: (ROOT / "data/donors" / name / "genes.txt").read_text().split("\n")
        for name in DONOR_ORDER
    }


def axis_for(bundle: Path, donors: dict[str, list[str]]) -> tuple[list[str], list[int]]:
    """`production_operator_columns` semantics, from gene lists alone."""
    official = pd.read_csv(bundle / "gene_names.csv")["gene_name"].to_numpy()
    panel = pd.read_csv(bundle / "pert_counts.csv")["target_gene"].to_numpy()
    tindex = {str(g): i for i, g in enumerate(official)}
    union: set[str] = set()
    for genes in donors.values():
        union |= {g for g in map(str, genes) if g in tindex}
    union -= set(map(str, panel))
    columns = sorted(tindex[g] for g in union)
    return [str(official[i]) for i in columns], columns


def main() -> None:
    donors = donor_symbols()
    result: dict[str, object] = {"donors": list(DONOR_ORDER)}
    axes: dict[str, list[str]] = {}
    for tag, bundle in BUNDLES.items():
        if not (bundle / "gene_names.csv").is_file():
            result[tag] = {"present": False, "path": str(bundle)}
            continue
        symbols, columns = axis_for(bundle, donors)
        axes[tag] = symbols
        panel = pd.read_csv(bundle / "pert_counts.csv")["target_gene"].astype(str).tolist()
        result[tag] = {
            "present": True,
            "path": str(bundle.relative_to(ROOT)),
            "panel_size": len(panel),
            "panel_sha256": axis_sha256(panel),
            "n_cols": len(symbols),
            "axis_sha256": axis_sha256(symbols),
            "columns_sha256_int64": hashlib.sha256(
                b"".join(int(c).to_bytes(8, "little") for c in columns)
            ).hexdigest(),
        }

    # --- the silent failure, quantified ------------------------------------------------
    abc = axes["abc"]
    for tag, symbols in axes.items():
        if tag == "abc":
            continue
        pairs = min(len(abc), len(symbols))
        wrong = sum(1 for i in range(pairs) if abc[i] != symbols[i])
        result[tag]["mispaired_vs_abc"] = wrong  # type: ignore[index]
        result[tag]["mispaired_fraction_vs_abc"] = wrong / pairs  # type: ignore[index]
        result[tag]["same_n_cols_as_abc"] = len(symbols) == len(abc)  # type: ignore[index]
        result[tag]["panel_overlap_with_abc"] = len(  # type: ignore[index]
            set(pd.read_csv(BUNDLES[tag] / "pert_counts.csv")["target_gene"].astype(str))
            & set(pd.read_csv(BUNDLES["abc"] / "pert_counts.csv")["target_gene"].astype(str))
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1, sort_keys=True))
    print(json.dumps(result, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()

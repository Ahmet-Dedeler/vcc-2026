"""534 Part B -- operator-axis bookkeeping (PREREGISTRATION-B section 5, third bullet). Light: set operations only.

The operator's axis is `union(the fitted donors' genes on the official axis) - panel`, ordered by official index;
its digest is build_union_whitener.py's `axis_sha256` (sha256 of the newline-joined ordered symbols; plolo's
`axis_sha256`). This recomputes that digest for any (bundle, donor set) exactly as the builder does, so that every R6
operator can be checked against:
  * its own json's `axis_sha256` (the fit is for the panel it claims);
  * A/B/C's `edf6ba88...` (and A/B/C's held-K562 axis `6eaba848...`);
  * every tuned-on panel's axis, for all four donor sets (the 9 panel files the preregistration excludes).
`plolo.Operator._check_axis` for the panel's own fold axis runs inside comp534 / full534 (op.axis_checked recorded).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import hashlib  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DONOR_SETS = {"held_k562_gwps": ("orion_hct116", "orion_hek293t"), "held_orion_hct116": ("k562_gwps", "orion_hek293t"),
              "held_orion_hek293t": ("k562_gwps", "orion_hct116"), "full3": ("k562_gwps", "orion_hct116", "orion_hek293t")}
ABC_AXES = {"edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f": "A/B/C 17,942-col axis (3-donor, held-HCT116, held-HEK293T)",
            "6eaba848c81b77dd1d4b83c0d61b7a8689c1e39762a34da7dcfe26ea1b5d765f": "A/B/C held-K562 axis"}
TUNED_ON = ("data/challenge", "experiments/180-def-readiness/work/def-bundle", "experiments/199-def-runbook/work/def-bundle",
            "experiments/199-def-runbook/work/def-bundle-r", "experiments/199-def-runbook/work/def-bundle-badaxis",
            "experiments/238-essential-free-panel/work/ef-bundle", "experiments/238-essential-free-panel/work/ef-bundle-r",
            "experiments/247-oct22-guards/work/census/collapsed-bundle", "experiments/323-current-family-validation/work/bundle")
_GENES: dict = {}


def donor_genes(name: str) -> list[str]:
    if name not in _GENES:
        _GENES[name] = (ROOT / "data/donors" / name / "genes.txt").read_text().split("\n")
    return _GENES[name]


def axis_sha(bundle: Path, donors) -> tuple[str, int]:
    """build_union_whitener.py's axis digest for (bundle, donors), recomputed."""
    bundle = Path(bundle)
    gfile = bundle / "gene_names.csv" if (bundle / "gene_names.csv").exists() else ROOT / "data/challenge/gene_names.csv"
    target_genes = pd.read_csv(gfile)["gene_name"].to_numpy()
    panel = set(str(t) for t in pd.read_csv(bundle / "pert_counts.csv")["target_gene"])
    tindex = {g: i for i, g in enumerate(target_genes)}
    union = set()
    for n in donors:
        union |= {g for g in donor_genes(n) if g in tindex}
    union -= panel
    wcols = sorted(tindex[g] for g in union)
    syms = [str(target_genes[i]) for i in wcols]
    return hashlib.sha256("\n".join(syms).encode("utf-8")).hexdigest(), len(wcols)


def tuned_on_axes() -> dict[str, str]:
    out = {}
    for b in TUNED_ON:
        for tag, ds in DONOR_SETS.items():
            sha, n = axis_sha(ROOT / b, ds)
            out[sha] = f"{b} x {tag} ({n} cols)"
    return out


def check_ops(panel: str, ops: dict[str, Path]) -> dict:
    """Assert every R6 operator's json axis == the recomputed axis for its own panel, and differs from A/B/C's and from
    every tuned-on panel's (all four donor sets)."""
    bundle = HERE / panel / "bundle"
    tuned = tuned_on_axes()
    rec = {}
    for tag, path in ops.items():
        meta = json.loads(Path(path).with_suffix(".json").read_text())
        sha, n = axis_sha(bundle, DONOR_SETS[tag])
        assert meta["axis_sha256"] == sha, f"{panel}/{tag}: json axis {meta['axis_sha256'][:12]} != recomputed {sha[:12]}"
        assert int(meta["n_cols"]) == n, f"{panel}/{tag}: n_cols {meta['n_cols']} != {n}"
        assert list(meta["donors"]) == list(DONOR_SETS[tag]), f"{panel}/{tag}: donors {meta['donors']}"
        assert sha not in ABC_AXES, f"{panel}/{tag}: axis == {ABC_AXES[sha]}"
        assert sha not in tuned, f"{panel}/{tag}: axis == tuned-on {tuned[sha]}"
        rec[tag] = dict(axis_sha256=sha, n_cols=n, op_sha256=meta["sha256"], differs_from_abc=True,
                        differs_from_tuned_on=f"all {len(tuned)} (panel x donor set) axes")
    return rec


if __name__ == "__main__":
    t = tuned_on_axes()
    print(f"{len(t)} distinct tuned-on axes (9 panel files x 4 donor sets)")
    for k, v in sorted(t.items(), key=lambda kv: kv[1]):
        print(f"  {k[:16]} {v}")

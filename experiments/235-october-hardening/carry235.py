"""235 item B — the COVERAGE-AWARE sizing rule, measured on every source in the registry.

222 put `solve_carrying_row_fraction` in for the H1 union.  233 then measured the four
October-only sources at +5.5% to +17.7% relative cosine each on a disjoint panel while their
**all-rows median is 0.0000 on both panels** -- so under the shipped `solve_cd4_fraction`
(median over EVERY panel row) each of them would be sized at a fraction of **exactly zero** and
deliver nothing, silently.  This measures, for all ten registered sources and all three
panels on disk, both rules side by side:

    carrying-row rule   g_per_unit = median(||added_i|| : ||added_i|| > 0) / base
    shipped all-rows    f_per_unit = median(||added_i|| over ALL rows)      / base

The base is a stand-in whose median row norm is EXACTLY the shipped `h1union` build's
(4.229295604454299) -- there is no D/E/F whitened base to measure today, and the base enters
both rules as the same positive denominator, so the comparison between them is exact and the
solved fractions are directly comparable to the shipped one.

KNOWN-ANSWER ROW (rule 10's corollary): the shipped pair on the A/B/C panel must reproduce
`submissions/222-h1union-g028.h5ad.meta.json` EXACTLY --
    carrying median 6.1654229367116375 | g_per_unit 1.4577895501601277
    solved 0.19207161964444322 at g = 0.28 | 77 carrying rows | all-rows median 0.0

    uv run python experiments/235-october-hardening/carry235.py
"""

from __future__ import annotations

import json
import resource
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import make_submission_pm_cd4 as M  # noqa: E402

PANELS = {
    "abc": ROOT / "data/challenge",
    "def199_pessimistic": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "def199_realistic": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
}
OCTOBER_FOUR = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
SHIPPED_META = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"
SHIPPED_BASE_NORM = 4.229295604454299
G = 0.28


def peak_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def stub_base(n_rows: int, n_genes: int, norm: float = SHIPPED_BASE_NORM) -> np.ndarray:
    """A base whose median row norm is EXACTLY `norm` (one nonzero per row)."""
    base = np.zeros((n_rows, n_genes), dtype=np.float64)
    base[:, 0] = norm
    got = float(np.median(np.linalg.norm(base, axis=1)))
    if got != norm:
        raise AssertionError(f"stub base median row norm {got!r} != {norm!r}")
    return base


def panel_and_genes(bundle: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    panel = tuple(map(str, pd.read_csv(bundle / "pert_counts.csv")["target_gene"]))
    genes = tuple(map(str, pd.read_csv(bundle / "gene_names.csv")["gene_name"]))
    return panel, genes


def coverage(name: str, panel: tuple[str, ...]) -> int:
    path = Path(str(M.UNION_SOURCE_REGISTRY[name]["path"]))
    perts = set((path / "perts.txt").read_text().split("\n"))
    return len(set(panel) & perts)


def size_sources(
    names: tuple[str, ...],
    panel: tuple[str, ...],
    genes: tuple[str, ...],
    *,
    target_g: float = G,
) -> dict[str, object]:
    """Build the union of `names` on `panel` and size it BOTH ways, production code only."""
    prepared, details = [], []
    for name in names:
        bank = M.EXP062.Bank.load(name, Path(str(M.UNION_SOURCE_REGISTRY[name]["path"])))
        source, detail = M.single_source_residual(
            bank, panel, genes, set(panel), chunk_rows=256, max_rss_gib=24.0
        )
        prepared.append(source)
        details.append(detail)
    union = M.union_residual(prepared) if len(prepared) > 1 else prepared[0]
    base = stub_base(len(panel), len(genes))

    # the shipped all-rows rule, through the production function itself
    _, all_rows_frame = M.solve_cd4_fraction(base, union)
    row: dict[str, object] = {
        "sources": list(names),
        "covered_targets": {d["condition"]: int(d["covered_targets"]) for d in details},
        "all_rows_f_per_unit": all_rows_frame["f_per_unit_fraction"],
        "all_rows_would_solve": (
            target_g / all_rows_frame["f_per_unit_fraction"]
            if all_rows_frame["f_per_unit_fraction"] > 0
            else 0.0
        ),
    }
    try:
        fraction, frame = M.solve_carrying_row_fraction(base, union, target_g)
        row["carrying_raised"] = False
        row["solved_fraction"] = fraction
        row["frame"] = frame
    except Exception as exc:  # noqa: BLE001
        row["carrying_raised"] = True
        row["carrying_message"] = str(exc)
        row["solved_fraction"] = 0.0
    return row


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {
        "registry": sorted(M.UNION_SOURCE_REGISTRY),
        "target_g": G,
        "base_median_row_norm": SHIPPED_BASE_NORM,
    }
    checks: list[dict[str, object]] = []

    # --- 1. coverage of every registered source on every panel we hold --------------------
    cov: dict[str, dict[str, int]] = {}
    for tag, bundle in PANELS.items():
        panel, _ = panel_and_genes(bundle)
        cov[tag] = {name: coverage(name, panel) for name in sorted(M.UNION_SOURCE_REGISTRY)}
    result["coverage"] = cov
    print(json.dumps(cov, indent=1))

    # --- 2. KNOWN-ANSWER ROW: the shipped pair on A/B/C ------------------------------------
    shipped = json.loads(SHIPPED_META.read_text())["effects"]["h1_union"]["residual_frame"]
    panel, genes = panel_and_genes(PANELS["abc"])
    ka = size_sources(M.UNION_SOURCES_DEFAULT, panel, genes)
    frame = ka["frame"]
    ka_rows = {
        "residual_median_row_norm_over_carrying_rows": (
            frame["residual_median_row_norm_over_carrying_rows"],
            shipped["residual_median_row_norm_over_carrying_rows"],
        ),
        "g_per_unit_fraction": (
            frame["g_per_unit_fraction"], shipped["g_per_unit_fraction"]
        ),
        "solved_fraction": (frame["solved_fraction"], shipped["solved_fraction"]),
        "carrying_rows": (frame["carrying_rows"], shipped["carrying_rows"]),
        "coverage_frac": (frame["coverage_frac"], shipped["coverage_frac"]),
        "residual_median_row_norm_over_all_rows": (
            frame["residual_median_row_norm_over_all_rows"],
            shipped["residual_median_row_norm_over_all_rows"],
        ),
        "shipped_median_rule_would_solve": (
            frame["shipped_median_rule_would_solve"],
            shipped["shipped_median_rule_would_solve"],
        ),
    }
    for key, (got, want) in ka_rows.items():
        checks.append({
            "check": f"KA h1union/{key}",
            "got": got,
            "want": want,
            "verdict": "PASS" if got == want else "FAIL",
        })
    result["known_answer_h1_pair_abc"] = {k: v[0] for k, v in ka_rows.items()}

    # --- 3. every registered source, sized on every panel it covers ------------------------
    sized: dict[str, dict[str, object]] = {}
    for tag, bundle in PANELS.items():
        panel, genes = panel_and_genes(bundle)
        for name in sorted(M.UNION_SOURCE_REGISTRY):
            if cov[tag][name] == 0:
                continue
            row = size_sources((name,), panel, genes)
            sized[f"{tag}/{name}"] = row
            print(
                f"[{tag}/{name}] covered={cov[tag][name]:3d}  "
                f"carrying_solved={row['solved_fraction']:.6f}  "
                f"all_rows_would_solve={row['all_rows_would_solve']:.6f}",
                flush=True,
            )
            checks.append({
                "check": f"{tag}/{name}: the carrying-row rule sizes it above zero",
                "got": row["solved_fraction"],
                "verdict": "PASS" if float(row["solved_fraction"]) > 0.0 else "FAIL",
            })
            if cov[tag][name] <= len(panel) // 2:
                checks.append({
                    "check": f"{tag}/{name}: the shipped all-rows rule would deliver ZERO",
                    "got": row["all_rows_would_solve"],
                    "verdict": (
                        "PASS" if float(row["all_rows_would_solve"]) == 0.0 else "FAIL"
                    ),
                })
    result["sized"] = sized

    # --- 4. the October four, together, on the October panels ------------------------------
    october: dict[str, object] = {}
    for tag in ("def199_pessimistic", "def199_realistic"):
        panel, genes = panel_and_genes(PANELS[tag])
        present = tuple(n for n in OCTOBER_FOUR if cov[tag][n] > 0)
        if not present:
            october[tag] = {"present": []}
            continue
        row = size_sources(present, panel, genes)
        october[tag] = row
        checks.append({
            "check": f"{tag}: the October four size above zero as a union",
            "got": row["solved_fraction"],
            "verdict": "PASS" if float(row["solved_fraction"]) > 0.0 else "FAIL",
        })
        checks.append({
            "check": f"{tag}: the shipped all-rows rule delivers ZERO for that union",
            "got": row["all_rows_would_solve"],
            "verdict": "PASS" if float(row["all_rows_would_solve"]) == 0.0 else "FAIL",
        })
    result["october_union"] = october

    result["checks"] = checks
    result["peak_rss_gib"] = round(peak_gib(), 3)
    result["all_pass"] = all(c["verdict"] == "PASS" for c in checks)
    (WORK / "carry235.json").write_text(json.dumps(result, indent=1, sort_keys=True))
    for check in checks:
        print(f"[{check['verdict']}] {check['check']}  got={check['got']!r}"
              + (f" want={check['want']!r}" if "want" in check else ""))
    print(f"\nITEM B {'PASS' if result['all_pass'] else 'FAIL'}  "
          f"{len(checks)} checks, peak {result['peak_rss_gib']} GiB")
    raise SystemExit(0 if result["all_pass"] else 1)


if __name__ == "__main__":
    main()

"""247 item C, step 1 — measure CD4's two sizing rules on the real banks, without a build.

`solve_cd4_fraction` sizes CD4 by the median row norm over ALL 300 panel rows.  At coverage
`c` the median lands at relative position `(c - 0.5) / c` of the carrying rows' norm
distribution, so it UNDERSTATES the per-unit norm and the solved fraction OVER-DOSES the
carrying rows by `median_carrying / median_all` (242 #2; 1.4222x measured by 229 at 58.7%).

Both rules depend on the whitened base only through its median row norm, which the shipped
meta records exactly (4.229295604454299), so a stub base with that exact median (235's
`carry235.stub_base`) makes every number below EXACTLY the number the builder would compute --
that is how 235's test suite reproduces the shipped `h1_union` frame to the last bit.

Known-answer rows (all must reproduce `submissions/222-h1union-g028.h5ad.meta.json`):
    residual_median_row_norm_at_unit_fraction  4.464717385709141
    f_per_unit_fraction                        1.0556645369046553
    solved_fraction at VCC_CD4_TARGET_F=0.34   0.32207201067578145
    fully_supported_targets                    276

Then: g_cd4 = 0.34 * median_carrying / median_all on A/B/C, and the carrying rule at g_cd4
must reproduce 0.32207201067578145 to <= 1e-9.  Then both rules on the pessimistic panels.

    uv run python experiments/247-oct22-guards/cd4dose247.py
"""

from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import make_submission_pm_cd4 as M  # noqa: E402

SHIPPED_META = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"
SHIPPED_BASE_NORM = 4.229295604454299
SHIPPED_F = 0.34
PANELS = {
    # tag: (bundle, cd4 banks extracted FOR that panel)
    "abc": (ROOT / "data/challenge", ROOT / "data/external/cd4t/cloud-compact/cd4_compact"),
    "def199_pessimistic": (ROOT / "experiments/199-def-runbook/work/def-bundle",
                           ROOT / "experiments/199-def-runbook/work/cd4_compact_def199"),
    "ef238_pessimistic_essential_free": (
        ROOT / "experiments/238-essential-free-panel/work/ef-bundle",
        ROOT / "experiments/238-essential-free-panel/work/cd4/efp/cd4_compact"),
}


def peak_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def stub_base(n_rows: int, n_genes: int, norm: float = SHIPPED_BASE_NORM) -> np.ndarray:
    base = np.zeros((n_rows, n_genes), dtype=np.float64)
    base[:, 0] = norm
    assert float(np.median(np.linalg.norm(base, axis=1))) == norm
    return base


def panel_and_genes(bundle: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    panel = tuple(map(str, pd.read_csv(bundle / "pert_counts.csv")["target_gene"]))
    genes = tuple(map(str, pd.read_csv(bundle / "gene_names.csv")["gene_name"]))
    return panel, genes


def measure(tag: str, bundle: Path, cd4_root: Path) -> dict[str, object]:
    panel, genes = panel_and_genes(bundle)
    if not (cd4_root / "cd4_rest" / "perts.txt").exists():
        return {"tag": tag, "present": False, "cd4_root": str(cd4_root)}
    t0 = time.time()
    _banks, cd4 = M.prepare_frozen_cd4(cd4_root, panel, genes, panel,
                                       chunk_rows=256, max_rss_gib=24.0)
    added = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    median_all = float(np.median(norms))
    median_carrying = float(np.median(live)) if live.size else 0.0
    base = stub_base(len(panel), len(genes))
    row: dict[str, object] = {
        "tag": tag,
        "present": True,
        "bundle": str(bundle.relative_to(ROOT)),
        "cd4_root": str(cd4_root.relative_to(ROOT)),
        "n_panel": len(panel),
        "fully_supported_targets": int(np.sum(np.asarray(cd4.target_condition_count)
                                              == len(M.CD4_NAMES))),
        "carrying_rows": int(live.size),
        "coverage_frac": float(live.size / norms.size),
        "median_all_rows": median_all,
        "median_carrying_rows": median_carrying,
        "overdose_factor_all_rows_rule": (median_carrying / median_all) if median_all > 0 else None,
        "all_median_percentile_of_carrying": (
            float(np.mean(live <= median_all) * 100.0) if live.size else None),
        "seconds": round(time.time() - t0, 1),
    }
    # the shipped rule, through the production function, at the shipped target
    saved = M.CD4_TARGET_F
    M.CD4_TARGET_F = SHIPPED_F
    try:
        try:
            legacy_fraction, legacy_frame = M.solve_cd4_fraction(base, cd4)
            row["legacy_rule"] = {"raised": False, "solved_fraction": legacy_fraction,
                                  "frame": legacy_frame}
        except Exception as exc:  # noqa: BLE001
            row["legacy_rule"] = {"raised": True, "message": str(exc)}
    finally:
        M.CD4_TARGET_F = saved
    return row | {"_cd4": cd4, "_base": base}


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    shipped = json.loads(SHIPPED_META.read_text())["effects"]["cd4"]
    frame = shipped["residual_frame"]
    checks: list[dict[str, object]] = []

    def chk(name: str, got, want, tol: float = 0.0) -> None:
        ok = (got == want) if tol == 0.0 else (abs(float(got) - float(want)) <= tol)
        checks.append({"check": name, "got": got, "want": want, "tol": tol,
                       "verdict": "PASS" if ok else "FAIL"})

    rows: dict[str, dict[str, object]] = {}
    for tag, (bundle, cd4_root) in PANELS.items():
        rows[tag] = measure(tag, bundle, cd4_root)
        r = rows[tag]
        if not r["present"]:
            print(f"[{tag}] NOT PRESENT: {r['cd4_root']}")
            continue
        print(f"[{tag}] complete {r['fully_supported_targets']}/{r['n_panel']}  "
              f"carrying {r['carrying_rows']} ({r['coverage_frac']:.4f})  "
              f"median_all {r['median_all_rows']!r}  median_carrying {r['median_carrying_rows']!r}  "
              f"overdose {r['overdose_factor_all_rows_rule']!r}  {r['seconds']}s", flush=True)

    # --- 1. KNOWN-ANSWER: the A/B/C frame, to the last bit -------------------------------
    abc = rows["abc"]
    lf = abc["legacy_rule"]["frame"]
    chk("KA abc: residual_median_row_norm_at_unit_fraction",
        lf["residual_median_row_norm_at_unit_fraction"],
        frame["residual_median_row_norm_at_unit_fraction"])
    chk("KA abc: whitened_base_median_row_norm (stub)", lf["whitened_base_median_row_norm"],
        frame["whitened_base_median_row_norm"])
    chk("KA abc: f_per_unit_fraction", lf["f_per_unit_fraction"], frame["f_per_unit_fraction"])
    chk("KA abc: solved_fraction at f=0.34", abc["legacy_rule"]["solved_fraction"],
        frame["solved_fraction"])
    chk("KA abc: fully_supported_targets", abc["fully_supported_targets"],
        shipped["fully_supported_targets"])
    chk("KA abc: carrying rows == fully supported targets", abc["carrying_rows"],
        shipped["fully_supported_targets"])

    # --- 2. g_cd4: the A/B/C DELIVERED carrying-row dose ----------------------------------
    ratio_abc = abc["overdose_factor_all_rows_rule"]
    g_cd4 = SHIPPED_F * ratio_abc
    print(f"\nA/B/C  median_carrying/median_all = {ratio_abc!r}")
    print(f"g_cd4 = 0.34 * {ratio_abc!r} = {g_cd4!r}")
    carrying_fraction, carrying_frame = M.solve_carrying_row_fraction(
        abc["_base"], abc["_cd4"], g_cd4)
    chk("carrying rule at g_cd4 reproduces effects.cd4.fraction (<= 1e-9, NUMERIC not byte)",
        carrying_fraction, frame["solved_fraction"], tol=1e-9)
    chk("carrying frame: residual_median_row_norm_over_all_rows == legacy median",
        carrying_frame["residual_median_row_norm_over_all_rows"],
        frame["residual_median_row_norm_at_unit_fraction"])
    exact = carrying_fraction == frame["solved_fraction"]
    print(f"carrying rule at g_cd4: solved {carrying_fraction!r}  shipped {frame['solved_fraction']!r}"
          f"  |diff| {abs(carrying_fraction - frame['solved_fraction']):.3e}  bit-identical={exact}")

    # --- 3. the pessimistic panels: delivered per-row dose under both rules ----------------
    # Delivered dose on the rows CD4 carries, as a multiple of the base's median row norm:
    #   legacy (all-rows median at F=0.34):  F * median_carrying / median_all = F * overdose
    #   carrying rule at g_cd4:              g_cd4 exactly
    # A/B/C delivered = F * ratio_abc = g_cd4, so "relative to A/B/C's" is the ratio to g_cd4.
    # This is base-independent: the base is the same positive denominator in both rules.
    delivered: dict[str, object] = {}
    for tag, r in rows.items():
        if not r["present"]:
            continue
        ratio = r["overdose_factor_all_rows_rule"]
        legacy_delivered = (SHIPPED_F * ratio) if ratio else None
        cf, cfr = M.solve_carrying_row_fraction(r["_base"], r["_cd4"], g_cd4)
        delivered[tag] = {
            "coverage_frac": r["coverage_frac"],
            "carrying_rows": r["carrying_rows"],
            "overdose_factor": ratio,
            "legacy_rule_raised": r["legacy_rule"]["raised"],
            "legacy_solved_fraction": r["legacy_rule"].get("solved_fraction"),
            "legacy_delivered_dose_x_base": legacy_delivered,
            "legacy_delivered_relative_to_abc": (legacy_delivered / g_cd4) if legacy_delivered else None,
            "carrying_solved_fraction": cf,
            "carrying_delivered_dose_x_base": g_cd4,
            "carrying_delivered_relative_to_abc": 1.0,
            "carrying_frame_shipped_median_rule_would_solve": cfr["shipped_median_rule_would_solve"],
        }
        d = delivered[tag]
        print(f"[{tag}] coverage {d['coverage_frac']:.4f}  legacy delivered {d['legacy_delivered_dose_x_base']!r}"
              f" (x{d['legacy_delivered_relative_to_abc']!r} of A/B/C)  carrying delivered {g_cd4!r} (x1.0)")

    # --- 4. a synthetic residual at 176/300 (229's coverage) as a second carrying-rule KA ---
    # 229's h1_union frame at 58.67%: ratio 1.4222 = 6.439850/4.528221.  A synthetic residual
    # with those exact medians must give that ratio through the production functions.
    rng = np.random.default_rng(247)
    n_rows, n_genes = 300, 64
    synth = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    carrying = rng.choice(n_rows, 176, replace=False)
    norms_target = rng.uniform(3.0, 10.0, size=176)
    for i, row_i in enumerate(carrying):
        v = rng.standard_normal(n_genes)
        synth[row_i] = v / np.linalg.norm(v) * norms_target[i]
        votes[row_i] = 1
    synth_res = M.EXP062.CD4Prepared(mode="synth", mean_effect=synth, gene_votes=votes,
                                     target_condition_count=votes.max(axis=1),
                                     condition_details=[])
    sb = stub_base(n_rows, n_genes)
    s_norms = np.linalg.norm(synth, axis=1)
    s_all, s_car = float(np.median(s_norms)), float(np.median(s_norms[s_norms > 0]))
    _, s_frame = M.solve_carrying_row_fraction(sb, synth_res, 0.35)
    saved = M.CD4_TARGET_F
    M.CD4_TARGET_F = 0.35
    try:
        s_legacy, _ = M.solve_cd4_fraction(sb, synth_res)
    finally:
        M.CD4_TARGET_F = saved
    chk("synthetic 176/300: legacy solve == carrying frame's shipped_median_rule_would_solve",
        s_legacy, s_frame["shipped_median_rule_would_solve"])
    chk("synthetic 176/300: legacy/carrying solved ratio == median_carrying/median_all",
        s_legacy / s_frame["solved_fraction"], s_car / s_all, tol=1e-12)
    print(f"[synthetic 176/300] overdose factor {s_car / s_all!r}, legacy/carrying solve ratio "
          f"{s_legacy / s_frame['solved_fraction']!r}")

    result = {
        "g_cd4": g_cd4,
        "g_cd4_repr": repr(g_cd4),
        "ratio_abc_median_carrying_over_median_all": ratio_abc,
        "abc_median_carrying_rows": abc["median_carrying_rows"],
        "abc_median_all_rows": abc["median_all_rows"],
        "carrying_rule_at_g_cd4_solved_fraction": carrying_fraction,
        "carrying_rule_at_g_cd4_frame": carrying_frame,
        "shipped_solved_fraction": frame["solved_fraction"],
        "abs_diff": abs(carrying_fraction - frame["solved_fraction"]),
        "bit_identical": exact,
        "panels": {tag: {k: v for k, v in r.items() if not k.startswith("_")}
                   for tag, r in rows.items()},
        "delivered": delivered,
        "checks": checks,
        "peak_rss_gib": round(peak_gib(), 3),
        "all_pass": all(c["verdict"] == "PASS" for c in checks),
    }
    (WORK / "cd4dose247.json").write_text(json.dumps(result, indent=1, sort_keys=True, default=str))
    print()
    for c in checks:
        print(f"[{c['verdict']}] {c['check']}  got={c['got']!r} want={c['want']!r}"
              + (f" tol={c['tol']}" if c["tol"] else ""))
    print(f"\nCD4 DOSE {'PASS' if result['all_pass'] else 'FAIL'}  {len(checks)} checks, "
          f"peak {result['peak_rss_gib']} GiB")
    raise SystemExit(0 if result["all_pass"] else 1)


if __name__ == "__main__":
    main()

"""223 Phase 0 — VERIFY THE SILENT-ZERO MECHANISM, OR ABORT THE MANDATE.

Three independent checks, in increasing strength:

  A. STRUCTURAL.  Read `solve_cd4_fraction` out of the production builder's source and
     assert the residual median is taken over EVERY row of `added` (`axis=1` then a bare
     `np.median` with no mask), and that the harness mirror `plolo.solve_fraction` is the
     same expression.  A textual check is not enough on its own, so it is only step one.

  B. BEHAVIOURAL, on the real objects.  Take the real CD4 residual for a real PLOLO fold
     and mask it down to 21% of the fold's targets -- a synthetic partial-coverage source
     with exactly the geometry of `telohaec` (21.7%) or `gse295214_h1` (21.3%).  Then call
     BOTH sizing functions and record what each does.

  C. END-TO-END, bit-identity.  Push the 21%-coverage residual through 162's PLOLO
     instrument at the board incumbent's dose f = 0.34 and compare all endpoints to the
     no-residual baseline.  The mandate's claim is that they are BIT-IDENTICAL.

  KNOWN-ANSWER ROW (rule 10).  The same instrument, same fold, same dose, with the FULL
  CD4 residual (92% coverage) must NOT be bit-identical -- otherwise the comparison is
  measuring a dead instrument rather than a silent zero.  And a 60%-coverage mask (just
  over the 50% line) must also be live.

    uv run python phase0_repro.py
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, Operator, PROD_OPS, build_fold, cd4_for,  # noqa: E402
                   deliver, endpoints, solve_fraction)

ALPHA = 0.982
GAMMA = 0.65
LAM = 0.0
FOLD = "orion_hct116"
DOSE = 0.34
SEED = 20260904
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"


# --------------------------------------------------------------- A. structural
def structural_check() -> dict:
    """Parse the builder and assert the median is over all rows, with no mask."""
    tree = ast.parse(BUILDER.read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "solve_cd4_fraction")
    src = ast.get_source_segment(BUILDER.read_text(), fn)

    # the residual norm line
    line = next(ln.strip() for ln in src.splitlines() if "residual_norm" in ln and "=" in ln)
    assert line == ("residual_norm = float(np.median(np.linalg.norm(added, axis=1)))"), line
    # `added` is the full-height masked matrix, not a subset of rows
    mk = next(ln.strip() for ln in src.splitlines() if ln.strip().startswith("added ="))
    assert mk == "added = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0)", mk
    # no row selection anywhere in the function
    for bad in ("> 0]", "carrying", "nonzero", "[mask", "live"):
        assert bad not in src.replace("gene_votes > 0", ""), f"unexpected row filter: {bad}"

    # the harness mirror
    ph = (ROOT / "experiments/162-retrodicting-validator/plolo.py").read_text()
    ptree = ast.parse(ph)
    pfn = next(n for n in ast.walk(ptree)
               if isinstance(n, ast.FunctionDef) and n.name == "solve_fraction")
    psrc = ast.get_source_segment(ph, pfn)
    assert "rn = float(np.median(np.linalg.norm(added, axis=1)))" in psrc
    assert "target_f / per_unit if per_unit > 0 else 0.0" in psrc

    return dict(
        builder_line=line, builder_mask=mk,
        builder_zero_behaviour="raises ValueError when VCC_CD4_TARGET_F is set",
        harness_zero_behaviour="returns fraction 0.0 SILENTLY",
        verdict="CONFIRMED: median over every panel row, no carrying-row filter")


# ------------------------------------------------------------------ residual ops
def mask_to_coverage(res, frac: float, rng) -> "E062.CD4Prepared":
    """Keep a random `frac` of the CARRYING rows; zero the rest. Coverage-only surgery."""
    votes = np.asarray(res.gene_votes).copy()
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    tcc = np.asarray(res.target_condition_count).copy()
    carry = np.flatnonzero(tcc > 0)
    keep_n = int(round(frac * len(tcc)))
    keep = rng.choice(carry, size=min(keep_n, carry.size), replace=False)
    drop = np.setdiff1d(np.arange(len(tcc)), keep)
    eff[drop] = 0.0
    votes[drop] = 0
    tcc[drop] = 0
    return E062.CD4Prepared(mode=res.mode + f"+cov{frac:.2f}", mean_effect=eff,
                            gene_votes=votes, target_condition_count=tcc,
                            condition_details=res.condition_details)


def carrying_solve(whitened, res, target_g):
    """218 section 7.1 / the builder's `solve_carrying_row_fraction`, mirrored."""
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    cm = float(np.median(live)) if live.size else 0.0
    bn = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    per_unit = cm / bn if bn else 0.0
    return (target_g / per_unit if per_unit > 0 else 0.0), per_unit


def score(final, fold, truth) -> dict:
    acc = {}
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        dl = deliver(final, b, gate, "minnorm")
        e = endpoints(dl[:, fold["scored"]], truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        acc[ctx] = e
    return {k: float(np.mean([acc[c][k] for c in "ABC"])) for k in acc["A"]}


def main() -> None:
    out: dict[str, object] = {}
    print("=" * 78)
    print("A. STRUCTURAL")
    out["structural"] = structural_check()
    for k, v in out["structural"].items():
        print(f"   {k}: {v}")

    print("=" * 78)
    print(f"B/C. BEHAVIOURAL + END-TO-END on PLOLO fold '{FOLD}' at f = {DOSE}")
    fold = build_fold(FOLD)
    targets = fold["targets"]
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(targets)

    op_path = (Path(FOLDS[FOLD]["op"]) if FOLDS[FOLD]["op_kind"] != "production"
               else PROD_OPS / "whitener_union_pruned.npz")
    whitened = Operator(op_path).apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

    base_final = whitened.copy()
    base_final[:, fold["panel_positions"]] = 0.0
    base = score(base_final, fold, truth)
    print(f"   baseline (no residual): cos_agg {base['cos_agg']:.12f} "
          f"pds {base['pds']:.12f}")

    rng = np.random.default_rng(SEED)
    rows = []
    for cov in (0.21, 0.40, 0.60, 1.00):
        res = cd4 if cov >= 1.0 else mask_to_coverage(cd4, cov, rng)
        added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
        norms = np.linalg.norm(added, axis=1)
        carrying = int((norms > 0).sum())
        med_all = float(np.median(norms))
        frac, per_unit = solve_fraction(whitened, res, DOSE)
        gfrac, gper = carrying_solve(whitened, res, DOSE)

        if frac > 0.0:
            final, preserve = E062.add_post_whiten_residual(
                whitened, res, frac, fold["wcols"])
        else:
            final = whitened.copy()
        final[:, fold["panel_positions"]] = 0.0
        e = score(final, fold, truth)

        bit_identical = bool(np.array_equal(final, base_final))
        max_abs = float(np.abs(final - base_final).max())
        d_cos = e["cos_agg"] - base["cos_agg"]
        d_pds = e["pds"] - base["pds"]
        rows.append(dict(
            coverage_requested=cov, carrying_rows=carrying, n_targets=len(targets),
            coverage_measured=carrying / len(targets),
            median_row_norm_all_rows=med_all,
            median_row_norm_carrying_rows=float(np.median(norms[norms > 0]))
            if carrying else 0.0,
            f_per_unit=per_unit, solved_fraction=frac,
            g_per_unit_carrying=gper, solved_fraction_carrying=gfrac,
            cos_agg=e["cos_agg"], pds=e["pds"],
            d_cos_agg=d_cos, d_pds=d_pds,
            effect_matrix_bit_identical_to_baseline=bit_identical,
            max_abs_effect_diff=max_abs))
        print(f"   cov {cov:4.0%} ({carrying:3d}/{len(targets)}) "
              f"med_all {med_all:.6e}  f_per_unit {per_unit:.6e}  frac {frac:.5f}  "
              f"| carrying-row g_per_unit {gper:.4f} frac {gfrac:.5f}  "
              f"| dcos {d_cos:+.3e} dpds {d_pds:+.3e}  "
              f"BIT-IDENTICAL={bit_identical}")
    out["fold"] = FOLD
    out["dose"] = DOSE
    out["baseline"] = base
    out["rows"] = rows

    # ---------------------------------------------------------------- verdicts
    r21 = rows[0]
    r40 = rows[1]
    r60 = rows[2]
    r100 = rows[3]
    checks = {
        "21pct_median_over_all_rows_is_exactly_zero":
            r21["median_row_norm_all_rows"] == 0.0,
        "21pct_f_per_unit_is_exactly_zero": r21["f_per_unit"] == 0.0,
        "21pct_solved_fraction_is_exactly_zero": r21["solved_fraction"] == 0.0,
        "21pct_effect_bit_identical_to_baseline":
            r21["effect_matrix_bit_identical_to_baseline"],
        "21pct_endpoints_bit_identical": (r21["d_cos_agg"] == 0.0
                                          and r21["d_pds"] == 0.0),
        "21pct_carrying_row_rule_is_NONVACUOUS": r21["g_per_unit_carrying"] > 0.0,
        "40pct_also_silently_zero": (r40["solved_fraction"] == 0.0
                                     and r40["effect_matrix_bit_identical_to_baseline"]),
        "KA_60pct_is_LIVE": (r60["solved_fraction"] > 0.0
                             and not r60["effect_matrix_bit_identical_to_baseline"]),
        "KA_full_CD4_is_LIVE": (r100["solved_fraction"] > 0.0
                                and not r100["effect_matrix_bit_identical_to_baseline"]),
    }
    out["checks"] = checks
    print("=" * 78)
    for k, v in checks.items():
        print(f"   {'PASS' if v else 'FAIL'}  {k}")
    out["verdict"] = ("MECHANISM CONFIRMED" if all(checks.values())
                      else "MECHANISM NOT REPRODUCED -- MANDATE VOID")
    print(f"\n   VERDICT: {out['verdict']}")
    (WORK / "phase0_repro.json").write_text(json.dumps(out, indent=2))
    print(f"-> {WORK / 'phase0_repro.json'}")
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

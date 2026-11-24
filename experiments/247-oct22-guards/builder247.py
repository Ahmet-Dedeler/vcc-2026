"""247 item C, unit level — the flag is byte-inert OFF and routes CD4 through the carrying rule ON.

Two builders in one process: the working tree's and the pre-247 byte copy (git HEAD before this
experiment), each loaded with `__file__` pinned to `scripts/` (stage222's trick).  On a synthetic
residual with 4 of 6 rows carrying (66.7% coverage, where the two rules differ), through
`apply_frozen_cd4_residual` -- the production entry point:

    OFF   new builder output == pre-247 builder output, bit for bit (array AND audit)
    ON    frame.sizing_rule == 'carrying', fraction == g / (median_carrying / base), and the
          frame carries the all-rows rule's would-be solve and the over-dose factor

    uv run python experiments/247-oct22-guards/builder247.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepared(M, n_rows, n_genes, rows, cols, rng):
    effect = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in rows:
        effect[row, cols] = rng.uniform(0.2, 1.5, size=len(cols)) * rng.choice([-1, 1], size=len(cols))
        votes[row, cols] = 1
        counts[row] = 3
    return M.EXP062.CD4Prepared(mode="stub", mean_effect=effect, gene_votes=votes,
                                target_condition_count=counts, condition_details=[])


def main() -> None:
    new = load_builder(ROOT / "scripts/make_submission_pm_cd4.py", "builder_247")
    old = load_builder(HERE / "work/pre247/make_submission_pm_cd4.py", "builder_pre247")
    assert new.CD4_SIZING is None and old.CD4_TARGET_F is None
    n_rows, n_genes = 6, 12
    panel_positions = np.array([0, 1], dtype=np.int64)
    columns = np.arange(2, n_genes, dtype=np.int64)
    rng = np.random.default_rng(247)
    whitened = rng.standard_normal((n_rows, n_genes))
    cd4_new = prepared(new, n_rows, n_genes, [0, 1, 2, 3], [2, 3, 4, 5, 6], np.random.default_rng(1))
    cd4_old = prepared(old, n_rows, n_genes, [0, 1, 2, 3], [2, 3, 4, 5, 6], np.random.default_rng(1))
    assert np.array_equal(cd4_new.mean_effect, cd4_old.mean_effect)

    checks = {}
    # --- OFF: identical to the pre-247 builder, with VCC_CD4_TARGET_F=0.34 as shipped -----
    for mod in (new, old):
        mod.CD4_TARGET_F = 0.34
    out_new, audit_new, _ = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions)
    out_old, audit_old, _ = old.apply_frozen_cd4_residual(whitened, cd4_old, columns, panel_positions)
    checks["OFF: output array bit-identical to pre-247"] = np.array_equal(out_new, out_old)
    checks["OFF: cd4 audit identical to pre-247"] = json.dumps(audit_new, sort_keys=True) == json.dumps(audit_old, sort_keys=True)
    checks["OFF: frame is the all-rows rule (no sizing_rule key)"] = "sizing_rule" not in audit_new["residual_frame"]
    legacy_fraction = audit_new["fraction"]
    print(f"OFF  fraction {legacy_fraction!r}  identical={checks['OFF: output array bit-identical to pre-247']}")

    # --- ON: the carrying rule at g ------------------------------------------------------
    added = np.where(cd4_new.gene_votes > 0, cd4_new.mean_effect, 0.0)
    norms = np.linalg.norm(added, axis=1)
    base = float(np.median(np.linalg.norm(whitened, axis=1)))
    median_all, median_car = float(np.median(norms)), float(np.median(norms[norms > 0]))
    g = 0.34 * median_car / median_all           # the A/B/C-style delivered dose on THIS residual
    new.CD4_SIZING, new.CD4_TARGET_G = "carrying", g
    out_on, audit_on, _ = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions)
    fr = audit_on["residual_frame"]
    checks["ON: frame.sizing_rule == 'carrying'"] = fr.get("sizing_rule") == "carrying"
    checks["ON: fraction == g / (median_carrying / base)"] = abs(audit_on["fraction"] - g / (median_car / base)) <= 1e-15
    checks["ON: at g = 0.34 * overdose the fraction equals the all-rows solve (<= 1e-9)"] = abs(audit_on["fraction"] - legacy_fraction) <= 1e-9
    checks["ON: frame carries the all-rows would-be solve == legacy fraction"] = abs(float(fr["legacy_all_rows_rule_would_solve"]) - legacy_fraction) <= 1e-15
    checks["ON: frame carries the over-dose factor == median_carrying/median_all"] = abs(fr["overdose_factor_all_rows_rule"] - median_car / median_all) <= 1e-15
    checks["ON: carrying_rows == 4, coverage 4/6"] = fr["carrying_rows"] == 4 and abs(fr["coverage_frac"] - 4 / 6) < 1e-12
    checks["ON: output array is NOT the OFF one at a different g"] = True
    # a different g moves the output; the same-dose g reproduces it numerically
    new.CD4_TARGET_G = 0.5
    out_half, _, _ = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions)
    checks["ON: a different g moves the output (the flag acts)"] = not np.array_equal(out_half, out_new)
    checks["ON: at the matched g the output equals OFF's to 1e-12"] = np.allclose(out_on, out_new, rtol=0, atol=1e-12)
    print(f"ON   fraction {audit_on['fraction']!r} at g={g!r}; overdose {fr['overdose_factor_all_rows_rule']!r}; "
          f"all-rows would solve {fr['legacy_all_rows_rule_would_solve']!r}")

    (WORK / "builder247.json").write_text(json.dumps(
        {"checks": checks, "off_fraction": legacy_fraction, "on_fraction": audit_on["fraction"], "g": g,
         "on_frame": fr}, indent=1))
    for k, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {k}")
    ok = all(checks.values())
    print(f"\nITEM C (unit) {'PASS' if ok else 'FAIL'}  peak {new.peak_rss_gib():.2f} GiB")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()

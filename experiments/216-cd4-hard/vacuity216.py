"""216 Phase 3 gate — non-vacuity, asserted on every surface before any member is read.

Rule 12 and the mandate both require that an unmeasured member come with the MECHANISM that
pins it, and 200/207 found the mechanism here: production sizes the CD4 residual by a MEDIAN
row norm over every row, and on `h1-test` 81 of 85 rows are zero, so the median is 0, so the
solved fraction is 0 and the residual is IDENTICALLY zero.  That is re-measured here rather
than inherited.

Three surfaces, three questions:
  (a) h1-test, RESIDUAL form  -- expected: exactly vacuous (reproduce 207)
  (b) h1-test, LIBRARY form   -- NOT sized by a median, so it is NOT vacuous; measured, since
      the library arm is one of the three forms this experiment prices and the exemption must
      be claimed per form, not per experiment
  (c) the four PLOLO folds at EVERY dose this experiment reports, including the new tail

    uv run python vacuity216.py
"""
from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)
from sweep216 import ALPHA, GAMMA, LAM  # noqa: E402

DOSES = [0.05, 0.10, 0.20, 0.34, 0.40, 0.55, 0.80, 1.00, 2.00]
H1 = ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz"


def main() -> None:
    out: dict[str, object] = {}

    # ------------------------------------------------- (a)/(b) h1-test
    E = np.load(H1, allow_pickle=True)
    h1 = tuple(str(t) for t in E["targets"])
    cd4h = cd4_for(h1)
    votes = np.asarray(cd4h.gene_votes)
    eff = np.asarray(cd4h.mean_effect)
    added = np.where(votes > 0, eff, 0.0)
    rn = np.linalg.norm(added, axis=1)
    carry = int((np.asarray(cd4h.target_condition_count) == 3).sum())
    med = float(np.median(rn))
    panel = set(line.split(",")[0] for line in
                (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
                if line.strip())
    out["h1_test_residual"] = dict(
        n_targets=len(h1), targets_in_panel=len(set(h1) & panel),
        carrying_all_three=carry, nonzero_rows=int((rn > 0).sum()),
        median_row_norm=med, mean_row_norm=float(rn.mean()), max_row_norm=float(rn.max()),
        f_per_unit_is_zero=bool(med == 0.0),
        solved_fraction=0.0 if med == 0.0 else None,
        vacuous=bool(med == 0.0),
        mechanism="production `solve_cd4_fraction` sizes the residual by the MEDIAN row "
                  "norm over all rows; with 81/85 rows zero the median is 0, f_per_unit is "
                  "0 and the solved fraction is 0 -- the residual is identically zero, not "
                  "small.  Production RAISES on this input; `plolo.solve_fraction` returns "
                  "0.0 silently (207 §6.1 handed that back).")
    out["h1_test_library"] = dict(
        rows_the_library_form_changes=int((rn > 0).sum()),
        fraction_of_rows=float((rn > 0).sum() / len(h1)),
        vacuous=bool((rn > 0).sum() == 0),
        note="the library form mixes CD4 into `pre` directly and is NOT sized by a median, "
             "so it is non-vacuous on h1-test -- but on 4 of 85 rows only, which is 4.7% of "
             "the surface.  Any h1-test member read for the library arm is therefore a "
             "4-row measurement and is reported as unusable rather than as a number.")
    print(f"(a) h1-test residual: {len(h1)} targets, {out['h1_test_residual']['targets_in_panel']} "
          f"in panel, carrying {carry}, nonzero rows {int((rn>0).sum())}, "
          f"median row norm {med!r} -> VACUOUS={med == 0.0}")
    print(f"(b) h1-test library: changes {int((rn>0).sum())}/85 rows "
          f"({100*(rn>0).sum()/len(h1):.1f}%)")

    # ------------------------------------------------- (c) PLOLO, every dose
    rows: dict[float, list] = {}
    per_fold = {}
    for held in FOLDS:
        fold = build_fold(held)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        a = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0)
        carry = int((np.asarray(cd4.target_condition_count) == 3).sum())
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        w = fold["wcols"]
        wh = op.apply(pre * ALPHA, w, GAMMA, LAM)
        del op
        gc.collect()
        base = wh.copy()
        base[:, fold["panel_positions"]] = 0.0
        nb = float(np.linalg.norm(base))
        med_p = float(np.median(np.linalg.norm(a, axis=1)))
        base_med = float(np.median(np.linalg.norm(wh, axis=1)))
        for v in DOSES:
            frac, per_unit = solve_fraction(wh, cd4, v)
            fin, _ = E062.add_post_whiten_residual(wh, cd4, frac, w)
            fin[:, fold["panel_positions"]] = 0.0
            rows.setdefault(v, []).append(
                (frac, float(np.linalg.norm(fin - base)) / nb))
            del fin
        per_fold[held] = dict(n_targets=len(fold["targets"]), carrying=carry,
                              median_residual_row_norm=med_p,
                              median_whitened_row_norm=base_med, nonvacuous=bool(med_p > 0))
        print(f"  {held:<15} targets {len(fold['targets']):3d} carrying {carry:3d} "
              f"median residual row norm {med_p:.4f} (base {base_med:.4f}) NONZERO={med_p>0}",
              flush=True)
        del wh, base, pre, cd4, fold, a
        gc.collect()
    out["plolo_folds"] = per_fold
    disp = {}
    print(f"\n  {'f':>6}{'solved fraction':>18}{'||delivered(f)-delivered(0)||/||delivered(0)||':>50}")
    for v in DOSES:
        fr = float(np.mean([x[0] for x in rows[v]]))
        dz = float(np.mean([x[1] for x in rows[v]]))
        disp[f"{v:g}"] = dict(fraction=fr, rel_displacement=dz)
        print(f"  {v:>6.2f}{fr:>18.5f}{dz:>50.5f}")
    out["displacement"] = disp
    out["all_plolo_nonvacuous"] = bool(all(v["nonvacuous"] for v in per_fold.values()))
    out["min_displacement"] = min(d["rel_displacement"] for d in disp.values())
    (WORK / "vacuity216.json").write_text(json.dumps(out, indent=2))
    print(f"\nall PLOLO folds non-vacuous: {out['all_plolo_nonvacuous']}; "
          f"smallest displacement reported anywhere: {out['min_displacement']:.5f}")


if __name__ == "__main__":
    main()

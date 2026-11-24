"""207 — non-vacuity, both directions (INBOX-1, verified rather than accepted).

(a) h1-test: reproduce 200's claim that production's median-based sizing rule makes the CD4
    residual IDENTICALLY zero there, not merely small.
(b) PLOLO: prove this experiment's own instrument is NOT in that state, at every dose, by a
    stated norm on the delivered effect.

    uv run python vacuity207.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, solve_fraction  # noqa: E402
from sweep207 import ALPHA, GAMMA, LAM  # noqa: E402

out = {}

# ---------------------------------------------------------------- (a) h1-test is vacuous
E = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
            allow_pickle=True)
h1 = tuple(str(t) for t in E["targets"])
cd4h = cd4_for(h1)
added = np.where(np.asarray(cd4h.gene_votes) > 0, np.asarray(cd4h.mean_effect), 0.0)
rn = np.linalg.norm(added, axis=1)
carry = int((np.asarray(cd4h.target_condition_count) == 3).sum())
med = float(np.median(rn))
print(f"(a) h1-test: {len(h1)} targets, CD4 carries all-3-conditions on {carry}, "
      f"nonzero row norms {int((rn>0).sum())}")
print(f"    MEDIAN row norm of the residual = {med!r}  ->  solve_cd4_fraction f_per_unit = "
      f"{med/1.0 if med else 0.0!r}  ->  solved fraction = "
      f"{'0.0 (SILENT NO-OP)' if med == 0 else 'nonzero'}")
print(f"    mean row norm (for contrast) = {float(rn.mean()):.6f}; "
      f"max = {float(rn.max()):.6f}")
out["h1_test"] = dict(n_targets=len(h1), carrying=carry, nonzero_rows=int((rn > 0).sum()),
                      median_row_norm=med, mean_row_norm=float(rn.mean()),
                      vacuous=bool(med == 0.0))

# ------------------------------------------------------------- (b) PLOLO is NOT vacuous
print("\n(b) PLOLO folds -- solved fraction and delivered-effect displacement vs f=0")
DOSES = [0.05, 0.10, 0.20, 0.34, 0.40, 0.80]
rows = {}
for held in FOLDS:
    fold = build_fold(held, limit_gib=7.0)
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(fold["targets"])
    a = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0)
    carry = int((np.asarray(cd4.target_condition_count) == 3).sum())
    op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                  if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
    w = fold["wcols"]
    wh = op.apply(pre * ALPHA, w, GAMMA, LAM)
    base = wh.copy(); base[:, fold["panel_positions"]] = 0.0
    nb = float(np.linalg.norm(base))
    med = float(np.median(np.linalg.norm(a, axis=1)))
    for v in DOSES:
        frac, per_unit = solve_fraction(wh, cd4, v)
        fin, _ = E062.add_post_whiten_residual(wh, cd4, frac, w)
        fin[:, fold["panel_positions"]] = 0.0
        rows.setdefault(v, []).append((frac, float(np.linalg.norm(fin - base)) / nb))
    print(f"  {held:<15} targets {len(fold['targets']):3d}  CD4 carries {carry:3d}  "
          f"median residual row norm {med:.4f}  (NONZERO: {med > 0})")
    out.setdefault("plolo", {})[held] = dict(n_targets=len(fold["targets"]), carrying=carry,
                                             median_row_norm=med)
    del op
print(f"\n  {'f':>6}{'solved fraction':>18}{'||delivered(f) - delivered(0)|| / ||delivered(0)||':>52}")
for v in DOSES:
    fr = float(np.mean([x[0] for x in rows[v]]))
    dz = float(np.mean([x[1] for x in rows[v]]))
    print(f"  {v:>6.2f}{fr:>18.5f}{dz:>52.5f}")
    out.setdefault("displacement", {})[f"{v:g}"] = dict(fraction=fr, rel_displacement=dz)
(WORK / "vacuity207.json").write_text(json.dumps(out, indent=2))
print(f"\nwrote {WORK/'vacuity207.json'}")

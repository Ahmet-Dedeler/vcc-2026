"""218 Phase 0d — the CONFOUND control. Alignment must be compared on MATCHED targets.

166 section 4.2 measured the trap this closes: between the SAME two cell lines, the choice of
target set moves the transfer cosine by **x3.98**, while the choice of line moves it by
**x1.36**.  So a source that carries 61 of a fold's 290 targets cannot have its
`cos_truth` compared to a source that carries 282 of them -- the target set, not the
biology, would be doing the work.  Without this control 218 would have reported that Pan's
H1 hESC screen aligns 3x better with our recipients than CD4 does.

This writes, per fold, the PER-TARGET cosine of every source to the held line's truth, so any
subset comparison is post hoc:

  cos_truth_matched[S] : source S's mean cosine on ITS OWN carrying rows
  cd4_on_S             : CD4's mean cosine ON THE SAME ROWS
  cons_on_S            : the production consensus's mean cosine on the same rows
  ratio_vs_cd4         : cos_truth_matched[S] / cd4_on_S       <- the only comparable number

    uv run python matched218.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import FOLDS, build_fold, cd4_for  # noqa: E402
from resid218 import bank_for, multi_source_residual  # noqa: E402

CANDIDATES = ("kolf2", "feng_ipsc", "k562_gwps", "orion_hct116", "orion_hek293t",
              "perturbmatch", "telohaec", "gse295214_h1", "h1_2025", "myllia_a549",
              "frangieh_melanoma")
SEED = 20260904


def per_target_cos(S, T):
    """cos(S_i, T_i) per row; NaN where S_i is all zero."""
    ns = np.linalg.norm(S, axis=1)
    nt = np.linalg.norm(T, axis=1)
    out = np.full(S.shape[0], np.nan)
    ok = (ns > 0) & (nt > 0)
    out[ok] = (S[ok] * T[ok]).sum(1) / (ns[ok] * nt[ok])
    return out


def main() -> None:
    rows = []
    for held in FOLDS:
        fold = build_fold(held)
        targets = list(fold["targets"])
        sc = fold["scored"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        pre = np.asarray(fold["pre"], dtype=np.float64)[:, sc]
        cd4 = cd4_for(targets)
        cd4_eff = np.where(np.asarray(cd4.gene_votes) > 0,
                           np.asarray(cd4.mean_effect), 0.0)[:, sc]
        cd4_pt = per_target_cos(cd4_eff, truth)
        cons_pt = per_target_cos(pre, truth)

        srcs = {"cd4": cd4_pt, "consensus": cons_pt}
        carry = {"cd4": np.asarray(cd4.target_condition_count) > 0,
                 "consensus": np.ones(len(targets), dtype=bool)}
        for name in CANDIDATES:
            if name == held:
                continue
            try:
                r = multi_source_residual([bank_for(name)], targets, "gene_norm")
            except Exception as exc:                      # noqa: BLE001
                print(f"  [{held}] skip {name}: {exc}")
                continue
            eff = np.where(np.asarray(r.gene_votes) > 0,
                           np.asarray(r.mean_effect), 0.0)[:, sc]
            srcs[name] = per_target_cos(eff, truth)
            carry[name] = np.asarray(r.target_condition_count) > 0

        print(f"[{held}] {len(targets)} targets")
        print(f"  {'source':16s} {'n':>4s} {'cos(S)':>9s} {'cd4 same rows':>14s} "
              f"{'cons same rows':>15s} {'S/cd4':>7s} {'S/cons':>7s}")
        for name, pt in srcs.items():
            m = carry[name] & np.isfinite(pt)
            if m.sum() < 5:
                continue
            s = float(np.nanmean(pt[m]))
            c = float(np.nanmean(cd4_pt[m & np.isfinite(cd4_pt)]))
            k = float(np.nanmean(cons_pt[m & np.isfinite(cons_pt)]))
            rng = np.random.default_rng(SEED)
            idx = np.flatnonzero(m)
            boots = [float(np.nanmean(pt[rng.choice(idx, len(idx))])) for _ in range(2000)]
            rows.append(dict(fold=held, source=name, n=int(m.sum()),
                             cos_source=s, cos_cd4_same_rows=c, cos_cons_same_rows=k,
                             ratio_vs_cd4=(s / c if c else None),
                             ratio_vs_consensus=(k and s / k),
                             ci95=[float(np.percentile(boots, 2.5)),
                                   float(np.percentile(boots, 97.5))],
                             per_target=np.where(np.isfinite(pt), pt, None).tolist(),
                             carrying=m.tolist()))
            print(f"  {name:16s} {int(m.sum()):4d} {s:+9.5f} {c:+14.5f} {k:+15.5f} "
                  f"{s/c if c else float('nan'):7.3f} {s/k if k else float('nan'):7.3f}")
        rows[-1]["targets"] = targets

    (WORK / "matched218.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {WORK/'matched218.json'}")


if __name__ == "__main__":
    main()

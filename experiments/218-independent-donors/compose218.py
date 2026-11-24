"""218 Phase 2 — do TWO independent residuals COMPOSE, on top of the shipped CD4 dose?

Phase 0 measures each source alone.  The acquisition question is different and it is the one
the `1/k` law actually asks: does adding a SECOND independent source on top of the incumbent
buy anything the incumbent does not already have?  If independence is the mechanism, a second
maximally-independent source should compose; if the mechanism is transfer alignment, a weaker
partner adds noise at any dose.

FORM.  Both residuals are solved against the SAME denominator -- the whitened base -- so `f`
means the same thing for both, and then added in ONE step followed by production's single
global preserve rescale.  That is exactly `add_post_whiten_residual` with a merged residual
whose `mean_effect` is already `frac_cd4 * cd4 + frac_x * x` and `fraction = 1.0`.
The CD4 dose is NEVER swept here (that is 216's and 217's territory); it is pinned at the
board incumbent `cd4f034`'s `VCC_CD4_TARGET_F = 0.34`, and also reported at the older
shipped 0.10 so the rider can be re-priced against either incumbent (208's rule).

KNOWN-ANSWER ROWS -- the run aborts if either fails:
  KA-A  `g = 0` must reproduce `resid218`'s `cd4_gn0.34` bit-for-bit (same object, two paths).
  KA-B  `g = 0` at cd4 f = 0.10 must reproduce 207's `f0.1`.
FALSIFIERS:
  * a target-shuffled second source at every `g` (whole rows permuted by a derangement);
  * the leak-free `kolf2` fold, and no source is ever used on the fold where it is the truth;
  * a per-arm non-vacuity assertion (median residual row norm and carrying-row count).

    uv run python compose218.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)
from resid218 import (ALPHA, GAMMA, LAM, bank_for, multi_source_residual,  # noqa: E402
                      score_effects, shuffled, source_residual)

CD4_F = (0.34, 0.10)                       # board incumbent, and the older shipped dose
# Sources whose fold coverage clears production's MEDIAN sizing rule (>50%) get the
# production solve; the two below it are re-sized on their CARRYING rows instead, because
# under the shipped rule they are identically zero (200's silent-zero, verified in 207 and
# reproduced here in `resid218`).  CD4's own sizing is never touched -- that is 216's.
SECOND = ("feng_ipsc", "kolf2")
SECOND_CARRY = ("gse295214_h1", "h1_2025")
G = (0.0, 0.02, 0.05, 0.10, 0.17, 0.28, 0.45)


def solve_fraction_carrying(whitened, res, target_f):
    """`f` as a ratio of the source's typical CARRYING-row norm to the base's typical norm.

    Production's `solve_cd4_fraction` takes the median over EVERY row, so a source covering
    <=50% of the panel has median 0 and the solve returns 0 -- the residual is not small, it
    is identically absent.  For a partial-coverage source the meaningful denominator is the
    rows it actually carries.  Reported separately and never applied to CD4.
    """
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    rn = float(np.median(live)) if live.size else 0.0
    bn = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    per_unit = rn / bn if bn else 0.0
    return (target_f / per_unit if per_unit > 0 else 0.0), per_unit


def union_rows(sources):
    """Per-gene mean over whichever sources carry a coordinate; union of carried rows."""
    tot = np.zeros_like(np.asarray(sources[0].mean_effect), dtype=np.float64)
    votes = np.zeros(tot.shape, dtype=np.uint16)
    counts = np.zeros(np.asarray(sources[0].target_condition_count).shape, dtype=np.uint16)
    for s in sources:
        v = np.asarray(s.gene_votes) > 0
        tot[v] += np.asarray(s.mean_effect)[v]
        votes[v] += 1
        counts += (np.asarray(s.target_condition_count) > 0).astype(np.uint16)
    eff = np.zeros_like(tot)
    live = votes > 0
    eff[live] = tot[live] / votes[live]
    return E062.CD4Prepared(mode="union:" + "+".join(s.mode for s in sources),
                            mean_effect=eff, gene_votes=np.minimum(votes, 255).astype(np.uint8),
                            target_condition_count=np.minimum(counts, 255).astype(np.uint8),
                            condition_details=[])


def restrict_rows(src, like):
    """`src`'s residual, kept only on the rows `like` carries (row-set-matched control)."""
    keep = np.asarray(like.target_condition_count) > 0
    eff = np.asarray(src.mean_effect).copy()
    votes = np.asarray(src.gene_votes).copy()
    counts = np.asarray(src.target_condition_count).copy()
    eff[~keep] = 0.0
    votes[~keep] = 0
    counts[~keep] = 0
    return E062.CD4Prepared(mode=f"{src.mode}|rows_of_{like.mode}", mean_effect=eff,
                            gene_votes=votes, target_condition_count=counts,
                            condition_details=[])


def merged(cd4, frac_cd4, x, frac_x):
    """One combined add: `frac_cd4 * cd4 + frac_x * x`, then production's single preserve."""
    vc = np.asarray(cd4.gene_votes) > 0
    ec = np.where(vc, np.asarray(cd4.mean_effect), 0.0) * frac_cd4
    if x is None or frac_x == 0.0:
        eff, votes, counts = ec, np.asarray(cd4.gene_votes).copy(), \
            np.asarray(cd4.target_condition_count).copy()
    else:
        vx = np.asarray(x.gene_votes) > 0
        ex = np.where(vx, np.asarray(x.mean_effect), 0.0) * frac_x
        eff = ec + ex
        votes = (vc | vx).astype(np.uint8)
        counts = np.maximum(np.asarray(cd4.target_condition_count),
                            np.asarray(x.target_condition_count)).astype(np.uint8)
    return E062.CD4Prepared(mode="merged", mean_effect=eff, gene_votes=votes,
                            target_condition_count=counts, condition_details=[])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "compose218.json"))
    args = ap.parse_args()

    ref = {}
    for r in json.loads((WORK / "resid218.json").read_text()):
        ref[(r["fold"], r["arm"])] = r
    ka207 = {}
    for r in json.loads((ROOT / "experiments/207-cd4-dose/work/sweep207.json").read_text()):
        ka207[(r["fold"], r["arm"])] = r

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held)
        targets = fold["targets"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

        seconds = {}
        sizer = {}
        for name in SECOND:
            if name == held:
                continue
            seconds[name] = source_residual(name, targets, "gene_norm")
            seconds[f"{name}_SHUF"] = shuffled(seconds[name])
            sizer[name] = sizer[f"{name}_SHUF"] = solve_fraction
        for name in SECOND_CARRY:
            if name == held:
                continue
            tag = f"{name}CARRY"
            seconds[tag] = source_residual(name, targets, "gene_norm")
            seconds[f"{tag}_SHUF"] = shuffled(seconds[tag])
            sizer[tag] = sizer[f"{tag}_SHUF"] = solve_fraction_carrying
            # THE CONTROL THAT MATTERS: the same rows, the same sizing, but CD4's OWN
            # effect instead of the new source's.  If a partial-coverage arm is really just
            # "more residual on a favourable subset of rows", this reproduces it and the new
            # source is worth nothing.  It also absorbs any row-selection effect, because the
            # row set is identical by construction.
            ctag = f"{name}CARRYCD4"
            seconds[ctag] = restrict_rows(cd4, seconds[tag])
            sizer[ctag] = solve_fraction_carrying
        # The two H1 sources are different screens of the same cell type and their panel
        # target sets are almost disjoint (64 and 20, sharing 7 -- 115's README), so their
        # UNION covers 77/300 rather than 64.  Averaged where both carry a target.
        if all(f"{n}CARRY" in seconds for n in SECOND_CARRY):
            seconds["h1PAIRCARRY"] = union_rows(
                [seconds[f"{n}CARRY"] for n in SECOND_CARRY])
            seconds["h1PAIRCARRY_SHUF"] = shuffled(seconds["h1PAIRCARRY"])
            seconds["h1PAIRCARRYCD4"] = restrict_rows(cd4, seconds["h1PAIRCARRY"])
            for t in ("h1PAIRCARRY", "h1PAIRCARRY_SHUF", "h1PAIRCARRYCD4"):
                sizer[t] = solve_fraction_carrying

        for cf in CD4_F:
            frac_cd4, _ = solve_fraction(whitened, cd4, cf)
            for sname, sres in [("none", None)] + sorted(seconds.items()):
                for g in G:
                    if sname == "none" and g != 0.0:
                        continue
                    if sname != "none" and g == 0.0:
                        continue
                    arm = f"cd4{cf:g}+{sname}{g:g}"
                    if (held, arm) in done:
                        continue
                    t0 = time.time()
                    if sres is None or g == 0.0:
                        frac_x = 0.0
                        med_x, carry_x = 0.0, 0
                    else:
                        frac_x, _ = sizer[sname](whitened, sres, g)
                        added = np.where(np.asarray(sres.gene_votes) > 0,
                                         np.asarray(sres.mean_effect), 0.0)
                        norms = np.linalg.norm(added, axis=1)
                        med_x = float(np.median(norms[norms > 0])) if (norms > 0).any() else 0.0
                        carry_x = int(np.sum(np.asarray(sres.target_condition_count) > 0))
                        assert med_x > 0.0 and frac_x > 0.0, \
                            f"VACUOUS second source {sname} on {held}"
                    res = merged(cd4, frac_cd4, sres, frac_x)
                    final, preserve = E062.add_post_whiten_residual(
                        whitened, res, 1.0, fold["wcols"])
                    final[:, fold["panel_positions"]] = 0.0
                    e = score_effects(final, fold, truth)
                    rows.append(dict(
                        fold=held, arm=arm, cd4_f=cf, second=sname, g=g,
                        n_targets=len(targets), **e,
                        audit=dict(frac_cd4=float(frac_cd4), frac_x=float(frac_x),
                                   median_second_row_norm=med_x, carrying_second=carry_x,
                                   preserve_scale=float(
                                       preserve["post_add_global_preserve_scale"]))))
                    out_path.write_text(json.dumps(rows, indent=2))
                    print(f"[{held}/{arm}] fx {frac_x:.5f} cos_agg {e['cos_agg']:.6f} "
                          f"pds {e['pds']:.6f} q_w {e['a15_sgn_w']:.6f} "
                          f"({time.time()-t0:.0f}s)", flush=True)
                    if sname == "none":
                        r = ref.get((held, f"cd4_gn{cf:g}"))
                        assert r is not None, f"no reference arm cd4_gn{cf:g}"
                        dc = abs(r["cos_agg"] - e["cos_agg"])
                        dp = abs(r["pds"] - e["pds"])
                        assert dc < 1e-9 and dp < 1e-9, \
                            f"KA-A FAILED {held} cd4 f={cf}: dcos {dc:.2e} dpds {dp:.2e}"
                        print(f"    KA-A PASS: g=0 == resid218 cd4_gn{cf:g}", flush=True)
                        if cf == 0.10:
                            r2 = ka207.get((held, "f0.1"))
                            if r2 is not None:
                                d2 = abs(r2["cos_agg"] - e["cos_agg"])
                                assert d2 < 1e-9, f"KA-B FAILED {held}: {d2:.2e}"
                                print("    KA-B PASS: == 207 f0.1", flush=True)
        del op, whitened, seconds


if __name__ == "__main__":
    main()

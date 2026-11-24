"""216 Phase 1d — what would the pseudobulk object's extra 7,658 genes actually buy?

Phase 1 found the one materially different processed object in the published release:
`GWCD4i.pseudobulk_merged.h5ad` carries **17,771** official genes against the DE table's
**10,113**, and `energy216.py` measured that the extra genes hold **11.94%** of the gated
truth energy our CD4 axis cannot reach (87.13% -> 99.07%).

Turning that into a cosine is an elasticity question, and assuming a sqrt law would be
exactly the kind of unmeasured step this project keeps punishing.  So it is measured: the
CD4 residual's gene support is progressively CUT DOWN from below and the cosine is
re-measured at fixed `f`, then the local slope is extrapolated upward by the +11.94% the
pseudobulk would add.

The cut is from the LOW-EXPRESSION end, because that is what the extra genes are -- they are
precisely the genes the atlas's own DESeq2 expression filter dropped.  Cutting at random
would measure the wrong slope.

`f` is held fixed, so `solve_fraction` re-solves `fraction` against the shrunken residual and
the DELIVERED dose is constant across the rows.  That is the control that makes the rows
comparable.

    uv run python coverage216.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
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
from sweep216 import ALPHA, GAMMA, LAM, score_effects  # noqa: E402

DROP_ENERGY = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45]
DOSES = [0.10, 0.34]


def main() -> None:
    t0 = time.time()
    out_path = WORK / "coverage216.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        votes0 = np.asarray(cd4.gene_votes)
        scored = np.asarray(fold["scored"])

        # gated truth energy per gene, averaged over the three challenge frames
        e_gene = np.zeros(len(scored))
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            g = ((cpm > plolo.MIN_CPM) & (b > 0))[scored]
            e_gene += np.where(g, (truth ** 2).sum(0), 0.0)
        e_gene /= 3.0

        # control CPM on the official axis, from challenge frame A (the binding context)
        bA = plolo.challenge_bctrl("A")
        exprA = np.expm1(bA)
        cpmA = 1e6 * exprA / max(float(exprA.sum()), 1e-12)

        # the CD4 axis, expressed on the SCORED columns; rank it low-expression first
        cd4_axis = votes0.max(axis=0) > 0                       # genes CD4 ever votes on
        idx_scored = np.arange(len(scored))
        on_axis = cd4_axis[scored]
        order = idx_scored[on_axis][np.argsort(cpmA[scored][on_axis])]
        e_axis_total = float(e_gene[on_axis].sum())
        e_all_total = float(e_gene.sum())
        cum = np.cumsum(e_gene[order]) / max(e_axis_total, 1e-300)

        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
        del op
        gc.collect()

        for de in DROP_ENERGY:
            k = int(np.searchsorted(cum, de))
            drop_cols = scored[order[:k]]
            votes = votes0.copy()
            votes[:, drop_cols] = 0
            eff = np.where(votes > 0, np.asarray(cd4.mean_effect), 0.0)
            resid = E062.CD4Prepared(
                mode=cd4.mode + f"+cut{de:g}", mean_effect=eff, gene_votes=votes,
                target_condition_count=np.asarray(cd4.target_condition_count).copy(),
                condition_details=cd4.condition_details)
            kept_energy = float(e_gene[on_axis].sum() - e_gene[order[:k]].sum())
            for f_val in DOSES:
                arm = f"cut{de:g}_f{f_val:g}"
                if (held, arm) in done:
                    continue
                frac, per_unit = solve_fraction(whitened, resid, f_val)
                final, _ = E062.add_post_whiten_residual(
                    whitened, resid, frac, fold["wcols"])
                final[:, fold["panel_positions"]] = 0.0
                e = score_effects(final, fold, truth)
                rows.append(dict(
                    fold=held, arm=arm, drop_energy_target=de, f=f_val,
                    genes_dropped=int(k),
                    genes_on_cd4_axis=int(on_axis.sum()),
                    energy_share_of_axis_kept=float(kept_energy / max(e_axis_total, 1e-300)),
                    gated_energy_share_kept=float(kept_energy / max(e_all_total, 1e-300)),
                    fraction=float(frac), f_per_unit=float(per_unit),
                    seconds=round(time.time() - t0, 1),
                    **{kk: vv for kk, vv in e.items() if not kk.startswith("_")}))
                out_path.write_text(json.dumps(rows))
                print(f"[{held}/{arm}] dropped {k} genes, axis energy kept "
                      f"{kept_energy/e_axis_total:.4f} -> cos_agg {e['cos_agg']:.6f} "
                      f"pds {e['pds']:.6f}", flush=True)
                del final
            del votes, eff, resid
            gc.collect()
        del whitened, pre, cd4, fold, truth
        gc.collect()
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()

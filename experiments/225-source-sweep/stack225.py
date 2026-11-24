"""225 Phase 4 — the stack curve.  Does adding sources one at a time accumulate, or saturate?

This is the experiment's whole point.  215 found the leaders' direction arrived in ~5 small
separately-resolved steps.  If five public sources give ~5x DLD-1, the "they added datasets
one at a time" reading is right and the remaining 48 days go to acquisition.  If they
saturate at ~1.5x, that reading is dead and 224's learned-transfer channel wins.

Three families of arm, all marginal over `base222` (= 222's `cd4@0.34 + h1pair@0.28`):

  CUMULATIVE  the ranked survivors added one at a time, k = 1 .. n, each priced at its own
              dose sweep.  The curve's shape is the answer.
  LOO         the full stack minus each member — the mandate's falsifier 4 in its strict
              form: a source that duplicates a direction the stack already holds is worth 0.
  ALLSMALL    the union of EVERY census source with >=1 panel target, including the ones
              that are individually dead.  The "many small steps" hypothesis in its purest
              form: five sources that are each worth nothing may still be worth something
              together, because 159's 1/k law is about independent estimates, not rows.

Every arm carries the same F1 (target-shuffled) and F3 (row-matched CD4) controls as the
per-source sweep, and `dld1_lowmoi` alone rides every table as the known-answer row.

    uv run python stack225.py --stack dld1_lowmoi igvf_cardiomyocyte ...
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from harness225 import (CANCER, WORK, Fold, restrict_rows, shuffled,  # noqa: E402
                        source_residual, union_rows)

ALL_SMALL = ("dld1_lowmoi", "myllia_a549", "igvf_cardiomyocyte", "tian_ipsc_neuron",
             "sunshine_calu3", "xucao_hek293", "jiang_6line", "frangieh_melanoma",
             "norman_k562", "tian_ipsc_neuron_crispra", "wessels_thp1", "gasperini_k562",
             "tian2019_ipsc")
G = (0.10, 0.17, 0.28, 0.35, 0.45, 0.60)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stack", nargs="+", required=True,
                    help="ranked survivors, best first — added cumulatively")
    ap.add_argument("--allsmall", nargs="*", default=list(ALL_SMALL))
    ap.add_argument("--folds", nargs="*", default=list(CANCER))
    ap.add_argument("--doses", nargs="*", type=float, default=list(G))
    ap.add_argument("--out", default=str(WORK / "stack225.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        t0 = time.time()
        F = Fold(held)
        print(f"[{held}] fold {time.time() - t0:.0f}s", flush=True)

        def res_of(names):
            live = [n for n in names if n != held]
            if not live:
                return None, live
            return union_rows([source_residual(n, F.targets, "gene_norm")
                               for n in live]), live

        specs: dict[str, list[str]] = {}
        for k in range(1, len(args.stack) + 1):
            specs[f"CUM{k}"] = list(args.stack[:k])
        if len(args.stack) > 1:
            for drop in args.stack:
                specs[f"LOO-{drop}"] = [s for s in args.stack if s != drop]
        specs["ALLSMALL"] = list(args.allsmall)
        # the ALLSMALL union minus the one source we already know is live, so the question
        # "is there anything in the tail at all?" gets its own number
        specs["ALLSMALL-nodld1"] = [s for s in args.allsmall if s != "dld1_lowmoi"]

        for tag, names in specs.items():
            R, live = res_of(names)
            if R is None:
                continue
            variants = {"": R, "SHUF": shuffled(R), "CD4ROWS": restrict_rows(F.cd4, R)}
            for vtag, res in variants.items():
                for g in args.doses:
                    arm = f"{tag}{vtag}_{g:g}"
                    if (held, arm) in done:
                        continue
                    t1 = time.time()
                    a = F.arm(res, g, base="base222")
                    rows.append(dict(fold=held, arm=arm, kind="stack", tag=tag,
                                     variant=vtag, members=live, n_members=len(live), **a))
                    out_path.write_text(json.dumps(rows, indent=2))
                    print(f"[{held}/{arm}] n={len(live)} carry "
                          f"{a['audit']['carrying_rows']}/{a['audit']['n_rows']} "
                          f"relcos {a['rel_cos_pct']:+.3f}% dpds {a['d_pds']:+.5f} "
                          f"({time.time() - t1:.0f}s)", flush=True)
        del F
        import resid218
        resid218._PREP_CACHE.clear()

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n{len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()

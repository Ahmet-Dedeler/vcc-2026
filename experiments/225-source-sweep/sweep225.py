"""225 Phase 2/3 — price every fetchable census member on PLOLO, marginal over 222's arm.

FRAME (208's rider rule).  The baseline is `cd4@0.34 + h1pair@0.28` — `base222`, the arm 222
is shipping — so every number is the MARGINAL value on top of the CURRENT stack, not on top
of the stack the source was originally measured against.  Each source is also run on
`cd4@0.34` alone so the two framings can be separated.

SIZING is always the CARRYING-ROW rule.  The mandate's trap: production's median-over-every-
row rule reads `f_per_unit = 0` for any source below 50% coverage — which is every source in
this sweep — and a live source then scores exactly zero and looks dead.  Every arm asserts
`median_row_norm_carrying > 0` and records both medians.

FALSIFIERS, run on every source at every dose:
  F1  target-shuffled control (218's derangement; geometry, norms, support, carrying count
      all preserved exactly).  Two seeds for anything that survives.
  F2  non-carrying rows must not move (218's bar: |Δ| <= 0.018 pp).
  F3  row-matched CD4 control — the same rows, the same sizing, CD4's information instead.
  F4  leak audit — donor <-> held-surface overlap recorded on every arm.
  F5  non-vacuity assertion (above).
  F6  marginal-not-standalone: the `base222` framing IS the marginal test; §Phase 4's
      leave-one-out over the surviving stack completes it.

    uv run python sweep225.py [--folds ...] [--sources ...] [--doses ...]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from harness225 import (CANCER, WORK, Fold, carrying_mask, restrict_rows,  # noqa: E402
                        shuffled, solve_fraction, source_residual, sources218)

# every census member with a donor-shaped table on disk and >=1 panel target.
# `dld1_lowmoi` is the KNOWN-ANSWER row and rides in every screen (mandate rule).
SRC = ("dld1_lowmoi",                                   # 223's live source — known answer
       "igvf_cardiomyocyte",                            # 87/300 — the headline candidate
       "tian_ipsc_neuron", "sunshine_calu3", "xucao_hek293",
       "tian_ipsc_neuron_crispra", "tian_ipsc_neuron_crispra_flip",
       "norman_k562", "norman_k562_flip",
       "jiang_6line", "wessels_thp1", "gasperini_k562", "tian2019_ipsc",
       "frangieh_melanoma", "myllia_a549")              # 223's two smallest, for continuity
G = (0.10, 0.17, 0.28, 0.35, 0.45)


def leak_row(sname: str, held: str, targets) -> dict:
    p = Path(sources218.SOURCES[sname]["path"])
    root = Path(__file__).resolve().parents[2]
    src_targets = {x.strip() for x in (root / p / "perts.txt").read_text().split()
                   if x.strip()}
    return dict(source=sname, held_line=held,
                source_is_held_line=bool(sname == held),
                source_targets_on_fold=len(src_targets & set(targets)),
                fold_targets=len(targets),
                note="h1_2025 covers 85/85 of h1-test's targets; h1-test is NOT used "
                     "anywhere in 225 — every measurement is on 162's PLOLO frame")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(CANCER))
    ap.add_argument("--sources", nargs="*", default=list(SRC))
    ap.add_argument("--doses", nargs="*", type=float, default=list(G))
    ap.add_argument("--bases", nargs="*", default=["base222", "cd4"])
    ap.add_argument("--out", default=str(WORK / "sweep225.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        t0 = time.time()
        F = Fold(held)
        print(f"[{held}] fold {time.time() - t0:.0f}s  "
              f"base222 cos {F.scored['base222']['cos']:.6f} "
              f"pds {F.scored['base222']['pds']:.6f}", flush=True)
        for name in ("f0", "cd4", "base222"):
            if (held, name) not in done:
                e = F.scored[name]
                rows.append(dict(fold=held, arm=name, kind="baseline",
                                 **{k: v for k, v in e.items() if not k.startswith("_")}))

        for sname in args.sources:
            if sname == held:
                continue
            S = source_residual(sname, F.targets, "gene_norm")
            leak = leak_row(sname, held, F.targets)
            variants = {"": S, "SHUF": shuffled(S), "SHUF2": shuffled(S, seed=771103),
                        "CD4ROWS": restrict_rows(F.cd4, S)}
            for vtag, res in variants.items():
                for base in args.bases:
                    if vtag and base != "base222":
                        continue                      # falsifiers only on the headline frame
                    for g in args.doses:
                        arm = f"{sname}{vtag}_{base}_{g:g}"
                        if (held, arm) in done:
                            continue
                        t1 = time.time()
                        a = F.arm(res, g, base=base)
                        # the OLD production rule, recorded so the void stays visible
                        fr_old, _ = solve_fraction(F.whitened, res, g)
                        a["old_rule_fraction"] = float(fr_old)
                        rows.append(dict(fold=held, arm=arm, kind="arm", source=sname,
                                         variant=vtag, leak=leak, **a))
                        out_path.write_text(json.dumps(rows, indent=2))
                        print(f"[{held}/{arm}] frac {a['audit']['fraction']:.5f} "
                              f"carry {a['audit']['carrying_rows']}/{a['audit']['n_rows']} "
                              f"relcos {a['rel_cos_pct']:+.3f}% dpds {a['d_pds']:+.5f} "
                              f"OLD {fr_old:.3f} ({time.time() - t1:.0f}s)", flush=True)
        del F
        # the residual cache is keyed on (source, mode, targets) — nothing in it survives
        # a fold change, and 15 sources x 18,533 genes is ~2 GB if it is left to accumulate.
        import resid218
        resid218._PREP_CACHE.clear()

    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n{len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()

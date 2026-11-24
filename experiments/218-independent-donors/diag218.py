"""218 Phase 2 diagnostics — is the partial-coverage arm real, or is it an artefact?

`cos_agg = sqrt(1 - sse)` is ENERGY-weighted, so an arm that touches only 21% of the rows
could in principle move it by inflating those rows rather than by pointing them better.
Four checks, and the arm only survives if all four hold:

  1. MACRO cosine (unweighted mean of per-target cosines) must move in the same direction
     as `cos_agg`.  It cannot be inflated by row energy.
  2. `pds` -- a pure rank statistic, blind to magnitude -- must not fall.
  3. The row-matched CD4 control (`...CARRYCD4`) must be MUCH smaller.  Same rows, same
     sizing rule, same dose, CD4's information instead of the new source's: this absorbs
     the row-selection effect and any "more residual is better" effect entirely.
  4. The target-shuffled control must be negative at every dose.

Plus leave-one-fold-out dose selection (207's falsifier 3), because a dose chosen on the
same folds it is scored on is in-sample.

    uv run python diag218.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
from price218 import GATE_SHAPE, price  # noqa: E402

CD4_F = 0.34
G = (0.02, 0.05, 0.10, 0.17, 0.28, 0.45)


def main() -> None:
    rows = json.loads((WORK / "compose218.json").read_text())
    folds = sorted({r["fold"] for r in rows})
    by = {(r["fold"], r["arm"]): r for r in rows}
    base = {f: by[(f, f"cd4{CD4_F:g}+none0")] for f in folds}

    all_fams = sorted({r["second"] for r in rows if r["second"] != "none"})
    families = [f for f in all_fams if not f.endswith("_SHUF")]

    print("=" * 118)
    print(f"CHECK 1+2 — macro cosine and `pds` must agree with `cos_agg` "
          f"(second source on top of cd4 f={CD4_F:g})")
    print("=" * 118)
    print(f"{'arm':32s} {'nf':>3s} {'rel cos_agg%':>13s} {'rel MACRO cos%':>15s} "
          f"{'d pds':>10s} {'d q_w':>10s}")
    macro = {}
    for fam in all_fams:
        for g in G:
            arm = f"cd4{CD4_F:g}+{fam}{g:g}"
            sel = [(f, by[(f, arm)]) for f in folds if (f, arm) in by]
            if len(sel) < 2:
                continue
            ra = float(np.mean([r["cos_agg"] / base[f]["cos_agg"] - 1 for f, r in sel]))
            rm = float(np.mean([r["a23_cos"] / base[f]["a23_cos"] - 1 for f, r in sel]))
            dp = float(np.mean([r["pds"] - base[f]["pds"] for f, r in sel]))
            dq = float(np.mean([r["a15_sgn_w"] - base[f]["a15_sgn_w"] for f, r in sel]))
            nk = [(f, r) for f, r in sel if f != "kolf2"]
            macro[(fam, g)] = dict(
                n=len(sel), rel_agg=ra, rel_macro=rm, d_pds=dp, d_qw=dq,
                rel_agg_nokolf=float(np.mean(
                    [r["cos_agg"] / base[f]["cos_agg"] - 1 for f, r in nk])) if nk else None,
                rel_macro_nokolf=float(np.mean(
                    [r["a23_cos"] / base[f]["a23_cos"] - 1 for f, r in nk])) if nk else None,
                d_pds_nokolf=float(np.mean(
                    [r["pds"] - base[f]["pds"] for f, r in nk])) if nk else None,
                per_fold={f: dict(rel_agg=r["cos_agg"] / base[f]["cos_agg"] - 1,
                                  rel_macro=r["a23_cos"] / base[f]["a23_cos"] - 1,
                                  d_pds=r["pds"] - base[f]["pds"]) for f, r in sel})
            if fam in families:
                print(f"{arm:32s} {len(sel):3d} {100*ra:+13.3f} {100*rm:+15.3f} "
                      f"{dp:+10.6f} {dq:+10.6f}")
        if fam in families:
            print()

    print("=" * 118)
    print("CHECK 3 — the ROW-MATCHED CD4 control: same rows, same sizing, CD4's information")
    print("=" * 118)
    print(f"{'source':22s} {'g':>5s} {'real rel%':>10s} {'rowCD4 rel%':>12s} "
          f"{'shuffled rel%':>14s} {'real/rowCD4':>12s}")
    ctrl = {}
    for fam in families:
        if not fam.endswith("CARRY"):
            continue
        for g in G:
            r = macro.get((fam, g))
            c = macro.get((f"{fam}CD4", g))
            s = macro.get((f"{fam}_SHUF", g))
            if r is None:
                continue
            ratio = (r["rel_agg"] / c["rel_agg"]) if (c and c["rel_agg"] != 0) else None
            ctrl[(fam, g)] = dict(real=r, rowcd4=c, shuf=s, ratio=ratio)
            print(f"{fam:22s} {g:5.2f} {100*r['rel_agg']:+10.3f} "
                  f"{(100*c['rel_agg']) if c else float('nan'):+12.3f} "
                  f"{(100*s['rel_agg']) if s else float('nan'):+14.3f} "
                  f"{ratio if ratio else float('nan'):12.2f}")
        print()

    print("=" * 118)
    print("LEAVE-ONE-FOLD-OUT dose selection (207 falsifier 3) — selection and reporting "
          "frames never overlap")
    print("=" * 118)
    lofo = {}
    for fam in families:
        if fam.endswith("_SHUF"):
            continue
        held_prices = []
        chosen = []
        for held in folds:
            best, bestv = None, -1e9
            for g in G:
                arm = f"cd4{CD4_F:g}+{fam}{g:g}"
                sel = [(f, by[(f, arm)]) for f in folds
                       if f != held and (f, arm) in by]
                if len(sel) < 2:
                    continue
                rel = float(np.mean([r["cos_agg"] / base[f]["cos_agg"] - 1 for f, r in sel]))
                dp = float(np.mean([r["pds"] - base[f]["pds"] for f, r in sel]))
                v = price(rel, dp, k=0)["total_pds_from_plolo"]
                if v > bestv:
                    best, bestv = g, v
            if best is None:
                continue
            arm = f"cd4{CD4_F:g}+{fam}{best:g}"
            if (held, arm) not in by:
                continue
            r = by[(held, arm)]
            rel = r["cos_agg"] / base[held]["cos_agg"] - 1
            dp = r["pds"] - base[held]["pds"]
            held_prices.append(price(rel, dp, k=0)["total_pds_from_plolo"])
            chosen.append((held, best, rel, dp, held_prices[-1]))
        if held_prices:
            lofo[fam] = dict(mean=float(np.mean(held_prices)),
                             chosen=[(h, g, rel, dp, p) for h, g, rel, dp, p in chosen])
            print(f"{fam:22s} out-of-sample price {np.mean(held_prices):+.5f}  "
                  f"(gate {GATE_SHAPE:+.4f})  per fold: "
                  + ", ".join(f"{h[:12]}@g={g:g}:{p:+.5f}" for h, g, _, _, p in chosen))

    (WORK / "diag218.json").write_text(json.dumps(
        dict(macro={f"{k[0]}@{k[1]}": v for k, v in macro.items()},
             control={f"{k[0]}@{k[1]}": v for k, v in ctrl.items()},
             lofo=lofo), indent=2))
    print(f"\nwrote {WORK/'diag218.json'}")


if __name__ == "__main__":
    main()

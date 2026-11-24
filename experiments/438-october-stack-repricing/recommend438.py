"""438: the re-priced stack table for one panel and the recommended weights + dump recipe (used by dayof438.py; also run
on the rehearsal panels).

Decision rule (fixed 2026-09-23 on the rehearsal panels, before any October data; RESULTS.md section 5):
  * price on the LEVELLED route (the re-level is mandatory on D/E/F, 432), k = 1, at the board's q (the q A/B/C ships),
    paired on the same three cancer folds, key K3 (the only wired key).
  * the INCUMBENT is H1 x2, CD4 x1: the lowest-regret point over A/B/C and the three rehearsal draws (section 4.5; CD4
    x1.5, 417e's, is the highest-regret CD4 dose). It is 417e's H1 x2 half alone, CD4 at the DEF_ENV default.
  * candidates: H1 x{1, 1.5, 2, 2.5, 3} x CD4 x{0.75 .. 2}. H1 is capped at x3: A/B/C's own levelled optimum, and the
    rehearsal panels' gains beyond it are not robust (efp: leave-top-3-out negative).
  * a candidate MOVES only if (i) its paired mean advantage over the incumbent >= MARGIN, (ii) it is positive on 3/3
    folds, and (iii) its advantage survives 236's leave-top-3-out (the union of each fold's top-3 contributing targets
    removed) with a positive mean. Among movers, the largest leave-top-3-out advantage wins.
  * MARGIN = 0.0010 overall.
  * q: the board's q ships; the table reports the q 1.5 step at the 417a-calibrated transfer (k 0.65, k_pds 0.42) and
    flags a VETO if it is <= 0 on any fold of the day's panel.

    uv run python recommend438.py <tag> [fold_ops]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402

import numpy as np  # noqa: E402

from price438 import CANCER, REF, SLOTS, WORK, Panel  # noqa: E402
from robust438 import top3  # noqa: E402

MARGIN = 0.0010
CD4_G, H1_G = 0.3443655278649641, 0.28
H1S, CD4S = (0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5), (0.75, 1, 1.25, 1.5, 1.75, 2)
H1_CAND = (1, 1.5, 2, 2.5, 3)


def pf(P, arm, ref=REF, level=True, k=1.0):
    return {f: P.price(f, arm, k, ref, level)["ov"] for f in P.folds}


def line(P, lab, arm, ref=REF, boot=False):
    out = []
    for level in (False, True):
        p = pf(P, arm, ref, level)
        m = float(np.mean(list(p.values())))
        b = P.boot(arm, ref, level) if boot else None
        out.append(" / ".join(f"{p[f]:+.4f}" for f in P.folds) + f" | **{m:+.5f}** {sum(v > 0 for v in p.values())}/3"
                   + (f" [{b['ci'][0]:+.4f}, {b['ci'][1]:+.4f}]" if b else ""))
    rho = "/".join(f"{P.E[arm]['rho'][c]:.3f}" for c in SLOTS)
    r = P.price(P.folds[0], arm, 1.0, ref, True)
    dd = {k: float(np.mean([P.price(f, arm, 1.0, ref, True)[k] for f in P.folds])) for k in ("d_rel", "dpds_bd", "dmse")}
    return (f"| {lab} | {dd['d_rel']:+.4f} | {dd['dpds_bd']:+.5f} | {dd['dmse']:+.4f} | {out[0]} | {out[1]} | {rho} |")


INCUMBENT = (2, 1)


def main(tag, fold_ops="held", q_board="1.25", **kw):
    q_inc = q_board
    P = Panel(tag, fold_ops, **kw)
    full = P.full
    ctxs = full.get("contexts", list(SLOTS))
    print(f"\n### {tag} ({fold_ops} fold operators), contexts {ctxs}, folds {P.folds}")
    cov = full.get("coverage")
    if cov:
        print(f"coverage: donors {cov['donors']}, H1 union rows {cov['h1_union_rows']} {cov['h1_sources']}, "
              f"CD4 complete {cov['cd4_complete_rows']}, partial {cov['cd4_partial_rows']}")
    print("retention ret_f " + ", ".join(f"{f} {P.ret[f]:.3f}" for f in P.folds) + "; claims "
          + "/".join(f"{P.claims[c]:.5f}" for c in SLOTS))
    for n in (REF, f"h2c1.5|K3+sh{q_inc}"):
        e = P.E[n]
        print(f"{n}: rho unlevelled " + "/".join(f"{e['rho'][c]:.4f}" for c in SLOTS)
              + ("; levelled " + "/".join(f"{e['rho_level'][c]:.4f}" for c in SLOTS) if "rho_level" in e else ""))
    hdr = ("| design (vs the 417a analogue) | Δrel cos (lev) | Δpds_bd (lev) | Δmse k=1 (lev) | unlevelled k=1 per fold "
           "\\| mean n+ [CI] | LEVELLED k=1 per fold \\| mean n+ [CI] | ρ unlevelled |\n|---|---:|---:|---:|---|---|---|")
    print(hdr)
    rows = [("H1 x2 half", "h2c1|K3+sh1.25", True), ("CD4 x1.5 half", "h1c1.5|K3+sh1.25", True),
            ("stack (H1 x2 + CD4 x1.5), q 1.25", "h2c1.5|K3+sh1.25", True),
            ("stack, q 1.5", "h2c1.5|K3+sh1.5", True), ("PP, q 1.5", "h1c1|K3+sh1.5", False),
            ("stack, LODO key (0.5, 3), q 1.25", "h2c1.5|V1_5_3+sh1.25", True),
            ("stack, LODO key (0.5, 3), q 1.5", "h2c1.5|V1_5_3+sh1.5", True)]
    rows += [(f"H1 x{m:g}, CD4 x1.5", f"h{m:g}c1.5|K3+sh1.25", False) for m in (0, 1, 1.5, 2.5, 3)]
    rows += [(f"H1 x2, CD4 x{m:g}", f"h2c{m:g}|K3+sh1.25", False) for m in (0.75, 1, 1.25, 1.75, 2)]
    for lab, n, b in rows:
        print(line(P, lab, n, boot=b))
    # ---------------------------------------------------------------- the surface and the decision
    rec = {}
    inc_k = INCUMBENT
    for q in ("1.25", "1.5"):
        inc = f"h{inc_k[0]:g}c{inc_k[1]:g}|K3+sh{q}"
        grid = {(mh, mc): pf(P, f"h{mh:g}c{mc:g}|K3+sh{q}") for mh in H1S for mc in CD4S}
        mean = {k: float(np.mean(list(v.values()))) for k, v in grid.items()}
        best = max(mean, key=mean.get)
        print(f"\nLEVELLED k=1 surface, K3 q {q} (cancer mean vs 417a analogue; rows H1 x, cols CD4 x):")
        print("| H1 \\ CD4 | " + " | ".join(f"x{c:g}" for c in CD4S) + " |\n|---|" + "---:|" * len(CD4S))
        for mh in H1S:
            print(f"| x{mh:g} | " + " | ".join(("**" if (mh, mc) == best else "") + f"{mean[(mh, mc)]:+.5f}"
                                               + ("**" if (mh, mc) == best else "") for mc in CD4S) + " |")
        movers = []
        for k in [(mh, mc) for mh in H1_CAND for mc in CD4S if (mh, mc) != inc_k]:
            n = f"h{k[0]:g}c{k[1]:g}|K3+sh{q}"
            d = pf(P, n, inc)
            dm = float(np.mean(list(d.values())))
            if dm < MARGIN or not all(v > 0 for v in d.values()):
                continue
            Pk = Panel(tag, fold_ops, exclude=top3(P, n, inc), **kw)
            lt3 = float(np.mean([Pk.price(f, n, 1.0, inc, True)["ov"] for f in Pk.folds]))
            movers.append((k, dm, lt3))
        movers = [m for m in movers if m[2] > 0]
        choice = max(movers, key=lambda m: m[2])[0] if movers else inc_k
        rec[q] = dict(best=best, best_mean=mean[best], abc_weights_mean=mean[(2, 1.5)], incumbent_mean=mean[inc_k],
                      choice=choice, movers=movers, choice_mean=mean[choice],
                      abc_weights_cost=mean[choice] - mean[(2, 1.5)])
        print(f"argmax H1 x{best[0]:g} CD4 x{best[1]:g} {mean[best]:+.5f}; incumbent (H1 x{inc_k[0]:g}, CD4 x{inc_k[1]:g}) "
              f"{mean[inc_k]:+.5f}; A/B/C weights (x2, x1.5) {mean[(2, 1.5)]:+.5f}")
        print("movers vs the incumbent (paired mean >= %.4f, 3/3, leave-top-3-out > 0): " % MARGIN
              + (", ".join(f"H1 x{m[0][0]:g} CD4 x{m[0][1]:g} {m[1]:+.5f} (lt3 {m[2]:+.5f})" for m in movers) or "none"))
        print(f"RULE CHOICE at q {q}{' (THE BOARD Q)' if q == q_board else ''}: H1 x{choice[0]:g}, CD4 x{choice[1]:g} "
              f"(vs A/B/C weights {mean[choice] - mean[(2, 1.5)]:+.5f})")
    c = rec[q_board]["choice"]
    n125, n15 = f"h{c[0]:g}c{c[1]:g}|K3+sh1.25", f"h{c[0]:g}c{c[1]:g}|K3+sh1.5"
    dq = {f: P.price(f, n15, 0.65, n125, True, k_pds=0.42)["ov"] for f in P.folds}
    dq1 = pf(P, n15, n125)
    veto = not all(v > 0 for v in dq.values())
    print(f"\nq 1.5 vs 1.25 on H1 x{c[0]:g} CD4 x{c[1]:g}, levelled: k=1 " + " / ".join(f"{dq1[f]:+.4f}" for f in P.folds)
          + f" (mean {np.mean(list(dq1.values())):+.5f}); 417a-calibrated " + " / ".join(f"{dq[f]:+.4f}" for f in P.folds)
          + f" (mean {np.mean(list(dq.values())):+.5f})" + (" -> VETO q 1.5 on this panel" if veto else " -> no veto"))
    mh, mc = c
    q_choice = f"{q_board} (the board's)" if not (veto and q_board == "1.5") else "1.25 (q 1.5 VETOED)"
    print(f"\nRECOMMENDATION ({tag}): H1 x{mh:g}, CD4 x{mc:g}, q {q_choice}, key K3 (410 re-fitted on the day's operator)")
    print(f"  dump recipe (432 section 8 step 6): VCC_H1_UNION_G={H1_G * mh:.10g}"
          + ("" if mc == 1 else f" VCC_CD4_TARGET_G={CD4_G * mc:.16g}") + " VCC_BULK_DUMP=$DUMP"
          + (" (CD4 at the DEF_ENV default: no CD4 override)" if mc == 1 else ""))
    out = dict(tag=tag, fold_ops=fold_ops, rec=rec, q_step_calibrated=dq, q_step_k1=dq1, q_veto=veto,
               choice=dict(h1=mh, cd4=mc, q=q_choice), env=dict(VCC_H1_UNION_G=H1_G * mh, VCC_CD4_TARGET_G=CD4_G * mc))
    (WORK / tag / f"recommend438-{fold_ops}.json").write_text(json.dumps(out, indent=1, default=lambda o: list(o) if isinstance(o, tuple) else float(o)))
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("fold_ops", nargs="?", default="held")
    ap.add_argument("--full", default=None)
    ap.add_argument("--fold-energy", action="store_true")
    ap.add_argument("--board-q", default="1.25", choices=("1.25", "1.5"))
    a = ap.parse_args()
    main(a.tag, a.fold_ops, a.board_q, full_tag=a.full, energy="fold" if a.fold_energy else "full")

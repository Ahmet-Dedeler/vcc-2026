"""216 — INBOX-3's four deliverables from the per-context sweep.

Known-answer row: the pooled column recomputed here must reproduce `sweep216.json` exactly,
since `ctx216.py` re-runs the identical grid and only changes what it stores.

    uv run python ctxanalyse216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CTX = ("A", "B", "C")
LINEAGE = {"A": "JURKAT — T-lymphoblastic (the donor's own lineage)",
           "B": "HELA — cervix", "C": "CAL33 — head-and-neck squamous"}
SIGMA_PDS = 0.001494
K_PDS_BOARD = 1.9869          # board216.py, measured on the cd4f034 board point


def turnover(fs, ys):
    i = int(np.argmax(ys))
    if i in (0, len(ys) - 1):
        return float(fs[i]), "grid endpoint"
    x = np.log(np.asarray(fs[i - 1:i + 2], float))
    y = np.asarray(ys[i - 1:i + 2], float)
    a, b, _ = np.polyfit(x, y, 2)
    if a >= 0:
        return float(fs[i]), "grid"
    return float(np.exp(-b / (2 * a))), "quadratic in log f"


def main() -> None:
    rows = json.loads((WORK / "ctx216.json").read_text())
    sw = {(r["fold"], r["arm"]): r for r in
          json.loads((WORK / "sweep216.json").read_text())}
    by = {(r["fold"], r["f"]): r for r in rows}
    doses = sorted({r["f"] for r in rows})
    have = [f for f in doses if all((h, f) in by for h in FOLDS)]

    # known-answer: the pooled column must reproduce sweep216
    worst = 0.0
    for f in have:
        for h in FOLDS:
            worst = max(worst, abs(by[(h, f)]["a23_cos_agg"]
                                   - sw[(h, "f%g" % f)]["cos_agg"]))
    print(f"[{'PASS' if worst < 1e-12 else 'FAIL'}] KA pooled cos_agg reproduces "
          f"sweep216.json: max|diff| = {worst:.3e}")

    def mc(f, ctx, key):
        return float(np.mean([by[(h, f)][f"a23_{ctx}_{key}"] for h in FOLDS]))

    print(f"\n{'f':>7}" + "".join(f"{'cos '+c:>12}" for c in CTX)
          + f"{'pooled':>12}" + "".join(f"{'pds '+c:>11}" for c in CTX))
    print("-" * 90)
    curves = {c: [] for c in CTX}
    pcurves = {c: [] for c in CTX}
    pooled = []
    for f in have:
        cs = [mc(f, c, "cos_agg") for c in CTX]
        ps = [mc(f, c, "pds") for c in CTX]
        for c, v, p in zip(CTX, cs, ps):
            curves[c].append(v)
            pcurves[c].append(p)
        pooled.append(float(np.mean([by[(h, f)]["a23_cos_agg"] for h in FOLDS])))
        print(f"{f:>7.3g}" + "".join(f"{v:>12.6f}" for v in cs) + f"{pooled[-1]:>12.6f}"
              + "".join(f"{p:>11.6f}" for p in ps))

    fs = [f for f in have if f > 0]
    idx = [have.index(f) for f in fs]
    out = {"known_answer_max_diff": worst, "per_context": {}}
    print(f"\n{'context':<44}{'cos argmax f':>14}{'value':>11}{'method':>22}")
    argmax = {}
    for c in CTX:
        ys = [curves[c][i] for i in idx]
        fstar, how = turnover(fs, ys)
        grid = fs[int(np.argmax(ys))]
        argmax[c] = fstar
        out["per_context"][c] = dict(lineage=LINEAGE[c], f_star=fstar,
                                     grid_argmax=grid, value=max(ys),
                                     method=how, curve=dict(zip(map(str, fs), ys)),
                                     pds_curve=dict(zip(map(str, fs),
                                                        [pcurves[c][i] for i in idx])))
        print(f"{c + ' (' + LINEAGE[c] + ')':<44}{fstar:>14.3f}{max(ys):>11.6f}{how:>22}")
    ys_p = [pooled[i] for i in idx]
    f_pool, how_p = turnover(fs, ys_p)
    out["pooled"] = dict(f_star=f_pool, method=how_p)
    print(f"{'POOLED (what we would ship)':<44}{f_pool:>14.3f}{max(ys_p):>11.6f}{how_p:>22}")

    spread = max(argmax.values()) / min(argmax.values())
    out["spread"] = dict(min_f=min(argmax.values()), max_f=max(argmax.values()),
                         ratio=spread, argmax_by_context=argmax,
                         A_is_highest=bool(argmax["A"] == max(argmax.values())),
                         A_is_lowest=bool(argmax["A"] == min(argmax.values())))
    print(f"\nSPREAD between the three maximisers: "
          f"{min(argmax.values()):.3f} to {max(argmax.values()):.3f} = "
          f"x{spread:.3f}.  A (the T-lineage context) is "
          f"{'the HIGHEST' if out['spread']['A_is_highest'] else ('the LOWEST' if out['spread']['A_is_lowest'] else 'in the middle')}.")

    # INBOX-3 item 4: the loss from shipping one context's optimum into another
    print(f"\nloss from shipping the POOLED optimum (f = {f_pool:.3f}) into each context, "
          f"and from shipping each other context's optimum:")
    print(f"{'context':<10}{'own opt':>10}{'at pooled':>12}{'rel loss':>11}"
          f"{'worst cross-context loss':>26}")
    losses = {}
    for c in CTX:
        ys = np.asarray([curves[c][i] for i in idx])
        own = float(ys.max())

        def at(fv):
            j = int(np.argmin([abs(np.log(x / fv)) for x in fs]))
            return float(ys[j])
        at_pool = at(f_pool)
        cross = [at(argmax[o]) for o in CTX if o != c]
        losses[c] = dict(own=own, at_pooled=at_pool,
                         rel_loss_pooled=at_pool / own - 1.0,
                         worst_cross=min(cross) / own - 1.0)
        print(f"{c:<10}{own:>10.6f}{at_pool:>12.6f}{100*(at_pool/own-1):>10.3f}%"
              f"{100*losses[c]['worst_cross']:>25.3f}%")
    out["transfer_loss"] = losses
    out["max_rel_loss_at_pooled"] = min(v["rel_loss_pooled"] for v in losses.values())

    # per-context pds cost at the candidate doses, in board-projected raw units
    print(f"\nper-context `pds` cost vs f = 0.34, PLOLO raw and board-projected "
          f"(x{K_PDS_BOARD:.3f}):")
    print(f"{'f':>7}" + "".join(f"{'PLOLO '+c:>12}" for c in CTX)
          + "".join(f"{'board '+c:>12}" for c in CTX))
    pds_tab = {}
    for f in (0.20, 0.28, 0.40, 0.55, 0.80, 1.00):
        if f not in have:
            continue
        d = [mc(f, c, "pds") - mc(0.34, c, "pds") for c in CTX]
        pds_tab[str(f)] = dict(plolo=dict(zip(CTX, d)),
                               board=dict(zip(CTX, [K_PDS_BOARD * x for x in d])))
        print(f"{f:>7.2f}" + "".join(f"{x:>+12.5f}" for x in d)
              + "".join(f"{K_PDS_BOARD*x:>+12.5f}" for x in d))
    out["pds_cost_per_context"] = pds_tab
    (WORK / "ctxanalyse216.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

"""214 — read `work/pds214.json` into the known-answer gates and the substitution table.

Predictor: 162's pre-committed one, reproduced from 200's `report_pds.py` -- the unweighted
mean over all four donor folds of the challenge-frame mean (A/B/C).  Scaling: raw -> scaled
by the board `pds` span 0.45208, scaled -> overall by rule 13's k = 1.011 and /6.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
SPAN = 0.95228 - 0.50020
K13 = 1.011
SIGMA = 0.001494        # 162's leave-one-arm-out RMS on scaled `pds`

R: dict = {}


def cell(fold, arm, rule, field="pds", frames=("chalA", "chalB", "chalC")):
    v = [R[(fold, arm, f, rule)][field] for f in frames if (fold, arm, f, rule) in R]
    return float(np.mean(v)) if len(v) == len(frames) else None


def delta(arm, base, rule="inf", field="pds", frames=("chalA", "chalB", "chalC")):
    per = {}
    for f in FOLDS:
        a, b = cell(f, arm, rule, field, frames), cell(f, base, rule, field, frames)
        if a is not None and b is not None:
            per[f] = a - b
    if not per:
        return None
    v = np.asarray(list(per.values()))
    return dict(per=per, mean=float(v.mean()), sd=float(v.std()), n=int(v.size),
                scaled=float(v.mean()) / SPAN,
                overall=float(v.mean()) / SPAN * K13 / 6,
                signs="".join("+" if per[f] > 0 else ("-" if per[f] < 0 else "0")
                              for f in FOLDS if f in per))


def main() -> None:
    rows = json.loads(open(os.path.join(WORK, "pds214.json")).read())
    for r in rows:
        R[(r["fold"], r["arm"], r["frame"], r["rule"])] = r
    fails, gates = [], []

    def gate(name, got, want, tol, unit=""):
        ok = abs(got - want) <= tol
        gates.append(dict(gate=name, got=got, want=want, tol=tol, pass_=bool(ok)))
        print(f"  {name:<52} {got:>12.6f}   expect {want:>11.6f}   "
              f"|d| {abs(got-want):.2e}   {'PASS' if ok else 'FAIL'}{unit}")
        if not ok:
            fails.append(name)

    print("=" * 118)
    print("PHASE 0 — KNOWN-ANSWER ROWS (rule 10).  Every one aborts.")
    print("=" * 118)

    print("\nKA1 — 200's `SHIP` (= a23) absolute chal-mean `pds`, per fold")
    KA1 = dict(orion_hct116=0.791823, orion_hek293t=0.777375,
               k562_gwps=0.789121, kolf2=0.610128)
    for f in FOLDS:
        gate(f"KA1 SHIP {f}", cell(f, "SHIP", "inf"), KA1[f], 5e-6)
    m = float(np.mean([cell(f, "SHIP", "inf") for f in FOLDS]))
    gate("KA1 SHIP 4-fold mean", m, 0.742112, 5e-6)

    print("\nKA1c — the mandate's own a15 / a23 pair, which is the `k562_gwps` fold")
    gate("KA1c a15 (mass fix OFF) k562_gwps", cell("k562_gwps", "noMassFix", "inf"),
         0.795784, 5e-6)
    gate("KA1c a23 (mass fix ON)  k562_gwps", cell("k562_gwps", "SHIP", "inf"),
         0.789121, 5e-6)

    print("\nKA1b — state-carry control: `noMassFix` must equal `g_pm050_f010` exactly "
          "(field-identical configs)")
    w = max(abs(cell(f, "noMassFix", "inf") - cell(f, "g_pm050_f010", "inf")) for f in FOLDS)
    gate("KA1b worst |d| over folds", w, 0.0, 0.0)

    print("\nKA2/3/4 — 200's leave-one-out rows against `SHIP`, on scaled `pds`")
    LOO = [("noMassFix", "KA2 mass fix (a15 - a23)", -0.004014, 3e-5),
           ("noPM_mf", "KA3 PerturbMatch (noPM - SHIP)", -0.011187, 3e-4),
           ("noWhiten_mf", "KA4a union whitening", -0.151892, 2e-3),
           ("gate20_mf", "KA4b CPM gate 10 vs 20", -0.004060, 3e-4),
           ("noCD4_mf", "KA4c CD4 residual f=0.10", -0.015496, 3e-4)]
    for a, label, want, tol in LOO:
        d = delta(a, "SHIP")
        gate(label, d["scaled"], want, tol)
        print(f"       raw {d['mean']:+.6f}  fold sd {d['sd']:.6f}  overall "
              f"{d['overall']:+.6f}  signs {d['signs']}")

    print("\nKA5 — 196/208's rounding rows, chalA, a23 renderer, f = 0.10, pm 0.50")
    for rule, want in (("floor", 0.739016), ("ham", 0.742224), ("inf", 0.742137)):
        v = float(np.mean([cell(f, "mf_pm050_f010", rule, frames=("chalA",))
                           for f in FOLDS]))
        gate(f"KA5 {rule} 4-fold mean chalA", v, want, 5e-6)

    print("\nKA6 — 169's cached PLOLO chalA rows at m = inf")
    KA6 = dict(k562_gwps=0.790712, orion_hct116=0.790896,
               orion_hek293t=0.771621, kolf2=0.615317)
    for f in FOLDS:
        gate(f"KA6 {f} chalA inf", cell(f, "mf_pm050_f010", "inf", frames=("chalA",)),
             KA6[f], 5e-6)

    if fails:
        json.dump(gates, open(os.path.join(WORK, "gates214.json"), "w"), indent=1)
        raise SystemExit(f"ABORT (rule 10): {fails}")
    print("\n  ALL KNOWN-ANSWER ROWS PASS")

    # ---------------------------------------------------------------- Phase 1: the grid
    DOSES = [(None, "0.00"), (0.05, "0.05"), (0.10, "0.10"), (0.20, "0.20"),
             (0.34, "0.34"), (0.50, "0.50"), (0.80, "0.80")]

    def tag(pref, pm, f):
        return f"{pref}_pm{int(round(pm*100)):03d}_f{0 if f is None else int(round(f*100)):03d}"

    print("\n" + "=" * 118)
    print("PHASE 1 — IS PERTURBMATCH A SUBSTITUTE FOR THE CD4 DOSE?")
    print("  value of PerturbMatch = pds(pm 0.50) - pds(pm 0.00), at each CD4 dose")
    print("  incumbent renderer (mass_fix None), gamma 0.65, alpha 0.982, 4 folds x A/B/C")
    print("=" * 118)
    print(f"  {'CD4 f':>7}{'pds(pm.5)':>12}{'pds(pm0)':>11}{'d raw':>10}{'fold sd':>9}"
          f"{'scaled':>10}{'overall':>10}{'/sigma':>8}  signs"
          f"{'cos(pm.5)':>12}{'cos(pm0)':>11}{'d cos %':>9}")
    grid = {}
    for f, lab in DOSES:
        a1, a0 = tag("g", 0.50, f), tag("g", 0.0, f)
        d = delta(a1, a0)
        dc = delta(a1, a0, field="cos")
        p1 = float(np.mean([cell(x, a1, "inf") for x in FOLDS]))
        p0 = float(np.mean([cell(x, a0, "inf") for x in FOLDS]))
        c1 = float(np.mean([cell(x, a1, "inf", "cos") for x in FOLDS]))
        c0 = float(np.mean([cell(x, a0, "inf", "cos") for x in FOLDS]))
        grid[lab] = dict(pds_pm50=p1, pds_pm00=p0, d=d, dcos=dc, cos_pm50=c1, cos_pm00=c0)
        print(f"  {lab:>7}{p1:>12.6f}{p0:>11.6f}{d['mean']:>+10.6f}{d['sd']:>9.6f}"
              f"{d['scaled']:>+10.6f}{d['overall']:>+10.6f}{d['scaled']/SIGMA:>8.1f}  "
              f"{d['signs']:<6}{c1:>12.6f}{c0:>11.6f}"
              f"{100*(c1/c0-1):>+9.3f}")

    print("\n  the mirror image — value of the CD4 DOSE at each PerturbMatch setting")
    print(f"  {'CD4 f':>7}{'d pds | pm=0.50':>18}{'d pds | pm=0.00':>18}"
          f"{'interaction':>14}{'overall interact':>18}")
    for f, lab in DOSES[1:]:
        d5 = delta(tag("g", 0.50, f), tag("g", 0.50, None))
        d0 = delta(tag("g", 0.0, f), tag("g", 0.0, None))
        inter = d5["mean"] - d0["mean"]
        print(f"  {lab:>7}{d5['mean']:>+18.6f}{d0['mean']:>+18.6f}{inter:>+14.6f}"
              f"{inter/SPAN*K13/6:>+18.6f}")

    print("\n  the blend fraction's shape (pm 0 -> 1), at both doses, raw `pds`")
    print(f"  {'pm':>6}" + "".join(f"{lab:>14}" for _, lab in DOSES if lab in ("0.10", "0.34")))
    for pm in (0.0, 0.25, 0.50, 0.75, 1.00):
        vals = []
        for f, lab in DOSES:
            if lab not in ("0.10", "0.34"):
                continue
            a = tag("g", pm, f)
            vals.append(float(np.mean([cell(x, a, "inf") for x in FOLDS]))
                        if (FOLDS[0], a, "chalA", "inf") in R else float("nan"))
        print(f"  {pm:>6.2f}" + "".join(f"{v:>14.6f}" for v in vals))

    print("\n  the OPERATOR: PerturbMatch's value with the whitener ON (gamma 0.65) and "
          "OFF (gamma 0)")
    print(f"  {'CD4 f':>7}{'gamma 0.65 d':>15}{'gamma 0 d':>13}{'ratio':>9}"
          f"{'pds(g0,pm.5)':>15}{'pds(g0,pm0)':>14}")
    for f, lab in DOSES:
        if lab not in ("0.00", "0.10", "0.34"):
            continue
        dw = delta(tag("g", 0.50, f), tag("g", 0.0, f))
        dn = delta(tag("nw", 0.50, f), tag("nw", 0.0, f))
        p1 = float(np.mean([cell(x, tag("nw", 0.50, f), "inf") for x in FOLDS]))
        p0 = float(np.mean([cell(x, tag("nw", 0.0, f), "inf") for x in FOLDS]))
        print(f"  {lab:>7}{dw['mean']:>+15.6f}{dn['mean']:>+13.6f}"
              f"{(dn['mean']/dw['mean'] if dw['mean'] else float('nan')):>9.2f}"
              f"{p1:>15.6f}{p0:>14.6f}")

    print("\n  the RENDERER: effect level (058/064's instrument) vs delivered, and the "
          "mass fix")
    print(f"  {'CD4 f':>7}{'d intended':>13}{'d delivered':>14}{'d deliv (a23 mf)':>19}")
    for f, lab in DOSES:
        if lab not in ("0.00", "0.10", "0.34"):
            continue
        di = delta(tag("g", 0.50, f), tag("g", 0.0, f), rule="intended")
        dd = delta(tag("g", 0.50, f), tag("g", 0.0, f))
        dm = (delta(tag("mf", 0.50, f), tag("mf", 0.0, f))
              if (FOLDS[0], tag("mf", 0.50, f), "chalA", "inf") in R else None)
        print(f"  {lab:>7}{di['mean']:>+13.6f}{dd['mean']:>+14.6f}"
              + (f"{dm['mean']:>+19.6f}" if dm else f"{'-':>19}"))

    # ---------------------------------------------------------------- Phase 2: composition
    print("\n" + "=" * 118)
    print("PHASE 2 — `pds` FOR ARMS A, B AND A∘B, MARGINAL TO `cd4f034` "
          "(f = 0.34, pm 0.50, mass_fix None)")
    print("=" * 118)
    INC = "g_pm050_f034"
    comp = {}
    print(f"  {'arm':<28}{'rule':>7}{'pds':>11}{'d raw':>10}{'fold sd':>9}{'scaled':>10}"
          f"{'overall':>10}{'/sigma':>8}  signs      cos      d cos %")
    for label, a, rule in (("incumbent cd4f034", INC, "floor"),
                           ("A  drop PerturbMatch", "g_pm000_f034", "floor"),
                           ("B  Hamilton rounding", INC, "ham"),
                           ("A∘B both", "g_pm000_f034", "ham"),
                           ("(reference) no rounding", INC, "inf")):
        d = delta(a, INC, rule=rule) if not (a == INC and rule == "floor") else None
        base = [cell(f, INC, "floor") for f in FOLDS]
        # deltas are always taken against the incumbent's OWN floor row
        per = {f: cell(f, a, rule) - cell(f, INC, "floor") for f in FOLDS}
        v = np.asarray(list(per.values()))
        pds = float(np.mean([cell(f, a, rule) for f in FOLDS]))
        cs = float(np.mean([cell(f, a, rule, "cos") for f in FOLDS]))
        cs0 = float(np.mean([cell(f, INC, "floor", "cos") for f in FOLDS]))
        signs = "".join("+" if per[f] > 0 else ("-" if per[f] < 0 else "0") for f in FOLDS)
        sc = float(v.mean()) / SPAN
        comp[label] = dict(pds=pds, raw=float(v.mean()), sd=float(v.std()), scaled=sc,
                           overall=sc * K13 / 6, per=per, cos=cs, dcos_pct=100*(cs/cs0-1),
                           signs=signs)
        print(f"  {label:<28}{rule:>7}{pds:>11.6f}{v.mean():>+10.6f}{v.std():>9.6f}"
              f"{sc:>+10.6f}{sc*K13/6:>+10.6f}{sc/SIGMA:>8.1f}  {signs:<6}"
              f"{cs:>11.6f}{100*(cs/cs0-1):>+11.4f}")
    a = comp["A  drop PerturbMatch"]["overall"]
    b = comp["B  Hamilton rounding"]["overall"]
    ab = comp["A∘B both"]["overall"]
    print(f"\n  COMPOSITION, measured not inferred (178's rule):  A + B = {a+b:+.6f}   "
          f"measured A∘B = {ab:+.6f}   difference {ab-(a+b):+.6f}")

    print("\n  per-fold `pds` for the four composition arms (leak-free fold `kolf2` "
          "included, falsifier 3)")
    print(f"  {'arm':<28}" + "".join(f"{f:>16}" for f in FOLDS))
    for label, a, rule in (("incumbent cd4f034", INC, "floor"),
                           ("A  drop PerturbMatch", "g_pm000_f034", "floor"),
                           ("B  Hamilton rounding", INC, "ham"),
                           ("A∘B both", "g_pm000_f034", "ham")):
        print(f"  {label:<28}" + "".join(f"{cell(f, a, rule):>16.6f}" for f in FOLDS))

    json.dump(dict(gates=gates, grid=grid, comp=comp, span=SPAN, k13=K13, sigma=SIGMA),
              open(os.path.join(WORK, "report214.json"), "w"), indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'report214.json')}")


if __name__ == "__main__":
    sys.exit(main())

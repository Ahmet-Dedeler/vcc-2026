"""231 sec3 — the h1-test deflator, RE-MEASURED on the currently shipped stack.

214 sec1.5 measured that h1-test overstates any PerturbMatch reading by **x1.533**, because
PerturbMatch enters the K562 slot only and h1-test carries TWO production donors where the
shipped stack carries THREE.  That deflator was measured against `cd4f034`.

The shipped stack is now `h1union`: three donors PLUS the `gse295214_h1`/`h1_2025` union
residual at `VCC_H1_UNION_G = 0.28`.  PerturbMatch's share of it is smaller again, so 208's
rider rule (*"re-price every rider against the current incumbent"*) requires the deflator to
be re-measured.  This does that with 214's own instrument.

Three surfaces, one statistic — the relative displacement of dropping PerturbMatch:
    ||E(pm=0.50) - E(pm=0.00)|| / ||E(pm=0.50)||
plus the per-row median and the per-row cosine.

KNOWN-ANSWER ROWS (rule 10; the run aborts on any failure), all from 214 sec1.5:
  KA1  challenge / cd4f034 (3 donors, 300 targets, 251 PM-covered):
       rel_frobenius 0.27059, row median 0.24704, row-cos median 0.97081
  KA2  h1-test (2 donors, 85 targets, 79 PM-covered):
       rel_frobenius 0.41486, row median 0.36445, row-cos median 0.93314
  KA3  the deflator itself: 0.41486 / 0.27059 = 1.533
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for p in ("experiments/162-retrodicting-validator", "experiments/200-production-ablation",
          "experiments/218-independent-donors", "experiments/225-source-sweep",
          "experiments/196-emitted-depth", "src"):
    sys.path.insert(0, str(ROOT / p))
sys.path.insert(0, str(ROOT / "experiments/214-cheap-positives"))

import harness225 as H                                                   # noqa: E402
import plolo                                                             # noqa: E402
import chalparts200                                                      # noqa: E402
import h1parts200                                                        # noqa: E402
import mse200                                                            # noqa: E402
from folds200 import consensus                                           # noqa: E402
from plolo import E062, PROD_OPS, Operator, solve_fraction               # noqa: E402
from resid218 import ALPHA, GAMMA, LAM, source_residual                  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows       # noqa: E402

plolo.WORK = WORK / "plolo"

CD4_F, UNION_G, UNION = H.CD4_F, H.UNION_G, H.UNION
FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<44} got {got:>11.5f} exp {exp:>10.5f}  {note}")


def rowstats(P, Q):
    """214's `disp214.rowstats`, verbatim."""
    d = P - Q
    rel = float(np.linalg.norm(d) / np.linalg.norm(P))
    per = np.linalg.norm(d, axis=1) / np.maximum(np.linalg.norm(P, axis=1), 1e-300)
    cs = ((P * Q).sum(1) / np.maximum(np.linalg.norm(P, axis=1) * np.linalg.norm(Q, axis=1),
                                      1e-300))
    return dict(rel_frobenius=rel, rel_row_median=float(np.median(per)),
                rel_row_max=float(per.max()), cos_median=float(np.median(cs)),
                cos_min=float(cs.min()), rows_moved=int((per > 1e-12).sum()),
                rows=int(P.shape[0]))


def chal_effects(p, c4, op, pm, h1=None):
    """The shipped construction on the challenge panel, at one `pm` and one stack."""
    pre = consensus(p, pm_blend=pm, weights=None) * ALPHA
    whitened = op.apply(pre, p["wcols"], GAMMA, LAM)
    del pre
    frac_cd4, _ = solve_fraction(whitened, c4, CD4_F)
    if h1 is None:
        m, fr_u = c4, 0.0
    else:
        fr_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        assert fr_u > 0.0, "the H1 union sized to zero on the challenge panel"
        m = merged(c4, frac_cd4, h1, fr_u)
        frac_cd4 = 1.0            # `merged` has already applied both fractions
    final, _ = E062.add_post_whiten_residual(whitened, m, frac_cd4, p["wcols"])
    del whitened
    final = np.array(final, dtype=np.float64, copy=True)
    final[:, p["panel_positions"]] = 0.0
    return final, float(fr_u)


def main():
    out = {}
    print("=== challenge panel, 3 donors ===")
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    targets = tuple(p["targets"]) if "targets" in p else None
    if targets is None:
        targets = tuple(plolo.PANEL)
    print(f"  donors {list(p['src'])}  targets {len(targets)}  "
          f"PM-available {int(np.asarray(p['pm_avail']).sum())}")

    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    carrying = int((np.asarray(h1.target_condition_count) > 0).sum())
    print(f"  H1 union ({'+'.join(UNION)}) carries {carrying}/{len(targets)} panel targets")

    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    E = {}
    for base, hh in (("cd4f034", None), ("h1union", h1)):
        for pm in (0.50, 0.00):
            E[(base, pm)], fr = chal_effects(p, c4, op, pm, hh)
            if pm == 0.50:
                out.setdefault("union_fraction", {})[base] = fr
    del op

    for base in ("cd4f034", "h1union"):
        out[f"challenge_{base}"] = rowstats(E[(base, 0.50)], E[(base, 0.00)])
        out[f"challenge_{base}"]["pm_available"] = int(np.asarray(p["pm_avail"]).sum())
        out[f"challenge_{base}"]["h1_union_carrying"] = carrying if base == "h1union" else 0

    # ---------------------------------------------------------------- h1-test, 2 donors
    print("\n=== h1-test, 2 donors (k562_gwps + orion_hek293t) ===")
    # 214's `disp214` instrument verbatim: `h1parts200.effects_for`, with every arm
    # renormalised to 169's artifact median row norm (h1parts200's own convention, which is
    # what makes this a pure SHAPE comparison on that surface).
    ph = h1parts200.parts()
    c4h = h1parts200.cd4(ph["targets"])
    oph = h1parts200.H1Operator()
    art = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                  allow_pickle=True)
    norm = float(np.median(np.linalg.norm(art["H"].astype(np.float64), axis=1)))
    Eh = {pm: h1parts200.effects_for(
        dict(pm=pm, gamma=GAMMA, lam=LAM, cd4_f=CD4_F, cd4_frac=None, weights=None),
        ph, c4h, oph, norm) for pm in (0.50, 0.00)}
    del oph, art
    out["h1test"] = rowstats(Eh[0.50], Eh[0.00])
    out["h1test"]["donors"] = list(ph["src"])
    out["h1test"]["pm_available"] = int(np.asarray(ph["pm_avail"]).sum())

    # ---------------------------------------------------------------- gates
    print("\n=== Phase 3 known-answer gates ===")
    gate("KA1 challenge/cd4f034 rel_frobenius",
         out["challenge_cd4f034"]["rel_frobenius"], 0.27059, 5e-4)
    gate("KA1 challenge/cd4f034 row median",
         out["challenge_cd4f034"]["rel_row_median"], 0.24704, 5e-4)
    gate("KA1 challenge/cd4f034 row-cos median",
         out["challenge_cd4f034"]["cos_median"], 0.97081, 5e-4)
    gate("KA1 challenge PM-available", out["challenge_cd4f034"]["pm_available"], 251, 0)
    gate("KA2 h1-test rel_frobenius", out["h1test"]["rel_frobenius"], 0.41486, 5e-4)
    gate("KA2 h1-test row median", out["h1test"]["rel_row_median"], 0.36445, 5e-4)
    gate("KA2 h1-test row-cos median", out["h1test"]["cos_median"], 0.93314, 5e-4)
    gate("KA2 h1-test PM-available", out["h1test"]["pm_available"], 79, 0)
    d214 = out["h1test"]["rel_frobenius"] / out["challenge_cd4f034"]["rel_frobenius"]
    gate("KA3 214's deflator 1.533", d214, 1.533, 3e-3)

    if FAILURES:
        print("\nPHASE 3 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # ---------------------------------------------------------------- the new number
    d231 = out["h1test"]["rel_frobenius"] / out["challenge_h1union"]["rel_frobenius"]
    out["deflator_cd4f034"] = d214
    out["deflator_h1union"] = d231
    out["deflator_ratio"] = d231 / d214

    print("\n" + "=" * 92)
    print("THE DEFLATOR, RE-MEASURED ON THE SHIPPED STACK")
    print("=" * 92)
    print(f"{'surface / stack':<28} {'||dE||/||E||':>12} {'row median':>11} "
          f"{'row-cos med':>12} {'rows moved':>11}")
    for k, lbl in (("h1test", "h1-test (2 donors)"),
                   ("challenge_cd4f034", "challenge / cd4f034"),
                   ("challenge_h1union", "challenge / h1union  <-- SHIPPED")):
        s = out[k]
        print(f"{lbl:<28} {s['rel_frobenius']:>12.5f} {s['rel_row_median']:>11.5f} "
              f"{s['cos_median']:>12.5f} {s['rows_moved']:>7}/{s['rows']}")
    print(f"\n  214's deflator (h1-test / cd4f034)   = {d214:.4f}   <- what 214 applied")
    print(f"  231's deflator (h1-test / h1union)   = {d231:.4f}   <- what the SHIPPED stack needs")
    print(f"  the extra dilution h1union adds      = {d231 / d214:.4f}x")
    print(f"\n  h1union's union residual sizes to fraction "
          f"{out['union_fraction']['h1union']:.6f} on the challenge panel, carrying "
          f"{carrying}/{len(targets)} targets.")

    json.dump(out, open(WORK / "disp231.json", "w"), indent=2)
    print(f"\nwrote {WORK / 'disp231.json'}")


if __name__ == "__main__":
    main()

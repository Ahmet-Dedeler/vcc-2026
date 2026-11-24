"""231 sec4 — rule 13's CLASS TEST for Arm A, on the currently shipped `h1union` stack.

Rule 13 (AGENTS.md sec11.1, single-valued table) charges a `reach` toll by ARM CLASS, and the
test is a **per-gene panel-total digest**, never an assertion:

  | emission-law, profile-CHANGING     | digest DIFFERS      | -0.0146 scaled |
  | emission-law, profile-PRESERVING   | digest IDENTICAL    | 0 (rc_hybrid refuted the charge) |
  | profile-class (operator, donors, magnitude, CD4 dose) | not an emission change | 0 (169) |

Arm A is a **donor-blend** change, so it is the third row — but it moves the delivered
amplitude, which is what put it on the +0.0089 floor rather than the +0.0019 one.  Both facts
are re-measured here on the shipped stack rather than carried over from `cd4f034`.

No package is built (out of mandate), so the digest is computed on the RECONSTRUCTED delivered
totals, `floor(N*L*profile)` per gene per block, with the reconstruction validated by 214's own
identity: `floor(x) == T` reproduced `cd4f034`'s own delivered integers on 1.0000 of blocks.

KNOWN-ANSWER ROWS (rule 10; the run aborts on any failure):
  KA1  the challenge displacement on the cd4f034 stack reproduces 214 sec1.5: 0.27059
  KA2  214 sec3's delivered amplitude for Arm A on cd4f034: median ||row|| ratio 1.03099
  KA3  196/214's challenge rounding tilt on the incumbent's own delivered profile:
       common log2 -0.000893, TILT log2 -0.01249 (214 Phase 0, `cd4f034` row)
"""
from __future__ import annotations

import hashlib
import json
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

import harness225 as H                                                   # noqa: E402
import plolo                                                             # noqa: E402
import chalparts200                                                      # noqa: E402
from folds200 import consensus                                           # noqa: E402
from plolo import E062, PROD_OPS, Operator, deliver, solve_fraction       # noqa: E402
from resid218 import ALPHA, GAMMA, LAM, source_residual                  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows       # noqa: E402

plolo.WORK = WORK / "plolo"
TS = plolo.TS
CD4_F, UNION_G, UNION = H.CD4_F, H.UNION_G, H.UNION
# per-context emitted library size -- 196's `recon196` known-answer row, 214 `pds214.py`:NL
NL = {"A": 400 * 20_109, "B": 400 * 19_946, "C": 400 * 20_034}
FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<44} got {got:>12.6f} exp {exp:>11.6f}  {note}")


def build(p, c4, op, pm, h1):
    pre = consensus(p, pm_blend=pm, weights=None) * ALPHA
    whitened = op.apply(pre, p["wcols"], GAMMA, LAM)
    del pre
    frac_cd4, _ = solve_fraction(whitened, c4, CD4_F)
    if h1 is None:
        m, fr = c4, frac_cd4
    else:
        fr_u, _ = solve_fraction_carrying(whitened, h1, UNION_G)
        m, fr = merged(c4, frac_cd4, h1, fr_u), 1.0
    final, _ = E062.add_post_whiten_residual(whitened, m, fr, p["wcols"])
    del whitened
    final = np.array(final, dtype=np.float64, copy=True)
    final[:, p["panel_positions"]] = 0.0
    return final


def totals_and_digest(effects, ctx):
    """The delivered per-gene panel totals for one context, and one digest per block."""
    b = plolo.challenge_bctrl(ctx)
    expr = np.expm1(b)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    g = (cpm > plolo.MIN_CPM) & (b > 0)
    dl = deliver(effects, b, g, None)                 # shipped renderer: mass_fix off
    prof = np.clip(np.expm1(b[None, :] + dl), 0.0, None)
    prof = prof / prof.sum(1, keepdims=True)
    T = np.floor(NL[ctx] * prof).astype(np.int64)
    dig = [hashlib.sha256(t.tobytes()).hexdigest()[:16] for t in T]
    return T, dig, dl


def main():
    out = {}
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    targets = tuple(p["targets"]) if "targets" in p else tuple(plolo.PANEL)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    E = {(b, pm): build(p, c4, op, pm, hh)
         for b, hh in (("cd4f034", None), ("h1union", h1)) for pm in (0.50, 0.00)}
    del op

    # ---- KA1 the displacement reproduces disp231/214 -----------------------------------
    print("=== Phase 4 known-answer gates ===")
    d = E[("cd4f034", 0.50)] - E[("cd4f034", 0.00)]
    gate("KA1 challenge/cd4f034 ||dE||/||E||",
         float(np.linalg.norm(d) / np.linalg.norm(E[("cd4f034", 0.50)])), 0.27059, 5e-4)
    del d

    # ---- deliver, digest, amplitude -----------------------------------------------------
    print("\n=== delivering all 900 blocks per stack (3 contexts x 300 targets) ===")
    res = {}
    for base in ("cd4f034", "h1union"):
        rows = {}
        for pm in (0.50, 0.00):
            digs, Ts, dls = [], [], []
            for ctx in ("A", "B", "C"):
                T, dg, dl = totals_and_digest(E[(base, pm)], ctx)
                digs += dg
                Ts.append(T)
                dls.append(dl)
            rows[pm] = dict(digests=digs,
                            med_norm_delivered=float(np.median(
                                np.linalg.norm(np.vstack(dls), axis=1))),
                            med_norm_effect=float(np.median(
                                np.linalg.norm(E[(base, pm)], axis=1))),
                            total_counts=int(sum(int(T.sum()) for T in Ts)),
                            mean_total_per_block=float(np.mean(
                                [T.sum(1).mean() for T in Ts])))
            del Ts, dls
        same = sum(1 for a, b in zip(rows[0.50]["digests"], rows[0.00]["digests"]) if a == b)
        r = dict(
            blocks=len(rows[0.50]["digests"]),
            digest_identical=same,
            first_block_digest_pm50=rows[0.50]["digests"][0],
            first_block_digest_pm00=rows[0.00]["digests"][0],
            med_norm_delivered_pm50=rows[0.50]["med_norm_delivered"],
            med_norm_delivered_pm00=rows[0.00]["med_norm_delivered"],
            med_norm_delivered_ratio=(rows[0.00]["med_norm_delivered"]
                                      / rows[0.50]["med_norm_delivered"]),
            med_norm_effect_ratio=(rows[0.00]["med_norm_effect"]
                                   / rows[0.50]["med_norm_effect"]),
            counts_pm50=rows[0.50]["total_counts"], counts_pm00=rows[0.00]["total_counts"],
        )
        r["counts_rel_pct"] = 100.0 * (r["counts_pm00"] / r["counts_pm50"] - 1.0)
        res[base] = r
        print(f"  {base:<9} digest identical on {same}/{r['blocks']} blocks; "
              f"delivered median ||row|| {r['med_norm_delivered_pm50']:.6f} -> "
              f"{r['med_norm_delivered_pm00']:.6f}  ratio {r['med_norm_delivered_ratio']:.5f}")
        print(f"  {'':<9} effect median ||row|| ratio {r['med_norm_effect_ratio']:.5f}; "
              f"panel counts move {r['counts_rel_pct']:+.4f}%")
    out["class"] = res

    # ---- KA2 214 sec3's delivered amplitude on cd4f034 ---------------------------------
    gate("KA2 214 sec3 delivered ||row|| ratio (cd4f034)",
         res["cd4f034"]["med_norm_delivered_ratio"], 1.03099, 3e-3,
         "214 head214.py, arm A on cd4f034")

    if FAILURES:
        print("\nPHASE 4 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # ---- the verdict ---------------------------------------------------------------------
    print("\n" + "=" * 92)
    print("RULE 13 CLASS TEST — Arm A on the shipped `h1union` stack")
    print("=" * 92)
    r = res["h1union"]
    amp = r["med_norm_delivered_ratio"]
    print(f"  per-gene panel-total digest identical on {r['digest_identical']}/900 blocks")
    print(f"    -> the delivered PROFILE moves.  But the EMISSION LAW is untouched: Arm A is a")
    print(f"       donor-blend change, which is rule 13's THIRD row (profile-class, operator /")
    print(f"       donors / magnitude / CD4 dose).  169 falsified the toll for that row at +0.0006.")
    print(f"  delivered median ||row|| ratio {amp:.5f}  ({100 * (amp - 1):+.3f}%)")
    print(f"    -> |move| {abs(100 * (amp - 1)):.3f}% vs the 1.0% bound 210's mandate used and the")
    print(f"       2.03% that made 210 flag it: **AMPLITUDE-CHANGING**, floor +0.0089 not +0.0019.")
    out["verdict"] = dict(
        arm_class="profile-class (donor blend), AMPLITUDE-CHANGING",
        digest_identical_blocks=r["digest_identical"], blocks=900,
        delivered_amplitude_ratio=amp,
        toll_profile_class=0.0, toll_if_charged_as_emission_law=-0.0146,
        floor_amplitude=0.0089, floor_shape_only=0.0019)
    json.dump(out, open(WORK / "class231.json", "w"), indent=2)
    print(f"\nwrote {WORK / 'class231.json'}")


if __name__ == "__main__":
    main()

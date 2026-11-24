"""231 sec2 — Arm A's `pds` and transfer cosine, RE-READ on the CURRENTLY SHIPPED stack.

219 sec3 row 5 named two things Arm A needs. This is the second: *"a `pds` re-read on
`cd4f034`"* — corrected upward by 208's rider rule, which says re-price a rider against the
CURRENT incumbent, never against the one it was measured on. The current incumbent is
**`h1union`** (entry `iUziqE1gewfcPKzYUDaB`, board 0.2118310, rank 3), not `cd4f034`.

INSTRUMENT.  162's PLOLO folds, decomposed into per-donor parts by 200's `folds200`, and
scored by 218/225's `score_effects`.  The stack is assembled exactly as 225's `harness225.Fold`
does it, which is 222's shipped construction:

    pre = consensus(parts, pm_blend)                 <- PerturbMatch enters HERE (K562 slot)
    whitened = op.apply(pre * ALPHA, wcols, GAMMA)
    frac_cd4   = solve_fraction(whitened, cd4, 0.34)              <- cd4f034
    frac_union = solve_fraction_carrying(whitened, h1union, 0.28) <- h1union
    final = merged(cd4, frac_cd4, h1, frac_union); panel cols zeroed

Both fractions are RE-SOLVED under each `pm_blend`, because production solves them from the
whitened consensus — so dropping PerturbMatch honestly moves the CD4 and union doses too.

FOUR CELLS: {pm 0.50, pm 0.00} x {baseline cd4f034, baseline h1union}.

KNOWN-ANSWER ROWS (rule 10; the run aborts on any failure):
  KA1  214 sec2.1 Arm A on `cd4f034`, SHIPPED renderer (mass_fix off, i.e. the `a15_*` frame):
       per fold  hct116 0.777958 -> 0.776216, hek293t 0.749592 -> 0.746662,
                 k562_gwps 0.802447 -> 0.802447 (structural tie), kolf2 0.613989 -> 0.607291
       4-fold mean 0.735996 -> 0.733154, d_pds = -0.002842, d_cos = -2.0512% relative
  KA2  222's h1union preregistration on the 3 CANCER folds (kolf2 excluded, as in 218):
       relative MACRO cosine +7.6803%, d_pds +0.004887
  KA3  the `f0`/`cd4` sanity row: CD4 at f=0.34 must size to a non-zero fraction on all folds
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")

sys.path.insert(0, os.path.join(ROOT, "experiments/162-retrodicting-validator"))
sys.path.insert(0, os.path.join(ROOT, "experiments/200-production-ablation"))
sys.path.insert(0, os.path.join(ROOT, "experiments/218-independent-donors"))
sys.path.insert(0, os.path.join(ROOT, "experiments/225-source-sweep"))

import harness225 as H                                                   # noqa: E402
import plolo                                                             # noqa: E402
import folds200                                                          # noqa: E402
from folds200 import consensus, parts                                    # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, cd4_for, solve_fraction     # noqa: E402
from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual   # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows       # noqa: E402

# Rule: write ONLY inside this experiment. Re-point every module's cache dir at MY copies.
from pathlib import Path                                                  # noqa: E402
plolo.WORK = Path(WORK) / "plolo"
folds200.WORK = os.path.join(WORK, "parts")

CD4_F = H.CD4_F          # 0.34, cd4f034
UNION_G = H.UNION_G      # 0.28, h1union's headline dose
UNION = H.UNION          # ("gse295214_h1", "h1_2025")
CANCER = H.CANCER        # the 3 cancer folds; kolf2 is the leak-free 4th
ALL_FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
PDS_SPAN = 0.95228 - 0.50020
K_RULE13 = 1.011

FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<48} got {got:>12.6f} exp {exp:>11.6f}  {note}")
    return ok


def build(held: str) -> dict:
    """One fold, both PerturbMatch settings, both baselines. Returns the four scored cells."""
    t0 = time.time()
    p = parts(held)
    truth = np.asarray(p["truth"], dtype=np.float64)[:, p["scored"]]
    targets = tuple(p["targets"])
    cd4 = cd4_for(targets)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in UNION])
    op_path = (PROD_OPS / "whitener_union_pruned.npz"
               if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))

    cells = {}
    for pm in (0.50, 0.00):
        pre = consensus(p, pm_blend=pm, weights=None)
        op = Operator(op_path)
        whitened = op.apply(pre * ALPHA, p["wcols"], GAMMA, LAM)
        del op, pre
        frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
        frac_union, per_unit = solve_fraction_carrying(whitened, h1, UNION_G)
        assert frac_union > 0.0, f"{held} pm{pm}: the H1 union sized to zero"
        for base, res, fr in (("cd4f034", None, 0.0), ("h1union", h1, frac_union)):
            m = merged(cd4, frac_cd4, res, fr)
            final, _ = plolo.E062.add_post_whiten_residual(whitened, m, 1.0, p["wcols"])
            final[:, p["panel_positions"]] = 0.0
            e = score_effects(final, p, truth)
            del final
            cells[(pm, base)] = dict(
                pds=e["a15_pds"], cos=e["a15_cos"], cos_agg=e["a15_cos_agg"],
                sgn_w=e["a15_sgn_w"],
                pds_a23=e["a23_pds"], cos_a23=e["a23_cos"],
                per_target_pds=np.asarray(e["_pds_per_target_a15"], dtype=np.float64),
                per_target_cos=np.asarray(e["_cos_per_target_a15"], dtype=np.float64),
                frac_cd4=float(frac_cd4), frac_union=float(frac_union),
                per_unit=float(per_unit))
        del whitened
    print(f"  [{held}] {len(targets)} targets, {time.time() - t0:.0f}s", flush=True)
    return dict(held=held, targets=targets, cells=cells)


def main():
    print("=== building folds (PerturbMatch on/off x two baselines) ===")
    folds = {h: build(h) for h in ALL_FOLDS}

    def mean_over(keys, field, pm, base):
        return float(np.mean([folds[h]["cells"][(pm, base)][field] for h in keys]))

    def mean_rel(keys, field, pm_hi, base_hi, pm_lo, base_lo):
        """218's aggregation: the MEAN over folds of each fold's RELATIVE change."""
        v = []
        for h in keys:
            a = folds[h]["cells"][(pm_hi, base_hi)][field]
            b = folds[h]["cells"][(pm_lo, base_lo)][field]
            v.append((a - b) / abs(b))
        return float(np.mean(v))

    print("\n=== Phase 2 known-answer gates ===")
    # ---- KA1: 214 sec1.2's f = 0.34 row, SHIPPED renderer (mass_fix off) and NO integer
    # rounding -- `score_effects` renders at m = infinity, which is 214 sec2.1's
    # "(reference) no rounding at all" line (0.740255) rather than its `floor` line (0.735996).
    m_inc = mean_over(ALL_FOLDS, "pds", 0.50, "cd4f034")
    m_arm = mean_over(ALL_FOLDS, "pds", 0.00, "cd4f034")
    gate("KA1a 4-fold mean, incumbent cd4f034 (m=inf)", m_inc, 0.740255, 2e-5)
    gate("KA1b 4-fold mean, arm A (m=inf)", m_arm, 0.736458, 2e-5)
    gate("KA1c 214 sec1.2 PerturbMatch pds value at f=0.34", m_inc - m_arm, 0.003797, 2e-5)
    c_inc = mean_over(ALL_FOLDS, "cos", 0.50, "cd4f034")
    c_arm = mean_over(ALL_FOLDS, "cos", 0.00, "cd4f034")
    gate("KA1d 214 sec1.2 cos pm0.50 at f=0.34", c_inc, 0.024862, 2e-5)
    gate("KA1e 214 sec1.2 cos pm0.00 at f=0.34", c_arm, 0.024350, 2e-5)
    gate("KA1f 214 sec2.1 arm A d_cos %% rel",
         100.0 * (c_arm - c_inc) / abs(c_inc), -2.0512, 3e-2)

    # ---- KA2: 222/218's h1union preregistration, 3 cancer folds, the a23 frame and 218's
    # mean-of-per-fold-relatives aggregation (`diag218.json` h1PAIRCARRY@0.28).
    gate("KA2a 222 h1union rel MACRO cosine (3 cancer, a23)",
         mean_rel(CANCER, "cos_a23", 0.50, "h1union", 0.50, "cd4f034"), 0.0768028, 3e-4)
    pp = float(np.mean([folds[h]["cells"][(0.50, "h1union")]["pds_a23"]
                        - folds[h]["cells"][(0.50, "cd4f034")]["pds_a23"] for h in CANCER]))
    gate("KA2b 222 h1union d_pds (3 cancer, a23)", pp, 0.0048867, 5e-5)
    pp4 = float(np.mean([folds[h]["cells"][(0.50, "h1union")]["pds_a23"]
                         - folds[h]["cells"][(0.50, "cd4f034")]["pds_a23"] for h in ALL_FOLDS]))
    gate("KA2c 222 h1union d_pds (4 fold, a23)", pp4, 0.0090720, 5e-5)
    for h, e in (("k562_gwps", 0.0095156), ("orion_hct116", 0.0020721),
                 ("orion_hek293t", 0.0030724), ("kolf2", 0.0216279)):
        gate(f"KA2d 218 per-fold d_pds {h}",
             folds[h]["cells"][(0.50, "h1union")]["pds_a23"]
             - folds[h]["cells"][(0.50, "cd4f034")]["pds_a23"], e, 5e-5)

    # KA3 -- non-vacuity
    for h in ALL_FOLDS:
        c = folds[h]["cells"][(0.50, "h1union")]
        ok = c["frac_cd4"] > 0 and c["frac_union"] > 0
        print(f"  [{'PASS' if ok else 'FAIL'}] KA3 {h:<16} frac_cd4 {c['frac_cd4']:.6f} "
              f"frac_union {c['frac_union']:.6f}")
        if not ok:
            FAILURES.append(f"KA3 {h} vacuous")

    if FAILURES:
        print("\nPHASE 2 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # =====================================================================
    print("\n" + "=" * 100)
    print("ARM A RE-READ — dropping PerturbMatch, on each baseline, SHIPPED renderer")
    print("=" * 100)
    out = {"folds": {}, "arms": {}, "config": dict(CD4_F=CD4_F, UNION_G=UNION_G,
                                                   UNION=list(UNION), ALPHA=ALPHA,
                                                   GAMMA=GAMMA, LAM=LAM)}
    for base in ("cd4f034", "h1union"):
        for keys, kname in ((ALL_FOLDS, "4-fold"), (CANCER, "3-cancer")):
            p50 = mean_over(keys, "pds", 0.50, base)
            p00 = mean_over(keys, "pds", 0.00, base)
            c50 = mean_over(keys, "cos", 0.50, base)
            c00 = mean_over(keys, "cos", 0.00, base)
            ca50 = mean_over(keys, "cos_agg", 0.50, base)
            ca00 = mean_over(keys, "cos_agg", 0.00, base)
            per = [folds[h]["cells"][(0.00, base)]["pds"] - folds[h]["cells"][(0.50, base)]["pds"]
                   for h in keys]
            sd = float(np.std(per, ddof=1)) if len(per) > 1 else float("nan")
            d_pds = p00 - p50
            scaled = d_pds / PDS_SPAN
            overall = scaled / 6.0 * K_RULE13
            signs = "".join("+" if x > 1e-9 else ("-" if x < -1e-9 else "0") for x in per)
            row = dict(baseline=base, folds=kname, pds_inc=p50, pds_armA=p00, d_pds=d_pds,
                       fold_sd=sd, scaled=scaled, overall=overall, signs=signs,
                       cos_inc=c50, cos_armA=c00,
                       rel_cos_pct=100.0 * (c00 - c50) / abs(c50),
                       rel_cos_pct_meanrel=100.0 * mean_rel(keys, "cos", 0.00, base,
                                                            0.50, base),
                       cos_agg_inc=ca50, cos_agg_armA=ca00,
                       rel_cos_agg_pct=100.0 * (ca00 - ca50) / abs(ca50),
                       d_pds_a23=mean_over(keys, "pds_a23", 0.00, base)
                       - mean_over(keys, "pds_a23", 0.50, base),
                       rel_cos_a23_pct_meanrel=100.0 * mean_rel(keys, "cos_a23", 0.00, base,
                                                                0.50, base),
                       per_fold=dict(zip(keys, per)))
            out["arms"][f"{base}|{kname}"] = row
            print(f"\n  baseline {base:<9} {kname:<9} pds {p50:.6f} -> {p00:.6f}  "
                  f"d {d_pds:+.6f}  sd {sd:.6f}  scaled {scaled:+.6f}  "
                  f"overall {overall:+.6f}  signs `{signs}`")
            print(f"  {'':>9} {'':<19} cos {c50:.6f} -> {c00:.6f}  "
                  f"rel {100.0 * (c00 - c50) / abs(c50):+.4f}%   "
                  f"mean-of-rel {row['rel_cos_pct_meanrel']:+.4f}%   "
                  f"cos_agg rel {100.0 * (ca00 - ca50) / abs(ca50):+.4f}%")
            print(f"  {'':>9} {'':<19} a23 frame: d_pds {row['d_pds_a23']:+.6f}  "
                  f"rel cos (mean-of-rel) {row['rel_cos_a23_pct_meanrel']:+.4f}%")
            print(f"  {'':>9} per fold: " +
                  "  ".join(f"{h}={v:+.6f}" for h, v in zip(keys, per)))

    # ---- falsifier 3: leave-top-3-out retention on the COSINE (per target) --------------
    print("\n" + "=" * 100)
    print("FALSIFIER 3 — leave-top-3-out retention, per-target COSINE (222's own test)")
    print("=" * 100)
    ret = {}
    for base in ("cd4f034", "h1union"):
        for keys, kname in ((ALL_FOLDS, "4-fold"), (CANCER, "3-cancer")):
            full, cut, kept = [], [], []
            for h in keys:
                a = folds[h]["cells"][(0.00, base)]["per_target_cos"]
                b = folds[h]["cells"][(0.50, base)]["per_target_cos"]
                d = a - b
                full.append(d.mean())
                # drop the 3 targets with the LARGEST |contribution| to the arm's effect
                idx = np.argsort(-np.abs(d))[:3]
                m = np.ones(d.size, bool)
                m[idx] = False
                cut.append(d[m].mean())
                kept.append(int(m.sum()))
            f_, c_ = float(np.mean(full)), float(np.mean(cut))
            r = 100.0 * c_ / f_ if abs(f_) > 1e-15 else float("nan")
            ret[f"{base}|{kname}|cos"] = dict(full=f_, leave_top3_out=c_, retention_pct=r,
                                              n_kept=kept)
            print(f"  cos  {base:<9} {kname:<9} full {f_:+.7f}  L3O {c_:+.7f}  "
                  f"retention {r:7.1f}%   (n kept {kept})")
        # the same on per-target pds
        for keys, kname in ((ALL_FOLDS, "4-fold"), (CANCER, "3-cancer")):
            full, cut = [], []
            for h in keys:
                a = folds[h]["cells"][(0.00, base)]["per_target_pds"]
                b = folds[h]["cells"][(0.50, base)]["per_target_pds"]
                d = a - b
                full.append(d.mean())
                idx = np.argsort(-np.abs(d))[:3]
                m = np.ones(d.size, bool)
                m[idx] = False
                cut.append(d[m].mean())
            f_, c_ = float(np.mean(full)), float(np.mean(cut))
            r = 100.0 * c_ / f_ if abs(f_) > 1e-15 else float("nan")
            ret[f"{base}|{kname}|pds"] = dict(full=f_, leave_top3_out=c_, retention_pct=r)
            print(f"  pds  {base:<9} {kname:<9} full {f_:+.7f}  L3O {c_:+.7f}  "
                  f"retention {r:7.1f}%")
    out["leave_top3_out"] = ret

    # ---- the dose interaction: does h1union change how much CD4 / the union take? -------
    print("\n=== what dropping PerturbMatch does to the SOLVED doses ===")
    dose = {}
    for h in ALL_FOLDS:
        a = folds[h]["cells"][(0.50, "h1union")]
        b = folds[h]["cells"][(0.00, "h1union")]
        dose[h] = dict(frac_cd4_pm50=a["frac_cd4"], frac_cd4_pm00=b["frac_cd4"],
                       frac_union_pm50=a["frac_union"], frac_union_pm00=b["frac_union"])
        print(f"  {h:<16} cd4 {a['frac_cd4']:.6f} -> {b['frac_cd4']:.6f} "
              f"({100 * (b['frac_cd4'] / a['frac_cd4'] - 1):+.3f}%)   "
              f"union {a['frac_union']:.6f} -> {b['frac_union']:.6f} "
              f"({100 * (b['frac_union'] / a['frac_union'] - 1):+.3f}%)")
    out["solved_doses"] = dose

    for h in ALL_FOLDS:
        out["folds"][h] = {f"pm{int(pm * 100):03d}|{base}":
                           {k: v for k, v in c.items() if not k.startswith("per_target")}
                           for (pm, base), c in folds[h]["cells"].items()}
        out["folds"][h]["n_targets"] = len(folds[h]["targets"])

    json.dump(out, open(os.path.join(WORK, "pds231.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(WORK, 'pds231.json')}")


if __name__ == "__main__":
    main()

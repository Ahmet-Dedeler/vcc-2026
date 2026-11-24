"""209 — assemble the deflation screen into board units, with the falsifiers.

Instruments and conversions, all borrowed unchanged from 204/162/197 so the numbers are
comparable to the arms this experiment is adjudicating against:

  `cos`   PRIMARY (mandate).  PLOLO macro cosine, paired target-cluster bootstrap on the
          per-target cosine, 4,000 draws.  204 §8 measured this instrument at +-0.20%
          relative = +-0.00006 overall and recommended it over `pds` for a direction arm.
  `pds`   PLOLO, 4 folds, a23 render.  overall = mean_fold d_pds / 0.452 / 6 * 1.011.
  `mse`   197's identity, priced FROM the cosine (096's `claim` route), both under the
          relative transfer (204's default) and the conservative absolute transfer.
  `nmae`  187 §5's exact identity: with alpha fixed, d raw_nmae = -2 * W_hi * d q_w,
          W_hi = 0.0153, and d overall = -d raw_nmae / (SPAN_NMAE * 6).
  `k_eff` 159 §5.5's law 1/cos^2(k) = 66.5 + 1920/k, on the board->surface anchor 197 §4.1
          uses for a15 (x2.078).  This is the mechanism claim.
  `q`     the ordering half of `reach`/`fid`.  Recorded, NOT priced (rule 13: `reach` is
          not locally priceable, 2/6 signs).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

SPAN_PDS = 0.452
CALIB = 1.011
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
LEAKFREE = ("kolf2",)                 # 162: the only leak-free fold
NB = 4000
SEED = 20260903

# 197's inversion at a15's own terms (204's `analyse204.py`, unchanged)
SUM_D = 0.22914103411323902
KC = 0.0054102686242422674
CLAIM = 0.24614814267726884
A_TERM = 0.2331583274832182
B_TERM = 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249                   # 197, rc_hybrid: our measured board direction

# 187 §5 / 204 §6
W_HI = 0.0153
SPAN_NMAE = 0.60                      # AGENTS.md §4: nmae span 0.57-0.63, midpoint
BOARD_SURFACE_ANCHOR = 2.078          # 197 §4.1, a15


def mse_at(c: float) -> float:
    q = CLAIM / SUM_D
    eps = A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw: float) -> float:
    return min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0)


def pds_overall(d_raw: float) -> float:
    return d_raw / SPAN_PDS / 6.0 * CALIB


def k_eff(board_cos: float) -> float:
    c = board_cos / BOARD_SURFACE_ANCHOR
    denom = 1.0 / (c * c) - 66.5
    return float(1920.0 / denom) if denom > 0 else float("inf")


def main() -> None:
    rows = json.loads((WORK / "scores209.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in rows}
    arms: list[str] = []
    for r in rows:
        if r["arm"] not in arms:
            arms.append(r["arm"])
    folds = [f for f in FOLDS if any(k[0] == f for k in by)]

    # ------------------------------------------- KA-4: the a15 -> a23 mass-fix price ----
    mf = [by[(f, "equal")]["a23_pds"] - by[(f, "equal")]["a15_pds"] for f in folds]
    ka = dict(per_fold_raw={f: v for f, v in zip(folds, mf)},
              mean_scaled=float(np.mean(mf) / SPAN_PDS), expected_scaled=0.004014,
              rel_err=float(np.mean(mf) / SPAN_PDS / 0.004014 - 1.0),
              a15_pds_k562_fold=by[("k562_gwps", "equal")]["a15_pds"],
              a23_pds_k562_fold=by[("k562_gwps", "equal")]["a23_pds"],
              q_k562_fold=by[("k562_gwps", "equal")]["a15_sgn"],
              q_w_k562_fold=by[("k562_gwps", "equal")]["a15_sgn_w"])
    print("=== KA-4: PLOLO known-answer rows (162 / 204 §0) ===")
    for label, got, want in (("a15 pds, k562_gwps fold", ka["a15_pds_k562_fold"], 0.79578409),
                             ("a23 pds, k562_gwps fold", ka["a23_pds_k562_fold"], 0.78912131),
                             ("a15 q,   k562_gwps fold", ka["q_k562_fold"], 0.50164311),
                             ("a15 q_w, k562_gwps fold", ka["q_w_k562_fold"], 0.50487994),
                             ("a15->a23 mass fix, 4-fold scaled", ka["mean_scaled"], 0.004014)):
        print(f"  [{'PASS' if abs(got - want) < max(2e-6, 2e-4 * abs(want)) else 'FAIL'}] "
              f"{label:36s} {got:.8f}  expected {want:.8f}")
    ka["identity_max_abs_dev_from_162_pre"] = {
        f: by[(f, "equal")]["audit"]["max_abs_dev_from_162_pre"] for f in folds}
    print(f"  identity arm vs 162's cached `pre`: max abs dev "
          + " ".join(f"{f} {v:.1e}" for f, v in
                     ka['identity_max_abs_dev_from_162_pre'].items()))

    # -------------------------------------------------------------- paired bootstrap ----
    base_p = {f: np.asarray(by[(f, "equal")]["_pds_per_target_a23"]) for f in folds}
    base_c = {f: np.asarray(by[(f, "equal")]["_cos_per_target_a23"]) for f in folds}
    rng = np.random.default_rng(SEED)
    idx = {f: rng.integers(0, len(base_p[f]), (NB, len(base_p[f]))) for f in folds}
    cos_base = float(np.mean([by[(f, "equal")]["cos"] for f in folds]))

    out = {"known_answer": ka, "folds": folds, "cos_base_4fold": cos_base, "arms": {}}
    for arm in arms:
        if any((f, arm) not in by for f in folds):
            continue
        per_p = {f: np.asarray(by[(f, arm)]["_pds_per_target_a23"]) for f in folds}
        per_c = {f: np.asarray(by[(f, arm)]["_cos_per_target_a23"]) for f in folds}
        d_p_fold = {f: float(per_p[f].mean() - base_p[f].mean()) for f in folds}
        d_c_fold = {f: float(per_c[f].mean() - base_c[f].mean()) for f in folds}
        d_p = float(np.mean(list(d_p_fold.values())))
        d_c = float(np.mean(list(d_c_fold.values())))
        boot_p = np.mean([[(per_p[f][idx[f][b]] - base_p[f][idx[f][b]]).mean()
                           for b in range(NB)] for f in folds], axis=0)
        boot_c = np.mean([[(per_c[f][idx[f][b]] - base_c[f][idx[f][b]]).mean()
                           for b in range(NB)] for f in folds], axis=0)
        rel_cos = d_c / cos_base
        board_rel = OURS_COS * (1.0 + rel_cos)
        board_abs = OURS_COS + d_c
        rec = dict(
            note=by[(folds[0], arm)].get("note", ""),
            audit={f: {k: v for k, v in by[(f, arm)]["audit"].items()
                       if k in ("lam", "dim", "restore_norm", "k562_max_abs_dev",
                                "orthonormality_err", "max_abs_dev_from_162_pre")}
                   for f in folds},
            removed_energy={f: {d: a["removed_energy_frac"] for d, a in
                                by[(f, arm)]["audit"].get("per_donor", {}).items()}
                            for f in folds},
            d_pds_per_fold=d_p_fold, d_pds=d_p,
            d_pds_ci=[float(np.percentile(boot_p, 2.5)),
                      float(np.percentile(boot_p, 97.5))],
            p_pds_gt0=float((boot_p > 0).mean()),
            d_pds_leakfree=float(np.mean([d_p_fold[f] for f in LEAKFREE])),
            d_cos_per_fold=d_c_fold, d_cos=d_c,
            d_cos_ci=[float(np.percentile(boot_c, 2.5)),
                      float(np.percentile(boot_c, 97.5))],
            p_cos_gt0=float((boot_c > 0).mean()),
            rel_cos=rel_cos,
            d_cos_leakfree=float(np.mean([d_c_fold[f] for f in LEAKFREE])),
            d_q=float(np.mean([by[(f, arm)]["a15_sgn"] - by[(f, "equal")]["a15_sgn"]
                               for f in folds])),
            d_q_w=float(np.mean([by[(f, arm)]["a15_sgn_w"] - by[(f, "equal")]["a15_sgn_w"]
                                 for f in folds])),
            row_norm_ratio={f: by[(f, arm)]["median_row_norm_final"]
                            / by[(f, "equal")]["median_row_norm_final"] for f in folds},
            implied_board_cos_rel=board_rel, implied_board_cos_abs=board_abs,
            k_eff_rel=k_eff(board_rel), k_eff_abs=k_eff(board_abs),
            overall_pds=pds_overall(d_p),
            overall_mse_rel=float((score_mse(mse_at(board_rel))
                                   - score_mse(mse_at(OURS_COS))) / 6.0),
            overall_mse_abs=float((score_mse(mse_at(board_abs))
                                   - score_mse(mse_at(OURS_COS))) / 6.0),
        )
        rec["overall_nmae"] = float(2.0 * W_HI * rec["d_q_w"] / (SPAN_NMAE * 6.0))
        rec["total_rel"] = rec["overall_pds"] + rec["overall_mse_rel"] + rec["overall_nmae"]
        rec["total_abs"] = rec["overall_pds"] + rec["overall_mse_abs"] + rec["overall_nmae"]
        out["arms"][arm] = rec

    out["k_eff_incumbent"] = k_eff(OURS_COS)

    hdr = (f"{'arm':20s} " + " ".join(f"{f[:7]:>9s}" for f in folds)
           + f" {'Δpds':>9s} {'95% CI':>21s} {'p>0':>5s} | "
             f"{'Δcos':>10s} {'95% CI':>23s} {'p>0':>5s} {'relcos':>8s} "
             f"{'Δq':>8s} {'Δq_w':>8s} {'k_eff':>6s} {'TOT':>9s} {'TOT|abs':>9s}")
    print("\n=== every arm, sorted by Δcos (the PRIMARY endpoint) ===")
    print(hdr)
    print("-" * len(hdr))
    for a in sorted(out["arms"], key=lambda x: -out["arms"][x]["d_cos"]):
        r = out["arms"][a]
        print(f"{a:20s} " + " ".join(f"{r['d_pds_per_fold'][f]:+9.5f}" for f in folds)
              + f" {r['d_pds']:+9.5f} [{r['d_pds_ci'][0]:+.5f},{r['d_pds_ci'][1]:+.5f}]"
                f" {r['p_pds_gt0']:5.3f} | {r['d_cos']:+10.6f} "
                f"[{r['d_cos_ci'][0]:+.6f},{r['d_cos_ci'][1]:+.6f}] {r['p_cos_gt0']:5.3f} "
                f"{100 * r['rel_cos']:+7.2f}% {r['d_q']:+8.5f} {r['d_q_w']:+8.5f} "
                f"{r['k_eff_rel']:6.2f} {r['total_rel']:+9.5f} {r['total_abs']:+9.5f}")

    # ------------------------------------------------------------------- the maximisers
    argmax_cos = max(out["arms"], key=lambda a: out["arms"][a]["d_cos"])
    argmax_pds = max(out["arms"], key=lambda a: out["arms"][a]["d_pds"])
    print(f"\nmaximiser of Δcos : {argmax_cos}  ({out['arms'][argmax_cos]['d_cos']:+.6f}, "
          f"{100 * out['arms'][argmax_cos]['rel_cos']:+.2f}% relative)")
    print(f"maximiser of Δpds : {argmax_pds}  ({out['arms'][argmax_pds]['d_pds']:+.5f})")
    out["argmax"] = dict(cos=argmax_cos, pds=argmax_pds)

    # ----------------------------------------------------------- the strength sweeps ---
    print("\n=== the strength sweeps (mandate Phase 2) ===")
    for fam, keys in (("subC (excess subspace, C, r=2)",
                       ["equal", "subC_l025", "subC_l05", "subC_l075", "subC_l1"]),
                      ("subW (excess subspace, WIDE, r=2)",
                       ["equal", "subW_l025", "subW_l05", "subW_l075", "subW_l1"]),
                      ("blob2 (the long-gene axis alone)",
                       ["equal", "blob2_l025", "blob2_l05", "blob2_l075", "blob2"])):
        print(f"  {fam}")
        print("    " + " ".join(f"{k:>12s}" for k in keys))
        for label, f in (("Δcos", lambda r: f"{r['d_cos']:+.6f}"),
                         ("Δpds (4-fold)", lambda r: f"{r['d_pds']:+.5f}"),
                         ("Δpds (kolf2)", lambda r: f"{r['d_pds_leakfree']:+.5f}"),
                         ("Δq_w", lambda r: f"{r['d_q_w']:+.5f}"),
                         ("k_eff", lambda r: f"{r['k_eff_rel']:.2f}")):
            print(f"    {label:14s} " + " ".join(
                f"{f(out['arms'][k]) if k in out['arms'] else '':>12s}" for k in keys))

    print(f"\nincumbent k_eff (board cos {OURS_COS}) = {out['k_eff_incumbent']:.2f}")
    (WORK / "analyse209.json").write_text(json.dumps(out, indent=2))
    print("wrote work/analyse209.json")


if __name__ == "__main__":
    main()

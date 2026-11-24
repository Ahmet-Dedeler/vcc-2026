"""204 — assemble the weighting screen into board units, with the falsifiers.

Instruments and conversions, all borrowed unchanged:
  `pds`     PLOLO, 4 folds, challenge frames A/B/C, a23 render (162, 194).
            overall = mean_fold Δpds_raw / SPAN_PDS(0.452) / 6 * 1.011      (rule 13)
  interval  paired target-cluster bootstrap on the per-target rank statistic, resampling
            targets within each fold, 4000 draws (194 §5.2's protocol).
  `mse`     through the transfer cosine and 197's inversion identity, NOT through a local
            `mse`: 162 measured that PLOLO cannot price `mse` (LOO RMS 2.5x worse than
            predicting zero).  A relative move in PLOLO's macro cosine is carried onto the
            board cosine proportionally and priced by 197's `mse_at(cos)`.
  `nmae`    `q_w`, the |t|-weighted per-gene sign accuracy (187 §5: `nmae` moves only
            through `q_w` once alpha is fixed, and alpha is fixed here), a15 render --
            194 §0.3 showed the a23 mass repair puts a 1e-3 numerical floor under the sign
            channel, so `q`/`q_w` are read through a15 for every arm.
  `reach`/`fid` ordering halves: `q` (unweighted sign accuracy), a15 render.
  `jac`     pinned by construction at fixed `keep_top` (the call set is not a function of
            the donor weights once the profile is delivered) -- stated, not measured.

KNOWN-ANSWER: the a15 -> a23 mass-fix price must come back at +0.004014 scaled `pds`
(162 §2), from the same four folds and the same pre-committed unweighted mean.
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

# 197's inversion, a15's own terms (work/invert197.json)
SUM_D = 0.22914103411323902
KC = 0.0054102686242422674
CLAIM = 0.24614814267726884
A_TERM = 0.2331583274832182
B_TERM = 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249            # our measured board direction (197, rc_hybrid)


def mse_at(c: float) -> float:
    q = CLAIM / SUM_D
    eps = A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw: float) -> float:
    return min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0)


def pds_overall(d_raw: float) -> float:
    return d_raw / SPAN_PDS / 6.0 * CALIB


def main() -> None:
    rows = json.loads((WORK / "scores204.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in rows}
    arms = []
    for r in rows:
        if r["arm"] not in arms:
            arms.append(r["arm"])
    folds = [f for f in FOLDS if any(k[0] == f for k in by)]

    # ---------------- known-answer: the a15 -> a23 mass-fix price ----------------------
    mf = [by[(f, "equal")]["a23_pds"] - by[(f, "equal")]["a15_pds"] for f in folds]
    ka = dict(per_fold_raw={f: v for f, v in zip(folds, mf)},
              mean_raw=float(np.mean(mf)), mean_scaled=float(np.mean(mf) / SPAN_PDS),
              expected_scaled=0.004014,
              rel_err=float(np.mean(mf) / SPAN_PDS / 0.004014 - 1.0))
    print(f"[known-answer] a15->a23 mass fix: {ka['mean_scaled']:+.6f} scaled "
          f"against 162's +0.004014  ({100*ka['rel_err']:+.2f}%)")

    base = {f: np.asarray(by[(f, "equal")]["_pds_per_target_a23"]) for f in folds}
    rng = np.random.default_rng(20260903)
    NB = 4000
    idx = {f: rng.integers(0, len(base[f]), (NB, len(base[f]))) for f in folds}

    out = {"known_answer_mass_fix": ka, "folds": folds, "arms": {}}
    for arm in arms:
        if any((f, arm) not in by for f in folds):
            continue
        per = {f: np.asarray(by[(f, arm)]["_pds_per_target_a23"]) for f in folds}
        d_fold = {f: float(per[f].mean() - base[f].mean()) for f in folds}
        d_mean = float(np.mean([d_fold[f] for f in folds]))
        boot = np.mean([[(per[f][idx[f][b]] - base[f][idx[f][b]]).mean()
                         for b in range(NB)] for f in folds], axis=0)
        d_cos = float(np.mean([by[(f, arm)]["cos"] - by[(f, "equal")]["cos"]
                               for f in folds]))
        cos_base = float(np.mean([by[(f, "equal")]["cos"] for f in folds]))
        rel_cos = d_cos / cos_base
        board_cos = OURS_COS * (1.0 + rel_cos)
        d_mse_ovr = (score_mse(mse_at(board_cos)) - score_mse(mse_at(OURS_COS))) / 6.0
        rec = dict(
            weights=by[(folds[0], arm)]["audit"]["weights"],
            note=by[(folds[0], arm)].get("note", ""),
            per_fold_pds={f: by[(f, arm)]["pds"] for f in folds},
            d_pds_per_fold=d_fold, d_pds=d_mean,
            d_pds_sd=float(np.std([d_fold[f] for f in folds], ddof=1)),
            ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            p_gt0=float((boot > 0).mean()),
            folds_positive=int(sum(v > 0 for v in d_fold.values())),
            d_cos=d_cos, rel_cos=rel_cos, implied_board_cos=board_cos,
            d_sse=float(np.mean([by[(f, arm)]["sse"] - by[(f, "equal")]["sse"]
                                 for f in folds])),
            d_q=float(np.mean([by[(f, arm)]["a15_sgn"] - by[(f, "equal")]["a15_sgn"]
                               for f in folds])),
            d_q_w=float(np.mean([by[(f, arm)]["a15_sgn_w"] - by[(f, "equal")]["a15_sgn_w"]
                                 for f in folds])),
            overall_pds=pds_overall(d_mean),
            overall_mse_via_cos=float(d_mse_ovr),
        )
        rec["overall_pds_plus_mse"] = rec["overall_pds"] + rec["overall_mse_via_cos"]
        out["arms"][arm] = rec

    (WORK / "analyse204.json").write_text(json.dumps(out, indent=2))

    order = sorted(out["arms"], key=lambda a: -out["arms"][a]["d_pds"])
    hdr = (f"{'arm':16s} {'w (k562,hct,hek)':26s} " + " ".join(f"{f[:7]:>9s}" for f in folds)
           + f" {'meanΔpds':>10s} {'sd':>8s} {'95% CI':>20s} {'p>0':>6s} "
             f"{'Δcos':>9s} {'Δq_w':>9s} {'ovr':>9s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for a in order:
        r = out["arms"][a]
        w = r["weights"]
        v = np.array([w["k562_gwps"], w["orion_hct116"], w["orion_hek293t"]])
        v = v / v.sum() * 3
        print(f"{a:16s} {'(' + ', '.join(f'{x:.3f}' for x in v) + ')':26s} " +
              " ".join(f"{r['d_pds_per_fold'][f]:+9.5f}" for f in folds) +
              f" {r['d_pds']:+10.5f} {r['d_pds_sd']:8.5f}"
              f"  [{r['ci95'][0]:+.5f},{r['ci95'][1]:+.5f}] {r['p_gt0']:6.3f}"
              f" {r['d_cos']:+9.5f} {r['d_q_w']:+9.5f} {r['overall_pds']:+9.5f}")
    print("\nabsolute `equal` baseline: " +
          "  ".join(f"{f} {by[(f, 'equal')]['pds']:.5f}" for f in folds))
    print("wrote work/analyse204.json")


if __name__ == "__main__":
    main()

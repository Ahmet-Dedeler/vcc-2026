"""204 Phase 2/3 — the six-member price, the cosine interval, and the mechanism regression.

THE SIX MEMBERS, and the instrument each is priced on (rule 12).

`pds`   PLOLO, 4 folds, challenge frames, a23 render; paired target bootstrap.
        overall = mean_fold Δpds_raw / 0.452 / 6 × 1.011.

`mse`   NOT priced on a local `mse` -- 162 measured PLOLO's `mse` LOO RMS at 2.5× worse
        than predicting zero.  Priced through 096/197's identity instead.  A donor-mixing
        arm holds the emission law, `alpha`, the spread policy and `keep_top` fixed, so
        `A`, `B`, `Kc` and `claim` are unchanged and the ONLY term that moves is `X`:

            MSE_c − 1 = (A_c + Kc_c − min(claim_c, B_c) − 2X) / ΣD ,  X = cos·√ΣD·mean√A_c
            d raw_mse / d cos = −2√(q(1+ε)) = −2.0826   →   +0.3653 overall per +1.0 cosine

        which is AGENTS.md §8.2's +0.00365 overall per +0.01 of cosine, re-derived.  That
        `A` really is unchanged is checked, not assumed: `median_row_norm_final` is
        reported for every arm.  PLOLO's macro cosine is carried onto the board cosine
        RELATIVELY (Δcos_board = 0.079249 · Δcos_PLOLO / cos_PLOLO) -- stated as the
        assumption it is, and the absolute Δcos is reported alongside.

`nmae`  187 §5's exact gene-wise identity `nmae = 1 − 2 q_w W_hi − (2 q_lo − 1) K_lo + K_hi`.
        With `alpha` fixed, `K_lo`, `K_hi` and `W_hi` are fixed and the only live term is
        `q_w`:  d raw_nmae = −2 W_hi · d q_w, W_hi = 0.0153 (187's measured value on the
        reference pool; the board's implied 0.012 is more favourable, so 0.0153 is the
        conservative choice).  Scaled = −Δraw / 0.5963, /6 for overall.
        The h1-test emitted route (188) is NOT run: see RESULTS.md §6.

`reach` / `fid`  ordering halves read off `q` (unweighted per-gene sign accuracy, a15
        render).  Yield/`n_real` halves unreadable on any local harness (rule 9).
        No `reach` toll: this is a PROFILE-class arm, and 169 falsified the toll there
        (+0.0006), the −0.0146 toll being an EMISSION-LAW toll (audit 181).

`jac`   the call set is `keep_top` applied to the delivered profile; the arm does not
        touch `keep_top`, and 192 owns the call set.  Reported as unmeasured with the
        mechanism, per rule 12's narrow exemption.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

SPAN_PDS, CALIB = 0.452, 1.011
SPAN_NMAE = 0.5963
W_HI = 0.0153
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")

SUM_D, KC, CLAIM = 0.22914103411323902, 0.0054102686242422674, 0.24614814267726884
A_TERM, B_TERM = 0.2331583274832182, 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249


def mse_at(c):
    q, eps = CLAIM / SUM_D, A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw):
    return min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0)


def main() -> None:
    rows = json.loads((WORK / "scores204.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in rows}
    arms = list(dict.fromkeys(r["arm"] for r in rows))
    folds = [f for f in FOLDS if (f, "equal") in by]

    NB = 4000
    rng = np.random.default_rng(20260903)
    base_p = {f: np.asarray(by[(f, "equal")]["_pds_per_target_a23"]) for f in folds}
    base_c = {f: np.asarray(by[(f, "equal")]["_cos_per_target_a23"]) for f in folds}
    idx = {f: rng.integers(0, len(base_p[f]), (NB, len(base_p[f]))) for f in folds}
    cos_base = float(np.mean([base_c[f].mean() for f in folds]))

    out = {"cos_base_plolo": cos_base, "arms": {}}
    for arm in arms:
        if any((f, arm) not in by for f in folds):
            continue
        p = {f: np.asarray(by[(f, arm)]["_pds_per_target_a23"]) for f in folds}
        c = {f: np.asarray(by[(f, arm)]["_cos_per_target_a23"]) for f in folds}
        d_p = float(np.mean([p[f].mean() - base_p[f].mean() for f in folds]))
        d_c = float(np.mean([c[f].mean() - base_c[f].mean() for f in folds]))
        bp = np.mean([[(p[f][idx[f][b]] - base_p[f][idx[f][b]]).mean() for b in range(NB)]
                      for f in folds], axis=0)
        bc = np.mean([[(c[f][idx[f][b]] - base_c[f][idx[f][b]]).mean() for b in range(NB)]
                      for f in folds], axis=0)
        d_qw = float(np.mean([by[(f, arm)]["a15_sgn_w"] - by[(f, "equal")]["a15_sgn_w"]
                              for f in folds]))
        d_q = float(np.mean([by[(f, arm)]["a15_sgn"] - by[(f, "equal")]["a15_sgn"]
                             for f in folds]))
        board_cos = OURS_COS * (1.0 + d_c / cos_base)
        o_pds = d_p / SPAN_PDS / 6.0 * CALIB
        o_mse = (score_mse(mse_at(board_cos)) - score_mse(mse_at(OURS_COS))) / 6.0
        o_nmae = -(-2.0 * W_HI * d_qw) / SPAN_NMAE / 6.0
        # the conservative alternative transfer: carry the ABSOLUTE PLOLO cosine move onto
        # the board rather than the relative one.  The two bracket the `mse` price.
        o_mse_abs = (score_mse(mse_at(OURS_COS + d_c))
                     - score_mse(mse_at(OURS_COS))) / 6.0
        norms = [by[(f, arm)]["median_row_norm_final"] /
                 by[(f, "equal")]["median_row_norm_final"] for f in folds]
        out["arms"][arm] = dict(
            weights=by[(folds[0], arm)]["audit"]["weights"],
            d_pds=d_p, d_pds_ci=[float(np.percentile(bp, 2.5)),
                                 float(np.percentile(bp, 97.5))],
            d_pds_p=float((bp > 0).mean()),
            d_pds_per_fold={f: float(p[f].mean() - base_p[f].mean()) for f in folds},
            d_cos=d_c, d_cos_ci=[float(np.percentile(bc, 2.5)),
                                 float(np.percentile(bc, 97.5))],
            d_cos_p=float((bc > 0).mean()),
            d_cos_per_fold={f: float(c[f].mean() - base_c[f].mean()) for f in folds},
            rel_cos=d_c / cos_base, implied_board_cos=board_cos,
            d_q=d_q, d_q_w=d_qw,
            row_norm_ratio=[float(x) for x in norms],
            overall=dict(pds=o_pds, mse=o_mse, mse_abs_transfer=o_mse_abs,
                         nmae=o_nmae, total_priced=o_pds + o_mse + o_nmae,
                         total_priced_abs=o_pds + o_mse_abs + o_nmae))

    # ---- mechanism: what degree of freedom does the score actually respond to? ---------
    A = []
    y_p, y_c = [], []
    names = []
    for arm, r in out["arms"].items():
        if arm.endswith("_pt"):
            continue                    # per-target arms are not a point in weight space
        w = np.array([r["weights"]["k562_gwps"], r["weights"]["orion_hct116"],
                      r["weights"]["orion_hek293t"]], dtype=float)
        w = w / w.sum()
        # coordinates: log K562 share, and log of the HCT116:HEK293T split at fixed K562
        lk = np.log(w[0] / np.sqrt(w[1] * w[2]))
        lo = np.log(w[1] / w[2])
        A.append([lk, lo, lk * lk, lo * lo])
        y_p.append(r["d_pds"])
        y_c.append(r["d_cos"])
        names.append(arm)
    A = np.column_stack([np.ones(len(A))] + list(np.array(A).T))
    for label, y in (("pds", np.array(y_p)), ("cos", np.array(y_c))):
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ beta
        ss = 1 - ((y - pred) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-30)
        # partial R2 of each coordinate
        part = {}
        for j, nm in ((1, "k562_share"), (2, "orion_split")):
            keep = [i for i in range(A.shape[1]) if i not in (j, j + 2)]
            b2, *_ = np.linalg.lstsq(A[:, keep], y, rcond=None)
            r2b = 1 - ((y - A[:, keep] @ b2) ** 2).sum() / max(
                ((y - y.mean()) ** 2).sum(), 1e-30)
            part[nm] = float(ss - r2b)
        out[f"regression_{label}"] = dict(
            n_arms=len(names), arms=names, intercept=float(beta[0]),
            slope_log_k562_share=float(beta[1]), slope_log_orion_split=float(beta[2]),
            quad_k562=float(beta[3]), quad_orion=float(beta[4]),
            r2=float(ss), partial_r2=part)
        print(f"[regression {label}] n={len(names)}  R2={ss:.4f}   "
              f"d/dlog(K562 share) {beta[1]:+.5f}   d/dlog(HCT:HEK) {beta[2]:+.5f}   "
              f"partial R2: k562 {part['k562_share']:.4f}  "
              f"orion_split {part['orion_split']:.4f}")

    (WORK / "six204.json").write_text(json.dumps(out, indent=2))

    order = sorted(out["arms"], key=lambda a: -out["arms"][a]["overall"]["total_priced"])
    print(f"\n{'arm':16s} {'Δpds':>9s} {'p>0':>6s} {'Δcos':>10s} "
          f"{'Δcos 95% CI':>22s} {'p>0':>6s} {'ovr pds':>9s} {'ovr mse':>9s} "
          f"{'ovr nmae':>9s} {'TOTAL':>9s} {'TOT|abs':>9s}")
    print("-" * 128)
    for a in order:
        r = out["arms"][a]
        o = r["overall"]
        print(f"{a:16s} {r['d_pds']:+9.5f} {r['d_pds_p']:6.3f} {r['d_cos']:+10.6f} "
              f"[{r['d_cos_ci'][0]:+.6f},{r['d_cos_ci'][1]:+.6f}] {r['d_cos_p']:6.3f} "
              f"{o['pds']:+9.5f} {o['mse']:+9.5f} {o['nmae']:+9.6f} "
              f"{o['total_priced']:+9.5f} {o['total_priced_abs']:+9.5f}")
    rn = np.array([out["arms"][a]["row_norm_ratio"] for a in out["arms"]])
    print(f"\ndelivered median row-norm ratio vs `equal`: min {rn.min():.6f} "
          f"max {rn.max():.6f}  (the `mse` identity holds A fixed; this is the check)")
    print("wrote work/six204.json")


if __name__ == "__main__":
    main()

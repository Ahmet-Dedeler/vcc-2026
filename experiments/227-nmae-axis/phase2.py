"""227 Phase 2 — what the per-row additive offset IS, and whether it is reachable.

Phase 1 found the only channel bigger than the Aginglab gap: an L1-optimal CONSTANT added to
our predicted lfc, one scalar per perturbation, worth **+0.0542 raw** on h1-test with no per-gene
information at all.  This script asks three questions in order and stops being interesting if
any of them answers no.

  Q1  WHAT IS IT?  A constant in lfc space is a multiplicative rescale of our predicted
      arithmetic-mean-CPM profile.  Under the scorer's per-cell CPM normalisation the profile
      sums to 1e6, so a shift of one gene set must be paid for by another: the offset is a
      COMPOSITIONAL shift.  Is `delta_p` a property of the gate only, or of the whole row?  Is
      it the truth's compositional shift that we under-predict because we shrink?

  Q2  IS IT REACHABLE?  Fit `delta_p` from observables only -- our own predicted profile, our
      own call set -- leave-one-perturbation-out, and price the honest version.  Shuffle null
      on every fitted arm.

  Q3  WHAT DOES IT COST?  `pds` and 197's energy-weighted aggregate cosine on the SAME run,
      every time.  An arm that buys `nmae` by losing cosine is the failure mode most of the
      board's `nmae` table is in and is reported as a failure, not a trade.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface           # noqa: E402
from phase1 import cos_agg, gate_of, nmae_of, rows_of   # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
TABLE = os.path.join(E169, "de169.q1.00_mf.npz")
EFFECT = os.path.join(E169, "work/h1test-effect.npz")
RNG = np.random.default_rng(227_000)


def add_on_gate(P0, G, d):
    """Add a per-row constant to the gate columns only."""
    P = P0.copy()
    for i, (p, s) in enumerate(G):
        P[p, s] = P0[p, s] + d[i]
    return P


def add_on_all(P0, G, d, cols):
    """Add a per-row constant to EVERY scored column -- the realisable shape, since the gate
    is the reference's and we do not know it."""
    P = P0.copy()
    for i, (p, s) in enumerate(G):
        P[p, cols] = P0[p, cols] + d[i]
    return P


def price(name, P, P0, G, T, S, cols, base, b):
    return dict(arm=name, nmae=nmae_of(G, P, T), d_nmae=base - nmae_of(G, P, T),
                pds=S.pds(pred_lfc=P), d_pds=S.pds(pred_lfc=P) - b["pds"],
                cos=cos_agg(P, T, cols),
                d_cos_pct=100 * (cos_agg(P, T, cols) / b["cos_agg"] - 1.0))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    S = Surface(TABLE)
    T = S.real_lfc.astype(np.float64)
    P0 = np.nan_to_num(S.pred_lfc.astype(np.float64))
    G = gate_of(S)
    cols = np.flatnonzero(np.isfinite(T).all(0) & np.isfinite(P0).all(0))
    base = nmae_of(G, P0, T)
    b = dict(nmae=base, pds=S.pds(), cos_agg=cos_agg(P0, T, cols))
    n = len(G)

    # ------------------------------------------------------------------ Q1  what is it
    d_orc = np.array([np.median(T[p, s] - P0[p, s]) for p, s in G])       # the oracle offset
    med_t_gate = np.array([np.median(T[p, s]) for p, s in G])
    med_h_gate = np.array([np.median(P0[p, s]) for p, s in G])
    med_t_all = np.array([np.median(T[p, cols]) for p, _ in G])           # whole-row shift
    med_h_all = np.array([np.median(P0[p, cols]) for p, _ in G])
    # mass-weighted mean lfc: pinned to ~0 by the CPM normalisation on both sides
    wm_t = np.array([float(np.average(T[p, cols])) for p, _ in G])
    wm_h = np.array([float(np.average(P0[p, cols])) for p, _ in G])

    def cc(a, c):
        return float(np.corrcoef(a, c)[0, 1])

    out["Q1_what_is_it"] = dict(
        n_rows=n,
        d_oracle=dict(mean=float(d_orc.mean()), sd=float(d_orc.std()),
                      q=[float(np.quantile(d_orc, x)) for x in (0.05, 0.5, 0.95)]),
        med_t_gate=dict(mean=float(med_t_gate.mean()), sd=float(med_t_gate.std())),
        med_h_gate=dict(mean=float(med_h_gate.mean()), sd=float(med_h_gate.std())),
        med_t_all=dict(mean=float(med_t_all.mean()), sd=float(med_t_all.std())),
        med_h_all=dict(mean=float(med_h_all.mean()), sd=float(med_h_all.std())),
        mean_lfc_t=dict(mean=float(wm_t.mean()), sd=float(wm_t.std())),
        mean_lfc_h=dict(mean=float(wm_h.mean()), sd=float(wm_h.std())),
        corr_d_vs_med_t_gate=cc(d_orc, med_t_gate),
        corr_d_vs_med_t_all=cc(d_orc, med_t_all),
        corr_med_t_gate_vs_med_t_all=cc(med_t_gate, med_t_all),
        corr_med_h_all_vs_med_t_all=cc(med_h_all, med_t_all),
        corr_med_h_gate_vs_med_t_gate=cc(med_h_gate, med_t_gate),
        # if we merely SHRINK the truth's compositional shift, this ratio is the shrinkage
        slope_t_on_h_all=float(np.polyfit(med_h_all, med_t_all, 1)[0]),
        slope_t_on_h_gate=float(np.polyfit(med_h_gate, med_t_gate, 1)[0]),
    )

    # ------------------------------------------------------------------ Q2/Q3  arms
    arms = []
    arms.append(price("oracle offset, gate only", add_on_gate(P0, G, d_orc),
                      P0, G, T, S, cols, base, b))
    arms.append(price("oracle offset, all genes", add_on_all(P0, G, d_orc, cols),
                      P0, G, T, S, cols, base, b))

    # -- the realisable family: rescale OUR OWN row-median shift by a single gain g ----------
    gains = {}
    for g in (1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0):
        d = g * med_h_all
        gains[f"g={g:g}"] = price(f"our median-shift x{g:g} (all genes)",
                                  add_on_all(P0, G, d, cols), P0, G, T, S, cols, base, b)
    out["gain_sweep"] = gains
    best_g = min(gains, key=lambda k: gains[k]["nmae"])
    arms.append(gains[best_g])

    # -- LOPO linear fit of delta on observables --------------------------------------------
    X = np.column_stack([np.ones(n), med_h_all, med_h_gate,
                         np.array([np.abs(P0[p, cols]).mean() for p, _ in G])])
    d_lopo = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        w, *_ = np.linalg.lstsq(X[m], d_orc[m], rcond=None)
        d_lopo[i] = X[i] @ w
    out["lopo_fit"] = dict(corr_pred_vs_oracle=cc(d_lopo, d_orc),
                           r2=float(1 - ((d_lopo - d_orc) ** 2).sum()
                                    / ((d_orc - d_orc.mean()) ** 2).sum()),
                           sd_pred=float(d_lopo.std()), sd_oracle=float(d_orc.std()))
    arms.append(price("LOPO-fitted offset (all genes)", add_on_all(P0, G, d_lopo, cols),
                      P0, G, T, S, cols, base, b))
    arms.append(price("LOPO-fitted offset (gate only)", add_on_gate(P0, G, d_lopo),
                      P0, G, T, S, cols, base, b))

    # -- shuffle nulls -----------------------------------------------------------------------
    nulls = {"oracle offset shuffled": [], "LOPO offset shuffled": []}
    for _ in range(20):
        pr = RNG.permutation(n)
        nulls["oracle offset shuffled"].append(
            nmae_of(G, add_on_all(P0, G, d_orc[pr], cols), T))
        nulls["LOPO offset shuffled"].append(
            nmae_of(G, add_on_all(P0, G, d_lopo[pr], cols), T))
    out["shuffle_null"] = {k: dict(mean_nmae=float(np.mean(v)), sd=float(np.std(v)),
                                   mean_d_nmae=float(base - np.mean(v)))
                           for k, v in nulls.items()}

    # -- the row-clip channel: 33% of rows score worse than predicting nothing ---------------
    r = rows_of(G, P0, T)
    key_norm = np.array([np.linalg.norm(P0[p, cols]) for p, _ in G])
    key_k = np.array([np.abs(P0[p, s]).sum() / np.abs(T[p, s]).sum() for p, s in G])  # oracle
    out["row_clip"] = dict(
        frac_rows_above_1=float((r > 1.0).mean()),
        oracle_clip_gain=float(base - np.minimum(r, 1.0).mean()),
        corr_row_nmae_vs_prednorm=cc(r, key_norm),
        corr_row_nmae_vs_K=cc(r, key_k),
    )

    out["base"] = b
    out["arms"] = arms
    json.dump(out, open(os.path.join(WORK, "phase2.json"), "w"), indent=1)

    # ------------------------------------------------------------------ report
    q = out["Q1_what_is_it"]
    print(f"baseline: nmae {b['nmae']:.6f}  pds {b['pds']:.6f}  cos_agg {b['cos_agg']:.6f}"
          f"   ({n} rows)")
    print("\n== Q1  what the per-row offset is ==")
    print(f"  oracle offset delta_p        mean {q['d_oracle']['mean']:+.4f}  "
          f"sd {q['d_oracle']['sd']:.4f}")
    print(f"  median TRUE lfc on the gate  mean {q['med_t_gate']['mean']:+.4f}  "
          f"sd {q['med_t_gate']['sd']:.4f}")
    print(f"  median OUR  lfc on the gate  mean {q['med_h_gate']['mean']:+.4f}  "
          f"sd {q['med_h_gate']['sd']:.4f}")
    print(f"  median TRUE lfc, ALL genes   mean {q['med_t_all']['mean']:+.4f}  "
          f"sd {q['med_t_all']['sd']:.4f}")
    print(f"  median OUR  lfc, ALL genes   mean {q['med_h_all']['mean']:+.4f}  "
          f"sd {q['med_h_all']['sd']:.4f}")
    print(f"  unweighted mean lfc  truth {q['mean_lfc_t']['mean']:+.4f}"
          f"/{q['mean_lfc_t']['sd']:.4f}   ours {q['mean_lfc_h']['mean']:+.4f}"
          f"/{q['mean_lfc_h']['sd']:.4f}")
    print(f"  corr(delta, median t on gate)      {q['corr_d_vs_med_t_gate']:+.4f}")
    print(f"  corr(delta, median t all genes)    {q['corr_d_vs_med_t_all']:+.4f}")
    print(f"  corr(median t gate, median t all)  {q['corr_med_t_gate_vs_med_t_all']:+.4f}")
    print(f"  corr(OUR median all, TRUE med all) {q['corr_med_h_all_vs_med_t_all']:+.4f}"
          f"   slope truth~ours {q['slope_t_on_h_all']:.3f}")
    print(f"  corr(OUR median gate,TRUE med gate){q['corr_med_h_gate_vs_med_t_gate']:+.4f}"
          f"   slope truth~ours {q['slope_t_on_h_gate']:.3f}")

    print("\n== gain sweep: our own median shift, multiplied by g, added to every gene ==")
    print(f"  {'g':>8s} {'nmae':>9s} {'d nmae':>9s} {'pds':>9s} {'d pds':>9s} "
          f"{'cos_agg':>9s} {'d cos %':>9s}")
    for k, v in gains.items():
        print(f"  {k:>8s} {v['nmae']:9.6f} {v['d_nmae']:+9.6f} {v['pds']:9.6f} "
              f"{v['d_pds']:+9.6f} {v['cos']:9.6f} {v['d_cos_pct']:+9.3f}")

    print("\n== arms, all six-member-relevant endpoints on the same run ==")
    print(f"  {'arm':34s} {'nmae':>9s} {'d nmae':>9s} {'pds':>9s} {'d pds':>9s} "
          f"{'cos_agg':>9s} {'d cos %':>9s}")
    for a in arms:
        print(f"  {a['arm']:34s} {a['nmae']:9.6f} {a['d_nmae']:+9.6f} {a['pds']:9.6f} "
              f"{a['d_pds']:+9.6f} {a['cos']:9.6f} {a['d_cos_pct']:+9.3f}")

    lf = out["lopo_fit"]
    print(f"\n  LOPO fit of delta on observables: corr {lf['corr_pred_vs_oracle']:+.4f}  "
          f"R2 {lf['r2']:+.4f}  sd(pred) {lf['sd_pred']:.4f} vs sd(oracle) "
          f"{lf['sd_oracle']:.4f}")
    print("\n== shuffle nulls (row labels permuted) ==")
    for k, v in out["shuffle_null"].items():
        print(f"  {k:28s} nmae {v['mean_nmae']:.6f} +- {v['sd']:.6f}   "
              f"d nmae {v['mean_d_nmae']:+.6f}")
    rc = out["row_clip"]
    print(f"\n== row-clip channel ==")
    print(f"  rows above 1.0 {100*rc['frac_rows_above_1']:.1f}%   oracle clip gain "
          f"{rc['oracle_clip_gain']:+.6f}")
    print(f"  corr(row nmae, ||pred||) {rc['corr_row_nmae_vs_prednorm']:+.4f}   "
          f"corr(row nmae, row K) {rc['corr_row_nmae_vs_K']:+.4f}")


if __name__ == "__main__":
    main()

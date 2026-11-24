"""227 Phase 0 — known-answer rows.  No instrument, no conclusions.

Five gates, in the order the mandate names them:

  KA-1  the board affine map: retrodict OUR OWN scaled `nmae` (+0.11074) from raw 0.93389,
        with the map fitted on the 657-entry board itself, and read Aginglab's row the same way.
  KA-2  188's exact `nmae` identity `q|k-1| + (1-q)(k+1)` with |t|-WEIGHTED q, and its two
        slopes -(2q-1) below the kink and +1 above it, reproduced gene-wise on h1-test.
  KA-3  187 KA-b: h1-test raw `nmae` at s = 1 on 170's table = 0.959881 (target removed).
  KA-4  188's identity188 gate: my own `nmae` implementation vs `debench.Surface.nmae()`.
  KA-5  187's h1-test amplitude response: `nmae(s)` turnover at s ~ 1.29, and the
        first-order slope d nmae / d ln alpha = -0.0211.

KA-6 is the one this experiment turns on and it is FREE, because 169 emitted the SAME arm at
five seeds: if the quasi emitter has exact gene margins and a fixed per-cell library size, the
arithmetic mean of per-cell CPM is a DETERMINISTIC function of the design and `nmae` must be
bit-identical across seeds.  If it is not, candidate 1 (emission spread) is alive and everything
below is wrong.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface  # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
BOARD = os.path.join(REPO, "data/leaderboard/2026-09-04T0433Z.json")
US, THEM, LEADER = "Ahmet Dedeler", "Aginglab.com", "cqawesome"


def nmae_parts(S, pred_lfc=None):
    """Gene-wise decomposition of `nmae` over the reference gate.

    Returns per-perturbation the pieces 188's identity needs, all |t|-weighted within a row and
    then averaged across rows the way the metric averages (mean of per-row ratios).
    """
    P = S.pred_lfc if pred_lfc is None else pred_lfc
    rows = []
    for p in range(len(S.perts)):
        s = np.flatnonzero(S.real_sig[p] & np.isfinite(S.real_lfc[p]))
        if s.size < 10:
            continue
        t = S.real_lfc[p, s].astype(np.float64)
        den = np.abs(t).sum()
        if not np.isfinite(den) or den == 0:
            continue
        h = np.where(np.isfinite(P[p, s]), P[p, s], 0.0).astype(np.float64)
        w = np.abs(t) / den                       # |t| weights, sum 1
        k = np.abs(h) / np.maximum(np.abs(t), 1e-300)
        ok = np.sign(h) == np.sign(t)
        lo, hi = k < 1, k >= 1
        K_lo, K_hi = float((w * k)[lo].sum()), float((w * k)[hi].sum())
        rows.append(dict(
            nmae=float(np.abs(h - t).sum() / den),
            n=int(s.size),
            q_w=float(w[ok].sum()),               # |t|-weighted sign accuracy (all genes)
            q_c=float(ok.mean()),                 # plain count accuracy
            K=float(np.abs(h).sum() / den),       # Sigma|h| / Sigma|t|
            K_lo=K_lo, K_hi=K_hi,
            W_hi=float(w[hi].sum()),
            # the two accuracies 187's identity actually needs:
            #   W_c_hi -- |t|-weighted mass of CORRECT-sign genes ABOVE the kink
            #   qK_lo  -- k-weighted sign accuracy BELOW the kink
            W_c_hi=float(w[ok & hi].sum()),
            qK_lo=float((w * k)[ok & lo].sum() / max(K_lo, 1e-300)),
            mean_abs_t=float(np.abs(t).mean()),
        ))
    return rows


def agg(rows, key):
    return float(np.mean([r[key] for r in rows]))


def identity_nmae(rows):
    """187 sec 1's gene-wise form, with the two weightings made explicit:

        nmae = 1 - 2 W_c_hi - (2 qK_lo - 1) K_lo + K_hi

    Derivation: split the gate by sign-correctness and by the kink at k = 1.  Correct/below
    contributes w(1-k), correct/above w(k-1), wrong contributes w(1+k); the w's sum to 1.
    `W_c_hi` is |t|-weighted over genes ABOVE the kink; `qK_lo` is k-weighted BELOW it.  This
    is an identity, not a fit: it must equal the metric to ~1e-12."""
    v = [1.0 - 2 * r["W_c_hi"] - (2 * r["qK_lo"] - 1) * r["K_lo"] + r["K_hi"]
         for r in rows]
    return float(np.mean(v))


def uniform_identity(q, k):
    """The 02:15 / CURRENT_STATE reading: p = +-k t, uniform k, |t|-weighted q."""
    return q * abs(k - 1.0) + (1.0 - q) * (k + 1.0)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}

    # ------------------------------------------------------------------ KA-1  board map
    # `score_nmae` is the mean over three contexts of (b_c - raw_c)/span_c while the reported
    # raw is the mean of raw_c, so the aggregate map is affine ONLY where the per-context spans
    # agree.  Fit it on the competitive band and quote the residual as the instrument's
    # precision -- 188 warns explicitly against doing this arithmetic naively.
    E = json.load(open(BOARD))["entries"]
    raw = np.array([e["de_wilcoxon_lfc_nmae"] for e in E], float)
    sc = np.array([e["score_nmae"] for e in E], float)
    m = (np.isfinite(raw) & np.isfinite(sc) & (raw > 0.75) & (raw < 1.10))
    A = np.polyfit(raw[m], sc[m], 1)
    span = -1.0 / A[0]
    b = -A[1] / A[0]                     # scaled = (b - raw)/span
    resid = np.abs(A[0] * raw[m] + A[1] - sc[m])

    def rows_for(name):
        r = [e for e in E if e["team_name"] == name]
        return max(r, key=lambda e: e["score_avg"]) if r else None
    us, them, ldr = rows_for(US), rows_for(THEM), rows_for(LEADER)
    out["KA1_board_map"] = dict(
        n=int(m.sum()), slope=float(A[0]), intercept=float(A[1]),
        baseline_b=float(b), span=float(span),
        max_abs_resid=float(resid.max()), med_abs_resid=float(np.median(resid)),
        p95_abs_resid=float(np.percentile(resid, 95)),
        us_raw=us["de_wilcoxon_lfc_nmae"], us_scaled_board=us["score_nmae"],
        us_scaled_model=float(A[0] * us["de_wilcoxon_lfc_nmae"] + A[1]),
        them_raw=them["de_wilcoxon_lfc_nmae"], them_scaled_board=them["score_nmae"],
        them_scaled_model=float(A[0] * them["de_wilcoxon_lfc_nmae"] + A[1]),
        leader_raw=ldr["de_wilcoxon_lfc_nmae"], leader_scaled_board=ldr["score_nmae"],
        us_overall=us["score_avg"], them_overall=them["score_avg"],
        leader_overall=ldr["score_avg"],
    )

    # ------------------------------------------------------------------ KA-3/4/2  h1-test
    tables = {
        "169.q1.00_mf (prod geometry, seed 90k)": f"{E169}/de169.q1.00_mf.npz",
        "169.q1.00_mf_s2": f"{E169}/de169.q1.00_mf_s2.npz",
        "169.q1.00_mf_s3": f"{E169}/de169.q1.00_mf_s3.npz",
        "169.q1.00_mf_s4": f"{E169}/de169.q1.00_mf_s4.npz",
        "169.q1.00_mf_s5": f"{E169}/de169.q1.00_mf_s5.npz",
        "169.KA_a15_raw": f"{E169}/de169.KA_a15_raw.npz",
        "170.de_g0.65": f"{REPO}/experiments/170-gamma-mse/work/de_g0.65.npz",
    }
    tab = {}
    for name, path in tables.items():
        if not os.path.exists(path):
            tab[name] = dict(missing=True)
            continue
        S = Surface(path)
        rows = nmae_parts(S)
        mine = float(np.mean([r["nmae"] for r in rows]))
        tab[name] = dict(
            n_perts=int(len(S.perts)), n_genes=int(S.genes.size), n_scored=len(rows),
            gate_median=float(np.median([r["n"] for r in rows])),
            nmae_debench=S.nmae(), nmae_mine=mine,
            delta_vs_debench=abs(mine - S.nmae()),
            nmae_identity=identity_nmae(rows),
            q_w=agg(rows, "q_w"), q_count=agg(rows, "q_c"),
            K=agg(rows, "K"), K_lo=agg(rows, "K_lo"), K_hi=agg(rows, "K_hi"),
            W_hi=agg(rows, "W_hi"), mean_abs_t=agg(rows, "mean_abs_t"),
            uniform_k_at_qw=uniform_identity(agg(rows, "q_w"), agg(rows, "K")),
        )
        del S
    out["h1_tables"] = tab

    # ------------------------------------------------------------------ KA-6  seed invariance
    seeds = [k for k in tab if k.startswith("169.q1.00_mf")]
    vals = [tab[k]["nmae_debench"] for k in seeds]
    out["KA6_seed_invariance"] = dict(
        seeds=seeds, nmae=vals,
        spread=float(max(vals) - min(vals)), sd=float(np.std(vals)),
    )

    # ------------------------------------------------------------------ KA-5  amplitude
    S = Surface(tables["170.de_g0.65"]) if os.path.exists(tables["170.de_g0.65"]) else None
    if S is not None:
        curve = {}
        for s in (0.5, 0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.29, 1.4, 1.6, 2.0, 3.0):
            curve[f"{s:.2f}"] = S.nmae(pred_lfc=S.pred_lfc * s)
        eps = 1e-3
        slope = (S.nmae(pred_lfc=S.pred_lfc * (1 + eps))
                 - S.nmae(pred_lfc=S.pred_lfc * (1 - eps))) / (2 * eps)
        ss = np.linspace(0.5, 3.0, 251)
        vv = [S.nmae(pred_lfc=S.pred_lfc * float(x)) for x in ss]
        out["KA5_amplitude"] = dict(
            curve=curve, dnmae_dlnalpha=float(slope),
            argmin_s=float(ss[int(np.argmin(vv))]), min_nmae=float(np.min(vv)),
            gain_raw=float(S.nmae() - np.min(vv)),
        )
        del S

    json.dump(out, open(os.path.join(WORK, "phase0.json"), "w"), indent=1)

    # ------------------------------------------------------------------ report
    k1 = out["KA1_board_map"]
    print(f"== KA-1  board affine map (fitted on {k1['n']} entries, raw in [0.75, 1.10]) ==")
    print(f"  scaled = ({k1['baseline_b']:.6f} - raw) / {k1['span']:.6f}"
          f"   |resid|: med {k1['med_abs_resid']:.2e}  p95 {k1['p95_abs_resid']:.2e}"
          f"  max {k1['max_abs_resid']:.2e}")
    for who, pre in (("us", "us"), ("Aginglab.com", "them"), ("cqawesome", "leader")):
        r, s = k1[f"{pre}_raw"], k1[f"{pre}_scaled_board"]
        mm = k1.get(f"{pre}_scaled_model")
        extra = f"  model {mm:+.5f}  err {mm-s:+.2e}" if mm is not None else ""
        print(f"  {who:14s} raw {r:.5f}  scaled {s:+.5f}{extra}")
    print()
    print("== KA-2/3/4  h1-test tables ==")
    print(f"  {'table':38s} {'gate':>5s} {'nmae':>9s} {'identity':>9s} {'|d|':>8s} "
          f"{'q_w':>6s} {'q_cnt':>6s} {'K':>7s} {'W_hi':>7s} {'E|t|':>6s} {'unif-k':>8s}")
    for name, d in tab.items():
        if d.get("missing"):
            print(f"  {name:38s}   MISSING")
            continue
        print(f"  {name:38s} {d['gate_median']:5.0f} {d['nmae_debench']:9.6f} "
              f"{d['nmae_identity']:9.6f} {abs(d['nmae_identity']-d['nmae_debench']):8.1e} "
              f"{d['q_w']:6.4f} {d['q_count']:6.4f} {d['K']:7.4f} {d['W_hi']:7.4f} "
              f"{d['mean_abs_t']:6.3f} {d['uniform_k_at_qw']:8.5f}")
    print()
    k6 = out["KA6_seed_invariance"]
    print("== KA-6  seed invariance of `nmae` under the quasi emitter (5 matched seeds) ==")
    for s, v in zip(k6["seeds"], k6["nmae"]):
        print(f"  {s:38s} {v:.12f}")
    print(f"  spread across seeds = {k6['spread']:.3e}   sd = {k6['sd']:.3e}")
    print()
    if "KA5_amplitude" in out:
        k5 = out["KA5_amplitude"]
        print("== KA-5  h1-test amplitude response (170 table, lfc scaled by s) ==")
        print("  " + "  ".join(f"{s}:{v:.6f}" for s, v in k5["curve"].items()))
        print(f"  d nmae/d ln alpha = {k5['dnmae_dlnalpha']:+.4f}   "
              f"argmin s = {k5['argmin_s']:.3f}   gain = {k5['gain_raw']:.5f} raw")


if __name__ == "__main__":
    main()

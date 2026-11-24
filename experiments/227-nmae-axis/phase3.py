"""227 Phase 3 — the honest arms, the falsifiers, and what Aginglab's number would take.

  A. The (q_w, K) frontier measured on OUR OWN magnitude distribution, not on the uniform-k
     identity: hold |h| fixed and re-sign genes to hit a target |t|-weighted sign accuracy, then
     hold the signs and rescale.  This converts "raw nmae 0.88902" into a statement about sign
     accuracy ON THE REFERENCE'S GATE, which is what the metric actually reads.
  B. The one remaining channel with an observable key: per-row amplitude.  Phase 2 measured
     corr(row nmae, ||pred||) = +0.42, so rows we predict loudly score worse.  Fitted
     leave-one-perturbation-out, priced with `pds` and cosine on the same run, with a shuffle null.
  C. Falsifier 3 -- leave-top-3-out retention (222 sec 9's protocol: delete the three targets
     carrying the most gain and report what survives).
  D. A correction to the board reading the mandate is built on: `mse` is CLAMPED to [0,1], so a
     team whose raw `mse` exceeds the baseline reports score_mse = 0 exactly and 197's inversion
     returns a cosine that was never measured.  Recompute the -0.394 on the teams where the
     inversion is actually defined.
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
from debench import Surface                                    # noqa: E402
from phase1 import cos_agg, gate_of, nmae_of, rows_of          # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
TABLE = os.path.join(E169, "de169.q1.00_mf.npz")
BOARD = os.path.join(REPO, "data/leaderboard/2026-09-04T0433Z.json")
RNG = np.random.default_rng(227_003)
MSE_B, MSE_SPAN = 0.987862879, 0.950162566
COS_A, COS_B = 1.023611, 2.074495


def qw_of(G, P, T):
    v = []
    for p, s in G:
        w = np.abs(T[p, s])
        v.append(w[np.sign(P[p, s]) == np.sign(T[p, s])].sum() / w.sum())
    return float(np.mean(v))


def K_of(G, P, T):
    return float(np.mean([np.abs(P[p, s]).sum() / np.abs(T[p, s]).sum() for p, s in G]))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    S = Surface(TABLE)
    T = S.real_lfc.astype(np.float64)
    P0 = np.nan_to_num(S.pred_lfc.astype(np.float64))
    G = gate_of(S)
    cols = np.flatnonzero(np.isfinite(T).all(0) & np.isfinite(P0).all(0))
    base = nmae_of(G, P0, T)
    b = dict(nmae=base, pds=S.pds(), cos_agg=cos_agg(P0, T, cols),
             q_w=qw_of(G, P0, T), K=K_of(G, P0, T), credit=1.0 - base)
    out["base"] = b
    n = len(G)

    # ------------------------------------------------------------------ A. the (q_w, K) frontier
    # ANCHORED at our own prediction: to raise q_w we flip a random subset of the genes we
    # currently get WRONG to correct (and vice versa to lower it), leaving every magnitude and
    # every other gene untouched.  Independent re-signing would also destroy the coupling
    # between |h| and correctness that our prediction actually has, and would not pass through
    # our own point -- which is the whole use of a frontier.
    front = []
    for f in (-0.30, -0.15, 0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00):
        P = P0.copy()
        for p, s in G:
            ok = np.sign(P0[p, s]) == np.sign(T[p, s])
            pool = np.flatnonzero(~ok if f >= 0 else ok)
            k = int(round(abs(f) * pool.size))
            if k:
                pick = s[RNG.permutation(pool)[:k]]
                P[p, pick] = -P0[p, pick]
        v = nmae_of(G, P, T)
        front.append(dict(kind="fix-wrong", param=f, q_w=qw_of(G, P, T), K=K_of(G, P, T),
                          nmae=v, credit=1 - v, credit_ratio=(1 - v) / b["credit"],
                          cos=cos_agg(P, T, cols), pds=S.pds(pred_lfc=P)))
    # THE GATE THE MANDATE ACTUALLY SETS: is `nmae` reachable at FIXED cosine?  `nmae` weights a
    # gene by |t|; the aggregate cosine weights it by |h||t| energy.  So if the axes separate at
    # all, the separation is maximal when we fix the wrong genes in ascending |h| (cheapest in
    # cosine) or in descending |t|/|h| (the nmae-per-unit-cosine optimum).  Both are ORACLE
    # orderings -- they read the truth -- so they are ceilings on the separation, not arms.
    for tag, keyf in (("fix-lo-|h|", lambda h, t: np.abs(h)),
                      ("fix-hi-|h|", lambda h, t: -np.abs(h)),
                      ("fix-hi-t/h", lambda h, t: -np.abs(t) / np.maximum(np.abs(h), 1e-12))):
        for f in (0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00):
            P = P0.copy()
            for p, s in G:
                ok = np.sign(P0[p, s]) == np.sign(T[p, s])
                pool = np.flatnonzero(~ok)
                if pool.size == 0:
                    continue
                order = pool[np.argsort(keyf(P0[p, s][pool], T[p, s][pool]), kind="stable")]
                k = int(round(f * pool.size))
                if k:
                    pick = s[order[:k]]
                    P[p, pick] = -P0[p, pick]
            v = nmae_of(G, P, T)
            front.append(dict(kind=tag, param=f, q_w=qw_of(G, P, T), K=K_of(G, P, T),
                              nmae=v, credit=1 - v, credit_ratio=(1 - v) / b["credit"],
                              cos=cos_agg(P, T, cols), pds=S.pds(pred_lfc=P)))

    # third axis: MAGNITUDE information.  Replace |h| on the gate by a noisy |t| at a controlled
    # rank correlation, keeping our own signs -- how much magnitude knowledge would it take?
    for lam in (0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0):
        P = P0.copy()
        for p, s in G:
            a = np.abs(P0[p, s])
            tgt = np.abs(T[p, s]) * (a.sum() / max(np.abs(T[p, s]).sum(), 1e-300))
            P[p, s] = np.sign(np.where(P0[p, s] == 0, 1.0, P0[p, s])) * ((1 - lam) * a
                                                                        + lam * tgt)
        v = nmae_of(G, P, T)
        front.append(dict(kind="magnitude", param=lam, q_w=qw_of(G, P, T), K=K_of(G, P, T),
                          nmae=v, credit=1 - v, credit_ratio=(1 - v) / b["credit"],
                          cos=cos_agg(P, T, cols), pds=S.pds(pred_lfc=P)))
    for sc in (0.5, 0.8, 1.0, 1.3, 1.6, 2.0, 3.0):
        P = P0 * sc
        v = nmae_of(G, P, T)
        front.append(dict(kind="rescale", param=sc, q_w=qw_of(G, P, T), K=K_of(G, P, T),
                          nmae=v, credit=1 - v, credit_ratio=(1 - v) / b["credit"],
                          cos=cos_agg(P, T, cols), pds=S.pds(pred_lfc=P)))
    out["frontier"] = front
    # what q_w reproduces Aginglab's credit RATIO of 1.6788 (0.11098 / 0.06611)?
    tgt = (1 - 0.88902) / (1 - 0.93389)
    out["q_w_required"] = dict(target_credit_ratio=float(tgt), q_w_ours=b["q_w"],
                               cos_ours=b["cos_agg"])
    for kind in ("fix-wrong", "fix-lo-|h|", "fix-hi-|h|", "fix-hi-t/h", "magnitude"):
        rs = [f for f in front if f["kind"] == kind]
        cr = np.array([f["credit_ratio"] for f in rs])
        o = np.argsort(cr)
        out["q_w_required"][kind] = dict(
            q_w_needed=float(np.interp(tgt, cr[o], np.array([f["q_w"] for f in rs])[o])),
            cos_at_target=float(np.interp(tgt, cr[o],
                                          np.array([f["cos"] for f in rs])[o])),
            pds_at_target=float(np.interp(tgt, cr[o],
                                          np.array([f["pds"] for f in rs])[o])),
            reachable=bool(cr.max() >= tgt))

    # ------------------------------------------------------------------ B. per-row amplitude
    r = rows_of(G, P0, T)
    key = np.array([np.linalg.norm(P0[p, cols]) for p, _ in G])
    ss = np.linspace(0.0, 2.0, 81)

    def fit_apply(idx_fit, idx_all, kv):
        """One scalar per tercile of the key, fitted on `idx_fit`, applied to `idx_all`."""
        ed = np.quantile(kv[idx_fit], [1 / 3, 2 / 3])
        binf, bina = np.searchsorted(ed, kv[idx_fit]), np.searchsorted(ed, kv[idx_all])
        sc = np.ones(3)
        for t in range(3):
            m = np.flatnonzero(binf == t)
            if m.size == 0:
                continue
            cand = [sum(np.abs(P0[G[idx_fit[i]][0], G[idx_fit[i]][1]] * x
                               - T[G[idx_fit[i]][0], G[idx_fit[i]][1]]).sum()
                        / np.abs(T[G[idx_fit[i]][0], G[idx_fit[i]][1]]).sum() for i in m)
                    for x in ss]
            sc[t] = ss[int(np.argmin(cand))]
        return sc[bina]

    idx = np.arange(n)
    s_lopo = np.array([fit_apply(np.delete(idx, i), np.array([i]), key)[0] for i in idx])
    P = P0.copy()
    for i, (p, s) in enumerate(G):
        P[p, :] = P0[p, :] * s_lopo[i]
    arm_amp = dict(arm="per-row amplitude, LOPO tercile of ||pred||",
                   nmae=nmae_of(G, P, T), d_nmae=base - nmae_of(G, P, T),
                   pds=S.pds(pred_lfc=P), d_pds=S.pds(pred_lfc=P) - b["pds"],
                   cos=cos_agg(P, T, cols),
                   d_cos_pct=100 * (cos_agg(P, T, cols) / b["cos_agg"] - 1),
                   scalars=[float(x) for x in np.unique(s_lopo)])
    nulls = []
    for _ in range(20):
        pr = RNG.permutation(n)
        Q = P0.copy()
        for i, (p, s) in enumerate(G):
            Q[p, :] = P0[p, :] * s_lopo[pr][i]
        nulls.append(nmae_of(G, Q, T))
    arm_amp["shuffle_null_nmae"] = float(np.mean(nulls))
    arm_amp["shuffle_null_d"] = float(base - np.mean(nulls))
    out["arm_per_row_amplitude"] = arm_amp

    # ------------------------------------------------------------------ C. leave-top-3-out
    def retention(P):
        per = np.array([np.abs(P0[p, s] - T[p, s]).sum() / np.abs(T[p, s]).sum()
                        - np.abs(P[p, s] - T[p, s]).sum() / np.abs(T[p, s]).sum()
                        for p, s in G])
        tot = float(per.mean())
        o = np.argsort(-per)
        rec = {}
        for k in (1, 3, 5):
            keep = np.setdiff1d(np.arange(n), o[:k])
            rec[f"k={k}"] = dict(retained=float(per[keep].sum() / n),
                                 retained_pct=float(100 * per[keep].sum() / n / tot)
                                 if tot != 0 else float("nan"),
                                 topk_share_pct=float(100 * per[o[:k]].sum() / per.sum())
                                 if per.sum() != 0 else float("nan"))
        return dict(total=tot, by_k=rec)

    Pamp = P0.copy()
    for i, (p, s) in enumerate(G):
        Pamp[p, :] = P0[p, :] * s_lopo[i]
    d_orc = np.array([np.median(T[p, s] - P0[p, s]) for p, s in G])
    Porc = P0.copy()
    for i, (p, s) in enumerate(G):
        Porc[p, s] = P0[p, s] + d_orc[i]
    out["retention"] = dict(per_row_amplitude=retention(Pamp),
                            oracle_per_row_offset=retention(Porc),
                            global_scalar_0p80=retention(P0 * 0.80))

    # ------------------------------------------------------------------ D. board clamp
    E = json.load(open(BOARD))["entries"]
    tm = {}
    for e in E:
        t = e["team_name"]
        if t not in tm or e["score_avg"] > tm[t]["score_avg"]:
            tm[t] = e
    top = [e for e in tm.values() if e["score_avg"] > 0.15]
    cosv = np.array([(COS_A - e["expr_mse_unbiased_capped_norm"]) / COS_B for e in top])
    nmv = np.array([e["score_nmae"] for e in top])
    clamped = np.array([e["score_mse"] <= 1e-12 for e in top])
    rawmse = np.array([e["expr_mse_unbiased_capped_norm"] for e in top])
    out["board_clamp"] = dict(
        n_top=len(top), n_clamped=int(clamped.sum()),
        corr_all=float(np.corrcoef(cosv, nmv)[0, 1]),
        corr_unclamped=float(np.corrcoef(cosv[~clamped], nmv[~clamped])[0, 1]),
        n_unclamped=int((~clamped).sum()),
        clamped_raw_mse_min=float(rawmse[clamped].min()) if clamped.any() else None,
        mse_baseline=MSE_B,
        clamped_mean_nmae_sc=float(nmv[clamped].mean()) if clamped.any() else None,
        unclamped_mean_nmae_sc=float(nmv[~clamped].mean()),
    )

    json.dump(out, open(os.path.join(WORK, "phase3.json"), "w"), indent=1)

    # ------------------------------------------------------------------ report
    print(f"h1-test base: nmae {b['nmae']:.6f} (credit {b['credit']:.5f})  q_w {b['q_w']:.4f} "
          f"K {b['K']:.4f}  pds {b['pds']:.5f}  cos_agg {b['cos_agg']:.5f}")
    print(f"\n== A. the (q_w, K) frontier on OUR OWN magnitude distribution ==")
    print(f"  {'kind':8s} {'param':>7s} {'q_w':>7s} {'K':>7s} {'nmae':>9s} {'credit':>8s} "
          f"{'x ours':>7s} {'pds':>8s} {'cos_agg':>8s}")
    for f in front:
        print(f"  {f['kind']:8s} {f['param']:7.4f} {f['q_w']:7.4f} {f['K']:7.4f} "
              f"{f['nmae']:9.6f} {f['credit']:8.5f} {f['credit_ratio']:7.3f} "
              f"{f['pds']:8.5f} {f['cos']:8.5f}")
    qr = out["q_w_required"]
    print(f"\n  Aginglab's credit is {qr['target_credit_ratio']:.4f}x ours.  Reaching it:")
    for kind in ("fix-wrong", "fix-lo-|h|", "fix-hi-|h|", "fix-hi-t/h", "magnitude"):
        d = qr[kind]
        print(f"    via {kind:10s}: q_w {d['q_w_needed']:.4f} (ours {qr['q_w_ours']:.4f}), "
              f"cos_agg {d['cos_at_target']:.5f} = {100*(d['cos_at_target']/qr['cos_ours']-1):+.1f}% "
              f"of ours, pds {d['pds_at_target']:.5f}   reachable={d['reachable']}")

    a = arm_amp
    print(f"\n== B. per-row amplitude, LOPO tercile of ||pred|| ==")
    print(f"  nmae {a['nmae']:.6f}  d {a['d_nmae']:+.6f}   pds {a['pds']:.6f} "
          f"d {a['d_pds']:+.6f}   cos_agg {a['cos']:.6f} d {a['d_cos_pct']:+.3f}%")
    print(f"  fitted scalars {['%.3f' % x for x in a['scalars']]}   "
          f"shuffle null d {a['shuffle_null_d']:+.6f}")

    print(f"\n== C. leave-top-3-out retention on `nmae` (222 sec 9 protocol) ==")
    for nm, d in out["retention"].items():
        s = "  ".join(f"k={k.split('=')[1]}: {v['retained_pct']:6.1f}% "
                      f"(top-{k.split('=')[1]} share {v['topk_share_pct']:.1f}%)"
                      for k, v in d["by_k"].items())
        print(f"  {nm:26s} total {d['total']:+.6f}   {s}")

    bc = out["board_clamp"]
    print(f"\n== D. the board's `mse` clamp ==")
    print(f"  of {bc['n_top']} teams above overall 0.15, {bc['n_clamped']} report score_mse "
          f"= 0 EXACTLY -- their raw `mse` is at or past the baseline {bc['mse_baseline']:.6f}\n"
          f"  (min raw among them {bc['clamped_raw_mse_min']:.5f}), so 197's inversion returns "
          f"a number that was never measured.")
    print(f"  corr(implied cosine, scaled nmae): all {bc['corr_all']:+.4f}   "
          f"UNCLAMPED ONLY (n={bc['n_unclamped']}) {bc['corr_unclamped']:+.4f}")
    print(f"  mean scaled nmae: clamped {bc['clamped_mean_nmae_sc']:+.5f}   "
          f"unclamped {bc['unclamped_mean_nmae_sc']:+.5f}")


if __name__ == "__main__":
    main()

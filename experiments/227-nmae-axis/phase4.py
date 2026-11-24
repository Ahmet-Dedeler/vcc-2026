"""227 Phase 4 — INBOX-1 and INBOX-2, answered together.

INBOX-1: do not lean on 197's `mse` inversion; `pds_cosine` is reported DIRECTLY on the board, and
the retrodiction target is a SIX-TUPLE, not a member.  So this script re-prices the Phase-2
falsification against `pds` (which needs no inversion and has no jackknife term) and fits the joint
shape: high `nmae` AND lower `pds` AND positive `fid` AND negative `jac`.

INBOX-2: `q` (|t|-weighted sign accuracy on the significance-selected pool) is direction over the
POOL; `pds` is direction over ALL genes; the claim is that they come apart.  Phase 3 measured that
they do not come apart along any sign-fixing route.  This script asks the question the other way
round: **is there a construction that raises `q` while LOWERING `pds`?**  Two candidates:

  A. CONCENTRATION -- keep only our top-N genes by |h| and predict no change elsewhere.  A sparse
     prediction keeps its confident signs (where `q` lives) and throws away the broad structure
     (where `pds` lives).  This is also the natural construction for "calls more, looser".
  B. GATE-SELECTIVE accuracy -- an oracle that fixes signs on the reference's gate and breaks them
     off it.  Not shippable; it is the existence proof for whether the two axes can be decoupled
     at all, and it bounds how extreme a model would have to be.

Every arm is priced on all five members `debench` can compute deterministically, in scaled units.
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
from debench import ANCHOR, Surface, scaled                    # noqa: E402
from phase1 import cos_agg, gate_of, nmae_of                   # noqa: E402
from phase3 import qw_of                                       # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
TABLE = os.path.join(E169, "de169.q1.00_mf.npz")
NMAE_SPAN = 0.603369                     # phase0 KA-1, board-fitted
RNG = np.random.default_rng(227_005)

# the six-tuple to hit, as Aginglab MINUS us, in scaled board units (INBOX-1 item 1)
TARGET = dict(nmae=+0.07323, pds=-0.03025, fid=+0.01970, jac=-0.03295, reach=-0.02213)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    S = Surface(TABLE)
    T = S.real_lfc.astype(np.float64)
    P0 = np.nan_to_num(S.pred_lfc.astype(np.float64))
    G = gate_of(S)
    cols = np.flatnonzero(np.isfinite(T).all(0) & np.isfinite(P0).all(0))
    call0 = S.pred_sig.copy()

    def six(P, call):
        nm = nmae_of(G, P, T)
        return dict(nmae_raw=nm, nmae_sc=(1.000222 - nm) / NMAE_SPAN,
                    pds_raw=S.pds(pred_lfc=P), pds_sc=scaled("pds", S.pds(pred_lfc=P)),
                    fid_sc=scaled("fid", S.fid(call=call)),
                    jac_sc=scaled("jac", S.jac(call=call)),
                    reach_sc=scaled("reach", S.reach(key=np.abs(P), call=call)),
                    cos_agg=cos_agg(P, T, cols), q_w=qw_of(G, P, T),
                    n_pred=float(np.mean([int(call[p].sum()) for p in range(len(S.perts))])))

    b = six(P0, call0)
    out = {"base": b, "target_vs_us": TARGET, "arms": {}}

    def add(name, P, call):
        v = six(P, call)
        out["arms"][name] = dict(
            **v, **{f"d_{k}": v[f"{k}_sc"] - b[f"{k}_sc"]
                    for k in ("nmae", "pds", "fid", "jac", "reach")},
            d_cos_pct=100 * (v["cos_agg"] / b["cos_agg"] - 1), d_q_w=v["q_w"] - b["q_w"])

    # ---- A. concentration: keep top-N by |h|, predict no change elsewhere ------------------
    order = np.argsort(-np.abs(P0), axis=1)
    for N in (200, 500, 1000, 2000, 4000, 8000):
        P = np.zeros_like(P0)
        call = np.zeros_like(call0)
        for p in range(len(S.perts)):
            k = order[p, :N]
            P[p, k] = P0[p, k]
            call[p, k] = True
        add(f"concentrate top-{N} by |h|", P, call)

    # ---- A'. the mirror: WIDEN the call set at our own profile (their fid/jac shape) -------
    for N in (2000, 4000, 6000, 9000):
        call = np.zeros_like(call0)
        for p in range(len(S.perts)):
            call[p, order[p, :N]] = True
        add(f"same profile, call top-{N}", P0, call)

    # ---- C. magnitude REALLOCATION at fixed budget toward the pool the Wilcoxon selects -----
    # 188 established the reference's BH set is reproduced 0.7219 by |lfc|*sqrt(cpm) -- i.e. the
    # gate is selected by POWER, which rises with expression.  So tilting our magnitude budget
    # toward high-CPM genes moves it onto the scored pool without changing the budget (row L2 is
    # restored) and without changing any sign.  This is a SHAPE change, not amplitude.
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    cpm = np.array([E["cpm"][pos[g]] for g in S.genes], dtype=np.float64)
    u = np.maximum(cpm / np.median(cpm), 1e-9)
    for beta in (-0.30, -0.15, 0.15, 0.30, 0.50, 1.00):
        w = u ** beta
        P = P0 * w[None, :]
        n0 = np.linalg.norm(P0, axis=1, keepdims=True)
        P = P * (n0 / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300))
        add(f"CPM tilt beta={beta:+.2f} (row L2 fixed)", P, call0)

    # ---- B. gate-selective accuracy (ORACLE existence proof) --------------------------------
    ing = np.zeros(P0.shape, bool)
    for p, s in G:
        ing[p, s] = True
    for f_in, f_out in ((0.25, 0.0), (0.25, 0.10), (0.25, 0.25), (0.50, 0.25), (0.50, 0.50)):
        P = P0.copy()
        for p in range(len(S.perts)):
            wrong_in = np.flatnonzero(ing[p] & (np.sign(P0[p]) != np.sign(T[p])))
            right_out = np.flatnonzero(~ing[p] & (np.sign(P0[p]) == np.sign(T[p]))
                                       & np.isfinite(T[p]))
            for pool, f in ((wrong_in, f_in), (right_out, f_out)):
                k = int(round(f * pool.size))
                if k:
                    pick = pool[RNG.permutation(pool.size)[:k]]
                    P[p, pick] = -P0[p, pick]
        add(f"ORACLE fix {f_in:.0%} on gate / break {f_out:.0%} off gate", P, call0)

    # ---- D. the per-row LEVEL channel, dialled to THEIR nmae, priced on `pds` --------------
    # Phase 2's oracle offset is the only channel whose sign is right on every member.  Dial it
    # to exactly their +0.07323 scaled `nmae` and read what it costs on the board's own direct
    # direction metric.  A model that predicts the per-perturbation RESPONSE LEVEL would buy this.
    d_orc = np.array([np.median(T[p, s] - P0[p, s]) for p, s in G])
    lvl = []
    for gam in (0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0):
        P = P0.copy()
        for i, (p, s) in enumerate(G):
            P[p, cols] = P0[p, cols] + gam * d_orc[i]
        v = six(P, call0)
        lvl.append(dict(gamma=gam, d_nmae=v["nmae_sc"] - b["nmae_sc"],
                        d_pds=v["pds_sc"] - b["pds_sc"],
                        d_fid=v["fid_sc"] - b["fid_sc"], d_jac=v["jac_sc"] - b["jac_sc"],
                        d_cos_pct=100 * (v["cos_agg"] / b["cos_agg"] - 1)))
    nn = np.array([r["d_nmae"] for r in lvl])
    o = np.argsort(nn)
    out["level_channel"] = dict(
        sweep=lvl,
        gamma_at_their_nmae=float(np.interp(TARGET["nmae"], nn[o],
                                            np.array([r["gamma"] for r in lvl])[o])),
        d_pds_at_their_nmae=float(np.interp(TARGET["nmae"], nn[o],
                                            np.array([r["d_pds"] for r in lvl])[o])),
        their_d_pds=TARGET["pds"])

    json.dump(out, open(os.path.join(WORK, "phase4.json"), "w"), indent=1)

    # ---------------------------------------------------------------- report
    print(f"h1-test base: nmae {b['nmae_raw']:.6f} (sc {b['nmae_sc']:+.5f})  "
          f"pds {b['pds_raw']:.5f} (sc {b['pds_sc']:+.5f})  q_w {b['q_w']:.4f}  "
          f"cos_agg {b['cos_agg']:.5f}  n_pred {b['n_pred']:.0f}")
    print(f"\nTARGET (Aginglab minus us, scaled board units): "
          + "  ".join(f"{k} {v:+.5f}" for k, v in TARGET.items()))
    print(f"\n  {'arm':40s} {'d nmae':>9s} {'d pds':>9s} {'d fid':>9s} {'d jac':>9s} "
          f"{'d reach':>9s} {'d q_w':>7s} {'d cos%':>8s} {'n_pred':>7s}")
    print(f"  {'TARGET':40s} {TARGET['nmae']:+9.5f} {TARGET['pds']:+9.5f} "
          f"{TARGET['fid']:+9.5f} {TARGET['jac']:+9.5f} {TARGET['reach']:+9.5f} "
          f"{'?':>7s} {'?':>8s} {'?':>7s}")
    for name, d in out["arms"].items():
        print(f"  {name:40s} {d['d_nmae']:+9.5f} {d['d_pds']:+9.5f} {d['d_fid']:+9.5f} "
              f"{d['d_jac']:+9.5f} {d['d_reach']:+9.5f} {d['d_q_w']:+7.4f} "
              f"{d['d_cos_pct']:+8.2f} {d['n_pred']:7.0f}")

    # ---- the INBOX-1 falsification, done on `pds` instead of the mse inversion -------------
    P3 = json.load(open(os.path.join(WORK, "phase3.json")))
    fr = [f for f in P3["frontier"] if f["kind"] == "fix-lo-|h|"]
    cr = np.array([f["credit_ratio"] for f in fr])
    pv = np.array([f["pds"] for f in fr])
    o = np.argsort(cr)
    tgt = (1 - 0.88902) / (1 - 0.93389)
    pds_at = float(np.interp(tgt, cr[o], pv[o]))
    d_pds_raw = pds_at - P3["base"]["pds"]
    PL = json.load(open(os.path.join(WORK, "plolo227.json")))
    a = PL["summary"]["3 cancer folds"]["arms"]["fix lo-|h| wrong signs, f=0.4"]
    plolo_rate = a["d_pds"] / a["d_nmae"]
    obs = 0.82276 - 0.83572
    out["pds_falsification"] = dict(
        h1_predicted_d_pds_raw=d_pds_raw,
        plolo_d_pds_per_raw_nmae=plolo_rate,
        plolo_predicted_d_pds_raw=float(plolo_rate * 0.04487),
        observed_d_pds_raw=obs, us_pds_raw=0.83572, them_pds_raw=0.82276)
    json.dump(out, open(os.path.join(WORK, "phase4.json"), "w"), indent=1)
    lc = out["level_channel"]
    print(f"\n== D. the per-row LEVEL channel (oracle offset x gamma), the only one whose sign "
          f"is right on every member ==")
    print(f"  {'gamma':>7s} {'d nmae':>9s} {'d pds':>9s} {'d fid':>9s} {'d jac':>9s} "
          f"{'d cos%':>8s}")
    for r in lc["sweep"]:
        print(f"  {r['gamma']:7.2f} {r['d_nmae']:+9.5f} {r['d_pds']:+9.5f} "
              f"{r['d_fid']:+9.5f} {r['d_jac']:+9.5f} {r['d_cos_pct']:+8.2f}")
    print(f"  at THEIR nmae (+{TARGET['nmae']:.5f}): gamma {lc['gamma_at_their_nmae']:.3f}, "
          f"d pds {lc['d_pds_at_their_nmae']:+.5f}  vs their actual "
          f"{lc['their_d_pds']:+.5f}  (ratio {lc['d_pds_at_their_nmae']/lc['their_d_pds']:.2f}x)")

    print(f"\n== INBOX-1: the falsification redone on `pds_cosine`, which the board reports "
          f"DIRECTLY ==")
    print(f"  h1-test, cheapest sign-fix route at their credit ratio: "
          f"d pds_cosine {d_pds_raw:+.5f}")
    print(f"  PLOLO 3 cancer folds: {plolo_rate:+.4f} raw pds per raw nmae  ->  "
          f"{plolo_rate*0.04487:+.5f} for their 0.04487")
    print(f"  OBSERVED: their raw pds_cosine {0.82276:.5f} vs ours {0.83572:.5f} = "
          f"{obs:+.5f}")
    print(f"  the direction hypothesis is wrong by {d_pds_raw - obs:+.5f} (h1) / "
          f"{plolo_rate*0.04487 - obs:+.5f} (PLOLO) raw `pds_cosine`.")


if __name__ == "__main__":
    main()

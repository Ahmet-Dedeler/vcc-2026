"""232 / INBOX-1 item 2 — what does the WHITENING OPERATOR cost in |t|-weighted SIGN accuracy?

Phase 1b's survival rows showed `acc pre 0.57163 -> whitened 0.55468` on `orion_hct116`.  The
operator has only ever been optimised against DIRECTION metrics (170 the vertex, 193 the ridge,
212 the exponent, 226 the fitting sample); none of them scored it on SIGN accuracy, and §8.1d
says `nmae` pays for signs.  So the trade is measured here, and ONLY the trade.

`VCC_WHITEN_GAMMA` is NOT re-opened as a knob -- 170 measured the vertex end-to-end (every
gamma in [0.35, 0.80] net negative on all six) and a19 lost -0.0135 on the board.  What is
measured is the fact that production's gamma was chosen without anyone knowing its `q` cost.

Three quantities per gamma, three cancer folds, everything else at the shipped stack:
  * `q` at four stages -- consensus `pre`, whitened, final (after CD4 + H1 union), delivered;
  * `cos_agg` and `pds`, which is what 170/193/212 optimised;
  * the implied `nmae` cost of the `q` change, by TWO routes: 188's `nmae = 1 + k(1 - 2q)`
    identity at our own fitted `k`, and 232's directly MEASURED exchange rate.

    uv run python gamma232.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common232 import (CANCER, POOL_N, WORK, FoldBundle, ctx_gate, de_pools,  # noqa: E402
                       deliver, to_lfc, umd, weighted_sign_accuracy)

GAMMAS = (0.0, 0.35, 0.5, 0.65, 0.8)          # 0.65 is production
PROD_GAMMA = 0.65
# 232 §5, measured on 14 oracle/antioracle doses at bit-identical magnitude
D_NMAE_PER_QW = -0.2573
COMPOSITE_PER_QW = 0.3990                      # overall points per unit q_w, k_cos 1.384


def main() -> None:
    t0 = time.time()
    bctx, gate = {}, {}
    for c in ("A", "B", "C"):
        bctx[c], gate[c] = ctx_gate(c)
    out: dict = {}
    for held in CANCER:
        b = FoldBundle(held)
        pools = de_pools(b)
        pool = pools["wald89"]
        sc, truth = b.sc, b.truth
        LR = to_lfc(b.truth_h, b.fold["b_ctrl"], gate["A"])[:, sc]
        pre_sign = np.sign(b.pre[:, sc]).astype(np.int8)
        q_pre = weighted_sign_accuracy(pre_sign, truth, pool)[0]
        for g in GAMMAS:
            wh = b.whitened(gamma=g)
            q_wh = weighted_sign_accuracy(np.sign(wh[:, sc]).astype(np.int8), truth, pool)[0]
            del wh
            fin = b.final_gamma(g)
            q_fin = weighted_sign_accuracy(np.sign(fin[:, sc]).astype(np.int8), truth, pool)[0]
            e = b.score(fin)
            H = deliver(fin, bctx["A"], gate["A"], None)
            q_del = weighted_sign_accuracy(np.sign(H[:, sc]).astype(np.int8), truth, pool)[0]
            LP = to_lfc(H, bctx["A"], gate["A"])[:, sc]
            U, M, D = umd(LP, LR, pool)
            nmae = float(np.mean(1 + (U - 2 * M) / D))
            del fin, H, LP
            out.setdefault(held, {})[f"{g:g}"] = dict(
                q_pre=q_pre, q_whitened=q_wh, q_final=q_fin, q_delivered=q_del,
                cos_agg=e["cos_agg"], pds=e["pds"], nmae=nmae,
                U_over_D=float(np.mean(U / D)), M_over_D=float(np.mean(M / D)))
            print(f"[{held}/g={g:g}] q pre {q_pre:.5f} wh {q_wh:.5f} fin {q_fin:.5f} "
                  f"del {q_del:.5f} | cos_agg {e['cos_agg']:.6f} pds {e['pds']:.6f} "
                  f"nmae {nmae:.5f} ({time.time()-t0:.0f}s)", flush=True)
            (WORK / "gamma232.json").write_text(json.dumps(out, indent=1))
        del b

    def mean(g, k):
        return float(np.mean([out[f][f"{g:g}"][k] for f in CANCER]))

    print("\nTHE TRADE — three cancer folds, shipped stack, only `gamma` varies. `wald89`.")
    hdr = (f"{'gamma':>7}{'q delivered':>13}{'dq vs g0':>11}{'q whitened':>12}"
           f"{'cos_agg':>10}{'pds':>10}{'nmae':>10}{'d nmae vs g0':>14}")
    print(hdr)
    print("-" * len(hdr))
    q0 = mean(0.0, "q_delivered")
    n0 = mean(0.0, "nmae")
    rows = {}
    for g in GAMMAS:
        qd, nm = mean(g, "q_delivered"), mean(g, "nmae")
        rows[g] = dict(q_delivered=qd, dq_vs_g0=qd - q0, q_whitened=mean(g, "q_whitened"),
                       q_pre=mean(g, "q_pre"), q_final=mean(g, "q_final"),
                       cos_agg=mean(g, "cos_agg"), pds=mean(g, "pds"), nmae=nm,
                       d_nmae_vs_g0=nm - n0)
        star = "  <- PRODUCTION" if g == PROD_GAMMA else ""
        print(f"{g:>7g}{qd:>13.5f}{qd - q0:>+11.5f}{rows[g]['q_whitened']:>12.5f}"
              f"{rows[g]['cos_agg']:>10.6f}{rows[g]['pds']:>10.6f}{nm:>10.5f}"
              f"{nm - n0:>+14.5f}{star}")

    # ---- the operator's own sign cost, isolated: consensus -> whitened, at production gamma
    qp = mean(PROD_GAMMA, "q_pre")
    qw = mean(PROD_GAMMA, "q_whitened")
    dq_op = qw - qp
    # k from 188's identity, fitted on THIS experiment's own oracle response: dnmae/dq = -2k
    k_fit = -D_NMAE_PER_QW / 2.0
    print(f"\nOPERATOR SIGN COST, isolated (consensus -> whitened, gamma = {PROD_GAMMA}):")
    print(f"  q  {qp:.5f} -> {qw:.5f}   d q = {dq_op:+.5f}  ({100*dq_op/qp:+.2f}% relative)")
    # NOTE: 188's identity gives `dnmae/dq = -2k`, so quoting the identity at a `k` FITTED
    # from this experiment's own oracle response is the SAME route, not a second one.  What
    # is independent is that the fitted k matches 187's separately measured |t|-weighted
    # median per-gene k of 0.12.
    print(f"  232's measured dnmae/dq = {D_NMAE_PER_QW:.4f} implies k = {k_fit:.4f} under "
          f"188's identity `nmae = 1 + k(1-2q)`  (187 measured median per-gene k = 0.12)")
    print(f"  implied d nmae from the q change alone: {D_NMAE_PER_QW*dq_op:+.5f}   "
          f"(directly measured d nmae, gamma {PROD_GAMMA} vs 0: {rows[PROD_GAMMA]['d_nmae_vs_g0']:+.5f}"
          f" -- the q route explains "
          f"{100*D_NMAE_PER_QW*dq_op/rows[PROD_GAMMA]['d_nmae_vs_g0']:.0f}% of it)")
    print(f"  implied overall points at 232's measured rate: {COMPOSITE_PER_QW*dq_op:+.5f}")
    print(f"  ... and the operator BUYS (production vs gamma=0): "
          f"cos_agg {mean(0.0,'cos_agg'):.6f} -> {mean(PROD_GAMMA,'cos_agg'):.6f} "
          f"({100*(mean(PROD_GAMMA,'cos_agg')/mean(0.0,'cos_agg')-1):+.1f}%), "
          f"pds {mean(0.0,'pds'):.5f} -> {mean(PROD_GAMMA,'pds'):.5f}")

    out["summary"] = dict(by_gamma=rows, prod_gamma=PROD_GAMMA,
                          operator_sign_cost=dict(
                              q_pre=qp, q_whitened=qw, dq=dq_op,
                              k_fitted_from_oracle=k_fit,
                              implied_d_nmae_identity=-2 * k_fit * dq_op,
                              implied_d_nmae_measured=D_NMAE_PER_QW * dq_op,
                              implied_overall=COMPOSITE_PER_QW * dq_op))
    (WORK / "gamma232.json").write_text(json.dumps(out, indent=1))
    print(f"\n-> {WORK / 'gamma232.json'}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()

"""232 Phase 1 — is the donors' MAJORITY sign actually better than the sign of their MEAN?

Measured before anything is built, exactly as the mandate requires.  Three questions, in
the order that decides whether Phase 2 happens at all:

  1. HOW OFTEN CAN THE TWO RULES EVEN DIFFER?  Where every covering donor agrees in sign,
     the mean, the majority, every weighted vote and the oracle all return the SAME sign.
     Only the DISAGREEMENT set is contestable.  So the |t|-weighted disagreement share is
     an upper bound on the mass any sign-combination rule can move, and

         prize_ceiling = mass_disagree * (1 - acc[sign(mean)] on the disagreement set)

     is an arithmetic ceiling on the |t|-weighted sign accuracy any rule can reach --
     including the oracle.  It is computed here before any arm exists.

  2. ON THAT SET, IS ANY IMPLEMENTABLE RULE BETTER THAN `sign(mean)`?
     `sign(mean)` breaks a disagreement by MAGNITUDE: the donor with the larger |effect|
     wins.  The alternatives break it by count (majority), by standardised response
     (`z`), or by donor identity (the per-donor-scalar family, which is what an arity-2
     vote degenerates to).  A per-GENE weight -- 189's `theta`, or a symmetric cross-donor
     reliability -- cancels exactly and is not a rule at all; see `common232.rho_weights`.

  3. DOES A SIGN FIXED UPSTREAM SURVIVE TO THE SCORER?  The consensus sign is two
     transformations away from what is scored: the whitening operator is a dense linear
     map on `wcols`, then two residual sources are added, then the renderer.  Phase 1b
     measures the |t|-weighted agreement between `sign(pre)`, `sign(final)` and
     `sign(delivered)` on the same pool.  If the consensus sign does not survive, an
     upstream vote cannot pay however good it is.

ARITY, STATED PLAINLY.  Production mixes THREE donor lines.  PLOLO holds one out, so the
three cancer folds have TWO covering donors and the unweighted majority is a TIE on every
contested pair -- it is not a rule there, and that is reported as a structural fact rather
than hidden.  Two things restore arity 3: the `kolf2` fold (all three production donors,
leak-free truth, reported apart), and a `+kolf2` voter set on the cancer folds (a
NON-production source used only where it is not the truth).

    uv run python sign232.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common232 import (ALL_FOLDS, WORK, FoldBundle, ctx_gate, de_pools,  # noqa: E402
                       deliver, load_sources, pool_weight_share, stack_sources,
                       weighted_sign_accuracy, z_weights, rho_weights)

RNG_SEED = 20260904


def sign_structure(E, avail, supp, sc):
    """Per (target, gene) on the scored axis: coverage count and the sign tallies."""
    Es = E[:, :, sc]
    cov = avail[:, :, None] & supp[:, None, sc]
    s = np.sign(Es) * cov
    pos = (s > 0).sum(0).astype(np.int16)
    neg = (s < 0).sum(0).astype(np.int16)
    ncov = cov.sum(0).astype(np.int16)
    return Es, cov, s, pos, neg, ncov


def rules(Es, cov, s, pos, neg, ncov, avail, supp, sc, s_mean):
    """Every sign rule, on the scored axis.  Ties resolved to `sign(mean)` (leave alone)."""
    def tb(x):
        return np.where(x != 0, x, s_mean).astype(np.int8)

    out = {"mean": s_mean.astype(np.int8)}
    out["maj_strict"] = np.sign(pos - neg).astype(np.int8)
    out["maj"] = tb(np.sign(pos - neg))
    zw = z_weights(np.ascontiguousarray(Es), avail, supp[:, sc])
    out["zvote"] = tb(np.sign((s * zw).sum(0)))
    del zw
    for k in range(Es.shape[0]):
        out[f"donor{k}"] = tb(s[k])
    return out


def analyse(bundle: FoldBundle, src: dict, use_extra: bool, tag: str) -> list[dict]:
    held = bundle.held
    names, E, avail, supp = stack_sources(src, use_extra=use_extra)
    sc = bundle.sc
    truth = bundle.truth
    pools = de_pools(bundle)

    Es, cov, s, pos, neg, ncov = sign_structure(E, avail, supp, sc)
    s_mean = np.sign(Es.sum(0)).astype(np.int8)
    # A target NO source carries takes 062's fallback row, whose sign is a property of the
    # donors' mean RESPONSE profile and not of any vote.  Those rows have `ncov == 0`, so
    # they are already outside every arm's reach; the identity is asserted on the rest.
    covered_t = src["votes"] > 0
    if not use_extra:
        ref = np.sign(np.asarray(bundle.fold["pre"], dtype=np.float64)[:, sc]).astype(np.int8)
        agree = float((ref[covered_t] == s_mean[covered_t]).mean())
        assert agree > 0.9999, (f"{held}: sign(sum over donors) != sign(pre) on "
                                f"{100*(1-agree):.4f}% of covered pairs")
        print(f"[{held}] sign(sum over donors) == sign(pre) on {100*agree:.4f}% of "
              f"covered-target pairs; {int((~covered_t).sum())} fallback target(s) excluded",
              flush=True)

    disagree = (pos > 0) & (neg > 0)
    unan = (ncov >= 2) & ~disagree & ((pos + neg) >= 2)
    single = ncov <= 1

    R = rules(Es, cov, s, pos, neg, ncov, avail, supp, sc, s_mean)
    R["oracle"] = np.sign(truth).astype(np.int8)

    rows = []
    for pool_name, pool in pools.items():
        base = dict(fold=held, voters=tag, n_voters=int(E.shape[0]),
                    names=names, pool=pool_name)
        # ---- coverage honesty and the contestable mass
        cov_share = {f"ncov{c}": pool_weight_share(truth, pool, ncov == c)
                     for c in range(0, int(ncov.max()) + 1)}
        cnt = {f"ncov{c}": float(np.mean([(ncov[i][pool[i]] == c).mean()
                                          for i in range(truth.shape[0])]))
               for c in range(0, int(ncov.max()) + 1)}
        mass_dis = pool_weight_share(truth, pool, disagree)
        mass_tie = pool_weight_share(truth, pool, disagree & (pos == neg))
        mass_unan = pool_weight_share(truth, pool, unan)
        mass_single = pool_weight_share(truth, pool, single)
        cnt_dis = float(np.mean([disagree[i][pool[i]].mean()
                                 for i in range(truth.shape[0])]))

        accs, accs_dis, accs_unan = {}, {}, {}
        for rn, rs in R.items():
            accs[rn] = weighted_sign_accuracy(rs, truth, pool)[0]
            accs_dis[rn] = weighted_sign_accuracy(rs, truth, pool, disagree)[0]
            accs_unan[rn] = weighted_sign_accuracy(rs, truth, pool, unan)[0]

        ceiling = mass_dis * (1.0 - accs_dis["mean"]) if np.isfinite(accs_dis["mean"]) else 0.0
        rows.append(dict(**base, acc=accs, acc_disagree=accs_dis, acc_unanimous=accs_unan,
                         mass_disagree=mass_dis, mass_tie=mass_tie, mass_unanimous=mass_unan,
                         mass_singleton=mass_single, count_disagree=cnt_dis,
                         mass_by_ncov=cov_share, count_by_ncov=cnt,
                         prize_ceiling_on_sign_accuracy=ceiling))
        print(f"[{held}/{tag}/{pool_name}] |t|-w disagree mass {mass_dis:.4f} "
              f"(count {cnt_dis:.4f})  acc mean {accs['mean']:.5f} maj {accs['maj']:.5f} "
              f"z {accs['zvote']:.5f} | on disagree: mean {accs_dis['mean']:.5f} "
              f"maj {accs_dis['maj']:.5f} z {accs_dis['zvote']:.5f} | ceiling "
              f"{ceiling:+.5f}", flush=True)
    return rows


def survival(bundle: FoldBundle, src: dict) -> dict:
    """Phase 1b — does the consensus sign survive whitening, the residuals, the renderer?"""
    sc, truth = bundle.sc, bundle.truth
    pools = de_pools(bundle)
    pre_sign = np.sign(np.asarray(bundle.fold["pre"], dtype=np.float64)[:, sc]).astype(np.int8)
    wh = bundle.whitened()
    wh_sign = np.sign(wh[:, sc]).astype(np.int8)
    del wh
    fin = bundle.final()
    fin_sign = np.sign(fin[:, sc]).astype(np.int8)
    bA, gA = ctx_gate("A")
    H = deliver(fin, bA, gA, "minnorm")
    del_sign = np.sign(H[:, sc]).astype(np.int8)
    del fin, H
    out = {}
    for pn, pool in pools.items():
        def share(a, b):
            vals = []
            for i in range(truth.shape[0]):
                w = np.abs(truth[i, pool[i]])
                t = float(w.sum())
                if t > 0:
                    vals.append(float(w[a[i, pool[i]] == b[i, pool[i]]].sum() / t))
            return float(np.mean(vals))
        out[pn] = dict(
            pre_vs_whitened=share(pre_sign, wh_sign),
            pre_vs_final=share(pre_sign, fin_sign),
            pre_vs_delivered=share(pre_sign, del_sign),
            final_vs_delivered=share(fin_sign, del_sign),
            acc_pre=weighted_sign_accuracy(pre_sign, truth, pool)[0],
            acc_whitened=weighted_sign_accuracy(wh_sign, truth, pool)[0],
            acc_final=weighted_sign_accuracy(fin_sign, truth, pool)[0],
            acc_delivered=weighted_sign_accuracy(del_sign, truth, pool)[0])
        print(f"[{bundle.held}/survival/{pn}] pre->whitened {out[pn]['pre_vs_whitened']:.4f} "
              f"pre->final {out[pn]['pre_vs_final']:.4f} "
              f"pre->delivered {out[pn]['pre_vs_delivered']:.4f} | acc pre "
              f"{out[pn]['acc_pre']:.5f} whitened {out[pn]['acc_whitened']:.5f} "
              f"final {out[pn]['acc_final']:.5f} delivered {out[pn]['acc_delivered']:.5f}",
              flush=True)
    return out


def main() -> None:
    t0 = time.time()
    out_path = WORK / "sign232.json"
    prev = json.loads(out_path.read_text()) if out_path.exists() else {}
    rows, surv = prev.get("rows", []), prev.get("survival", {})
    done = {(r["fold"], r["voters"], r["pool"]) for r in rows}
    folds = [a for a in sys.argv[1:] if not a.startswith("-")] or list(ALL_FOLDS)
    for held in folds:
        if (held, "production", "wald89") in done:
            print(f"[{held}] done", flush=True)
            continue
        b = FoldBundle(held)
        src = load_sources(held)
        rows += analyse(b, src, use_extra=False, tag="production")
        if "extra_E" in src:
            rows += analyse(b, src, use_extra=True, tag="production+kolf2")
        surv[held] = survival(b, src)
        out_path.write_text(json.dumps(
            dict(rows=rows, survival=surv, seconds=time.time() - t0), indent=1))
        del b, src
    print(f"\n-> {out_path}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()

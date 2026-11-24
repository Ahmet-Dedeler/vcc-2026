"""227 — reconciling Aginglab's row: how precise is 197's cosine inversion, really?

The Phase-1 frontier says reaching their `nmae` credit costs at least +5.0% relative aggregate
cosine (h1-test) / +1.8% per +0.0065 raw (PLOLO).  Their board `mse` says their cosine is 0.6%
BELOW ours.  Those cannot both be true, so one of them is being read through an instrument that
is not precise enough.

197's inversion `cos = (1.023611 - raw_mse) / 2.074495` is calibrated on OUR arms.  The scored
numerator (AGENTS.md 4.1b) is

    N_p = (1/G_p) ( ||b_hat - b||^2  -  rho * min(C_hat_p, C_p)  -  C_p )

and `rho * min(C_hat, C)` -- the jackknife credit for our OWN emitted noise, `claim` in 187's
notation -- is a property of the EMITTER, not of the direction.  A submission whose emission is
tighter than ours claims less of it and prints a higher raw `mse` at the same cosine.  This
script prices that sensitivity from 187's own cached constants and asks how large an emission
difference would reconcile the two readings.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E187 = os.path.join(REPO, "experiments/187-joint-amplitude/work/mse187.json")

SIGMA_D = 0.229141          # 187 sec 1, C_p-free LS solve on a8+a10+a11
COS_B = 2.074495            # 197's inversion slope
US_MSE, THEM_MSE = 0.81726, 0.81851
US_COS = (1.023611 - US_MSE) / COS_B
US_NMAE, THEM_NMAE = 0.93389, 0.88902


def main() -> None:
    J = json.load(open(E187))
    claim = {c: J["a23"][c]["claim180"] for c in ("A", "B", "C")}
    claim_mean = float(np.mean(list(claim.values())))

    # d raw_mse / d claim  = -1/SigmaD ; then read through 197's inversion as a cosine
    d_mse_per_1pct = claim_mean * 0.01 / SIGMA_D
    d_cos_per_1pct = d_mse_per_1pct / COS_B
    pct_cos_per_1pct_claim = 100.0 * d_cos_per_1pct / US_COS

    # what the two readings each require
    obs_d_mse = THEM_MSE - US_MSE
    obs_d_cos_pct = 100.0 * (-obs_d_mse / COS_B) / US_COS

    # PLOLO exchange rate, 3 cancer folds: +0.006497 raw nmae per +1.775% relative cos_agg
    P = json.load(open(os.path.join(WORK, "plolo227.json")))
    a = P["summary"]["3 cancer folds"]["arms"]["fix lo-|h| wrong signs, f=0.2"]
    rate = a["d_nmae"] / a["d_cos_pct"]                 # raw nmae per 1% relative cosine
    need_cos_pct = (THEM_NMAE - US_NMAE) / -rate        # their gap is negative (better)
    need_claim_pct = need_cos_pct / pct_cos_per_1pct_claim

    # h1-test's own version of the same exchange rate
    H = json.load(open(os.path.join(WORK, "phase3.json")))
    h = H["q_w_required"]["fix-lo-|h|"]
    h_cos_pct = 100.0 * (h["cos_at_target"] / H["q_w_required"]["cos_ours"] - 1.0)

    out = dict(
        claim_per_context=claim, claim_mean=claim_mean, sigma_d=SIGMA_D,
        d_raw_mse_per_1pct_claim=d_mse_per_1pct,
        cos_pct_misread_per_1pct_claim=pct_cos_per_1pct_claim,
        observed_d_raw_mse=obs_d_mse, observed_d_cos_pct=obs_d_cos_pct,
        plolo_rate_nmae_per_1pct_cos=rate,
        cos_pct_needed_for_their_nmae=need_cos_pct,
        claim_pct_that_reconciles=need_claim_pct,
        h1_cos_pct_at_their_credit=h_cos_pct,
        implied_mse_if_direction=US_MSE - need_cos_pct / 100.0 * US_COS * COS_B,
    )
    json.dump(out, open(os.path.join(WORK, "recon227.json"), "w"), indent=1)

    print("== 197's cosine inversion: how much of it is emitter, not direction? ==")
    print(f"  187's a23 jackknife claim per context: "
          + "  ".join(f"{c} {v:.6f}" for c, v in claim.items())
          + f"   mean {claim_mean:.6f}")
    print(f"  SigmaD {SIGMA_D:.6f}  ->  a 1% relative change in the claim moves raw `mse` by "
          f"{d_mse_per_1pct:.6f},")
    print(f"  which 197's inversion reads as {pct_cos_per_1pct_claim:.2f}% of our transfer "
          f"cosine ({US_COS:.5f}).")
    print(f"\n  Aginglab's raw `mse` is {obs_d_mse:+.5f} vs ours -> read as "
          f"{obs_d_cos_pct:+.2f}% cosine.")
    print(f"  That is {abs(obs_d_cos_pct)/pct_cos_per_1pct_claim:.2f}% of a jackknife-claim "
          f"difference.  The board cannot resolve their cosine below a few percent.")
    print(f"\n== what their `nmae` would require, if it came from direction ==")
    print(f"  PLOLO (3 cancer folds), the cheapest sign-fix ordering there is: "
          f"{rate:+.6f} raw nmae per 1% relative cos_agg")
    print(f"  their nmae edge {THEM_NMAE - US_NMAE:+.5f} raw  ->  needs "
          f"{need_cos_pct:+.1f}% relative cosine")
    print(f"  h1-test cross-check at their credit ratio: {h_cos_pct:+.1f}% relative cos_agg")
    print(f"  that cosine would print raw `mse` {out['implied_mse_if_direction']:.5f}; "
          f"they print {THEM_MSE:.5f}")
    print(f"\n  RECONCILIATION: a {need_claim_pct:.1f}% smaller jackknife claim than ours -- "
          f"i.e. a tighter emission --\n  makes those two readings the SAME submission.  Their "
          f"`jac` {-0.03295:+.5f} and `fid` {+0.01970:+.5f} against ours independently say their "
          f"call set\n  is wider and looser than ours, which is what a tighter emission "
          f"produces.")


if __name__ == "__main__":
    main()

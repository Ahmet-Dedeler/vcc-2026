"""215 — the CHANNEL DECOMPOSITION instrument.

Four channels, their board signatures, and the arithmetic that separates them.

    DIRECTION   pds up, mse down, nmae down together        (a14 -> a15, the CD4 residual)
    ORDERING    pds, mse, nmae EXACTLY unchanged; DE moves   (a15 -> rc_hybrid -> rc_callset)
    AMPLITUDE   mse moves, pds pinned by scale invariance    (a10 -> a11, alpha 1.0796 -> 0.838)
    GEOMETRY    mse moves with NO pds corroboration and nmae moving the OTHER way  (a15 -> a19)

The load-bearing asymmetry is that `pds` is the only member that is
  * EXACTLY invariant to the emitted order and call set (board-proven twice, 16 s.f.), and
  * invariant to a global amplitude to within a measured +-0.0015 raw,
so `pds` is the amplitude-free, ordering-free readout of a real change in the prediction, and
`mse` is the amplitude-contaminated one.  Their RATIO is the discriminant.
"""
import numpy as np

# ---------------------------------------------------------------- the scaling map
# fitted in panel215.py on 2,208 board observations, agrees with docs/BOARD-LEDGER.md to <=0.11%
SLOPE = {"pds": 2.21433, "mse": -1.05211, "nmae": -1.65443,
         "fid": 3.36305, "reach": 1.12921, "jac": 2.72607}
ZERO = {"pds": 0.50026, "mse": 0.98792, "nmae": 1.00034,
        "fid": 0.51241, "reach": 0.07902, "jac": 0.03037}
# the board clamps two members (per context, then averages -> convex near the kink)
FLOOR = {"mse": 0.0, "nmae": -6.0}
MEMBERS = ["pds", "mse", "nmae", "fid", "reach", "jac"]


def scale(m, raw):
    v = SLOPE[m] * (raw - ZERO[m])
    if m in FLOOR:
        v = max(FLOOR[m], v)
    return v


# ---------------------------------------------------------------- 197's inversion
COS_A, COS_B = 1.023611, 2.074495
def cos_floor(raw_mse):
    return (COS_A - raw_mse) / COS_B
def is_informative(raw_mse):
    return raw_mse < COS_A


# ---------------------------------------------------------------- channel response vectors
# DIRECTION, per +1.000 of transfer cosine.  Measured on a14 -> a15 (the CD4 residual), the only
# single-variable direction arm this project has put on the board.  197 measured its cosine move
# at +0.00507 (0.07418 -> 0.07925) by an inversion with eight passing known-answer gates.
DCOS_A15 = 0.00507
E_DIR = {          # raw member delta per unit cosine
    "pds":   (0.835987 - 0.832702) / DCOS_A15,     # +0.6479
    "mse":   (0.873805 - 0.882775) / DCOS_A15,     # -1.7692
    "nmae":  (0.944251 - 0.945764) / DCOS_A15,     # -0.2984
    "fid":   (0.503536 - 0.502881) / DCOS_A15,     # +0.1292
    "reach": (0.181173 - 0.180615) / DCOS_A15,     # +0.1101
    "jac":   (0.030080 - 0.030081) / DCOS_A15,     # -0.0002
}
# 197's identity says d(raw mse)/d(cos) = -2.074495 exactly at the rho cliff with claim at the
# cap.  The measured -1.7692 is 15% shallower because a15's own `claim` and `A` moved too.  Both
# are carried; the identity value is the one used for the field, where only `raw` is observable.
E_DIR_MSE_IDENTITY = -COS_B

# the direction channel's member RATIOS to pds -- these are what the decomposition actually uses,
# because pds is the amplitude-free member.
RATIO_DIR = {m: E_DIR[m] / E_DIR["pds"] for m in MEMBERS}
RATIO_DIR_MSE_IDENTITY = E_DIR_MSE_IDENTITY / E_DIR["pds"]     # -3.2019

# ---------------------------------------------------------------- noise floors (raw units)
# sigma_pds: our one pure-amplitude board arm (a10 -> a11, alpha -23%) moved raw pds by +0.001503
#            while 197 measured its cosine flat to -0.00004.  That is the quantisation/seed floor.
# sigma_mse/nmae: rc_hybrid and rc_callset reproduced a15's to SIXTEEN significant figures at a
#            pinned profile, so the floor here is set by re-seeded twins, not by the scorer.
# sigma_reach/fid/jac: 170's re-seeded twin at fixed config moves raw reach by up to 0.0096.
SIGMA = {"pds": 0.0004, "mse": 0.0020, "nmae": 0.0005,
         "fid": 0.0015, "reach": 0.0050, "jac": 0.0005}
K_SIG = 3.0                      # a member "moved" if |delta| > K_SIG * sigma

# `pds` is scale-invariant in exact arithmetic but NOT in integer counts: a10 -> a11 cut alpha by
# 23% and moved raw pds by +0.001503 against a cosine 197 measured FLAT (-0.00004).  That leakage
# is proportional to the size of the amplitude move, for which |delta raw mse| is the observable.
# So the pds significance threshold is adaptive:  sigma_eff = max(sigma_seed, LEAK * |dmse|).
PDS_LEAK = 0.001503 / 0.442426          # = 0.003397, measured on a10 -> a11


def sigma_pds_eff(d_mse):
    return max(SIGMA["pds"], PDS_LEAK * abs(d_mse))


def classify(before, after, ratio_mse=None):
    """Decompose one arm-to-arm move.  `before`/`after` are dicts of raw members.

    Returns a dict with the four channel magnitudes, the overall-units budget, and a label.
    The budget is an EXACT partition of the observed delta-overall.
    """
    if ratio_mse is None:
        ratio_mse = RATIO_DIR_MSE_IDENTITY
    d = {m: after[m] - before[m] for m in MEMBERS}
    s_pds = sigma_pds_eff(d["mse"])
    sig = {m: abs(d[m]) > K_SIG * SIGMA[m] for m in MEMBERS}
    sig["pds"] = abs(d["pds"]) > K_SIG * s_pds

    # --- the two direction readouts
    dcos_pds = d["pds"] / E_DIR["pds"]                 # amplitude-free, ordering-free
    dcos_mse = -d["mse"] / COS_B                       # = delta of 197's floor; amp-contaminated
    inf_b, inf_a = is_informative(before["mse"]), is_informative(after["mse"])
    both_informative = inf_b and inf_a

    # --- corroboration ratio: +1 pure direction, 0 pure amplitude, negative geometry
    rho = dcos_pds / dcos_mse if abs(dcos_mse) > 1e-9 else float("nan")

    # --- counterfactual: what the six members would read if ONLY direction had moved, by dcos_pds
    cf = {"pds": after["pds"]}
    cf["mse"] = before["mse"] + ratio_mse * d["pds"]
    for m in ["nmae", "fid", "reach", "jac"]:
        cf[m] = before[m] + RATIO_DIR[m] * d["pds"]

    sca_obs = {m: scale(m, after[m]) - scale(m, before[m]) for m in MEMBERS}
    sca_dir = {m: scale(m, cf[m]) - scale(m, before[m]) for m in MEMBERS}
    sca_res = {m: sca_obs[m] - sca_dir[m] for m in MEMBERS}

    ovr = lambda dd, ms: sum(dd[m] for m in ms) / 6.0
    b_direction = ovr(sca_dir, MEMBERS)
    b_ordering = ovr(sca_res, ["fid", "reach", "jac"])
    # split the (mse, nmae) residual between AMPLITUDE and GEOMETRY.
    #   amplitude   moves mse hard and nmae essentially not at all  (a10->a11: |dnmae/dmse| 0.0036)
    #   geometry    moves nmae the OPPOSITE way to mse              (a15->a19: dnmae/dmse = -0.056)
    res_mse, res_nmae = d["mse"] - (cf["mse"] - before["mse"]), d["nmae"] - (cf["nmae"] - before["nmae"])
    if abs(res_mse) > K_SIG * SIGMA["mse"] and res_mse * res_nmae < 0 \
            and abs(res_nmae) > K_SIG * SIGMA["nmae"]:
        geom_share = min(1.0, abs(res_nmae / res_mse) / 0.056)
    else:
        geom_share = 0.0
    pool = ovr(sca_res, ["mse", "nmae"])
    b_geometry = geom_share * pool
    b_amplitude = (1.0 - geom_share) * pool

    # --- the label
    # `mse` only enters the OVERALL score while it is unclamped.  When both endpoints sit above
    # the clamp point the mse move is real but UNPRICED, and neither AMPLITUDE nor GEOMETRY can
    # be told apart there -- both are invisible to the board.  In that regime the only member
    # that carries information is `pds`, so the label comes from `pds` and the DE members alone.
    priced = (scale("mse", before["mse"]) > 0.0) or (scale("mse", after["mse"]) > 0.0)
    moved_de = sig["fid"] or sig["reach"] or sig["jac"]
    # ORDERING is the STRICT signature of our own rc_hybrid / rc_callset arms: all three
    # direction members pinned at the SEED floor (not the amplitude-inflated one), DE moving.
    pinned = (abs(d["pds"]) <= K_SIG * SIGMA["pds"]
              and abs(d["mse"]) <= K_SIG * SIGMA["mse"]
              and abs(d["nmae"]) <= K_SIG * SIGMA["nmae"])
    if pinned and moved_de:
        label = "ORDERING"
    elif pinned:
        label = "NULL"
    elif not priced:
        # the mse move is invisible to `overall`, so AMPLITUDE and GEOMETRY cannot be told
        # apart and the pds move may be amplitude leakage; only a resolved pds move is readable
        label = "DIRECTION(pds)" if sig["pds"] else "UNRESOLVED"
    elif not sig["pds"] and sig["mse"]:
        label = "GEOMETRY" if geom_share > 0.5 else "AMPLITUDE"
    elif sig["pds"] and not sig["mse"]:
        label = "DIRECTION"
    else:
        if rho > 0.15:
            label = "DIRECTION" if rho > 0.60 else "DIRECTION+AMPLITUDE"
        elif rho < -0.10:
            label = "GEOMETRY"
        else:
            label = "GEOMETRY" if geom_share > 0.5 else "AMPLITUDE"

    return {
        "d_raw": d, "sig": sig,
        "dcos_pds": dcos_pds, "dcos_mse": dcos_mse, "rho": rho,
        "informative_before": inf_b, "informative_after": inf_a,
        "both_informative": both_informative,
        "cos_floor_before": cos_floor(before["mse"]), "cos_floor_after": cos_floor(after["mse"]),
        "sca_obs": sca_obs, "sca_dir": sca_dir, "sca_res": sca_res,
        "ovr_obs": ovr(sca_obs, MEMBERS),
        "b_direction": b_direction, "b_ordering": b_ordering,
        "b_amplitude": b_amplitude, "b_geometry": b_geometry,
        "geom_share": geom_share, "label": label, "mse_priced": bool(priced),
        # uncertainty on the direction magnitude, from the adaptive pds floor
        "sigma_pds_eff": s_pds,
        "dcos_pds_sd": s_pds / E_DIR["pds"],
    }

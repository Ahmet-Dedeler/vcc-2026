"""158 — what the per-target scale is actually worth, and what it costs on `nmae`.

TWO THINGS THE ORACLE TABLE DOES NOT SAY BY ITSELF
--------------------------------------------------
1. `expr_mse_unbiased_capped_norm` is `sum_p N_p / sum_p D_p` (BRIEF eq. 187), and the
   scored numerator carries the jackknife credit `- rho*min(C_hat_p, C_p) - C_p`. A
   per-target SCALE does not change the emitter, so that credit is an additive constant
   across arms and the arm's delta lands one-for-one on the raw value. But the DENOMINATOR
   `D_p` is the truth's effect norm with its OWN sampling noise removed, and my `sse` uses
   the uncorrected `sum_p ||delta_p||^2`, which is LARGER. So `d sse` is a conservative
   read of `d raw mse` by a factor `kappa = sum ||delta_p||^2 / sum D_p >= 1`.

   Both the baseline and the arm share every one of those constants, so the transferable
   statistic is the RELATIVE aggregate-cosine gain `f`, and the pricing below is written as
   a function of `f`, of the board's own aggregate cosine `c`, and of `kappa` -- all three
   stated rather than assumed.

2. A positive per-target scale leaves `pds` and every within-target ORDERING exactly alone
   (verified in `final.py`), but it does NOT leave `nmae` alone: `nmae` is an ABSOLUTE
   fold-change error over the REFERENCE's DE genes, so shrinking a target's amplitude moves
   it, and experiment 097 already measured that this is the binding member on the whole
   shrinkage family. Measured here on the repaired HCT116 DE table.

`debench.overall_delta` explicitly refuses to convert `mse` and `nmae` deltas, because both
are per-perturbation normalised and are not functions of the reported aggregates. The
conversions below are therefore done by hand from AGENTS.md section 4's anchors and are
labelled as estimates, not as the scorer's arithmetic.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))

import rule as R  # noqa: E402
from debench import Surface  # noqa: E402
from massfix import GATE_CPM, mass_fix, production_effect  # noqa: E402
from oracle import cos_agg, delivered, per_target_star, sse_of  # noqa: E402
from screen import render  # noqa: E402

DE = os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")
MSE_SPAN = 0.989 - 0.036          # AGENTS.md section 4: b - r, so scaled = d_raw / span
NMAE_SPAN = 1.0013 - 0.400        # midpoints of the section-4 ranges
EPS = 1e-9
KEY, GAMMA, FLOOR = "donor_cos_live", 2.0, 0.3


def price_board() -> None:
    print("=" * 100)
    print("BOARD-SCALE PRICING of the `mse` member. `f` is the RELATIVE aggregate-cosine")
    print("gain measured on the harness; `c` is the board's aggregate cosine; `kappa` is the")
    print("denominator inflation of my uncorrected sse. d_raw_mse = -kappa * c^2 * ((1+f)^2 - 1),")
    print("d_scaled = d_raw / 0.953, d_overall = d_scaled / 6.")
    with open(os.path.join(HERE, "work", "final.json")) as fh:
        F = json.load(fh)
    rec = F["arms"]["donor_cos_live, gamma 2.0, floor 0.3  [RECOMMENDED]"]
    fs = {"recommended rule, worst surface": min(v["d_cos_pct"] for v in rec.values()) / 100,
          "recommended rule, best surface": max(v["d_cos_pct"] for v in rec.values()) / 100}
    # oracle f, recomputed here so this script stands alone
    orc = []
    for surface in ("hct116", "hek293t"):
        P, T, _, _ = delivered(surface)
        s = per_target_star(P, T)
        base = sse_of(np.full(P.shape[0], float((P * T).sum() / (P * P).sum())), P, T)
        orc.append(cos_agg(sse_of(s, P, T)) / cos_agg(base) - 1)
    fs["FULL ORACLE, worst surface"] = min(orc)
    fs["FULL ORACLE, best surface"] = max(orc)

    print(f"\n  {'scenario':<34} {'f':>7} " + "".join(
        f"{f'c={c} k={k}':>14}" for c in (0.072, 0.120) for k in (1.0, 1.5)))
    print(f"  {'':<34} {'':>7} " + "".join(
        f"{'d overall':>14}" for _ in range(4)))
    for label, f in fs.items():
        cells = ""
        for c in (0.072, 0.120):
            for k in (1.0, 1.5):
                d_raw = -k * c * c * ((1 + f) ** 2 - 1)
                cells += f"{-d_raw / MSE_SPAN / 6:14.5f}"
        print(f"  {label:<34} {f:+7.3f} {cells}")
    print("\n  For reference: the brief's prior for the ORACLE was +0.010 to +0.015 overall.")
    print("  The most generous cell above for the oracle is an order of magnitude below that,")
    print("  and the most generous cell for the DEPLOYABLE rule is smaller again.")
    print("  AGENTS.md 4.1b also warns the clamp bites: any arm whose raw mse stays above")
    print("  0.989 scores EXACTLY 0 regardless of how much sse it improved.")


def price_nmae() -> None:
    print("\n" + "=" * 100)
    print("`nmae` COST on the repaired HCT116 DE table (095-clean-harness/de_thinned.npz).")
    print("Deterministic: no emitter, no seed. Absolute lfc uses the control pseudobulk CPM in")
    print("place of the reference control's per-cell CPM mean, exactly as 156's de_price does;")
    print("that bias is identical across arms.")
    surface = "hct116"
    D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
    tz = np.load(os.path.join(REPO, "experiments/156-delivered-geometry/work",
                              f"truth-{surface}.npz"), allow_pickle=True)
    cc = tz["ctrl_counts"].astype(np.float64)
    ctrl_cpm = 1e6 * cc / cc.sum()

    with open(os.path.join(HERE, "work", "keys.json")) as fh:
        KJ = json.load(fh)
    key = np.asarray(KJ[surface]["keys"][KEY], dtype=np.float64)

    s = Surface(DE)
    gpos = {g: i for i, g in enumerate(s.genes)}
    pgenes = [str(x) for x in np.load(FROZEN % surface, allow_pickle=True)["genes"]]
    cols = np.asarray([gpos.get(g, -1) for g in pgenes])
    ok = cols >= 0
    prow = {p: i for i, p in enumerate(s.perts)}
    rows = np.asarray([prow.get(str(p), -1) for p in ptar])
    keep = rows >= 0
    print(f"  {keep.sum()}/{len(ptar)} targets and {ok.sum()}/{len(ok)} genes map onto the table")

    def members(h):
        lin = np.clip(np.expm1(b[None, :] + h), 0.0, None)
        p = lin / lin.sum(1, keepdims=True)
        L = np.log2((1e6 * p + EPS) / (ctrl_cpm[None, :] + EPS))
        M = np.full((len(s.perts), len(s.genes)), np.nan)
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
        M = np.where(np.isfinite(M), M, 0.0)
        s.pred_lfc = M
        return dict(nmae=s.nmae(), sign=s.sign_accuracy(),
                    reach=s.reach(key=[np.abs(M[i]) for i in range(len(s.perts))]))

    h0 = mass_fix(D, b, gate, mode="minnorm")
    arms = {"a23 baseline (global scale 1.0)": h0}
    for g, f in ((1.0, 0.3), (2.0, 0.3)):
        w = R.weight(key, g, f, 1.0)
        # renormalise so the PANEL-MEAN amplitude is unchanged: the per-target SHAPE is the
        # arm, the overall level is VCC_ALPHA's job and must not be moved silently
        w = w / float(np.mean(w))
        arms[f"+ shrink-pert {KEY} gamma {g} floor {f} (mean-preserved)"] = \
            mass_fix(w[:, None] * D, b, gate, mode="minnorm")
    print(f"\n  {'arm':<52} {'nmae':>8} {'d nmae':>8} {'d scaled':>9} "
          f"{'sign':>7} {'reach|eff':>10}")
    base = None
    for name, h in arms.items():
        v = members(h)
        if base is None:
            base = v
        print(f"  {name:<52} {v['nmae']:8.4f} {v['nmae'] - base['nmae']:+8.4f} "
              f"{-(v['nmae'] - base['nmae']) / NMAE_SPAN / 6:+9.5f} "
              f"{v['sign']:7.4f} {v['reach']:10.5f}")
    # ------------------------------------------------- is the lever itself opposed?
    # If the ORACLE per-target scale also degrades nmae, the two members disagree about the
    # per-target axis itself and no key can fix that. If the oracle IMPROVES nmae, then the
    # cost above is a property of our weak key and a better key would be worth chasing.
    print("\n  IS THE LEVER OPPOSED? nmae under the ORACLE per-target scale (needs the truth,")
    print("  so this is a diagnostic ceiling, not an arm). Mean-normalised, same as above.")
    P, T, ptar_o, _ = delivered(surface)
    s_star = np.maximum(per_target_star(P, T), 0.0)
    s_norm = s_star / float(np.mean(s_star))
    print(f"  {'arm':<52} {'nmae':>8} {'d nmae':>8} {'d scaled':>9}")
    v0 = members(h0)["nmae"]
    for label, ww in (("ORACLE max(s_p*,0), mean-normalised", s_norm),
                      ("ORACLE, re-levelled at its own nmae optimum", None)):
        if ww is None:
            best = min(members(m * mass_fix(s_norm[:, None] * D, b, gate, mode="minnorm"))["nmae"]
                       for m in np.linspace(0.5, 1.8, 27))
            val = best
        else:
            val = members(mass_fix(ww[:, None] * D, b, gate, mode="minnorm"))["nmae"]
        print(f"  {label:<52} {val:8.4f} {val - v0:+8.4f} "
              f"{-(val - v0) / NMAE_SPAN / 6:+9.5f}")

    # ---------------------------------------------------------------- the falsifier
    # The nmae cost could be an artefact of holding the LEVEL fixed: production already has a
    # free global amplitude (VCC_ALPHA), so each shape must be priced at its OWN best level,
    # not at the incumbent's. If the cost vanishes there, the verdict flips.
    print("\n  FALSIFIER — each shape re-levelled at its own nmae-optimal global amplitude")
    print(f"  {'arm':<52} {'m*':>6} {'nmae at m*':>11} {'d nmae':>8} {'d scaled':>9}")
    grid = np.concatenate([np.linspace(0.5, 1.6, 23)])
    ref = None
    for name, h in arms.items():
        vals = []
        for m in grid:
            vals.append(members(m * h)["nmae"])
        j = int(np.argmin(vals))
        if ref is None:
            ref = vals[j]
        print(f"  {name:<52} {grid[j]:6.3f} {vals[j]:11.4f} {vals[j] - ref:+8.4f} "
              f"{-(vals[j] - ref) / NMAE_SPAN / 6:+9.5f}")

    print("\n  `d scaled` is the estimated contribution to OVERALL (nmae is lower-is-better,")
    print("  so a positive nmae delta is a LOSS). n_real on this table is 68-92 against the")
    print("  challenge's regime, so per AGENTS.md 11.1 rule 9 the ABSOLUTE fid/jac/reach")
    print("  values here are unreadable; nmae's gate is the reference's own significant set")
    print("  and its RATIO is the readable part, but it is n = 1 surface.")


if __name__ == "__main__":
    price_board()
    price_nmae()

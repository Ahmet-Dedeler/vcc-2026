"""158 — the deployable rule: stability grid, gain decomposition, bootstrap, null, six members.

Everything reported here is either
  (a) evaluated on the surface whose hyperparameters came from the OTHER surface, or
  (b) a ZERO-degree-of-freedom configuration -- production's already-shipped
      `--shrink-pert-key donor-cos` at its already-frozen `PERT_SHRINK_GAMMA = 1.0`,
      `PERT_SHRINK_FLOOR = 0.3`, `PERT_SHRINK_MISSING = damp` -- which spends no fitting
      budget on either surface and is therefore readable on both.
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
from massfix import mass_fix, production_effect  # noqa: E402
from oracle import (SURFACES, cos_agg, delivered, global_star, pds_of,  # noqa: E402
                    per_target_star, sse_of)
from screen import render  # noqa: E402

WORK = os.path.join(HERE, "work")
# vcc2026 mse anchors, AGENTS.md section 4: baseline b, replicate r. Lower raw is better and
# the metric is clamped to [0, 1] after scaling, so raw > b scales to exactly 0.
MSE_B, MSE_R = 0.989, 0.036


def scaled_mse(raw: float) -> float:
    return float(np.clip((raw - MSE_B) / (MSE_R - MSE_B), 0.0, 1.0))


def main() -> None:
    with open(os.path.join(WORK, "keys.json")) as fh:
        KJ = json.load(fh)
    S = {}
    for surface in SURFACES:
        P, T, ptar, aux = delivered(surface)
        S[surface] = dict(P=P, T=T, ptar=ptar, aux=aux,
                          keys={k: np.asarray(v, dtype=np.float64)
                                for k, v in KJ[surface]["keys"].items()},
                          base=R.score(np.ones(P.shape[0]), P, T),
                          s=per_target_star(P, T), lab=R.clusters(P))

    report = {}

    # -------------------------------------------------------- 1. stability of the grid
    print("=" * 100)
    print("GRID STABILITY for the two donor-cosine keys — d sse on BOTH surfaces at every")
    print("(gamma, floor). A key whose optimum is a spike is a key that will not transfer.")
    for k in ("donor_cos", "donor_cos_live"):
        print(f"\n  {k}")
        print("    gamma  " + "".join(f"{f'floor {f}':>19s}" for f in R.FLOORS))
        print("           " + "".join(f"{'hct':>9s}{'hek':>10s}" for _ in R.FLOORS))
        for g in R.GAMMAS:
            cells = ""
            for f in R.FLOORS:
                for surf in SURFACES:
                    d = S[surf]
                    v = d["base"]["sse"] - R.score(R.weight(d["keys"][k], g, f, 1.0),
                                                   d["P"], d["T"])["sse"]
                    cells += f"{1e6 * v:9.1f}" if surf == "hct116" else f"{1e6 * v:10.1f}"
            print(f"    {g:5.2f}  {cells}     (units: 1e-6 of sse)")

    # ------------------------------- 2. the zero-DOF production default, and the recommendation
    print("\n" + "=" * 100)
    print("ARMS. `d sse` positive = mse improves. `capture` = fraction of the POSITIVE-")
    print("constrained oracle gain (sse_global - sse_{s*>=0}) that the arm actually collects.")
    arms = {
        "production default: donor_cos, gamma 1.0, floor 0.3 (0 DOF)":
            ("donor_cos", 1.0, 0.3, "cont"),
        "donor_cos, gamma 2.0, floor 0.3":
            ("donor_cos", 2.0, 0.3, "cont"),
        "donor_cos_live, gamma 1.0, floor 0.3 (0 DOF, key swapped)":
            ("donor_cos_live", 1.0, 0.3, "cont"),
        "donor_cos_live, gamma 2.0, floor 0.3  [RECOMMENDED]":
            ("donor_cos_live", 2.0, 0.3, "cont"),
        "donor_cos_live, gamma 2.0, floor 0.5":
            ("donor_cos_live", 2.0, 0.5, "cont"),
        "donor_cos_live, gamma 4.0, floor 0.3":
            ("donor_cos_live", 4.0, 0.3, "cont"),
        "DISCRETE donor_cos_live, worst 50% -> 0.50":
            ("donor_cos_live", 0.50, 0.50, "disc"),
        "DISCRETE donor_cos_live, worst 50% -> 0.15":
            ("donor_cos_live", 0.50, 0.15, "disc"),
        "DISCRETE donor_cos_live, worst 25% -> 0.15":
            ("donor_cos_live", 0.25, 0.15, "disc"),
    }
    print(f"  {'arm':<58} {'surface':>8} {'sse':>10} {'d sse':>10} {'cos':>8} "
          f"{'d cos %':>8} {'capture':>8} {'d pds':>9}")
    out_arms = {}
    for label, (k, a, b, form) in arms.items():
        rec = {}
        for surf in SURFACES:
            d = S[surf]
            w = (R.weight(d["keys"][k], a, b, 1.0) if form == "cont"
                 else R.quantile_shrink(d["keys"][k], a, b, 1.0))
            r = R.score(w, d["P"], d["T"])
            pos = sse_of(np.maximum(d["s"], 0.0), d["P"], d["T"])
            ceil = d["base"]["sse"] - pos
            dsse = d["base"]["sse"] - r["sse"]
            g = global_star(w[:, None] * d["P"], d["T"])
            dp = pds_of((g * w)[:, None] * d["P"], d["T"]) - pds_of(d["P"], d["T"])
            rec[surf] = dict(sse=r["sse"], d_sse=dsse, cos=r["cos"],
                             d_cos_pct=100 * (r["cos"] / cos_agg(d["base"]["sse"]) - 1),
                             capture=dsse / ceil, d_pds=dp)
            print(f"  {label if surf == SURFACES[0] else '':<58} {surf:>8} "
                  f"{r['sse']:10.6f} {dsse:+10.6f} {r['cos']:8.5f} "
                  f"{100 * (r['cos'] / cos_agg(d['base']['sse']) - 1):+7.1f}% "
                  f"{dsse / ceil:8.3f} {dp:+9.1e}")
        out_arms[label] = rec
    report["arms"] = out_arms

    # ---------------------------------------------- 3. bootstrap and permutation, recommended
    print("\n" + "=" * 100)
    print("BOOTSTRAP (target-cluster, 20 Ward clusters of our own predicted directions,")
    print("2000 resamples, paired, global scale re-optimised inside each resample) and")
    print("PERMUTATION NULL (2000 shuffles of the key across targets, form and")
    print("hyperparameters held fixed).")
    boots = {}
    for label, (k, a, b, form) in (
            ("donor_cos_live g2 f0.3", ("donor_cos_live", 2.0, 0.3, "cont")),
            ("donor_cos g1 f0.3 (0 DOF)", ("donor_cos", 1.0, 0.3, "cont")),
            ("DISCRETE live worst50->0.15", ("donor_cos_live", 0.50, 0.15, "disc"))):
        rec = {}
        for surf in SURFACES:
            d = S[surf]
            w = (R.weight(d["keys"][k], a, b, 1.0) if form == "cont"
                 else R.quantile_shrink(d["keys"][k], a, b, 1.0))
            lo, med, hi = R.boot_gain(w, d["P"], d["T"], d["lab"])
            p, null = R.perm_null(d["keys"][k], a, b, 1.0, d["P"], d["T"],
                                  discrete=(form == "disc"))
            obs = d["base"]["sse"] - R.score(w, d["P"], d["T"])["sse"]
            rec[surf] = dict(obs=obs, lo=lo, med=med, hi=hi, p=p,
                             null_med=float(np.median(null)),
                             null_p95=float(np.percentile(null, 95)))
            print(f"  {label:30s} {surf:>8}  obs {obs:+.6f}  "
                  f"boot95 [{lo:+.6f}, {hi:+.6f}]  perm p {p:.4f}  "
                  f"null med {np.median(null):+.6f} p95 {np.percentile(null, 95):+.6f}")
        boots[label] = rec
    report["bootstrap"] = boots

    # --------------------------------------- 4. what this is worth on the SCALED mse member
    print("\n" + "=" * 100)
    print("PRICING THE `mse` MEMBER HONESTLY (AGENTS.md section 4 anchors b=0.989, r=0.036,")
    print("clamped to [0,1]). Harness cosines are 3-5x below the board-implied 0.072 because")
    print("the harness truth carries 150 cells/target against the challenge's 400 (156 caveat),")
    print("so the RELATIVE cosine gain is extrapolated to the board and the absolute sse is not.")
    board_cos = 0.072
    print(f"  {'scenario':<44} {'cos_agg':>9} {'min raw mse':>12} {'scaled mse':>11} "
          f"{'overall/6':>10}")
    for label, mult in (("board today (implied)", 1.0),
                        ("+ recommended rule, worst surface", None),
                        ("+ recommended rule, best surface", None),
                        ("+ FULL ORACLE per-target scale, worst surface", None),
                        ("+ FULL ORACLE per-target scale, best surface", None)):
        if label == "board today (implied)":
            c = board_cos
        elif "recommended" in label:
            rs = [out_arms["donor_cos_live, gamma 2.0, floor 0.3  [RECOMMENDED]"][s]["d_cos_pct"]
                  for s in SURFACES]
            c = board_cos * (1 + (min(rs) if "worst" in label else max(rs)) / 100)
        else:
            rs = []
            for s in SURFACES:
                d = S[s]
                rs.append(100 * (cos_agg(sse_of(d["s"], d["P"], d["T"]))
                                 / cos_agg(d["base"]["sse"]) - 1))
            c = board_cos * (1 + (min(rs) if "worst" in label else max(rs)) / 100)
        raw = 1.0 - c * c
        print(f"  {label:<44} {c:9.5f} {raw:12.6f} {scaled_mse(raw):11.4f} "
              f"{scaled_mse(raw) / 6:10.5f}")
    print("  NOTE the challenge `mse` also carries the jackknife/rho correction; this is the")
    print("  UNCORRECTED core at its own best global amplitude, i.e. an upper bound on the")
    print("  reachable raw value. AGENTS.md 4.1b: raw stays above the 0.989 baseline, hence")
    print("  scaled 0, until aggregate cosine clears 0.105.")

    # ------------------------------------------ 5. the other five members, rule 12
    print("\n" + "=" * 100)
    print("ALL SIX MEMBERS for the recommended rule (AGENTS.md 11.1 rule 12)")
    k, a, b = "donor_cos_live", 2.0, 0.3
    for surf in SURFACES:
        d = S[surf]
        P, T = d["P"], d["T"]
        w = R.weight(d["keys"][k], a, b, 1.0)
        # pds under the idealised post-render scale, and under the production pre-render path
        D, Tf, bc, b_ref, gate, mask, ptar = production_effect(surf)
        pre = render(mass_fix(w[:, None] * D, bc, gate, mode="minnorm"), bc)[:, mask]
        print(f"  {surf}: pds post-render scale {pds_of(w[:, None] * P, T):.10f} "
              f"(d {pds_of(w[:, None] * P, T) - pds_of(P, T):+.2e})  |  "
              f"pds production pre-render path {pds_of(pre, T):.10f} "
              f"(d {pds_of(pre, T) - pds_of(P, T):+.2e})")
        # within-target orderings are invariant under a positive row scale -> reach is exact
        o0 = [np.argsort(np.argsort(-np.abs(P[i]))) for i in range(P.shape[0])]
        o1 = [np.argsort(np.argsort(-np.abs((w[:, None] * P)[i]))) for i in range(P.shape[0])]
        same = all(np.array_equal(x, y) for x, y in zip(o0, o1))
        print(f"          within-target |effect| ordering identical on all "
              f"{P.shape[0]} targets: {same}  -> `reach` keyed on |effect| is EXACTLY invariant")

    with open(os.path.join(WORK, "final.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'final.json')}")


if __name__ == "__main__":
    main()

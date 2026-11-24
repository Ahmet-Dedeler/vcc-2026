"""189 Phase 0 — reproduce the coordinator's pilot table.

If this does not reproduce, the experiment aborts (mandate: "Phase 0 ... If you cannot, abort
and say so").  Every downstream script reads `work/pilot189.npz`, so a drift here is caught
once rather than propagating into three priced uses.

Coordinator's pilot, 2026-09-03 02:55 UTC:
    common genes / three-way targets          7,976 / 9,073
    genes with support                        7,916
    r(HCT116, HEK293T)  mean 0.5379 sd 0.0467 p5 0.499 p50 0.530 p95 0.594
    spearman(r_HE, r_HK) +0.3833  (r_HE, r_EK) +0.3359  (r_HK, r_EK) +0.4763
    held-out decile spread of r_HE read on r_HK   0.5515 vs 0.5113 = +0.0402
    spearman(r_HE, b_ctrl) +0.5215  (mean |effect|) +0.1125  (across-target sd) +0.1578
    partialled: resid r_HE vs resid r_HK +0.2133,  vs resid r_EK +0.2251
    held-out decile spread of the residual prior +0.0173,  of expression alone +0.0505
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from common189 import (WORK, bvec, centre, decile_spread, gate_mask, load_table, pair_agreement,
                       panel_targets, partial_rank, spearman, submatrix)

LINES = ["orion_hct116", "orion_hek293t", "k562_gwps"]
SHORT = {"orion_hct116": "HCT", "orion_hek293t": "HEK", "k562_gwps": "K562"}
MIN_SUPPORT = 50            # "genes with >50 supported pairs"


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    tabs = {t: load_table(t) for t in LINES}
    assert all(tabs.values()), "missing donor table"

    genes = sorted(set.intersection(*[set(t["genes"]) for t in tabs.values()]))
    targets = sorted(set.intersection(*[set(t["perts"]) for t in tabs.values()]))
    print(f"common genes {len(genes)}   three-way targets {len(targets)}")
    assert len(genes) == 7976, f"gene axis {len(genes)} != 7976"
    assert len(targets) == 9073, f"target axis {len(targets)} != 9073"

    panel = panel_targets()
    print(f"challenge-panel targets inside the three-way pool: {len(panel & set(targets))}"
          f" of {len(panel)}")

    D, G, B = {}, {}, {}
    for t in LINES:
        A = centre(submatrix(tabs[t], targets, genes).astype(np.float32))
        D[t] = A
        G[t] = gate_mask(A)
        B[t] = bvec(tabs[t], genes)
        print(f"  {SHORT[t]:<5} block {A.shape} gated {G[t].mean():.4f} "
              f"|effect| mean {np.abs(A).mean():.5f}")

    pairs = {}
    for a, b in [("orion_hct116", "orion_hek293t"), ("orion_hct116", "k562_gwps"),
                 ("orion_hek293t", "k562_gwps")]:
        k, n = pair_agreement(D[a], D[b], G[a], G[b])
        pairs[f"{SHORT[a]}{SHORT[b]}"] = (k, n)

    kHE, nHE = pairs["HCTHEK"]
    kHK, nHK = pairs["HCTK562"]
    kEK, nEK = pairs["HEKK562"]

    sup = (nHE > MIN_SUPPORT) & (nHK > MIN_SUPPORT) & (nEK > MIN_SUPPORT)
    print(f"\ngenes with support (>{MIN_SUPPORT} pairs in all three): {int(sup.sum())}")

    rHE = kHE[sup] / nHE[sup]
    rHK = kHK[sup] / nHK[sup]
    rEK = kEK[sup] / nEK[sup]

    print(f"r(HCT,HEK)  mean {rHE.mean():.4f} sd {rHE.std(ddof=1):.4f} "
          f"p5 {np.percentile(rHE, 5):.3f} p50 {np.percentile(rHE, 50):.3f} "
          f"p95 {np.percentile(rHE, 95):.3f}")
    s_he_hk = spearman(rHE, rHK)
    s_he_ek = spearman(rHE, rEK)
    s_hk_ek = spearman(rHK, rEK)
    print(f"spearman(r_HE, r_HK) {s_he_hk:+.4f}   (r_HE, r_EK) {s_he_ek:+.4f}   "
          f"(r_HK, r_EK) {s_hk_ek:+.4f}")

    hi, lo, sp = decile_spread(rHE, rHK)   # unweighted: the pilot's convention
    print(f"held-out decile of r_HE on r_HK: top {hi:.4f} bottom {lo:.4f} spread {sp:+.4f}")

    # ---- confounders -----------------------------------------------------------------
    # Computed on the PAIR whose agreement is being explained (HCT116, HEK293T), mean over the
    # two lines.  `work/probe_confounder.py` enumerated 32 conventions; this one lands on the
    # pilot's +0.1125 / +0.1578 to four decimals and the three-line variant does not.
    # `b_ctrl` is HCT116's alone (+0.5215 to four decimals against 11 alternatives in
    # `work/probe_bctrl.py`) -- HCT116 is the line common to r_HE and the held-out r_HK.
    PAIR = ["orion_hct116", "orion_hek293t"]
    bmean = B["orion_hct116"][sup]
    magn = np.mean([np.abs(D[t]).mean(axis=0) for t in PAIR], axis=0)[sup]
    sdev = np.mean([D[t].std(axis=0) for t in PAIR], axis=0)[sup]
    print(f"\nspearman(r_HE, b_ctrl)        {spearman(rHE, bmean):+.4f}")
    print(f"spearman(r_HE, mean|effect|)  {spearman(rHE, magn):+.4f}")
    print(f"spearman(r_HE, across-tgt sd) {spearman(rHE, sdev):+.4f}")

    X = np.column_stack([bmean, magn, sdev])
    eHE = partial_rank(rHE, X)
    eHK = partial_rank(rHK, X)
    eEK = partial_rank(rEK, X)
    print(f"partialled: spearman(resid r_HE, resid r_HK) {spearman(eHE, eHK):+.4f}   "
          f"(resid r_EK) {spearman(eHE, eEK):+.4f}")

    _, _, sp_res = decile_spread(eHE, rHK)
    _, _, sp_expr = decile_spread(bmean, rHK)
    print(f"held-out decile spread: residual prior {sp_res:+.4f}   expression alone "
          f"{sp_expr:+.4f}")

    out = dict(genes=np.array(genes), targets=np.array(targets), support=sup,
               kHE=kHE, nHE=nHE, kHK=kHK, nHK=nHK, kEK=kEK, nEK=nEK,
               b_hct=B["orion_hct116"], b_hek=B["orion_hek293t"], b_k562=B["k562_gwps"],
               magn_all=np.mean([np.abs(D[t]).mean(axis=0) for t in LINES], axis=0),
               sd_all=np.mean([D[t].std(axis=0) for t in LINES], axis=0))
    np.savez_compressed(os.path.join(WORK, "pilot189.npz"), **out)

    got = dict(n_genes=len(genes), n_targets=len(targets), n_support=int(sup.sum()),
               rHE_mean=float(rHE.mean()), rHE_sd=float(rHE.std(ddof=1)),
               s_he_hk=s_he_hk, s_he_ek=s_he_ek, s_hk_ek=s_hk_ek,
               decile_hi=hi, decile_lo=lo, decile_spread=sp,
               s_expr=spearman(rHE, bmean), s_mag=spearman(rHE, magn),
               s_sd=spearman(rHE, sdev),
               partial_hk=spearman(eHE, eHK), partial_ek=spearman(eHE, eEK),
               spread_resid=sp_res, spread_expr=sp_expr)
    json.dump(got, open(os.path.join(WORK, "pilot189.json"), "w"), indent=2)

    EXP = dict(n_genes=7976, n_targets=9073, n_support=7916, rHE_mean=0.5379, rHE_sd=0.0467,
               s_he_hk=0.3833, s_he_ek=0.3359, s_hk_ek=0.4763, decile_hi=0.5515,
               decile_lo=0.5113, decile_spread=0.0402, s_expr=0.5215, s_mag=0.1125,
               s_sd=0.1578, partial_hk=0.2133, partial_ek=0.2251, spread_resid=0.0173,
               spread_expr=0.0505)
    TOL = dict(n_genes=0, n_targets=0, n_support=0)
    print(f"\n{'row':<16}{'here':>10}{'pilot':>10}{'delta':>10}")
    bad = []
    for k_, v in EXP.items():
        g = got[k_]
        tol = TOL.get(k_, 0.006)
        ok = abs(g - v) <= tol
        bad += [] if ok else [k_]
        print(f"{'OK  ' if ok else 'FAIL'}{k_:<12}{g:>10.4f}{v:>10.4f}{g - v:>+10.4f}")
    if bad:
        print(f"\nPILOT DID NOT REPRODUCE: {bad}")
        sys.exit(3)
    print("\npilot reproduces; the statistic is the coordinator's")


if __name__ == "__main__":
    main()

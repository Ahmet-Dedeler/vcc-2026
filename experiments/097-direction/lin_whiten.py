"""Isolate the whitener: ONE analytic array, gamma swept, everything else held fixed.

`experiments/031-somatic/run.py` builds the shipped profile as

    E = MultiDonorModel(..., gene_norm=False, eb=False).effect(p, alpha=1, beta=0, uniform)
    E[:, wcols] = Whitener.apply(E[:, wcols], gamma=0.65, lam=0, renorm=False), median-renormed
    E *= alpha * effect_scale / median||E||          # a global scalar; kappa cannot see it

With `gene_norm=False, eb=False` that first line is EXACTLY a per-donor-renormalised mean of the
donors' gene-specific effects -- i.e. `lin_null.py`'s "plain renorm-mean" up to a global scalar.
So the whole analytic difference between the incumbent's profile and that arm is the whitener,
and sweeping gamma on the SAME array measures its cost on the L2 objective directly.

`gamma = 0.65` was fitted for `pds_cosine`, a RETRIEVAL metric, in experiment 017 -- before
`expr_mse` left its dead zone.  A lever closed under one objective is not closed under a new one
(RESEARCH-METHOD rule 9).

`wcols` is rebuilt from `data/challenge/gene_names.csv` + the two donors' gene lists minus the
300 panel targets; it reproduces the shipped operator's 17,942-gene axis exactly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import lin_kit as K
from lin_null import donor_S_on_targets

sys_path_note = "whitener: data/donors/whitener_union.npz (17942 x 17942)"
GAMMAS = (0.0, 0.15, 0.3, 0.45, 0.65, 0.85, 1.0)
DONORS = ("k562_gwps", "orion_hek293t")


def wcols_axis():
    tg = pd.read_csv(f"{K.ROOT}/data/challenge/gene_names.csv")["gene_name"].to_numpy()
    panel = set(pd.read_csv(f"{K.ROOT}/data/challenge/pert_counts.csv")["target_gene"]
                .astype(str))
    ti = {str(g): i for i, g in enumerate(tg)}
    union = set()
    for d in DONORS:
        union |= {g for g in open(f"{K.ROOT}/data/donors/{d}/genes.txt").read().split("\n")
                  if g in ti}
    union -= panel
    return tg, np.array(sorted(ti[g] for g in union), dtype=np.int64)


def main():
    from vcc2026.whiten import Whitener
    tg, wcols = wcols_axis()
    print(f"  union axis rebuilt: {wcols.size} genes (shipped operator is 17942) -- "
          f"{'MATCH' if wcols.size == 17942 else 'MISMATCH'}")
    W = Whitener.load(f"{K.ROOT}/data/donors/whitener_union.npz")
    print(f"  {sys_path_note}; eval range [{W.evals.min():.3e}, {W.evals.max():.3e}], "
          f"condition {W.evals.max()/W.evals.min():.3e}")
    cpos = {str(g): i for i, g in enumerate(tg)}
    wpos = {int(c): j for j, c in enumerate(wcols)}

    for name in ("h1-test", "orion"):
        t = K.prep(name)
        P = t["dh"].shape[0]
        r0 = K.shape_full(t["dh"], t)
        k0 = r0["kappa"]
        corr = np.sqrt(r0["B"] / (r0["B"] - float((t["Chat"] * t["w"]).sum()) * (1 - 1 / P)))
        # the incumbent's own analytic array, before the whitener
        acc = np.zeros((P, t["genes"].size))
        cnt = np.zeros(P)
        for dn in DONORS:
            M, cov, med, _ = donor_S_on_targets(t, dn)
            acc[cov] += M[cov] / med
            cnt[cov] += 1
        acc[cnt > 0] /= cnt[cnt > 0, None]
        # harness gene -> position within the whitener's own axis
        col = np.full(t["genes"].size, -1, np.int64)
        for j, g in enumerate(t["genes"]):
            c = cpos.get(str(g))
            if c is not None and c in wpos:
                col[j] = wpos[c]
        hit = col >= 0
        print(f"\n################ {name} (P={P}) ################")
        print(f"  incumbent (whitened + emitted) kappa {k0:.5f}, emission-corrected "
              f"{k0*corr:.5f}, eps {r0['eps']:+.5f}")
        print(f"  {int(hit.sum())}/{t['genes'].size} harness genes sit on the whitener axis; "
              f"the rest are left untouched, exactly as run.py does")
        # N is the emitter's own sampling variance inside the budget: sum_p w_p Chat_p (1-1/n).
        # Adding it to an ANALYTIC arm's budget projects that arm to what the scorer would read
        # after it had been emitted by the same policy, which is the only fair comparison to a
        # measured, emitted incumbent.
        N = float((t["Chat"] * t["w"]).sum()) * (1 - 1 / P)
        print(f"    emitter noise inside the budget  N = {N:.5f}   "
              f"(= {N/r0['B']*100:.1f}% of the incumbent's B)")
        print(f"    {'gamma':>7}{'kappa':>10}{'M vs e0':>9}{'M/corr':>9}{'eps':>10}"
              f"{'cos':>9}{'mse_min':>9}{'X':>10}{'B':>9}{'k emit':>9}{'M emit':>8}")
        Ew = np.zeros((P, wcols.size))
        Ew[:, col[hit]] = acc[:, hit]
        med0 = float(np.median(np.linalg.norm(Ew, axis=1))) or 1.0
        rowsg = []
        for g in GAMMAS:
            E = acc.copy()
            if g > 0:
                out = W.apply(Ew, g, 0.0, renorm=False)
                out *= med0 / (float(np.median(np.linalg.norm(out, axis=1))) or 1.0)
                E[:, hit] = out[:, col[hit]]
            r = K.shape_full(E, t)
            r["gamma"] = g
            r["k_emit"] = r["X"] / np.sqrt((r["B"] + N) * r["sumD"])
            rowsg.append(r)
        # The analytic array is at median row norm ~0.8; the SHIPPED profile is emitted at
        # alpha * effect_scale = 1.08 * 6.817.  B scales as amplitude^2 while N does not, so the
        # emitted projection is only meaningful once the amplitude is calibrated.  Calibrate it
        # by requiring the gamma = 0.65 row (which IS the shipped model) to reproduce the
        # incumbent's own SIGNAL budget, B_incumbent - N.
        ref = [r for r in rowsg if r["gamma"] == 0.65][0]
        s2 = (r0["B"] - N) / ref["B"]
        nrel = N / s2
        for r in rowsg:
            r["k_emit"] = (np.sqrt(s2) * r["X"]) / np.sqrt((s2 * r["B"] + N) * r["sumD"])
        for r in rowsg:
            print(f"    {r['gamma']:>7.2f}{r['kappa']:>10.5f}{r['kappa']/k0:>9.3f}"
                  f"{r['kappa']/(k0*corr):>9.3f}{r['eps']:>+10.5f}{r['cos']:>9.5f}"
                  f"{r['mse_min']:>9.4f}{r['X']:>10.5f}{r['B']:>9.5f}"
                  f"{r['k_emit']:>9.5f}{r['k_emit']/ref['k_emit']:>8.3f}")
        med_an = float(np.median(np.linalg.norm(acc, axis=1)))
        print(f"    amplitude calibration: s = {np.sqrt(s2):.3f} (median effect norm "
              f"{med_an:.3f} -> {med_an*np.sqrt(s2):.3f}); the shipped profile is built at "
              f"alpha*effect_scale = {1.08*6.817:.3f} before the min_cpm gate, so the gate and "
              f"the log1p round trip compress it by {med_an*np.sqrt(s2)/(1.08*6.817):.2f}x")
        print(f"    SELF-CHECK: gamma=0.65 projected-to-emitted kappa {ref['k_emit']:.5f} vs "
              f"the MEASURED emitted incumbent {k0:.5f} (ratio {ref['k_emit']/k0:.3f}) -- the "
              f"gap is X only, i.e. the gate + fallback, since B was calibrated to match.")
        print(f"    WHITENER COST, decomposed: X falls {ref['X']/rowsg[0]['X']:.3f}x from "
              f"gamma 0 to 0.65 while B moves only {ref['B']/rowsg[0]['B']:.3f}x -- the "
              f"operator destroys alignment, it does not buy budget.")
        K._CACHE.clear()


if __name__ == "__main__":
    main()

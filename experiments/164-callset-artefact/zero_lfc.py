"""164 — the sub-gate band's SCORED direction, exactly, on the real challenge controls.

`blocks.py` established that in every challenge context the scorer's DE universe
(mean-per-cell CPM > 5) strictly contains production's spread gate (PSEUDOBULK CPM > 10, via
`--min-cpm 10.0`), and that 794 / 958 / 997 genes sit in the band between them.  For those
genes `profile_from_effect` sets the predicted effect to EXACTLY zero
(make_submission_pm_cd4.py:1281), yet the scorer still sees them and still scores their
direction.

So their scored `lfc` carries no model content at all.  It is exactly

    lfc_0(g) = log2( (pb_cpm_g + eps) / (cell_cpm_g + eps) )

  * numerator   -- our emitted cells all carry the same library, and the quasi emitter
                   repairs each gene's panel total exactly, so our per-cell mean CPM is the
                   control PSEUDOBULK cpm: ratio-of-sums.
  * denominator -- the scorer's control mean is the arithmetic mean of per-cell CPM:
                   mean-of-ratios.

Same cells, two normalisations.  That difference IS our predicted log fold change for 8-10%
of every scored gene set.  This script measures it, and measures how often such a gene is
CALLED, with the exact integer emission law and exact Mann-Whitney ties.

KNOWN-ANSWER GATES, both of which could fail:
  1. the block census must reproduce `blocks.log`: 794 / 958 / 997 band genes, 0 genes
     prod-gated but not DE;
  2. the zero-effect call rate at `spread 0.25 / phi 4` must reproduce 157 `zero_effect.log`
     (19.0 / 22.8 / 16.2 %) and 153 §15 (19.0 / 21.6 / 16.0 %).

Controls only -- the challenge ships no reference truth, so PRECISION is not computable here
and is not claimed.  Direction, size and call rate are.  No cell-eval2, no build, no submit.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = "~/Code/virtual-cell-challenge-2026"
N, EPS = 400, 1e-9
SCORER_CPM, PROD_CPM, KEEP_TOP = 5.0, 10.0, 8000
rng = np.random.default_rng(0)


def fbar_for(values, counts, nz_sorted, n_zero, m):
    """Exact Mann-Whitney Fbar of an emitted count vector against the control column."""
    tot = 0.0
    for e, c in zip(values, counts):
        if c == 0:
            continue
        if e == 0.0:
            below, eq = 0, n_zero
        else:
            lo = np.searchsorted(nz_sorted, e, "left")
            hi = np.searchsorted(nz_sorted, e, "right")
            below, eq = n_zero + lo, hi - lo
        tot += c * (below + 0.5 * eq) / m
    return tot / N


def main() -> None:
    print("=" * 108)
    print("(1) BLOCK CENSUS + the normalisation artefact lfc_0, real challenge controls")
    print("=" * 108)
    store = {}
    for ctx in "ABC":
        A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
        X = A.X.tocsc().astype(np.float64)
        lib = np.asarray(X.sum(1)).ravel()
        lib[lib == 0] = 1.0
        L = int(np.median(lib))
        D = (sp.diags(1e6 / lib) @ X).tocsc()          # per-cell CPM
        m = D.shape[0]
        cell_cpm = np.asarray(D.mean(0)).ravel()        # mean-of-ratios (the scorer's control)
        colsum = np.asarray(X.sum(0)).ravel()
        p = colsum / colsum.sum()                       # zero-effect profile = control pseudobulk
        pb_cpm = 1e6 * p                                # ratio-of-sums (our emitted cells)

        de = cell_cpm > SCORER_CPM                      # the scorer's DE universe
        gate = (pb_cpm > PROD_CPM) & (np.log1p(pb_cpm) > 0)   # production's spread gate
        band = de & ~gate
        lfc0 = np.log2((pb_cpm + EPS) / (cell_cpm + EPS))

        print(f"\nctx {ctx}:  {m} control cells, median library {L}")
        print(f"  scorer DE universe (cell-CPM>5) {int(de.sum()):>6}   "
              f"prod gate (pb-CPM>10) {int(gate.sum()):>6}   "
              f"band (DE & not gated) {int(band.sum()):>5}   "
              f"gated & not DE {int((gate & ~de).sum()):>3}")
        print(f"  {'block':<28}{'n':>7}{'median lfc0':>13}{'mean lfc0':>11}"
              f"{'DOWN':>8}{'|lfc0|>0.1':>12}{'|lfc0|>0.5':>12}")
        for nm, msk in (("sub-gate band (5-10 CPM)", band),
                        ("gated (pb-CPM>10)", gate & de),
                        ("whole DE universe", de)):
            v = lfc0[msk]
            print(f"  {nm:<28}{v.size:>7}{np.median(v):>13.4f}{v.mean():>11.4f}"
                  f"{float((v < 0).mean()):>8.1%}{float((np.abs(v) > 0.1).mean()):>12.1%}"
                  f"{float((np.abs(v) > 0.5).mean()):>12.1%}")
        store[ctx] = dict(D=D, m=m, L=L, p=p, de=de, gate=gate, band=band, lfc0=lfc0,
                          cell_cpm=cell_cpm, pb_cpm=pb_cpm)
        del A, X

    print("\n" + "=" * 108)
    print("(2) ZERO-EFFECT CALL RATE and the SCORED direction of what gets called")
    print("    'z-up' is the Wilcoxon sign (157's number).  'lfc-up' is what fid ACTUALLY scores.")
    print("=" * 108)
    print(f"{'ctx':>4}{'block':>26}{'spread/phi':>12}{'Fbar med':>10}{'|z| med':>9}"
          f"{'sig@0':>8}{'rate':>8}{'z-up':>8}{'lfc-up':>8}")
    rates = {}
    for ctx in "ABC":
        s = store[ctx]
        D, m, L, p, lfc0 = s["D"], s["m"], s["L"], s["p"], s["lfc0"]
        fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
        step = 1e6 / L
        # Production's real blocks: the band always sits at tail_spread; the gated genes split
        # keep_top / tail by |effect| -- but at zero effect there is no ranking, so the two
        # spreads are evaluated over the whole DE universe and combined by block size below.
        for label, msk in (("sub-gate band", s["band"]), ("gated (pb>10)", s["gate"] & s["de"])):
            idx = np.flatnonzero(msk)
            for spread, phi in ((0.25, 4.0), (1.0, 4.0), (1.0, 0.0), (0.5, 4.0)):
                fb = np.empty(idx.size)
                for k, g in enumerate(idx):
                    nz = np.sort(D.data[D.indptr[g]:D.indptr[g + 1]])
                    n_zero = m - nz.size
                    T = int(np.floor(N * L * p[g]))
                    t_rand = int(np.floor(T * spread))
                    t_even = T - t_rand
                    b, r = divmod(t_even, N)
                    even = np.full(N, b)
                    even[:r] += 1
                    if t_rand > 0 and phi > 0:
                        draw = rng.poisson(rng.gamma(1.0 / phi, phi * (t_rand / N), N))
                    elif t_rand > 0:
                        draw = rng.poisson(np.full(N, t_rand / N))
                    else:
                        draw = np.zeros(N, dtype=int)
                    vals, cts = np.unique((even + draw) * step, return_counts=True)
                    fb[k] = fbar_for(vals, cts, nz, n_zero, m)
                z = (fb - 0.5) * fac
                sig = np.abs(z) > 2.8
                rates[(ctx, label, spread, phi)] = (float(sig.mean()), idx.size)
                lu = float((lfc0[idx][sig] > 0).mean()) if sig.any() else float("nan")
                print(f"{ctx:>4}{label:>26}{f'{spread:.2f}/phi{phi:g}':>12}"
                      f"{np.median(fb):>10.4f}{np.median(np.abs(z)):>9.1f}"
                      f"{int(sig.sum()):>8d}{sig.mean():>8.1%}"
                      f"{((z[sig] > 0).mean() if sig.any() else 0):>8.1%}{lu:>8.1%}")

    print("\n" + "=" * 108)
    print("(3) THE SHIPPED MIXTURE, with the block composition blocks.py established")
    print("=" * 108)
    print(f"{'ctx':>4}{'policy':>26}{'band':>8}{'keep8000':>10}{'gated tail':>12}"
          f"{'sig@0':>9}{'rate of DE universe':>22}")
    for ctx in "ABC":
        s = store[ctx]
        nband = int(s["band"].sum())
        ngate = int((s["gate"] & s["de"]).sum())
        nkeep = min(KEEP_TOP, ngate)
        ntail = ngate - nkeep
        nde = int(s["de"].sum())
        for nm, (sp_keep, ph_keep), (sp_tail, ph_tail) in (
                ("a15 as shipped", (0.25, 4.0), (1.0, 4.0)),
                ("armA tail phi 0", (0.25, 4.0), (1.0, 0.0)),
                ("armB tail spread 0.5", (0.25, 4.0), (0.5, 4.0)),
                ("armC a8 flat", (0.25, 4.0), (0.25, 4.0))):
            r_band = rates[(ctx, "sub-gate band", sp_tail, ph_tail)][0]
            r_keep = rates[(ctx, "gated (pb>10)", sp_keep, ph_keep)][0]
            r_tail = rates[(ctx, "gated (pb>10)", sp_tail, ph_tail)][0]
            tot = nband * r_band + nkeep * r_keep + ntail * r_tail
            print(f"{ctx:>4}{nm:>26}{nband:>8}{nkeep:>10}{ntail:>12}"
                  f"{tot:>9.0f}{tot / nde:>22.1%}")

    np.savez_compressed(f"{ROOT}/experiments/164-callset-artefact/zero_lfc.npz",
                        **{f"{k}_{c}": store[c][k] for c in "ABC"
                           for k in ("lfc0", "de", "gate", "band", "cell_cpm", "pb_cpm")})


if __name__ == "__main__":
    main()

"""171 — the six-member ledger, board-scaled, with twins averaged.

Three corrections the raw arm table needs before any of it can be believed:

1. **Twins.** h1-test `reach` has a large same-policy seed spread -- measured here at 0.0041
   (ks000 pair), 0.0107 (armC pair) and 0.0123 (kt20000 pair), against arm effects of ~0.020.
   Single-seed rows are not resolvable.  Every reported delta below is a twin MEAN.

2. **The board calibration.** 157 section 5c and 164 section 6 established that h1-test
   OVERSTATES this family's `reach` moves.  The board's own a8 <-> a15 pair reads
   0.1890 vs 0.1811 = -0.0079 raw; the same pair on h1-test reads the twin-mean difference
   below.  Their ratio is the shrink factor, and armC is therefore the known-answer row for
   the calibration itself.

3. **`fid` and `jac` are unreadable raw (rule 9).**  h1-test's n_real is 1-5,484 (median 115,
   mean 793) against the challenge's ~240 with n_pred ~4,424, so `max(n_pred, n_real)` and the
   jaccard union are in a different regime.  Those two members are carried over from 164's
   board-anchored projections for the arm with the matching n_pred, and that adoption is
   stated rather than hidden.

`mse` from `claim171b.py` on the real challenge panel, correct branch (`msebranch171.py`).
`pds` from the pseudobulk digest, which is an identity, not an estimate.
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
from diag171 import Walk, kstar, spec_order  # noqa: E402

W = f"{ROOT}/experiments/171-reach-prefix/work"
D164 = f"{ROOT}/experiments/164-callset-artefact"
REACH_SPAN = 0.96414 - 0.07902          # debench ANCHOR, fitted on 780 board entries
MSE_SPAN = 0.950162566
BOARD_A8_A15 = 0.1890 - 0.1811          # +0.0079 raw: a8's flat law beats a15's on the board


def reach(path):
    Wk = Walk(path)
    out = []
    for p in range(Wk.perts.size):
        pool = np.flatnonzero(Wk.real_sig[p])
        if pool.size == 0:
            continue
        o = spec_order(pool, Wk.genes, Wk.pred_sig[p], Wk.pred_padj[p], Wk.pred_p[p],
                       Wk.pred_lfc[p])
        w = pool[o]
        w = w[Wk.adj[p, w]]
        if w.size == 0:
            out.append(0.0)
            continue
        ok = np.sign(Wk.pred_lfc[p, w]) == np.sign(Wk.real_lfc[p, w])
        out.append(kstar(ok)[0] / pool.size)
    return float(np.mean(out))


# arm -> (seed-A table, seed-B table).  a15's twin lives in 164's directory (same seeds,
# same code path, bit-identical n_pred -- verified by this experiment's gates).
PAIRS = {
    "a15 (incumbent)":      (f"{W}/de171.a15.npz",          f"{D164}/de164.a15_b.npz"),
    "armC a8 flat 0.25":    (f"{W}/de171.armC_a8_flat.npz", f"{W}/de171.armC_a8_flat_b.npz"),
    "kt20000 (all gated)":  (f"{W}/de171.kt20000.npz",      f"{W}/de171.kt20000_b.npz"),
    "ks000 keep_spread 0":  (f"{W}/de171.ks000.npz",        f"{W}/de171.ks000_b.npz"),
    "kt10000 = KNEE arm":   (f"{W}/de171.kt10000.npz",      f"{W}/de171.kt10000_b.npz"),
    "kt9000":               (f"{W}/de171.kt9000.npz",       f"{W}/de171.kt9000_b.npz"),
}
SINGLES = {
    "ks010 keep_spread .10": f"{W}/de171.ks010.npz",
    "kt2000_ks000":          f"{W}/de171.kt2000_ks000.npz",
    "kt2000_ks025":          f"{W}/de171.kt2000_ks025.npz",
    "flat100_phi4":          f"{W}/de171.flat100_phi4.npz",
}
# `claim` from claim171.py / claim171b.py, and the resulting d overall on the correct branch
CLAIM = {
    "a15 (incumbent)": (0.24581, 0.000000),
    "armC a8 flat 0.25": (0.16652, -0.018520),
    "kt20000 (all gated)": (0.18049, -0.018520),
    "ks000 keep_spread 0": (0.19180, -0.018520),
    "ks010 keep_spread .10": (0.21191, -0.015943),
    "kt10000 = KNEE arm": (0.23580, 0.000000),
    "kt9000": (0.24148, 0.000000),
    "kt2000_ks000": (0.24653, 0.000000),
    "kt2000_ks025": (0.24653, 0.000000),
    "flat100_phi4": (0.24653, 0.000000),
}


def main():
    print("TWIN-AVERAGED h1-test raw reach, and the seed spread each average rests on")
    print(f"  {'arm':<24}{'seed A':>9}{'seed B':>9}{'mean':>9}{'|A-B|':>9}"
          f"{'d vs a15':>10}")
    means = {}
    for nm, (a, b) in PAIRS.items():
        ra, rb = reach(a), reach(b)
        means[nm] = (ra + rb) / 2
        print(f"  {nm:<24}{ra:>9.5f}{rb:>9.5f}{means[nm]:>9.5f}{abs(ra-rb):>9.5f}"
              f"{means[nm]-means['a15 (incumbent)']:>+10.5f}")
    for nm, p in SINGLES.items():
        r = reach(p)
        means[nm] = r
        print(f"  {nm:<24}{r:>9.5f}{'--':>9}{r:>9.5f}{'--':>9}"
              f"{r-means['a15 (incumbent)']:>+10.5f}   (single seed)")
    print()

    d_armC = means["armC a8 flat 0.25"] - means["a15 (incumbent)"]
    shrink = BOARD_A8_A15 / d_armC
    print("CALIBRATION (the known-answer row for the transfer itself)")
    print(f"  board a8 - a15 raw reach            {BOARD_A8_A15:+.5f}")
    print(f"  h1-test armC - a15, twin-averaged   {d_armC:+.5f}")
    print(f"  => h1-test OVERSTATES by {1/shrink:.2f}x ; shrink factor {shrink:.4f}")
    print(f"     (164 section 6 derived 0.428 from the same pair with one seed each)")
    print()

    print("BOARD-SCALED SIX-MEMBER LEDGER (one variable per arm, a15's VCC_ALPHA held)")
    print("  pds  = 0.000000 EXACTLY: the delivered per-gene panel total is sha256-identical")
    print("         across all 12 arms on 85/85 targets, and BRIEF section 1 defines pds on")
    print("         that pseudobulk.  This is an identity, not a small measurement.")
    print("  nmae = measured on h1-test (BRIEF section 7 reads the DE lfc, so it CAN move);")
    print("         every arm lands inside the 0.00008 re-seed spread.")
    print("  fid / jac = 164's board-anchored projections for the matching n_pred; raw")
    print("         h1-test values are unreadable under rule 9 and are NOT used.")
    print()
    print(f"  {'arm':<24}{'claim':>9}{'d mse':>10}{'d reach raw':>12}{'d reach':>10}"
          f"{'d fid':>10}{'d jac':>10}{'NET overall':>13}")
    # 164 section 6, board-anchored, for the n_pred-matched arm
    FIDJAC = {
        "armC a8 flat 0.25": (+0.000852, +0.000570),
        "kt20000 (all gated)": (+0.000852, +0.000570),   # n_pred 3588 vs armC's 3589
        "ks000 keep_spread 0": (-0.000200, -0.000150),   # n_pred UP 6034->7051: sign only
        "ks010 keep_spread .10": (-0.000100, -0.000080),
        "kt10000 = KNEE arm": (+0.000600, +0.000400),
        "kt9000": (+0.000400, +0.000260),
        "kt2000_ks000": (-0.000400, -0.000300),
        "kt2000_ks025": (-0.000400, -0.000300),
        "flat100_phi4": (-0.000500, -0.000350),
    }
    for nm in list(PAIRS) + list(SINGLES):
        if nm == "a15 (incumbent)":
            print(f"  {nm:<24}{0.24581:>9.5f}{0.0:>+10.6f}{0.0:>+12.5f}{0.0:>+10.6f}"
                  f"{0.0:>+10.6f}{0.0:>+10.6f}{0.0:>+13.6f}")
            continue
        draw = (means[nm] - means["a15 (incumbent)"]) * shrink
        dre = draw / REACH_SPAN / 6
        cl, dmse = CLAIM.get(nm, (None, None))
        df, dj = FIDJAC.get(nm, (0.0, 0.0))
        if dmse is None:
            print(f"  {nm:<24}{'unpriced':>9}{'?':>10}{draw:>+12.5f}{dre:>+10.6f}"
                  f"{df:>+10.6f}{dj:>+10.6f}{'not priced':>13}")
            continue
        net = dmse + dre + df + dj
        print(f"  {nm:<24}{cl:>9.5f}{dmse:>+10.6f}{draw:>+12.5f}{dre:>+10.6f}"
              f"{df:>+10.6f}{dj:>+10.6f}{net:>+13.6f}")
    print()

    print("THE ONLY mse-FREE POINT ON THE keep_top LADDER")
    print("  claim(8000) = 0.24581, claim(8400) = 0.23006, knee = 0.23268.  Linear in")
    print("  keep_top between those two, the knee is crossed at keep_top ~ 8333, i.e. only")
    print("  333 of the 976 gated tail genes (34%) can be promoted for free.")
    frac = (0.24581 - 0.23268) / (0.24581 - 0.23006)
    kt = 8000 + 400 * frac
    share = (kt - 8000) / 976
    draw_full = (means["kt20000 (all gated)"] - means["a15 (incumbent)"]) * shrink
    print(f"    knee at keep_top = {kt:.0f}   ({share*100:.1f}% of the gated tail promoted)")
    print(f"    full promotion is worth d raw reach {draw_full:+.5f} on the board scale;")
    print(f"    {share*100:.1f}% of it, if the response is linear in the block size, is "
          f"{draw_full*share:+.5f}")
    print(f"    => d overall {draw_full*share/REACH_SPAN/6:+.6f}, against rule 13's")
    print(f"       shape-only ship gate of +0.0019 and an h1-test re-seed spread of 0.0107")
    print(f"       raw reach ({0.0107*shrink/REACH_SPAN/6:+.6f} overall) on the SAME family.")


if __name__ == "__main__":
    main()

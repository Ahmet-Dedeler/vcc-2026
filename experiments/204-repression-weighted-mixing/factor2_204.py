"""204 Phase 1d — is the rank-1 rejection real, and is the between-block disagreement
bigger than sampling?

`factor204.py` fitted `rho_ij = c_i c_j` on three overlapping line blocks and got
`a_rel(HCT116)` from 0.68 to 1.69.  Two things have to be separated before that can be
called a model failure rather than noise:

  1. within a block, is the rank-1 residual bigger than the target-resampling error?
  2. across blocks, is the disagreement in `a_rel` bigger than each block's own interval?

Both are answered by resampling TARGETS (the unit the metric averages over, and the unit
194's own bootstrap used), recomputing the whole Gram, refitting, and carrying the
tetrads through.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))

from factor204 import BLOCKS, CELL_LINE, LINES, PLATFORM, fit_rank1  # noqa: E402

BLOCKS = dict(BLOCKS)
BLOCKS["h1_block"] = ("k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025")
BLOCKS["panel_h1_block"] = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2",
                            "h1_2025")
from prep204 import DONORS, prepared, shared_targets  # noqa: E402
from plolo import DONOR_PATHS, PANEL_SET  # noqa: E402

NB = 2000


def main() -> None:
    tg = shared_targets()
    P = {}
    for name, path in LINES.items():
        DONOR_PATHS[name] = path
        p = prepared(name, tg, "shared") if name in DONORS else prepared(
            name, tg, "allshared")
        P[name] = (np.asarray(p["effects"], dtype=np.float32),
                   np.asarray(p["target_available"], bool))
    support = {n: (np.abs(P[n][0]).sum(0) > 0) for n in LINES}

    out: dict = {}
    for bname, block in BLOCKS.items():
        rows = np.flatnonzero(np.logical_and.reduce([P[n][1] for n in block]))
        cols = np.flatnonzero(np.logical_and.reduce([support[n] for n in block]))
        if len(rows) < 30:
            continue
        k = len(block)
        M = [P[n][0][np.ix_(rows, cols)].astype(np.float64) for n in block]
        # per-target pair sums: pair[t, i, j] = <e_i(t), e_j(t)>
        pair = np.zeros((len(rows), k, k))
        for i in range(k):
            for j in range(i, k):
                pair[:, i, j] = pair[:, j, i] = (M[i] * M[j]).sum(1)
        del M

        def stat(G):
            """Two different objects come out of the same fit and must not be confused:
            `c_i` is donor i's RELIABILITY (the fraction of its row that is shared
            signal), and `A_i = c_i * ||d_i||` is its shared-signal AMPLITUDE, which is
            the coefficient on `s` when the consensus averages the production-frame rows.
            The mandate's `a_i` is the amplitude."""
            d = np.sqrt(np.diag(G))
            R = G / np.outer(d, d)
            c = fit_rank1(R)
            resid = R - np.outer(c, c)
            off = ~np.eye(k, dtype=bool)
            return R, c, resid, float(np.sqrt((resid[off] ** 2).mean())), c * d

        R0, c0, res0, rms0, A0 = stat(pair.sum(0))
        rng = np.random.default_rng(20260903)
        bs_c, bs_rms, bs_res, bs_a, bs_A = [], [], [], [], []
        i0 = block.index(DONORS[0])
        for _ in range(NB):
            idx = rng.integers(0, len(rows), len(rows))
            R, c, resid, rms, A = stat(pair[idx].sum(0))
            bs_c.append(c)
            bs_rms.append(rms)
            bs_res.append(resid[np.triu_indices(k, 1)])
            bs_a.append(c / c[i0])
            bs_A.append(A / A[i0])
        bs_c = np.array(bs_c)
        bs_res = np.array(bs_res)
        bs_a = np.array(bs_a)
        bs_A = np.array(bs_A)

        pairs = []
        for n, (i, j) in enumerate(itertools.combinations(range(k), 2)):
            v = bs_res[:, n]
            pairs.append(dict(
                pair=[block[i], block[j]], rho=float(R0[i, j]),
                resid=float(res0[i, j]),
                ci=[float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
                z=float(res0[i, j] / max(v.std(), 1e-12)),
                same_cell_line=CELL_LINE[block[i]] == CELL_LINE[block[j]],
                same_platform=PLATFORM[block[i]] == PLATFORM[block[j]]))

        rec = dict(
            lines=list(block), n_targets=int(len(rows)), n_genes=int(len(cols)),
            n_panel=int(sum(1 for r in rows if tg[r] in PANEL_SET)),
            rank1_rms_resid=rms0,
            rank1_rms_ci=[float(np.percentile(bs_rms, 2.5)),
                          float(np.percentile(bs_rms, 97.5))],
            mean_abs_rho=float(np.abs(R0[~np.eye(k, dtype=bool)]).mean()),
            pairs=pairs,
            c_rel={d: float(c0[block.index(d)] / c0[block.index(DONORS[0])])
                   for d in DONORS if d in block},
            c_rel_ci={d: [float(np.percentile(bs_a[:, block.index(d)], 2.5)),
                          float(np.percentile(bs_a[:, block.index(d)], 97.5))]
                      for d in DONORS if d in block},
            a_rel={d: float(A0[block.index(d)] / A0[block.index(DONORS[0])])
                   for d in DONORS if d in block},
            a_rel_ci={d: [float(np.percentile(bs_A[:, block.index(d)], 2.5)),
                          float(np.percentile(bs_A[:, block.index(d)], 97.5))]
                      for d in DONORS if d in block},
            row_norm={d: float(np.sqrt(np.diag(pair.sum(0))[block.index(d)]))
                      for d in DONORS if d in block})
        out[bname] = rec
        print(f"\n[{bname}] n_t={len(rows)} n_panel={rec['n_panel']} n_g={len(cols)}")
        print(f"  rank-1 RMS resid {rms0:.5f}  CI {rec['rank1_rms_ci']}  "
              f"on mean|rho| {rec['mean_abs_rho']:.5f} "
              f"({100*rms0/rec['mean_abs_rho']:.1f}%)")
        for lab in ("c_rel", "a_rel"):
            print(f"  {lab}: " + "  ".join(
                f"{d.split('_')[-1]} {rec[lab][d]:.3f} "
                f"[{rec[lab + '_ci'][d][0]:.3f},{rec[lab + '_ci'][d][1]:.3f}]"
                for d in rec[lab]))
        for p in sorted(pairs, key=lambda x: -abs(x["z"]))[:6]:
            flag = ("SAME LINE" if p["same_cell_line"]
                    else ("same platform" if p["same_platform"] else ""))
            print(f"    {p['pair'][0]:14s} {p['pair'][1]:14s} resid {p['resid']:+.5f} "
                  f"CI [{p['ci'][0]:+.5f},{p['ci'][1]:+.5f}] z {p['z']:+6.2f}  {flag}")

    (WORK / "factor2_204.json").write_text(json.dumps(out, indent=2))
    print("\nwrote work/factor2_204.json")


if __name__ == "__main__":
    main()

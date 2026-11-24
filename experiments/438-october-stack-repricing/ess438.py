"""438 diagnostic: does a step's PLOLO gain live on common-essential targets (238's lesson)?

For a panel, a step (arm vs ref) and each cancer fold: the relative aggregate-cosine change and the board-emulated pds
change, on all targets and with the panel's K562-essential targets removed (the list 238 used:
data/donors/k562_essential/perts.txt, n 2,021). A panel drawn the way Arc draws (essential-free) has none; 199's
realistic draw has 34.

    uv run python ess438.py <panel> <arm> <ref> [level]
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
ESS = set(l.strip() for l in (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines() if l.strip())


def cos_sub(x, a, d, keep):
    return float(np.mean([x[c, keep].sum() / np.sqrt(a[c, keep].sum() * d[keep].sum()) for c in range(3)]))


def main(panel, arm, ref, level=False):
    L = ":L" if level else ""
    out = {}
    for f in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        Z = np.load(WORK / panel / f"arrays438-{f}-held.npz")
        tg = [str(t) for t in np.load(WORK / panel / f"fold-{f}-held.npz", allow_pickle=True)["targets"]]
        ess = np.asarray([t in ESS for t in tg])
        d = Z["d"]
        g = lambda n, k: Z[f"{n}{L}@R16__{k}"]
        xa, aa, pa, xr, ar, pr = g(arm, "x"), g(arm, "a"), g(arm, "pdsl"), g(ref, "x"), g(ref, "a"), g(ref, "pdsl")
        allk = np.ones(len(tg), bool)
        r_all = cos_sub(xa, aa, d, allk) / cos_sub(xr, ar, d, allk) - 1
        r_ne = cos_sub(xa, aa, d, ~ess) / cos_sub(xr, ar, d, ~ess) - 1
        dp_all = float(pa.mean() - pr.mean())
        dp_ne = float(pa[:, ~ess].mean() - pr[:, ~ess].mean())
        dx = (xa - xr).sum(0)
        share = float(dx[ess].sum() / dx.sum()) if ess.any() else 0.0
        top = np.argsort(-dx)[:3]                                  # 236's leave-top-3-out falsifier
        keep = np.ones(len(tg), bool)
        keep[top] = False
        r_lt3 = cos_sub(xa, aa, d, keep) / cos_sub(xr, ar, d, keep) - 1
        out[f] = dict(n=len(tg), n_ess=int(ess.sum()), rel_all=r_all, rel_nonessential=r_ne, dpds_all=dp_all,
                      dpds_nonessential=dp_ne, share_of_dx_on_essential=share, rel_minus_top3=r_lt3,
                      top3=[tg[i] for i in top], top3_share=float(dx[top].sum() / dx.sum()))
        print(f"{panel} {f}: n {len(tg)} (essential {int(ess.sum())}) | rel cos all {r_all:+.4f} non-ess {r_ne:+.4f} "
              f"minus top-3 {r_lt3:+.4f} (top-3 {','.join(tg[i] for i in top)}: {100 * dx[top].sum() / dx.sum():.0f}% of sum dx) | "
              f"dpds_bd all {dp_all:+.5f} non-ess {dp_ne:+.5f} | essential share of sum dx {share:.2f}")
    return out


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], len(sys.argv) > 4)

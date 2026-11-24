"""204 Phase 1c — is `a_i` even identified?  A one-factor (rank-1) test on the cross-line
Gram, with a positive control.

Phase 1 got two answers for the same quantity: the self-gene residual says
K562 > HCT116 > HEK293T, the 3-donor shared-target Gram says the opposite.  Both cannot be
right, and before averaging them one has to ask whether the object they are estimating
exists at all.

The mandate's model is `d_i = a_i s + e_i` with the `e_i` INDEPENDENT ACROSS LINES.  That
model makes a hard, falsifiable prediction on more than three lines: the off-diagonal
cross-line correlation matrix is **rank one**,

    rho_ij = c_i c_j,     c_i = a_i |s| / ||d_i||   (donor i's reliability)

so every tetrad vanishes: `rho_ij rho_kl = rho_ik rho_jl`.  With N lines there are
N(N-1)/2 equations and N unknowns, so from N = 4 the model is over-identified and can be
rejected.  If it is rejected, no scalar `a_i` exists and the disagreement between routes
A and B is not a precision problem -- it is the model failing.

POSITIVE CONTROL, built in: `k562_essential` is the SAME CELL LINE as `k562_gwps`
(different library).  Its line-specific component is shared with `k562_gwps` by
construction, so that one pair MUST show a large positive residual against the rank-1 fit.
A test that does not flag it is not sensitive enough to flag anything.
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

from prep204 import DONORS, prepared, shared_targets  # noqa: E402
from plolo import DONOR_PATHS, PANEL_SET  # noqa: E402

LINES = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    "kolf2": ROOT / "data/external/kolf2",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "h1_2025": ROOT / "data/donors/h1_2025",
}
# which lines share a cell line / a platform, for reading the residual structure
CELL_LINE = dict(k562_gwps="K562", k562_essential="K562", orion_hct116="HCT116",
                 orion_hek293t="HEK293T", kolf2="KOLF2", rpe1="RPE1",
                 nadig_jurkat="JURKAT", h1_2025="H1")
PLATFORM = dict(k562_gwps="replogle", k562_essential="replogle", rpe1="replogle",
                orion_hct116="xatlas", orion_hek293t="xatlas",
                nadig_jurkat="nadig", kolf2="kolf2", h1_2025="vcc2025")

BLOCKS = {
    "panel_block": ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"),
    "replogle_block": ("k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
                       "nadig_jurkat", "k562_essential"),
    "five_no_sameline": ("k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
                         "nadig_jurkat"),
}


def fit_rank1(R: np.ndarray, iters: int = 500) -> np.ndarray:
    """Least squares c minimising sum_{i<j} (R_ij - c_i c_j)^2 (the diagonal is ignored)."""
    n = R.shape[0]
    off = ~np.eye(n, dtype=bool)
    c = np.sqrt(np.maximum(np.abs(R[off]).mean(), 1e-6)) * np.ones(n)
    for _ in range(iters):
        for i in range(n):
            j = [k for k in range(n) if k != i]
            num = float(R[i, j] @ c[j])
            den = float(c[j] @ c[j])
            c[i] = num / den if den > 0 else 0.0
    return c


def main() -> None:
    tg = shared_targets()
    P = {}
    for name, path in LINES.items():
        DONOR_PATHS[name] = path
        p = prepared(name, tg, "shared") if name in DONORS else prepared(
            name, tg, "allshared")
        P[name] = (np.asarray(p["effects"], dtype=np.float32),
                   np.asarray(p["target_available"], bool))
        print(f"[factor] {name}: {int(P[name][1].sum())} of {len(tg)} targets", flush=True)

    support = {n: (np.abs(P[n][0]).sum(0) > 0) for n in LINES}
    out: dict = {"blocks": {}}

    for bname, block in BLOCKS.items():
        rows = np.flatnonzero(np.logical_and.reduce([P[n][1] for n in block]))
        cols = np.flatnonzero(np.logical_and.reduce([support[n] for n in block]))
        if len(rows) < 30:
            out["blocks"][bname] = dict(skipped=f"only {len(rows)} common targets")
            continue
        M = {n: P[n][0][np.ix_(rows, cols)].astype(np.float64) for n in block}
        k = len(block)
        G = np.zeros((k, k))
        for i in range(k):
            for j in range(i, k):
                G[i, j] = G[j, i] = float((M[block[i]] * M[block[j]]).sum())
        d = np.sqrt(np.diag(G))
        R = G / np.outer(d, d)
        c = fit_rank1(R)
        fit = np.outer(c, c)
        off = ~np.eye(k, dtype=bool)
        resid = R - fit
        pairs = []
        for i, j in itertools.combinations(range(k), 2):
            pairs.append(dict(pair=[block[i], block[j]], rho=float(R[i, j]),
                              fitted=float(fit[i, j]), resid=float(resid[i, j]),
                              rel_resid=float(resid[i, j] / R[i, j]),
                              same_cell_line=CELL_LINE[block[i]] == CELL_LINE[block[j]],
                              same_platform=PLATFORM[block[i]] == PLATFORM[block[j]]))
        rms = float(np.sqrt((resid[off] ** 2).mean()))
        # tetrads on every 4-subset
        tet = []
        for quad in itertools.combinations(range(k), 4):
            a, b, cc, dd = quad
            for (i, j, m, n) in ((a, b, cc, dd), (a, cc, b, dd), (a, dd, b, cc)):
                tet.append(dict(quad=[block[q] for q in (i, j, m, n)],
                                value=float(R[i, j] * R[m, n] - R[i, m] * R[j, n])))
        out["blocks"][bname] = dict(
            lines=list(block), n_targets=int(len(rows)), n_genes=int(len(cols)),
            n_panel=int(sum(1 for r in rows if tg[r] in PANEL_SET)),
            rho=R.tolist(), c=dict(zip(block, c.tolist())),
            rank1_rms_resid=rms,
            mean_abs_rho=float(np.abs(R[off]).mean()),
            rel_rms=float(rms / np.abs(R[off]).mean()),
            pairs=pairs, tetrads=tet,
            # the production trio's attenuation implied by the fitted reliabilities
            a_rel_from_c={d_: float(c[block.index(d_)] / c[block.index(DONORS[0])])
                          for d_ in DONORS if d_ in block})
        print(f"\n[{bname}] n_targets={len(rows)} n_genes={len(cols)} "
              f"rank1 RMS resid {rms:.5f} on mean|rho| {np.abs(R[off]).mean():.5f} "
              f"({100*rms/np.abs(R[off]).mean():.1f}%)", flush=True)
        print("   c = " + "  ".join(f"{n}:{c[i]:.4f}" for i, n in enumerate(block)))
        for p in sorted(pairs, key=lambda x: -abs(x["resid"])):
            flag = ("SAME LINE" if p["same_cell_line"] else
                    ("same platform" if p["same_platform"] else ""))
            print(f"   {p['pair'][0]:15s} {p['pair'][1]:15s} rho {p['rho']:+.5f} "
                  f"fit {p['fitted']:+.5f} resid {p['resid']:+.5f} "
                  f"({100*p['rel_resid']:+6.1f}%)  {flag}")

    (WORK / "factor204.json").write_text(json.dumps(out, indent=2))
    print("\nwrote work/factor204.json")


if __name__ == "__main__":
    main()

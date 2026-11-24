"""Amplitude response of every COMPUTABLE term, on the real challenge panel.

The pipeline's `VCC_ALPHA` multiplies the effect in log space before `expm1`
(`make_submission_pm_cd4.py:348`, `multi_donor.py:169`), and the quasi emitter repairs each
gene's panel total exactly, so the metric's profile is

    b_hat_p(alpha) = log1p( TS * expm1(b_ctrl + alpha*e_p) / Z(alpha) )   ~=  b_ctrl + alpha*e_p

To vary alpha from a SHIPPED submission we invert that: read `e_p = b_hat_p - b_ctrl` off the
submitted group sums, rescale it by `s = alpha_new / alpha_old`, and re-render. `Z` drift makes
this exact only to the extent `sum_g expm1(b_ctrl + s*e)` tracks; it is checked at `s = 1`,
where the reconstruction must reproduce the shipped sums.

    A(s), B(s)   exact from the rendered sums
    C_hat_p(s)   needs the CELLS -- re-emitted through the production `emit_counts` for a
                 subsample and applied as a ratio to the shipped C_hat_p
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.sparse import csr_matrix

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights  # noqa: E402
from cell_eval2.moments import jackknife_correction  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402
from vcc2026.donor_transfer import emit_counts  # noqa: E402

TS = 50_000.0
CTX = ("A", "B", "C")
NCELL = 400


def render(bctrl, dh, S, s, floor=False):
    """Group count sums for the profile `b_ctrl + s*dh`.

    `floor=False` by default and that is deliberate. `dh` is read off ALREADY-FLOORED shipped
    sums, so at `s = 1` the unfloored render is an exact identity (`expm1(b_hat) = TS*P/S`,
    renormalised, times `S`, is `P`) and every term reproduces the shipped one to machine
    precision. Flooring a second time adds the emitter's truncation bias twice -- measured
    +1.8 to +2.1 % on `A` at `s = 1`, i.e. it would corrupt the very anchor the sweep is read
    against. The shipped submission already carries ONE flooring; the counterfactual carries
    the same one through `dh`.
    """
    prof = np.expm1(bctrl[None, :] + s * dh)
    np.clip(prof, 0.0, None, out=prof)
    prof /= prof.sum(axis=1, keepdims=True)
    if not floor:
        return S[:, None] * prof
    L = np.floor(S / NCELL)
    return np.floor(NCELL * L[:, None] * prof).astype(np.int64)


def load(stats_path):
    z = np.load(stats_path, allow_pickle=True)
    ctrl = np.load(f"{ROOT}/experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    perts = z["perts"].astype(str)
    genes = z["genes"].astype(str)
    gpos = {g: i for i, g in enumerate(genes)}
    out = {}
    for c in CTX:
        sel = np.flatnonzero(np.char.startswith(perts, f"{c}|"))
        tgt = np.asarray([p.split("|", 1)[1] for p in perts[sel]])
        bctrl = bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
        pm = bulk_lognorm_means(z["sums"][sel].astype(np.float64), TS)
        out[c] = dict(sel=sel, tgt=tgt, bctrl=bctrl, dh=pm - bctrl[None, :],
                      S=z["sums"][sel].sum(1).astype(np.float64), jk=z["jk"][sel],
                      Cctrl=float(ctrl[f"jk_{c}"]),
                      cols={k: gpos[t] for k, t in enumerate(tgt) if t in gpos},
                      G=genes.size)
    return out


def terms(d, sums):
    G, cols = d["G"], d["cols"]
    w = _row_weights(len(d["tgt"]), G, cols)
    pm = bulk_lognorm_means(np.asarray(sums, dtype=np.float64), TS)
    dh = pm - d["bctrl"][None, :]
    a = np.einsum("ij,ij->i", dh, dh)
    for k, col in cols.items():
        a[k] -= dh[k, col] ** 2
    return float((a * w).sum()), _across_pert_budget(pm, cols, w), w


def chat_ratio(d, s, n_sub=24, spread=0.25, phi=0.0, seed=0):
    """`mean C_hat_p(s) / mean C_hat_p(1)` over a subsample, through the real emitter."""
    rng = np.random.default_rng(seed)
    idx = np.linspace(0, len(d["tgt"]) - 1, n_sub).astype(int)
    got = {}
    for tag, ss in (("one", 1.0), ("s", s)):
        sums = render(d["bctrl"], d["dh"][idx], d["S"][idx], ss)
        vals = []
        for j in range(idx.size):
            prof = sums[j] / max(sums[j].sum(), 1)
            libs = np.full(NCELL, d["S"][idx[j]] / NCELL)
            v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(seed + j),
                                     mode="quasi", spread=spread, overdispersion=phi)
            X = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                           shape=(NCELL, d["G"]))
            vals.append(jackknife_correction(X, np.zeros(NCELL, np.intp), 1, TS, chunk=64)[0])
        got[tag] = np.asarray(vals)
    return got["s"], got["one"], idx


if __name__ == "__main__":
    stats = sys.argv[1]
    grid = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [
        0.2, 0.3, 0.37050, 0.45, 0.55, 0.65, 0.8, 1.0, 1.2]
    D = load(stats)
    print(f"{'s':>8}{'A':>12}{'B':>12}{'A-B':>12}   (mean over the three contexts)")
    rec = {}
    for s in grid:
        As, Bs = [], []
        for c in CTX:
            d = D[c]
            sums = render(d["bctrl"], d["dh"], d["S"], s)
            A, B, w = terms(d, sums)
            As.append(A); Bs.append(B)
        rec[s] = (float(np.mean(As)), float(np.mean(Bs)), As, Bs)
        print(f"{s:>8.4f}{np.mean(As):>12.6f}{np.mean(Bs):>12.6f}{np.mean(As)-np.mean(Bs):>12.6f}")
    np.savez(f"{ROOT}/experiments/096-mse-floor/work/alpha_terms_"
             f"{stats.split('/')[-1].split('-')[0]}.npz",
             grid=np.asarray(grid),
             A=np.asarray([rec[s][2] for s in grid]),
             B=np.asarray([rec[s][3] for s in grid]))

"""188 — the same maps, but on a SCALE-FREE x-axis, so they can be carried to another surface.

`budget188.py` fits `m` against the raw `|p|` of the PLOLO effect matrix.  That is fine inside
PLOLO, where every fold's effect matrix has the same median row norm (4.90 / 4.90 / 4.95 /
4.14) on the same 7,022-gene challenge gate.  It is NOT fine on h1-test, whose artifact is
gated at `--min-cpm 5.0` over 10,783 of 18,080 genes: the same map would be read on a different
part of its own domain and the comparison would be meaningless.

So both axes are made dimensionless per row before pooling:

    u = |p_g| / rms_gated(p)        x-axis, dimensionless, ~1 at a typical gene
    v = r_g*sign(p_g) / rms_gated(r)   y-axis, in units of the row's own truth scale

`m(u)` is then a pure shape and the row-L2 restore removes the y scale entirely.  The maps are
fitted leave-one-line-out for the PLOLO read-out and POOLED over all four donor lines for the
h1-test read-out -- h1-test's truth is H1 hESC, which is not one of the four folds, so the
pooled map cannot see its answer.

Also dumps `maps_norm.json` for `realize188.py`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "experiments/169-sharpening-six-member"))
from mapfit import Map, fit_map, fit_map_budget, to_lfc                    # noqa: E402
from budget188 import FOLDS, N_DE, PDS_SPAN, deliver, load, nmae_star, pds_of  # noqa: E402
from sharp import sharpen                                                  # noqa: E402


def rms(A, mask):
    return np.sqrt(np.maximum((A ** 2 * mask[None, :]).sum(1) / max(mask.sum(), 1), 1e-300))


def apply_norm(P, m, gate, mode="row"):
    """`p -> sign(p)*s*m(|p|/s)` with `s` the row's own gated rms, then row L2 restored."""
    A = np.atleast_2d(np.asarray(P, dtype=np.float64))
    g = np.asarray(gate, dtype=bool)
    G = np.where(g[None, :], A, 0.0)
    s = rms(G, g)[:, None]
    out = np.where(g[None, :], np.sign(G) * m(np.abs(G) / s), 0.0)
    if mode == "row":
        n0 = np.linalg.norm(G, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        out = out * (n0 / np.maximum(n1, 1e-300))[:, None]
    return out[0] if P.ndim == 1 else out


def main() -> None:
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    bA, gA = z["bA"], z["gA"]
    F = {f: load(f) for f in FOLDS}
    cols, D = {}, {}
    for f in FOLDS:
        mk = np.zeros(cg.size, bool)
        mk[F[f]["scored"]] = True
        cols[f] = np.flatnonzero(mk & cg)
        c = cols[f]
        P, R = F[f]["pred"][:, c], F[f]["truth"][:, c]
        one = np.ones(c.size, bool)
        sP, sR = rms(P, one)[:, None], rms(R, one)[:, None]
        live = (P != 0) & (R != 0)
        de = np.zeros(R.shape, bool)
        k = min(N_DE, R.shape[1])
        np.put_along_axis(de, np.argpartition(-np.abs(R), k - 1, axis=1)[:, :k], True, 1)
        D[f] = ((np.abs(P) / sP)[live], (R * np.sign(P) / sR)[live], de[live])
        print(f"  {f:<14} pred rms/row {float(np.median(sP)):.5f}  truth rms/row "
              f"{float(np.median(sR)):.5f}  {c.size} cols")

    maps = {}
    for held in list(FOLDS) + ["POOLED"]:
        src = [f for f in FOLDS if f != held]
        u = np.concatenate([D[f][0] for f in src])
        v = np.concatenate([D[f][1] for f in src])
        d = np.concatenate([D[f][2] for f in src])
        maps[("median", held)] = fit_map(u, v)
        maps[("budget", held)], mu = fit_map_budget(u, v)
        maps[("de", held)] = fit_map(u[d], v[d])
        maps[("debudget", held)], mu2 = fit_map_budget(u[d], v[d])
        print(f"  [{held}] budget mu = {mu:.4f}   DE-budget mu = {mu2:.4f}")

    with open(os.path.join(WORK, "maps_norm.json"), "w") as fh:
        json.dump({f"{k}|{h}": maps[(k, h)].to_dict() for (k, h) in maps}, fh, indent=1)

    print(f"\n{'=' * 92}\nPOOLED normalised maps (u = |p| / rms, v in truth-rms units)")
    M = {k: maps[(k, "POOLED")] for k in ("median", "budget", "de", "debudget")}
    print(f"  {'pctile':>7}{'u':>9}{'m_median':>10}{'m/u':>8}{'m_budget':>10}{'m/u':>8}"
          f"{'m_de':>9}{'m_debud':>10}")
    xs = M["median"].xs
    for i in list(range(0, len(xs), 4)) + [len(xs) - 1]:
        u = xs[i]
        print(f"  {100.0*(i+0.5)/len(xs):>7.1f}{u:>9.4f}{M['median'](u):>10.5f}"
              f"{M['median'](u)/max(u,1e-9):>8.4f}{M['budget'](u):>10.5f}"
              f"{M['budget'](u)/max(u,1e-9):>8.4f}{M['de'](u):>9.5f}"
              f"{M['debudget'](u):>10.5f}")

    print(f"\n{'=' * 92}\nHELD-OUT, challenge frame A, delivered — normalised maps against the "
          f"power family")
    print(f"  {'arm':<12}{'d pds raw':>11}{'d pds scaled':>14}{'d nmae*':>10}"
          f"{'DE L1 over k=1':>16}{'rA':>8}")
    rows = {}
    for name in ("identity", "power1.15", "median", "budget", "de", "debudget"):
        dp, dn, ov, ra = [], [], [], []
        for held in FOLDS:
            P0, c = F[held]["pred"], cols[held]
            LR = to_lfc(F[held]["truth"], F[held]["b"], gA)
            de = np.zeros((P0.shape[0], c.size), bool)
            k = min(N_DE, c.size)
            np.put_along_axis(de, np.argpartition(-np.abs(LR[:, c]), k - 1, axis=1)[:, :k],
                              True, 1)
            if name == "identity":
                Q = np.where(gA[None, :], P0, 0.0)
            elif name == "power1.15":
                Q = sharpen(P0, 1.15, gA, mode="row")
            else:
                Q = apply_norm(P0, maps[(name, held)], gA)
            Dl = deliver(Q, bA, gA)
            LP = to_lfc(Dl, bA, gA)
            base = deliver(np.where(gA[None, :], P0, 0.0), bA, gA)
            LB = to_lfc(base, bA, gA)
            dp.append(pds_of(Dl[:, F[held]["scored"]], F[held]["truth"][:, F[held]["scored"]])
                      - pds_of(base[:, F[held]["scored"]],
                               F[held]["truth"][:, F[held]["scored"]]))
            dn.append(nmae_star(LP[:, c], LR[:, c], de) - nmae_star(LB[:, c], LR[:, c], de))
            a = np.abs(LP[:, c])[de]
            ov.append(float(a[a > np.abs(LR[:, c])[de]].sum() / max(a.sum(), 1e-12)))
            ra.append(float(np.median(np.linalg.norm(Q, axis=1)) /
                            np.median(np.linalg.norm(np.where(gA[None, :], P0, 0.0), axis=1))))
        rows[name] = dict(dpds=float(np.mean(dp)), dpds_scaled=float(np.mean(dp)) / PDS_SPAN,
                          dnmae=float(np.mean(dn)), over=float(np.mean(ov)),
                          rA=float(np.mean(ra)), per_fold=dp)
        print(f"  {name:<12}{np.mean(dp):>+11.6f}{np.mean(dp)/PDS_SPAN:>+14.6f}"
              f"{np.mean(dn):>+10.5f}{np.mean(ov):>16.4f}{np.mean(ra):>8.5f}")
    with open(os.path.join(WORK, "normfit188.json"), "w") as fh:
        json.dump(rows, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/maps_norm.json and normfit188.json")


if __name__ == "__main__":
    sys.exit(main())

"""518 — the strength key on H1 (third surface, exact mse legs; h1arms518's emulated recipe, each arm re-levelled).

z_p = centred percentile rank of the robust response strength of each H1 target in the three production sources
(K562-GWPS, Orion HCT116, Orion HEK293T; strength518's statistic, the per-gene MAD scale taken over the 300 A/B/C panel
rows of each bank, applied to the H1 targets' rows). No CD4 / VIPerturb / KOLF2 term: they carry the A/B/C panel only.
Arms: E0 rows x exp(theta z), theta in {0.5, 1, 1.25, 1.5}; control = z permuted. On the H1 recipe the key is 410's
(clip (0.3, 3)), so the price is the strength key over 410's key, not over the LODO key.

    uv run python experiments/518-bulk-mse-gap/h1skey518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import json  # noqa: E402

import numpy as np  # noqa: E402

import h1arms518 as H  # noqa: E402
from lofo518 import crank  # noqa: E402
from render518 import render, sharpen_on  # noqa: E402

WORK = HERE / "work"
PANEL = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
PSET = set(PANEL)
BANKS = {"k562_gwps": ROOT / "data/donors/k562_gwps", "orion_hct116": ROOT / "data/donors/orion_hct116",
         "orion_hek293t": ROOT / "data/donors/orion_hek293t"}


def strength_for(path, targets, thr=4.0):
    perts = [l.strip() for l in (path / "perts.txt").read_text().splitlines() if l.strip()]
    genes = [l.strip() for l in (path / "genes.txt").read_text().splitlines() if l.strip()]
    b = np.load(path / "b_ctrl.npy").astype(np.float64)
    cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
    keep = (cpm > 10.0) & np.asarray([g not in PSET and g not in set(targets) for g in genes])
    pi = {p: i for i, p in enumerate(perts)}
    delta = np.load(path / "delta.npy", mmap_mode="r")
    prow = sorted(pi[t] for t in PANEL if t in pi)
    P = np.asarray(delta[np.asarray(prow)], np.float64)[:, keep]
    sc = 1.4826 * np.median(np.abs(P - np.median(P, 0)), 0)
    ok = sc > 0
    cov = [t for t in targets if t in pi]
    rows = [pi[t] for t in cov]
    Rr = np.asarray(delta[np.asarray(sorted(rows))], np.float64)[np.argsort(np.argsort(rows))][:, keep][:, ok] / sc[ok]
    nz = np.median(np.abs(Rr), 1) / 0.6745
    cnt = (np.abs(Rr) > thr * nz[:, None]).sum(1)
    o = dict(zip(cov, np.log1p(cnt).tolist()))
    return np.asarray([o.get(t, np.nan) for t in targets])


def main():
    D = H.load()
    M = H.MSE(D)
    b = D["b"]
    E0 = np.load(WORK / "h1arms518_recipe.npz")["E0"]
    gpos = {g: i for i, g in enumerate(D["genes"])}
    own = np.zeros(b.size, bool)
    own[[gpos[t] for t in D["targets"] if t in gpos]] = True
    cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
    g10 = (cpm > 10) & (b > 0)
    cols = []
    for k, p in BANKS.items():
        v = strength_for(p, D["targets"])
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
        print(f"{k}: covers {int(np.isfinite(v).sum())}/{len(D['targets'])} H1 targets")
    z = np.nanmean(np.vstack(cols), 0)
    z = np.where(np.isfinite(z), z, 0.0)

    def lev(E):
        Eq = sharpen_on(E, 1.75, g10 & ~own)
        best = None
        for f in np.round(np.arange(0.7, 1.41, 0.01), 2):
            h = render(f * Eq, b, "add", g10)
            r_, rho = M.raw(h)
            if best is None or r_ < best[0]:
                best = (r_, rho, float(f), h)
        return best
    base = lev(E0)
    pc = np.flatnonzero(own)
    pb = H.pds(base[3], D["dl"], pc)
    out = dict(base=dict(raw=base[0], pds=pb), arms={})
    rng = np.random.default_rng(518)
    zp = z[rng.permutation(z.size)]
    for nm, th, zz in (("th0.5", 0.5, z), ("th1", 1.0, z), ("th1.25", 1.25, z), ("th1.5", 1.5, z), ("CTRL th1.25 permuted", 1.25, zp)):
        r = lev(E0 * np.exp(th * zz)[:, None])
        pa = H.pds(r[3], D["dl"], pc)
        out["arms"][nm] = dict(raw=r[0], draw=r[0] - base[0], pds=pa, dpds=pa - pb, level=r[2])
        print(f"{nm:22s} raw {r[0]:.5f} ({r[0] - base[0]:+.5f})  pds {pa:.5f} ({pa - pb:+.5f})  level {r[2]:.2f}")
    (WORK / "h1skey518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()


def null(nperm=20):
    """20 permuted-strength controls at theta 1 and 1.25: the null distribution of a per-target reweighting on H1."""
    D = H.load()
    M = H.MSE(D)
    b = D["b"]
    E0 = np.load(WORK / "h1arms518_recipe.npz")["E0"]
    gpos = {g: i for i, g in enumerate(D["genes"])}
    own = np.zeros(b.size, bool)
    own[[gpos[t] for t in D["targets"] if t in gpos]] = True
    cpm = 1e6 * np.expm1(b) / np.expm1(b).sum()
    g10 = (cpm > 10) & (b > 0)
    cols = []
    for k, p in BANKS.items():
        v = strength_for(p, D["targets"])
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    z = np.nanmean(np.vstack(cols), 0)
    z = np.where(np.isfinite(z), z, 0.0)

    def lev(E):
        Eq = sharpen_on(E, 1.75, g10 & ~own)
        best = None
        for f in np.round(np.arange(0.6, 1.41, 0.01), 2):
            r_, _ = M.raw(render(f * Eq, b, "add", g10))
            if best is None or r_ < best:
                best = r_
        return best
    base = lev(E0)
    res = {}
    for th in (1.0, 1.25):
        arm = lev(E0 * np.exp(th * z)[:, None]) - base
        nulls = [lev(E0 * np.exp(th * z[np.random.default_rng(1000 + s).permutation(z.size)])[:, None]) - base for s in range(nperm)]
        p = float(np.mean(np.asarray(nulls) <= arm))
        res[f"th{th:g}"] = dict(arm=arm, null_mean=float(np.mean(nulls)), null_sd=float(np.std(nulls)), p_null_le_arm=p, nulls=nulls)
        print(f"theta {th:g}: strength key d raw {arm:+.5f}; permuted null mean {np.mean(nulls):+.5f} sd {np.std(nulls):.5f}; "
              f"P(null <= arm) {p:.2f}")
    js = json.loads((WORK / "h1skey518.json").read_text())
    js["null"] = res
    (WORK / "h1skey518.json").write_text(json.dumps(js, indent=1, default=float))

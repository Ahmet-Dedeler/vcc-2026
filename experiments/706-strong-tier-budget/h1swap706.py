#!/usr/bin/env python3
"""706c: the strong-tier direction swap on h1-test, the Arc-lab / Flex-chemistry surface (A/B/C's lab and chemistry).

dirswap706 found, on 522's 9 October-like folds, that blending the source donors' raw unit-mean direction into the
strong-z tier raises the aggregate cosine by +5 … +31 %, with or without the across-target mean. But every fold's truth is
a public 10x 3' line (K562 GWPS, Orion HCT116 / HEK293T): the same family as the raw donor rows. Our stack also carries
Arc's H1 and the CD4 banks. If the fold gain is lab/chemistry affinity, it vanishes on an Arc-lab truth.

Surface: h1-test (85 targets, the 2025 Arc H1 Flex data; truth = pseudobulk log1p(5e4) deltas vs its controls, as
eval700_h1 / overread702 build it). Base: 435's production key rows (lib521 state 0, the H1 harness's bulk direction; the
h1_2025 donor never enters). Raw direction: the three production donors' unit-mean U (K562 GWPS, Orion HCT116, Orion
HEK293T) on the state's gene axis, panel genes excluded.
Tiers: tertiles of the donors' mean robust strength (here: mean row norm over the donors carrying the target, a stand-in
for z_src on a panel 522 never ranked). Blend as dirswap706: norm-preserving, β 0.5 and 1, on the top tier / top two /
all, with and without centring U. ew-cos = Σ<P,T> / sqrt(Σ|P|² Σ|T|²) over targets (the cos_agg analogue, no render) and
mini-pds. Target bootstrap 4,000 draws, seed 706.
usage: uv run python experiments/706-strong-tier-budget/h1swap706.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402

rng = np.random.default_rng(706)


def donor(nm):
    D = ROOT / "data/donors" / nm
    return (D / "genes.txt").read_text().split(), (D / "perts.txt").read_text().split(), np.load(D / "delta.npy", mmap_mode="r")


def ewcos(P, T, idx):
    return float((P[idx] * T[idx]).sum() / np.sqrt((P[idx] ** 2).sum() * (T[idx] ** 2).sum()))


def minipds(P, T):
    Pu = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    Tu = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = Pu @ Tu.T
    own = np.diag(S)
    n = len(S)
    return float(np.mean(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)))


def main():
    import anndata as ad
    S = L.load_state(0)
    targets, genes = list(S["targets"]), list(S["genes"])
    a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad")
    vp = {g: i for i, g in enumerate(a.var_names)}
    cols = np.array([vp[g] for g in genes])
    tg = a.obs["target_gene"].to_numpy()
    X = a.X.tocsr()

    def pb(m):
        P = np.asarray(X[m].sum(0)).ravel().astype(float)
        return np.log1p(5e4 * P / P.sum())[cols]
    bc = pb(tg == "non-targeting")
    T = np.stack([pb(tg == x) - bc for x in targets])
    keep = ~np.isin(np.array(genes), np.array(targets))
    B = np.asarray(S["key"], float).copy()
    B[:, ~keep] = 0.0
    T[:, ~keep] = 0.0
    U = np.zeros_like(B)
    norms = []
    for nm in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        dg, dp, dd = donor(nm)
        gpos = {g: i for i, g in enumerate(dg)}
        ppos = {p: i for i, p in enumerate(dp)}
        gc = np.array([gpos.get(g, -1) for g in genes])
        ok = gc >= 0
        nr = np.full(len(targets), np.nan)
        for i, t in enumerate(targets):
            if t in ppos:
                r = np.zeros(len(genes))
                r[ok] = np.asarray(dd[ppos[t]], float)[gc[ok]]
                r[~keep] = 0.0
                n_ = np.linalg.norm(r)
                if n_ > 0:
                    U[i] += r / n_
                    nr[i] = n_
        norms.append(nr)
    with np.errstate(all="ignore"):
        strength = np.nanmean(np.vstack(norms), 0)
    have = np.linalg.norm(U, axis=1) > 0
    st = np.where(np.isfinite(strength), strength, -np.inf)
    r = np.argsort(np.argsort(st[have], kind="stable"), kind="stable")
    tier = np.full(len(targets), -1)
    tier[np.flatnonzero(have)] = np.minimum((r * 3) // have.sum(), 2)
    idx = np.arange(len(targets))
    base_ew, base_pds = ewcos(B, T, idx), minipds(B, T)
    out = {"n_targets": len(targets), "n_have": int(have.sum()), "base": {"ew": base_ew, "minipds": base_pds},
           "shared_truth": float(len(T) * (T.mean(0) @ T.mean(0)) / (T * T).sum())}
    for centre in (False, True):
        Uc = U.copy()
        if centre:
            Uc[have] -= Uc[have].mean(0, keepdims=True)
            Uc[:, ~keep] = 0.0
        for sel_name, sel in (("top", tier == 2), ("top2", tier >= 1), ("all", have)):
            for beta in (0.5, 1.0):
                P = B.copy()
                for i in np.flatnonzero(sel):
                    nb, nu = np.linalg.norm(B[i]), np.linalg.norm(Uc[i])
                    if nb == 0 or nu == 0:
                        continue
                    m = (1 - beta) * B[i] / nb + beta * Uc[i] / nu
                    P[i] = nb * m / np.linalg.norm(m)
                d = []
                for _ in range(4000):
                    s = rng.choice(idx, len(idx))
                    d.append(ewcos(P, T, s) / ewcos(B, T, s) - 1)
                key = f"{'c_' if centre else ''}{sel_name}_b{beta:g}"
                out[key] = {"dew_pct": 100 * (ewcos(P, T, idx) / base_ew - 1),
                            "ci95": [100 * float(np.quantile(d, .025)), 100 * float(np.quantile(d, .975))],
                            "dminipds": minipds(P, T) - base_pds, "n_sel": int(sel.sum())}
    print(json.dumps(out, indent=1))
    (HERE / "work/h1swap706.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

"""276 — POST-HOC robustness: the same statistic on the WIDEST gene set each pair shares.

The frozen construction uses one common 2,994-gene axis for every pair, which is the cleanest
comparability property but throws away 84 % of the challenge's genes and costs 0.097 of the
donor-donor baseline's rho (RESULTS.md §2.3).  Here each pair is scored on the intersection of
its two sources' own gene symbols with the challenge's 18,533, after centring each source once
on its own across-target mean over its own full official-gene set.  Everything else — the pair
rules, the leave-the-held-line-out sets, the outcomes, the bootstrap — is unchanged.

Also: a split-panel confirmation (the 300 targets cut in two by a deterministic hash) so the
prediction is reported on two disjoint halves of the panel rather than once.

    PYTHONDONTWRITEBYTECODE=1 uv run python wide276.py > work/wide276.md
"""
from __future__ import annotations

import hashlib
import resource

import numpy as np
from scipy.stats import rankdata, spearmanr

import lib276 as L

SEED = 20260906
B = 2000
E262 = L.ROOT / "experiments/262-where-pds-is-lost/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
HELD = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
DD = {"k562_gwps": ("orion_hct116", "orion_hek293t"),
      "orion_hct116": ("k562_gwps", "orion_hek293t"),
      "orion_hek293t": ("k562_gwps", "orion_hct116")}
IN_PRED = ("k562_gwps", "perturbmatch", "orion_hct116", "orion_hek293t",
           "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "gse295214_h1", "h1_2025")


def boot(x, y, b=B, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    rng = np.random.default_rng(seed)
    d = np.array([spearmanr(x[s := rng.integers(0, x.size, x.size)], y[s]).statistic
                  for _ in range(b)])
    return (float(spearmanr(x, y).statistic),
            [float(v) for v in np.percentile(d, [2.5, 97.5])], int(x.size))


def paired(x1, x2, y, b=B, seed=SEED):
    ok = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(y)
    x1, x2, y = np.asarray(x1)[ok], np.asarray(x2)[ok], np.asarray(y)[ok]
    d0 = spearmanr(x1, y).statistic - spearmanr(x2, y).statistic
    rng = np.random.default_rng(seed)
    d = np.array([spearmanr(x1[s := rng.integers(0, x1.size, x1.size)], y[s]).statistic
                  - spearmanr(x2[s], y[s]).statistic for _ in range(b)])
    return float(d0), [float(v) for v in np.percentile(d, [2.5, 97.5])], int(x1.size)


def main() -> None:
    keys = list(L.SOURCES)
    off = set(L.OFFICIAL)
    rows, avail, gsets = {}, {}, {}
    for k in keys:
        have = off & set(L.read_bank(L.SOURCES[k][0])[1])
        ax = [g for g in L.OFFICIAL if g in have]
        r, a = L.load_source(k, ax)
        rows[k] = r.astype(np.float32); avail[k] = a
        gsets[k] = {g: i for i, g in enumerate(ax)}
    pairs = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys))
             if L.SOURCES[keys[i]][1] != L.SOURCES[keys[j]][1]]
    C = np.full((len(pairs), len(L.PANEL)), np.nan)
    ngene = np.zeros(len(pairs), dtype=int)
    for p, (i, j) in enumerate(pairs):
        ki, kj = keys[i], keys[j]
        shared = [g for g in gsets[ki] if g in gsets[kj]]
        ngene[p] = len(shared)
        if not shared:
            continue
        ci = np.asarray([gsets[ki][g] for g in shared])
        cj = np.asarray([gsets[kj][g] for g in shared])
        m = avail[ki] & avail[kj]
        if not m.any():
            continue
        C[p, m] = L.cosine_rows(rows[ki][m][:, ci].astype(np.float64),
                                rows[kj][m][:, cj].astype(np.float64))
    print("## POST-HOC — the same map on each pair's own shared genes\n")
    print(f"- {len(pairs)} cross-line pairs; shared genes per pair: median {int(np.median(ngene))}, "
          f"min {int(ngene.min())}, max {int(ngene.max())} (the frozen map uses 2,994 for every pair)")
    with np.errstate(invalid="ignore"):
        c_raw_w = np.nanmean(C, axis=0)
        npair = np.isfinite(C).sum(axis=0)
    zc = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    c_raw_f = np.asarray(zc["c_raw"])
    print(f"- `C_raw` wide: mean {np.nanmean(c_raw_w):+.4f}, median {np.nanmedian(c_raw_w):+.4f}, "
          f"share <= 0 {np.nanmean(c_raw_w[np.isfinite(c_raw_w)] <= 0)*100:.1f} %; "
          f"rho with the frozen `C_raw` **{spearmanr(c_raw_w, c_raw_f).statistic:+.3f}**\n")

    def score(drop):
        d = {keys.index(k) for k in drop if k in keys}
        keep = [p for p, (i, j) in enumerate(pairs) if i not in d and j not in d]
        with np.errstate(invalid="ignore"):
            return np.nanmean(C[keep], axis=0)

    idx = {t: i for i, t in enumerate(L.PANEL)}
    cols: dict[str, list[np.ndarray]] = {k: [] for k in ("c_lolo", "c_out", "dd", "pds", "anti")}
    tnames = []
    for fold in CANCER:
        z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]
        sel = np.asarray([idx[t] for t in tg])
        drop = L.HELD_LINE_SOURCES[HELD[fold]]
        a, b_ = DD[fold]
        sh = [g for g in gsets[a] if g in gsets[b_]]
        ca = np.asarray([gsets[a][g] for g in sh]); cb = np.asarray([gsets[b_][g] for g in sh])
        m = avail[a] & avail[b_]
        ddv = np.full(len(L.PANEL), np.nan)
        ddv[m] = L.cosine_rows(rows[a][m][:, ca].astype(np.float64),
                               rows[b_][m][:, cb].astype(np.float64))
        vals = dict(c_lolo=score(drop)[sel], dd=ddv[sel],
                    c_out=score(tuple(set(drop) | set(IN_PRED)))[sel],
                    pds=np.asarray(z["cd4svd5|pds"], dtype=np.float64),
                    anti=(np.asarray(z["cov|base_cos_truth"], dtype=np.float64) <= 0).astype(float))
        for k, v in vals.items():
            r = np.full(v.size, np.nan)
            ok = np.isfinite(v)
            r[ok] = rankdata(v[ok]) / ok.sum() if k not in ("anti",) else v[ok]
            cols[k].append(r)
        tnames += tg
    P = {k: np.concatenate(v) for k, v in cols.items()}
    print("| predictor (wide axis) | rho with `pds`, pooled | frozen-axis value |")
    print("|---|---|---|")
    froz = {"c_lolo": "+0.387 [+0.326, +0.444]", "dd": "+0.307 [+0.239, +0.372]",
            "c_out": "+0.198 [+0.127, +0.263]"}
    lab = {"c_lolo": "`C_lolo` (held line out)", "dd": "baseline: donor-donor agreement",
           "c_out": "`C_out` (held line out AND every incumbent source)"}
    for k in ("c_lolo", "dd", "c_out"):
        r, ci, n = boot(P[k], P["pds"])
        print(f"| {lab[k]} | **{r:+.3f}** [{ci[0]:+.3f}, {ci[1]:+.3f}] (n {n}) | {froz[k]} |")
    d0, ci, n = paired(P["c_lolo"], P["dd"], P["pds"])
    print(f"\npaired dRho (`C_lolo` - donor-donor), wide axis: **{d0:+.3f}** "
          f"[{ci[0]:+.3f}, {ci[1]:+.3f}] on {n} rows "
          f"(frozen axis +0.066 [+0.004, +0.133])\n")

    # ---- split-panel confirmation ----------------------------------------------------------
    print("## POST-HOC — the prediction on two disjoint halves of the panel\n")
    half = np.asarray([int(hashlib.sha256(t.encode()).hexdigest(), 16) % 2 for t in tnames])
    print("| half | n rows | `C_lolo` (frozen axis) | `C_lolo` (wide axis) | donor-donor (wide) |")
    print("|---|---:|---|---|---|")
    zp = np.load(L.WORK / "predict276.npz", allow_pickle=True)
    def _rank(v):
        r = np.full(v.size, np.nan)
        ok = np.isfinite(v)
        r[ok] = rankdata(v[ok]) / ok.sum()
        return r
    froz_lolo = np.concatenate([_rank(np.asarray(zp[f"{f}|c_lolo"])) for f in CANCER])
    for h in (0, 1):
        m = half == h
        r1, c1, n1 = boot(froz_lolo[m], P["pds"][m])
        r2, c2, _ = boot(P["c_lolo"][m], P["pds"][m])
        r3, c3, _ = boot(P["dd"][m], P["pds"][m])
        print(f"| {h} | {n1} | {r1:+.3f} [{c1[0]:+.3f}, {c1[1]:+.3f}] | "
              f"{r2:+.3f} [{c2[0]:+.3f}, {c2[1]:+.3f}] | {r3:+.3f} [{c3[0]:+.3f}, {c3[1]:+.3f}] |")
    print(f"\npeak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.2f} GiB")
    np.savez_compressed(L.WORK / "wide276.npz", c_raw_wide=c_raw_w, n_pairs_wide=npair,
                        targets=np.asarray(L.PANEL, dtype=object))


if __name__ == "__main__":
    main()

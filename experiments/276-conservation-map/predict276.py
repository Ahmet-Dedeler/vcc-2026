"""276 item 4 — does conservation predict transfer?

Outcomes (pre-existing, untouched):
  * `cd4svd5` per-target `pds` and `cov|base_cos_truth` on the four PLOLO folds
    (`experiments/262-where-pds-is-lost/work/diag262-<fold>.npz`);
  * `INCUMBENT|pds_eff` on h1-test (`experiments/273-.../work/instrument273.npz`).

Predictors: `C_raw` (all sources), `C_lolo` (every source from the held line removed),
`C_rel` (reliability-corrected, four sources).  Baselines: 262's donor-donor agreement
recomputed in this residual space, the CD4 condition count, the target gene's own expression
in the challenge control, and the number of pairs (a coverage control).

    PYTHONDONTWRITEBYTECODE=1 uv run python predict276.py
"""
from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import lib276 as L

SEED = 20260906
B = 2000
E262 = L.ROOT / "experiments/262-where-pds-is-lost/work"

# fold key -> (held cell line, the two donor rows that make the fold's base)
FOLDS = {
    "k562_gwps": ("K562", ("orion_hct116", "orion_hek293t")),
    "orion_hct116": ("HCT116", ("k562_gwps", "orion_hek293t")),
    "orion_hek293t": ("HEK293T", ("k562_gwps", "orion_hct116")),
    "kolf2": ("KOLF2.1J", ("k562_gwps", "orion_hct116", "orion_hek293t")),
}
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def boot_spearman(x, y, b=B, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    if x.size < 5:
        return dict(rho=None, ci=[None, None], n=int(x.size))
    rng = np.random.default_rng(seed)
    r0 = float(spearmanr(x, y).statistic)
    d = np.empty(b)
    for i in range(b):
        s = rng.integers(0, x.size, x.size)
        d[i] = spearmanr(x[s], y[s]).statistic
    return dict(rho=r0, ci=[float(v) for v in np.percentile(d, [2.5, 97.5])], n=int(x.size))


def auc_rest_vs_group(x, group, b=B, seed=SEED):
    """P(x of a NON-group target > x of a group target); 0.5 = chance, > 0.5 = separates."""
    ok = np.isfinite(x)
    x, group = np.asarray(x)[ok], np.asarray(group)[ok].astype(bool)
    if group.sum() < 3 or (~group).sum() < 3:
        return dict(auc=None, ci=[None, None], n_group=int(group.sum()))
    a, c = x[~group], x[group]
    def _auc(a, c):
        from scipy.stats import rankdata
        r = rankdata(np.concatenate([a, c]))
        return float((r[:a.size].sum() - a.size * (a.size + 1) / 2) / (a.size * c.size))
    rng = np.random.default_rng(seed)
    d = np.empty(b)
    for i in range(b):
        d[i] = _auc(a[rng.integers(0, a.size, a.size)], c[rng.integers(0, c.size, c.size)])
    return dict(auc=_auc(a, c), ci=[float(v) for v in np.percentile(d, [2.5, 97.5])],
                n_group=int(group.sum()), n_rest=int((~group).sum()))


def partial_spearman(x, y, z, b=B, seed=SEED):
    """Spearman of x and y after linearly removing rank(z) from both (rank-partial)."""
    from scipy.stats import rankdata
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = np.asarray(x)[ok], np.asarray(y)[ok], np.asarray(z)[ok]
    if x.size < 10:
        return dict(rho=None, ci=[None, None], n=int(x.size))

    def _pr(x, y, z):
        rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
        Z = np.column_stack([np.ones_like(rz), rz])
        bx = np.linalg.lstsq(Z, rx, rcond=None)[0]
        by = np.linalg.lstsq(Z, ry, rcond=None)[0]
        ex, ey = rx - Z @ bx, ry - Z @ by
        return float(np.corrcoef(ex, ey)[0, 1])
    rng = np.random.default_rng(seed)
    d = np.empty(b)
    for i in range(b):
        s = rng.integers(0, x.size, x.size)
        d[i] = _pr(x[s], y[s], z[s])
    return dict(rho=_pr(x, y, z), ci=[float(v) for v in np.percentile(d, [2.5, 97.5])],
                n=int(x.size))


def paired_delta_spearman(x1, x2, y, b=B, seed=SEED):
    """rho(x1, y) - rho(x2, y) with a paired target bootstrap over the shared targets."""
    ok = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(y)
    x1, x2, y = np.asarray(x1)[ok], np.asarray(x2)[ok], np.asarray(y)[ok]
    if x1.size < 10:
        return dict(delta=None, ci=[None, None], n=int(x1.size))
    d0 = spearmanr(x1, y).statistic - spearmanr(x2, y).statistic
    rng = np.random.default_rng(seed)
    d = np.empty(b)
    for i in range(b):
        s = rng.integers(0, x1.size, x1.size)
        d[i] = spearmanr(x1[s], y[s]).statistic - spearmanr(x2[s], y[s]).statistic
    return dict(delta=float(d0), ci=[float(v) for v in np.percentile(d, [2.5, 97.5])],
                p_le0=float((d <= 0).mean()), n=int(x1.size))


def lolo_score(keys, Ccross, pairs, drop_keys):
    """mean cross-line pair cosine after removing every pair that touches a dropped source."""
    drop = {keys.index(k) for k in drop_keys if k in keys}
    keep = np.asarray([p for p, (i, j) in enumerate(pairs) if i not in drop and j not in drop])
    with np.errstate(invalid="ignore"):
        s = np.nanmean(Ccross[keep], axis=0)
        n = np.isfinite(Ccross[keep]).sum(axis=0)
    return s, n


def main() -> None:
    zc = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    keys = [str(k) for k in zc["keys"]]
    Ccross = np.asarray(zc["Ccross"], dtype=np.float64)
    pairs = [tuple(p) for p in np.asarray(zc["pairs_cross"])]
    c_raw = np.asarray(zc["c_raw"]); n_pairs = np.asarray(zc["n_pairs"])
    zr = np.load(L.WORK / "rel276.npz", allow_pickle=True)
    c_rel = np.asarray(zr["c_rel"])
    hkeys = [str(k) for k in zr["hkeys"]]
    Crel = np.asarray(zr["Crel"], dtype=np.float64)
    prs = [tuple(p) for p in np.asarray(zr["prs"])]
    idx = {t: i for i, t in enumerate(L.PANEL)}

    # --- donor-donor agreement baseline, recomputed in THIS residual space -------------------
    zrow = np.load(L.WORK / "rows276.npz", allow_pickle=True)
    R = np.asarray(zrow["R"]); A = np.asarray(zrow["A"])
    def pcos(a, b):
        i, j = keys.index(a), keys.index(b)
        m = A[i] & A[j]
        out = np.full(len(L.PANEL), np.nan)
        out[m] = L.cosine_rows(R[i][m].astype(np.float64), R[j][m].astype(np.float64))
        return out
    dd = {
        "k562_gwps": pcos("orion_hct116", "orion_hek293t"),
        "orion_hct116": pcos("k562_gwps", "orion_hek293t"),
        "orion_hek293t": pcos("k562_gwps", "orion_hct116"),
    }
    trio = np.nanmean(np.vstack([pcos("k562_gwps", "orion_hct116"),
                                 pcos("k562_gwps", "orion_hek293t"),
                                 pcos("orion_hct116", "orion_hek293t")]), axis=0)
    dd["kolf2"] = trio

    # --- per-source scale-free response magnitude, and the measured reliabilities ------------
    nrm: dict[str, np.ndarray] = {}
    for i, k in enumerate(keys):
        v = np.full(len(L.PANEL), np.nan)
        m = A[i]
        n = np.linalg.norm(R[i][m].astype(np.float64), axis=1)
        v[m] = n / np.median(n)
        nrm[k] = v
    rfull = {str(k): np.asarray(zr[f"r_full|{k}"]) for k in hkeys}

    # --- target expression in the challenge controls (mean over A/B/C) ----------------------
    gi = {g: i for i, g in enumerate(L.OFFICIAL)}
    bc = np.mean([np.load(E262 / f"plolo/bctrl-challenge-{c}.npy") for c in "ABC"], axis=0)
    expr = np.array([bc[gi[t]] if t in gi else np.nan for t in L.PANEL])

    ka, rows, out = [], [], {}
    per_fold = {}
    for fold, (line, _) in FOLDS.items():
        z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]
        sel = np.asarray([idx[t] for t in tg])
        pds = np.asarray(z["cd4svd5|pds"], dtype=np.float64)
        base = np.asarray(z["cov|base_cos_truth"], dtype=np.float64)
        cd4c = np.asarray(z["cov|cd4_count"], dtype=np.float64)
        anti = base <= 0
        ka.append(dict(fold=fold, pds_mean=float(pds.mean()),
                       rho_base_pds=float(spearmanr(base, pds).statistic), n=len(tg)))
        drop = L.HELD_LINE_SOURCES.get(line, ())
        c_lolo, n_lolo = lolo_score(keys, Ccross, pairs, drop)
        c_rel_lolo, _ = lolo_score(hkeys, Crel, prs, drop)
        # "is it just a measurable response?" controls, both LOLO
        keep_h = [k for k in hkeys if k not in drop]
        with np.errstate(invalid="ignore"):
            rel_mean = np.nanmean(np.vstack([rfull[k] for k in keep_h]), axis=0)
            nrm_mean = np.nanmean(np.vstack([nrm[k] for k in keys if k not in drop]), axis=0)
        per_fold[fold] = dict(
            rel_mean=rel_mean[sel], nrm_mean=nrm_mean[sel],
            targets=tg, sel=sel, pds=pds, base=base, anti=anti, cd4c=cd4c,
            c_raw=c_raw[sel], c_lolo=c_lolo[sel], n_lolo=n_lolo[sel], c_rel=c_rel[sel],
            c_rel_lolo=c_rel_lolo[sel],
            dd=dd[fold][sel], expr=expr[sel], n_pairs=n_pairs[sel], dropped=list(drop))

    print("KA-1/KA-2:", json.dumps(ka, indent=1))

    # --- the prediction table ---------------------------------------------------------------
    PRED = [("C_raw (all 26 sources)", "c_raw"), ("C_lolo (held line removed)", "c_lolo"),
            ("C_rel (reliability-corrected, 4 sources)", "c_rel"),
            ("C_rel_lolo (reliability-corrected, held line out)", "c_rel_lolo"),
            ("baseline: donor-donor agreement", "dd"),
            ("baseline: CD4 condition count", "cd4c"),
            ("baseline: target expression (control log1p)", "expr"),
            ("control: measured reliability of the non-held sources", "rel_mean"),
            ("control: response magnitude of the non-held sources", "nrm_mean"),
            ("control: number of pairs", "n_pairs")]
    from scipy.stats import rankdata

    def pool(kk: str) -> np.ndarray:
        """within-fold rank of a predictor, concatenated over the three cancer folds; NaN kept."""
        out = []
        for fold in CANCER:
            v = np.asarray(per_fold[fold][kk], dtype=np.float64)
            r = np.full(v.size, np.nan)
            ok = np.isfinite(v)
            r[ok] = rankdata(v[ok]) / ok.sum()
            out.append(r)
        return np.concatenate(out)

    py = pool("pds")
    pa = np.concatenate([per_fold[f]["anti"] for f in CANCER])
    pooled_cols = {kk: pool(kk) for _, kk in PRED}
    table = []
    for label, kk in PRED:
        row = dict(predictor=label, key=kk)
        for fold in CANCER + ("kolf2",):
            row[f"rho|{fold}"] = boot_spearman(per_fold[fold][kk], per_fold[fold]["pds"])
        row["rho|pooled3"] = boot_spearman(pooled_cols[kk], py)
        row["auc|pooled3_anti"] = auc_rest_vs_group(pooled_cols[kk], pa)
        table.append(row)
        print(f"{label:52s} pooled rho {row['rho|pooled3']['rho']:+.3f} "
              f"[{row['rho|pooled3']['ci'][0]:+.3f},{row['rho|pooled3']['ci'][1]:+.3f}]  "
              f"AUC(rest>anti) {row['auc|pooled3_anti']['auc']:.3f}")

    # --- is C_lolo more than "the response is measurable"? ----------------------------------
    paired, partial = {}, {}
    for other in ("dd", "rel_mean", "nrm_mean", "n_pairs", "c_rel_lolo"):
        paired[f"c_lolo - {other}"] = paired_delta_spearman(
            pooled_cols["c_lolo"], pooled_cols[other], py)
        partial[f"c_lolo | {other}"] = partial_spearman(
            pooled_cols["c_lolo"], py, pooled_cols[other])
        print(f"  c_lolo vs {other:10s} paired dRho {paired[f'c_lolo - {other}']['delta']:+.3f} "
              f"{[round(v,3) for v in paired[f'c_lolo - {other}']['ci']]}   "
              f"partial rho {partial[f'c_lolo | {other}']['rho']:+.3f} "
              f"{[round(v,3) for v in partial[f'c_lolo | {other}']['ci']]}")

    # --- h1-test: only 5 of its 85 targets are on the 2026 panel, so the map is rebuilt on ITS
    #     targets, with the same frozen construction (same sources, same axis, same centring
    #     rule applied to the rows the source contributes to THIS map). --------------------------
    zi = np.load(L.ROOT / "experiments/273-flex-reference-instrument/work/instrument273.npz",
                 allow_pickle=True)
    h1t = [str(t) for t in zi["targets"]]
    hpds = np.asarray(zi["INCUMBENT|pds_eff"], dtype=np.float64)
    hcos = np.asarray(zi["base_cos_eff"], dtype=np.float64)
    axis = [str(g) for g in zrow["axis"]]
    Rh = np.zeros((len(keys), len(h1t), len(axis)), dtype=np.float32)
    Ah = np.zeros((len(keys), len(h1t)), dtype=bool)
    for i, k in enumerate(keys):
        Rh[i], Ah[i] = L.load_source(k, axis, targets=h1t)
    Ch = np.full((len(pairs), len(h1t)), np.nan)
    for p, (i, j) in enumerate(pairs):
        m = Ah[i] & Ah[j]
        if m.any():
            Ch[p, m] = L.cosine_rows(Rh[i][m].astype(np.float64), Rh[j][m].astype(np.float64))
    with np.errstate(invalid="ignore"):
        h_raw = np.nanmean(Ch, axis=0)
        h_np = np.isfinite(Ch).sum(axis=0)
    h_lolo, h_nlolo = lolo_score(keys, Ch, pairs, L.HELD_LINE_SOURCES["H1"])
    # the two-donor baseline on h1's targets
    def pcos_h(a, b):
        i, j = keys.index(a), keys.index(b)
        m = Ah[i] & Ah[j]
        o = np.full(len(h1t), np.nan)
        o[m] = L.cosine_rows(Rh[i][m].astype(np.float64), Rh[j][m].astype(np.float64))
        return o
    with np.errstate(invalid="ignore"):
        trio_h = np.nanmean(np.vstack([pcos_h("k562_gwps", "orion_hct116"),
                                       pcos_h("k562_gwps", "orion_hek293t"),
                                       pcos_h("orion_hct116", "orion_hek293t")]), axis=0)
    expr_h = np.array([bc[gi[t]] if t in gi else np.nan for t in h1t])
    h1 = {}
    for label, vec in (("C_raw (all sources incl. H1)", h_raw),
                       ("C_lolo (H1 sources removed)", h_lolo),
                       ("baseline: donor-donor agreement", trio_h),
                       ("baseline: target expression", expr_h),
                       ("control: number of pairs", h_np.astype(float))):
        h1[label] = dict(rho_pds=boot_spearman(vec, hpds), rho_cos=boot_spearman(vec, hcos))
        print(f"h1-test  {label:34s} rho(pds) {h1[label]['rho_pds']['rho']:+.3f} "
              f"{[round(v,3) for v in h1[label]['rho_pds']['ci']]}  n={h1[label]['rho_pds']['n']}")
    n_lolo_h1 = h_nlolo

    out = dict(ka=ka, table=table, h1_test=h1, paired=paired, partial=partial,
               h1_pairs_per_target=dict(median_all=float(np.median(h_np)),
                                        median_lolo=float(np.median(h_nlolo))),
               lolo_dropped={f: per_fold[f]["dropped"] for f in FOLDS},
               n_lolo_pairs_median={f: float(np.median(per_fold[f]["n_lolo"])) for f in FOLDS},
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (L.WORK / "predict276.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(L.WORK / "predict276.npz",
                        **{f"{f}|{k}": np.asarray(v) for f, d in per_fold.items()
                           for k, v in d.items() if k not in ("dropped",)},
                        h1_targets=np.asarray(h1t, dtype=object), h1_pds=hpds, h1_cos=hcos,
                        h1_c_lolo=h_lolo, h1_c_raw=h_raw, h1_n_pairs=h_np,
                        h1_dd=trio_h, h1_expr=expr_h)
    print(f"peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()

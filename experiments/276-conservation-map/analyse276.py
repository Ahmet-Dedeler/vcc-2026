"""276 items 3b + 5 + the deliverable table.

  (a) a second, independent reproducibility test of the score: split the CELL LINES into two
      disjoint groups, build the map inside each, and correlate the two per-target scores —
      "does conservation measured on one set of lines predict conservation on a disjoint set?"
  (b) the LINCS/CMap version of the same statistic (declared secondary, DEFINITION-276 §3),
      cross-line cosine on the 978 measured landmark genes, and its prediction of transfer;
  (c) the top / bottom decile by name and what they are made of (expression breadth across
      CCLE lines, essentiality, share of the response energy in the production whitener's
      top-decile — the shared programme);
  (d) work/conservation276.csv, the per-target deliverable.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse276.py
"""
from __future__ import annotations

import csv
import json
import resource
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

import lib276 as L

SEED = 20260906
B = 2000
E262 = L.ROOT / "experiments/262-where-pds-is-lost/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
HELD = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}


def boot_spearman(x, y, b=B, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    if x.size < 5:
        return dict(rho=None, ci=[None, None], n=int(x.size))
    rng = np.random.default_rng(seed)
    d = np.empty(b)
    for i in range(b):
        s = rng.integers(0, x.size, x.size)
        d[i] = spearmanr(x[s], y[s]).statistic
    return dict(rho=float(spearmanr(x, y).statistic),
                ci=[float(v) for v in np.percentile(d, [2.5, 97.5])], n=int(x.size))


def main() -> None:
    zc = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    keys = [str(k) for k in zc["keys"]]
    Ccross = np.asarray(zc["Ccross"], dtype=np.float64)
    pairs = [tuple(p) for p in np.asarray(zc["pairs_cross"])]
    c_raw = np.asarray(zc["c_raw"]); n_pairs = np.asarray(zc["n_pairs"])
    null_mean = np.asarray(zc["null_mean"])
    boot_lo, boot_hi = np.asarray(zc["boot_lo"]), np.asarray(zc["boot_hi"])
    c_within = np.asarray(zc["c_within"])
    zr = np.load(L.WORK / "rel276.npz", allow_pickle=True)
    c_rel = np.asarray(zr["c_rel"]); n_rel = np.asarray(zr["n_rel"])
    zrow = np.load(L.WORK / "rows276.npz", allow_pickle=True)
    R = np.asarray(zrow["R"]); A = np.asarray(zrow["A"])
    axis = [str(g) for g in zrow["axis"]]
    lines = [L.SOURCES[k][1] for k in keys]
    uniq_lines = sorted(set(lines))
    out: dict = {}

    # ---- (a) disjoint-LINE split reproducibility -------------------------------------------
    def score_from_lines(group: set[str]) -> tuple[np.ndarray, np.ndarray]:
        keep = np.asarray([p for p, (i, j) in enumerate(pairs)
                           if lines[i] in group and lines[j] in group])
        if keep.size == 0:
            return np.full(len(L.PANEL), np.nan), np.zeros(len(L.PANEL), int)
        with np.errstate(invalid="ignore"):
            return np.nanmean(Ccross[keep], axis=0), np.isfinite(Ccross[keep]).sum(axis=0)

    declared_G1 = {"K562", "HEK293T", "HEK293", "KOLF2.1J", "cardiomyocyte", "A549",
                   "HAP1", "MCF7", "melanoma", "Calu-3"}
    declared_G2 = set(uniq_lines) - declared_G1
    s1, n1 = score_from_lines(declared_G1)
    s2, n2 = score_from_lines(declared_G2)
    dec = boot_spearman(s1, s2)
    print(f"(a) declared disjoint-line split: rho {dec['rho']:+.3f} "
          f"[{dec['ci'][0]:+.3f},{dec['ci'][1]:+.3f}] on n={dec['n']}  "
          f"(G1 {len(declared_G1)} lines, G2 {len(declared_G2)})")

    rng = np.random.default_rng(SEED)
    rand = []
    for _ in range(50):
        perm = list(uniq_lines)
        rng.shuffle(perm)
        g1 = set(perm[: len(perm) // 2])
        a, _ = score_from_lines(g1)
        b_, _ = score_from_lines(set(uniq_lines) - g1)
        ok = np.isfinite(a) & np.isfinite(b_)
        if ok.sum() >= 30:
            rand.append(float(spearmanr(a[ok], b_[ok]).statistic))
    rand = np.asarray(rand)
    print(f"    50 random line splits: median rho {np.median(rand):+.3f} "
          f"[{np.percentile(rand,10):+.3f}, {np.percentile(rand,90):+.3f}], "
          f"share > 0 {float((rand>0).mean()):.2f}")
    out["disjoint_line_split"] = dict(
        declared=dict(G1=sorted(declared_G1), G2=sorted(declared_G2), **dec),
        random=dict(n=int(rand.size), median=float(np.median(rand)),
                    p10=float(np.percentile(rand, 10)), p90=float(np.percentile(rand, 90)),
                    share_positive=float((rand > 0).mean()), values=[float(v) for v in rand]))

    # ---- (b) the LINCS/CMap version of the same statistic -----------------------------------
    gi_txt = (L.ROOT / "data/external/cmap/geneinfo_beta.txt").read_text().splitlines()
    hdr = gi_txt[0].split("\t")
    ci_id, ci_sym, ci_space = hdr.index("gene_id"), hdr.index("gene_symbol"), hdr.index("feature_space")
    ent2sym, landmark = {}, set()
    for ln in gi_txt[1:]:
        f = ln.split("\t")
        ent2sym[f[ci_id]] = f[ci_sym]
        if f[ci_space] == "landmark":
            landmark.add(f[ci_id])
    cmap = {}
    for tag, fn in (("cmap_sh", "cmap270_trt_sh_cgs.npz"), ("cmap_xpr", "cmap270_trt_xpr.npz")):
        z = np.load(L.ROOT / "experiments/270-data-reconnaissance/work" / fn, allow_pickle=True)
        genes = [str(g) for g in z["genes"]]
        lm = np.asarray([i for i, g in enumerate(genes) if g in landmark], dtype=np.int64)
        syms = [ent2sym.get(genes[i], "") for i in lm]
        tgt = np.asarray([str(t) for t in z["target"]], dtype=object)
        lin = np.asarray([str(x) for x in z["line"]], dtype=object)
        Z = np.asarray(z["Z"], dtype=np.float64)[:, lm]
        L_lines = sorted(set(lin))
        Rc = np.full((len(L_lines), len(L.PANEL), lm.size), np.nan)
        for li, ln_ in enumerate(L_lines):
            m = lin == ln_
            rows = {t: Z[k] for k, t in zip(np.flatnonzero(m), tgt[m])}
            have = [t for t in L.PANEL if t in rows]
            if not have:
                continue
            M = np.vstack([rows[t] for t in have])
            M -= M.mean(axis=0, keepdims=True)      # same centring rule
            for q, t in enumerate(have):
                Rc[li, L.PANEL.index(t)] = M[q]
        cmap[tag] = dict(lines=L_lines, R=Rc, syms=syms)
        prs = [(i, j) for i in range(len(L_lines)) for j in range(i + 1, len(L_lines))]
        C = np.full((len(prs), len(L.PANEL)), np.nan)
        for p, (i, j) in enumerate(prs):
            m = np.isfinite(Rc[i][:, 0]) & np.isfinite(Rc[j][:, 0])
            if m.any():
                C[p, m] = L.cosine_rows(Rc[i][m], Rc[j][m])
        with np.errstate(invalid="ignore"):
            cmap[tag]["score"] = np.nanmean(C, axis=0)
            cmap[tag]["n"] = np.isfinite(C).sum(axis=0)
            cmap[tag]["C"] = C
            cmap[tag]["prs"] = prs
        print(f"(b) {tag}: {len(L_lines)} lines, {lm.size} landmark genes, "
              f"score median {np.nanmedian(cmap[tag]['score']):+.4f} on "
              f"{int(np.isfinite(cmap[tag]['score']).sum())} targets")

    # its prediction of transfer, pooled over the three cancer folds, with LOLO
    idx = {t: i for i, t in enumerate(L.PANEL)}
    pooled_y, pooled = [], {k: [] for k in ("cmap_sh", "cmap_xpr", "cmap_sh_lolo", "cmap_xpr_lolo")}
    for fold in CANCER:
        z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]
        sel = np.asarray([idx[t] for t in tg])
        pds = np.asarray(z["cd4svd5|pds"], dtype=np.float64)
        r = np.full(pds.size, np.nan); r[:] = rankdata(pds) / pds.size
        pooled_y.append(r)
        held = L.CMAP_HELD_LINE[HELD[fold]]
        for tag in ("cmap_sh", "cmap_xpr"):
            for suffix, drop in (("", ()), ("_lolo", held)):
                Ls = cmap[tag]["lines"]
                di = {Ls.index(h) for h in drop if h in Ls}
                keep = [p for p, (i, j) in enumerate(cmap[tag]["prs"])
                        if i not in di and j not in di]
                with np.errstate(invalid="ignore"):
                    s = np.nanmean(cmap[tag]["C"][keep], axis=0)[sel]
                v = np.full(s.size, np.nan)
                ok = np.isfinite(s)
                v[ok] = rankdata(s[ok]) / ok.sum()
                pooled[tag + suffix].append(v)
    py = np.concatenate(pooled_y)
    cmap_pred = {k: boot_spearman(np.concatenate(v), py) for k, v in pooled.items()}
    for k, v in cmap_pred.items():
        print(f"    {k:16s} pooled rho {v['rho']:+.3f} "
              f"[{v['ci'][0]:+.3f},{v['ci'][1]:+.3f}] n={v['n']}")
    out["cmap"] = dict(
        summary={t: dict(n_lines=len(cmap[t]["lines"]), lines=cmap[t]["lines"],
                         n_landmark=len(cmap[t]["syms"]),
                         n_targets=int(np.isfinite(cmap[t]["score"]).sum()),
                         median=float(np.nanmedian(cmap[t]["score"])))
                 for t in ("cmap_sh", "cmap_xpr")},
        prediction=cmap_pred,
        rho_with_C_raw={t: boot_spearman(cmap[t]["score"], c_raw) for t in ("cmap_sh", "cmap_xpr")})
    print(f"    rho(C_cmap_sh, C_raw)  {out['cmap']['rho_with_C_raw']['cmap_sh']['rho']:+.3f}; "
          f"rho(C_cmap_xpr, C_raw) {out['cmap']['rho_with_C_raw']['cmap_xpr']['rho']:+.3f}")

    # ---- (c) what the deciles are made of ---------------------------------------------------
    # expression breadth across CCLE lines
    zz = np.load(L.ROOT / "data/external/ccle_broad_2025/ccle_cache.npz", allow_pickle=True)
    cg = {str(g): i for i, g in enumerate(zz["genes"])}
    X = zz["X"]
    breadth = np.full(len(L.PANEL), np.nan); ccle_med = np.full(len(L.PANEL), np.nan)
    ccle_sd = np.full(len(L.PANEL), np.nan)
    for i, t in enumerate(L.PANEL):
        if t in cg:
            v = np.asarray(X[:, cg[t]], dtype=np.float64)
            breadth[i] = float((v > 1.0).mean())
            ccle_med[i] = float(np.median(v)); ccle_sd[i] = float(np.std(v))
    ess = set((L.ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
    essential = np.asarray([1.0 if t in ess else 0.0 for t in L.PANEL])

    # the production whitener's top decile = the shared programme
    shared_share = np.full(len(L.PANEL), np.nan)
    wmeta = json.loads((L.ROOT / "data/donors/whitener_union_pruned.json").read_text())
    union_genes: set[str] = set()
    for k in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        _, g, _ = L.read_bank(L.SOURCES[k][0])
        union_genes |= set(g) & set(L.OFFICIAL)
    wcols = [g for g in L.OFFICIAL if g in union_genes and g not in set(L.PANEL)]
    print(f"(c) reconstructed whitener axis: {len(wcols)} columns (meta n_cols {wmeta['n_cols']})")
    if len(wcols) == int(wmeta["n_cols"]):
        w = np.load(L.ROOT / "data/donors/whitener_union_pruned.npz")
        ev, V = np.asarray(w["evals"]), np.asarray(w["evecs"], dtype=np.float32)
        order = np.argsort(ev)[::-1]
        top = order[: max(1, len(order) // 10)]
        Vt = np.ascontiguousarray(V[:, top])
        del V
        wi = {g: i for i, g in enumerate(wcols)}
        cols_in_axis = [(wi[g], q) for q, g in enumerate(axis) if g in wi]
        wsel = np.asarray([a for a, _ in cols_in_axis]); asel = np.asarray([b for _, b in cols_in_axis])
        Vs = Vt[wsel]                                    # (|axis ∩ wcols|, top)
        with np.errstate(invalid="ignore"):
            M = np.nanmean(np.where(A[:, :, None], R, np.nan), axis=0)   # mean response across sources
        Ma = np.asarray(M[:, asel], dtype=np.float32)
        proj = Ma @ Vs
        num = (proj ** 2).sum(axis=1)
        den = (Ma.astype(np.float64) ** 2).sum(axis=1)
        shared_share = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
        del Vt, Vs
        print(f"    whitener top decile: {len(top)} directions, "
              f"{len(wsel)} of the map's {len(axis)} genes on its axis; "
              f"median energy share {np.nanmedian(shared_share):.3f}")
    else:
        print("    whitener axis could not be reconstructed exactly; shared-programme share skipped")

    dec_cols = dict(breadth=breadth, ccle_median_tpm=ccle_med, ccle_sd=ccle_sd,
                    essential=essential, shared_share=shared_share,
                    n_pairs=n_pairs.astype(float), c_within=c_within)
    q = np.nanpercentile(c_raw, [10, 90])
    lo_m, hi_m = c_raw <= q[0], c_raw >= q[1]
    deciles = {}
    for name, v in dec_cols.items():
        deciles[name] = dict(bottom=float(np.nanmedian(v[lo_m])), top=float(np.nanmedian(v[hi_m])),
                             all=float(np.nanmedian(v)),
                             rho_with_C=boot_spearman(v, c_raw)["rho"])
        print(f"    {name:16s} bottom {deciles[name]['bottom']:+.3f}  top {deciles[name]['top']:+.3f}"
              f"  rho with C {deciles[name]['rho_with_C']:+.3f}")
    out["deciles"] = dict(
        cut=[float(q[0]), float(q[1])], stats=deciles,
        top_names=[L.PANEL[i] for i in np.argsort(-np.nan_to_num(c_raw, nan=-9))[:30]],
        bottom_names=[L.PANEL[i] for i in np.argsort(np.nan_to_num(c_raw, nan=9))[:30]])
    print("    TOP 15:", ", ".join(out["deciles"]["top_names"][:15]))
    print("    BOTTOM 15:", ", ".join(out["deciles"]["bottom_names"][:15]))

    # ---- (d) the deliverable table ----------------------------------------------------------
    zp = np.load(L.WORK / "predict276.npz", allow_pickle=True)
    per = {}
    for fold in CANCER + ("kolf2",):
        tg = [str(t) for t in zp[f"{fold}|targets"]]
        per[fold] = dict(zip(tg, zp[f"{fold}|pds"]))
    lolo = {}
    for fold in CANCER + ("kolf2",):
        tg = [str(t) for t in zp[f"{fold}|targets"]]
        lolo[fold] = dict(zip(tg, zp[f"{fold}|c_lolo"]))
    rows = []
    for i, t in enumerate(L.PANEL):
        r = dict(target=t, C_raw=c_raw[i], n_pairs=int(n_pairs[i]),
                 boot_lo=boot_lo[i], boot_hi=boot_hi[i], null_mean=null_mean[i],
                 C_within=c_within[i], n_within=int(np.asarray(zc["n_within"])[i]),
                 C_rel=c_rel[i], n_rel_pairs=int(n_rel[i]),
                 C_cmap_sh=cmap["cmap_sh"]["score"][i], C_cmap_xpr=cmap["cmap_xpr"]["score"][i],
                 half_A=np.asarray(zr["sA"])[i], half_B=np.asarray(zr["sB"])[i],
                 ccle_breadth=breadth[i], ccle_median_tpm=ccle_med[i],
                 essential=int(essential[i]), shared_programme_share=shared_share[i])
        for fold in CANCER + ("kolf2",):
            r[f"pds_{fold}"] = per[fold].get(t, np.nan)
            r[f"C_lolo_{fold}"] = lolo[fold].get(t, np.nan)
        rows.append(r)
    with (L.WORK / "conservation276.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, (float, np.floating)) and np.isfinite(v)
                            else ("" if isinstance(v, (float, np.floating)) else v))
                        for k, v in r.items()})
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    (L.WORK / "analyse276.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"wrote work/conservation276.csv ({len(rows)} rows); peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()

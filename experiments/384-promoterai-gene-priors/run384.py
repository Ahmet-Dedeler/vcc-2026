#!/usr/bin/env python3
"""384 — PromoterAI-gene-prior proxies: build per-gene features, test on H1.

Modes: ka (known answers, in-memory, writes nothing outside work/),
build (parse downloads -> work/features384.npz), stats, reach.
Reads 323/189/348 inputs read-only; never writes outside this experiment.
"""
from __future__ import annotations
import argparse, gzip, json, sys, time
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "data"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/348-reach-ordering"))
import reach348 as R348  # noqa: E402 (read-only reuse of walk conventions)

BASE323 = ROOT / "experiments/323-current-family-validation/work"
BASE189 = ROOT / "experiments/189-per-gene-reliability/work"
ALPHA = 0.05
EXP_O1 = 0.25434726355678233
EXP_O4 = 0.28074543623823345
EXP_TAU = 0.3417639845890673


def load_h1():
    d = R348.load()
    n, g = len(d["targets"]), len(d["rgenes"])
    qrm, lrm = R348.build_long(d)
    pp, qp, lp = d["pp"], d["qp"], d["lp"]
    ref_sig = qrm < ALPHA
    own = np.equal.outer(d["targets"], d["rgenes"])
    pool = ref_sig & ~own
    n_conf = pool.sum(1)
    adjud = np.isfinite(lrm) & (lrm != 0)
    commit = np.isfinite(lp) & (lp != 0)
    match = adjud & commit & (np.sign(lp) == np.sign(lrm))
    z = np.load(BASE323 / "reference.npz", allow_pickle=True)
    full_genes, keep = z["genes"].astype(str), z["keep"]
    assert np.array_equal(full_genes[keep], d["rgenes"])
    ref_lin = np.asarray(z["ref_lin"], dtype=float)
    z.close()
    return d, dict(n=n, g=g, qrm=qrm, lrm=lrm, pool=pool, n_conf=n_conf,
                   adjud=adjud, match=match, ref_lin=ref_lin)


def lex_order(feat_names, *key_mats):
    """348 lex() verbatim: keys primary-last, gene name least significant."""
    n = key_mats[0].shape[0]
    idx = np.empty_like(key_mats[0], dtype=int)
    for i in range(n):
        idx[i] = np.lexsort((feat_names[i],) + tuple(k[i] for k in key_mats))
    return idx


def mode_ka():
    d, h = load_h1()
    n, g = h["n"], h["g"]
    feat = np.tile(d["rgenes"], (n, 1))
    qm = np.where(np.isfinite(d["qp"]), d["qp"], np.inf)
    pm = np.where(np.isfinite(d["pp"]), d["pp"], np.inf)
    ala_m = np.where(np.isfinite(np.abs(d["lp"])), -np.abs(d["lp"]), np.inf)
    sig = np.isfinite(d["qp"]) & (d["qp"] < ALPHA)
    cur = lex_order(feat, ala_m, pm, qm, ~sig)
    R1, _ = R348.walk_reach(cur, h["match"] & h["pool"], h["adjud"] & h["pool"], h["n_conf"])
    D = R348.design_matrix(d)
    ad_m = np.where(np.isfinite(np.abs(D)), -np.abs(D), np.inf)
    o4 = lex_order(feat, ad_m)
    R4, _ = R348.walk_reach(o4, h["match"] & h["pool"], h["adjud"] & h["pool"], h["n_conf"])
    m1, m4 = float(np.nanmean(R1)), float(np.nanmean(R4))
    z = np.load(BASE189 / "prior189.npz", allow_pickle=True)
    tau = float(z["tau"])
    th = (z["genes"].astype(str).shape, z["theta"].shape)
    z.close()
    ok = (abs(m1 - EXP_O1) < 1e-9 and abs(m4 - EXP_O4) < 1e-9
          and abs(tau - EXP_TAU) < 1e-12 and th[0] == (10725,) and th[1] == (10725,))
    print(f"O1={m1:.8f} (exp {EXP_O1:.8f}) diff={m1-EXP_O1:.2e}")
    print(f"O4={m4:.8f} (exp {EXP_O4:.8f}) diff={m4-EXP_O4:.2e}")
    print(f"tau={tau:.10f} theta-genes={th}")
    print("KA:", "PASS" if ok else "FAIL")
    return ok


def parse_gtf_tss(gtf_path):
    """gene_name -> (chrom, tss). First gene record wins on duplicates."""
    tss = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            attrs = {}
            for kv in c[8].strip().strip(";").split(";"):
                kv = kv.strip()
                if not kv or " " not in kv:
                    continue
                k, v = kv.split(" ", 1)
                attrs[k] = v.strip().strip('"')
            name = attrs.get("gene_name")
            if not name or name in tss:
                continue
            s, e = int(c[3]), int(c[4])
            tss[name] = (c[0], s if c[6] == "+" else e)
    return tss


def mode_build():
    import pandas as pd
    d, h = load_h1()
    uni = list(d["rgenes"])
    U = {g: i for i, g in enumerate(uni)}
    G = len(uni)
    feats = {}
    cov = {}
    # --- theta_g from 189 ---
    z = np.load(BASE189 / "prior189.npz", allow_pickle=True)
    pg, th = z["genes"].astype(str), np.asarray(z["theta"], float)
    z.close()
    v = np.full(G, np.nan)
    hit = 0
    for gname, t in zip(pg, th):
        if gname in U:
            v[U[gname]] = t
            hit += 1
    feats["theta"] = v
    cov["theta"] = hit
    # --- control expression ---
    feats["log_expr"] = np.log1p(np.maximum(h["ref_lin"], 0))
    cov["log_expr"] = G
    # --- gnomAD constraint (canonical transcript per gene) ---
    # no canonical flag in v2.1.1 by_gene: one row per transcript;
    # per gene symbol keep the max-exp_lof (most powered) transcript.
    loeuf = {}
    gfh = gzip.open(DATA / "gnomad.lof_by_gene.txt.bgz", "rt")
    gdf = pd.read_csv(gfh, sep="\t",
                      usecols=["gene", "gene_id", "exp_lof", "oe_lof_upper",
                               "pLI", "lof_z", "mis_z"], low_memory=False)
    gfh.close()
    gdf["exp_lof"] = pd.to_numeric(gdf["exp_lof"], errors="coerce")
    gdf = gdf.sort_values("exp_lof").drop_duplicates("gene", keep="last")
    for row in gdf.itertuples():
        gid = row.gene_id.split(".")[0] if isinstance(row.gene_id, str) else None
        loeuf[(row.gene, gid)] = (row.oe_lof_upper, row.pLI, row.lof_z, row.mis_z)
    # map: prefer symbol, else ENSG->symbol via GENCODE
    tss = parse_gtf_tss(DATA / "gencode.v44.basic.annotation.gtf.gz")
    print(f"TSS parsed: {len(tss)} genes")
    ensg2sym = {}
    with gzip.open(DATA / "gencode.v44.basic.annotation.gtf.gz", "rt") as f:
        for line in f:
            if not line or line[0] == "#" or "\tgene\t" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            gid, gname = None, None
            for kv in c[8].strip().strip(";").split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id "):
                    gid = kv.split('"')[1].split(".")[0]
                elif kv.startswith("gene_name "):
                    gname = kv.split('"')[1]
            if gid and gname and gid not in ensg2sym:
                ensg2sym[gid] = gname
    vl = np.full(G, np.nan)
    vp = np.full(G, np.nan)
    vz = np.full(G, np.nan)
    vmz = np.full(G, np.nan)
    hl = 0
    for (key, gid), (lu, pliv, lz, mz) in loeuf.items():
        sym = key if (key in U) else ensg2sym.get(gid)
        if sym in U:
            vl[U[sym]] = lu
            vp[U[sym]] = pliv
            vz[U[sym]] = lz
            vmz[U[sym]] = mz
            hl += 1
    feats["loeuf"] = vl  # lower = more constrained
    feats["pli"] = vp
    feats["lof_z"] = vz
    feats["mis_z"] = vmz
    cov["loeuf"] = hl
    cov["pli"] = int(np.isfinite(vp).sum())
    cov["lof_z"] = int(np.isfinite(vz).sum())
    cov["mis_z"] = int(np.isfinite(vmz).sum())
    # --- GTEx breadth ---
    with gzip.open(DATA / "gtex_v8_median_tpm.gct.gz", "rt") as f:
        f.readline()
        dims = f.readline().rstrip("\n").split("\t")
        hdr = f.readline().rstrip("\n").split("\t")
        tissues = hdr[2:]
        rows = {}
        for line in f:
            c = line.rstrip("\n").split("\t")
            s = c[1]
            x = np.array([float(v) for v in c[2:]], dtype=np.float32)
            rows.setdefault(s, []).append(x)
    breadth = {}
    for s, arrs in rows.items():
        m = np.median(np.stack(arrs), axis=0)
        breadth[s] = (int((m >= 1.0).sum()), float(np.log1p(m).mean()))
    vb = np.full(G, np.nan)
    vm = np.full(G, np.nan)
    hb = 0
    for s, (b, m) in breadth.items():
        if s in U:
            vb[U[s]] = b
            vm[U[s]] = m
            hb += 1
    feats["gtex_breadth"] = vb
    feats["gtex_meanlog"] = vm
    cov["gtex_breadth"] = hb
    cov["gtex_meanlog"] = hb
    # --- CpG-island promoter flag ---
    cpg = pd.read_csv(DATA / "cpgIslandExt.hg38.txt.gz", sep="\t", header=None,
                      usecols=[1, 2, 3])
    cpg.columns = ["chrom", "start", "end"]
    by_chr = {c: df for c, df in cpg.groupby("chrom", sort=False)}
    vc = np.full(G, np.nan)
    vn = np.full(G, np.nan)
    hc = 0
    for sym, (chrom, t) in tss.items():
        if sym not in U or chrom not in by_chr:
            continue
        df = by_chr[chrom]
        lo, hi = t - 2000, t + 500
        n = int(((df["start"] < hi) & (df["end"] > lo)).sum())
        vc[U[sym]] = 1.0 if n > 0 else 0.0
        vn[U[sym]] = n
        hc += 1
    feats["cpg_promoter"] = vc
    feats["cpg_count"] = vn
    cov["cpg_promoter"] = hc
    cov["cpg_count"] = hc
    WORK.mkdir(parents=True, exist_ok=True)
    keys = list(feats.keys())
    M = np.stack([feats[k] for k in keys])
    np.savez_compressed(WORK / "features384.npz", genes=np.array(uni),
                        keys=np.array(keys), M=M)
    with open(WORK / "features384_coverage.json", "w") as f:
        json.dump({"universe": G, "coverage": cov,
                   "cpg_islands": int(len(cpg)),
                   "gtex_tissues": len(tissues)}, f, indent=2)
    print("coverage:", json.dumps(cov, indent=2))


def half_acc(h, scored_idx, sel):
    idx = scored_idx[sel]
    P = h["pool"][idx]
    M = (h["match"] & h["pool"])[idx]
    sup = P.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        a = M.sum(0) / np.maximum(sup, 1)
    a[sup == 0] = np.nan
    return a, sup


def spear(a, b):
    from scipy.stats import spearmanr
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10:
        return np.nan, 0
    return spearmanr(a[m], b[m]).statistic, int(m.sum())


def mode_stats():
    from scipy.stats import spearmanr  # noqa: F401
    d, h = load_h1()
    z = np.load(WORK / "features384.npz", allow_pickle=True)
    keys = [str(k) for k in z["keys"]]
    M = z["M"]
    z.close()
    F = {k: M[i] for i, k in enumerate(keys)}
    scored = np.flatnonzero(h["n_conf"] > 0)
    n_sc = len(scored)
    ev, od = np.arange(0, n_sc, 2), np.arange(1, n_sc, 2)
    aA, sA = half_acc(h, scored, ev)
    aB, sB = half_acc(h, scored, od)
    aF, sF = half_acc(h, scored, np.arange(n_sc))
    print(f"scored targets={n_sc} genes={h['g']} "
          f"support>=10: full={int((sF>=10).sum())} A={int((sA>=10).sum())} B={int((sB>=10).sum())}")
    rng = np.random.default_rng(38400)
    NPERM = 500
    out = {"support_ge10": {k: int(v) for k, v in
           [("full", (sF >= 10).sum()), ("A", (sA >= 10).sum()), ("B", (sB >= 10).sum())]},
           "features": {}}
    for k in keys:
        f = F[k]
        r_full, n1 = spear(f[sF >= 10], aF[sF >= 10])
        r_th, n2 = spear(f, F["theta"])
        r_ex, n3 = spear(f, F["log_expr"])
        # held-out spread on half B (deciles of f, accuracy on B;
        # binary features: 1-group vs 0-group)
        def groups(ff):
            u = np.unique(ff[np.isfinite(ff)])
            if len(u) <= 3:
                return ff == u.max(), ff == u.min(), True
            lo, hi = np.nanquantile(ff, [0.1, 0.9])
            return ff >= hi, ff <= lo, False
        vb = np.isfinite(f) & (sB >= 10)
        fb = f[vb]
        ab = aB[vb]
        topB, botB, isbin = groups(fb)
        spreadB = float(np.nanmean(ab[topB]) - np.nanmean(ab[botB]))
        va = np.isfinite(f) & (sA >= 10)
        fa = f[va]
        topA, botA, _ = groups(fa)
        spreadA = float(np.nanmean(aA[va][topA]) - np.nanmean(aA[va][botA]))
        null = np.empty(NPERM)
        perm = fb.copy()
        for s in range(NPERM):
            rng.shuffle(perm)
            t2, b2, _ = groups(perm)
            null[s] = np.nanmean(ab[t2]) - np.nanmean(ab[b2])
        nm, nsd = float(np.mean(null)), float(np.std(null))
        zsc = (spreadB - nm) / nsd if nsd > 0 else 0.0
        p1 = float((null >= spreadB).mean())
        # partial correlation with aF given theta + log_expr
        base = np.stack([F["theta"], F["log_expr"]]).T
        mall = np.isfinite(f) & np.isfinite(aF) & (sF >= 10) & np.all(np.isfinite(base), 1)
        Xb = np.column_stack([np.ones(mall.sum()), base[mall]])
        y = aF[mall]
        beta_b, *_ = np.linalg.lstsq(Xb, y, rcond=None)
        r2b = 1 - np.var(y - Xb @ beta_b) / np.var(y)
        Xf = np.column_stack([Xb, f[mall]])
        beta_f, *_ = np.linalg.lstsq(Xf, y, rcond=None)
        r2f = 1 - np.var(y - Xf @ beta_f) / np.var(y)
        res_y = y - Xb @ beta_b
        res_f = f[mall] - Xb @ np.linalg.lstsq(Xb, f[mall], rcond=None)[0]
        pr = float(np.corrcoef(res_y, res_f)[0, 1])
        out["features"][k] = {
            "binary": bool(isbin),
            "n_full": n1, "spear_acc": None if np.isnan(r_full) else float(r_full),
            "n_th": n2, "spear_theta": None if np.isnan(r_th) else float(r_th),
            "n_ex": n3, "spear_expr": None if np.isnan(r_ex) else float(r_ex),
            "spreadA": spreadA, "spreadB": spreadB,
            "null_mean": nm, "null_sd": nsd, "z": float(zsc), "p_one_sided": p1,
            "r2_base": float(r2b), "r2_plus": float(r2f),
            "dR2": float(r2f - r2b), "partial_r": pr, "n_ols": int(mall.sum())}
        print(f"{k:12s} sp_acc={r_full:+.4f} sp_th={r_th:+.4f} sp_ex={r_ex:+.4f} "
              f"sprA={spreadA:+.4f} sprB={spreadB:+.4f} null={nm:+.4f}±{nsd:.4f} "
              f"z={zsc:+.2f} p={p1:.3f} dR2={r2f-r2b:+.5f} pr={pr:+.4f}")
    with open(WORK / "stats384.json", "w") as fo:
        json.dump(out, fo, indent=2)


def mode_reach():
    d, h = load_h1()
    n = h["n"]
    feat = np.tile(d["rgenes"], (n, 1))
    qm = np.where(np.isfinite(d["qp"]), d["qp"], np.inf)
    pm = np.where(np.isfinite(d["pp"]), d["pp"], np.inf)
    ala_m = np.where(np.isfinite(np.abs(d["lp"])), -np.abs(d["lp"]), np.inf)
    sig = np.isfinite(d["qp"]) & (d["qp"] < ALPHA)
    D = R348.design_matrix(d)
    ad = np.abs(D)
    ad_m = np.where(np.isfinite(ad), -ad, np.inf)
    z = np.load(WORK / "features384.npz", allow_pickle=True)
    keys = [str(k) for k in z["keys"]]
    M = z["M"]
    z.close()
    F = {k: M[i] for i, k in enumerate(keys)}
    st = json.load(open(WORK / "stats384.json"))["features"]
    # best non-trivial feature by held-out spreadB (orient: higher accuracy first)
    cand = [k for k in keys if k not in ("theta", "log_expr")]
    best = max(cand, key=lambda k: st[k]["spreadB"])
    sgn = 1.0 if (st[best]["spear_acc"] or 0) >= 0 else -1.0
    fb = np.tile(np.where(np.isfinite(F[best]), sgn * F[best], np.inf), (n, 1))
    schemes = {
        "O1_current": lex_order(feat, ala_m, pm, qm, ~sig),
        "O4_abs_design": lex_order(feat, ad_m),
        "theta_desc": lex_order(feat, np.tile(np.where(
            np.isfinite(F["theta"]), -F["theta"], np.inf), (n, 1))),
        "expr_desc": lex_order(feat, np.tile(np.where(
            np.isfinite(F["log_expr"]), -F["log_expr"], np.inf), (n, 1))),
        f"FEAT_{best}": lex_order(feat, fb),
        f"FEAT_{best}_then_design": lex_order(feat, ad_m, fb),
    }
    res = {}
    for name, idx in schemes.items():
        Rr, _ = R348.walk_reach(idx, h["match"] & h["pool"],
                                 h["adjud"] & h["pool"], h["n_conf"])
        assert set(np.flatnonzero(np.isnan(Rr))) == set(np.flatnonzero(h["n_conf"] == 0))
        res[name] = Rr
        print(f"{name}: mean={np.nanmean(Rr):.6f} dO1={np.nanmean(Rr)-np.nanmean(res['O1_current']):+.6f}")
    # shuffled-feature control (3 seeds)
    rng = np.random.default_rng(38401)
    for s in range(3):
        sh = np.tile(rng.permutation(np.where(np.isfinite(F[best]), F[best], 0)), (n, 1))
        Rr, _ = R348.walk_reach(lex_order(feat, sh), h["match"] & h["pool"],
                                 h["adjud"] & h["pool"], h["n_conf"])
        res[f"FEAT_shuf{s}"] = Rr
        print(f"FEAT_shuf{s}: mean={np.nanmean(Rr):.6f}")
    # paired-target bootstrap (10k) for best-feature deltas
    scored = np.flatnonzero(h["n_conf"] > 0)
    B = 10000
    rb = np.random.default_rng(38402).integers(0, len(scored), size=(B, len(scored)))
    out = {"best": best, "sign": sgn, "schemes": {},
           "n_scored": len(scored)}
    base1, base4 = res["O1_current"][scored], res["O4_abs_design"][scored]
    for name in schemes:
        v = res[name][scored]
        d1, d4 = v - base1, v - base4
        out["schemes"][name] = {"mean": float(np.nanmean(v)),
                                "dO1": float(np.nanmean(d1)),
                                "dO4": float(np.nanmean(d4))}
    for tgt, ref in (("dO1", base1), ("dO4", base4)):
        v = res[f"FEAT_{best}"][scored] - ref
        boots = v[rb].mean(1)
        lo, hi = np.quantile(boots, [0.025, 0.975])
        out["schemes"][f"FEAT_{best}"][f"boot_{tgt}"] = {
            "mean": float(boots.mean()), "lo": float(lo), "hi": float(hi),
            "p_gt0": float((boots > 0).mean())}
        v2 = res[f"FEAT_{best}_then_design"][scored] - ref
        b2 = v2[rb].mean(1)
        lo2, hi2 = np.quantile(b2, [0.025, 0.975])
        out["schemes"][f"FEAT_{best}_then_design"][f"boot_{tgt}"] = {
            "mean": float(b2.mean()), "lo": float(lo2), "hi": float(hi2),
            "p_gt0": float((b2 > 0).mean())}
        print(f"FEAT {tgt}: mean={boots.mean():+.5f} CI=[{lo:+.5f},{hi:+.5f}] P(>0)={(boots>0).mean():.3f}")
        print(f"FEAT+design {tgt}: mean={b2.mean():+.5f} CI=[{lo2:+.5f},{hi2:+.5f}] P(>0)={(b2>0).mean():.3f}")
    with open(WORK / "reach384.json", "w") as fo:
        json.dump(out, fo, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["ka", "build", "stats", "reach"])
    a = ap.parse_args()
    t0 = time.monotonic()
    ok = True
    if a.mode == "ka":
        ok = mode_ka()
    elif a.mode == "build":
        mode_build()
    elif a.mode == "stats":
        mode_stats()
    elif a.mode == "reach":
        mode_reach()
    print(f"elapsed {time.monotonic()-t0:.1f}s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

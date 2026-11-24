"""302 — claim C6, computed from the sums `access302.py` wrote.  `X` is NOT opened here.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse302.py [--tag mock]

Writes work/stats302[_tag].json, work/symmetry302[_tag].csv and prints the verdict line.
"""
from __future__ import annotations

import argparse
import csv
import json
import resource
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib302 as L  # noqa: E402

LOG: list[str] = []


def P(m: str = "") -> None:
    print(m, flush=True)
    LOG.append(m)


def cos_pair(i: int, j: int, G, Aown, Bown, na2, nb2, own, owncol_pos) -> float:
    """cos(a_i, b_j) with the own columns of BOTH TFs zeroed in both vectors."""
    cols = {own[i], own[j]} - {-1}
    num = G[i, j]
    ai2, bj2 = na2[i], nb2[j]
    for k in cols:
        c = owncol_pos[k]
        num -= Aown[i, c] * Bown[j, c]
        ai2 -= Aown[i, c] ** 2
        bj2 -= Bown[j, c] ** 2
    d = np.sqrt(max(ai2, 0.0) * max(bj2, 0.0))
    return float(num / d) if d > 0 else np.nan


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    tag = f"_{args.tag}" if args.tag else ""
    out: dict = {}

    z = np.load(L.WORK / f"sums302{tag}.npz", allow_pickle=True)
    axis = np.asarray(z["axis"]).astype(str).tolist()
    scored = np.asarray(z["scored"]).astype(str).tolist()
    primary = np.asarray(z["primary"]).astype(str).tolist()
    secondary = np.asarray(z["secondary"]).astype(str).tolist()
    used_codes = np.asarray(z["used_codes"]).astype(int).tolist()
    acc_sym, acc_code = z["acc_sym"], z["acc_code"]
    lib_sym, lib_code = z["lib_sym"], z["lib_code"]
    cells_sym, cells_code = z["cells_sym"], z["cells_code"]
    draws = np.asarray(z["draws"])
    g_gfp, g_mch = int(z["g_gfp"]), int(z["g_mch"])
    c_gfp, c_mch = int(z["c_gfp"]), int(z["c_mch"])
    axis_h1 = np.asarray(z["axis_h1"])
    nA, nT = len(axis), len(scored)
    P(f"sums: axis {nA}, scored {nT} (primary {len(primary)}, secondary {len(secondary)}), "
      f"{len(used_codes)} ORF codes")

    # ------------------------------------------------------------------ ORF side
    S_sym = acc_sym.sum(axis=1)                       # (nT+2, nA)
    Lsym = lib_sym.sum(axis=1)                        # (nT+2,)
    Csym = cells_sym.sum(axis=1)
    Pr = L.profile(S_sym, Lsym[:, None])              # log1p CPM profiles
    p_ctrl = L.profile(S_sym[g_gfp] + S_sym[g_mch], Lsym[g_gfp] + Lsym[g_mch])
    mu_T = Pr[:nT].mean(axis=0)
    A = Pr[:nT] - mu_T                                # ORF residual (control cancels exactly)
    R = Pr[:nT] - p_ctrl                              # ORF response (for magnitude only)

    # ------------------------------------------------------------------ CRISPRi side
    perts, h1_genes, _ = L.read_h1_bank()
    h1_row = {p: i for i, p in enumerate(perts)}
    Dm = np.load(L.H1BANK / "delta.npy", mmap_mode="r")
    rows = np.asarray([h1_row[s] for s in scored])
    lo, hi = int(rows.min()), int(rows.max())
    Dsub = np.asarray(Dm[lo:hi + 1])[rows - lo][:, axis_h1].astype(np.float64)
    del Dm
    B = Dsub - Dsub.mean(axis=0, keepdims=True)
    P(f"CRISPRi rows {B.shape}; ORF residual {A.shape}")

    # ------------------------------------------------------------------ own columns
    apos = {g: i for i, g in enumerate(axis)}
    own = np.asarray([apos.get(s, -1) for s in scored])
    own_cols = sorted({int(k) for k in own if k >= 0})
    owncol_pos = {k: i for i, k in enumerate(own_cols)}
    Aown = A[:, own_cols].copy()
    Bown = B[:, own_cols].copy()
    na2 = (A ** 2).sum(axis=1)
    nb2 = (B ** 2).sum(axis=1)
    G = A @ B.T
    P(f"own gene on the axis for {len(own_cols)} of {nT} scored symbols")

    Sg = np.asarray([cos_pair(i, i, G, Aown, Bown, na2, nb2, own, owncol_pos)
                     for i in range(nT)])

    # ------------------------------------------------------------------ KA-G, KA-F
    hv = own >= 0
    orf_own = np.asarray([A[i, own[i]] if own[i] >= 0 else np.nan for i in range(nT)])
    kd_own = np.asarray([B[i, own[i]] if own[i] >= 0 else np.nan for i in range(nT)])
    pri = np.asarray([s in set(primary) for s in scored])
    m = hv & pri
    ka_g = dict(n=int(m.sum()),
                orf_median=float(np.median(orf_own[m])), orf_frac_pos=float((orf_own[m] > 0).mean()),
                kd_median=float(np.median(kd_own[m])), kd_frac_neg=float((kd_own[m] < 0).mean()))
    ka_g["pass_"] = bool(ka_g["orf_median"] > 0 and ka_g["orf_frac_pos"] >= 0.8
                         and ka_g["kd_median"] < 0 and ka_g["kd_frac_neg"] >= 0.8)
    P(f"KA-G  own-gene coordinate on {ka_g['n']} primary TFs: ORF residual median "
      f"{ka_g['orf_median']:+.4f} ({ka_g['orf_frac_pos']:.1%} positive); CRISPRi residual median "
      f"{ka_g['kd_median']:+.4f} ({ka_g['kd_frac_neg']:.1%} negative) "
      f"-> {'PASS' if ka_g['pass_'] else 'FAIL'}")

    v = L.profile(S_sym[g_gfp], Lsym[g_gfp]) - L.profile(S_sym[g_mch], Lsym[g_mch])
    gv = B @ v
    vn2 = float(v @ v)
    ctrl_cos = np.empty(nT)
    for i in range(nT):
        num, a2, b2 = gv[i], vn2, nb2[i]
        if own[i] >= 0:
            k = own[i]
            num -= v[k] * B[i, k]; a2 -= v[k] ** 2; b2 -= B[i, k] ** 2
        ctrl_cos[i] = num / np.sqrt(max(a2, 0) * max(b2, 0))
    ka_f = dict(median_abs=float(np.median(np.abs(ctrl_cos[pri]))),
                max_abs=float(np.abs(ctrl_cos[pri]).max()),
                median=float(np.median(ctrl_cos[pri])),
                share_gt05=float((np.abs(ctrl_cos[pri]) > 0.05).mean()))
    ka_f["pass_"] = bool(ka_f["median_abs"] < 0.05)
    P(f"KA-F  GFP-vs-mCherry no-signal row vs the knockdown rows: median |cos| "
      f"{ka_f['median_abs']:.4f} (median {ka_f['median']:+.4f}, max |cos| {ka_f['max_abs']:.4f}, "
      f"{ka_f['share_gt05']:.1%} above 0.05) -> {'PASS' if ka_f['pass_'] else 'FAIL'}")

    # ------------------------------------------------------------------ ORF reliability
    rel = np.zeros(nT)
    for (a_sh, b_sh) in [(d[:L.N_SHARDS // 2], d[L.N_SHARDS // 2:]) for d in draws]:
        SA, LA = acc_sym[:, a_sh].sum(axis=1), lib_sym[:, a_sh].sum(axis=1)
        SB, LB = acc_sym[:, b_sh].sum(axis=1), lib_sym[:, b_sh].sum(axis=1)
        PA, PB = L.profile(SA, LA[:, None]), L.profile(SB, LB[:, None])
        AA = PA[:nT] - PA[:nT].mean(axis=0)
        BB = PB[:nT] - PB[:nT].mean(axis=0)
        for i in range(nT):
            x, y = AA[i], BB[i]
            num, a2, b2 = float(x @ y), float(x @ x), float(y @ y)
            if own[i] >= 0:
                k = own[i]
                num -= x[k] * y[k]; a2 -= x[k] ** 2; b2 -= y[k] ** 2
            d = np.sqrt(max(a2, 0) * max(b2, 0))
            rel[i] += (num / d if d > 0 else 0.0) / len(draws)
    P(f"ORF split-half reliability ({len(draws)} draws, 16 shards): primary median "
      f"{np.median(rel[pri]):+.4f}, q10 {np.percentile(rel[pri], 10):+.4f}, "
      f"q90 {np.percentile(rel[pri], 90):+.4f}")

    # ------------------------------------------------------------------ covariates
    cells = Csym[:nT].astype(float)
    umi = Lsym[:nT] / np.maximum(cells, 1)
    n_iso = np.asarray([sum(1 for c in used_codes if True) for _ in range(nT)])  # replaced below
    code_sym = {}
    orfmap = L.read_orf_map()
    for c in used_codes:
        code_sym.setdefault(orfmap[c], []).append(c)
    n_iso = np.asarray([len(code_sym.get(s, [])) for s in scored])
    n_iso50 = np.asarray([sum(1 for c in code_sym.get(s, [])
                              if cells_code[used_codes.index(c)] >= L.MIN_CELLS)
                          for s in scored])
    orf_mag = np.linalg.norm(R, axis=1)
    kd_mag = np.linalg.norm(Dsub, axis=1)
    # batch-1 fraction per scored symbol (obs metadata, free) — a declared covariate
    import h5py
    with h5py.File(L.ATLAS, "r") as _f:
        _codes = _f["obs/perturbation/codes"][:].astype(np.int64)
        _batch = _f["obs/batch/codes"][:].astype(np.float64)
    _c2s = L.read_orf_map()
    _sym_of_code = np.full(int(_codes.max()) + 1, -1, np.int64)
    _sidx = {s_: i for i, s_ in enumerate(scored)}
    for _c, _s in _c2s.items():
        if _s in _sidx and _c < _sym_of_code.size:
            _sym_of_code[_c] = _sidx[_s]
    _g = _sym_of_code[_codes]
    _ok = _g >= 0
    batch1 = (np.bincount(_g[_ok], weights=_batch[_ok], minlength=nT)
              / np.maximum(np.bincount(_g[_ok], minlength=nT), 1))[:nT]
    del _codes, _batch, _g, _ok

    # 279 / 290 features
    zn = {}
    with open(L.ZNF279) as fh:
        for r in csv.DictReader(fh):
            zn[r["gene"]] = r

    def col279(name):
        return np.asarray([float(zn[s][name]) if s in zn and zn[s][name] not in ("", "nan")
                           else np.nan for s in scored])
    C_held, r_bar = col279("C_held"), col279("r_bar")
    d290 = np.load(L.PERSRC290, allow_pickle=True)
    keys290 = np.asarray(d290["keys"]).astype(str).tolist()
    g290 = {g: i for i, g in enumerate(np.asarray(d290["genes"]).astype(str).tolist())}
    qi = keys290.index("gse295214_h1")

    def col290(name, agg=False):
        M = d290[name]
        vals = []
        for s in scored:
            if s not in g290:
                vals.append(np.nan); continue
            j = g290[s]
            with np.errstate(invalid="ignore"):
                vals.append(float(np.nanmean(M[:, j])) if agg else float(M[qi, j]))
        return np.asarray(vals)
    gen_pc20 = col290("pc20")
    gen_pc20_all = col290("pc20", agg=True)
    breadth_h1 = col290("breadth")
    breadth_all = col290("breadth", agg=True)
    mag290 = col290("mag")

    # ------------------------------------------------------------------ C6a
    idx_p = np.flatnonzero(pri)
    Sp = Sg[idx_p]
    med = float(np.median(Sp))
    lo_, hi_, _ = L.boot_ci(lambda ii: float(np.median(Sp[ii])), len(Sp), 20000, L.SEED)
    sh_neg = float((Sp < -0.05).mean()); sh_pos = float((Sp > 0.05).mean())
    rng = np.random.default_rng(L.SEED)
    null = np.empty(200)
    for t in range(200):
        while True:
            perm = rng.permutation(len(idx_p))
            if not (perm == np.arange(len(idx_p))).any():
                break
        null[t] = np.median([cos_pair(idx_p[a], idx_p[perm[a]], G, Aown, Bown, na2, nb2,
                                      own, owncol_pos) for a in range(len(idx_p))])
    c6a = dict(n=len(Sp), median=med, ci=[lo_, hi_],
               share_below_m05=sh_neg, share_above_p05=sh_pos,
               n_below_m05=int((Sp < -0.05).sum()), n_above_p05=int((Sp > 0.05).sum()),
               share_negative=float((Sp < 0).mean()),
               null_median=float(np.median(null)), null_sd=float(null.std()),
               null_frac_le_obs=float((null <= med).mean()),
               z=float((med - null.mean()) / null.std()) if null.std() > 0 else np.nan)
    c6a["pass_"] = bool(med < 0 and hi_ < 0 and sh_neg > sh_pos)
    P(f"\nC6a  n {c6a['n']}: median S {med:+.4f} [{lo_:+.4f}, {hi_:+.4f}]; "
      f"share < −0.05 {sh_neg:.1%} ({c6a['n_below_m05']}) vs > +0.05 {sh_pos:.1%} "
      f"({c6a['n_above_p05']}); {c6a['share_negative']:.1%} negative; null median "
      f"{c6a['null_median']:+.4f} ± {c6a['null_sd']:.4f} (z {c6a['z']:+.1f}, "
      f"{c6a['null_frac_le_obs']:.0%} of 200 draws ≤ observed) "
      f"-> {'PASS' if c6a['pass_'] else 'FAIL'}")

    # ------------------------------------------------------------------ C6b / C6c
    def rho_cell(feat, name, minus=True):
        y = -Sg if minus else Sg
        ok = pri & np.isfinite(feat) & np.isfinite(y)
        ii = np.flatnonzero(ok)
        r = L.spearman(feat[ii], y[ii])
        lo2, hi2, ple = L.boot_ci(lambda k: L.spearman(feat[ii[k]], y[ii[k]]),
                                  len(ii), 10000, L.SEED)
        return dict(name=name, n=len(ii), rho=r, ci=[lo2, hi2], p_le0=ple,
                    pass_=bool(np.isfinite(r) and np.isfinite(lo2) and r > 0 and lo2 > 0))
    c6b = rho_cell(C_held, "C_held")
    c6c = rho_cell(gen_pc20, "gen_pc20 (gse295214_h1)")
    P(f"C6b  n {c6b['n']}: Spearman(C_held, −S) {c6b['rho']:+.4f} "
      f"[{c6b['ci'][0]:+.4f}, {c6b['ci'][1]:+.4f}] -> {'PASS' if c6b['pass_'] else 'FAIL'}")
    P(f"C6c  n {c6c['n']}: Spearman(gen_pc20, −S) {c6c['rho']:+.4f} "
      f"[{c6c['ci'][0]:+.4f}, {c6c['ci'][1]:+.4f}] -> {'PASS' if c6c['pass_'] else 'FAIL'}")

    base = {}
    for feat, name in [(kd_mag, "kd_mag"), (mag290, "mag290_h1"), (orf_mag, "orf_mag"),
                       (rel, "orf_reliability"), (r_bar, "r_bar"), (cells, "orf_cells"),
                       (umi, "orf_umi_per_cell"), (n_iso.astype(float), "n_isoforms"),
                       (gen_pc20_all, "gen_pc20_mapwide"), (breadth_h1, "breadth_h1"),
                       (breadth_all, "breadth_mapwide"), (batch1, "orf_batch1_fraction")]:
        base[name] = rho_cell(feat, name)
        P(f"   baseline {name:22s} rho(−S) {base[name]['rho']:+.4f} "
          f"[{base[name]['ci'][0]:+.4f}, {base[name]['ci'][1]:+.4f}]  n {base[name]['n']}")
    Z = np.column_stack([kd_mag, orf_mag, rel])
    okp = pri & np.isfinite(C_held) & np.isfinite(Z).all(axis=1)
    part_b = L.partial_spearman(C_held[okp], -Sg[okp], Z[okp])
    okc = pri & np.isfinite(gen_pc20) & np.isfinite(Z).all(axis=1)
    part_c = L.partial_spearman(gen_pc20[okc], -Sg[okc], Z[okc])
    P(f"   partial (given kd_mag, orf_mag, reliability): C6b {part_b:+.4f} (n {int(okp.sum())}), "
      f"C6c {part_c:+.4f} (n {int(okc.sum())})")

    # ------------------------------------------------------------------ secondaries
    sec: dict = {}
    absS = np.abs(Sg)
    sec["abs_median"] = float(np.median(absS[pri]))
    okA = pri & np.isfinite(C_held)
    iiA = np.flatnonzero(okA)
    sec["abs_C_held"] = dict(rho=L.spearman(C_held[iiA], absS[iiA]), n=len(iiA))
    okB = pri & np.isfinite(gen_pc20)
    iiB = np.flatnonzero(okB)
    sec["abs_gen_pc20"] = dict(rho=L.spearman(gen_pc20[iiB], absS[iiB]), n=len(iiB))
    sc = np.flatnonzero(~pri)
    if sc.size:
        Ss = Sg[sc]
        lo3, hi3, _ = L.boot_ci(lambda ii: float(np.median(Ss[ii])), len(Ss), 20000, L.SEED)
        sec["cohort_50_99"] = dict(n=len(Ss), median=float(np.median(Ss)), ci=[lo3, hi3],
                                   share_neg=float((Ss < 0).mean()))
        P(f"   secondary cohort (50–99 cells) n {len(Ss)}: median S "
          f"{np.median(Ss):+.4f} [{lo3:+.4f}, {hi3:+.4f}]")
    sec["disattenuated_median"] = float(np.median(
        Sp / np.sqrt(np.clip(rel[idx_p], 1e-3, None))))
    sec["frac_positive"] = float((Sp > 0).mean())
    sec["median_minus_null"] = float(c6a["median"] - c6a["null_median"])

    # per-isoform
    cpos = {c: i for i, c in enumerate(used_codes)}
    P_code = L.profile(acc_code, lib_code[:, None])
    iso_rows, iso_pairs = [], []
    for i, s in enumerate(scored):
        cs = [c for c in code_sym.get(s, []) if cells_code[cpos[c]] >= L.MIN_CELLS]
        if len(cs) < 2 or not pri[i]:
            continue
        vs = []
        for c in cs:
            a = P_code[cpos[c]] - mu_T
            x, y = a, B[i]
            num, a2, b2 = float(x @ y), float(x @ x), float(y @ y)
            if own[i] >= 0:
                k = own[i]
                num -= x[k] * y[k]; a2 -= x[k] ** 2; b2 -= y[k] ** 2
            d = np.sqrt(max(a2, 0) * max(b2, 0))
            sc_ = num / d if d > 0 else np.nan
            iso_rows.append(dict(symbol=s, code=int(c), cells=int(cells_code[cpos[c]]), S=sc_))
            vs.append((c, a))
        for u in range(len(vs)):
            for w in range(u + 1, len(vs)):
                x, y = vs[u][1], vs[w][1]
                k = own[i]
                num, a2, b2 = float(x @ y), float(x @ x), float(y @ y)
                if k >= 0:
                    num -= x[k] * y[k]; a2 -= x[k] ** 2; b2 -= y[k] ** 2
                iso_pairs.append(num / np.sqrt(max(a2, 0) * max(b2, 0)))
    if iso_rows:
        iso_S = np.asarray([r["S"] for r in iso_rows])
        sec["per_isoform"] = dict(n_rows=len(iso_rows),
                                  n_tfs=len({r["symbol"] for r in iso_rows}),
                                  median_S=float(np.median(iso_S)),
                                  share_neg=float((iso_S < 0).mean()),
                                  isoform_agreement_median=float(np.median(iso_pairs)),
                                  n_pairs=len(iso_pairs))
        P(f"   per-isoform: {len(iso_rows)} ORF rows over "
          f"{sec['per_isoform']['n_tfs']} TFs, median S {np.median(iso_S):+.4f}; "
          f"isoform agreement median cos {np.median(iso_pairs):+.4f} ({len(iso_pairs)} pairs)")

    order = np.argsort(np.where(pri, Sg, np.nan))
    top_neg = [(scored[i], float(Sg[i]), int(cells[i])) for i in order[:20]]
    top_pos = [(scored[i], float(Sg[i]), int(cells[i]))
               for i in order[np.isfinite(Sg[order]) & pri[order]][::-1][:20]]
    sec["top_mirrored"] = top_neg
    sec["top_aligned"] = top_pos

    # ------------------------------------------------------------------ outputs
    with open(L.WORK / f"symmetry302{tag}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["symbol", "cohort", "orf_cells", "orf_umi_per_cell", "n_isoforms",
                    "n_isoforms_ge50", "orf_reliability", "S_g", "C_held", "kd_mag", "mag290_h1",
                    "gen_pc20_h1", "gen_pc20_mapwide", "breadth_h1", "orf_mag", "r_bar",
                    "batch1_fraction", "own_gene_on_axis", "orf_own_gene", "kd_own_gene"])
        for i, s in enumerate(scored):
            w.writerow([s, "primary" if pri[i] else "secondary", int(cells[i]),
                        f"{umi[i]:.1f}", int(n_iso[i]), int(n_iso50[i]), f"{rel[i]:.6f}",
                        f"{Sg[i]:.6f}",
                        "" if not np.isfinite(C_held[i]) else f"{C_held[i]:.6f}",
                        f"{kd_mag[i]:.6f}",
                        "" if not np.isfinite(mag290[i]) else f"{mag290[i]:.6f}",
                        "" if not np.isfinite(gen_pc20[i]) else f"{gen_pc20[i]:.6f}",
                        "" if not np.isfinite(gen_pc20_all[i]) else f"{gen_pc20_all[i]:.6f}",
                        "" if not np.isfinite(breadth_h1[i]) else f"{breadth_h1[i]:.6f}",
                        f"{orf_mag[i]:.6f}",
                        "" if not np.isfinite(r_bar[i]) else f"{r_bar[i]:.6f}",
                        f"{batch1[i]:.4f}",
                        int(own[i] >= 0),
                        "" if own[i] < 0 else f"{orf_own[i]:.6f}",
                        "" if own[i] < 0 else f"{kd_own[i]:.6f}"])
    if iso_rows:
        with open(L.WORK / f"isoform302{tag}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["symbol", "code", "cells", "S"])
            for r in iso_rows:
                w.writerow([r["symbol"], r["code"], r["cells"], f"{r['S']:.6f}"])

    out.update(n_axis=nA, n_scored=nT, KA_F=ka_f, KA_G=ka_g, C6a=c6a, C6b=c6b, C6c=c6c,
               baselines=base, partial_C6b=part_b, partial_C6c=part_c, secondaries=sec,
               depth=dict(cells_median=float(np.median(cells[pri])),
                          cells_q10=float(np.percentile(cells[pri], 10)),
                          umi_median=float(np.median(umi[pri])),
                          rel_median=float(np.median(rel[pri]))),
               rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, log=LOG)
    (L.WORK / f"stats302{tag}.json").write_text(json.dumps(out, indent=1, default=float))

    def verdict(c):
        return "PASS" if c["pass_"] else "FAIL"
    P(f"\n**C6a (overexpression reverses knockdown, hESC): {verdict(c6a)} "
      f"(median S {c6a['median']:+.4f} [{c6a['ci'][0]:+.4f}, {c6a['ci'][1]:+.4f}]; "
      f"share < −0.05 {c6a['share_below_m05']:.1%} vs > +0.05 {c6a['share_above_p05']:.1%}; "
      f"null {c6a['null_median']:+.4f} ± {c6a['null_sd']:.4f}); "
      f"C6b (conservation predicts mirroring): {verdict(c6b)} (ρ {c6b['rho']:+.4f} "
      f"[{c6b['ci'][0]:+.4f}, {c6b['ci'][1]:+.4f}]; partial {part_b:+.4f}); "
      f"C6c (geometry sign): {verdict(c6c)} (ρ {c6c['rho']:+.4f} "
      f"[{c6c['ci'][0]:+.4f}, {c6c['ci'][1]:+.4f}])**")
    P(f"peak RSS {out['rss_gib']:.2f} GiB")


if __name__ == "__main__":
    sys.exit(main())

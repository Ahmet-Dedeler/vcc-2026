"""295 item 4 — the statistics of claim C5.  Reads only artefacts already written: the access
table (`access295.csv`), the map (`held279.npz`, `persrc290.npz`, `pred284.csv`) and 292's
metadata.  No reserve file is opened again.

Statistical code (spearman, block bootstrap, partial Spearman) is COPIED from
`experiments/291-f13-independent-surfaces/lib291.py` lines 73-78, 198-248, per that file's own
"copied, never imported across experiments" convention.

Writes work/reserve295.csv (the per-gene table the mandate asks for) and work/stats295.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/295-reserve-r3-claim"
WORK = HERE / "work"

HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
PERSRC = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"
PRED284 = ROOT / "experiments/284-conservation-family-split/work/pred284.csv"

SEED = 20260906
N_BOOT = 2000
N_NULL = 200

CRISPRI = {"AdamsonWeissman2016_GSM2406681_10X010", "DixitRegev2016",
           "TianKampmann2019_day7neuron"}
TIAN = "TianKampmann2019_day7neuron"


# ---- lib291.py:73-78 --------------------------------------------------------------------------
def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[m], b[m]).statistic)


# ---- lib291.py:198-207 ------------------------------------------------------------------------
def boot_rho(a, b, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    vals = [spearman(a[sel], b[sel]) for sel in
            (np.concatenate([idx_of[uniq[p]] for p in rng.integers(0, len(uniq), size=len(uniq))])
             for _ in range(n))]
    return dict(point=spearman(a, b), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)),
                n=int((np.isfinite(a) & np.isfinite(b)).sum()))


def boot_rho_rows(a, b, rng, n: int = N_BOOT) -> dict:
    """plain row (gene) bootstrap."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[m], np.asarray(b)[m]
    vals = [spearman(a[i], b[i]) for i in
            (rng.integers(0, a.size, size=a.size) for _ in range(n))]
    return dict(point=spearman(a, b), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)), n=int(a.size))


# ---- lib291.py:210-232 ------------------------------------------------------------------------
def _rank(x):
    x = np.asarray(x, float)
    r = np.full(len(x), np.nan)
    m = np.isfinite(x)
    r[m] = stats.rankdata(x[m])
    return r


def partial_spearman(a, b, Z) -> float:
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    m = np.isfinite(a) & np.isfinite(b) & np.isfinite(Z).all(axis=1)
    if m.sum() < 5:
        return float("nan")
    ra, rb = _rank(np.asarray(a)[m]), _rank(np.asarray(b)[m])
    RZ = np.column_stack([_rank(Z[m, j]) for j in range(Z.shape[1])] + [np.ones(int(m.sum()))])
    ea = ra - RZ @ np.linalg.lstsq(RZ, ra, rcond=None)[0]
    eb = rb - RZ @ np.linalg.lstsq(RZ, rb, rcond=None)[0]
    if ea.std() < 1e-12 or eb.std() < 1e-12:
        return float("nan")
    return float(np.corrcoef(ea, eb)[0, 1])


def boot_partial(a, b, Z, groups, rng, n: int = N_BOOT) -> dict:
    groups = np.asarray(groups)
    uniq = np.array(sorted(set(groups.tolist())))
    idx_of = {g: np.flatnonzero(groups == g) for g in uniq}
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    vals = []
    for _ in range(n):
        pick = rng.integers(0, len(uniq), size=len(uniq))
        sel = np.concatenate([idx_of[uniq[p]] for p in pick])
        vals.append(partial_spearman(a[sel], b[sel], Z[sel]))
    return dict(point=partial_spearman(a, b, Z), lo=float(np.nanpercentile(vals, 2.5)),
                hi=float(np.nanpercentile(vals, 97.5)), sd=float(np.nanstd(vals)))


def fmt(d) -> str:
    return f"{d['point']:+.4f} [{d['lo']:+.4f}, {d['hi']:+.4f}]"


# ------------------------------------------------------------------------------------- the table
def build_table() -> pd.DataFrame:
    acc = pd.read_csv(WORK / "access295.csv")
    meta = json.loads((WORK / "meta295.json").read_text())
    z = np.load(HELD279, allow_pickle=True)
    genes = [str(g) for g in z["genes"]]
    gi = {g: i for i, g in enumerate(genes)}
    keys = [str(k) for k in z["keys"]]
    lines = np.asarray([str(l) for l in z["lines"]])
    Ccross, pairs = z["Ccross"], z["pairs_cross"]
    ps = np.load(PERSRC, allow_pickle=True)
    assert [str(g) for g in ps["genes"]] == genes and [str(k) for k in ps["keys"]] == keys
    pc20, breadth, magsrc = ps["pc20"], ps["breadth"], ps["mag"]
    p284 = pd.read_csv(PRED284).set_index("gene")

    out = []
    for _, r in acc.iterrows():
        g, stem = r["gene"], r["file"]
        i = gi[g]
        excl = set(meta["files"][stem]["excluded_banks"])
        rows_ok = np.asarray([k not in excl for k in keys])
        fl = meta["files"][stem]["line"]
        pm = np.asarray([(lines[a] != fl) and (lines[b] != fl) for a, b in pairs])
        with np.errstate(invalid="ignore"):
            c_excl = float(np.nanmean(Ccross[pm, i])) if pm.any() else np.nan
            gp = float(np.nanmean(pc20[rows_ok, i]))
            br = float(np.nanmean(breadth[rows_ok, i]))
            mg = np.nanmean(magsrc[rows_ok, i])
            lm = float(np.log(mg)) if np.isfinite(mg) and mg > 0 else np.nan
        rec = dict(r)
        rec.update(C_held=float(z["c_held"][i]), C_held_excl=c_excl,
                   n_map_pairs=int(z["n_pairs"][i]), mag=float(z["mag"][i]),
                   gen_pc20=gp, breadth=br, logmag=lm,
                   M_ann=float(p284.loc[g, "pred_oof_raw"]),
                   r_bar=float(p284.loc[g, "r_bar"]),
                   modality=("CRISPRi" if stem in CRISPRI else "Cas9-KO"))
        out.append(rec)
    df = pd.DataFrame(out)
    # primary unit: one row per gene, the file with the MOST cells (ties -> first file name)
    df = df.sort_values(["gene", "cells", "file"], ascending=[True, False, True])
    df["primary"] = ~df.duplicated("gene", keep="first")
    return df.reset_index(drop=True)


def cell(df: pd.DataFrame, xcol: str, rng, label: str) -> dict:
    a = df[xcol].to_numpy(float)
    b = df["A"].to_numpy(float)
    fb = boot_rho(a, b, df["file"].to_numpy(), rng)
    gb = boot_rho_rows(a, b, rng)
    return dict(label=label, x=xcol, n=int(np.isfinite(a * b).sum()),
                rho=fb["point"], file_block=fb, gene_boot=gb)


def main() -> None:
    rng = np.random.default_rng(SEED)
    df = build_table()
    df.to_csv(WORK / "reserve295.csv", index=False)
    P = df[df["primary"]].reset_index(drop=True)
    print(f"rows {len(df)} ((file,gene) pairs); primary genes {len(P)}", flush=True)

    res: dict = {"n_pairs": int(len(df)), "n_primary": int(len(P))}

    # ---------------- per-file table -----------------------------------------------------------
    tab = []
    for stem, sub in df.groupby("file"):
        tab.append(dict(file=stem, modality=sub["modality"].iloc[0], line=sub["line"].iloc[0],
                        n=len(sub), median_cells=float(sub["cells"].median()),
                        umi_per_cell=float(sub["umi_file"].iloc[0]),
                        median_r_split=float(sub["r_split"].median()),
                        median_de=float(sub["de"].median()),
                        median_A=float(sub["A"].median()),
                        share_A_le0=float((sub["A"] <= 0).mean()),
                        median_C_held=float(sub["C_held"].median()),
                        n_primary=int(sub["primary"].sum())))
    res["per_file"] = tab
    for t in tab:
        print(f"  {t['file']:<40} n {t['n']:3d} cells {t['median_cells']:7.0f} "
              f"umi {t['umi_per_cell']:7.0f} r {t['median_r_split']:+.3f} DE {t['median_de']:6.0f} "
              f"A {t['median_A']:+.4f}", flush=True)

    # ---------------- C5a ------------------------------------------------------------------------
    print("\n=== C5a: rho(C_held, A_g) ===", flush=True)
    c5a = cell(P, "C_held", rng, "C5a primary (174 genes, pooled, Tian included)")
    c5a["pass"] = bool(c5a["rho"] >= 0.15 and c5a["file_block"]["lo"] > 0)
    print(f"  POOLED   rho {c5a['rho']:+.4f}  file-block {fmt(c5a['file_block'])}  "
          f"gene-boot {fmt(c5a['gene_boot'])}  n={c5a['n']}  -> "
          f"{'PASS' if c5a['pass'] else 'FAIL'}", flush=True)
    res["C5a"] = c5a

    base = {}
    for col in ["M_ann", "mag", "r_split", "de", "r_bar", "C_held_excl", "logmag", "cells"]:
        base[col] = cell(P, col, rng, f"baseline {col}")
        print(f"  baseline {col:<12} rho {base[col]['rho']:+.4f}  "
              f"file-block {fmt(base[col]['file_block'])}", flush=True)
    res["C5a_baselines"] = base

    pa = boot_partial(P["C_held"].to_numpy(float), P["A"].to_numpy(float),
                      np.column_stack([P["r_split"].to_numpy(float), P["mag"].to_numpy(float)]),
                      P["file"].to_numpy(), rng)
    res["C5a_partial_given_rel_mag"] = pa
    print(f"  partial(C_held, A | r_split, mag) {fmt(pa)}", flush=True)
    pa2 = boot_partial(P["C_held"].to_numpy(float), P["A"].to_numpy(float),
                       np.column_stack([P["r_split"].to_numpy(float), P["mag"].to_numpy(float),
                                        np.log1p(P["cells"].to_numpy(float)),
                                        np.log1p(P["umi_file"].to_numpy(float))]),
                       P["file"].to_numpy(), rng)
    res["C5a_partial_plus_depth"] = pa2
    print(f"  partial(+ log cells, log UMI/cell)  {fmt(pa2)}", flush=True)

    # null
    nrng = np.random.default_rng(SEED + 1)
    a = P["C_held"].to_numpy(float); b = P["A"].to_numpy(float)
    nulls = np.asarray([spearman(a, nrng.permutation(b)) for _ in range(N_NULL)])
    res["C5a_null"] = dict(mean=float(nulls.mean()), sd=float(nulls.std()),
                           n_ge_obs=int((nulls >= c5a["rho"]).sum()), n=N_NULL,
                           z=float((c5a["rho"] - nulls.mean()) / nulls.std()))
    print(f"  null {nulls.mean():+.4f} +- {nulls.std():.4f}, "
          f"{res['C5a_null']['n_ge_obs']}/{N_NULL} reach the observed, "
          f"z {res['C5a_null']['z']:+.1f}", flush=True)

    # ---------------- C5b ------------------------------------------------------------------------
    print("\n=== C5b: rho(gen_pc20, A_g) ===", flush=True)
    c5b = cell(P, "gen_pc20", rng, "C5b primary")
    c5b["pass"] = bool(c5b["file_block"]["lo"] > 0)
    print(f"  gen_pc20 rho {c5b['rho']:+.4f}  file-block {fmt(c5b['file_block'])}  "
          f"gene-boot {fmt(c5b['gene_boot'])}  -> {'PASS' if c5b['pass'] else 'FAIL'}", flush=True)
    br = cell(P, "breadth", rng, "breadth (reported)")
    print(f"  breadth  rho {br['rho']:+.4f}  file-block {fmt(br['file_block'])}", flush=True)
    res["C5b"] = c5b
    res["C5b_breadth"] = br

    # ---------------- strata ----------------------------------------------------------------------
    print("\n=== strata (declared) ===", flush=True)
    strata = {}
    for name, sub in [("CRISPRi", P[P.modality == "CRISPRi"]),
                      ("Cas9-KO", P[P.modality == "Cas9-KO"]),
                      ("pooled_no_Tian", P[P.file != TIAN]),
                      ("core7_only", P[P.file != TIAN])]:
        if name == "core7_only":
            continue
        s = {}
        for col, tag in [("C_held", "C5a"), ("gen_pc20", "C5b")]:
            s[tag] = cell(sub, col, rng, f"{name} {tag}")
        strata[name] = dict(n=int(len(sub)), **s)
        print(f"  {name:<16} n {len(sub):3d}  C5a {s['C5a']['rho']:+.4f} "
              f"{fmt(s['C5a']['file_block'])}   C5b {s['C5b']['rho']:+.4f} "
              f"{fmt(s['C5b']['file_block'])}", flush=True)
    res["strata"] = strata

    # ---------------- labelled secondaries ---------------------------------------------------------
    print("\n=== labelled secondaries ===", flush=True)
    sec = {}
    sec["all_205_pairs_C5a"] = cell(df, "C_held", rng, "all (file,gene) pairs")
    sec["all_205_pairs_C5b"] = cell(df, "gen_pc20", rng, "all (file,gene) pairs")
    # on-target column KEPT
    P2 = P.copy(); P2["A"] = P2["A_onT"]
    sec["ontarget_kept_C5a"] = cell(P2, "C_held", rng, "on-target column kept")
    sec["ontarget_kept_C5b"] = cell(P2, "gen_pc20", rng, "on-target column kept")
    # disattenuated outcome
    P3 = P.copy()
    P3["A"] = P3["A"] / np.sqrt(np.maximum(P3["r_split"], 1e-6))
    P3.loc[P3["r_split"] <= 0, "A"] = np.nan
    sec["disattenuated_C5a"] = cell(P3, "C_held", rng, "A / sqrt(r_split), r>0 only")
    sec["disattenuated_C5b"] = cell(P3, "gen_pc20", rng, "A / sqrt(r_split), r>0 only")
    # within-file ranks (removes any between-file offset)
    P4 = P.copy()
    P4["A"] = P4.groupby("file")["A"].transform(lambda s: stats.rankdata(s) / len(s))
    P4["C_held_r"] = P4.groupby("file")["C_held"].transform(lambda s: stats.rankdata(s) / len(s))
    P4["gen_pc20_r"] = P4.groupby("file")["gen_pc20"].transform(lambda s: stats.rankdata(s) / len(s))
    sec["within_file_ranks_C5a"] = cell(P4, "C_held_r", rng, "within-file ranks")
    sec["within_file_ranks_C5b"] = cell(P4, "gen_pc20_r", rng, "within-file ranks")
    # C_held_excl as the predictor for C5a
    sec["C_held_excl_C5a"] = base["C_held_excl"]
    for k, v in sec.items():
        print(f"  {k:<26} rho {v['rho']:+.4f}  file-block {fmt(v['file_block'])}  n {v['n']}",
              flush=True)
    res["secondaries"] = sec

    # ---------------- descriptive -------------------------------------------------------------------
    print("\n=== descriptive (post-hoc, labelled) ===", flush=True)
    d = {}
    top = P.sort_values("A", ascending=False)
    d["best"] = top.head(10)[["gene", "file", "A", "r_split", "cells", "C_held"]].to_dict("records")
    d["worst"] = top.tail(10)[["gene", "file", "A", "r_split", "cells", "C_held"]].to_dict("records")
    for pair, tag in [(("DatlingerBock2017", "DatlingerBock2021"), "datlinger_2017_vs_2021"),
                      (("PapalexiSatija2021_eccite_RNA",
                        "PapalexiSatija2021_eccite_arrayed_RNA"), "papalexi_pooled_vs_arrayed")]:
        x = df[df.file == pair[0]].set_index("gene")
        y = df[df.file == pair[1]].set_index("gene")
        common = sorted(set(x.index) & set(y.index))
        if not common:
            continue
        d[tag] = dict(n=len(common),
                      rho_A=spearman(x.loc[common, "A"], y.loc[common, "A"]),
                      median_A_x=float(x.loc[common, "A"].median()),
                      median_A_y=float(y.loc[common, "A"].median()),
                      median_r_x=float(x.loc[common, "r_split"].median()),
                      median_r_y=float(y.loc[common, "r_split"].median()),
                      genes=common)
        print(f"  {tag}: n {len(common)}, rho(A_x, A_y) {d[tag]['rho_A']:+.4f}, "
              f"median A {d[tag]['median_A_x']:+.4f} vs {d[tag]['median_A_y']:+.4f}", flush=True)
    cross = df[df.gene.isin(["EGR1", "IRF1", "SOCS1"])][["gene", "file", "A", "r_split", "cells"]]
    d["cross_study_genes"] = cross.to_dict("records")
    # depth associations (rule 9, reported never conditioned)
    d["rho_A_cells"] = spearman(P["cells"], P["A"])
    d["rho_A_rsplit"] = spearman(P["r_split"], P["A"])
    d["rho_A_umi"] = spearman(P["umi_file"], P["A"])
    d["rho_A_de"] = spearman(P["de"], P["A"])
    d["median_A"] = float(P["A"].median())
    d["share_A_le0"] = float((P["A"] <= 0).mean())
    print(f"  median A {d['median_A']:+.4f}, share <=0 {d['share_A_le0']:.3f}; "
          f"rho(A, cells) {d['rho_A_cells']:+.3f}, rho(A, r_split) {d['rho_A_rsplit']:+.3f}, "
          f"rho(A, UMI/cell) {d['rho_A_umi']:+.3f}, rho(A, DE) {d['rho_A_de']:+.3f}", flush=True)
    res["descriptive"] = d

    (WORK / "stats295.json").write_text(json.dumps(res, indent=1, default=float))
    print("\nHEADLINE: C5a " + ("PASS" if c5a["pass"] else "FAIL")
          + f" (rho {c5a['rho']:+.4f} {fmt(c5a['file_block'])}); C5b "
          + ("PASS" if c5b["pass"] else "FAIL")
          + f" (gen_pc20 {c5b['rho']:+.4f} {fmt(c5b['file_block'])})", flush=True)


if __name__ == "__main__":
    main()

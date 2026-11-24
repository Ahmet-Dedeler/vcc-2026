"""281 — assemble the gene-level FEATURE table for the 17,106 held genes of 279 (+ the 300 panel targets).

Nothing here touches the target beyond copying it out of 279's deliverable; no model is fitted.
Every feature is build-time observable: none is derived from a measured perturbation response.

Outputs (all inside this experiment):
  work/feat281.csv     one row per gene, features + target + split keys
  work/feat281.json    provenance, coverage, known-answer rows
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

T0 = time.time()
LOG: list[str] = []


def P(msg: str) -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    print(line, flush=True)
    LOG.append(line)


out: dict = {"sources": {}, "ka": {}, "coverage": {}}

# ---------------------------------------------------------------------------------------------
# 1. the target: 279's held conservation map, and 276's panel map
# ---------------------------------------------------------------------------------------------
Z279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
Z276 = ROOT / "experiments/276-conservation-map/work/conservation276.csv"
z = pd.read_csv(Z279)
P(f"279 held map: {z.shape[0]} genes x {z.shape[1]} cols  ({Z279.relative_to(ROOT)})")

# KA-1: reproduce 279 RESULTS §1.1 / §3 headline numbers straight off the deliverable.
ka1 = dict(
    n=int(len(z)),
    median_C_held=float(np.median(z["C_held"])),
    share_le0=float((z["C_held"] <= 0).mean()),
    znf_delta_median=float(np.median(z.loc[z["class"] == "ZNF", "C_held"])
                          - np.median(z.loc[z["class"] == "NONZNF", "C_held"])),
    n_znf=int((z["class"] == "ZNF").sum()),
)
ka1["pass"] = bool(ka1["n"] == 17106 and abs(ka1["median_C_held"] - 0.0186) < 5e-4
                   and abs(ka1["share_le0"] - 0.306) < 1e-3
                   and abs(ka1["znf_delta_median"] + 0.01047) < 1e-4 and ka1["n_znf"] == 696)
out["ka"]["KA-1_279_headline"] = ka1
P(f"KA-1 (279 headline): n {ka1['n']}, median C_held {ka1['median_C_held']:+.6f} (req +0.0186), "
  f"share<=0 {ka1['share_le0']:.3f} (req 0.306), ZNF dmed {ka1['znf_delta_median']:+.5f} "
  f"(req -0.01047), n_ZNF {ka1['n_znf']} (req 696) -> {'PASS' if ka1['pass'] else 'FAIL'}")
if not ka1["pass"]:
    sys.exit("KA-1 failed — aborting (AGENTS.md §known-answer rule)")

p276 = pd.read_csv(Z276)
P(f"276 panel map: {p276.shape[0]} targets x {p276.shape[1]} cols")

df = pd.DataFrame({
    "gene": z["gene"].astype(str),
    "set": "held",
    "y_raw": z["C_held"].astype(float),
    "y_rel": z["C_rel"].astype(float),
    "r_bar": z["r_bar"].astype(float),
    "n_pairs": z["n_pairs"].astype(float),
    "n_lines": z["n_lines"].astype(float),
    "znf_class": z["class"].astype(str),
})
panel = pd.DataFrame({
    "gene": p276["target"].astype(str),
    "set": "panel",
    "y_raw": p276["C_raw"].astype(float),
    "y_rel": p276["C_rel"].astype(float),
    "r_bar": np.nan,                       # 276 carries no per-target r_bar column
    "n_pairs": p276["n_pairs"].astype(float),
    "n_lines": np.nan,
    "znf_class": "PANEL",
})
# The panel's per-target reference reliability comes from 271's `pertarget271-<fold>.npz` `ref_rel`
# (the held line's own split-half residual reliability, R271 §3). PROGRAMME RULE 2: the
# `orion_hek293t` fold's split-half B TRUTHS on the 300 panel targets are reserved, so that fold is
# NOT read here — the panel's r_bar is the mean of `k562_gwps` and `orion_hct116` only.
rel_acc: dict[str, list[float]] = {}
for fold in ("k562_gwps", "orion_hct116"):
    zz271 = np.load(ROOT / f"experiments/271-donor-reliability-estimation/work/pertarget271-{fold}.npz",
                    allow_pickle=True)
    for t, r in zip(np.asarray(zz271["targets"]).astype(str), np.asarray(zz271["ref_rel"])):
        if np.isfinite(r):
            rel_acc.setdefault(str(t), []).append(float(r))
panel["r_bar"] = panel["gene"].map({k: float(np.mean(v)) for k, v in rel_acc.items()})
P(f"panel r_bar from 271 (k562_gwps + orion_hct116 only; hek293t reserved): "
  f"{panel['r_bar'].notna().sum()}/300 targets, median {np.nanmedian(panel['r_bar']):+.4f}")

df = pd.concat([df, panel], ignore_index=True)
P(f"combined table: {len(df)} rows ({(df['set'] == 'held').sum()} held, {(df['set'] == 'panel').sum()} panel)")
overlap = set(df.loc[df["set"] == "held", "gene"]) & set(df.loc[df["set"] == "panel", "gene"])
P(f"held/panel gene overlap: {len(overlap)} (must be 0 by 279's construction)")
out["coverage"]["held_panel_overlap"] = len(overlap)

# ---------------------------------------------------------------------------------------------
# 2. HGNC: family class, cytoband, ensembl id, uniprot ids
#    source: data/external/xcell/hgnc_protein_coding.txt (HGNC complete protein-coding set,
#    genenames.org; HGNC data are freely available, no licence restriction)
# ---------------------------------------------------------------------------------------------
HG = ROOT / "data/external/xcell/hgnc_protein_coding.txt"
hg = pd.read_csv(HG, sep="\t", low_memory=False,
                 usecols=["symbol", "gene_group", "location", "ensembl_gene_id", "uniprot_ids",
                          "locus_group"])
hg = hg.drop_duplicates("symbol").set_index("symbol")
out["sources"]["hgnc"] = dict(path=str(HG.relative_to(ROOT)), rows=int(len(hg)))
df["gene_group_all"] = df["gene"].map(hg["gene_group"]).fillna("")
df["gene_group1"] = df["gene_group_all"].str.split("|").str[0].fillna("")
df["cytoband"] = df["gene"].map(hg["location"]).fillna("")
df["ensembl_gene_id"] = df["gene"].map(hg["ensembl_gene_id"]).fillna("")
df["uniprot_ids"] = df["gene"].map(hg["uniprot_ids"]).fillna("")
P(f"HGNC: gene_group on {(df['gene_group1'] != '').sum()}/{len(df)}, "
  f"cytoband on {(df['cytoband'] != '').sum()}, ensembl on {(df['ensembl_gene_id'] != '').sum()}")

# ---------------------------------------------------------------------------------------------
# 3. CCLE expression breadth / level  (276's convention, analyse276.py:193-202)
#    source: data/external/ccle_broad_2025/ccle_cache.npz  (Broad CCLE / DepMap 2025 mRNA TPM,
#    1,600 samples x 19,209 genes, log2(TPM+1))
# ---------------------------------------------------------------------------------------------
CC = ROOT / "data/external/ccle_broad_2025/ccle_cache.npz"
zz = np.load(CC, allow_pickle=True)
cg = {str(g): i for i, g in enumerate(zz["genes"])}
X = zz["X"]
P(f"CCLE cache: {X.shape[0]} samples x {X.shape[1]} genes")
out["sources"]["ccle"] = dict(path=str(CC.relative_to(ROOT)), shape=list(X.shape))
breadth = np.full(len(df), np.nan)
ccle_med = np.full(len(df), np.nan)
ccle_sd = np.full(len(df), np.nan)
for i, g in enumerate(df["gene"].to_numpy()):
    j = cg.get(g)
    if j is not None:
        v = np.asarray(X[:, j], dtype=np.float64)
        breadth[i] = float((v > 1.0).mean())
        ccle_med[i] = float(np.median(v))
        ccle_sd[i] = float(np.std(v))
df["ccle_breadth"] = breadth
df["ccle_median_tpm"] = ccle_med
df["ccle_sd"] = ccle_sd
del X, zz

# KA-2: my CCLE reader must reproduce 276's own columns on the 300 panel targets.
pm = p276.set_index("target")
sel = df["set"] == "panel"
a = df.loc[sel, ["gene", "ccle_breadth", "ccle_median_tpm"]].set_index("gene")
b = pm.loc[a.index, ["ccle_breadth", "ccle_median_tpm"]]
d_breadth = float(np.nanmax(np.abs(a["ccle_breadth"].to_numpy() - b["ccle_breadth"].to_numpy())))
d_med = float(np.nanmax(np.abs(a["ccle_median_tpm"].to_numpy() - b["ccle_median_tpm"].to_numpy())))
ka2 = dict(max_abs_diff_breadth=d_breadth, max_abs_diff_median=d_med,
           pass_=bool(d_breadth < 1e-6 and d_med < 1e-6))
out["ka"]["KA-2_ccle_matches_276"] = ka2
P(f"KA-2 (CCLE reader vs conservation276.csv on 300 panel targets): "
  f"max|Δbreadth| {d_breadth:.2e}, max|Δmedian| {d_med:.2e} -> {'PASS' if ka2['pass_'] else 'FAIL'}")
if not ka2["pass_"]:
    sys.exit("KA-2 failed — aborting")

# ---------------------------------------------------------------------------------------------
# 4. DepMap 24Q2 Public gene effect (Chronos), downloaded by 269 (figshare 25880521)
#    source: experiments/269-utility-router-feasibility/work/depmap/CRISPRGeneEffect.csv
#    Terms: DepMap Public data, CC BY 4.0 (depmap.org).
# ---------------------------------------------------------------------------------------------
DM = ROOT / "experiments/269-utility-router-feasibility/work/depmap/CRISPRGeneEffect.csv"
P(f"reading DepMap gene effect ({DM.stat().st_size / 1e6:.0f} MB) ...")
dm = pd.read_csv(DM, index_col=0)
sym = pd.Index([re.sub(r"\s*\(\d+\)$", "", c) for c in dm.columns])
P(f"DepMap: {dm.shape[0]} models x {dm.shape[1]} genes; symbols unique {sym.is_unique}")
out["sources"]["depmap"] = dict(path=str(DM.relative_to(ROOT)), shape=list(dm.shape),
                                bytes=int(DM.stat().st_size))
E = dm.to_numpy(dtype=np.float32)
del dm
dm_mean = np.nanmean(E, axis=0)
dm_sd = np.nanstd(E, axis=0)
with np.errstate(invalid="ignore"):
    dm_frac_dep = np.nanmean(E < -0.5, axis=0)
    dm_frac_strong = np.nanmean(E < -1.0, axis=0)
dm_min = np.nanmin(E, axis=0)
del E
tab = pd.DataFrame(dict(dm_mean_effect=dm_mean, dm_sd_effect=dm_sd,
                        dm_frac_dependent=dm_frac_dep, dm_frac_strong=dm_frac_strong,
                        dm_min_effect=dm_min), index=sym)
tab = tab[~tab.index.duplicated()]
for c in tab.columns:
    df[c] = df["gene"].map(tab[c])
# DepMap-style "common essential": dependent (effect < -0.5) in >= 90 % of the models.
df["dm_common_essential"] = (df["dm_frac_dependent"] >= 0.90).astype(float)
df.loc[df["dm_mean_effect"].isna(), "dm_common_essential"] = np.nan
P(f"DepMap coverage: {df['dm_mean_effect'].notna().sum()}/{len(df)} genes; "
  f"common-essential {int(np.nansum(df['dm_common_essential']))}")

# KA-3: the file must be byte-for-byte the artifact 269 recorded — sha256 `d15514918130…`,
# 419,330,454 B, 1,150 models, all 300 panel targets covered. (R269's prose also says "17,931
# genes"; the file's header carries 18,443 gene columns. The digest and the byte count decide,
# and they match: R269's gene count is the number that is wrong, recorded here, not papered over.)
DM_SHA = "d155149181308f8e16b6a4677cde3e5af1ddef68ada1947486ae8c0bcb231452"
_h = hashlib.sha256()
with DM.open("rb") as fh:
    for chunk in iter(lambda: fh.read(1 << 22), b""):
        _h.update(chunk)
dm_sha = _h.hexdigest()
n_panel_cov = int(df.loc[df["set"] == "panel", "dm_mean_effect"].notna().sum())
ka3 = dict(sha256=dm_sha, sha256_required=DM_SHA,
           bytes=int(DM.stat().st_size), bytes_required=419330454,
           models=int(out["sources"]["depmap"]["shape"][0]),
           genes=int(out["sources"]["depmap"]["shape"][1]),
           genes_claimed_by_R269=17931,
           panel_targets_covered=n_panel_cov,
           pass_=bool(dm_sha == DM_SHA and DM.stat().st_size == 419330454
                      and out["sources"]["depmap"]["shape"][0] == 1150 and n_panel_cov == 300))
out["sources"]["depmap"]["sha256"] = dm_sha
out["ka"]["KA-3_depmap_matches_269"] = ka3
P(f"KA-3 (DepMap artifact vs R269): sha256 {dm_sha[:12]}… (req {DM_SHA[:12]}…), "
  f"{ka3['bytes']} B, {ka3['models']} models x {ka3['genes']} genes (R269's prose says 17,931 — "
  f"the header carries {ka3['genes']}), panel covered {n_panel_cov}/300 -> "
  f"{'PASS' if ka3['pass_'] else 'FAIL'}")
if not ka3["pass_"]:
    sys.exit("KA-3 failed — aborting")

# Replogle's essential library membership — 276's `essential` convention (analyse276.py:204-205).
ess = set((ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
df["replogle_essential"] = df["gene"].isin(ess).astype(float)
out["sources"]["replogle_essential"] = dict(path="data/donors/k562_essential/perts.txt", n=len(ess))
# KA-2b: 276's own `essential` column on the panel must agree exactly.
agree = int((df.loc[sel, "replogle_essential"].to_numpy()
             == pm.loc[a.index, "essential"].to_numpy().astype(float)).sum())
out["ka"]["KA-2b_essential_matches_276"] = dict(agree=agree, n=300, pass_=bool(agree == 300))
P(f"KA-2b (Replogle-essential flag vs 276's `essential`): {agree}/300 -> "
  f"{'PASS' if agree == 300 else 'FAIL'}")

# ---------------------------------------------------------------------------------------------
# 5. Regulons: CollecTRI (OmniPath) and DoRothEA A-C
#    sources: data/external/regulon/collectri_omnipath.tsv, data/external/omnipath-dorothea/dorothea-abc.tsv
# ---------------------------------------------------------------------------------------------
CT = ROOT / "data/external/regulon/collectri_omnipath.tsv"
ct = pd.read_csv(CT, sep="\t", usecols=["source_genesymbol", "target_genesymbol"])
ct_out = ct.groupby("source_genesymbol")["target_genesymbol"].nunique()
ct_in = ct.groupby("target_genesymbol")["source_genesymbol"].nunique()
df["collectri_regulon_size"] = df["gene"].map(ct_out).fillna(0.0)
df["collectri_indegree"] = df["gene"].map(ct_in).fillna(0.0)
df["is_tf_collectri"] = (df["collectri_regulon_size"] > 0).astype(float)
out["sources"]["collectri"] = dict(path=str(CT.relative_to(ROOT)), edges=int(len(ct)),
                                   tfs=int(ct_out.shape[0]))
P(f"CollecTRI: {len(ct)} edges, {ct_out.shape[0]} TFs; "
  f"{int(df['is_tf_collectri'].sum())} of our genes are TFs")

DO = ROOT / "data/external/omnipath-dorothea/dorothea-abc.tsv"
do = pd.read_csv(DO, sep="\t", usecols=["source_genesymbol", "target_genesymbol", "dorothea_level"])
do_out = do.groupby("source_genesymbol")["target_genesymbol"].nunique()
df["dorothea_regulon_size"] = df["gene"].map(do_out).fillna(0.0)
df["is_tf_dorothea"] = (df["dorothea_regulon_size"] > 0).astype(float)
out["sources"]["dorothea"] = dict(path=str(DO.relative_to(ROOT)), edges=int(len(do)),
                                  tfs=int(do_out.shape[0]))
P(f"DoRothEA A-C: {len(do)} edges, {do_out.shape[0]} TFs")

# ---------------------------------------------------------------------------------------------
# 6. Protein-complex membership
#    (a) hu.MAP 2.0  humap2_complexes_20200809.txt (Drew et al. 2021; humap2.proteincomplexes.org,
#        freely downloadable, no registration)  -> work/humap2_complexes.txt
#    (b) Reactome ComplexParticipantsPubMedIdentifiers_human.txt (reactome.org, CC0) -> work/
# ---------------------------------------------------------------------------------------------
HM = WORK / "humap2_complexes.txt"
hm = pd.read_csv(HM)
n_cx = {}
mx_cx = {}
for names in hm["genenames"].fillna(""):
    gs = [g for g in str(names).split() if g]
    for g in set(gs):
        n_cx[g] = n_cx.get(g, 0) + 1
        mx_cx[g] = max(mx_cx.get(g, 0), len(gs))
df["humap_n_complexes"] = df["gene"].map(n_cx).fillna(0.0)
df["humap_max_complex_size"] = df["gene"].map(mx_cx).fillna(0.0)
out["sources"]["humap2"] = dict(path="work/humap2_complexes.txt", complexes=int(len(hm)),
                                url="https://humap2.proteincomplexes.org/static/downloads/humap2/"
                                    "humap2_complexes_20200809.txt")
P(f"hu.MAP 2.0: {len(hm)} complexes; {(df['humap_n_complexes'] > 0).sum()} of our genes in >=1")

RX = WORK / "reactome_complexes_human.txt"
rx = pd.read_csv(RX, sep="\t")
up2n: dict[str, int] = {}
up2mx: dict[str, int] = {}
for parts in rx["participants"].fillna(""):
    accs = set(re.findall(r"UniProt:([A-Za-z0-9]+)", str(parts)))
    for u in accs:
        up2n[u] = up2n.get(u, 0) + 1
        up2mx[u] = max(up2mx.get(u, 0), len(accs))
def _first_up(s: str) -> str:
    s = str(s)
    return s.split("|")[0] if s else ""
ups = df["uniprot_ids"].map(_first_up)
df["reactome_n_complexes"] = ups.map(up2n).fillna(0.0)
df["reactome_max_complex_size"] = ups.map(up2mx).fillna(0.0)
out["sources"]["reactome"] = dict(path="work/reactome_complexes_human.txt", complexes=int(len(rx)),
                                  url="https://reactome.org/download/current/"
                                      "ComplexParticipantsPubMedIdentifiers_human.txt")
P(f"Reactome: {len(rx)} human complexes; {(df['reactome_n_complexes'] > 0).sum()} of our genes in >=1")

# ---------------------------------------------------------------------------------------------
# 7. Ensembl BioMart: paralogue count and ortholog depth (a gene-age proxy)
#    source: useast.ensembl.org BioMart, current release, homologs page. Optional — if the fetch
#    did not land, the two columns are all-NaN and are dropped from the frozen feature set.
# ---------------------------------------------------------------------------------------------
para_p = WORK / "bm_paralog.tsv"
if para_p.exists() and para_p.stat().st_size > 1000:
    pa = pd.read_csv(para_p, sep="\t")
    col = [c for c in pa.columns if "aralog" in c][0]
    cnt = pa.dropna(subset=[col]).groupby(pa.columns[0])[col].nunique()
    df["n_paralogues"] = df["ensembl_gene_id"].map(cnt).fillna(0.0)
    df.loc[df["ensembl_gene_id"] == "", "n_paralogues"] = np.nan
    out["sources"]["ensembl_paralogues"] = dict(path="work/bm_paralog.tsv", rows=int(len(pa)))
    P(f"BioMart paralogues: {len(pa)} pairs; median count {np.nanmedian(df['n_paralogues']):.0f}")
else:
    df["n_paralogues"] = np.nan
    P("BioMart paralogues: NOT AVAILABLE")

SPECIES = ["mmusculus", "drerio", "dmelanogaster", "celegans", "scerevisiae"]
depth = np.zeros(len(df))
got = []
for sp in SPECIES:
    p = WORK / f"bm_orth_{sp}.tsv"
    if p.exists() and p.stat().st_size > 1000:
        s = set(pd.read_csv(p, sep="\t").iloc[:, 0].astype(str))
        col = df["ensembl_gene_id"].isin(s).astype(float)
        df[f"orth_{sp}"] = col
        depth += col.to_numpy()
        got.append(sp)
        P(f"BioMart orthologs {sp}: {len(s)} human genes with a homolog; ours {int(col.sum())}")
if got:
    df["ortholog_depth"] = depth
    df.loc[df["ensembl_gene_id"] == "", "ortholog_depth"] = np.nan
    out["sources"]["ensembl_orthologs"] = dict(species=got)
else:
    df["ortholog_depth"] = np.nan
    P("BioMart orthologs: NOT AVAILABLE")

# ---------------------------------------------------------------------------------------------
# 8. write
# ---------------------------------------------------------------------------------------------
df["log_n_pairs"] = np.log1p(df["n_pairs"])
df["log_collectri_regulon_size"] = np.log1p(df["collectri_regulon_size"])
df["log_dorothea_regulon_size"] = np.log1p(df["dorothea_regulon_size"])
df["log_n_paralogues"] = np.log1p(df["n_paralogues"])

FEATS = ["dm_mean_effect", "dm_sd_effect", "dm_frac_dependent", "dm_frac_strong", "dm_min_effect",
         "dm_common_essential", "replogle_essential",
         "ccle_breadth", "ccle_median_tpm", "ccle_sd",
         "collectri_regulon_size", "log_collectri_regulon_size", "collectri_indegree",
         "is_tf_collectri", "dorothea_regulon_size", "log_dorothea_regulon_size", "is_tf_dorothea",
         "humap_n_complexes", "humap_max_complex_size",
         "reactome_n_complexes", "reactome_max_complex_size",
         "n_paralogues", "log_n_paralogues", "ortholog_depth"]
cov = {f: float(df.loc[df["set"] == "held", f].notna().mean()) for f in FEATS if f in df}
out["coverage"]["features_held"] = cov
P("feature coverage on the 17,106 held genes:")
for f, v in cov.items():
    P(f"    {f:30s} {v:6.3f}")

df.to_csv(WORK / "feat281.csv", index=False)
out["n_rows"] = int(len(df))
out["log"] = LOG
(WORK / "feat281.json").write_text(json.dumps(out, indent=1, default=float))
P(f"wrote work/feat281.csv ({len(df)} rows x {df.shape[1]} cols) and work/feat281.json")

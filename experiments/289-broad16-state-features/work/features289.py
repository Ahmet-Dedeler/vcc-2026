"""289 step 2: build the frozen state features of PREREGISTRATION.md §2 for the 12 working Broad-16 lines.

Writes work/state289.csv (per line: every feature + the nuisance depth columns) and work/features289.json
(the calls, the KA-2 instrument check, the CN region definitions, agreement with the mutation annotation).
Nothing here touches a perturbation response except the basal (OR-control) log-profile, which is the recipient's
build-time state and is exactly what 268's basal PCs are built from.
"""
import json, sys, numpy as np, pandas as pd
from pathlib import Path

ROOT = Path("~/Code/virtual-cell-challenge-2026")
E268 = ROOT / "experiments/268-broad16-context-residuals/work"
W = ROOT / "experiments/289-broad16-state-features/work"
sys.path.insert(0, str(E268))
from load268 import WORKING, RESERVED, META, load_halves, build_responses, control_profiles, usable_targets, CONTROLS  # noqa

DEP = ROOT / "experiments/269-utility-router-feasibility/work/depmap"
DM = W / "depmap"
GMT = ROOT / "data/external/broad16line/metadata/hallmark_symbols.gmt"
GLOC = ROOT / "data/external/broad16line/data/gene_locations.csv"
QUAL = ROOT / "data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv"

lines = list(WORKING)
ids = {l: META.loc[l, "arxspan_id"] for l in lines}
out = {"lines": lines, "arxspan": ids}

# ---------------------------------------------------------------- hallmark sets (KA-4)
def read_gmt(path):
    S = {}
    for ln in open(path):
        p = ln.rstrip("\n").split("\t")
        S[p[0]] = [g for g in p[2:] if g]
    return S
HM = read_gmt(GMT)
P53SET = set(HM["HALLMARK_P53_PATHWAY"]); E2FSET = set(HM["HALLMARK_E2F_TARGETS"]) | set(HM["HALLMARK_G2M_CHECKPOINT"])
out["KA4_hallmark_sizes"] = {k: len(HM[k]) for k in ("HALLMARK_P53_PATHWAY", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")}
out["KA4_pass"] = all(len(HM[k]) == 200 for k in out["KA4_hallmark_sizes"])
print("KA-4 hallmark sizes", out["KA4_hallmark_sizes"], "PASS" if out["KA4_pass"] else "FAIL", flush=True)

# ---------------------------------------------------------------- DepMap gene effects (F-p53) + KA-2
hdr = pd.read_csv(DEP / "CRISPRGeneEffect.csv", nrows=0)
cols = list(hdr.columns)
def find(sym):
    c = [x for x in cols if x.split(" (")[0] == sym]
    assert len(c) == 1, (sym, c)
    return c[0]
want = {g: find(g) for g in ("TP53", "MDM2", "IER3IP1")}
ge = pd.read_csv(DEP / "CRISPRGeneEffect.csv", usecols=[cols[0]] + list(want.values()), index_col=0)
ge.columns = [c.split(" (")[0] for c in ge.columns]
# KA-2: instrument check across all 1,150 models
sub = ge[["TP53", "MDM2"]].dropna()
r = float(np.corrcoef(sub.TP53, sub.MDM2)[0, 1])
dep = sub.MDM2 <= -0.5
ka2 = {"n_models": int(len(sub)), "pearson_TP53_MDM2": r,
       "mean_TP53_if_MDM2_dep": float(sub.TP53[dep].mean()), "mean_TP53_if_not": float(sub.TP53[~dep].mean()),
       "n_MDM2_dep": int(dep.sum())}
ka2["pass"] = bool(r <= -0.3 and ka2["mean_TP53_if_MDM2_dep"] > ka2["mean_TP53_if_not"])
out["KA2"] = ka2
print("KA-2", json.dumps(ka2), flush=True)

eff = ge.loc[[ids[l] for l in lines]]
eff.index = lines
p53_wt = ((eff.MDM2 <= -0.5) & (eff.TP53 >= 0.0)).astype(int)
p53_score = (eff.TP53 - eff.MDM2)

# mutation annotation (agreement + S4 feature)
def mut_col(path):
    h = pd.read_csv(path, nrows=0); c = [x for x in h.columns if x.split(" (")[0] == "TP53"]
    assert len(c) == 1
    return pd.read_csv(path, usecols=[h.columns[0], c[0]], index_col=0).iloc[:, 0]
dam = mut_col(DM / "OmicsSomaticMutationsMatrixDamaging.csv")
hot = mut_col(DM / "OmicsSomaticMutationsMatrixHotspot.csv")
tp53_damaging = pd.Series({l: float(dam.get(ids[l], np.nan)) for l in lines})
tp53_hotspot = pd.Series({l: float(hot.get(ids[l], np.nan)) for l in lines})
p53_wt_mut = ((tp53_damaging == 0) & (tp53_hotspot == 0)).astype(int)  # no damaging and no hotspot TP53 mutation

# ---------------------------------------------------------------- F-CN region-level copy number
gl = pd.read_csv(GLOC)
gl["chr"] = gl["Chromosome/scaffold name"].astype(str)
gl["band"] = gl["Karyotype band"].astype(str)
def region(chrom, prefixes):
    m = (gl.chr == chrom) & gl.band.str.startswith(tuple(prefixes))
    assert m.sum() > 0, (chrom, prefixes)
    return float(gl.loc[m, "Gene start (bp)"].min()), float(gl.loc[m, "Gene end (bp)"].max()), int(m.sum())
REGIONS = {"cn_18q21": ("18", ["q21"]), "cn_11q23_24": ("11", ["q23", "q24"]), "cn_18q": ("18", ["q"]), "cn_11q": ("11", ["q"])}
RCOORD = {k: region(*v) for k, v in REGIONS.items()}
out["cn_regions"] = {k: {"chrom": REGIONS[k][0], "bands": REGIONS[k][1], "start": RCOORD[k][0], "end": RCOORD[k][1], "n_genes": RCOORD[k][2]} for k in REGIONS}
print("CN regions", json.dumps(out["cn_regions"]), flush=True)

prof = pd.read_csv(DM / "OmicsProfiles.csv")
def region_cn(segfile, required=True):
    seg = pd.read_csv(segfile)
    if "SegmentMean" not in seg.columns: seg = seg.rename(columns={"SegmentAbsoluteCN": "SegmentMean"})
    seg["Chromosome"] = seg.Chromosome.astype(str)
    have = set(seg.ProfileID.unique())
    pick = {}
    for l in lines:
        p = prof[(prof.ModelID == ids[l]) & (prof.ProfileID.isin(have))]
        if not len(p):
            assert not required, l
            pick[l] = None; continue
        wgs = p[p.Datatype == "wgs"]
        pick[l] = (wgs if len(wgs) else p).sort_values("ProfileID").ProfileID.iloc[0]
    vals = {}
    for l in lines:
        if pick[l] is None:
            vals[l] = {k: np.nan for k in REGIONS} | {k + "_nseg": 0 for k in REGIONS}; continue
        s = seg[seg.ProfileID == pick[l]]
        row = {}
        for k, (chrom, _) in REGIONS.items():
            st, en, _ = RCOORD[k]
            ss = s[(s.Chromosome == chrom) & (s.End >= st) & (s.Start <= en)].copy()
            ss["ov"] = np.minimum(ss.End, en) - np.maximum(ss.Start, st)
            ss = ss[ss.ov > 0]
            row[k] = float(np.average(ss.SegmentMean, weights=ss.ov)) if len(ss) else np.nan
            row[k + "_nseg"] = int(len(ss))
        vals[l] = row
    return pd.DataFrame(vals).T, pick
cn_rel, pick_rel = region_cn(DM / "OmicsCNSegmentsProfile.csv")
cn_abs, pick_abs = region_cn(DM / "OmicsAbsoluteCNSegmentsProfile.csv", required=False)
out["cn_profiles_relative"] = pick_rel; out["cn_profiles_absolute"] = pick_abs

codel = ((cn_rel.cn_18q21 < 0.85) & (cn_rel.cn_11q23_24 < 0.85)).astype(int)

# ---------------------------------------------------------------- basal hallmark scores (12-line descriptive version)
H, genes = load_halves(lines)
T82 = usable_targets(H, lines, 6)
R, ok, avail = build_responses(H, genes, lines, T82, keys=("full",), min_cpm=5.0)
C = control_profiles(H, genes, lines, ok)  # (12, G)
gk = genes[ok]
mP = np.isin(gk, list(P53SET)); mE = np.isin(gk, list(E2FSET))
e2f_g2m = pd.Series({l: float(C[j, mE].mean() - C[j].mean()) for j, l in enumerate(lines)})
p53_path = pd.Series({l: float(C[j, mP].mean() - C[j].mean()) for j, l in enumerate(lines)})
out["gene_axis_all12"] = {"G": int(ok.sum()), "n_p53_pathway": int(mP.sum()), "n_e2f_g2m": int(mE.sum()), "n_targets": len(T82)}
print("gene axis", out["gene_axis_all12"], flush=True)

# ---------------------------------------------------------------- depth (nuisance) + lineage
q = pd.read_csv(QUAL).set_index("cell_line")
inv = json.load(open(E268 / "inventory.json"))

df = pd.DataFrame({
    "line": lines,
    "arxspan": [ids[l] for l in lines],
    "batch": [META.loc[l, "batch"] for l in lines],
    "lineage": [META.loc[l, "OncotreeLineage"] for l in lines],
    "hpv": [1 if l in ("SKGII", "C4I") else 0 for l in lines],
    "TP53_effect": [float(eff.TP53[l]) for l in lines],
    "MDM2_effect": [float(eff.MDM2[l]) for l in lines],
    "IER3IP1_effect": [float(eff.IER3IP1[l]) for l in lines],
    "p53_wt": [int(p53_wt[l]) for l in lines],
    "p53_score": [float(p53_score[l]) for l in lines],
    "TP53_damaging": [tp53_damaging[l] for l in lines],
    "TP53_hotspot": [tp53_hotspot[l] for l in lines],
    "p53_wt_mut": [int(p53_wt_mut[l]) for l in lines],
    "e2f_g2m": [float(e2f_g2m[l]) for l in lines],
    "p53_path": [float(p53_path[l]) for l in lines],
    "cn_18q21": [float(np.log2(cn_rel.cn_18q21[l])) for l in lines],
    "cn_11q23_24": [float(np.log2(cn_rel.cn_11q23_24[l])) for l in lines],
    "cn_18q": [float(np.log2(cn_rel.cn_18q[l])) for l in lines],
    "cn_11q": [float(np.log2(cn_rel.cn_11q[l])) for l in lines],
    "cn_18q21_rel": [float(cn_rel.cn_18q21[l]) for l in lines],
    "cn_11q23_24_rel": [float(cn_rel.cn_11q23_24[l]) for l in lines],
    "cn_18q21_abs": [float(cn_abs.cn_18q21[l]) for l in lines],
    "cn_11q23_24_abs": [float(cn_abs.cn_11q23_24[l]) for l in lines],
    "codel": [int(codel[l]) for l in lines],
    "avg_umis_per_ko": [float(q.avg_umis_per_ko[l]) for l in lines],
    "total_umis_per_ko": [float(q.total_umis_per_ko[l]) for l in lines],
    "n_cells_per_ko": [float(q.n_cells_per_ko[l]) for l in lines],
}).set_index("line")
df.to_csv(W / "state289.csv")

# degeneracy contingency (prereg §2 F-p53)
out["p53_wt_counts"] = {"wt": int(df.p53_wt.sum()), "def": int((1 - df.p53_wt).sum())}
out["p53_binary_degenerate"] = bool(min(out["p53_wt_counts"].values()) < 2)
out["p53_mut_counts"] = {"wt": int(df.p53_wt_mut.sum()), "def": int((1 - df.p53_wt_mut).sum())}
out["p53_dependency_vs_mutation_agreement"] = float((df.p53_wt == df.p53_wt_mut).mean())
out["codel_count"] = int(df.codel.sum())
json.dump(out, open(W / "features289.json", "w"), indent=1, default=float)
print(df.round(3).to_string(), flush=True)
print("p53 binary counts", out["p53_wt_counts"], "degenerate", out["p53_binary_degenerate"],
      "| mutation counts", out["p53_mut_counts"], "agreement", out["p53_dependency_vs_mutation_agreement"],
      "| codel lines", out["codel_count"], flush=True)

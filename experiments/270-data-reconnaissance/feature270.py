"""270 — compute the FROZEN sign-stability feature (work/FEATURE-DEFINITION.md, sha f8a6bd39…) from GPSAdb 2.0.

Inputs : work/gpsadb_metadata_detail.csv (per-group metadata, downloaded from gpsadb.com Download tab),
         the figshare archive GPSAdiffResults.zip read member-by-member over HTTP range requests (only the
         groups whose perturbed gene is on the 300-target panel; compressed-byte cap CAP).
Outputs: work/gpsadb_members/<Dxxxxx>_diffresults.csv (raw tables kept), work/feature270.csv (per target),
         work/feature270_pairs.csv (per pair), work/feature270.json (coverage counts, bytes fetched).
"""
import io, json, os, re, sys, time, zipfile, numpy as np, pandas as pd, urllib.request
from types import SimpleNamespace
W = "~/Code/virtual-cell-challenge-2026/experiments/270-data-reconnaissance/work"
PANEL = set(pd.read_csv("~/Code/virtual-cell-challenge-2026/data/challenge/pert_counts.csv")["target_gene"])
URL = "https://ndownloader.figshare.com/files/55769468"
CAP = int(os.environ.get("CAP270", 1_500_000_000))
N_TOP = {"primary": 100, "n50": 50, "n200": 200}
MIN_UNION = 20

sys.path.insert(0, W)
src = open(f"{W}/remotezip270.py").read().split("url, out = sys.argv[1], sys.argv[2]")[0]
ns = {}; exec(src, ns); RangeFile = ns["RangeFile"]

meta_path = sys.argv[1] if len(sys.argv) > 1 else f"{W}/gpsadb_metadata_detail.csv"
meta = pd.read_csv(meta_path, low_memory=False)
meta.columns = [c.strip() for c in meta.columns]
print("metadata columns:", list(meta.columns)[:40], "rows", len(meta))
# --- column resolution (recorded in the json; no other freedom)
def pick(cands):
    for c in cands:
        for col in meta.columns:
            if col.lower() == c.lower(): return col
    for c in cands:
        for col in meta.columns:
            if c.lower() in col.lower(): return col
    return None
if {"Dataset", "gene", "cell_line", "GSE", "method"} <= set(meta.columns):  # work/gpsadb_metadata_merged.csv (explicit names)
    COL = {"id": "Dataset", "gene": "gene", "line": "cell_line", "study": "GSE", "method": "method"}
else:
    COL = {"id": pick(["Dataset", "dataset_id", "GPSA_ID", "ID"]), "gene": pick(["perturbation_gene", "perturbed_gene", "gene_symbol", "Gene", "gene"]),
           "line": pick(["cell_line", "cellline", "cell"]), "study": pick(["GSE", "gse", "study", "series"]),
           "method": pick(["perturbation_method", "method", "perturb_type", "type"])}
print("column map:", COL)
assert COL["id"] and COL["gene"], "cannot resolve id/gene columns"
meta["_id"] = meta[COL["id"]].astype(str).str.extract(r"(D\d+)", expand=False).fillna(meta[COL["id"]].astype(str))
meta["_gene"] = meta[COL["gene"]].astype(str).str.strip()
meta["_line"] = meta[COL["line"]].astype(str).str.strip() if COL["line"] else "NA"
meta["_study"] = meta[COL["study"]].astype(str).str.extract(r"(GSE\d+)", expand=False) if COL["study"] else None
if COL["study"]: meta["_study"] = meta["_study"].fillna(meta[COL["study"]].astype(str))
meta["_method"] = meta[COL["method"]].astype(str).str.strip() if COL["method"] else "NA"

# --- coverage counts (metadata only)
cov = meta[meta["_gene"].isin(PANEL)].copy()
per_t = cov.groupby("_gene").agg(n_groups=("_id", "nunique"), n_lines=("_line", "nunique"), n_studies=("_study", "nunique"),
                                methods=("_method", lambda s: ";".join(sorted(set(s)))))
summary = {"panel_targets_covered": int(per_t.shape[0]), "targets_ge2_lines": int((per_t.n_lines >= 2).sum()),
           "targets_ge2_studies": int((per_t.n_studies >= 2).sum()), "targets_ge2_groups": int((per_t.n_groups >= 2).sum()),
           "groups_on_panel": int(cov["_id"].nunique()), "method_mix_groups_on_panel": cov["_method"].value_counts().to_dict(),
           "method_mix_all_groups": meta["_method"].value_counts().head(12).to_dict(), "column_map": COL,
           "all_groups": int(meta["_id"].nunique()), "all_genes": int(meta["_gene"].nunique()), "all_lines": int(meta["_line"].nunique()),
           "all_studies": int(meta["_study"].nunique()) if COL["study"] else None}
print(json.dumps(summary, indent=1)[:2500])
per_t.to_csv(f"{W}/gpsadb_panel_coverage.csv")

# --- fetch the panel groups' DE tables (only those with >=2 groups; singletons cannot form a pair)
need = cov[cov["_gene"].isin(per_t.index[per_t.n_groups >= 2])]["_id"].unique().tolist()
listing = pd.read_csv(f"{W}/gpsadb_zip_listing.tsv", sep="\t"); sizes = dict(zip(listing.name, listing.compress_size))
members = [f"{i}_diffresults.csv" for i in need if f"{i}_diffresults.csv" in sizes]
print(f"groups needed {len(need)}, members present in the zip {len(members)}, compressed bytes {sum(sizes[m] for m in members):,}")
os.makedirs(f"{W}/gpsadb_members", exist_ok=True)
rf = RangeFile(URL); zf = zipfile.ZipFile(rf); got = 0; tot = 0; t0 = time.time()
for m in members:
    p = f"{W}/gpsadb_members/{m}"
    if os.path.exists(p): got += 1; continue
    if tot + sizes[m] > CAP: print("CAP reached", tot); break
    open(p, "wb").write(zf.read(m)); tot += sizes[m]; got += 1
    if got % 25 == 0: print(f"  {got}/{len(members)} members, {rf.fetched/1e6:.0f} MB fetched, {time.time()-t0:.0f} s", flush=True)
print(f"members on disk {got}/{len(members)}; bytes fetched this run {rf.fetched:,}")
summary["members_downloaded"] = got; summary["bytes_fetched"] = int(rf.fetched); summary["compressed_bytes_selected"] = int(tot)

# --- per-group top sets (frozen: rank by |Wald z| = |logFC/lfcSE|; exclude the target; drop lfc==0/NA)
def load_group(gid, gene):
    p = f"{W}/gpsadb_members/{gid}_diffresults.csv"
    if not os.path.exists(p) or os.path.getsize(p) < 1000: return None
    try: d = pd.read_csv(p, usecols=["logFC", "lfcSE", "gene"]).dropna(subset=["logFC", "gene"])
    except Exception as e: print("unreadable", p, e); return None
    d = d[(d.gene != gene) & (d.logFC != 0)]
    d = d.groupby("gene", as_index=False).agg(logFC=("logFC", "mean"), lfcSE=("lfcSE", "mean"))  # duplicate symbols averaged
    d["z"] = d.logFC / d.lfcSE.replace(0, np.nan)
    d["z"] = d["z"].fillna(0.0)
    return d.set_index("gene")
pairs = []; feats = []
for gene, grp in cov.groupby("_gene"):
    ids = grp.drop_duplicates("_id")
    recs = [SimpleNamespace(**r) for r in ids[["_id", "_gene", "_line", "_study", "_method"]].to_dict("records")]
    tabs = {r._id: load_group(r._id, gene) for r in recs}
    tabs = {k: v for k, v in tabs.items() if v is not None and len(v) >= 200}
    if len(tabs) < 2: continue
    info = {r._id: r for r in recs}
    tops = {}
    for k, d in tabs.items():
        tops[k] = {name: set(d["z"].abs().sort_values(ascending=False).index[:n]) for name, n in N_TOP.items()}
        tops[k]["lfc100"] = set(d["logFC"].abs().sort_values(ascending=False).index[:100])
    keys = list(tabs); rec = {"target": gene, "n_groups": len(keys), "n_studies": len({info[k]._study for k in keys}),
                              "n_lines": len({info[k]._line for k in keys}), "methods": ";".join(sorted({info[k]._method for k in keys}))}
    acc = {v: [] for v in ["SS_cross", "SS_cross_n50", "SS_cross_n200", "SS_cross_lfc100", "RHO_cross", "SS_line", "SS_within_study"]}
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]; da, db = tabs[a], tabs[b]
            common = da.index.intersection(db.index)
            cross_study = info[a]._study != info[b]._study; cross_line = info[a]._line != info[b]._line
            def agree(setname):
                U = [g for g in (tops[a][setname] | tops[b][setname]) if g in common]
                if len(U) < MIN_UNION: return np.nan, len(U)
                sa = np.sign(da.loc[U, "logFC"].values); sb = np.sign(db.loc[U, "logFC"].values)
                return float(np.mean((sa == sb) & (sa != 0))), len(U)
            c100, n100 = agree("primary"); c50, _ = agree("n50"); c200, _ = agree("n200"); cl, _ = agree("lfc100")
            U = [g for g in (tops[a]["primary"] | tops[b]["primary"]) if g in common]
            rho = float(pd.Series(da.loc[U, "logFC"].values).corr(pd.Series(db.loc[U, "logFC"].values), method="spearman")) if len(U) >= MIN_UNION else np.nan
            pairs.append(dict(target=gene, g1=a, g2=b, study1=info[a]._study, study2=info[b]._study, line1=info[a]._line, line2=info[b]._line,
                              cross_study=cross_study, cross_line=cross_line, n_union=n100, c100=c100, c50=c50, c200=c200, c_lfc100=cl, rho=rho))
            if cross_study:
                acc["SS_cross"].append(c100); acc["SS_cross_n50"].append(c50); acc["SS_cross_n200"].append(c200)
                acc["SS_cross_lfc100"].append(cl); acc["RHO_cross"].append(rho)
            else:
                acc["SS_within_study"].append(c100)
            if cross_line: acc["SS_line"].append(c100)
    for v, L in acc.items():
        L = [x for x in L if not (isinstance(x, float) and np.isnan(x))]
        rec[v] = float(np.mean(L)) if L else np.nan; rec["n_pairs_" + v] = len(L)
    feats.append(rec)
F = pd.DataFrame(feats); P = pd.DataFrame(pairs)
F.to_csv(f"{W}/feature270.csv", index=False); P.to_csv(f"{W}/feature270_pairs.csv", index=False)
summary["targets_with_SS_cross"] = int(F["SS_cross"].notna().sum()) if len(F) else 0
summary["targets_with_any_pair"] = int(len(F))
if len(F):
    summary["SS_cross_describe"] = F["SS_cross"].describe().to_dict(); summary["SS_within_study_describe"] = F["SS_within_study"].describe().to_dict()
json.dump(summary, open(f"{W}/feature270.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in summary.items() if k.startswith(("targets_", "members", "bytes", "SS_"))}, indent=1, default=str))

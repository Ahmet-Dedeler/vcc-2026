"""270 — extract the panel targets' CMap signatures (sh.cgs and xpr) from the level-5 gctx files and compute the
declared SECONDARY sign-stability features (work/FEATURE-DEFINITION-CMAP.md). Read-only on data/; writes work/ only.
Memory: reads the xpr matrix in chunk-aligned row blocks (262 rows × 12,328 × 4 B ≈ 13 MB each)."""
import h5py, json, time, numpy as np, pandas as pd, pyarrow.parquet as pq
W = "~/Code/virtual-cell-challenge-2026/experiments/270-data-reconnaissance/work"
C = "~/Code/virtual-cell-challenge-2026/data/external/cmap"
PANEL = set(pd.read_csv("~/Code/virtual-cell-challenge-2026/data/challenge/pert_counts.csv")["target_gene"])
N_TOP = 100; MIN_UNION = 20
t0 = time.time()
sig = pq.read_table(f"{C}/genetic_sigs.parquet").to_pandas()
sig["pert_type"] = sig.pert_type.astype(str); sig["cmap_name"] = sig.cmap_name.astype(str); sig["cell_iname"] = sig.cell_iname.astype(str)
sig = sig[sig.cmap_name.isin(PANEL)]
gi = pd.read_csv(f"{C}/geneinfo_beta.txt", sep="\t", dtype=str)
lm_ids = set(gi.loc[gi.feature_space.str.contains("landmark", case=False, na=False), "gene_id"].astype(str))
sym = dict(zip(gi.gene_id.astype(str), gi.gene_symbol.astype(str)))

def pull(gctx, wanted_ids):
    with h5py.File(gctx, "r") as h:
        col = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["0/META/COL/id"][()]])
        row = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in h["0/META/ROW/id"][()]])
        pos = {s: i for i, s in enumerate(col)}
        idx = np.array(sorted(pos[s] for s in wanted_ids if s in pos))
        d = h["0/DATA/0/matrix"]; out = np.empty((len(idx), d.shape[1]), np.float32)
        if d.chunks is None:  # contiguous layout: read the needed rows one by one (49 KB each)
            for k, r in enumerate(idx): out[k] = d[r]
            return col[idx], row, out
        ch = d.chunks[0]
        # chunked layout: group needed rows by chunk-aligned block (262 rows x 12,328 cols = 13 MB per block)
        blocks = {}
        for k, r in enumerate(idx): blocks.setdefault(r // ch, []).append((k, r))
        for b, items in sorted(blocks.items()):
            lo, hi = b * ch, min((b + 1) * ch, d.shape[0]); blk = d[lo:hi]
            for k, r in items: out[k] = blk[r - lo]
        return col[idx], row, out

feat_rows = []; summary = {}
sets = {}
for pt, path in [("trt_sh.cgs", f"{C}/level5_beta_trt_sh.gctx"), ("trt_xpr", f"{C}/level5_beta_trt_xpr.gctx")]:
    s = sig[sig.pert_type == pt]
    ids, genes, Z = pull(path, s.sig_id.tolist())
    print(f"{pt}: {len(ids)} signatures pulled in {time.time()-t0:.0f} s; genes {len(genes)}", flush=True)
    meta = s.set_index("sig_id").loc[ids]
    # per (gene, line) consensus: sh.cgs is one per pair already; xpr -> mean over its signatures
    key = pd.MultiIndex.from_arrays([meta.cmap_name.values, meta.cell_iname.values])
    df = pd.DataFrame(Z, index=key); per = df.groupby(level=[0, 1]).mean()
    sets[pt] = (per, genes)
    summary[pt] = {"signatures": int(len(ids)), "gene_line_pairs": int(len(per)), "targets": int(per.index.get_level_values(0).nunique()),
                   "lines": int(per.index.get_level_values(1).nunique())}
    np.savez_compressed(f"{W}/cmap270_{pt.replace('.', '_')}.npz", Z=per.values.astype(np.float32), target=per.index.get_level_values(0).values.astype(str),
                        line=per.index.get_level_values(1).values.astype(str), genes=np.array(genes).astype(str))
genes = list(sets["trt_sh.cgs"][1]); gx = list(sets["trt_xpr"][1])
assert set(genes) == set(gx), "gene sets differ between the two gctx files"
if genes != gx:  # align the xpr matrix to the sh gene order by gene id
    per_x_df = sets["trt_xpr"][0]; per_x_df.columns = gx; per_x_df = per_x_df[genes]; per_x_df.columns = range(len(genes))
    sets["trt_xpr"] = (per_x_df, genes); print("xpr columns re-ordered to the sh gene order")
lm_mask = np.array([g in lm_ids for g in genes]); gsym = np.array([sym.get(g, g) for g in genes])
print("landmark genes on axis:", int(lm_mask.sum()))

def topset(z, tgt, mask, n):
    zz = z.copy(); zz[~mask] = 0; zz[gsym == tgt] = 0
    return set(np.argsort(-np.abs(zz))[:n])
def agree(za, zb, U):
    U = np.array(sorted(U)); sa, sb = np.sign(za[U]), np.sign(zb[U])
    return float(np.mean((sa == sb) & (sa != 0))), len(U)
def pair_stat(za, zb, tgt, mask, n):
    U = topset(za, tgt, mask, n) | topset(zb, tgt, mask, n)
    if len(U) < MIN_UNION: return np.nan
    return agree(za, zb, U)[0]

per_sh, per_x = sets["trt_sh.cgs"][0], sets["trt_xpr"][0]
allmask = np.ones(len(genes), bool)
targets = sorted(set(per_sh.index.get_level_values(0)) | set(per_x.index.get_level_values(0)))
pairs = []
for tgt in targets:
    rec = {"target": tgt}
    for name, per in [("sh", per_sh), ("xpr", per_x)]:
        if tgt not in per.index.get_level_values(0): rec[f"n_lines_{name}"] = 0; continue
        sub = per.loc[tgt]; lines = list(sub.index); rec[f"n_lines_{name}"] = len(lines)
        vals = {"lm100": [], "all100": [], "lm50": []}
        for i in range(len(lines)):
            for j in range(i + 1, len(lines)):
                za, zb = sub.iloc[i].values, sub.iloc[j].values
                c1 = pair_stat(za, zb, tgt, lm_mask, 100); c2 = pair_stat(za, zb, tgt, allmask, 100); c3 = pair_stat(za, zb, tgt, lm_mask, 50)
                vals["lm100"].append(c1); vals["all100"].append(c2); vals["lm50"].append(c3)
                pairs.append(dict(target=tgt, kind=name, line1=lines[i], line2=lines[j], c_lm100=c1, c_all100=c2, c_lm50=c3))
        for k, L in vals.items():
            L = [x for x in L if not np.isnan(x)]
            rec[f"SS_cmap_{name}" + ("" if k == "lm100" else "_" + k)] = float(np.mean(L)) if L else np.nan
        rec[f"n_pairs_{name}"] = len([x for x in vals["lm100"] if not np.isnan(x)])
    # cross-method within line
    if tgt in per_sh.index.get_level_values(0) and tgt in per_x.index.get_level_values(0):
        a, b = per_sh.loc[tgt], per_x.loc[tgt]; common = [l for l in a.index if l in b.index]; L = []
        for l in common:
            c = pair_stat(a.loc[l].values, b.loc[l].values, tgt, lm_mask, 100); L.append(c)
            pairs.append(dict(target=tgt, kind="method", line1=l, line2=l, c_lm100=c, c_all100=np.nan, c_lm50=np.nan))
        L = [x for x in L if not np.isnan(x)]
        rec["SS_cmap_method"] = float(np.mean(L)) if L else np.nan; rec["n_lines_method"] = len(L)
    feat_rows.append(rec)
F = pd.DataFrame(feat_rows); F.to_csv(f"{W}/feature270_cmap.csv", index=False); pd.DataFrame(pairs).to_csv(f"{W}/feature270_cmap_pairs.csv", index=False)
summary["targets_total"] = int(len(F)); summary["targets_SS_sh"] = int(F["SS_cmap_sh"].notna().sum()) if "SS_cmap_sh" in F else 0
summary["targets_SS_xpr"] = int(F["SS_cmap_xpr"].notna().sum()) if "SS_cmap_xpr" in F else 0
summary["targets_SS_method"] = int(F["SS_cmap_method"].notna().sum()) if "SS_cmap_method" in F else 0
for c in [c for c in F.columns if c.startswith("SS_")]: summary[c + "_describe"] = F[c].describe().round(4).to_dict()
summary["seconds"] = round(time.time() - t0, 1)
json.dump(summary, open(f"{W}/feature270_cmap.json", "w"), indent=1, default=str); print(json.dumps(summary, indent=1, default=str))

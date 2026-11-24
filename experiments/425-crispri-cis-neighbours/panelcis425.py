"""425: the A/B/C panel's CRISPRi cis-neighbours, and whether 417a's bulk and carriers already predict them.

For each of the 300 panel targets: expressed non-target genes whose TSS lies within 5 kb of the target's TSS (GENCODE v44
gene-level TSS). Reports, per neighbour: distance, strand relation, control CPM in A/B/C, K562-GWPS donor delta (if the
target and gene are on the donor axis), 417a's dumped bulk and carrier values and their within-row rank.
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
tss = pd.read_csv(HERE / "work/tss_v44.tsv", sep="\t").set_index("gene_name")
z = np.load(ROOT / "experiments/417-bulk-sharpen/work/dump_ka417a_g065.npz", allow_pickle=True)
bulk, carr, perts, genes = z["bulk"], z["carriers"], z["perts"].astype(str), z["genes"].astype(str)
gpos = {g: i for i, g in enumerate(genes)}
panel = set(perts)


def ctrl_cpm(ctx):
    with h5py.File(ROOT / f"data/challenge/context_{ctx}.h5ad", "r") as h:
        ind = np.asarray(h["X"]["indices"][:], dtype=np.int64)
        dat = np.asarray(h["X"]["data"][:], dtype=np.float64)
        ip = np.asarray(h["X"]["indptr"][:])
        lib = np.add.reduceat(dat, ip[:-1])
        rowid = np.repeat(np.arange(ip.size - 1), np.diff(ip))
        cpm = np.bincount(ind, weights=dat / lib[rowid] * 1e6, minlength=genes.size) / (ip.size - 1)
    return cpm


cpm = {c: ctrl_cpm(c) for c in "ABC"}
k562 = np.load(ROOT / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
kg = {g: i for i, g in enumerate((ROOT / "data/donors/k562_gwps/genes.txt").read_text().split())}
kp = {g: i for i, g in enumerate((ROOT / "data/donors/k562_gwps/perts.txt").read_text().split())}

chrom_groups = {c: d for c, d in tss.groupby("chrom")}
rows = []
for pi, p in enumerate(perts):
    if p not in tss.index:
        continue
    t = tss.loc[p]
    d = chrom_groups[t["chrom"]]
    near = d[(np.abs(d["tss"] - t["tss"]) <= 5000) & (d.index != p)]
    for g, r in near.iterrows():
        if g not in gpos:
            continue
        j = gpos[g]
        row_b, row_c = bulk[pi], carr[pi]
        rank_b = int((row_b < row_b[j]).sum())          # 0 = most negative in the row
        kd = np.nan
        if p in kp and g in kg:
            kd = float(k562[kp[p], kg[g]])
        rows.append(dict(target=p, gene=g, dist=int(abs(r["tss"] - t["tss"])),
                         orient="same" if r["strand"] == t["strand"] else "opposite",
                         in_panel=g in panel, cpmA=cpm["A"][j], cpmB=cpm["B"][j], cpmC=cpm["C"][j],
                         k562_delta=kd, bulk=row_b[j], bulk_rank=rank_b, bulk_row_min=row_b.min(),
                         carrier=row_c[j]))
out = pd.DataFrame(rows).sort_values("dist")
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 400)
print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
out.to_csv(HERE / "work/panelcis425.csv", index=False)
expressed = out[(out[["cpmA", "cpmB", "cpmC"]].max(axis=1) > 5) & ~out["in_panel"]]
print(f"\npanel targets with >=1 expressed non-panel neighbour within 5 kb: {expressed['target'].nunique()} "
      f"({len(expressed)} pairs); within 1 kb: {expressed[expressed['dist'] <= 1000]['target'].nunique()}")

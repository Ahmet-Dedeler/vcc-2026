import sys
from pathlib import Path
import anndata as ad
sys.path.insert(0,"~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import DONOR_PATH, ROOT
chal=set(l.strip() for l in open(ROOT/"data/challenge/pert_counts.csv"))-{"target_gene",""}
axis=set(g.strip() for g in open(ROOT/"data/challenge/gene_names.csv") if g.strip())
pan={}
for hn,hp in [("h1-test","data/harness/h1-test.h5ad"),("orion","data/harness/orion-hct116-thinned.h5ad")]:
    a=ad.read_h5ad(ROOT/hp,backed="r"); pan[hn]=set(a.obs["target_gene"].astype(str))-{"non-targeting"}; del a
print(f"{'donor':16s}{'perts':>7}{'genes':>7}{'chal300':>9}{'axis':>7}{'h1_85':>7}{'ori_132':>9}")
for n,p in DONOR_PATH.items():
    if not p.exists(): print(f"{n:16s}  MISSING {p}"); continue
    pt=set(l.strip() for l in open(p/"perts.txt") if l.strip())
    gn=set(l.strip() for l in open(p/"genes.txt") if l.strip())
    print(f"{n:16s}{len(pt):>7}{len(gn):>7}{len(pt&chal):>9}{len(gn&axis):>7}"
          f"{len(pt&pan['h1-test']):>7}{len(pt&pan['orion']):>9}")

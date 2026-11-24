import numpy as np, pandas as pd, h5py, re
W="~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
C="~/Code/virtual-cell-challenge-2026/data/challenge/"
genes=[l.strip() for l in open(C+"gene_names.csv") if l.strip()][1:]
G=np.array(genes); gidx={g:i for i,g in enumerate(genes)}
targets=[l.strip() for l in open(C+"pert_counts.csv") if l.strip()][1:]
print("challenge genes:",len(genes),"targets:",len(targets))
def load_wide(path):
    df=pd.read_csv(path, index_col=0)
    syms=[re.sub(r"\s*\(\d+\)$","",c) for c in df.columns]
    df.columns=syms
    df=df.loc[:,~df.columns.duplicated(keep="first")]
    return df
expr=load_wide(W+"Expression_subset.csv"); print("expr:",expr.shape)
dep=load_wide(W+"CRISPRGeneEffect_subset.csv"); print("dep:",dep.shape)
inter_e=[g for g in genes if g in expr.columns]
inter_d=[g for g in genes if g in dep.columns]
print("expression coverage: %d/%d (loss %d)"%(len(inter_e),len(genes),len(genes)-len(inter_e)))
print("dependency coverage: %d/%d (loss %d)"%(len(inter_d),len(genes),len(genes)-len(inter_d)))
print("missing expr sample:",[g for g in genes if g not in expr.columns][:10])
E=np.full((len(expr),len(genes)),np.nan); D=np.full((len(dep),len(genes)),np.nan)
for j,g in enumerate(inter_e): E[:,gidx[g]]=expr[g].values
for j,g in enumerate(inter_d): D[:,gidx[g]]=dep[g].values
model_ids=list(expr.index)
dep_ids=list(dep.index)
np.savez_compressed(W+"three_lines.npz", genes=G, targets=np.array(targets),
  expr_ids=np.array(model_ids), expr_logtpm1p=E, dep_ids=np.array(dep_ids), dep_chronos=D,
  expr_missing=np.array([g for g in genes if g not in expr.columns]),
  dep_missing=np.array([g for g in genes if g not in dep.columns]))
print("wrote three_lines.npz")
# challenge control pseudobulk (total-count CPM, chunked backed read)
import anndata
cpm={}
for ctx in "ABC":
    a=anndata.read_h5ad(C+f"context_{ctx}.h5ad", backed="r")
    tot=np.zeros(a.shape[1]); L=0.0; n=a.shape[0]
    for s in range(0,n,2000):
        X=a.X[s:s+2000]; X=np.asarray(X.todense() if not isinstance(X,np.ndarray) else X)
        tot+=X.sum(axis=0); L+=X.sum()
    cpm[ctx]=tot/L*1e6
    print(ctx,"libtotal %.3g"%L,"medianCPMtop",np.sort(cpm[ctx])[-3:])
np.savez_compressed(W+"challenge_cpm.npz", **{"cpm_"+k:v for k,v in cpm.items()}, genes=G)

import numpy as np, pandas as pd
W="~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
B=np.load(W+"three_lines.npz", allow_pickle=True)
genes=list(B["genes"]); gidx={g:i for i,g in enumerate(genes)}
targets=list(B["targets"])
expr_ids=list(B["expr_ids"]); E=B["expr_logtpm1p"]
dep_ids=list(B["dep_ids"]); D=B["dep_chronos"]
ID={"A":"ACH-000995","B":"ACH-001086","C":"ACH-000518"}
NAMES={"ACH-000995":"JURKAT","ACH-001086":"HELA","ACH-000518":"CAL33","ACH-000551":"K562",
 "ACH-000971":"HCT116","ACH-000739":"HEPG2","ACH-000997":"HCT15","ACH-000823":"KYSE140",
 "ACH-000865":"KYSE450","ACH-000348":"RPMI7951","ACH-001652":"SKGII","ACH-000036":"U343",
 "ACH-001334":"C4I","ACH-001992":"ONE58","ACH-002042":"T3M5","ACH-000250":"KMRC20",
 "ACH-000460":"SNU8","ACH-001852":"HKGZCC"}
# (c) panel dependency table
donors=[m for m in dep_ids if m not in ID.values()]
rows=[]
for t in targets:
    j=gidx.get(t)
    r={"target":t}
    for ctx,mid in ID.items():
        r[f"chronos_{ctx}"]=float(D[dep_ids.index(mid),j]) if (mid in dep_ids and j is not None and np.isfinite(D[dep_ids.index(mid),j])) else np.nan
    dv=np.array([D[dep_ids.index(m),j] for m in donors if j is not None])
    dv=dv[np.isfinite(dv)]
    r["donor_mean"]=float(dv.mean()) if len(dv) else np.nan
    r["donor_sd"]=float(dv.std()) if len(dv) else np.nan
    r["n_donors_dep"]=int((dv<-0.5).sum())
    rows.append(r)
pt=pd.DataFrame(rows)
pt.to_csv(W+"panel_dependency.csv",index=False)
print("panel rows:",len(pt))
for ctx in "ABC":
    v=pt[f"chronos_{ctx}"].dropna()
    print(f"context {ctx}: n={len(v)} strong(<-1):{(v<-1).sum()} dep(<-0.5):{(v<-0.5).sum()} median:{v.median():.3f}")
    print("  most dependent:", pt.loc[v.nsmallest(5).index,["target",f"chronos_{ctx}"]].values.tolist())
# context-essential in line but not donors (divergent)
for ctx in "AC":
    v=pt[f"chronos_{ctx}"]
    div=pt[(v<-0.5)&(pt["donor_mean"]>-0.3)]
    print(ctx,"dependent-here/avg-not-in-donors:",div["target"].tolist())
# cross-check 24Q4 expr vs held 25Q3 CCLE cache on the 3 lines
C=np.load("~/Code/virtual-cell-challenge-2026/data/external/ccle_broad_2025/ccle_cache.npz")
smp=[str(x) for x in C["samples"]]; cg=[str(x) for x in C["genes"]]
ci={s:i for i,s in enumerate(smp)}
cm={"ACH-000995":"JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","ACH-001086":"HELA_CERVIX","ACH-000518":"CAL33_UPPER_AERODIGESTIVE_TRACT"}
X=C["X"]
for mid,lab in [("ACH-000995","A"),("ACH-001086","B"),("ACH-000518","C")]:
    i=expr_ids.index(mid); e=E[i]
    j=ci[cm[mid]]; t=np.asarray(X[j]).ravel()
    common=[g for g in genes if g in cg]
    ii=[gidx[g] for g in common]; jj=[cg.index(g) for g in common]
    a=np.log1p(np.maximum(e[ii],0)); b=np.log1p(np.maximum(t[jj],0))
    m=np.isfinite(a)&np.isfinite(b)
    print(f"24Q4-vs-25Q3mirror {lab}: pearson {np.corrcoef(a[m],b[m])[0,1]:.4f} n={m.sum()}")
# mutations summary
mu=pd.read_csv(W+"Mutations_3lines.csv", usecols=["HugoSymbol","ModelID","VariantType","VepImpact","ProteinChange"])
print("mut rows per line:"); print(mu["ModelID"].value_counts().to_string())
for mid,lab in [("ACH-000995","A/Jurkat"),("ACH-001086","B/HeLa"),("ACH-000518","C/CAL33")]:
    s=mu[(mu.ModelID==mid)&(mu.HugoSymbol.isin(targets))]
    print(lab,"panel-target mutations:",len(s))
    if len(s): print(s[["HugoSymbol","VariantType","VepImpact","ProteinChange"]].to_string())

import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
W="~/Code/virtual-cell-challenge-2026/experiments/360-depmap-three-lines/work/"
B=np.load(W+"three_lines.npz"); Pm=np.load(W+"challenge_cpm.npz")
genes=list(B["genes"]); gidx={g:i for i,g in enumerate(genes)}
expr_ids=list(B["expr_ids"]); E=B["expr_logtpm1p"]
dep_ids=list(B["dep_ids"]); D=B["dep_chronos"]
ID={"A":"ACH-000995","B":"ACH-001086","C":"ACH-000518"}
NAMES={"ACH-000995":"JURKAT","ACH-001086":"HELA","ACH-000518":"CAL33","ACH-000551":"K562",
 "ACH-000971":"HCT116","ACH-000739":"HEPG2","ACH-000997":"HCT15","ACH-000823":"KYSE140",
 "ACH-000865":"KYSE450","ACH-000348":"RPMI7951","ACH-001652":"SKGII","ACH-000036":"U343",
 "ACH-001334":"C4I","ACH-001992":"ONE58","ACH-002042":"T3M5","ACH-000250":"KMRC20",
 "ACH-000460":"SNU8","ACH-001852":"HKGZCC"}
print("=== (a) identity validation: DepMap log(TPM+1) vs challenge log(CPM+1) ===")
for ctx in "ABC":
    c=np.log1p(Pm["cpm_"+ctx])
    print(f"-- context {ctx} (claimed {NAMES[ID[ctx]]})")
    for mid in [ID["A"],ID["B"],ID["C"]]:
        i=expr_ids.index(mid); e=E[i]
        m=np.isfinite(e)&np.isfinite(c)
        print(f"   vs {NAMES[mid]:8s} pearson {pearsonr(e[m],c[m])[0]:.4f}  spearman {spearmanr(e[m],c[m])[0]:.4f}  n={m.sum()}")
print("=== (b) nearest donors by expression (pearson, common finite genes) ===")
donors=[m for m in expr_ids if m not in ID.values()]
for ctx in "ABC":
    c=np.log1p(Pm["cpm_"+ctx])
    rs=[]
    for mid in expr_ids:
        e=E[expr_ids.index(mid)]; m=np.isfinite(e)&np.isfinite(c)
        rs.append((NAMES.get(mid,mid),pearsonr(e[m],c[m])[0]))
    rs.sort(key=lambda t:-t[1])
    print(ctx,["%s %.4f"%t for t in rs])
print("=== (b2) nearest donors by dependency (pearson on common finite genes) ===")
for ctx,mid in ID.items():
    if mid not in dep_ids: print(ctx,"context line has NO DepMap CRISPR row -> skipped"); continue
print("(HeLa absent from 24Q4 CRISPRGeneEffect; Jurkat/Cal33 present)")
for ctx in "AC":
    mid=ID[ctx]; d=D[dep_ids.index(mid)]
    rs=[]
    for m2 in dep_ids:
        e=D[dep_ids.index(m2)]; m=np.isfinite(d)&np.isfinite(e)
        rs.append((NAMES.get(m2,m2),float(np.corrcoef(d[m],e[m])[0,1])))
    rs.sort(key=lambda t:-t[1])
    print(ctx,["%s %.4f"%t for t in rs])

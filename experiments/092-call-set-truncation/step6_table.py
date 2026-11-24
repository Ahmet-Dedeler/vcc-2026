import os, numpy as np
S = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]; pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
A=0.05; G=len(genes); gpos={g:i for i,g in enumerate(genes)}
fr=np.empty(G,np.int64); fr[np.argsort(genes)]=np.arange(G)
rows=[]; marg=[]
for i,t in enumerate(perts):
    ok=np.ones(G,bool); j=gpos.get(t)
    if j is not None: ok[j]=False
    rs=ok&(rq[i]<A); ps=ok&(pq[i]<A); ind=ok&(rl[i]!=0)&~np.isnan(rl[i])
    idx=np.flatnonzero(ps); key=np.lexsort((fr[idx],-np.abs(pl[i][idx]),pp[i][idx],pq[i][idx])); r=idx[key]
    hit=ind[r]&(pl[i][r]!=0)&(np.sign(pl[i][r])==np.sign(rl[i][r]))
    nr=int(rs.sum()); npd=int((ps&ind).sum()); k=int(hit.sum())
    inter=int((rs&ps).sum()); union=int((rs|ps).sum())
    rows.append((t,nr,npd,npd/nr,k,k/npd,k/max(npd,nr),inter/union))
    ch=np.cumsum(hit); n=r.size
    marg.append([ (int(ch[min(int(b*nr),n)-1]) if int(b*nr)>=1 else 0) for b in (0,) ])
rows.sort(key=lambda x:-x[1])
print(f"{'target':<10}{'n_real':>7}{'n_pred':>8}{'c':>8}{'k':>7}{'k/n_pred':>10}{'fid':>8}{'jac':>8}")
for r in rows: print(f"{r[0]:<10}{r[1]:>7}{r[2]:>8}{r[3]:>8.1f}{r[4]:>7}{r[5]:>10.4f}{r[6]:>8.4f}{r[7]:>8.4f}")

# marginal directional accuracy by rank band (pooled over perts, bands in units of n_real)
bands=[(0,0.25),(0.25,0.5),(0.5,1.0),(1.0,2.0),(2.0,5.0),(5.0,10.0),(10.0,1e9)]
tot=np.zeros(len(bands)); hits=np.zeros(len(bands))
for i,t in enumerate(perts):
    ok=np.ones(G,bool); j=gpos.get(t)
    if j is not None: ok[j]=False
    rs=ok&(rq[i]<A); ps=ok&(pq[i]<A); ind=ok&(rl[i]!=0)&~np.isnan(rl[i])
    idx=np.flatnonzero(ps); key=np.lexsort((fr[idx],-np.abs(pl[i][idx]),pp[i][idx],pq[i][idx])); r=idx[key]
    hit=(ind[r]&(pl[i][r]!=0)&(np.sign(pl[i][r])==np.sign(rl[i][r]))).astype(float)
    nr=int(rs.sum())
    for b,(lo,hi) in enumerate(bands):
        a=int(min(lo*nr,r.size)); bb=int(min(hi*nr,r.size))
        if bb>a: tot[b]+=bb-a; hits[b]+=hit[a:bb].sum()
print("\nPOOLED marginal directional accuracy by rank band (band edges in units of n_real):")
for b,(lo,hi) in enumerate(bands):
    lab=f"{lo:g}-{hi:g}" if hi<1e8 else f">{lo:g}"
    print(f"  rank {lab:>10} x n_real : n={int(tot[b]):>7}  acc={hits[b]/max(tot[b],1):.4f}")

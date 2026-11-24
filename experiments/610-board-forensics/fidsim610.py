# 610 (board forensics, 2026-10-02): fid on the h1-test truth (85 targets x 400 cells, 18,400 controls: the board shape) for a
# fixed call-propensity block of K genes (control-CPM order) under four sign sources: all-down, CLEAN4 (391 prior), the
# context own common response (leave-one-target-out P_down over the other 84 targets; an oracle for A/B/C), and the truth.
# Uses experiments/213-transport-law/work/de_real_h1.npz (the real DE table) and data/harness/h1-test.h5ad (controls only).
# Usage: uv run python experiments/610-board-forensics/fidsim610.py
import numpy as np, anndata as ad, scipy.sparse as sp
d=np.load('experiments/213-transport-law/work/de_real_h1.npz', allow_pickle=True)
genes=list(d['genes']); perts=list(d['perts']); padj=d['real_padj'].copy(); lfc=d['real_lfc'].copy()
P,G=lfc.shape; gi={g:i for i,g in enumerate(genes)}
for p,t in enumerate(perts):
    if t in gi: padj[p,gi[t]]=1.0; lfc[p,gi[t]]=np.nan
sig=padj<0.05; adj=np.isfinite(lfc)&(lfc!=0); nreal=sig.sum(1)
# control mean CPM
a=ad.read_h5ad('data/harness/h1-test.h5ad', backed='r')
ctrl=np.where(a.obs['target_gene'].values=='non-targeting')[0]
vn=list(a.var_names); col=np.array([vn.index(g) for g in genes])
X=a.X[ctrl.min():ctrl.max()+1] if False else a[ctrl].to_memory().X
X=sp.csr_matrix(X); L=np.asarray(X.sum(1)).ravel()
cpm=np.asarray(sp.diags(1e6/L)@X[:,col].tocsc().toarray() if False else (sp.diags(1e6/L)@X).tocsc()[:,col].mean(0)).ravel()
order=np.argsort(-cpm)
# sign sources (+1 = up, -1 = down)
c4=np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True)
c4g={g:i for i,g in enumerate(c4['genes'])}
pd_c4=np.array([c4['P_down'][c4g[g]] if g in c4g else np.nan for g in genes])
s_c4=np.where(np.isfinite(pd_c4), np.where(pd_c4>0.5,-1,1), -1)
down=(lfc<0)&adj; cnt=adj.sum(0)
def fid(calls_mask, signs):  # signs: P x G (+1/-1)
    F=[]
    for p in range(P):
        c=calls_mask[p]&adj[p]; npred=calls_mask[p].sum()
        k=((np.sign(lfc[p,c])==signs[p,c])).sum()
        if nreal[p]==0 and npred==0: continue
        if nreal[p]==0: F.append(k/npred); continue
        F.append(k/max(npred,nreal[p]))
    return np.mean(F)
truth=np.where(np.nan_to_num(lfc)<0,-1,1)
loo=np.zeros((P,G)); 
for p in range(P): loo[p]=(down.sum(0)-down[p])/np.maximum(cnt-adj[p],1)
s_loo=np.where(loo>0.5,-1,1)
S_all_down=-np.ones((P,G),int); S_c4=np.tile(s_c4,(P,1))
print('n_real median %d, share of targets with n_real>2200: %.2f'%(np.median(nreal),(nreal>2200).mean()))
print('%-6s %8s %8s %8s %8s'%('K','alldown','CLEAN4','ctxLOO','oracle'))
for K in [300,1000,2200,4000,6000,10779]:
    blk=np.zeros(G,bool); blk[order[:K]]=True; M=np.tile(blk,(P,1))
    print('%-6d %8.4f %8.4f %8.4f %8.4f'%(K,fid(M,S_all_down),fid(M,S_c4),fid(M,s_loo),fid(M,truth)))
# per-target oracle K = n_real with CLEAN4 signs on block-ordered genes (size matched)
M=np.zeros((P,G),bool)
for p in range(P): M[p,order[:max(nreal[p],300)]]=True
print('size-matched block K=max(n_real,300): CLEAN4 %.4f ctxLOO %.4f oracle %.4f'%(fid(M,S_c4),fid(M,s_loo),fid(M,truth)))

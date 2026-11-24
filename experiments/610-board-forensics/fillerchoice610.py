# 610 (2026-10-02): filler CHOICE on h1-test truth. Block of K genes chosen by (a) control CPM (call propensity, ours),
# (b) CLEAN4 sign confidence |P_down-0.5|, (c) a 50/50 mix; signs = CLEAN4. Reports fid and jac of the call set
# (jac = |R & calls| / |R | calls| per target, mean over 85; the target gene excluded), plus the ctx-LOO oracle signs.
import numpy as np
z=np.load('experiments/610-board-forensics/work/nullskew_h1test.npz', allow_pickle=True)
d=np.load('experiments/213-transport-law/work/de_real_h1.npz', allow_pickle=True)
genes=list(d['genes']); perts=list(d['perts']); padj=d['real_padj'].copy(); lfc=d['real_lfc'].copy()
assert list(z['genes'])==genes
P,G=lfc.shape; gi={g:i for i,g in enumerate(genes)}
for p,t in enumerate(perts):
    if t in gi: padj[p,gi[t]]=1.0; lfc[p,gi[t]]=np.nan
sig=padj<0.05; adj=np.isfinite(lfc)&(lfc!=0); nreal=sig.sum(1); cpm=z['cpm']
c4=np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True); c4g={g:i for i,g in enumerate(c4['genes'])}
pd4=np.array([c4['P_down'][c4g[g]] if g in c4g else np.nan for g in genes])
S4=np.tile(np.where(np.isfinite(pd4)&(pd4<=0.5),1,-1),(P,1))
down=(lfc<0)&adj; cnt=adj.sum(0)
loo=np.array([(down.sum(0)-down[p])/np.maximum(cnt-adj[p],1) for p in range(P)]); SL=np.where(loo>0.5,-1,1)
def fid(M,S):
    F=[]
    for p in range(P):
        c=M[p]&adj[p]; npred=M[p].sum(); k=(np.sign(lfc[p,c])==S[p,c]).sum()
        if nreal[p]==0 and npred==0: continue
        F.append(k/npred if nreal[p]==0 else k/max(npred,nreal[p]))
    return np.mean(F)
def jac(M):
    return np.mean([ ( (M[p]&sig[p]).sum() / max((M[p]|sig[p]).sum(),1) ) if (M[p]|sig[p]).sum()>0 else 1.0 for p in range(P)])
cpm_rank=np.argsort(-cpm); conf4=np.where(np.isfinite(pd4),np.abs(pd4-0.5),-1); c4_rank=np.argsort(-conf4)
print('%-5s %-28s %8s %8s %8s'%('K','block choice','fid(C4)','fid(LOO)','jac'))
for K in [1000,2200,4000]:
    for name,rank in [('CPM (call propensity)',cpm_rank),('CLEAN4 confidence',c4_rank)]:
        b=np.zeros(G,bool); b[rank[:K]]=True; M=np.tile(b,(P,1))
        print('%-5d %-28s %8.4f %8.4f %8.4f'%(K,name,fid(M,S4),fid(M,SL),jac(M)))
    b=np.zeros(G,bool); b[cpm_rank[:K//2]]=True
    extra=[g for g in c4_rank if not b[g]][:K-b.sum()]; b[extra]=True; M=np.tile(b,(P,1))
    print('%-5d %-28s %8.4f %8.4f %8.4f'%(K,'50/50 CPM + CLEAN4-conf',fid(M,S4),fid(M,SL),jac(M)))

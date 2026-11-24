# 610 (2026-10-02): does the NULL SIGN SKEW of a 400-cell group, estimated from a context's CONTROLS ALONE, predict the sign
# of each target's non-significant reference lfc?  Mechanism: the scorer's lfc = log2(mean CPM of 400 cells / mean CPM of
# controls); for a right-skewed per-cell distribution the 400-cell mean sits below the population mean more often than above,
# so unaffected genes read lfc < 0 more than half the time (and, since CPM rows sum to 1e6, abundant genes lean up).
# Instrument: h1-test truth (85 targets x 400 cells vs 18,400 controls = the board shape), real DE table from 213.
# Compares filler-sign sources on a fixed block of K genes: all-down, CLEAN4 (391), NULL-SKEW (controls only), ctx-LOO common
# response (oracle for A/B/C), truth.  Usage: uv run python experiments/610-board-forensics/nullskew610.py
import numpy as np, anndata as ad, scipy.sparse as sp, sys
rng=np.random.default_rng(610)
d=np.load('experiments/213-transport-law/work/de_real_h1.npz', allow_pickle=True)
genes=list(d['genes']); perts=list(d['perts']); padj=d['real_padj'].copy(); lfc=d['real_lfc'].copy()
P,G=lfc.shape; gi={g:i for i,g in enumerate(genes)}
for p,t in enumerate(perts):
    if t in gi: padj[p,gi[t]]=1.0; lfc[p,gi[t]]=np.nan
sig=padj<0.05; adj=np.isfinite(lfc)&(lfc!=0); nreal=sig.sum(1)
a=ad.read_h5ad('data/harness/h1-test.h5ad', backed='r')
ctrl=np.where(a.obs['target_gene'].values=='non-targeting')[0]
X=sp.csr_matrix(a[ctrl].to_memory().X).astype(np.float64)
L=np.asarray(X.sum(1)).ravel()
vn={g:i for i,g in enumerate(a.var_names)}; col=np.array([vn[g] for g in genes])
C=(sp.diags(1e6/L)@X).tocsc()[:,col].tocsr()   # per-cell CPM (normalised over ALL genes), restricted to the DE genes
N=C.shape[0]; tot=np.asarray(C.sum(0)).ravel()
B=int(sys.argv[1]) if len(sys.argv)>1 else 300
neg=np.zeros(G); pos=np.zeros(G); zsum=np.zeros(G)
for b in range(B):
    idx=rng.choice(N,400,replace=False)
    s=np.asarray(C[idx].sum(0)).ravel()
    m_g=s/400; m_c=(tot-s)/(N-400)          # leave-group-out control mean, as the scorer's pert-vs-control
    l=np.log2((m_g+1e-9)/(m_c+1e-9))
    neg+=(l<0); pos+=(l>0)
pnull=neg/np.maximum(neg+pos,1)
cpm=tot/N; order=np.argsort(-cpm)
print('null-skew P(lfc<0) over %d random 400-cell control groups: mean %.4f; genes with P>=0.55: %d, <=0.45: %d'%(B,pnull.mean(),(pnull>=0.55).sum(),(pnull<=0.45).sum()))
for lo,hi in [(0,300),(300,1000),(1000,2200),(2200,5000),(5000,G)]:
    print('  CPM rank %5d-%5d: mean P_null_down %.4f'%(lo,hi,pnull[order[lo:hi]].mean()))
c4=np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True); c4g={g:i for i,g in enumerate(c4['genes'])}
pd_c4=np.array([c4['P_down'][c4g[g]] if g in c4g else np.nan for g in genes])
down=(lfc<0)&adj; cnt=adj.sum(0)
loo=np.array([(down.sum(0)-down[p])/np.maximum(cnt-adj[p],1) for p in range(P)])
srcs={'alldown':-np.ones((P,G)),'CLEAN4':np.tile(np.where(np.isfinite(pd_c4)&(pd_c4<=0.5),1,-1),(P,1)),
      'NULLSKEW':np.tile(np.where(pnull>0.5,-1,1),(P,1)),'ctxLOO':np.where(loo>0.5,-1,1),'truth':np.where(np.nan_to_num(lfc)<0,-1,1)}
def fid(M,S):
    F=[]
    for p in range(P):
        c=M[p]&adj[p]; npred=M[p].sum(); k=(np.sign(lfc[p,c])==S[p,c]).sum()
        if nreal[p]==0 and npred==0: continue
        F.append(k/npred if nreal[p]==0 else k/max(npred,nreal[p]))
    return np.mean(F)
def prec_nonDE(S, M):
    return np.mean([ (np.sign(lfc[p,M[p]&adj[p]&~sig[p]])==S[p,M[p]&adj[p]&~sig[p]]).mean() for p in range(P)])
print('\nfid on h1-test truth, fixed block of the K highest-CPM genes; [precision on the block genes that are NOT DE]')
print('%-6s'%'K'+''.join('%20s'%k for k in srcs))
for K in [1000,2200,4000,10779]:
    blk=np.zeros(G,bool); blk[order[:K]]=True; M=np.tile(blk,(P,1))
    print('%-6d'%K+''.join('%11.4f [%.3f]'%(fid(M,S),prec_nonDE(S,M)) for S in srcs.values()))
conf=np.abs(pnull-0.5); oc=np.argsort(-conf)
for K in [1000,2200,4000]:
    blk=np.zeros(G,bool); blk[oc[:K]]=True; M=np.tile(blk,(P,1))
    print('block = top-%d by NULL-SKEW confidence: NULLSKEW fid %.4f [%.3f]  ctxLOO %.4f  CLEAN4 %.4f'%(K,fid(M,srcs['NULLSKEW']),prec_nonDE(srcs['NULLSKEW'],M),fid(M,srcs['ctxLOO']),fid(M,srcs['CLEAN4'])))
np.savez_compressed('experiments/610-board-forensics/work/nullskew_h1test.npz', genes=np.array(genes), pnull=pnull, cpm=cpm, B=B)

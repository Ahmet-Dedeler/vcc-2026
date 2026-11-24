# Per-target amplitude allocation under the ratio-of-sums MSE, on the 323 H1 harness outputs.
import numpy as np, json, h5py, scipy.sparse as sp, time
t0=time.time()
W='experiments/323-current-family-validation/work/'
z=np.load(W+'effects.npz', allow_pickle=True)
genes=z['genes']; targets=z['targets']; D=z['cap15']            # design effects (delivered log-delta space, per 318)
m=json.load(open(W+'metrics-cap15.json')); num=m['mse_legs']['numerator']; den=m['mse_legs']['denominator']
# ---- truth from h1-test cells
with h5py.File('data/harness/h1-test.h5ad') as f:
    X=f['X']; data=X['data'][:]; indices=X['indices'][:]; indptr=X['indptr'][:]
    n_cells,n_genes=X.attrs['shape']
    hg=f["var"]["_index"]["values"][:].astype(str)
    tg=f['obs']['target_gene']; cats=tg['categories'][:].astype(str); codes=tg['codes'][:]
Xs=sp.csr_matrix((data,indices,indptr),shape=(n_cells,n_genes)); del data,indices,indptr
TS=5e4
def profile(mask):
    P=np.asarray(Xs[mask].sum(axis=0)).ravel(); return np.log1p(TS*P/P.sum())
ctrl_name=[c for c in cats if 'non-targeting' in c or c=='control' or c=='NTC'][0]
b_ctrl=profile(codes==np.where(cats==ctrl_name)[0][0])
# gene alignment: design genes -> h1 genes
hidx={g:i for i,g in enumerate(hg)}
keep=[i for i,g in enumerate(genes) if g in hidx]; hcol=np.array([hidx[genes[i]] for i in keep]); keep=np.array(keep)
print(f"genes aligned {len(keep)} of {len(genes)} design / {len(hg)} h1; ctrl cells {(codes==np.where(cats==ctrl_name)[0][0]).sum()}; load {time.time()-t0:.0f}s")
rows=[]
for k,t in enumerate(targets):
    if t not in cats: continue
    ci=np.where(cats==t)[0][0]; delta=profile(codes==ci)-b_ctrl
    d=D[k,keep]; dl=delta[hcol]
    # exclude own target gene
    own=np.where(genes[keep]==t)[0]
    if len(own): d=d.copy(); dl=dl.copy(); d[own]=0; dl[own]=0
    rows.append(dict(t=t, dd=float(d@d), tt=float(dl@dl), dt=float(d@dl), N=num.get(t,np.nan), Dn=den.get(t,np.nan)))
R={k:np.array([r[k] for r in rows]) for k in ('dd','tt','dt','N','Dn')}; T=[r['t'] for r in rows]
cos=R['dt']/np.sqrt(R['dd']*R['tt']); nd=np.sqrt(R['dd']); nt=np.sqrt(R['tt'])
print(f"\nn targets {len(T)}; harness raw mse (sum N / sum D) = {np.nansum(R['N'])/np.nansum(R['Dn']):.5f}")
# ---- effect-size concentration
Dn=R['Dn']; order=np.argsort(-Dn); cum=np.cumsum(Dn[order])/Dn.sum()
print("share of sum D_p held by top 5/10/20/40 targets:", [round(float(cum[i-1]),3) for i in (5,10,20,40)])
print("per-target N_p/D_p: median %.3f, D-weighted %.3f; top-10 targets' N/D: %s"%(np.nanmedian(R['N']/Dn), np.nansum(R['N'][order[:10]])/np.nansum(Dn[order[:10]]), np.round((R['N']/Dn)[order[:10]],2).tolist()))
print("cos_p: median %.3f mean %.3f min %.3f max %.3f; D-weighted mean cos %.3f"%(np.median(cos),cos.mean(),cos.min(),cos.max(), (cos*Dn).sum()/Dn.sum()))
# ---- scale structure
from scipy.stats import spearmanr, pearsonr
print("design norm vs truth norm: spearman %.3f pearson %.3f | design norm vs sqrt(D_p): spearman %.3f"%(spearmanr(nd,nt)[0], pearsonr(nd,nt)[0], spearmanr(nd,np.sqrt(np.maximum(Dn,1e-12)))[0]))
print("truth norm range: %.2f-%.2f (median %.2f); design norm range %.2f-%.2f (median %.2f)"%(nt.min(),nt.max(),np.median(nt),nd.min(),nd.max(),np.median(nd)))
# ---- allocation counterfactuals in the signal part: S = sum_p |a_p d_p - delta_p|^2 ; delta MSE = dS / sum D
def S(a): return float(np.sum(a*a*R['dd'] - 2*a*R['dt'] + R['tt']))
a_glob=R['dt'].sum()/R['dd'].sum()
a_or=R['dt']/R['dd']                       # oracle per-target
sumD=np.nansum(Dn)
base=S(np.full(len(T),a_glob))
print(f"\nglobal alpha* = {a_glob:.4f} (design already delivered at alpha~1)")
print(f"signal MSE term S/sumD: global-alpha {base/sumD:.4f} | alpha=1 {S(np.ones(len(T)))/sumD:.4f}")
print(f"ORACLE per-target alpha: S/sumD {S(a_or)/sumD:.4f}  -> delta raw mse {(S(a_or)-base)/sumD:+.4f}  (scaled ~ x{1/0.955:.2f}, overall /6 -> {-(S(a_or)-base)/sumD/0.955/6:+.4f})")
# proxy allocations: a_p = k * f_p with f from predictable quantities, k fit globally
def proxy(f,label):
    k=(f*R['dt']).sum()/((f*f)*R['dd']).sum(); a=k*f; print(f"proxy {label:38s}: S/sumD {S(a)/sumD:.4f} -> delta raw mse {(S(a)-base)/sumD:+.4f}")
proxy(np.ones(len(T)),'uniform (=global)')
proxy(nt/nd,'ORACLE norm match (|truth|/|design|)')
proxy(1/nd,'equalize design norms (a_p ∝ 1/|d_p|)')
proxy(nd/nd.mean(),'amplify by design norm (a_p ∝ |d_p|)')
proxy(nd**0.5/np.sqrt(nd.mean()),'a_p ∝ |d_p|^0.5')
# the truth's norm rank as a mild proxy check
proxy(np.sqrt(np.maximum(Dn,1e-12))/nd,'a_p ∝ sqrt(D_p)/|d_p| (oracle-ish)')
# how far are we from optimal per target: distribution of a_or
print("\noracle a_p quantiles (10/25/50/75/90):", np.round(np.quantile(a_or,[.1,.25,.5,.75,.9]),3).tolist())
print("oracle a_p vs design norm spearman %.3f; oracle a_p vs truth norm spearman %.3f; oracle a_p vs cos spearman %.3f"%(spearmanr(a_or,nd)[0],spearmanr(a_or,nt)[0],spearmanr(a_or,cos)[0]))
np.savez('/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/0acea69f-98f2-413e-8a3f-27adb48f049a/scratchpad/mse_alloc.npz', T=np.array(T), **R, cos=cos, nd=nd, nt=nt, a_or=a_or)
print(f"done {time.time()-t0:.0f}s")

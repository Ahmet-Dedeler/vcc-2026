import sys, time, numpy as np, anndata as ad
from pathlib import Path
ROOT=Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0,str(ROOT/"experiments/196-emitted-depth")); sys.path.insert(0,str(ROOT/"src"))
import recon196, vcc2026.rc_emit as RC
from vcc2026.donor_transfer import emit_counts
t0=time.time()
D=recon196.load_chal("A")
A=ad.read_h5ad(ROOT/"data/challenge/context_A.h5ad"); X=A.X.tocsr(); del A
cb=RC.control_block_from_counts(X, min_cpm=5.0); del X
scorer, builder = cb.keep, D["gate"]; G=D["genes"].size; n=400
for i in (0,15,30):
    Rec=recon196.Recon(D["sums"][i], D["r"], builder, n=n, seed=1000+i)
    M1=D["sums"][i].astype(np.int64); pcpm=1e6*M1/M1.sum()
    eff=np.abs(np.log2((pcpm+1e-9)/(D["cpm"]+1e-9)))
    order=np.lexsort((np.arange(G), -np.where(scorer, eff, -np.inf)))
    sv=np.full(G,1.0); sv[order[:8000]]=0.25
    for m in (1,2,4,8,16):
        prof=Rec.profile(m); keyv=np.where(builder, np.log1p(50000.0*prof)-D["b"], 0.0)
        try:
            B,dg=RC.rc_hybrid_block(prof, np.full(n,Rec.emit_L(m)), 0.1,
                np.random.default_rng(90000+i), pass2_rng=np.random.default_rng(90000+i),
                spread=sv, overdispersion=4.0, control=cb, key=keyv, jitter_seed=1000+i,
                emit_counts=emit_counts)
            print(f"  t{i} m{m}: OK  nnz={int((B>0).sum())} maxcell={int(B.sum(1).max())} "
                  f"forced={dg['forced']} resid={dg['max_abs_resid']} moves={dg['moves']} "
                  f"unres={dg['unresolved']} subgate/cell={dg['sub_gate_counts_per_cell']:.1f} "
                  f"clip={dg['clip_frac']:.4f} [{time.time()-t0:.0f}s]", flush=True)
            del B
        except Exception as e:
            print(f"  t{i} m{m}: FAIL {type(e).__name__}: {e} [{time.time()-t0:.0f}s]", flush=True)

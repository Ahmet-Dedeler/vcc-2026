"""196 falsifier 1 — the DEEPEST LEGAL depth, measured rather than extrapolated."""
import sys, time, numpy as np
from pathlib import Path
ROOT=Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0,str(ROOT/"experiments/196-emitted-depth")); sys.path.insert(0,str(ROOT/"src"))
import recon196
from vcc2026.donor_transfer import emit_counts
CAP, CELL, n = 4.75e9, 1_000_000, 400
t0=time.time()
D=recon196.load_chal("A"); G=D["genes"].size
base=None
print(f"{'depth':>7}{'emit L':>10}{'nnz/block':>12}{'ratio':>8}{'package':>16}{'% cap':>8}"
      f"{'max cell':>11}{'% cell cap':>12}")
for m in (1,16,32,45,64,96):
    nn=[]; mx=0; Ls=[]
    for i in (0,15,30):
        R=recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000+i)
        M1=D["sums"][i].astype(np.int64); pcpm=1e6*M1/M1.sum()
        eff=np.abs(np.log2((pcpm+1e-9)/(D["cpm"]+1e-9)))
        order=np.lexsort((np.arange(G), -np.where(D["gate"], eff, -np.inf)))
        sv=np.full(G,1.0); sv[order[:8000]]=0.25
        Lm=R.emit_L(m); Ls.append(Lm)
        v,c,ip=emit_counts(R.profile(m), np.full(n,Lm), 0.1, np.random.default_rng(90000+i),
                           mode="quasi", spread=sv, overdispersion=4.0)
        rows=np.add.reduceat(v.astype(np.float64), np.asarray(ip)[:-1])
        nn.append(v.size); mx=max(mx,int(rows.max()))
    a=float(np.mean(nn)); base=base or a
    pkg=2_377_814_190*a/base
    print(f"{m:>7}{int(np.mean(Ls)):>10,}{a:>12,.0f}{a/base:>8.4f}{pkg:>16,.0f}"
          f"{100*pkg/CAP:>7.2f}%{mx:>11,}{100*mx/CELL:>11.2f}%", flush=True)
print(f"done [{time.time()-t0:.0f}s]")

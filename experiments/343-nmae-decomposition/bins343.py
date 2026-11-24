"""343c -- sign accuracy on the reference gate by |lfc| bin and by recipient basal CPM bin (cap15, H1)."""
import numpy as np, polars as pl, json
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]; W = REPO/"experiments/323-current-family-validation/work"; OUT = Path(__file__).resolve().parent/"work"
s = np.load(W/"stats-cap15.npz"); r = np.load(W/"reference.npz")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]; targets = s["targets"].astype(str); gi = {g:i for i,g in enumerate(genes)}
basal = np.asarray(r["ref_lin"], dtype=float) if "ref_lin" in r.files else None
print("reference.npz arrays:", r.files)
ref = pl.read_parquet(W/"reference_de.parquet"); lhat = np.nan_to_num(np.asarray(s["lfc"], float)); qhat = np.asarray(s["q"], float)
L, LH, B, T, C = [], [], [], [], []
for ti,t in enumerate(targets):
    sub = ref.filter(pl.col("target")==t); g=sub["feature"].to_numpy().astype(str); q=sub["p_adj"].to_numpy(); l=sub["log2_fold_change"].to_numpy().astype(float)
    gate=(q<.05)&np.isfinite(l)&(l!=0)&(g!=t)
    if gate.sum()<10: continue
    idx=np.array([gi[x] for x in g[gate]]); L.append(l[gate]); LH.append(lhat[ti,idx]); C.append(qhat[ti,idx]<.05); T.append(np.full(gate.sum(),ti))
    if basal is not None: B.append(basal[idx])
L=np.concatenate(L); LH=np.concatenate(LH); C=np.concatenate(C); T=np.concatenate(T); ok=np.sign(LH)==np.sign(L); w=np.abs(L)
print(f"pooled gate genes {len(L):,} over {len(set(T))} targets; unweighted acc {ok.mean():.3f}; weighted {(w*ok).sum()/w.sum():.3f}; up-share of reference {np.mean(L>0):.3f}; our up-share {np.mean(LH>0):.3f}")
def table(x, name, nb=5):
    qs=np.quantile(x, np.linspace(0,1,nb+1)); rows=[]
    for i in range(nb):
        m=(x>=qs[i])&(x<=qs[i+1]) if i==nb-1 else (x>=qs[i])&(x<qs[i+1])
        rows.append((f"{qs[i]:.3g}–{qs[i+1]:.3g}", int(m.sum()), ok[m].mean(), (w[m]*ok[m]).sum()/w[m].sum(), w[m].sum()/w.sum(), C[m].mean(), np.abs(LH[m]).sum()/w[m].sum()))
    print(f"\n{name} quintile | n | acc | w-acc | share of NMAE mass | our-called frac | magnitude ratio")
    for r_ in rows: print(f"{r_[0]:>14s} | {r_[1]:7d} | {r_[2]:.3f} | {r_[3]:.3f} | {r_[4]:.3f} | {r_[5]:.3f} | {r_[6]:.3f}")
table(np.abs(L), "|ref lfc|")
if basal is not None: table(np.log10(B_:=np.concatenate(B)+1e-3), "log10 recipient basal CPM")
# sign accuracy split by reference direction
for d,name in ((L>0,"ref UP"),(L<0,"ref DOWN")):
    print(f"{name}: n={d.sum():,} acc={ok[d].mean():.3f} w-acc={(w[d]*ok[d]).sum()/w[d].sum():.3f} our mean lfc={LH[d].mean():+.3f}")

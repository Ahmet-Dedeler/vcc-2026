"""343b -- where the gate-gene signs come from: designed effect vs delivered lfc, called vs uncalled."""
import json, sys, numpy as np, polars as pl
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]; W = REPO/"experiments/323-current-family-validation/work"; OUT = Path(__file__).resolve().parent/"work"
ARM = sys.argv[1] if len(sys.argv) > 1 else "cap15"
s = np.load(W/f"stats-{ARM}.npz"); r = np.load(W/"reference.npz"); e = np.load(W/"effects.npz")
genes_all = r["genes"].astype(str); keep = r["keep"].astype(bool); genes = genes_all[keep]; targets = s["targets"].astype(str)
assert (e["targets"].astype(str) == targets).all()
eg = {g:i for i,g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes]); miss = col < 0; print("DE-universe genes absent from effects axis:", int(miss.sum()))
eff = np.asarray(e[ARM], dtype=np.float64)[:, np.where(miss, 0, col)]; eff[:, miss] = 0.   # sign 0 never matches: counted as design miss
lhat = np.asarray(s["lfc"], dtype=np.float64); qhat = np.asarray(s["q"], dtype=np.float64); phat = np.asarray(s["p"], dtype=np.float64)
ref = pl.read_parquet(W/"reference_de.parquet"); gi = {g:i for i,g in enumerate(genes)}
acc = {k: [] for k in ("w_called","w_uncalled","w_all","design_w_all","design_w_called","design_w_uncalled","agree_design_delivered_w","n_gate","mass_called",
                       "u_called","u_uncalled","design_u_all","nmae_best_scalar_grid","best_a","reach_like")}
grid = np.concatenate([np.linspace(0.25,1,4), np.linspace(1.25,6,20)])
for ti,t in enumerate(targets):
    sub = ref.filter(pl.col("target")==t); g=sub["feature"].to_numpy().astype(str); q=sub["p_adj"].to_numpy(); l=sub["log2_fold_change"].to_numpy().astype(float)
    gate=(q<.05)&np.isfinite(l)&(l!=0)&(g!=t)
    if gate.sum()<10: continue
    idx=np.array([gi[x] for x in g[gate]]); l=l[gate]; w=np.abs(l); D=w.sum()
    lh=np.nan_to_num(lhat[ti,idx]); ef=eff[ti,idx]; called=qhat[ti,idx]<.05
    sd=np.sign(lh)==np.sign(l); se=np.sign(ef)==np.sign(l); ag=np.sign(ef)==np.sign(lh)
    acc["n_gate"].append(int(gate.sum())); acc["mass_called"].append(w[called].sum()/D)
    acc["w_all"].append((w*sd).sum()/D); acc["design_w_all"].append((w*se).sum()/D); acc["agree_design_delivered_w"].append((w*ag).sum()/D)
    acc["u_called"].append(sd[called].mean() if called.any() else np.nan); acc["u_uncalled"].append(sd[~called].mean() if (~called).any() else np.nan)
    acc["design_u_all"].append(se.mean())
    acc["w_called"].append((w[called]*sd[called]).sum()/max(w[called].sum(),1e-12) if called.any() else np.nan)
    acc["w_uncalled"].append((w[~called]*sd[~called]).sum()/max(w[~called].sum(),1e-12) if (~called).any() else np.nan)
    acc["design_w_called"].append((w[called]*se[called]).sum()/max(w[called].sum(),1e-12) if called.any() else np.nan)
    acc["design_w_uncalled"].append((w[~called]*se[~called]).sum()/max(w[~called].sum(),1e-12) if (~called).any() else np.nan)
    vals=[np.abs(a*lh-l).sum()/D for a in grid]; j=int(np.argmin(vals)); acc["nmae_best_scalar_grid"].append(vals[j]); acc["best_a"].append(float(grid[j]))
    # reach-like: order gate genes by our confidence (called first, then p asc, then |lfc| desc); deepest prefix with >=.9 sign accuracy
    order=np.lexsort((-np.abs(lh), phat[ti,idx], ~called)); c=np.cumsum(sd[order]); k=np.arange(1,len(order)+1); ok=np.where(c/k>=.9)[0]
    acc["reach_like"].append((ok[-1]+1)/len(order) if ok.size else 0.)
summ={k: float(np.nanmean(v)) for k,v in acc.items()}
summ["n_targets"]=len(acc["n_gate"]); summ["arm"]=ARM
(OUT/f"signs-{ARM}.json").write_text(json.dumps(summ,indent=1)); print(json.dumps(summ,indent=1))

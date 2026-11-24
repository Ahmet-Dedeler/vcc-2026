"""KA-1b: distinct guide ids per (target, condition) over the 4 donors vs DE_stats.n_guides.
KA-6 : my across-donor cosine vs the authors' own donor_correlation_all_mean (independent estimator)."""
import numpy as np, pandas as pd, h5py, json
from scipy.stats import spearmanr
ROOT = "~/Code/virtual-cell-challenge-2026"
z = np.load("obs280.npz", allow_pickle=True)
obs = pd.DataFrame({k: z[k] for k in ("cond","gtype","pert","guide","keep_for_DE","keep_eff","keep_min")})
for c in ("cond","gtype","pert","guide"): obs[c]=obs[c].astype(str)
panel = sorted(set(l.split(",")[0] for l in open(f"{ROOT}/data/challenge/pert_counts.csv").read().strip().split("\n")[1:]))
de = pd.read_csv("de_stats_panel280.csv")
out={}
for lab, m in (("all guides", np.ones(len(obs),bool)), ("keep_for_DE", obs.keep_for_DE.values),
               ("keep_eff", obs.keep_eff.values), ("keep_min", obs.keep_min.values)):
    tg = obs[(obs.gtype=="targeting") & obs.pert.isin(panel) & m]
    g = tg.groupby(["pert","cond"]).guide.nunique().reset_index(name="ng")
    j = de.merge(g, left_on=["target","cond"], right_on=["pert","cond"])
    d = np.abs(j.n_guides.values - j.ng.values)
    print(f"  KA-1b [{lab:<12}] distinct guides: exact {int((d==0).sum())}/{len(j)}  max|d| {d.max():.0f}")
    out[lab]=dict(n=int(len(j)), exact=int((d==0).sum()), max_abs=float(d.max()))
print()
res={}
for cond in ("Rest","Stim8hr","Stim48hr"):
    zz = np.load(f"q1_pertarget_{cond}_keepDE.npz", allow_pickle=True)
    df = pd.DataFrame(dict(target=[str(x) for x in zz["targets"]], cos=zz["cos"]))
    j = df.merge(de[de.cond==cond][["target","donor_correlation_all_mean","donor_correlation_hits_mean"]], on="target")
    for col in ("donor_correlation_all_mean","donor_correlation_hits_mean"):
        v=j[col].values.astype(float); m=np.isfinite(v)&np.isfinite(j.cos.values)
        rho,p = spearmanr(v[m], j.cos.values[m])
        print(f"  KA-6 {cond:<9} Spearman(my across-donor cosine, authors' {col}) = {rho:+.3f} (p={p:.2g}, n={m.sum()})")
        res[f"{cond}:{col}"]=dict(rho=float(rho), p=float(p), n=int(m.sum()))
json.dump(dict(KA1b=out, KA6=res), open("ka6_280.json","w"), indent=1)

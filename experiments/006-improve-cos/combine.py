"""Follow-up: sweep the target-expression gate and combine it with per-gene normalisation."""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TS
CONTROL, PC = "non-targeting", "target_gene"

def load(path):
    a = ad.read_h5ad(path); X = a.X.tocsr()
    lab = a.obs[PC].astype(str).to_numpy(); genes = np.asarray(a.var_names)
    Pc = np.asarray(X[lab==CONTROL].sum(0)).ravel().astype(np.float64)
    b = np.log1p(TS*Pc/Pc.sum())
    ps, ds = [], []
    for g in sorted(set(lab)-{CONTROL}):
        P = np.asarray(X[lab==g].sum(0)).ravel().astype(np.float64)
        ds.append(np.log1p(TS*P/P.sum())-b); ps.append(g)
    return np.asarray(ps), np.vstack(ds), b, genes

donor = Donor.load("data/donors/k562_gwps")
m, s = donor.decompose()
didx = {p:i for i,p in enumerate(donor.perts)}

def cosm(P,Q):
    Pc=P-P.mean(0,keepdims=True); Qc=Q-Q.mean(0,keepdims=True)
    na,nb=np.linalg.norm(Pc,axis=1),np.linalg.norm(Qc,axis=1); ok=(na>0)&(nb>0)
    return float(np.mean(np.sum(Pc[ok]*Qc[ok],axis=1)/(na[ok]*nb[ok])))

for name, path in [("h1-test", "data/harness/h1-test.h5ad"),
                   ("h1-trainval", "data/harness/h1-trainval.h5ad")]:
    ps, true, b_ctrl, genes = load(path)
    gi={g:i for i,g in enumerate(genes)}
    tc,dc=[],[]
    for j,g in enumerate(donor.genes):
        t=gi.get(g)
        if t is not None: tc.append(t); dc.append(j)
    tc,dc=np.asarray(tc),np.asarray(dc)
    cov=np.asarray([p in didx for p in ps]); rows=[didx[p] for p in ps[cov]]
    T=true[cov][:,tc]; S=s[np.ix_(rows,dc)]; S_all=s[:,dc]
    cpm = 1e6*np.expm1(b_ctrl)/np.expm1(b_ctrl).sum()
    cpm_t = cpm[tc]
    mu = np.abs(S_all).mean(0)+1e-9
    print(f"\n=== {name} ({len(rows)} perts) ===")
    print(f"  {'gate':>14}  {'plain':>8} {'/mean|e|':>9} {'/sd':>8}")
    sd = S_all.std(0)+1e-9
    base = cosm(S, T)
    for thr in (0, 1, 5, 10, 20, 50, 100, 200):
        g = (cpm_t > thr)
        if g.sum() < 200: break
        a1 = cosm(S*g, T); a2 = cosm(S/mu*g, T); a3 = cosm(S/sd*g, T)
        star = "  <--" if a2 == max(a1,a2,a3) and thr>0 else ""
        print(f"  CPM > {thr:<8} {a1:8.4f} {a2:9.4f} {a3:8.4f}  (n={g.sum()}){star}")
    print(f"  baseline (no gate, plain): {base:.4f}")

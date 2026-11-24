"""Experiment 006 — can any cheap transform raise the transfer cosine?

Two attempts failed in exp 004 (low-rank truncation, supervised per-gene weights). Both
were reweightings of a fixed donor estimate. This tries transforms of a different kind:
changing which genes are allowed to dominate the direction.

Cosine is driven by whichever genes carry the largest values. A gene that is noisy in the
donor inflates ||pred|| without contributing signal. Standardising per gene, or restricting
to genes that respond reproducibly, changes which genes set the direction.

Offline against the h1-test harness -- seconds per variant, no cell-eval2.
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TS
CONTROL, PC = "non-targeting", "target_gene"

a = ad.read_h5ad("data/harness/h1-test.h5ad"); X = a.X.tocsr()
lab = a.obs[PC].astype(str).to_numpy(); genes = np.asarray(a.var_names)
Pc = np.asarray(X[lab == CONTROL].sum(0)).ravel().astype(np.float64)
b_ctrl = np.log1p(TS * Pc / Pc.sum())
ps, ds = [], []
for g in sorted(set(lab) - {CONTROL}):
    P = np.asarray(X[lab == g].sum(0)).ravel().astype(np.float64)
    ds.append(np.log1p(TS * P / P.sum()) - b_ctrl); ps.append(g)
true = np.vstack(ds); ps = np.asarray(ps)

donor = Donor.load("data/donors/k562_gwps")
gi = {g: i for i, g in enumerate(genes)}
tc, dc = [], []
for j, g in enumerate(donor.genes):
    t = gi.get(g)
    if t is not None: tc.append(t); dc.append(j)
tc, dc = np.asarray(tc), np.asarray(dc)
didx = {p: i for i, p in enumerate(donor.perts)}
cov = np.asarray([p in didx for p in ps])
rows = [didx[p] for p in ps[cov]]
T = true[cov][:, tc]
m, s = donor.decompose()
S = s[np.ix_(rows, dc)]
S_all = s[:, dc]
print(f"{len(rows)} perts x {len(tc)} shared genes\n")

def cosm(P, Q):
    Pc_ = P - P.mean(0, keepdims=True); Qc = Q - Q.mean(0, keepdims=True)
    na, nb = np.linalg.norm(Pc_, axis=1), np.linalg.norm(Qc, axis=1)
    ok = (na > 0) & (nb > 0)
    return float(np.mean(np.sum(Pc_[ok]*Qc[ok], axis=1)/(na[ok]*nb[ok])))

res = {}
def run(name, P):
    c = cosm(P, T); res[name] = c
    print(f"  {name:44} {c:+.4f}   {c - res.get('baseline', c):+.4f}")

run("baseline", S)
sd = S_all.std(0) + 1e-9
run("gene z-score (divide by donor sd)", S / sd)
run("gene sd^-0.5 (half standardisation)", S / np.sqrt(sd))
run("gene sd^+0.5 (anti-standardisation)", S * np.sqrt(sd))
run("divide by mean |effect| per gene", S / (np.abs(S_all).mean(0) + 1e-9))

def ranknorm(M):
    r = np.argsort(np.argsort(M, axis=1), axis=1).astype(np.float64)
    return (r - r.mean(1, keepdims=True)) / (r.std(1, keepdims=True) + 1e-9)
run("rank-normalise per perturbation", ranknorm(S))

for q in (0.5, 0.75, 0.9):
    run(f"keep top {int((1-q)*100)}% most-variable genes", S * (sd > np.quantile(sd, q)))

cpm = 1e6 * np.expm1(b_ctrl) / np.expm1(b_ctrl).sum()
for thr in (1.0, 5.0, 20.0):
    run(f"keep target CPM > {thr}", S * (cpm[tc] > thr))

run("clip at +-3 gene sd", np.clip(S, -3*sd, 3*sd))
run("sign only (binarised)", np.sign(S))

best = max(res, key=res.get)
print(f"\nbest: {best}\n  {res[best]:+.4f} vs baseline {res['baseline']:+.4f} "
      f"({100*(res[best]-res['baseline'])/abs(res['baseline']):+.1f}%)")
Path("experiments/006-improve-cos/variants.json").write_text(json.dumps(res, indent=2))

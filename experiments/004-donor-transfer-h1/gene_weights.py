"""Can we learn WHICH genes transfer, and does that generalise to an unseen line?

Low-rank denoising failed (see lowrank.py), so the transferable signal is not a few
shared programs. The next hypothesis is per-GENE: some genes' responses are portable
across cell contexts (core machinery) and others are context-specific. If that is a
property of the gene rather than of the target line, a weight vector `w` fitted on one
donor->target pair should improve a *different* donor->target pair. That is the only
form we could actually use in 2026, where no paired data for A/B/C exists.

Fit: per gene j, the least-squares slope of true-target effect on donor effect,
     w_j = sum_p (delta_target[p,j] * delta_donor[p,j]) / sum_p delta_donor[p,j]^2
Apply: pred = w * delta_donor.

Three tests, in increasing honesty:
  1. w fitted on H1 train+val, evaluated on H1 test      (same line, held-out perts)
  2. w fitted on K562->RPE1, evaluated on K562->H1       (UNSEEN LINE -- the real test)
  3. shuffled-gene control                                (null)
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad, scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TransferModel, TS
CONTROL, PERT_COL = "non-targeting", "target_gene"

def target_deltas(path):
    a = ad.read_h5ad(path); X = a.X.tocsr()
    lab = a.obs[PERT_COL].astype(str).to_numpy(); genes = np.asarray(a.var_names)
    Pc = np.asarray(X[lab==CONTROL].sum(0)).ravel().astype(np.float64)
    b = np.log1p(TS*Pc/Pc.sum())
    ps, ds = [], []
    for g in sorted(set(lab)-{CONTROL}):
        P = np.asarray(X[lab==g].sum(0)).ravel().astype(np.float64)
        ds.append(np.log1p(TS*P/P.sum())-b); ps.append(g)
    return np.asarray(ps), np.vstack(ds), b, genes

donor = Donor.load("data/donors/k562_gwps")
m_d, s_d = donor.decompose()
didx = {p:i for i,p in enumerate(donor.perts)}

def align(genes):
    gi = {g:i for i,g in enumerate(genes)}
    tc, dc = [], []
    for j,g in enumerate(donor.genes):
        t = gi.get(g)
        if t is not None: tc.append(t); dc.append(j)
    return np.asarray(tc), np.asarray(dc)

def cos_mean(pred, true):
    pc = pred - pred.mean(0,keepdims=True); tc = true - true.mean(0,keepdims=True)
    na, nb = np.linalg.norm(pc,axis=1), np.linalg.norm(tc,axis=1)
    ok = (na>0)&(nb>0)
    return float(np.mean(np.sum(pc[ok]*tc[ok],axis=1)/(na[ok]*nb[ok])))

# ---- H1 test target (evaluation set) ----
pT, dT, bT, gT = target_deltas("data/harness/h1-test.h5ad")
tcT, dcT = align(gT)
covT = np.array([p in didx for p in pT]); pT, dT = pT[covT], dT[covT]
D_T = s_d[np.ix_([didx[p] for p in pT], dcT)]        # donor effects, donor-covered genes
T_T = dT[:, tcT]                                      # true H1 effects, same genes
print(f"eval: H1 test, {len(pT)} perts x {len(tcT)} shared genes")
base = cos_mean(D_T, T_T)
print(f"  baseline (w=1, unweighted donor transfer): cos = {base:.4f}\n")

def fit_w(D, T, ridge=1.0):
    num = (D*T).sum(0); den = (D*D).sum(0) + ridge
    return num/den

results = {"baseline": base}

# ---- Test 1: fit on H1 train+val (same line, different perturbations) ----
pA, dA, bA, gA = target_deltas("data/harness/h1-trainval.h5ad")
tcA, dcA = align(gA)
covA = np.array([p in didx for p in pA]); pA, dA = pA[covA], dA[covA]
# use the intersection of gene columns with the eval set
common = np.intersect1d(gA[tcA], gT[tcT])
sel_A = np.array([list(gA[tcA]).index(g) for g in common])
sel_T = np.array([list(gT[tcT]).index(g) for g in common])
D_A = s_d[np.ix_([didx[p] for p in pA], dcA)][:, sel_A]
T_A = dA[:, tcA][:, sel_A]
w1 = fit_w(D_A, T_A)
c1 = cos_mean(D_T[:, sel_T]*w1, T_T[:, sel_T])
c1b = cos_mean(D_T[:, sel_T], T_T[:, sel_T])
print(f"Test 1  w fitted on H1 train+val ({len(pA)} perts), eval on H1 test:")
print(f"  unweighted on same genes {c1b:.4f}  ->  weighted {c1:.4f}   ({c1-c1b:+.4f})")
results["h1_trainval_fit"] = {"unweighted": c1b, "weighted": c1}

# ---- Test 2: fit on K562->RPE1, eval on K562->H1 (UNSEEN LINE) ----
rpe = Donor.load("data/donors/rpe1"); m_r, s_r = rpe.decompose()
shared_p = np.intersect1d(donor.perts, rpe.perts)
shared_g = np.intersect1d(donor.genes, rpe.genes)
gi_d = {g:i for i,g in enumerate(donor.genes)}; gi_r = {g:i for i,g in enumerate(rpe.genes)}
ridx = {p:i for i,p in enumerate(rpe.perts)}
D_R = s_d[np.ix_([didx[p] for p in shared_p], [gi_d[g] for g in shared_g])]
T_R = s_r[np.ix_([ridx[p] for p in shared_p], [gi_r[g] for g in shared_g])]
w2_full = fit_w(D_R, T_R)
gpos = {g:i for i,g in enumerate(shared_g)}
inboth = np.array([g in gpos for g in gT[tcT]])
w2 = np.ones(len(tcT)); w2[inboth] = w2_full[[gpos[g] for g in gT[tcT][inboth]]]
c2 = cos_mean(D_T*w2, T_T)
print(f"\nTest 2  w fitted on K562->RPE1 ({len(shared_p)} perts), eval on K562->H1 (UNSEEN LINE):")
print(f"  unweighted {base:.4f}  ->  weighted {c2:.4f}   ({c2-base:+.4f})")
results["rpe1_fit_to_h1"] = {"unweighted": base, "weighted": c2}

# ---- Test 3: shuffled control ----
rng = np.random.default_rng(0)
c3 = cos_mean(D_T*rng.permutation(w2), T_T)
print(f"\nTest 3  same weights, gene labels shuffled (null): {c3:.4f}   ({c3-base:+.4f})")
results["shuffled_null"] = c3

# binary variant: keep only genes with reliably positive transfer
for q in (0.25, 0.5, 0.75):
    thr = np.quantile(w2_full, q)
    mask = np.ones(len(tcT)); mask[inboth] = (w2_full[[gpos[g] for g in gT[tcT][inboth]]] > thr).astype(float)
    print(f"  binary gate, keep genes with w > q{int(q*100)} ({thr:+.3f}): "
          f"cos {cos_mean(D_T*mask, T_T):.4f}")

Path("experiments/004-donor-transfer-h1/gene_weights.json").write_text(json.dumps(results, indent=2))
np.save("data/donors/k562_gwps/transfer_w_rpe1.npy", w2_full.astype(np.float32))
np.save("data/donors/k562_gwps/transfer_w_genes.npy", shared_g)
print("\nsaved transfer weights fitted on K562->RPE1")

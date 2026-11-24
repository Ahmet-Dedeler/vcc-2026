"""KA-2 diagnosis: is the gap to the cd4_compact bank the ESTIMATOR or a bug?
(i)  on-target sign: my pooled delta at the target's OWN gene must be strongly negative (CRISPRi knockdown)
(ii) restrict to well-expressed genes, where log1p(CPM) and log-fold-change nearly coincide -> the cosine
     to the DESeq2 rows should rise if the gap is the functional, not an error."""
import numpy as np, pandas as pd, sys, h5py
sys.path.insert(0, ".")
from lib280 import DORDER, Bank, cosines, residualise
ROOT = "~/Code/virtual-cell-challenge-2026"
BK = f"{ROOT}/data/external/cd4t/cloud-compact/cd4_compact"
bank = Bank("."); SC, CC, CO = bank.scored, bank.cols_cd4, bank.cols_off
o = bank.obs
official = [l.strip() for l in open(f"{ROOT}/data/challenge/gene_names.csv").read().strip().split("\n")[1:]]
offpos = {g: i for i, g in enumerate(official)}
off2map = -np.ones(len(official), dtype=np.int64); off2map[CO] = np.arange(len(CO))

for cond, sub in (("Rest","cd4_rest"), ("Stim48hr","cd4_stim48hr")):
    perts = open(f"{BK}/{sub}/perts.txt").read().split("\n")
    genes = open(f"{BK}/{sub}/genes.txt").read().split("\n")
    delta = np.load(f"{BK}/{sub}/delta.npy", mmap_mode="r")
    base  = np.load(f"{BK}/{sub}/base_mean.npy", mmap_mode="r")
    tg = o[(o.gtype=="targeting") & o.pert.isin(perts)]
    tt = sorted(set(tg.pert))
    b = bank.ctrl_profile(cond, DORDER)[CC]
    M = np.zeros((len(tt), len(CC)))
    for k,t in enumerate(tt): M[k] = bank.profile(bank.counts(tg[tg.pert==t].i.values))[CC] - b
    # (i) on-target
    ont = []
    for k,t in enumerate(tt):
        j = offpos.get(t)
        if j is None: continue
        c = off2map[j]
        if c >= 0: ont.append(M[k, c])
    ont = np.array(ont)
    print(f"{cond}: ON-TARGET own-gene delta over {len(ont)} panel targets: median {np.median(ont):+.4f}, "
          f"frac < 0 = {np.mean(ont<0):.3f}, q90 {np.quantile(ont,.9):+.4f}")
    # (ii) expression-stratified cosine to the bank
    bpos = np.asarray([offpos.get(g,-1) for g in genes])
    take = (bpos>=0) & (off2map[np.clip(bpos,0,None)]>=0)
    mapcols = off2map[bpos[take]]; keep = SC[mapcols]
    mapcols, tidx = mapcols[keep], np.where(take)[0][keep]
    R = residualise(M); R[:,~SC]=0.0
    mine = R[:, mapcols]
    th_raw = np.asarray(delta[[perts.index(t) for t in tt]][:, tidx], dtype=np.float64)
    th = th_raw - th_raw.mean(0, keepdims=True)
    bm = np.asarray(base[[perts.index(t) for t in tt]][:, tidx], dtype=np.float64)
    med_bm = np.nanmedian(bm, axis=0)
    for q in (0.0, 0.5, 0.8, 0.9):
        s = med_bm >= np.nanquantile(med_bm, q)
        cs = cosines(mine[:, s], th[:, s])
        print(f"   genes with median baseMean >= q{q:.2f} ({s.sum():5d} genes): median cosine {np.nanmedian(cs):+.4f}")

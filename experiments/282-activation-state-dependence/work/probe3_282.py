"""282 POST-HOC probe 3: the best-powered deconfounded state contrast, the deconfounded decomposition,
the deconfounded S4 test, and the descriptive tail.

Replicates that share NEITHER the donor NOR the guide:  r1 = (pair P, guide 0),  r2 = (pair Pbar, guide 1).
Everything below is on the 74 targets that carry the same two guides in all 4 donors x 3 states.

    PYTHONDONTWRITEBYTECODE=1 uv run python probe3_282.py
"""
from __future__ import annotations

import csv
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, ".")
from lib282 import (CONDS, DORDER, PAIRINGS, PAIRS, ROOT, SEED, STATE_PAIRS, W280, Bank, boot_median,  # noqa: E402
                    boot_median_diff, boot_spearman, build_index, cosines, sd_index, thin)

P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
OUT: dict = {}
bank = Bank(W280)
SC, CC = bank.scored, bank.cols_cd4
T, IDX = build_index(bank, "common")
TG = [t for t in T if all(len(IDX[(c, d, t)]) == 2 for c in CONDS for d in DORDER)]
rng = np.random.default_rng(SEED)
bud = {t: min(float(IDX[(c, d, t)].total.values[k]) + float(IDX[(c, d2, t)].total.values[k])
              for c in CONDS for d, d2 in PAIRS for k in (0, 1)) for t in TG}


def pair_guide_rows(matched):
    out = {}
    for c in CONDS:
        for p in PAIRS:
            b = bank.ctrl_profile(c, p)[CC]
            for k in (0, 1):
                M = np.zeros((len(TG), len(CC)))
                for j, t in enumerate(TG):
                    S = bank.counts([int(IDX[(c, d, t)].i.values[k]) for d in p])
                    if matched:
                        S = thin(S, bud[t], rng)
                    M[j] = bank.profile(S)[CC] - b
                M -= M.mean(0, keepdims=True)
                M[:, ~SC] = 0.0
                out[("".join(p), k, c)] = M[:, SC]
    return out


# ==================================================================== 1. pooled state contrast
P("=== the pooled DECONFOUNDED state contrast (neither donor nor guide shared on the diff-donor cells) ===")
pool = []
for matched in (False, True):
    Rg = pair_guide_rows(matched)
    same_s, diff_s = [], []
    for a, b_ in PAIRINGS:                       # different donors on the two sides, guide swapped
        ak, bk = "".join(a), "".join(b_)
        for c in CONDS:
            same_s.append(cosines(Rg[(ak, 0, c)], Rg[(bk, 1, c)]))
            same_s.append(cosines(Rg[(ak, 1, c)], Rg[(bk, 0, c)]))
        for si, sj in STATE_PAIRS:
            diff_s.append(cosines(Rg[(ak, 0, si)], Rg[(bk, 1, sj)]))
            diff_s.append(cosines(Rg[(ak, 1, si)], Rg[(bk, 0, sj)]))
    X = np.nanmean(np.vstack(same_s), axis=0)    # SAME state,  diff donors, diff guide
    Y = np.nanmean(np.vstack(diff_s), axis=0)    # DIFF state,  diff donors, diff guide
    d, lo, hi, n = boot_median_diff(X, Y)
    sd = sd_index(Y, X)
    m, slo, shi, sn = boot_median(sd)
    pooled = 1 - np.nanmedian(Y) / np.nanmedian(X)
    tag = "depth-matched" if matched else "full depth"
    P(f"  {tag:<14} same-state {np.nanmedian(X):+.4f}  diff-state {np.nanmedian(Y):+.4f}   "
      f"cost {d:+.4f} [{lo:+.4f},{hi:+.4f}]   SD median {m:+.4f} [{slo:+.4f},{shi:+.4f}]  pooled {pooled:+.4f}")
    pool.append(dict(arm=tag, same_state=float(np.nanmedian(X)), diff_state=float(np.nanmedian(Y)),
                     cost=d, lo=lo, hi=hi, sd_median=m, sd_lo=slo, sd_hi=shi, sd_pooled=float(pooled), n=n))
    if not matched:
        SAME_S, DIFF_S, Rfull = X, Y, Rg
    else:
        del Rg
pd.DataFrame(pool).to_csv("pooled_state282.csv", index=False)
OUT["pooled_state"] = pool

# ==================================================================== 2. deconfounded decomposition
P("\n=== the DECONFOUNDED shared / state-specific / noise decomposition (74 targets) ===")
nT = len(TG)
acc = {k: np.zeros(nT) for k in ("T", "R", "Sh")}
nrep = 0
for a, b_ in PAIRINGS:
    ak, bk = "".join(a), "".join(b_)
    for (g1, g2) in ((0, 1), (1, 0)):
        r1 = np.stack([Rfull[(ak, g1, s)] for s in CONDS])
        r2 = np.stack([Rfull[(bk, g2, s)] for s in CONDS])
        m1, m2 = r1.mean(0), r2.mean(0)
        acc["T"] += 0.5 * ((r1 ** 2).sum(2) + (r2 ** 2).sum(2)).mean(0)
        acc["R"] += (r1 * r2).sum(2).mean(0)
        acc["Sh"] += (m1 * m2).sum(1)
        nrep += 1
        del r1, r2
for k in acc:
    acc[k] /= nrep
acc["Sp"] = acc["R"] - acc["Sh"]
acc["Nz"] = acc["T"] - acc["R"]
tot = acc["T"].sum()
P(f"  pooled (energy-weighted, {nT} targets): shared {acc['Sh'].sum()/tot:.4f}   "
  f"state-specific & reproducible {acc['Sp'].sum()/tot:.4f}   noise {acc['Nz'].sum()/tot:.4f}   "
  f"(reproducible {acc['R'].sum()/tot:.4f});  state share of the reproducible part "
  f"{acc['Sp'].sum()/(acc['Sh'].sum()+acc['Sp'].sum()):.4f}")
OUT["decomp_deconf"] = {k: float(acc[k].sum() / tot) for k in ("Sh", "Sp", "Nz", "R")}
OUT["decomp_deconf"]["state_share_of_reproducible"] = float(
    acc["Sp"].sum() / (acc["Sh"].sum() + acc["Sp"].sum()))
# effective rank of the deconfounded reproducible parts
a, b_ = PAIRINGS[0]
ak, bk = "".join(a), "".join(b_)
r1 = np.stack([Rfull[(ak, 0, s)] for s in CONDS])
r2 = np.stack([Rfull[(bk, 1, s)] for s in CONDS])
m1, m2 = r1.mean(0), r2.mean(0)
D1 = (r1 - m1).reshape(3 * nT, -1)
D2 = (r2 - m2).reshape(3 * nT, -1)
for nm, X1, X2 in (("state-specific", D1, D2), ("shared", m1, m2)):
    G = 0.5 * (X1 @ X2.T + X2 @ X1.T)
    ev = np.linalg.eigvalsh(G)
    pos = ev[ev > 0][::-1]
    pr = float(pos.sum() ** 2 / (pos ** 2).sum())
    P(f"  effective rank ({nm}, {ak}|{bk} guide-swapped): PR {pr:.1f} of {X1.shape[0]} rows; "
      f"top1 {pos[0]/pos.sum():.3f} top5 {pos[:5].sum()/pos.sum():.3f} top20 {pos[:20].sum()/pos.sum():.3f}")
    OUT.setdefault("rank_deconf", {})[nm] = dict(PR=pr, n_rows=int(X1.shape[0]),
                                                 top1=float(pos[0] / pos.sum()),
                                                 top5=float(pos[:5].sum() / pos.sum()),
                                                 top20=float(pos[:20].sum() / pos.sum()))
    del G

# ==================================================================== 3. deconfounded S4
P("\n=== S4 re-run on the DECONFOUNDED index (POST-HOC) ===")
cz = np.load("cnocd4_282.npz", allow_pickle=True)
CONS = pd.DataFrame(dict(target=[str(x) for x in cz["targets"]], C_noCD4=cz["C_noCD4"],
                         n_pairs_noCD4=cz["n_pairs_noCD4"], C_raw=cz["C_raw"])).set_index("target")
c = CONS.reindex(TG)
sdx = sd_index(DIFF_S, SAME_S)
s4 = []
for pred in ("C_noCD4", "C_raw"):
    v = c[pred].values.astype(float)
    if pred == "C_noCD4":
        v = np.where(c["n_pairs_noCD4"].values >= 3, v, np.nan)
    for nm, x in (("SDx (deconfounded)", sdx), ("diff-state cos", DIFF_S), ("same-state cos", SAME_S),
                  ("state cost (same-diff)", SAME_S - DIFF_S)):
        rho, lo, hi, n = boot_spearman(x, v)
        P(f"    rho({nm:<24}, {pred:<8}) {rho:+.3f} [{lo:+.3f},{hi:+.3f}]  n={n}")
        s4.append(dict(pred=pred, stat=nm, rho=rho, lo=lo, hi=hi, n=n))
pd.DataFrame(s4).to_csv("s4_deconf282.csv", index=False)

# ==================================================================== 4. descriptive tail
P("\n=== descriptive tail (POST-HOC, and NOISY: per-target cosines are ~0.01 with ~0.05 spread) ===")
FAM, NAME = {}, {}
with open(f"{ROOT}/data/external/xcell/hgnc_protein_coding.txt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        FAM[r["symbol"]] = ";".join(x.strip() for x in (r["gene_group"] or "").split("|") if x.strip())
        NAME[r["symbol"]] = r["name"]
IMM = ("interleukin", "chemokine", "cytokine", "interferon", "T cell receptor", "tumor necrosis",
       "janus kinase", "signal transducer and activator", "nuclear factor kappa", "mitogen-activated",
       "cyclin", "proliferat", "CD4", "CD8", "CD2", "CD3", "CD28")


def immune(t):
    s = (NAME.get(t, "") + " " + FAM.get(t, "")).lower()
    return any(k.lower() in s for k in IMM)


tail = pd.DataFrame(dict(target=TG, same_state=SAME_S, diff_state=DIFF_S,
                         cost=SAME_S - DIFF_S, SDx=sdx,
                         family=[FAM.get(t, "") for t in TG], name=[NAME.get(t, "") for t in TG],
                         immune=[immune(t) for t in TG]))
tail = tail.merge(c[["C_noCD4", "C_raw"]], left_on="target", right_index=True, how="left")
tail.to_csv("tail282.csv", index=False)
for lab, df in (("MOST state-dependent (largest same-state minus diff-state cosine)",
                 tail.nlargest(20, "cost")),
                ("LEAST state-dependent (smallest / negative cost)", tail.nsmallest(20, "cost"))):
    P(f"  {lab}:")
    P("    " + ", ".join(f"{r.target} {r.cost:+.3f}" for r in df.itertuples()))
    P(f"    immune-programme labelled: {int(df.immune.sum())}/20   (background {tail.immune.mean():.2f})")
OUT["tail"] = dict(n_immune_top20=int(tail.nlargest(20, 'cost').immune.sum()),
                   n_immune_bot20=int(tail.nsmallest(20, 'cost').immune.sum()),
                   background=float(tail.immune.mean()), n=len(tail))
# is the ranking itself reproducible?  split the 6 pairing/guide replicates in half
half1, half2 = [], []
for i, (a, b_) in enumerate(PAIRINGS):
    ak, bk = "".join(a), "".join(b_)
    for j, (g1, g2) in enumerate(((0, 1), (1, 0))):
        acc_s, acc_d = [], []
        for cnd in CONDS:
            acc_s.append(cosines(Rfull[(ak, g1, cnd)], Rfull[(bk, g2, cnd)]))
        for si, sj in STATE_PAIRS:
            acc_d.append(cosines(Rfull[(ak, g1, si)], Rfull[(bk, g2, sj)]))
        cost = np.nanmean(np.vstack(acc_s), 0) - np.nanmean(np.vstack(acc_d), 0)
        (half1 if (2 * i + j) % 2 == 0 else half2).append(cost)
h1 = np.nanmean(np.vstack(half1), 0)
h2 = np.nanmean(np.vstack(half2), 0)
rho, lo, hi, n = boot_spearman(h1, h2)
P(f"  split-half reproducibility of the per-target state cost: Spearman {rho:+.3f} [{lo:+.3f},{hi:+.3f}] n={n}")
OUT["tail"]["split_half_rho"] = dict(rho=rho, lo=lo, hi=hi, n=n)
json.dump(OUT, open("probe3_282.json", "w"), indent=1, default=float)
P(f"[{time.time()-t0:.0f}s]")

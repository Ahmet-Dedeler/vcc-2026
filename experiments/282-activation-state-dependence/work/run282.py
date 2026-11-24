"""282 -- the across-state statistics, the nulls, the SD index, the S4 family/tercile contrast,
the shared-vs-state-specific decomposition, and the post-hoc descriptive tail.

Preregistered in ../PREREGISTRATION.md (sha256 7d75116b…, hashed 2026-09-06T08:03:50Z).
    PYTHONDONTWRITEBYTECODE=1 uv run python run282.py
"""
from __future__ import annotations

import csv
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, ".")
from lib282 import (CONDS, DORDER, PAIRINGS, PAIRS, ROOT, SEED, STATE_PAIRS, W280, Bank,  # noqa: E402
                    boot_median, boot_median_diff, boot_median_diff_unpaired, boot_spearman, build_index,
                    build_rows, cosines,
                    depth_budget, sd_index, stat_S1, stat_S12, stat_S2, stat_W)

P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
OUT: dict = dict(seed=SEED)
bank = Bank(W280)
P(f"bank {bank.X.shape}, scored {bank.scored.sum()}  [{time.time()-t0:.0f}s]")

# ---------------------------------------------------------------- conservation predictors
cz = np.load("cnocd4_282.npz", allow_pickle=True)
CONS = pd.DataFrame(dict(target=[str(x) for x in cz["targets"]], C_noCD4=cz["C_noCD4"],
                         n_pairs_noCD4=cz["n_pairs_noCD4"], C_raw=cz["C_raw"])).set_index("target")
# HGNC families + names (the table 279 used)
FAM, NAME, LOC = {}, {}, {}
with open(f"{ROOT}/data/external/xcell/hgnc_protein_coding.txt") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        FAM[r["symbol"]] = [x.strip() for x in (r["gene_group"] or "").split("|") if x.strip()]
        NAME[r["symbol"]] = r["name"]
        LOC[r["symbol"]] = r["location"]
DE = pd.read_csv(f"{W280}/de_stats_panel280.csv")


def arm(filt: str, matched: bool, tag: str, keep_rows=False):
    """Build the 18 (donor pair, state) residual matrices and read the four cells of the 2x2."""
    targets, IDX = build_index(bank, filt)
    budget = depth_budget(bank, IDX, targets) if matched else None
    R = build_rows(bank, IDX, targets, budget=budget)
    s1, c1 = stat_S1(R, targets)
    s2, c2 = stat_S2(R, targets)
    s12, c12 = stat_S12(R, targets)
    res = dict(tag=tag, filt=filt, matched=matched, targets=targets, S1=s1, S2=s2, S12=s12,
               cells1=c1, cells2=c2, cells12=c12)
    if keep_rows:
        res["R"], res["IDX"], res["budget"] = R, IDX, budget
    else:
        del R
    return res


def row(lab, x):
    m, lo, hi, n = boot_median(x)
    P(f"    {lab:<42} n={n:3d}  median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  mean {np.nanmean(x):+.4f}"
      f"  <=0 {np.nanmean(x <= 0):.3f}")
    return dict(stat=lab, n=n, median=m, lo=lo, hi=hi, mean=float(np.nanmean(x)),
                frac_le0=float(np.nanmean(x <= 0)))


# ==================================================================== the four arms
ARMS = [("common", False, "PRIMARY  state-common guides"),
        ("keepDE", False, "SENS-a   280 per-state keep_for_DE"),
        ("all",    False, "SENS-b   all guides"),
        ("common", True,  "N3       depth-matched (min over the 18 cells)")]
main_rows, res_all = [], {}
for i, (filt, matched, tag) in enumerate(ARMS):
    P(f"\n=== {tag}  [{time.time()-t0:.0f}s] ===")
    A = arm(filt, matched, tag, keep_rows=(i == 0))
    res_all[tag] = A
    P(f"  n targets = {len(A['targets'])}")
    rr = [row("S1  across STATE, same donors", A["S1"]),
          row("S2  across DONORS, same state  [ceiling]", A["S2"]),
          row("S12 across state AND donors", A["S12"])]
    for r in rr:
        r.update(arm_=tag)
    main_rows += rr
    for lab, x, y in (("S2 - S1  (state costs more than donors)", A["S2"], A["S1"]),
                      ("S2 - S12 (value of sharing the STATE)", A["S2"], A["S12"]),
                      ("S1 - S12 (value of sharing the DONORS)", A["S1"], A["S12"])):
        d, lo, hi, n = boot_median_diff(x, y)
        P(f"    {lab:<42} {d:+.4f} [{lo:+.4f},{hi:+.4f}]  n={n}")
        main_rows.append(dict(arm_=tag, stat=lab, n=n, median=d, lo=lo, hi=hi, mean=np.nan, frac_le0=np.nan))
    sd = sd_index(A["S1"], A["S2"])
    A["SD"] = sd
    m, lo, hi, n = boot_median(sd)
    sdp = 1 - np.nanmedian(A["S1"]) / np.nanmedian(A["S2"])
    P(f"    SD = 1 - S1/S2  (S2>0)                     n={n:3d}  median {m:+.4f} [{lo:+.4f},{hi:+.4f}]"
      f"   SD_pooled {sdp:+.4f}   SD>0 {np.nanmean(sd > 0):.3f}")
    main_rows.append(dict(arm_=tag, stat="SD = 1 - S1/S2 (S2>0)", n=n, median=m, lo=lo, hi=hi,
                          mean=float(np.nanmean(sd[np.isfinite(sd)])), frac_le0=float(np.nanmean(sd <= 0))))
    main_rows.append(dict(arm_=tag, stat="SD_pooled = 1 - med(S1)/med(S2)", n=n, median=float(sdp),
                          lo=np.nan, hi=np.nan, mean=np.nan, frac_le0=np.nan))
    strong = A["S2"] >= 0.05
    if strong.sum() >= 20:
        m2, lo2, hi2, n2 = boot_median(sd[strong])
        P(f"    SD on S2 >= 0.05                           n={n2:3d}  median {m2:+.4f} [{lo2:+.4f},{hi2:+.4f}]")
        main_rows.append(dict(arm_=tag, stat="SD on S2>=0.05", n=n2, median=m2, lo=lo2, hi=hi2,
                              mean=np.nan, frac_le0=np.nan))
pd.DataFrame(main_rows).to_csv("main282.csv", index=False)

A = res_all["PRIMARY  state-common guides"]
T, R = A["targets"], A["R"]

# ==================================================================== per state pair / per state
P("\n=== S1 by state pair, S2 by state (primary arm) ===")
brk = []
for (si, sj) in STATE_PAIRS:
    v = np.nanmean(np.vstack([A["cells1"][("".join(p), f"{si}|{sj}")] for p in PAIRS]), axis=0)
    m, lo, hi, n = boot_median(v)
    P(f"    S1 {si:<9}·{sj:<9} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(v<=0):.3f}")
    brk.append(dict(kind="S1", key=f"{si}|{sj}", n=n, median=m, lo=lo, hi=hi, frac_le0=float(np.nanmean(v <= 0))))
for s in CONDS:
    v = np.nanmean(np.vstack([A["cells2"][(f"{''.join(a)}|{''.join(b)}", s)] for a, b in PAIRINGS]), axis=0)
    m, lo, hi, n = boot_median(v)
    P(f"    S2 {s:<19} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(v<=0):.3f}")
    brk.append(dict(kind="S2", key=s, n=n, median=m, lo=lo, hi=hi, frac_le0=float(np.nanmean(v <= 0))))
# N4 pairing symmetry: S1 per donor pair
P("  N4 pairing symmetry -- S1 per donor pair:")
for p in PAIRS:
    v = np.nanmean(np.vstack([A["cells1"][("".join(p), f"{si}|{sj}")] for si, sj in STATE_PAIRS]), axis=0)
    m, lo, hi, n = boot_median(v)
    P(f"    S1 {'+'.join(p):<19} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]")
    brk.append(dict(kind="S1_by_pair", key="+".join(p), n=n, median=m, lo=lo, hi=hi, frac_le0=np.nan))
pd.DataFrame(brk).to_csv("breakdown282.csv", index=False)

# ==================================================================== nulls
P("\n=== nulls ===")
rng = np.random.default_rng(SEED)
perm = rng.permutation(len(T))
while np.any(perm == np.arange(len(T))):
    perm = rng.permutation(len(T))
n1 = np.nanmean(np.vstack([cosines(R[("".join(p), si)], R[("".join(p), sj)][perm])
                           for p in PAIRS for si, sj in STATE_PAIRS]), axis=0)
P(f"  N1 target-shuffled state pairing: mean {np.nanmean(n1):+.5f}  median {np.nanmedian(n1):+.5f}"
  f"  sd {np.nanstd(n1):.4f}   (real S1 median {np.nanmedian(A['S1']):+.4f})")
W = stat_W(bank, A["IDX"], T)
mW, loW, hiW, nW = boot_median(W)
P(f"  W  counting floor (disjoint count halves, HALF depth): median {mW:+.4f} [{loW:+.4f},{hiW:+.4f}]")
OUT["nulls"] = dict(N1_mean=float(np.nanmean(n1)), N1_median=float(np.nanmedian(n1)),
                    N1_sd=float(np.nanstd(n1)), N1_pass=bool(abs(np.nanmean(n1)) < 0.02),
                    W_median=mW, W_lo=loW, W_hi=hiW)

# ==================================================================== S4
P("\n=== S4: does cross-line conservation predict within-individual state dependence? ===")
cons = CONS.reindex(T)
s4 = []
for pred in ("C_noCD4", "C_raw"):
    v = cons[pred].values.astype(float)
    if pred == "C_noCD4":
        v = np.where(cons["n_pairs_noCD4"].values >= 3, v, np.nan)
    for arm_ in ARMS:
        tag = arm_[2]
        B = res_all[tag]
        if B["targets"] != T:
            vv = CONS.reindex(B["targets"])[pred].values.astype(float)
            if pred == "C_noCD4":
                vv = np.where(CONS.reindex(B["targets"])["n_pairs_noCD4"].values >= 3, vv, np.nan)
        else:
            vv = v
        rho, lo, hi, n = boot_spearman(B["SD"], vv)
        star = "   <-- PRIMARY" if (pred == "C_noCD4" and tag.startswith("PRIMARY")) else ""
        P(f"    rho(SD, {pred:<8}) {tag:<34} {rho:+.3f} [{lo:+.3f},{hi:+.3f}]  n={n}{star}")
        s4.append(dict(pred=pred, arm_=tag, stat="rho(SD,C)", rho=rho, lo=lo, hi=hi, n=n))
    # the same on S1 and S2 separately (interpretation, preregistered as "reported alongside")
    for nm, x in (("S1", A["S1"]), ("S2", A["S2"]), ("S2-S1", A["S2"] - A["S1"])):
        rho, lo, hi, n = boot_spearman(x, v)
        P(f"      rho({nm:<5}, {pred:<8})                                   {rho:+.3f} [{lo:+.3f},{hi:+.3f}]  n={n}")
        s4.append(dict(pred=pred, arm_="PRIMARY", stat=f"rho({nm},C)", rho=rho, lo=lo, hi=hi, n=n))
# terciles of C_noCD4
v = np.where(cons["n_pairs_noCD4"].values >= 3, cons["C_noCD4"].values.astype(float), np.nan)
q = np.nanquantile(v, [1 / 3, 2 / 3])
ter = np.where(np.isfinite(v), np.digitize(v, q), -1)
P("  SD / S1 / S2 by C_noCD4 tercile (low -> high conservation):")
for k in range(3):
    m = ter == k
    P(f"    tercile {k}: n={m.sum():3d}  SD {np.nanmedian(A['SD'][m]):+.4f}   S1 {np.nanmedian(A['S1'][m]):+.4f}"
      f"   S2 {np.nanmedian(A['S2'][m]):+.4f}   S2-S1 {np.nanmedian(A['S2'][m]-A['S1'][m]):+.4f}")
    s4.append(dict(pred="C_noCD4", arm_="PRIMARY", stat=f"tercile{k}", rho=float(np.nanmedian(A["SD"][m])),
                   lo=float(np.nanmedian(A["S1"][m])), hi=float(np.nanmedian(A["S2"][m])), n=int(m.sum())))
pd.DataFrame(s4).to_csv("s4_282.csv", index=False)

# families with >= 8 analysis targets (preregistered)
P("  SD by HGNC family (>= 8 analysis targets; preregistered):")
famrows = []
cnt: dict[str, list[int]] = {}
for k, t in enumerate(T):
    for f in FAM.get(t, []):
        cnt.setdefault(f, []).append(k)
for f, ix in sorted(cnt.items(), key=lambda kv: -len(kv[1])):
    if len(ix) < 8:
        continue
    ix = np.asarray(ix)
    oth = np.setdiff1d(np.arange(len(T)), ix)
    d, lo, hi, n = boot_median_diff_unpaired(A["SD"][ix], A["SD"][oth])
    P(f"    {f:<38} n={len(ix):3d}  SD {np.nanmedian(A['SD'][ix]):+.4f} vs rest {np.nanmedian(A['SD'][oth]):+.4f}"
      f"   Δ {d:+.4f} [{lo:+.4f},{hi:+.4f}]   S1 {np.nanmedian(A['S1'][ix]):+.4f} S2 {np.nanmedian(A['S2'][ix]):+.4f}")
    famrows.append(dict(family=f, n=len(ix), sd=float(np.nanmedian(A["SD"][ix])),
                        sd_rest=float(np.nanmedian(A["SD"][oth])), d=d, lo=lo, hi=hi,
                        s1=float(np.nanmedian(A["S1"][ix])), s2=float(np.nanmedian(A["S2"][ix]))))
pd.DataFrame(famrows).to_csv("families282.csv", index=False)
OUT["S4_primary"] = s4[0]

# ==================================================================== decomposition
P(f"\n=== the shared / state-specific / noise decomposition  [{time.time()-t0:.0f}s] ===")
nT, nG = len(T), R[("D1D2", "Rest")].shape[1]
acc = {k: np.zeros(nT) for k in ("T", "R", "Sh", "Sp", "Nz")}
prs = []
for a, b in PAIRINGS:
    r1 = np.stack([R[("".join(a), s)] for s in CONDS])          # (3, nT, nG)
    r2 = np.stack([R[("".join(b), s)] for s in CONDS])
    m1, m2 = r1.mean(0), r2.mean(0)
    Tt = 0.5 * ((r1 ** 2).sum(2) + (r2 ** 2).sum(2)).mean(0)
    Rr = (r1 * r2).sum(2).mean(0)
    Sh = (m1 * m2).sum(1)
    acc["T"] += Tt / 3
    acc["R"] += Rr / 3
    acc["Sh"] += Sh / 3
    # effective ranks, cross-replicate Gram in the (target, state) space
    D1 = (r1 - m1).reshape(3 * nT, nG)
    D2 = (r2 - m2).reshape(3 * nT, nG)
    G = 0.5 * (D1 @ D2.T + D2 @ D1.T)
    ev = np.linalg.eigvalsh(G)
    pos = ev[ev > 0][::-1]
    pr_sp = float(pos.sum() ** 2 / (pos ** 2).sum())
    Gs = 0.5 * (m1 @ m2.T + m2 @ m1.T)
    evs = np.linalg.eigvalsh(Gs)
    poss = evs[evs > 0][::-1]
    pr_sh = float(poss.sum() ** 2 / (poss ** 2).sum())
    prs.append(dict(pairing=f"{''.join(a)}|{''.join(b)}", pr_state_specific=pr_sp, pr_shared=pr_sh,
                    top1_sp=float(pos[0] / pos.sum()), top5_sp=float(pos[:5].sum() / pos.sum()),
                    top20_sp=float(pos[:20].sum() / pos.sum()),
                    top1_sh=float(poss[0] / poss.sum()), top5_sh=float(poss[:5].sum() / poss.sum()),
                    top20_sh=float(poss[:20].sum() / poss.sum()),
                    n_pos_sp=int(len(pos)), n_pos_sh=int(len(poss)), n_rows_sp=3 * nT, n_rows_sh=nT,
                    trace_sp=float(pos.sum()), trace_sh=float(poss.sum())))
    del r1, r2, D1, D2, G, Gs
acc["Sp"] = acc["R"] - acc["Sh"]
acc["Nz"] = acc["T"] - acc["R"]
tot = acc["T"].sum()
P(f"  pooled (energy-weighted over {nT} targets):  shared {acc['Sh'].sum()/tot:.4f}   "
  f"state-specific & reproducible {acc['Sp'].sum()/tot:.4f}   noise {acc['Nz'].sum()/tot:.4f}"
  f"   (reproducible total {acc['R'].sum()/tot:.4f})")
fr = {k: acc[k] / acc["T"] for k in ("Sh", "Sp", "Nz", "R")}
for k, lab in (("Sh", "shared across states"), ("Sp", "state-specific & reproducible"), ("Nz", "noise")):
    m, lo, hi, n = boot_median(fr[k])
    P(f"    per-target fraction {lab:<32} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]")
P("  effective rank (participation ratio) of the reproducible part:")
for r in prs:
    P(f"    {r['pairing']}: state-specific PR {r['pr_state_specific']:.1f} of {r['n_rows_sp']} rows"
      f"  (top1 {r['top1_sp']:.3f}, top5 {r['top5_sp']:.3f}, top20 {r['top20_sp']:.3f});"
      f"  shared PR {r['pr_shared']:.1f} of {r['n_rows_sh']}"
      f"  (top1 {r['top1_sh']:.3f}, top5 {r['top5_sh']:.3f}, top20 {r['top20_sh']:.3f})")
pd.DataFrame(prs).to_csv("rank282.csv", index=False)
OUT["decomp"] = dict(pooled_shared=float(acc["Sh"].sum() / tot), pooled_state=float(acc["Sp"].sum() / tot),
                     pooled_noise=float(acc["Nz"].sum() / tot), pooled_repro=float(acc["R"].sum() / tot),
                     n_targets=nT)

# decomposition by tercile and by family
P("  decomposition by C_noCD4 tercile (of the REPRODUCIBLE part: shared vs state-specific):")
dec_rows = []
for k in range(3):
    m = ter == k
    sh, sp = acc["Sh"][m].sum(), acc["Sp"][m].sum()
    P(f"    tercile {k}: n={m.sum():3d}  shared/T {sh/acc['T'][m].sum():.4f}  state/T {sp/acc['T'][m].sum():.4f}"
      f"  state/(shared+state) {sp/(sh+sp):.4f}")
    dec_rows.append(dict(group=f"C_noCD4 tercile {k}", n=int(m.sum()), shared_T=float(sh / acc["T"][m].sum()),
                         state_T=float(sp / acc["T"][m].sum()), state_share=float(sp / (sh + sp))))
for f, ix in sorted(cnt.items(), key=lambda kv: -len(kv[1])):
    if len(ix) < 8:
        continue
    ix = np.asarray(ix)
    sh, sp = acc["Sh"][ix].sum(), acc["Sp"][ix].sum()
    P(f"    {f:<38} n={len(ix):3d}  shared/T {sh/acc['T'][ix].sum():.4f}  state/T {sp/acc['T'][ix].sum():.4f}"
      f"  state/(shared+state) {sp/(sh+sp):.4f}")
    dec_rows.append(dict(group=f, n=len(ix), shared_T=float(sh / acc["T"][ix].sum()),
                         state_T=float(sp / acc["T"][ix].sum()), state_share=float(sp / (sh + sp))))
sh_all, sp_all = acc["Sh"].sum(), acc["Sp"].sum()
dec_rows.append(dict(group="ALL", n=nT, shared_T=float(sh_all / tot), state_T=float(sp_all / tot),
                     state_share=float(sp_all / (sh_all + sp_all))))
pd.DataFrame(dec_rows).to_csv("decomp282.csv", index=False)

# ==================================================================== the per-target table
de48 = DE[DE.cond == "Stim48hr"].set_index("target")
tab = pd.DataFrame(dict(target=T, S1=A["S1"], S2=A["S2"], S12=A["S12"], W=W, SD=A["SD"],
                        S2_minus_S1=A["S2"] - A["S1"]))
for (si, sj) in STATE_PAIRS:
    tab[f"S1_{si}_{sj}"] = np.nanmean(np.vstack([A["cells1"][("".join(p), f"{si}|{sj}")] for p in PAIRS]), axis=0)
for s in CONDS:
    tab[f"S2_{s}"] = np.nanmean(np.vstack([A["cells2"][(f"{''.join(a)}|{''.join(b)}", s)] for a, b in PAIRINGS]),
                                axis=0)
for k, lab in (("T", "energy_total"), ("R", "energy_repro"), ("Sh", "energy_shared"),
               ("Sp", "energy_state"), ("Nz", "energy_noise")):
    tab[lab] = acc[k]
tab["frac_shared"] = acc["Sh"] / acc["T"]
tab["frac_state"] = acc["Sp"] / acc["T"]
tab["frac_noise"] = acc["Nz"] / acc["T"]
tab["C_noCD4"] = cons["C_noCD4"].values
tab["n_pairs_noCD4"] = cons["n_pairs_noCD4"].values
tab["C_raw"] = cons["C_raw"].values
tab["C_tercile"] = ter
tab["hgnc_families"] = [";".join(FAM.get(t, [])) for t in T]
tab["hgnc_name"] = [NAME.get(t, "") for t in T]
tab["n_de_genes_Stim48hr"] = de48.reindex(T)["n_total_de_genes"].values
tab["n_cells_Stim48hr"] = de48.reindex(T)["n_cells_target"].values
for tag in ("SENS-a   280 per-state keep_for_DE", "SENS-b   all guides",
            "N3       depth-matched (min over the 18 cells)"):
    B = res_all[tag]
    idx = {t: i for i, t in enumerate(B["targets"])}
    short = tag.split()[0]
    for nm in ("S1", "S2", "SD"):
        tab[f"{nm}_{short}"] = [B[nm][idx[t]] if t in idx else np.nan for t in T]
tab.to_csv("state282.csv", index=False)
P(f"\n  wrote state282.csv  ({len(tab)} targets, {len(tab.columns)} columns)")

OUT["elapsed_s"] = time.time() - t0
json.dump(OUT, open("run282.json", "w"), indent=1, default=float)
P(f"[{time.time()-t0:.0f}s]")

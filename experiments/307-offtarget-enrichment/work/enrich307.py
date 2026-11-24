"""307 step 4 -- THE PRIMARY STATISTIC and its joint selection-matched null.

Implements PREREGISTRATION.md sections 5-8 exactly.  Reads only 305's saved integer sums and 301's
saved sequence products (via `seq307.py` / `pool307.py` outputs).  No h5ad is opened.

    PYTHONDONTWRITEBYTECODE=1 uv run python enrich307.py
"""
from __future__ import annotations

import collections
import csv
import json
import pickle
import time
from pathlib import Path

import numpy as np
from scipy import stats

import lib307 as X

B_NULL = 100_000
B_BOOT = 20_000
N_DERANGE = 200
N_CALIB = 200
M_POOL = 200
MIN_CELL = 50          # below this, use the whole hard-stratum cell
MIN_WIDEN = 20         # below this, widen to adjacent basal deciles
LS = (12, 11, 10, 9)
FLOOR_CP50K = 0.5      # 305's measurement floor, on RAW CP50k

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


rng = np.random.default_rng(X.SEED307)

# ------------------------------------------------------------------ inputs
SP = {s: X.Split(s) for s in X.SPLITS}
GENES = SP["train"].genes
GIX = SP["train"].gix
NG = GENES.size
for s in X.SPLITS:
    assert np.array_equal(SP[s].genes, GENES)

cons_rows = X.load_constructs()
gene_meta = X.load_genes()
D = pickle.load(open(X.HERE / "poolhits307.pkl", "rb"))
pi_off = D["pi_off"]
lookup = pickle.load(open(X.HERE / "seq307_lookup.pkl", "rb"))
LIB_TARGETS = set(lookup["lib_targets"])
CAND_L = lookup["cand_L"]          # construct -> {L: {gene: [offsets]}}
P(f"{len(cons_rows)} constructs, {NG} genes, {len(LIB_TARGETS)} library-targeted genes excluded")

# covariate matrix over all measured genes (0 where the gene has no promoter window)
COV_NAMES = ["gc", "log1p_tss_density", "map_mult", "lowcomp", "mean_offset_100",
             "log1p_n_windows", "log1p_n_sites"]
COV = np.zeros((NG, len(COV_NAMES)))
PI = {L: np.zeros(NG) for L in LS}
HAS_WIN = np.zeros(NG, dtype=bool)
for g, r in gene_meta.items():
    i = GIX[g]
    HAS_WIN[i] = True
    COV[i] = [float(r["gc"]), np.log1p(float(r["tss_density_10kb"])), float(r["map_mult"]),
              float(r["lowcomp"]), float(r["mean_offset"]) / 100.0,
              np.log1p(float(r["n_windows"])), np.log1p(float(r["n_sites"]))]
for L in LS:
    for g, v in pi_off[L].items():
        PI[L][GIX[g]] = v
IS_LIB = np.zeros(NG, dtype=bool)
for g in LIB_TARGETS:
    IS_LIB[GIX[g]] = True

# ---------------------------------------------------- base eligible universe per seed length
BASE = {}
PIQ = {}
ZCOV = {}
for L in LS:
    ix = np.flatnonzero(HAS_WIN & (PI[L] > 0) & ~IS_LIB)
    BASE[L] = ix
    q = np.quantile(PI[L][ix], [0.2, 0.4, 0.6, 0.8])
    PIQ[L] = np.searchsorted(q, PI[L], side="right")            # 0..4 over all genes
    mu, sd = COV[ix].mean(0), COV[ix].std(0)
    sd[sd == 0] = 1.0
    ZCOV[L] = (COV - mu) / sd
    P(f"  eligible universe at L={L}: {ix.size:,} genes")

# basal decile bins per (split, L)
BASDEC = {}
for s in X.SPLITS:
    raw = 50000.0 * SP[s].ctrl_pooled / SP[s].ctrl_pooled.sum()
    for L in LS:
        q = np.quantile(SP[s].basal[BASE[L]], np.arange(1, 10) / 10.0)
        BASDEC[(s, L)] = np.searchsorted(q, SP[s].basal, side="right")   # 0..9
    SP[s].raw_cp50k = raw

# ------------------------------------------------------------------ per-construct objects
EFF = {}
SPLIT_OF = {}
for r in cons_rows:
    c = r["construct"]
    if r["unit"] == "targeting":
        s = r["split"]
    else:
        s = "train"                              # NTC sums are identical across splits (KA-3)
    SPLIT_OF[c] = s
    EFF[c] = SP[s].effect(c, loo_control=(r["unit"] != "targeting"))
P(f"effect vectors built for {len(EFF)} constructs, peak {X.peak_gib():.2f} GiB")


def exclusions(c: str, row: dict) -> np.ndarray:
    """Genes ineligible as pseudo-candidates for construct c (PREREGISTRATION sec 5.2)."""
    bad = set()
    tgt = row["intended_target"]
    if tgt:
        bad.add(tgt)
    for L in LS:
        bad |= set(CAND_L[c][L])
    return np.asarray(sorted(GIX[g] for g in bad if g in GIX), dtype=np.int64)


def build_pools(c: str, row: dict):
    """One matched pool per candidate.  Returns list of dicts."""
    L = int(row["L_star"]) if row["L_star"] else 0
    if L == 0:
        return []
    s = SPLIT_OF[c]
    ix = BASE[L]
    cell = BASDEC[(s, L)] * 5 + PIQ[L]
    bad = set(exclusions(c, row).tolist())
    out = []
    for g in row["candidates"].split(","):
        gi = GIX[g]
        widened = 0
        members = ix[(cell[ix] == cell[gi])]
        members = np.asarray([m for m in members if m not in bad])
        if members.size < MIN_WIDEN:                        # widen to adjacent basal deciles
            dec, piq = cell[gi] // 5, cell[gi] % 5
            keep = (PIQ[L][ix] == piq) & (np.abs(BASDEC[(s, L)][ix] - dec) <= 1)
            members = np.asarray([m for m in ix[keep] if m not in bad])
            widened = 1
        if members.size > MIN_CELL:
            d = np.linalg.norm(ZCOV[L][members] - ZCOV[L][gi], axis=1)
            members = members[np.argsort(d, kind="stable")[:M_POOL]]
        out.append(dict(gene=g, gi=gi, L=L, members=members, widened=widened))
    return out


def percentiles(e: np.ndarray, members: np.ndarray):
    """Mid-rank percentile of every pool member within the pool (lower = more repressed)."""
    v = e[members]
    order = np.argsort(v, kind="stable")
    sv = v[order]
    lo = np.searchsorted(sv, v, side="left")
    hi = np.searchsorted(sv, v, side="right")
    return (lo + 0.5 * (hi - lo)) / v.size


def statistic(c: str, pools: list, e: np.ndarray | None = None):
    """Observed z_c and, for the null, the per-slot percentile vectors."""
    e = EFF[c] if e is None else e
    if not pools:
        return 0.0, 0.0, [], 0
    qs, qpools = [], []
    for p in pools:
        qp = percentiles(e, p["members"])
        qpools.append(qp)
        v = e[p["gi"]]
        pv = e[p["members"]]
        qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
    k = len(qs)
    r = float(np.mean(qs))
    return float(np.sqrt(12 * k) * (0.5 - r)), r, qpools, k


# ------------------------------------------------------------------ observed
POOLS, ZOBS, ROBS, QPOOL, KC = {}, {}, {}, {}, {}
for r in cons_rows:
    c = r["construct"]
    POOLS[c] = build_pools(c, r)
    ZOBS[c], ROBS[c], QPOOL[c], KC[c] = statistic(c, POOLS[c])

TARG = [r["construct"] for r in cons_rows if r["unit"] == "targeting"]
NTCC = [r["construct"] for r in cons_rows if r["unit"] == "non-targeting"]
S_obs = float(np.mean([ZOBS[c] for c in TARG]))
P(f"\nPRIMARY  S_obs = {S_obs:+.5f}   (n = {len(TARG)} targeting constructs, "
  f"{sum(KC[c] for c in TARG)} candidate genes)")


# ------------------------------------------------------------------ the joint null
def joint_null(constructs, B, rng):
    """Redraw every construct's whole candidate set from its matched pools; return (B, n) z."""
    Z = np.zeros((B, len(constructs)), dtype=np.float64)
    dup_iters = 0
    for j, c in enumerate(constructs):
        k = KC[c]
        if k == 0:
            continue
        qp = QPOOL[c]
        mem = [POOLS[c][i]["members"] for i in range(k)]
        idx = np.stack([rng.integers(0, qp[i].size, B) for i in range(k)], axis=1)
        if k > 1:                                   # without replacement across the k slots
            for _ in range(50):
                gid = np.stack([mem[i][idx[:, i]] for i in range(k)], axis=1)
                gs = np.sort(gid, axis=1)
                bad = np.flatnonzero((gs[:, 1:] == gs[:, :-1]).any(axis=1))
                if bad.size == 0:
                    break
                dup_iters += 1
                for i in range(k):
                    idx[bad, i] = rng.integers(0, qp[i].size, bad.size)
        q = np.stack([qp[i][idx[:, i]] for i in range(k)], axis=1)
        Z[:, j] = np.sqrt(12 * k) * (0.5 - q.mean(axis=1))
    return Z, dup_iters


Zn, dup = joint_null(TARG, B_NULL, rng)
Sn = Zn.mean(axis=1)
p_joint = (1 + int(np.sum(Sn >= S_obs))) / (B_NULL + 1)
P(f"joint null: B = {B_NULL:,} draws, floor p = {1 / (B_NULL + 1):.2e}; "
  f"null S mean {Sn.mean():+.5f} sd {Sn.std():.5f}, "
  f"95% [{np.quantile(Sn, 0.025):+.5f}, {np.quantile(Sn, 0.975):+.5f}]")
P(f"PRIMARY p_joint (one-sided, repression) = {p_joint:.6g}   "
  f"[{'<=' if p_joint <= 0.05 else '>'} alpha = 0.05]")

# leave-one-construct-out
tot = Zn.sum(axis=1)
tot_obs = float(np.sum([ZOBS[c] for c in TARG]))
loo = {}
for j, c in enumerate(TARG):
    s_obs_j = (tot_obs - ZOBS[c]) / (len(TARG) - 1)
    s_null_j = (tot - Zn[:, j]) / (len(TARG) - 1)
    loo[c] = dict(S=float(s_obs_j),
                  p=(1 + int(np.sum(s_null_j >= s_obs_j))) / (B_NULL + 1))
loo_max = max(v["p"] for v in loo.values())
loo_worst = max(loo, key=lambda c: loo[c]["p"])
P(f"leave-one-construct-out: max p = {loo_max:.6g} (deleting {loo_worst}); "
  f"min p = {min(v['p'] for v in loo.values()):.6g}")

VERDICT = "PASS" if (p_joint <= 0.05 and loo_max <= 0.05) else "FAIL"
P(f"\n>>> C8 VERDICT: {VERDICT}  (rule: p_joint <= 0.05 AND max leave-one-out p <= 0.05)")

# construct-level bootstrap
zt = np.asarray([ZOBS[c] for c in TARG])
bs = zt[rng.integers(0, zt.size, (B_BOOT, zt.size))].mean(axis=1)
ci = (float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)))
P(f"construct-level bootstrap (n = {len(TARG)} targeting constructs, {B_BOOT:,} resamples): "
  f"S = {S_obs:+.5f}, 95% CI [{ci[0]:+.5f}, {ci[1]:+.5f}]")

# interpretable effect size, log1p-CP50k units
dlt = []
for c in TARG:
    if KC[c] == 0:
        dlt.append(0.0)
        continue
    e = EFF[c]
    dlt.append(float(np.mean([e[p["gi"]] - np.median(e[p["members"]]) for p in POOLS[c]])))
dlt = np.asarray(dlt)
dbs = dlt[rng.integers(0, dlt.size, (B_BOOT, dlt.size))].mean(axis=1)
dci = (float(np.quantile(dbs, 0.025)), float(np.quantile(dbs, 0.975)))
P(f"effect size Delta = {dlt.mean():+.6f} log1p-CP50k per targeting construct, "
  f"95% CI [{dci[0]:+.6f}, {dci[1]:+.6f}]  (negative = candidates repressed)")

# ------------------------------------------------------------------ KA-4  no-signal row
ka4 = []
by_split = collections.defaultdict(list)
for c in TARG:
    by_split[SPLIT_OF[c]].append(c)
rowmap = {r["construct"]: r for r in cons_rows}
for _ in range(N_DERANGE):
    perm = {}
    for s, cs in by_split.items():
        while True:
            q = list(rng.permutation(cs))
            if all(a != b for a, b in zip(cs, q)):
                break
        perm.update(dict(zip(cs, q)))
    zs = []
    for c in TARG:
        d = perm[c]                      # score construct c's expression on d's candidate set
        if KC[d] == 0:
            zs.append(0.0)
            continue
        bad = set(exclusions(c, rowmap[c]).tolist())
        e, qs = EFF[c], []
        for p in POOLS[d]:
            if p["gi"] in bad:
                continue
            mem = p["members"][~np.isin(p["members"], list(bad))] if bad else p["members"]
            if mem.size < 20:
                continue
            v, pv = e[p["gi"]], e[mem]
            qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
        zs.append(float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))) if qs else 0.0)
    ka4.append(float(np.mean(zs)))
ka4 = np.asarray(ka4)
ka4_ok = bool(abs(ka4.mean()) <= 0.05 and np.quantile(ka4, 0.025) <= 0 <= np.quantile(ka4, 0.975))
P(f"\nKA-4 no-signal row ({N_DERANGE} within-split derangements, construct scored against another "
  f"construct's candidate set): mean S {ka4.mean():+.5f}, "
  f"95% [{np.quantile(ka4, 0.025):+.5f}, {np.quantile(ka4, 0.975):+.5f}] "
  f"-> {'PASS' if ka4_ok else 'FAIL'}")

# ------------------------------------------------------------------ KA-5  null calibration
pick = rng.choice(B_NULL, N_CALIB, replace=False)
pcal = np.asarray([(1 + int(np.sum(np.delete(Sn, i) >= Sn[i]))) / B_NULL for i in pick])
ks = stats.kstest(pcal, "uniform")
frac05 = float(np.mean(pcal <= 0.05))
ka5_ok = bool(ks.pvalue >= 0.05 and 0.02 <= frac05 <= 0.10)
P(f"KA-5 null calibration ({N_CALIB} null draws each p-valued against the rest): "
  f"KS uniformity p = {ks.pvalue:.4f}, fraction p<=0.05 = {frac05:.4f} "
  f"-> {'PASS' if ka5_ok else 'FAIL'}")

# ------------------------------------------------------------------ KA-6  matching balance
bal = {}
for name, vec in [("pi_off_Lstar", None), ("basal_log1p_cp50k", None)] + \
                 [(n, COV[:, i]) for i, n in enumerate(COV_NAMES)]:
    a, b = [], []
    for c in TARG:
        for p in POOLS[c]:
            if vec is None and name == "pi_off_Lstar":
                v, w = np.log10(PI[p["L"]][p["gi"]]), np.log10(PI[p["L"]][p["members"]])
            elif vec is None:
                s = SPLIT_OF[c]
                v, w = SP[s].basal[p["gi"]], SP[s].basal[p["members"]]
            else:
                v, w = vec[p["gi"]], vec[p["members"]]
            a.append(float(v))
            b.append(np.asarray(w, dtype=float))
    a = np.asarray(a)
    bm = np.asarray([x.mean() for x in b])
    sd = np.sqrt((a.var() + np.mean([x.var() for x in b])) / 2)
    bal[name] = dict(cand_mean=float(a.mean()), pool_mean=float(bm.mean()),
                     smd=float((a.mean() - bm.mean()) / sd) if sd > 0 else 0.0)
ka6_ok = all(abs(v["smd"]) <= 0.25 for v in bal.values())
P("KA-6 matching balance (real candidates vs their matched pools, standardized mean difference):")
for n, v in bal.items():
    P(f"     {n:22s} candidate {v['cand_mean']:+9.4f}  pool {v['pool_mean']:+9.4f}  "
      f"SMD {v['smd']:+.4f}")
P(f"     -> {'PASS' if ka6_ok else 'FAIL'} (rule |SMD| <= 0.25 on every axis)")

# ------------------------------------------------------------------ secondaries
sec: dict = {}

# S1 failing vs working
zf = np.asarray([ZOBS[c] for c in TARG if rowmap[c]["group"] == "failing"])
zw = np.asarray([ZOBS[c] for c in TARG if rowmap[c]["group"] == "working"])
db = np.asarray([zf[rng.integers(0, zf.size, zf.size)].mean()
                 - zw[rng.integers(0, zw.size, zw.size)].mean() for _ in range(B_BOOT)])
mw = stats.mannwhitneyu(zf, zw, alternative="two-sided")
sec["S1_failing_vs_working"] = dict(
    n_failing=int(zf.size), n_working=int(zw.size),
    mean_failing=float(zf.mean()), mean_working=float(zw.mean()),
    diff=float(zf.mean() - zw.mean()),
    ci95=[float(np.quantile(db, 0.025)), float(np.quantile(db, 0.975))],
    mannwhitney_p=float(mw.pvalue))
P(f"\nS1 failing ({zf.size}) vs working ({zw.size}) constructs: z {zf.mean():+.4f} vs "
  f"{zw.mean():+.4f}, diff {zf.mean() - zw.mean():+.4f} "
  f"[{np.quantile(db, 0.025):+.4f}, {np.quantile(db, 0.975):+.4f}], MWU p {mw.pvalue:.4f}")

# S2 the 31 non-targeting constructs
S_ntc = float(np.mean([ZOBS[c] for c in NTCC]))
Zn2, _ = joint_null(NTCC, B_NULL, np.random.default_rng(X.SEED307 + 1))
Sn2 = Zn2.mean(axis=1)
p_ntc = (1 + int(np.sum(Sn2 >= S_ntc))) / (B_NULL + 1)
zn_ = np.asarray([ZOBS[c] for c in NTCC])
bn = zn_[rng.integers(0, zn_.size, (B_BOOT, zn_.size))].mean(axis=1)
sec["S2_non_targeting"] = dict(n=len(NTCC), S=S_ntc, p_joint=p_ntc,
                               ci95=[float(np.quantile(bn, 0.025)), float(np.quantile(bn, 0.975))])
P(f"S2 the 31 non-targeting constructs (candidate-bearing observations, NOT negative controls): "
  f"S = {S_ntc:+.5f} [{np.quantile(bn, 0.025):+.5f}, {np.quantile(bn, 0.975):+.5f}], "
  f"joint-null p = {p_ntc:.6g}")

# S3 GRB10_P2 -> SNX8 and its partner
s3 = {}
for c in ("GRB10_P2", "GRB10_P1"):
    e = EFF[c]
    s3[c] = dict(L_star=int(rowmap[c]["L_star"]), z=ZOBS[c], k=KC[c],
                 candidates={p["gene"]: dict(
                     e=float(e[p["gi"]]),
                     q=float((np.sum(e[p["members"]] < e[p["gi"]])
                              + 0.5 * np.sum(e[p["members"]] == e[p["gi"]])) / p["members"].size),
                     pool=int(p["members"].size)) for p in POOLS[c]})
    for g in ("SNX8", "NEUROG3"):
        s3[c][f"e_{g}"] = float(e[GIX[g]])
sec["S3_GRB10"] = s3
P(f"S3 GRB10_P2 -> SNX8: e = {s3['GRB10_P2']['candidates']['SNX8']['e']:+.6f}, matched-pool "
  f"percentile {s3['GRB10_P2']['candidates']['SNX8']['q']:.4f}; partner GRB10_P1 e(SNX8) = "
  f"{s3['GRB10_P1']['e_SNX8']:+.6f}")

# S4 measurement floor
below = []
for c in TARG:
    s = SPLIT_OF[c]
    for p in POOLS[c]:
        if SP[s].raw_cp50k[p["gi"]] < FLOOR_CP50K:
            below.append((c, p["gene"], float(SP[s].raw_cp50k[p["gi"]])))
z_above = []
for c in TARG:
    s = SPLIT_OF[c]
    keep = [p for p in POOLS[c] if SP[s].raw_cp50k[p["gi"]] >= FLOOR_CP50K]
    if not keep:
        z_above.append(0.0)
        continue
    e = EFF[c]
    qs = [(np.sum(e[p["members"]] < e[p["gi"]]) + 0.5 * np.sum(e[p["members"]] == e[p["gi"]]))
          / p["members"].size for p in keep]
    z_above.append(float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))))
sec["S4_floor"] = dict(n_candidates=int(sum(KC[c] for c in TARG)), n_below_floor=len(below),
                       below=[[a, b, round(v, 4)] for a, b, v in below],
                       S_above_floor_only=float(np.mean(z_above)))
P(f"S4 measurement floor (RAW CP50k >= {FLOOR_CP50K}): {len(below)} of "
  f"{sum(KC[c] for c in TARG)} candidates below it; S on above-floor candidates only = "
  f"{np.mean(z_above):+.5f}")

# S6 level and population restrictions
sub = {}
for lab, cs in [("L*=12 constructs", [c for c in TARG if int(rowmap[c]["L_star"]) == 12]),
                ("L*=11 constructs", [c for c in TARG if int(rowmap[c]["L_star"]) == 11]),
                ("L*=10 constructs", [c for c in TARG if int(rowmap[c]["L_star"]) == 10]),
                ("L*=9 constructs", [c for c in TARG if int(rowmap[c]["L_star"]) == 9])]:
    if not cs:
        continue
    s_ = float(np.mean([ZOBS[c] for c in cs]))
    sn_ = Zn[:, [TARG.index(c) for c in cs]].mean(axis=1)
    sub[lab] = dict(n=len(cs), S=s_, p_joint=(1 + int(np.sum(sn_ >= s_))) / (B_NULL + 1),
                    constructs=cs)
    P(f"S6 {lab:18s} n={len(cs):2d}  S = {s_:+.5f}  joint-null p = {sub[lab]['p_joint']:.6g}")
sec["S6_by_level"] = sub

# S7 promiscuity of the 2025 library vs the pool
pm = json.load(open(X.HERE / "pool307_meta.json"))["pseudo_construct"]
obs_l = collections.Counter(int(rowmap[c]["L_star"]) for c in TARG)
exp_l = np.asarray([pm[f"L*={L}"]["frac"] for L in LS])
exp_l = exp_l / exp_l.sum() * len(TARG)
chi = stats.chisquare([obs_l[L] for L in LS], exp_l)
sec["S7_promiscuity"] = dict(observed={str(L): obs_l[L] for L in LS},
                             expected={str(L): float(e) for L, e in zip(LS, exp_l)},
                             chi2=float(chi.statistic), p=float(chi.pvalue))
P(f"S7 L* of the 34 real constructs {[obs_l[L] for L in LS]} vs the hCRISPRi-v2 pseudo-construct "
  f"expectation {[round(float(e), 1) for e in exp_l]}: chi2 = {chi.statistic:.2f}, p = {chi.pvalue:.3f}")

# ------------------------------------------------------------------ per-construct CSV
with open(X.HERE / "enrichment307.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["construct", "unit", "group", "intended_target", "split", "kd_lfc", "cells",
                "L_star", "k_candidates", "candidates", "pool_size_mean", "pools_widened",
                "null_draws", "r_c", "z_c", "delta_log1p_cp50k",
                "null_z_mean", "null_z_sd", "null_z_p2.5", "null_z_p97.5", "p_construct",
                "n_below_floor", "min_cand_cp50k"])
    for r in cons_rows:
        c = r["construct"]
        s = SPLIT_OF[c]
        j = TARG.index(c) if c in TARG else None
        zn2 = Zn[:, j] if j is not None else Zn2[:, NTCC.index(c)]
        e = EFF[c]
        d = float(np.mean([e[p["gi"]] - np.median(e[p["members"]]) for p in POOLS[c]])) if POOLS[c] else 0.0
        cp = [float(SP[s].raw_cp50k[p["gi"]]) for p in POOLS[c]]
        w.writerow([c, r["unit"], r["group"], r["intended_target"], r["split"], r["kd_lfc"],
                    r["cells_full"], r["L_star"], KC[c], r["candidates"],
                    round(float(np.mean([p["members"].size for p in POOLS[c]])), 1) if POOLS[c] else 0,
                    sum(p["widened"] for p in POOLS[c]), B_NULL,
                    round(ROBS[c], 5), round(ZOBS[c], 5), round(d, 6),
                    round(float(zn2.mean()), 5), round(float(zn2.std()), 5),
                    round(float(np.quantile(zn2, 0.025)), 5), round(float(np.quantile(zn2, 0.975)), 5),
                    round((1 + int(np.sum(zn2 >= ZOBS[c]))) / (B_NULL + 1), 6),
                    sum(1 for v in cp if v < FLOOR_CP50K), round(min(cp), 4) if cp else ""])
P(f"\nwrote enrichment307.csv ({len(cons_rows)} constructs)")

res = dict(
    preregistration_sha256=open(X.HERE / "PREREGISTRATION.sha256").read().split()[0],
    primary=dict(statistic="S = mean over the 34 targeting constructs of sqrt(12 k_c) (0.5 - r_c)",
                 unit="targeting construct", n=len(TARG), S=S_obs, ci95_construct_bootstrap=list(ci),
                 p_joint=p_joint, B=B_NULL, attainable_floor_p=1 / (B_NULL + 1),
                 null_mean=float(Sn.mean()), null_sd=float(Sn.std()),
                 null_ci=[float(np.quantile(Sn, 0.025)), float(np.quantile(Sn, 0.975))],
                 leave_one_out_max_p=loo_max, leave_one_out_worst=loo_worst,
                 leave_one_out={c: v for c, v in loo.items()},
                 verdict=VERDICT,
                 delta_log1p_cp50k=float(dlt.mean()), delta_ci95=list(dci),
                 n_candidate_genes=int(sum(KC[c] for c in TARG)),
                 duplicate_redraw_iterations=int(dup)),
    ka=dict(KA4=dict(mean=float(ka4.mean()),
                     ci=[float(np.quantile(ka4, 0.025)), float(np.quantile(ka4, 0.975))],
                     n=N_DERANGE, passed=ka4_ok),
            KA5=dict(ks_p=float(ks.pvalue), frac_le_05=frac05, n=N_CALIB, passed=ka5_ok),
            KA6=dict(balance=bal, passed=ka6_ok)),
    secondaries=sec,
    per_construct={c: dict(z=ZOBS[c], r=ROBS[c], k=KC[c], L_star=rowmap[c]["L_star"],
                           unit=rowmap[c]["unit"], group=rowmap[c]["group"]) for c in ZOBS},
    peak_rss_gib=round(X.peak_gib(), 3), seconds=round(time.time() - t0, 1))
(X.HERE / "results307.json").write_text(json.dumps(res, indent=1, default=float))
(X.HERE / "enrich307.log").write_text("\n".join(log) + "\n")
P(f"done {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")

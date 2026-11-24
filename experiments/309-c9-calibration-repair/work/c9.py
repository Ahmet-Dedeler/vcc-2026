"""309 step 2 -- THE REPAIRED PRIMARY, its null, the remaining known-answer rows, the secondaries.

Implements PREREGISTRATION.md sections 6-10 exactly.  Reads only 305's saved integer sums and
301's / 307's saved sequence products.  No h5ad is opened.

    PYTHONDONTWRITEBYTECODE=1 uv run python c9.py
"""
from __future__ import annotations

import collections
import csv
import json
import pickle
import time

import numpy as np
from scipy import stats

import lib309 as X
import pools309 as PL

B_NULL = 100_000
B_BOOT = 20_000
N_DERANGE = 5_000
N_CALIB = 200
FLOOR_CP50K = 0.5
LS = X.LS

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


# ------------------------------------------------------------------ inputs
SP = {s: X.Split(s) for s in X.SPLITS}
GIX = SP["train"].gix
cons_rows = X.load_constructs()
rowmap = {r["construct"]: r for r in cons_rows}
ALL = [r["construct"] for r in cons_rows]
TARG = [r["construct"] for r in cons_rows if r["unit"] == "targeting"]
NTCC = [r["construct"] for r in cons_rows if r["unit"] != "targeting"]
gene_meta = X.load_genes()
lookup = pickle.load(open(X.E307 / "seq307_lookup.pkl", "rb"))
pi_off = pickle.load(open(X.E307 / "poolhits307.pkl", "rb"))["pi_off"]
SPLIT_OF = {r["construct"]: (r["split"] if r["unit"] == "targeting" else "train")
            for r in cons_rows}
EFF = {c: SP[SPLIT_OF[c]].effect(c, loo_control=(rowmap[c]["unit"] != "targeting")) for c in ALL}
NG = SP["train"].genes.size
PI = {L: np.zeros(NG) for L in LS}
for L in LS:
    for g, v in pi_off[L].items():
        PI[L][GIX[g]] = v
P(f"{len(ALL)} constructs ({len(TARG)} targeting, {len(NTCC)} NTC), {NG} genes; "
  f"effect vectors built, peak {X.peak_gib():.2f} GiB")

# recipient groups: the other targeting constructs of the same split; the other NTC constructs
BY_SPLIT = collections.defaultdict(list)
for c in TARG:
    BY_SPLIT[SPLIT_OF[c]].append(c)
GROUP = {c: BY_SPLIT[SPLIT_OF[c]] for c in TARG}
for c in NTCC:
    GROUP[c] = NTCC


# ================================================================== cross-recipient machinery
def build(constructs, hard_windows, global_excl, recips=None):
    """Pools + the (recipient x pool) percentile tensors PREREGISTRATION sec 6.1 needs."""
    PL.setup(SP, cons_rows, gene_meta, lookup, pi_off,
             hard_windows=hard_windows, global_excl=global_excl)
    POOLS = {c: PL.build_pools(c, rowmap[c]) for c in constructs}
    R = {c: (recips[c] if recips else [d for d in GROUP[c] if d != c]) for c in constructs}
    obj = {}
    for c in constructs:
        pl = POOLS[c]
        k = len(pl)
        if k == 0:
            obj[c] = dict(k=0, pools=[], recips=R[c])
            continue
        rec = R[c]
        QG = np.zeros((len(rec), k))                       # candidate percentile per recipient
        Qm = [np.zeros((len(rec), p["members"].size)) for p in pl]   # member percentiles
        DG = np.zeros((len(rec), k))                       # candidate - pool median, log1p-CP50k
        for j, d in enumerate(rec):
            e = EFF[d]
            for i, p in enumerate(pl):
                QG[j, i] = PL.q_of(e, p)
                Qm[i][j] = PL.percentiles(e, p["members"])
                DG[j, i] = e[p["gi"]] - np.median(e[p["members"]])
        e = EFF[c]
        qG = np.asarray([PL.q_of(e, p) for p in pl])
        qm = [PL.percentiles(e, p["members"]) for p in pl]
        dG = np.asarray([e[p["gi"]] - np.median(e[p["members"]]) for p in pl])
        obj[c] = dict(k=k, pools=pl, recips=rec, qG=qG, qm=qm, dG=dG,
                      QG=QG, Qm=Qm, DG=DG,
                      qbarG=QG.mean(0), qbarm=[Qm[i].mean(0) for i in range(k)],
                      dbarG=DG.mean(0))
    return POOLS, obj


def z_unc(o) -> float:
    return 0.0 if o["k"] == 0 else float(np.sqrt(12 * o["k"]) * (0.5 - o["qG"].mean()))


def b_lvl(o) -> float:
    """The construct's wrong-recipient level b_c."""
    return 0.0 if o["k"] == 0 else float(np.sqrt(12 * o["k"]) * (0.5 - o["qbarG"].mean()))


def z_rep(o) -> float:
    return 0.0 if o["k"] == 0 else float(np.sqrt(12 * o["k"]) * (o["qbarG"].mean() - o["qG"].mean()))


def z_in(o, j: int) -> float:
    """z of this construct's set scored in its j-th wrong recipient's cells (uncentred)."""
    return float(np.sqrt(12 * o["k"]) * (0.5 - o["QG"][j].mean()))


def joint_null_rep(constructs, obj, B, rng, centred=True):
    """Redraw every construct's whole candidate set from its matched pools; recompute both terms."""
    Z = np.zeros((B, len(constructs)))
    for j, c in enumerate(constructs):
        o = obj[c]
        k = o["k"]
        if k == 0:
            continue
        mem = [p["members"] for p in o["pools"]]
        idx = np.stack([rng.integers(0, mem[i].size, B) for i in range(k)], axis=1)
        if k > 1:                                   # without replacement across the k slots
            for _ in range(50):
                gid = np.stack([mem[i][idx[:, i]] for i in range(k)], axis=1)
                gs = np.sort(gid, axis=1)
                bad = np.flatnonzero((gs[:, 1:] == gs[:, :-1]).any(axis=1))
                if bad.size == 0:
                    break
                for i in range(k):
                    idx[bad, i] = rng.integers(0, mem[i].size, bad.size)
        own = np.stack([o["qm"][i][idx[:, i]] for i in range(k)], axis=1).mean(axis=1)
        if centred:
            bar = np.stack([o["qbarm"][i][idx[:, i]] for i in range(k)], axis=1).mean(axis=1)
            Z[:, j] = np.sqrt(12 * k) * (bar - own)
        else:
            Z[:, j] = np.sqrt(12 * k) * (0.5 - own)
    return Z


def verdict_block(constructs, obj, label, seed, centred=True, B=B_NULL):
    """S, joint null p, leave-one-out, construct bootstrap, Delta -- the primary machinery."""
    z = np.asarray([(z_rep(obj[c]) if centred else z_unc(obj[c])) for c in constructs])
    S = float(z.mean())
    rng = np.random.default_rng(seed)
    Zn = joint_null_rep(constructs, obj, B, rng, centred=centred)
    Sn = Zn.mean(axis=1)
    p = (1 + int(np.sum(Sn >= S))) / (B + 1)
    tot, tot_obs = Zn.sum(axis=1), float(z.sum())
    n = len(constructs)
    loo = {}
    for j, c in enumerate(constructs):
        sj = (tot_obs - z[j]) / (n - 1)
        loo[c] = dict(S=float(sj), p=(1 + int(np.sum((tot - Zn[:, j]) / (n - 1) >= sj))) / (B + 1))
    bs = z[rng.integers(0, n, (B_BOOT, n))].mean(axis=1)
    ci = [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))]
    d = np.asarray([0.0 if obj[c]["k"] == 0 else
                    float(obj[c]["dG"].mean() - (obj[c]["dbarG"].mean() if centred else 0.0))
                    for c in constructs])
    dbs = d[rng.integers(0, n, (B_BOOT, n))].mean(axis=1)
    lm = max(v["p"] for v in loo.values())
    lw = max(loo, key=lambda c: loo[c]["p"])
    return dict(label=label, n=n, centred=centred, S=S, ci95=ci, p_joint=p, B=B,
                floor=1 / (B + 1), null_mean=float(Sn.mean()), null_sd=float(Sn.std()),
                null_ci=[float(np.quantile(Sn, 0.025)), float(np.quantile(Sn, 0.975))],
                loo_max_p=lm, loo_worst=lw, loo={c: v for c, v in loo.items()},
                delta=float(d.mean()),
                delta_ci=[float(np.quantile(dbs, 0.025)), float(np.quantile(dbs, 0.975))],
                z=z.tolist(), Sn=Sn)


# ================================================================== C8's machinery, reproduced
P("\n--- C8's machinery, reproduced (KA-4 and KA-6 regression rows) ---")
POOLS8, OBJ8 = build(ALL, hard_windows=False, global_excl=False)
S8 = float(np.mean([z_unc(OBJ8[c]) for c in TARG]))
r307 = json.load(open(X.E307 / "results307.json"))
P(f"C8 uncentred S reproduced from C8's pools: {S8:+.5f} vs published "
  f"{r307['primary']['S']:+.5f}  (|diff| {abs(S8 - r307['primary']['S']):.2e})")

# KA-4 (unchanged machinery, C8's drop convention): construct scored against another's set
CANDIX = {c: np.asarray(sorted(PL.exclusions(c, rowmap[c]))) for c in ALL}


def ka4_unchanged(POOLS, n_der, rng):
    der = []
    for _ in range(n_der):
        perm = {}
        for s, cs in BY_SPLIT.items():
            while True:
                q = list(rng.permutation(cs))
                if all(a != b for a, b in zip(cs, q)):
                    break
            perm.update(dict(zip(cs, q)))
        zs = []
        for c in TARG:                       # score c's expression on perm[c]'s candidate set
            d = perm[c]
            bad = CANDIX[c]
            e, qs = EFF[c], []
            for p in POOLS[d]:
                if p["gi"] in bad:
                    continue
                mem = p["members"][~np.isin(p["members"], bad)]
                if mem.size < 20:
                    continue
                v, pv = e[p["gi"]], e[mem]
                qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
            zs.append(float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))) if qs else 0.0)
        der.append(float(np.mean(zs)))
    return np.asarray(der)


ka4 = ka4_unchanged(POOLS8, N_DERANGE, np.random.default_rng(X.SEED309 + 11))
ka4_ref = 0.0717
ka4_ok = bool(abs(ka4.mean() - ka4_ref) <= 0.010)
P(f"KA-4 (UNCHANGED machinery, {N_DERANGE:,} within-split derangements): mean {ka4.mean():+.5f} "
  f"(SE {ka4.std() / np.sqrt(ka4.size):.5f}), 95% [{np.quantile(ka4, 0.025):+.5f}, "
  f"{np.quantile(ka4, 0.975):+.5f}]; C8 published {ka4_ref:+.4f} -> "
  f"{'PASS' if ka4_ok else 'FAIL'} (rule: reproduction within +-0.010)")
P(f"     for the record, against C8's ORIGINAL |mean| <= 0.05 rule this row still "
  f"{'FAILS' if abs(ka4.mean()) > 0.05 else 'passes'} -- that is the defect C9 repairs")

# KA-6 (unchanged machinery): the nine SMDs must reproduce C8's
bal8 = PL.balance(TARG, POOLS8, SP, SPLIT_OF, PI)
ref8 = r307["ka"]["KA6"]["balance"]
d6 = {n: abs(bal8[n]["smd"] - ref8[n]["smd"]) for n in bal8}
ka6_ok = bool(max(d6.values()) <= 0.002)
P(f"KA-6 (UNCHANGED machinery): nine-axis SMD table reproduces C8 to max |diff| "
  f"{max(d6.values()):.2e} -> {'PASS' if ka6_ok else 'FAIL'} (rule: within 0.002)")
for n in bal8:
    P(f"     {n:22s} SMD {bal8[n]['smd']:+.4f}   (C8 {ref8[n]['smd']:+.4f})")

# ================================================================== C9's REPAIRED machinery
P("\n--- C9's REPAIRED machinery (150-cell strata + global L*-candidate exclusion + centring) ---")
POOLS9, OBJ9 = build(ALL, hard_windows=True, global_excl=True)
sizes = [p["members"].size for c in TARG for p in POOLS9[c]]
wid = sum(p["widened"] for c in TARG for p in POOLS9[c])
small = [(c, p["gene"], p["cell_n"]) for c in TARG for p in POOLS9[c] if p["cell_n"] < 20]
P(f"repaired pools over the 34 targeting constructs: {len(sizes)} pools, size mean "
  f"{np.mean(sizes):.1f}, min {min(sizes)}, max {max(sizes)}; {wid} needed the declared widening; "
  f"{len(small)} strata cells held < 20 eligible genes")
if small:
    P("     widened candidates: " + ", ".join(f"{c}/{g} (cell n={n})" for c, g, n in small))
nrec = {c: len(OBJ9[c]["recips"]) for c in TARG}
P("wrong recipients per construct, by split: "
  + ", ".join(f"{s} {len(cs) - 1} (n={len(cs)})" for s, cs in sorted(BY_SPLIT.items())))

# KA-6R: the nine axes on the repaired pools, C8's ORIGINAL rule
bal9 = PL.balance(TARG, POOLS9, SP, SPLIT_OF, PI)
ka6r_ok = all(abs(v["smd"]) <= 0.25 for v in bal9.values())
P(f"KA-6R (REPAIRED pools, C8's original |SMD| <= 0.25 rule on every axis) -> "
  f"{'PASS' if ka6r_ok else 'FAIL'}")
for n in bal9:
    P(f"     {n:22s} candidate {bal9[n]['cand_mean']:+9.4f}  pool {bal9[n]['pool_mean']:+9.4f}  "
      f"SMD {bal9[n]['smd']:+.4f}   (C8 {ref8[n]['smd']:+.4f})")

# the wrong-recipient levels
bvec = np.asarray([b_lvl(OBJ9[c]) for c in TARG])
zunc9 = np.asarray([z_unc(OBJ9[c]) for c in TARG])
P(f"\nwrong-recipient levels b_c over the 34 constructs: mean {bvec.mean():+.5f}, "
  f"sd {bvec.std(ddof=1):.5f}, range [{bvec.min():+.4f}, {bvec.max():+.4f}]; "
  f"uncentred S on the repaired pools {zunc9.mean():+.5f}")

# ---------------------------------------------------------------- THE PRIMARY
PRI = verdict_block(TARG, OBJ9, "C9 primary (repaired, centred)", X.SEED309, centred=True)
P(f"\nPRIMARY  S = {PRI['S']:+.5f}   (n = 34 targeting constructs, "
  f"{sum(OBJ9[c]['k'] for c in TARG)} candidate genes)")
P(f"joint null: B = {B_NULL:,}, floor p = {PRI['floor']:.2e}; null S mean "
  f"{PRI['null_mean']:+.5f} sd {PRI['null_sd']:.5f}, 95% "
  f"[{PRI['null_ci'][0]:+.5f}, {PRI['null_ci'][1]:+.5f}]   (C8's null: mean +0.00005 sd 0.17173)")
P(f"PRIMARY p_joint (one-sided, repression) = {PRI['p_joint']:.6g}  "
  f"[{'<=' if PRI['p_joint'] <= 0.05 else '>'} alpha = 0.05]")
P(f"leave-one-construct-out: max p = {PRI['loo_max_p']:.6g} (deleting {PRI['loo_worst']}), "
  f"min p = {min(v['p'] for v in PRI['loo'].values()):.6g}")
P(f"construct-level bootstrap (n = 34, {B_BOOT:,} resamples): 95% CI "
  f"[{PRI['ci95'][0]:+.5f}, {PRI['ci95'][1]:+.5f}]")
P(f"effect size Delta = {PRI['delta']:+.6f} log1p-CP50k per construct, 95% "
  f"[{PRI['delta_ci'][0]:+.6f}, {PRI['delta_ci'][1]:+.6f}]  (negative = candidates repressed)")
RULE = bool(PRI["p_joint"] <= 0.05 and PRI["loo_max_p"] <= 0.05)
P(f">>> pass rule (p_joint <= 0.05 AND max leave-one-out p <= 0.05): "
  f"{'BOTH CONDITIONS MET' if RULE else 'NOT MET'}")

# ---------------------------------------------------------------- KA-4R
ZC = {c: {d: z_in(OBJ9[c], j) for j, d in enumerate(OBJ9[c]["recips"])} for c in TARG}
BC = {c: b_lvl(OBJ9[c]) for c in TARG}
rng4 = np.random.default_rng(X.SEED309 + 13)
der = []
for _ in range(N_DERANGE):
    perm = {}
    for s, cs in BY_SPLIT.items():
        while True:
            q = list(rng4.permutation(cs))
            if all(a != b for a, b in zip(cs, q)):
                break
        perm.update(dict(zip(cs, q)))
    der.append(float(np.mean([ZC[perm[c]][c] - BC[perm[c]] for c in TARG])))
der = np.asarray(der)
ka4r_ok = bool(abs(der.mean()) <= 0.05
               and np.quantile(der, 0.025) <= 0 <= np.quantile(der, 0.975))
p_perm = (1 + int(np.sum(der >= PRI["S"]))) / (der.size + 1)
P(f"\nKA-4R (REPAIRED statistic, {N_DERANGE:,} within-split derangements): mean {der.mean():+.5f} "
  f"(SE {der.std() / np.sqrt(der.size):.5f}), 95% [{np.quantile(der, 0.025):+.5f}, "
  f"{np.quantile(der, 0.975):+.5f}] -> {'PASS' if ka4r_ok else 'FAIL'} "
  f"(C8's original rule: |mean| <= 0.05 and the interval brackets 0)")
P(f"     as declared in sec 6.1 this row is zero in EXPECTATION BY CONSTRUCTION; its content is "
  f"correct implementation + Monte-Carlo bracketing.  Its draw distribution also gives the "
  f"wrong-recipient permutation p for the observed S: {p_perm:.6g} (floor {1 / (der.size + 1):.2e})")

# ---------------------------------------------------------------- KA-5 on both nulls
def ka5(Sn, seed):
    r = np.random.default_rng(seed)
    pick = r.choice(Sn.size, N_CALIB, replace=False)
    pc = np.asarray([(1 + int(np.sum(np.delete(Sn, i) >= Sn[i]))) / Sn.size for i in pick])
    ks = stats.kstest(pc, "uniform")
    f5 = float(np.mean(pc <= 0.05))
    return float(ks.pvalue), f5, bool(ks.pvalue >= 0.05 and 0.02 <= f5 <= 0.10)


UNC8 = verdict_block(TARG, OBJ8, "C8 uncentred, reproduced", X.SEED309 + 3, centred=False)
k5a = ka5(PRI["Sn"], X.SEED309 + 17)
k5b = ka5(UNC8["Sn"], X.SEED309 + 19)
ka5_ok = bool(k5a[2] and k5b[2])
P(f"KA-5 null calibration ({N_CALIB} draws each p-valued against the rest): REPAIRED null "
  f"KS p = {k5a[0]:.4f}, frac<=0.05 = {k5a[1]:.4f}; C8's null KS p = {k5b[0]:.4f}, "
  f"frac<=0.05 = {k5b[1]:.4f} -> {'PASS' if ka5_ok else 'FAIL'} (both)")
P(f"     (C8's own reproduction check: uncentred S on C8 pools {UNC8['S']:+.5f}, "
  f"p_joint {UNC8['p_joint']:.6g}, published +0.78567 / 9.9999e-06)")

KA = dict(KA0="reported in ka309.json", KA1="ka309.json", KA2="ka309.json", KA3="ka309.json",
          KA7="ka309.json",
          KA4=dict(mean=float(ka4.mean()), se=float(ka4.std() / np.sqrt(ka4.size)),
                   ci=[float(np.quantile(ka4, 0.025)), float(np.quantile(ka4, 0.975))],
                   c8_reference=ka4_ref, n=N_DERANGE, rule="reproduction within +-0.010",
                   passed=ka4_ok, fails_c8_original_rule=bool(abs(ka4.mean()) > 0.05)),
          KA5=dict(repaired=dict(ks_p=k5a[0], frac=k5a[1], passed=k5a[2]),
                   c8=dict(ks_p=k5b[0], frac=k5b[1], passed=k5b[2]), n=N_CALIB, passed=ka5_ok),
          KA6=dict(balance=bal8, c8_reference=ref8, max_abs_diff=float(max(d6.values())),
                   rule="reproduction within 0.002", passed=ka6_ok),
          KA4R=dict(mean=float(der.mean()), se=float(der.std() / np.sqrt(der.size)),
                    ci=[float(np.quantile(der, 0.025)), float(np.quantile(der, 0.975))],
                    n=N_DERANGE, rule="|mean| <= 0.05 and interval brackets 0", passed=ka4r_ok,
                    permutation_p_for_S=p_perm),
          KA6R=dict(balance=bal9, rule="|SMD| <= 0.25 on every axis", passed=ka6r_ok))

# ================================================================== secondaries
sec: dict = {}

# S1 failing vs working, on the repaired statistic
zrep = {c: z_rep(OBJ9[c]) for c in TARG}
rngS = np.random.default_rng(X.SEED309 + 23)
zf = np.asarray([zrep[c] for c in TARG if rowmap[c]["group"] == "failing"])
zw = np.asarray([zrep[c] for c in TARG if rowmap[c]["group"] == "working"])
db = np.asarray([zf[rngS.integers(0, zf.size, zf.size)].mean()
                 - zw[rngS.integers(0, zw.size, zw.size)].mean() for _ in range(B_BOOT)])
mw = stats.mannwhitneyu(zf, zw, alternative="two-sided")
sec["S1_failing_vs_working"] = dict(n_failing=int(zf.size), n_working=int(zw.size),
                                    mean_failing=float(zf.mean()), mean_working=float(zw.mean()),
                                    diff=float(zf.mean() - zw.mean()),
                                    ci95=[float(np.quantile(db, 0.025)),
                                          float(np.quantile(db, 0.975))],
                                    mannwhitney_p=float(mw.pvalue))
P(f"\nS1 failing ({zf.size}) vs working ({zw.size}): z_rep {zf.mean():+.4f} vs {zw.mean():+.4f}, "
  f"diff {zf.mean() - zw.mean():+.4f} [{np.quantile(db, 0.025):+.4f}, "
  f"{np.quantile(db, 0.975):+.4f}], MWU p {mw.pvalue:.4f}   (C8: -0.077 [-0.869, +0.699], p 0.595)")

# S2 the 31 non-targeting constructs, identical pipeline
NTC = verdict_block(NTCC, OBJ9, "31 non-targeting constructs (repaired)", X.SEED309 + 5)
sec["S2_non_targeting"] = {k: v for k, v in NTC.items() if k not in ("Sn", "loo", "z")}
P(f"S2 the 31 non-targeting constructs (candidate-bearing observations, NOT negative controls): "
  f"S = {NTC['S']:+.5f} [{NTC['ci95'][0]:+.5f}, {NTC['ci95'][1]:+.5f}], "
  f"joint-null p = {NTC['p_joint']:.6g}   (C8: +0.2857 [-0.0691, +0.6358], p 0.0537)")

# S3 partner-construct specificity, UNCENTRED, on the repaired pools
by_target = collections.defaultdict(list)
for c in TARG:
    by_target[rowmap[c]["intended_target"]].append(c)
partner = {}
for t, cs in by_target.items():
    if len(cs) == 2:
        partner[cs[0]], partner[cs[1]] = cs[1], cs[0]
PART = {}
for c in TARG:
    if c not in partner:
        continue
    d = partner[c]
    o = OBJ9[c]
    if d in o["recips"]:
        PART[c] = z_in(o, o["recips"].index(d))
    else:                       # partner in another split: score it directly
        PART[c] = float(np.sqrt(12 * o["k"])
                        * (0.5 - np.mean([PL.q_of(EFF[d], p) for p in o["pools"]])))
pc = [c for c in TARG if c in PART]
zo = np.asarray([z_unc(OBJ9[c]) for c in pc])
zp = np.asarray([PART[c] for c in pc])
bp = zp[rngS.integers(0, zp.size, (B_BOOT, zp.size))].mean(axis=1)
dfp = np.asarray([zo[rngS.integers(0, zo.size, zo.size)].mean()
                  - zp[rngS.integers(0, zp.size, zp.size)].mean() for _ in range(B_BOOT)])
sec["S3_partner"] = dict(n=len(pc), z_own=float(zo.mean()), z_partner=float(zp.mean()),
                         partner_ci=[float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975))],
                         diff=float(zo.mean() - zp.mean()),
                         diff_ci=[float(np.quantile(dfp, 0.025)), float(np.quantile(dfp, 0.975))],
                         wilcoxon_p=float(stats.wilcoxon(zo, zp).pvalue),
                         same_split=int(sum(1 for c in pc if partner[c] in OBJ9[c]["recips"])),
                         per_construct={c: dict(own=float(z_unc(OBJ9[c])), partner=float(PART[c]))
                                        for c in pc})
P(f"S3 partner-construct specificity on the REPAIRED pools, uncentred (n = {len(pc)}): own "
  f"{zo.mean():+.5f} vs PARTNER (same intended target, different protospacers) {zp.mean():+.5f} "
  f"[{np.quantile(bp, 0.025):+.5f}, {np.quantile(bp, 0.975):+.5f}]; own - partner "
  f"{zo.mean() - zp.mean():+.5f} [{np.quantile(dfp, 0.025):+.5f}, {np.quantile(dfp, 0.975):+.5f}], "
  f"paired Wilcoxon p {stats.wilcoxon(zo, zp).pvalue:.5f}   (C8: +0.041 vs +0.786, diff +0.745)")

# S6 by seed level
sub = {}
zt = np.asarray([zrep[c] for c in TARG])
for lab in (12, 11, 10, 9):
    cs = [c for c in TARG if int(rowmap[c]["L_star"]) == lab]
    if not cs:
        continue
    blk = verdict_block(cs, OBJ9, f"L*={lab}", X.SEED309 + 100 + lab)
    sub[f"L*={lab}"] = dict(n=blk["n"], S=blk["S"], p_joint=blk["p_joint"], ci95=blk["ci95"],
                            constructs=cs)
    P(f"S6 L*={lab:2d}  n={blk['n']:2d}  S = {blk['S']:+.5f}  joint-null p = {blk['p_joint']:.6g}")
sec["S6_by_level"] = sub

# sign count + measurement floor
npos = int(np.sum(zt > 0))
sgn = stats.binomtest(npos, zt.size, 0.5, alternative="greater")
below = [(c, p["gene"], float(SP[SPLIT_OF[c]].raw_cp50k[p["gi"]]))
         for c in TARG for p in POOLS9[c] if SP[SPLIT_OF[c]].raw_cp50k[p["gi"]] < FLOOR_CP50K]
zab = []
for c in TARG:
    o = OBJ9[c]
    keep = [i for i, p in enumerate(o["pools"])
            if SP[SPLIT_OF[c]].raw_cp50k[p["gi"]] >= FLOOR_CP50K]
    if not keep:
        zab.append(0.0)
        continue
    zab.append(float(np.sqrt(12 * len(keep))
                     * (o["qbarG"][keep].mean() - o["qG"][keep].mean())))
sec["S11_sign_and_floor"] = dict(n_pos=npos, n=int(zt.size), binom_p=float(sgn.pvalue),
                                 median_z=float(np.median(zt)),
                                 n_candidates=int(sum(OBJ9[c]["k"] for c in TARG)),
                                 n_below_floor=len(below),
                                 S_above_floor_only=float(np.mean(zab)))
P(f"S11 sign count {npos}/34 constructs with z_rep > 0 (exact binomial p = {sgn.pvalue:.6g}), "
  f"median {np.median(zt):+.5f}; {len(below)}/{sum(OBJ9[c]['k'] for c in TARG)} candidates below "
  f"the CP50k {FLOOR_CP50K} floor, S above-floor only {np.mean(zab):+.5f}   (C8: 24/34, +0.6907)")

# ================================================================== the per-construct CSV
with open(X.HERE / "enrichment309.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["construct", "unit", "group", "intended_target", "split", "kd_lfc", "cells",
                "L_star", "k_candidates", "candidates", "pool_size_mean", "pools_widened",
                "n_wrong_recipients", "r_c", "z_c_uncentred", "b_c_wrong_recipient",
                "z_c_repaired", "delta_own", "delta_wrong_recipient", "delta_repaired",
                "null_z_mean", "null_z_sd", "null_z_p2.5", "null_z_p97.5", "p_construct",
                "z_partner_uncentred", "n_below_floor", "min_cand_cp50k"])
    for grp, blk in ((TARG, PRI), (NTCC, NTC)):
        Zn = joint_null_rep(grp, OBJ9, B_NULL, np.random.default_rng(X.SEED309 + 31))
        for j, c in enumerate(grp):
            r, o = rowmap[c], OBJ9[c]
            zn = Zn[:, j]
            cp = [float(SP[SPLIT_OF[c]].raw_cp50k[p["gi"]]) for p in o["pools"]]
            w.writerow([c, r["unit"], r["group"], r["intended_target"], r["split"], r["kd_lfc"],
                        r["cells_full"], r["L_star"], o["k"], r["candidates"],
                        round(float(np.mean([p["members"].size for p in o["pools"]])), 1)
                        if o["pools"] else 0,
                        sum(p["widened"] for p in o["pools"]), len(o["recips"]),
                        round(float(o["qG"].mean()), 5) if o["k"] else "",
                        round(z_unc(o), 5), round(b_lvl(o), 5), round(z_rep(o), 5),
                        round(float(o["dG"].mean()), 6) if o["k"] else 0.0,
                        round(float(o["dbarG"].mean()), 6) if o["k"] else 0.0,
                        round(float(o["dG"].mean() - o["dbarG"].mean()), 6) if o["k"] else 0.0,
                        round(float(zn.mean()), 5), round(float(zn.std()), 5),
                        round(float(np.quantile(zn, 0.025)), 5),
                        round(float(np.quantile(zn, 0.975)), 5),
                        round((1 + int(np.sum(zn >= z_rep(o)))) / (zn.size + 1), 6),
                        round(PART[c], 5) if c in PART else "",
                        sum(1 for v in cp if v < FLOOR_CP50K), round(min(cp), 4) if cp else ""])
P(f"\nwrote enrichment309.csv ({len(ALL)} constructs)")

res = dict(
    preregistration_sha256=open(X.HERE / "PREREGISTRATION.sha256").read().split()[0],
    primary={k: v for k, v in PRI.items() if k != "Sn"},
    pass_rule_conditions_met=RULE,
    c8_reproduction=dict(S_uncentred_c8_pools=S8, published=r307["primary"]["S"],
                         abs_diff=abs(S8 - r307["primary"]["S"]),
                         p_joint=UNC8["p_joint"], published_p=r307["primary"]["p_joint"]),
    repaired_pools=dict(n_pools=len(sizes), size_mean=float(np.mean(sizes)), size_min=int(min(sizes)),
                        size_max=int(max(sizes)), n_widened=int(wid),
                        n_cells_below_20=len(small),
                        small=[[c, g, n] for c, g, n in small],
                        recipients_per_construct=nrec),
    wrong_recipient_levels=dict(mean=float(bvec.mean()), sd=float(bvec.std(ddof=1)),
                                min=float(bvec.min()), max=float(bvec.max()),
                                per_construct={c: float(BC[c]) for c in TARG},
                                uncentred_S_repaired_pools=float(zunc9.mean())),
    ka=KA, secondaries=sec,
    peak_rss_gib=round(X.peak_gib(), 3), seconds=round(time.time() - t0, 1))
(X.HERE / "results309.json").write_text(json.dumps(res, indent=1, default=float))
np.savez_compressed(X.HERE / "state309.npz", z_rep=np.asarray([zrep[c] for c in TARG]),
                    b=bvec, z_unc=zunc9, targ=np.asarray(TARG), Sn=PRI["Sn"], der=der, ka4=ka4)
(X.HERE / "c9.log").write_text("\n".join(log) + "\n")
P(f"done {res['seconds']}s, peak {res['peak_rss_gib']:.2f} GiB")

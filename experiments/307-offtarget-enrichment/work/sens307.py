"""307 step 5 -- the declared second null (PREREGISTRATION sec 8.5), the extended wrong-recipient
permutation, and a POST-HOC repaired-matching sensitivity.

  S5   real-spacer resampling null.  A pseudo-construct is a random PAIR of stringency-0
       hCRISPRi-v2 protospacers rescanned through 301's exact rule; a draw is kept only if it
       yields L* = L*(c) and, after c's exclusions, exactly k_c candidates.  Scoring is
       matching-free: each candidate's percentile among the WHOLE eligible universe at L*(c).
       This null preserves the sequence rule jointly and does NOT match expression -- it is the
       complement of the primary null, not a replacement.
  P1   the KA-4 wrong-recipient derangement null, extended from 200 to 5,000 draws (post-hoc
       Monte-Carlo extension of a preregistered statistic, no new statistic).
  P2   POST-HOC repaired matching: promoter-window count added to the hard strata, which is what
       KA-6 failed on.  Reported without a verdict.

    PYTHONDONTWRITEBYTECODE=1 uv run python sens307.py
"""
from __future__ import annotations

import collections
import json
import pickle
import time
from pathlib import Path

import numpy as np
from scipy import stats

import lib307 as X

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


LS = (12, 11, 10, 9)
N_PAIRS = 2_000_000
MAX_ACC = 5_000
MIN_ACC = 200
N_DERANGE2 = 5_000
rng = np.random.default_rng(X.SEED307 + 7)

SP = {s: X.Split(s) for s in X.SPLITS}
GENES = SP["train"].genes
GIX = SP["train"].gix
NG = GENES.size
cons_rows = X.load_constructs()
rowmap = {r["construct"]: r for r in cons_rows}
TARG = [r["construct"] for r in cons_rows if r["unit"] == "targeting"]
gene_meta = X.load_genes()
lookup = pickle.load(open(X.HERE / "seq307_lookup.pkl", "rb"))
LIB_TARGETS = set(lookup["lib_targets"])
CAND_L = lookup["cand_L"]
D = pickle.load(open(X.HERE / "poolhits307.pkl", "rb"))
pool, hits, pi_off = D["pool"], D["hits"], D["pi_off"]

PI = {L: np.zeros(NG) for L in LS}
for L in LS:
    for g, v in pi_off[L].items():
        PI[L][GIX[g]] = v
IS_LIB = np.zeros(NG, dtype=bool)
for g in LIB_TARGETS:
    IS_LIB[GIX[g]] = True
HAS_WIN = np.zeros(NG, dtype=bool)
for g in gene_meta:
    HAS_WIN[GIX[g]] = True
BASE = {L: np.flatnonzero(HAS_WIN & (PI[L] > 0) & ~IS_LIB) for L in LS}

SPLIT_OF = {r["construct"]: (r["split"] if r["unit"] == "targeting" else "train")
            for r in cons_rows}
EFF = {c: SP[SPLIT_OF[c]].effect(c, loo_control=(rowmap[c]["unit"] != "targeting"))
       for c in [r["construct"] for r in cons_rows]}


def excl_set(c: str) -> set:
    bad = set()
    if rowmap[c]["intended_target"]:
        bad.add(rowmap[c]["intended_target"])
    for L in LS:
        bad |= set(CAND_L[c][L])
    return bad


# ======================================================== S5 real-spacer resampling null
hv = pickle.load(open(X.E301 / "horlbeck_supp3.pkl", "rb"))["hCRISPRi-v2"]
ixh = {h: i for i, h in enumerate(hv[0])}


def clean(v) -> str:
    v = str(v)
    return v[2:-1] if v.startswith('="') and v.endswith('"') else v


strg = np.asarray([clean(r[ixh["off-target stringency"]]) for r in hv[1]
                   if len(clean(r[ixh["protospacer sequence"]])) >= 12
                   and not (set(clean(r[ixh["protospacer sequence"]]).upper()) - set("ACGT"))])
assert strg.size == len(pool)
keep = np.flatnonzero(strg == "0")          # all 68 real 2025 spacers carry stringency 0
P(f"S5 pseudo-spacer pool: {keep.size:,} stringency-0 hCRISPRi-v2 protospacers "
  f"of {len(pool):,}")

firstL = np.zeros(len(pool), dtype=np.int8)
for L in LS:                                 # hits are nested (a 12-mer match implies an 11-mer one)
    for i, h in enumerate(hits[L]):
        if h and firstL[i] == 0:
            firstL[i] = L
ia = keep[rng.integers(0, keep.size, N_PAIRS)]
ib = keep[rng.integers(0, keep.size, N_PAIRS)]
lstar = np.maximum(firstL[ia], firstL[ib])
buckets = {L: np.flatnonzero(lstar == L) for L in LS}
P("S5 pseudo-construct pairs available: "
  + ", ".join(f"L*={L}: {buckets[L].size:,}" for L in LS))

UNIV_Q = {}          # (construct) -> percentile of every eligible gene, in that construct's cells


def universe_percentile(c: str, L: int, bad: set):
    ix = BASE[L]
    if bad:
        badix = np.asarray([GIX[g] for g in bad if g in GIX])
        ix = ix[~np.isin(ix, badix)]
    v = EFF[c][ix]
    order = np.argsort(v, kind="stable")
    sv = v[order]
    return ix, v, sv


s5 = {}
for c in TARG:
    L = int(rowmap[c]["L_star"])
    k = int(rowmap[c]["k"])
    bad = excl_set(c)
    ix, v, sv = universe_percentile(c, L, bad)
    n = v.size
    pos = {int(g): i for i, g in enumerate(ix)}

    def perc(gene_idxs):
        vals = EFF[c][gene_idxs]
        lo = np.searchsorted(sv, vals, side="left")
        hi = np.searchsorted(sv, vals, side="right")
        return float(np.mean((lo + 0.5 * (hi - lo)) / n)) if len(vals) else np.nan

    obs = perc(np.asarray([GIX[g] for g in rowmap[c]["candidates"].split(",")]))
    W_obs = float(np.sqrt(12 * k) * (0.5 - obs))
    acc, tried = [], 0
    for j in buckets[L]:
        tried += 1
        u = (set(hits[L][ia[j]]) | set(hits[L][ib[j]])) - bad
        u = {g for g in u if g in GIX and not IS_LIB[GIX[g]]}
        if len(u) != k:
            continue
        acc.append(np.sqrt(12 * k) * (0.5 - perc(np.asarray([GIX[g] for g in u]))))
        if len(acc) >= MAX_ACC:
            break
    a = np.asarray(acc)
    s5[c] = dict(L=L, k=k, W_obs=W_obs, n_accepted=int(a.size), pairs_tried=int(tried),
                 null_mean=float(a.mean()) if a.size else None,
                 null_sd=float(a.std()) if a.size else None,
                 p=float((1 + int(np.sum(a >= W_obs))) / (a.size + 1)) if a.size else None)
ok = [c for c in TARG if s5[c]["n_accepted"] >= MIN_ACC]
P(f"S5 constructs with >= {MIN_ACC} accepted size- and level-matched pseudo-constructs: "
  f"{len(ok)}/{len(TARG)}  (short: {[c for c in TARG if c not in ok]})")
if ok:
    W_obs_agg = float(np.mean([s5[c]["W_obs"] for c in ok]))
    nb = min(s5[c]["n_accepted"] for c in ok)
    # joint draws: the b-th accepted pseudo-construct of every construct
    Wn = np.zeros((nb, len(ok)))
    for j, c in enumerate(ok):
        L, k, bad = s5[c]["L"], s5[c]["k"], excl_set(c)
        ix, v, sv = universe_percentile(c, L, bad)
        n = v.size
        got = 0
        for jj in buckets[L]:
            u = (set(hits[L][ia[jj]]) | set(hits[L][ib[jj]])) - bad
            u = {g for g in u if g in GIX and not IS_LIB[GIX[g]]}
            if len(u) != k:
                continue
            vals = EFF[c][np.asarray([GIX[g] for g in u])]
            lo = np.searchsorted(sv, vals, side="left")
            hi = np.searchsorted(sv, vals, side="right")
            Wn[got, j] = np.sqrt(12 * k) * (0.5 - np.mean((lo + 0.5 * (hi - lo)) / n))
            got += 1
            if got >= nb:
                break
    Sn5 = Wn.mean(axis=1)
    p5 = (1 + int(np.sum(Sn5 >= W_obs_agg))) / (Sn5.size + 1)
    P(f"S5 aggregate over the {len(ok)} resolvable constructs: W_obs = {W_obs_agg:+.5f}, "
      f"null mean {Sn5.mean():+.5f} sd {Sn5.std():.5f} "
      f"95% [{np.quantile(Sn5, 0.025):+.5f}, {np.quantile(Sn5, 0.975):+.5f}], "
      f"p = {p5:.6g} ({Sn5.size:,} joint draws, floor {1 / (Sn5.size + 1):.2e})")
else:
    W_obs_agg = p5 = None
    Sn5 = np.asarray([])

# ======================================================== P1 extended wrong-recipient permutation
P("\nP1 extended wrong-recipient derangement null: rebuilding pools ...")
import enrich_pools307 as EP          # noqa: E402  -- pool builder shared with enrich307

EP.setup(SP, cons_rows, gene_meta, lookup, pi_off)
POOLS = {c: EP.build_pools(c, rowmap[c]) for c in [r["construct"] for r in cons_rows]}
ZOBS = {}
for c in POOLS:
    ZOBS[c] = EP.z_of(EFF[c], POOLS[c])
S_obs = float(np.mean([ZOBS[c] for c in TARG]))
_ref = json.load(open(X.HERE / "results307.json"))["primary"]["S"]
assert abs(S_obs - _ref) < 1e-12, (S_obs, _ref)
P(f"   pool rebuild reproduces enrich307.py's S: {S_obs:+.5f} vs {_ref:+.5f} "
  f"(|diff| {abs(S_obs - _ref):.2e})")
by_split = collections.defaultdict(list)
for c in TARG:
    by_split[SPLIT_OF[c]].append(c)
der = []
rng2 = np.random.default_rng(X.SEED307 + 11)
for _ in range(N_DERANGE2):
    perm = {}
    for s, cs in by_split.items():
        while True:
            q = list(rng2.permutation(cs))
            if all(a != b for a, b in zip(cs, q)):
                break
        perm.update(dict(zip(cs, q)))
    zs = []
    for c in TARG:
        d = perm[c]
        bad = excl_set(c)
        badix = np.asarray([GIX[g] for g in bad if g in GIX])
        e, qs = EFF[c], []
        for p in POOLS[d]:
            if p["gi"] in badix:
                continue
            mem = p["members"][~np.isin(p["members"], badix)]
            if mem.size < 20:
                continue
            v, pv = e[p["gi"]], e[mem]
            qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
        zs.append(float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))) if qs else 0.0)
    der.append(float(np.mean(zs)))
der = np.asarray(der)
p_perm = (1 + int(np.sum(der >= S_obs))) / (der.size + 1)
P(f"P1 {N_DERANGE2:,} within-split derangements: mean {der.mean():+.5f} "
  f"(SE {der.std() / np.sqrt(der.size):.5f}), sd {der.std():.5f}, "
  f"95% [{np.quantile(der, 0.025):+.5f}, {np.quantile(der, 0.975):+.5f}]")
P(f"   S_obs = {S_obs:+.5f}  ->  wrong-recipient permutation p = {p_perm:.6g} "
  f"(floor {1 / (der.size + 1):.2e}); this reading is immune to the matched-pool "
  f"construction because any matching bias enters observed and permuted alike")

# ======================================================== P2 post-hoc repaired matching
P("\nP2 POST-HOC repaired matching (promoter-window count added to the hard strata) -- "
  "reported without a verdict:")
EP.setup(SP, cons_rows, gene_meta, lookup, pi_off, hard_windows=True)
POOLS2 = {c: EP.build_pools(c, rowmap[c]) for c in [r["construct"] for r in cons_rows]}
Z2 = {c: EP.z_of(EFF[c], POOLS2[c]) for c in POOLS2}
S2 = float(np.mean([Z2[c] for c in TARG]))
Zn2 = EP.joint_null(TARG, POOLS2, EFF, 100_000, np.random.default_rng(X.SEED307 + 13))
Sn2 = Zn2.mean(axis=1)
p2 = (1 + int(np.sum(Sn2 >= S2))) / (Sn2.size + 1)
bal2 = EP.balance(TARG, POOLS2, SP, SPLIT_OF, PI)
der2 = []
for _ in range(500):
    perm = {}
    for s, cs in by_split.items():
        while True:
            q = list(rng2.permutation(cs))
            if all(a != b for a, b in zip(cs, q)):
                break
        perm.update(dict(zip(cs, q)))
    zs = []
    for c in TARG:
        d = perm[c]
        bad = excl_set(c)
        badix = np.asarray([GIX[g] for g in bad if g in GIX])
        e, qs = EFF[c], []
        for p in POOLS2[d]:
            if p["gi"] in badix:
                continue
            mem = p["members"][~np.isin(p["members"], badix)]
            if mem.size < 20:
                continue
            v, pv = e[p["gi"]], e[mem]
            qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
        zs.append(float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs))) if qs else 0.0)
    der2.append(float(np.mean(zs)))
der2 = np.asarray(der2)
P(f"   S = {S2:+.5f}, joint-null p = {p2:.6g}; KA-4-style bias {der2.mean():+.5f} "
  f"(SE {der2.std() / np.sqrt(der2.size):.5f}); worst |SMD| "
  f"{max(abs(v['smd']) for v in bal2.values()):.4f} on "
  f"{max(bal2, key=lambda k: abs(bal2[k]['smd']))}")
for n, v in bal2.items():
    P(f"     {n:22s} SMD {v['smd']:+.4f}")

out = dict(S5=dict(per_construct=s5, resolvable=ok,
                   aggregate=dict(W_obs=W_obs_agg, p=p5, n_draws=int(Sn5.size),
                                  null_mean=float(Sn5.mean()) if Sn5.size else None,
                                  null_sd=float(Sn5.std()) if Sn5.size else None)),
           P1=dict(n=int(der.size), mean=float(der.mean()), sd=float(der.std()),
                   se=float(der.std() / np.sqrt(der.size)),
                   ci=[float(np.quantile(der, 0.025)), float(np.quantile(der, 0.975))],
                   S_obs=S_obs, p_perm=p_perm, floor=1 / (der.size + 1)),
           P2=dict(S=S2, p_joint=p2, bias=float(der2.mean()),
                   bias_se=float(der2.std() / np.sqrt(der2.size)), balance=bal2),
           peak_rss_gib=round(X.peak_gib(), 3), seconds=round(time.time() - t0, 1))
(X.HERE / "sens307.json").write_text(json.dumps(out, indent=1, default=float))
(X.HERE / "sens307.log").write_text("\n".join(log) + "\n")
P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")

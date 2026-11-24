"""309 step 3 -- the declared secondaries that need their own machinery (PREREGISTRATION sec 10).

  Sec-1  the real-spacer resampling null (C8's S5), uncentred AND centred by the C9 construction
  Sens-1 repair variant (b): scalar subtraction with a propagated interval
  Sens-2 centring over all 33 other targeting constructs instead of within-split only
  Sens-3 the global L*-candidate exclusion switched off
  Sens-4 each repair alone: the KA-6 repair uncentred, and the KA-4 repair on C8's 50-cell pools

    PYTHONDONTWRITEBYTECODE=1 uv run python sens309.py
"""
from __future__ import annotations

import json
import pickle
import time

import numpy as np

import lib309 as X
import pools309 as PL
import c9 as C          # re-runs the primary in-process; its assertions guard this script

t0 = time.time()
log: list[str] = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    log.append(s)


LS = X.LS
N_PAIRS = 2_000_000
MAX_ACC = 5_000
MIN_ACC = 200
rng = np.random.default_rng(X.SEED309 + 7)
TARG, NTCC, ALL = C.TARG, C.NTCC, C.ALL
EFF, rowmap, GIX = C.EFF, C.rowmap, C.GIX
OBJ9, OBJ8, PI = C.OBJ9, C.OBJ8, C.PI
out: dict = {}

# ============================================================ Sec-1  real-spacer resampling null
D = pickle.load(open(X.E307 / "poolhits307.pkl", "rb"))
pool, hits = D["pool"], D["hits"]
lookup = C.lookup
IS_LIB = {g for g in lookup["lib_targets"]}
hv = pickle.load(open(X.E301 / "horlbeck_supp3.pkl", "rb"))["hCRISPRi-v2"]
ixh = {h: i for i, h in enumerate(hv[0])}


def clean(v) -> str:
    v = str(v)
    return v[2:-1] if v.startswith('="') and v.endswith('"') else v


strg = np.asarray([clean(r[ixh["off-target stringency"]]) for r in hv[1]
                   if len(clean(r[ixh["protospacer sequence"]])) >= 12
                   and not (set(clean(r[ixh["protospacer sequence"]]).upper()) - set("ACGT"))])
assert strg.size == len(pool)
keep = np.flatnonzero(strg == "0")
P(f"Sec-1 pseudo-spacer pool: {keep.size:,} stringency-0 hCRISPRi-v2 protospacers of {len(pool):,}")
firstL = np.zeros(len(pool), dtype=np.int8)
for L in LS:
    for i, h in enumerate(hits[L]):
        if h and firstL[i] == 0:
            firstL[i] = L
ia = keep[rng.integers(0, keep.size, N_PAIRS)]
ib = keep[rng.integers(0, keep.size, N_PAIRS)]
lstar = np.maximum(firstL[ia], firstL[ib])
buckets = {L: np.flatnonzero(lstar == L) for L in LS}

# the eligible universe at each L, as C9's repaired pools define it (global exclusion applied)
PL.setup(C.SP, C.cons_rows, C.gene_meta, lookup, C.pi_off, hard_windows=True, global_excl=True)
BASE = {L: PL._S["BASE"][L] for L in LS}


def excl_set(c):
    bad = {rowmap[c]["intended_target"]} if rowmap[c]["intended_target"] else set()
    for L in LS:
        bad |= set(lookup["cand_L"][c][L])
    return bad


def universe(c, L, bad):
    ix = BASE[L]
    badix = np.asarray([GIX[g] for g in bad if g in GIX], dtype=np.int64)
    if badix.size:
        ix = ix[~np.isin(ix, badix)]
    return ix


sec1 = {}
for c in TARG:
    L, k = int(rowmap[c]["L_star"]), int(rowmap[c]["k"])
    bad = excl_set(c)
    ix = universe(c, L, bad)
    rec = OBJ9[c]["recips"]
    sv = {d: np.sort(EFF[d][ix]) for d in [c] + list(rec)}
    n = ix.size

    def perc(gi_arr, d):
        vals = EFF[d][gi_arr]
        lo = np.searchsorted(sv[d], vals, side="left")
        hi = np.searchsorted(sv[d], vals, side="right")
        return float(np.mean((lo + 0.5 * (hi - lo)) / n))

    gi_obs = np.asarray([GIX[g] for g in rowmap[c]["candidates"].split(",")], dtype=np.int64)
    s = np.sqrt(12 * k)
    W_own = s * (0.5 - perc(gi_obs, c))
    W_wr = float(np.mean([s * (0.5 - perc(gi_obs, d)) for d in rec]))
    acc_own, acc_rep, tried = [], [], 0
    for j in buckets[L]:
        tried += 1
        u = (set(hits[L][ia[j]]) | set(hits[L][ib[j]])) - bad
        u = {g for g in u if g in GIX and g not in IS_LIB}      # C8's S5 filter, unchanged
        if len(u) != k:
            continue
        gi = np.asarray([GIX[g] for g in u], dtype=np.int64)
        w_own = s * (0.5 - perc(gi, c))
        w_wr = float(np.mean([s * (0.5 - perc(gi, d)) for d in rec]))
        acc_own.append(w_own)
        acc_rep.append(w_own - w_wr)
        if len(acc_own) >= MAX_ACC:
            break
    sec1[c] = dict(L=L, k=k, W_own=W_own, W_wr=W_wr, W_rep=W_own - W_wr,
                   n_accepted=len(acc_own), pairs_tried=tried,
                   own=np.asarray(acc_own), rep=np.asarray(acc_rep))
ok = [c for c in TARG if sec1[c]["n_accepted"] >= MIN_ACC]
P(f"Sec-1 constructs with >= {MIN_ACC} size- and level-matched pseudo-constructs: "
  f"{len(ok)}/{len(TARG)}" + ("" if len(ok) == len(TARG) else f"  short: {set(TARG) - set(ok)}"))
nb = min(sec1[c]["n_accepted"] for c in ok)
res1 = {}
for tag in ("own", "rep"):
    W = np.stack([sec1[c][tag][:nb] for c in ok], axis=1)
    obs = float(np.mean([sec1[c]["W_own" if tag == "own" else "W_rep"] for c in ok]))
    Sn = W.mean(axis=1)
    p = (1 + int(np.sum(Sn >= obs))) / (Sn.size + 1)
    res1[tag] = dict(W_obs=obs, p=p, n_draws=int(Sn.size), null_mean=float(Sn.mean()),
                     null_sd=float(Sn.std()),
                     null_ci=[float(np.quantile(Sn, 0.025)), float(np.quantile(Sn, 0.975))])
    P(f"Sec-1 {'UNCENTRED' if tag == 'own' else 'CENTRED  '} real-spacer null over {len(ok)} "
      f"constructs: W_obs = {obs:+.5f}, null mean {Sn.mean():+.5f} sd {Sn.std():.5f} "
      f"95% [{np.quantile(Sn, 0.025):+.5f}, {np.quantile(Sn, 0.975):+.5f}], p = {p:.6g} "
      f"({Sn.size:,} joint draws, floor {1 / (Sn.size + 1):.2e})")
P(f"      (C8's uncentred reading: W_obs +0.74470, p 0.00149)")
out["Sec1_real_spacer_null"] = dict(
    resolvable=len(ok), n_joint_draws=nb, uncentred=res1["own"], centred=res1["rep"],
    per_construct={c: {k: v for k, v in sec1[c].items() if k not in ("own", "rep")} for c in TARG})

# ============================================================ Sens-1  repair variant (b)
z8 = np.asarray([C.z_unc(OBJ9[c]) for c in TARG])
b = np.asarray([C.b_lvl(OBJ9[c]) for c in TARG])
rngb = np.random.default_rng(X.SEED309 + 41)
bs = z8[rngb.integers(0, 34, (20_000, 34))].mean(axis=1)
se_b = float(b.std(ddof=1) / np.sqrt(34))
lo = float(np.quantile(bs, 0.025)) - b.mean() - 1.96 * se_b
hi = float(np.quantile(bs, 0.975)) - b.mean() + 1.96 * se_b
P(f"\nSens-1 repair variant (b), scalar subtraction: S = {z8.mean():+.5f} - {b.mean():+.5f} = "
  f"{z8.mean() - b.mean():+.5f}, propagated 95% [{lo:+.5f}, {hi:+.5f}] "
  f"(z-bootstrap CI shifted by b-bar, widened by 1.96 x SE(b-bar) = {1.96 * se_b:.5f})")
P(f"       the two variants have the SAME point estimate by linearity "
  f"(mean(z) - mean(b) = mean(z - b) = {C.PRI['S']:+.5f}); they differ only in the interval and "
  f"the null.  Variant (a)'s interval is [{C.PRI['ci95'][0]:+.5f}, {C.PRI['ci95'][1]:+.5f}]")
out["Sens1_variant_b"] = dict(S=float(z8.mean() - b.mean()), ci95=[lo, hi], b_bar=float(b.mean()),
                              se_b_bar=se_b, variant_a_S=C.PRI["S"], variant_a_ci=C.PRI["ci95"],
                              b_sd=float(b.std(ddof=1)),
                              b_range=[float(b.min()), float(b.max())])

# ============================================================ Sens-2  centring over all 33
rec33 = {c: [d for d in TARG if d != c] for c in TARG}
rec33.update({c: [d for d in NTCC if d != c] for c in NTCC})
_, OBJ33 = C.build(ALL, hard_windows=True, global_excl=True, recips=rec33)
B33 = C.verdict_block(TARG, OBJ33, "Sens-2 all-33 centring", X.SEED309 + 43)
P(f"Sens-2 centring over all 33 other targeting constructs (cross-split): S = {B33['S']:+.5f} "
  f"[{B33['ci95'][0]:+.5f}, {B33['ci95'][1]:+.5f}], p = {B33['p_joint']:.6g}, "
  f"max leave-one-out p = {B33['loo_max_p']:.6g}")
out["Sens2_all33"] = {k: v for k, v in B33.items() if k not in ("Sn", "loo", "z")}

# ============================================================ Sens-3  global exclusion off
_, OBJnx = C.build(ALL, hard_windows=True, global_excl=False)
Bnx = C.verdict_block(TARG, OBJnx, "Sens-3 global exclusion off", X.SEED309 + 45)
P(f"Sens-3 the global L*-candidate exclusion switched OFF: S = {Bnx['S']:+.5f} "
  f"[{Bnx['ci95'][0]:+.5f}, {Bnx['ci95'][1]:+.5f}], p = {Bnx['p_joint']:.6g} "
  f"(primary keeps it ON: {C.PRI['S']:+.5f}; difference {C.PRI['S'] - Bnx['S']:+.5f})")
out["Sens3_no_global_exclusion"] = {k: v for k, v in Bnx.items() if k not in ("Sn", "loo", "z")}

# ============================================================ Sens-4  each repair alone
B6 = C.verdict_block(TARG, OBJ9, "Sens-4a KA-6 repair only (uncentred)", X.SEED309 + 47,
                     centred=False)
B4 = C.verdict_block(TARG, OBJ8, "Sens-4b KA-4 repair only (C8 pools, centred)", X.SEED309 + 49)
P(f"Sens-4a KA-6 repair alone (150-cell pools, UNCENTRED): S = {B6['S']:+.5f} "
  f"[{B6['ci95'][0]:+.5f}, {B6['ci95'][1]:+.5f}], p = {B6['p_joint']:.6g}   "
  f"(C8's post-hoc P2 read +0.78813)")
P(f"Sens-4b KA-4 repair alone (C8's 50-cell pools, CENTRED): S = {B4['S']:+.5f} "
  f"[{B4['ci95'][0]:+.5f}, {B4['ci95'][1]:+.5f}], p = {B4['p_joint']:.6g}   "
  f"(C8's post-hoc P5 read +0.71399)")
P(f"Sens-4  both repairs (the C9 primary): S = {C.PRI['S']:+.5f} "
  f"[{C.PRI['ci95'][0]:+.5f}, {C.PRI['ci95'][1]:+.5f}], p = {C.PRI['p_joint']:.6g}")
out["Sens4a_ka6_only"] = {k: v for k, v in B6.items() if k not in ("Sn", "loo", "z")}
out["Sens4b_ka4_only"] = {k: v for k, v in B4.items() if k not in ("Sn", "loo", "z")}

out["peak_rss_gib"] = round(X.peak_gib(), 3)
out["seconds"] = round(time.time() - t0, 1)
(X.HERE / "sens309.json").write_text(json.dumps(out, indent=1, default=float))
(X.HERE / "sens309.log").write_text("\n".join(log) + "\n")
P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")

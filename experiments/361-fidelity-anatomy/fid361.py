#!/usr/bin/env python3
"""361 — fid anatomy + CPM-block sign counterfactuals. Exact arithmetic on saved arrays only.

Reads (read-only): 323 work stats-cap15.npz, reference.npz, effects.npz, per-target-cap15.csv,
metrics-cap15.json, h1-anchors.json. Writes only into experiments/361-fidelity-anatomy/work/.
"""
import json
import os
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
BASE = "~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work"
ALPHA = 0.05
DRAWS, SEED = 10000, 20260909

os.makedirs(WORK, exist_ok=True)
t0 = time.monotonic()

# ---------- load ----------
s = np.load(os.path.join(BASE, "stats-cap15.npz"))
targets = s["targets"].astype(str)
pp, qp, lp = s["p"], s["q"], s["lfc"]
s.close()
r = np.load(os.path.join(BASE, "reference.npz"), allow_pickle=False)
genes_full = r["genes"].astype(str)
rtargets = r["targets"].astype(str)
keep = r["keep"]
r_q, r_lfc = r["q"], r["lfc"]
ref_lin = r["ref_lin"]
r.close()
ugenes = genes_full[keep]
assert np.array_equal(targets, rtargets), "target order drift"
n, g = len(targets), len(ugenes)
assert pp.shape == (n, g) and r_q.shape == (n, g), (pp.shape, r_q.shape)
ze = np.load(os.path.join(BASE, "effects.npz"), allow_pickle=False)
egenes = ze["genes"].astype(str)
etargets = ze["targets"].astype(str)
cap15 = ze["cap15"]
ze.close()
print(f"targets={n} universe={g} full={len(genes_full)}", flush=True)

# design matrix on universe by name (2 missing -> 0, per 347)
pos = {gg: i for i, gg in enumerate(egenes)}
epos = {t: i for i, t in enumerate(etargets)}
D = np.zeros((n, g))
missing = [gg for gg in ugenes if gg not in pos]
for j, gg in enumerate(ugenes):
    if gg in pos:
        D[:, j] = cap15[[epos[t] for t in targets], pos[gg]]
print(f"design coverage {g - len(missing)}/{g} missing={missing}", flush=True)
assert all(t in ugenes for t in targets), "own target outside universe!"

own = np.equal.outer(targets, ugenes)  # (n,g) True at own target gene
adjud = np.isfinite(r_lfc) & (r_lfc != 0)
commit = np.isfinite(lp) & (lp != 0)
print(f"adjudicable frac={adjud.mean():.4f} commit frac={commit.mean():.4f}", flush=True)

ref_sig = r_q < ALPHA
pred_sig = np.isfinite(qp) & (qp < ALPHA)

# ---------- metric implementations ----------
def fid_per_target(qp_, lp_):
    ps = np.isfinite(qp_) & (qp_ < ALPHA)
    scored = ps & adjud & ~own
    n_pred = scored.sum(1).astype(float)
    n_real = (ref_sig & ~own).sum(1).astype(float)
    match = scored & np.isfinite(lp_) & (lp_ != 0) & (np.sign(lp_) == np.sign(r_lfc))
    k = match.sum(1).astype(float)
    den = np.maximum(n_pred, n_real)
    with np.errstate(invalid="ignore", divide="ignore"):
        F = np.where(den > 0, k / den, np.nan)
    return F, n_pred, n_real, k

def nmae_per_target(lp_):
    gate = ref_sig & np.isfinite(r_lfc) & ~own
    cnt = gate.sum(1)
    num = np.where(gate, np.abs(lp_ - r_lfc), 0.0).sum(1)
    den = np.where(gate, np.abs(r_lfc), 0.0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where((cnt >= 10) & np.isfinite(den) & (den != 0), num / den, np.nan)
    return v, cnt

def reach_per_target(lp_, qp_, pp_):
    pool = ref_sig & ~own
    n_conf = pool.sum(1).astype(float)
    R = np.full(n, np.nan)
    # ordering keys per target: sig-first, q asc, p asc, |lfc| desc, gene asc
    qm = np.where(np.isfinite(qp_), qp_, np.inf)
    pm = np.where(np.isfinite(pp_), pp_, np.inf)
    ala = np.abs(lp_)
    sig = np.isfinite(qp_) & (qp_ < ALPHA)
    feat = ugenes
    for i in range(n):
        if n_conf[i] == 0:
            continue
        idx = np.lexsort((feat, -np.where(np.isfinite(ala[i]), ala[i], -np.inf),
                          pm[i], qm[i], ~sig[i]))
        # walk adjudicable pool genes in that order
        m = ((np.isfinite(lp_[i]) & (lp_[i] != 0)
              & (np.sign(lp_[i]) == np.sign(r_lfc[i]))) & pool[i])[idx]
        a = (adjud[i] & pool[i])[idx]
        ma = m[a]
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / np.arange(1, ma.size + 1)
        ok = np.flatnonzero(pur >= 0.9)
        kstar = int(ok.max() + 1) if ok.size else 0
        R[i] = kstar / n_conf[i]
    return R, n_conf

# ---------- known answers ----------
rec = pl.read_csv(os.path.join(BASE, "per-target-cap15.csv"))
def recorded(metric):
    d = rec.filter(pl.col("metric") == metric).sort("perturbation")
    return dict(zip(d["perturbation"].to_list(), d["value"].to_list()))

F0, n_pred0, n_real0, k0 = fid_per_target(qp, lp)
N0, cnt0 = nmae_per_target(lp)
R0, nc0 = reach_per_target(lp, qp, pp)
f_rec = recorded("de_wilcoxon_direction_fidelity_yield_raw")
n_rec = recorded("de_wilcoxon_lfc_nmae")
r_rec = recorded("de_wilcoxon_direction_reach_raw")
df = max(abs(float(F0[i]) - f_rec[t]) for i, t in enumerate(targets))
print(f"fid: repro mean={np.nanmean(F0):.8f} recorded=0.49755452 maxdiff={df:.3e}", flush=True)
# nmae recorded only for 58 eligible
dn = max(abs(float(N0[i]) - n_rec[t]) for i, t in enumerate(targets) if t in n_rec)
print(f"nmae: repro mean={np.nanmean(N0):.8f} recorded=0.95083646 maxdiff={dn:.3e} eligible={(cnt0>=10).sum()}", flush=True)
rr = [(t, float(R0[i]), r_rec[t]) for i, t in enumerate(targets) if t in r_rec]
dr = max(abs(a - b) for _, a, b in rr if np.isfinite(a) and np.isfinite(b))
nanmatch = sum(1 for _, a, b in rr if (not np.isfinite(a)) and (not np.isfinite(b)))
print(f"reach: repro mean={np.nanmean(R0):.8f} recorded=0.25434726 maxdiff={dr:.3e} nanmatch={nanmatch}/85", flush=True)
assert df <= 1e-9, "fid known-answer FAILED"
assert dn <= 1e-9, "nmae known-answer FAILED"
assert dr <= 1e-9 and nanmatch == 3, "reach known-answer FAILED"
print("KNOWN-ANSWER GATE: PASS (fid, nmae, reach all <=1e-9)", flush=True)

# ---------- decomposition ----------
calls = pred_sig & adjud & ~own  # == pred_sig & ~own here
is_refsig = ref_sig
n_calls = calls.sum(1).astype(float)
share_refsig = np.where(n_calls > 0, (calls & is_refsig).sum(1) / n_calls, np.nan)
match0 = calls & (np.sign(lp) == np.sign(r_lfc))
acc_all = np.where(n_calls > 0, match0.sum(1) / n_calls, np.nan)
cs = calls & is_refsig
cn = calls & ~is_refsig
acc_s = np.where(cs.sum(1) > 0, (match0 & is_refsig).sum(1) / np.maximum(cs.sum(1), 1), np.nan)
acc_n = np.where(cn.sum(1) > 0, (match0 & ~is_refsig).sum(1) / np.maximum(cn.sum(1), 1), np.nan)

# head / block partition
head = np.zeros((n, g), bool)
block = np.zeros((n, g), bool)
rank_ref = np.argsort(-ref_lin)  # desc basal CPM
for i in range(n):
    h = np.argsort(-np.abs(D[i]), kind="stable")[:400]
    head[i, h] = True
    rest = [j for j in rank_ref if not head[i, j]][:2200]
    block[i, rest] = True
other = ~(head | block)
print(f"head calls/target med={np.median((calls&head).sum(1)):.0f} "
      f"block calls med={np.median((calls&block).sum(1)):.0f} "
      f"other-called med={np.median((calls&other).sum(1)):.0f}", flush=True)

def part_stats(mask, name):
    c = calls & mask
    mm = match0 & mask
    cp = c.sum()
    acc = mm.sum() / cp if cp else np.nan
    # contribution to mean F: sum over targets of k_part/max(n_pred,n_real)
    den = np.maximum(n_pred0, n_real0)
    contrib = float(np.sum(mm.sum(1) / np.where(den > 0, den, np.nan)))
    contrib = contrib / n  # mean over 85
    print(f"{name}: pooled calls={cp} acc={acc:.4f} contrib_to_meanF={contrib:.4f}", flush=True)
    return {"calls": int(cp), "acc": float(acc), "contrib": float(contrib)}

parts = {}
for m_, nm in ((head, "head400"), (block, "cpmblock2200"), (other, "other")):
    parts[nm] = part_stats(m_, nm)

# reference sign tilt
neg = r_lfc < 0
tilt_block_calls = float((neg & calls & block).sum() / (calls & block).sum())
tilt_all = float(neg[adjud].mean())
tilt_refsig = float(neg[ref_sig & ~own].mean()) if (ref_sig & ~own).sum() else np.nan
tilt_block_per_target = np.where((calls & block).sum(1) > 0,
                                 (neg & calls & block).sum(1) / np.maximum((calls & block).sum(1), 1),
                                 np.nan)
print(f"tilt frac-down: block-calls pooled={tilt_block_calls:.4f} "
      f"per-target mean={np.nanmean(tilt_block_per_target):.4f} all-adjud={tilt_all:.4f} "
      f"refsig={tilt_refsig:.4f}", flush=True)

decomp = {
    "fid_mean": float(np.nanmean(F0)),
    "n_pred_med": float(np.median(n_pred0)), "n_real_med": float(np.median(n_real0)),
    "n_real_mean": float(np.mean(n_real0)), "k_med": float(np.median(k0)),
    "share_calls_refsig_pooled": float((calls & is_refsig).sum() / calls.sum()),
    "share_calls_refsig_per_target_mean": float(np.nanmean(share_refsig)),
    "kacc_pooled_all": float(match0.sum() / calls.sum()),
    "kacc_per_target_mean": float(np.nanmean(acc_all)),
    "kacc_on_refsig_pooled": float((match0 & is_refsig).sum() / max((calls & is_refsig).sum(), 1)),
    "kacc_on_nonsig_pooled": float((match0 & ~is_refsig).sum() / max((calls & ~is_refsig).sum(), 1)),
    "kacc_on_refsig_per_target_mean": float(np.nanmean(acc_s)),
    "kacc_on_nonsig_per_target_mean": float(np.nanmean(acc_n)),
    "parts": parts,
    "tilt_block_calls_pooled": tilt_block_calls,
    "tilt_block_calls_per_target_mean": float(np.nanmean(tilt_block_per_target)),
    "tilt_all_adjudicable": tilt_all,
    "tilt_refsignificant": tilt_refsig,
}
with open(os.path.join(WORK, "decomposition.json"), "w") as f:
    json.dump(decomp, f, indent=2)

# ---------- counterfactuals ----------
# LOO common vectors
ref_common_sign = np.zeros((n, g))
des_common_sign = np.zeros((n, g))
ref_zero_ties = 0
des_zero_ties = 0
ref_tot = r_lfc.sum(0)
des_tot = D.sum(0)
for i in range(n):
    rc = (ref_tot - r_lfc[i]) / (n - 1)
    dc = (des_tot - D[i]) / (n - 1)
    ref_common_sign[i] = np.sign(rc)
    des_common_sign[i] = np.sign(dc)
ref_zero_ties = int((((ref_common_sign == 0) & block).sum()))
des_zero_ties = int((((des_common_sign == 0) & block).sum()))
print(f"LOO ties on block genes: ref {ref_zero_ties}, design {des_zero_ties}", flush=True)

top_half = np.zeros((n, g), bool)  # top half of BLOCK by basal CPM
for i in range(n):
    bgenes = np.flatnonzero(block[i])
    order = bgenes[np.argsort(-ref_lin[bgenes], kind="stable")]
    top_half[i, order[: len(order) // 2]] = True

abs_lp = np.abs(lp)
scenarios = {}
# e: current
scenarios["e_current"] = lp.copy()
# a: all-down on block
la = lp.copy()
la[block] = -abs_lp[block]
scenarios["a_alldown"] = la
# b: LOO ref-common sign on block (ties keep current)
lb = lp.copy()
use = block & (ref_common_sign != 0)
lb[use] = ref_common_sign[use] * abs_lp[use]
scenarios["b_looref"] = lb
# c: LOO design-common sign on block
lc = lp.copy()
use = block & (des_common_sign != 0)
lc[use] = des_common_sign[use] * abs_lp[use]
scenarios["c_loodesign"] = lc
# d: down on top-half block, design sign elsewhere
ld = lp.copy()
ld[block & top_half] = -abs_lp[block & top_half]
scenarios["d_tophalfdown"] = ld

res = {}
for name, L in scenarios.items():
    F, _, _, _ = fid_per_target(qp, L)
    N, _ = nmae_per_target(L)
    R, _ = reach_per_target(L, qp, pp)
    res[name] = {"fid": F, "nmae": N, "reach": R}
    print(f"{name}: fid={np.nanmean(F):.6f} nmae={np.nanmean(N):.6f} "
          f"reach={np.nanmean(R):.6f} (elig nmae={np.isfinite(N).sum()} reach={np.isfinite(R).sum()})",
          flush=True)

# flip-rate diagnostics: what fraction of block calls each rule flips
for name in ("a_alldown", "b_looref", "c_loodesign", "d_tophalfdown"):
    L = scenarios[name]
    flipped = ((np.sign(L) != np.sign(lp)) & calls & block).sum()
    tot = (calls & block).sum()
    print(f"{name}: flips {flipped}/{tot} = {flipped/tot:.3f} of block calls", flush=True)

# ---------- paired bootstrap vs current ----------
rng = np.random.default_rng(SEED)
draws = rng.integers(0, n, (DRAWS, n))
counts = np.zeros((DRAWS, n), dtype=np.int64)
np.add.at(counts, (np.repeat(np.arange(DRAWS), n), draws.ravel()), 1)
counts = np.vstack([np.ones(n, dtype=np.int64), counts])  # row0 = point

def boot(values):
    v = np.asarray(values, float)
    valid = np.isfinite(v)
    den = counts[:, valid].sum(1)
    m = (counts[:, valid] @ v[valid]) / den
    assert abs(m[0] - np.nanmean(v)) < 1e-12
    return m

def ci(v):
    v = np.asarray(v, float)
    f = v[np.isfinite(v)]
    return {"lo": float(np.percentile(f, 2.5)), "hi": float(np.percentile(f, 97.5))}

anchors = json.load(open(os.path.join(BASE, "h1-anchors.json")))
B, Rp = anchors["baseline"], anchors["replicate"]
table = {}
for name in scenarios:
    if name == "e_current":
        continue
    entry = {}
    for m in ("fid", "nmae", "reach"):
        a = boot(res[name][m])
        b = boot(res["e_current"][m])
        d = a - b
        span = Rp[m] - B[m]
        dsp = {"point_raw": float(d[0]), **ci(d[1:]),
               "point_scaled": float(d[0] / span)}
        entry[m] = dsp
    # H1-local overall: mean of scaled member deltas / 6 (mse/pds/jac carried at 0)
    entry["overall_H1local"] = float(sum(entry[m]["point_scaled"] for m in ("fid", "nmae", "reach")) / 6)
    table[name] = entry
    print(name, json.dumps({m: {k: round(v, 6) for k, v in entry[m].items()} for m in ("fid", "nmae", "reach")},
                           indent=1), flush=True)
    print(f"  overall_H1local Δ = {entry['overall_H1local']:+.6f}", flush=True)

with open(os.path.join(WORK, "counterfactuals.json"), "w") as f:
    json.dump({"point_means": {k: {"fid": float(np.nanmean(v["fid"])),
                                   "nmae": float(np.nanmean(v["nmae"])),
                                   "reach": float(np.nanmean(v["reach"]))} for k, v in res.items()},
               "contrasts_vs_current": table,
               "anchors": {"baseline": {m: B[m] for m in ("fid", "nmae", "reach")},
                           "replicate": {m: Rp[m] for m in ("fid", "nmae", "reach")}},
               "loo_ties_on_block": {"ref": ref_zero_ties, "design": des_zero_ties},
               "draws": DRAWS, "seed": SEED}, f, indent=2)

# per-target detail for the record
pl.DataFrame({"target": targets,
              "fid_current": res["e_current"]["fid"], "fid_alldown": res["a_alldown"]["fid"],
              "fid_looref": res["b_looref"]["fid"], "fid_loodesign": res["c_loodesign"]["fid"],
              "fid_tophalf": res["d_tophalfdown"]["fid"],
              "nmae_current": res["e_current"]["nmae"], "nmae_alldown": res["a_alldown"]["nmae"],
              "nmae_looref": res["b_looref"]["nmae"], "nmae_loodesign": res["c_loodesign"]["nmae"],
              "nmae_tophalf": res["d_tophalfdown"]["nmae"],
              "reach_current": res["e_current"]["reach"], "reach_alldown": res["a_alldown"]["reach"],
              "reach_looref": res["b_looref"]["reach"], "reach_loodesign": res["c_loodesign"]["reach"],
              "reach_tophalf": res["d_tophalfdown"]["reach"],
              "n_pred": n_pred0, "n_real": n_real0, "k": k0,
              "tilt_block": tilt_block_per_target}).write_csv(os.path.join(WORK, "per-target-361.csv"))

import resource
print(f"done in {time.monotonic()-t0:.1f}s peakRSS={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f}GB", flush=True)

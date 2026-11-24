#!/usr/bin/env python3
"""364 — amplitude break-even vs gate sign accuracy. Exact arithmetic on saved
arrays only. Reads 323/317 work files (read-only); writes only into this dir's work/.
"""
import json, os, time
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import polars as pl

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)
BASE = "~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work"
S317 = "~/Code/virtual-cell-challenge-2026/experiments/317-operating-point/work/sweep317.json"
ALPHA = 0.05
t0 = time.monotonic()

def mem():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

# ---------- load ----------
s = np.load(os.path.join(BASE, "stats-cap15.npz"))
targets = s["targets"].astype(str)
pp, qp, lp0 = s["p"].astype(float), s["q"].astype(float), s["lfc"].astype(float)
s.close()
r = np.load(os.path.join(BASE, "reference.npz"), allow_pickle=False)
genes_full = r["genes"].astype(str); rtargets = r["targets"].astype(str)
keep = r["keep"].astype(bool)
r_q, r_lfc = r["q"].astype(float), r["lfc"].astype(float)
r.close()
ugenes = genes_full[keep]
assert np.array_equal(targets, rtargets)
n, g = len(targets), len(ugenes)
assert pp.shape == (n, g) and r_q.shape == (n, g), (pp.shape, r_q.shape, n, g)

own = np.equal.outer(targets, ugenes)
adjud = np.isfinite(r_lfc) & (r_lfc != 0)
ref_sig = (r_q < ALPHA) & adjud
pred_sig = np.isfinite(qp) & (qp < ALPHA)
gate = ref_sig & ~own
elig = gate.sum(1) >= 10
print(f"targets={n} universe={g} eligible={elig.sum()} (expect 58)", flush=True)
print(f"rss={mem():.2f}GB", flush=True)

# ---------- metric functionals (361 conventions) ----------
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
    return F

def nmae_per_target(lp_):
    gt = ref_sig & ~own
    cnt = gt.sum(1)
    num = np.where(gt, np.abs(lp_ - r_lfc), 0.0).sum(1)
    den = np.where(gt, np.abs(r_lfc), 0.0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = np.where((cnt >= 10) & np.isfinite(den) & (den != 0), num / den, np.nan)
    return v

def reach_per_target(lp_, qp_, pp_):
    pool = ref_sig & ~own
    n_conf = pool.sum(1).astype(float)
    R = np.full(n, np.nan)
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
        m = ((np.isfinite(lp_[i]) & (lp_[i] != 0)
              & (np.sign(lp_[i]) == np.sign(r_lfc[i]))) & pool[i])[idx]
        a = (adjud[i] & pool[i])[idx]
        ma = m[a]
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / np.arange(1, ma.size + 1)
        ok = np.flatnonzero(pur >= 0.9)
        R[i] = (int(ok.max() + 1) / n_conf[i]) if ok.size else 0.0
    return R

# ---------- 1. known answers ----------
F0 = fid_per_target(qp, lp0); N0 = nmae_per_target(lp0); R0 = reach_per_target(lp0, qp, pp)
rec = pl.read_csv(os.path.join(BASE, "per-target-cap15.csv"))
def recorded(metric):
    d = rec.filter(pl.col("metric") == metric).sort("perturbation")
    return dict(zip(d["perturbation"].to_list(), d["value"].to_list()))
f_rec = recorded("de_wilcoxon_direction_fidelity_yield_raw")
n_rec = recorded("de_wilcoxon_lfc_nmae")
r_rec = recorded("de_wilcoxon_direction_reach_raw")
df = max(abs(float(F0[i]) - f_rec[t]) for i, t in enumerate(targets))
dn = max(abs(float(N0[i]) - n_rec[t]) for i, t in enumerate(targets) if t in n_rec)
rr = [(t, float(R0[i]), r_rec[t]) for i, t in enumerate(targets) if t in r_rec]
dr = max(abs(a - b) for _, a, b in rr if np.isfinite(a) and np.isfinite(b))
nanmatch = sum(1 for _, a, b in rr if (not np.isfinite(a)) and (not np.isfinite(b)))
print(f"KNOWN fid mean={np.nanmean(F0):.8f} (want 0.49755452) maxdiff={df:.2e}", flush=True)
print(f"KNOWN nmae mean={np.nanmean(N0):.8f} (want 0.95083646) maxdiff={dn:.2e}", flush=True)
print(f"KNOWN reach mean={np.nanmean(R0):.8f} (want 0.25434726) maxdiff={dr:.2e} nan={nanmatch}/85", flush=True)
ka_pass = (df <= 1e-9 and dn <= 1e-9 and dr <= 1e-9 and nanmatch == 3
           and abs(np.nanmean(F0) - 0.49755452) < 5e-9
           and abs(np.nanmean(N0) - 0.95083646) < 5e-9
           and abs(np.nanmean(R0) - 0.25434726) < 5e-9)
print(f"KNOWN-ANSWER GATE: {'PASS' if ka_pass else 'FAIL'}", flush=True)
if not ka_pass:
    raise SystemExit("known-answer FAILED — stopping per rule 3")

# delivered sign accuracy (pooled + mean-over-targets, |l|-weighted, eligible only)
w = np.where(gate & elig[:, None], np.abs(r_lfc), 0.0)
Dtot = w.sum()
correct0 = (np.sign(lp0) == np.sign(r_lfc)) & (lp0 != 0)
p_pooled_0 = float((w * correct0).sum() / Dtot)
pw_t = np.array([float((w[i] * correct0[i]).sum() / w[i].sum()) for i in range(n) if elig[i]])
p_mean_0 = float(pw_t.mean())
m_ratio = float(np.abs(np.where(gate & elig[:, None], lp0, 0.0)).sum() / Dtot)
print(f"delivered pooled w-acc={p_pooled_0:.6f} mean-over-targets={p_mean_0:.6f} mag-ratio={m_ratio:.4f}", flush=True)
print(f"eligible gate pairs={int((gate & elig[:, None]).sum())}", flush=True)

# ---------- 2+3. sign-flip model x amplitude grid ----------
P_GRID = [0.60, 0.62, 0.65, 0.70, 0.75]
A_GRID = [1.0, 1.2, 1.45, 2.0]
DRAWS = 20
abs_lp0 = np.abs(lp0)

# wrong-sign pool (eligible targets only): flat index lists
wrong_pairs = [(i, j) for i in np.flatnonzero(elig) for j in np.flatnonzero(gate[i] & ~correct0[i])]
wrong_mass = np.array([abs(r_lfc[i, j]) for (i, j) in wrong_pairs])
order_mass = wrong_mass  # fixed masses; shuffle order per draw
need_base_correct = (w * correct0).sum()
print(f"wrong-sign pairs={len(wrong_pairs)} mass={wrong_mass.sum():.1f} total D={Dtot:.1f}", flush=True)

results = {}   # results[p][a] = dict(mean nmae/fid/reach over draws + sd + achieved p)
for p in P_GRID:
    need = p * Dtot - need_base_correct
    print(f"--- p={p} need flip mass={need:.1f} ({100*need/Dtot:.2f}% of D)", flush=True)
    if need <= 0:
        print("  already above target — using delivered signs (no flips)", flush=True)
    per_draw = []
    for d in range(DRAWS):
        rng = np.random.default_rng(1000 + d)
        perm = rng.permutation(len(wrong_pairs))
        s_new = np.sign(lp0).copy()
        if need > 0:
            cum = np.cumsum(wrong_mass[perm])
            k = int(np.searchsorted(cum, need, side="left")) + 1
            k = min(k, len(perm))
            for t in perm[:k]:
                i, j = wrong_pairs[t]
                s_new[i, j] = np.sign(r_lfc[i, j])
        # achieved pooled accuracy
        p_ach = float((w * (s_new == np.sign(r_lfc))).sum() / Dtot)
        row = {"p_ach": p_ach}
        for a in A_GRID:
            L = lp0.copy()
            gm = gate  # scale on ALL gate genes (eligible or not — ineligible contribute nothing to nmae mean but keep fid exact)
            L[gm] = a * s_new[gm] * abs_lp0[gm]
            N = nmae_per_target(L); F = fid_per_target(qp, L); R = reach_per_target(L, qp, pp)
            row[a] = (float(np.nanmean(N[elig])), float(np.nanmean(F)), float(np.nanmean(R)))
        per_draw.append(row)
    agg = {"p_ach_mean": float(np.mean([r["p_ach"] for r in per_draw]))}
    for a in A_GRID:
        for m, k in (("nmae", 0), ("fid", 1), ("reach", 2)):
            v = np.array([r[a][k] for r in per_draw])
            agg[f"{m}_a{a}_mean"] = float(v.mean())
            agg[f"{m}_a{a}_sd"] = float(v.std(ddof=1)) if len(v) > 1 else 0.0
    results[str(p)] = agg
    print(f"  p_ach={agg['p_ach_mean']:.5f} " +
          " ".join(f"nmae@{a}={agg[f'nmae_a{a}_mean']:.5f}±{agg[f'nmae_a{a}_sd']:.5f}" for a in A_GRID), flush=True)

# invariance check: fid/reach across a within each (p, draw) — recompute max |diff| on draw 0
print("--- invariance across a (draw 0, p=0.65 flips) ---", flush=True)
rng = np.random.default_rng(1000)
perm = rng.permutation(len(wrong_pairs))
need = 0.65 * Dtot - need_base_correct
s_chk = np.sign(lp0).copy()
cum = np.cumsum(wrong_mass[perm]); k = min(int(np.searchsorted(cum, need, side="left")) + 1, len(perm))
for t in perm[:k]:
    i, j = wrong_pairs[t]
    s_chk[i, j] = np.sign(r_lfc[i, j])
Ls = {}
for a in A_GRID:
    L = lp0.copy(); L[gate] = a * s_chk[gate] * abs_lp0[gate]; Ls[a] = L
Fx = [fid_per_target(qp, Ls[a]) for a in A_GRID]
Rx = [reach_per_target(Ls[a], qp, pp) for a in A_GRID]
print(f"  max|fid(a)-fid(1)| = {max(np.nanmax(np.abs(Fx[k]-Fx[0])) for k in range(4)):.3e}", flush=True)
print(f"  max|reach(a)-reach(1)| = {max(np.nanmax(np.abs(Rx[k]-Rx[0])) for k in range(4)):.3e}", flush=True)

json.dump({"p_grid": P_GRID, "a_grid": A_GRID, "draws": DRAWS,
           "p_pooled_0": p_pooled_0, "p_mean_0": p_mean_0, "m_ratio": m_ratio,
           "need_mass": {str(p): float(p*Dtot-need_base_correct) for p in P_GRID},
           "table": results,
           "seconds": time.monotonic()-t0, "peakRSS_GB": mem()},
          open(os.path.join(WORK, "breakeven.json"), "w"), indent=1)
print(f"done in {time.monotonic()-t0:.1f}s peakRSS={mem():.2f}GB", flush=True)

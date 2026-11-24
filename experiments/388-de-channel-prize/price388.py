#!/usr/bin/env python3
"""388 pricing — every estimator as a signed (85 x 10779) field; gate/nmae/reach proxies.

E0 delivered lfc | E1 designed cap15 | E2 un-whitened cons (172, 2 donors:
k562_gwps+orion_hek293t) | E3 flip (173 analogue: E1 mags, E2 signs where they
disagree) | E4 N/A (donorstats carry gene-level stats only, no signed rows) |
E5 deranged E2 | E6 random signs x |E0|.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
import json
import numpy as np
import polars as pl
from pathlib import Path

REPO = Path("~/Code/virtual-cell-challenge-2026")
W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/388-de-channel-prize/work"
SEED_BOOT = 20260913
DRAWS = 10000
GRID343 = np.concatenate([np.linspace(0.25, 1, 4), np.linspace(1.25, 6, 20)])
GRIDWIDE = np.concatenate([np.linspace(0.05, 1, 8), np.linspace(1.25, 24, 48)])

# ---------- load ----------
s = np.load(W / "stats-cap15.npz")
r = np.load(W / "reference.npz")
e = np.load(W / "effects.npz")
f172 = np.load(REPO / "experiments/172-sign-confidence-key/work/h1test-feats.npz")
genes = r["genes"].astype(str)[r["keep"].astype(bool)]
targets = s["targets"].astype(str)
assert (targets == r["targets"].astype(str)).all()
ref = pl.read_parquet(W / "reference_de.parquet")
tcol = ref["target"].to_numpy().astype(str)
gcol = ref["feature"].to_numpy().astype(str)
gi = {g: i for i, g in enumerate(genes)}
E0 = np.where(np.isfinite(np.asarray(s["lfc"], dtype=float)),
              np.asarray(s["lfc"], dtype=float), 0.)
QHAT = np.asarray(s["q"], dtype=float)
PHAT = np.asarray(s["p"], dtype=float)

eg = {g: i for i, g in enumerate(e["genes"].astype(str))}
col = np.array([eg.get(g, -1) for g in genes])
miss = col < 0
E1 = np.asarray(e["cap15"], dtype=float)[:, np.where(miss, 0, col)]
E1[:, miss] = 0.0
print(f"E1: {int(miss.sum())} universe genes missing from effects axis -> 0")

E2 = np.asarray(f172["cons"], dtype=float)
assert E2.shape == E0.shape, (E2.shape, E0.shape)
print(f"E2: nonzero fraction {float((E2 != 0).mean()):.3f}")

E3 = E1.copy()
dis = (np.sign(E1) != np.sign(E2)) & (E2 != 0)
E3[dis] = np.sign(E2[dis]) * np.abs(E1[dis])
print(f"E3: flipped {float(dis.mean()):.4f} of entries")

rng = np.random.default_rng(388)
perm = rng.permutation(len(targets))
while (perm == np.arange(len(targets))).any():
    perm = rng.permutation(len(targets))
E5 = E2[perm]
E6 = np.sign(rng.random(E0.shape) - 0.5) * np.abs(E0)
E6[E6 == 0] = 1.0
EST = {"E0": E0, "E1": E1, "E2": E2, "E3": E3, "E5": E5, "E6": E6}

# ---------- gate (343 convention) ----------
elig, GIDX, LREF = [], [], []
for k, t in enumerate(targets):
    m = (tcol == t)
    g = gcol[m]
    q = ref["p_adj"].to_numpy().astype(float)[m]
    l = ref["log2_fold_change"].to_numpy().astype(float)[m]
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10:
        continue
    elig.append(k)
    GIDX.append(np.array([gi[x] for x in g[gate]]))
    LREF.append(l[gate])
elig = np.array(elig)
print(f"eligible targets: {len(elig)}")
np.savez_compressed(OUT / "est388.npz",
                    **{k: v.astype(np.float32) for k, v in EST.items()},
                    elig=elig)

# ---------- 1+2. sign accuracy + analytic nmae ----------
res = {}
for name, X in EST.items():
    pw_g, pw_c, up_c, mass_c = [], [], [], []
    nm_own, nm_e0m = [], []
    per_t_scalar, per_t_wide = [], []
    for e_i, k in enumerate(elig):
        idx = GIDX[e_i]
        l = LREF[e_i]
        D = np.abs(l).sum()
        x = X[k, idx]
        sd = np.sign(x) == np.sign(l)
        pw_g.append((np.abs(l) * sd).sum() / D)
        called = QHAT[k, idx] < .05
        mass_c.append(np.abs(l[called]).sum() / D)
        if called.any():
            pw_c.append((np.abs(l[called]) * sd[called]).sum() / np.abs(l[called]).sum())
            up_c.append(sd[called].mean())
        else:
            pw_c.append(np.nan)
            up_c.append(np.nan)
        nm_own.append(np.abs(x - l).sum() / D)
        xe0 = np.sign(x) * np.abs(E0[k, idx])
        nm_e0m.append(np.abs(xe0 - l).sum() / D)
        vals = [np.abs(a * x - l).sum() / D for a in GRID343]
        per_t_scalar.append(min(vals))
        valsW = [np.abs(a * x - l).sum() / D for a in GRIDWIDE]
        per_t_wide.append(min(valsW))
    res[name] = dict(pw_gate=np.array(pw_g), pw_call=np.array(pw_c),
                     up_call=np.array(up_c), mass_called=np.array(mass_c),
                     nmae_own=np.array(nm_own), nmae_e0mag=np.array(nm_e0m),
                     nmae_pt=np.array(per_t_scalar), nmae_ptwide=np.array(per_t_wide))

# global-scalar + decile maps (fit pooled over gate pairs, score mean over targets)
pool_x = {n: np.concatenate([X[k, GIDX[i]] for i, k in enumerate(elig)]) for n, X in EST.items()}
pool_l = np.concatenate(LREF)
for name in EST:
    x = pool_x[name]
    a_star = min(GRID343, key=lambda a: np.abs(a * x - pool_l).sum())
    nm = []
    for e_i, k in enumerate(elig):
        l = LREF[e_i]
        nm.append(np.abs(a_star * EST[name][k, GIDX[e_i]] - l).sum() / np.abs(l).sum())
    res[name]["nmae_global"] = np.array(nm)
    res[name]["a_global"] = float(a_star)
    ax = np.abs(x)
    qs = np.quantile(ax, np.linspace(0, 1, 11))
    dec = np.clip(np.searchsorted(qs[1:-1], ax), 0, 9)
    amap = np.zeros(10)
    for d in range(10):
        sel = dec == d
        if sel.sum() == 0 or (x[sel] == 0).all():
            amap[d] = 0.0
        else:
            amap[d] = min(GRID343, key=lambda a: np.abs(a * x[sel] - pool_l[sel]).sum())
    nm2 = []
    for e_i, k in enumerate(elig):
        xx = EST[name][k, GIDX[e_i]]
        dd = np.clip(np.searchsorted(qs[1:-1], np.abs(xx)), 0, 9)
        nm2.append(np.abs(amap[dd] * xx - LREF[e_i]).sum() / np.abs(LREF[e_i]).sum())
    res[name]["nmae_dec"] = np.array(nm2)

summary = {}
for name in EST:
    v = res[name]
    summary[name] = dict(pw_gate=float(v["pw_gate"].mean()),
                         pw_call=float(np.nanmean(v["pw_call"])),
                         up_call=float(np.nanmean(v["up_call"])),
                         nmae_own=float(v["nmae_own"].mean()),
                         nmae_e0mag=float(v["nmae_e0mag"].mean()),
                         nmae_global=float(v["nmae_global"].mean()),
                         nmae_dec=float(v["nmae_dec"].mean()),
                         nmae_pt=float(v["nmae_pt"].mean()),
                         nmae_ptwide=float(v["nmae_ptwide"].mean()),
                         a_global=v["a_global"])
    print(f"{name}: pw_gate={summary[name]['pw_gate']:.4f} "
          f"fid(up_call)={summary[name]['up_call']:.4f} nmae_own={summary[name]['nmae_own']:.4f} "
          f"nmae_e0mag={summary[name]['nmae_e0mag']:.4f} nmae_pt={summary[name]['nmae_pt']:.4f} "
          f"nmae_ptwide={summary[name]['nmae_ptwide']:.4f}")

# ---------- 3. reach proxy (walk_reach VERBATIM from reach348.py:125-145) ----------
def walk_reach(order_idx, match, adjud, n_conf):
    """order_idx: per-target ascending rank positions over pool genes."""
    n = match.shape[0]
    R = np.full(n, np.nan)
    Ks = np.zeros(n, dtype=int)
    for i in range(n):
        if n_conf[i] == 0:
            continue
        m = match[i][order_idx[i]]
        a = adjud[i][order_idx[i]]
        ma = m[a]
        denom = np.arange(1, ma.size + 1)
        if ma.size == 0:
            R[i] = 0.0
            continue
        pur = np.cumsum(ma) / denom
        ok = np.flatnonzero(pur >= 0.9)
        k = int(ok.max() + 1) if ok.size else 0
        Ks[i] = k
        R[i] = k / n_conf[i]
    return R, Ks

nT, nG = len(targets), len(genes)
qrm = np.asarray(r["q"], dtype=float)
lrm = np.asarray(r["lfc"], dtype=float)
ref_sig = qrm < 0.05
own = np.equal.outer(targets, genes)
pool = ref_sig & ~own
n_conf = pool.sum(1)
adjud = np.isfinite(lrm) & (lrm != 0)
feat = np.tile(genes, (nT, 1))
qm = np.where(np.isfinite(QHAT), QHAT, np.inf)
pm = np.where(np.isfinite(PHAT), PHAT, np.inf)
ala = np.abs(E0)
ala_m = np.where(np.isfinite(ala), -ala, np.inf)
sig = np.isfinite(QHAT) & (QHAT < 0.05)
cur = np.empty((nT, nG), dtype=int)
for i in range(nT):
    cur[i] = np.lexsort((feat[i], ala_m[i], pm[i], qm[i], ~sig[i]))
scored = n_conf > 0
print(f"reach: {nT} targets, {int(scored.sum())} scored")

for name, X in EST.items():
    commit = np.isfinite(X) & (X != 0)
    match = adjud & commit & (np.sign(X) == np.sign(lrm))
    Rcur, Kcur = walk_reach(cur, match & pool, adjud & pool, n_conf)
    ax = np.abs(X)
    ax_m = np.where(np.isfinite(ax), -ax, np.inf)
    o_abs = np.empty((nT, nG), dtype=int)
    for i in range(nT):
        o_abs[i] = np.lexsort((feat[i], ax_m[i]))
    Rabs, Kabs = walk_reach(o_abs, match & pool, adjud & pool, n_conf)
    res[name]["R_cur"] = Rcur
    res[name]["K_cur"] = Kcur
    res[name]["R_abs"] = Rabs
    res[name]["K_abs"] = Kabs
    summary[name]["reach_cur"] = float(np.nanmean(Rcur))
    summary[name]["k_cur"] = float(np.nanmean(Kcur[scored]))
    summary[name]["reach_abs"] = float(np.nanmean(Rabs))
    summary[name]["k_abs"] = float(np.nanmean(Kabs[scored]))
    print(f"{name}: reach_cur={summary[name]['reach_cur']:.5f} k*={summary[name]['k_cur']:.1f} "
          f"reach_abs={summary[name]['reach_abs']:.5f} k*={summary[name]['k_abs']:.1f}")

# ---------- 4. paired-target bootstrap vs E0 (seed 20260913, 342 convention) ----------
brng = np.random.default_rng(SEED_BOOT)
idx58 = np.arange(len(elig))
draws58 = brng.integers(0, len(elig), (DRAWS, len(elig)))
sc_idx = np.flatnonzero(scored)
draws_sc = brng.integers(0, len(sc_idx), (DRAWS, len(sc_idx)))
boot = {}
base_pw = res["E0"]["pw_gate"]
base_nm = res["E0"]["nmae_own"]
base_rc = res["E0"]["R_cur"][scored]
for name in EST:
    if name == "E0":
        continue
    d_pw = res[name]["pw_gate"] - base_pw
    b_pw = d_pw[draws58].mean(1)
    d_nm = base_nm - res[name]["nmae_own"]  # + means X better (lower nmae)
    b_nm = d_nm[draws58].mean(1)
    d_rc = res[name]["R_cur"][scored] - base_rc
    b_rc = d_rc[draws_sc].mean(1)
    boot[name] = dict(
        d_pw_mean=float(d_pw.mean()), d_pw_se=float(b_pw.std()),
        d_pw_p=float((b_pw > 0).mean()),
        d_nm_mean=float(d_nm.mean()), d_nm_se=float(b_nm.std()),
        d_nm_p=float((b_nm > 0).mean()),
        d_rc_mean=float(d_rc.mean()), d_rc_se=float(b_rc.std()),
        d_rc_p=float((b_rc > 0).mean()))
    print(f"{name} vs E0: d_pw={boot[name]['d_pw_mean']:+.4f} (se {boot[name]['d_pw_se']:.4f}, "
          f"P={boot[name]['d_pw_p']:.3f}) | d_nm={boot[name]['d_nm_mean']:+.4f} "
          f"(se {boot[name]['d_nm_se']:.4f}, P={boot[name]['d_nm_p']:.3f}) | "
          f"d_reach={boot[name]['d_rc_mean']:+.5f} (se {boot[name]['d_rc_se']:.5f}, "
          f"P={boot[name]['d_rc_p']:.3f})")

(OUT / "summary388.json").write_text(json.dumps(
    {"summary": summary, "bootstrap_vs_E0": boot,
     "notes": {"elig": len(elig), "scored_reach": int(scored.sum()),
               "E2_provenance": "h1test-feats.npz cons, SOURCES=(k562_gwps, orion_hek293t)",
               "E4": "N/A: donorstats.*.npz carry gene-level stats only (m/marg_abs/marg_sd/resp_frac/scale), no signed per-(target,gene) rows",
               "grids": "343 grid to 6x; wide grid to 24x diagnostic"}}, indent=1))
pl.DataFrame({"target": [targets[k] for k in elig],
              **{f"{n}_{c}": res[n][c].tolist()
                 for n in EST for c in ("pw_gate", "nmae_own", "nmae_e0mag", "nmae_pt")}}).write_csv(
    OUT / "per-target-388.csv")
print("wrote work/summary388.json work/per-target-388.csv work/est388.npz")

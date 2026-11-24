"""Realised reach: re-emit under each spread policy and score what the scorer would see.

Everything in experiment 095 so far is a COUNTERFACTUAL -- it reorders a fixed DE table by
fiat. The mechanisms we can actually ship set a per-gene `spread`, and the realised p-value
order is whatever the Wilcoxon then produces, which is only an approximation of the key.
Experiment 093's board result is what happens when you skip this step.

This drives the UNMODIFIED production `emit_counts` with each policy's spread vector, runs
the exact DE path, and reports realised reach. The trick that makes it cheap is exp093's:
in the quasi emitter the only thing carried from the profile into the block is the per-gene
panel total, and the emitter preserves that integer exactly, so the totals can be read off
the existing prediction and fed back.
"""
from __future__ import annotations
import os, sys
import numpy as np
import anndata as ad
from scipy.stats import false_discovery_control, norm

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from vcc2026.donor_transfer import emit_counts            # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z            # noqa: E402

# Harness and prediction are configurable so the same measurement can be repeated on a
# second, independent surface -- which is the only way to tell a real emission effect from
# one that is specific to a harness's geometry.
#
#   REAL=data/harness/h1-test.h5ad \
#   PRED=experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad \
#   COMBOS="production:0.25:0,production_b:0.25:0,phi4:0.25:4.0,phi4_b:0.25:4.0" \
#   uv run python realize.py
REAL = os.environ.get("REAL", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad")
PRED = os.environ.get(
    "PRED",
    f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
)
if not os.path.isabs(REAL): REAL = os.path.join(ROOT, REAL)
if not os.path.isabs(PRED): PRED = os.path.join(ROOT, PRED)
CTRL, COL, ALPHA, MIN_CPM = "non-targeting", "target_gene", 0.05, 5.0
N_EMIT, EMIT_L, SPREAD = 400, 100_000, 0.25
B_R, R_R = 0.079022, 0.964145
NP = int(os.environ.get("N_PERTS", "40"))

def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))

R = ad.read_h5ad(REAL)
rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
X = R.X; X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)
ci = np.flatnonzero(rlab == CTRL)
C = norm_log(X[ci]); cpm = np.expm1(C).mean(0); keep = cpm > MIN_CPM
genes_full = np.asarray([str(g) for g in R.var_names])

# --- gene-axis alignment: harness panel vs prediction panel --------------------------
# h1-test carries 18,080 genes against the production panel's 18,533, and everything
# downstream indexes on the HARNESS axis. Map the prediction onto it once, and drop from
# the gate any gene the prediction does not carry -- emitting zeros for those would feed
# a fake all-zero column into the rank test and score it as a real gene.
P = ad.read_h5ad(PRED)
plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
_pg = np.asarray([str(g) for g in P.var_names])
if np.array_equal(_pg, genes_full):
    PIDX = None
else:
    _pp = {g: i for i, g in enumerate(_pg)}
    PIDX = np.array([_pp.get(g, -1) for g in genes_full], dtype=np.int64)
    _miss = int((PIDX < 0).sum())
    keep &= PIDX >= 0
    print(f"gene-axis align: prediction {_pg.size} genes vs harness {genes_full.size}; "
          f"{_miss} harness genes absent from the prediction, dropped from the gate")
cs = np.ascontiguousarray(np.sort(C[:, keep].T, axis=1)); tie = ctrl_tiesums(cs)
ref_lin = np.expm1(C[:, keep]).mean(0); del C
G_full, G = genes_full.size, int(keep.sum())
gk = genes_full[keep]
fr = np.empty(G, np.int64); fr[np.argsort(gk)] = np.arange(G)
gpos = {g: i for i, g in enumerate(gk)}
print(f"harness {REAL.split('/')[-1]}: {ci.size} controls, {G}/{G_full} genes past the gate")

# Per-gene target dispersion, from the control cells' own RAW counts.
#
# CAREFUL WITH THE ALGEBRA -- the first version of this got it wrong and assigned phi=0 to
# 97.9% of genes. The emitted per-cell count is (even part) + NB(mu_s, phi) where
# `mu_s = spread * T/n` is the SCATTERED mean, but the count's own mean is the full `T/n`.
# Only the scattered part carries variance, so
#     Var = mu_s + phi * mu_s^2,      CV^2 = Var / (T/n)^2.
# Matching a gene's control CV therefore means
#     phi_g = (CV_ctrl^2 * (T/n)^2 - mu_s) / mu_s^2,
# NOT `CV^2 - 1/mu_s`, which is the formula for a gene whose whole total is scattered.
# Dividing by the wrong mean makes the target look already-exceeded and zeroes it out.
_ctrl_raw = X[ci]
_cm = _ctrl_raw.mean(0)
_cv2_full = np.zeros(G_full)
_nzc = _cm > 0
_cv2_full[_nzc] = (_ctrl_raw.std(0)[_nzc] / _cm[_nzc]) ** 2
print(f"control CV (gated genes): median {np.sqrt(_cv2_full[keep]).mean():.3f}")
del _ctrl_raw

# reference DE, once
targets = [t for t in sorted(set(rlab) - {CTRL}) if (plab == t).sum() > 0][:NP]
print(f"{len(targets)} perturbations\n")

# donor agreement key on the gated axis
DON = {}
for name in ("k562_gwps", "orion_hek293t"):
    root = f"{ROOT}/data/donors/{name}"
    dg = np.array(open(f"{root}/genes.txt").read().split())
    dp = np.array(open(f"{root}/perts.txt").read().split())
    gi = {g: i for i, g in enumerate(dg)}
    DON[name] = (np.load(f"{root}/delta.npy", mmap_mode="r"),
                 np.array([gi.get(g, -1) for g in gk]), {p: i for i, p in enumerate(dp)})
def agree(t):
    vs = []
    for n in DON:
        delta, cols, pi = DON[n]
        if t not in pi: return np.zeros(G)
        row = np.asarray(delta[pi[t]], float); v = np.zeros(G); ok = cols >= 0
        v[ok] = row[cols[ok]]; vs.append(v)
    a, b = vs
    return np.sqrt(np.abs(a)*np.abs(b)) * ((np.sign(a) == np.sign(b)) & (a != 0) & (b != 0))

def record(name, ps, rs, rl_, pl_, j, nc, gs, rlv, plv, idn):
    """fid / jac / nmae for one perturbation under one policy, scorer-exact."""
    ok = np.ones(rl_.size, bool)
    if j is not None: ok[j] = False
    ind = ok & (rl_ != 0) & ~np.isnan(rl_)
    scored = ps & ind
    match = scored & (pl_ != 0) & ~np.isnan(pl_) & (np.sign(pl_) == np.sign(rl_))
    n_real, n_pred_, k = int(rs.sum()), int(scored.sum()), int(match.sum())
    if max(n_pred_, n_real) > 0:
        FID[name].append(k / max(n_pred_, n_real))
    union = int((rs | (ps & ok)).sum())
    JAC[name].append(1.0 if union == 0 else int((rs & ps & ok).sum()) / union)
    gate = np.isfinite(rlv)
    if gate.sum() >= 10:
        den = np.abs(rlv[gate]).mean()
        if den > 0 and np.isfinite(den):
            NMAE[name].append(np.abs(np.where(np.isfinite(plv), plv, 0.0)[gate]
                                     - rlv[gate]).mean() / den)


def reach_of(order, idn, mt, nc):
    nd, nmt = np.cumsum(idn[order]), np.cumsum(mt[order])
    with np.errstate(invalid="ignore", divide="ignore"):
        pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
    g = np.flatnonzero((nd > 0) & (pur >= 0.9))
    return (nd[g[-1]] / nc) if g.size else 0.0

# (spread, phi). spread is the fraction of each gene's total that is scattered at all;
# phi is the gamma-Poisson dispersion of that scatter. Over-dispersion at full spread moves
# reach hugely but collapses n_pred, and the challenge's coverage term punishes that far
# harder than this harness can show. Lowering spread scatters less of the total, so the
# question is whether any combination buys reach while leaving n_pred near production.
# Policies as (name, spread, phi). Override without editing this file:
#
#   COMBOS="production:0.25:0,s0.25_phi3:0.25:3,s0.25_phi3_b:0.25:3" \
#   N_PERTS=60 uv run python realize.py
#
# ALWAYS include at least one duplicated policy. Over-dispersion carries an order of
# magnitude more emission noise than the Poisson emitter -- two phi=4 draws differ by 0.014
# raw reach against production's 0.0015 -- so a single-seed result is unreliable at the
# +-0.002 overall level. A name ending in `_b` is re-seeded, which is what makes it a
# control rather than a duplicate.
_DEFAULT = ("production:0.25:0,production_seed2:0.25:0,"
            "s0.25_phi3.0:0.25:3.0,s0.25_phi4.0:0.25:4.0,s0.25_phi4.0_b:0.25:4.0,"
            "s0.25_phi5.0:0.25:5.0")
COMBOS = [(n, float(sp), float(ph)) for n, sp, ph in
          (c.split(":") for c in os.environ.get("COMBOS", _DEFAULT).split(","))]
POLICIES = [c[0] for c in COMBOS]


def emit_nb(totals, phi, rng, n=N_EMIT):
    """Gamma-Poisson scatter with the SAME exact-total repair the quasi emitter uses.

    `spread` is capped at Poisson, which experiment 095 showed is far too weak to reorder
    p-values -- it moves them by a few orders of magnitude against a spread of hundreds.
    Over-dispersion is the one emission lever that can move them enough. Experiment 093 §7
    probed it on the depth-confounded harness (phi 0.15 gave the best fid it measured, phi
    0.40 sent reach to 0.101) and never validated it. This is that mechanism, on a repaired
    harness, through the same realisation path as every other policy here.

    The exact-total repair is what keeps `pds`, `mse` and `nmae` untouched, so it is
    preserved verbatim: whatever the scatter, each gene's panel total is restored to the
    integer the profile implies.
    """
    G_ = totals.size
    out = np.zeros((n, G_), dtype=np.int64)
    nz = np.flatnonzero(totals > 0)
    mu = totals[nz] / n
    if phi <= 0:
        draw = rng.poisson(np.broadcast_to(mu, (n, nz.size)))
    else:
        # NB as a gamma mixture of Poissons: Var = mu + phi*mu^2
        lam = rng.gamma(shape=1.0 / phi, scale=phi * np.broadcast_to(mu, (n, nz.size)))
        draw = rng.poisson(lam)
    # exact-total repair, same policy as emit_counts: add at random cells, remove only
    # from cells that still hold a count so nothing goes negative
    for c, g in enumerate(nz):
        col = draw[:, c]
        resid = int(totals[g]) - int(col.sum())
        if resid > 0:
            idx = rng.integers(0, n, resid)
            np.add.at(col, idx, 1)
        while resid < 0:
            hot = np.flatnonzero(col > 0)
            take = min(hot.size, -resid)
            pick = rng.choice(hot, take, replace=False)
            col[pick] -= 1
            resid += take
        draw[:, c] = col
    out[:, nz] = draw
    return out.astype(np.float64)
acc = {k: [] for k in POLICIES}; npred = {k: [] for k in POLICIES}
FID = {k: [] for k in POLICIES}; JAC = {k: [] for k in POLICIES}
NMAE = {k: [] for k in POLICIES}
for i, t in enumerate(targets):
    rrows = np.flatnonzero(rlab == t)
    Br = norm_log(X[rrows])
    z = wilcoxon_z(cs, tie, np.ascontiguousarray(Br[:, keep].T))
    rp = 2*norm.sf(np.abs(z)); rq = false_discovery_control(rp, method="bh")
    rl = np.log2((np.expm1(Br[:, keep]).mean(0)+1e-9)/(ref_lin+1e-9)); del Br
    rs = rq < ALPHA
    j = gpos.get(t)
    if j is not None: rs[j] = False
    gs = np.flatnonzero(rs); nc = gs.size
    if nc == 0: continue

    prows = np.flatnonzero(plab == t)
    blk = P[prows].to_memory().X
    blk = np.asarray(blk.todense() if hasattr(blk, "todense") else blk, dtype=np.float64)
    totals = blk.sum(0).astype(np.int64); del blk
    if PIDX is not None:                      # prediction axis -> harness axis
        _t = np.zeros(G_full, np.int64); _m = PIDX >= 0
        _t[_m] = totals[PIDX[_m]]; totals = _t
    profile = (totals + 0.5) / (N_EMIT * EMIT_L)
    eff = np.zeros(G_full); eff[keep] = np.abs(rl)          # |lfc| stands in for |effect|
    ag = np.zeros(G_full); ag[keep] = agree(t)
    gated = np.zeros(G_full, bool); gated[keep] = True

    def graded(score, lo, hi):
        v = np.full(G_full, SPREAD); gi_ = np.flatnonzero(gated)
        o = np.lexsort((gi_, -score[gi_])); rk = np.empty(gi_.size); rk[o] = np.arange(gi_.size)
        v[gi_] = lo + (hi-lo) * (rk / max(gi_.size-1, 1)); return v
    # Buildable |predicted effect| on the full axis: our own emitted profile against the
    # control's expression. NOTHING from the reference enters this -- an earlier version
    # used the reference's own lfc here, which was an oracle leak (see RESULTS section 8a).
    _pc = np.zeros(G_full); _tot_sum = max(totals.sum(), 1)
    _pc[:] = 1e6 * totals / _tot_sum
    _cc = np.zeros(G_full); _cc[keep] = cpm[keep]
    eff_pred = np.abs(np.log2((_pc + 1e-9) / (_cc + 1e-9)))
    gated_f = np.zeros(G_full, bool); gated_f[keep] = True
    eff_pred = np.where(gated_f, eff_pred, 0.0)

    order_eff = np.lexsort((np.arange(G_full), -np.where(gated_f, eff_pred, -np.inf)))
    def block(k, tail=SPREAD, kept=0.0):
        v = np.full(G_full, tail); v[order_eff[:k]] = kept; return v
    vecs = {}
    for name, sp, phi in COMBOS:
        if name.startswith("keeptight"):
            # Over-disperse only the genes we do NOT trust: phi=0 on the top-K by
            # |predicted effect| so they stay tight and significant, phi=PHI on the rest so
            # they slide down the p-order and out of the call set. The intent is a8's reach
            # gain WITHOUT the call-set collapse that cost it `jac`.
            K = int(name.split("K")[1].split("_")[0])
            order = np.lexsort((np.arange(G_full), -np.where(gated_f, eff_pred, -np.inf)))
            per = np.full(G_full, float(phi)); per[order[:K]] = 0.0
            vecs[name] = (np.full(G_full, sp), per)
        elif name.startswith("pergene"):
            mean_cell = np.maximum(totals / N_EMIT, 1e-12)
            mu_s = np.maximum(mean_cell * sp, 1e-12)
            target_var = _cv2_full * mean_cell ** 2
            # cap so a handful of very tight genes cannot demand an absurd phi
            per = np.clip((target_var - mu_s) / mu_s ** 2, 0.0,
                          float(phi) if phi > 0 else 50.0)
            vecs[name] = (np.full(G_full, sp), per)
        else:
            vecs[name] = (np.full(G_full, sp), phi)

    for name, (sv, phi_) in vecs.items():
        if False:
            phi = float(name.split("phi")[1])
            rng = np.random.default_rng(90_000 + i)
            blk2 = emit_nb(totals, phi, rng)
            assert np.array_equal(blk2.sum(0).astype(np.int64), totals), name
            nb = norm_log(blk2); del blk2
            z2 = wilcoxon_z(cs, tie, np.ascontiguousarray(nb[:, keep].T))
            pp2 = 2*norm.sf(np.abs(z2)); pq2 = false_discovery_control(pp2, method="bh")
            pl2 = np.log2((np.expm1(nb[:, keep]).mean(0)+1e-9)/(ref_lin+1e-9)); del nb
            ps = pq2 < ALPHA
            if j is not None: ps[j] = False
            npred[name].append(int(ps.sum()))
            rlv, plv = rl[gs], pl2[gs]
            idn = (~np.isnan(rlv)) & (rlv != 0)
            mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
            absl = np.abs(np.nan_to_num(plv, nan=-np.inf))
            acc[name].append(reach_of(
                np.lexsort((fr[gs], -absl, pp2[gs], pq2[gs], ~ps[gs])), idn, mt, nc))
            record(name, ps, rs, rl, pl2, j, nc, gs, rlv, plv, idn)
            continue
        # production_seed2 is an IDENTICAL policy on a different emission seed: it measures
        # how much of any delta below is just emission noise. Without it this table cannot
        # distinguish a real effect from a re-roll.
        _seed = 90_000
        if name == "production_seed2": _seed = 5_000_000
        elif name.endswith("_b"): _seed = 7_000_000   # second seed on the leading arm
        rng = np.random.default_rng(_seed + i)
        _od = phi_
        vals, cols_, indptr = emit_counts(profile, np.full(N_EMIT, EMIT_L), 0.1, rng,
                                          mode="quasi", spread=sv, overdispersion=_od)
        blk2 = np.zeros((N_EMIT, G_full))
        for r in range(N_EMIT):
            sl = slice(indptr[r], indptr[r+1]); blk2[r, cols_[sl]] = vals[sl]
        assert np.array_equal(blk2.sum(0).astype(np.int64), totals), name
        nb = norm_log(blk2); del blk2
        z2 = wilcoxon_z(cs, tie, np.ascontiguousarray(nb[:, keep].T))
        pp2 = 2*norm.sf(np.abs(z2)); pq2 = false_discovery_control(pp2, method="bh")
        pl2 = np.log2((np.expm1(nb[:, keep]).mean(0)+1e-9)/(ref_lin+1e-9)); del nb
        ps = pq2 < ALPHA
        if j is not None: ps[j] = False
        npred[name].append(int(ps.sum()))
        rlv, plv = rl[gs], pl2[gs]
        idn = (~np.isnan(rlv)) & (rlv != 0)
        mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
        absl = np.abs(np.nan_to_num(plv, nan=-np.inf))
        acc[name].append(reach_of(
            np.lexsort((fr[gs], -absl, pp2[gs], pq2[gs], ~ps[gs])), idn, mt, nc))
        record(name, ps, rs, rl, pl2, j, nc, gs, rlv, plv, idn)
    if (i+1) % 10 == 0: print(f"  {i+1}/{len(targets)}", flush=True)

B_F, R_F = 0.51239, 0.80964
B_J, R_J = 0.030376, 0.397270
b_r, b_f, b_j, b_n = (np.mean(acc["production"]), np.nanmean(FID["production"]),
                      np.mean(JAC["production"]), np.mean(NMAE["production"]))
print(f"\n{'policy':<18}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}"
      f"{'Δreach':>10}{'Δfid':>10}{'Δjac':>10}{'ΔTOTAL':>10}")
for k in POLICIES:
    r, f, jj, nm = (np.mean(acc[k]), np.nanmean(FID[k]), np.mean(JAC[k]), np.mean(NMAE[k]))
    dr = (r-b_r)/(R_R-B_R)/6; df = (f-b_f)/(R_F-B_F)/6; dj = (jj-b_j)/(R_J-B_J)/6
    print(f"{k:<18}{np.median(npred[k]):>8.0f}{r:>9.5f}{f:>9.4f}{jj:>9.4f}{nm:>9.4f}"
          f"{dr:>+10.5f}{df:>+10.5f}{dj:>+10.5f}{dr+df+dj:>+10.5f}")
print("\n(nmae is printed as a guard: the exact-total repair should hold it constant.)")

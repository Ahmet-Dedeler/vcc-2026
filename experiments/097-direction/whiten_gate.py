"""GATE 3: price the WHITENER across all six members, on top of a13, on both harnesses.

097 measured that dropping the whitening operator (`gamma` 0.65 -> 0) is worth 1.495 / 2.136 on
directional accuracy -- the largest single-lever number in the experiment, monotone on both
surfaces -- and declined to recommend it because `WHITEN_GAMMA` was fitted for `pds_cosine` in
experiment 017 and `pds` is the member we lead. This prices that trade.

Reproduces `experiments/031-somatic/run.py`'s build path exactly:

    E   = MultiDonorModel(donors, gene_norm=False, eb=False).effect(p, alpha=1, beta=0, uniform)
    E[:, wcols] = W.apply(E[:, wcols], gamma, 0.0, renorm=False), then median-renormalised
    E  *= alpha * effect_scale / median||E||
    profile = model.profile_from_effect(E[i], b_ctrl, min_cpm=5, gate_unexpressed=True)

with the CACHED harness whitener (never refitted), and `wcols` rebuilt and sha-checked against
the cached meta. `gamma = 0.65` therefore reproduces the shipped prediction and is the control.

On top of that, a13's shrinkage exactly as the builder applies it:
    alpha_g = exp(-gamma_gene * descending_midrank(donor_agreement)_g / 2000), target gene protected
    z_p     = max(rank01(mean pairwise donor cosine), 0.3) ** gamma_pert, missing -> damp
and each arm is then re-optimised onto its own rho cliff before emission.

    uv run python whiten_gate.py <h1|orion>
"""
from __future__ import annotations
import hashlib, json, os, sys, time
import numpy as np
import pandas as pd
from scipy.stats import rankdata, false_discovery_control, norm
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
import anndata as ad
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget
from cell_eval2.metrics.discrimination import discrimination_score
from decompose import terms
from msekit import CTRL, TS, group_stats
from vcc2026.donor_transfer import Donor, emit_counts
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.whiten import Whitener
from de_kernel import ctrl_tiesums, wilcoxon_z

WHICH = sys.argv[1]
REAL, PRED, OUTDIR = {
 "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test"),
 "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out")}[WHICH]
KEYDONS = {"h1": ["k562_gwps", "orion_hek293t", "orion_hct116"],
           "orion": ["k562_gwps", "orion_hek293t"]}[WHICH]     # leak rule: no self-donor
MODELDONS = ["k562_gwps", "orion_hek293t"]                     # 031-somatic's model, both harnesses
GAMMAS = [float(x) for x in os.environ.get("GAMMAS", "0.0,0.15,0.30,0.45,0.65,0.85").split(",")]
NP_ = int(os.environ.get("N_PERTS", "40"))
ALPHA_DE, MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL_SPREAD, PHI = 0.05, 5.0, 8000, 0.25, 1.0, 4.0
G_GENE, G_PERT, FLOOR, ALPHA_V, EFFSCALE = 0.25, 1.0, 0.3, 1.08, 6.817

rp, rs, rc, rj, rg = group_stats(REAL); pp, ps, pc, pj, pg = group_stats(PRED)
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
for k, c in t["cols"].items():
    t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
w, cols, sD, claim = t["w"], t["cols"], t["sumD"], t["claim"]
keep_i = t["keep"]; tg = np.array([str(x) for x in pp[keep_i]])
P, G = t["dh"].shape
ci_r = int(np.flatnonzero(rp == CTRL)[0])
bctrl = bulk_lognorm_means(rs[ci_r][None, :].astype(np.float64), TS)[0]
rmeans = bulk_lognorm_means(rs.astype(np.float64), TS)
S = ps.sum(1).astype(np.float64); NC = int(round(float(np.median(pc))))
print(f"# {WHICH}: P={P} scored rows, G={G}, NCELL={NC}", flush=True)

# ---- rebuild 031-somatic's model and wcols, sha-checked -------------------------------------
donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in MODELDONS}
model = MultiDonorModel(donors, rg, gene_norm=False, eb=False)
shared = set()
for n in model.names:
    shared |= set(model.target_genes[model.tcols[n]])
chal = set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str))
real_t = set(str(x) for x in rp if str(x) != CTRL)
shared -= (chal if real_t <= chal else real_t)
ti = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
meta = json.load(open(f"{OUTDIR}/whitener_union.json"))
sha = hashlib.sha256(wcols.tobytes()).hexdigest()
print(f"# wcols {wcols.size} vs cached {meta['n_cols']}; sha "
      f"{'MATCH' if sha == meta['wcols_sha256'] else 'MISMATCH -- STOP'}", flush=True)
assert sha == meta["wcols_sha256"], "wcols does not reproduce the cached operator"
E0 = np.vstack([model.effect(str(x), 1.0, 0.0, weights={n: 1.0 for n in model.names})
                for x in tg]).astype(np.float64)
W = Whitener.load(f"{OUTDIR}/whitener_union.npz")
print(f"# whitener loaded, {W.evals.size} evals, condition {W.evals.max()/W.evals.min():.3e}",
      flush=True)


def effect_at(gam):
    E = E0.copy()
    if gam > 0:
        before = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
        E[:, wcols] = W.apply(E[:, wcols], gam, 0.0, renorm=False)
        E[:, wcols] *= before / (float(np.median(np.linalg.norm(E[:, wcols], axis=1))) or 1.0)
    E *= ALPHA_V * EFFSCALE / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
    return E


def prof_of(eff_row):
    return model.profile_from_effect(eff_row, bctrl, min_cpm=MIN_CPM, gate_unexpressed=True)


# ---- MECHANISM CHECK: gamma=0.65 must reproduce the shipped prediction ----------------------
E65 = effect_at(0.65)
dh_rec = np.vstack([np.log1p(TS * prof_of(E65[i])) - bctrl for i in range(P)])
for k, c in cols.items():
    dh_rec[k, c] = 0.0
dh_ship = t["dh"]
csim = np.einsum("ij,ij->i", dh_rec, dh_ship) / np.maximum(
    np.linalg.norm(dh_rec, axis=1) * np.linalg.norm(dh_ship, axis=1), 1e-30)
kap_rec = float((np.einsum("ij,ij->i", dh_rec, t["d"]) * w).sum()) / np.sqrt(
    float(_across_pert_budget(dh_rec, cols, w)) * sD)
kap_ship = t["X"] / np.sqrt(t["budget"] * sD)
print(f"# RECONSTRUCTION CHECK at gamma=0.65: per-row cosine vs the shipped dh "
      f"median {np.median(csim):.5f} min {csim.min():.5f};  kappa {kap_rec:.5f} vs shipped "
      f"{kap_ship:.5f} (ratio {kap_rec/kap_ship:.4f})", flush=True)

# ---- keys (unchanged by gamma: they come from the donors, not from our effect) --------------
M_, H_ = {}, {}
for nd in KEYDONS:
    dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m_, s_ = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cc = np.array([gi.get(str(g), -1) for g in rg]); ok = cc >= 0
    pi = {p_: i for i, p_ in enumerate(dn.perts)}
    Mx = np.zeros((P, G)); h = np.zeros(P, bool)
    for r, x in enumerate(tg):
        if x in pi: Mx[r, ok] = s_[pi[x]][cc[ok]]; h[r] = True
    M_[nd] = Mx / (float(np.median(np.linalg.norm(s_[:, cc[ok]], axis=1))) or 1.0)
    H_[nd] = h; del dn, s_
ag = np.zeros((P, G)); dc = np.zeros(P); ndon = np.zeros(P, int)
for i in range(P):
    have = [n for n in KEYDONS if H_[n][i]]; ndon[i] = len(have)
    if len(have) < 2: continue
    st = np.vstack([M_[n][i] for n in have]); sg = np.sign(st)
    un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
    ag[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
    tot, cnt = 0.0, 0
    for a in range(len(have)):
        for b in range(a + 1, len(have)):
            x, y = M_[have[a]][i], M_[have[b]][i]
            tot += float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12); cnt += 1
    dc[i] = tot / cnt
zp = np.maximum(np.where(ndon >= 2, (rankdata(dc, method="average") - 0.5) / P, 0.0), FLOOR) ** G_PERT


def desc_midrank(x):
    n = x.size; o = np.argsort(-x, kind="stable"); od = x[o]
    e = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = e[:-1], e[1:]
    out = np.empty(n); out[o] = np.repeat((st + en - 1) / 2.0, en - st)
    return out


AL = np.ones((P, G))
for i in range(P):
    if np.unique(ag[i]).size >= 2:
        AL[i] = np.exp(-G_GENE * desc_midrank(ag[i]) / 2000.0)
    if i in cols: AL[i, cols[i]] = 1.0          # the builder protects the target's own gene
print(f"# a13 keys: mean per-gene alpha {AL.mean():.4f}, mean z {zp.mean():.4f}, "
      f"rows with <2 donors {int((ndon<2).sum())}", flush=True)

# ---- arms ----------------------------------------------------------------------------------
ARMS = [("a12_g0.65", 0.65, False)] + [(f"a13_g{g:g}", g, True) for g in GAMMAS]
PRE = {}
for nm, gam, shrink in ARMS:
    E = effect_at(gam)
    dh = np.vstack([np.log1p(TS * prof_of(E[i])) - bctrl for i in range(P)])
    for k, c in cols.items():
        dh[k, c] = 0.0
    if shrink:
        dh = zp[:, None] * AL * dh
    B = float(_across_pert_budget(dh, cols, w))
    s_c = float(np.sqrt(claim / B))
    X = float((np.einsum("ij,ij->i", dh, t["d"]) * w).sum())
    kap = X / np.sqrt(B * sD)
    pm = np.tile(bctrl, (P + 1, 1)); pm[:P] += s_c * dh
    pv = discrimination_score(pred_bulk=(np.append(tg, CTRL), pm), real_bulk=(rp, rmeans),
                              pert_col="target_gene", control=CTRL, distance="cosine",
                              exclude_target_gene=True, exclusion_scope="panel",
                              control_source="pred", genes=rg)
    pds = float(np.mean(list(pv.values()))) if isinstance(pv, dict) else float(pv)
    PRE[nm] = dict(E=E, shrink=shrink, s=s_c, kap=kap, pds=pds)
    print(f"#   {nm:<12} gamma={gam:<5} shrink={int(shrink)}  cliff s={s_c:.4f}  "
          f"kappa={kap:.5f}  pds={pds:.4f}", flush=True)

# ---- DE path -------------------------------------------------------------------------------
R = ad.read_h5ad(REAL)
rlab = np.asarray(R.obs["target_gene"].astype(str).values)
Xr = R.X; Xr = np.asarray(Xr.todense() if hasattr(Xr, "todense") else Xr, dtype=np.float64)
def norm_log(x):
    sm = x.sum(1); sm[sm == 0] = 1.0
    return np.log1p(x * (1e6 / sm[:, None]))
cidx = np.flatnonzero(rlab == CTRL)
Cn = norm_log(Xr[cidx]); cpm = np.expm1(Cn).mean(0); gate = cpm > MIN_CPM
cs = np.ascontiguousarray(np.sort(Cn[:, gate].T, axis=1)); tie = ctrl_tiesums(cs)
ref_lin = np.expm1(Cn[:, gate]).mean(0); del Cn
Gk = int(gate.sum()); gk = rg[gate]
fr = np.empty(Gk, np.int64); fr[np.argsort(gk)] = np.arange(Gk)
gpos = {g: i for i, g in enumerate(gk)}
expr = np.expm1(bctrl); cpm_b = 1e6 * expr / max(expr.sum(), 1e-12)
gated_b = (cpm_b > MIN_CPM) & (bctrl > 0)
POL = [(nm, sd) for nm, _, _ in ARMS for sd in (0, 1)]
RES = {f"{nm}{'_b' if sd else ''}": {k: [] for k in ("reach", "fid", "jac", "nmae", "npred")}
       for nm, sd in POL}
t0 = time.time()
for i_t, tname in enumerate(tg[:NP_]):
    ridx = np.flatnonzero(rlab == tname)
    if ridx.size == 0: continue
    Br = norm_log(Xr[ridx])
    z_ = wilcoxon_z(cs, tie, np.ascontiguousarray(Br[:, gate].T))
    rq = false_discovery_control(2 * norm.sf(np.abs(z_)), method="bh")
    rl = np.log2((np.expm1(Br[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del Br
    rsig = rq < ALPHA_DE
    j = gpos.get(tname)
    if j is not None: rsig[j] = False
    gs = np.flatnonzero(rsig); nc = gs.size
    if nc == 0: continue
    row = int(np.flatnonzero(tg == tname)[0])
    for nm, sd in POL:
        pr = PRE[nm]
        eff = pr["E"][row].copy()
        if pr["shrink"]:
            # shrinkage acts in log-effect space, exactly as the builder applies it
            eff = zp[row] * AL[row] * eff
        prof = prof_of(pr["s"] / (ALPHA_V * 1.0) * eff * (ALPHA_V * 1.0))
        prof = prof_of(pr["s"] * eff)
        sc_ = np.where(gated_b, np.abs(pr["s"] * eff), -np.inf)
        order = np.lexsort((np.arange(G), -sc_))
        sv = np.full(G, TAIL_SPREAD); sv[order[:min(KEEP_TOP, int(gated_b.sum()))]] = KEEP_SPREAD
        v, col, ip = emit_counts(prof, np.full(NC, S[row] / NC), 0.0,
                                 np.random.default_rng(7919 + i_t + 10_000 * sd),
                                 mode="quasi", spread=sv, overdispersion=PHI)
        blk = np.zeros((NC, G))
        for r_ in range(NC):
            blk[r_, col[ip[r_]:ip[r_+1]]] = v[ip[r_]:ip[r_+1]]
        nb = norm_log(blk); del blk
        z2 = wilcoxon_z(cs, tie, np.ascontiguousarray(nb[:, gate].T))
        p2 = 2 * norm.sf(np.abs(z2)); q2 = false_discovery_control(p2, method="bh")
        pl = np.log2((np.expm1(nb[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del nb
        psig = q2 < ALPHA_DE
        if j is not None: psig[j] = False
        key = f"{nm}{'_b' if sd else ''}"
        RES[key]["npred"].append(int(psig.sum()))
        ok = np.ones(Gk, bool)
        if j is not None: ok[j] = False
        ind = ok & (rl != 0) & ~np.isnan(rl)
        scored = psig & ind
        match = scored & (pl != 0) & ~np.isnan(pl) & (np.sign(pl) == np.sign(rl))
        nr_, np_, kk = int(rsig.sum()), int(scored.sum()), int(match.sum())
        if max(np_, nr_) > 0: RES[key]["fid"].append(kk / max(np_, nr_))
        union = int((rsig | (psig & ok)).sum())
        RES[key]["jac"].append(1.0 if union == 0 else int((rsig & psig & ok).sum()) / union)
        rlv, plv = rl[gs], pl[gs]
        idn = (~np.isnan(rlv)) & (rlv != 0)
        mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
        den = np.abs(rlv[np.isfinite(rlv)]).mean()
        if den > 0: RES[key]["nmae"].append(np.abs(np.nan_to_num(plv) - rlv).mean() / den)
        o = np.lexsort((fr[gs], -np.abs(np.nan_to_num(plv, nan=-np.inf)), p2[gs], q2[gs], ~psig[gs]))
        nd_, nmt = np.cumsum(idn[o]), np.cumsum(mt[o])
        with np.errstate(invalid="ignore", divide="ignore"):
            pur = np.where(nd_ > 0, nmt / np.maximum(nd_, 1), np.nan)
        gg_ = np.flatnonzero((nd_ > 0) & (pur >= 0.9))
        RES[key]["reach"].append((nd_[gg_[-1]] / nc) if gg_.size else 0.0)
    if (i_t + 1) % 10 == 0:
        print(f"#   {i_t+1} perts [{time.time()-t0:.0f}s]", flush=True)

print(f"\n# {WHICH} WHITENER GATE  ({len(RES['a12_g0.65']['reach'])} perturbations)")
print(f"  {'arm':<14}{'kappa':>9}{'pds':>9}{'n_pred':>8}{'reach':>8}{'reach_b':>9}"
      f"{'fid':>9}{'jac':>8}{'nmae':>9}")
for nm, _, _ in ARMS:
    a, b = RES[nm], RES[nm + "_b"]
    if not a["reach"]: continue
    f = lambda d, k: float(np.mean(d[k]))
    print(f"  {nm:<14}{PRE[nm]['kap']:>9.5f}{PRE[nm]['pds']:>9.4f}"
          f"{(f(a,'npred')+f(b,'npred'))/2:>8.0f}{f(a,'reach'):>8.4f}{f(b,'reach'):>9.4f}"
          f"{(f(a,'fid')+f(b,'fid'))/2:>9.4f}{(f(a,'jac')+f(b,'jac'))/2:>8.4f}"
          f"{(f(a,'nmae')+f(b,'nmae'))/2:>9.4f}")

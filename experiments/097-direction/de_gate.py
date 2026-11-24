"""GATE 1: price the four DE members for a13's exact configuration, on BOTH harnesses.

Everything 097 and 096 measured on these arms is mse-only. This runs the real DE path --
production `emit_counts` with a13's tiered spread vector, the scanpy-identical rank-sum against
the harness's own control cells, BH at 0.05 -- and reports reach / fid / jac / nmae / n_pred.

STANDING WARNINGS, per RESEARCH-METHOD rule 1, reproduced in the output:
  * `n_real` here is 85/132-panel-scale, not the challenge's ~2,000, so the LEVELS of `fid`,
    `jac` and `reach` do not transfer. Only `fid`'s PRECISION leg k/n_pred and the
    ordering-driven part of `reach` do. `jac` does not factorise at all -- read its SIGN only.
  * `n_pred` is reported because `fid`'s coverage leg min(1, n_pred/n_real) is the term that
    has cost this project three board predictions (v14, v15, a8).
  * Every arm carries a re-seeded twin (`_b`) as the noise floor. Over-dispersion carries an
    order of magnitude more emission noise than the Poisson emitter, so a single-seed
    difference smaller than the twin gap is not a result.

    ARMS="a12,a13,a13_gp2,rec" N_PERTS=40 uv run python de_gate.py <h1|orion>
"""
from __future__ import annotations
import os, sys, time
import numpy as np
from scipy.stats import rankdata, false_discovery_control, norm
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
import anndata as ad
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget
from decompose import terms
from msekit import CTRL, TS, group_stats
from vcc2026.donor_transfer import emit_counts, Donor
from de_kernel import ctrl_tiesums, wilcoxon_z

WHICH = sys.argv[1]
REAL, PRED, DONS, THIRD = {
 "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], "orion_hct116"),
 "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], None)}[WHICH]
NP_ = int(os.environ.get("N_PERTS", "40"))
ALPHA_DE, MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL_SPREAD, PHI = 0.05, 5.0, 8000, 0.25, 1.0, 4.0

rp, rs, rc, rj, rg = group_stats(REAL); pp, ps, pc, pj, pg = group_stats(PRED)
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
for k, c in t["cols"].items():
    t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
dh0, w, cols, claim = t["dh"], t["w"], t["cols"], t["claim"]
keep_i = t["keep"]; tg = np.array([str(x) for x in pp[keep_i]])
P, G = dh0.shape
ci_r = int(np.flatnonzero(rp == CTRL)[0])
bctrl = bulk_lognorm_means(rs[ci_r][None, :].astype(np.float64), TS)[0]
S = ps.sum(1).astype(np.float64); NC = int(round(float(np.median(pc))))

# ---- keys (production's all-donor-unanimity agreement + mean pairwise donor cosine) --------
dons = DONS + ([THIRD] if THIRD else [])
M, H = {}, {}
for nd in dons:
    dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s_ = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cc = np.array([gi.get(str(g), -1) for g in rg]); ok = cc >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    Mx = np.zeros((P, G)); h = np.zeros(P, bool)
    for r, x in enumerate(tg):
        if x in pi: Mx[r, ok] = s_[pi[x]][cc[ok]]; h[r] = True
    M[nd] = Mx / (float(np.median(np.linalg.norm(s_[:, cc[ok]], axis=1))) or 1.0)
    H[nd] = h; del dn, s_
ag = np.zeros((P, G)); dc = np.zeros(P); ndon = np.zeros(P, int)
for i in range(P):
    have = [n for n in dons if H[n][i]]; ndon[i] = len(have)
    if len(have) < 2: continue
    st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
    un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
    ag[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
    tot, cnt = 0.0, 0
    for a in range(len(have)):
        for b in range(a + 1, len(have)):
            x, y = M[have[a]][i], M[have[b]][i]
            tot += float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12); cnt += 1
    dc[i] = tot / cnt
al = np.empty((P, G))
for i in range(P):
    al[i] = rankdata(ag[i], method="average") / G
    if np.unique(ag[i]).size < 2: al[i] = 1.0
zr = (rankdata(dc, method="average") - 0.5) / P; zr[ndon < 2] = 0.0
zr = np.maximum(zr, 0.3)

def _desc_midrank(x):
    n = x.size
    o = np.argsort(-x, kind="stable"); od = x[o]
    e = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = e[:-1], e[1:]
    out = np.empty(n); out[o] = np.repeat((st + en - 1) / 2.0, en - st)
    return out


def alpha_of(form, g):
    """097's power form vs the BUILDER's exponential form. They are different families -- see
    scale_check.py. `pow` is (ascending_midrank/G)**g; `exp` is
    exp(-g * descending_midrank / 2000), which is `agreement_shrink_vector` verbatim including
    its <2-distinct-values guard."""
    A = np.empty((P, G))
    for i in range(P):
        if np.unique(ag[i]).size < 2:
            A[i] = 1.0
        elif form == "pow":
            A[i] = (rankdata(ag[i], method="average") / G) ** g
        else:
            A[i] = np.exp(-g * _desc_midrank(ag[i]) / 2000.0)
    return A


# (form, gamma_gene, gamma_pert)
ARMS = {"a12":       ("pow", 0.0,  0.0),    # no shrinkage at all: a12's base
        "a13":       ("pow", 0.25, 1.0),    # what 097 GATED -- 097's own parameterisation
        "a13_gp2":   ("pow", 0.25, 2.0),
        "rec":       ("pow", 1.0,  2.0),
        "bldr_g1":   ("exp", 1.0,  1.0),    # the BUILDER's default = a13 AS BUILT
        "bldr_g0.25":("exp", 0.25, 1.0),
        "bldr_g0.05":("exp", 0.05, 1.0)}    # ~ 097's gg=0.25 by kappa
SEL = os.environ.get("ARMS", "a12,a13,a13_gp2,rec").split(",")
POL = []
for nm in SEL:
    form, gg, gp = ARMS[nm]
    W = (zr ** gp)[:, None] * alpha_of(form, gg)
    B = float(_across_pert_budget(W * dh0, cols, w))
    POL.append((nm, W, float(np.sqrt(claim / B)), 0))
    POL.append((nm + "_b", W, float(np.sqrt(claim / B)), 1))   # re-seeded twin
print(f"# {WHICH}: P={P} scored rows, NCELL={NC}, {len(dons)}-donor keys, "
      f"rows with <2 donors {int((ndon<2).sum())}", flush=True)
for nm, W, s_c, sd in POL:
    if sd == 0:
        print(f"#   {nm:<12} form={ARMS[nm][0]} gamma=(gene {ARMS[nm][1]}, pert {ARMS[nm][2]})"
              f"  cliff s={s_c:.4f}  mean weight {W.mean():.4f}", flush=True)

# ---- harness DE setup ----------------------------------------------------------------------
R = ad.read_h5ad(REAL)
rlab = np.asarray(R.obs["target_gene"].astype(str).values)
Xr = R.X; Xr = np.asarray(Xr.todense() if hasattr(Xr, "todense") else Xr, dtype=np.float64)
def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))
cidx = np.flatnonzero(rlab == CTRL)
Cn = norm_log(Xr[cidx]); cpm = np.expm1(Cn).mean(0); gate = cpm > MIN_CPM
cs = np.ascontiguousarray(np.sort(Cn[:, gate].T, axis=1)); tie = ctrl_tiesums(cs)
ref_lin = np.expm1(Cn[:, gate]).mean(0); del Cn
Gk = int(gate.sum()); gk = rg[gate]
fr = np.empty(Gk, np.int64); fr[np.argsort(gk)] = np.arange(Gk)
gpos = {g: i for i, g in enumerate(gk)}
print(f"# {cidx.size} controls, {Gk}/{G} genes past the {MIN_CPM} CPM gate", flush=True)

expression = np.expm1(bctrl); cpm_b = 1e6 * expression / max(expression.sum(), 1e-12)
gated_b = (cpm_b > MIN_CPM) & (bctrl > 0)
RES = {nm: {k: [] for k in ("reach", "fid", "jac", "nmae", "npred")} for nm, _, _, _ in POL}
tgt = [x for x in tg[:NP_]]
t0 = time.time()
for i_t, tname in enumerate(tgt):
    ridx = np.flatnonzero(rlab == tname)
    if ridx.size == 0: continue
    Br = norm_log(Xr[ridx])
    z = wilcoxon_z(cs, tie, np.ascontiguousarray(Br[:, gate].T))
    rq = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh")
    rl = np.log2((np.expm1(Br[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del Br
    rsig = rq < ALPHA_DE
    j = gpos.get(tname)
    if j is not None: rsig[j] = False
    gs = np.flatnonzero(rsig); nc = gs.size
    if nc == 0: continue
    row = int(np.flatnonzero(tg == tname)[0])
    for nm, W, s_c, sd in POL:
        eff = s_c * W[row] * dh0[row]
        prof = np.expm1(bctrl + eff); np.clip(prof, 0.0, None, out=prof)
        tot = prof.sum()
        if tot <= 0: continue
        prof = prof / tot
        # a13's tiered spread: top-8000 gated genes by |effect| keep 0.25, the tail goes to 1.0
        sc_ = np.where(gated_b, np.abs(eff), -np.inf)
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
        RES[nm]["npred"].append(int(psig.sum()))
        ok = np.ones(Gk, bool)
        if j is not None: ok[j] = False
        ind = ok & (rl != 0) & ~np.isnan(rl)
        scored = psig & ind
        match = scored & (pl != 0) & ~np.isnan(pl) & (np.sign(pl) == np.sign(rl))
        n_real_, n_pred_, kk = int(rsig.sum()), int(scored.sum()), int(match.sum())
        if max(n_pred_, n_real_) > 0: RES[nm]["fid"].append(kk / max(n_pred_, n_real_))
        union = int((rsig | (psig & ok)).sum())
        RES[nm]["jac"].append(1.0 if union == 0 else int((rsig & psig & ok).sum()) / union)
        rlv, plv = rl[gs], pl[gs]
        idn = (~np.isnan(rlv)) & (rlv != 0)
        mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
        den = np.abs(rlv[np.isfinite(rlv)]).mean()
        if den > 0: RES[nm]["nmae"].append(np.abs(np.nan_to_num(plv) - rlv).mean() / den)
        o = np.lexsort((fr[gs], -np.abs(np.nan_to_num(plv, nan=-np.inf)), p2[gs], q2[gs], ~psig[gs]))
        nd_, nmt = np.cumsum(idn[o]), np.cumsum(mt[o])
        with np.errstate(invalid="ignore", divide="ignore"):
            pur = np.where(nd_ > 0, nmt / np.maximum(nd_, 1), np.nan)
        gg_ = np.flatnonzero((nd_ > 0) & (pur >= 0.9))
        RES[nm]["reach"].append((nd_[gg_[-1]] / nc) if gg_.size else 0.0)
    if (i_t + 1) % 10 == 0:
        print(f"#   {i_t+1}/{len(tgt)} [{time.time()-t0:.0f}s]", flush=True)

print(f"\n# {WHICH} DE GATE  ({len(RES['a12']['reach'])} perturbations, "
      f"tiered spread keep-top {KEEP_TOP}@{KEEP_SPREAD} tail {TAIL_SPREAD}, phi {PHI})")
print(f"  {'arm':<11}{'n_pred':>9}{'reach':>9}{'fid':>9}{'k/n_pred':>10}{'jac':>9}{'nmae':>9}")
base = None
for nm, _, _, _ in POL:
    r = RES[nm]
    if not r["reach"]: continue
    v = tuple(float(np.mean(r[k])) for k in ("npred", "reach", "fid", "jac", "nmae"))
    if base is None: base = v
    print(f"  {nm:<11}{v[0]:>9.0f}{v[1]:>9.4f}{v[2]:>9.4f}{'':>10}{v[3]:>9.4f}{v[4]:>9.4f}")

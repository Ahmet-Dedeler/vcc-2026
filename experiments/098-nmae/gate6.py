"""Full six-member gate for the arms 098 actually wants priced, with re-seeded twins.

Arms are direction-only (the whitening exponent), amplitude held at a12's `A` so nothing here
is an amplitude comparison.  Emission is production `emit_counts` with a12's tiered spread
(top-8000 by |effect| at 0.25, tail 1.0, phi 4), scanpy-identical rank-sum, BH 0.05.

`nmae` is reported BOTH realised (from the emitted cells) and analytic, at three reference-gate
depths, because gateladder.py shows the BH-0.05 cohort ranks amplitude arms with the opposite
sign to the board.  The realised-vs-analytic gap is the mechanism check (rule 3).

    N_PERTS=50 uv run python gate6.py <h1|orion>
"""
from __future__ import annotations
import hashlib, json, os, sys, time
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
import nmaekit as K
from msekit import TS, group_stats
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor, emit_counts
from vcc2026.multi_donor import MultiDonorModel
from vcc2026.whiten import Whitener
from de_kernel import ctrl_tiesums, wilcoxon_z
import anndata as ad

WHICH = sys.argv[1]; NP_ = int(os.environ.get("N_PERTS", "50"))
OUTDIR = {"h1": f"{ROOT}/experiments/095-clean-harness/out-h1test",
          "orion": f"{ROOT}/experiments/095-clean-harness/out"}[WHICH]
DEPTHS = (0.05, 0.35, 0.60)
ALPHA_V, EFFSCALE, MIN_CPM, ALPHA_DE = 1.08, 6.817, 5.0, 0.05
KEEP_TOP, KEEP_SPREAD, TAIL_SPREAD, PHI = 8000, 0.25, 1.0, 4.0

H = K.Harness(WHICH)
z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{WHICH}.npz", allow_pickle=True)
QQ, RL = z["q"], z["lfc"].astype(np.float64)
donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in ("k562_gwps", "orion_hek293t")}
model = MultiDonorModel(donors, H.genes, gene_norm=False, eb=False)
shared = set()
for n in model.names: shared |= set(model.target_genes[model.tcols[n]])
chal = set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str))
real_t = set(str(x) for x in H.rperts)
shared -= (chal if real_t <= chal else real_t)
ti = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
assert hashlib.sha256(wcols.tobytes()).hexdigest() == \
       json.load(open(f"{OUTDIR}/whitener_union.json"))["wcols_sha256"]
E0 = np.vstack([model.effect(str(x), 1.0, 0.0, weights={n: 1.0 for n in model.names})
                for x in H.tg]).astype(np.float64)
W = Whitener.load(f"{OUTDIR}/whitener_union.npz")


def dh_at(gam):
    E = E0.copy()
    if gam > 0:
        b0 = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
        E[:, wcols] = W.apply(E[:, wcols], gam, 0.0, renorm=False)
        E[:, wcols] *= b0 / (float(np.median(np.linalg.norm(E[:, wcols], axis=1))) or 1.0)
    E *= ALPHA_V * EFFSCALE / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
    dh = np.vstack([np.log1p(TS * model.profile_from_effect(E[i], H.bctrl, min_cpm=MIN_CPM,
                                                            gate_unexpressed=True)) - H.bctrl
                    for i in range(H.P)])
    for k, c in H.cols.items(): dh[k, c] = 0.0
    return dh


A0, _ = K.AXR(H, H.dh0)
cpm = np.zeros(H.G); cpm[H.gate] = H.ref_lin
d65, d00 = dh_at(0.65), dh_at(0.0)
ARMS = {"a12 (g0.65)": d65, "g0.45": dh_at(0.45), "g0.00": d00}
m200 = cpm > 200
dd = d65.copy(); dd[:, m200] = d00[:, m200]; ARMS["dewhiten CPM>200"] = dd
for nm in list(ARMS):
    A1, _ = K.AXR(H, ARMS[nm]); ARMS[nm] = ARMS[nm] * float(np.sqrt(A0 / A1))

# ---- harness DE setup -------------------------------------------------------------------
R = ad.read_h5ad({"h1": f"{ROOT}/data/harness/h1-test.h5ad",
                  "orion": f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"}[WHICH])
rlab = np.asarray(R.obs["target_gene"].astype(str).values)
Xr = R.X; Xr = np.asarray(Xr.todense() if hasattr(Xr, "todense") else Xr, dtype=np.float64)
def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))
cidx = np.flatnonzero(rlab == "non-targeting")
Cn = norm_log(Xr[cidx]); gate = H.gate
cs = np.ascontiguousarray(np.sort(Cn[:, gate].T, axis=1)); tie = ctrl_tiesums(cs)
ref_lin = np.expm1(Cn[:, gate]).mean(0); del Cn
Gk = int(gate.sum()); gk = H.genes[gate]
fr = np.empty(Gk, np.int64); fr[np.argsort(gk)] = np.arange(Gk)
gpos = {g: i for i, g in enumerate(gk)}
rp, rs_, rc, rj, rg = group_stats({"h1": f"{ROOT}/data/harness/h1-test.h5ad",
                                   "orion": f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"}[WHICH])
pp, ps, pc, pj, pg = group_stats(K.PRED[WHICH])
S = ps.sum(1).astype(np.float64); NC = int(round(float(np.median(pc))))
expr = np.expm1(H.bctrl); gated_b = (1e6 * expr / expr.sum() > MIN_CPM) & (H.bctrl > 0)
G = H.G

POL = [(nm, sd) for nm in ARMS for sd in (0, 1)]
RES = {(nm, sd): {k: [] for k in ("reach", "fid", "jac", "npred")} for nm, sd in POL}
NM = {(nm, sd, a): [] for nm, sd in POL for a in DEPTHS}
AN = {(nm, a): [] for nm in ARMS for a in DEPTHS}
LAN = {nm: H.lfc_hat(ARMS[nm][H.rows]) for nm in ARMS}
sel = list(zip(range(len(H.rows)), H.rows, H.rrows))[:NP_]
t0 = time.time()
for kk, (k, i, j) in enumerate(sel):
    tname = H.tg[i]
    ridx = np.flatnonzero(rlab == tname)
    if ridx.size == 0: continue
    Br = norm_log(Xr[ridx])
    zz = wilcoxon_z(cs, tie, np.ascontiguousarray(Br[:, gate].T))
    rq = false_discovery_control(2 * norm.sf(np.abs(zz)), method="bh")
    rl = np.log2((np.expm1(Br[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del Br
    jt = gpos.get(tname)
    gates = {}
    for a in DEPTHS:
        s = (rq < a) & np.isfinite(rl) & (rl != 0)
        if jt is not None: s[jt] = False
        gs = np.flatnonzero(s)
        gates[a] = gs if gs.size >= 10 else None
    rsig = rq < ALPHA_DE
    if jt is not None: rsig[jt] = False
    gs05 = np.flatnonzero(rsig)
    if gs05.size == 0: continue
    for nm in ARMS:
        for a in DEPTHS:
            if gates[a] is None: continue
            r_ = rl[gates[a]]
            AN[(nm, a)].append(np.abs(LAN[nm][k][gates[a]] - r_).sum() / np.abs(r_).sum())
    for nm, sd in POL:
        eff = ARMS[nm][i]
        prof = np.expm1(H.bctrl + eff); np.clip(prof, 0.0, None, out=prof)
        tot = prof.sum()
        if tot <= 0: continue
        prof = prof / tot
        sc_ = np.where(gated_b, np.abs(eff), -np.inf)
        order = np.lexsort((np.arange(G), -sc_))
        sv = np.full(G, TAIL_SPREAD); sv[order[:min(KEEP_TOP, int(gated_b.sum()))]] = KEEP_SPREAD
        v, col, ip = emit_counts(prof, np.full(NC, S[i] / NC), 0.0,
                                 np.random.default_rng(7919 + kk + 10_000 * sd),
                                 mode="quasi", spread=sv, overdispersion=PHI)
        blk = np.zeros((NC, G))
        for r2 in range(NC): blk[r2, col[ip[r2]:ip[r2+1]]] = v[ip[r2]:ip[r2+1]]
        nb = norm_log(blk); del blk
        z2 = wilcoxon_z(cs, tie, np.ascontiguousarray(nb[:, gate].T))
        p2 = 2 * norm.sf(np.abs(z2)); q2 = false_discovery_control(p2, method="bh")
        pl = np.log2((np.expm1(nb[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del nb
        psig = q2 < ALPHA_DE
        if jt is not None: psig[jt] = False
        RES[(nm, sd)]["npred"].append(int(psig.sum()))
        ok = np.ones(Gk, bool)
        if jt is not None: ok[jt] = False
        ind = ok & (rl != 0) & ~np.isnan(rl)
        scored = psig & ind
        match = scored & (pl != 0) & ~np.isnan(pl) & (np.sign(pl) == np.sign(rl))
        n_real_, n_pred_, kmt = int(rsig.sum()), int(scored.sum()), int(match.sum())
        if max(n_pred_, n_real_) > 0: RES[(nm, sd)]["fid"].append(kmt / max(n_pred_, n_real_))
        union = int((rsig | (psig & ok)).sum())
        RES[(nm, sd)]["jac"].append(1.0 if union == 0 else int((rsig & psig & ok).sum()) / union)
        for a in DEPTHS:
            if gates[a] is None: continue
            r_ = rl[gates[a]]
            NM[(nm, sd, a)].append(np.abs(np.nan_to_num(pl[gates[a]]) - r_).sum() / np.abs(r_).sum())
        o = np.lexsort((fr[gs05], -np.abs(np.nan_to_num(pl[gs05], nan=-np.inf)),
                        p2[gs05], q2[gs05], ~psig[gs05]))
        ndc, nmt = np.cumsum(ind[gs05][o]), np.cumsum(match[gs05][o])
        with np.errstate(invalid="ignore", divide="ignore"):
            pur = np.where(ndc > 0, nmt / np.maximum(ndc, 1), np.nan)
        gg_ = np.flatnonzero((ndc > 0) & (pur >= 0.9))
        RES[(nm, sd)]["reach"].append((ndc[gg_[-1]] / gs05.size) if gg_.size else 0.0)
    if (kk + 1) % 10 == 0: print(f"#  {kk+1}/{len(sel)} [{time.time()-t0:.0f}s]", flush=True)

print(f"\n===== {WHICH} SIX-MEMBER GATE ({len(RES[('a12 (g0.65)',0)]['reach'])} perturbations, "
      f"amplitude held at a12) =====")
hdr = f"  {'arm':<20}{'seed':>5}{'n_pred':>8}{'reach':>8}{'fid':>8}{'jac':>8}" + \
      "".join(f"{'nmae@'+str(a):>11}" for a in DEPTHS) + f"{'pds':>9}{'kappa':>9}"
print(hdr)
for nm in ARMS:
    A, X = K.AXR(H, ARMS[nm]); B = float(_across_pert_budget(ARMS[nm], H.cols, H.w))
    for sd in (0, 1):
        r = RES[(nm, sd)]
        print(f"  {nm:<20}{sd:>5}{np.mean(r['npred']):>8.0f}{np.mean(r['reach']):>8.4f}"
              f"{np.mean(r['fid']):>8.4f}{np.mean(r['jac']):>8.4f}"
              + "".join(f"{np.mean(NM[(nm,sd,a)]):>11.4f}" for a in DEPTHS)
              + (f"{K.pds(H, ARMS[nm]):>9.4f}{X/np.sqrt(B):>9.5f}" if sd == 0 else ""))
    print(f"  {'':<20}{'anal':>5}{'':>8}{'':>8}{'':>8}{'':>8}"
          + "".join(f"{np.mean(AN[(nm,a)]):>11.4f}" for a in DEPTHS))

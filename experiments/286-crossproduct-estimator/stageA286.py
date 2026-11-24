"""286 STAGE A -- HELD VIEWS on 280's GWCD4i bank (the cleanest test on disk).

    PYTHONDONTWRITEBYTECODE=1 uv run python stageA286.py

Fit.   S-hat, N-hat on the 500 NON-PANEL Stim48hr targets only (`extra280_X.npy` / `extra280_rows.npy`), views = the
       DONOR split: for each of the three disjoint 2+2 pairings, A = the mean of one pair's per-donor residual rows,
       B = the mean of the other pair's; the moments are averaged over the three pairings.  Held targets (the panel)
       never enter the fit.  A GUIDE split (guide-rank 0 vs guide-rank 1, pooled over donors) is fitted as the second,
       independent view definition.
       n = 500 << m = 17,471, so S-hat and N-hat live in a subspace of dimension p <= 2n: everything is computed
       EXACTLY in an orthonormal basis P of span([U; D]) and the estimator is the identity-with-floor outside it
       (G_phi = P G~_phi P^T + phi (I - P P^T), exact -- see RESULTS).

Test.  The 297 panel targets, all three disjoint 2+2 pairings x both directions: the prediction is built from TWO
       donors (E0 = the mean of their residual rows = 280's k = 2 object) and scored against the OTHER two donors'
       consensus, in 280's residual cosine and in the scorer's per-target rank statistic with all panel targets as
       distractors.  Primary state Stim48hr; Rest and Stim8hr secondary; HELD STATE = S-hat fitted on Stim48hr,
       scored on the held donors' Rest views.

KA.    (1) 280's Q2 `full_k2` per-target array reproduced at 0.0e+00 by this file's own row builder.
       (2) KA-M, the moment identity C = S + 2N at machine precision.
       (3) KA-E1id: E1 with N = 0 and phi = 1 equals E0 exactly.
Controls per cell: Cs shuffled pairing, Cw S-hat from another state, Cr the same spectrum on RANDOM directions.
"""
from __future__ import annotations

import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg as sla

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/286-crossproduct-estimator"
WORK = HERE / "work"
E280 = ROOT / "experiments/280-gwcd4i-donor-replicates/work"
sys.path.insert(0, str(E280))
sys.path.insert(0, str(HERE))
import est286 as E  # noqa: E402
from lib280 import CONDS, DORDER, DONORS, PAIRINGS, Bank, cosines, residualise  # noqa: E402

B_BOOT = 20_000
SEED = 20260906
PHI_MAIN = 0.10
LOG: list[str] = []
OUT: dict = dict(est286_sha256=E.sha256_self(), seed=SEED, b_boot=B_BOOT, ka=[], failures=[])
rng_boot = np.random.default_rng(SEED)


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def check(label: str, got: float, want: float, tol: float) -> bool:
    d = abs(float(got) - float(want))
    ok = d <= tol
    OUT["ka"].append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
    if not ok:
        OUT["failures"].append(f"{label}: |d| {d:.3e}")
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
    return ok


def boot(diff: np.ndarray) -> dict:
    d = diff[np.isfinite(diff)]
    idx = rng_boot.integers(0, d.size, size=(B_BOOT, d.size))
    mu = d[idx].mean(1)
    return dict(point=float(d.mean()), lo=float(np.percentile(mu, 2.5)), hi=float(np.percentile(mu, 97.5)),
                p_le0=float((mu <= 0).mean()), n=int(d.size), median=float(np.median(d)))


def fb(b: dict) -> str:
    return f"{b['point']:+.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}]"


def pds_rank(Pred: np.ndarray, Truth: np.ndarray) -> np.ndarray:
    """The scorer's per-target rank statistic: the fraction of the OTHER truths farther (lower cosine) than the own
    truth, ties at half; denominator n - 1."""
    Pn = Pred / np.maximum(np.linalg.norm(Pred, axis=1, keepdims=True), 1e-300)
    Tn = Truth / np.maximum(np.linalg.norm(Truth, axis=1, keepdims=True), 1e-300)
    Cm = Pn @ Tn.T
    own = np.diag(Cm).copy()
    n = Cm.shape[0]
    out = np.empty(n)
    for i in range(n):
        c = np.delete(Cm[i], i)
        out[i] = ((c < own[i]).sum() + 0.5 * (c == own[i]).sum()) / (n - 1)
    return out


# ============================================================================================= the bank
t0 = time.time()
bank = Bank(str(E280))
SC, CC = bank.scored, bank.cols_cd4
M = int(SC.sum())
o = bank.obs
tg = o[(o.gtype == "targeting") & o.keep_for_DE]
IDX = {(c, d, t): g.sort_values("guide") for (c, d, t), g in tg.groupby(["cond", "D", "pert"])}
P(f"bank: X {bank.X.shape}, mapped {len(CC)}, scored m = {M}; rss {rss():.2f} GiB, {time.time()-t0:.0f}s")


def donor_rows(cond: str, donor: str, targets, index=IDX, counts=None) -> np.ndarray:
    """Per-donor residual rows for `targets` (n, M): profile(sum of the target's guide pseudobulks) - the donor's
    control profile, centred on the column mean over `targets`, panel columns zeroed, restricted to `scored`."""
    b = bank.ctrl_profile(cond, [donor])[CC]
    X = np.zeros((len(targets), len(CC)))
    for k, t in enumerate(targets):
        S = counts(cond, donor, t) if counts is not None else bank.counts(index[(cond, donor, t)].i.values)
        X[k] = bank.profile(S)[CC] - b
    R = residualise(X)
    R[:, ~SC] = 0.0
    return R[:, SC]


# ============================================================================================= KA-1: 280's Q2 k=2
cond0 = "Stim48hr"
panel_t = [t for t in sorted(set(tg.pert)) if all((cond0, d, t) in IDX for d in DORDER)]
full = {d: donor_rows(cond0, d, panel_t) for d in DORDER}
import itertools  # noqa: E402
acc = []
for h in DORDER:
    for sub in itertools.combinations([d for d in DORDER if d != h], 2):
        acc.append(cosines(np.mean([full[d] for d in sub], axis=0), full[h]))
q2_mine = np.nanmean(np.vstack(acc), axis=0)
z = np.load(E280 / f"q2_{cond0}_keepDE.npz", allow_pickle=True)
assert [str(x) for x in z["targets"]] == panel_t, "Stage A: target order differs from 280's Q2"
check("KA-1 280's Q2 full_k2 per-target array reproduced (Stim48hr)", float(np.abs(q2_mine - z["full_k2"]).max()), 0.0, 0.0)
P(f"  Q2 k=2 (full) median cosine {np.median(q2_mine):+.5f} over {len(panel_t)} panel targets")
del acc, q2_mine
gc.collect()


# ============================================================================================= the fit set
zx = np.load(E280 / "obs280.npz", allow_pickle=True)
sel = np.load(E280 / "extra280_rows.npy")
XE = np.load(E280 / "extra280_X.npy", mmap_mode="r")
ext = pd.DataFrame({k: zx[k][sel] for k in ("donor", "cond", "gtype", "pert", "guide", "n_cells", "total", "keep_for_DE")})
for c in ("donor", "cond", "pert", "guide", "gtype"):
    ext[c] = ext[c].astype(str)
ext["D"] = ext.donor.map(DONORS)
ext["j"] = np.arange(len(ext))
ext = ext[ext.keep_for_DE.astype(bool)]
EIDX = {(d, t): g.sort_values("guide") for (d, t), g in ext.groupby(["D", "pert"])}
fit_t = sorted({t for t in set(ext.pert) if all((d, t) in EIDX for d in DORDER)})
P(f"fit set: {len(fit_t)} non-panel Stim48hr targets carried by all four donors "
  f"(of {len(set(ext.pert))} sampled); rss {rss():.2f} GiB")


def ext_counts(_cond, donor, t):
    r = np.sort(EIDX[(donor, t)].j.values)
    return np.asarray(XE[r], dtype=np.float64).sum(0)


def ext_counts_guide(donor, t, rank):
    g = EIDX[(donor, t)]
    if len(g) <= rank:
        return None
    r = np.asarray([np.sort(g.j.values)[rank]])
    return np.asarray(XE[r], dtype=np.float64).sum(0)


fit_full = {d: donor_rows(cond0, d, fit_t, counts=ext_counts) for d in DORDER}
P(f"  fit rows built: 4 x ({len(fit_t)}, {M}); rss {rss():.2f} GiB, {time.time()-t0:.0f}s")


# ============================================================================================= moments in a basis
def moments_from_views(A_list, B_list):
    """S-hat, N-hat averaged over the given (A, B) view pairs, and an orthonormal basis P of span([U; D]).
    Returns (P, S_red, N_red, audit) with S-hat = P S_red P^T and N-hat = P N_red P^T EXACTLY."""
    U = np.vstack([0.5 * (a + b) for a, b in zip(A_list, B_list)])
    D = np.vstack([a - b for a, b in zip(A_list, B_list)])
    n = A_list[0].shape[0]
    K = len(A_list)
    Z = np.vstack([U, D])
    # a rank-revealing orthonormal basis of span([U; D]) (the three pairings share one U, so the rank is < 2nK)
    _, sv, Vt = np.linalg.svd(Z, full_matrices=False)
    keep = sv > 1e-9 * sv.max()
    Pb = np.ascontiguousarray(Vt[keep].T)                      # (m, p)
    Ur, Dr = U @ Pb, D @ Pb
    G_U = Ur.T @ Ur
    G_D = Dr.T @ Dr
    S_red = (G_U - 0.25 * G_D) / (n * K)
    N_red = G_D / (4.0 * n * K)
    C_red = (G_U + 0.25 * G_D) / (n * K)
    aud = dict(n=n, n_pairs=K, p=int(Pb.shape[1]), m=int(Pb.shape[0]),
               ka_exact=float(np.linalg.norm(C_red - S_red - 2 * N_red) / np.linalg.norm(C_red)),
               literal=float(np.linalg.norm(C_red - S_red - N_red) / np.linalg.norm(C_red)),
               trace_S=float(np.trace(S_red)), trace_N=float(np.trace(N_red)), trace_C=float(np.trace(C_red)),
               signal_share_of_trace=float(np.trace(S_red) / np.trace(C_red)))
    return Pb, S_red, N_red, aud


def donor_views(rows: dict, pairings=PAIRINGS, perm=None):
    A, B = [], []
    for (p1, p2) in pairings:
        a = np.mean([rows[d] for d in p1], axis=0)
        b = np.mean([rows[d] for d in p2], axis=0)
        A.append(a)
        B.append(b if perm is None else b[perm])
    return A, B


# the primary fit (donor split, three pairings) and its shuffled-pairing control
PB, S_red, N_red, aud_fit = moments_from_views(*donor_views(fit_full))
check("KA-M exact C = S + 2N (Stage A fit)", aud_fit["ka_exact"], 0.0, 1e-9)
P(f"  fit moments: p = {aud_fit['p']} (<= 2 n K = {6*len(fit_t)}), tr S / tr C {aud_fit['signal_share_of_trace']:+.4f}, "
  f"literal ||C-S-N||/||C|| {aud_fit['literal']:.4f}")
perm = E.derangement(len(fit_t))
PBs, S_reds, N_reds, aud_shuf = moments_from_views(*donor_views(fit_full, perm=perm))
P(f"  shuffled-pairing control: p = {aud_shuf['p']}, tr S / tr C {aud_shuf['signal_share_of_trace']:+.4f}")
lead_real = float(np.linalg.eigvalsh(S_red).max())
lead_shuf = float(np.linalg.eigvalsh(S_reds).max())
P(f"  KA-shuf leading eigenvalue: real {lead_real:.5e}, shuffled {lead_shuf:.5e}, ratio {lead_real/lead_shuf:.3f}")
OUT["fit"] = dict(donor=aud_fit, shuffled=aud_shuf, lead_real=lead_real, lead_shuf=lead_shuf,
                  lead_ratio=lead_real / lead_shuf, n_fit_targets=len(fit_t), perm_seed=E.DERANGE_SEED)

# the guide split (the second view definition)
guide_ok = [t for t in fit_t if all(ext_counts_guide(d, t, 1) is not None for d in DORDER)]
if len(guide_ok) >= 100:
    gA = {d: donor_rows(cond0, d, guide_ok, counts=lambda c, d_, t, r=0: ext_counts_guide(d_, t, r)) for d in DORDER}
    gB = {d: donor_rows(cond0, d, guide_ok, counts=lambda c, d_, t: ext_counts_guide(d_, t, 1)) for d in DORDER}
    A_g = [np.mean([gA[d] for d in DORDER], axis=0)]
    B_g = [np.mean([gB[d] for d in DORDER], axis=0)]
    PBg, S_redg, N_redg, aud_g = moments_from_views(A_g, B_g)
    P(f"  guide split ({len(guide_ok)} targets with >= 2 keep_for_DE guides in all four donors): p = {aud_g['p']}, "
      f"tr S / tr C {aud_g['signal_share_of_trace']:+.4f}")
    OUT["fit"]["guide"] = aud_g
    del gA, gB
else:
    PBg = None
    P(f"  guide split: only {len(guide_ok)} targets carry two keep_for_DE guides in all four donors -- SKIPPED")
    OUT["fit"]["guide"] = dict(skipped=True, n=len(guide_ok))
gc.collect()


# ============================================================================================= the estimator
def make_filter(Pb, S_red, N_red, J: float, m_full: int):
    """The reduced factors and the eps of the FULL-space definition (eps = EPS_REL * tr(N) / m)."""
    p = Pb.shape[1]
    eps_rel_red = E.EPS_REL * p / m_full            # keeps eps = EPS_REL * tr(N) / m_full, the full-space definition
    lam, V = np.linalg.eigh(S_red)
    Sp = (V * np.clip(lam, 0.0, None)) @ V.T
    L, Q, k, aud = E.wiener_factors(Sp, N_red, J, eps_rel=eps_rel_red)
    aud["clipped_energy"] = float(np.abs(lam[lam < 0]).sum() / np.abs(lam).sum())
    aud["share_negative"] = float((lam < 0).mean())
    aud["p"] = p
    return dict(Pb=Pb, L=L, Q=Q, k=k, lam=lam, V=V, aud=aud, p=p)


def apply_filter(F, X: np.ndarray, phi: float, *, mirror: bool = False, rank: int | None = None) -> np.ndarray:
    """G_phi X^T with G_phi = P G~_phi P^T + phi (I - P P^T) (exact; see the module docstring)."""
    Pb = F["Pb"]
    Xin = X @ Pb                                            # (n, p) coordinates
    if rank is None:
        Y = E.apply_e1(Xin, F["L"], F["Q"], F["k"], phi, mirror=mirror)
    else:
        Y = E.apply_e2(Xin, F["L"], F["lam"], F["V"], rank)
        phi = 0.0
    out = Y @ Pb.T
    if phi:
        out = out + phi * (X - Xin @ Pb.T)
    return out


F_main = make_filter(PB, S_red, N_red, 0.5, M)      # J = 0.5: the row filtered is a TWO-donor mean, N-hat is
P(f"  E1 factors: {F_main['aud']}")                 # calibrated for the FOUR-donor mean -> N/J = 2 N
F_shuf = make_filter(PBs, S_reds, N_reds, 0.5, M)
rng_r = np.random.default_rng(E.RAND_SEED)
Pr, _ = np.linalg.qr(rng_r.standard_normal((M, PB.shape[1])))
F_rand = make_filter(np.ascontiguousarray(Pr), S_red, N_red, 0.5, M)   # same spectrum, RANDOM directions
del Pr
gc.collect()

# KA-E1id: phi = 1 raises every gain to 1, so G_1 = I exactly (on span(P) and on its complement) -- the identity row
Xt = np.mean([full[d] for d in ("D1", "D2")], axis=0)[:40]
rel = float(np.abs(apply_filter(F_main, Xt, 1.0) - Xt).max() / np.abs(Xt).max())
check("KA-E1id E1 at phi = 1 (every gain floored to 1) equals E0, relative", rel, 0.0, 1e-4)
del Xt
gc.collect()


# ============================================================================================= the held-view test
def evaluate(cond_pred: str, cond_truth: str, filters: dict, tag: str) -> dict:
    rows_p = full if cond_pred == cond0 else {d: donor_rows(cond_pred, d, panel_c[cond_pred]) for d in DORDER}
    rows_t = rows_p if cond_truth == cond_pred else {d: donor_rows(cond_truth, d, panel_c[cond_pred]) for d in DORDER}
    cells = {}
    for name, (F, phi, rank, mirror) in filters.items():
        cs, pd_ = [], []
        for (p1, p2) in PAIRINGS:
            for pred_pair, truth_pair in ((p1, p2), (p2, p1)):
                X0 = np.mean([rows_p[d] for d in pred_pair], axis=0)
                T = np.mean([rows_t[d] for d in truth_pair], axis=0)
                Xf = X0 if F is None else apply_filter(F, X0, phi, mirror=mirror, rank=rank)
                cs.append(cosines(Xf, T))
                pd_.append(pds_rank(Xf, T))
        cells[name] = dict(cos=np.nanmean(np.vstack(cs), axis=0), pds=np.nanmean(np.vstack(pd_), axis=0))
    base = cells["E0"]
    rec = {}
    for name, v in cells.items():
        rec[name] = dict(cos_median=float(np.nanmedian(v["cos"])), cos_mean=float(np.nanmean(v["cos"])),
                         pds_mean=float(np.nanmean(v["pds"])),
                         d_cos=boot(v["cos"] - base["cos"]), d_pds=boot(v["pds"] - base["pds"]))
        P(f"    {tag:26s} {name:14s} cos median {rec[name]['cos_median']:+.5f}  d cos {fb(rec[name]['d_cos'])}  "
          f"pds {rec[name]['pds_mean']:.5f}  d pds {fb(rec[name]['d_pds'])}")
    np.savez_compressed(WORK / f"stageA-{tag}.npz", targets=np.asarray(panel_c[cond_pred]),
                        **{f"{k}|cos": v["cos"] for k, v in cells.items()},
                        **{f"{k}|pds": v["pds"] for k, v in cells.items()})
    return rec


panel_c = {c: [t for t in sorted(set(tg.pert)) if all((c, d, t) in IDX for d in DORDER)] for c in CONDS}
FILTERS = {"E0": (None, 0.0, None, False)}
for phi in E.PHIS:
    FILTERS[f"E1:phi{phi:g}"] = (F_main, float(phi), None, False)
for r in E.RANKS:
    FILTERS[f"E2:r{r}"] = (F_main, 0.0, int(r), False)
FILTERS[f"E1:phi{PHI_MAIN:g}:mirror"] = (F_main, PHI_MAIN, None, True)
FILTERS[f"E1:phi{PHI_MAIN:g}:shufpair"] = (F_shuf, PHI_MAIN, None, False)
FILTERS[f"E1:phi{PHI_MAIN:g}:randdir"] = (F_rand, PHI_MAIN, None, False)

P(f"\n== held VIEWS: predict from two donors, score against the other two ({len(panel_c[cond0])} panel targets) ==")
OUT["heldview"] = {cond0: evaluate(cond0, cond0, FILTERS, f"heldview-{cond0}")}
P("\n== held STATE: S-hat is fitted on Stim48hr views only, so Rest and Stim8hr ARE the held-state readings ==")
for c in ("Rest", "Stim8hr"):
    P(f"  -- held state {c} ({len(panel_c[c])} targets)")
    OUT["heldview"][c] = evaluate(c, c, FILTERS, f"heldstate-{c}")
OUT["heldstate"] = OUT["heldview"]["Rest"]
OUT["cw_note"] = ("280's extra280 extraction covers Stim48hr only, so a wrong-STATE S-hat cannot be fitted on "
                  "non-panel targets; Cw is realised as (i) the cross-state application (S-hat fitted on Stim48hr, "
                  "scored on the held donors' Rest / Stim8hr views) and (ii) `:randdir`, the SAME spectrum on random "
                  "directions -- a stronger control than a wrong-source covariance, since it matches the energy and "
                  "the eigenvalue profile exactly and differs only in the directions.")

OUT["seconds"] = round(time.time() - t0, 1)
OUT["peak_rss_gib"] = rss()
(WORK / "stageA286.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "stageA286.log").write_text("\n".join(LOG) + "\n")
P(f"\n[Stage A] done {(time.time()-t0)/60:.1f} min, peak rss {rss():.2f} GiB, KA "
  f"{sum(1 for k in OUT['ka'] if k['ok'])}/{len(OUT['ka'])}, failures {OUT['failures']}")

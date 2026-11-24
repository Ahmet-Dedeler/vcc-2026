"""177 — `mse` pricing for the RESAMPLING emission law, on the REAL challenge panel.

Verbatim reuse of `experiments/171-reach-prefix/claim171.py` (which is 164's, which is 096's
`mech_claim_grid2`), with 164's three known-answer gates kept unchanged and the resampling arms
substituted for the spread-schedule ones.

WHY THIS IS NEEDED even though every real-depth arm holds the pseudobulk bit-identical
(`realize177.py` records one sha256 per (arm, target) and all of them match): `claim` is the
delete-one jackknife of OUR OWN emitted noise, which is exactly the quantity an emission law
moves.  164 put a15 at `claim` 0.24581 against the reference cap K = 0.24653 -- 99.7% of the
ceiling -- and 171 section 5 measured the knee at `s^2 B = 0.23268`, i.e. a15 sits **+0.01313
ABOVE** it, and `mse` is exactly flat in `claim` above the knee.  The `resample_mix` docstring
claims `C_hat` rises to ~1.03 of the real-cell value; if so the arm is `mse`-free.  Measure it,
do not assert it (AGENTS.md 11.1 rule 12, the a19 failure mode).

The resample arms need the challenge's own CONTROL CELLS, which `mech_alpha_chal` does not
carry, so `data/challenge/context_{A,B,C}.h5ad` is read directly; its gene axis is verified
identical to `v11-stats.npz`'s 18,533.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import os
import sys

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, HERE)
from signal177 import fracpair  # noqa: E402
import vcc2026.donor_transfer as _dt  # noqa: E402
from realize177 import _fast_repair_totals  # noqa: E402

_dt._repair_totals = _fast_repair_totals
from cell_eval2.moments import jackknife_correction  # noqa: E402
import mech_alpha_chal as M  # noqa: E402

D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
NSUB = 16
MIN_CPM, KEEP_TOP = 10.0, 8000
NPOOL = 800
CV_Q, CPM_FLOOR, SCHED_REPS = 60.0, 15.0, 2

# name, law, tail_spread(None=flat at keep_spread), keep_top, keep_spread, phi_k, phi_t, jitter
POLICIES = [
    ("KNOWN flat 0.25/phi0", "quasi", None, 0, 0.25, 0.0, 0.0, False),
    ("KNOWN flat 1.00/phi0", "quasi", None, 0, 1.00, 0.0, 0.0, False),
    ("KNOWN flat 0.25/phi4 (a8)", "quasi", None, 0, 0.25, 4.0, 4.0, False),
    ("a15 AS SHIPPED", "quasi", 1.0, KEEP_TOP, 0.25, 4.0, 4.0, False),
    ("rs_all  (resample, all genes)", "rs", 1.0, KEEP_TOP, 0.25, 4.0, 4.0, False),
    ("rs_sched (+ fracpair x2)", "rs", 1.0, KEEP_TOP, 0.25, 4.0, 4.0, True),
    ("hy_kt2000 (top2000 parametric)", "hy", 1.0, 2000, 0.25, 4.0, 4.0, False),
]


def spread_phi(e, bctrl, min_cpm, tail, ktop, kspread, phi_keep, phi_tail):
    """Exactly `emission_spread_vector` + `emit_counts`'s per-gene phi (164's function)."""
    n = e.size
    if tail is None:                                   # flat law
        return np.full(n, float(kspread)), np.full(n, float(phi_keep)), np.zeros(n, bool)
    sv = np.full(n, float(tail))                       # trap: initialised to tail_spread
    pv = np.full(n, float(phi_tail))
    top = np.zeros(n, bool)
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > min_cpm) & (bctrl > 0)
    score = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k > 0:
        order = np.lexsort((np.arange(n), -score))[:k]
        sv[order] = float(kspread)
        pv[order] = float(phi_keep)
        top[order] = True
    return sv, pv, top


def load_controls():
    """400-cell donor pools and the pooled control profile, per context."""
    out = {}
    for c in M.CTX:
        Aa = ad.read_h5ad(f"{ROOT}/data/challenge/context_{c}.h5ad")
        Xc = Aa.X.tocsr()
        rng = np.random.default_rng(77)
        idx = np.sort(rng.choice(Xc.shape[0], min(NPOOL, Xc.shape[0]), replace=False))
        pool = np.asarray(Xc[idx].todense(), dtype=np.float64)
        e = np.expm1(D[c]["bctrl"])
        out[c] = (pool, np.clip(e, 0.0, None) / np.clip(e, 0.0, None).sum())
        del Aa, Xc
        print(f"  control pool {c}: {pool.shape}, median lib {np.median(pool.sum(1)):.0f}",
              flush=True)
    return out


def main() -> None:
    print(f"reference cap K = Sum_p w_p C_p = {K:.5f}   sumD = {sD:.5f}   "
          f"A = {A:.5f}   B = {B:.5f}   X = {X:.6f}   Kc = {Kc:.6f}")
    print(f"gate: pseudobulk CPM > {MIN_CPM}\n")
    pools = load_controls()
    # the expression-compensating schedule, per context, exactly as `realize177.py` builds it
    sched = {}
    for c in M.CTX:
        pool, crel = pools[c]
        lib = pool.sum(1)
        lib[lib == 0] = 1.0
        dcpm = pool * (1e6 / lib[:, None])
        e1, e2 = dcpm.mean(0), (dcpm ** 2).mean(0)
        cv = np.sqrt(np.maximum(e2 - e1 ** 2, 0.0)) / np.maximum(e1, 1e-12)
        ok = np.flatnonzero((e1 > CPM_FLOOR) & (D[c]["bctrl"] > 0))
        Tq = float(np.percentile(cv[ok], CV_Q))
        sg = np.zeros(e1.size)
        sg[ok] = np.clip(np.sqrt(np.maximum(Tq ** 2 - cv[ok] ** 2, 0.0)), 0.0, 1.0)
        sched[c] = sg
        del dcpm
    print()
    print(f"{'policy':<32}{'gated':>7}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'cap?':>6}"
          f"{'rho':>6}{'s*':>8}{'mse@s*':>9}{'score':>8}{'d overall':>11}")
    out = {}
    for name, law, tail, ktop, kspread, phi_k, phi_t, jit in POLICIES:
        rows, ng = [], []
        for c in M.CTX:
            d = D[c]
            pool, crel = pools[c]
            idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
            vals = []
            for j in range(idx.size):
                e = d["dh"][idx[j]]
                prof = np.expm1(d["bctrl"] + e)
                np.clip(prof, 0, None, out=prof)
                prof /= prof.sum()
                sv, pv, top = spread_phi(e, d["bctrl"], MIN_CPM, tail, ktop, kspread,
                                         phi_k, phi_t)
                if j == 0:
                    expr = np.expm1(d["bctrl"])
                    cc = 1e6 * expr / expr.sum()
                    ng.append(int(((cc > MIN_CPM) & (d["bctrl"] > 0)).sum()))
                rng = np.random.default_rng(3000 + j)
                S = float(d["S"][idx[j]])
                if law == "quasi":
                    libs = np.full(M.NCELL, S / M.NCELL)
                    v, col, ip = M.emit_counts(prof, libs, 0.0, rng, mode="quasi",
                                               spread=sv, overdispersion=pv)
                else:
                    sel = rng.choice(pool.shape[0], M.NCELL, replace=False)
                    Cc = pool[sel]
                    if jit:
                        Cc = fracpair(Cc, sched[c], rng, reps=SCHED_REPS)
                    dlib = Cc.sum(1)
                    dlib[dlib == 0] = 1.0
                    Ld = float(int(min(np.median(dlib), 1_000_000)))
                    pmask = top if law == "hy" else np.zeros(d["G"], bool)
                    v, col, ip = M.emit_counts(
                        prof * (S / (M.NCELL * Ld)), dlib, 0.0, rng, mode="resample_mix",
                        spread=sv, overdispersion=pv, control_counts=Cc,
                        control_profile=crel, param_mask=pmask)
                    del Cc
                Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64),
                                 ip.astype(np.int64)), shape=(M.NCELL, d["G"]))
                vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1,
                                                 M.TS, chunk=64)[0])
            cp = float(cpn[c].mean())
            rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals, cp)))))
        ch = float(np.mean([r[0] for r in rows]))
        cp = float(np.mean([r[1] for r in rows]))
        claim = float(np.mean([300 * r[2] / 18532 for r in rows]))
        rho = min(1.0, B / claim)
        s = np.sqrt(claim / B)
        mse_opt = 1 + (s * s * (A - B) + Kc - 2 * s * X) / sD
        sc = min(max((0.987862879 - mse_opt) / 0.950162566, 0), 1)
        out[name] = dict(claim=claim, mse_opt=mse_opt, sc=sc, s=s, ch=ch)
        base = out.get("a15 AS SHIPPED")
        dov = (sc - base["sc"]) / 6 if base else 0.0
        print(f"{name:<32}{int(np.mean(ng)):>7}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
              f"{'YES' if claim >= K - 1e-9 else 'no':>6}{rho:>6.2f}{s:>8.4f}"
              f"{mse_opt:>9.4f}{sc:>8.4f}{dov:>+11.6f}", flush=True)

    print("\n--- known-answer gates (164 section 1, from 096's own logs) ---")
    G = [("KNOWN flat 0.25/phi0", 0.11860, 0.02),
         ("KNOWN flat 1.00/phi0", 0.19216, 0.02),
         ("KNOWN flat 0.25/phi4 (a8)", 0.16478, 0.03)]
    bad = []
    for nm, want, tol in G:
        got = out[nm]["claim"]
        ok = abs(got - want) / want < tol
        print(f"  {nm:<32} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(nm)
    a15 = out["a15 AS SHIPPED"]
    ok15 = abs(a15["claim"] - 0.24581) / 0.24581 < 0.02
    print(f"  {'a15 AS SHIPPED (164 section 5)':<32} claim {a15['claim']:.5f}  "
          f"expect 0.24581  ({100*(a15['claim']-0.24581)/0.24581:+.2f}%)  "
          f"{'PASS' if ok15 else 'FAIL'}")
    if not ok15:
        bad.append("a15")
    if bad:
        sys.exit(f"ABORT: claim known-answer gate failed for {bad}")

    print("\n--- held at a15's SHIPPED VCC_ALPHA 0.982 (171 section 5: the knee is s^2 B =")
    print("    0.23268 and a15 sits +0.01313 above it, so `mse` is flat in `claim` up there)")
    print(f"{'arm':<32}{'claim':>9}{'d claim':>10}{'d raw mse':>11}"
          f"{'raw mse':>10}{'score_mse':>11}{'d overall':>11}")
    s15 = 0.982 / 1.0796
    mse15 = 1 + (s15 * s15 * (A - B) + Kc - 2 * s15 * X) / sD
    knee = s15 * s15 * B
    for nm in [p[0] for p in POLICIES]:
        cl = out[nm]["claim"]
        dc = min(cl, knee) - min(a15["claim"], knee)     # `min(s^2 B, claim)` saturates
        raw = mse15 - dc / sD
        sc = min(max((0.987862879 - raw) / 0.950162566, 0), 1)
        sc15 = min(max((0.987862879 - mse15) / 0.950162566, 0), 1)
        print(f"{nm:<32}{cl:>9.5f}{dc:>+10.5f}{-dc/sD:>+11.5f}"
              f"{raw:>10.5f}{sc:>11.5f}{(sc-sc15)/6:>+11.6f}")
    print(f"\nknee s^2 B = {knee:.5f}; every arm at or above it is exactly `mse`-free.")


if __name__ == "__main__":
    main()

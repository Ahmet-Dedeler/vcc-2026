"""164 — `mse` pricing for every call-set arm.  AGENTS.md 11.1 rule 12 (the a19 failure).

`claim = Sum_p w_p min(C_hat_p, C_p)` is the ONLY term in the `mse` numerator that a
(spread, phi) change can move: the quasi emitter repairs each gene's panel total exactly, so
`A`, `B`, `X` -- and therefore `pds` and `nmae` -- are bit-invariant across these arms
(`donor_transfer.py:290-300`).  `claim` enters with a MINUS sign, so losing claim RAISES raw
`mse` (worse).

    mse - 1 = ( s^2*(A - B) + Kc - 2*s*X ) / sumD   at the claim-matched amplitude s
    (the shipped form; `rho = min(1, B/claim) = 1` throughout, which is re-checked below)

Machinery is `experiments/096-mse-floor/mech_claim_grid2.py` + `mech_split_policy.py`,
reused verbatim except that the per-gene (spread, phi) vector is built the way PRODUCTION
builds it -- `emission_spread_vector`: gate on PSEUDOBULK CPM > --min-cpm, rank the gated
genes by |effect|, top `keep_top` get `keep_spread`, and THE WHOLE VECTOR IS INITIALISED TO
`tail_spread` so sub-gate genes sit in the tail too.

KNOWN-ANSWER GATES, all three from 096's own logs, any of which could fail:
    (0.25, phi 0)  claim 0.11860   work/claim_grid2.log
    (1.00, phi 0)  claim 0.19216   work/claim_grid.log
    (0.25, phi 4)  claim 0.16478 / 0.16676 (a8 as shipped, work/board_solve.npz)

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.sparse import csr_matrix

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction  # noqa: E402
import mech_alpha_chal as M  # noqa: E402

D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
NSUB = 16
MIN_CPM, KEEP_TOP = 10.0, 8000

# name, min_cpm, tail_spread(None=flat), keep_top, keep_spread, phi_keep, phi_tail
POLICIES = [
    ("KNOWN flat 0.25/phi0", MIN_CPM, None, 0, 0.25, 0.0, 0.0),
    ("KNOWN flat 1.00/phi0", MIN_CPM, None, 0, 1.00, 0.0, 0.0),
    ("KNOWN flat 0.25/phi4 (a8)", MIN_CPM, None, 0, 0.25, 4.0, 4.0),
    ("a15 AS SHIPPED", MIN_CPM, 1.0, KEEP_TOP, 0.25, 4.0, 4.0),
    ("armA tail phi 0", MIN_CPM, 1.0, KEEP_TOP, 0.25, 4.0, 0.0),
    ("armB tail spread 0.5", MIN_CPM, 0.5, KEEP_TOP, 0.25, 4.0, 4.0),
    ("armC a8 flat 0.25/phi4", MIN_CPM, None, 0, 0.25, 4.0, 4.0),
    ("armB' tail 0.75", MIN_CPM, 0.75, KEEP_TOP, 0.25, 4.0, 4.0),
    ("armA' tail phi 1", MIN_CPM, 1.0, KEEP_TOP, 0.25, 4.0, 1.0),
    ("a15 @ --min-cpm 5 (SPREAD ONLY)", 5.0, 1.0, KEEP_TOP, 0.25, 4.0, 4.0),
]


def spread_phi(e, bctrl, min_cpm, tail, ktop, kspread, phi_keep, phi_tail):
    """Exactly `emission_spread_vector` + `emit_counts`'s per-gene phi."""
    n = e.size
    if tail is None:                                   # flat law
        return np.full(n, float(kspread)), np.full(n, float(phi_keep))
    sv = np.full(n, float(tail))                       # trap: initialised to tail_spread
    pv = np.full(n, float(phi_tail))
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > min_cpm) & (bctrl > 0)
    score = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k > 0:
        order = np.lexsort((np.arange(n), -score))[:k]
        sv[order] = float(kspread)
        pv[order] = float(phi_keep)
    return sv, pv


def main() -> None:
    print(f"reference cap K = Sum_p w_p C_p = {K:.5f}   sumD = {sD:.5f}   "
          f"A = {A:.5f}   B = {B:.5f}   X = {X:.6f}   Kc = {Kc:.6f}")
    print(f"gate: pseudobulk CPM > min_cpm; keep_top = {KEEP_TOP}\n")
    print(f"{'policy':<34}{'gated':>7}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'cap?':>6}"
          f"{'rho':>6}{'s*':>8}{'ALPHA*':>9}{'mse@s*':>9}{'mse@a15':>9}{'score':>8}{'overall':>9}")
    out = {}
    for name, mcpm, tail, ktop, kspread, phi_k, phi_t in POLICIES:
        rows, ng = [], []
        for c in M.CTX:
            d = D[c]
            idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
            vals = []
            for j in range(idx.size):
                e = d["dh"][idx[j]]
                prof = np.expm1(d["bctrl"] + e)
                np.clip(prof, 0, None, out=prof)
                prof /= prof.sum()
                sv, pv = spread_phi(e, d["bctrl"], mcpm, tail, ktop, kspread, phi_k, phi_t)
                if j == 0:
                    expr = np.expm1(d["bctrl"])
                    cc = 1e6 * expr / expr.sum()
                    ng.append(int(((cc > mcpm) & (d["bctrl"] > 0)).sum()))
                libs = np.full(M.NCELL, d["S"][idx[j]] / M.NCELL)
                v, col, ip = M.emit_counts(
                    prof, libs, 0.0, np.random.default_rng(3000 + j), mode="quasi",
                    spread=sv, overdispersion=pv)
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
        s = np.sqrt(claim / B)                       # claim-matched optimal amplitude
        mse_opt = 1 + (s * s * (A - B) + Kc - 2 * s * X) / sD
        # a15 ships VCC_ALPHA 0.982 against v11's 1.0796 -> s_a15 = 0.982/1.0796.
        # 096 (b) established that this IS the per-context argmin at a15's claim, i.e. a15
        # sits ON the cliff, so an arm that LOWERS claim cannot compensate with amplitude.
        s15 = 0.982 / 1.0796
        mse_a15 = 1 + (s15 * s15 * (A - B) + Kc - 2 * s15 * X) / sD
        sc = min(max((0.987862879 - mse_opt) / 0.950162566, 0), 1)
        out[name] = dict(claim=claim, mse_opt=mse_opt, sc=sc, s=s, ch=ch)
        print(f"{name:<34}{int(np.mean(ng)):>7}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
              f"{'YES' if claim >= K - 1e-9 else 'no':>6}{rho:>6.2f}{s:>8.4f}"
              f"{1.0796*s:>9.4f}{mse_opt:>9.4f}{mse_a15:>9.4f}{sc:>8.4f}{sc/6:>9.5f}",
              flush=True)

    # -------- gates --------
    print("\n--- known-answer gates ---")
    G = [("KNOWN flat 0.25/phi0", 0.11860, 0.02),
         ("KNOWN flat 1.00/phi0", 0.19216, 0.02),
         ("KNOWN flat 0.25/phi4 (a8)", 0.16478, 0.03)]
    bad = []
    for nm, want, tol in G:
        got = out[nm]["claim"]
        ok = abs(got - want) / want < tol
        print(f"  {nm:<30} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(nm)
    if bad:
        sys.exit(f"ABORT: claim known-answer gate failed for {bad}")

    # -------- the cliff, stated as a number --------
    print("\n--- the rho cliff: a15 sits ON the per-context argmin, so no amplitude rescue ---")
    base = out["a15 AS SHIPPED"]
    print(f"{'arm':<34}{'claim':>9}{'d claim':>10}{'raw mse @ its own s*':>22}"
          f"{'d raw mse':>11}{'d score_mse':>13}{'d overall':>11}")
    for nm in ("a15 AS SHIPPED", "armA tail phi 0", "armB tail spread 0.5",
               "armC a8 flat 0.25/phi4", "armB' tail 0.75", "armA' tail phi 1",
               "a15 @ --min-cpm 5 (SPREAD ONLY)"):
        o = out[nm]
        print(f"{nm:<34}{o['claim']:>9.5f}{o['claim']-base['claim']:>+10.5f}"
              f"{o['mse_opt']:>22.5f}{o['mse_opt']-base['mse_opt']:>+11.5f}"
              f"{o['sc']-base['sc']:>+13.5f}{(o['sc']-base['sc'])/6:>+11.6f}")
    print("\n  Note: `d overall` above already grants each arm its OWN optimal amplitude.")
    print("  Held at a15's shipped VCC_ALPHA the loss is larger, because a15's alpha was")
    print("  chosen for a15's claim.  Delta claim / sumD is the first-order cost:")
    for nm in ("armA tail phi 0", "armB tail spread 0.5", "armC a8 flat 0.25/phi4"):
        dc = out[nm]["claim"] - base["claim"]
        print(f"    {nm:<32} d claim {dc:+.5f} -> d raw mse {-dc/sD:+.5f} "
              f"(no compensation) -> d overall {(-(-dc/sD)/0.950162566)/6:+.6f}")


if __name__ == "__main__":
    main()

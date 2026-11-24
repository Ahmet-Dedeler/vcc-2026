"""171 — `mse` pricing for the reach-ordering arms, on the REAL challenge panel.

Verbatim reuse of `experiments/164-callset-artefact/claim164.py`'s machinery (which is itself
096's `mech_claim_grid2` + `mech_split_policy`), with 164's three known-answer gates kept
unchanged and this experiment's arms substituted.

WHY THIS IS NEEDED even though every arm here holds the pseudobulk bit-identical: `claim`
is the delete-one jackknife of OUR OWN emitted noise, which is exactly the quantity a
(spread, phi) policy moves.  164 measured a15 at `claim` 0.24581 against the reference cap
K = 0.24653 -- 99.7% of the ceiling, with `score_mse` clamping to zero the moment raw `mse`
passes 0.987862879.  So `mse` is the binding constraint on this whole family and it must be
computed, not asserted (AGENTS.md 11.1 rule 12, the a19 failure mode).

The arms below differ from 164's in DIRECTION: 164 only lowered the tail's dispersion, which
loses claim.  These lower `keep_spread` (which loses the keep block's claim) or shrink
`keep_top` (which moves genes INTO the high-dispersion tail and gains claim).  Whether the
family can be made claim-neutral is the question.

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

# name, min_cpm, tail_spread(None=flat at keep_spread), keep_top, keep_spread, phi_k, phi_t
POLICIES = [
    ("KNOWN flat 0.25/phi0", MIN_CPM, None, 0, 0.25, 0.0, 0.0),
    ("KNOWN flat 1.00/phi0", MIN_CPM, None, 0, 1.00, 0.0, 0.0),
    ("KNOWN flat 0.25/phi4 (a8)", MIN_CPM, None, 0, 0.25, 4.0, 4.0),
    ("a15 AS SHIPPED", MIN_CPM, 1.0, KEEP_TOP, 0.25, 4.0, 4.0),
    ("ks000  keep_spread 0.00", MIN_CPM, 1.0, KEEP_TOP, 0.00, 4.0, 4.0),
    ("ks010  keep_spread 0.10", MIN_CPM, 1.0, KEEP_TOP, 0.10, 4.0, 4.0),
    ("kt2000_ks000", MIN_CPM, 1.0, 2000, 0.00, 4.0, 4.0),
    ("kt2000_ks025", MIN_CPM, 1.0, 2000, 0.25, 4.0, 4.0),
    ("kt4000_ks000", MIN_CPM, 1.0, 4000, 0.00, 4.0, 4.0),
    ("flat100_phi4", MIN_CPM, None, 0, 1.00, 4.0, 4.0),
]


def spread_phi(e, bctrl, min_cpm, tail, ktop, kspread, phi_keep, phi_tail):
    """Exactly `emission_spread_vector` + `emit_counts`'s per-gene phi (164's function)."""
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
    print(f"gate: pseudobulk CPM > {MIN_CPM}\n")
    print(f"{'policy':<30}{'gated':>7}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'cap?':>6}"
          f"{'rho':>6}{'s*':>8}{'mse@s*':>9}{'score':>8}{'d overall':>11}")
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
        s = np.sqrt(claim / B)
        mse_opt = 1 + (s * s * (A - B) + Kc - 2 * s * X) / sD
        sc = min(max((0.987862879 - mse_opt) / 0.950162566, 0), 1)
        out[name] = dict(claim=claim, mse_opt=mse_opt, sc=sc, s=s, ch=ch)
        base = out.get("a15 AS SHIPPED")
        dov = (sc - base["sc"]) / 6 if base else 0.0
        print(f"{name:<30}{int(np.mean(ng)):>7}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
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
        print(f"  {nm:<30} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(nm)
    a15 = out["a15 AS SHIPPED"]
    ok15 = abs(a15["claim"] - 0.24581) / 0.24581 < 0.02
    print(f"  {'a15 AS SHIPPED (164 section 5)':<30} claim {a15['claim']:.5f}  "
          f"expect 0.24581  ({100*(a15['claim']-0.24581)/0.24581:+.2f}%)  "
          f"{'PASS' if ok15 else 'FAIL'}")
    if not ok15:
        bad.append("a15")
    if bad:
        sys.exit(f"ABORT: claim known-answer gate failed for {bad}")

    print("\n--- held at a15's SHIPPED VCC_ALPHA 0.982 (a15 is ON the per-context cliff, so")
    print("    an arm that loses claim gets no amplitude rescue: d raw mse = -d claim / sumD)")
    print(f"{'arm':<30}{'claim':>9}{'d claim':>10}{'d raw mse':>11}"
          f"{'raw mse':>10}{'score_mse':>11}{'d overall':>11}")
    s15 = 0.982 / 1.0796
    mse15 = 1 + (s15 * s15 * (A - B) + Kc - 2 * s15 * X) / sD
    for nm in [p[0] for p in POLICIES]:
        dc = out[nm]["claim"] - a15["claim"]
        raw = mse15 - dc / sD
        sc = min(max((0.987862879 - raw) / 0.950162566, 0), 1)
        sc15 = min(max((0.987862879 - mse15) / 0.950162566, 0), 1)
        print(f"{nm:<30}{out[nm]['claim']:>9.5f}{dc:>+10.5f}{-dc/sD:>+11.5f}"
              f"{raw:>10.5f}{sc:>11.5f}{(sc-sc15)/6:>+11.6f}")


if __name__ == "__main__":
    main()

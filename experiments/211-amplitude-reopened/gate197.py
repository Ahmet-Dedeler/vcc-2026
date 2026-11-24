"""197 Phase 0 — the known-answer gates that license the inversion.

Every number here is already on the record.  If any one of them fails, the instrument is not
the instrument 096/164/169/171/187 used and the whole inversion is void (rule 10's corollary).

Gates, and where each comes from:

  G1  claim(spread 0.25, phi 0)   = 0.11860      096 `claim_grid2.log`   (tol 3%)
  G2  claim(spread 1.00, phi 0)   = 0.19216      096 `claim_grid.log`    (tol 3%)
  G3  claim(spread 0.25, phi 4)   = 0.16478      096 `board_solve` / a8  (tol 3%)
  G4  claim(a15 as shipped)       = 0.24581      164 §5                  (tol 3%)
  G5  a15 raw mse / score_mse     = 0.88228 / 0.11112     164 §5 / 169 / 171
  G6  the C_p-free LS solve on a8+a10+a11  ->  X 0.019240, sumD 0.229141   187 [KA-i]
  G7  OUT-OF-SAMPLE, both predicted before the arm was built:
        a10 raw mse predicted 1.3567   (board 1.3524)
        a11 raw mse predicted 0.9021   (board 0.9100)

The 3% tolerance on G1-G4 is 096's own documented spread between its two `claim` grids on one
policy (`claim_grid` 0.116800 vs `claim_grid2` 0.11860 = 1.5%).

    uv run python gate197.py

Writes work/gate197.json.  No cell-eval2, no build, no submission, nothing written outside
experiments/197-cosine-inversion/.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
W096 = f"{REPO}/experiments/096-mse-floor/work"
sys.path.insert(0, f"{REPO}/src")
sys.path.insert(0, f"{REPO}/experiments/096-mse-floor")

from cell_eval2.moments import jackknife_correction          # noqa: E402
import mech_alpha_chal as M                                  # noqa: E402

CTX = ("A", "B", "C")
NSUB = 16
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI = 10.0, 8000, 0.25, 1.0, 4.0
ALPHA_PROD = 1.0796
S15 = 0.982 / ALPHA_PROD
B_ANCHOR, SPAN = 0.987862879, 0.950162566

BOARD_RAW = {"a8": 1.362666960949, "a10": 1.3524148, "a11": 0.9099893, "a12": 0.9013492,
             "a14": 0.882775, "a15": 0.873805, "a23": 0.872353, "knee": 0.874491}
ALPHA = {"a8": 1.0796, "a10": 1.0796, "a11": 0.838, "a12": 1.0108, "a14": 0.982}


def score_mse(raw):
    return min(max((B_ANCHOR - raw) / SPAN, 0.0), 1.0)


def spread_phi(e, bctrl, tail, ktop, kspread, phi_keep, phi_tail, min_cpm=MIN_CPM):
    """`emission_spread_vector` verbatim (164's transcription; the whole vector starts at
    `tail_spread`, so sub-gate genes sit in the tail)."""
    n = e.size
    if tail is None:
        return np.full(n, float(kspread)), np.full(n, float(phi_keep)), None
    sv = np.full(n, float(tail))
    pv = np.full(n, float(phi_tail))
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > min_cpm) & (bctrl > 0)
    score = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    order = np.lexsort((np.arange(n), -score))[:k] if k > 0 else np.array([], int)
    sv[order] = float(kspread)
    pv[order] = float(phi_keep)
    return sv, pv, order


def claim_of(D, cpn, tail, ktop, kspread, phi_k, phi_t):
    """`claim = sum_p w_p min(C_hat_p, C_p)`, context-mean, on the real challenge panel."""
    rows = []
    for c in CTX:
        d = D[c]
        idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
        vals = []
        for j in range(idx.size):
            e = d["dh"][idx[j]]
            prof = np.expm1(d["bctrl"] + e)
            np.clip(prof, 0, None, out=prof)
            prof /= prof.sum()
            sv, pv, _ = spread_phi(e, d["bctrl"], tail, ktop, kspread, phi_k, phi_t)
            libs = np.full(M.NCELL, d["S"][idx[j]] / M.NCELL)
            v, col, ip = M.emit_counts(prof, libs, 0.0, np.random.default_rng(3000 + j),
                                       mode="quasi", spread=sv, overdispersion=pv)
            Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                            shape=(M.NCELL, d["G"]))
            vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1, M.TS,
                                             chunk=64)[0])
        cp = float(cpn[c].mean())
        rows.append((float(np.mean(vals)), cp,
                     float(np.mean(np.minimum(vals, cp)))))
    return (float(np.mean([300 * r[2] / 18532 for r in rows])),
            float(np.mean([r[0] for r in rows])))


def main():
    os.makedirs(WORK, exist_ok=True)
    out, lines = {}, []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))

    bs = np.load(f"{W096}/board_solve.npz")
    sD, X0 = float(bs["sumD"]), float(bs["X"])
    A_v11, B_v11 = float(bs["v11_A"]), float(bs["v11_B"])
    A_a8, B_a8 = float(bs["a8_A"]), float(bs["a8_B"])
    Kc, K = float(bs["v11_Kc"]), float(bs["K"])
    cl_a8, cl_a10 = float(bs["a8_claim"]), 0.169139
    P(f"096 board_solve : X {X0:.6f}  sumD {sD:.6f}  Kc {Kc:.6f}  cap K {K:.5f}")
    P(f"                  v11 A {A_v11:.6f} B {B_v11:.6f} | a8/a10 A {A_a8:.6f} B {B_a8:.6f}")

    # ---------------------------------------------------------------- G1-G4, the emitter gates
    D = M.load(f"{W096}/v11-stats.npz")
    cpn = dict(np.load(f"{W096}/cp-proxy-ntc.npz").items())
    P("\n--- G1-G4  challenge-panel `claim` (16 targets x 3 contexts, production emitter) ---")
    GATES = [("G1 flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860),
             ("G2 flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216),
             ("G3 flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478),
             ("G4 a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581)]
    bad, claims = [], {}
    for name, tail, kt, ks, pk, pt, want in GATES:
        got, chat = claim_of(D, cpn, tail, kt, ks, pk, pt)
        rel = (got - want) / want
        ok = abs(rel) < 0.03
        claims[name.split()[0]] = got
        P(f"  {name:<26} claim {got:.5f}  expect {want:.5f}  ({100*rel:+.2f}%)  "
          f"C_hat {chat:.2f}   {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    out["claims"] = claims
    if bad:
        (open(f"{WORK}/gate197.log", "w")).write("\n".join(lines))
        sys.exit(f"ABORT (rule 10): claim gate failed for {bad}")

    # ---------------------------------------------------------------- G5, a15 raw / score
    a15_claim = claims["G4"]
    a15_raw = 1 + (S15**2 * A_v11 - 2 * S15 * X0 + Kc
                   - min(S15**2 * B_v11, a15_claim)) / sD
    P(f"\n--- G5  a15 through the identity at its shipped VCC_ALPHA 0.982 ---")
    P(f"  s = {S15:.5f}   s^2 B = {S15**2*B_v11:.5f}   claim = {a15_claim:.5f}   "
      f"branch = {'s^2 B' if S15**2*B_v11 <= a15_claim else 'claim'}")
    P(f"  raw mse {a15_raw:.5f} (expect 0.88228)   score_mse {score_mse(a15_raw):.5f} "
      f"(expect 0.11112)   {'PASS' if abs(a15_raw-0.88228) < 5e-4 else 'FAIL'}")
    out["G5"] = dict(raw=a15_raw, score=score_mse(a15_raw))
    assert abs(a15_raw - 0.88228) < 5e-4, "G5 FAILED"

    # ---------------------------------------------------------------- G6, the C_p-free solve
    from mech_a12_audit import cp_proxies, measure
    cp16 = {c: cp_proxies()["ntc"][c][0] for c in CTX}
    arms = {}
    for arm in ("a8", "a10", "a11", "a12"):
        rows = measure(f"{W096}/{arm}-stats.npz", cp16)
        arms[arm] = {r["ctx"]: r for r in rows}
    P("\n--- G6  the C_p-free least-squares solve on a8 + a10 + a11 ---")
    Mmat, bvec = [], []
    for a in ("a8", "a10", "a11"):
        s = ALPHA[a] / ALPHA_PROD
        num = np.mean([arms[a][c]["A"] + arms[a][c]["Kc"]
                       - min(arms[a][c]["claim"], arms[a][c]["B"]) for c in CTX])
        Mmat.append([-2.0 * s, -(BOARD_RAW[a] - 1.0)])
        bvec.append(-num)
    sol, *_ = np.linalg.lstsq(np.asarray(Mmat), np.asarray(bvec), rcond=None)
    Xp, sDp = float(sol[0]), float(sol[1])
    okG6 = abs(Xp - 0.019240) < 2e-5 and abs(sDp - 0.229141) < 2e-4
    P(f"  X = {Xp:.6f} (187 KA-i: 0.019240)   sumD = {sDp:.6f} (0.229141)   "
      f"{'PASS' if okG6 else 'FAIL'}")
    for a in ("a8", "a10", "a11"):
        s = ALPHA[a] / ALPHA_PROD
        r = 1 + np.mean([arms[a][c]["A"] + arms[a][c]["Kc"]
                         - min(arms[a][c]["claim"], arms[a][c]["B"]) for c in CTX]
                        ) / sDp - 2 * s * Xp / sDp
        P(f"    {a}: model {r:.6f}  board {BOARD_RAW[a]:.6f}  residual {r-BOARD_RAW[a]:+.6f}")
    out["G6"] = dict(X=Xp, sumD=sDp)
    assert okG6, "G6 FAILED"

    # ---------------------------------------------------------------- G7, the out-of-sample record
    P("\n--- G7  the two OUT-OF-SAMPLE predictions 096 made before either arm existed ---")
    P("  (096's own v11/v12 two-point solve: X 0.018846, sumD 0.228572; a8's pruned-whitener A)")
    a10_pred = 1 + (A_a8 - 2 * X0 + Kc - min(B_a8, cl_a10)) / sD
    s11 = ALPHA["a11"] / ALPHA_PROD
    a11_pred = 1 + (s11**2 * A_a8 - 2 * s11 * X0 + Kc
                    - min(s11**2 * B_a8, cl_a10)) / sD
    for tag, pred, want, board in (("a10", a10_pred, 1.3567, BOARD_RAW["a10"]),
                                   ("a11", a11_pred, 0.9021, BOARD_RAW["a11"])):
        P(f"  {tag}: reproduced prediction {pred:.4f} (096 published {want:.4f})  "
          f"board {board:.4f}  error {pred-board:+.4f}  "
          f"{'PASS' if abs(pred-want) < 5e-4 else 'FAIL'}")
        out[f"G7_{tag}"] = dict(pred=pred, published=want, board=board)
    assert abs(a10_pred - 1.3567) < 5e-4 and abs(a11_pred - 0.9021) < 5e-4, "G7 FAILED"

    P("\nALL PHASE-0 GATES PASS.")
    json.dump(out, open(f"{WORK}/gate197.json", "w"), indent=1, default=float)
    open(f"{WORK}/gate197.log", "w").write("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

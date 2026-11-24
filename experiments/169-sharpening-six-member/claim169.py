"""169 — price `mse` for the sharpening arm on the REAL challenge panel, via 096's `claim` route.

WHY THIS AND NOT A HARNESS
--------------------------
162 §3 measured that a cross-line validator cannot price `mse` at all (LOO RMS 0.112 against
a predict-zero null of 0.044 -- 2.5x WORSE than predicting nothing), because its truth is a
different cell line and the amplitude it wants is 17-46x smaller than production ships.  The
only instrument that has ever retrodicted a board `mse` move is 096's decomposition,

    mse - 1  =  ( s^2*A  -  2*s*X  +  Kc  -  min(s^2*B, claim) ) / sumD                (*)

evaluated on `data/challenge/context_{A,B,C}.h5ad`, where `A`, `B`, `Kc` and `claim` are all
exactly computable from our own side and only `X` and `sumD` come from the board solve.

164 established the state this arm has to be priced against: **a15's `claim` is 0.24581
against the reference cap 0.24653 (99.7%), and `s^2*B = 0.23268`, so `min(s^2*B, claim)` is
on the `s^2*B` branch with only 5.3% of headroom.**  An arm that drops `claim` below 0.23268
falls off a cliff whose slope is `-1/sumD = -4.375` raw `mse` per unit `claim`, with no
amplitude compensation available.  Sharpening does not change the emission policy, but it
DOES change `A`, `X`, `B` and the profile the emitter samples, so `claim` has to be
re-derived rather than reused.

WHAT MOVES, AND WHAT IS TRANSFERRED
-----------------------------------
`A`, `B`, `claim`  : computed here, on the challenge panel, no reference.
`X`                : the one term that needs the withheld truth.  Priced at THREE settings --
                     `rX = 1` (the null), and the two per-surface ratios `xratio169.py`
                     measures -- and the GO/NO-GO uses the null, because 156/162 warn
                     explicitly against pricing a profile arm by scaling a harness cosine.

KNOWN-ANSWER GATES (rule 10's corollary), all from 164 §1 / 096's own logs:
    flat 0.25/phi0   claim 0.11860        flat 1.00/phi0   claim 0.19216-0.19438
    flat 0.25/phi4   claim 0.16478-0.16676 (a8 as shipped)
    a15 as shipped   claim 0.24581, raw mse 0.88228, score_mse 0.11112
Tolerance 3%, the documented spread between 096's own two grids on one policy.

No cell-eval2, no build, no submission.
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
sys.path.insert(0, f"{REPO}/src")
sys.path.insert(0, f"{REPO}/experiments/096-mse-floor")
sys.path.insert(0, HERE)
from cell_eval2.moments import jackknife_correction          # noqa: E402
import mech_alpha_chal as M                                  # noqa: E402
from sharp import sharpen, sharpen_lfc                       # noqa: E402

# a15 ships `VCC_MASS_FIX=minnorm`, so the delivered effect IS the mass-fixed one and the
# renderer's final division is the identity.  `MASSFIX=0` reproduces 164's exact convention
# (which used the raw v11 `dh`) and is kept so the two can be compared: without the fix the
# renderer's own normaliser drifts with `q` and moves `A` and `B` TOGETHER, which is why the
# common mode `A - B` is the quantity that actually enters (*).
MASSFIX = os.environ.get("MASSFIX", "1") != "0"
SPACE = os.environ.get("SPACE", "h")
_PROD = None


def massfix(h, bctrl, gate):
    global _PROD
    if not MASSFIX:
        return h
    if _PROD is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "prod_builder_169c", f"{REPO}/scripts/make_submission_pm_cd4.py")
        _PROD = importlib.util.module_from_spec(spec)
        sys.modules["prod_builder_169c"] = _PROD
        spec.loader.exec_module(_PROD)
    out = np.empty_like(h)
    for i in range(h.shape[0]):
        out[i] = _PROD.mass_consistent_effect(h[i], bctrl, gate, mode="minnorm")
    return out

STATS = f"{REPO}/experiments/096-mse-floor/work/v11-stats.npz"
NSUB = 16
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI = 10.0, 8000, 0.25, 1.0, 4.0
S15 = 0.982 / 1.0796                    # a15's amplitude relative to v11's, 096 (b)
SCORE = lambda raw: min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)  # noqa: E731
QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)


def spread_phi(e, bctrl, tail, ktop, kspread, phi_keep, phi_tail, min_cpm=MIN_CPM):
    """`emission_spread_vector` verbatim (164's transcription, including the trap that the
    whole vector is initialised to `tail_spread` so sub-gate genes sit in the tail)."""
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


def claim_of(D, effects, tail, ktop, kspread, phi_k, phi_t, cpn, keeptop_ref=None):
    """`claim = Sum_p w_p min(C_hat_p, C_p)` for a given per-target effect matrix."""
    rows, ov = [], []
    for c in M.CTX:
        d = D[c]
        E = effects[c]
        idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
        vals = []
        for j in range(idx.size):
            e = E[idx[j]]
            prof = np.expm1(d["bctrl"] + e)
            np.clip(prof, 0, None, out=prof)
            prof /= prof.sum()
            sv, pv, order = spread_phi(e, d["bctrl"], tail, ktop, kspread, phi_k, phi_t)
            if keeptop_ref is not None and order is not None:
                ref = keeptop_ref[c][j]
                ov.append(len(set(order.tolist()) & set(ref.tolist())) / max(len(ref), 1))
            libs = np.full(M.NCELL, d["S"][idx[j]] / M.NCELL)
            v, col, ip = M.emit_counts(prof, libs, 0.0, np.random.default_rng(3000 + j),
                                       mode="quasi", spread=sv, overdispersion=pv)
            Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                            shape=(M.NCELL, d["G"]))
            vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1, M.TS,
                                             chunk=64)[0])
        cp = float(cpn[c].mean())
        rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals, cp)))))
    ch = float(np.mean([r[0] for r in rows]))
    cp = float(np.mean([r[1] for r in rows]))
    claim = float(np.mean([300 * r[2] / 18532 for r in rows]))
    return claim, ch, cp, (float(np.mean(ov)) if ov else None)


def keeptop_of(D, effects):
    out = {}
    for c in M.CTX:
        d, E = D[c], effects[c]
        idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
        out[c] = [spread_phi(E[i], d["bctrl"], TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI)[2]
                  for i in idx]
    return out


def AB(D, effects, s):
    """`A` and `B` (the across-perturbation budget) at amplitude `s`, per 096's `terms`."""
    As, Bs = [], []
    for c in M.CTX:
        d = D[c]
        sums = M.render(d["bctrl"], effects[c], d["S"], s)
        a, b, _ = M.terms(d, sums)
        As.append(a)
        Bs.append(b)
    return float(np.mean(As)), float(np.mean(Bs))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    D = M.load(STATS)
    cpn = dict(np.load(f"{REPO}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(f"{REPO}/experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
    print(f"board solve: sumD {sD:.6f}  X {X0:.6f}  Kc {Kc:.6f}  cap K {K:.5f}")
    print(f"v11 A {A0:.6f}  B {B0:.6f} | a15 amplitude s = {S15:.5f}, "
          f"s^2 B = {S15**2 * B0:.5f}\n")

    gate = {}
    for c in M.CTX:
        b = D[c]["bctrl"]
        expr = np.expm1(b)
        gate[c] = (1e6 * expr / expr.sum() > MIN_CPM) & (b > 0)
    print("production gate (pb CPM > 10) per context: "
          + ", ".join(f"{c}={int(gate[c].sum())}" for c in M.CTX))

    # ---------------------------------------------------------------- known answers
    raw = {c: D[c]["dh"] for c in M.CTX}
    print(f"\n--- known-answer gates (164 §1) ---")
    GATES = [("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860, 0.03),
             ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216, 0.03),
             ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478, 0.03),
             ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581, 0.03)]
    bad = []
    for name, tail, kt, ks, pk, pt, want, tol in GATES:
        got, ch, cp, _ = claim_of(D, raw, tail, kt, ks, pk, pt, cpn)
        ok = abs(got - want) / want < tol
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f}  "
              f"({100 * (got - want) / want:+.2f}%)  C_hat {ch:.2f}  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    if bad:
        sys.exit(f"ABORT: claim known-answer gate failed for {bad}")

    A_chk, B_chk = AB(D, raw, 1.0)
    print(f"  A/B round trip at s=1: A {A_chk:.6f} (board_solve {A0:.6f}, "
          f"{100*(A_chk/A0-1):+.3f}%)  B {B_chk:.6f} ({100*(B_chk/B0-1):+.3f}%)")
    assert abs(A_chk / A0 - 1) < 5e-3 and abs(B_chk / B0 - 1) < 5e-3, "A/B round trip failed"

    a15_claim, _, _, _ = claim_of(D, raw, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
    a15_raw = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc
                   - min(S15**2 * B0, a15_claim)) / sD
    print(f"  a15 raw mse {a15_raw:.5f} (164: 0.88228)  score_mse {SCORE(a15_raw):.5f} "
          f"(164: 0.11112)")

    # ---------------------------------------------------------------- the arms
    rx = json.loads(open(os.path.join(WORK, "xratio169.json")).read())
    ktref = keeptop_of(D, raw)
    print(f"\n--- sharpening on the challenge panel, `row` renormalisation, "
          f"mass fix {'ON' if MASSFIX else 'OFF'}, space {SPACE} ---")
    print(f"{'q':>5}{'A':>10}{'rA':>9}{'B':>10}{'claim':>9}{'d claim':>10}{'C_hat':>8}"
          f"{'keeptop ov':>11}{'rho branch':>12}")
    rec = {}
    for q in QS:
        def _t(c):
            if SPACE == "lfc":
                return sharpen_lfc(raw[c], q, D[c]["bctrl"], gate[c], mode="row")
            return sharpen(raw[c], q, gate[c], mode="row")
        eff = {c: massfix(_t(c), D[c]["bctrl"], gate[c]) for c in M.CTX}
        cl, ch, cp, ov = claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn,
                                  keeptop_ref=ktref)
        Aq, Bq = AB(D, eff, 1.0)
        branch = "s^2 B" if S15**2 * Bq <= cl else "CLAIM (cliff)"
        rec[q] = dict(A=Aq, B=Bq, claim=cl, chat=ch, ov=ov)
        print(f"{q:>5.2f}{Aq:>10.5f}{Aq/rec[1.00]['A']:>9.5f}{Bq:>10.5f}{cl:>9.5f}"
              f"{cl - rec[1.00]['claim']:>+10.5f}{ch:>8.2f}"
              f"{(ov if ov is not None else float('nan')):>11.5f}{branch:>12}")

    print(f"\n--- raw `mse` and `score_mse` at a15's shipped VCC_ALPHA 0.982 ---")
    for tag, key in (("rX = 1 (NULL, used for the GO/NO-GO)", None),
                     ("rX = hct116 (measured)", "hct116|row"),
                     ("rX = hek293t (measured)", "hek293t|row")):
        print(f"\n  {tag}")
        print(f"  {'q':>5}{'X':>11}{'raw mse':>11}{'d raw':>10}{'score_mse':>11}"
              f"{'d score':>10}{'d overall':>11}")
        base = None
        for q in QS:
            r = rec[q]
            if key is None:
                Xq = X0
            else:
                row = [z for z in rx[key] if abs(z["q"] - q) < 1e-9][0]
                zero = [z for z in rx[key] if abs(z["q"] - 1.0) < 1e-9][0]
                Xq = X0 * row["X"] / zero["X"]
            mse = 1 + (S15**2 * r["A"] - 2 * S15 * Xq + Kc
                       - min(S15**2 * r["B"], r["claim"])) / sD
            sc = SCORE(mse)
            if base is None:
                base = (mse, sc)
            print(f"  {q:>5.2f}{Xq:>11.6f}{mse:>11.5f}{mse - base[0]:>+10.5f}"
                  f"{sc:>11.5f}{sc - base[1]:>+10.5f}{(sc - base[1]) / 6:>+11.6f}")

    with open(os.path.join(WORK, "claim169%s%s.json" % ("" if MASSFIX else "-raw", "" if SPACE == "h" else "-lfc")), "w") as fh:
        json.dump({str(k): v for k, v in rec.items()}, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'claim169.json')}")


if __name__ == "__main__":
    sys.exit(main())

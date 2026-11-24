"""188 — `mse` for this experiment's arms, via 096's `claim` route on the real challenge panel.

Same instrument, same four known-answer gates and the same A/B round trip as 169 §3, because
162 §3 measured that a cross-line validator cannot price `mse` at all (LOO RMS 0.112 against a
predict-zero null of 0.044) and 096's decomposition is the only thing that has retrodicted a
board `mse` move:

    mse - 1 = ( s^2 A - 2 s X + Kc - min(s^2 B, claim) ) / sumD

`X` is the one term needing the withheld truth; it is held at the null (`rX = 1`), which is what
rule 13 credits a shape-only arm.

NOTE ON WHICH MAPS ARE PRICEABLE HERE.  096's cached object is **v11's delivered** effect
(`alpha = 1.0796`, post-render), not a23's pre-render effect matrix the maps were fitted on.
A power transform with the row norm restored is EXACTLY invariant to a global rescale of the
row -- `(c|h|)^q` renormalised is `|h|^q` renormalised -- so `q` arms transfer to this object
unchanged.  A fitted map on a RAW `|p|` axis does not, which is precisely why `normfit188.py`
exists; only the scale-free (`nmap_*`) maps are priced here, and the raw-axis ones are recorded
as NOT PRICEABLE on this instrument rather than priced wrongly.
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
sys.path.insert(0, f"{REPO}/experiments/169-sharpening-six-member")
sys.path.insert(0, HERE)
from cell_eval2.moments import jackknife_correction    # noqa: E402
import mech_alpha_chal as M                            # noqa: E402
from sharp import sharpen, sharpen_lfc                 # noqa: E402
from mapfit import Map                                 # noqa: E402
from normfit188 import apply_norm                      # noqa: E402
from tilt188 import tilt                               # noqa: E402

STATS = f"{REPO}/experiments/096-mse-floor/work/v11-stats.npz"
NSUB = 16
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI = 10.0, 8000, 0.25, 1.0, 4.0
S15 = 0.982 / 1.0796
SCORE = lambda raw: min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)  # noqa: E731
_PROD = None


def massfix(h, b, gate):
    global _PROD
    if _PROD is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "prod_builder_188c", f"{REPO}/scripts/make_submission_pm_cd4.py")
        _PROD = importlib.util.module_from_spec(spec)
        sys.modules["prod_builder_188c"] = _PROD
        spec.loader.exec_module(_PROD)
    out = np.empty_like(h)
    for i in range(h.shape[0]):
        out[i] = _PROD.mass_consistent_effect(h[i], b, gate, mode="minnorm")
    return out


def spread_phi(e, bctrl, tail, ktop, kspread, phi_keep, phi_tail, min_cpm=MIN_CPM):
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


def claim_of(D, effects, tail, ktop, kspread, phi_k, phi_t, cpn):
    rows = []
    for c in M.CTX:
        d, E = D[c], effects[c]
        idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
        vals = []
        for j in range(idx.size):
            e = E[idx[j]]
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
        rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals, cp)))))
    return (float(np.mean([300 * r[2] / 18532 for r in rows])),
            float(np.mean([r[0] for r in rows])))


def AB(D, effects, s):
    As, Bs = [], []
    for c in M.CTX:
        d = D[c]
        a, b, _ = M.terms(d, M.render(d["bctrl"], effects[c], d["S"], s))
        As.append(a); Bs.append(b)
    return float(np.mean(As)), float(np.mean(Bs))


def main() -> None:
    D = M.load(STATS)
    cpn = dict(np.load(f"{REPO}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(f"{REPO}/experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0 = float(bs["v11_Kc"]), float(bs["X"])
    gate, cpm = {}, {}
    for c in M.CTX:
        b = D[c]["bctrl"]
        expr = np.expm1(b)
        cpm[c] = 1e6 * expr / expr.sum()
        gate[c] = (cpm[c] > MIN_CPM) & (b > 0)
    raw = {c: D[c]["dh"] for c in M.CTX}

    print("--- known-answer gates (164 §1 / 096) ---")
    bad = []
    for name, tail, kt, ks, pk, pt, want in (
            ("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860),
            ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216),
            ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478),
            ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581)):
        got, ch = claim_of(D, raw, tail, kt, ks, pk, pt, cpn)
        ok = abs(got - want) / want < 0.03
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(name)
    if bad:
        sys.exit(f"ABORT: {bad}")
    A_chk, B_chk = AB(D, raw, 1.0)
    print(f"  A/B round trip at s=1: A {A_chk:.6f} vs {A0:.6f} "
          f"({100*(A_chk/A0-1):+.3f}%)  B {B_chk:.6f} vs {B0:.6f} "
          f"({100*(B_chk/B0-1):+.3f}%)")
    a15_claim, _ = claim_of(D, raw, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
    a15_raw = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, a15_claim)) / sD
    print(f"  a15 raw mse {a15_raw:.5f} (164: 0.88228)  score_mse {SCORE(a15_raw):.5f} "
          f"(164: 0.11112)")

    NRM = {k: Map.from_dict(v) for k, v in
           json.loads(open(os.path.join(WORK, "maps_norm.json")).read()).items()}
    # The BASELINE must be the GATED identity, because every arm gates: `sharpen(.,1.0,gate)`
    # is `where(gate, raw, 0)`.  Pricing against the UNGATED v11 effect instead charges the arm
    # for zeroing the renderer's off-gate residue -- worth d raw mse -0.0025, 10x the arm.
    # 169 §3 used the gated q=1.00 row for exactly this reason; the ungated row is kept below
    # as the known-answer row, since 164's `claim` 0.24581 is defined on it.
    ARMS = {
        "identity_GATED":  lambda c: np.where(gate[c], raw[c], 0.0),
        "identity_ungated (KA)": lambda c: raw[c],
        "q1.15":           lambda c: sharpen(raw[c], 1.15, gate[c], mode="row"),
        "q1.15L":          lambda c: sharpen_lfc(raw[c], 1.15, D[c]["bctrl"], gate[c],
                                                 mode="row"),
        "tilt-0.10":       lambda c: tilt(raw[c], gate[c], cpm[c], -0.10, 1.0),
        "q1.15_tilt-0.10": lambda c: tilt(raw[c], gate[c], cpm[c], -0.10, 1.15),
        "nmap_median":     lambda c: apply_norm(raw[c], NRM["median|POOLED"], gate[c]),
        "nmap_budget":     lambda c: apply_norm(raw[c], NRM["budget|POOLED"], gate[c]),
    }
    print(f"\n--- the arms, mass fix ON, `rX = 1` null, a15's VCC_ALPHA 0.982 ---")
    print(f"{'arm':<18}{'A':>10}{'rA':>9}{'B':>10}{'claim':>9}{'d claim':>10}"
          f"{'raw mse':>10}{'d raw':>9}{'score':>9}{'d overall':>11}{'branch':>13}")
    rec, base = {}, None
    for name, fn in ARMS.items():
        eff = {c: massfix(fn(c), D[c]["bctrl"], gate[c]) for c in M.CTX}
        cl, ch = claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
        Aq, Bq = AB(D, eff, 1.0)
        mse = 1 + (S15**2 * Aq - 2 * S15 * X0 + Kc - min(S15**2 * Bq, cl)) / sD
        sc = SCORE(mse)
        if base is None:
            base = (Aq, cl, mse, sc)
        rec[name] = dict(A=Aq, rA=Aq / base[0], B=Bq, claim=cl, chat=ch, mse=mse, score=sc,
                         d_overall=(sc - base[3]) / 6)
        print(f"{name:<18}{Aq:>10.5f}{Aq/base[0]:>9.5f}{Bq:>10.5f}{cl:>9.5f}"
              f"{cl-base[1]:>+10.5f}{mse:>10.5f}{mse-base[2]:>+9.5f}{sc:>9.5f}"
              f"{(sc-base[3])/6:>+11.6f}"
              f"{('s^2 B' if S15**2*Bq <= cl else 'CLAIM CLIFF'):>13}")
    with open(os.path.join(WORK, "claim188.json"), "w") as fh:
        json.dump(rec, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/claim188.json")


if __name__ == "__main__":
    sys.exit(main())

"""178 — `mse` for the composed-rider arms, on the REAL challenge panel, through 096's route.

Rule 12 exists because a19 was selected on a screen that never computed `mse` and lost
**-0.01352** there: a whitener change reshapes the effect vector and `mse` reads the profile
directly. Every arm here changes the whitening exponent, the donor mixture or the sign pattern,
so none of them gets the rule-12 exemption.

THE ROUTE, and what is measured versus transferred
--------------------------------------------------
096 solved the board's own `mse` identity for the shipped submission:

    mse_raw = 1 + ( s^2 (A - B) + Kc - 2 s X ) / sumD ,     s* = sqrt(claim / B)

with `X = 0.018846`, `sumD = 0.228572`, `Kc = 0.005410`, `A = 0.283614`, `B = 0.281236`
(`experiments/096-mse-floor/work/board_solve.npz`). `A`, `B` and `claim` are properties of OUR
prediction and are re-measured here for every arm on contexts A/B/C. `X = sum_p w_p <dh_p, d_p>`
needs the challenge truth, which we do not have, so only its RATIO is transferred, from
h1-test's delivered aggregate cosine (`kappa = X / sqrt(B sumD)`, so
`X_arm = X_a15 * (kappa_arm/kappa_a15) * sqrt(B_arm/B_a15)`). That is 097 §9h's discipline
verbatim: a ratio inside one frame, never an absolute across frames.

KNOWN-ANSWER ROWS (rule 10's corollary)
---------------------------------------
  claim at three flat emission policies   0.11860 / 0.19216 / 0.16478   (096's own logs, via
                                          164/171; the mandate quotes 0.19438 / 0.16652 from a
                                          later run, and both are printed)
  claim for a15 as shipped                0.24581 (164 section 5)
  mse_raw at the a15 arm                  must land near the board's own raw `mse` 0.872

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
from cell_eval2.moments import jackknife_correction  # noqa: E402
import mech_alpha_chal as M  # noqa: E402

TS = 50_000.0
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI = 10.0, 8000, 0.25, 1.0, 4.0
NSUB = 16
ARMS = ("a15", "g045", "g030", "flip", "w211", "g045_flip", "g045_w211",
        "all3_g045", "all3_g030")
# name, min_cpm, tail_spread(None = flat at keep_spread), keep_top, keep_spread, phi_k, phi_t
FLAT_POLICIES = [
    ("KNOWN flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860),
    ("KNOWN flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216),
    ("KNOWN flat 0.25/phi4", None, 0, 0.25, 4.0, 4.0, 0.16478),
]


def load_builder():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "prod_builder_178_mse", ROOT / "scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def spread_phi(e, bctrl, tail, ktop, kspread, phi_keep, phi_tail):
    """`emission_spread_vector` + `emit_counts`'s per-gene phi (164/171's function)."""
    n = e.size
    if tail is None:
        return np.full(n, float(kspread)), np.full(n, float(phi_keep))
    sv = np.full(n, float(tail))
    pv = np.full(n, float(phi_tail))
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > MIN_CPM) & (bctrl > 0)
    score = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k > 0:
        order = np.lexsort((np.arange(n), -score))[:k]
        sv[order] = float(kspread)
        pv[order] = float(phi_keep)
    return sv, pv


def claim_for(D, dh_by_ctx, cpn, tail, ktop, kspread, phi_k, phi_t):
    rows = []
    for c in M.CTX:
        d = D[c]
        dh = dh_by_ctx[c]
        idx = np.linspace(0, dh.shape[0] - 1, NSUB).astype(int)
        vals = []
        for j in range(idx.size):
            e = dh[idx[j]]
            prof = np.expm1(d["bctrl"] + e)
            np.clip(prof, 0, None, out=prof)
            prof /= prof.sum()
            sv, pv = spread_phi(e, d["bctrl"], tail, ktop, kspread, phi_k, phi_t)
            libs = np.full(M.NCELL, d["S"][idx[j]] / M.NCELL)
            v, col, ip = M.emit_counts(prof, libs, 0.0, np.random.default_rng(3000 + j),
                                       mode="quasi", spread=sv, overdispersion=pv)
            Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64),
                             ip.astype(np.int64)), shape=(M.NCELL, d["G"]))
            vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1, TS,
                                             chunk=64)[0])
        cp = float(cpn[c].mean())
        rows.append(float(np.mean(np.minimum(vals, cp))))
    return float(np.mean([300 * r / 18532 for r in rows]))


def main() -> None:
    prod = load_builder()
    D = M.load(str(ROOT / "experiments/096-mse-floor/work/v11-stats.npz"))
    cpn = dict(np.load(ROOT / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(ROOT / "experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0, X0, Kc = (float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"]),
                          float(bs["X"]), float(bs["v11_Kc"]))
    print(f"board solve: X {X0:.6f}  sumD {sD:.5f}  A {A0:.5f}  B {B0:.5f}  Kc {Kc:.6f}")

    panel = [line.split(",")[0] for line in
             (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
             if line.strip()]
    official = [line.strip() for line in
                (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                if line.strip()]
    assert len(panel) == 300 and len(official) == 18533

    # ---- known-answer rows first: v11's own dh through the flat policies ----------------
    v11_dh = {c: D[c]["dh"] for c in M.CTX}
    print("\n=== known-answer claim rows (v11 profile, flat emission policies) ===")
    ok = True
    for name, tail, ktop, ks, pk, pt, want in FLAT_POLICIES:
        got = claim_for(D, v11_dh, cpn, tail, ktop, ks, pk, pt)
        good = abs(got - want) / want < 0.02
        ok &= good
        print(f"  {name:<24} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if good else 'FAIL'}")
    if not ok:
        print("  *** a known-answer row failed -- everything below is void ***")

    # ---- per-arm A, B, claim on the challenge frame -------------------------------------
    h1 = json.loads((WORK / "h1rebuild178.json").read_text())["arms"]
    gates, bctrls = {}, {}
    for c in M.CTX:
        b = D[c]["bctrl"]
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gates[c] = (cpm > MIN_CPM) & (b > 0)
        bctrls[c] = b
    pos = {c: {t: i for i, t in enumerate(D[c]["tgt"])} for c in M.CTX}

    print(f"\n{'arm':<12}{'A/A15':>8}{'B/B15':>8}{'claim':>9}{'kappa/k15':>11}"
          f"{'s*':>8}{'mse raw':>10}{'mse sc':>9}{'d ovr mse':>11}")
    res: dict[str, dict] = {}
    for arm in ARMS:
        eff = np.load(WORK / f"eff.{arm}.npy")
        assert eff.shape == (300, 18533)
        dh_by_ctx, Ac, Bc = {}, [], []
        for c in M.CTX:
            rows = np.asarray([pos[c][t] for t in panel if t in pos[c]], dtype=np.int64)
            src = np.asarray([i for i, t in enumerate(panel) if t in pos[c]], dtype=np.int64)
            dh = np.zeros((len(D[c]["tgt"]), 18533))
            g = gates[c]
            b = bctrls[c]
            for r, s in zip(rows, src):
                gated = np.where(g, eff[s], 0.0)
                gated = prod.mass_consistent_effect(gated, b, g, mode="minnorm")
                p = np.clip(np.expm1(b + gated), 0.0, None)
                p /= p.sum()
                dh[r] = np.log1p(TS * p) - b
            dh_by_ctx[c] = dh
            sums = M.render(b, dh, D[c]["S"], 1.0)
            a, bb, _ = M.terms(D[c], sums)
            Ac.append(a)
            Bc.append(bb)
        claim = claim_for(D, dh_by_ctx, cpn, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI)
        res[arm] = dict(A=float(np.mean(Ac)), B=float(np.mean(Bc)), claim=claim)
        del eff, dh_by_ctx

    base = res["a15"]
    kap0 = float(np.sqrt(max(1.0 - h1["a15"]["sse_c"], 0.0)))
    out = {}
    for arm in ARMS:
        r = res[arm]
        rA, rB = r["A"] / base["A"], r["B"] / base["B"]
        kap = float(np.sqrt(max(1.0 - h1[arm]["sse_c"], 0.0)))
        kr = kap / kap0
        A, B = A0 * rA, B0 * rB
        X = X0 * kr * np.sqrt(rB)
        s = float(np.sqrt(r["claim"] / B))
        mse = 1.0 + (s * s * (A - B) + Kc - 2 * s * X) / sD
        sc = min(max((0.987862879 - mse) / 0.950162566, 0.0), 1.0)
        out[arm] = dict(**r, ratio_A=rA, ratio_B=rB, kappa_ratio=kr, s=s,
                        mse_raw=mse, mse_scaled=sc)
        d = (sc - out["a15"]["mse_scaled"]) / 6 if arm != "a15" else 0.0
        out[arm]["d_overall"] = d
        print(f"{arm:<12}{rA:>8.4f}{rB:>8.4f}{r['claim']:>9.5f}{kr:>11.4f}"
              f"{s:>8.4f}{mse:>10.5f}{sc:>9.5f}{d:>+11.6f}")
    print(f"\n  a15 arm raw mse {out['a15']['mse_raw']:.5f} against the board's own 0.872 "
          f"({100*(out['a15']['mse_raw']/0.872-1):+.1f}%) -- the anchor's absolute level is "
          f"096's, so only the DELTAS are read")
    (WORK / "mse178.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {WORK/'mse178.json'}")


if __name__ == "__main__":
    main()

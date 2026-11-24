"""200 — `mse` for every ablated component, through 096's `claim` route on the real panel.

WHY THIS ROUTE AND NOTHING ELSE.  Rule 13: "`mse` is priceable ONLY through 096's `claim`
route and ONLY for emission/dispersion arms."  162 §3 measured a cross-line validator's `mse`
LOO RMS at 0.112 against a predict-zero null of 0.044 -- 2.5x worse than predicting zero -- so
PLOLO is disqualified on this member by measurement, not by taste.

    mse - 1 = ( s^2 A - 2 s X + Kc - min(s^2 B, claim) ) / sumD                          (*)

`A`, `B`, `claim` are ours; `X`, `Kc`, `sumD` come from 096's two-point board solve.

TWO ARM CLASSES, PRICED DIFFERENTLY, AND THE DIFFERENCE IS THE WHOLE POINT.

  EMISSION arms (the tiered spread policy, `--overdispersion`) leave the per-gene panel totals
  and therefore `A`, `B` and `X` bit-identical -- v14 and a8 confirmed `d pds = 0.000000` on
  the board -- and move `mse` through `claim` alone, which is exactly the case rule 13 says
  this route is valid for.  They are priced ABSOLUTELY, on 096's own `dh`, and the four
  known-answer gates below ARE those arms' baselines.

  GEOMETRY arms move `A`, `B`, `claim` and `X`.  `X` needs the withheld truth, so it is held
  at the null (rX = 1) exactly as 169 §3 did, and `A`/`B`/`claim` are transferred as RATIOS
  measured on the reconstruction of the panel's effect matrix against its own SHIP row --
  never as levels, because the reconstruction zeroes the 300 panel coordinates and production
  does not.  A geometry arm's `mse` here is therefore a PARTIAL price and is labelled as one.

KNOWN-ANSWER GATES (164 §1 / 096's own logs), each aborting:
    flat 0.25/phi0   0.11860      flat 1.00/phi0   0.19216
    flat 0.25/phi4   0.16478      a15 as shipped   0.24581
    A/B round trip at s = 1, and a15 raw mse 0.88228 / score_mse 0.11112.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/162-retrodicting-validator")
sys.path.insert(0, HERE)

from cell_eval2.moments import jackknife_correction          # noqa: E402
import mech_alpha_chal as M                                  # noqa: E402

STATS = f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz"
NSUB = 16
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI = 10.0, 8000, 0.25, 1.0, 4.0
S15 = 0.982 / 1.0796                    # a15's amplitude relative to v11's, 096 (b)
SCORE = lambda raw: min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)  # noqa: E731

# --- the emission-law ablations.  (label, tail_spread, keep_top, keep_spread, phi) --------
# a15 ships tail 1.0 / keep_top 8000 / keep_spread 0.25 / phi 4.0.  `tail=None` is a FLAT
# spread over every gene, which is what a8 and every build before a12 used.
EMISSION_ARMS = [
    ("SHIP  tier8000 phi4",       TAIL, KEEP_TOP, KEEP_SPREAD, PHI),
    ("-tier (flat 0.25) phi4",    None, 0,        0.25,        PHI),   # == a8's law
    ("-phi  tier8000 phi0",       TAIL, KEEP_TOP, KEEP_SPREAD, 0.0),
    ("-tier -phi (flat 0.25 phi0)", None, 0,      0.25,        0.0),   # == v11's law
    ("tier 4000 phi4",            TAIL, 4000,     KEEP_SPREAD, PHI),
    ("tier 14000 phi4",           TAIL, 14000,    KEEP_SPREAD, PHI),
    ("tier8000 tail 0.50 phi4",   0.50, KEEP_TOP, KEEP_SPREAD, PHI),
    ("flat 1.00 phi4",            None, 0,        1.00,        PHI),
    ("flat 1.00 phi0",            None, 0,        1.00,        0.0),
    ("tier8000 phi 2",            TAIL, KEEP_TOP, KEEP_SPREAD, 2.0),
    ("tier8000 phi 8",            TAIL, KEEP_TOP, KEEP_SPREAD, 8.0),
]


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


def claim_of(D, effects, tail, ktop, kspread, phi_k, phi_t, cpn, min_cpm=MIN_CPM):
    """`claim = Sum_p w_p min(C_hat_p, C_p)` for a given per-target effect matrix."""
    rows = []
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
            sv, pv, _ = spread_phi(e, d["bctrl"], tail, ktop, kspread, phi_k, phi_t,
                                   min_cpm=min_cpm)
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
    return claim, ch, cp


def AB(D, effects, s):
    As, Bs = [], []
    for c in M.CTX:
        d = D[c]
        sums = M.render(d["bctrl"], effects[c], d["S"], s)
        a, b, _ = M.terms(d, sums)
        As.append(a)
        Bs.append(b)
    return float(np.mean(As)), float(np.mean(Bs))


# ------------------------------------------------------------------ geometry reconstruction
_PROD = None


def _builder():
    """The PRODUCTION `mass_consistent_effect`, not a copy of it (156's lesson, a19's cost)."""
    global _PROD
    if _PROD is None:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "prod_builder_200", f"{ROOT}/scripts/make_submission_pm_cd4.py")
        _PROD = importlib.util.module_from_spec(spec)
        sys.modules["prod_builder_200"] = _PROD
        spec.loader.exec_module(_PROD)
    return _PROD


def mass_fix_row(h, b, gate, mode):
    return _builder().mass_consistent_effect(h, b, gate, mode=mode)


def gate_of(bctrl, min_cpm):
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    return (cpm > min_cpm) & (bctrl > 0)


def geometry_effects(cfg, p, c4, op):
    """One ablated challenge-panel effect matrix, context-independent by construction."""
    from folds200 import consensus
    from plolo import E062, solve_fraction
    pre = consensus(p, pm_blend=cfg["pm"], weights=cfg["weights"]) * cfg["alpha"]
    if cfg["gamma"] is None:
        whitened = np.asarray(pre, dtype=np.float64).copy()
    else:
        whitened = op.apply(pre, p["wcols"], cfg["gamma"], cfg["lam"])
    if cfg["cd4_f"] is not None:
        frac, _ = solve_fraction(whitened, c4, cfg["cd4_f"])
        final, _ = E062.add_post_whiten_residual(whitened, c4, frac, p["wcols"])
    elif cfg["cd4_frac"] is not None:
        final, _ = E062.add_post_whiten_residual(whitened, c4, float(cfg["cd4_frac"]),
                                                 p["wcols"])
    else:
        final = whitened
    final = np.array(final, dtype=np.float64, copy=True)
    final[:, p["panel_positions"]] = 0.0
    return final


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    D = M.load(STATS)
    cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
    sD, A0, B0 = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
    Kc, X0, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
    print(f"board solve: sumD {sD:.6f}  X {X0:.6f}  Kc {Kc:.6f}  cap K {K:.5f}")
    print(f"v11 A {A0:.6f}  B {B0:.6f} | a15 amplitude s = {S15:.5f}, "
          f"s^2 B = {S15**2 * B0:.5f}")

    raw = {c: D[c]["dh"] for c in M.CTX}

    # ------------------------------------------------------------------- known-answer gates
    print(f"\n--- known-answer gates (164 §1, 169 §3) ---")
    GATES = [("flat 0.25/phi0", None, 0, 0.25, 0.0, 0.0, 0.11860, 0.03),
             ("flat 1.00/phi0", None, 0, 1.00, 0.0, 0.0, 0.19216, 0.03),
             ("flat 0.25/phi4 (a8)", None, 0, 0.25, 4.0, 4.0, 0.16478, 0.03),
             ("a15 as shipped", TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, 0.24581, 0.03)]
    bad = []
    for name, tail, kt, ks, pk, pt, want, tol in GATES:
        got, ch, cp = claim_of(D, raw, tail, kt, ks, pk, pt, cpn)
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
    a15_claim, _, _ = claim_of(D, raw, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn)
    a15_raw = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, a15_claim)) / sD
    print(f"  a15 raw mse {a15_raw:.5f} (164: 0.88228)  score_mse {SCORE(a15_raw):.5f} "
          f"(164: 0.11112)")
    assert abs(a15_raw - 0.88228) < 3e-4 and abs(SCORE(a15_raw) - 0.11112) < 3e-4

    out: dict[str, object] = dict(anchor=dict(sumD=sD, A0=A0, B0=B0, Kc=Kc, X0=X0, K=K,
                                              S15=S15, a15_claim=a15_claim,
                                              a15_raw=a15_raw, a15_score=SCORE(a15_raw)))

    # ------------------------------------------------------------------- EMISSION ablations
    print(f"\n--- EMISSION-law ablations (profile bit-identical: A, B, X pinned) ---")
    print(f"  {'arm':<30}{'claim':>10}{'C_hat':>9}{'branch':>14}{'raw mse':>10}"
          f"{'score':>9}{'d score':>10}{'d overall':>11}")
    erec = {}
    base = None
    for label, tail, kt, ks, phi in EMISSION_ARMS:
        cl, ch, cp = claim_of(D, raw, tail, kt, ks, phi, phi, cpn)
        branch = "s^2 B" if S15**2 * B0 <= cl else "CLAIM (cliff)"
        mse = 1 + (S15**2 * A0 - 2 * S15 * X0 + Kc - min(S15**2 * B0, cl)) / sD
        sc = SCORE(mse)
        if base is None:
            base = (mse, sc)
        erec[label] = dict(claim=cl, chat=ch, cp=cp, raw=mse, score=sc,
                           d_score=sc - base[1], d_overall=(sc - base[1]) / 6, branch=branch)
        print(f"  {label:<30}{cl:>10.5f}{ch:>9.2f}{branch:>14}{mse:>10.5f}{sc:>9.5f}"
              f"{sc - base[1]:>+10.5f}{(sc - base[1]) / 6:>+11.6f}")
    out["emission"] = erec
    with open(os.path.join(WORK, "mse200.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ------------------------------------------------------------------- GEOMETRY ablations
    if os.environ.get("SKIP_GEOMETRY"):
        print("\n(geometry arms skipped)")
        return
    import chalparts200
    from pds200 import ARMS
    from plolo import PROD_OPS, Operator
    p = chalparts200.parts()
    c4 = chalparts200.cd4()
    tgt_index = {t: i for i, t in enumerate(p["targets"])}
    order = {c: np.asarray([tgt_index[t] for t in D[c]["tgt"]], dtype=np.int64)
             for c in M.CTX}
    for c in M.CTX:
        assert len(order[c]) == len(D[c]["tgt"])

    skip = {"SHIP-dup", "noOperator", "a11", "a12", "BARE+whiten_unpruned"}
    names = [a for a in ARMS if a not in skip]
    if os.environ.get("GEOM_ONLY"):
        keep = set(os.environ["GEOM_ONLY"].split(",")) | {"SHIP"}
        names = [a for a in names if a in keep]
    by_op: dict[str, list[str]] = {}
    for a in names:
        by_op.setdefault(ARMS[a]["op"], []).append(a)
    grec = {}
    print(f"\n--- GEOMETRY ablations (rX = 1 null; A/B/claim as RATIOS vs SHIP) ---")
    print(f"  {'arm':<18}{'A':>11}{'rA':>9}{'B':>11}{'rB':>9}{'claim':>10}{'rClaim':>9}"
          f"{'raw mse':>10}{'score':>9}{'d overall':>11}")
    for opkey, arms in by_op.items():
        op = Operator(PROD_OPS / f"whitener_{opkey}.npz")
        for a in arms:
            cfg = ARMS[a]
            E = geometry_effects(cfg, p, c4, op)
            # Gate and (optionally) mass-fix per context, exactly as `profile_from_effect`
            # does at render time -- both are ablated flags, so neither can be assumed.
            eff = {}
            for c in M.CTX:
                b = D[c]["bctrl"]
                g = gate_of(b, cfg["min_cpm"])
                h = np.where(g[None, :], E[order[c]], 0.0)
                if cfg["mass"] is not None:
                    h = np.vstack([mass_fix_row(h[i], b, g, cfg["mass"])
                                   for i in range(h.shape[0])])
                eff[c] = h
            Aq, Bq = AB(D, eff, 1.0)
            cl, ch, cp = claim_of(D, eff, TAIL, KEEP_TOP, KEEP_SPREAD, PHI, PHI, cpn,
                                  min_cpm=cfg["min_cpm"])
            grec[a] = dict(A=Aq, B=Bq, claim=cl, chat=ch)
            print(f"    built {a:<18} A {Aq:.5f} B {Bq:.5f} claim {cl:.5f}", flush=True)
            with open(os.path.join(WORK, "mse200-geom.json"), "w") as fh:
                json.dump(grec, fh, indent=1, default=float)
        del op
    ship = grec["SHIP"]
    for a in names:
        r = grec[a]
        rA, rB, rC = r["A"] / ship["A"], r["B"] / ship["B"], r["claim"] / ship["claim"]
        mse = 1 + (S15**2 * A0 * rA - 2 * S15 * X0 + Kc
                   - min(S15**2 * B0 * rB, a15_claim * rC)) / sD
        sc = SCORE(mse)
        r.update(rA=rA, rB=rB, rClaim=rC, raw=mse, score=sc,
                 d_score=sc - SCORE(a15_raw), d_overall=(sc - SCORE(a15_raw)) / 6)
        print(f"  {a:<18}{r['A']:>11.5f}{rA:>9.5f}{r['B']:>11.5f}{rB:>9.5f}"
              f"{r['claim']:>10.5f}{rC:>9.5f}{mse:>10.5f}{sc:>9.5f}{r['d_overall']:>+11.6f}")
    out["geometry"] = grec
    with open(os.path.join(WORK, "mse200.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'mse200.json')}")


if __name__ == "__main__":
    sys.exit(main())

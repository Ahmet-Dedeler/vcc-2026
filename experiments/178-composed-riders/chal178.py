"""178 — `mse` for every arm on the CHALLENGE panel, through 170's per-context identity.

THIS SUPERSEDES `mse178.py` FOR THE VERDICT, AND WHY
----------------------------------------------------
`mse178.py` priced every arm at its own rho cliff (`s* = sqrt(claim/B)`), i.e. it silently
re-tuned `VCC_ALPHA` for each arm. 170 measured, production-staged on the real panel, that this
is the wrong convention for a rider that moves `B`: a15 sits at `rho = 0.938`, and the members
are read at the SHIPPED amplitude, so the identity that governs the board is

    MSE_c - 1 = ( A_c - 2 X_c + Kc_c - min(claim_c, B_c) ) / sumD

applied **per context**, with `min` evaluated at each arm (170's own correction: linearising it
around the shipped point turns -0.0039 into +0.0021 and is the whole result). `A`, `B` and
`claim` are properties of our own submission and are measured here; `X` is solved from a23's
own board point at the a15 configuration and carried across arms by the ratio measured on
h1-test.

`X` on h1-test is recovered from `h1rebuild178.json` without a second pipeline: `endpoints`
reports `sse = 1 - x^2/(a d)` and `gstar = x/a`, so `x = (1 - sse) d / gstar` and, `d` being
the same truth for every arm, `X_arm/X_a15 = [(1-sse)/gstar]_arm / [(1-sse)/gstar]_a15`.
**Cross-check against 170, which measured the same ratio through its own reconstruction:
gamma 0.45 / gamma 0.65 = 1.0889 there, 1.0901 here (+0.1%).**

Both amplitude conventions are reported:
  `mse_ship`  at the shipped alpha (170's convention, the one the board reads);
  `mse_knee`  after re-tuning alpha onto the arm's own cliff, `alpha_new = 0.982 sqrt(B15/B)`,
              which is legitimate only if the arm ships with that alpha AND `nmae`/`fid`/
              `reach` are re-read there too (they are not, so it is a sensitivity, not the
              verdict).

KNOWN-ANSWER ROWS
-----------------
  * `claim` at the a15 arm must land on 164's 0.24581 / 170's 0.24610.
  * `rho` at the a15 arm must land on 170's 0.938 / 164's 0.9466.
  * `A`, `B` at the a15 arm must land on 170's 0.23152 / 0.23031, and the `g045` arm on
    170's gamma-0.45 row (A 0.25669, B 0.25539) -- an independent staging of the same object.
  * with X solved at the a15 arm, its raw `mse` reproduces a23's board point BY CONSTRUCTION;
    the `g045` arm must then reproduce 170's measured raw `mse` 0.90965.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
sys.path.insert(0, str(ROOT / "experiments/170-gamma-mse"))

CTX = ("A", "B", "C")
NSUB = 24
TS = 50_000.0
MIN_CPM, KEEP_TOP, KEEP_SPREAD, TAIL, PHI, DISPERSION = 10.0, 8000, 0.25, 1.0, 4.0, 0.1
NCELL = 400
SUMD = 0.229141                              # 170 / MSE-MECHANISM's corrected solve
MSE_MAP = (0.987862879, 0.950162566)
A23_SCORE_MSE = 0.121581                     # a23's board point; every arm sits on a23
ARMS = ("a15", "g045", "g030", "flip", "w211", "flip_w211", "g045_flip", "g045_w211",
        "all3_g045", "all3_g030")


def load_builder():
    import importlib.util
    import os
    os.environ.setdefault("VCC_ALPHA", "0.982")
    os.environ.setdefault("VCC_CD4_TARGET_F", "0.1")
    os.environ["VCC_MASS_FIX"] = "minnorm"       # every arm sits on a23
    spec = importlib.util.spec_from_file_location(
        "prod_builder_178_chal", ROOT / "scripts/make_submission_pm_cd4.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    assert mod.MASS_FIX == "minnorm"
    return mod


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts
    from chal import context_library_medians          # 170's, cached

    B = load_builder()
    ctrl = np.load(ROOT / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cpn = dict(np.load(ROOT / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    genes = np.asarray([line.strip() for line in
                        (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                        if line.strip()])
    perts = np.asarray([line.split(",")[0] for line in
                        (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
                        if line.strip()])
    assert np.array_equal(genes, ctrl["genes"].astype(str))
    bctrl = {c: bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
             for c in CTX}
    Kc = {c: 300 * float(ctrl[f"jk_{c}"]) / 18532 for c in CTX}
    Cp = {c: float(np.mean(cpn[c])) for c in CTX}
    Lmed = context_library_medians()
    gpos = {g: i for i, g in enumerate(genes)}
    cols = {i: gpos[p] for i, p in enumerate(perts) if p in gpos}
    w = _row_weights(perts.size, genes.size, cols)
    sub = np.linspace(0, perts.size - 1, NSUB).astype(int)
    print(f"[chal178] Kc {[round(Kc[c], 6) for c in CTX]}  median control library "
          f"{[int(Lmed[c]) for c in CTX]}", flush=True)

    out: dict[str, dict] = {}
    for arm in ARMS:
        t0 = time.time()
        eff = np.load(WORK / f"eff.{arm}.npy")
        per = {}
        for c in CTX:
            b = bctrl[c]
            P = np.empty_like(eff)
            for i in range(eff.shape[0]):
                P[i] = B.profile_from_effect(eff[i], b, min_cpm=MIN_CPM)
            pm = np.log1p(TS * P)
            dh = pm - b[None, :]
            a = np.einsum("ij,ij->i", dh, dh)
            for k, col in cols.items():
                a[k] -= dh[k, col] ** 2
            A = float((a * w).sum())
            bud = float(_across_pert_budget(pm, cols, w))
            libs = np.full(NCELL, Lmed[c])
            vals = []
            for j in sub:
                spread = B.emission_spread_vector(
                    eff[j], b, min_cpm=MIN_CPM, tail_spread=TAIL, keep_top=KEEP_TOP,
                    keep_spread=KEEP_SPREAD, spread_map="block")
                v, cl, ip = emit_counts(P[j], libs, DISPERSION,
                                        np.random.default_rng(3000 + int(j)), mode="quasi",
                                        overdispersion=PHI, spread=spread)
                Xm = csr_matrix((v.astype(np.float64), cl.astype(np.int64),
                                 ip.astype(np.int64)), shape=(NCELL, genes.size))
                vals.append(jackknife_correction(Xm, np.zeros(NCELL, np.intp), 1, TS,
                                                 chunk=64)[0])
            claim = float(perts.size * np.mean(
                np.minimum(vals, PRED_TRACE_CAP_K * Cp[c])) / (genes.size - 1))
            per[c] = dict(A=A, B=bud, claim=claim, Kc=Kc[c],
                          rho=min(1.0, bud / claim) if claim else 1.0)
            del P, pm, dh
        out[arm] = dict(per=per,
                        A=float(np.mean([per[c]["A"] for c in CTX])),
                        B=float(np.mean([per[c]["B"] for c in CTX])),
                        claim=float(np.mean([per[c]["claim"] for c in CTX])))
        print(f"[chal178] {arm:<11} A {out[arm]['A']:.6f}  B {out[arm]['B']:.6f}  "
              f"claim {out[arm]['claim']:.6f}  rho "
              f"{[round(per[c]['rho'], 4) for c in CTX]}  [{time.time()-t0:.0f}s]", flush=True)
        del eff
        (WORK / "chal178.json").write_text(json.dumps(out, indent=1))

    # ---------------------------------------------------------------- X, from h1-test
    h1 = json.loads((WORK / "h1rebuild178.json").read_text())["arms"]
    xr = {a: ((1 - h1[a]["sse"]) / h1[a]["gstar"]) / ((1 - h1["a15"]["sse"]) / h1["a15"]["gstar"])
          for a in ARMS}
    target_raw = MSE_MAP[0] - A23_SCORE_MSE * MSE_MAP[1]
    base = out["a15"]
    num = float(np.mean([base["per"][c]["A"] + base["per"][c]["Kc"]
                         - min(base["per"][c]["claim"], base["per"][c]["B"]) for c in CTX]))
    X0 = (num - SUMD * (target_raw - 1.0)) / 2.0
    print(f"\n[chal178] a23 board anchor score_mse {A23_SCORE_MSE} -> raw {target_raw:.6f}; "
          f"solved X at the a15 arm = {X0:.6f} (170 solved 0.017766)")

    print("\n=== known-answer rows against 170's own staging ===")
    for label, got, want, tol in (
            ("a15 claim", base["claim"], 0.24610, 0.004),
            ("a15 A", base["A"], 0.23152, 0.01),
            ("a15 B", base["B"], 0.23031, 0.01),
            ("a15 rho (mean)", float(np.mean([base["per"][c]["rho"] for c in CTX])),
             0.9381, 0.02),
            ("g045 A", out["g045"]["A"], 0.25669, 0.01),
            ("g045 B", out["g045"]["B"], 0.25539, 0.01),
            ("g045 X ratio (h1)", xr["g045"], 1.0889, 0.01)):
        rel = abs(got - want) / abs(want)
        print(f"  {label:<20} got {got:.6f}  170 {want:.6f}  ({100*rel:+.2f}%)  "
              f"{'PASS' if rel < tol else 'FAIL'}")

    print(f"\n{'arm':<12}{'A':>10}{'B':>10}{'claim':>9}{'rho A/B/C':>22}{'X/X15':>8}"
          f"{'raw mse':>10}{'d ovr ship':>12}{'alpha*':>8}{'d ovr knee':>12}")
    res = {}
    for arm in ARMS:
        r = out[arm]
        X = X0 * xr[arm]
        raw = 1.0 + float(np.mean([
            r["per"][c]["A"] + r["per"][c]["Kc"]
            - min(r["per"][c]["claim"], r["per"][c]["B"]) - 2 * X for c in CTX])) / SUMD
        sc = min(max((MSE_MAP[0] - raw) / MSE_MAP[1], 0.0), 1.0)
        # knee: rescale the effect by s so that s^2 B = claim, per context, then re-evaluate
        knee_raw = 1.0 + float(np.mean([
            (lambda s: s * s * r["per"][c]["A"] + r["per"][c]["Kc"]
             - min(s * s * r["per"][c]["B"], r["per"][c]["claim"]) - 2 * s * X)(
                float(np.sqrt(r["per"][c]["claim"] / r["per"][c]["B"])))
            for c in CTX])) / SUMD
        knee_sc = min(max((MSE_MAP[0] - knee_raw) / MSE_MAP[1], 0.0), 1.0)
        alpha_star = 0.982 * float(np.mean(
            [np.sqrt(r["per"][c]["claim"] / r["per"][c]["B"]) for c in CTX]))
        res[arm] = dict(A=r["A"], B=r["B"], claim=r["claim"], X=X, x_ratio=xr[arm],
                        raw=raw, scaled=sc, knee_raw=knee_raw, knee_scaled=knee_sc,
                        alpha_star=alpha_star,
                        rho=[r["per"][c]["rho"] for c in CTX])
    b0 = res["a15"]
    for arm in ARMS:
        r = res[arm]
        print(f"{arm:<12}{r['A']:>10.5f}{r['B']:>10.5f}{r['claim']:>9.5f}"
              f"{'/'.join(f'{v:.3f}' for v in r['rho']):>22}{r['x_ratio']:>8.4f}"
              f"{r['raw']:>10.5f}{(r['scaled']-b0['scaled'])/6:>+12.6f}"
              f"{r['alpha_star']:>8.4f}{(r['knee_scaled']-b0['knee_scaled'])/6:>+12.6f}")
    print(f"\n  g045 raw mse {res['g045']['raw']:.5f} against 170's measured 0.90965 "
          f"({100*(res['g045']['raw']/0.90965-1):+.2f}%)")
    (WORK / "chal178_solve.json").write_text(json.dumps(res, indent=1))
    print(f"wrote {WORK/'chal178_solve.json'}")


if __name__ == "__main__":
    main()

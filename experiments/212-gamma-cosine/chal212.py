"""212 — the challenge panel's `A`, `B`, `claim`, rho and Ĉ/C_p at every gamma, in the
INCUMBENT's own geometry.

This is `experiments/170-gamma-mse/chal.py` with three differences, all recorded:

  1. `WORK` is this experiment's. Effects come from `work/eff/` (built by `gate212.py`,
     gated 20/20 bit-identical to 170's for the nine gammas 170 ran).
  2. The default geometry is **a15's, i.e. NO mass fix** — because that is what the
     incumbent `rc_callset` ships (`experiments/202-callset-rider/build202.sh`: "no
     mass_fix"). 170 ran only the `minnorm` (a23) geometry and noted the a15 sweep as not
     done. The board anchor is therefore a15's own `score_mse` 0.1200566151966586, which
     `rc_hybrid` and `rc_callset` both reproduced to sixteen significant figures.
  3. The grid is extended to 15 points and `Ĉ_p`, `C_p` and `Ĉ/C_p` are reported per
     context, because Phase 1's crux is WHICH cliff 208's "2.6x further away" refers to.

Two cliffs live in this table and they are not the same object:

  (a) the **rho = 1 cliff**: `s^2 B_c` overtaking `claim_c` in `min(claim_c, s^2 B_c)`.
      This is the one gamma moves, because gamma moves `A` and `B`.
  (b) the **claim-saturation cliff**: `Ĉ_p` falling back through `K * C_p` in
      `min(Ĉ_p, K C_p)`. This is the one 208 measured at 1.125 (a15) -> 2.940
      (`rc_callset`). It is an EMISSION property; gamma cannot move it.

KNOWN-ANSWER ROWS (rule 10):
  * `claim` at gamma 0.65 must land on 164's 0.24581 / 170's 0.24610 (+-0.5%);
  * per-context `Kc` must be 096 (e)'s 0.005084 / 0.005228 / 0.005919 exactly;
  * rho at gamma 0.65 must land on 170's 0.938 / 164's 0.9466;
  * the whole minnorm column must reproduce 170's `chal_minnorm.json` to 1e-12 on the nine
    shared gammas (`--geometry minnorm --ka170`).

    uv run python chal212.py                      # a15 geometry (the incumbent's)
    VCC_MASS_FIX=minnorm uv run python chal212.py  # a23 geometry, 170's known-answer column
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(REPO / "experiments/170-gamma-mse"))
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments/096-mse-floor"))

# a15's shipped emission policy (`submissions/a15-cd4f010.vcc.meta.json`), copied from
# `170/h1kit.py`'s module constants rather than imported: importing h1kit sets
# `VCC_MASS_FIX=minnorm` at module scope and its `builder()` RAISES if the mass fix is off,
# so it cannot express the incumbent's own geometry. The production builder itself refuses
# `VCC_MASS_FIX` together with `--emitter rc_hybrid|rc_callset`
# (`scripts/make_submission_pm_cd4.py:2755`), which is the decisive proof that the incumbent
# ships without it.
TS = 50_000.0
NCELL = 400
MIN_CPM = 10.0
TAIL_SPREAD, KEEP_TOP, KEEP_SPREAD = 1.0, 8000, 0.25
OVERDISPERSION, DISPERSION = 4.0, 0.1


def builder():
    """The PRODUCTION builder module, with whatever `VCC_MASS_FIX` the caller set."""
    import importlib.util
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_212c", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_212c"] = mod
    spec.loader.exec_module(mod)
    return mod

CTX = ("A", "B", "C")
NSUB = 24
GRID = [0.00, 0.30, 0.35, 0.40, 0.45, 0.50, 0.525, 0.55, 0.575, 0.60, 0.625,
        0.65, 0.70, 0.75, 0.80]

SUMD = 0.229141
MSE_MAP = (0.987862879, 0.950162566)
A15_SCORE_MSE, A23_SCORE_MSE = 0.1200566151966586, 0.121581
KC_KA = (0.005084, 0.005228, 0.005919)          # 096 (e), via 170 gate 7


def raw_from_score(s):
    return MSE_MAP[0] - s * MSE_MAP[1]


def context_library_medians():
    """Median control library per context. Cached by 170; recomputed here if absent."""
    import h5py
    cache = WORK / "chal_libs.json"
    if cache.exists():
        return {k: float(v) for k, v in json.loads(cache.read_text()).items()}
    src = REPO / "experiments/170-gamma-mse/work/chal_libs.json"
    if src.exists():                      # read 170's, write our own copy
        d = {k: float(v) for k, v in json.loads(src.read_text()).items()}
        cache.write_text(json.dumps(d, indent=1))
        return d
    out = {}
    for c in CTX:
        with h5py.File(REPO / f"data/challenge/context_{c}.h5ad", "r") as h:
            ip = h["X/indptr"][:]
            data = h["X/data"]
            tot = np.empty(ip.size - 1, dtype=np.float64)
            for lo in range(0, ip.size - 1, 2000):
                hi = min(lo + 2000, ip.size - 1)
                blk = np.asarray(data[int(ip[lo]):int(ip[hi])], dtype=np.float64)
                tot[lo:hi] = np.add.reduceat(blk, (ip[lo:hi] - ip[lo]).astype(np.int64))
        out[c] = float(np.median(tot))
    cache.write_text(json.dumps(out, indent=1))
    return out


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts

    B = builder()
    massfix = B.MASS_FIX
    tag = "minnorm" if massfix else "a15"
    print(f"[chal] MASS_FIX={tag}  ALPHA={B.ALPHA}  CD4_TARGET_F={B.CD4_TARGET_F}  "
          f"PRED_TRACE_CAP_K={PRED_TRACE_CAP_K}", flush=True)

    ctrl = np.load(REPO / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cpn = dict(np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    pre = np.load(REPO / "experiments/170-gamma-mse/work/pre.npz", allow_pickle=True)
    genes = pre["genes"].astype(str)
    perts = pre["perts"].astype(str)
    assert np.array_equal(genes, ctrl["genes"].astype(str))

    bctrl = {c: bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
             for c in CTX}
    Kc = {c: 300 * float(ctrl[f"jk_{c}"]) / 18532 for c in CTX}
    Cp = {c: float(np.mean(cpn[c])) for c in CTX}
    Lmed = context_library_medians()
    for i, c in enumerate(CTX):
        rel = abs(Kc[c] / KC_KA[i] - 1.0)
        print(f"[chal] Kc[{c}] {Kc[c]:.6f} vs 096 (e) {KC_KA[i]}  rel {rel:.2e}  "
              f"{'PASS' if rel < 1e-4 else 'FAIL'}")
    cap = float(np.mean([perts.size * Cp[c] / (genes.size - 1) for c in CTX]))
    print(f"[chal] C_p proxy {[round(Cp[c], 4) for c in CTX]}  cap {cap:.5f} vs 0.24653; "
          f"median control library {[int(Lmed[c]) for c in CTX]}", flush=True)

    gpos = {g: i for i, g in enumerate(genes)}
    cols = {i: gpos[p] for i, p in enumerate(perts) if p in gpos}
    w = _row_weights(perts.size, genes.size, cols)
    sub = np.linspace(0, perts.size - 1, NSUB).astype(int)

    out_path = WORK / f"chal212_{tag}.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {r["gamma"] for r in out}
    for gamma in GRID:
        if gamma in done:
            continue
        t0 = time.time()
        eff = np.load(WORK / "eff" / f"chal_g{gamma:.3f}.npz")["effect"]
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
                    eff[j], b, min_cpm=MIN_CPM, tail_spread=TAIL_SPREAD,
                    keep_top=KEEP_TOP, keep_spread=KEEP_SPREAD, spread_map="block")
                v, cl, ip = emit_counts(P[j], libs, DISPERSION,
                                        np.random.default_rng(3000 + int(j)), mode="quasi",
                                        overdispersion=OVERDISPERSION, spread=spread)
                Xm = csr_matrix((v.astype(np.float64), cl.astype(np.int64),
                                 ip.astype(np.int64)), shape=(NCELL, genes.size))
                vals.append(jackknife_correction(Xm, np.zeros(NCELL, np.intp), 1, TS,
                                                 chunk=64)[0])
            chat = float(np.mean(vals))
            capped = PRED_TRACE_CAP_K * Cp[c]
            claim = float(perts.size * np.mean(np.minimum(vals, capped)) / (genes.size - 1))
            per[c] = dict(A=A, B=bud, claim=claim, chat=chat, Kc=Kc[c], Cp=Cp[c],
                          chat_over_cp=chat / Cp[c],
                          rows_under_cap=int(np.sum(np.asarray(vals) < capped)),
                          claim_cap=float(perts.size * capped / (genes.size - 1)),
                          # the rho = 1 cliff, in units of the shipped amplitude:
                          s_cliff=float(np.sqrt(claim / bud)) if bud > 0 else np.nan,
                          rho=min(1.0, bud / claim) if claim else 1.0)
            del P, pm, dh
        row = dict(gamma=gamma, per={c: per[c] for c in CTX},
                   A=float(np.mean([per[c]["A"] for c in CTX])),
                   B=float(np.mean([per[c]["B"] for c in CTX])),
                   claim=float(np.mean([per[c]["claim"] for c in CTX])),
                   chat=float(np.mean([per[c]["chat"] for c in CTX])),
                   mean_sqrtA=float(np.mean([np.sqrt(per[c]["A"]) for c in CTX])),
                   secs=round(time.time() - t0, 1))
        out.append(row)
        out.sort(key=lambda r: r["gamma"])
        print(f"[chal] gamma {gamma:.3f}  A {row['A']:.6f}  B {row['B']:.6f}  "
              f"claim {row['claim']:.6f}  C_hat {row['chat']:.3f}  "
              f"rho {[round(per[c]['rho'], 4) for c in CTX]}  "
              f"s_cliff {[round(per[c]['s_cliff'], 4) for c in CTX]}  [{row['secs']:.0f}s]",
              flush=True)
        out_path.write_text(json.dumps(out, indent=1))

    base = [r for r in out if r["gamma"] == 0.65][0]
    print(f"\n[chal] CLAIM KNOWN-ANSWER: gamma 0.65 claim {base['claim']:.5f} vs 164's "
          f"0.24581 / 170's 0.24610 (cap {cap:.5f} vs 0.24653); rho "
          f"{[round(base['per'][c]['rho'], 4) for c in CTX]}, mean "
          f"{np.mean([base['per'][c]['rho'] for c in CTX]):.4f} vs 170's 0.938 / "
          f"164's 0.9466")

    target_raw = raw_from_score(A23_SCORE_MSE if massfix else A15_SCORE_MSE)
    num = float(np.mean([base["per"][c]["A"] + base["per"][c]["Kc"]
                         - min(base["per"][c]["claim"], base["per"][c]["B"]) for c in CTX]))
    X = (num - SUMD * (target_raw - 1.0)) / 2.0
    cos = X / (np.sqrt(SUMD) * base["mean_sqrtA"])
    print(f"[chal] board anchor: {'a23' if massfix else 'a15'} score_mse "
          f"{A23_SCORE_MSE if massfix else A15_SCORE_MSE} -> raw {target_raw:.6f}")
    print(f"[chal] solved X at gamma 0.65 = {X:.6f}; implied transfer cosine {cos:.6f} "
          f"(197's a15 measurement: 0.079249)")
    with open(WORK / f"chal212_{tag}_solve.json", "w") as fh:
        json.dump({"X_at_065": X, "cos_at_065": float(cos), "sumD": SUMD,
                   "target_raw": target_raw, "cap": cap,
                   "anchor": "a23" if massfix else "a15"}, fh, indent=1)


if __name__ == "__main__":
    main()

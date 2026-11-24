"""213 — the challenge panel's `A`, `B`, `claim`, rho and the RE-SOLVED alpha, per law.

`experiments/212-gamma-cosine/chal212.py`, keyed by law instead of gamma, with the RENDER
replaced by the law's own (`law213.deliver`-equivalent, a15 geometry -- no mass fix, which is
what the incumbent ships).  Everything else -- 096's control statistics, the `C_p` proxy,
`cell_eval2`'s own `_across_pert_budget` / `_row_weights` / `jackknife_correction`, a15's
emission policy, the 24-target subsample and its seeds -- is 212's, unchanged.

KNOWN-ANSWER ROWS (rule 10), all at `logc20`:
  * per-context `Kc` == 096 (e)'s 0.005084 / 0.005228 / 0.005919
  * the `C_p` cap == 0.24653  (212 gate 7)
  * `claim` == 212's 0.24608 / 170's 0.24610 / 164's 0.24581   (+-0.5%)
  * mean rho == 170's 0.938 / 164's 0.9466, per context 0.9926 / 0.9407 / 0.8821 (212 2.3)
  * the solved `X` implies transfer cosine 0.0773 (212) / 0.079249 (197) -- 212 documents the
    -2.45% level offset between the two and it is a level, not a ratio.

    uv run python chal213.py --laws logc20 logc40
"""
from __future__ import annotations

import argparse
import gc
import importlib.util
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
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments/096-mse-floor"))
sys.path.insert(0, str(HERE))
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")

import law213 as L  # noqa: E402

TS = 50_000.0
NCELL = 400
MIN_CPM = 10.0
TAIL_SPREAD, KEEP_TOP, KEEP_SPREAD = 1.0, 8000, 0.25
OVERDISPERSION, DISPERSION = 4.0, 0.1
CTX = ("A", "B", "C")
NSUB = 24
SUMD = 0.229141                       # 187's `C_p`-free LS solve (096/197/212 all use it)
MSE_MAP = (0.987862879, 0.950162566)  # raw = MSE_MAP[0] - score*MSE_MAP[1]
A15_SCORE_MSE = 0.1200566151966586
KC_KA = (0.005084, 0.005228, 0.005919)
SPAN_MSE = MSE_MAP[1]


def builder():
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_213c", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_213c"] = mod
    spec.loader.exec_module(mod)
    return mod


def context_library_medians():
    cache = WORK / "chal_libs.json"
    if cache.exists():
        return {k: float(v) for k, v in json.loads(cache.read_text()).items()}
    src = REPO / "experiments/170-gamma-mse/work/chal_libs.json"
    d = {k: float(v) for k, v in json.loads(src.read_text()).items()}
    cache.write_text(json.dumps(d, indent=1))
    return d


def solve_alpha(per, X, sumd=SUMD, lo=0.5, hi=1.6, n=8001):
    """Minimise the exact per-context numerator over a global amplitude `s`.

    `min(claim_c, s^2 B_c)` is evaluated per context at every `s` -- 170's own correction,
    that linearising the `min` around the shipped point inverts the result.
    """
    s = np.linspace(lo, hi, n)
    num = np.zeros_like(s)
    for c in CTX:
        A, Bc, claim, Kc = per[c]["A"], per[c]["B"], per[c]["claim"], per[c]["Kc"]
        num += (s * s * A + Kc - np.minimum(claim, s * s * Bc) - 2.0 * s * X) / len(CTX)
    raw = 1.0 + num / sumd
    i = int(np.argmin(raw))
    return float(s[i]), float(raw[i]), float(1.0 + num[np.argmin(np.abs(s - 0.982))] / sumd)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["logc20", "logc40"])
    ap.add_argument("--cos", nargs="*", type=float, default=None,
                    help="board transfer cosine per law; default 197's 0.079249 x the "
                         "PLOLO ratio read from work/plolo213_main.json")
    args = ap.parse_args()

    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction
    from cell_eval2.prep import bulk_lognorm_means
    from vcc2026.donor_transfer import emit_counts

    B = builder()
    if B.MASS_FIX is not None:
        raise SystemExit("this measures the INCUMBENT's geometry: leave VCC_MASS_FIX unset")
    laws = L.all_laws()

    ctrl = np.load(REPO / "experiments/096-mse-floor/work/ctrl-stats.npz", allow_pickle=True)
    cpn = dict(np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    genes = np.load(WORK / "pre_logc20.npz", allow_pickle=True)["genes"].astype(str)
    perts = np.load(WORK / "pre_logc20.npz", allow_pickle=True)["perts"].astype(str)
    assert np.array_equal(genes, ctrl["genes"].astype(str))
    bctrl = {c: bulk_lognorm_means(ctrl[f"sum_{c}"][None, :].astype(np.float64), TS)[0]
             for c in CTX}
    Kc = {c: 300 * float(ctrl[f"jk_{c}"]) / 18532 for c in CTX}
    Cp = {c: float(np.mean(cpn[c])) for c in CTX}
    Lmed = context_library_medians()
    for i, c in enumerate(CTX):
        rel = abs(Kc[c] / KC_KA[i] - 1.0)
        print(f"[chal213] Kc[{c}] {Kc[c]:.6f} vs 096 (e) {KC_KA[i]}  rel {rel:.2e}  "
              f"{'PASS' if rel < 1e-4 else 'FAIL'}")
    cap = float(np.mean([perts.size * Cp[c] / (genes.size - 1) for c in CTX]))
    print(f"[chal213] C_p cap {cap:.5f} vs 212's 0.24653", flush=True)

    gpos = {g: i for i, g in enumerate(genes)}
    cols = {i: gpos[p] for i, p in enumerate(perts) if p in gpos}
    w = _row_weights(perts.size, genes.size, cols)
    sub = np.linspace(0, perts.size - 1, NSUB).astype(int)

    # PLOLO cosine ratios -> board cosine level, 197's 0.079249 at the incumbent
    plolo = json.loads((WORK / "plolo213_main.json").read_text())
    F = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
    macro = {}
    for k in set(r["law"] for r in plolo):
        v = {r["fold"]: r["a15_cos"] for r in plolo if r["law"] == k}
        if all(f in v for f in F):
            macro[k] = float(np.mean([v[f] for f in F]))

    out_path = WORK / "chal213.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {r["law"] for r in out}

    for key in args.laws:
        if key in done:
            continue
        law = laws[key]
        eff = np.load(WORK / "eff" / f"chal_{key}.npz")["effect"]
        t0 = time.time()
        per = {}
        for c in CTX:
            b = bctrl[c]
            p_ctrl = np.expm1(b) / TS
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > MIN_CPM) & (b > 0)
            P = np.empty_like(eff)
            for i in range(eff.shape[0]):
                g = np.where(gate, eff[i], 0.0)
                pp = np.clip(law.render(g, p_ctrl, b), 0.0, None)
                P[i] = pp / float(pp.sum())
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
                          s_cliff=float(np.sqrt(claim / bud)) if bud > 0 else float("nan"),
                          rho=min(1.0, bud / claim) if claim else 1.0,
                          median_row_norm=float(np.median(np.linalg.norm(dh, axis=1))))
            del P, pm, dh
            gc.collect()
        mean_sqrtA = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
        row = dict(law=key, per=per,
                   A=float(np.mean([per[c]["A"] for c in CTX])),
                   B=float(np.mean([per[c]["B"] for c in CTX])),
                   claim=float(np.mean([per[c]["claim"] for c in CTX])),
                   chat=float(np.mean([per[c]["chat"] for c in CTX])),
                   mean_sqrtA=mean_sqrtA, cap=cap,
                   plolo_cos_macro=macro.get(key), secs=round(time.time() - t0, 1))
        out.append(row)
        out_path.write_text(json.dumps(out, indent=1))
        print(f"[chal213] {key}  A {row['A']:.6f}  B {row['B']:.6f}  claim "
              f"{row['claim']:.6f}  C_hat {row['chat']:.3f}  rho "
              f"{[round(per[c]['rho'], 4) for c in CTX]}  s_cliff "
              f"{[round(per[c]['s_cliff'], 4) for c in CTX]}  [{row['secs']:.0f}s]",
              flush=True)
        del eff
        gc.collect()

    # ------------------------------------------------------------------ solve and price
    base = [r for r in out if r["law"] == "logc20"]
    if not base:
        return
    base = base[0]
    target_raw = MSE_MAP[0] - A15_SCORE_MSE * MSE_MAP[1]
    num = float(np.mean([base["per"][c]["A"] + base["per"][c]["Kc"]
                         - min(base["per"][c]["claim"], base["per"][c]["B"]) for c in CTX]))
    X0 = (num - SUMD * (target_raw - 1.0)) / 2.0
    cos0 = X0 / (np.sqrt(SUMD) * base["mean_sqrtA"])
    print(f"\n[chal213] board anchor a15 score_mse {A15_SCORE_MSE} -> raw {target_raw:.6f}")
    print(f"[chal213] solved X {X0:.6f}; implied transfer cosine {cos0:.6f} "
          f"(212: 0.077305, 197: 0.079249)")
    print(f"[chal213] claim {base['claim']:.5f} (212 0.24608 / 170 0.24610 / 164 0.24581); "
          f"rho {[round(base['per'][c]['rho'], 4) for c in CTX]} "
          f"(212 0.9926/0.9407/0.8821)")

    priced = []
    for r in out:
        # the law's board cosine, transferred by the SCALE-INVARIANT PLOLO ratio (212 3.5)
        ratio = (r["plolo_cos_macro"] / base["plolo_cos_macro"]) if base["plolo_cos_macro"] \
            else 1.0
        cos = cos0 * ratio
        X = cos * np.sqrt(SUMD) * r["mean_sqrtA"]
        s_star, raw_star, raw_at_982 = solve_alpha(r["per"], X)
        s0, raw0, raw0_982 = solve_alpha(base["per"], X0)
        priced.append(dict(
            law=r["law"], cos_ratio=ratio, board_cos=cos, X=X,
            alpha_star=s_star, raw_mse_at_alpha_star=raw_star,
            raw_mse_at_982=raw_at_982,
            d_overall_mse_at_star=(raw0 - raw_star) / SPAN_MSE / 6.0,
            d_overall_mse_at_982=(raw0_982 - raw_at_982) / SPAN_MSE / 6.0,
            rA=r["A"] / base["A"], rNorm=(r["per"]["A"]["median_row_norm"]
                                          / base["per"]["A"]["median_row_norm"]),
        ))
    (WORK / "chal213_price.json").write_text(json.dumps(
        dict(base_cos=cos0, base_X=X0, sumD=SUMD, span_mse=SPAN_MSE, rows=priced), indent=1))
    print(f"\n{'law':<10}{'cos ratio':>10}{'board cos':>11}{'alpha*':>9}"
          f"{'raw mse*':>10}{'d ovr mse*':>12}{'rA':>8}")
    for p in priced:
        print(f"{p['law']:<10}{p['cos_ratio']:10.4f}{p['board_cos']:11.6f}"
              f"{p['alpha_star']:9.4f}{p['raw_mse_at_alpha_star']:10.6f}"
              f"{p['d_overall_mse_at_star']:+12.5f}{p['rA']:8.4f}")


if __name__ == "__main__":
    main()

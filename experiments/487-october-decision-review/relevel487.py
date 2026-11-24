"""487 decision 4: sensitivity of the final score to an error in the D/E/F mse-cliff estimate, in each direction.

Model (ntclevel417 / 187): near context c's cliff s*, raw mse changes by -0.31 per unit scale below the cliff and +1.78
above (A/B/C slopes). The shipped rule levels to rho* = exp(2 sd Phi^-1(q)), q = 0.31/2.09, sd = 0.01 in ln scale (an
unfitted context), i.e. a margin m = sd Phi^-1(q) in ln scale; the secant loop stops within 1 % of rho* in rho
(+-0.5 % in scale). One raw mse unit in one context = 1 / (0.9502 * 18) overall (scaled mse span 0.950162566, 187;
6 members x 3 contexts).
Known answer: level_rho(0.01) = 0.97933 (DEF-RUNBOOK's 0.9793; lib438.level_rho).
usage: uv run python experiments/487-october-decision-review/relevel487.py -> work/relevel487.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
A, B = 0.31, 1.78                 # raw mse per unit scale, below / above the cliff
PER_RAW = 1.0 / (0.950162566 * 18)  # overall per raw-mse unit in one context
q = A / (A + B)


def margin(sd):
    return sd * norm.ppf(q)       # ln-scale position of the shipped point relative to the (proxy) cliff


def loss(pos):
    """raw mse lost vs sitting exactly on the true cliff; pos = ln(s_shipped / s*_true)."""
    return np.where(pos > 0, B * pos, -A * pos)


def main():
    out = {"q": q, "ka_rho_star_sd0.01": float(np.exp(2 * margin(0.01)))}
    assert abs(out["ka_rho_star_sd0.01"] - 0.9793) < 5e-5, out
    m = margin(0.01)
    print(f"KA: rho* = {out['ka_rho_star_sd0.01']:.5f} (margin {100 * m:+.3f} % in scale)")
    # deterministic proxy error: the true cliff sits at delta (ln scale) relative to the proxy cliff
    rows = []
    for d in (-0.30, -0.10, -0.05, -0.03, -0.02, -0.0088, 0.0, 0.0088, 0.02, 0.03, 0.05, 0.10, 0.293):
        for tol in (-0.005, 0.0, 0.005):     # where inside the secant loop's +-1 % rho band the build landed
            pos = m + tol - d
            L = float(loss(pos))
            rows.append({"delta_pct": round(100 * d, 2), "tol_pct": round(100 * tol, 2), "pos_pct": round(100 * pos, 2),
                         "raw_mse_per_context": round(L, 5), "overall_per_context": round(L * PER_RAW, 5)})
    out["deterministic"] = rows
    print("true cliff vs proxy (delta, % scale; negative = cliff LOWER than the proxy says), shipped at the loop's centre:")
    for r in rows:
        if r["tol_pct"] == 0.0:
            print(f"  delta {r['delta_pct']:+7.2f} %: shipped {r['pos_pct']:+6.2f} % from the true cliff, raw mse +{r['raw_mse_per_context']:.4f}, "
                  f"overall -{r['overall_per_context']:.5f} per context (x3 if all three share it: -{3 * r['overall_per_context']:.4f})")
    # expected loss when the proxy error is random, delta ~ N(bias, sd_true), vs the shipped margin and the optimal margin
    ev = []
    for bias in (0.0, -0.0088):
        for sd_true in (0.0025, 0.005, 0.01, 0.02, 0.03):
            z = np.random.default_rng(0).normal(bias, sd_true, 400000)
            shipped = float(loss(m - z).mean())
            grid = np.linspace(-0.08, 0.02, 2001)
            best = min(((float(loss(g - z).mean()), g) for g in grid))
            ev.append({"bias_pct": round(100 * bias, 2), "sd_true_pct": round(100 * sd_true, 2),
                       "E_raw_shipped": round(shipped, 5), "E_overall_shipped_per_ctx": round(shipped * PER_RAW, 6),
                       "best_margin_pct": round(100 * best[1], 3), "E_overall_best_per_ctx": round(best[0] * PER_RAW, 6),
                       "regret_overall_per_ctx": round((shipped - best[0]) * PER_RAW, 6)})
    out["expected"] = ev
    print("random proxy error delta ~ N(bias, sd): expected loss per context (overall units) and regret vs the best margin")
    for r in ev:
        print(f"  bias {r['bias_pct']:+.2f} % sd {r['sd_true_pct']:.2f} %: shipped -{r['E_overall_shipped_per_ctx']:.5f}, "
              f"best margin {r['best_margin_pct']:+.2f} % -> -{r['E_overall_best_per_ctx']:.5f}, regret {r['regret_overall_per_ctx']:.5f}")
    # slope uncertainty: if the above-cliff slope on D/E/F were x times A/B/C's, what rho* is optimal and what is the regret?
    sl = []
    for f in (0.5, 1.0, 2.0, 4.0):
        qq = A / (A + f * B)
        mm = 0.01 * norm.ppf(qq)
        z = np.random.default_rng(1).normal(0, 0.01, 400000)
        Lf = lambda pos: np.where(pos > 0, f * B * pos, -A * pos)
        sl.append({"above_slope_x": f, "rho_star_opt": round(float(np.exp(2 * mm)), 4),
                   "regret_overall_per_ctx_of_shipping_0.9793": round(float((Lf(m - z).mean() - Lf(mm - z).mean()) * PER_RAW), 6)})
    out["slope_sensitivity"] = sl
    print("above-cliff slope sensitivity:", json.dumps(sl))
    (HERE / "work/relevel487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

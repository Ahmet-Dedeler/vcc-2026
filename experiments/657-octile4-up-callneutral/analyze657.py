#!/usr/bin/env python3
"""657 — read the octile-4 UP routes on 618 against PREREGISTRATION.md (ed36f45f), with 661's scoring code path:
lib521 on 618's state (REF / EFF / ST435 -> ref_tv, the tv key, g2mde), price528.rows_at (realised order) at 521's regimes,
policy521.SGN (192's INV) for overall-eq; calls per target (q < 0.05, own gene excluded).

Base = c4h hook + 630 (stats-rnup_c4h_s<s>_b_k100_q20r3p10). Arms, paired per seed:
  A  = c4h_o4 emission + 630            (route A: the hook)
  B  = c4h emission + 630 with _o4      (route B: 630's selection)
  AB = c4h_o4 emission + 630 with _o4
  A0 = c4h_o4 emission vs c4h emission  (route A without 630; a component)
and, when Phase 2 exists, the same with DN (_D200) against base_dn = c4h + 630 + DN.
Gates (prereg): |dcalls| <= 10 per target (every seed); |dfid|, |djac| <= 0.0005 raw at f 0.35 and 0.5; downside: overall
>= -0.0005 at f 0.35 AND f 0.5 (mean over seeds) -> ELIGIBLE, else CLOSED.
Octile-4 diagnostics: on uncalled octile-4 gate cells of the truth, the delivered sign's agreement with the truth and the
mean delivered lfc, per arm.
usage: uv run python experiments/657-octile4-up-callneutral/analyze657.py"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W618 = ROOT / "experiments/618-h1-trainval-surface/work"
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
import lib521 as L  # noqa: E402

L.REF = W618 / "ref_tv/reference.npz"
L.EFF = W618 / "tv/effects.npz"
L.ST435 = W618 / "g2mde/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
L.load_state.cache_clear()
import price528 as P5  # noqa: E402
from policy521 import SGN, means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REG = ("1.0", "0.5", "0.35", "cmp")
MEM = ("fid", "jac", "reach", "nmae")
SEEDS = (0, 1, 2)
EM = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
ST = HERE / "work/stacks"
V8 = np.load(ROOT / "experiments/625-module-readout/work/modules625_v8.npz", allow_pickle=True)
OCT = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))


def load(path):
    if not Path(path).exists():
        return None
    z = np.load(path, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float), sums=np.asarray(z["sums"]))


def arms(s):
    st = lambda nm, v: ST / f"stats-rnup_{nm}_s{s}_{v}.npz"  # noqa: E731
    return {
        "base": st("c4h", "b_k100_q20r3p10"),
        "A": st("o4", "b_k100_q20r3p10"),
        "B": st("c4h", "b_k100_q20r3p10_o4"),
        "AB": st("o4", "b_k100_q20r3p10_o4"),
        "em_c4h": HERE / "work/c4h" / EM.format(s=s),
        "A0": HERE / "work/c4h_o4" / EM.format(s=s),
        "base_dn": st("c4h", "b_k100_q20r3p10_D200"),
        "A_dn": st("o4", "b_k100_q20r3p10_D200"),
        "B_dn": st("c4h", "b_k100_q20r3p10_D200_o4"),
        "AB_dn": st("o4", "b_k100_q20r3p10_D200_o4"),
    }


COMPARE = (("A", "base"), ("B", "base"), ("AB", "base"), ("A0", "em_c4h"),
           ("A_dn", "base_dn"), ("B_dn", "base_dn"), ("AB_dn", "base_dn"))


def main():
    states = {s: L.load_state(s) for s in SEEDS}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in SEEDS for r in REG}
    own = states[0]["own"]
    genes = np.asarray(states[0]["genes"]).astype(str)
    o4 = np.array([OCT.get(g, -1) == 4 for g in genes])
    rl = states[0]["rl"]
    runs = {s: {k: load(p) for k, p in arms(s).items()} for s in SEEDS}
    out, lines = {}, []
    say = lambda x="": (print(x), lines.append(x))  # noqa: E731
    say("arm vs ref | seeds | dcalls/target (per seed) | sums==ref | per regime: overall-eq [fid jac reach nmae] (mean over seeds)")
    for arm, ref in COMPARE:
        seeds = [s for s in SEEDS if runs[s][arm] is not None and runs[s][ref] is not None]
        if not seeds:
            continue
        rec = {"seeds": seeds}
        dc = [float(((runs[s][arm]["q"] < 0.05) & ~own).sum(1).mean() - ((runs[s][ref]["q"] < 0.05) & ~own).sum(1).mean()) for s in seeds]
        rec["d_calls_per_target"] = dc
        rec["sums_equal_ref"] = all(np.array_equal(runs[s][arm]["sums"], runs[s][ref]["sums"]) for s in seeds)
        for r in REG:
            dms = np.stack([means(P5.rows_at(states[s], runs[s][arm], pools[(s, r)])) - means(P5.rows_at(states[s], runs[s][ref], pools[(s, r)]))
                            for s in seeds])
            rec[r] = dict(zip(MEM, map(float, dms.mean(0))), overall=float(np.mean(dms @ SGN)),
                          overall_seeds=[float(x) for x in dms @ SGN],
                          sd=dict(zip(MEM, map(float, dms.std(0, ddof=1) if len(seeds) > 1 else [np.nan] * 4))))
        # octile-4 diagnostics on uncalled octile-4 gate cells of the truth (f 1 pool), arm vs ref, seed-mean
        ag, mg = [], []
        for s in seeds:
            for which in (arm, ref):
                A = runs[s][which]
                m = pools[(s, "1.0")] & o4[None, :] & ~((A["q"] < 0.05)) & ~own & np.isfinite(rl) & (rl != 0)
                ag.append((which, float((np.sign(np.nan_to_num(A["lfc"])) == np.sign(rl))[m].mean()) if m.any() else np.nan))
                mg.append((which, float(np.nan_to_num(A["lfc"])[m].mean()) if m.any() else np.nan))
        rec["o4_uncalled_sign_agree"] = {w: float(np.nanmean([v for x, v in ag if x == w])) for w in (arm, ref)}
        rec["o4_uncalled_mean_lfc"] = {w: float(np.nanmean([v for x, v in mg if x == w])) for w in (arm, ref)}
        gates = dict(calls=all(abs(x) <= 10 for x in dc),
                     fid_jac=all(abs(rec[r][m]) <= 0.0005 for r in ("0.35", "0.5") for m in ("fid", "jac")),
                     downside=all(rec[r]["overall"] >= -0.0005 for r in ("0.35", "0.5")))
        rec["gates"] = gates
        rec["verdict"] = "ELIGIBLE" if all(gates.values()) else "CLOSED"
        out[f"{arm}|{ref}"] = rec
        say(f"{arm:6s} vs {ref:8s} | {seeds} | dcalls {np.mean(dc):+.1f} {[round(x, 1) for x in dc]} | sums==ref {rec['sums_equal_ref']} | "
            + " | ".join(f"f {r}: {rec[r]['overall']:+.5f} [{rec[r]['fid']:+.5f} {rec[r]['jac']:+.5f} {rec[r]['reach']:+.5f} {rec[r]['nmae']:+.5f}]" for r in REG))
        say(f"        o4 uncalled: sign agreement {rec['o4_uncalled_sign_agree'][arm]:.3f} (ref {rec['o4_uncalled_sign_agree'][ref]:.3f}), mean delivered lfc "
            f"{rec['o4_uncalled_mean_lfc'][arm]:+.4f} (ref {rec['o4_uncalled_mean_lfc'][ref]:+.4f}) | gates {gates} -> {rec['verdict']}")
    (HERE / "work/analyze657.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/analyze657.out").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""706: where does the bulk's remaining error sit by target strength, and how much is AMPLITUDE vs DIRECTION?

Motivation: the board's leader (Illumina AI) sits at mse 0.79 scaled with a field-level pds. mse is Σ N_p / Σ D_p, so the
strong targets carry it. The question is whether our per-target amplitude allocation (the LODO key, percentile cap) leaves
the strong tier under- or over-weighted on held-out lines. If it does not, the remaining strong-tier gap is direction.

Surfaces: 522's exact fold surfaces (fold522.load, its own known answers asserted on load):
  - 9 October-like folds (efp / efr / real × held-out K562-GWPS / Orion HCT116 / Orion HEK293T),
  - 4 A/B/C PLOLO folds (held-out k562 / hct116 / hek293t / kolf2).
Arm: production's key, the percentile cap P 96.64697522769474 (q 1.75 row sharpening, render, ballast KMAX 16), and the
same with 522's θ 0.25 floor-preserving rider (705's arm), for reference.
Per frame (A/B/C slots) the realised rows give per-target x_p = <P_p, T_p>, a_p = |P_p|², d_p = |T_p|² on the scored genes.
cos_agg X = Σx / sqrt(Σa Σd) (the instrument's transfer quantity; 518: +1 % X ≈ +0.0007 overall).
Reported, frame-averaged:
  - X; the per-tier cosine Σx_t / sqrt(Σa_t Σd_t); each tier's share of Σd, Σa and Σx;
  - k_t: tier t's optimal common scale relative to the panel's, (Σx_t/Σa_t) / (Σx/Σa). k > 1 means under-weighted;
  - X with a per-tier optimal scale (3 numbers, Σa fixed): the size of any tier-amplitude lever;
  - X with an oracle per-target scale (s_p ∝ x_p/a_p for x_p > 0, else 0; Σa fixed): the amplitude ceiling.
Tiers: tertiles of z_src (the donors' strength rank: predictable on the day), of the key weight w, and of truth energy d
(an oracle, for orientation only).
usage: uv run python experiments/706-strong-tier-budget/budget706.py [panel ...]     (~1-2 GB; a few min per fold)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402

OUT = HERE / "work"
OUT.mkdir(parents=True, exist_ok=True)


def tiers(v, n=3):
    """tertile labels 0 (low) .. n-1 (high) by rank; ties broken by order."""
    r = np.argsort(np.argsort(v, kind="stable"), kind="stable")
    return np.minimum((r * n) // len(v), n - 1)


def budget(x, a, d, lab):
    """frame-level budget for one labelling; x, a are (frames, n), d is (n,)."""
    out = {"X": [], "X_tier_oracle": [], "X_target_oracle": []}
    per = {t: {"cos": [], "k": [], "share_d": [], "share_a": [], "share_x": []} for t in range(3)}
    for c in range(x.shape[0]):
        xc, ac = x[c], a[c]
        Sx, Sa, Sd = xc.sum(), ac.sum(), d.sum()
        out["X"].append(Sx / np.sqrt(Sa * Sd))
        g = Sx / Sa
        num = 0.0
        for t in range(3):
            m = lab == t
            xt, at, dt = xc[m].sum(), ac[m].sum(), d[m].sum()
            per[t]["cos"].append(xt / np.sqrt(at * dt))
            per[t]["k"].append((xt / at) / g)
            per[t]["share_d"].append(dt / Sd)
            per[t]["share_a"].append(at / Sa)
            per[t]["share_x"].append(xt / Sx)
            num += max(xt, 0.0) ** 2 / at
        out["X_tier_oracle"].append(np.sqrt(num / Sd))
        pos = xc > 0
        out["X_target_oracle"].append(np.sqrt((xc[pos] ** 2 / ac[pos]).sum() / Sd))
    res = {k: float(np.mean(v)) for k, v in out.items()}
    res["tiers"] = {str(t): {k: float(np.mean(v)) for k, v in per[t].items()} for t in range(3)}
    return res


def main():
    panels = sys.argv[1:] or ["efp", "efr", "real", "abc"]
    allres = {}
    for panel, held in FD.surfaces(tuple(panels)):
        t0 = time.time()
        S = FD.load(panel, held, check=True)
        kv, car = S.kv, S.carried
        wp = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_PROD))
        z = FD.zsrc(S)
        wr = np.where(wp <= 0.5 + 1e-12, wp, np.maximum(0.5, wp * np.exp(0.25 * z)))
        rec = {"n": int(S.n), "ka": {k: (float(v) if not isinstance(v, str) else v) for k, v in S.ka.items()}}
        for arm, w in (("pct", wp), ("pct_rider025", wr)):
            s, arr = S.score(w)
            x, a, d = np.asarray(arr["x"], float), np.asarray(arr["a"], float), np.asarray(arr["d"], float)
            live = d > 0
            x, a, d = x[:, live], a[:, live], d[live]
            r = {"cos_agg": float(s["cos_agg"]), "pds": float(s["pds"]), "n_live": int(live.sum())}
            for name, v in (("z_src", z[live]), ("w", w[live]), ("truth_d", d)):
                r[f"by_{name}"] = budget(x, a, d, tiers(v))
            rec[arm] = r
        rec["seconds"] = round(time.time() - t0, 1)
        allres[f"{panel}/{held}"] = rec
        p = rec["pct"]
        print(f"{panel}/{held}: X {p['by_z_src']['X']:.4f} | tier-oracle(z) {p['by_z_src']['X_tier_oracle']:.4f} "
              f"| target-oracle {p['by_z_src']['X_target_oracle']:.4f} | k by z (lo/mid/hi) "
              f"{[round(p['by_z_src']['tiers'][t]['k'], 3) for t in '012']} | k by truth "
              f"{[round(p['by_truth_d']['tiers'][t]['k'], 3) for t in '012']} | rider X "
              f"{rec['pct_rider025']['by_z_src']['X']:.4f} | {rec['seconds']} s", flush=True)
        (OUT / "budget706.json").write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()

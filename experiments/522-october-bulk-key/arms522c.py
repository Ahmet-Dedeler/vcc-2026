"""522c (EXPLORATORY) — the rider with 479's P2g partner. P2g (2.5 / 2.0 / 1.75 by the key's tercile, every target at the
key floor -> LOW) reads the keyed/unkeyed row ratio, so a rider that multiplies floor targets by exp(0.25 z) in
[0.88, 1.13] breaks the floor block (the builder's guard catches only the single minimum ratio). The floor-preserving
rider (fp) multiplies only targets above the floor (w > 0.5) and leaves the floor block intact; fpc (the form
def_strength_key522.py writes) also clips at the floor, so no target drops below the floor block (the builder's P2g floor is
the MINIMUM ratio, and its tie rule would send a floor block that contains a cut to MID).
Arms: pct; pct x S0.25; pct x S0.25fp / fpc; with P2g: pct + P2g, pct x S0.25 + P2g (naive), fp / fpc + P2g; with S_P2:
pct x S0.25fp / fpc + S_P2 (S_P2 reads the unkeyed rows, so it is untouched by any rider).
KNOWN ANSWERS (abort unless |d| < 1e-12): Oct folds: `p96.6 + P2g` (P 96.647, comp546o's tercile rule) == 546o's row;
A/B/C folds: `cap 6 + P2g` (arms540c's terc rule) == 540c's `V1_0.5_6 + P2g` row; pct and pct x S0.25 == arms522's rows.

    uv run python experiments/522-october-bulk-key/arms522c.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402
import lib438 as L  # noqa: E402
from sp2_522 import P2, score_q, terc  # noqa: E402

WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
FLOOR = 0.5
TOL = 1e-12


def p2g_oct(W, carried):
    """comp546o's P2g: carried-only cut points, floor -> LOW, non-carried -> MID."""
    fl = W <= FLOOR + 1e-12
    ec = np.quantile(W[carried], [1 / 3, 2 / 3])
    tk = np.where(fl, 0, np.where(~carried, 1, np.digitize(W, ec)))
    return np.asarray([P2[k] for k in tk], float)


def p2g_abc(W):
    """arms540c's P2g: terc(W, P2, floor 0.5)."""
    t = np.digitize(W, np.quantile(W, [1 / 3, 2 / 3]))
    t = np.where(W <= FLOOR + 1e-12, 0, t)
    return np.asarray([P2[k] for k in t], float)


def main():
    out = {}
    for p, f in FD.surfaces(("abc", "efp", "efr", "real")):
        if f == "kolf2":
            continue
        t0 = time.time()
        S = FD.load(p, f)
        z = FD.zsrc(S)
        kv, car = S.kv, S.carried
        wp = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_PROD))
        m = np.exp(0.25 * z)
        wr = wp * m
        wfp = np.where(wp > FLOOR + 1e-12, wp * m, wp)
        wfpc = np.where(wp > FLOOR + 1e-12, np.maximum(FLOOR, wp * m), wp)    # also never below the floor (P2g's vmin)
        p2g = (lambda W: p2g_oct(W, car)) if S.kind == "oct" else p2g_abc
        strength = np.linalg.norm(S.fin[:, S.nonpanel], axis=1) if S.kind == "oct" else np.linalg.norm(S.F.fST[:, S.F.nonpanel], axis=1)
        qs = terc(strength, P2)
        q1 = np.full(S.n, 1.75)
        ka = {}
        if S.kind == "oct":
            w546 = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_546))
            s, _ = score_q(S, w546, p2g(w546))
            ref = json.loads((FD.W546 / f"comp546o-{p}-{f}.json").read_text())["arms"]["p96.6 + P2g"]["realised"]["16"]
            d = L.cmp(s, ref)
            assert d < TOL, f"KA p96.6 + P2g {p}/{f}: {d}"
            ka["p96.6 + P2g == 546o"] = d
        else:
            w6 = FD.wclip(kv, car, 0.5, 6.0)
            s, _ = score_q(S, w6, p2g(w6))
            ref = json.loads((W540 / f"arms540c-{f}.json").read_text())["arms"]["V1_0.5_6 + P2g"]["realised"]["16"]
            d = max(abs(s["cos_agg"] - ref["cos_agg"]), abs(s["pds"] - ref["pds"]))
            assert d < TOL, f"KA cap6 + P2g {f}: {d}"
            ka["cap6 + P2g == 540c"] = d
        arms = {"pct": (wp, q1), "pctxS0.25": (wr, q1), "pctxS0.25fp": (wfp, q1),
                "pct + P2g": (wp, p2g(wp)), "pctxS0.25 + P2g (naive)": (wr, p2g(wr)), "pctxS0.25fp + P2g": (wfp, p2g(wfp)),
                "pct + S_P2": (wp, qs), "pctxS0.25fp + S_P2": (wfp, qs),
                "pctxS0.25fpc": (wfpc, q1), "pctxS0.25fpc + P2g": (wfpc, p2g(wfpc)), "pctxS0.25fpc + S_P2": (wfpc, qs)}
        ref522 = json.loads((WORK / f"arms522-{p}-{f}.json").read_text())["arms"]
        rec = {}
        for nm, (w, qv) in arms.items():
            s, _ = score_q(S, w, qv)
            rec[nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"])
        for nm, rn in (("pct", "pctPROD"), ("pctxS0.25", "pctPRODxS0.25")):
            d = abs(rec[nm]["cos_agg"] - ref522[rn]["cos_agg"])
            assert d < TOL, f"KA {nm} vs arms522 {p}/{f}: {d}"
        ka["pct, pct x S0.25 == arms522"] = 0.0
        c = {k: v["cos_agg"] for k, v in rec.items()}
        res = dict(ka=ka, arms=rec,
                   fp_vs_pct=c["pctxS0.25fp"] / c["pct"] - 1, rider_vs_pct=c["pctxS0.25"] / c["pct"] - 1,
                   P2g_on_pct=c["pct + P2g"] / c["pct"] - 1,
                   naive_rider_under_P2g=c["pctxS0.25 + P2g (naive)"] / c["pct + P2g"] - 1,
                   fp_rider_under_P2g=c["pctxS0.25fp + P2g"] / c["pct + P2g"] - 1,
                   fp_rider_under_S_P2=c["pctxS0.25fp + S_P2"] / c["pct + S_P2"] - 1,
                   fpc_vs_pct=c["pctxS0.25fpc"] / c["pct"] - 1,
                   fpc_rider_under_P2g=c["pctxS0.25fpc + P2g"] / c["pct + P2g"] - 1,
                   fpc_rider_under_S_P2=c["pctxS0.25fpc + S_P2"] / c["pct + S_P2"] - 1,
                   n_clipped_to_floor=int(np.sum((wp > FLOOR + 1e-12) & (wp * m < FLOOR))),
                   n_floor=int(np.sum(wp <= FLOOR + 1e-12)))
        out[f"{p}/{f}"] = res
        print(f"[{p}/{f}] KA {ka}; floor {res['n_floor']}; rider {100 * res['rider_vs_pct']:+.2f} fp {100 * res['fp_vs_pct']:+.2f}; "
              f"P2g on pct {100 * res['P2g_on_pct']:+.2f}; under P2g: naive {100 * res['naive_rider_under_P2g']:+.2f} "
              f"fp {100 * res['fp_rider_under_P2g']:+.2f}; fp under S_P2 {100 * res['fp_rider_under_S_P2']:+.2f}; "
              f"fpc {100 * res['fpc_vs_pct']:+.2f} |P2g {100 * res['fpc_rider_under_P2g']:+.2f} |S_P2 {100 * res['fpc_rider_under_S_P2']:+.2f} "
              f"(clipped {res['n_clipped_to_floor']}) ({time.time() - t0:.0f}s)",
              flush=True)
        (WORK / "arms522c.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

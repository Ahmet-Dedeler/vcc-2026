#!/usr/bin/env python3
"""537 — the called-gene magnitude scale s (513f; FINAL-PHASE-FEEDBACK-PLAN rows 8-9), priced analytically on H1.

Question: if a builder knob multiplied the DELIVERED lfc of every called gene by s (calls, signs and order held; pds/mse
pinned by the ballast), what would nmae do, and where is the optimum as the called genes' sign accuracy rises?
The board transferred every pure-magnitude lever at 1.6-2.3x H1 (docs/H1-BOARD-CALIBRATION.md), and 513f read that as a
sign that A/B/C's called-gene sign accuracy exceeds H1's. This script prices s on H1 (435 + zs = 446bj2z's DE design,
3 seeds, 521's regimes) as is, and under an emulation of higher sign accuracy: a fixed random fraction pi of the
wrong-signed called gate genes (the regime's pool) is given the right sign (the same draw for every s).
Known answer: s = 1, pi = 0 reproduces 435's zs nmae per seed (metrics json) at f 1.
Members other than nmae do not move with s (calls, signs, order held); with pi > 0 fid / reach move too (reported).
Writes work/scale537.json. Light (< 1 GB).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402
from table521 import regime_pool  # noqa: E402

SGRID = (0.5, 0.67, 0.8, 0.9, 1.0, 1.1, 1.25, 1.4, 1.6, 2.0)
PIS = (0.0, 0.1, 0.2, 0.35, 0.5)
REG = ("1.0", "0.5", "0.35", "cmp")
INV_NMAE = 1.64787


def main():
    out = {"KA": {}, "rows": []}
    for seed in (0, 1, 2):
        S = L.load_state(seed)
        A = P5.load_arm("base", seed, True)
        call = (A["q"] < L.ALPHA) & ~S["own"]
        rl = S["rl"]
        # KA: nmae at s 1 on the f 1 pool == 435's zs metrics
        nm = float(np.nanmean(L.nmae_rows(A["lfc"], L.pools(S, 1.0), rl)))
        ok = abs(nm - A["raw"]["nmae"]) < 1e-12
        out["KA"][f"s{seed}"] = [nm, A["raw"]["nmae"], ok]
        print(f"KA seed {seed}: nmae at s 1 {nm:.12f} vs 435 zs metrics {A['raw']['nmae']:.12f} -> {ok}", flush=True)
        assert ok
        u = np.random.default_rng(5370 + seed).random(call.shape)
        for reg in REG:
            pool = regime_pool(S, reg)[0]
            gate = pool & np.isfinite(rl) & (rl != 0)
            cg = call & gate
            right = np.sign(A["lfc"]) == np.sign(rl)
            a0 = float(right[cg].mean())
            for pi in PIS:
                fix = cg & ~right & (u < pi)
                lfc_pi = np.where(fix, -A["lfc"], A["lfc"])
                a = float((np.sign(lfc_pi) == np.sign(rl))[cg].mean())
                base = np.nanmean(L.nmae_rows(lfc_pi, pool, rl))
                for s in SGRID:
                    lf = np.where(call, lfc_pi * s, lfc_pi)
                    v = float(np.nanmean(L.nmae_rows(lf, pool, rl)))
                    out["rows"].append(dict(seed=seed, reg=reg, pi=pi, acc_called_gate=a, acc0=a0, s=s, nmae=v,
                                            d_vs_s1=v - base, overall_vs_s1=-(v - base) * INV_NMAE / 6))
        # where the gate's |truth| mass sits (f 1): called vs uncalled, and the delivered/true magnitude ratio when right
        g1 = L.pools(S, 1.0) & np.isfinite(rl) & (rl != 0)
        cg1 = call & g1
        mass_c = float(np.abs(rl[cg1]).sum() / np.abs(rl[g1]).sum())
        rr = (np.sign(A["lfc"]) == np.sign(rl))
        ratio = np.abs(A["lfc"][cg1 & rr]) / np.abs(rl[cg1 & rr])
        out.setdefault("anatomy", {})[f"s{seed}"] = dict(called_share_of_gate_mass=mass_c, acc_called_gate=float(rr[cg1].mean()),
                                                         acc_uncalled_gate=float(rr[g1 & ~call].mean()),
                                                         ratio_when_right_q25_50_75=np.quantile(ratio, [.25, .5, .75]).tolist())
        print(f"anatomy seed {seed}: {out['anatomy'][f's{seed}']}", flush=True)
    # summary: seed-mean overall vs s = 1 per (reg, pi, s), and the argmin s
    import collections
    agg = collections.defaultdict(list)
    for r in out["rows"]:
        agg[(r["reg"], r["pi"], r["s"])].append(r)
    summ = {}
    for reg in REG:
        for pi in PIS:
            vals = {s: float(np.mean([x["overall_vs_s1"] for x in agg[(reg, pi, s)]])) for s in SGRID}
            acc = float(np.mean([x["acc_called_gate"] for x in agg[(reg, pi, 1.0)]]))
            sbest = max(vals, key=vals.get)
            summ[f"{reg}|{pi}"] = dict(acc_called_gate=acc, best_s=sbest, best_overall=vals[sbest], curve=vals)
            print(f"f {reg:4s} pi {pi:.2f} (called-gate acc {acc:.3f}): best s {sbest:4.2f} ({vals[sbest]:+.5f} overall vs s 1) | "
                  + " ".join(f"s{s:g}:{vals[s]:+.5f}" for s in SGRID), flush=True)
    out["summary"] = summ
    (HERE / "work/scale537.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

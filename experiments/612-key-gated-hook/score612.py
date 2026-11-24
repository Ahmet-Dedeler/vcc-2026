#!/usr/bin/env python3
"""612: price the real-emitter H1 runs. Each run scored with its OWN realized order (scorer's rule) by 521's exact member
functions at its four regimes (f 1 / 0.5 / 0.35 / cmp, table521.regime_pool); 3 seeds; overall-equivalent with 192's
board-fitted INV (lib521.overall_delta). Comparisons: plain hook (537 cpm+c4h) vs base (435 = 446bj2z's construction),
key-gated (c4kg) vs base, identity key-gated (idkg) vs base, c4kg vs c4h, idkg vs c4kg. KA: c4kgall s0 == c4h s0 bit for
bit. Board expectation per 611 (RESULTS-611.md): nmae x1.2-1.6; reach read on the board (provisional for delivery arms).
usage: uv run python experiments/612-key-gated-hook/score612.py [modes ...]   (missing runs are skipped)"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L   # noqa: E402
import table521 as TB  # noqa: E402
TAG = "stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
SRC = {"base": ROOT / "experiments/435-zs-on-coarse/work/h1", "c4h": ROOT / "experiments/537-arm-pipeline/work/h1_cpm+c4h",
       "c4kg": ROOT / "experiments/612-key-gated-hook/work/h1_c4kg", "idkg": ROOT / "experiments/612-key-gated-hook/work/h1_idkg",
       "c4kgall": ROOT / "experiments/612-key-gated-hook/work/h1_c4kgall",
       "idkg3": ROOT / "experiments/612-key-gated-hook/work/h1_idkg3",
       "idall3": ROOT / "experiments/612-key-gated-hook/work/h1_idall3"}
S = L.load_state(0)
own, rl = S["own"], S["rl"]
adj = L.adjud(S)
POOLS = {reg: TB.regime_pool(S, reg)[0] for reg in TB.REGIMES}


def load(mode, s):
    p = SRC[mode] / TAG.format(s=s)
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    return dict(q=np.asarray(z["q"], float), p=np.asarray(z["p"], float), lfc=np.asarray(z["lfc"], float), sums=z["sums"])


def members(run, reg):
    pool = POOLS[reg]
    call = (run["q"] < 0.05) & ~own
    St = dict(S); St["oq"], St["op"], St["ol"] = run["q"], run["p"], run["lfc"]
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], adj[i], run["lfc"][i], rl[i]) for i in range(len(S["targets"]))])
    return dict(fid=float(np.nanmean(L.fid_rows(call, run["lfc"], pool, adj, rl))), jac=float(np.mean(L.jac_rows(call, pool))),
                reach=float(np.nanmean(rr)), nmae=float(np.nanmean(L.nmae_rows(run["lfc"], pool, rl))))


out = {}
ka = (load("c4kgall", 0), load("c4h", 0))
if ka[0] is not None:
    eq = all(np.array_equal(ka[0][k], ka[1][k], equal_nan=np.issubdtype(np.asarray(ka[0][k]).dtype, np.floating))
             for k in ("q", "p", "lfc", "sums"))
    print(f"KA c4kgall s0 == c4h s0 (q, p, lfc, sums bit-equal): {'PASS' if eq else 'FAIL'}")
    out["ka"] = bool(eq)
pairs = [("c4h", "base"), ("c4kg", "base"), ("idkg", "base"), ("c4kg", "c4h"), ("idkg", "c4kg"), ("idkg", "c4h"),
         ("idkg3", "base"), ("idkg3", "c4h"), ("idkg3", "c4kg"), ("idkg3", "idkg"),   # idkg3: 604 v3 rule (no extension)
         ("idall3", "base"), ("idall3", "c4h"), ("idall3", "idkg3")]                  # idall3: v3 source, no key gate
for arm, ref in pairs:
    rows = []
    for s in (0, 1, 2):
        a, b = load(arm, s), load(ref, s)
        if a is None or b is None:
            continue
        assert np.array_equal(a["sums"], b["sums"]), f"{arm} vs {ref} s{s}: bulk sums differ (not a DE-only change)"
        rows.append({reg: (members(a, reg), members(b, reg)) for reg in TB.REGIMES})
    if not rows:
        continue
    out[f"{arm}-{ref}"] = {}
    print(f"\n{arm} - {ref} ({len(rows)} seeds; raw member deltas, overall-equivalent at 192's INV)")
    for reg in TB.REGIMES:
        d = {m: float(np.mean([r[reg][0][m] - r[reg][1][m] for r in rows])) for m in ("fid", "jac", "reach", "nmae")}
        ov = float(np.mean([L.overall_delta(r[reg][0], r[reg][1]) for r in rows]))
        sd = float(np.std([L.overall_delta(r[reg][0], r[reg][1]) for r in rows])) if len(rows) > 1 else float("nan")
        per = {m: [float(r[reg][0][m] - r[reg][1][m]) for r in rows] for m in ("nmae", "reach")}   # per-seed (added for 612 A1)
        out[f"{arm}-{ref}"][reg] = dict(d, overall=ov, overall_sd=sd, per_seed=per, seeds=len(rows))
        board_nmae = (d["nmae"] * 1.2, d["nmae"] * 1.6)
        print(f"  {reg:4s} fid {d['fid']:+.5f} jac {d['jac']:+.5f} reach {d['reach']:+.5f} nmae {d['nmae']:+.5f} | overall "
              f"{ov:+.5f} (sd {sd:.5f}) | board nmae x1.2-1.6: {board_nmae[0]:+.5f} .. {board_nmae[1]:+.5f}")
(ROOT / "experiments/612-key-gated-hook/work/score612.json").write_text(json.dumps(out, indent=1))

#!/usr/bin/env python3
"""625d — price the real-emitter H1 octile flips (block prior + hook), paired against `none` per seed.

Members by 521's exact functions at its four regimes (f 1 / 0.5 / 0.35 / cmp; table521.regime_pool), each run with its own
realised order; overall-equivalent with 192's board-fitted INV (lib521.overall_delta) — score612's code path verbatim.
KA: `none` s0 (and s1) == 612's c4kg s0 (s1) bit for bit (q, p, lfc, sums).
Decode validation on H1 (what the readout's logic assumes): for octile k, CLEAN4's sign agreement with H1's TRUTH on the
gate genes of that octile (|l|-weighted, 521's pools). Agreement < 0.5 predicts flipping HELPS (overall delta > 0).
usage: uv run python experiments/625-module-readout/score625d.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import table521 as TB  # noqa: E402

TAG = "stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
SRC = {"c4kg": ROOT / "experiments/612-key-gated-hook/work/h1_c4kg", "none": HERE / "octnone/work/h1"}
for k in range(8):
    SRC[f"oct{k}"] = HERE / f"oct{k}/work/h1"
S = L.load_state(0)
own, rl = S["own"], S["rl"]
adj = L.adjud(S)
POOLS = {reg: TB.regime_pool(S, reg)[0] for reg in TB.REGIMES}
V8 = np.load(HERE / "work/modules625_v8.npz", allow_pickle=True)
OCTMAP = dict(zip(V8["genes"].astype(str), V8["octile"].astype(int)))
GOCT = np.array([OCTMAP.get(str(g), -1) for g in S["genes"]])
CP = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
c4pos = {g: i for i, g in enumerate(CP["genes"].astype(str))}
C4SIGN = np.array([(-1.0 if CP["P_down"][c4pos[g]] > 0.5 else 1.0) if g in c4pos else np.nan for g in S["genes"]])


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


def truth_agreement(k, reg="1.0"):
    pool = POOLS[reg]
    m = pool & (GOCT[None, :] == k) & np.isfinite(C4SIGN)[None, :] & np.isfinite(rl) & (rl != 0) & ~own
    agree = (np.sign(rl) == C4SIGN[None, :])
    w = np.abs(rl)
    return float((agree * w * m).sum() / max((w * m).sum(), 1e-12)), int(m.sum())


def main() -> None:
    out = {"ka": {}}
    for s in (0, 1):
        a, b = load("none", s), load("c4kg", s)
        if a is None or b is None:
            continue
        eq = all(np.array_equal(a[k], b[k], equal_nan=np.issubdtype(np.asarray(a[k]).dtype, np.floating)) for k in ("q", "p", "lfc", "sums"))
        out["ka"][f"s{s}"] = bool(eq)
        print(f"KA none s{s} == 612 c4kg s{s} (q, p, lfc, sums bit-equal): {'PASS' if eq else 'FAIL'}")
    out["octiles"] = {}
    for k in range(8):
        rows = []
        for s in (0, 1, 2):
            a, b = load(f"oct{k}", s), load("none", s)
            if a is None or b is None:
                continue
            assert np.array_equal(a["sums"], b["sums"]), f"oct{k} s{s}: bulk sums differ (not a prior-only change)"
            rows.append({reg: (members(a, reg), members(b, reg)) for reg in TB.REGIMES})
        if not rows:
            continue
        ag, n = truth_agreement(k)
        rec = {"seeds": len(rows), "clean4_truth_agreement_f1.0": ag, "n_gate_genes_f1.0": n, "regimes": {}}
        print(f"\noct{k} - none ({len(rows)} seeds); CLEAN4 vs H1 truth agreement on this octile's gate genes {ag:.3f} (n {n}) "
              f"-> flip predicted to {'HELP' if ag < 0.5 else 'HURT'}")
        for reg in TB.REGIMES:
            d = {m: float(np.mean([r[reg][0][m] - r[reg][1][m] for r in rows])) for m in ("fid", "jac", "reach", "nmae")}
            ovs = [L.overall_delta(r[reg][0], r[reg][1]) for r in rows]
            # reach on H1 is driven by a few tiny-n_real targets' discrete jumps (oct2 s0: 5 targets with n_real 6-25 carry
            # it), which do not transfer across panels; the transfer rule reads nmae + fid, so report overall without reach too
            ovn = [L.overall_delta({**r[reg][0], "reach": r[reg][1]["reach"]}, r[reg][1]) for r in rows]
            per_seed = {m: [float(r[reg][0][m] - r[reg][1][m]) for r in rows] for m in ("fid", "jac", "reach", "nmae")}
            rec["regimes"][reg] = dict(d, overall=float(np.mean(ovs)), overall_sd=float(np.std(ovs)) if len(ovs) > 1 else float("nan"),
                                       overall_noreach=float(np.mean(ovn)),
                                       overall_noreach_sd=float(np.std(ovn)) if len(ovn) > 1 else float("nan"), per_seed=per_seed)
            print(f"  {reg:4s} fid {d['fid']:+.5f} jac {d['jac']:+.5f} reach {d['reach']:+.5f} nmae {d['nmae']:+.5f} | overall "
                  f"{np.mean(ovs):+.5f} (sd {rec['regimes'][reg]['overall_sd']:.5f}) | w/o reach {np.mean(ovn):+.5f} "
                  f"(sd {rec['regimes'][reg]['overall_noreach_sd']:.5f}) | nmae per seed {' '.join(f'{v:+.5f}' for v in per_seed['nmae'])}")
        out["octiles"][str(k)] = rec
    (HERE / "work/score625d.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""625e: how much of a sign SOURCE's change reaches the DELIVERED lfc through the real production emitter (H1, 625d runs).

An octile flip changes the source (CLEAN4's P) on every covered gene of the octile. The emitter applies the source only
where it acts: the CPM block (the block prior) and the uncalled weak-key half with |P - 1/2| >= tau (the 612 hook). The
head's signs follow the key, and the delivered magnitudes on those routes are small. Per octile and seed, on the gate
genes of 521's pools (f 1.0 and 0.35):
  source change  = -(right - wrong) of CLEAN4 on the octile's gate genes / all gate genes. This is the change in gate-sign
                   accuracy if every gate gene's delivered sign followed the flipped source.
  delivered      = the change in delivered-lfc sign accuracy over all gate genes (none -> oct<k>, same seed).
  delivery ratio = delivered / source.
  where          = delivered sign flips on the octile's gate genes, split by called (q < 0.05) and uncalled.
  magnitude      = mean |delivered lfc| on the flipped gate genes vs mean |real lfc| there.
  exchange       = overall-equivalent without reach (fid + jac + nmae, 192's INV) per unit of SOURCE and DELIVERED change.
Reach is left out, because on H1 it moves through a few tiny-n_real targets' discrete jumps (625d note).
usage: uv run python experiments/625-module-readout/delivery625.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.argv = sys.argv[:1]
_spec = importlib.util.spec_from_file_location("score625d", HERE / "score625d.py")
s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s)
L, S, rl, own = s.L, s.S, s.rl, s.own


def oracle_mask(k):
    """625f's flip set for oracle<k>, on the lib521 kept axis (oracle625.py's definition, recomputed here)."""
    pool = s.TB.regime_pool(S, "1.0")[0] & ~own
    npool = pool.sum(0)
    pdn = np.where(npool >= 3, ((rl < 0) & pool).sum(0) / np.maximum(npool, 1), np.nan)
    c4 = s.C4SIGN                                   # -1 = CLEAN4 says DOWN
    return (s.GOCT == k) & np.isfinite(c4) & np.isfinite(pdn) & ((pdn > 0.5) != (c4 < 0))


def rows_for(label, k, a, b, flipmask, seed, out):
    for reg in ("1.0", "0.35"):
        pool = s.POOLS[reg]
        gate = pool & np.isfinite(rl) & (rl != 0) & ~own
        m = gate & flipmask[None, :]
        okb = (np.sign(b["lfc"]) == np.sign(rl)) & (b["lfc"] != 0)
        oka = (np.sign(a["lfc"]) == np.sign(rl)) & (a["lfc"] != 0)
        c4ok = np.sign(rl) == s.C4SIGN[None, :]
        src = float(-((2.0 * c4ok - 1.0) * m).sum() / gate.sum())
        dlv = float(((oka & gate).sum() - (okb & gate).sum()) / gate.sum())
        ma, mb = s.members(a, reg), s.members(b, reg)
        parts = {mm: float((ma[mm] - mb[mm]) * L.INV[mm] / 6.0 * (-1 if mm == "nmae" else 1)) for mm in ("fid", "jac", "nmae")}
        dcalls = float((((a["q"] < 0.05) & ~own).sum(1) - ((b["q"] < 0.05) & ~own).sum(1)).mean())
        out.append(dict(mode=label, octile=k, seed=seed, regime=reg, source_change=src, delivered_change=dlv,
                        delivery_ratio=dlv / src if src else float("nan"), d_calls_per_target=dcalls, **{f"part_{x}": v for x, v in parts.items()},
                        nmae_only_per_001=parts["nmae"] / src * 0.01 if src else float("nan"),
                        nmae_fid_per_001=(parts["nmae"] + parts["fid"]) / src * 0.01 if src else float("nan")))


def main() -> None:
    out = []
    for k in range(8):
        for seed in (0, 1, 2):
            a, b = s.load(f"oct{k}", seed), s.load("none", seed)
            if a is None or b is None:
                continue
            for reg in ("1.0", "0.35"):
                pool = s.POOLS[reg]
                gate = pool & np.isfinite(rl) & (rl != 0) & ~own
                m = gate & (s.GOCT[None, :] == k)
                okb = (np.sign(b["lfc"]) == np.sign(rl)) & (b["lfc"] != 0)
                oka = (np.sign(a["lfc"]) == np.sign(rl)) & (a["lfc"] != 0)
                c4ok = np.sign(rl) == s.C4SIGN[None, :]
                src = float(-((2.0 * c4ok - 1.0) * m).sum() / gate.sum())
                dlv = float(((oka & gate).sum() - (okb & gate).sum()) / gate.sum())
                flip = m & (np.sign(a["lfc"]) != np.sign(b["lfc"])) & (a["lfc"] != 0) & (b["lfc"] != 0)
                called = b["q"] < 0.05
                ma, mb = s.members(a, reg), s.members(b, reg)
                ov_noreach = L.overall_delta({**ma, "reach": mb["reach"]}, mb)
                out.append(dict(octile=k, seed=seed, regime=reg, gate_pairs=int(m.sum()), share_of_gate=float(m.sum() / gate.sum()),
                                source_change=src, delivered_change=dlv, delivery_ratio=dlv / src if src else float("nan"),
                                flips=int(flip.sum()), flips_called=int((flip & called).sum()), flips_uncalled=int((flip & ~called).sum()),
                                mean_abs_delivered_on_flips=float(np.mean(np.abs(a["lfc"][flip]))) if flip.any() else float("nan"),
                                mean_abs_real_on_flips=float(np.mean(np.abs(rl[flip]))) if flip.any() else float("nan"),
                                d_nmae=float(ma["nmae"] - mb["nmae"]), d_fid=float(ma["fid"] - mb["fid"]), overall_noreach=float(ov_noreach),
                                per_001_source=float(ov_noreach / src * 0.01) if src else float("nan"),
                                per_001_delivered=float(ov_noreach / dlv * 0.01) if dlv else float("nan")))
    print("oct seed reg  | gate share | source dacc | delivered dacc | ratio | flips called/uncalled | |lfc| dlv/real | dnmae    dfid     | overall w/o reach | per +0.01 source / delivered")
    for r in out:
        print(f" {r['octile']}   {r['seed']}  {r['regime']:4s} | {r['share_of_gate']:.3f}      | {r['source_change']:+.5f}    | {r['delivered_change']:+.5f}       | "
              f"{r['delivery_ratio']:.2f}  | {r['flips_called']:5d}/{r['flips_uncalled']:<5d}         | {r['mean_abs_delivered_on_flips']:.3f}/{r['mean_abs_real_on_flips']:.3f}    | "
              f"{r['d_nmae']:+.5f} {r['d_fid']:+.5f} | {r['overall_noreach']:+.5f}          | {r['per_001_source']:+.5f} / {r['per_001_delivered']:+.5f}")
    for reg in ("1.0", "0.35"):
        rr = [r for r in out if r["regime"] == reg]
        if not rr:
            continue
        src = np.array([r["source_change"] for r in rr]); ov = np.array([r["overall_noreach"] for r in rr])
        dl = np.array([r["delivered_change"] for r in rr])
        slope_src = float((src @ ov) / (src @ src))           # through the origin: overall per unit source change
        slope_dlv = float((dl @ ov) / (dl @ dl))
        print(f"f {reg}: pooled exchange, overall w/o reach per +0.01 gate-sign: SOURCE {slope_src * 0.01:+.5f}, DELIVERED {slope_dlv * 0.01:+.5f}; "
              f"mean delivery ratio {np.mean([r['delivery_ratio'] for r in rr]):.2f} (n {len(rr)})")
    # 625f: oracle-direction arms (IMPROVE signs on the octile genes where H1's own programme disagrees with CLEAN4)
    orc = []
    for k in range(8):
        s.SRC[f"oracle{k}"] = HERE / f"oracle{k}/work/h1"
        a, b = s.load(f"oracle{k}", 0), s.load("none", 0)
        if a is None or b is None:
            continue
        assert np.array_equal(a["sums"], b["sums"]), f"oracle{k}: bulk sums differ"
        fm = oracle_mask(k)
        rows_for(f"oracle{k}", k, a, b, fm, 0, orc)
        # the same octile's whole-octile flip, decomposed identically, for the side-by-side
        w = s.load(f"oct{k}", 0)
        if w is not None:
            rows_for(f"oct{k}", k, w, b, s.GOCT == k, 0, orc)
    if orc:
        print("\n625f: oracle (improve) vs whole-octile (degrade) flips, seed 0")
        print("mode     reg  | source dacc | delivered | ratio | dcalls | parts nmae / fid / jac | nmae-only per +0.01 | nmae+fid per +0.01")
        for r in orc:
            print(f"{r['mode']:8s} {r['regime']:4s} | {r['source_change']:+.5f}   | {r['delivered_change']:+.5f} | {r['delivery_ratio']:.2f}  | {r['d_calls_per_target']:+6.0f} | "
                  f"{r['part_nmae']:+.5f} / {r['part_fid']:+.5f} / {r['part_jac']:+.5f} | {r['nmae_only_per_001']:+.5f} | {r['nmae_fid_per_001']:+.5f}")
    (HERE / "work/delivery625.json").write_text(json.dumps({"octile_flips": out, "oracle_625f": orc}, indent=1))


if __name__ == "__main__":
    main()

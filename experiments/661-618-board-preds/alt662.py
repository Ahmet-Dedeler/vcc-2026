#!/usr/bin/env python3
"""661 + 662: the I-662 ALTERNATIVE for 661's registered 618 predictions (630 alone, DN on 630, idall3 − c4h, 625's octile
flips). It re-scores the same 618 runs on the h1-trainval state at 662's level f* = 0.61 (band 0.58–0.65 shown) beside
the f 1 / 0.5 / 0.35 / cmp values 661 used.
Why: K20's board read (10-04) put the board's level for call-set changes at f ≈ 0.61, not 619's f 0.35 (662 RESULTS). The
delivery pairs cannot tell f 0.35 from f 0.61 for nmae, so 661's f 0.35 lines STAND as registered. This adds a second,
labelled prediction before each read, so that the level-dependent levers (DN above all) can discriminate the two levels.
Known answers (asserted, 1e-12): 630 alone and (b) 630 + DN == 618/work/stacks618b.json; idall3 − c4h == 618/work/score618.json;
octile k vs none == work/score661oct.json; all at f 1 / 0.5 / 0.35 / cmp.
-> work/alt662.json.   PYTHONDONTWRITEBYTECODE=1 uv run python experiments/661-618-board-preds/alt662.py"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W618 = ROOT / "experiments/618-h1-trainval-surface/work"
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
import lib521 as L  # noqa: E402

L.REF = W618 / "ref_tv/reference.npz"
L.EFF = W618 / "tv/effects.npz"
L.ST435 = W618 / "g2mde/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
L.load_state.cache_clear()
import price528 as P5  # noqa: E402
from policy521 import means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REGS = ("1.0", "0.5", "0.35", "cmp", "0.58", "0.6", "0.61", "0.62", "0.63", "0.65", "0.7")
MEM = ("fid", "jac", "reach", "nmae")
SEEDS = (0, 1, 2)
TZ = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
RUNS = {"stk_base": W618 / "stk_base" / TZ, "630u": W618 / "stacks/stats-rnup_t618u_s{s}_b_k100_q20r3p10.npz",
        "630b_DN": W618 / "stacks/stats-rnup_t618b_s{s}_b_k100_q20r3p10_D200.npz",
        "g1": W618 / "g1" / TZ, "c4h": W618 / "hk_c4h" / TZ, "idall3": W618 / "hk_id" / TZ,
        **{f"oct{k}": HERE / f"work/tv_oct{k}" / TZ for k in ["none"] + list(range(8))}}
PAIRS = {"630 alone": ("630u", "stk_base"), "(b) 630 + DN": ("630b_DN", "stk_base"), "DN on 630": ("630b_DN", "630u"),
         "idall3 - c4h": ("idall3", "c4h"), "c4h - base": ("c4h", "g1"),
         **{f"octile {k}": (f"oct{k}", "octnone") for k in range(8)}}


def load(mode, s):
    z = np.load(str(RUNS[mode]).format(s=s), allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float))


def main():
    states = {s: L.load_state(s) for s in SEEDS}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in SEEDS for r in REGS}
    means_c = {}
    for mode in RUNS:
        for s in SEEDS:
            A = load(mode, s)
            for r in REGS:
                means_c[(mode, s, r)] = means(P5.rows_at(states[s], A, pools[(s, r)]))
    out = {}
    for name, (a, b) in PAIRS.items():
        out[name] = {}
        for r in REGS:
            d = np.stack([means_c[(a, s, r)] - means_c[(b, s, r)] for s in SEEDS])
            out[name][r] = dict(d=dict(zip(MEM, map(float, d.mean(0)))), sd=dict(zip(MEM, map(float, d.std(0, ddof=1)))))
    # known answers
    st = json.load(open(W618 / "stacks618b.json"))
    sc = json.load(open(W618 / "score618.json"))
    oc = json.load(open(HERE / "work/score661oct.json"))
    err = 0.0
    for r in ("1.0", "0.5", "0.35", "cmp"):
        err = max(err, max(abs(st["630 alone (x1)"][r]["d"][m] - out["630 alone"][r]["d"][m]) for m in MEM),
                  max(abs(st["(b) 630 + DN"][r]["d"][m] - out["(b) 630 + DN"][r]["d"][m]) for m in MEM),
                  max(abs(sc["idall3-c4h"][r][m] - out["idall3 - c4h"][r]["d"][m]) for m in MEM),
                  max(abs(oc[str(k)][r][m] - out[f"octile {k}"][r]["d"][m]) for k in range(8) for m in MEM))
    assert err < 1e-12, err
    out["KA_max_err"] = err
    (HERE / "work/alt662.json").write_text(json.dumps(out, indent=1))
    for name in PAIRS:
        print(name)
        for r in REGS:
            d = out[name][r]["d"]
            print(f"  f {r:4s} " + " ".join(f"{m} {d[m]:+.5f}" for m in MEM) + f"  (nmae sd {out[name][r]['sd']['nmae']:.5f})")
    print(f"KA (stacks618b, score618, score661oct at the four stored levels): max |err| {err:.1e}")


if __name__ == "__main__":
    main()

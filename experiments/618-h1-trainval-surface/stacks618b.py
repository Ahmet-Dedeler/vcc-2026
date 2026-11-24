#!/usr/bin/env python3
"""618 stacks, second pass (10-03): price every arm vs plain with price630's OWN functions (528's rows_at, 521's means /
SGN, the same pools and seeds) but WITHOUT its "per-gene sums pinned" assert.
Why: x1.25 is NOT bulk-neutral in this harness. Its per-gene block totals differ from plain in ~1.35M cells over all 126
targets (|d| up to ~5e4 counts), so price630 refused it (rc 1 for t618a and t618x). 630 and DN preserve sums, so (b) and
630-alone priced normally.
- KA: (b) and 630-alone must reproduce price618's value_mean exactly at every regime.
- Reported per arm and regime: the registered overall-equivalent (SGN over fid / jac / reach / nmae, as price630) AND the
  nmae-only overall-equivalent (619's amendment: reach is unpriced for delivery levers).
- The decision is mechanical on the REGISTERED statistic (PRED-618-stacks.md / the coordinator's rule). The nmae-only
  version is shown beside it.
-> work/stacks618b.json"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work"
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
import lib521 as L  # noqa: E402

L.REF = W / "ref_tv/reference.npz"
L.EFF = W / "tv/effects.npz"
L.ST435 = W / "g2mde/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
L.load_state.cache_clear()
import price528 as P5  # noqa: E402
from policy521 import SGN, means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REG = ("1.0", "0.5", "0.35", "cmp")
TB = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
TX = "stats-x125_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
ARMS = {"(a) x1.25 + 630": W / "stacks/stats-rnup_t618a_s{s}_b_k100_q20r3p10.npz",
        "(b) 630 + DN": W / "stacks/stats-rnup_t618b_s{s}_b_k100_q20r3p10_D200.npz",
        "630 alone (x1)": W / "stacks/stats-rnup_t618u_s{s}_b_k100_q20r3p10.npz",
        "x1.25 alone": W / "stk_x125" / TX}
KA = {"(b) 630 + DN": ("t618b", "b_k100_q20r3p10_D200"), "630 alone (x1)": ("t618u", "b_k100_q20r3p10")}
SG_NMAE = np.array([0.0, 0.0, 0.0, SGN[3]])   # nmae-only overall-equivalent (fid, jac, reach zeroed)


def load(path):
    z = np.load(path, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float))


def main():
    seeds = (0, 1, 2)
    states = {s: L.load_state(s) for s in seeds}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in seeds for r in REG}
    base = {s: load(W / "stk_base" / TB.format(s=s)) for s in seeds}
    Rb = {(s, r): P5.rows_at(states[s], base[s], pools[(s, r)]) for s in seeds for r in REG}
    out, lines = {}, []
    for arm, path in ARMS.items():
        rec = {}
        for r in REG:
            dms = []
            for s in seeds:
                A = load(str(path).format(s=s))
                dms.append(means(P5.rows_at(states[s], A, pools[(s, r)])) - means(Rb[(s, r)]))
            dm = np.mean(np.stack(dms), 0)
            rec[r] = dict(value=float(np.mean([SGN @ d for d in dms])), value_nmae_only=float(np.mean([SG_NMAE @ d for d in dms])),
                          d=dict(zip(("fid", "jac", "reach", "nmae"), map(float, dm))),
                          per_seed=[float(SGN @ d) for d in dms])
        out[arm] = rec
        lines.append(f"{arm:18s} " + "  ".join(f"f {r} {rec[r]['value']:+.5f} (nmae-only {rec[r]['value_nmae_only']:+.5f})" for r in REG))
    # KA vs price618 for the sums-preserving arms
    ka = {}
    for arm, (name, v) in KA.items():
        pj = json.loads((W / "stacks" / f"price630_{name}.json").read_text())[v]
        ka[arm] = max(abs(pj[r]["value_mean"] - out[arm][r]["value"]) for r in REG)
    ka_ok = all(x < 1e-12 for x in ka.values())
    lines.append(f"KA vs price618 (max |diff| over regimes): {ka} -> {'PASS' if ka_ok else 'FAIL'}")
    mins = {arm: min(out[arm]["0.5"]["value"], out[arm]["0.35"]["value"]) for arm in out}
    mins_n = {arm: min(out[arm]["0.5"]["value_nmae_only"], out[arm]["0.35"]["value_nmae_only"]) for arm in out}
    d = mins["(a) x1.25 + 630"] - mins["(b) 630 + DN"]
    dn = mins_n["(a) x1.25 + 630"] - mins_n["(b) 630 + DN"]
    pick = "(a) x1.25 + 630" if d > 0.0005 else "(b) 630 + DN"
    pick_n = "(a) x1.25 + 630" if dn > 0.0005 else "(b) 630 + DN"
    h1 = "TRANSFERS" if 0.0004 <= d <= 0.0024 else "does NOT transfer"
    lines.append(f"REGISTERED statistic: MIN(a) {mins['(a) x1.25 + 630']:+.5f} - MIN(b) {mins['(b) 630 + DN']:+.5f} = {d:+.5f} -> PICK {pick} "
                 f"(H1 predicted +0.0014: H1's ranking {h1})")
    lines.append(f"nmae-only (619): MIN(a) {mins_n['(a) x1.25 + 630']:+.5f} - MIN(b) {mins_n['(b) 630 + DN']:+.5f} = {dn:+.5f} -> {pick_n}")
    out["decision"] = dict(registered=dict(min_a=mins["(a) x1.25 + 630"], min_b=mins["(b) 630 + DN"], diff=d, pick=pick, h1_ranking=h1),
                           nmae_only=dict(min_a=mins_n["(a) x1.25 + 630"], min_b=mins_n["(b) 630 + DN"], diff=dn, pick=pick_n),
                           ka_vs_price618=ka, ka_pass=ka_ok,
                           caveat="x1.25 changes per-gene block totals (bulk) in this harness; its pds / mse effect is not in "
                                  "these DE-member values. (b) is bulk-neutral (sums pinned).")
    (W / "stacks618b.json").write_text(json.dumps(out, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

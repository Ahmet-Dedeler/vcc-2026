#!/usr/bin/env python3
"""537 — price the real-emitter H1 runs (run537.py; plus 528's / 435's recorded runs where a seed exists there) on the four DE
members at 521's regimes, with target bootstraps, exactly as 528's price528.py does (its functions are imported read-only).

A comparison is (arm, base): per seed, the arm's and the base's stats (zs form = 446bj2z's DE construction) are scored with
521's exact member functions (realized order) on the truth's pools at f 1 / 0.5 / 0.35 (PowerSurface) and cmp; the value is
192's board-fitted overall (policy521.SGN); board columns are 523's convention (fid, jac at 0.72x; reach x1 | x0 | x-0.76).
CIs: 4,000 target resamples of the seed-averaged rows. pds/mse: the per-gene sums must be equal (a DE-only change).
Stats are looked up in this directory's work/h1_<mode>/ first, then 528's, then 435's (mode `base`).
Known answers: `cpm` seed 0 == 435 (ka537.py); price528's own KA (435's members from stats == its metrics json) runs first.
Writes work/price537.json and prints one line per (comparison, regime).
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
import price528 as P5  # noqa: E402  (read-only: rows_at, boot, boot_member, SGN, board_cols, means ...)
from global521 import board_cols  # noqa: E402
from policy521 import SGN, means, per_target_contrib, top3_share  # noqa: E402
from table521 import regime_pool  # noqa: E402
from vcc2026.rc_emit_k import block_from_key  # noqa: E402

TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
DIRS = [HERE / "work", ROOT / "experiments/528-callprop-block/work"]
ST435 = ROOT / "experiments/435-zs-on-coarse/work/h1"
REG = ("1.0", "0.5", "0.35", "cmp")
MEMBERS = ("fid", "jac", "reach", "nmae")
COMPARISONS = [
    # (label, arm, base)
    ("REF cp (528, kept) on 446bj2z", "lcap-kept", "base"),
    ("P1 CLEAN4 hook on 446bj2z (= 496a's knob on the zs base)", "cpm+c4h", "base"),
    ("P1 cp on the CLEAN4-hook base", "lcap-kept+c4h", "cpm+c4h"),
    ("P1 CLEAN4 hook on the cp base", "lcap-kept+c4h", "lcap-kept"),
    ("P1 cp + CLEAN4 hook, total vs 446bj2z", "lcap-kept+c4h", "base"),
    ("P2 K 2,000 alone (446bj2zk20's construction)", "cpm-n1600", "base"),
    ("P2 cp kept at K 2,000", "lcap-kept-n1600", "cpm-n1600"),
    ("P2 cp coupled at K 2,000", "lcap-n1600", "cpm-n1600"),
    ("P2 cp kept + K 2,000, total vs 446bj2z", "lcap-kept-n1600", "base"),
    ("P3 K 1,600 alone (446bj2zk16's construction)", "cpm-n1200", "base"),
    ("P3 cp kept at K 1,600", "lcap-kept-n1200", "cpm-n1200"),
    ("P3 cp kept + K 1,600, total vs 446bj2z", "lcap-kept-n1200", "base"),
    ("P2/3 K 2,000 vs K 1,600 (523's dose pair, K20 - k16)", "cpm-n1600", "cpm-n1200"),
    ("REF Arc hook r 0 (446bj2zmz's construction)", "cpm+mzg", "base"),
    ("REF cp on the Arc-hook base r 0 (528)", "lcap-kept+mzg", "cpm+mzg"),
    ("P4 Arc hook at r 0.37 (the 458d read)", "cpm+mzg-r0.37", "base"),
    ("P4 cp on the Arc-hook base at r 0.37", "lcap-kept+mzg-r0.37", "cpm+mzg-r0.37"),
    ("P4 cp + Arc hook at r 0.37, total vs 446bj2z", "lcap-kept+mzg-r0.37", "base"),
    ("P4 Arc hook r 0.37 vs CLEAN4 hook (source effect at fixed delivery)", "cpm+mzg-r0.37", "cpm+c4h"),
    ("P4 cp+Arc hook r 0.37 vs cp+CLEAN4 hook (the source at fixed block)", "lcap-kept+mzg-r0.37", "lcap-kept+c4h"),
    ("P5 CLEAN4 hook at K 2,000", "cpm-n1600+c4h", "cpm-n1600"),
    ("P5 cp on the hook base at K 2,000", "lcap-kept-n1600+c4h", "cpm-n1600+c4h"),
    ("P5 cp + hook + K 2,000, total vs 446bj2z", "lcap-kept-n1600+c4h", "base"),
    ("P6 CLEAN4 hook at K 1,600", "cpm-n1200+c4h", "cpm-n1200"),
    ("P6 cp on the hook base at K 1,600", "lcap-kept-n1200+c4h", "cpm-n1200+c4h"),
    ("P6 cp + hook + K 1,600, total vs 446bj2z", "lcap-kept-n1200+c4h", "base"),
    ("P7 CLEAN4 hook tau 0.03 on 446bj2z", "cpm+c4h0.03", "base"),
    ("P7 CLEAN4 hook tau 0.03 on the cp base", "lcap-kept+c4h0.03", "lcap-kept"),
    ("P7 cp + hook tau 0.03 vs cp + hook tau 0.05 (the hook's dose on the block)", "lcap-kept+c4h0.03", "lcap-kept+c4h"),
    ("P7 cp + hook tau 0.03, total vs 446bj2z", "lcap-kept+c4h0.03", "base"),
    ("P8 CLEAN4 hook attenuated r 0.17 on 446bj2z", "cpm+c4h-r0.17", "base"),
    ("P8 CLEAN4 hook r 0.17 on the cp base", "lcap-kept+c4h-r0.17", "lcap-kept"),
    ("P8 cp + hook r 0.17, total vs 446bj2z", "lcap-kept+c4h-r0.17", "base"),
    ("P9 sign-shuffled CLEAN4 hook on 446bj2z (zero-information floor)", "cpm+c4h-s", "base"),
    ("P9 sign-shuffled CLEAN4 hook on the cp base", "lcap-kept+c4h-s", "lcap-kept"),
    ("P9 cp + sign-shuffled hook, total vs 446bj2z", "lcap-kept+c4h-s", "base"),
    ("P9 CLEAN4 signs vs shuffled signs on the cp base (the sign information)", "lcap-kept+c4h", "lcap-kept+c4h-s"),
    ("P10 sign-shuffled hook on the cp base at K 2,000", "lcap-kept-n1600+c4h-s", "lcap-kept-n1600"),
    ("P10 cp + sign-shuffled hook + K 2,000, total vs 446bj2z", "lcap-kept-n1600+c4h-s", "base"),
    ("P10 CLEAN4 vs shuffled signs on the cp base at K 2,000", "lcap-kept-n1600+c4h", "lcap-kept-n1600+c4h-s"),
    ("P10 CPHK20 vs CPH (block + hook, K 2,000 vs 2,600)", "lcap-kept-n1600+c4h", "lcap-kept+c4h"),
    ("P11 block at swap 0.5 (kept), vs 446bj2z", "lcap-k1100-kept", "base"),
    ("P11 CLEAN4 hook on the swap-0.5 block", "lcap-k1100-kept+c4h", "lcap-k1100-kept"),
    ("P11 swap-0.5 block + hook, total vs 446bj2z", "lcap-k1100-kept+c4h", "base"),
    ("P11 swap-0.5 block + hook vs CPH (the dose partner vs the full swap)", "lcap-k1100-kept+c4h", "lcap-kept+c4h"),
    ("P11 CLEAN4 hook tau 0.10 on 446bj2z", "cpm+c4h0.10", "base"),
    ("P11 CLEAN4 hook tau 0.10 on the cp base", "lcap-kept+c4h0.10", "lcap-kept"),
    ("P11 cp + hook tau 0.10 vs cp + hook tau 0.05 (the hook's dose partner on the block)", "lcap-kept+c4h0.10", "lcap-kept+c4h"),
    ("P12 sign-shuffled Arc hook on 446bj2z (MZ's zero-information floor)", "cpm+mzg-s", "base"),
    ("P12 sign-shuffled Arc hook on the cp base (CPMZ's floor increment)", "lcap-kept+mzg-s", "lcap-kept"),
    ("P12 cp + sign-shuffled Arc hook, total vs 446bj2z", "lcap-kept+mzg-s", "base"),
    ("P12 shuffled Arc vs shuffled CLEAN4 hook, on the block (the sources' structure without signs)", "lcap-kept+mzg-s", "lcap-kept+c4h-s"),
    ("P12 Arc (r 0) vs CLEAN4 hook, on the block (CPMZ vs CPH at H1's own information)", "lcap-kept+mzg", "lcap-kept+c4h"),
    ("P12 Arc hook r 0 on the cp base", "lcap-kept+mzg", "lcap-kept"),
]


def stats_path(mode, seed, zs=True):
    t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
    if mode == "base":
        return ST435 / f"stats-{t}.npz", ST435 / f"metrics-{t}.json"
    for d in DIRS:
        p = d / f"h1_{mode}" / f"stats-{t}.npz"
        m = d / f"h1_{mode}" / f"metrics-{t}.json"
        if p.exists() and m.exists():
            return p, m
    return None, None


def load(mode, seed, zs=True):
    p, m = stats_path(mode, seed, zs)
    if p is None:
        return None
    z = np.load(p, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float),
                sums=np.asarray(z["sums"]), raw=json.load(open(m))["raw"], src=str(p.relative_to(ROOT)))


def parse(mode):
    import re
    if mode == "base":
        return "cpm", 0, 2200
    m = re.fullmatch(r"(cpm|lcap)(?:-k([0-9]+))?(?:-kept)?(?:-n([0-9]+))?(?:\+c4h|\+mzg)?(?:[0-9]\.[0-9]+)?(?:-r[0-9.]+)?(?:-s)?(?:-bs[0-9.]+)?", mode)
    return m.group(1), int(m.group(2) or 0), int(m.group(3) or 2200)


def prescribed(mode, S, keys):
    key, keep, n = parse(mode)
    T = S["head"].shape[0]
    if key == "cpm" and n == 2200:
        return S["head"] | S["block"]
    kk = S["ref_lin"] if key == "cpm" else None
    out = np.zeros_like(S["head"])
    for i in range(T):
        bk = S["ref_lin"] if kk is not None else keys[key][i]
        out[i] = S["head"][i] | block_from_key(S["head"][i], S["ref_lin"], n, bk, n if key == "cpm" else keep)
    return out


def main():
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in (0, 1, 2) for r in REG}
    keys = dict(np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz", allow_pickle=False))
    out = {"KA": {}, "comparisons": {}}
    for s in (0, 1, 2):          # price528's KA: members recomputed from 435's stats == its metrics json
        A = load("base", s, False)
        R = P5.rows_at(states[s], A, pools[(s, "1.0")])
        got = dict(zip(MEMBERS, (float(np.nanmean(R[m])) for m in range(4))))
        ok = all(abs(got[m] - A["raw"][m]) < 1e-12 for m in MEMBERS)
        out["KA"][f"base_s{s}"] = ok
        assert ok, got
    print("KA: 435's members from stats == its metrics json on seeds 0/1/2", flush=True)
    for label, arm, b0 in COMPARISONS:
        seeds = [s for s in (0, 1, 2) if load(arm, s) is not None and load(b0, s) is not None]
        if not seeds:
            print(f"{label:70s} -- no runs yet ({arm} vs {b0})", flush=True)
            continue
        rec = {"arm": arm, "base": b0, "seeds": seeds,
               "sources": {s: [load(arm, s)["src"], load(b0, s)["src"]] for s in seeds}}
        rec["pinned"] = {s: dict(sums_equal=bool(np.array_equal(load(arm, s)["sums"], load(b0, s)["sums"])),
                                 dpds=load(arm, s)["raw"]["pds"] - load(b0, s)["raw"]["pds"],
                                 dmse=load(arm, s)["raw"]["mse"] - load(b0, s)["raw"]["mse"]) for s in seeds}
        for reg in REG:
            per, RA, RB = [], [], []
            for s in seeds:
                Ra = P5.rows_at(states[s], load(arm, s), pools[(s, reg)])
                Rb = P5.rows_at(states[s], load(b0, s), pools[(s, reg)])
                dm = means(Ra) - means(Rb)
                per.append(dict(value=float(SGN @ dm), d=dict(zip(MEMBERS, dm.tolist())), board=list(board_cols(dm))))
                RA.append(Ra)
                RB.append(Rb)
            Ra, Rb = np.nanmean(np.stack(RA), 0), np.nanmean(np.stack(RB), 0)
            dm = np.array([np.mean([p["d"][m] for p in per]) for m in MEMBERS])
            lo, hi, p0 = P5.boot(Ra, Rb)
            contrib = per_target_contrib(Ra, Rb)
            rec[reg] = dict(per_seed=per, value_mean=float(np.mean([p["value"] for p in per])), ci=[lo, hi], p_le0=p0,
                            d_mean=dict(zip(MEMBERS, dm.tolist())), fid_jac_overall=float(SGN[0] * dm[0] + SGN[1] * dm[1]),
                            board=list(board_cols(dm)), top3_share=top3_share(contrib), n_gain=int((contrib > 0).sum()))
            r = rec[reg]
            print(f"{label[:70]:70s} f {reg:4s} n{len(seeds)} ov {r['value_mean']:+.5f} [{lo:+.5f}, {hi:+.5f}] "
                  f"(seeds {' '.join(f'{p['value']:+.5f}' for p in per)}) | fid {dm[0]:+.5f} jac {dm[1]:+.5f} reach {dm[2]:+.5f} "
                  f"nmae {dm[3]:+.5f} | fid+jac {r['fid_jac_overall']:+.5f} | board [{r['board'][0]:+.5f} | {r['board'][1]:+.5f} | "
                  f"{r['board'][2]:+.5f}]", flush=True)
        # deliverability, seed 0 (or the first seed): prescribed call set vs realized calls
        s = seeds[0]
        A, B = load(arm, s), load(b0, s)
        own = states[s]["own"]
        C_a, C_b = prescribed(arm, states[s], keys) & ~own, prescribed(b0, states[s], keys) & ~own
        D_a, D_b = (A["q"] < L.ALPHA) & ~own, (B["q"] < L.ALPHA) & ~own
        rec["deliverability"] = dict(
            prescribed_per_target=[float(C_a.sum(1).mean()), float(C_b.sum(1).mean())],
            realized_per_target=[float(D_a.sum(1).mean()), float(D_b.sum(1).mean())],
            unsuppressible_per_target=[float((D_a & ~C_a).sum(1).mean()), float((D_b & ~C_b).sum(1).mean())],
            missed_per_target=[float((C_a & ~D_a).sum(1).mean()), float((C_b & ~D_b).sum(1).mean())])
        dl = rec["deliverability"]
        print(f"   deliverability s{s} (arm | base): prescribed {dl['prescribed_per_target'][0]:.0f} | {dl['prescribed_per_target'][1]:.0f}; "
              f"realized {dl['realized_per_target'][0]:.0f} | {dl['realized_per_target'][1]:.0f}; unsuppressible "
              f"{dl['unsuppressible_per_target'][0]:.1f} | {dl['unsuppressible_per_target'][1]:.1f}; missed "
              f"{dl['missed_per_target'][0]:.1f} | {dl['missed_per_target'][1]:.1f}; sums equal "
              f"{all(v['sums_equal'] for v in rec['pinned'].values())}", flush=True)
        out["comparisons"][label] = rec
    (HERE / "work/price537.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

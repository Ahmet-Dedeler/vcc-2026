"""215 — sensitivity of every conclusion to the instrument's one free parameter, and the
cqawesome window in detail.

The decomposition has exactly ONE fitted constant: the direction channel's pds response,
E_DIR["pds"] = 0.6479 raw pds per unit transfer cosine, taken from a14 -> a15.  Everything
else is either a board measurement or 197's identity.  This sweeps it.
"""
import json, os, importlib
import numpy as np
import chan215 as C
from phase1215 import series, usable

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
M = C.MEMBERS
TEAMS = ["cqawesome", "Aginglab.com", "Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng",
         "Garrrbage", "Slurm Shadi", "Haystack", "Ahmet Dedeler"]
BASE = C.E_DIR["pds"]


def budget(panel, team):
    rs = [r for r in series(panel, team)]
    tot = dict(d=0.0, o=0.0, a=0.0, g=0.0, t=0.0)
    for b, a in zip(rs[:-1], rs[1:]):
        if not (usable(b) and usable(a)):
            continue
        r = C.classify({m: b["raw_" + m] for m in M}, {m: a["raw_" + m] for m in M})
        tot["d"] += r["b_direction"]; tot["o"] += r["b_ordering"]
        tot["a"] += r["b_amplitude"]; tot["g"] += r["b_geometry"]; tot["t"] += r["ovr_obs"]
    return tot


def main():
    panel = json.load(open(os.path.join(WORK, "panel215.json")))

    print("=" * 118)
    print("S1 — sensitivity of the per-team DIRECTION share to the one fitted constant")
    print(f"     E_DIR['pds'] = d(raw pds)/d(cos).  Base {BASE:.4f} from a14 -> a15 (197's +0.00507).")
    print("     0.686 is the board's own median on 13 DIRECTION-labelled strangers' arms "
          "(phase2 V1, rho median 0.944).")
    print("=" * 118)
    scales = [("x0.70", 0.70 * BASE), ("x1.00 (base)", BASE),
              ("board median 0.686", 0.686), ("x1.30", 1.30 * BASE)]
    print(f"{'team':<16}" + "".join(f"{lab:>21}" for lab, _ in scales))
    out = {}
    for t in TEAMS:
        cells = []
        for lab, v in scales:
            C.E_DIR["pds"] = v
            for m in M:
                C.RATIO_DIR[m] = C.E_DIR[m] / v
            C.RATIO_DIR_MSE_IDENTITY = C.E_DIR_MSE_IDENTITY / v
            b = budget(panel, t)
            pct = 100 * b["d"] / b["t"] if abs(b["t"]) > 1e-9 else float("nan")
            cells.append((b["d"], pct))
        out[t] = cells
        print(f"{t[:15]:<16}" + "".join(f"{d:>+12.4f} ({p:>4.0f}%)" for d, p in cells))
    C.E_DIR["pds"] = BASE
    for m in M:
        C.RATIO_DIR[m] = C.E_DIR[m] / BASE
    C.RATIO_DIR_MSE_IDENTITY = C.E_DIR_MSE_IDENTITY / BASE

    # ---------------------------------------------------------------- cqawesome in detail
    print("\n" + "=" * 118)
    print("S2 — cqawesome, the only top-of-board team whose direction climb is MEASURABLE, "
          "row by row")
    print("=" * 118)
    rs = [r for r in series(panel, "cqawesome") if usable(r)]
    print(f"{'date':<18}{'#':>3} {'model':<16}{'pds':>9}{'raw mse':>10}{'cos floor':>11}"
          f"{'reach':>9}{'nmae':>9}{'overall':>9}   note")
    prev = None
    for r in rs:
        inf = C.is_informative(r["raw_mse"])
        cf = C.cos_floor(r["raw_mse"]) if inf else float("nan")
        note = "" if inf else "vacuous (raw mse >= 1.023611): no cosine information"
        print(f"{r['date'][:16]:<18}{r['n']:>3} {str(r['model'])[:15]:<16}{r['raw_pds']:>9.4f}"
              f"{r['raw_mse']:>10.4f}{cf:>11.5f}{r['raw_reach']:>9.4f}{r['raw_nmae']:>9.4f}"
              f"{r['overall']:>9.4f}   {note}")
    inf_rows = [r for r in rs if C.is_informative(r["raw_mse"])]
    a, b = inf_rows[0], inf_rows[-1]
    print(f"\n  first informative row  {a['date'][:10]} #{a['n']} {a['model']}: "
          f"cos floor {C.cos_floor(a['raw_mse']):.5f}, pds {a['raw_pds']:.4f}")
    print(f"  latest row             {b['date'][:10]} #{b['n']} {b['model']}: "
          f"cos floor {C.cos_floor(b['raw_mse']):.5f}, pds {b['raw_pds']:.4f}")
    print(f"  MEASURED cosine gain inside the visible window: "
          f"{C.cos_floor(b['raw_mse']) - C.cos_floor(a['raw_mse']):+.5f}  "
          f"({C.cos_floor(b['raw_mse'])/C.cos_floor(a['raw_mse']):.2f}x)")
    print(f"  pds gain over the same window: {b['raw_pds'] - a['raw_pds']:+.4f} raw = "
          f"{(b['raw_pds'] - a['raw_pds'])/BASE:+.5f} in cosine units by the amplitude-free route")
    print(f"  the two routes agree to "
          f"{abs((b['raw_pds']-a['raw_pds'])/BASE / (C.cos_floor(b['raw_mse'])-C.cos_floor(a['raw_mse'])) - 1):.1%}")

    # what fraction of their FINAL cosine had already arrived before the window opened?
    print(f"\n  Before 2026-08-30 their raw mse never fell below 1.0236, so the board carries "
          f"NO cosine reading at all for the first 18 of their 27 submissions.")
    print(f"  What it does carry is pds, which is amplitude-free: "
          f"{rs[0]['raw_pds']:.4f} (#{rs[0]['n']}, {rs[0]['date'][:10]}) -> "
          f"{a['raw_pds']:.4f} (#{a['n']}, {a['date'][:10]}) -> {b['raw_pds']:.4f} "
          f"(#{b['n']}, {b['date'][:10]}).")

    json.dump({"sensitivity": {t: [[float(d), float(p)] for d, p in v] for t, v in out.items()},
               "scales": [[l, float(v)] for l, v in scales]},
              open(os.path.join(WORK, "sens215.json"), "w"), indent=1)
    print("\nwrote work/sens215.json")


if __name__ == "__main__":
    main()

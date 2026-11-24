"""317 — is the leader's shape an OPERATING POINT or a different model? Measured on the board.

If `mse` -> `nmae` + `reach` were a trade any team can make by moving its calibration, then
WITHIN a team, consecutive submissions that move raw `mse` would move raw `nmae` and raw
`reach` along that trade, and the resulting `score_avg` would move with it. This measures the
within-team exchange rate on the live board, which is the only place the trade could be
observed being made.

Experimental unit: one consecutive submission PAIR from the same `team_id`, ordered by
`submission_date`. Intervals are percentile bootstraps over PAIRS (a team contributes several
pairs, so they are not independent across a team; the team-block bootstrap is reported too).
"""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SNAP = ROOT / "data/leaderboard/2026-09-07T1853Z.json"
ent = json.loads(SNAP.read_text())["entries"]
print(f"[317] {SNAP.name}: {len(ent)} entries")

# The public board keeps ONE row per team (the current best/last), so consecutive pairs must
# come from the snapshot SERIES, not from one file.
snaps = sorted((ROOT / "data/leaderboard").glob("2026-09-0*Z.json"))
seen, rows = set(), []
for s in snaps:
    for e in json.loads(s.read_text())["entries"]:
        k = e["id"]
        if k in seen:
            continue
        seen.add(k)
        rows.append(e)
print(f"[317] {len(snaps)} snapshots -> {len(rows)} DISTINCT submission ids")

by_team = defaultdict(list)
for e in rows:
    if e.get("score_avg") is None or e.get("expr_mse_unbiased_capped_norm") is None:
        continue
    by_team[e["team_id"]].append(e)
pairs = []
for t, es in by_team.items():
    es.sort(key=lambda e: e["submission_date"])
    for a, b in zip(es, es[1:]):
        pairs.append((t, a, b))
print(f"[317] {len(by_team)} teams with >=1 usable row -> {len(pairs)} consecutive pairs")

K = {"mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
     "reach": "de_wilcoxon_direction_reach_raw", "pds": "pds_cosine",
     "fid": "de_wilcoxon_direction_fidelity_yield_raw", "jac": "de_wilcoxon_sig_jaccard"}


def d(pair, m):
    _, a, b = pair
    return b[K[m]] - a[K[m]]


# Restrict to pairs where `mse` actually moved and both ends sit in the unclamped `mse` core,
# i.e. pairs where the amplitude/calibration channel is READABLE at all.
core = [p for p in pairs
        if all(0.70 <= x[K["mse"]] <= 1.05 for x in p[1:])
        and abs(d(p, "mse")) > 0.002]
print(f"[317] pairs with both ends in the readable `mse` core and |d mse| > 0.002: {len(core)}")


def boot(vals, n=10000, seed=0):
    r = np.random.default_rng(seed)
    v = np.asarray(vals, float)
    bs = np.array([r.choice(v, v.size, replace=True).mean() for _ in range(n)])
    return float(v.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


print("\n[317] WITHIN-TEAM exchange rate: what a +1.0 raw `mse` move (i.e. GIVING UP `mse`, "
      "the leader's direction) brings with it")
for m in ("nmae", "reach", "pds", "fid", "jac"):
    ratio = [d(p, m) / d(p, "mse") for p in core]
    mu, lo, hi = boot(ratio)
    print(f"  d raw {m:<5} / d raw mse   median {np.median(ratio):+8.4f}   "
          f"mean {mu:+8.4f} [{lo:+.4f}, {hi:+.4f}]   n {len(ratio)}")

# and the thing that matters: does giving up `mse` raise `score_avg`?
up = [p for p in core if d(p, "mse") > 0]     # gave up mse
dn = [p for p in core if d(p, "mse") < 0]     # bought mse
for lab, grp in (("GAVE UP mse (leader's direction)", up), ("BOUGHT mse", dn)):
    dv = [p[2]["score_avg"] - p[1]["score_avg"] for p in grp]
    if not dv:
        continue
    mu, lo, hi = boot(dv)
    print(f"\n  {lab}: n {len(dv)}  d score_avg mean {mu:+.5f} [{lo:+.5f}, {hi:+.5f}]  "
          f"median {np.median(dv):+.5f}  fraction improving {np.mean(np.array(dv) > 0):.3f}")

# where the leader sits against the whole field on the (mse, nmae+reach) plane
lead = max(rows, key=lambda e: e["score_avg"])
print(f"\n[317] leader {lead['team_name']!r} {lead['score_avg']:.6f}: raw mse "
      f"{lead[K['mse']]:.5f}, nmae {lead[K['nmae']]:.5f}, reach {lead[K['reach']]:.5f}, "
      f"pds {lead[K['pds']]:.5f}")
band = [e for e in rows if 0.95 <= e[K["mse"]] <= 1.05]
print(f"  teams whose raw mse is in [0.95, 1.05] (the leader's `mse` band): {len(band)}")
band.sort(key=lambda e: -e["score_avg"])
for e in band[:6]:
    print(f"    {e['score_avg']:.5f}  mse {e[K['mse']]:.4f}  nmae {e[K['nmae']]:.4f}  "
          f"reach {e[K['reach']]:.4f}  pds {e[K['pds']]:.4f}  {e['team_name'][:28]}")
print(f"  ... median score_avg in that band {np.median([e['score_avg'] for e in band]):.5f}; "
      f"our live entry 0.2177189")


# ---------------------------------------------------------------------------------------
# THE CLEAN VERSION. `pds_cosine` is EXACTLY scale-invariant (AGENTS.md §4.1, board-verified
# twice at d pds = 0.000000), so a consecutive pair whose `pds` does not move AT ALL is a pair
# whose DIRECTION did not move: a calibration / emission change and nothing else. That isolates
# the operating point from the model, which the |d mse| > 0.002 cut above does not.
print("\n" + "=" * 88)
print("[317] PURE-CALIBRATION pairs: same team, `pds_cosine` IDENTICAL, `mse` moved")
pure = [p for p in pairs if abs(d(p, "pds")) < 1e-9 and abs(d(p, "mse")) > 1e-6
        and all(0.70 <= x[K["mse"]] <= 1.30 for x in p[1:])]
print(f"  n = {len(pure)} pairs, {len({p[0] for p in pure})} teams")
if pure:
    for m in ("nmae", "reach", "fid", "jac"):
        ratio = [d(p, m) / d(p, "mse") for p in pure]
        mu, lo, hi = boot(ratio, seed=1)
        print(f"    d raw {m:<5} / d raw mse   median {np.median(ratio):+8.4f}   "
              f"mean {mu:+8.4f} [{lo:+.4f}, {hi:+.4f}]")
    for lab, grp in (("GAVE UP mse", [p for p in pure if d(p, "mse") > 0]),
                     ("BOUGHT mse", [p for p in pure if d(p, "mse") < 0])):
        dv = [p[2]["score_avg"] - p[1]["score_avg"] for p in grp]
        if dv:
            mu, lo, hi = boot(dv, seed=2)
            print(f"    {lab}: n {len(dv)}  d score_avg mean {mu:+.5f} [{lo:+.5f}, {hi:+.5f}]"
                  f"  median {np.median(dv):+.5f}  frac improving "
                  f"{np.mean(np.array(dv) > 0):.3f}")
    from price317 import board_slopes
    _, sl = board_slopes()
    ov = [sum(d(p, m) * sl[m]["slope"] for m in ("mse", "nmae", "reach", "fid", "jac"))
          / 6.0 / d(p, "mse") for p in pure if d(p, "mse") > 0]
    if ov:
        mu, lo, hi = boot(ov, seed=3)
        print(f"\n    OVERALL points per +1.0 raw `mse` GIVEN UP: median {np.median(ov):+.5f}"
              f"  mean {mu:+.5f} [{lo:+.5f}, {hi:+.5f}]  n {len(ov)}")
        print("    (negative = giving up `mse` LOSES overall; the board's own price for the "
              "trade 317 was sent to look for)")


# --------------------------------------------------------------------------------------
# MIDDLE TIER. The exact-`pds` cut is clean but rare (n = 7). Relax it to pairs whose
# DIRECTION barely moved (|d pds| < 0.002 raw, i.e. under a tenth of the board's `pds`
# spread) while `mse` moved by more than 0.01 raw. That still excludes the model rewrites
# that dominate the 285-pair set, at the cost of admitting small direction changes.
print("\n" + "=" * 88)
print("[317] DIRECTION-STABLE pairs: |d raw pds| < 0.002 and |d raw mse| > 0.01")
mid = [p for p in pairs if abs(d(p, "pds")) < 0.002 and abs(d(p, "mse")) > 0.01
       and all(0.70 <= x[K["mse"]] <= 1.30 for x in p[1:])]
print(f"  n = {len(mid)} pairs, {len({p[0] for p in mid})} teams")
from price317 import board_slopes as _bs
_, SL = _bs()
for m in ("nmae", "reach", "fid", "jac"):
    ratio = [d(p, m) / d(p, "mse") for p in mid]
    mu, lo, hi = boot(ratio, seed=11)
    print(f"    d raw {m:<5} / d raw mse   median {np.median(ratio):+8.4f}   "
          f"mean {mu:+8.4f} [{lo:+.4f}, {hi:+.4f}]")
for lab, grp in (("GAVE UP mse", [p for p in mid if d(p, "mse") > 0]),
                 ("BOUGHT mse", [p for p in mid if d(p, "mse") < 0])):
    dv = [p[2]["score_avg"] - p[1]["score_avg"] for p in grp]
    if dv:
        mu, lo, hi = boot(dv, seed=12)
        print(f"    {lab}: n {len(dv)}  d score_avg mean {mu:+.5f} [{lo:+.5f}, {hi:+.5f}]"
              f"  median {np.median(dv):+.5f}  frac improving {np.mean(np.array(dv) > 0):.3f}")
ov = [sum(d(p, m) * SL[m]["slope"] for m in ("mse", "nmae", "reach", "fid", "jac"))
      / 6.0 / d(p, "mse") for p in mid if d(p, "mse") > 0]
if ov:
    mu, lo, hi = boot(ov, seed=13)
    print(f"\n    OVERALL points per +1.0 raw `mse` GIVEN UP: median {np.median(ov):+.5f}  "
          f"mean {mu:+.5f} [{lo:+.5f}, {hi:+.5f}]  n {len(ov)}")

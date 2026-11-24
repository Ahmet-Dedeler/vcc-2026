#!/usr/bin/env python3
"""Early warning: is anyone closing our ORDERING cushion?

`docs/SLOT-PROGRAM.md` and 222's decision tree both carry the same caution — our +0.0157
ordering lead is an *emission-side* lever (`rc_hybrid` + `rc_callset` bought it in two slots
without touching the profile), and any leader who notices it takes most of it in one slot.
The leader sits at `fid` rank ~390 and `jac` rank ~249 while holding `pds`/`mse` rank 1, so
they have not noticed yet. This script says when that changes.

Reads the snapshots in data/leaderboard/ and reports, per team, the move in the ordering sum
(score_fid + score_reach + score_jac) over the window. Run it after any board snapshot.

    uv run python scripts/cushion_watch.py [--days 2]
"""
import argparse, glob, json, os, sys

ORD = ("score_fid", "score_reach", "score_jac")
US = "zwIq3oGOWaXPhxAfioUl"


def load(path):
    d = json.load(open(path))
    rows = d["entries"] if isinstance(d, dict) else d
    return [r for r in rows if isinstance(r, dict) and r.get("score_avg") is not None]


def osum(r):
    return sum(r.get(m) or 0 for m in ORD)


def best(rows, key, val):
    c = [r for r in rows if r.get(key) == val]
    return max(c, key=lambda r: r["score_avg"]) if c else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", type=float, default=2.0)
    ap.add_argument("--min-overall", type=float, default=0.17)
    a = ap.parse_args()

    snaps = sorted(glob.glob("data/leaderboard/2026-*.json"), key=os.path.getmtime)
    if len(snaps) < 2:
        sys.exit("need at least two snapshots")
    cutoff = os.path.getmtime(snaps[-1]) - a.days * 86400
    window = [s for s in snaps if os.path.getmtime(s) >= cutoff] or snaps[-2:]
    first, last = load(window[0]), load(window[-1])
    print(f"window {os.path.basename(window[0])} -> {os.path.basename(window[-1])} "
          f"({len(window)} snapshots)")

    us = best(last, "team_id", US)
    moves = []
    for t in {r["team_name"] for r in last if r["score_avg"] >= a.min_overall}:
        p, q = best(first, "team_name", t), best(last, "team_name", t)
        if p and q:
            moves.append((osum(q) - osum(p), osum(q), q["score_avg"], t))
    moves.sort(reverse=True)

    print(f"\n{'d(ord sum)':>11} {'ord sum':>9} {'overall':>9}  team")
    for d, o, s, t in moves[:12]:
        flag = ""
        if us is not None and t != us.get("team_name"):
            # a rival both moving on ordering AND already ahead of us overall
            if d > 0.02 and s > us["score_avg"]:
                flag = "   <-- RIVAL CLOSING, ahead of us"
            elif o > osum(us):
                flag = "   <-- higher ordering sum than ours"
        print(f"{d:>+11.4f} {o:>9.4f} {s:>9.5f}  {t[:26]}{flag}")

    if us is not None:
        print(f"\nours: ordering sum {osum(us):.4f}, overall {us['score_avg']:.5f} "
              f"({us.get('model_name')})")
        print("NOTE: the displayed row lags the newest entry until the next snapshot.")


if __name__ == "__main__":
    main()

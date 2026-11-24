"""197 — convert the answer into what it demands: score at the leader's cosine, and the
board's own empirical price of a unit of direction.

Three things:
  1. What OUR `mse` would read at the leader's cosine, holding everything else at a15's.
  2. What a unit of implied cosine is worth in OVERALL on this board, measured across the
     informative cohort rather than modelled — this is the number allocation needs.
  3. The same inversion on the LIVE board (fetched into work/, read-only), because the
     mandated snapshot is ~12 h old and the top of the board moves ~0.005/day.

    uv run python target197.py
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
CTX = ("A", "B", "C")
B_ANCHOR, SPAN = 0.987862879, 0.950162566
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def main():
    inv = json.load(open(f"{WORK}/invert197.json"))
    ph3 = json.load(open(f"{WORK}/phase3197.json"))
    sumD, K, floor0, coeff = inv["sumD"], inv["K"], inv["floor0"], inv["coeff"]
    ours, exact = inv["ours"], ph3["exact"]
    ours_cos = ph3["ours_cos"]
    lead = ph3["leader"]

    P("=" * 100)
    P("197 — WHAT THE ANSWER DEMANDS")
    P("=" * 100)

    # ---- 1. our mse at the leader's cosine -------------------------------------------------
    per = ours["a15"]["per"]
    Kc = float(np.mean([per[c]["Kc"] for c in CTX]))
    claim = float(np.mean([per[c]["claim"] for c in CTX]))
    A = float(np.mean([per[c]["A"] for c in CTX]))
    B = float(np.mean([per[c]["B"] for c in CTX]))
    eps = A / B - 1.0
    q = claim / sumD

    def mse_at(c, qq=q, ee=eps):
        return 1.0 + Kc / sumD + qq * ee - 2.0 * c * np.sqrt(qq * (1.0 + ee))

    def sc(raw):
        return min(max((B_ANCHOR - raw) / SPAN, 0.0), 1.0)
    P("\n1. OUR `mse` at a given aggregate cosine, at our own claim and our own optimal alpha")
    P(f"   (a15's terms: claim {claim:.5f} = {100*claim/K:.1f}% of the cap, eps {eps:.4f}):")
    P(f"   {'cos':>8}{'raw mse':>11}{'score_mse':>11}{'overall vs a15':>16}{'what it is':<34}")
    for tag, c in (("a23 / today", exact["a23"]), ("a15", exact["a15"]),
                   ("Aginglab floor", 0.09478), ("leader floor", lead["floor"]),
                   ("leader floor x1.09 (tightness)", lead["floor"] * 1.09),
                   ("4-donor consensus", 0.08732),
                   ("k->inf ceiling (board)", 0.12263 * exact["a23"] / 0.03814)):
        raw = mse_at(c)
        P(f"   {c:>8.5f}{raw:>11.5f}{sc(raw):>11.5f}"
          f"{(sc(raw)-0.1200566)/6:>16.5f}   {tag:<34}")
    P("   (a15's own board score_mse is 0.1200566; the `overall` column is the `mse` member")
    P("    alone, /6, and does NOT include what the same cosine would buy on the other five.)")

    # ---- 2. the board's own price of a unit of direction ------------------------------------
    field = inv["field"]
    x = np.array([t["floor"] for t in field])
    P("\n2. THE BOARD'S OWN PRICE OF DIRECTION, binned rather than fitted.  A straight OLS over")
    P("   the informative cohort extrapolates badly (it predicts overall 0.343 at the leader's")
    P("   cosine against an observed 0.218), so the honest object is the conditional mean.")
    order = np.argsort(-x)
    P(f"   {'quartile of cos floor':<26}{'n':>4}{'cos':>9}{'overall':>10}{'raw pds':>10}"
      f"{'nmae':>9}{'reach':>9}{'fid':>9}")
    bins = np.array_split(order, 4)
    binned = []
    for i, b in enumerate(bins):
        g = [field[j] for j in b]
        row = dict(n=len(g), cos=float(np.mean([t["floor"] for t in g])),
                   overall=float(np.mean([t["overall"] for t in g])),
                   pds=float(np.mean([t["pds"] for t in g])),
                   nmae=float(np.mean([t["nmae"] for t in g])),
                   reach=float(np.mean([t["reach"] for t in g])),
                   fid=float(np.mean([t["fid"] for t in g])))
        binned.append(row)
        P(f"   {'Q%d (highest cos first)' % (i+1):<26}{row['n']:>4}{row['cos']:>9.5f}"
          f"{row['overall']:>10.5f}{row['pds']:>10.5f}{row['nmae']:>9.5f}"
          f"{row['reach']:>9.4f}{row['fid']:>9.4f}")
    us_row = next(t for t in field if "Ahmet" in t["team"])
    P(f"   {'US (abc-knee)':<26}{1:>4}{ours_cos:>9.5f}{us_row['overall']:>10.5f}"
      f"{us_row['pds']:>10.5f}{us_row['nmae']:>9.5f}{us_row['reach']:>9.4f}"
      f"{us_row['fid']:>9.4f}")
    slope = (binned[0]["overall"] - binned[1]["overall"]) / (binned[0]["cos"] - binned[1]["cos"])
    P(f"\n   Q4 is contaminated (mean `nmae` 1.278 — teams that are broken, not merely worse),")
    P(f"   so the usable comparison is Q1 vs Q2: {slope*0.01:+.5f} overall per +0.01 cosine,")
    P(f"   i.e. {slope*(lead['floor']-ours_cos):+.5f} over our own gap.  Cross-sectional and")
    P("   therefore confounded — teams with more direction also have better ORDERING.")
    P("   Note our own row against the same table: our `pds` (0.8368) beats even Q1's mean")
    P("   (0.8254) while our `reach` (0.1773) is below Q3's (0.1706 is close, Q2 is 0.2029) at a")
    P("   cosine between Q1's and Q2's.  That is 8.1's ordering deficit, seen from outside.")

    P("\n   The exact decomposition of the leader's lead over us, member by member (scaled):")
    lb = json.load(open(f"{REPO}/data/leaderboard/2026-09-03T0030Z.json"))["entries"]
    top = max(lb, key=lambda t: t["score_avg"])
    me = next(t for t in lb if "Ahmet" in t["team_name"])
    P(f"   {'member':<8}{'us':>11}{'%s' % top['team_name'][:12]:>14}{'delta':>11}"
      f"{'share of gap':>14}")
    gap = top["score_avg"] - me["score_avg"]
    dec = {}
    for m in ("pds", "mse", "nmae", "fid", "reach", "jac"):
        dm = (top[f"score_{m}"] - me[f"score_{m}"]) / 6.0
        dec[m] = dm
        P(f"   {m:<8}{me['score_'+m]:>11.5f}{top['score_'+m]:>14.5f}{dm:>11.5f}"
          f"{100*dm/gap:>13.1f}%")
    P(f"   {'TOTAL':<8}{me['score_avg']:>11.5f}{top['score_avg']:>14.5f}{gap:>11.5f}"
      f"{100.0:>13.1f}%")
    dirsum = dec["mse"] + dec["pds"] + dec["nmae"]
    P(f"\n   `mse` + `pds` + `nmae` — the three members 8.1 calls PURE DIRECTION — account for")
    P(f"   {dirsum:+.5f} of the {gap:+.5f} gap, i.e. {100*dirsum/gap:.0f}%.  `reach` adds "
      f"{dec['reach']:+.5f} more and is")
    P("   split between the two channels; `fid` and `jac` run the other way.")
    P(f"   Bracket for closing the direction gap: [{(sc(mse_at(lead['floor']))-0.1200566)/6:+.5f}"
      f" (`mse` alone, exact through the identity),")
    P(f"                                           {dirsum:+.5f} (`mse`+`pds`+`nmae`, the "
      f"members that read direction and nothing else)].")
    rows = dict(binned=binned, slope=slope, decomposition=dec, gap=gap)
    d = slope * (lead["floor"] - ours_cos)

    # ---- 3. the live board ------------------------------------------------------------------
    P("\n3. THE SAME INVERSION ON THE LIVE BOARD (fetched 2026-09-03, read-only, into work/).")
    live = f"{WORK}/lb-live.json"
    if os.path.exists(live):
        d2 = json.load(open(live))
        e = d2["entries"] if isinstance(d2, dict) else d2
        for t in e:
            t["_floor"] = (floor0 - t["expr_mse_unbiased_capped_norm"]) / coeff
        inf2 = sorted([t for t in e if t["_floor"] > 0], key=lambda t: -t["_floor"])
        P(f"   {len(e)} teams, {len(inf2)} informative (was {len(field)} of 612 at 00:30Z).")
        P(f"   {'#':>3} {'team':<24}{'overall':>9}{'raw mse':>9}{'cos floor':>11}{'x ours':>8}")
        for i, t in enumerate(inf2[:8]):
            P(f"   {i+1:>3} {t['team_name'][:23]:<24}{t['score_avg']:>9.5f}"
              f"{t['expr_mse_unbiased_capped_norm']:>9.5f}{t['_floor']:>11.5f}"
              f"{t['_floor']/ours_cos:>8.3f}")
        best = inf2[0]
        P(f"   the board's leading direction has moved to {best['_floor']:.5f} = "
          f"{best['_floor']/ours_cos:.2f}x ours in the ~12 h since the mandated snapshot.")
    out = dict(price=rows, gap_overall=d, mse_only=(sc(mse_at(lead["floor"])) - 0.1200566) / 6)
    json.dump(out, open(f"{WORK}/target197.json", "w"), indent=1, default=float)
    open(f"{WORK}/target197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/target197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

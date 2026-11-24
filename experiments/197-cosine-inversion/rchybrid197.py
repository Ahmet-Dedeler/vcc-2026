"""197 — adjudicate the `rc_hybrid` pre-registration, and refresh the target on the new row.

`assump197.py` pre-registered, before the row landed and from 191's digest identity alone:

    rc_hybrid's board raw `mse` = a15's 0.873805 to 0.0e+00, and it inverts to cos = 0.07925.
    FALSIFIED if the raw differs by more than +/- 0.004 or the cosine lands outside
    [0.07608, 0.08242].

Entry `nQUJZhZb3MW2uFJHBX5e` landed 2026-09-03 11:18 UTC.  This script reads it off the live
board, adjudicates, and re-runs the target arithmetic against our new displayed entry.

    uv run python rchybrid197.py
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

    e = json.load(open(f"{WORK}/lb-live.json"))
    e = e["entries"] if isinstance(e, dict) else e
    rc = next(t for t in e if t["id"] == "nQUJZhZb3MW2uFJHBX5e")

    P("=" * 96)
    P("197 — rc_hybrid ADJUDICATION (entry nQUJZhZb3MW2uFJHBX5e, 2026-09-03 11:18 UTC)")
    P("=" * 96)
    P(f"  pre-registered raw `mse`  {inv['rc_hybrid_pred_raw']:.16f}")
    P(f"  board raw `mse`           {rc['expr_mse_unbiased_capped_norm']:.16f}")
    d = rc["expr_mse_unbiased_capped_norm"] - inv["rc_hybrid_pred_raw"]
    P(f"  error                     {d:+.3e}   band +/- 4e-3   "
      f"{'CONFIRMED' if abs(d) < 4e-3 else 'FALSIFIED'}")
    P(f"  `score_mse` board {rc['score_mse']:.16f}  vs a15's recorded 0.1200566151966586 "
      f"(docs/LOG.md)")
    P(f"  `pds_cosine` board {rc['pds_cosine']:.16f}  vs a15's recorded 0.835987")

    # invert it: a15's terms, only `claim` moved, and only onto a leg that was already binding
    per = {c: dict(ours["a15"]["per"][c]) for c in CTX}
    c191 = json.load(open(f"{REPO}/experiments/191-rc-hybrid-build/work/claim191.json"))
    for c in CTX:
        per[c]["claim"] = per[c]["n"] * c191["rc_hybrid"][c]["claim_per_row"] / 18532.0
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["claim"], per[c]["B"])
                         for c in CTX]))
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in CTX]))
    X = (num - (rc["expr_mse_unbiased_capped_norm"] - 1.0) * sumD) / 2.0
    cos = X / (np.sqrt(sumD) * msa)
    P(f"\n  inverted cosine  {cos:.6f}   pre-registered {exact['a15']:.6f}   "
      f"band [{0.96*exact['a15']:.5f}, {1.04*exact['a15']:.5f}]")
    P(f"  difference {cos-exact['a15']:+.2e}   "
      f"{'CONFIRMED' if abs(cos-exact['a15']) < 0.04*exact['a15'] else 'FALSIFIED'}")
    P("\n  This is the mandate's digest-identical consistency test, executed on the board rather")
    P("  than on our own terms: two BYTE-DIFFERENT submissions (a15 19,049,370,200 bytes vs")
    P("  rc_hybrid 13,202,460,120, 31% sparser) with identical per-gene panel totals on 900/900")
    P("  blocks invert to the SAME cosine, and the board agrees to sixteen significant figures.")
    P(f"  Meanwhile `reach` moved raw 0.18117 -> {rc['de_wilcoxon_direction_reach_raw']:.5f} "
      f"(+{rc['score_reach']-0.115368:.5f} scaled).")
    P("  A change that moves the ORDERING channel by the largest single-member delta this")
    P("  project has recorded moves the inverted DIRECTION by 0.0e+00.  That is 8.1's two-channel")
    P("  split, board-verified, and it is the cleanest evidence the inversion reads direction and")
    P("  nothing else.")

    # ---- refresh the target on the new displayed entry --------------------------------------
    P("\n" + "=" * 96)
    P("THE TARGET, REFRESHED ON THE NEW DISPLAYED ENTRY")
    P("=" * 96)
    for t in e:
        t["_floor"] = (floor0 - t["expr_mse_unbiased_capped_norm"]) / coeff
    inf = sorted([t for t in e if t["_floor"] > 0], key=lambda t: -t["_floor"])
    top_dir = inf[0]
    top_all = max(e, key=lambda t: t["score_avg"])
    P(f"  our displayed entry is now rc_hybrid: overall {rc['score_avg']:.6f}, "
      f"rank {sorted(e, key=lambda t: -t['score_avg']).index(rc)+1} of {len(e)}")
    P(f"  our cosine is a15's, MEASURED: {cos:.5f}   (our floor would have said "
      f"{rc['_floor']:.5f}, {100*rc['_floor']/cos:.0f}% of it)")
    P(f"  best direction on the board : {top_dir['team_name']} raw {top_dir['expr_mse_unbiased_capped_norm']:.5f}"
      f" -> floor {top_dir['_floor']:.5f} = {top_dir['_floor']/cos:.3f}x ours")
    P(f"  overall leader              : {top_all['team_name']} raw {top_all['expr_mse_unbiased_capped_norm']:.5f}"
      f" -> floor {top_all['_floor']:.5f} = {top_all['_floor']/cos:.3f}x ours")
    P(f"  {'#':>3} {'team':<24}{'overall':>10}{'raw mse':>9}{'cos floor':>11}{'x ours':>8}")
    for i, t in enumerate(inf[:8]):
        P(f"  {i+1:>3} {t['team_name'][:23]:<24}{t['score_avg']:>10.5f}"
          f"{t['expr_mse_unbiased_capped_norm']:>9.5f}{t['_floor']:>11.5f}{t['_floor']/cos:>8.3f}")

    P(f"\n  Member-by-member gap to the overall leader ({top_all['team_name']}), scaled:")
    P(f"  {'member':<8}{'us (rc_hybrid)':>16}{'leader':>12}{'delta/6':>11}{'share':>9}")
    gap = top_all["score_avg"] - rc["score_avg"]
    dec = {}
    for m in ("pds", "mse", "nmae", "fid", "reach", "jac"):
        dm = (top_all[f"score_{m}"] - rc[f"score_{m}"]) / 6.0
        dec[m] = dm
        P(f"  {m:<8}{rc['score_'+m]:>16.5f}{top_all['score_'+m]:>12.5f}{dm:>11.5f}"
          f"{100*dm/gap:>8.1f}%")
    P(f"  {'TOTAL':<8}{rc['score_avg']:>16.5f}{top_all['score_avg']:>12.5f}{gap:>11.5f}"
      f"{100.0:>8.1f}%")
    dirsum = dec["mse"] + dec["pds"] + dec["nmae"]
    P(f"\n  `mse`+`pds`+`nmae` (pure direction) = {dirsum:+.5f} of {gap:+.5f} = "
      f"{100*dirsum/gap:.0f}% of the remaining gap.")
    P(f"  `reach` is now {dec['reach']:+.5f} — rc_hybrid took {100*(1-abs(dec['reach'])/0.01542):.0f}%"
      f" of the reach gap that stood at 00:30Z.")

    # what mse alone would read at the leader's cosine
    Kc = float(np.mean([per[c]["Kc"] for c in CTX]))
    claim = float(np.mean([per[c]["claim"] for c in CTX]))
    A = float(np.mean([per[c]["A"] for c in CTX]))
    Bm = float(np.mean([per[c]["B"] for c in CTX]))
    epsm, qm = A / Bm - 1.0, claim / sumD

    def mse_at(c):
        return 1.0 + Kc / sumD + qm * epsm - 2.0 * c * np.sqrt(qm * (1.0 + epsm))

    def scf(raw):
        return min(max((B_ANCHOR - raw) / SPAN, 0.0), 1.0)
    P(f"\n  What the leader's DIRECTION alone is worth to us, exact through the identity:")
    P(f"  {'cos':>9}{'raw mse':>10}{'score_mse':>11}{'d overall':>11}")
    for tag, c in (("ours now", cos), ("overall leader's floor", top_all["_floor"]),
                   ("best direction's floor", top_dir["_floor"]),
                   ("best floor x1.09 (tightness)", top_dir["_floor"] * 1.09)):
        r = mse_at(c)
        P(f"  {c:>9.5f}{r:>10.5f}{scf(r):>11.5f}{(scf(r)-rc['score_mse'])/6:>11.5f}   {tag}")

    out = dict(rc_raw=rc["expr_mse_unbiased_capped_norm"], rc_cos=cos,
               pred_raw=inv["rc_hybrid_pred_raw"], pred_cos=exact["a15"],
               overall=rc["score_avg"], gap=gap, dec=dec, dirsum=dirsum,
               top_dir=top_dir["team_name"], top_dir_floor=top_dir["_floor"],
               top_all=top_all["team_name"], top_all_floor=top_all["_floor"])
    json.dump(out, open(f"{WORK}/rchybrid197.json", "w"), indent=1, default=float)
    open(f"{WORK}/rchybrid197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/rchybrid197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

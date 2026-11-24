"""215 Phase 2 — sweep the WHOLE board for direction arms, and validate the instrument on
strangers' rows.

The validation is the important half.  For every arm-to-arm pair where BOTH endpoints are
informative (raw mse < 1.023611) the two channel readouts are independent:

    dcos_pds  = d(raw pds)  / 0.6479        <- amplitude-free, ordering-free
    dcos_mse  = -d(raw mse) / 2.074495      <- 197's floor difference, amplitude-CONTAMINATED

A pure-direction arm must satisfy dcos_pds = dcos_mse.  A pure-amplitude arm must satisfy
dcos_pds = 0 with dcos_mse != 0.  Nothing in the instrument forces either.
"""
import json, os
import numpy as np
import chan215 as C
from phase1215 import series, usable

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
M = C.MEMBERS
CAL_LO, CAL_HI = 0.79, 0.88


def all_pairs(panel):
    teams = sorted({r["team"] for r in panel})
    out = []
    for t in teams:
        rs = series(panel, t)
        for b, a in zip(rs[:-1], rs[1:]):
            if not (usable(b) and usable(a)):
                continue
            rb = {m: b["raw_" + m] for m in M}
            ra = {m: a["raw_" + m] for m in M}
            r = C.classify(rb, ra)
            out.append({
                "team": t, "from": b["model"], "to": a["model"],
                "date_from": b["date"][:16], "date_to": a["date"][:16],
                "n_from": b["n"], "n_to": a["n"], "dn": (a["n"] or 0) - (b["n"] or 0),
                "snaps_from": b["n_snaps"], "snaps_to": a["n_snaps"],
                "pds_from": b["raw_pds"], "pds_to": a["raw_pds"],
                "mse_from": b["raw_mse"], "mse_to": a["raw_mse"],
                "ovr_from": b["overall"], "ovr_to": a["overall"],
                "in_cal_band": bool(CAL_LO <= b["raw_pds"] <= CAL_HI
                                    and CAL_LO <= a["raw_pds"] <= CAL_HI),
                "d_pds": r["d_raw"]["pds"], "d_mse": r["d_raw"]["mse"],
                "d_nmae": r["d_raw"]["nmae"], "d_reach": r["d_raw"]["reach"],
                "d_fid": r["d_raw"]["fid"], "d_jac": r["d_raw"]["jac"],
                **{k: r[k] for k in ("label", "rho", "dcos_pds", "dcos_mse", "dcos_pds_sd",
                                     "sigma_pds_eff", "b_direction", "b_ordering",
                                     "b_amplitude", "b_geometry", "ovr_obs",
                                     "informative_before", "informative_after",
                                     "both_informative", "cos_floor_before", "cos_floor_after",
                                     "mse_priced")},
            })
    return out


def main():
    panel = json.load(open(os.path.join(WORK, "panel215.json")))
    pairs = all_pairs(panel)
    print(f"board-wide: {len({p['team'] for p in pairs})} teams with >=2 usable observations, "
          f"{len(pairs)} arm-to-arm pairs")

    # ---------------------------------------------------------------- VALIDATION
    v = [p for p in pairs if p["both_informative"] and p["in_cal_band"]
         and abs(p["dcos_mse"]) > 0.005]
    print("\n" + "=" * 116)
    print("V1 — the two independent readouts on STRANGERS' rows "
          f"(both endpoints informative, pds in [{CAL_LO},{CAL_HI}], |dcos_mse| > 0.005)")
    print("=" * 116)
    print(f"{'team':<18}{'from -> to':<30}{'dn':>3}{'dcos_pds':>10}{'dcos_mse':>10}"
          f"{'rho':>8}  {'label':<22}")
    for p in sorted(v, key=lambda q: -abs(q["dcos_mse"])):
        nm = f"{(p['from'] or '?')[:12]} -> {(p['to'] or '?')[:12]}"
        print(f"{p['team'][:17]:<18}{nm:<30}{p['dn']:>3}{p['dcos_pds']:>+10.4f}"
              f"{p['dcos_mse']:>+10.4f}{p['rho']:>8.2f}  {p['label']:<22}")
    if len(v) >= 4:
        a = np.array([p["dcos_pds"] for p in v]); b = np.array([p["dcos_mse"] for p in v])
        sl = float(np.dot(a, b) / np.dot(b, b))
        rr = float(np.corrcoef(a, b)[0, 1])
        print(f"\n  regression through the origin: dcos_pds = {sl:.3f} * dcos_mse "
              f"(a pure-direction field would give 1.000);  Pearson r = {rr:+.3f}, n = {len(v)}")

    # the direction-labelled subset only
    vd = [p for p in v if p["label"].startswith("DIRECTION") and p["dcos_mse"] > 0]
    va = [p for p in v if p["label"] == "AMPLITUDE"]
    for nm, sub in (("DIRECTION-labelled", vd), ("AMPLITUDE-labelled", va)):
        if sub:
            r = np.array([p["rho"] for p in sub])
            print(f"  {nm:<20} n={len(sub):>2}  rho median {np.median(r):+.3f} "
                  f"[{r.min():+.3f}, {r.max():+.3f}]   (direction -> 1, amplitude -> 0)")

    # ---------------------------------------------------------------- the direction arms
    # A pair is a DIRECTION arm if the amplitude-free readout moved and the other three buckets
    # are each under half of it.
    is_dir = lambda p: (p["dcos_pds"] > 3 * p["dcos_pds_sd"]
                        and abs(p["b_ordering"]) < 0.5 * abs(p["b_direction"])
                        and abs(p["b_amplitude"]) < 0.5 * abs(p["b_direction"])
                        and abs(p["b_geometry"]) < 0.5 * abs(p["b_direction"]))
    dirs = [p for p in pairs if is_dir(p)]

    # (A) MEASURED: both endpoints informative -> Delta cosine is a difference of two 197 floors,
    #     no pds->cosine conversion needed at all.
    meas = sorted([p for p in dirs if p["both_informative"]],
                  key=lambda p: -p["dcos_mse"])
    print("\n" + "=" * 130)
    print("PHASE 2A — DIRECTION arms whose cosine move is MEASURED "
          "(both endpoints raw mse < 1.023611, so both floors are informative)")
    print("=" * 130)
    print(f"{'#':>3} {'team':<17}{'from -> to':<32}{'date_to':<12}{'dn':>3}"
          f"{'d cos (floor)':>14}{'d cos (pds)':>12}{'rho':>7}{'cos_from':>10}{'cos_to':>9}{'dovr':>9}")
    for i, p in enumerate(meas, 1):
        nm = f"{(p['from'] or '?')[:13]} -> {(p['to'] or '?')[:13]}"
        print(f"{i:>3} {p['team'][:16]:<17}{nm:<32}{p['date_to'][:10]:<12}{p['dn']:>3}"
              f"{p['dcos_mse']:>+14.4f}{p['dcos_pds']:>+12.4f}{p['rho']:>7.2f}"
              f"{p['cos_floor_before']:>10.5f}{p['cos_floor_after']:>9.5f}{p['ovr_obs']:>+9.4f}")

    # (B) BOUNDED: the arm ENDS informative but starts vacuous -> the endpoint cosine is a floor
    #     and the move is bounded below by (floor_after - 0), reported as a bound not a point.
    part = sorted([p for p in dirs if p["informative_after"] and not p["informative_before"]],
                  key=lambda p: -p["cos_floor_after"])
    print("\n" + "=" * 130)
    print("PHASE 2B — DIRECTION arms that ARRIVE in the informative region "
          "(start vacuous, end informative): the endpoint is a BOUND, the move is not measurable")
    print("=" * 130)
    print(f"{'#':>3} {'team':<17}{'from -> to':<32}{'date_to':<12}{'dn':>3}"
          f"{'raw mse from':>13}{'raw mse to':>12}{'cos_to >=':>11}{'d pds':>9}{'dovr':>9}")
    for i, p in enumerate(part[:20], 1):
        nm = f"{(p['from'] or '?')[:13]} -> {(p['to'] or '?')[:13]}"
        print(f"{i:>3} {p['team'][:16]:<17}{nm:<32}{p['date_to'][:10]:<12}{p['dn']:>3}"
              f"{p['mse_from']:>13.3f}{p['mse_to']:>12.4f}{p['cos_floor_after']:>11.5f}"
              f"{p['d_pds']:>+9.4f}{p['ovr_obs']:>+9.4f}")

    # (C) everything else: a pds move with NO cosine reading available at all.
    rest = [p for p in dirs if not p["informative_after"]]
    inband_rest = [p for p in rest if p["in_cal_band"]]
    print(f"\nPHASE 2C — {len(rest)} further DIRECTION-labelled arms end with raw mse >= 1.023611.")
    print("  For these the board carries NO cosine information: the mse floor is vacuous and the")
    print(f"  pds -> cosine rate is calibrated only inside pds [{CAL_LO},{CAL_HI}] "
          f"({len(inband_rest)} of the {len(rest)} are in band).")
    print("  They are reported as raw pds moves and NOT converted to cosine.")
    print(f"\n{'#':>3} {'team':<17}{'from -> to':<32}{'date_to':<12}{'dn':>3}"
          f"{'pds_from':>10}{'pds_to':>9}{'d pds':>9}{'dovr':>9}{'band':>6}")
    for i, p in enumerate(sorted(rest, key=lambda q: -q["d_pds"])[:18], 1):
        nm = f"{(p['from'] or '?')[:13]} -> {(p['to'] or '?')[:13]}"
        print(f"{i:>3} {p['team'][:16]:<17}{nm:<32}{p['date_to'][:10]:<12}{p['dn']:>3}"
              f"{p['pds_from']:>10.4f}{p['pds_to']:>9.4f}{p['d_pds']:>+9.4f}"
              f"{p['ovr_obs']:>+9.4f}{'in' if p['in_cal_band'] else 'OUT':>6}")

    # ---------------------------------------------------------------- the informative frontier
    print("\n" + "=" * 130)
    print("PHASE 2D — WHEN each of the board's top cosines first became visible, and where "
          "the team's pds was then")
    print("=" * 130)
    lb = json.load(open(os.path.join(WORK, "lb-live.json")))["entries"]
    top = sorted([e for e in lb if e["expr_mse_unbiased_capped_norm"] is not None
                  and C.is_informative(e["expr_mse_unbiased_capped_norm"])],
                 key=lambda e: e["expr_mse_unbiased_capped_norm"])[:14]
    print(f"{'team':<19}{'model':<16}{'cos floor':>10}{'raw pds':>9}   "
          f"{'first informative row':<34}{'pds then':>9}{'pds now':>9}")
    for e in top:
        t = e["team_name"]
        rs = [r for r in series(panel, t) if usable(r)]
        firstinf = next((r for r in rs if C.is_informative(r["raw_mse"])), None)
        fi = (f"{firstinf['date'][:10]} #{firstinf['n']} {str(firstinf['model'])[:14]}"
              if firstinf else "-")
        print(f"{t[:18]:<19}{str(e['model_name'])[:15]:<16}"
              f"{C.cos_floor(e['expr_mse_unbiased_capped_norm']):>10.5f}"
              f"{e['pds_cosine']:>9.4f}   {fi:<34}"
              f"{(firstinf['raw_pds'] if firstinf else float('nan')):>9.4f}"
              f"{e['pds_cosine']:>9.4f}")

    json.dump({"pairs": pairs, "validation": v, "direction_arms": dirs,
               "measured": meas, "bounded": part},
              open(os.path.join(WORK, "phase2215.json"), "w"), indent=1, default=float)
    print(f"\nwrote work/phase2215.json  ({len(pairs)} pairs, {len(dirs)} direction arms, "
          f"{len(meas)} with a MEASURED cosine move)")


if __name__ == "__main__":
    main()

"""215 Phase 1 — decompose every leader trajectory, arm by arm, and report a channel budget.

Reads work/panel215.json (built and gated by panel215.py).  Uses the instrument in chan215.py,
which passed all five known-answer arms in phase0215.py.
"""
import json, os, sys
import numpy as np
import chan215 as C

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
M = C.MEMBERS

TEAMS = ["cqawesome", "Aginglab.com", "Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng",
         "Garrrbage", "Slurm Shadi", "Haystack", "Ahmet Dedeler"]
# the pds band inside which the direction rate E_DIR["pds"] = 0.6479 was calibrated (a14/a15)
CAL_LO, CAL_HI = 0.79, 0.88


def series(panel, team):
    rs = sorted([r for r in panel if r["team"] == team], key=lambda r: r["date"])
    # drop exact re-observations of the same entry id (already deduped) and identical relaunches
    out = []
    for r in rs:
        if out and all(abs(r["raw_" + m] - out[-1]["raw_" + m]) < 1e-12 for m in M):
            out[-1]["dup_of_next"] = True      # e.g. cqawesome mp4 / mp4_relaunch1
            continue
        out.append(r)
    return out


def usable(r):
    """A row whose DE stage actually ran.  Broken/shakedown submissions read fid at chance-or-worse
    (the board's broken rows sit at fid 0.23-0.45 against a working 0.49-0.52) or nmae above 1.05."""
    return r["raw_nmae"] <= 1.05 and r["raw_fid"] >= 0.45


def pair_rows(panel, team):
    rs = series(panel, team)
    rows = []
    for b, a in zip(rs[:-1], rs[1:]):
        rb = {m: b["raw_" + m] for m in M}
        ra = {m: a["raw_" + m] for m in M}
        r = C.classify(rb, ra)
        dn = (a["n"] or 0) - (b["n"] or 0)
        in_band = CAL_LO <= b["raw_pds"] <= CAL_HI and CAL_LO <= a["raw_pds"] <= CAL_HI
        rows.append({
            "team": team, "from": b["model"], "to": a["model"],
            "date_from": b["date"][:16], "date_to": a["date"][:16],
            "n_from": b["n"], "n_to": a["n"], "dn": dn,
            "snaps_from": b["n_snaps"], "snaps_to": a["n_snaps"],
            "pds_from": b["raw_pds"], "pds_to": a["raw_pds"],
            "mse_from": b["raw_mse"], "mse_to": a["raw_mse"],
            "ovr_from": b["overall"], "ovr_to": a["overall"],
            "in_cal_band": bool(in_band),
            "usable": bool(usable(b) and usable(a)),
            **{k: r[k] for k in ("label", "rho", "dcos_pds", "dcos_mse", "dcos_pds_sd",
                                 "sigma_pds_eff", "geom_share", "b_direction", "b_ordering",
                                 "b_amplitude", "b_geometry", "ovr_obs",
                                 "informative_before", "informative_after", "both_informative",
                                 "cos_floor_before", "cos_floor_after")},
            "d_raw": r["d_raw"], "sig": r["sig"],
        })
    return rows, rs


def main():
    panel = json.load(open(os.path.join(WORK, "panel215.json")))
    all_rows, summary = [], []

    for t in TEAMS:
        rows, rs = pair_rows(panel, t)
        if not rows:
            continue
        all_rows += rows
        print("\n" + "=" * 132)
        print(f"### {t}   —  {len(rs)} distinct observations, latest submission #{rs[-1]['n']}, "
              f"overall {rs[0]['overall']:+.4f} -> {rs[-1]['overall']:+.4f}")
        print("=" * 132)
        print(f"{'from -> to':<34}{'dn':>3}{'dpds':>9}{'dmse':>10}{'dnmae':>9}{'dreach':>9}"
              f"{'rho':>8}{'dcosP':>9}{'label':<20}{'dovr':>9}  {'dir':>8}{'ord':>8}{'amp':>8}{'geo':>8}")
        for r in rows:
            tag = ""
            if not r["in_cal_band"]:
                tag += "*"
            if r["dn"] > 1:
                tag += "+"
            if not r["both_informative"]:
                tag += "v"
            if not r["usable"]:
                tag += "B"
            nm = f"{(r['from'] or '?')[:14]} -> {(r['to'] or '?')[:14]}"
            print(f"{nm:<34}{r['dn']:>3}{r['d_raw']['pds']:>+9.4f}{r['d_raw']['mse']:>+10.4f}"
                  f"{r['d_raw']['nmae']:>+9.4f}{r['d_raw']['reach']:>+9.4f}"
                  f"{r['rho']:>8.2f}{r['dcos_pds']:>+9.4f}{r['label']+tag:<20}"
                  f"{r['ovr_obs']:>+9.4f}  {r['b_direction']:>+8.4f}{r['b_ordering']:>+8.4f}"
                  f"{r['b_amplitude']:>+8.4f}{r['b_geometry']:>+8.4f}")
        clean = [r for r in rows if r["usable"]]
        tot = {k: sum(r[k] for r in clean) for k in
               ("b_direction", "b_ordering", "b_amplitude", "b_geometry", "ovr_obs")}
        cs = [r for r in rs if usable(r)]
        dpds_tot = (cs[-1]["raw_pds"] - cs[0]["raw_pds"]) if len(cs) > 1 else 0.0
        # --- direction-arrival concentration: how much of the pds climb is in the biggest step?
        dps = [r["d_raw"]["pds"] for r in clean]
        pos = sum(x for x in dps if x > 0)
        big = max(dps) if dps else 0.0
        n_res = sum(1 for r in clean if abs(r["d_raw"]["pds"]) > C.K_SIG * r["sigma_pds_eff"])
        # when did the prediction reach 95% of its final rank quality (measured from the first
        # usable row's pds as the origin)?
        arrive = None
        if len(cs) > 1 and dpds_tot > 1e-9:
            for r in cs:
                if (r["raw_pds"] - cs[0]["raw_pds"]) >= 0.95 * dpds_tot:
                    arrive = r
                    break
        summary.append({"team": t, "n_obs": len(rs), "n_pairs": len(rows),
                        "n_clean_pairs": len(clean),
                        "first_date": rs[0]["date"][:10], "last_date": rs[-1]["date"][:10],
                        "clean_first_date": cs[0]["date"][:10] if cs else None,
                        "clean_first_model": cs[0]["model"] if cs else None,
                        "clean_first_n": cs[0]["n"] if cs else None,
                        "n_last": rs[-1]["n"],
                        "ovr_first": rs[0]["overall"], "ovr_last": rs[-1]["overall"],
                        "pds_first": rs[0]["raw_pds"], "pds_last": rs[-1]["raw_pds"],
                        "pds_clean_first": cs[0]["raw_pds"] if cs else None,
                        "dpds_total": dpds_tot,
                        "dpds_biggest_step": big,
                        "dpds_biggest_share": (big / pos) if pos > 1e-9 else float("nan"),
                        "n_resolvable_pds_steps": n_res,
                        "arrive_date": arrive["date"][:10] if arrive else None,
                        "arrive_model": arrive["model"] if arrive else None,
                        "arrive_n": arrive["n"] if arrive else None,
                        "mse_first": rs[0]["raw_mse"], "mse_last": rs[-1]["raw_mse"],
                        "cos_floor_last": C.cos_floor(rs[-1]["raw_mse"]),
                        "informative_last": C.is_informative(rs[-1]["raw_mse"]),
                        **tot})
        print(f"{'TOTAL':<34}{'':>3}{dpds_tot:>+9.4f}{'':>10}{'':>9}{'':>9}{'':>8}{'':>9}{'':<20}"
              f"{tot['ovr_obs']:>+9.4f}  {tot['b_direction']:>+8.4f}{tot['b_ordering']:>+8.4f}"
              f"{tot['b_amplitude']:>+8.4f}{tot['b_geometry']:>+8.4f}")

    print("\n\n" + "=" * 132)
    print("PER-TEAM CHANNEL BUDGET  (overall-score units; the four buckets sum to the observed climb)")
    print("=" * 132)
    print(f"{'team':<16}{'obs':>4}{'first':>11}{'last':>11}{'climb':>9}"
          f"{'direction':>11}{'ordering':>10}{'amplitude':>11}{'geometry':>10}"
          f"{'  %dir':>7}{'%ord':>6}{'%amp':>6}{'%geo':>6}{'dpds':>9}{'cosfloor':>10}")
    for s in summary:
        c = s["ovr_obs"]
        pct = lambda x: 100 * x / c if abs(c) > 1e-9 else float("nan")
        cf = s["cos_floor_last"] if s["informative_last"] else float("nan")
        print(f"{s['team'][:15]:<16}{s['n_obs']:>4}{s['first_date']:>11}{s['last_date']:>11}"
              f"{c:>+9.4f}{s['b_direction']:>+11.4f}{s['b_ordering']:>+10.4f}"
              f"{s['b_amplitude']:>+11.4f}{s['b_geometry']:>+10.4f}"
              f"{pct(s['b_direction']):>7.0f}{pct(s['b_ordering']):>6.0f}"
              f"{pct(s['b_amplitude']):>6.0f}{pct(s['b_geometry']):>6.0f}"
              f"{s['dpds_total']:>+9.4f}{cf:>10.5f}")

    print("\nflags:  * pds outside the [0.79,0.88] band where the direction rate is calibrated")
    print("        + more than one submission elapsed between the two observations (aggregate)")
    print("        v at least one endpoint has raw mse >= 1.023611, so its cosine FLOOR is vacuous")
    print("        B at least one endpoint is a BROKEN/shakedown row (fid < 0.45 or nmae > 1.05);")
    print("          B pairs are excluded from the per-team budget below")

    print("\n" + "=" * 132)
    print("DIRECTION-ARRIVAL PROFILE  (pds is the amplitude-free, ordering-free readout)")
    print("=" * 132)
    print(f"{'team':<16}{'clean start':>13}{'pds0':>8}{'pdsN':>8}{'d pds':>9}"
          f"{'biggest step':>14}{'share':>7}{'nsteps':>8}   {'95% of the climb reached at':<34}")
    for s in summary:
        print(f"{s['team'][:15]:<16}{str(s['clean_first_date']):>13}"
              f"{(s['pds_clean_first'] or 0):>8.4f}{s['pds_last']:>8.4f}{s['dpds_total']:>+9.4f}"
              f"{s['dpds_biggest_step']:>+14.4f}{s['dpds_biggest_share']:>7.0%}"
              f"{s['n_resolvable_pds_steps']:>8d}   "
              f"{str(s['arrive_date'])} #{s['arrive_n']} {str(s['arrive_model'])[:20]:<22}")

    json.dump({"pairs": all_rows, "summary": summary},
              open(os.path.join(WORK, "phase1215.json"), "w"), indent=1, default=float)
    print(f"\nwrote work/phase1215.json  ({len(all_rows)} pairs, {len(summary)} teams)")


if __name__ == "__main__":
    main()

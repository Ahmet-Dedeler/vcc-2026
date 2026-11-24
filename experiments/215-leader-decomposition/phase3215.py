"""215 Phase 3 — did the leaders' direction arrive in one step, accumulate in many, or was it
there from the start?

Three tests, all on `pds`, the amplitude-free ordering-free readout:
  T1  did ANY team's earliest observed submission already carry a high direction?
  T2  how concentrated is each leader's pds climb -- one step, or many?
  T3  Haystack: ordering structural from day one against a direction that climbed steadily.
"""
import json, os
import numpy as np
import chan215 as C
from phase1215 import series, usable

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
M = C.MEMBERS
LEADERS = ["cqawesome", "jw_13288", "Hongjian Zhou", "Andrew Liu", "MT Zeng",
           "Garrrbage", "Aginglab.com", "Slurm Shadi", "Haystack", "Ahmet Dedeler"]


def main():
    panel = json.load(open(os.path.join(WORK, "panel215.json")))
    ph1 = json.load(open(os.path.join(WORK, "phase1215.json")))
    teams = sorted({r["team"] for r in panel})

    # ---------------------------------------------------------------- T1
    print("=" * 122)
    print("T1 — DID ANYONE START HIGH?  Every team's EARLIEST observed usable row, "
          "ranked by the direction it already carried")
    print("=" * 122)
    rows = []
    for t in teams:
        rs = [r for r in series(panel, t) if usable(r)]
        if not rs:
            continue
        f, l = rs[0], rs[-1]
        rows.append({"team": t, "n_first": f["n"], "date_first": f["date"][:10],
                     "model_first": f["model"], "pds_first": f["raw_pds"],
                     "mse_first": f["raw_mse"],
                     "cos_first": C.cos_floor(f["raw_mse"]) if C.is_informative(f["raw_mse"]) else None,
                     "pds_last": l["raw_pds"], "mse_last": l["raw_mse"],
                     "cos_last": C.cos_floor(l["raw_mse"]) if C.is_informative(l["raw_mse"]) else None,
                     "n_last": l["n"], "n_obs": len(rs)})
    # who was already informative on their earliest OBSERVED row?
    early = [r for r in rows if r["cos_first"] is not None and r["cos_first"] > 0.03]
    early.sort(key=lambda r: -r["cos_first"])
    print(f"teams with >=1 usable observation: {len(rows)}")
    print(f"of those, already informative (raw mse < 1.023611) on their earliest OBSERVED row "
          f"with cos floor > 0.03: {len(early)}")
    print(f"\n{'team':<20}{'earliest observed row':<34}{'#':>4}{'cos floor then':>15}"
          f"{'pds then':>10}{'cos floor now':>14}{'pds now':>9}")
    for r in early[:16]:
        cl = r["cos_last"] if r["cos_last"] is not None else float("nan")
        print(f"{r['team'][:19]:<20}"
              f"{r['date_first']+' '+str(r['model_first'])[:20]:<34}{r['n_first']:>4}"
              f"{r['cos_first']:>15.5f}{r['pds_first']:>10.4f}{cl:>14.5f}{r['pds_last']:>9.4f}")
    best_start = max((r["cos_first"] for r in early), default=float("nan"))
    print(f"\n  HIGHEST cosine anyone has EVER carried on their earliest observed row: {best_start:.5f}")
    print(f"  HIGHEST cosine on the board today: 0.13211 (cqawesome m0903c)")
    # and the same restricted to genuinely early submissions (n <= 6)
    veryearly = [r for r in early if (r["n_first"] or 99) <= 6]
    veryearly.sort(key=lambda r: -r["cos_first"])
    print(f"\n  restricted to submissions #1-#6 (a team's first days): {len(veryearly)} teams, "
          f"best cos floor {max((r['cos_first'] for r in veryearly), default=float('nan')):.5f}")
    for r in veryearly[:8]:
        print(f"    {r['team'][:19]:<20}#{r['n_first']:<3} {r['date_first']}  "
              f"cos {r['cos_first']:.5f}  pds {r['pds_first']:.4f}")

    # ---------------------------------------------------------------- T2
    print("\n" + "=" * 122)
    print("T2 — CONCENTRATION OF THE CLIMB.  Gini-style: what share of each leader's total pds "
          "gain sits in its 1 / 2 / 3 largest steps?")
    print("=" * 122)
    print(f"{'team':<16}{'clean pairs':>12}{'resolvable':>11}{'total dpds':>12}"
          f"{'top1':>7}{'top2':>7}{'top3':>7}   {'reading':<30}")
    conc = []
    for t in LEADERS:
        prs = [p for p in ph1["pairs"] if p["team"] == t and p["usable"]]
        if not prs:
            continue
        dps = sorted([p["d_raw"]["pds"] for p in prs], reverse=True)
        pos = sum(x for x in dps if x > 0)
        nres = sum(1 for p in prs
                   if abs(p["d_raw"]["pds"]) > C.K_SIG * p["sigma_pds_eff"])
        sh = [sum(dps[:k]) / pos if pos > 1e-9 else float("nan") for k in (1, 2, 3)]
        if pos < 0.005:
            read = "FLAT since first observation"
        elif sh[0] > 0.70:
            read = "ONE step"
        elif sh[1] > 0.80:
            read = "TWO steps"
        else:
            read = "MANY steps (incremental)"
        conc.append({"team": t, "n_pairs": len(prs), "n_res": nres, "pos": pos,
                     "top1": sh[0], "top2": sh[1], "top3": sh[2], "reading": read})
        print(f"{t[:15]:<16}{len(prs):>12}{nres:>11}{pos:>+12.4f}"
              f"{sh[0]:>7.0%}{sh[1]:>7.0%}{sh[2]:>7.0%}   {read:<30}")

    # ---------------------------------------------------------------- T3
    print("\n" + "=" * 122)
    print("T3 — HAYSTACK, the natural experiment: ORDERING structural from day one, "
          "DIRECTION climbing steadily")
    print("=" * 122)
    hs = [r for r in series(panel, "Haystack") if usable(r)]
    us = [r for r in series(panel, "Ahmet Dedeler") if usable(r)]
    print(f"{'date':<12}{'#':>3} {'model':<14}{'pds':>9}{'raw mse':>10}{'cos floor':>11}"
          f"{'reach':>9}{'fid':>8}{'jac':>8}{'overall':>9}")
    for r in hs:
        cf = C.cos_floor(r["raw_mse"]) if C.is_informative(r["raw_mse"]) else float("nan")
        print(f"{r['date'][:10]:<12}{r['n']:>3} {str(r['model'])[:13]:<14}{r['raw_pds']:>9.4f}"
              f"{r['raw_mse']:>10.4f}{cf:>11.5f}{r['raw_reach']:>9.4f}{r['raw_fid']:>8.4f}"
              f"{r['raw_jac']:>8.4f}{r['overall']:>9.4f}")
    d_reach = hs[-1]["raw_reach"] - hs[0]["raw_reach"]
    d_pds = hs[-1]["raw_pds"] - hs[0]["raw_pds"]
    print(f"\n  Haystack over its whole visible history:  raw reach {hs[0]['raw_reach']:.4f} -> "
          f"{hs[-1]['raw_reach']:.4f}  ({d_reach:+.4f}),  raw pds {hs[0]['raw_pds']:.4f} -> "
          f"{hs[-1]['raw_pds']:.4f}  ({d_pds:+.4f})")
    print(f"  Its THIRD submission already read reach {hs[0]['raw_reach']:.4f} at pds "
          f"{hs[0]['raw_pds']:.4f}.")
    print(f"  Ours read reach {us[0]['raw_reach']:.4f} at submission #{us[0]['n']} and only reached "
          f"{us[-1]['raw_reach']:.4f} at #{us[-1]['n']}, after two arms designed for it.")
    # reach conditional on pds: is Haystack's reach explained by its pds?
    live = json.load(open(os.path.join(WORK, "lb-live.json")))["entries"]
    ok = [e for e in live if e["pds_cosine"] is not None and e["de_wilcoxon_direction_reach_raw"]
          is not None and e["de_wilcoxon_lfc_nmae"] <= 1.05 and e["de_wilcoxon_direction_fidelity_yield_raw"] >= 0.45]
    x = np.array([e["pds_cosine"] for e in ok]); y = np.array([e["de_wilcoxon_direction_reach_raw"] for e in ok])
    A = np.vstack([x, np.ones(len(x))]).T
    co = np.linalg.lstsq(A, y, rcond=None)[0]
    resid = y - A @ co
    sd = resid.std(ddof=2)
    def z_of(team):
        e = [q for q in ok if q["team_name"] == team]
        if not e:
            return float("nan")
        e = e[0]
        return (e["de_wilcoxon_direction_reach_raw"] - (co[0] * e["pds_cosine"] + co[1])) / sd
    print(f"\n  Board-wide reach ~ pds regression on {len(ok)} live usable rows: "
          f"reach = {co[0]:+.4f}*pds {co[1]:+.4f},  residual sd {sd:.4f}")
    print(f"  {'team':<18}{'reach residual z':>18}")
    for t in ["Haystack", "Ahmet Dedeler", "cqawesome", "Aginglab.com", "Hongjian Zhou",
              "Garrrbage", "Slurm Shadi"]:
        print(f"  {t:<18}{z_of(t):>+18.2f}")

    json.dump({"start_high": rows, "concentration": conc,
               "reach_vs_pds": {"slope": float(co[0]), "intercept": float(co[1]),
                                "resid_sd": float(sd), "n": len(ok)}},
              open(os.path.join(WORK, "phase3215.json"), "w"), indent=1, default=float)
    print("\nwrote work/phase3215.json")


if __name__ == "__main__":
    main()

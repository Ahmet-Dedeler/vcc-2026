"""215 — out-of-sample validation of the instrument on STRANGERS' rows, and the decisive
amplitude-exclusion test for the three leaders whose cosine move is measurable.

V1 (in phase2215.py)  the two independent readouts agree on DIRECTION-labelled arms (rho -> 1)
                      and disagree on AMPLITUDE-labelled arms (rho -> 0)
V2 (here)             the amplitude channel's leakage into pds, measured on the field
V3 (here)             can each leader's pds climb be explained by amplitude leakage?  No.
"""
import json, os
import numpy as np
import chan215 as C

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
M = C.MEMBERS


def main():
    d = json.load(open(os.path.join(WORK, "phase2215.json")))
    P = [p for p in d["pairs"] if p["label"] == "AMPLITUDE" and abs(p["d_mse"]) > 0.02]
    r = np.array([abs(p["d_pds"]) / abs(p["d_mse"]) for p in P])
    p90 = float(np.percentile(r, 90))
    print("=" * 104)
    print(f"V2 — the amplitude channel's pds leakage, on {len(P)} AMPLITUDE-labelled strangers' arms")
    print("=" * 104)
    print(f"  |d pds| / |d mse| :  median {np.median(r):.5f}   p90 {p90:.5f}   max {r.max():.5f}")
    print(f"  our own calibration (a10 -> a11, alpha 1.0796 -> 0.838): {C.PDS_LEAK:.5f}")
    print("  The field distribution is an UPPER envelope, not a measurement of pure leakage: an arm")
    print("  that mixes a little direction into an amplitude move lands in this bucket too.  The p90")
    print(f"  {p90:.4f} is therefore used below as the CONSERVATIVE leak bound, {p90/C.PDS_LEAK:.1f}x our own.")

    print("\n" + "=" * 104)
    print("V3 — could each leader's pds climb be amplitude leakage instead of direction?")
    print("=" * 104)
    panel = json.load(open(os.path.join(WORK, "panel215.json")))

    def row(team, model):
        q = [x for x in panel if x["team"] == team and x["model"] == model]
        assert len(q) == 1, (team, model, len(q))
        return q[0]

    cases = [
        ("cqawesome", "m0830b", "m0903c", "the whole informative window"),
        ("jw_13288", "willow_34", "juniper_46", "the cluster's own informative window"),
        ("Hongjian Zhou", "pm-r14b", "pm-r14f", "after the cluster became informative"),
        ("Andrew Liu", "AL-4746", "AL-1466", "after the cluster became informative"),
        ("MT Zeng", "Kestrel4", "Kestrel2", "their own Kestrel line, skipping Petrel1"),
        ("Aginglab.com", "GeroAI_v11", "GeroAI_v19", "from their first informative row"),
        ("Haystack", "GSA-015", "GSA-023", "from their first informative row"),
        ("Ahmet Dedeler", "The Houseplant Has Notes", "Peer Review Rejected This", "a14-era -> now"),
    ]
    print(f"{'team':<16}{'window':<26}{'d cos floor':>12}{'d pds':>9}"
          f"{'leak bound':>11}{'x bound':>9}{'rho':>7}   verdict")
    out = []
    for team, m1, m2, note in cases:
        b, a = row(team, m1), row(team, m2)
        rb = {m: b["raw_" + m] for m in M}
        ra = {m: a["raw_" + m] for m in M}
        cl = C.classify(rb, ra)
        dmse = abs(ra["mse"] - rb["mse"])
        bound = p90 * dmse
        x = abs(cl["d_raw"]["pds"]) / bound if bound > 0 else float("inf")
        both = C.is_informative(rb["mse"]) and C.is_informative(ra["mse"])
        dcf = (C.cos_floor(ra["mse"]) - C.cos_floor(rb["mse"])) if both else float("nan")
        verdict = ("DIRECTION, amplitude excluded" if x > 3 else
                   "not resolvable against amplitude leakage")
        print(f"{team[:15]:<16}{(m1[:11]+'->'+m2[:11]):<26}{dcf:>+12.5f}"
              f"{cl['d_raw']['pds']:>+9.4f}{bound:>11.5f}{x:>9.1f}{cl['rho']:>7.2f}   {verdict}")
        out.append({"team": team, "from": m1, "to": m2, "note": note,
                    "d_cos_floor": dcf, "d_pds": cl["d_raw"]["pds"],
                    "leak_bound": bound, "x_bound": x, "rho": cl["rho"],
                    "both_informative": both, "verdict": verdict,
                    "cos_from": C.cos_floor(rb["mse"]) if C.is_informative(rb["mse"]) else None,
                    "cos_to": C.cos_floor(ra["mse"]) if C.is_informative(ra["mse"]) else None})

    json.dump({"leak_p90": p90, "leak_median": float(np.median(r)),
               "leak_calibrated": C.PDS_LEAK, "n_amp_arms": len(P), "cases": out},
              open(os.path.join(WORK, "valid215.json"), "w"), indent=1, default=float)
    print("\nwrote work/valid215.json")


if __name__ == "__main__":
    main()

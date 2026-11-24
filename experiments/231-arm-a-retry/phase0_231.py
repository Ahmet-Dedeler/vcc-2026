"""231 Phase 0 — reconcile 214 Arm A's `reach` claim against 220's outside bound.

220 bounded the ENTIRE ordering channel from outside: board raw `reach` ours 0.24983,
field best (luka/PertMind) 0.2908, so matching the best team in the world on `reach` is
worth 0.0410 * 1.12927 / 6 = +0.0077 overall.

219 says every point of Arm A is in `reach` (and `nmae`). So Arm A's `reach` claim can be
converted to an IMPLIED board raw `reach` and compared against the field maximum. If the
implied value exceeds the field best, the arm is refuted on arithmetic alone.

Known-answer rows carried here:
  KA0-a  220's quoted ours/Aginglab/field-best raw reach reproduce from the snapshot 220 read
  KA0-b  220's own conversion 0.0410 * 1.12927 / 6 = +0.0077
  KA0-c  the reach SPAN implied by our own board row (raw, scaled) equals 1/1.12927
  KA0-d  AGENTS 4.1c's board table (cqawesome 0.22431 / us) reproduces from the snapshot
Aborts on any failure.
"""

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)

FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<52} got {got:>14.7f}  exp {exp:>12.7f}  {note}")
    return ok


def load(path):
    d = json.load(open(path))
    return d["entries"] if isinstance(d, dict) and "entries" in d else d


def latest_by_team(entries):
    """Board rank uses a team's LATEST submission (AGENTS 8)."""
    best = {}
    for e in entries:
        t = e.get("team_name")
        d = e.get("submission_date") or ""
        if t not in best or d > best[t]["submission_date"]:
            best[t] = e
    return best


def main():
    snaps = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/2026-09-04*.json")))
    older = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/lb_*.json")))
    print(f"snapshots today: {[os.path.basename(s) for s in snaps]}")
    newest = snaps[-1]
    ent = load(newest)
    print(f"\nUsing {os.path.basename(newest)}  ({len(ent)} rows)")

    lat = latest_by_team(ent)
    rows = sorted(lat.values(), key=lambda e: -(e.get("score_avg") or -9))

    us = lat.get("Ahmet Dedeler")
    assert us is not None, "our team row not found"

    # --- our current shipped row -------------------------------------------------
    print("\n=== our current board row ===")
    for k in ("model_name", "submission_date", "score_avg", "de_wilcoxon_direction_reach_raw",
              "score_reach", "pds_cosine", "score_pds", "expr_mse_unbiased_capped_norm",
              "score_mse", "de_wilcoxon_lfc_nmae", "score_nmae",
              "de_wilcoxon_direction_fidelity_yield_raw", "score_fid",
              "de_wilcoxon_sig_jaccard", "score_jac"):
        print(f"  {k:<44} {us.get(k)}")
    our_rank = 1 + sum(1 for r in rows if (r.get("score_avg") or -9) > us["score_avg"])
    print(f"  rank (latest-per-team) = {our_rank} of {len(rows)}")

    # --- the reach field ---------------------------------------------------------
    reach_rows = [(r["team_name"], r["de_wilcoxon_direction_reach_raw"], r["score_reach"],
                   r.get("score_avg"))
                  for r in rows if r.get("de_wilcoxon_direction_reach_raw") is not None]
    reach_rows.sort(key=lambda x: -x[1])
    print("\n=== top 12 board raw `reach`, latest-per-team ===")
    for n, raw, sc, ov in reach_rows[:12]:
        print(f"  {n[:34]:<34} raw {raw:.6f}  scaled {sc:+.6f}  overall {ov:.6f}")
    field_best_name, field_best_raw, field_best_sc, _ = reach_rows[0]

    ours_raw = us["de_wilcoxon_direction_reach_raw"]
    ours_sc = us["score_reach"]
    our_reach_rank = 1 + sum(1 for x in reach_rows if x[1] > ours_raw)
    print(f"\n  ours     raw {ours_raw:.6f} scaled {ours_sc:+.6f}  (reach rank {our_reach_rank})")
    print(f"  field best {field_best_name}: raw {field_best_raw:.6f} scaled {field_best_sc:+.6f}")

    # span implied by our own row: scaled = (raw - b)/(r - b); the RATIO d_scaled/d_raw = 1/span
    # 220 used 1.12927 as the raw->scaled multiplier. Recover it from two board rows.
    print("\n=== Phase 0 known-answer gates ===")
    # KA0-c: recover the raw->scaled slope from the field, least squares over rows with both
    xs = [(r["de_wilcoxon_direction_reach_raw"], r["score_reach"]) for r in rows
          if r.get("de_wilcoxon_direction_reach_raw") is not None and r.get("score_reach") is not None]
    n = len(xs)
    mx = sum(a for a, _ in xs) / n
    my = sum(b for _, b in xs) / n
    num = sum((a - mx) * (b - my) for a, b in xs)
    den = sum((a - mx) ** 2 for a, _ in xs)
    slope = num / den
    intercept = my - slope * mx
    resid = max(abs(b - (slope * a + intercept)) for a, b in xs)
    print(f"  raw->scaled fit over {n} rows: slope {slope:.6f}  intercept {intercept:.6f}  "
          f"max resid {resid:.2e}")
    gate("KA0-c reach raw->scaled slope == 220's 1.12927", slope, 1.12927, 5e-4)

    # KA0-a: 220's quoted numbers. 220 wrote them on 2026-09-04; find the snapshot that carries
    # ours ~ 0.24983 and field best 0.2908, and record how far our own row has since moved.
    print("\n  --- KA0-a: locate the snapshot 220 quoted, and OUR OWN reach history ---")
    hit = None
    our_hist = []
    for s in snaps + older:
        e2 = load(s)
        l2 = latest_by_team(e2)
        u2 = l2.get("Ahmet Dedeler")
        if u2 is None or u2.get("de_wilcoxon_direction_reach_raw") is None:
            continue
        rr = sorted((x for x in l2.values() if x.get("de_wilcoxon_direction_reach_raw") is not None),
                    key=lambda x: -x["de_wilcoxon_direction_reach_raw"])
        ag = l2.get("Aginglab.com")
        print(f"    {os.path.basename(s):<28} ours raw {u2['de_wilcoxon_direction_reach_raw']:.5f} "
              f"scaled {u2['score_reach']:+.5f} `{str(u2.get('model_name'))[:26]:<26}` | "
              f"best {rr[0]['de_wilcoxon_direction_reach_raw']:.5f} ({rr[0]['team_name'][:14]})"
              + (f" aging {ag['de_wilcoxon_direction_reach_raw']:.5f}" if ag else ""))
        our_hist.append({"snapshot": os.path.basename(s),
                         "model_name": u2.get("model_name"),
                         "raw": u2["de_wilcoxon_direction_reach_raw"],
                         "scaled": u2["score_reach"],
                         "field_best_raw": rr[0]["de_wilcoxon_direction_reach_raw"],
                         "field_best_team": rr[0]["team_name"]})
        if abs(u2["de_wilcoxon_direction_reach_raw"] - 0.24983) < 1.0e-4:
            hit = (s, u2, rr)
    assert hit is not None, "could not locate the snapshot carrying 220's ours~0.24983"
    s220, u220, rr220 = hit
    print(f"  220's snapshot = {os.path.basename(s220)}  (model `{u220.get('model_name')}`)")
    gate("KA0-a1 220's ours raw reach", u220["de_wilcoxon_direction_reach_raw"], 0.24983, 1.0e-4,
         "220 rounded/read live; snapshot is 0.24988")
    gate("KA0-a2 220's field best raw reach", rr220[0]["de_wilcoxon_direction_reach_raw"],
         0.2908, 5e-4, f"({rr220[0]['team_name']})")

    # KA0-e: the EXACT raw->scaled slope on OUR OWN rows (same 3 contexts, same anchors), which
    # is the slope a delta on our own submission actually travels along.
    pts = sorted({(h["raw"], h["scaled"]) for h in our_hist})
    own_slopes = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            dr = pts[j][0] - pts[i][0]
            if abs(dr) > 1e-4:
                own_slopes.append((pts[j][1] - pts[i][1]) / dr)
    own_slope = sum(own_slopes) / len(own_slopes)
    print(f"  own-row raw->scaled slope over {len(own_slopes)} pairs: "
          f"{own_slope:.6f}  (min {min(own_slopes):.6f} max {max(own_slopes):.6f})")
    gate("KA0-e own-row slope == 220's 1.12927", own_slope, 1.12927, 2e-3)

    # KA0-b: 220's own conversion
    conv = (0.2908 - 0.24983) * 1.12927 / 6.0
    gate("KA0-b 220's 0.0410 x 1.12927 / 6", conv, 0.0077, 5e-5)

    # KA0-d: AGENTS 4.1c board table
    e313 = load(os.path.join(ROOT, "data/leaderboard/2026-09-04T0313Z.json"))
    l313 = latest_by_team(e313)
    gate("KA0-d1 4.1c cqawesome overall", l313["cqawesome"]["score_avg"], 0.22431, 5e-5)
    gate("KA0-d2 4.1c ours overall", l313["Ahmet Dedeler"]["score_avg"], 0.20492, 5e-5)

    if FAILURES:
        print("\nPHASE 0 FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # =========================================================================
    # THE RECONCILIATION
    # =========================================================================
    print("\n" + "=" * 78)
    print("PHASE 0 RECONCILIATION — Arm A's reach claim against 220's outside bound")
    print("=" * 78)

    SPAN_MULT = slope  # raw -> scaled

    # 214's arm A reach member, overall board units, at each transfer factor it reported.
    # 214 Phase 2.4: `reach` h1-test x 1.290 -> +0.007621 (wald89) / and the native pool row.
    # Recompute from the h1-test raw deltas 214 measured so the chain is explicit.
    h1_raw = {"native": 0.02605, "wald89": 0.03124}   # 214 sec 2.3 deltas (A - incumbent)
    # 214 sec 2.3 exact: native incumbent 0.19437 -> A 0.22042 = +0.02605
    #                    wald89 incumbent 0.18118 -> A 0.21242 = +0.03124
    deflator = 1.533   # 214 sec 1.5, measured
    factors = {"1.290 (rc_hybrid, order-changing)": 1.290,
               "0.4277 (164, spread-policy)": 0.4277,
               "1.000 (unity)": 1.000}

    out = {
        "snapshot": os.path.basename(newest),
        "our_row": {k: us.get(k) for k in
                    ("model_name", "submission_date", "score_avg",
                     "de_wilcoxon_direction_reach_raw", "score_reach")},
        "our_rank": our_rank,
        "our_reach_rank": our_reach_rank,
        "field_best": {"team": field_best_name, "raw": field_best_raw, "scaled": field_best_sc},
        "raw_to_scaled_slope": slope,
        "raw_to_scaled_intercept": intercept,
        "raw_to_scaled_max_resid": resid,
        "own_row_slope": own_slope,
        "our_reach_history": our_hist,
        "gap_to_field_best_raw": field_best_raw - ours_raw,
        "gap_to_field_best_overall": (field_best_raw - ours_raw) * slope / 6.0,
        "bound_220_gap_raw": 0.2908 - 0.24983,
        "bound_220_overall": (0.2908 - 0.24983) * 1.12927 / 6.0,
        "claims": [],
    }
    print(f"\n220's bound was written on raw reach {0.24983:.5f}; ours is now {ours_raw:.5f}. "
          f"The remaining ordering prize has shrunk from "
          f"{(0.2908 - 0.24983) * 1.12927 / 6.0:+.6f} to "
          f"{(field_best_raw - ours_raw) * slope / 6.0:+.6f} overall.")

    print(f"\nour CURRENT board reach : raw {ours_raw:.6f}  scaled {ours_sc:+.6f}  "
          f"(model `{us.get('model_name')}`, overall {us['score_avg']:.7f}, rank {our_rank})")
    print(f"field best              : raw {field_best_raw:.6f}  ({field_best_name})")
    print(f"gap to field best       : raw {field_best_raw - ours_raw:+.6f}  = "
          f"{(field_best_raw - ours_raw) * slope / 6.0:+.6f} overall")

    print(f"\n{'variant':<46} {'d_overall':>10} {'d_raw':>9} {'implied':>9} {'% of gap':>9} {'verdict':>10}")
    print("-" * 98)
    for pool, hr in h1_raw.items():
        for fname, f in factors.items():
            for defl, dname in ((1.0, "no deflator"), (deflator, "/1.533")):
                d_scaled_reach = hr * f / defl          # scaled board reach delta
                d_overall = d_scaled_reach / 6.0
                d_raw = d_scaled_reach / slope
                implied = ours_raw + d_raw
                pct = 100.0 * d_raw / (field_best_raw - ours_raw)
                verdict = "REFUTED" if implied > field_best_raw else "ok"
                label = f"{pool:<7} x{f:<6.4f} {dname}"
                print(f"{label:<46} {d_overall:>+10.6f} {d_raw:>+9.5f} {implied:>9.5f} "
                      f"{pct:>8.1f}% {verdict:>10}")
                out["claims"].append({
                    "pool": pool, "h1_raw_delta": hr, "transfer_factor": f,
                    "deflator": defl, "deflator_name": dname,
                    "d_scaled_reach": d_scaled_reach, "d_overall": d_overall,
                    "d_board_raw_reach": d_raw, "implied_board_raw_reach": implied,
                    "pct_of_field_gap": pct,
                    "exceeds_field_best": bool(implied > field_best_raw),
                })

    # where would each land us on the reach ladder
    print("\n=== where the headline variants would place us on the board's reach ladder ===")
    ladder = [(x[0], x[1]) for x in reach_rows]
    for c in out["claims"]:
        if c["transfer_factor"] in (1.290, 0.4277) and c["pool"] == "wald89":
            imp = c["implied_board_raw_reach"]
            rank = 1 + sum(1 for _, r in ladder if r > imp)
            c["implied_reach_rank"] = rank
            print(f"  wald89 x{c['transfer_factor']:.4f} {c['deflator_name']:<12} -> raw {imp:.5f}"
                  f"  reach rank {rank} (from {our_reach_rank})")

    json.dump(out, open(os.path.join(WORK, "phase0_231.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(WORK, 'phase0_231.json')}")


if __name__ == "__main__":
    main()

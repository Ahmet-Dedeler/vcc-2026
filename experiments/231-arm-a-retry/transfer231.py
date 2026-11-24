"""231 §1 — re-derive the h1-test -> board `reach` transfer factor on TODAY's instruments.

219 §3 row 5: Arm A needs "its `reach` transfer factor, which is the whole arm". 214 used
164's 0.4277 as a sensitivity and 185's 1.290 for its headline. Neither was re-derived after
`rc_callset`, `cd4f034` and `h1union` landed. This script derives it from every consecutive
single-variable board pair we own, with the h1-test delta measured on the SAME instrument.

Known-answer rows (abort on failure):
  KA1  185's 1.290 reproduces from its own stated inputs (rc_hybrid 0.21681 - a15 twin 0.17910;
       board +0.04865)
  KA2  the ledger's a8/a15 0.390 reproduces from the board rows
  KA3  every board raw `reach` used here matches the live snapshot for that entry id
  KA4  the ledger's scaling check: a8 raw 0.1890 -> scaled 0.124349 at slope 1.1292
"""

import glob
import json
import math
import os
import statistics
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
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<50} got {got:>12.6f} exp {exp:>11.6f}  {note}")


# --- board rows, by entry id (the id is the primary key; the model name is a joke string) ----
# arm identification from docs/BOARD-LEDGER.md sec 3b + the *BOARD-RESULT.md files.
ARMS = {
    "a8":         dict(entry="lfnx5c9vp1gVhlNWcemT", overall=0.1543607),
    "a14":        dict(entry="0NvkHnMPzLFaJl8phAxS", overall=0.1704242),
    "a23":        dict(entry="IzQVfDpeDWufkPT8Lmd2", overall=0.1733597),
    "abc-knee":   dict(entry="N3ui4hCkKx9bxExG2pxE", overall=0.1740437),
    "rc_hybrid":  dict(entry="nQUJZhZb3MW2uFJHBX5e", overall=0.1834190),
    "rc_callset": dict(entry="31oB9yWy8CSNBDEkY4iW", overall=0.1918638),
    "cd4f034":    dict(entry="QSe3JlNUYZ2BzvlcsL14", overall=0.2049196),
    "h1union":    dict(entry="iUziqE1gewfcPKzYUDaB", overall=0.2118310),
}
# a15 never appears in a saved snapshot (superseded before the next pull); its members come from
# docs/BOARD-LEDGER.md's raw table, which was built from experiments/096-mse-floor/A15-BOARD-RESULT.md
A15_BOARD_RAW_REACH = 0.181173
A15_BOARD_SCALED_REACH = 0.115368

# --- h1-test DELIVERED `reach`, same surface (85 targets, native BH pool), same statistic -----
# sources, all quoted with their provenance:
H1 = {
    "a15":        dict(reach=0.18283, src="172 / 220 sec0 KA (`reach(production order)`)"),
    "a15_twin":   dict(reach=0.17910, src="185's re-seeded twin MEAN, the baseline it used"),
    "a8":         dict(reach=0.20291, src="docs/BOARD-LEDGER.md sec4 (a8-flat)"),
    "rc_hybrid":  dict(reach=0.21681, src="185 BOARD-RESULT-rc_hybrid.md; BOARD-LEDGER sec5"),
    "rc_callset": dict(reach=0.223341, src="214 / 220 sec0 KA, bit-identical to de202"),
    "cd4f034":    dict(reach=0.223341, src="214: h1-test effect matrices IDENTICAL to rc_callset, "
                                           "max|diff| 0.000e+00 -- the dose is invisible here"),
}

SLOPE = 1.1292  # board raw -> scaled for reach (ledger, fitted on 612 live entries)


def load_board():
    seen = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json"))):
        d = json.load(open(f))
        e = d["entries"] if isinstance(d, dict) and "entries" in d else d
        for r in e:
            if r.get("team_name") == "Ahmet Dedeler":
                seen[r["id"]] = r
    return seen


def main():
    board = load_board()
    print("=== Phase 1 known-answer gates ===")

    for name, a in ARMS.items():
        r = board.get(a["entry"])
        assert r is not None, f"entry {a['entry']} ({name}) not in any snapshot"
        gate(f"KA3 {name} overall", r["score_avg"], a["overall"], 5e-6, f"[{a['entry']}]")
        a["reach_raw"] = r["de_wilcoxon_direction_reach_raw"]
        a["reach_scaled"] = r["score_reach"]
        a["pds_raw"] = r["pds_cosine"]
        a["mse_raw"] = r["expr_mse_unbiased_capped_norm"]
        a["date"] = r["submission_date"][:19]

    ARMS["a15"] = dict(entry="(superseded, from A15-BOARD-RESULT via ledger)",
                       overall=0.1740865, reach_raw=A15_BOARD_RAW_REACH,
                       reach_scaled=A15_BOARD_SCALED_REACH, date="2026-08-30")

    # KA4 the ledger's own scaling check
    gate("KA4 a8 raw->scaled at slope 1.1292",
         (ARMS["a8"]["reach_raw"] - 0.0790) * SLOPE, 0.124349, 3e-4)

    # KA1 reproduce 185's 1.290 from its own stated inputs
    ka1_h1 = H1["rc_hybrid"]["reach"] - H1["a15_twin"]["reach"]
    ka1_board = ARMS["rc_hybrid"]["reach_raw"] - A15_BOARD_RAW_REACH
    gate("KA1a 185's h1-test delta", ka1_h1, 0.03771, 5e-5)
    gate("KA1b 185's board delta", ka1_board, 0.04865, 1.5e-4, "185 rounded a15 raw reach")
    gate("KA1c 185's transfer factor 1.290", ka1_board / ka1_h1, 1.290, 5e-3)

    # KA2 reproduce the ledger's a8/a15 0.390
    ka2_h1 = H1["a8"]["reach"] - H1["a15"]["reach"]
    ka2_board = ARMS["a8"]["reach_raw"] - A15_BOARD_RAW_REACH
    gate("KA2a ledger a8/a15 h1-test delta", ka2_h1, 0.02008, 5e-5)
    gate("KA2b ledger a8/a15 board delta", ka2_board, 0.007827, 5e-5)
    gate("KA2c ledger a8/a15 factor 0.390", ka2_board / ka2_h1, 0.390, 3e-3)

    if FAILURES:
        print("\nPHASE 1 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # =========================================================================
    print("\n" + "=" * 96)
    print("THE BOARD `reach` TRAJECTORY — every entry, and the h1-test value paired with it")
    print("=" * 96)
    order = ["a8", "a14", "a23", "abc-knee", "a15", "rc_hybrid", "rc_callset", "cd4f034", "h1union"]
    print(f"{'arm':<12} {'date':<12} {'overall':>10} {'board raw':>10} {'board scal':>11} {'h1-test':>10}  class")
    CLS = {"a8": "emission (spread policy)", "a14": "profile (per-context cap)",
           "a23": "renderer (mass fix)", "abc-knee": "profile (per-context alpha)",
           "a15": "profile (CD4 residual)", "rc_hybrid": "emission-law, profile-PRESERVING",
           "rc_callset": "emission-law, profile-PRESERVING (call set)",
           "cd4f034": "PROFILE (CD4 dose 0.10->0.34)",
           "h1union": "PROFILE (new donor source, g=0.28)"}
    for k in order:
        a = ARMS[k]
        h = H1.get(k, {}).get("reach")
        hs = f"{h:.5f}" if h is not None else "     --"
        print(f"{k:<12} {a['date'][:10]:<12} {a['overall']:>10.7f} {a['reach_raw']:>10.6f} "
              f"{a['reach_scaled']:>+11.6f} {hs:>10}  {CLS[k]}")

    # =========================================================================
    print("\n" + "=" * 96)
    print("THE TRANSFER FACTOR, on every pair where BOTH sides are measured")
    print("=" * 96)
    PAIRS = [
        ("a8", "a15", "emission / spread policy", "the ledger sec4 pair (confounded: alpha, CD4, cap also move)"),
        ("rc_hybrid", "a15_twin", "emission-law, profile-PRESERVING", "185's own pair -- the source of 1.290"),
        ("rc_hybrid", "a15", "emission-law, profile-PRESERVING", "same pair, a15's SINGLE run instead of the twin mean"),
        ("rc_callset", "rc_hybrid", "emission-law, profile-PRESERVING", "the call-set arm, single-variable"),
        ("cd4f034", "rc_callset", "PROFILE (donor blend / dose)", "**the class Arm A belongs to**"),
    ]
    rows = []
    print(f"{'pair':<26} {'h1 delta':>10} {'board delta':>12} {'factor':>10}   class / note")
    print("-" * 110)
    for hi, lo, cls, note in PAIRS:
        h_hi = H1[hi]["reach"]
        h_lo = H1[lo]["reach"]
        b_hi = ARMS[hi]["reach_raw"] if hi in ARMS else None
        b_lo = ARMS[lo]["reach_raw"] if lo in ARMS else A15_BOARD_RAW_REACH
        dh = h_hi - h_lo
        db = b_hi - b_lo
        f = (db / dh) if abs(dh) > 1e-9 else float("inf")
        fs = f"{f:>10.4f}" if math.isfinite(f) else f"{'INFINITE':>10}"
        print(f"{hi + ' - ' + lo:<26} {dh:>+10.5f} {db:>+12.6f} {fs}   {cls}")
        print(f"{'':<26} {'':<10} {'':<12} {'':<10}   {note}")
        rows.append(dict(pair=f"{hi}-{lo}", h1_delta=dh, board_delta=db,
                         factor=(f if math.isfinite(f) else None), arm_class=cls, note=note))

    # the cd4f034 pair, spelled out
    dh0 = H1["cd4f034"]["reach"] - H1["rc_callset"]["reach"]
    db0 = ARMS["cd4f034"]["reach_raw"] - ARMS["rc_callset"]["reach_raw"]
    print("\n" + "!" * 96)
    print(f"cd4f034 - rc_callset: h1-test delta = {dh0:+.9f}  (EXACTLY ZERO -- 214 measured the two")
    print(f"  h1-test effect matrices identical, max|diff| 0.000e+00), while the board delivered")
    print(f"  {db0:+.6f} raw reach = {db0 * SLOPE / 6.0:+.6f} overall.")
    print(f"  The transfer factor for the PROFILE class is therefore UNDEFINED (division by zero),")
    print(f"  not large or small. h1-test cannot see this arm class's `reach` at all.")
    print("!" * 96)

    # h1union: board move with no h1-test read
    dbh = ARMS["h1union"]["reach_raw"] - ARMS["cd4f034"]["reach_raw"]
    print(f"\nh1union - cd4f034: board {dbh:+.6f} raw reach = {dbh * SLOPE / 6.0:+.6f} overall; "
          f"h1-test UNMEASURABLE\n  (h1_2025 is 85/85 oracle leakage on h1-test -- AGENTS sec6.4 trap 1 -- "
          f"so the H1-union\n  source cannot be scored there at all).")

    # --- the error bar ------------------------------------------------------
    fin = [r["factor"] for r in rows if r["factor"] is not None]
    print("\n" + "=" * 96)
    print("THE ERROR BAR")
    print("=" * 96)
    print(f"  measured finite factors: {[round(x, 4) for x in fin]}")
    lg = [math.log(x) for x in fin]
    gm = math.exp(statistics.mean(lg))
    gsd = math.exp(statistics.stdev(lg))
    print(f"  geometric mean {gm:.4f}, geometric sd x/{gsd:.3f}  -> 1-sigma band "
          f"[{gm / gsd:.4f}, {gm * gsd:.4f}]")
    print(f"  range {min(fin):.4f} .. {max(fin):.4f}  = a spread of {max(fin) / min(fin):.2f}x, "
          f"plus one UNDEFINED point")
    print(f"  164's 0.4277 and 185's 1.290 are BOTH inside this set and {1.290 / 0.4277:.2f}x apart.")

    # class-restricted: what is the honest factor for a PROFILE arm?
    print("\n  Restricted to the arm class Arm A belongs to (PROFILE / donor-blend):")
    print(f"    the only measured pair is cd4f034-rc_callset, and it is 0/0-undefined.")
    print(f"    n = 0 usable pairs. **There is no measured transfer factor for this arm class.**")

    out = dict(arms={k: {kk: vv for kk, vv in v.items()} for k, v in ARMS.items()},
               h1=H1, pairs=rows, slope=SLOPE,
               finite_factors=fin, geometric_mean=gm, geometric_sd=gsd,
               profile_class_pairs_usable=0,
               cd4f034_h1_delta=dh0, cd4f034_board_delta=db0,
               h1union_board_delta=dbh)
    json.dump(out, open(os.path.join(WORK, "transfer231.json"), "w"), indent=2, default=str)
    print(f"\nwrote {os.path.join(WORK, 'transfer231.json')}")


if __name__ == "__main__":
    main()

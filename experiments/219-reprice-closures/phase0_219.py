"""219 Phase 0 — re-derive AGENTS.md 8.4's correction from the artifacts, not from its prose,
and put a RANGE on it.

Everything here is arithmetic on numbers already on disk.  No new measurement.

Outputs work/phase0_219.json and prints a compact report.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

out: dict = {}
fail: list[str] = []


def gate(name, expected, measured, tol=1e-4, rel=True):
    e, m = float(expected), float(measured)
    d = abs(m / e - 1.0) if (rel and e) else abs(m - e)
    ok = d <= tol
    out.setdefault("gates", []).append(
        dict(name=name, expected=e, measured=m, dev=d, rel=rel, pass_=bool(ok)))
    if not ok:
        fail.append(f"{name}: expected {e!r} got {m!r} (dev {d:.3e})")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<62s} exp {e:>14.8g}  got {m:>14.8g}")
    return ok


# =======================================================================================
# 0.  The board rows this audit uses.  docs/BOARD-LEDGER.md + 210's BOARD-RESULT.md.
# =======================================================================================
SCALED = {                     # scaled members, which average to Overall
    "a14":         dict(pds=0.738614, mse=0.110615, nmae=0.091972, fid=-0.032586,
                        reach=0.114542, jac=-0.000612, overall=0.1704242),
    "a15":         dict(pds=0.745871, mse=0.120057, nmae=0.094239, fid=-0.030397,
                        reach=0.115368, jac=-0.000619, overall=0.1740865),
    "rc_callset":  dict(pds=0.7458708290005253, mse=0.12005661519665856,
                        nmae=0.09429042916917090, fid=-0.029268,
                        reach=0.192898, jac=0.027335, overall=0.1918638367),
    "cd4f034":     dict(pds=0.737800, mse=0.159276, nmae=0.104861, fid=-0.020127,
                        reach=0.219300, jac=0.028408, overall=0.2049195830),
}
RAW = {   # raw members, where the ledger records them
    "a14": dict(mse=0.882775), "a15": dict(mse=0.873805),
    "rc_callset": dict(mse=0.8738045),      # identical to a15/rc_hybrid to 16 s.f. on scaled
    "cd4f034": dict(mse=0.83635),
}
MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")

print("\n=== KA block A — the board arithmetic (181 gate 1's construction) ===")
for arm, row in SCALED.items():
    avg = sum(row[m] for m in MEMBERS) / 6.0
    gate(f"A1 {arm}: six scaled members average to Overall", row["overall"], avg, tol=2e-4)

# the mse scaling map, solved from the a14/a15 pair, checked out of sample on abc-knee
slope_mse = (SCALED["a15"]["mse"] - SCALED["a14"]["mse"]) / (RAW["a14"]["mse"] - RAW["a15"]["mse"])
b_mse = RAW["a15"]["mse"] + SCALED["a15"]["mse"] / slope_mse
gate("A2 mse slope from a14/a15 vs the ledger's published −1.051", 1.051, slope_mse, tol=2e-3)
gate("A3 abc-knee scaled mse out of sample (recorded 0.119322)",
     0.119322, slope_mse * (b_mse - 0.874491), tol=2e-3)
out["mse_map"] = dict(slope=slope_mse, baseline=b_mse)

# rc_callset / cd4f034 raw mse implied by the map, against the transcribed values
raw_rc_from_map = b_mse - SCALED["rc_callset"]["mse"] / slope_mse
raw_cd_from_map = b_mse - SCALED["cd4f034"]["mse"] / slope_mse
out["raw_mse_from_map"] = dict(rc_callset=raw_rc_from_map, cd4f034=raw_cd_from_map,
                               rc_transcribed=0.87367, cd_transcribed=RAW["cd4f034"]["mse"])
print(f"  raw mse implied by the map: rc_callset {raw_rc_from_map:.6f} "
      f"(210's table transcribes 0.87367), cd4f034 {raw_cd_from_map:.6f} "
      f"(210's table transcribes {RAW['cd4f034']['mse']})")

# =======================================================================================
# 1.  197's inversion, re-derived, and the cosine ratio the board actually delivered
# =======================================================================================
print("\n=== Phase 0.1 — 197's one-sided inversion, re-derived from its own constants ===")
Kc, SUMD, K = 0.005410, 0.229141, 0.246529
NUM_C = 1.0 + Kc / SUMD
DEN_C = 2.0 * np.sqrt(K / SUMD)
gate("B1 197's numerator constant 1.023611", 1.023611, NUM_C, tol=1e-5)
gate("B2 197's denominator constant 2.074495", 2.074495, DEN_C, tol=1e-5)


def cos_floor(raw_mse: float) -> float:
    return (NUM_C - raw_mse) / DEN_C


# 197 §2.2's tightness table, verbatim, is the only calibration of floor -> truth
TIGHT = {"a11": (0.90999, 0.05477, 0.07550), "a12": (0.90135, 0.05894, 0.07420),
         "a14": (0.88277, 0.06789, 0.07418), "a15": (0.87381, 0.07221, 0.07925),
         "a23": (0.87235, 0.07291, 0.07785), "abc-knee": (0.87449, 0.07188, 0.07452)}
for name, (raw, floor, meas) in TIGHT.items():
    gate(f"B3 197 floor reproduces for {name}", floor, cos_floor(raw), tol=3e-3)
near_cliff = [floor / meas for n, (r, floor, meas) in TIGHT.items() if n in
              ("a14", "a15", "a23", "abc-knee")]
out["floor_tightness_near_cliff"] = dict(values=near_cliff,
                                         lo=min(near_cliff), hi=max(near_cliff))
print(f"  floor/measured on the four near-cliff rows: "
      f"{min(near_cliff):.3f}–{max(near_cliff):.3f}")

f_rc_t, f_cd = cos_floor(0.87367), cos_floor(RAW["cd4f034"]["mse"])
f_rc_m, f_cd_m = cos_floor(raw_rc_from_map), cos_floor(raw_cd_from_map)
gate("B4 §8.4's 0.07228 -> 0.09027 pair", 0.07228, f_rc_t, tol=1e-3)
gate("B5 §8.4's cd4f034 floor 0.09027", 0.09027, f_cd, tol=1e-3)
ratio_transcribed = f_cd / f_rc_t
ratio_from_map = f_cd_m / f_rc_m
gate("B6 210's '+24.9% board cosine'", 1.249, ratio_transcribed, tol=3e-3)
out["board_cos_ratio_cd4f034"] = dict(from_transcribed_raw=ratio_transcribed,
                                      from_scaling_map=ratio_from_map)
print(f"  board cosine ratio cd4f034/rc_callset: {ratio_transcribed:.4f} "
      f"(transcribed raw) / {ratio_from_map:.4f} (via the scaling map)")

# the same ratio for the a14 -> a15 dose step, both estimators
a1415_floor_ratio = cos_floor(RAW["a15"]["mse"]) / cos_floor(RAW["a14"]["mse"])
a1415_meas_ratio = TIGHT["a15"][2] / TIGHT["a14"][2]
out["board_cos_ratio_a14_a15"] = dict(floor=a1415_floor_ratio, measured=a1415_meas_ratio)
gate("B7 197 §1.2's a14->a15 +6.8% (measured cosines)", 1.0683, a1415_meas_ratio, tol=3e-3)
print(f"  a14->a15 board cosine ratio: floors {a1415_floor_ratio:.4f}, "
      f"197's measured {a1415_meas_ratio:.4f}  -> the two estimators agree to "
      f"{abs(a1415_floor_ratio - a1415_meas_ratio)*100:.2f} pp")

# =======================================================================================
# 2.  The PLOLO forward measurements, in BOTH cosine functionals
# =======================================================================================
print("\n=== Phase 0.2 — the PLOLO forward cosine, in both functionals ===")
rows212 = json.loads((ROOT / "experiments/212-gamma-cosine/work/cd4f212_summary.json").read_text())
g65 = [r for r in rows212 if abs(r["gamma"] - 0.65) < 1e-12][0]
macro_010, macro_034 = g65["cos_f010"], g65["cos_f034"]
gate("C1 212's macro cosine at f=0.10, gamma 0.65", 0.0218418, macro_010, tol=1e-3)
gate("C2 212's macro cosine at f=0.34", 0.0248622, macro_034, tol=1e-3)
gate("C3 212 §4.3's +13.83%", 1.1383, macro_034 / macro_010, tol=2e-3)

sw = json.loads((ROOT / "experiments/207-cd4-dose/work/sweep207.json").read_text())


def agg(arm: str) -> float:
    v = [r["cos_agg"] for r in sw if r["arm"] == arm]
    assert len(v) == 4, (arm, len(v))
    return float(np.mean(v))


agg_f0, agg_a14 = agg("f0"), agg("a14")
agg_010, agg_034 = agg("f0.1"), agg("f0.34")
gate("C4 207's cos_agg at f=0.10 (its own table: 0.020086)", 0.020086, agg_010, tol=1e-3)
gate("C5 207's cos_agg at f=0.34 (its own table: 0.023698)", 0.023698, agg_034, tol=1e-3)
gate("C6 207's a14 cos_agg (its k_cos derivation used 0.019053)", 0.019053, agg_a14, tol=2e-3)
out["plolo_cos"] = dict(macro_f010=macro_010, macro_f034=macro_034,
                        agg_f0=agg_f0, agg_a14=agg_a14, agg_f010=agg_010, agg_f034=agg_034)

k_cos_a1415_meas = a1415_meas_ratio_pct = (a1415_meas_ratio - 1) / (agg_010 / agg_a14 - 1)
k_cos_a1415_floor = (a1415_floor_ratio - 1) / (agg_010 / agg_a14 - 1)
k_cos_cd_agg = (ratio_transcribed - 1) / (agg_034 / agg_010 - 1)
k_cos_cd_macro = (ratio_transcribed - 1) / (macro_034 / macro_010 - 1)
gate("C7 207's own published k_cos = 1.2598 (a14->a15, cos_agg, measured cosines)",
     1.2598, k_cos_a1415_meas, tol=2e-2)
out["k_cos"] = dict(a14_a15_measured_cos=k_cos_a1415_meas,
                    a14_a15_floor_cos=k_cos_a1415_floor,
                    cd4f034_cos_agg=k_cos_cd_agg,
                    cd4f034_macro_212=k_cos_cd_macro)
print(f"\n  k_cos (board cosine gain / PLOLO cosine gain), FOUR readings:")
print(f"    a14->a15, 207's cos_agg, 197 measured cosines : {k_cos_a1415_meas:.3f}")
print(f"    a14->a15, 207's cos_agg, floor-ratio          : {k_cos_a1415_floor:.3f}")
print(f"    cd4f034,  207's cos_agg (197's functional)    : {k_cos_cd_agg:.3f}")
print(f"    cd4f034,  212's macro-per-fold  <- §8.4's 1.80: {k_cos_cd_macro:.3f}")

# =======================================================================================
# 3.  The propagation factor:  overall / mse-only, on BOTH board points
# =======================================================================================
print("\n=== Phase 0.3 — how a direction gain propagates, on TWO board points ===")


def deltas(a, b):
    d = {m: SCALED[b][m] - SCALED[a][m] for m in MEMBERS}
    d["overall_scaled"] = sum(d.values())
    d["overall"] = d["overall_scaled"] / 6.0
    return d


d_a1415 = deltas("a14", "a15")
d_cd = deltas("rc_callset", "cd4f034")
gate("D1 a14->a15 overall reconstructs from its six members",
     SCALED["a15"]["overall"] - SCALED["a14"]["overall"], d_a1415["overall"], tol=2e-3)
gate("D2 rc_callset->cd4f034 overall = 210's +0.0130557",
     0.0130557, d_cd["overall"], tol=2e-3)
gate("D3 210's realised mse leg +0.006537 overall", 0.006537, d_cd["mse"] / 6, tol=2e-3)

prop = {}
for name, d, dcos in (("a14->a15 (f 0.053->0.10)", d_a1415,
                       cos_floor(RAW["a15"]["mse"]) - cos_floor(RAW["a14"]["mse"])),
                      ("rc_callset->cd4f034 (f 0.10->0.34)", d_cd, f_cd - f_rc_t)):
    mse_ovr = d["mse"] / 6
    de4 = sum(d[m] for m in ("nmae", "fid", "reach", "jac")) / 6
    pds = d["pds"] / 6
    prop[name] = dict(
        d_cos_floor=dcos, mse_ovr=mse_ovr, pds_ovr=pds, de4_ovr=de4,
        total_ovr=d["overall"], P_total_over_mse=d["overall"] / mse_ovr,
        rate_mse_per_cos=mse_ovr / dcos, rate_pds_per_cos=pds / dcos,
        rate_de4_per_cos=de4 / dcos, rate_total_per_cos=d["overall"] / dcos)
    p = prop[name]
    print(f"\n  {name}")
    print(f"    Δcos (floor)  {dcos:+.5f}      total {d['overall']:+.6f}")
    print(f"    mse {mse_ovr:+.6f}   pds {pds:+.6f}   DE-four {de4:+.6f}")
    print(f"    P = total/mse = {p['P_total_over_mse']:.3f}")
    print(f"    per unit cosine: mse {p['rate_mse_per_cos']:+.4f}  "
          f"pds {p['rate_pds_per_cos']:+.4f}  DE4 {p['rate_de4_per_cos']:+.4f}  "
          f"TOTAL {p['rate_total_per_cos']:+.4f}")
out["propagation"] = prop
gate("D4 the mse exchange rate reproduces 8.2/212's +0.3655 per unit cosine, point 1",
     0.3655, prop["a14->a15 (f 0.053->0.10)"]["rate_mse_per_cos"], tol=0.03)
gate("D5 the same rate on the second board point",
     0.3655, prop["rc_callset->cd4f034 (f 0.10->0.34)"]["rate_mse_per_cos"], tol=0.03)

# =======================================================================================
# 4.  §8.4's own multipliers, re-derived, and their decomposition
# =======================================================================================
print("\n=== Phase 0.4 — §8.4's 2.51 / 1.80 / 4.5, re-derived and decomposed ===")
PRED = dict(mse=+0.005206, pds=-0.000685, nmae=-0.000152, fid=0.0, reach=0.0, jac=0.0)
PRED["total"] = sum(PRED[m] for m in MEMBERS)
gate("E1 207/210's preregistered total +0.00437", 0.004369, PRED["total"], tol=2e-3)
M_84 = d_cd["overall"] / PRED["mse"]
gate("E2 §8.4's 2.51x = realised total / PREDICTED mse", 2.51, M_84, tol=5e-3)
E_mse = (d_cd["mse"] / 6) / PRED["mse"]
P_cd = prop["rc_callset->cd4f034 (f 0.10->0.34)"]["P_total_over_mse"]
print(f"  2.51 decomposes EXACTLY as  P x E  =  {P_cd:.4f} x {E_mse:.4f} = {P_cd*E_mse:.4f}")
gate("E3 the decomposition is exact", M_84, P_cd * E_mse, tol=1e-6)
out["multiplier_84"] = dict(M=M_84, P=P_cd, E=E_mse,
                            P_other_board_point=prop["a14->a15 (f 0.053->0.10)"]["P_total_over_mse"])

# the double-count check
print(f"\n  §8.4 states 4.5x = 2.51 x 1.80.  But 207's price ALREADY carried a PLOLO->board")
print(f"  cosine correction of k_cos = 1.2598, and E = {E_mse:.3f} is what is LEFT of the")
print(f"  transfer error after it.  Implied total transfer factor = 1.2598 x {E_mse:.3f} "
      f"= {1.2598*E_mse:.3f},")
print(f"  against 207's cos_agg reading of {k_cos_cd_agg:.3f} — agreement "
      f"{abs(1.2598*E_mse/k_cos_cd_agg-1)*100:.0f}%.")
out["double_count"] = dict(k_cos_applied_by_207=1.2598, residual_E=E_mse,
                           implied_total_transfer=1.2598 * E_mse,
                           direct_cos_agg_reading=k_cos_cd_agg)

# =======================================================================================
# 5.  The counter-rate: operator-geometry arms.  212 §3.4, re-expressed with DE on.
# =======================================================================================
print("\n=== Phase 0.5 — the same rule applied to OPERATOR-GEOMETRY arms ===")
R_MSE, R_PDS = 0.3655, -0.3517               # 212 §3.4, gamma axis, measured
de4_lo = prop["a14->a15 (f 0.053->0.10)"]["rate_de4_per_cos"]
de4_hi = prop["rc_callset->cd4f034 (f 0.10->0.34)"]["rate_de4_per_cos"]
M_geom = [(R_MSE + R_PDS + d) / R_MSE for d in (de4_lo, de4_hi)]
out["operator_geometry_multiplier"] = dict(rate_mse=R_MSE, rate_pds=R_PDS,
                                           de4_lo=de4_lo, de4_hi=de4_hi,
                                           M_lo=min(M_geom), M_hi=max(M_geom))
print(f"  d(mse)/d(cos) {R_MSE:+.4f}, d(pds)/d(cos) {R_PDS:+.4f} (212 §3.4, gamma axis)")
print(f"  d(DE4)/d(cos) measured on the two board points: {de4_lo:+.4f} .. {de4_hi:+.4f}")
print(f"  => multiplier for an operator-geometry arm = {min(M_geom):.2f} .. {max(M_geom):.2f}"
      f"   (against {M_84:.2f} for a geometry-preserving one)")

# =======================================================================================
# 6.  KA block B — 198's own reconstruction rows, on this audit's arithmetic
# =======================================================================================
print("\n=== KA block B — 198's known-answer rows (150 / 151 six-member reconstruction) ===")
E150 = dict(pds=0.0, mse=0.0, nmae=+0.006929, fid=-0.006765, reach=-0.045128, jac=+0.014694)
E151 = dict(pds=0.0, mse=+0.000341, nmae=+0.006315, fid=-0.003736, reach=-0.043628,
            jac=+0.014645)
gate("F1 150's six-member average (198 R1)", -0.005045, sum(E150.values()) / 6, tol=1e-3)
gate("F2 151's six-member average (198 R2)", -0.004344, sum(E151.values()) / 6, tol=1e-3)
gate("F3 150 ex-reach (198: +0.002476 overall)", 0.002476,
     sum(v for k, v in E150.items() if k != "reach") / 6, tol=2e-3)

print("\n=== KA block C — 181's board arithmetic on the snapshot it used ===")
lb = json.loads((ROOT / "data/leaderboard/2026-09-02T0630Z.json").read_text())
entries = lb if isinstance(lb, list) else lb.get("entries", lb.get("data", []))
key = {}
if entries:
    e0 = entries[0]
    key = {k: k for k in e0}
bad, n = 0, 0
for e in entries:
    try:
        vals = [float(e[f"score_{m}"]) for m in
                ("pds", "mse", "nmae", "fid", "reach", "jac")]
        tot = float(e["score_avg"])
    except Exception:
        continue
    n += 1
    if abs(np.mean(vals) - tot) > 1e-6:
        bad += 1
if n:
    gate(f"G1 six scaled members average to score_avg on all {n} board entries",
         0.0, bad, tol=0.5, rel=False)
else:
    print("  [skip] board snapshot schema not matched; keys:", sorted(key)[:14])
    out["board_snapshot_keys"] = sorted(key)

print("\n=== KA block D — PLOLO, from 212's cached fold scores (gate-10 verified vs 204) ===")
plolo = json.loads((ROOT / "experiments/212-gamma-cosine/work/plolo212.json").read_text())
p65 = {r["fold"]: r for r in plolo if abs(r["gamma"] - 0.65) < 1e-12 and r["arm"] == "equal"}
gate("H1 a15 pds, k562_gwps fold, equal weights (204 row 2)",
     0.79578409, p65["k562_gwps"]["a15_pds"], tol=1e-7)
gate("H2 a23 pds, same fold (204 row 2b)", 0.78912131, p65["k562_gwps"]["a23_pds"], tol=1e-7)
SPAN_PDS = 0.452
massfix = float(np.mean([p65[f]["a23_pds"] - p65[f]["a15_pds"] for f in p65])) / SPAN_PDS
gate("H3 a15->a23 mass-fix price, 4 folds, scaled (162/204 row 3)", 0.004014, massfix, tol=1e-3)
out["plolo_ka"] = dict(a15_pds=p65["k562_gwps"]["a15_pds"],
                       a23_pds=p65["k562_gwps"]["a23_pds"], massfix=massfix,
                       provenance="the k562_gwps leg is REBUILT FROM SCRATCH by "
                                  "ka_plolo219.py (max_abs_dev_from_162_pre 4.13e-08); the "
                                  "other three folds enter only the 4-fold mass-fix mean and "
                                  "come from 212 work/plolo212.json, whose gate 10 reproduced "
                                  "204's 11 published scalars at rel 0.00e+00.")

# =======================================================================================
print("\n=== SUMMARY: the multiplier, with its range ===")
band = dict(
    P_lo=min(p["P_total_over_mse"] for p in prop.values()),
    P_hi=max(p["P_total_over_mse"] for p in prop.values()),
    E_point=E_mse,
    M_point=M_84,
)
band["M_lo"] = band["P_lo"] * 1.0
band["M_hi"] = band["P_hi"] * (k_cos_cd_agg / 1.2598)
out["multiplier_band"] = band
print(f"  propagation P            : {band['P_lo']:.2f} .. {band['P_hi']:.2f}  (2 board points)")
print(f"  transfer residual E      : {band['E_point']:.2f}   (1 board point)")
print(f"  END-TO-END M (geometry-preserving cosine arm, priced 207's way):")
print(f"      point {band['M_point']:.2f},  range {band['M_lo']:.2f} .. {band['M_hi']:.2f}")
print(f"  END-TO-END M (operator-geometry arm): {min(M_geom):.2f} .. {max(M_geom):.2f}")
print(f"  §8.4's headline 4.5x is NOT reproduced: it multiplies 2.51 (which already "
      f"contains\n      the transfer error) by 1.80 (the transfer error, on a mismatched "
      f"cosine functional).")

out["failures"] = fail
(WORK / "phase0_219.json").write_text(json.dumps(out, indent=1, default=float))
print(f"\n{len(out['gates'])} gates, {len(fail)} FAIL")
for f in fail:
    print("  FAIL:", f)

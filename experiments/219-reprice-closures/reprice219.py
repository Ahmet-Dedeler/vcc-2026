"""219 Phase 2 — re-price the Class-A arms with the DE members ON, using the rates Phase 0
measured on the two board points rather than 8.4's prose multipliers.

Every input is transcribed from a named file and line; nothing is re-measured here.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
P0 = json.loads((WORK / "phase0_219.json").read_text())

COS_NOW = 0.090267669              # board cosine floor after cd4f034 (Phase 0 B5)
RATE_MSE = 0.3655                  # 8.2 / 212 §3.4, reproduced on both board points (D4/D5)
prop = P0["propagation"]
RATE_TOT_LO = min(p["rate_total_per_cos"] for p in prop.values())
RATE_TOT_HI = max(p["rate_total_per_cos"] for p in prop.values())
FLOOR_SHAPE, FLOOR_AMP = 0.0019, 0.0089

out: dict = {}
print("=== 1. The §8.2 single-arm bar, re-derived on the total rate ===")
print(f"  total overall per unit board cosine, measured on the two board points: "
      f"{RATE_TOT_LO:.4f} .. {RATE_TOT_HI:.4f}")
print(f"  mse-only rate §8.2 used                                             : {RATE_MSE:.4f}")
bar = {}
for name, floor in (("shape-only +0.0019", FLOOR_SHAPE), ("amplitude +0.0089", FLOOR_AMP)):
    old = floor / RATE_MSE
    new_hi, new_lo = floor / RATE_TOT_LO, floor / RATE_TOT_HI
    bar[name] = dict(abs_old=old, abs_lo=new_lo, abs_hi=new_hi,
                     rel_old=old / COS_NOW, rel_lo=new_lo / COS_NOW, rel_hi=new_hi / COS_NOW)
    print(f"  {name:22s}: needs Δcos {old:.5f} (={old/COS_NOW*100:5.1f}% of today's cosine) "
          f"on the OLD rate  ->  {new_lo:.5f}..{new_hi:.5f} "
          f"(={new_lo/COS_NOW*100:4.1f}%..{new_hi/COS_NOW*100:4.1f}%) corrected")
out["bar"] = bar
print("  §8.2's published '+31%' was computed at the pre-cd4f034 cosine 0.079249 on the")
print(f"  mse-only rate; at today's cosine the same construction reads "
      f"{FLOOR_AMP/RATE_MSE/COS_NOW*100:.0f}%.")

print("\n=== 2. Arms re-priced with the DE members ON ===")
rows = []


def price(name, dcos_rel, note, pds_ovr=None, nmae_ovr=None, source=""):
    """dcos_rel = relative cosine gain in BOARD units."""
    dcos = dcos_rel * COS_NOW
    mse = RATE_MSE * dcos
    lo = RATE_TOT_LO * dcos
    hi = RATE_TOT_HI * dcos
    if pds_ovr is not None or nmae_ovr is not None:
        # measured members override the board-rate composite for their own leg
        de4_lo = prop["a14->a15 (f 0.053->0.10)"]["rate_de4_per_cos"] * dcos
        de4_hi = prop["rc_callset->cd4f034 (f 0.10->0.34)"]["rate_de4_per_cos"] * dcos
        base_lo = mse + (pds_ovr or 0.0) + (nmae_ovr if nmae_ovr is not None else de4_lo)
        base_hi = mse + (pds_ovr or 0.0) + (nmae_ovr if nmae_ovr is not None else de4_hi)
        lo, hi = min(base_lo, base_hi), max(base_lo, base_hi)
    r = dict(arm=name, dcos_rel=dcos_rel, dcos_abs=dcos, mse_only=mse,
             six_lo=lo, six_hi=hi, note=note, source=source)
    rows.append(r)
    print(f"\n  {name}")
    print(f"    Δcos {dcos_rel:+.1%} rel = {dcos:+.5f} abs")
    print(f"    mse-only (the old price)  {mse:+.5f}")
    print(f"    six-member, corrected     {lo:+.5f} .. {hi:+.5f}")
    print(f"    {note}")
    return r


price("205 — one extra effective donor line (X-Atlas/Pisces Jurkat, if it existed)",
      0.213,
      f"205 priced it at +0.0062 through the mse-only rate and called it 'a fifth of the gap'. "
      f"Corrected it clears the +0.0089 amplitude floor. 166 proved no such line is downloadable; "
      f"this changes the VALUE of the search, not its feasibility. -> 218.",
      source="experiments/205-xcell-pisces/RESULTS.md:44, :384")

price("205 — two fully independent new lines (k = 6.0)", 0.310,
      "205's own top row, which it recorded as 'clears, but no such pair exists'.",
      source="experiments/205-xcell-pisces/RESULTS.md:388")

price("158 — deployable per-target confidence scale, AS 158 MEASURED IT (Orion pair)",
      0.162, "REFUTED as an input: 165 §7 measured the Orion pair overstating a profile-geometry "
             "cosine gain by 12x and 48x against the board's a23 point. Deflated, f_board is "
             "+0.3%..+2.5% and the arm is dead. Its nmae cost (-0.0029 re-levelled) is measured "
             "and does not deflate.",
      pds_ovr=0.0, nmae_ovr=-0.0029,
      source="experiments/158-target-confidence/RESULTS.md:20, :266; 165 RESULTS.md:333")

price("158 — the same arm at 165's measured Orion deflator (12x, the kinder end)",
      0.162 / 12.0, "pds exactly 0 (measured), reach exactly invariant (measured), fid/jac pinned "
                    "by rc_callset's prescribed call set (208 §2). nmae -0.0029 is the binding leg.",
      pds_ovr=0.0, nmae_ovr=-0.0029,
      source="same")

print("\n=== 3. Arms whose margin turns on the RENDERER, not on the direction ===")
ham = dict(
    arm="208 rider A — Hamilton rounding, on the renderer that actually ships",
    priced_208=+0.001197, pds_208_raw=0.001124, pds_214_raw=0.002069,
    six_member_208=+0.0029, de_208=+0.00170)
ham["renderer_factor"] = ham["pds_214_raw"] / ham["pds_208_raw"]
ham["priced_corrected"] = ham["priced_208"] * ham["renderer_factor"]
ham["six_corrected"] = ham["de_208"] + (ham["six_member_208"] - ham["de_208"]) * ham["renderer_factor"]
print(f"  208 priced it at {ham['priced_208']:+.6f} vs the +0.0019 shape floor -> VETO.")
print(f"  214 §4 re-measured its PLOLO pds on the SHIPPED renderer (VCC_MASS_FIX unset):")
print(f"      {ham['pds_208_raw']:.6f} -> {ham['pds_214_raw']:.6f} raw = "
      f"{ham['renderer_factor']:.2f}x")
print(f"  rule-13 priced composite corrected : {ham['priced_corrected']:+.6f}  "
      f"({'CLEARS' if ham['priced_corrected'] > FLOOR_SHAPE else 'fails'} +0.0019)")
print(f"  six-member corrected              : {ham['six_corrected']:+.6f}")
print("  BUT: 208 proved the arm is NOT pseudobulk-preserving (900/900 blocks differ), so the")
print("  emission-law toll -0.00243 is chargeable; with it the arm is "
      f"{ham['priced_corrected'] - 0.00243:+.6f}. The toll's two calibration points (a21, a22)")
print("  moved the profile ENORMOUSLY; Hamilton moves the challenge tilt -0.01253 -> +0.00030 log2.")
out["hamilton"] = ham

print("\n=== 4. 200's ablation: the rows priced on fewer than six members ===")
two_member = {
    "VCC_ALPHA 0.982 (vs v11's 1.0796)": +0.018542,
    "post-whitening CD4 residual (f 0.10)": +0.002635,
    "pruned basis (union_pruned)": +0.000373,
}
board_cd4_half_step = P0["propagation"]["a14->a15 (f 0.053->0.10)"]["total_ovr"]
print("  three of ten rows are priced on pds+mse only (n=2 in 200's own table):")
for k, v in two_member.items():
    print(f"    {k:42s} {v:+.6f}")
print(f"\n  The board has since measured HALF of the CD4 row directly: the a14->a15 step")
print(f"  (f 0.053 -> 0.10) realised {board_cd4_half_step:+.6f} overall, against 200's")
print(f"  {two_member['post-whitening CD4 residual (f 0.10)']:+.6f} for the WHOLE residual "
      f"(f 0 -> 0.10).")
print(f"  => 200's CD4 row understates by at least "
      f"{board_cd4_half_step/two_member['post-whitening CD4 residual (f 0.10)']:.2f}x, and the")
print("     under-statement is exactly the four DE members it carried at zero.")
out["ablation200"] = dict(two_member_rows=two_member,
                          board_half_step=board_cd4_half_step,
                          understatement=board_cd4_half_step /
                          two_member["post-whitening CD4 residual (f 0.10)"])

out["rows"] = rows
(WORK / "reprice219.json").write_text(json.dumps(out, indent=1, default=float))
print("\nwrote work/reprice219.json")

"""382 item 2: recompute H1-local scaled overalls + deltas from raw members + anchors.

scaled = (raw-b)/(r-b); mse scaled clamped to [0,1]; nmae floored at -6 (no-op
here); overall = mean of six. Compares against 346 results.json full-precision
contrasts AND the 4-decimal assignment headlines.
"""
import json

W346 = "experiments/346-same-line-calibration-h1/work"
anchors = json.load(open("experiments/323-current-family-validation/work/h1-anchors.json"))
b, r = anchors["baseline"], anchors["replicate"]
MEMBERS = ["pds", "mse", "fid", "reach", "jac", "nmae"]

results = json.load(open(f"{W346}/results.json"))
rec_arms = {a["arm"]: a for a in results["arms"]}
rec_contrasts = results["contrasts"]

out = {"arms": {}, "deltas": {}}
for arm in ["cap15", "same_ridge_b05", "same_ridge_b1", "diag_cross_b1"]:
    raw = json.load(open(f"{W346}/metrics-{arm}.json"))["raw"]
    scaled = {}
    for m in MEMBERS:
        s = (raw[m] - b[m]) / (r[m] - b[m])
        if m == "mse":
            s = min(max(s, 0.0), 1.0)
        if m == "nmae":
            s = max(s, -6.0)
        scaled[m] = s
    overall = float(sum(scaled.values()) / 6)
    rec = rec_arms[arm]
    dscale = {m: abs(scaled[m] - rec["local_scaled"][m]) for m in MEMBERS}
    out["arms"][arm] = {
        "scaled_recomputed": scaled,
        "max_absdiff_vs_recorded_local_scaled": max(dscale.values()),
        "overall_recomputed": overall,
        "overall_recorded": rec["local_h1_overall"],
        "overall_absdiff": abs(overall - rec["local_h1_overall"]),
    }

base = out["arms"]["cap15"]["overall_recomputed"]
checks = {"same_ridge_b05": ("+0.0071", "same_ridge_b05-minus-cap15"),
          "same_ridge_b1": ("-0.0149", "same_ridge_b1-minus-cap15"),
          "diag_cross_b1": ("-0.0123", "diag_cross_b1-minus-cap15")}
for arm, (headline, ckey) in checks.items():
    d = out["arms"][arm]["overall_recomputed"] - base
    full = rec_contrasts[ckey]["overall"]
    out["deltas"][arm] = {
        "delta_recomputed": d,
        "contrast_recorded_full_precision": full,
        "absdiff_vs_contrast": abs(d - full),
        "match_contrast_tol1e6": bool(abs(d - full) <= 1e-6),
        "assignment_headline": headline,
        "absdiff_vs_headline": abs(d - float(headline)),
    }

print(json.dumps(out, indent=1))
json.dump(out, open("experiments/382-night-audit/audit2.json", "w"), indent=1)

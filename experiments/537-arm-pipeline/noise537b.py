#!/usr/bin/env python3
"""537 (amendment 1, item 2) — the DE channel's board noise, made auditable.

1. KA-230: 230 §6's "sd 0.027011 scaled = 0.004502 overall" is recomputed from 230's own table (work/fit230.json): the six
   digest-identical arms, each a DIFFERENT emission law against its incumbent. Abort unless it reproduces to 1e-12.
2. Board pairs of ours, by entry id, raw members from data/board-entries/<id>.json (nothing else is read):
   B1 bulk-only changes with the same DE design and the same builder seed (SEED 0): the ballast count K changes per block,
      so the carriers are re-drawn on other cells. Their DE deltas are a re-roll of one design (+ the ballast's nmae term).
   B2 DE changes that share every emission draw (468: prior-only; zs is a post-pass on the same blocks) or touch few blocks
      (the repair floor): their deltas are the design's effect, with no re-roll.
3. H1 (real production emitter, 85 targets): for every DE-knob class with >= 2 emission seeds on record, the seed-to-seed sd
   of the PAIRED effect (arm - its base, same seed) per member, and for 435's base the per-target re-roll. Scaled to the
   board's 900 target-rows by sqrt(85 / 900) (a panel mean of independent per-target terms). This is the read noise of an
   arm of that class: it answers whether a knob that re-draws the emission (call set, hook) inherits B1's sigma.
Units: raw members; overall = raw x (1/span) / 6 with 192's board-fitted INV (fid 3.32482, jac 2.72081, reach 1.12927,
nmae 1.64787). Writes work/noise537b.json. Light (< 0.5 GB).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INV = dict(fid=3.32482, jac=2.72081, reach=1.12927, nmae=1.64787)
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
ID = {  # arm -> board entry id (docs/BOARD-LEDGER.md; submissions/SUBMITTED_LEDGER.json)
    "410b": "4WgVi6FIE4iVw9PmRUP3", "417a": "T88GsCjp53Bmk8xQcJMS", "417e": "wPP8U70kj4mOIvtA83CO",
    "433a": "FCbBcvGxF1HOFPxPISUC", "433h": "wMoDAi1mca08FA2Bfqkd", "433k": "T76JWD3tOajk5IYvzwUW",
    "446b": "UDgjOKL9gI1udAexq6ue", "446bja": "AsHteaDAXCWQo57fLYti", "446bj2": "YvVru3DAF89Jiyr34olR",
    "446br2": "wE7usCaR6TAlPp16TR0V", "446bj2z": "be90smZUDcIV7j5MgQ03", "446bd2": "R5uE5Rve3R9nWXGbeSLM",
    "453dr": "q3VuTImh4FRGKXH2waCp", "457dr": "1FfirZsGsAZvaXYTd39J", "458d": "b7UI3NRy1rdVNRRyZwmU",
    "458g": "JL8llNLHaTyDQzm798TV",
}
B1 = [("410b", "417a", "sharpen q 1.25 (bulk through the ballast)"),
      ("417a", "417e", "bulk stack H1 x2 + CD4 x1.5 via VCC_BULK_LOAD"),
      ("433h", "433k", "sharpen q 1.5 + re-level"),
      ("433k", "446b", "sharpen q 1.75 + re-level"),
      ("446b", "446bj2", "LODO bulk key + re-level"),
      ("446b", "446bd2", "DLD-1 as the second union residual"),
      ("453dr", "457dr", "October path: CD4 x1.5 in the stack")]
B2 = [("446bj2", "446bj2z", "zs post-pass (same blocks)"),
      ("446b", "446bja", "prior-only: A's block prior (468: shares draws)"),
      ("446b", "458d", "prior-only: Arc programme block prior"),
      ("458d", "458g", "prior-only: the q33 gate"),
      ("446b", "446br2", "zero-repair floor -5 -> -2 (a few blocks)")]


def entry(eid):
    d = json.loads((ROOT / f"data/board-entries/{eid}.json").read_text())
    e = d.get("entry", d)
    assert e["entry_id"] == eid
    return {m: float(e[k]) for m, k in RAW.items()} | {"overall": float(e["score_avg"])}


def ov(m, v):
    return v * INV[m] / 6.0 * (-1.0 if m == "nmae" else 1.0)


def ka230():
    t = json.loads((ROOT / "experiments/230-toll-magnitude/work/fit230.json").read_text())
    rows = [r for r in t["table"] if r["digest_same"] == 900]
    z = np.array([r["d_reach"] for r in rows])
    sd = float(z.std(ddof=1))
    assert abs(sd - 0.027010964677330575) < 1e-12 and abs(sd / 6 - 0.004501827446221763) < 1e-12, sd
    assert abs(sd - t["falsifier2"]["zero_sd"]) < 1e-15
    return {"pairs": [{k: r[k] for k in ("arm", "incumbent", "what", "d_reach", "lib_cv", "frac_resampled")} for r in rows],
            "sd_scaled": sd, "sd_overall": sd / 6}


def board():
    out = {}
    for name, pairs in (("B1", B1), ("B2", B2)):
        rows = []
        for a, b, what in pairs:
            ea, eb = entry(ID[a]), entry(ID[b])
            d = {m: eb[m] - ea[m] for m in RAW}
            rows.append({"pair": f"{b} - {a}", "ids": [ID[b], ID[a]], "what": what, **{m: d[m] for m in RAW},
                         "overall": eb["overall"] - ea["overall"],
                         "de_overall": sum(ov(m, d[m]) for m in ("fid", "jac", "reach", "nmae"))})
        out[name] = rows
    r = np.array([x["reach"] for x in out["B1"]])
    f = np.array([x["fid"] for x in out["B1"]])
    j = np.array([x["jac"] for x in out["B1"]])
    n = np.array([x["nmae"] for x in out["B1"]])
    n_ex = np.array([x["nmae"] for x in out["B1"] if x["pair"] not in ("433k - 433h", "446b - 433k")])
    k = r.size
    # the rms of k paired differences estimates sigma_diff; chi-square 95 % interval for it
    from scipy.stats import chi2
    lo, hi = np.sqrt(k / chi2.ppf(0.975, k)), np.sqrt(k / chi2.ppf(0.025, k))
    rms = lambda v: float(np.sqrt(np.mean(v ** 2)))  # noqa: E731
    out["B1_summary"] = {
        "n_pairs": int(k),
        "reach_rms_raw": rms(r), "reach_rms_raw_ci95": [rms(r) * lo, rms(r) * hi],
        "reach_rms_overall": ov("reach", rms(r)), "reach_rms_overall_ci95": [ov("reach", rms(r) * lo), ov("reach", rms(r) * hi)],
        "fid_rms_raw": rms(f), "jac_rms_raw": rms(j), "nmae_rms_raw": rms(n), "nmae_rms_raw_excl_q_steps": rms(n_ex),
        "de_overall_rms": rms(np.array([x["de_overall"] for x in out["B1"]])),
        "pds_mse_moved": "yes (bulk arms)"}
    return out


def h1():
    """seed sd of paired effects (arm - base, same seed) on H1, per member, and the board-scaled read noise."""
    p528 = json.loads((ROOT / "experiments/528-callprop-block/work/price528.json").read_text())["arms"]
    r517 = json.loads((ROOT / "experiments/517-flex-programme/work/real517.json").read_text())
    scale = np.sqrt(85.0 / 900.0)
    classes = {}
    # 528: the call-set block (zs form), per seed at f 1 (the same targets / pools as the board read's f is unknown)
    for arm in ("lcap-kept", "lcap", "lcap-dec", "shuf", "ltv", "lcap-kept+mzg", "lcap+mzg", "lcap-kept-n1200"):
        a = p528.get(f"{arm} +zs")
        if a is None or len(a["seeds"]) < 2:
            continue
        ps = a["1.0"]["per_seed"]
        classes[f"528 {arm} (vs {a['base']})"] = {m: [p["d"][m] for p in ps] for m in ("fid", "jac", "reach", "nmae")}
    # 517: the uncalled-gene hook (real emitter), per seed
    for arm in ("c4-mzg0.05", "c4-arcg0.05", "c4-arcg0.05-r0.25", "g33-arcg0.05"):
        ps = r517.get(arm, [])
        if len(ps) < 2:
            continue
        classes[f"517 {arm} (vs its no-hook base)"] = {m: [p[m] for p in ps] for m in ("fid", "jac", "reach", "nmae")}
    res = {}
    for name, d in classes.items():
        row = {}
        for m, v in d.items():
            v = np.asarray(v, float)
            sd = float(v.std(ddof=1))
            row[m] = dict(values=v.tolist(), h1_seed_sd=sd, board_sd_raw=sd * scale, board_sd_overall=abs(ov(m, sd * scale)))
        row["all_four_overall"] = float(np.sqrt(sum(row[m]["board_sd_overall"] ** 2 for m in ("fid", "jac", "reach", "nmae"))))
        row["n_seeds"] = len(next(iter(d.values())))
        res[name] = row
    # pooled over classes (variance-weighted by df)
    pooled = {}
    for m in ("fid", "jac", "reach", "nmae"):
        num = sum((res[c]["n_seeds"] - 1) * res[c][m]["h1_seed_sd"] ** 2 for c in res)
        den = sum(res[c]["n_seeds"] - 1 for c in res)
        sd = float(np.sqrt(num / den))
        pooled[m] = dict(h1_seed_sd=sd, board_sd_raw=sd * scale, board_sd_overall=abs(ov(m, sd * scale)), df=int(den))
    pooled["all_four_overall"] = float(np.sqrt(sum(pooled[m]["board_sd_overall"] ** 2 for m in ("fid", "jac", "reach", "nmae"))))
    base = json.loads((HERE / "work/noise537.json").read_text())["h1"]
    return {"classes": res, "pooled_paired": pooled, "base_reroll_independent": base}


def main():
    k = ka230()
    print(f"KA-230 PASS: six digest-identical arms, sd {k['sd_scaled']:.12f} scaled = {k['sd_overall']:.12f} overall")
    for r in k["pairs"]:
        print(f"   {r['arm']:10s} vs {r['incumbent']:10s} d_reach {r['d_reach']:+.6f}  lib CV {r['lib_cv']:.4f}  {r['what']}")
    b = board()
    for name in ("B1", "B2"):
        print(f"\n{name}: pair (ids) | raw d fid, jac, reach, nmae, pds, mse | DE overall")
        for x in b[name]:
            print(f"  {x['pair']:18s} {x['ids'][0]} - {x['ids'][1]} | {x['fid']:+.6f} {x['jac']:+.6f} {x['reach']:+.6f} "
                  f"{x['nmae']:+.6f} {x['pds']:+.6f} {x['mse']:+.6f} | {x['de_overall']:+.6f}  ({x['what']})")
    s = b["B1_summary"]
    print(f"\nB1 (n {s['n_pairs']}): reach rms {s['reach_rms_raw']:.5f} raw [{s['reach_rms_raw_ci95'][0]:.5f}, "
          f"{s['reach_rms_raw_ci95'][1]:.5f}] = {s['reach_rms_overall']:.5f} overall [{s['reach_rms_overall_ci95'][0]:.5f}, "
          f"{s['reach_rms_overall_ci95'][1]:.5f}]; fid {s['fid_rms_raw']:.6f}, jac {s['jac_rms_raw']:.6f}, nmae {s['nmae_rms_raw']:.6f} "
          f"({s['nmae_rms_raw_excl_q_steps']:.6f} without the two q steps' ballast term); DE overall rms {s['de_overall_rms']:.5f}")
    h = h1()
    print("\nH1 paired-effect seed sd (arm - base, same seed), board-scaled x sqrt(85/900):")
    for c, r in h["classes"].items():
        print(f"  {c:48s} n {r['n_seeds']} | " + "  ".join(
            f"{m} {r[m]['h1_seed_sd']:.5f} -> {r[m]['board_sd_raw']:.5f} raw" for m in ("fid", "jac", "reach", "nmae"))
            + f" | overall {r['all_four_overall']:.5f}")
    p = h["pooled_paired"]
    print("  POOLED: " + "  ".join(f"{m} board {p[m]['board_sd_raw']:.6f} raw ({p[m]['board_sd_overall']:.6f} ov, df {p[m]['df']})"
                                  for m in ("fid", "jac", "reach", "nmae")) + f" | all four {p['all_four_overall']:.5f} overall")
    (HERE / "work/noise537b.json").write_text(json.dumps({"ka230": k, "board": b, "h1": h}, indent=1))


if __name__ == "__main__":
    main()

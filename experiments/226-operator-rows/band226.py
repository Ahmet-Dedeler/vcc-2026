"""226 — merge every run, put a SEED BAND on the curve, and price the result.

Phase 1 produced a curve with a huge rank cliff below `n = p` and a small, non-monotone scatter
above it.  A scatter is not a slope, so the whole verdict turns on one number: how big is the
draw-to-draw spread of a refit at FIXED `n`?  This merges

    work/rows226.json          the 7-point curve, seed 1
    work/rows226_band.json     5 independent draws at n = 22,045, plus the permuted nulls
    work/rows226_caplift.json  the same donors with the 12,000-row cap LIFTED -- the only way
                               to extend the axis ABOVE production without changing composition

and reports the full-rank points against the band, not against each other.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = HERE / "work"
SPAN_PDS = 0.452
PDS_OVR = 1.011 / 6.0
FRAMES = ("chalA", "chalB", "chalC")
SOURCES = ("rows226.json", "rows226_band.json", "rows226_caplift.json")


def per_target(row, field):
    return np.mean([np.asarray(row["per_frame"][f][field], float) for f in FRAMES], axis=0)


def main() -> None:
    runs, rows = {}, []
    for name in SOURCES:
        p = W / name
        if not p.exists():
            print(f"  (missing {name})")
            continue
        d = json.loads(p.read_text())
        runs[name] = d["pool"]
        for r in d["rows"]:
            r = dict(r)
            r["source"] = name
            r["pool_total"] = d["pool"]["total"]
            r["cap"] = d["pool"]["cap"]
            rows.append(r)
    p_dim = next(iter(runs.values()))["p"]

    real = [r for r in rows if not r["permuted"]]
    perm = [r for r in rows if r["permuted"]]

    # ---- the seed band: every (n) with more than one independent draw -------------------
    print("=" * 100)
    print("THE SEED BAND — independent subsamples at the SAME row count")
    print("=" * 100)
    band = {}
    by_n: dict[int, list] = {}
    for r in real:
        by_n.setdefault(r["n_rows"], []).append(r)
    for n, rs in sorted(by_n.items()):
        if len(rs) < 2:
            continue
        c = np.array([x["chal_cos"] for x in rs])
        ca = np.array([x["chal_cos_agg"] for x in rs])
        pd_ = np.array([x["chal_pds"] for x in rs])
        band[n] = dict(n_draws=len(rs), seeds=[x["seed"] for x in rs],
                       cos_mean=float(c.mean()), cos_sd=float(c.std(ddof=1)),
                       cos_min=float(c.min()), cos_max=float(c.max()),
                       cos_agg_mean=float(ca.mean()), cos_agg_sd=float(ca.std(ddof=1)),
                       pds_mean=float(pd_.mean()), pds_sd=float(pd_.std(ddof=1)),
                       pds_min=float(pd_.min()), pds_max=float(pd_.max()))
        b = band[n]
        print(f"  n = {n}  ({len(rs)} draws, seeds {b['seeds']})")
        print(f"    cos      mean {b['cos_mean']:.7f}  sd {b['cos_sd']:.7f} "
              f"({100*b['cos_sd']/b['cos_mean']:.2f}%)  range [{b['cos_min']:.7f}, {b['cos_max']:.7f}]")
        print(f"    cos_agg  mean {b['cos_agg_mean']:.7f}  sd {b['cos_agg_sd']:.7f}")
        print(f"    pds      mean {b['pds_mean']:.6f}  sd {b['pds_sd']:.6f} "
              f"  range [{b['pds_min']:.6f}, {b['pds_max']:.6f}]")

    # ---- every full-rank point, expressed in band sigmas --------------------------------
    # `f = 1.0` consumes the whole pool at every seed, so its "draws" are the same operator
    # and its sd is exactly 0 by construction -- it measures reproducibility, not subsample
    # noise.  The reference band must be one with genuine draw-to-draw variation.
    varying = {n: b for n, b in band.items() if b["cos_sd"] > 0}
    ref_n = max(varying, key=lambda n: varying[n]["n_draws"]) if varying else None
    print("\n" + "=" * 100)
    print("THE FULL-RANK CURVE, IN UNITS OF THE BAND'S OWN SIGMA")
    if ref_n:
        print(f"reference band: n = {ref_n}, {band[ref_n]['n_draws']} draws  "
              f"(cos sd {band[ref_n]['cos_sd']:.7f}, pds sd {band[ref_n]['pds_sd']:.6f})")
    print("=" * 100)
    print(f"{'n rows':>8} {'n/p':>6} {'cap':>6} {'seed':>5} {'floored':>8} "
          f"{'cos':>11} {'z(cos)':>8} {'cos_agg':>11} {'pds':>9} {'z(pds)':>8}")
    table = []
    fr = sorted([r for r in real if r["n_rows"] >= p_dim], key=lambda r: (r["n_rows"], r["seed"]))
    for r in fr:
        zc = zp = None
        if ref_n:
            b = band[ref_n]
            zc = (r["chal_cos"] - b["cos_mean"]) / b["cos_sd"]
            zp = (r["chal_pds"] - b["pds_mean"]) / b["pds_sd"]
        table.append(dict(n_rows=r["n_rows"], n_over_p=r["n_rows"] / p_dim, cap=r["cap"],
                          seed=r["seed"], cos=r["chal_cos"], cos_agg=r["chal_cos_agg"],
                          pds=r["chal_pds"], n_floored=r["n_floored"],
                          z_cos=zc, z_pds=zp, source=r["source"]))
        print(f"{r['n_rows']:>8d} {r['n_rows']/p_dim:>6.3f} {r['cap']:>6d} {r['seed']:>5d} "
              f"{r['n_floored']:>8d} {r['chal_cos']:>11.7f} "
              f"{(zc if zc is not None else float('nan')):>+8.2f} "
              f"{r['chal_cos_agg']:>11.7f} {r['chal_pds']:>9.6f} "
              f"{(zp if zp is not None else float('nan')):>+8.2f}")

    # ---- is the curve still rising at production's n? -----------------------------------
    verdict = {}
    # A SECOND, independent read on the draw noise, this time AT production's own row count:
    # the cap-lifted pool subsampled to 75% is n = 29,941 -- the same size as production's
    # 29,394 but a completely different draw of rows.  Two operators, same n, same donors.
    twin = [r for r in real if 29000 <= r["n_rows"] <= 30500]
    if len(twin) >= 2:
        print("\n  independent draws AT production's row count (second noise read):")
        for r in sorted(twin, key=lambda x: x["n_rows"]):
            print(f"    n={r['n_rows']:6d} cap={r['cap']:5d}  cos {r['chal_cos']:.7f}  "
                  f"pds {r['chal_pds']:.6f}")
        cs = [r["chal_cos"] for r in twin]
        ps = [r["chal_pds"] for r in twin]
        verdict["twin_at_production_n"] = dict(
            n=[r["n_rows"] for r in twin], cos=cs, pds=ps,
            cos_spread_pct=100 * (max(cs) - min(cs)) / np.mean(cs),
            pds_spread=max(ps) - min(ps))
        print(f"    spread: cos {100*(max(cs)-min(cs))/np.mean(cs):.2f}%   "
              f"pds {max(ps)-min(ps):.6f}")

    if ref_n and band:
        b = band[ref_n]
        prod = [r for r in real if r["n_rows"] == 29394 and r["cap"] == 12000]
        if prod:
            pr = prod[0]
            dc = pr["chal_cos"] - b["cos_mean"]
            dp = pr["chal_pds"] - b["pds_mean"]
            se_c = b["cos_sd"] / np.sqrt(b["n_draws"])
            se_p = b["pds_sd"] / np.sqrt(b["n_draws"])
            verdict["prod_vs_band"] = dict(
                n_band=ref_n, n_prod=pr["n_rows"],
                d_cos=dc, d_cos_pct=100 * dc / b["cos_mean"],
                z_cos_single=dc / b["cos_sd"], z_cos_mean=dc / se_c,
                d_pds=dp, d_pds_pct=100 * dp / b["pds_mean"],
                z_pds_single=dp / b["pds_sd"], z_pds_mean=dp / se_p,
                d_pds_overall=dp / SPAN_PDS * PDS_OVR)
            v = verdict["prod_vs_band"]
            print(f"\n  n {ref_n} (band mean) -> {pr['n_rows']} (production, exact):")
            print(f"    Δcos  {dc:+.7f}  ({v['d_cos_pct']:+.2f}%)   "
                  f"{v['z_cos_single']:+.2f} σ of a single draw, {v['z_cos_mean']:+.2f} σ of the band mean")
            print(f"    Δpds  {dp:+.6f}  ({v['d_pds_pct']:+.3f}%)   "
                  f"{v['z_pds_single']:+.2f} σ of a single draw, {v['z_pds_mean']:+.2f} σ of the band mean")
            print(f"    Δpds in overall units: {v['d_pds_overall']:+.6f}")
        # and the point ABOVE production, if the cap-lift ran
        above = sorted([r for r in real if r["n_rows"] > 29394], key=lambda r: r["n_rows"])
        if above and prod:
            hi = above[-1]
            dc = hi["chal_cos"] - prod[0]["chal_cos"]
            dp = hi["chal_pds"] - prod[0]["chal_pds"]
            verdict["above_production"] = dict(
                n_lo=prod[0]["n_rows"], n_hi=hi["n_rows"], cap_hi=hi["cap"],
                d_cos=dc, d_cos_pct=100 * dc / prod[0]["chal_cos"],
                z_cos=dc / b["cos_sd"], d_pds=dp, z_pds=dp / b["pds_sd"],
                d_pds_overall=dp / SPAN_PDS * PDS_OVR,
                cos_lo=prod[0]["chal_cos"], cos_hi=hi["chal_cos"],
                pds_lo=prod[0]["chal_pds"], pds_hi=hi["chal_pds"])
            v = verdict["above_production"]
            print(f"\n  ABOVE production: n {v['n_lo']} -> {v['n_hi']} (cap lifted to {v['cap_hi']}):")
            print(f"    Δcos  {dc:+.7f}  ({v['d_cos_pct']:+.2f}%)   {v['z_cos']:+.2f} σ")
            print(f"    Δpds  {dp:+.6f}   {v['z_pds']:+.2f} σ   "
                  f"overall {v['d_pds_overall']:+.6f}")

    # ---- falsifier 2 --------------------------------------------------------------------
    if perm:
        print("\n" + "=" * 100)
        print("FALSIFIER 2 — gene-label-permuted refit (same n, structure destroyed)")
        print("=" * 100)
        f2 = []
        for r in sorted(perm, key=lambda x: x["n_rows"]):
            m = [q for q in real if q["n_rows"] == r["n_rows"] and q["cap"] == r["cap"]]
            rec = dict(n_rows=r["n_rows"], cos=r["chal_cos"], cos_agg=r["chal_cos_agg"],
                       pds=r["chal_pds"], n_floored=r["n_floored"],
                       real_cos=m[0]["chal_cos"] if m else None,
                       real_pds=m[0]["chal_pds"] if m else None)
            f2.append(rec)
            print(f"  n={r['n_rows']:6d}  permuted: cos {r['chal_cos']:.7f} pds "
                  f"{r['chal_pds']:.6f} floored {r['n_floored']}")
            if m:
                print(f"           real:     cos {m[0]['chal_cos']:.7f} pds "
                      f"{m[0]['chal_pds']:.6f} floored {m[0]['n_floored']}")
        verdict["permuted_null"] = f2
        print("  212's kolf2 gamma=0 row (NO whitening at all): cos 0.0092404  pds 0.5885361")
        print("  -> permuting the gene labels collapses C toward isotropy, so the permuted")
        print("     operator IS approximately gamma=0.  A cosine gain that looks like this is")
        print("     a weaker operator, not a better-sampled one.")

    out = dict(pools=runs, band=band, full_rank_table=table, verdict=verdict,
               p=p_dim, span_pds=SPAN_PDS, pds_overall_factor=PDS_OVR)
    (W / "band226.json").write_text(json.dumps(out, indent=1))
    print("\nwrote work/band226.json")


if __name__ == "__main__":
    main()

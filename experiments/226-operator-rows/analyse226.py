"""226 analysis — the cosine-vs-fitting-rows curve, its trade rate against the gamma axis,
and falsifiers 2/3/4.  Pure post-processing: reads `work/*.json`, fits nothing, writes a table.

The pricing question this has to answer
---------------------------------------
212 closed gamma with a RATIO: along the gamma axis the operator buys transfer cosine and `pds`
pays for it, and the two cancel to 3.8%.  A ratio is instrument-independent, which is why that
closure survived the 2026-09-04 discovery that the PLOLO->board cosine factor is not a constant.

The row axis has to be judged the same way, on the same fold, the same render and the same code
path.  So the number here is `d(pds)/d(cos)` along ROWS, quoted against `d(pds)/d(cos)` along
GAMMA.  If they match, "more rows" is gamma in disguise and 212 already closed it.  If the row
axis buys cosine without the `pds` charge, it is a new lever.

Nothing here multiplies a PLOLO cosine by a scalar to reach the board (AGENTS.md 2026-09-04).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = HERE / "work"

SPAN_PDS = 0.452           # 162's board span for pds
PDS_OVR = 1.011 / 6.0      # 193 §5: overall contribution = scaled/6 x 1.011
FRAMES = ("chalA", "chalB", "chalC")


def load(name):
    p = W / name
    return json.loads(p.read_text()) if p.exists() else None


def per_target(row, field):
    """Mean across the three challenge frames of a per-target vector."""
    return np.mean([np.asarray(row["per_frame"][f][field], float) for f in FRAMES], axis=0)


def retention(row_lo, row_hi, field, k=3):
    """Falsifier 4 — leave-top-`k`-out retention of the lo->hi effect on ONE channel.

    The `k` targets with the largest individual contribution to the effect are removed and the
    macro mean is recomputed on the remaining targets.  The retrieval pool is NOT re-run, so
    this measures how concentrated the effect is, which is what 222 §9 compares across channels.
    """
    a, b = per_target(row_lo, field), per_target(row_hi, field)
    d = b - a
    full = float(d.mean())
    top = np.argsort(-np.abs(d))[:k]
    keep = np.ones(len(d), bool)
    keep[top] = False
    trimmed = float(d[keep].mean())
    return dict(full=full, trimmed=trimmed,
                retention=(trimmed / full if full else None),
                top_idx=top.tolist(), top_contrib=d[top].tolist(),
                n_targets=int(len(d)), n_kept=int(keep.sum()))


def main() -> None:
    cur = load("rows226.json")
    gam = load("gamma226.json")
    if cur is None:
        raise SystemExit("work/rows226.json missing")

    real = sorted([r for r in cur["rows"] if not r["permuted"]], key=lambda r: r["n_rows"])
    perm = sorted([r for r in cur["rows"] if r["permuted"]], key=lambda r: r["n_rows"])
    p = cur["pool"]["p"]
    base = [r for r in real if r["frac"] == 1.0][0]

    out = {"pool": cur["pool"], "fold": cur["fold"], "config": cur["config"],
           "n_targets": cur["n_targets"], "curve": [], "gamma_trade": (gam or {}).get("gamma_trade"),
           "falsifiers": {}}

    print("=" * 104)
    print(f"PHASE 1 — cosine vs fitting rows.  fold {cur['fold']}, config {cur['config']}, "
          f"p = {p}, {cur['n_targets']} targets")
    print("=" * 104)
    print(f"{'frac':>6} {'n rows':>8} {'n/p':>6} {'rank-def':>9} "
          f"{'cos(macro)':>11} {'Δ% vs 100%':>11} {'cos_agg':>10} {'Δ%':>8} "
          f"{'pds':>9} {'Δ raw':>10} {'cos_wrong':>10}")
    for r in real:
        dcos = 100 * (r["chal_cos"] - base["chal_cos"]) / base["chal_cos"]
        dagg = 100 * (r["chal_cos_agg"] - base["chal_cos_agg"]) / base["chal_cos_agg"]
        dpds = r["chal_pds"] - base["chal_pds"]
        rec = dict(frac=r["frac"], n_rows=r["n_rows"], n_over_p=r["n_rows"] / p,
                   rank_deficient=bool(r["n_rows"] < p), cos=r["chal_cos"],
                   cos_agg=r["chal_cos_agg"], pds=r["chal_pds"],
                   cos_wrong=r.get("chal_cos_wrong"), n_floored=r["n_floored"],
                   d_cos_pct=dcos, d_cos_agg_pct=dagg, d_pds_raw=dpds,
                   d_pds_scaled=dpds / SPAN_PDS, d_pds_overall=dpds / SPAN_PDS * PDS_OVR)
        out["curve"].append(rec)
        print(f"{r['frac']:>6.3f} {r['n_rows']:>8d} {r['n_rows']/p:>6.3f} "
              f"{'YES' if r['n_rows'] < p else 'no':>9} "
              f"{r['chal_cos']:>11.7f} {dcos:>+11.3f} {r['chal_cos_agg']:>10.7f} "
              f"{dagg:>+8.3f} {r['chal_pds']:>9.6f} {dpds:>+10.6f} "
              f"{(r.get('chal_cos_wrong') or float('nan')):>10.7f}")

    # ---------------- the saturation test the mandate asks for -----------------------
    def get(f):
        m = [r for r in real if abs(r["frac"] - f) < 1e-9]
        return m[0] if m else None

    half, full = get(0.5), get(1.0)
    tests = {}
    if half and full:
        tests["50_to_100"] = {
            "d_cos_pct": 100 * (full["chal_cos"] - half["chal_cos"]) / half["chal_cos"],
            "d_cos_agg_pct": 100 * (full["chal_cos_agg"] - half["chal_cos_agg"]) / half["chal_cos_agg"],
            "d_pds_raw": full["chal_pds"] - half["chal_pds"],
            "n_lo": half["n_rows"], "n_hi": full["n_rows"]}
    # the top of the curve: is it still rising where it matters?
    top = [r for r in real if r["n_rows"] >= p]
    if len(top) >= 2:
        lo, hi = top[0], top[-1]
        tests["full_rank_segment"] = {
            "n_lo": lo["n_rows"], "n_hi": hi["n_rows"],
            "d_cos_pct": 100 * (hi["chal_cos"] - lo["chal_cos"]) / lo["chal_cos"],
            "d_cos_agg_pct": 100 * (hi["chal_cos_agg"] - lo["chal_cos_agg"]) / lo["chal_cos_agg"],
            "d_pds_raw": hi["chal_pds"] - lo["chal_pds"]}
        # the LAST leg -- the local slope at production's own row count
        prev = top[-2]
        dcos = hi["chal_cos"] - prev["chal_cos"]
        dpds = hi["chal_pds"] - prev["chal_pds"]
        tests["local_slope_at_production"] = {
            "n_lo": prev["n_rows"], "n_hi": hi["n_rows"], "d_cos": dcos, "d_pds": dpds,
            "d_pds_raw": dpds,
            "d_cos_pct": 100 * dcos / prev["chal_cos"],
            "d_cos_agg_pct": 100 * (hi["chal_cos_agg"] - prev["chal_cos_agg"]) / prev["chal_cos_agg"],
            "d_pds_per_d_cos": (dpds / dcos if dcos else None)}
    out["saturation"] = tests

    print("\n" + "=" * 104)
    print("SATURATION TESTS")
    print("=" * 104)
    for k, v in tests.items():
        print(f"  {k}: n {v['n_lo']} -> {v['n_hi']}  Δcos {v['d_cos_pct']:+.3f}%  "
              f"Δcos_agg {v.get('d_cos_agg_pct', float('nan')):+.3f}%  "
              f"Δpds {v['d_pds_raw' if 'd_pds_raw' in v else 'd_pds']:+.6f} raw")

    # ---------------- the trade rate, against gamma's ---------------------------------
    print("\n" + "=" * 104)
    print("THE TRADE RATE — is the row axis the gamma axis in disguise?")
    print("=" * 104)
    if len(real) >= 2:
        # regress pds on cos across the whole curve, and across the full-rank part only
        for label, rows in (("all points", real), ("full-rank only", top)):
            if len(rows) < 2:
                continue
            x = np.array([r["chal_cos"] for r in rows])
            y = np.array([r["chal_pds"] for r in rows])
            slope = float(np.polyfit(x, y, 1)[0])
            rr = float(np.corrcoef(x, y)[0, 1])
            out.setdefault("row_trade", {})[label] = dict(d_pds_per_d_cos=slope, r=rr,
                                                          n_points=len(rows))
            print(f"  ROW axis  ({label:15s}): d(pds)/d(cos) = {slope:+9.3f}   r = {rr:+.4f}")
    if gam:
        for k, v in gam["gamma_trade"].items():
            if k.startswith("a15"):
                print(f"  GAMMA axis({k:15s}): d(pds)/d(cos) = {v['d_pds_per_d_cos']:+9.3f}   "
                      f"(Δcos {v['rel_d_cos_pct']:+.3f}%)")

    # ---------------- falsifier 2: the gene-label-permuted refit ----------------------
    if perm:
        print("\nFALSIFIER 2 — gene-label-permuted refit (same n, no structure)")
        f2 = []
        for r in perm:
            m = [q for q in real if q["n_rows"] == r["n_rows"]]
            rec = dict(n_rows=r["n_rows"], cos=r["chal_cos"], pds=r["chal_pds"],
                       real_cos=m[0]["chal_cos"] if m else None,
                       real_pds=m[0]["chal_pds"] if m else None)
            f2.append(rec)
            print(f"  n={r['n_rows']:6d}  permuted cos {r['chal_cos']:.7f} pds {r['chal_pds']:.6f}"
                  + (f"   |  real cos {m[0]['chal_cos']:.7f} pds {m[0]['chal_pds']:.6f}"
                     if m else ""))
        if len(f2) >= 2:
            a, b = f2[0], f2[-1]
            out["falsifiers"]["permuted_null"] = dict(
                rows=f2,
                d_cos_pct=100 * (b["cos"] - a["cos"]) / a["cos"],
                d_pds_raw=b["pds"] - a["pds"])
            print(f"  permuted-null slope over n: Δcos "
                  f"{100*(b['cos']-a['cos'])/a['cos']:+.3f}%  Δpds {b['pds']-a['pds']:+.6f}")
        else:
            out["falsifiers"]["permuted_null"] = dict(rows=f2)

    # ---------------- falsifier 3: the wrong recipient --------------------------------
    if real[0].get("chal_cos_wrong") is not None:
        print("\nFALSIFIER 3 — wrong-recipient control (fixed derangement, same operator)")
        lo, hi = real[0], real[-1]
        dw = 100 * (hi["chal_cos_wrong"] - lo["chal_cos_wrong"]) / abs(lo["chal_cos_wrong"])
        do = 100 * (hi["chal_cos"] - lo["chal_cos"]) / lo["chal_cos"]
        out["falsifiers"]["wrong_recipient"] = dict(
            n_lo=lo["n_rows"], n_hi=hi["n_rows"],
            cos_own_lo=lo["chal_cos"], cos_own_hi=hi["chal_cos"], d_own_pct=do,
            cos_wrong_lo=lo["chal_cos_wrong"], cos_wrong_hi=hi["chal_cos_wrong"],
            d_wrong_pct=dw, separation=do - dw,
            note="cos_wrong is NEGATIVE at every n; growing n makes it MORE negative, i.e. "
                 "mis-paired predictions get less similar while own pairs get more similar")
        print(f"  n {lo['n_rows']} -> {hi['n_rows']}:  own-pair cos {lo['chal_cos']:.7f} -> "
              f"{hi['chal_cos']:.7f} ({do:+.1f}%)   wrong-pair cos "
              f"{lo['chal_cos_wrong']:.7f} -> {hi['chal_cos_wrong']:.7f}")
        print("  every point, own vs wrong:")
        for r in real:
            print(f"    n={r['n_rows']:6d}  own {r['chal_cos']:+.7f}   "
                  f"wrong {r['chal_cos_wrong']:+.7f}   "
                  f"ratio {r['chal_cos']/abs(r['chal_cos_wrong']):6.2f}x")

    # ---------------- falsifier 4: leave-top-3-out ------------------------------------
    if "cos_per_target" in real[0]["per_frame"]["chalA"]:
        print("\nFALSIFIER 4 — leave-top-3-out retention (cosine and pds, separately)")
        lo, hi = real[0], real[-1]
        f4 = {"n_lo": lo["n_rows"], "n_hi": hi["n_rows"],
              "cos": retention(lo, hi, "cos_per_target"),
              "pds": retention(lo, hi, "pds_per_target")}
        out["falsifiers"]["leave_top3_out"] = f4
        for ch in ("cos", "pds"):
            v = f4[ch]
            ret = v["retention"]
            print(f"  {ch:4s}: full Δ {v['full']:+.7f}  minus top 3 {v['trimmed']:+.7f}  "
                  f"retention {100*ret:6.1f}%" if ret is not None else f"  {ch}: no effect")
        # and the same on the full-rank segment, which is the one the verdict rests on
        if len(top) >= 2:
            f4b = {"n_lo": top[0]["n_rows"], "n_hi": top[-1]["n_rows"],
                   "cos": retention(top[0], top[-1], "cos_per_target"),
                   "pds": retention(top[0], top[-1], "pds_per_target")}
            out["falsifiers"]["leave_top3_out_full_rank"] = f4b
            for ch in ("cos", "pds"):
                v = f4b[ch]
                ret = v["retention"]
                print(f"  [full-rank {top[0]['n_rows']}->{top[-1]['n_rows']}] {ch:4s}: "
                      f"full Δ {v['full']:+.7f}  minus top 3 {v['trimmed']:+.7f}  "
                      f"retention {100*ret:6.1f}%" if ret is not None else "")

    (W / "analyse226.json").write_text(json.dumps(out, indent=1))
    print("\nwrote work/analyse226.json")


if __name__ == "__main__":
    main()

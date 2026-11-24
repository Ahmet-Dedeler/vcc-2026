#!/usr/bin/env python3
"""553 — KA-S0 (light): lib553.score_stats on 323's stored H1 reference (its truth arrays and reference_de.parquet, not
recomputed) reproduces recorded runs' metrics json exactly, and times one scoring. Also checks lib553.scaled (cell-eval2's
score_one at the replicate anchor) against report323's simple formula with 323's h1-anchors.json.
usage: uv run python experiments/553-r7b-harness/kascore553.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import lib553 as S  # noqa: E402

RUNS = [ROOT / "experiments/528-callprop-block/work/h1_cpm/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz",
        ROOT / "experiments/528-callprop-block/work/h1_lcap-kept/stats-base_s1_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz",
        ROOT / "experiments/537-arm-pipeline/work/h1_lcap-kept+c4h/stats-base_s2_bps_prodem_mde0.1cap0.1min0.3coarse.npz"]


def main():
    with np.load(S.E323 / "work/reference.npz", allow_pickle=False) as z:
        ref = {k: z[k] for k in ("genes", "targets", "perts", "sums", "counts", "jk", "keep")}
    de = pl.read_parquet(S.E323 / "work/reference_de.parquet")
    anc = json.loads((S.E323 / "work/h1-anchors.json").read_text())
    out = {}
    for p in RUNS:
        z = np.load(p, allow_pickle=True)
        t0 = time.monotonic()
        raw, rows, _ = S.score_stats(ref, ref["keep"], de, z["sums"], z["jk"], z["p"], z["q"], z["lfc"])
        dt = time.monotonic() - t0
        rec = json.loads(Path(str(p).replace("stats-", "metrics-").replace(".npz", ".json")).read_text())["raw"]
        d = {m: abs(raw[m] - rec[m]) for m in S.MEMBERS}
        sc, ss = S.scaled(raw, anc["baseline"], anc["replicate"]), S.scaled_simple(raw, anc["baseline"], anc["replicate"])
        out[str(p.relative_to(ROOT))] = dict(abs_d=d, seconds=round(dt, 1), scaled=sc,
                                             scaled_vs_simple_max=max(abs(sc[k] - ss[k]) for k in sc))
        print(f"{p.parent.name}/{p.name}: max |d| {max(d.values()):.3e} ({dt:.1f} s); scaled overall {sc['overall']:.6f}, "
              f"|scaled - simple| max {out[str(p.relative_to(ROOT))]['scaled_vs_simple_max']:.2e}", flush=True)
    # cell-eval2's nmae dispatch is non-deterministic at the last ulp (repeated identical calls: 3 values in 6)
    ok = all(max(v["abs_d"].values()) <= 1e-12 for v in out.values())
    S.write_json(HERE / "work/h1ka/kascore553.json", dict(pass_=ok, runs=out, peak_rss_gb=S.rss_gb()))
    print(f"KA-S0 {'PASS' if ok else 'FAIL'}; peak RSS {S.rss_gb():.2f} GB")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

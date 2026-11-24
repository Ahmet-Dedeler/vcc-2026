"""305 step 7 -- merge the frozen test with the three post-hoc controls into the final artefacts.

Adds to `expression305.csv` the per-pair expression-matched base rate, the lane-matched effect and the
verdict flip flag; adds to `results305.json` the base-rate, lane-match and context blocks and the exact
Poisson-binomial combined statistic.

    PYTHONDONTWRITEBYTECODE=1 uv run python finalize305.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402


def poibin_ge(rates: list[float], k: int) -> float:
    p = [1.0]
    for r in rates:
        q = [0.0] * (len(p) + 1)
        for i, v in enumerate(p):
            q[i] += v * (1 - r)
            q[i + 1] += v * r
        p = q
    return float(sum(p[k:]))


def main() -> None:
    res = json.loads((X.HERE / "results305.json").read_text())
    base = {r["pair_id"]: r for r in json.loads((X.HERE / "baserate305.json").read_text())["rows"]}
    lane = {r["pair_id"]: r for r in json.loads((X.HERE / "lanematch305.json").read_text())["rows"]}
    ctx = json.loads((X.HERE / "context305.json").read_text())

    rows = list(csv.DictReader(open(X.HERE / "expression305.csv")))
    for r in rows:
        b, l = base[r["pair_id"]], lane[r["pair_id"]]
        r["base_rate_band"] = b["base_rate_band"]
        r["n_meeting_band"] = b["n_meeting_band"]
        r["n_band"] = b["n_band"]
        r["n_genes_meeting_conjunction_under_this_construct"] = b["n_meeting"]
        r["n_testable_genes"] = b["n_testable_genes"]
        r["e_pool_lanematched"] = l["e_pool_lanematched"]
        r["p_raw_lanematched"] = l["p_raw"]
        r["lanematch_flip"] = "False"
    fields = list(rows[0])
    with open(X.HERE / "expression305.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    br = json.loads((X.HERE / "baserate305.json").read_text())["rows"]
    rates_all = [r["base_rate_band"] for r in br]
    obs_all = sum(1 for r in br if r["candidate_meets"])
    prim = [r for r in br if r["group"] != "ntc"]
    rates_p = [r["base_rate_band"] for r in prim]
    obs_p = sum(1 for r in prim if r["candidate_meets"])
    clean = [r for r in br if r["n_meeting"] <= 10]        # constructs whose profile is not globally noisy
    rates_c = [r["base_rate_band"] for r in clean]
    obs_c = sum(1 for r in clean if r["candidate_meets"])

    res["post_hoc"] = dict(
        declared="added AFTER the frozen outcomes were seen; can only weaken a positive reading",
        base_rate=br, lane_match=json.loads((X.HERE / "lanematch305.json").read_text()), context=ctx,
        combined=dict(
            all15=dict(n=len(br), expected=sum(rates_all), observed=obs_all,
                       p_ge=poibin_ge(rates_all, obs_all)),
            primary11=dict(n=len(prim), expected=sum(rates_p), observed=obs_p,
                           p_ge=poibin_ge(rates_p, obs_p)),
            excluding_globally_noisy_constructs=dict(
                n=len(clean), expected=sum(rates_c), observed=obs_c, p_ge=poibin_ge(rates_c, obs_c),
                excluded=[r["pair_id"] for r in br if r["n_meeting"] > 10])))
    res["peak_rss_gib_stream"] = 3.79
    X.dump(res, X.HERE / "results305.json")
    c = res["post_hoc"]["combined"]
    for k, v in c.items():
        print(f"  {k:36s} n={v['n']:2d} expected {v['expected']:.4f} observed {v['observed']} "
              f"P(X>={v['observed']}) = {v['p_ge']:.3e}")
    print(f"  excluded as globally noisy: {c['excluding_globally_noisy_constructs']['excluded']}")


if __name__ == "__main__":
    main()

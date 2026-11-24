#!/usr/bin/env python3
"""553 — side-by-side of two or more access553.json reports (e.g. h1-test vs the R7b-shaped rehearsal h1r7): per arm the
overall-equivalent paired effect vs `cpm`, its seed sd, split-half sd (full-truth estimate), band sd, the target-bootstrap
95 % CI, the effect in band-sd units, and the raw member effects that move (nmae / fid / reach / jac); then the anchors and
the access's cost. Reads only access553.json files (outputs of an access already run).

usage: uv run python experiments/553-r7b-harness/cmp553.py work/access553_h1 work/access553_h1r7 [...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(d: str) -> tuple[str, dict]:
    p = Path(d)
    if not p.is_absolute():
        p = HERE / p
    return p.name.replace("access553_", ""), json.loads((p / "access553.json").read_text())


def main():
    reps = [load(d) for d in sys.argv[1:]]
    arms = []
    for _, r in reps:
        arms += [a for a in r["arms"] if a != "cpm" and a not in arms]
    print("| arm | surface | overall-equiv effect | seed sd | split-half sd (full-est) | band sd | effect / band sd | "
          "target-bootstrap 95 % | raw nmae / fid / reach / jac | band sd, cells + controls |")
    print("|---|---|---:|---:|---:|---:|---:|---|---|---:|")
    for arm in arms:
        for name, r in reps:
            row = r["arms"].get(arm)
            if row is None:
                continue
            o = row["effect_scaled"]["overall"]
            raw = row["effect"]
            ci = row.get("overall_boot_ci95", [float("nan")] * 2)
            mem = " / ".join(f"{raw[m]['mean']:+.4f}" for m in ("nmae", "fid", "reach", "jac"))
            bc = o.get("band_sd_ctrl")
            print(f"| `{arm}` | {name} | {o['mean']:+.5f} | {o['seed_sd']:.5f} | {o['split_half_sd_full_truth_est']:.5f} | "
                  f"{o['band_sd']:.5f} | {o['mean'] / o['band_sd']:+.1f} | [{ci[0]:+.4f}, {ci[1]:+.4f}] | {mem} | "
                  f"{'' if bc is None else f'{bc:.5f}'} |")
    print()
    print("| surface | targets | controls | truth cells / target (min / median / max) | anchors baseline (pds / mse / nmae / fid / reach / jac) | "
          "anchors replicate | `cpm` overall (scaled) | access seconds | peak RSS GB (max phase) |")
    print("|---|---:|---:|---|---|---|---:|---:|---:|")
    for name, r in reps:
        ref = r["reference"]
        c = sorted(ref["truth_cells_per_target"].values())
        b = r["anchors"]["baseline"]
        rp = r["anchors"]["replicate"]
        fmt = lambda a: " / ".join(f"{a[m]:.3f}" for m in ("pds", "mse", "nmae", "fid", "reach", "jac"))  # noqa: E731
        secs = sum(float(p.get("seconds", 0.0)) for p in r["phases"].values())
        peak = max(float(p.get("peak_rss_gb", 0.0)) for p in r["phases"].values())
        print(f"| {name} | {ref['n_targets']} | {ref['n_controls']} | {c[0]} / {c[len(c) // 2]} / {c[-1]} | {fmt(b)} | {fmt(rp)} | "
              f"{r['arms']['cpm']['scaled_mean']['overall']:.4f} | {secs:.0f} | {peak:.2f} |")


if __name__ == "__main__":
    main()

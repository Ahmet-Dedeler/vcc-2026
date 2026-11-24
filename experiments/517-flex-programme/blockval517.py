#!/usr/bin/env python3
"""517 — validate the ANALYTIC block-prior swap (price517.Ctx.block: intent change added to the delivered means, calls /
p / q / order held) against REAL emitter runs that changed only the block prior (462l, 467, 470). If the analytic path
reproduces their fid / nmae deltas, it can price block-prior variants (programme dose w, protect_q) before any heavy run.

Real runs (same seeds, same harness, only the prior differs; deltas vs 435's c4 metrics):
  arcfull (462l)  Arc on every target        gated33/50/67 (467)  Arc on strong, CLEAN4 on weak (q 0.33/0.50/0.67)
  wdown (470)     Arc strong, all-down weak  sarcs (470)          Arc strong-half programme on strong, CLEAN4 weak
usage: uv run python experiments/517-flex-programme/blockval517.py      (< 2 GB)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P  # noqa: E402
RC = P.RC

TAG = P.TAG
REAL = {"arcfull": ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull",
        "gated33": ROOT / "experiments/467-strength-gated-programme/work/h1_gated33",
        "gated50": ROOT / "experiments/467-strength-gated-programme/work/h1_gated50",
        "gated67": ROOT / "experiments/467-strength-gated-programme/work/h1_gated67",
        "wdown": ROOT / "experiments/470-weak-target-prior/work/h1_wdown"}


def main():
    t0 = time.monotonic()
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    SRC, _ = P.load_sources(genes_full, keep)
    c4 = RC.prior_signs_full(genes_full)[keep]
    arc = P.sign_from_P(SRC["arc"])
    st = json.loads((ROOT / "experiments/467-strength-gated-programme/work/strength467.json").read_text())
    order = st["targets_in_order"]
    vals = np.array([st["per_target"][str(i)][st["chosen_measure"]] for i in range(len(order))], float)
    name2v = {order[i]: vals[i] for i in range(len(order))}
    sv = np.array([name2v[t] for t in T])

    def gated(q):
        thr = float(np.quantile(vals, q))
        return lambda i, t: arc if sv[i] >= thr else c4
    arms = {"arcfull": lambda i, t: arc, "gated33": gated(0.33), "gated50": gated(0.50), "gated67": gated(0.67),
            "wdown": lambda i, t, _g=gated(0.33): _g(i, t) if sv[i] >= float(np.quantile(vals, 0.33)) else -np.ones_like(c4)}
    out = {}
    for seed in (0, 1, 2):
        ctx = P.Ctx("c4", seed, lambda i, t: c4)
        base_raw = ctx.raw_o
        rng = np.random.default_rng(5170 + seed)
        for name, fn in arms.items():
            lfc_new, rec = ctx.block(fn)
            sc, _ = ctx.score(lfc_new, rng, n_boot=200)
            f = REAL[name] / f"metrics-{TAG.format(s=seed)}.json"
            real = json.loads(f.read_text())["raw"] if f.exists() else None
            row = {"pred": {m: sc[m] for m in P.SP}, "pred_ov": sc["overall_equiv"],
                   "real": ({m: real[m] - base_raw[m] for m in P.SP} if real else None)}
            if real:
                row["real_ov"] = sum(row["real"][m] / P.SP[m] for m in P.SP) / 6
            out.setdefault(name, {})[str(seed)] = row
            print(f"s{seed} {name:8s} pred " + " ".join(f"{m} {sc[m]:+.5f}" for m in P.SP) + f" ov {sc['overall_equiv']:+.5f}"
                  + ("" if not real else " | real " + " ".join(f"{m} {row['real'][m]:+.5f}" for m in P.SP)
                     + f" ov {row['real_ov']:+.5f}") + f" | flips vs base {np.median(rec['flips_vs_base']):.0f}", flush=True)
    (HERE / "work/blockval517.json").write_text(json.dumps(out, indent=1))
    print(f"done {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()

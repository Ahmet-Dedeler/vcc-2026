#!/usr/bin/env python3
"""490 F-B: what do the K ballast cells do to the DE table the scorer reads, relative to the carriers' intent?

The K-class emitter designs every DE quantity on the n_car = 400 - K carriers: the prescribed call set and the z_cap
complement (Wilcoxon on n_car cells), the arithmetic-mean targets (prior sign x magnitude, MDE add). The scorer reads
all 400 cells, K of them ballast at ~1e6 counts carrying the BULK composition. This script compares, on 435's H1 runs:

  delivered   stats-<tag>.npz lfc (the 400-cell block, what the harness scored)
  intended    the carriers' own arithmetic mean at the recorded doses (recon_h1.py, exact)

with the call set held at the delivered one (p, q from the 400-cell block), so only the lfc channel differs.
Known-answer row: `delivered` must reproduce the harness's per-target CSV and metrics JSON.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score_de  # noqa: E402

ROOT = HERE.parents[1]
H1 = ROOT / "experiments/435-zs-on-coarse/work/h1"


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse"
    z = np.load(HERE / "work" / f"recon_{tag}.npz", allow_pickle=False)
    d = {k: z[k] for k in z.files}
    T, G = d["targets"].astype(str), d["genes_gated"].astype(str)
    # ---- known answer: the delivered tables reproduce the harness's own numbers -------------------------------------
    raw_del, per_del = score_de.score(T, G, d["p_del"], d["q_del"], d["lfc_del"])
    want = json.loads((H1 / f"metrics-{tag}.json").read_text())["raw"]
    ka = {m: abs(raw_del[m] - want[m]) for m in raw_del}
    print("KA delivered vs harness metrics JSON |diff|:", ka, flush=True)
    if max(ka.values()) > 1e-9:
        raise AssertionError("known answer failed")
    # ---- intended carriers' lfc, same calls/order --------------------------------------------------------------------
    raw_int, per_int = score_de.score(T, G, d["p_del"], d["q_del"], d["lfc_score_int"])
    print("delivered:", {k: round(v, 6) for k, v in raw_del.items()})
    print("intended :", {k: round(v, 6) for k, v in raw_int.items()})
    print("intended - delivered:", {k: round(raw_int[k] - raw_del[k], 6) for k in raw_del})
    # ---- where the delivered sign differs from the intended one -----------------------------------------------------
    own = np.equal.outer(T, G)
    ld, li = d["lfc_del"], d["lfc_score_int"]
    called = d["q_del"] < 0.05
    head, blk, prot, flip = d["head"], d["blk"], d["protected"], d["flipped"]
    sig_real = (d["q_real"] < 0.05) & ~own & np.isfinite(d["lfc_real"]) & (d["lfc_real"] != 0)
    diff = (np.sign(ld) != np.sign(li)) & ~own
    classes = {"head": head, "block protected": blk & prot, "block flipped": blk & flip & ~prot,
               "block kept": blk & ~flip & ~prot, "complement": ~head & ~blk}
    rows = {}
    for name, m in classes.items():
        mm = m & ~own
        rows[name] = dict(genes_per_target=float(mm.sum(1).mean()),
                          sign_changed_per_target=float((diff & mm).sum(1).mean()),
                          sign_changed_frac=float((diff & mm).sum() / max(mm.sum(), 1)),
                          on_truth_sig_per_target=float((diff & mm & sig_real).sum(1).mean()),
                          acc_intended_on_truth_sig=float(((np.sign(li) == np.sign(d["lfc_real"])) & mm & sig_real).sum()
                                                          / max((mm & sig_real).sum(), 1)),
                          acc_delivered_on_truth_sig=float(((np.sign(ld) == np.sign(d["lfc_real"])) & mm & sig_real).sum()
                                                           / max((mm & sig_real).sum(), 1)),
                          median_abs_dlfc=float(np.median(np.abs(ld - li)[mm])),
                          q90_abs_dlfc=float(np.quantile(np.abs(ld - li)[mm], 0.9)))
        print(name, {k: round(v, 5) for k, v in rows[name].items()}, flush=True)
    # called-set fid-style precision: delivered vs intended on the delivered call set
    cs = called & ~own & np.isfinite(d["lfc_real"]) & (d["lfc_real"] != 0)
    prec_d = ((np.sign(ld) == np.sign(d["lfc_real"])) & cs).sum() / cs.sum()
    prec_i = ((np.sign(li) == np.sign(d["lfc_real"])) & cs).sum() / cs.sum()
    print(f"pooled called-set sign precision: delivered {prec_d:.5f} intended {prec_i:.5f}")
    out = dict(tag=tag, known_answer_abs_diff=ka, raw_delivered=raw_del, raw_intended=raw_int,
               delta_intended_minus_delivered={k: raw_int[k] - raw_del[k] for k in raw_del}, classes=rows,
               pooled_called_precision=dict(delivered=float(prec_d), intended=float(prec_i)))
    (HERE / "work" / f"ballast_footprint_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

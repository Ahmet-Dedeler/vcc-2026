"""250 Phase 3b — cross-panel PINNED bases: deflate panel P's CD4 residual (fold-level) with
the PRODUCTION-shaped V_k of every OTHER panel (from basis250-production.npz), resolved sizing.
Answers: is the A/B/C basis special, or does any source-level basis transfer?

    /usr/bin/time -l uv run python cross250.py --panel efr --fold k562_gwps
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gate250 import wait_for_memory  # noqa: E402
import lib250 as L  # noqa: E402

PANELS = ("abc", "real", "pess", "efp", "efr")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, choices=PANELS)
    ap.add_argument("--fold", required=True, choices=("k562_gwps", "orion_hct116", "orion_hek293t"))
    args = ap.parse_args()
    panel, held = args.panel, args.fold
    sizing = "allrows" if panel == "abc" else "carrying"
    t0 = time.time()
    print(f"[{panel}/{held}] gate clear: {wait_for_memory(tag=f' cross {panel}/{held}')}", flush=True)
    st = L.stage(panel, held)
    cd4, whitened = st["cd4"], st["whitened"]
    bases = np.load(L.WORK / "basis250-production.npz")
    ref = {r["arm"]: r for r in json.loads(
        (L.WORK / f"sweep250-{panel}-{held}.json").read_text())["rows"]}
    rows = []
    per_target = {}

    def run(name, cd4_res, frac, audit, ka=None):
        final, preserve = L.compose(st, cd4_res, frac)
        e = L.score(st, final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        row = dict(panel=panel, fold=held, arm=name, n_targets=len(st["targets"]), sizing=sizing,
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   audit=dict(frac_cd4=float(frac), preserve=preserve, **audit))
        msg = f"[{panel}/{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f}"
        if ka is not None:
            r = ref[ka]
            d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["pds"] - e["pds"]))
            assert d < 1e-12, f"KA FAILED {name} vs sweep {ka}: {d:.3e}"
            msg += f"   KA PASS vs sweep250 {ka} ({d:.1e})"
        rows.append(row)
        print(msg, flush=True)

    frac0, _ = L.solve_cd4(whitened, cd4, sizing)
    run("BASE", cd4, frac0, dict(family="base"), ka="BASE")
    # the pinned A/B/C row must reproduce the sweep's (same function, same basis)
    for k in (5, 10):
        for src in PANELS:
            V = bases[f"V{k}_{src}"]
            d, aud = L.deflate_gene_basis(cd4, V, src)
            fr, pu = L.solve_cd4(whitened, d, sizing)
            run(f"pin{src}{k}", d, fr, dict(family="cross_pin", basis_panel=src,
                                            sizing_mode="resolved", per_unit=float(pu), **aud),
                ka=(f"pinABC{k}" if src == "abc" else None))
    (L.WORK / f"cross250-{panel}-{held}.json").write_text(json.dumps(
        dict(meta=dict(panel=panel, fold=held, sizing=sizing, seconds=round(time.time() - t0, 1),
                       peak_rss_gib=L.peak_gib()), rows=rows), indent=1))
    np.savez_compressed(L.WORK / f"pertarget250x-{panel}-{held}.npz", **per_target)
    print(f"[{panel}/{held}] cross done: {len(rows)} rows in {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()

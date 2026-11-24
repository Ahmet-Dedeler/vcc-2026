#!/usr/bin/env python3
"""537 — are the preregistered noise pairs (PREREGISTRATION-noise.md) really bulk-only? Diff each pair's package config
(submissions/<stem>.h5ad.meta.json, read-only) and flag any key outside the bulk set. `sameline` is allowed only when no
context carries a programme and no row changed (the file sha may differ; 548 §2.4: that pair shape moves DE by ~5e-5).
Run before the reads land; the result goes to RESULTS §2.5. < 0.2 GB.
usage: uv run python experiments/537-arm-pipeline/pairmeta537.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
S = ROOT / "submissions"
BULK = {"bulk_ctx_scale", "bulk_load", "bulk_load_sha256", "bulk_sharpen_q", "bulk_sharpen_pt", "bulk_sharpen_pt_sha256",
        "bulk_lodo", "bulk_cap", "pert_scale"}
PAIRS = [
    ("460dr − 453dr", "460dr-october-path-lodo-abc-q1.75", "453dr-october-path-abc-q1.75"),
    ("479p2g − 446bj2", "479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5",
     "446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5"),
    ("479sp2 − 446bj2", "479sp2-cap15-partial50-clean4q80-pertscale3-ballast-sharpensp2-stacklodo05load-ctxlevelr-mde010coarse-zr5",
     "446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5"),
    ("495 − 446bj2", "495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5",
     "446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5"),
    ("471c − 446b", "471c-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevel-mde010coarse-zr5",
     "446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5"),
    ("498p2g − 495", "498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5",
     "495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5"),
    ("498sp2 − 495", "498sp2-cap15-partial50-clean4q80-lodo6-ballast-sharpensp2-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5",
     "495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5"),
    ("497dr − 460dr", "497dr-october-path-lodopct-spt-abc", "460dr-october-path-lodo-abc-q1.75"),
]


def sameline_inert(v):
    if not isinstance(v, dict):
        return v in (None, {}, False)
    ctx = v.get("contexts", {})
    return all(not c.get("programme_default") and not c.get("rows_changed") for c in ctx.values())


def main():
    for lab, a, b in PAIRS:
        fa, fb = S / f"{a}.h5ad.meta.json", S / f"{b}.h5ad.meta.json"
        if not (fa.exists() and fb.exists()):
            print(f"{lab:16s} meta missing: {[str(f.name) for f in (fa, fb) if not f.exists()]}")
            continue
        ca, cb = json.loads(fa.read_text())["config"], json.loads(fb.read_text())["config"]
        diff = sorted(k for k in set(ca) | set(cb) if json.dumps(ca.get(k), sort_keys=True) != json.dumps(cb.get(k), sort_keys=True))
        bad = [k for k in diff if k not in BULK and not (k == "sameline" and sameline_inert(ca.get(k)) and sameline_inert(cb.get(k)))]
        print(f"{lab:16s} differing keys {diff} -> {'BULK-ONLY' if not bad else 'NOT bulk-only: ' + str(bad)}")
        for k in bad:
            print(f"      {k}: {json.dumps(cb.get(k))[:160]} -> {json.dumps(ca.get(k))[:160]}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""554 -- calibrate the dose guard on everything on disk: control-library spread per context (the preflight floor), the
per-context steering dose from every per-block rc_audit dump (the postbuild floors, with the dump matched to its meta), and
the pooled medians of every listed meta (the postbuild's no-dump fallback). Prints each value's margin to its floor and
writes work/calib554.json. Reads count matrices of CONTROL cells only (and h1-test's controls); light (< 3 GB).

    uv run python experiments/554-dose-guard/calib554.py
"""
from __future__ import annotations

import hashlib
import json
import statistics
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPREAD_MAX, PRIOR_MIN, MDE_MIN, MISSED_MAX = 3.0, 0.9, 0.9, 100.0

BUNDLES = {                                                    # name -> bundle dir (contexts from its manifest)
    "abc (data/challenge)": "data/challenge",
    "180 def-bundle (A/B/C relabelled)": "experiments/180-def-readiness/work/def-bundle",
    "199 def-bundle-r": "experiments/199-def-runbook/work/def-bundle-r",
    "238 ef-bundle (efp)": "experiments/238-essential-free-panel/work/ef-bundle",
    "238 ef-bundle-r (efr)": "experiments/238-essential-free-panel/work/ef-bundle-r",
    "534 R6 S1": "experiments/534-generalization-audit/work/partB/S1/bundle",
    "552 R7a T1": "experiments/552-r7-reserve/work/T1/bundle",
}
EXTRA_CONTROLS = {                                             # single control files outside a 3-context bundle
    "H1 (h1-test controls)": ("data/harness/h1-test.h5ad", "non-targeting"),
    "R7b as built (Calu-3)": ("experiments/552-r7-reserve/work/r7b/bundle/context_R.h5ad", None),
    "R7b floored 5,000": ("experiments/552-r7-reserve/work/r7b_floor5000/bundle/context_R.h5ad", None),
}
DUMPS = {                                                      # dump -> (meta, what it is)
    "experiments/429-zs-production/work/rc_audit_429a.jsonl":
        ("submissions/429a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel-zs250.h5ad.meta.json",
         "429a: plain MDE (the KNOWN FAILURE: A's MDE dose 0)"),
    "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl":
        ("submissions/433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse.h5ad.meta.json",
         "433a: coarse MDE (the production design since)"),
    "experiments/432-def-rehearsal-new-knobs/work/b0/rc_audit_b0.jsonl":
        ("experiments/432-def-rehearsal-new-knobs/work/b0/b0.h5ad.meta.json", "432 b0: October rehearsal (D/E/F = A/B/C)"),
    "experiments/432-def-rehearsal-new-knobs/work/b1/rc_audit_b1.jsonl":
        ("experiments/432-def-rehearsal-new-knobs/work/b1/b1.h5ad.meta.json", "432 b1"),
    "experiments/432-def-rehearsal-new-knobs/work/b2/rc_audit_b2.jsonl":
        ("experiments/432-def-rehearsal-new-knobs/work/b2/b2.h5ad.meta.json", "432 b2"),
    "experiments/432-def-rehearsal-new-knobs/work/ef_b1/rc_audit_ef_b1.jsonl":
        ("experiments/432-def-rehearsal-new-knobs/work/ef_b1/ef_b1.h5ad.meta.json", "432 ef_b1 (efp panel)"),
    "experiments/448-def-rehearsal-q175/work/efr/rc_audit_lvl1.jsonl": (None, "448 efr lvl1 (efr panel, q 1.75)"),
}
METAS = [
    "submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json",
    "submissions/446bj2zcp-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cp.h5ad.meta.json",
    "submissions/446bj2zcph-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph.h5ad.meta.json",
    "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json",
    "submissions/457dr-october-path-abc-cd4x15-q1.75.h5ad.meta.json",
    "submissions/460dr-october-path-lodo-abc-q1.75.h5ad.meta.json",
    "submissions/497dr-october-path-lodopct-spt-abc.h5ad.meta.json",
    "experiments/526-cloud-october-review2/work/efp_inert/474efp_final.h5ad.meta.json",
    "experiments/497-october-switches/work/efr497/497efr_final.h5ad.meta.json",
    "experiments/510-october-recipe-rehearsal/work/efr510/510efr_final.h5ad.meta.json",
    "experiments/452-def-day-rehearsal/work/efr/def452_final.h5ad.meta.json",
    "experiments/474-october-hardening/work/efr_retrace/465efr_final.h5ad.meta.json",
]


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def labels(h, path):
    g = h[path]
    if isinstance(g, h5py.Group):
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]], dtype=object)
        return cats[g["codes"][:]]
    v = g[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object)


def libraries(path: Path, ctrl_label: str | None) -> np.ndarray:
    """Per-cell library sizes (CSR row sums, streamed); with ctrl_label, only the rows carrying that target_gene."""
    with h5py.File(path, "r") as h:
        X = h["X"]
        ip = X["indptr"][:].astype(np.int64)
        n = ip.size - 1
        keep = np.ones(n, bool) if ctrl_label is None else (labels(h, "obs/target_gene") == ctrl_label)
        lib = np.zeros(n)
        d = X["data"]
        step = 1 << 24
        for s in range(0, int(ip[-1]), step):
            blk = d[s:s + step].astype(np.float64)
            rows = np.searchsorted(ip, np.arange(s, s + blk.size), "right") - 1
            np.add.at(lib, rows, blk)
        return lib[keep]


def spread(lib: np.ndarray) -> dict:
    q25, q50, q75 = (float(v) for v in np.quantile(lib, [0.25, 0.5, 0.75]))
    s = q75 / q25
    return dict(n=int(lib.size), q25=round(q25), q50=round(q50), q75=round(q75), spread=round(s, 4),
                margin=round(SPREAD_MAX - s, 4), trips=bool(s > SPREAD_MAX))


def dose_from_dump(rows: list[dict]) -> dict:
    out = {}
    for c in sorted({str(r["context"]) for r in rows}):
        sub = [r for r in rows if str(r["context"]) == c]
        v = {k: (float(np.median([float(r[k]) for r in sub if k in r])) if any(k in r for r in sub) else None)
             for k in ("prior_lambda", "mde_lambda", "call_set_missed")}
        m = {"prior_lambda": None if v["prior_lambda"] is None else round(v["prior_lambda"] - PRIOR_MIN, 4),
             "mde_lambda": None if v["mde_lambda"] is None else round(v["mde_lambda"] - MDE_MIN, 4),
             "call_set_missed": None if v["call_set_missed"] is None else round(MISSED_MAX - v["call_set_missed"], 2)}
        out[c] = dict(blocks=len(sub), **{k: (round(x, 4) if x is not None else None) for k, x in v.items()}, margin=m,
                      trips=any(x is not None and x < 0 for x in m.values()))
    return out


def main() -> None:
    res = {"floors": dict(spread=SPREAD_MAX, prior=PRIOR_MIN, mde=MDE_MIN, missed=MISSED_MAX), "spread": {}, "dumps": {},
           "metas": {}}
    print(f"== control-library q75/q25 per context (floor {SPREAD_MAX}) ==")
    seen = {}
    for name, d in BUNDLES.items():
        b = ROOT / d
        man = json.loads((b / "manifest.json").read_text())
        res["spread"][name] = {}
        for c in man["contexts"]:
            f = b / f"context_{c}.h5ad"
            key = sha(f)
            if key not in seen:
                seen[key] = spread(libraries(f, None))
            res["spread"][name][c] = dict(seen[key], sha256=key[:16])
            r = res["spread"][name][c]
            print(f"  {name:36s} {c}  n {r['n']:>6}  q25 {r['q25']:>6}  q50 {r['q50']:>6}  q75 {r['q75']:>6}  "
                  f"q75/q25 {r['spread']:.3f}  margin {r['margin']:+.3f}  {'TRIPS' if r['trips'] else 'ok'}  sha {r['sha256']}")
    for name, (f, lab) in EXTRA_CONTROLS.items():
        p = ROOT / f
        if not p.exists():
            print(f"  {name:36s} (missing: {f})")
            continue
        r = spread(libraries(p, lab))
        res["spread"][name] = {"controls": r}
        print(f"  {name:36s} -  n {r['n']:>6}  q25 {r['q25']:>6}  q50 {r['q50']:>6}  q75 {r['q75']:>6}  "
              f"q75/q25 {r['spread']:.3f}  margin {r['margin']:+.3f}  {'TRIPS' if r['trips'] else 'ok'}")
    print(f"\n== per-context steering dose from per-block dumps (floors: prior >= {PRIOR_MIN}, mde >= {MDE_MIN}, "
          f"missed <= {MISSED_MAX:g}) ==")
    for dpath, (mpath, what) in DUMPS.items():
        rows = [json.loads(ln) for ln in (ROOT / dpath).read_text().splitlines() if ln.strip()]
        match = None
        if mpath and (ROOT / mpath).exists():
            rc = json.loads((ROOT / mpath).read_text())["output"].get("rc_callset", {})
            match = len(rows) == rc.get("blocks") and all(
                float(np.median([float(r[k]) for r in rows if k in r])) == float(rc[k]["median"])
                for k in ("prior_lambda", "mde_lambda", "call_set_missed") if k in rc)
        per = dose_from_dump(rows)
        res["dumps"][dpath] = dict(what=what, meta=mpath, dump_matches_meta=match, per_context=per)
        print(f"  {what}  [dump == meta pooled: {match}]")
        for c, v in per.items():
            print(f"     {c}: prior {v['prior_lambda']}  mde {v['mde_lambda']}  missed {v['call_set_missed']}  "
                  f"margins {v['margin']}  {'TRIPS' if v['trips'] else 'ok'}")
    print("\n== pooled medians of the listed metas (the no-dump fallback) ==")
    for m in METAS:
        p = ROOT / m
        if not p.exists():
            print(f"  (missing {m})")
            continue
        rc = json.loads(p.read_text())["output"].get("rc_callset", {})
        v = {k: (rc[k]["median"] if k in rc else None) for k in ("prior_lambda", "mde_lambda", "call_set_missed")}
        v["min_prior"], v["min_mde"] = rc.get("prior_lambda", {}).get("min"), rc.get("mde_lambda", {}).get("min")
        trips = (v["prior_lambda"] is not None and v["prior_lambda"] < PRIOR_MIN) or \
                (v["mde_lambda"] is not None and v["mde_lambda"] < MDE_MIN) or \
                (v["call_set_missed"] is not None and v["call_set_missed"] > MISSED_MAX)
        res["metas"][m] = dict(blocks=rc.get("blocks"), **v, trips=bool(trips))
        print(f"  {Path(m).name[:60]:60s} blocks {rc.get('blocks')}  prior {v['prior_lambda']}  mde {v['mde_lambda']}  "
              f"missed {v['call_set_missed']}  (pooled min prior {v['min_prior']}, mde {v['min_mde']})  "
              f"{'TRIPS' if trips else 'ok'}")
    out = HERE / "work" / "calib554.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=float))
    print(f"\n-> {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

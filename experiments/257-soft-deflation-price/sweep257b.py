"""257b -- EXPLORATORY, post hoc, NOT part of the pre-fixed grid and NOT used for the verdict.

The per-direction cells (sweep257) said direction 1 carries half the `pds` at no aggregate-cosine
cost while direction 4 alone costs -0.41 % and buys nothing.  The directions are not additive, so
the complementary question is measured directly: the k = 5 projector with ONE direction LEFT IN
(`k5_drop{j}`: dirs = top-5 minus j, s = 1), resolved sizing.  Same known-answer rows as sweep257.

    uv run python sweep257b.py --fold <held>
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib257 as X  # noqa: E402

L = X.L
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    ref_a = {r["arm"]: r for r in json.loads((X.E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    ref_b = {r["arm"]: r for r in json.loads((X.E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}
    st = L.stage(held)
    ident = X.Identity(st)
    rows, per_target, identity, failures = [], {}, {}, []

    def run(name, parts, extra=None, ka=None, ka_ref=None):
        final = L.compose(st, parts)
        e = L.score(st, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        row = dict(fold=held, arm=name, **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                   delivered_norm=m["delivered_norm"], audit=extra or {})
        msg = f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f}"
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(d < TOL))
            msg += f"   KA {'PASS' if d < TOL else 'FAIL'} vs {ka} ({d:.1e})"
            if d >= TOL:
                failures.append(f"{name} vs {ka}")
        rows.append(row)
        print(msg, flush=True)

    run("STACK", L.stack_parts(st), dict(frac_cd4=st["frac_cd4"], frac_x=st["frac_x"]), "sweep246 STACK", ref_a["STACK"])
    d, aud = X.deflate_shrink(st["cd4"], 5, 1.0)
    parts, fr, pu = X.resolved_parts(st, d)
    run("k5_s1", parts, dict(**aud, fraction=fr, per_unit=pu), "B2_cd4_svd5_resolved", ref_b["B2_cd4_svd5_resolved"])
    for j in range(5):
        dirs = [i for i in range(5) if i != j]
        d, aud = X.deflate_shrink(st["cd4"], 5, 1.0, dirs=dirs)
        parts, fr, pu = X.resolved_parts(st, d)
        run(f"k5_drop{j+1}", parts, dict(**aud, fraction=fr, per_unit=pu, sizing="resolved"))
    meta = dict(fold=held, failures=failures, seconds=round(time.time() - t0, 1),
                peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3)
    (X.WORK / f"sweep257b-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=float))
    np.savez_compressed(X.WORK / f"pertarget257b-{held}.npz", **per_target)
    np.savez_compressed(X.WORK / f"identity257b-{held}.npz", **identity)
    print(f"[{held}] {len(rows)} rows, {len(failures)} KA failures, {meta['seconds']}s, peak {meta['peak_rss_gib']:.2f} GiB", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

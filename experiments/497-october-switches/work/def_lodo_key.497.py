"""460 — 415's LODO-V1 per-target bulk key at clip (0.5, 3) (`lodo415_lo0.5_hi3`, 420's artifact `8c82bef3`, the key
446bj2 carried into its stack dump on A/B/C) re-fitted for ANY panel: the October day-of fit.

    key_p = mean_j <s (b^{-j}_p + radd_p), d_{j,p}> / |F_p|^2      over the donors j carrying p (>= 2 carried)
    w_p   = clip(key_p / median_carried(key), lo, hi) on carried targets, 1 elsewhere

F is the h1c1 composition the key was designed on (415pp: normrel cap -1.5 + PM50 base, whitening at gamma 0.65, CD4
residual, H1 union g 0.28, CD4 partial g 0.5, one preserve), b^{-j} the whitened base rebuilt without donor j, d_j donor
j's raw bank row. The arithmetic is `experiments/438-october-stack-repricing/full438.py`'s (lib438: the panel, operator,
CD4 banks and composition RULES are arguments), which already asserts on A/B/C that this LODO-V1 equals 415pp's key to
1e-12. This script only writes it in 410's artifact format (targets, w, key, carried, n_donors) so
`def_pin_panel_artifacts.py --pert-scale <npz> --pert-scale-id lodo_lo0.5_hi3_def` can register it on the day.

Rules: `def` on a D/E/F bundle (CD4 carrying-row sizing at g 0.3443655, deflation by the pinned V5_abc basis: what the
October builds run); `abc` reproduces the A/B/C artifact.

KNOWN ANSWER (--ka): on data/challenge with the production operator (e2f65ecf) and rules abc, the arrays equal
`experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz` (8c82bef3): targets, carried and n_donors exactly,
key and w to <= 1e-12. Any failure exits 1.

usage: uv run python scripts/def_lodo_key.py --bundle <dir> --op <whitener.npz> --cd4-root <cd4_compact> --out <npz>
           [--rules def|abc] [--lo 0.5] [--hi 3 | --hi-pct 96.64697522769474] [--tag T] [--ka | --ka-pct]
497: --hi-pct P sets hi to the P-th percentile of key / median over the carried targets (543: P 96.64697522769474 == cap 6
on A/B/C to 4e-15; the rounded 96.647 gives 6.000055;
on October-like panels it caps at 8-16). --ka-pct checks that against 495's cap-6 artifact dc374281.
Peak RSS ~8 GiB (full438's composition), a few minutes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))

import numpy as np  # noqa: E402

import lib438 as L  # noqa: E402

ABC_OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_ART = (ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz",
           "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b")
ABC_ART_PCT = (ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",     # 497: cap 6 == p96.64697522769474
               "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--op", type=Path, required=True, help="the panel's 3-donor whitener (the day's re-fit)")
    ap.add_argument("--cd4-root", type=Path, required=True, help="the panel's cd4_compact (the day's extraction)")
    ap.add_argument("--rules", default="def", choices=("def", "abc"))
    ap.add_argument("--lo", type=float, default=0.5)
    ap.add_argument("--hi", type=float, default=3.0)
    ap.add_argument("--hi-pct", type=float, default=None,
                    help="497: hi = this percentile of key/median over carried targets (overrides --hi; 96.64697522769474 == cap 6 on A/B/C)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tag", default=None, help="lib438 work subdir for the frame caches (default: lodo_<bundle name>)")
    ap.add_argument("--ka", action="store_true")
    ap.add_argument("--ka-pct", action="store_true", help="497: KA of --hi-pct 96.64697522769474 against 495's cap-6 artifact")
    a = ap.parse_args()
    t0 = time.time()
    log = lambda *x: print("[def_lodo_key]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    bundle, op_path = a.bundle.resolve(), a.op.resolve()
    if a.ka:
        assert bundle == (ROOT / "data/challenge").resolve() and a.rules == "abc", "KA runs on data/challenge, rules abc"
        assert (a.lo, a.hi) == (0.5, 3.0) and a.hi_pct is None, "KA is the (0.5, 3) artifact"
    if a.ka_pct:
        assert bundle == (ROOT / "data/challenge").resolve() and a.rules == "abc", "KA-pct runs on data/challenge, rules abc"
        assert a.lo == 0.5 and a.hi_pct == 96.64697522769474 and not a.ka, "KA-pct is (0.5, p96.64697522769474) alone"
    tag = a.tag or f"lodo_{bundle.name}_{a.rules}"
    spec = dict(bundle=str(bundle), rules=a.rules, cd4_root=str(a.cd4_root.resolve()), full_op=str(op_path),
                folds=list(L.CANCER), synthetic=False, meta=None, fold_ops={})
    ctx = L.bind(tag, spec=spec, tag=tag)
    targets = list(ctx.panel)
    pp_pos = ctx.pp_pos
    wcols = L.E064.operator_columns([L.bank(n) for n in L.MAIN], L.OFFICIAL, ctx.panel_set)
    op = L.plolo.Operator(op_path)
    op._check_axis(wcols)                     # refuses unless OFFICIAL[wcols] hashes to the operator's axis record
    assert not (set(targets) & {L.OFFICIAL[i] for i in wcols}), "panel genes on the operator axis"
    if a.ka or a.ka_pct:
        assert op.sha256 == ABC_OP_SHA, f"KA needs the production operator, got {op.sha256}"
    log(f"panel {bundle} ({ctx.contexts}); rules {a.rules}; op {op.sha256[:16]} ({op.axis_checked}); "
        f"cd4 {a.cd4_root}")

    # composition (full438's lines, h1c1 only)
    dn = L.Donors438(targets, L.MAIN)
    ncar = sum(dn.A[s].astype(int) for s in L.MAIN)
    fb = None
    if (ncar == 0).any():
        src = [L.prepared(n, targets, pp_pos) for n in L.MAIN]
        fb, _ = L.E062.build_baseline(src)
    pre = L.pre_normrel(dn, L.MAIN, "cap", -1.5, fallback_rows=fb)
    cd4, partial, h1, daud = L.residuals(ctx, targets)
    whit = op.apply(pre * L.ALPHA, wcols, L.GAMMA, L.LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos, pinned=("h1c1",), size=1)
    carried = ncar >= 2
    bases = L.lodo_bases(op, dn, L.MAIN, "cap", -1.5, pre, wcols)
    d_rows = {s: L.donor_rows(s, targets)[0] for s in L.MAIN}
    kv1 = L.lodo_key(FC.final("h1c1"), FC.parts("h1c1"), bases, d_rows, {s: dn.A[s] for s in L.MAIN}, carried,
                     L.MAIN, pp_pos)
    key = np.nan_to_num(kv1, nan=0.0)
    hi = a.hi
    if a.hi_pct is not None:
        # 543: the cap adapts to the panel's key tail (numpy's default linear percentile of key / median, carried only)
        hi = float(np.percentile(key[carried] / float(np.median(key[carried])), a.hi_pct))
        assert hi > a.lo, f"percentile cap {hi} <= floor {a.lo}"
        a.hi = hi
    w, med = L.weights(key, carried, a.lo, hi)
    log(f"key fitted: carried {int(carried.sum())}/300, median {med:.6g}, composition {json.dumps(FC.meta['h1c1'])}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.out, targets=np.asarray(targets), w=w.astype(np.float64), key=key, carried=carried,
             n_donors=ncar.astype(np.int64))
    sha = hashlib.sha256(a.out.read_bytes()).hexdigest()
    q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
    side = dict(sha256=sha, lo=a.lo, hi=a.hi, hi_pct=a.hi_pct, key="LODO-V1 (415pp) on the h1c1 composition, lib438 arithmetic",
                rules=a.rules, donors=list(L.MAIN), bundle=str(bundle), contexts=list(ctx.contexts),
                operator=str(op_path), operator_sha256=op.sha256, operator_axis=op.axis_checked,
                cd4_root=str(a.cd4_root.resolve()), n_wcols=int(wcols.size), composition=FC.meta["h1c1"],
                deflation=daud, key_median_carried=med, carried=int(carried.sum()),
                n_donor_hist={int(k): int((ncar == k).sum()) for k in np.unique(ncar)},
                key_quantiles_carried=q(key[carried]), frac_key_le0=float((key[carried] <= 0).mean()),
                w_quantiles=q(w), frac_at_floor=float((w[carried] <= a.lo + 1e-12).mean()),
                frac_at_cap=float((w[carried] >= a.hi - 1e-12).mean()),
                top10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(-key)[:10]],
                bottom10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(key)[:10]])
    ref_path, ref_sha = ABC_ART
    if sorted(targets) == sorted(str(t) for t in np.load(ref_path)["targets"]):
        # on the A/B/C panel, report the distance to the board-tested artifact whatever the rules
        z = np.load(ref_path)
        side["vs_abc_artifact"] = dict(path=str(ref_path.relative_to(ROOT)), sha256=ref_sha,
                                       max_abs_key=float(np.abs(z["key"] - key).max()),
                                       max_abs_w=float(np.abs(z["w"] - w).max()),
                                       pearson_w=float(np.corrcoef(z["w"], w)[0, 1]),
                                       bytes_identical=sha == ref_sha)
    a.out.with_suffix(".json").write_text(json.dumps(side, indent=1))
    log(f"wrote {a.out} sha256 {sha}")
    print(json.dumps({k: side[k] for k in ("w_quantiles", "frac_at_floor", "frac_at_cap", "carried", "n_donor_hist")
                      + (("vs_abc_artifact",) if "vs_abc_artifact" in side else ())}, indent=1), flush=True)
    if a.ka_pct:
        pp, ps = ABC_ART_PCT
        assert L.sha256(pp) == ps, "KA-pct reference artifact sha"
        z = np.load(pp)
        okp = ([str(t) for t in z["targets"]] == targets and np.array_equal(z["carried"], carried)
               and np.array_equal(z["n_donors"], ncar.astype(np.int64))
               and float(np.abs(z["key"] - key).max()) <= 1e-12 and float(np.abs(z["w"] - w).max()) <= 1e-12)
        print(f"KA-pct {'PASS' if okp else 'FAIL'}: hi {hi!r} vs {pp.name} ({ps[:12]}): max|dw| "
              f"{float(np.abs(z['w'] - w).max()):.2e} max|dkey| {float(np.abs(z['key'] - key).max()):.2e}", flush=True)
        if not okp:
            raise SystemExit(1)
    if a.ka:
        assert L.sha256(ref_path) == ref_sha, "KA reference artifact sha"
        z = np.load(ref_path)
        ok = ([str(t) for t in z["targets"]] == targets and np.array_equal(z["carried"], carried)
              and np.array_equal(z["n_donors"], ncar.astype(np.int64))
              and float(np.abs(z["key"] - key).max()) <= 1e-12 and float(np.abs(z["w"] - w).max()) <= 1e-12)
        print(f"KA {'PASS' if ok else 'FAIL'}: vs {ref_path.name} ({ref_sha[:12]}): "
              f"{json.dumps(side['vs_abc_artifact'])}", flush=True)
        if not ok:
            raise SystemExit(1)


if __name__ == "__main__":
    main()

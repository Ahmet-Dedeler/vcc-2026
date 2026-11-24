"""497: `scripts/def_lodo_key.py --hi-pct P` -- the October form of the LODO key's cap (543): hi = the P-th percentile
(numpy default, linear) of key / median over the panel's carried targets, instead of a fixed number. On A/B/C, P 96.647
is cap 6 (495's artifact dc374281) only at the unrounded P 96.64697522769474 (weights to 4e-15); the rounded 96.647
lands at 6.000055 (10 weights off by 5.5e-5). New KA `--ka-pct`: data/challenge, rules abc, lo 0.5,
--hi-pct 96.64697522769474 must reproduce dc374281's targets / carried / n_donors exactly and key / w to <= 1e-12.
Without --hi-pct nothing changes (the (0.5, 3) KA still holds).
usage: uv run python experiments/497-october-switches/patch_lodokey497.py <def_lodo_key.py> [--check]
"""
import sys
from pathlib import Path

path = Path(sys.argv[1]); src = path.read_text()
edits = [
    ('           [--rules def|abc] [--lo 0.5] [--hi 3] [--tag T] [--ka]\n',
     '           [--rules def|abc] [--lo 0.5] [--hi 3 | --hi-pct 96.64697522769474] [--tag T] [--ka | --ka-pct]\n'
     '497: --hi-pct P sets hi to the P-th percentile of key / median over the carried targets (543: P 96.64697522769474 == cap 6\n'
     'on A/B/C to 4e-15; the rounded 96.647 gives 6.000055;\n'
     'on October-like panels it caps at 8-16). --ka-pct checks that against 495\'s cap-6 artifact dc374281.\n'),
    ('ABC_ART = (ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz",\n'
     '           "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b")\n',
     'ABC_ART = (ROOT / "experiments/420-bulk-stack/work/pert_scale_415pp_lo0.5_hi3.npz",\n'
     '           "8c82bef37d0a59b3ed1f4917ccad0d8fd24618a89fd9e2700c058b9aeff4339b")\n'
     'ABC_ART_PCT = (ROOT / "experiments/495-lodo-cap6/work/pert_scale_415pp_lo0.5_hi6.npz",     # 497: cap 6 == p96.64697522769474\n'
     '               "dc374281452ed4f7ae3a25b392dbfa4ff0d5cc591bd7b5b71e2c0a2dca86a4c0")\n'),
    ('    ap.add_argument("--hi", type=float, default=3.0)\n',
     '    ap.add_argument("--hi", type=float, default=3.0)\n'
     '    ap.add_argument("--hi-pct", type=float, default=None,\n'
     '                    help="497: hi = this percentile of key/median over carried targets (overrides --hi; 96.64697522769474 == cap 6 on A/B/C)")\n'),
    ('    ap.add_argument("--ka", action="store_true")\n',
     '    ap.add_argument("--ka", action="store_true")\n'
     '    ap.add_argument("--ka-pct", action="store_true", help="497: KA of --hi-pct 96.64697522769474 against 495\'s cap-6 artifact")\n'),
    ('        assert (a.lo, a.hi) == (0.5, 3.0), "KA is the (0.5, 3) artifact"\n',
     '        assert (a.lo, a.hi) == (0.5, 3.0) and a.hi_pct is None, "KA is the (0.5, 3) artifact"\n'
     '    if a.ka_pct:\n'
     '        assert bundle == (ROOT / "data/challenge").resolve() and a.rules == "abc", "KA-pct runs on data/challenge, rules abc"\n'
     '        assert a.lo == 0.5 and a.hi_pct == 96.64697522769474 and not a.ka, "KA-pct is (0.5, p96.64697522769474) alone"\n'),
    ('    if a.ka:\n        assert op.sha256 == ABC_OP_SHA, f"KA needs the production operator, got {op.sha256}"\n',
     '    if a.ka or a.ka_pct:\n        assert op.sha256 == ABC_OP_SHA, f"KA needs the production operator, got {op.sha256}"\n'),
    ('    w, med = L.weights(key, carried, a.lo, a.hi)\n',
     '    hi = a.hi\n'
     '    if a.hi_pct is not None:\n'
     '        # 543: the cap adapts to the panel\'s key tail (numpy\'s default linear percentile of key / median, carried only)\n'
     '        hi = float(np.percentile(key[carried] / float(np.median(key[carried])), a.hi_pct))\n'
     '        assert hi > a.lo, f"percentile cap {hi} <= floor {a.lo}"\n'
     '        a.hi = hi\n'
     '    w, med = L.weights(key, carried, a.lo, hi)\n'),
    ('    side = dict(sha256=sha, lo=a.lo, hi=a.hi, key=',
     '    side = dict(sha256=sha, lo=a.lo, hi=a.hi, hi_pct=a.hi_pct, key='),
    ('    if a.ka:\n        assert L.sha256(ref_path) == ref_sha, "KA reference artifact sha"\n',
     '    if a.ka_pct:\n'
     '        pp, ps = ABC_ART_PCT\n'
     '        assert L.sha256(pp) == ps, "KA-pct reference artifact sha"\n'
     '        z = np.load(pp)\n'
     '        okp = ([str(t) for t in z["targets"]] == targets and np.array_equal(z["carried"], carried)\n'
     '               and np.array_equal(z["n_donors"], ncar.astype(np.int64))\n'
     '               and float(np.abs(z["key"] - key).max()) <= 1e-12 and float(np.abs(z["w"] - w).max()) <= 1e-12)\n'
     '        print(f"KA-pct {\'PASS\' if okp else \'FAIL\'}: hi {hi!r} vs {pp.name} ({ps[:12]}): max|dw| "\n'
     '              f"{float(np.abs(z[\'w\'] - w).max()):.2e} max|dkey| {float(np.abs(z[\'key\'] - key).max()):.2e}", flush=True)\n'
     '        if not okp:\n'
     '            raise SystemExit(1)\n'
     '    if a.ka:\n        assert L.sha256(ref_path) == ref_sha, "KA reference artifact sha"\n'),
]
for old, new in edits:
    if src.count(old) != 1:
        raise SystemExit(f"anchor count {src.count(old)} != 1: {old[:70]!r}")
if "--hi-pct" in src:
    raise SystemExit("already patched")
for old, new in edits:
    src = src.replace(old, new)
compile(src, str(path), "exec")
if "--check" in sys.argv:
    print("patch applies cleanly (not written)")
else:
    path.write_text(src)
    print(f"patched {path}")

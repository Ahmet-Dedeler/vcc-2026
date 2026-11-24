"""429 tests for the production zs pass (`VCC_DE_MAG`, src/vcc2026/rank_mag.py) and its builder patch.

Dry-copy mode (before the PI applies anything):
  TEST_MODULE=<patched module copy>  TEST_BUILDER=<patched builder copy>  uv run python test_zs429.py
The builder copy is copied to scripts/.dry429_<pid>.py for the parse test (the builder resolves ROOT from its own path)
and deleted afterwards. Defaults: the production files (src/vcc2026/rank_mag.py, scripts/make_submission_pm_cd4.py).

Sections:
  1. module invariants on real H1 blocks (429's dumped production-bulk blocks, base and MDE carriers, K = 1..16; 423's
     cap15 blocks, K = 1): k = 1 is a bit-identical no-op; dn_k250_zs keeps per-gene totals, libraries, the ballast rows,
     non-negativity, every gene's rank sum R1 and the calls (scorer arithmetic, rp423.r1_tie), no sign flips, untouched
     genes' tie term; something moves; a K = 0 block works;
  2. known answer: the module equals rp423.postpass(mode zs) bit-for-bit on the same blocks (block and shared
     diagnostics); block_lfc equals the evaluator's lfc bit-for-bit; apply() == postpass(block_lfc) and its diagnostics
     are flat floats (the builder's rc_diagnostics_summary calls float() on every value);
  3. the scorer-order control block rebuilt from the h1-test controls equals 323's reference.npz keep / cs / ref_lin
     (sha256), and matches the emitter gate (skip with --no-control; the whole test peaks at 7.7 GB RSS);
  4. builder: parse + refusals on the dry copy; the pinned module sha; run order (zero repair < zs < to_csr); and the
     patch is additions only: deleting its snippets restores the unpatched builder byte-for-byte (TEST_BUILDER_ORIG).
usage: uv run python test_zs429.py [--no-control]
"""
from __future__ import annotations

import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(ROOT / "experiments/323-current-family-validation"), str(ROOT / "src"),
           str(ROOT / "experiments/423-rank-preserving-magnitude"), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import evaluate323 as E  # noqa: E402  (read-only)
import rp423  # noqa: E402  (423's reference implementation, read-only)
import patch_zs429 as P  # noqa: E402

MOD = Path(os.environ.get("TEST_MODULE", ROOT / "src/vcc2026/rank_mag.py"))
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BUILDER_ORIG = os.environ.get("TEST_BUILDER_ORIG")
spec_ = importlib.util.spec_from_file_location("rank_mag_under_test", MOD)
RM = importlib.util.module_from_spec(spec_)
sys.modules[spec_.name] = RM   # the dataclass needs its module registered
spec_.loader.exec_module(RM)
FAILS: list[str] = []


def check(ok, what):
    print(f"[{'PASS' if ok else 'FAIL'}] {what}", flush=True)
    if not ok:
        FAILS.append(what)


def load(path):
    z = np.load(path)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    return X, int(z["n_car"]), str(z["target"]) if "target" in z.files else path.stem


def blocks():
    out = sorted((HERE / "work/blocks").glob("base_s0_bps_*.npz"))
    out += sorted((ROOT / "experiments/423-rank-preserving-magnitude/work/explore").glob("block_s0_0[16].npz"))
    return out


def section_module(cc):
    paths = blocks()
    check(len(paths) >= 8, f"{len(paths)} real H1 blocks to test (429 dumps + 423 cap15)")
    shared = ("n_touch", "n_moved", "counts_moved", "frac_mean", "frac_massw", "frac_massw_dn", "frac_full", "frac_zero",
              "sign_guard", "dU_car_abs_sum", "sub_used", "sub_budget_exhausted_cells")
    for path in paths:
        X, n_car, tgt = load(path)
        name = f"{path.parent.name}/{path.stem} ({tgt}, K={X.shape[0] - n_car})"
        s, j, p, q, l = E.block_stats(X, cc)
        R0, T0 = rp423.r1_tie(X, cc)
        # k = 1 no-op
        X1, _ = RM.postpass(X, n_car, cc, l, RM.parse_spec("dn_k100_zs"))
        check(np.array_equal(X1, X), f"{name}: dn_k100_zs is a bit-identical no-op")
        # the shipped rule
        X2, dg = RM.postpass(X, n_car, cc, l, RM.parse_spec("dn_k250_zs"))
        touched = dg["_touched"]
        s2, j2, p2, q2, l2 = E.block_stats(X2, cc)
        R1, T1 = rp423.r1_tie(X2, cc)
        fin = np.isfinite(l) & np.isfinite(l2)
        check(np.array_equal(X2.sum(0), X.sum(0)) and np.array_equal(s2, s), f"{name}: per-gene block totals")
        check(np.array_equal(X2.sum(1), X.sum(1)), f"{name}: libraries")
        check(np.array_equal(X2[n_car:], X[n_car:]), f"{name}: ballast rows untouched")
        check(X2.min() >= 0, f"{name}: non-negative")
        check(np.array_equal(R1, R0), f"{name}: rank sum R1 identical for every DE-universe gene")
        check(np.array_equal(q < .05, q2 < .05), f"{name}: calls identical")
        check(not np.any(np.sign(l[fin]) != np.sign(l2[fin])), f"{name}: no lfc sign flip")
        check(np.array_equal(T1[~touched], T0[~touched]), f"{name}: untouched genes' tie term exact")
        moved_dn = bool(touched.any()) and bool(np.all(l[touched] < 0)) and bool(np.all(l2[touched] <= l[touched]))
        check(moved_dn, f"{name}: {int(touched.sum())} genes moved, all down, none raised")
        # known answer vs 423's implementation
        Xr, dr = rp423.postpass(X, n_car, cc, l, rp423.parse_variant("dn_k250_zs"))
        check(np.array_equal(Xr, X2), f"{name}: module == rp423.postpass (zs) bit-for-bit")
        check(all(dg[k] == dr[k] for k in shared), f"{name}: shared diagnostics equal rp423's")
        # production entry point
        lb = RM.block_lfc(X, cc)
        check(np.array_equal(lb, l, equal_nan=True), f"{name}: block_lfc == evaluator lfc bit-for-bit")
        X3, d3 = RM.apply(X, n_car, cc, "dn_k250_zs")
        flat = all(isinstance(v, float) for v in d3.values()) and all(k.startswith("zs_") for k in d3)
        check(np.array_equal(X3, X2) and flat,
              f"{name}: apply() == postpass(block_lfc); {len(d3)} flat float diagnostics")
    # K = 0: the carriers alone (the pass never needs the ballast)
    X, n_car, tgt = load(paths[0])
    Xc = X[:n_car]
    X4, d4 = RM.apply(Xc, Xc.shape[0], cc, "dn_k250_zs")
    R0, _ = rp423.r1_tie(Xc, cc)
    R1, _ = rp423.r1_tie(X4, cc)
    check(np.array_equal(X4.sum(0), Xc.sum(0)) and np.array_equal(X4.sum(1), Xc.sum(1)) and np.array_equal(R1, R0)
          and d4["zs_n_moved"] > 0, f"K = 0 block ({tgt} carriers only): runs, totals/libraries/R1 exact, "
                                    f"{int(d4['zs_n_moved'])} genes moved")
    # parse grammar
    for good in ("dn_k250_zs", "all_k140_zs", "q80dn_k200_zs", "dn_k1000_zs"):
        try:
            RM.parse_spec(good)
            ok = True
        except ValueError:
            ok = False
        check(ok, f"parse_spec accepts {good}")
    for bad in ("dn_k250", "dn_k250_pc", "dn_k250_ex", "all_k140u120_zs", "dn_k050_zs", "dn_k2000_zs", "xx_k250_zs", ""):
        try:
            RM.parse_spec(bad)
            ok = False
        except ValueError:
            ok = True
        check(ok, f"parse_spec refuses {bad!r}")


def section_control(cc, emit_keep):
    import h5py  # noqa: PLC0415
    # stream the h1-test rows in chunks and keep only the controls (the full matrix is 2.7 GB)
    parts = []
    with h5py.File(ROOT / "data/harness/h1-test.h5ad", "r") as h:
        lab = E.labels(h, "obs/target_gene")
        ip = h["X/indptr"][:]
        shape = tuple(h["X"].attrs["shape"])
        for lo in range(0, shape[0], 4000):
            hi = min(lo + 4000, shape[0])
            a, b = int(ip[lo]), int(ip[hi])
            x = sp.csr_matrix((h["X/data"][a:b], h["X/indices"][a:b], ip[lo:hi + 1] - a), shape=(hi - lo, shape[1]))
            parts.append(x[lab[lo:hi] == "non-targeting"])
    xc = sp.vstack(parts, format="csr")
    del parts
    sc = RM.scorer_control_block(xc)
    del xc

    def h_(a):
        # hash the buffer in place (tobytes() would copy the 1.6 GB cs block)
        return hashlib.sha256(np.ascontiguousarray(a).reshape(-1).view(np.uint8)).hexdigest()
    for k in ("keep", "cs", "ref_lin"):
        a, b = getattr(sc, k), getattr(cc, k)
        check(a.dtype == b.dtype and a.shape == b.shape and h_(a) == h_(b),
              f"scorer_control_block(h1-test controls).{k} == 323 reference.npz {k} (sha256 {h_(a)[:16]})")
    check(np.array_equal(sc.keep, emit_keep), "scorer gate == emitter gate on the h1-test controls")
    check(sc.n_controls == 18400, "18,400 controls")


def run_builder(copy, env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    code = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(copy)!r}); m=u.module_from_spec(s); "
            "sys.modules['b']=m; s.loader.exec_module(m); print('DM', m.DE_MAG, m.DE_MAG_MODULE_SHA256)")
    return subprocess.run([sys.executable, "-c", code], env=e, capture_output=True, text=True, cwd=str(ROOT))


def section_builder():
    text = BUILDER.read_text()
    mod_sha = hashlib.sha256(MOD.read_bytes()).hexdigest()
    check(f'DE_MAG_MODULE_SHA256 = "{mod_sha}"' in text, f"builder pins the module under test (sha {mod_sha[:16]})")
    check(mod_sha == P.MODULE_SHA256, "module under test == the tested src_rank_mag.py (patch_zs429.MODULE_SHA256)")
    i_zr = text.find("_zr.zero_repair(")
    i_zs = text.find("_rank_mag.apply(")
    i_csr = text.find("values, columns, pointers = rc_emit.to_csr(block)")
    check(0 < i_zr < i_zs < i_csr, "run order in the block loop: ballast_rows -> zero repair -> zs pass -> to_csr")
    i_ctrl = text.find("_rank_mag.scorer_control_block(control_matrix")
    i_del = text.find("                del control_matrix\n")
    check(0 < i_ctrl < i_del, "the scorer-order control block is built before the count matrix is released")
    if BUILDER_ORIG:
        orig = Path(BUILDER_ORIG).read_text()
        undo = text
        for snip in (P.KNOB.replace("__MODULE_SHA256__", P.MODULE_SHA256), P.CTRL, P.PASS, P.CFG,
                     'if DE_MAG is not None and PERT_SCALE_BALLAST is None:\n'
                     '    raise SystemExit("VCC_DE_MAG needs VCC_PERT_SCALE_BALLAST (it runs after ballast_rows and the zero repair)")\n',
                     "            rc_scorer_control = None   # 429: VCC_DE_MAG's evaluator-order control block\n",
                     "            rc_scorer_control = None\n"):
            undo = undo.replace(snip, "", 1)
        check(undo == orig, f"additions only: removing the patch's snippets restores {Path(BUILDER_ORIG).name} byte-for-byte")
    dry = ROOT / f"scripts/.dry429_{os.getpid()}.py"
    shutil.copyfile(BUILDER, dry)
    try:
        base = dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
        r = run_builder(dry, dict(base))
        check(r.returncode == 0 and "DM None" in r.stdout, f"unset: DE_MAG None {r.stderr[-300:]}")
        r = run_builder(dry, dict(base, VCC_DE_ZERO_REPAIR="-5", VCC_DE_MDE="0.1,0.1,0.3", VCC_DE_MAG="dn_k250_zs"))
        check(r.returncode == 0 and f"DM dn_k250_zs {mod_sha}" in r.stdout,
              f"VCC_DE_MAG=dn_k250_zs parses (with MDE and zero repair set) {r.stderr[-300:]}")
        for env, want in (
                (dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_DE_MAG="dn_k250_zs"),
                 "needs VCC_PERT_SCALE_BALLAST"),
                (dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1", VCC_DE_MAG="dn_k250_zs"),
                 "needs VCC_DE_SIGN_PRIOR"),
                (dict(base, VCC_DE_MAG="dn_k250_pc"), "must look like"),
                (dict(base, VCC_DE_MAG="dn_k250"), "must look like"),
                (dict(base, VCC_DE_MAG="dn_k050_zs"), "k must be in"),
                (dict(base, VCC_DE_MAG="x"), "must look like")):
            r = run_builder(dry, env)
            out = r.stderr + r.stdout
            check(r.returncode != 0 and want in out, f"refuses {env.get('VCC_DE_MAG')!r} "
                                                     f"({'no ballast' if 'VCC_PERT_SCALE_BALLAST' not in env else 'no sign prior' if 'VCC_DE_SIGN_PRIOR' not in env else 'grammar'}): {want!r}")
    finally:
        dry.unlink(missing_ok=True)
    check(not dry.exists(), f"dry copy {dry.name} deleted")


def main():
    print(f"module under test {MOD} sha {hashlib.sha256(MOD.read_bytes()).hexdigest()}")
    print(f"builder under test {BUILDER} sha {hashlib.sha256(BUILDER.read_bytes()).hexdigest()}")
    r, cc, emit_c = E.load_reference()
    emit_keep = np.asarray(emit_c.keep, dtype=bool).copy()
    del r, emit_c                  # only the scorer-order block is needed from here (the emitter's cs is 1.6 GB)
    section_module(cc)
    if "--no-control" not in sys.argv:
        section_control(cc, emit_keep)
    section_builder()
    print(f"{len(FAILS)} failed" + (": " + "; ".join(FAILS) if FAILS else ""))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()

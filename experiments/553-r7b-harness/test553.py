#!/usr/bin/env python3
"""553 — fixture tests (light, seconds): the group-by-group functions of lib553 against the full-AnnData functions they
replace, on synthetic data written to work/test553/.

  T1 read_rows(path, rows) == ad.read_h5ad(path).X[rows]  (and rows=None == the whole X)
  T2 anchors(...) == anchors323: full stats, the dispersed seed-0 baseline scored with fast_score, the 5-split replicate
     (A.replicate on the in-memory AnnData), bit for bit, on (a) 323's own 3-target fixture and (b) a larger fixture with
     unequal group sizes (odd counts), 150 controls and a categorical label column
  T3 controlref(...) on the fixture's controls == evaluate323.reference()'s control half computed inline
  T4 the sealing guard: open_truth refuses r7b's truth when the access is not armed
usage: uv run python experiments/553-r7b-harness/test553.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib553 as S  # noqa: E402


def fixture_big(seed=553):
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    G = 60
    arrays = [rng.poisson(6, (151, G)).astype(np.float32)]
    labels = [S.CTRL] * 151
    sizes = [37, 52, 41, 60, 45]
    for i, n in enumerate(sizes):
        x = rng.poisson(6, (n, G)).astype(np.float32)
        blk = x[:, 2 + 9 * i: 2 + 9 * (i + 1)]
        x[:, 2 + 9 * i: 2 + 9 * (i + 1)] = blk * 3 if i % 2 else rng.binomial(blk.astype(np.int64), 0.2)   # integer counts
        arrays.append(x)
        labels += [f"G{i}"] * n
    X = np.vstack(arrays)
    perm = rng.permutation(X.shape[0])            # interleave the groups, as h1-test's rows are
    X, labels = X[perm], [labels[i] for i in perm]
    return ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"target_gene": pd.Categorical(labels)},
                                                            index=[f"c{i}" for i in range(len(labels))]),
                      var=pd.DataFrame(index=[f"G{i}" for i in range(G)]))


def check_anchors(a, name, path):
    import anchors323 as A
    from cell_eval2.baseline import build_baseline_prediction, generic_response_profile
    a.write_h5ad(path)
    spec = dict(name=f"fixture-{name}", sealed=False, truth_h5ad=path)
    T = S.truth_rows(spec)
    rec, full = S.anchors(spec, T)
    # the access's path: anchors fed with truth_reference's control block and control moments
    T2 = S.truth_rows(spec)
    cm = S.control_moments(T2)
    ref, c, de = S.truth_reference(spec, T2, cm)
    rec2, full2 = S.anchors(spec, T2, ctrl_block=c, ctrl_moments=cm)
    import anndata as ad
    import evaluate323 as E
    import exact_de323
    from cell_eval2.moments import jackknife_correction
    b = ad.read_h5ad(path)
    # T5: truth_reference == evaluate323.reference()'s arithmetic, inline on the AnnData
    labs = b.obs.target_gene.to_numpy().astype(str)
    xc = b.X[labs == S.CTRL]
    cc = exact_de323.control_block(xc)
    tg = sorted(set(labs) - {S.CTRL})
    st = [E.block_stats(b.X[labs == t], cc) for t in tg]
    cp = np.asarray(xc.astype(np.float64).sum(0)).ravel().astype(np.int64)
    cj = jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, S.TS, chunk=32)[0]
    t5 = (np.array_equal(ref["sums"], np.vstack([*[x[0] for x in st], cp])) and np.array_equal(ref["jk"], np.array([*[x[1] for x in st], cj]))
          and all(np.array_equal(ref[k], np.asarray([x[i] for x in st])) for k, i in (("p", 2), ("q", 3), ("lfc", 4)))
          and np.array_equal(c.cs, cc.cs) and de.equals(E.de_frame(np.array(tg), np.asarray(b.var_names)[cc.keep],
                                                                   [x[2] for x in st], [x[3] for x in st], [x[4] for x in st])))
    print(f"T5 truth_reference [{name}]: {'PASS' if t5 else 'FAIL'}", flush=True)
    same_paths = (json.dumps(rec2["baseline"]) == json.dumps(rec["baseline"]) or
                  max(abs(rec2["baseline"][m] - rec["baseline"][m]) for m in S.MEMBERS) <= 1e-12) and \
        all(np.array_equal(full2[k], full[k]) for k in ("perts", "sums", "counts", "jk"))
    print(f"T2b anchors via the access's path (truth_reference's block + control moments) == standalone: "
          f"{'PASS' if same_paths else 'FAIL'}", flush=True)
    ref_full = A.stats(b)
    p = build_baseline_prediction(generic_response_profile(b, pert_col="target_gene", control=S.CTRL), b,
                                  pert_col="target_gene", control=S.CTRL, emit="dispersed", seed=0)
    ref_base = A.fast_score(A.stats(p), ref_full)
    splits, ref_rep = A.replicate(b, ref_full)
    d = {f"base.{m}": abs(rec["baseline"][m] - ref_base[m]) for m in S.MEMBERS}
    d.update({f"rep.{m}": abs(rec["replicate"][m] - ref_rep[m]) for m in S.MEMBERS})
    for s1, s2 in zip(rec["replicate_splits"], splits):
        d.update({f"split{s1['split']}.{m}": abs(s1["raw"][m] - s2["raw"][m]) for m in S.MEMBERS})
    same_full = all(np.array_equal(full[k], ref_full[k]) for k in ("perts", "sums", "counts", "jk")) and full["de"].equals(ref_full["de"])
    # the access's per-half path (each half read from the file on its own) == the in-process replicate
    for sp_ in rec["replicate_splits"]:
        ga, gb = S.replicate_halves(spec, sp_["split"])
        raw_h, _ = S.replicate_score(S.replicate_half_stats(spec, ga), S.replicate_half_stats(spec, gb), full["de"])
        d.update({f"perhalf{sp_['split']}.{m}": abs(raw_h[m] - sp_["raw"][m]) for m in S.MEMBERS})
    mx = max(d.values())
    ok = mx <= 1e-12 and same_full and t5 and same_paths   # cell-eval2's nmae varies at the last ulp between calls
    print(f"T2 anchors [{name}]: {'PASS' if ok else 'FAIL'} max |d| {mx:.3e}, full stats equal {same_full}; "
          f"baseline {rec['baseline']}", flush=True)
    return ok


def main():
    import anndata as ad
    import anchors323 as A
    out = HERE / "work/test553"
    out.mkdir(parents=True, exist_ok=True)
    ok = True
    # T1
    a = fixture_big()
    p = out / "fixture_big.h5ad"
    a.write_h5ad(p)
    b = ad.read_h5ad(p)
    rows = np.array([0, 3, 4, 17, 100, 250, 300, b.n_obs - 1])
    x1 = S.read_rows(p, rows, chunk_rows=64)
    t1 = (x1 != b.X[rows]).nnz == 0 and x1.dtype == b.X.dtype
    t1b = (S.read_rows(p, None) != b.X).nnz == 0
    print(f"T1 read_rows: {'PASS' if t1 and t1b else 'FAIL'}", flush=True)
    ok &= t1 and t1b
    # T2
    ok &= check_anchors(A.fixture(), "323-fixture", out / "fixture_323.h5ad")
    ok &= check_anchors(fixture_big(), "big", out / "fixture_big2.h5ad")
    # T3
    import exact_de323
    from cell_eval2.moments import jackknife_correction
    from vcc2026 import rc_emit
    cr = S.controlref(p)
    lab = b.obs.target_gene.to_numpy().astype(str)
    xc = b.X[lab == S.CTRL]
    ec = rc_emit.control_block_from_counts(xc)
    c = exact_de323.control_block(xc)
    cp = np.asarray(xc.astype(np.float64).sum(0)).ravel().astype(np.int64)
    cj = jackknife_correction(xc, np.zeros(xc.shape[0], dtype=np.int64), 1, S.TS, chunk=32)[0]
    t3 = all([np.array_equal(cr["keep"], c.keep), np.array_equal(cr["cs"], c.cs), np.array_equal(cr["tie"], c.tie),
              np.array_equal(cr["ref_lin"], c.ref_lin), np.array_equal(cr["emit_cs"], ec.cs),
              np.array_equal(cr["emit_tie"], ec.tie), np.array_equal(cr["emit_ref_lin"], ec.ref_lin),
              np.array_equal(cr["ctrl_sums"], cp), cr["ctrl_jk"] == cj,
              np.array_equal(cr["libraries"], np.asarray(xc.sum(1)).ravel().astype(float)),
              np.array_equal(cr["b_ctrl"], np.log1p(S.TS * cp / cp.sum()))])
    print(f"T3 controlref: {'PASS' if t3 else 'FAIL'}", flush=True)
    ok &= t3
    # T6: target exclusion == the target's cells absent from the file (reference, anchors incl. the replicate's RNG stream)
    out6 = HERE / "work/test553"
    a6 = fixture_big(seed=5536)
    a6.write_h5ad(out6 / "fx6_full.h5ad")
    keep6 = a6.obs.target_gene.to_numpy().astype(str) != "G2"
    a6[keep6].copy().write_h5ad(out6 / "fx6_noG2.h5ad")
    sp_full = dict(name="fx6-full", sealed=False, truth_h5ad=out6 / "fx6_full.h5ad")
    sp_no = dict(name="fx6-noG2", sealed=False, truth_h5ad=out6 / "fx6_noG2.h5ad")
    Ta, Tb = S.truth_rows(sp_full, exclude=("G2",)), S.truth_rows(sp_no)
    ra, ca, dea = S.truth_reference(sp_full, Ta, S.control_moments(Ta))
    rb, cb, deb = S.truth_reference(sp_no, Tb, S.control_moments(Tb))
    same_ref = all(np.array_equal(ra[k], rb[k]) for k in ("targets", "sums", "counts", "jk", "p", "q", "lfc")) and dea.equals(deb)
    rec_a, _ = S.anchors(sp_full, Ta)
    rec_b, _ = S.anchors(sp_no, Tb)
    dd = max(abs(rec_a[k][m] - rec_b[k][m]) for k in ("baseline", "replicate") for m in S.MEMBERS)
    for s1, s2 in zip(rec_a["replicate_splits"], rec_b["replicate_splits"]):
        dd = max(dd, max(abs(s1["raw"][m] - s2["raw"][m]) for m in S.MEMBERS))
    t6 = same_ref and dd <= 1e-12
    print(f"T6 exclusion == absent: {'PASS' if t6 else 'FAIL'} (reference equal {same_ref}, anchors max |d| {dd:.1e})", flush=True)
    ok &= t6
    # T4
    try:
        S.open_truth(S.surface_spec("r7b"))
        t4 = False
    except PermissionError:
        t4 = True
    print(f"T4 sealing guard: {'PASS' if t4 else 'FAIL'}", flush=True)
    ok &= t4
    # T8: control_block_keep with the block's own keep == exact_de323.control_block (the csplit phase's control halves)
    import exact_de323
    from scipy import sparse as _sp
    rng8 = np.random.default_rng(5538)
    x8 = _sp.csr_matrix(rng8.poisson(rng8.gamma(0.3, 3.0, size=(1, 400)), size=(150, 400)).astype(np.float32))
    c0 = exact_de323.control_block(x8)
    c1 = S.control_block_keep(x8, c0.keep)
    t8 = (np.array_equal(c0.keep, c1.keep) and np.array_equal(c0.cs, c1.cs) and np.array_equal(c0.tie, c1.tie)
          and np.array_equal(c0.ref_lin, c1.ref_lin) and c0.n_controls == c1.n_controls)
    print(f"T8 control_block_keep == control_block: {'PASS' if t8 else 'FAIL'} ({int(c0.keep.sum())} gated genes)", flush=True)
    ok &= t8
    # T7: a floored prediction side (r7bd<floor>) without 552's floored truth stops in the access's orchestrator before any
    # phase, and truth_rows refuses it too. A stand-in spec (unsealed, nonexistent truth path): no R7b truth is touched.
    import access553 as A
    real_spec, real_run = S.surface_spec, A.subprocess.run
    phases = []

    def stand_in(name):
        sp = real_spec(name)
        if name.startswith("r7bd"):
            sp = dict(sp, sealed=False, truth_h5ad=Path("/nonexistent/truth.h5ad"))
            sp.pop("floored_build_json", None)
        return sp
    S.surface_spec = stand_in
    A.subprocess.run = lambda *a_, **k_: phases.append(a_)
    argv0 = sys.argv
    sys.argv = ["access553.py", "r7bd5000", "--arms", "lcap-kept", "--seeds", "0", "--exclude-targets", "BST2",
                "--out", str(HERE / "work/test553/t7")]
    try:
        A.main()
        t7a = False
    except SystemExit as e:
        t7a = "has not built the floored truth" in str(e) and not phases
    finally:
        S.surface_spec, A.subprocess.run, sys.argv = real_spec, real_run, argv0
    try:
        S.truth_rows(stand_in("r7bd5000"))
        t7b = False
    except PermissionError:
        t7b = True
    t7 = t7a and t7b
    print(f"T7 floored prediction side without a floored truth refused: {'PASS' if t7 else 'FAIL'} (orchestrator {t7a}, "
          f"truth_rows {t7b})", flush=True)
    ok &= t7
    print("ALL PASS" if ok else "SOME FAILED", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

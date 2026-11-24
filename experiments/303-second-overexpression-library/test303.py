"""303 — synthetic and refusal tests.  Reads no reserve value.

Synthetic:  S1 zero / S2 exact mirror / S3 exact alignment / S4 partial mirror
            recovered at the declared size; S5 the shuffled-pairing null is
            centred where the construction puts it, not assumed at 0.
Refusal:    R1 wrong axis, R2 non-binary mask, R3 nonfinite input, R4 CPM sum
            drift, R5 a value read of a reserve dataset, R6 no release,
            R7 tampered release, R8 second access after a marker.
Hygiene:    H1 determinism under a fixed seed, H2 finite/dtype validation,
            H3 the metadata/preflight path opens no R5 slice of X or control.

    PYTHONDONTWRITEBYTECODE=1 uv run python test303.py
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib303 as L  # noqa: E402
import stats303 as S  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []
SEED = 20260907


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name:38s} {detail}")


# ------------------------------------------------------------------ synthetic
def synth(n_tf=100, n_gene=400, seed=SEED):
    rng = np.random.default_rng(seed)
    KD = rng.normal(size=(n_tf, n_gene))
    own = rng.permutation(n_gene)[:n_tf].astype(np.int64)
    return rng, KD, own


def run_synthetic() -> None:
    rng, KD, own = synth()
    n = KD.shape[0]

    # S1 — independent OE: median S ~ 0, no tail asymmetry
    OE = rng.normal(size=KD.shape)
    s1 = S.cell_stats(S.s_vector(OE, KD, own))
    check("S1 zero (independent draws)", abs(s1["median"]) < 0.03 and abs(s1["tail_asym"]) < 0.10,
          f"median {s1['median']:+.4f} tail_asym {s1['tail_asym']:+.3f}")

    # S2 — exact mirror: every S = -1
    s2 = S.cell_stats(S.s_vector(-KD, KD, own))
    check("S2 exact mirror -> S = -1", np.isclose(s2["median"], -1.0, atol=1e-9)
          and s2["share_below"] == 1.0, f"median {s2['median']:+.6f}")

    # S3 — exact alignment: every S = +1
    s3 = S.cell_stats(S.s_vector(KD.copy(), KD, own))
    check("S3 exact alignment -> S = +1", np.isclose(s3["median"], 1.0, atol=1e-9)
          and s3["share_above"] == 1.0, f"median {s3['median']:+.6f}")

    # S4 — planted partial mirror of size rho recovered within tolerance
    rho = 0.30
    noise = rng.normal(size=KD.shape)
    KDu = KD / np.linalg.norm(KD, axis=1, keepdims=True)
    nu = noise - (noise * KDu).sum(1, keepdims=True) * KDu
    nu /= np.linalg.norm(nu, axis=1, keepdims=True)
    OE4 = -rho * KDu + np.sqrt(1 - rho ** 2) * nu
    s4 = S.cell_stats(S.s_vector(OE4, KD, own))
    check("S4 planted mirror rho=-0.30", abs(s4["median"] + rho) < 0.02,
          f"median {s4['median']:+.4f} (target {-rho:+.2f})")

    # S5 — the null is measured, not assumed to be 0
    nd = S.null_distribution(OE4, KD, own, draws=60, seed=SEED)
    check("S5 shuffled null measured", np.isfinite(nd["median"]).all()
          and abs(float(np.median(nd["median"]))) < 0.05,
          f"null median {np.median(nd['median']):+.5f} sd {np.std(nd['median']):.5f}; "
          f"observed {s4['median']:+.4f} -> p="
          f"{L.perm_p(s4['median'], nd['median'], 'less'):.4g}")

    # S6 — the own-column removal is symmetric between matched and shuffled pairs
    perm = S.derangements(n, 1, np.random.default_rng(SEED))[0]
    a = S.s_vector(OE4, KD, own, None)
    b = S.s_vector(OE4, KD, own, perm)
    check("S6 own-column removal both sides", np.isfinite(a).all() and np.isfinite(b).all(),
          f"matched n={np.isfinite(a).sum()} shuffled n={np.isfinite(b).sum()}")

    # H1 — determinism
    n1 = S.null_distribution(OE4, KD, own, draws=20, seed=7)["median"]
    n2 = S.null_distribution(OE4, KD, own, draws=20, seed=7)["median"]
    check("H1 determinism (fixed seed)", np.array_equal(n1, n2), "identical null draws")

    # H2 — Holm control over four cells
    h = L.holm(dict(C7a=0.004, C7b=0.20, C7c=0.03, C7d=0.0001))
    check("H2 Holm familywise 0.05", h["C7d"]["reject"] and not h["C7b"]["reject"]
          and h["C7a"]["p_holm"] >= h["C7d"]["p_holm"],
          " ".join(f"{k}:{v['p_holm']:.4f}{'*' if v['reject'] else ''}" for k, v in
                   sorted(h.items())))


# -------------------------------------------------------------------- refusal
def run_refusal() -> None:
    m = L.read_db_meta()
    proj = L.build_projection(m.var_symbol, L.canonical_genes())
    lab = ["row0"]

    # R1 wrong axis / shape mismatch
    try:
        L.bank_transform(np.zeros((2, 5)), np.zeros((2, 6)), np.zeros((2, 5)), proj, lab)
        check("R1 wrong-axis refusal", False, "no exception")
    except L.Refuse as e:
        check("R1 wrong-axis refusal", True, str(e)[:60])

    # R2 non-binary mask
    x = np.zeros((1, m.n_var)); c = np.zeros((1, m.n_var))
    mm = np.full((1, m.n_var), 0.5)
    try:
        L.bank_transform(x, c, mm, proj, lab)
        check("R2 non-binary mask refusal", False, "no exception")
    except L.Refuse as e:
        check("R2 non-binary mask refusal", True, str(e)[:60])

    # R3 nonfinite input
    x2 = np.zeros((1, m.n_var)); x2[0, 0] = np.nan
    try:
        L.bank_transform(x2, c, np.ones((1, m.n_var)), proj, lab)
        check("R3 nonfinite refusal", False, "no exception")
    except L.Refuse as e:
        check("R3 nonfinite refusal", True, str(e)[:60])

    # R4 CPM10k sum drift (all-zero control sums to 0, not 10_000)
    try:
        L.bank_transform(np.zeros((1, m.n_var)), np.zeros((1, m.n_var)),
                         np.ones((1, m.n_var)), proj, lab, strict_sums=True)
        check("R4 CPM sum-drift refusal", False, "no exception")
    except L.Refuse as e:
        check("R4 CPM sum-drift refusal", True, str(e)[:60])

    # R5 a value read of a reserve dataset through the safe reader
    r5 = np.flatnonzero(m.dataset == L.RESERVE_DATASET).astype(np.int64)[:3]
    try:
        L.read_source_rows(L.DB, r5, L.FREE_DATASETS, m.dataset)
        check("R5 reserve value-read refusal", False, "READ THE RESERVE")
    except L.Refuse as e:
        check("R5 reserve value-read refusal", True, str(e)[:70])
    # ... but the mask alone is free
    try:
        got = L.read_source_rows(L.DB, r5, L.FREE_DATASETS, m.dataset,
                                 what=("layers/measure_mask",))
        check("R5b mask-only read allowed", got["layers/measure_mask"].shape == (3, m.n_var),
              "measure_mask is free per the R5 declaration")
    except L.Refuse as e:
        check("R5b mask-only read allowed", False, str(e)[:70])


def run_gate() -> None:
    """Current strict-schema guard tests replace the incompatible draft fixture.

    The original test303.json is a preflight record, not proof about revised
    code. Every refusal and successful synthetic access now runs in temporary
    fixtures owned by test_access_guard303; no real release is fabricated.
    """
    import io
    import unittest
    import test_access_guard303 as G
    suite = unittest.TestSuite([
        unittest.defaultTestLoader.loadTestsFromTestCase(G.AccessGuardTests),
        unittest.defaultTestLoader.loadTestsFromTestCase(G.RealInputSchemaTests),
    ])
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=1).run(suite)
    check("R6-R8 current strict-schema suite", result.wasSuccessful(),
          f"{result.testsRun} fixture/schema tests; "
          f"{len(result.failures)} failures, {len(result.errors)} errors")
    if not result.wasSuccessful():
        print(stream.getvalue())


# ----------------------------------------------- H3: preflight touches no R5 X
def run_no_reserve_read() -> None:
    m = L.read_db_meta()
    r5 = np.flatnonzero(m.dataset == L.RESERVE_DATASET)
    lo, hi = int(r5.min()), int(r5.max())
    touched: list[tuple[str, int, int]] = []

    orig_getitem = h5py.Dataset.__getitem__

    def spy(self, key):
        name = self.name
        if name in ("/X", "/layers/control"):
            idx = key[0] if isinstance(key, tuple) else key
            rows = np.atleast_1d(np.asarray(idx))
            if rows.dtype != bool and rows.size:
                inside = rows[(rows >= lo) & (rows <= hi)]
                if inside.size:
                    touched.append((name, int(inside.min()), int(inside.max())))
        return orig_getitem(self, key)

    h5py.Dataset.__getitem__ = spy
    try:
        import importlib
        import meta303
        importlib.reload(meta303)
        meta303.main()
    finally:
        h5py.Dataset.__getitem__ = orig_getitem
    check("H3 preflight reads no R5 X/control", not touched,
          f"{len(touched)} reserve value slices" if touched else "0 reserve value slices")


def main() -> int:
    t0 = time.time()
    run_synthetic()
    run_refusal()
    run_gate()
    run_no_reserve_read()
    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    L.dump_json(L.WORK / "test303.json",
                dict(results=[dict(name=n, pass_=o, detail=d) for n, o, d in RESULTS],
                     n=len(RESULTS), n_fail=n_fail, seconds=round(time.time() - t0, 1)))
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} tests pass "
          f"({time.time() - t0:.1f}s)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

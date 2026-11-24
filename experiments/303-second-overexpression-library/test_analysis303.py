"""Synthetic-only readout tests; no R5/source access, no real output writes.

    PYTHONDONTWRITEBYTECODE=1 uv run python test_analysis303.py
    ... test_analysis303.py --manifest /path/to/frozen-core-manifest.json

The second form also verifies the readout's USED inputs against the coordinator
manifest and emits a frozen_analysis303 receipt to stdout. The coordinator,
not this suite, saves that receipt. It never creates a release or access marker.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import spearmanr

import analyse303 as A
import stats303 as S


def toy(seed=S.SEED, n=40, g=160):
    rng = np.random.default_rng(seed)
    return rng, rng.normal(size=(n, g)), np.arange(n, dtype=np.int64)


class StatisticsTests(unittest.TestCase):
    def test_exact_mirror_and_alignment(self):
        rng, kd, own = toy()
        x = np.arange(len(kd), dtype=float)
        for sign in (-1, 1):
            oe, kk = S.cohort_residuals(sign * kd, kd)
            C = S.cosine_matrix(oe, kk, own)
            np.testing.assert_allclose(np.diag(C), sign, atol=1e-14)
            r = S.evaluate_c7(C, x, x, draws=2000, boot_draws=100)
            self.assertEqual(r["cells"]["C7a"]["verdict"], "PASS" if sign == -1 else "FAIL")
            if sign == -1:
                self.assertEqual(r["cells"]["C7a"]["p_raw"], 1/2001)

    def test_independent_zero(self):
        rng, kd, own = toy(n=100, g=600)
        oe, kk = S.cohort_residuals(rng.normal(size=kd.shape), kd)
        r = S.evaluate_c7(S.cosine_matrix(oe, kk, own), np.arange(100.),
                          np.arange(100.), draws=2000, boot_draws=100)
        self.assertLess(abs(r["observed_summary"]["median"]), .03)
        self.assertNotEqual(r["cells"]["C7a"]["verdict"], "PASS")

    def test_zero_norm_is_unestimable_not_zero_signal(self):
        _, kd, own = toy()
        oe, kk = S.cohort_residuals(np.zeros_like(kd), kd)
        r = S.evaluate_c7(S.cosine_matrix(oe, kk, own), np.arange(40.), np.arange(40.),
                          draws=200, boot_draws=100)
        self.assertTrue(np.isnan(r["S"]).all())
        for cell in r["cells"].values():
            self.assertEqual(cell["verdict"], "UNESTIMABLE")
            self.assertEqual(cell["p_raw"], 1)

    def test_positive_but_below_null_cannot_pass_a(self):
        # Both one-sided null tests reject, but observed median is positive.
        obs = np.array([-.2]*40 + [.01]*60)
        null = np.full((2000, 100), .8)
        r = S.evaluate_from_null(obs, null, np.arange(100.), np.arange(100.), boot_draws=100)
        a = r["cells"]["C7a"]
        self.assertLess(a["p_holm"], .05)
        self.assertGreater(a["observed"]["tail_asym"], 0)
        self.assertFalse(a["direction_ok"])
        self.assertEqual(a["verdict"], "FAIL")

    def test_negative_rho_cannot_pass_relative_null(self):
        obs = np.linspace(-.7, .7, 40)
        obs[19:22] = obs[19:22][::-1]
        null = np.tile(np.linspace(-.7, .7, 40), (2000, 1))
        r = S.evaluate_from_null(obs, null, np.arange(40.), np.arange(40.), boot_draws=100)
        for key in ("C7b", "C7c"):
            self.assertLess(r["cells"][key]["observed"], 0)
            self.assertFalse(r["cells"][key]["direction_ok"])
            self.assertNotEqual(r["cells"][key]["verdict"], "PASS")

    def test_generic_kd_alignability_is_preserved_in_null(self):
        # KD geometry predicts generic overlap with ANY OE row. Permuting the
        # feature instead of the OE matching would incorrectly credit this.
        feature = np.linspace(.1, .8, 100)
        C = np.broadcast_to(feature[:, None], (100, 100)).copy()
        r = S.evaluate_c7(C, feature, feature, draws=2000, boot_draws=100)
        d = r["cells"]["C7d"]
        self.assertAlmostEqual(d["observed"], 1)
        self.assertAlmostEqual(d["null"]["mean"], 1)
        self.assertEqual(d["p_raw"], 1)
        self.assertEqual(d["verdict"], "FAIL")
        np.testing.assert_allclose(r["null_S"], np.broadcast_to(feature, (2000, 100)))

    def test_missing_features_keep_four_cell_family(self):
        _, kd, own = toy()
        oe, kk = S.cohort_residuals(-kd, kd)
        r = S.evaluate_c7(S.cosine_matrix(oe, kk, own), np.full(40, np.nan),
                          np.ones(40), draws=2000, boot_draws=100)
        self.assertEqual(tuple(r["cells"]), S.CELLS)
        for k in S.CELLS[1:]:
            self.assertEqual(r["cells"][k]["verdict"], "UNESTIMABLE")
            self.assertEqual(r["cells"][k]["p_raw"], 1)
        self.assertEqual(r["cells"]["C7a"]["p_holm"], 4/2001)

    def test_uniform_permutations_include_fixed_points_and_identity(self):
        p = S.uniform_permutations(4, 2000, np.random.default_rng(S.SEED))
        self.assertEqual(len({tuple(row) for row in p}), 24)
        self.assertTrue(np.any(np.all(p == np.arange(4), axis=1)))
        np.testing.assert_array_equal(p, S.uniform_permutations(4, 2000, np.random.default_rng(S.SEED)))
        with self.assertRaises(ValueError):
            S.validate_permutations([[0, 0, 1, 2]], 4)
        with self.assertRaises(ValueError):
            S.derangements(1, 1, np.random.default_rng(1))

    def test_vectorized_vs_literal_all_pairings(self):
        rng, kd, own = toy(n=9, g=45)
        own[-1] = -1
        oe = rng.normal(size=kd.shape)
        for exclusion in ("cohort", "pairwise", "none"):
            C = S.cosine_matrix(oe, kd, own, exclusion=exclusion)
            for i in range(9):
                for j in range(9):
                    drop = (own[own >= 0] if exclusion == "cohort" else
                            np.unique([x for x in (own[i], own[j]) if x >= 0]).astype(int)
                            if exclusion == "pairwise" else np.array([], int))
                    self.assertAlmostEqual(C[i, j], S.pair_cosine(kd[i], oe[j], drop), places=12)

    def test_all_cohort_targets_excluded_even_under_shuffle(self):
        rng, kd, own = toy(n=9, g=45)
        oe = rng.normal(size=kd.shape)
        before = S.cosine_matrix(oe, kd, own)
        oe[:, own] = rng.normal(size=(9, 9))*1e6
        kd[:, own] = rng.normal(size=(9, 9))*1e6
        np.testing.assert_array_equal(before, S.cosine_matrix(oe, kd, own))

    def test_nan_strict_before_exclusion_and_typed_indices(self):
        rng, kd, own = toy()
        for value in (np.nan, np.inf):
            bad = kd.copy()
            bad[0, own[0]] = value
            with self.assertRaises(ValueError):
                S.cosine_matrix(bad, kd, own)
        for bad in (own.astype(float), own.astype(bool), own-1000):
            with self.assertRaises(ValueError):
                S.cosine_matrix(kd, kd, bad)
        with self.assertRaises(ValueError):
            S.residualise(kd, np.full(kd.shape[1], np.nan))
        self.assertEqual(S.permutation_p(-1, np.array([np.nan, 0.]), "less"), 1)

    def test_ranks_and_bootstrap_rerank_ties(self):
        x, y = np.array([1., 1, 3, 5, 7, 9]), np.array([5., 4, 2, 2, 8, 9])
        b = S._bootstrap(y, x, draws=100, seed=9)["rho"]
        rng = np.random.default_rng(9)
        ids = rng.integers(0, len(x), (100, len(x)))
        literal = np.array([spearmanr(x[i], y[i]).statistic for i in ids])
        np.testing.assert_allclose(b, literal, atol=1e-14, equal_nan=True)

    def test_holm_conjunction_and_complete_family(self):
        h = S.holm_four(dict(C7a=.01, C7b=.03, C7c=.04, C7d=1))
        self.assertAlmostEqual(h["C7a"]["p_holm"], .04)
        self.assertAlmostEqual(h["C7b"]["p_holm"], .09)
        self.assertAlmostEqual(h["C7c"]["p_holm"], .09)
        with self.assertRaises(ValueError):
            S.holm_four(dict(C7a=.01))

    def test_full_size_full_budget(self):
        rng, kd, own = toy(n=100, g=15309)
        own[-2:] = -1
        oe, kk = S.cohort_residuals(-kd, kd)
        self.assertEqual(int(S.fixed_axis_mask(own, 100, 15309).sum()), 15211)
        cons = np.linspace(0, 1, 100)
        cons[:3] = np.nan
        r = S.evaluate_c7(S.cosine_matrix(oe, kk, own), cons, rng.normal(size=100))
        self.assertEqual(r["null_S"].shape, (2000, 100))
        self.assertEqual(r["bootstrap"]["C7a_median"].shape, (20000,))
        self.assertEqual(r["cells"]["C7b"]["n"], 97)
        self.assertEqual(r["cells"]["C7a"]["verdict"], "PASS")
        self.assertIn("not between-study", r["ci_scope"])


def fixture():
    rng = np.random.default_rng(8)
    n, g = 8, 30
    names = tuple(f"TF{i}" for i in range(n))
    axis = names + tuple(f"G{i}" for i in range(g-n))
    d = A.Design(axis, names, names, tuple(f"sample_{x}" for x in names),
                 np.arange(n), np.arange(n), np.arange(n), np.arange(n),
                 np.ones(n, bool))
    z = dict(rows=d.r5_rows, symbols=np.asarray(names), obs_id=np.asarray(d.r5_ids),
             axis=np.asarray(axis), delta=rng.normal(size=(n, g)).astype(np.float32),
             mask=np.ones((n, g), bool), b_ctrl=np.ones(g, np.float32),
             control_cpm_sum=np.full(n, 10000.), perturbed_cpm_sum=np.full(n, 10000.),
             control_rows_identical=np.array([True]))
    acc = dict(cpm_gate_pass=True, control_rows_identical=True, mask_row_constant=True,
               all_finite=True, delta_sha256=A.sha_array(z["delta"]), mask_sha256=A.sha_array(z["mask"]),
               n_rows=n, n_axis=g)
    return d, z, acc


class DriverTests(unittest.TestCase):
    def test_kd_row_and_gene_alignment_before_values(self):
        d, _, _ = fixture()
        perts = d.cohort + tuple(f"P{i}" for i in range(2354-len(d.cohort)))
        for genes, labels in ((d.axis, perts[::-1]), (d.axis[::-1], perts)):
            with patch.object(A, "text_labels", side_effect=[genes, labels]), \
                    patch.object(np, "load", side_effect=AssertionError("values loaded")):
                with self.assertRaisesRegex(ValueError, "row labels|gene-axis"):
                    A.load_used(d)

    def test_product_exact_alignment_and_schema(self):
        d, z, acc = fixture()
        self.assertEqual(A.validate_product(z, d, acc)["validity"], "VALID")
        for key in ("axis", "symbols", "rows", "obs_id"):
            bad = dict(z, **{key: z[key][::-1]})
            with self.assertRaises(ValueError, msg=key):
                A.validate_product(bad, d, acc)
        for mask in (z["mask"].astype(int), np.full(z["mask"].shape, .5), np.zeros_like(z["mask"])):
            with self.assertRaises(ValueError):
                A.validate_product(dict(z, mask=mask), d, acc)
        bad = dict(z)
        del bad["mask"]
        with self.assertRaises(ValueError):
            A.validate_product(bad, d, acc)

    def test_product_nonfinite_and_wrong_boolean_refused(self):
        d, z, acc = fixture()
        for key in ("delta", "b_ctrl", "control_cpm_sum", "perturbed_cpm_sum"):
            bad = copy.deepcopy(z)
            bad[key].flat[0] = np.nan
            with self.assertRaises(ValueError, msg=key):
                A.validate_product(bad, d, acc)
        with self.assertRaises(ValueError):
            A.validate_product(z, d, dict(acc, cpm_gate_pass=1))
        with self.assertRaises(ValueError):
            A.validate_product(dict(z, control_rows_identical=np.array([1])), d, acc)

    def test_cpm_failure_voids_inference_and_full_report(self):
        d, z, acc = fixture()
        z["control_cpm_sum"][0] = 9998
        acc["cpm_gate_pass"] = False
        validity = A.validate_product(z, d, acc)
        self.assertEqual(validity["validity"], "VOID")
        kd = -z["delta"].astype(float)
        features = {k: np.arange(8.) for k in ("C_excl", "C_held", "gen_pc20", "r_bar", "breadth", "mag290")}
        used = dict(KD=kd, bank_mean=kd.mean(0), features=features, feature_provenance={},
                    Sg302=np.arange(8.), all302=np.ones(8, bool), primary302=np.ones(8, bool))
        report, rows, arrays = A.analyse_arrays(z, d, used, validity, draws=500, boot_draws=100)
        for cell in report["cells"].values():
            self.assertEqual(cell["verdict"], "VOID")
        self.assertEqual(report["cells"]["C7a"]["statistical_verdict"], "PASS")
        text = A.render_results(report)
        self.assertIn("VOID", text)
        self.assertIn("p_median_less", text)
        self.assertIn("p_tail_greater", text)
        self.assertIn("ALL58", text)
        self.assertIn("primary45", text)
        self.assertEqual(len(rows), 8)
        self.assertEqual(arrays["null_S"].shape, (500, 8))
        json.dumps(A.json_clean(report), allow_nan=False)

    def test_missing_access_refuses_without_array_load(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(A, "WORK", Path(tmp)):
            with patch.object(np, "load", side_effect=AssertionError("should not load")):
                self.assertEqual(A.main(["--analyse"]), 2)

    def test_product_hash_mismatch_before_npz_decode(self):
        # Temporary synthetic artifacts ONLY; no real marker/release/product.
        d, z, acc = fixture()
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            prereg = work / "PREREGISTRATION.md"
            prereg.write_text("synthetic test")
            man = dict(files={A.SOURCE_KEY: "a"*64, "PREREGISTRATION.md": A.sha_file(prereg)})
            rel = dict(schema_version=1, claim="C7", reserve="R5",
                       authorization="ONE_SHOT_R5_EXPRESSION", manifest=man)
            (work / "RELEASE.json").write_text(json.dumps(rel))
            rh = A.sha_file(work / "RELEASE.json")
            mh = A.manifest_sha256(man)
            marker = dict(claim="C7", release=rel, release_sha256=rh,
                          live=dict(manifest_sha256=mh, source_sha256="a"*64))
            (work / "ACCESS-MARKER.json").write_text(json.dumps(marker))
            acc.update(manifest_sha256=mh, release_sha256=rh, product_sha256="0"*64)
            (work / "access303.json").write_text(json.dumps(acc))
            (work / "nak303.npz").write_bytes(b"not a product")
            with patch.object(A, "verify_manifest", return_value={}), patch.object(A, "HERE", work), \
                    patch.object(A, "ROOT", work), patch.object(np, "load", side_effect=AssertionError("decoded")):
                with self.assertRaisesRegex(ValueError, "product hash mismatch"):
                    A.load_completed_product(d, work)

    def test_manifest_missing_readout_binding_refused(self):
        with patch.object(A, "input_hashes", return_value={"used.npy": "a"*64}):
            with self.assertRaisesRegex(ValueError, "unbound"):
                A.verify_manifest({"files": {}})

    def test_saved_product_success_and_postaccess_mismatch(self):
        d, z, acc = fixture()
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            prereg = work / "PREREGISTRATION.md"
            prereg.write_text("synthetic test")
            man = dict(files={A.SOURCE_KEY: "a"*64, "PREREGISTRATION.md": A.sha_file(prereg)})
            rel = dict(schema_version=1, claim="C7", reserve="R5",
                       authorization="ONE_SHOT_R5_EXPRESSION", manifest=man)
            (work / "RELEASE.json").write_text(json.dumps(rel))
            rh, mh = A.sha_file(work / "RELEASE.json"), A.manifest_sha256(man)
            marker = dict(claim="C7", release=rel, release_sha256=rh,
                          live=dict(manifest_sha256=mh, source_sha256="a"*64))
            (work / "ACCESS-MARKER.json").write_text(json.dumps(marker))
            np.savez_compressed(work / "nak303.npz", **z)
            acc.update(manifest_sha256=mh, release_sha256=rh,
                       product_sha256=A.sha_file(work / "nak303.npz"))
            (work / "access303.json").write_text(json.dumps(acc))
            with patch.object(A, "verify_manifest", return_value={}), patch.object(A, "HERE", work), \
                    patch.object(A, "ROOT", work):
                product, validity, provenance = A.load_completed_product(d, work)
                np.testing.assert_array_equal(product["delta"], z["delta"])
                self.assertEqual(validity["validity"], "VALID")
                self.assertEqual(provenance["product_sha256"], acc["product_sha256"])
                acc["release_sha256"] = "f"*64
                (work / "access303.json").write_text(json.dumps(acc))
                with self.assertRaisesRegex(ValueError, "release mismatch"):
                    A.load_completed_product(d, work)

    def test_full_driver_writes_only_temporary_outputs(self):
        d, z, acc = fixture()
        validity = A.validate_product(z, d, acc)
        kd = -z["delta"].astype(float)
        features = {k: np.arange(8.) for k in ("C_excl", "C_held", "gen_pc20", "r_bar", "breadth", "mag290")}
        used = dict(KD=kd, bank_mean=kd.mean(0), features=features, feature_provenance={},
                    Sg302=np.arange(8.), all302=np.ones(8, bool), primary302=np.ones(8, bool))
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            (work / "access303.json").write_text("{}")
            with patch.object(A, "WORK", work), patch.object(A, "HERE", work), \
                    patch.object(A, "load_design", return_value=d), \
                    patch.object(A, "load_used", return_value=used), \
                    patch.object(A, "input_hashes", return_value={}), \
                    patch.object(A, "load_completed_product", return_value=(z, validity, {"input_sha256": {}})), \
                    patch("sys.stdout", new_callable=io.StringIO):
                self.assertEqual(A.main(["--analyse"]), 0)
            result = A.read_json(work / "stats303.json")
            self.assertEqual(result["cells"]["C7a"]["verdict"], "PASS")
            self.assertEqual(len(A.read_table(work / "symmetry303.csv")), 8)
            with np.load(work / "analysis303.npz", allow_pickle=False) as arrays:
                self.assertEqual(arrays["null_S"].shape, (2000, 8))
            self.assertIn("Full readout", (work / "RESULTS.md").read_text())
            self.assertFalse((work / "RELEASE.json").exists())
            self.assertFalse((work / "ACCESS-MARKER.json").exists())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    result = unittest.TextTestRunner(verbosity=2, stream=sys.stderr).run(suite)
    if not result.wasSuccessful():
        return 1
    if args.manifest:
        manifest = A.read_json(args.manifest)
        A.verify_manifest(manifest)
        # Frozen features are tested on USED data only, not on R5 arrays.
        design = A.load_design()
        used = A.load_used(design)
        if np.isfinite(used["features"]["C_excl"]).sum() != 97:
            raise AssertionError("frozen C_excl coverage mismatch")
        print(json.dumps(dict(schema_version=1, manifest_sha256=A.manifest_sha256(manifest),
                              suite="frozen_analysis303", passed=True,
                              tests_run=result.testsRun, failures=0, errors=0)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

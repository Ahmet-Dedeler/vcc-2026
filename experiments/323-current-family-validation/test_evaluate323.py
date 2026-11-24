import tempfile
from pathlib import Path
import unittest

import h5py
import numpy as np

from evaluate323 import labels


class LabelTests(unittest.TestCase):
    def test_h5ad_string_encodings(self):
        with tempfile.TemporaryDirectory() as tmp:
            with h5py.File(Path(tmp) / "labels.h5", "w") as h:
                dt = h5py.string_dtype("utf-8")
                h.create_dataset("plain", data=["A", "B"], dtype=dt)
                cat = h.create_group("categorical")
                cat.create_dataset("categories", data=["B", "A"], dtype=dt)
                cat.create_dataset("codes", data=[1, 0])
                nullable = h.create_group("nullable")
                nullable.create_dataset("values", data=["A", "B"], dtype=dt)
                nullable.create_dataset("mask", data=[False, False])
                for name in ("plain", "categorical", "nullable"):
                    np.testing.assert_array_equal(labels(h, name), ["A", "B"])
                    self.assertEqual(labels(h, name).dtype.kind, "U")
                nullable["mask"][0] = True
                with self.assertRaisesRegex(ValueError, "missing"):
                    labels(h, "nullable")
                cat["codes"][0] = -1
                with self.assertRaisesRegex(ValueError, "missing"):
                    labels(h, "categorical")

    def test_exact_de_matches_scanpy(self):
        import anndata as ad
        import pandas as pd
        import polars as pl
        import scipy.sparse as sp
        from cell_eval2.de_compute import compute_de
        from exact_de323 import control_block, de_from_counts

        rng = np.random.default_rng(323)
        counts = rng.poisson(3, (80, 20)).astype(np.float32)
        counts[60:, :4] *= 3
        genes = np.array([f"G{i}" for i in range(20)])
        a = ad.AnnData(
            X=sp.csr_matrix(counts),
            obs=pd.DataFrame({"target_gene": ["non-targeting"]*60 + ["G0"]*20},
                             index=[str(i) for i in range(80)]),
            var=pd.DataFrame(index=genes))
        official = compute_de(
            a, backend="scanpy", groupby="target_gene", reference="non-targeting",
            mean_calc="arithmetic", epsilon=1e-9, input_type="counts",
            target_sum=1e6, clip_value=None, fdr_scope="per_pert",
            filter_gene_min_cpm_cell=5., threads=1)
        c = control_block(a.X[:60])
        _, p, q, lfc = de_from_counts(a.X[60:], c)
        fast = pl.DataFrame({"feature": genes[c.keep], "p": p, "q": q, "lfc": lfc})
        joined = official.join(fast, on="feature")
        for a, b in (("p_value", "p"), ("p_adj", "q"), ("log2_fold_change", "lfc")):
            self.assertLess(float((joined[a]-joined[b]).abs().max()), 1e-12)


if __name__ == "__main__":
    unittest.main()

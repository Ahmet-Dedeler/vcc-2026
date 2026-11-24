import io
import pickle
import unittest
import zipfile
import numpy as np
import acquire333 as A


class AcquisitionTests(unittest.TestCase):
    def test_numpy_only_pickle_roundtrip(self):
        mapping = {"ABC": np.arange(12, dtype=np.float32)}
        result = A.RestrictedNumpyUnpickler(io.BytesIO(pickle.dumps(mapping, protocol=4))).load()
        np.testing.assert_array_equal(result["ABC"], mapping["ABC"])

    def test_unknown_global_rejected(self):
        class Probe:
            def __reduce__(self):
                return (eval, ("1+1",))
        with self.assertRaises(pickle.UnpicklingError):
            A.RestrictedNumpyUnpickler(io.BytesIO(pickle.dumps(Probe()))).load()

    def test_schema_and_nonfinite(self):
        genes, x = A.validate_mapping({"BBB": [1., 2.], "AAA": [3., 4.]})
        self.assertEqual(genes.tolist(), ["AAA", "BBB"])
        self.assertEqual(x.dtype, np.float32)
        with self.assertRaises(ValueError):
            A.validate_mapping({"ABC": [1., np.nan]})
        with self.assertRaises(ValueError):
            A.validate_mapping({"ABC": [0., 0.]})

    def test_zip_path_and_exact_member(self):
        good = zipfile.ZipInfo("GenePT/" + A.MEMBER)
        self.assertEqual(A.select_member([good]), good)
        with self.assertRaises(ValueError):
            A.select_member([good, zipfile.ZipInfo("../escape")])
        with self.assertRaises(ValueError):
            A.select_member([good, good])


if __name__ == "__main__":
    unittest.main()

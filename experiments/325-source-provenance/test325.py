import unittest

from fetch325 import fields_encode, fields_soft, publication_and_samples, source_url


class MetadataTests(unittest.TestCase):
    def test_soft_preserves_repeated_fields(self):
        text = ("^SERIES = GSE1\n!Series_geo_accession = GSE1\n"
                "!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n")
        self.assertEqual(fields_soft(text, "GSE1")["Series_sample_id"], ["GSM1", "GSM2"])

    def test_soft_id_mismatch_rejected(self):
        with self.assertRaises(ValueError):
            fields_soft("^SERIES = GSE2\n!Series_geo_accession = GSE2\n", "GSE1")

    def test_encode_id_mismatch_rejected(self):
        with self.assertRaises(ValueError):
            fields_encode('{"accession":"ENCSR123ABC"}', "ENCSR999XYZ")

    def test_no_arbitrary_url_or_sample_download(self):
        for value in ("https://example.com", "GSM123", "../../secret", "GSE1;evil"):
            with self.assertRaises(ValueError):
                source_url(value)

    def test_same_publication_is_not_same_samples(self):
        a = {"kind": "GEO", "fields": {"Series_pubmed_id": ["42"], "Series_sample_id": ["GSM1"]}}
        b = {"kind": "GEO", "fields": {"Series_pubmed_id": ["42"], "Series_sample_id": ["GSM2"]}}
        pa, sa = publication_and_samples(a)
        pb, sb = publication_and_samples(b)
        self.assertEqual(pa & pb, {"42"})
        self.assertEqual(sa & sb, set())


if __name__ == "__main__":
    unittest.main()

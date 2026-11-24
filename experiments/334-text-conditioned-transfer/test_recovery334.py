from run_recovery334 import sources


def test_primary_priority():
    docs, source = sources({"A":"Gene Symbol A original"},
                           {"A":"Gene Symbol A alternative"}, ["A"])
    assert docs == ["original"] and source == ["NCBI"]


def test_header_only_fallback():
    docs, source = sources({"A":"Gene Symbol A"},
                           {"A":"Gene Symbol A Protein summary: function"}, ["A"])
    assert docs == ["Protein summary: function"]
    assert source == ["NCBI_UniProt_fallback"]


def test_no_alias_or_wrong_header():
    docs, source = sources({}, {"A":"Gene Symbol B wrong",
                               "alias":"Gene Symbol alias irrelevant"}, ["A"])
    assert docs == [""] and source == ["missing"]

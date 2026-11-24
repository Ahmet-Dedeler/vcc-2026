from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


extract = load_module("extract_058", "extract.py")
evaluate = load_module("evaluate_058", "evaluate_lodo.py")


def compressed_fixture() -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        text = io.TextIOWrapper(handle, encoding="utf-8", newline="")
        writer = csv.writer(text)
        writer.writerow(["gene", "A_K5", "B_K5", "C_K5"])
        writer.writerow(["G1", "1", "2", "3"])
        writer.writerow(["G2", "-1", "0", "1"])
        text.flush()
        text.detach()
    return buffer.getvalue()


def test_streaming_extractor_layout_and_hash(tmp_path: Path):
    source = compressed_fixture()
    output = tmp_path / "compact"
    metadata = extract.extract_stream(
        io.BytesIO(source),
        output,
        ("A", "B", "X"),
        expected_genes=2,
        expected_bytes=len(source),
        expected_md5=hashlib.md5(source).hexdigest(),
    )
    matrix = np.load(output / "delta.npy", mmap_mode="r")
    assert matrix.shape == (3, 2)
    np.testing.assert_allclose(matrix, [[1, -1], [2, 0], [3, 1]])
    assert matrix.flags.f_contiguous
    assert (output / "perts.txt").read_text().splitlines() == ["A", "B", "C"]
    assert metadata["vcc_target_coverage"] == 2
    assert (
        json.loads((output / "meta.json").read_text())["source_md5"]
        == hashlib.md5(source).hexdigest()
    )


def test_pds_uses_complete_gallery():
    truth = np.eye(4)
    prediction = truth.copy()
    assert np.all(evaluate.pds_per_query(prediction, truth) == 1)
    prediction[[0, 1]] = prediction[[1, 0]]
    scores = evaluate.pds_per_query(prediction, truth)
    np.testing.assert_allclose(scores, [1 / 3, 1 / 3, 1, 1])


def test_candidate_algebra():
    prepared = {}
    for name, value in {
        "k562_gwps": 1.0,
        "perturbmatch": 3.0,
        "orion_hct116": 5.0,
        "orion_hek293t": 7.0,
    }.items():
        matrix = np.full((2, 3), value)
        prepared[name] = type(
            "PreparedStub",
            (),
            {"source_raw": matrix, "source_gene_norm": matrix},
        )()
    found = evaluate.candidates(prepared, "orion_hct116", "raw")
    np.testing.assert_allclose(found["baseline_old_k562_plus_orion"], 4)
    np.testing.assert_allclose(found["replacement_pm_plus_orion"], 5)
    np.testing.assert_allclose(found["k562_ensemble_f0.25"], 1.5)
    np.testing.assert_allclose(found["full_ensemble_f0.25"], 4.25)


def test_panel_coordinates_are_excluded_by_construction(tmp_path: Path):
    official = ("PANEL", "G1", "G2")
    panel = {"PANEL"}
    features = tuple(gene for gene in official if gene not in panel)
    assert features == ("G1", "G2")
    assert not panel.intersection(features)

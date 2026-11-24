from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


EXTRACT = load_module("extract_kolf", "extract_kolf.py")
EVALUATE = load_module("evaluate_kolf_lodo", "evaluate_lodo.py")


def test_contract_accepts_exact_official_shape_and_code_markers() -> None:
    meta = {
        "perts": [f"P{i}" for i in range(10_167)],
        "genes": [f"G{i}" for i in range(38_606)],
        "scale": 0.04,
        "chunk_size": 2_500,
        "n_chunks": 5,
        "n_perts": 10_167,
        "n_genes": 38_606,
        "row_bytes": 38_606,
    }
    builder = "\n".join(
        (
            "int8 row-major",
            "q = np.clip(np.round(M / SCALE), -127, 127).astype(np.int8)",
            "chunk.tobytes(order='C')",
            "'row_bytes':  int(n_genes)",
        )
    ).encode()
    loader = "\n".join(
        (
            "Rows are fetched on demand via HTTP Range",
            "const offset   = (pi % cs) * rb;",
            "Range: `bytes=${offset}-${offset + rb - 1}`",
            "row = new Int8Array(buf);",
            "row[gi] * scale",
        )
    ).encode()
    assert (
        EXTRACT.validate_official_contract(json.dumps(meta).encode(), builder, loader)[
            "n_perts"
        ]
        == 10_167
    )
    meta["row_bytes"] = 38_605
    try:
        EXTRACT.validate_official_contract(json.dumps(meta).encode(), builder, loader)
    except ValueError as exc:
        assert "row_bytes" in str(exc)
    else:
        raise AssertionError("malformed row contract was accepted")


def test_range_client_rejects_full_response_before_read() -> None:
    class Response:
        status = 200

        def getheaders(self):
            return [("Content-Length", "96515000")]

        def read(self, _size):
            raise AssertionError("full response body must not be read")

        def close(self):
            pass

    class Connection:
        def request(self, *_args, **_kwargs):
            pass

        def getresponse(self):
            return Response()

        def close(self):
            pass

    client = EXTRACT.StrictRangeClient(
        "https://example.test/root", timeout=1, retries=0
    )
    client.connection = Connection()
    try:
        client.get("chunk_0.bin", 0, 9, 100)
    except RuntimeError as exc:
        assert "HTTP 200" in str(exc)
    else:
        raise AssertionError("full-body response was accepted")


def test_decode_scale_and_exact_axis_mapping() -> None:
    source = np.asarray([[-127, -1, 0, 2, 127]], dtype=np.int8)
    genes = ["A", "B", "C", "D", "E"]
    wanted = ["E", "B", "D"]
    columns = np.asarray([genes.index(gene) for gene in wanted])
    decoded = source[:, columns].astype(np.float32) * 0.04
    np.testing.assert_allclose(decoded, [[5.08, -0.04, 0.08]])


def write_bank(
    path: Path,
    genes: list[str],
    perts: list[str],
    delta: np.ndarray,
    *,
    kolf: bool,
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", delta.astype(np.float32))
    np.save(path / "b_ctrl.npy", np.ones(len(genes), dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes) + "\n")
    (path / "perts.txt").write_text("\n".join(perts) + "\n")
    if kolf:
        (path / "meta.json").write_text(
            json.dumps(
                {
                    "effect_definition": "delta = int8 source value * 0.04",
                    "vcc_target_coverage": len(perts),
                }
            )
        )


def test_lodo_excludes_panel_and_freezes_before_k562(tmp_path: Path) -> None:
    panel = ["P1", "P2", "P3", "P4"]
    genes = [*panel, "G1", "G2", "G3", "G4"]
    pd.DataFrame({"target_gene": panel}).to_csv(tmp_path / "panel.csv", index=False)
    pd.DataFrame({"gene_name": genes}).to_csv(tmp_path / "genes.csv", index=False)
    base = np.zeros((4, 8), dtype=np.float32)
    base[:, 4:] = np.asarray([[3, 0, 0, 0], [0, 3, 0, 0], [0, 0, 3, 0], [0, 0, 0, 3]])
    main_paths = {}
    for index, name in enumerate(("k562_gwps", "orion_hct116", "orion_hek293t")):
        path = tmp_path / name
        write_bank(path, genes, panel, base + index * 0.01, kolf=False)
        main_paths[name] = path
    kolf = tmp_path / "kolf"
    write_bank(kolf, genes, panel, base, kolf=True)

    old_defaults = EVALUATE.MAIN_DEFAULTS
    EVALUATE.MAIN_DEFAULTS = main_paths
    try:
        result = EVALUATE.run(
            argparse.Namespace(
                panel_csv=str(tmp_path / "panel.csv"),
                genes_csv=str(tmp_path / "genes.csv"),
                kolf=str(kolf),
                kolf_fractions=[0.1, 0.25, 0.5],
                chunk_rows=2,
                max_rss_gib=0,
                min_shared_targets=4,
                min_shared_genes=4,
                output=str(tmp_path / "results.json"),
            )
        )
    finally:
        EVALUATE.MAIN_DEFAULTS = old_defaults
    assert result["protocol"]["genes"] == 4
    assert result["protocol"]["retrieval_gallery_size"] == 4
    assert result["protocol"]["k562_candidate_evaluations_after_selection"] == 1
    assert len(result["development_rows"]) == 12
    assert result["frozen_k562"]["fraction"] in {0.1, 0.25, 0.5}


def test_development_selection_never_reads_k562_rows() -> None:
    rows = []
    for mode in ("raw", "gene_norm"):
        for held in EVALUATE.DEVELOPMENT_HELD:
            for fraction in (0.1, 0.25, 0.5):
                rows.append(
                    {
                        "mode": mode,
                        "held": held,
                        "fraction": fraction,
                        "delta_vs_baseline": (
                            1.0 if (mode, fraction) == ("gene_norm", 0.25) else 0.0
                        ),
                    }
                )
    selected = EVALUATE.choose_development_candidate(rows)
    assert selected["mode"] == "gene_norm"
    assert selected["fraction"] == 0.25


if __name__ == "__main__":
    test_contract_accepts_exact_official_shape_and_code_markers()
    test_range_client_rejects_full_response_before_read()
    test_decode_scale_and_exact_axis_mapping()
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        test_lodo_excludes_panel_and_freezes_before_k562(Path(directory))
    test_development_selection_never_reads_k562_rows()
    print("5 synthetic checks passed")

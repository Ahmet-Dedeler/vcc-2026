from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse as sp


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp081_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_hash_split_is_deterministic_disjoint_and_complete() -> None:
    targets = tuple(f"G{i}" for i in range(500))
    first = RUN.split_targets(targets)
    second = RUN.split_targets(tuple(reversed(targets)))
    assert first == second
    sets = [set(first[name]) for name in ("fit", "tune", "test")]
    assert not (sets[0] & sets[1] or sets[0] & sets[2] or sets[1] & sets[2])
    assert set.union(*sets) == set(targets)


def test_streamed_npz_operator_matches_dense() -> None:
    work = HERE / ".test-tmp"
    work.mkdir(exist_ok=True)
    path = work / "operator.npz"
    metadata = work / "operator.json"
    rng = np.random.default_rng(81)
    q, _ = np.linalg.qr(rng.normal(size=(9, 9)))
    evecs = np.asarray(q, dtype=np.float32)
    evals = np.linspace(0.2, 2.0, 9, dtype=np.float64)
    np.savez(path, evals=evals, evecs=evecs)
    wcols = np.arange(9, dtype=np.int64)
    metadata.write_text(
        json.dumps(
            {
                "donors": ["a", "b"],
                "gamma": RUN.WHITEN_GAMMA,
                "n_cols": 9,
                "wcols_sha256": RUN.hashlib.sha256(wcols.tobytes()).hexdigest(),
            }
        )
    )
    operator = RUN.FrozenOperator(
        path, metadata, RUN.sha256_file(path), ["a", "b"], wcols
    )
    matrix = rng.normal(size=(7, 9)).astype(np.float32)
    observed, details = operator.apply(
        matrix, block_rows=3, max_rss_gib=0.99
    )
    weights = evals ** (-RUN.WHITEN_GAMMA / 2.0)
    expected = ((matrix @ evecs) * weights) @ evecs.T
    expected *= np.median(np.linalg.norm(matrix, axis=1)) / np.median(
        np.linalg.norm(expected, axis=1)
    )
    assert np.allclose(observed, expected, rtol=2e-5, atol=2e-5)
    assert details["applied_scale"] > 0
    for item in work.iterdir():
        item.unlink()
    work.rmdir()


def test_igvf_parser_aggregates_duplicate_pairs_and_reads_tsv() -> None:
    work = HERE / ".test-tmp"
    work.mkdir(exist_ok=True)
    de = work / "de.csv.gz"
    guide = work / "guide.csv.gz"
    with gzip.open(de, "wt") as handle:
        handle.write(
            "gene_names,intended_target_name,num_cell,log2fc,pval-empirical\n"
            "R1,T1,10,1.0,0\n"
            "R1,T1,30,3.0,0\n"
            "R2,T1,40,-2.0,0\n"
        )
    with gzip.open(guide, "wt") as handle:
        handle.write("guide_id\ttargeting\ttype\tgene_name\n")
        handle.write("g1\tTRUE\ttarget\tT1\n")
    parsed = RUN.load_igvf(de, guide, ("R1", "R2"))
    assert parsed.de_rows == 3
    assert parsed.guide_rows == 1
    assert parsed.duplicate_pairs == 1
    assert parsed.effects["T1"][0] == 2.5
    assert parsed.effects["T1"][1] == -2.0
    for item in work.iterdir():
        item.unlink()
    work.rmdir()


def test_sparse_adapter_preserves_missing_targets_and_genes() -> None:
    baseline = np.arange(24, dtype=np.float32).reshape(4, 6) / 10
    truth = np.eye(4, 6, dtype=np.float32)
    truth /= np.linalg.norm(truth, axis=1, keepdims=True)
    residual = sp.csr_matrix(
        (
            np.asarray([1.0, -2.0], dtype=np.float32),
            (np.asarray([0, 2]), np.asarray([1, 4])),
        ),
        shape=(4, 6),
    )
    base_norm_sq = np.sum(baseline.astype(np.float64) ** 2, axis=1)
    pack = RUN.FoldPack(
        held="x",
        targets=("A", "B", "C", "D"),
        wcols=np.arange(6),
        retrieval_positions=np.arange(6),
        retrieval_official_positions=np.arange(6),
        baseline_pre_norm=np.linalg.norm(baseline, axis=1),
        baseline_retrieval=baseline,
        truth_unit=truth,
        base_dot=baseline @ truth.T,
        base_norm_sq=base_norm_sq,
        baseline_pds=RUN.pds_from_similarities(
            (baseline @ truth.T) / np.sqrt(base_norm_sq)[:, None]
        ),
        source_similarity=np.eye(4, dtype=np.float32),
        raw_by_representation={"log2fc": residual},
        whitening={},
        source_names=("a", "b"),
        whitener_sha256="x",
    )
    candidate_pds = RUN.evaluate_residual(pack, residual, 0.1)
    assert candidate_pds[1] == pack.baseline_pds[1]
    assert candidate_pds[3] == pack.baseline_pds[3]
    dense_candidate = baseline + 0.1 * residual.toarray()
    assert np.array_equal(
        dense_candidate[residual.toarray() == 0],
        baseline[residual.toarray() == 0],
    )


def test_scale_modes_leave_zero_rows_zero() -> None:
    raw = sp.csr_matrix(
        np.asarray([[1.0, 0.0], [0.0, 0.0], [0.0, -2.0]], dtype=np.float32)
    )
    baseline = np.ones((3, 2), dtype=np.float32)
    for mode in RUN.SCALE_MODES:
        params = RUN.fit_scale(
            sp.vstack([raw] * 50).tocsr(),
            np.linalg.norm(np.vstack([baseline] * 50), axis=1),
            mode,
        )
        scaled = RUN.scale_sparse(
            raw, np.linalg.norm(baseline, axis=1), params
        )
        assert scaled.getrow(1).nnz == 0


def test_nulls_preserve_sparsity_contract() -> None:
    raw = sp.csr_matrix(
        np.asarray(
            [[1.0, 0.0, -2.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.0]],
            dtype=np.float32,
        )
    )
    randomized = RUN.sign_randomize(raw, np.random.default_rng(3))
    assert np.array_equal(raw.indptr, randomized.indptr)
    assert np.array_equal(raw.indices, randomized.indices)
    assert np.array_equal(np.abs(raw.data), np.abs(randomized.data))
    strata = np.asarray([111, 111, 0])
    permutation = RUN.matched_permutation(strata, np.random.default_rng(4))
    assert permutation[2] == 2
    assert set(permutation[:2]) == {0, 1}


def test_panel_coordinates_are_removed_from_sparse_matrix() -> None:
    igvf = RUN.IGVFData(
        effects={"T": {0: 1.0, 1: -2.0}},
        num_cells={"T": 20},
        row_counts={"T": 2},
        de_targets=frozenset({"T"}),
        guide_targets=frozenset({"T"}),
        de_rows=2,
        guide_rows=1,
        response_genes=frozenset({"PANEL", "OK"}),
        duplicate_pairs=0,
    )
    matrix = RUN.igvf_matrix(
        igvf,
        ("T",),
        np.asarray([0, 1], dtype=np.int64),
        {0},
        "log2fc",
        10.0,
    )
    assert matrix[0, 0] == 0
    assert matrix[0, 1] == -2


def test_output_scope_rejects_other_paths() -> None:
    RUN.assert_output_scope(HERE / "ok.json")
    try:
        RUN.assert_output_scope(HERE.parent / "not-ok.json")
    except ValueError:
        pass
    else:
        raise AssertionError("outside output path was accepted")

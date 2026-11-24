"""235's regression tests — the Oct-22 panel guard, and the coverage-aware sizing rule.

Run them before Oct 22 and before any build on a new panel:

    uv run python experiments/235-october-hardening/tests/test_october_hardening.py

Same shape as `scripts/test_make_submission_pm_cd4.py` (a `TESTS` list and a `__main__`
runner) because `pytest` is not installed in this repo's `uv` environment.  They live here
rather than in that suite only because 235's mandate authorises edits to
`scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py` and nothing else
under `scripts/`; folding them in is a one-move follow-up for an agent with that permission.

The first six are synthetic and take about a second.  `test_every_registered_source...` reads
the real banks and takes a couple of minutes; it is the one that has to pass on Oct 22.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
ROOT = EXP.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import make_submission_pm_cd4 as M  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


CARRY = _load(EXP / "carry235.py", "carry235")
FITTER = _load(ROOT / "scripts/build_union_whitener.py", "build_union_whitener_under_test")

PRODUCTION_OPERATOR = ROOT / "data/donors/whitener_union_pruned.npz"
PRODUCTION_OPERATOR_SHA = (
    "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
)
PRODUCTION_AXIS_SHA = (
    "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
)


# --------------------------------------------------------------------------------- item A
class StubModel:
    """The three attributes `production_operator_columns` and the guard actually read."""

    def __init__(self, genes):
        self.target_genes = np.asarray(list(map(str, genes)))
        self.names = ["stub_donor"]
        self.tcols = {"stub_donor": np.arange(len(self.target_genes), dtype=np.int64)}


def _write_operator(directory: Path, n_cols: int, axis_sha: str | None) -> tuple[Path, Path]:
    rng = np.random.default_rng(0)
    evecs = np.linalg.qr(rng.standard_normal((n_cols, n_cols)))[0]
    evals = np.linspace(1.0, 2.0, n_cols)
    npz = directory / "op.npz"
    Whitener(evals, evecs).save(npz)
    meta: dict[str, object] = {
        "n_cols": n_cols,
        "basis": [10 * n_cols, n_cols],
        "donors": ["stub_donor"],
    }
    if axis_sha is not None:
        meta["axis_sha256"] = axis_sha
    js = directory / "op.json"
    js.write_text(json.dumps(meta))
    return npz, js


def test_the_two_axis_digest_implementations_are_the_same_function() -> None:
    """`build_union_whitener.axis_sha256` writes what the builder asserts. One function."""
    cases = [
        [],
        ["A"],
        ["TSPAN6", "TNMD", "DPM1"],
        ["DPM1", "TNMD", "TSPAN6"],           # order matters
        ["AÅ", "B", "C"],                # non-ASCII symbols
        [f"GENE{i}" for i in range(5000)],
    ]
    for case in cases:
        assert FITTER.axis_sha256(case) == M.operator_axis_sha256(case), case
    # and it is genuinely order-sensitive, which is the whole point
    assert M.operator_axis_sha256(["A", "B"]) != M.operator_axis_sha256(["B", "A"])


def test_the_axis_guard_fires_on_a_panel_the_operator_was_not_fitted_for() -> None:
    """Same donor, same `n_cols`, different panel -> must raise, not build."""
    genes = [f"G{i:02d}" for i in range(8)]
    model = StubModel(genes)
    panel_fitted = ("G00", "G01")
    panel_other = ("G06", "G07")        # same size, so the axis is the same LENGTH
    fitted_axis = [g for g in genes if g not in panel_fitted]
    other_axis = [g for g in genes if g not in panel_other]
    assert len(fitted_axis) == len(other_axis)          # the blind spot, reproduced
    assert fitted_axis != other_axis

    with tempfile.TemporaryDirectory() as tmp:
        npz, js = _write_operator(
            Path(tmp), len(fitted_axis), M.operator_axis_sha256(fitted_axis)
        )
        loaded = M.load_production_operator(
            model, panel_fitted, npz, js, expected_sha256=None, whitener_id="stub"
        )
        assert loaded.axis_sha256 == M.operator_axis_sha256(fitted_axis)
        assert loaded.axis_sha256_source.startswith("metadata:")

        try:
            M.load_production_operator(
                model, panel_other, npz, js, expected_sha256=None, whitener_id="stub"
            )
        except ValueError as exc:
            assert "DIFFERENT PANEL" in str(exc)
            assert M.operator_axis_sha256(other_axis) in str(exc)
        else:
            raise AssertionError(
                "the axis guard did NOT fire on a panel the operator was not fitted for"
            )


def test_an_operator_with_no_recorded_axis_is_refused_rather_than_trusted() -> None:
    genes = [f"G{i:02d}" for i in range(8)]
    model = StubModel(genes)
    panel = ("G00", "G01")
    with tempfile.TemporaryDirectory() as tmp:
        npz, js = _write_operator(Path(tmp), len(genes) - len(panel), None)
        assert M.file_sha256(npz) not in M.LEGACY_OPERATOR_AXIS_SHA256
        try:
            M.load_production_operator(
                model, panel, npz, js, expected_sha256=None, whitener_id="stub"
            )
        except ValueError as exc:
            assert "no axis digest" in str(exc)
        else:
            raise AssertionError("an unvalidatable operator was accepted")


def test_the_production_operator_is_pinned_by_path_and_by_sha_rule_11() -> None:
    assert PRODUCTION_OPERATOR.is_file()
    assert M.file_sha256(PRODUCTION_OPERATOR) == PRODUCTION_OPERATOR_SHA
    assert M.LEGACY_OPERATOR_AXIS_SHA256[PRODUCTION_OPERATOR_SHA] == PRODUCTION_AXIS_SHA
    assert M.WHITENER_REGISTRY["union_pruned"][0] == PRODUCTION_OPERATOR_SHA


def test_the_fitter_writes_an_axis_digest_into_every_operator_it_makes() -> None:
    """A refit on the day must produce a metadata json the builder can check."""
    source = (ROOT / "scripts/build_union_whitener.py").read_text()
    assert '"axis_sha256": axis_digest' in source
    assert "def axis_sha256(" in source


# --------------------------------------------------------------------------------- item B
def _prepared(n_rows: int, n_genes: int, rows, cols, value: float):
    effect = np.zeros((n_rows, n_genes), dtype=np.float64)
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in rows:
        effect[row, cols] = value
        votes[row, cols] = 1
        counts[row] = 1
    return M.EXP062.CD4Prepared(
        mode="stub",
        mean_effect=effect,
        gene_votes=votes,
        target_condition_count=counts,
        condition_details=[],
    )


def test_the_second_residual_is_sized_by_the_carrying_row_rule_and_cd4_is_not() -> None:
    """The routing item B has to hold for every source, tested behaviourally.

    The second residual covers 2 of 6 targets -- under 50%, where the shipped all-rows rule
    returns a fraction of EXACTLY ZERO and the source is absent with nothing to read it back.
    """
    n_rows, n_genes = 6, 10
    panel_positions = np.array([0, 1], dtype=np.int64)
    columns = np.arange(2, n_genes, dtype=np.int64)
    rng = np.random.default_rng(11)
    whitened = rng.standard_normal((n_rows, n_genes))
    cd4 = _prepared(n_rows, n_genes, range(n_rows), [2, 3, 4, 5], 0.5)
    extra = _prepared(n_rows, n_genes, [0, 1], [6, 7, 8, 9], 0.5)

    saved = M.H1_UNION_G
    M.H1_UNION_G = 0.28
    try:
        _out, cd4_audit, union_audit = M.apply_frozen_cd4_residual(
            whitened, cd4, columns, panel_positions, extra=extra, extra_audit={"mode": "x"}
        )
    finally:
        M.H1_UNION_G = saved

    assert "medians over all rows" in cd4_audit["residual_frame"]["definition"]
    frame = union_audit["residual_frame"]
    assert "over the rows it CARRIES" in frame["definition"]
    assert frame["carrying_rows"] == 2
    assert frame["residual_median_row_norm_over_all_rows"] == 0.0
    assert frame["shipped_median_rule_would_solve"] == 0.0
    assert frame["solved_fraction"] > 0.0
    assert frame["applied_fraction"] == frame["solved_fraction"]


def test_a_selected_source_that_carries_no_panel_target_now_raises() -> None:
    """225: `kolf2` 273 -> 0 and `feng_ipsc` 182 -> 0 in October. Silent before 235."""
    bundle = ROOT / "data/challenge"
    panel, genes = CARRY.panel_and_genes(bundle)
    name = "nadig_hepg2"                      # 0/300 of today's panel, by 225's census
    assert CARRY.coverage(name, panel) == 0
    saved = M.H1_UNION_SOURCES
    M.H1_UNION_SOURCES = ((name, Path(str(M.UNION_SOURCE_REGISTRY[name]["path"]))),)
    try:
        M.prepare_h1_union(panel, genes, panel, max_rss_gib=24.0)
    except ValueError as exc:
        assert "carries 0 of 300 panel targets" in str(exc)
    else:
        raise AssertionError("a zero-coverage source was accepted silently")
    finally:
        M.H1_UNION_SOURCES = saved


def test_the_shipped_h1_pair_reproduces_the_shipped_build_exactly() -> None:
    """Known-answer row: `submissions/222-h1union-g028.h5ad.meta.json`, to the last bit."""
    shipped = json.loads(CARRY.SHIPPED_META.read_text())
    frame = shipped["effects"]["h1_union"]["residual_frame"]
    panel, genes = CARRY.panel_and_genes(CARRY.PANELS["abc"])
    row = CARRY.size_sources(M.UNION_SOURCES_DEFAULT, panel, genes)
    got = row["frame"]
    for key in (
        "residual_median_row_norm_over_carrying_rows",
        "g_per_unit_fraction",
        "solved_fraction",
        "carrying_rows",
        "coverage_frac",
        "residual_median_row_norm_over_all_rows",
        "shipped_median_rule_would_solve",
    ):
        assert got[key] == frame[key], (key, got[key], frame[key])


def test_every_registered_source_sizes_above_zero_on_a_panel_it_covers() -> None:
    """Item B's claim, for all ten names -- including the October four.

    Each of the October four carries 0/300 today and ~30/300 of `DEF-RUNBOOK` section 5's
    pessimistic draw, which is the panel they are registered for.  If any registered source
    can be selected and still deliver a fraction of zero, this fails.
    """
    checked = 0
    for tag, bundle in CARRY.PANELS.items():
        panel, genes = CARRY.panel_and_genes(bundle)
        for name in sorted(M.UNION_SOURCE_REGISTRY):
            covered = CARRY.coverage(name, panel)
            if covered == 0:
                continue
            row = CARRY.size_sources((name,), panel, genes)
            assert not row["carrying_raised"], (tag, name, row.get("carrying_message"))
            assert float(row["solved_fraction"]) > 0.0, (tag, name, row)
            if covered <= len(panel) // 2:
                assert float(row["all_rows_would_solve"]) == 0.0, (tag, name, row)
            checked += 1
    assert checked >= len(M.UNION_SOURCE_REGISTRY), checked
    for name in CARRY.OCTOBER_FOUR:
        panel, _ = CARRY.panel_and_genes(CARRY.PANELS["def199_pessimistic"])
        assert CARRY.coverage(name, panel) > 0, name


TESTS = [
    test_the_two_axis_digest_implementations_are_the_same_function,
    test_the_axis_guard_fires_on_a_panel_the_operator_was_not_fitted_for,
    test_an_operator_with_no_recorded_axis_is_refused_rather_than_trusted,
    test_the_production_operator_is_pinned_by_path_and_by_sha_rule_11,
    test_the_fitter_writes_an_axis_digest_into_every_operator_it_makes,
    test_the_second_residual_is_sized_by_the_carrying_row_rule_and_cd4_is_not,
    test_a_selected_source_that_carries_no_panel_target_now_raises,
    test_the_shipped_h1_pair_reproduces_the_shipped_build_exactly,
    test_every_registered_source_sizes_above_zero_on_a_panel_it_covers,
]


if __name__ == "__main__":
    failures = 0
    for test in TESTS:
        try:
            test()
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"FAIL {test.__name__}: {type(exc).__name__}: {exc}")
        else:
            print(f"PASS {test.__name__}")
    print(f"{len(TESTS) - failures}/{len(TESTS)} passed")
    raise SystemExit(1 if failures else 0)

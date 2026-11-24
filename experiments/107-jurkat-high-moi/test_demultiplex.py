from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import demultiplex  # noqa: E402


EXPECTED_KERNSMOOTH_VERSION = "2.23.26"


def synthetic_hash_counts() -> tuple[np.ndarray, tuple[str, ...], tuple[str, ...]]:
    """Create deterministic negative, singlet, and doublet populations."""

    rng = np.random.default_rng(107)
    cells = 160
    hashes = 4
    counts = rng.integers(0, 4, size=(cells, hashes), dtype=np.int64)
    for cell in range(80):
        counts[cell, cell % hashes] = rng.integers(48, 150)
    for cell in range(80, 100):
        first = cell % hashes
        counts[cell, first] = rng.integers(50, 140)
        counts[cell, (first + 1) % hashes] = rng.integers(50, 140)
    return (
        counts,
        tuple(f"cell-{index:03d}" for index in range(cells)),
        tuple(f"hash-{index + 1}" for index in range(hashes)),
    )


R_ORACLE = r"""
suppressPackageStartupMessages({
  library(jsonlite)
  library(KernSmooth)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: oracle.R COUNTS.csv OUTPUT.json")
expected <- "2.23.26"
observed <- as.character(packageVersion("KernSmooth"))
if (!identical(observed, expected)) {
  stop(sprintf(
    "KernSmooth version mismatch: observed %s expected %s",
    observed,
    expected
  ))
}

localMaxima <- function(x) {
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) y <- y[-1]
  y
}

normalizeHashes <- function(barTable) {
  normalized <- as.data.frame(log2(barTable))
  for (i in seq_len(ncol(barTable))) {
    bad <- which(is.finite(normalized[, i]) == FALSE)
    normalized[bad, i] <- 0
    normalized[, i] <- normalized[, i] - mean(normalized[, i])
  }
  as.matrix(normalized)
}

analyzeRound <- function(barTable) {
  normalized <- normalizeHashes(barTable)
  q.range <- seq(0.01, 0.99, by=0.02)
  models <- vector("list", ncol(normalized))
  thresholds <- matrix(NA_real_, nrow=length(q.range), ncol=ncol(normalized))

  for (i in seq_len(ncol(normalized))) {
    model <- approxfun(bkde(normalized[, i], kernel="normal"))
    x <- seq(
      from=quantile(normalized[, i], 0.001, type=7),
      to=quantile(normalized[, i], 0.999, type=7),
      length.out=100
    )
    density <- model(x)
    extrema <- localMaxima(density)
    if (length(extrema) <= 1L) stop("undefined hash threshold")
    low <- extrema[which.max(density[extrema])]
    high <- max(extrema)
    if (low == high) stop("identical low/high extrema")
    for (q.index in seq_along(q.range)) {
      thresholds[q.index, i] <- quantile(
        c(x[high], x[low]), q.range[[q.index]], type=7
      )
    }
    models[[i]] <- list(
      evaluation_points=unname(x),
      density=unname(density),
      extrema=unname(extrema - 1L),
      low_extreme=unname(low - 1L),
      high_extreme=unname(high - 1L)
    )
  }

  calls.by.q <- vector("list", length(q.range))
  singlet.fraction <- numeric(length(q.range))
  for (q.index in seq_along(q.range)) {
    hits <- sweep(
      normalized,
      2L,
      thresholds[q.index, ],
      FUN=">="
    )
    hit.count <- rowSums(hits)
    calls <- rep("Negative", nrow(normalized))
    calls[hit.count > 1L] <- "Doublet"
    singlet <- which(hit.count == 1L)
    if (length(singlet)) {
      calls[singlet] <- colnames(normalized)[max.col(hits[singlet, , drop=FALSE])]
    }
    calls.by.q[[q.index]] <- unname(calls)
    singlet.fraction[[q.index]] <- mean(hit.count == 1L)
  }

  q.maxima <- localMaxima(singlet.fraction)
  if (!length(q.maxima)) stop("undefined q")
  selected.index <- q.maxima[which.max(singlet.fraction[q.maxima])]
  list(
    normalized=unname(split(normalized, row(normalized))),
    models=models,
    q_values=unname(q.range),
    thresholds=unname(split(thresholds, row(thresholds))),
    calls_by_q=calls.by.q,
    singlet_fraction=unname(singlet.fraction),
    q_local_maxima=unname(q.maxima - 1L),
    selected_q_index=unname(selected.index - 1L),
    selected_q=unname(q.range[[selected.index]]),
    selected_calls=calls.by.q[[selected.index]]
  )
}

barTable <- as.matrix(read.csv(args[[1L]], row.names=1L, check.names=FALSE))
storage.mode(barTable) <- "double"
active <- seq_len(nrow(barTable))
final <- rep(NA_character_, nrow(barTable))
rounds <- list()

for (round.number in 1:3) {
  result <- analyzeRound(barTable[active, , drop=FALSE])
  result$round_number <- round.number
  result$active_cell_ids <- rownames(barTable)[active]
  rounds[[round.number]] <- result
  negative <- result$selected_calls == "Negative"
  final[active[negative]] <- "Negative"
  if (!any(negative)) {
    final[active] <- result$selected_calls
    break
  }
  survivors <- active[!negative]
  if (!length(survivors)) break
  if (round.number == 3L) {
    final[survivors] <- result$selected_calls[!negative]
    break
  }
  active <- survivors
}

payload <- list(
  KernSmooth=observed,
  calls=unname(final),
  rounds=rounds
)
write_json(payload, args[[2L]], auto_unbox=TRUE, digits=NA, pretty=TRUE)
"""


class DemultiplexTests(unittest.TestCase):
    def test_normalization_replaces_nonfinite_before_centering(self) -> None:
        counts = np.array([[0, 1], [2, 4], [8, 0]], dtype=np.int64)
        observed = demultiplex.normalize_hash_counts(counts)
        expected_before_centering = np.array([[0.0, 0.0], [1.0, 2.0], [3.0, 0.0]])
        expected = expected_before_centering - expected_before_centering.mean(
            axis=0, keepdims=True
        )
        np.testing.assert_allclose(observed, expected, rtol=0, atol=0)
        np.testing.assert_allclose(observed.mean(axis=0), 0.0, rtol=0, atol=1e-16)

    def test_local_maxima_matches_source_plateau_and_endpoint_rules(self) -> None:
        np.testing.assert_array_equal(
            demultiplex.local_maxima([0, 1, 1, 0, 2]),
            [1, 4],
        )
        np.testing.assert_array_equal(
            demultiplex.local_maxima([1, 1, 0]),
            [],
        )
        np.testing.assert_array_equal(
            demultiplex.local_maxima([0, 0, 1]),
            [2],
        )

    def test_type7_quantile_and_ge_classification(self) -> None:
        self.assertEqual(
            demultiplex.type7_quantile(np.array([10.0, 2.0]), 0.25),
            4.0,
        )
        calls = demultiplex.classify_normalized(
            np.array([[1.0, 0.0], [1.0, 2.0], [0.0, 0.0]]),
            ("h1", "h2"),
            np.array([1.0, 2.0]),
        )
        self.assertEqual(calls, ("h1", "Doublet", "Negative"))

    def test_q_selection_uses_first_best_local_maximum(self) -> None:
        selected, maxima = demultiplex.select_quantile(
            np.array([0.1, 0.4, 0.2, 0.4, 0.1]),
            np.array([0.01, 0.03, 0.05, 0.07, 0.09]),
        )
        np.testing.assert_array_equal(maxima, [1, 3])
        self.assertEqual(selected, 1)

    def test_undefined_hash_threshold_stops_without_fallback(self) -> None:
        with self.assertRaisesRegex(demultiplex.DemultiplexError, "bandwidth"):
            demultiplex.classify_round(
                np.ones((20, 2), dtype=np.int64),
                tuple(f"cell-{index}" for index in range(20)),
                ("h1", "h2"),
            )

    def test_three_round_peeling_removes_only_negatives_and_reclassifies(self) -> None:
        counts = np.ones((3, 2), dtype=np.int64)
        cell_ids = ("a", "b", "c")
        seen: list[tuple[int, tuple[str, ...]]] = []

        def fake_round(
            raw_counts: np.ndarray,
            active_cell_ids: tuple[str, ...],
            hash_labels: tuple[str, ...],
            *,
            round_number: int,
        ) -> demultiplex.RoundResult:
            del raw_counts, hash_labels
            seen.append((round_number, active_cell_ids))
            calls = {
                1: ("Negative", "Doublet", "hash-1"),
                2: ("Negative", "hash-1"),
                3: ("hash-2",),
            }[round_number]
            return demultiplex.RoundResult(
                round_number=round_number,
                active_cell_ids=active_cell_ids,
                normalized_counts=np.empty((len(calls), 2)),
                hash_models=(),
                thresholds=np.empty((0, 2)),
                calls_by_q=np.empty((0, len(calls)), dtype=object),
                singlet_fraction=np.empty(0),
                q_local_maxima=np.empty(0, dtype=np.int64),
                selected_q_index=0,
                selected_q=0.01,
                selected_calls=calls,
            )

        with mock.patch.object(demultiplex, "classify_round", side_effect=fake_round):
            result = demultiplex.demultiplex(
                counts,
                cell_ids,
                ("hash-1", "hash-2"),
            )

        self.assertEqual(
            seen,
            [(1, ("a", "b", "c")), (2, ("b", "c")), (3, ("c",))],
        )
        self.assertEqual(result.calls, ("Negative", "Negative", "hash-2"))

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_full_core_matches_pinned_r_synthetic_oracle(self) -> None:
        counts, cell_ids, hash_labels = synthetic_hash_counts()
        with tempfile.TemporaryDirectory(prefix="vcc107-demultiplex-oracle-") as raw:
            directory = Path(raw)
            counts_path = directory / "counts.csv"
            script_path = directory / "oracle.R"
            output_path = directory / "oracle.json"
            header = "," + ",".join(hash_labels)
            np.savetxt(
                counts_path,
                counts,
                fmt="%d",
                delimiter=",",
                header=header,
                comments="",
            )
            lines = counts_path.read_text(encoding="utf-8").splitlines()
            lines[1:] = [
                f"{cell_id},{line}" for cell_id, line in zip(cell_ids, lines[1:])
            ]
            counts_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            script_path.write_text(R_ORACLE, encoding="utf-8")
            subprocess.run(
                ["Rscript", str(script_path), str(counts_path), str(output_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            oracle = json.loads(output_path.read_text(encoding="utf-8"))

        self.assertEqual(oracle["KernSmooth"], EXPECTED_KERNSMOOTH_VERSION)
        observed = demultiplex.demultiplex(counts, cell_ids, hash_labels)
        self.assertEqual(observed.calls, tuple(oracle["calls"]))
        self.assertEqual(len(observed.rounds), len(oracle["rounds"]))

        for actual, expected in zip(
            observed.rounds,
            oracle["rounds"],
            strict=True,
        ):
            self.assertEqual(actual.round_number, expected["round_number"])
            self.assertEqual(actual.active_cell_ids, tuple(expected["active_cell_ids"]))
            np.testing.assert_allclose(
                actual.normalized_counts,
                np.asarray(expected["normalized"]),
                rtol=0,
                atol=1e-14,
            )
            np.testing.assert_allclose(
                demultiplex.Q_VALUES,
                expected["q_values"],
                rtol=0,
                atol=0,
            )
            np.testing.assert_allclose(
                actual.thresholds,
                np.asarray(expected["thresholds"]),
                rtol=0,
                atol=2e-13,
            )
            np.testing.assert_array_equal(
                actual.calls_by_q,
                np.asarray(expected["calls_by_q"], dtype=object),
            )
            np.testing.assert_allclose(
                actual.singlet_fraction,
                expected["singlet_fraction"],
                rtol=0,
                atol=0,
            )
            np.testing.assert_array_equal(
                actual.q_local_maxima,
                expected["q_local_maxima"],
            )
            self.assertEqual(actual.selected_q_index, expected["selected_q_index"])
            self.assertEqual(actual.selected_q, expected["selected_q"])
            self.assertEqual(actual.selected_calls, tuple(expected["selected_calls"]))

            for model, expected_model in zip(
                actual.hash_models,
                expected["models"],
                strict=True,
            ):
                np.testing.assert_allclose(
                    model.evaluation_points,
                    expected_model["evaluation_points"],
                    rtol=0,
                    atol=1e-14,
                )
                np.testing.assert_allclose(
                    model.density,
                    expected_model["density"],
                    rtol=0,
                    atol=2e-13,
                )
                np.testing.assert_array_equal(model.extrema, expected_model["extrema"])
                self.assertEqual(model.low_extreme, expected_model["low_extreme"])
                self.assertEqual(model.high_extreme, expected_model["high_extreme"])


if __name__ == "__main__":
    unittest.main()

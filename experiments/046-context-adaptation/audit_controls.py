"""Negative controls for experiment 046's context-adaptation result.

The main result is only context adaptation if paired target labels matter.
This audit compares the fixed rank-32 PCA-ridge map against:

* the same fit after permuting target perturbation rows; and
* maps trained for each of the two wrong target contexts.

It uses no new hyperparameter search. The 259-panel outer scores are diagnostic
checks of the already-selected method, not another selection surface.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUN_PATH = HERE / "run.py"
SPEC = importlib.util.spec_from_file_location("context_adaptation_run", RUN_PATH)
RUN = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def main() -> None:
    genes, train_perts, outer_perts, train_mask = RUN.validate_cache()
    local_feat = np.flatnonzero(train_mask).astype(np.int64)
    all_feat = np.arange(len(genes), dtype=np.int64)
    banks = RUN.Banks()
    candidate = RUN.Candidate(
        "pca_ridge",
        ridge=0.3,
        rank=32,
        eta=1.0,
    )

    train_pairs = {}
    outer_pairs = {}
    correct_models = {}
    permuted_models = {}
    for target in RUN.OUTER_DONORS:
        sources = tuple(name for name in RUN.OUTER_DONORS if name != target)
        train_pairs[target] = RUN.make_train_pair(
            banks, target, sources, "raw", local_feat
        )
        outer_pairs[target] = RUN.make_outer_pair(
            banks, target, sources, "raw", local_feat, all_feat
        )
        seed = RUN.stable_seed(("audit", target))
        correct_models[target] = RUN.fit_models(
            train_pairs[target].local_query,
            train_pairs[target].local_truth,
            [candidate],
            seed,
        )[candidate.key]
        rng = np.random.default_rng(seed + 1)
        permuted_models[target] = RUN.fit_models(
            train_pairs[target].local_query,
            train_pairs[target].local_truth[rng.permutation(len(train_perts))],
            [candidate],
            seed,
        )[candidate.key]

    rows = {}
    for target in RUN.OUTER_DONORS:
        pair = outer_pairs[target]
        baseline = RUN.pds_from_scores(pair.full_query @ pair.full_truth.T)
        correct = RUN.pds_from_scores(
            RUN.model_scores(pair, correct_models[target])
        )
        permuted = RUN.pds_from_scores(
            RUN.model_scores(pair, permuted_models[target])
        )
        wrong = {}
        for source_target in RUN.OUTER_DONORS:
            if source_target == target:
                continue
            score = RUN.pds_from_scores(
                RUN.model_scores(pair, correct_models[source_target])
            )
            wrong[source_target] = {
                "score": score,
                "delta": score - baseline,
            }
        rows[target] = {
            "baseline": baseline,
            "paired": correct,
            "paired_delta": correct - baseline,
            "permuted": permuted,
            "permuted_delta": permuted - baseline,
            "wrong_target_maps": wrong,
        }
        print(
            f"{target:16} baseline={baseline:.5f} "
            f"paired={correct:.5f} ({correct-baseline:+.5f}) "
            f"permuted={permuted:.5f} ({permuted-baseline:+.5f}) "
            f"wrong_best={max(v['delta'] for v in wrong.values()):+.5f}",
            flush=True,
        )

    paired_beats_permuted = {
        target: rows[target]["paired"] > rows[target]["permuted"]
        for target in RUN.OUTER_DONORS
    }
    paired_beats_wrong = {
        target: all(
            rows[target]["paired"] > value["score"]
            for value in rows[target]["wrong_target_maps"].values()
        )
        for target in RUN.OUTER_DONORS
    }
    result = {
        "candidate": candidate.key,
        "training_perturbations": len(train_perts),
        "outer_vcc_perturbations": len(outer_perts),
        "rows": rows,
        "macro": {
            "baseline": float(np.mean([
                rows[target]["baseline"] for target in RUN.OUTER_DONORS
            ])),
            "paired_delta": float(np.mean([
                rows[target]["paired_delta"] for target in RUN.OUTER_DONORS
            ])),
            "permuted_delta": float(np.mean([
                rows[target]["permuted_delta"] for target in RUN.OUTER_DONORS
            ])),
            "wrong_target_best_delta": float(np.mean([
                max(
                    value["delta"]
                    for value in rows[target]["wrong_target_maps"].values()
                )
                for target in RUN.OUTER_DONORS
            ])),
        },
        "context_specificity_gate": {
            "criterion": (
                "The correctly paired target map must beat both the permuted-label "
                "map and every wrong-target map in each outer donor."
            ),
            "paired_beats_permuted": paired_beats_permuted,
            "paired_beats_wrong_target_maps": paired_beats_wrong,
            "passed": (
                all(paired_beats_permuted.values())
                and all(paired_beats_wrong.values())
            ),
        },
        "peak_rss_gib": RUN.rss_gib(),
    }
    (HERE / "results_controls.json").write_text(json.dumps(result, indent=2))
    print(
        f"wrote results_controls.json; peak RSS={RUN.rss_gib():.2f} GiB",
        flush=True,
    )


if __name__ == "__main__":
    main()

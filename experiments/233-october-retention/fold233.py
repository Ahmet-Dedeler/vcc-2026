"""233 — PLOLO fold construction, made panel-parameterised and component-decomposable.

`experiments/162-retrodicting-validator/plolo.py:build_fold` builds one object -- the
CONSENSUS `pre` -- and throws the per-source pieces away.  Phase 2 needs them: the attribution
question is *which* component's contribution shrinks on a disjoint panel, and that cannot be
asked of a pre-averaged matrix.  So this is `build_fold` re-expressed, keeping:

    per-source effects + availability, the PerturbMatch blend as a separable step,
    the truth, the control profile, the CPM gate, the scored columns, the operator axis

and it is checked against the original: on today's panel my consensus must equal
`plolo.build_fold`'s `pre` EXACTLY (`KA-fold`).  Everything else in this experiment is a ratio
against that object, so if the two disagree the retention numbers are measuring my
re-implementation rather than the panel.

The operator is loaded through `panels233.load_frozen_operator_233`, which verifies the AXIS
HASH and not merely `n_cols` -- `n_cols` is 17,942 on ALL THREE panels (measured, `--plan`),
which is exactly why the production check cannot see the Oct-22 trap.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import panels233 as P  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, E064, OFFICIAL  # noqa: E402

PM_BLEND = 0.5
DONOR_PATHS = plolo.DONOR_PATHS
PM_PATH = ROOT / "data/external/perturbmatch"


def _prepared(name, effects, avail, mean_on_axis, raw_median_norm, mapped, dups):
    return E062.PreparedSource(name=name, effects=effects,
                               target_available=np.asarray(avail, bool),
                               mean_on_axis=np.asarray(mean_on_axis, float),
                               raw_median_norm=float(raw_median_norm),
                               mapped_gene_count=int(mapped),
                               duplicate_symbols=tuple(dups))


def build(panel: str, held: str, limit_gib: float = 10.0) -> dict:
    """Build (or load) the fold for `panel`, holding out donor `held`."""
    cache = P.WORK / "folds" / f"{panel}__{held}.npz"
    cache.parent.mkdir(parents=True, exist_ok=True)
    sources = P.FOLD_SOURCES[held]
    panel_t = P.read_panel(panel)
    panel_set = set(panel_t)

    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in sources}
    held_bank = E062.Bank.load(held, DONOR_PATHS[held])
    targets = tuple(g for g in panel_t if g in held_bank.pert_pos)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in panel_t if g in official_pos],
                                 dtype=np.int64)

    # The operator AXIS is a property of (donor coverage, panel) alone, so it is computed
    # here; the operator FILE is verified against it in `verify_op` before it is applied.
    op_p = P.op_path(panel, held)
    donor_union: set[str] = set()
    for n in sources:
        donor_union |= set(banks[n].genes)
    wcols = np.asarray([i for i, g in enumerate(OFFICIAL)
                        if g in donor_union and g not in panel_set], dtype=np.int64)

    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if tuple(str(t) for t in z["targets"]) == targets:
            out = {k: z[k] for k in z.files if k != "audit"}
            out["audit"] = json.loads(str(z["audit"]))
            out.update(panel=panel, held=held, targets=targets, sources=sources,
                       panel_positions=panel_positions, wcols=wcols, op_path=op_p)
            print(f"[{panel}/{held}] fold cache hit ({len(targets)} targets)", flush=True)
            return out

    prep = {n: E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                        chunk_rows=256, max_rss_gib=limit_gib)
            for n, b in banks.items()}
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prep = E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions,
                                       chunk_rows=256, max_rss_gib=limit_gib)
    shared = (E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, panel_set)
              if "k562_gwps" in banks else np.zeros(0, dtype=np.int64))

    held_prep = E062.prepare_main_source(held_bank, targets, OFFICIAL, panel_positions,
                                         chunk_rows=256, max_rss_gib=limit_gib)
    truth = E062.raw_truth(held_prep)
    assert held_prep.target_available.all(), "a held target has no truth row"

    b_ctrl, held_cols = plolo._b_ctrl_on_axis(held_bank)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > plolo.MIN_CPM) & (b_ctrl > 0)
    scored_mask = np.zeros(len(OFFICIAL), dtype=bool)
    scored_mask[held_cols] = True
    scored_mask[panel_positions] = False
    scored = np.flatnonzero(scored_mask)

    pm_both = np.zeros(len(targets), dtype=bool)
    if "k562_gwps" in banks:
        pm_both = prep["k562_gwps"].target_available & pm_prep.target_available

    audit = dict(
        panel=panel, held=held, sources=list(sources), n_targets=len(targets),
        panel_coverage=f"{len(targets)}/300",
        source_coverage={n: int(prep[n].target_available.sum()) for n in sources},
        carriers={str(k): int(v) for k, v in zip(*np.unique(
            sum(prep[n].target_available.astype(int) for n in sources),
            return_counts=True))},
        pm_supported_targets=int(pm_both.sum()),
        op=str(op_p.relative_to(ROOT)),
        op_sha256=(E062.sha256_file(op_p)[:12] if op_p.exists() else "not fitted yet"),
        op_cols=int(len(wcols)), scored_genes=int(len(scored)),
        gated_genes=int(gate.sum()))

    payload = dict(
        truth=truth.astype(np.float64), b_ctrl=b_ctrl, gate=gate, scored=scored,
        pm_effects=np.asarray(pm_prep.effects, dtype=np.float64),
        pm_avail=np.asarray(pm_prep.target_available, bool), pm_shared=shared,
        pm_both=pm_both, targets=np.asarray(targets), audit=json.dumps(audit))
    for n in sources:
        payload[f"eff__{n}"] = np.asarray(prep[n].effects, dtype=np.float64)
        payload[f"avail__{n}"] = np.asarray(prep[n].target_available, bool)
        payload[f"mean__{n}"] = np.asarray(prep[n].mean_on_axis, dtype=np.float64)
        payload[f"rmn__{n}"] = np.float64(prep[n].raw_median_norm)
    np.savez_compressed(cache, **payload)
    print(f"[{panel}/{held}] built: {json.dumps(audit)}", flush=True)
    out = {k: v for k, v in payload.items() if k != "audit"}
    out["audit"] = audit
    out.update(panel=panel, held=held, targets=targets, sources=sources,
               panel_positions=panel_positions, wcols=wcols, op_path=op_p)
    return out


def verify_op(panel: str, held: str) -> dict:
    """Refuse an operator whose axis was fitted for a different panel (AGENTS.md §8.3)."""
    banks = [E062.Bank.load(n, DONOR_PATHS[n]) for n in P.FOLD_SOURCES[held]]
    op_p = P.op_path(panel, held)
    op, wcols = P.load_frozen_operator_233(op_p, banks, OFFICIAL, set(P.read_panel(panel)),
                                           expected_sha256=E062.sha256_file(op_p))
    return dict(path=str(op_p.relative_to(ROOT)), n_cols=int(len(wcols)),
                axis_sha256=json.loads(op_p.with_suffix(".json").read_text())["axis_sha256"])


def consensus(fold: dict, *, pm: bool = True, only: str | None = None) -> np.ndarray:
    """`E062.build_baseline` over the fold's sources, with the PerturbMatch blend optional.

    `only` restricts the consensus to a single source -- the ablation that prices the SECOND
    donor, i.e. how much of the score is bought by having more than one carrier.
    """
    names = [only] if only else list(fold["sources"])
    prepared = []
    for n in names:
        eff = np.asarray(fold[f"eff__{n}"], dtype=np.float64).copy()
        if n == "k562_gwps" and pm:
            rows = np.flatnonzero(np.asarray(fold["pm_both"], bool))
            sh = np.asarray(fold["pm_shared"], dtype=np.int64)
            if len(rows) and len(sh):
                pm_eff = np.asarray(fold["pm_effects"], dtype=np.float64)
                eff[np.ix_(rows, sh)] = ((1 - PM_BLEND) * eff[np.ix_(rows, sh)]
                                         + PM_BLEND * pm_eff[np.ix_(rows, sh)])
            name = "k562_gwps_pm50"
        else:
            name = n
        prepared.append(_prepared(name, eff, fold[f"avail__{n}"], fold[f"mean__{n}"],
                                  float(fold[f"rmn__{n}"]), 0, ()))
    pre, votes = E062.build_baseline(prepared)
    return pre, votes

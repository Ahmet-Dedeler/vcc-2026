"""162 — PLOLO: panel-restricted leave-one-library-out, scored after the renderer.

WHAT THIS IS
------------
The design a competitor (GISL) describes for a validator that retrodicts board deltas:
hold each donor library out as the TRUTH, over **the real 300 challenge panel targets**
rather than a harness's 85/132, rebuild the production effect matrix from the remaining
libraries, and score.  Three structural differences from every incumbent harness here:

  1. the decoy set is ~290, matching the challenge's 300, not 85 or 132.  `pds` is a
     retrieval metric and its value depends on how many decoys there are;
  2. the targets are the actual panel -- a deliberately non-common-essential gene set --
     so the response distribution is the shipped one;
  3. the truth is a deep pseudobulk library, not 150 harness cells, so the truth-side
     sampling noise that depresses every absolute on the Orion surfaces is much smaller.

MANDATORY: the prediction is scored AFTER `profile_from_effect` + the scorer's read-back
(experiment 156).  `identity/intended` and `identity/delivered` differ by +0.0032/+0.0038
raw pds, so an effect-level score measures an object the board never sees.

LEAKAGE
-------
Every fold uses a whitening operator fitted WITHOUT the held library:
  held k562_gwps    -> 062 held_k562_gwps      sources orion_hct116 + orion_hek293t
  held orion_hct116 -> 062 held_orion_hct116   sources k562_gwps + orion_hek293t
  held orion_hek293t-> 062 held_orion_hek293t  sources k562_gwps + orion_hct116
  held kolf2        -> PRODUCTION union_pruned sources all three
`kolf2` is not a production donor and contributed no row to any production operator, so
that fold is leak-free while using the exact shipped operator and the exact shipped
three-source mixture.  It is the closest available analogue of the shipped build, and it
is the only fold on which `a2-prune` (a union vs union_pruned comparison) is expressible.

CD4 is a production INPUT on every fold and is never a truth, so no fold can leak it.
`h1_2025` is excluded as a truth: it covers 20/300 panel targets.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))

TS = 50_000.0
EFFECT_SCALE = 6.817           # 062/064/production
MIN_CPM = 10.0                 # a15 shipped
PM_BLEND = 0.5
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
PM_PATH = ROOT / "data/external/perturbmatch"

PANEL = tuple(
    line.split(",")[0]
    for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
    if line.strip()
)
OFFICIAL = tuple(
    line.strip()
    for line in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
    if line.strip()
)
PANEL_SET = set(PANEL)
assert len(PANEL) == 300 and len(OFFICIAL) == 18533


def _load(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E062 = _load("e162_062", ROOT / "experiments/062-cd4-post-whiten/run.py")
E064 = _load("e162_064", ROOT / "experiments/064-pm-cd4-composition/run.py")

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    "kolf2": ROOT / "data/external/kolf2",
}
HELD_OPS = ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
PROD_OPS = ROOT / "data/donors"

FOLDS = {
    "k562_gwps": dict(sources=("orion_hct116", "orion_hek293t"),
                      op=HELD_OPS / "held_k562_gwps.npz", op_kind="held"),
    "orion_hct116": dict(sources=("k562_gwps", "orion_hek293t"),
                         op=HELD_OPS / "held_orion_hct116.npz", op_kind="held"),
    "orion_hek293t": dict(sources=("k562_gwps", "orion_hct116"),
                          op=HELD_OPS / "held_orion_hek293t.npz", op_kind="held"),
    "kolf2": dict(sources=("k562_gwps", "orion_hct116", "orion_hek293t"),
                  op=PROD_OPS / "whitener_union_pruned.npz", op_kind="production"),
}


# --------------------------------------------------------------------- renderer
def profile_from_effect(effect, b_ctrl, gate, mass_fix=None):
    """Production `profile_from_effect` for one row, gate precomputed."""
    gated = np.where(gate, effect, 0.0)
    if mass_fix is not None:
        gated = mass_consistent_effect(gated, b_ctrl, gate, mode=mass_fix)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS
    return profile / float(profile.sum())


def mass_consistent_effect(gated, b, on, mode="minnorm", iters=64, rtol=1e-12):
    """Verbatim port of production `mass_consistent_effect` (156)."""
    h = np.asarray(gated, dtype=np.float64).copy()
    for _ in range(iters):
        live = on & (b + h > 0.0)
        residual = float(np.clip(np.expm1(b + h), 0.0, None).sum() - TS)
        if abs(residual) <= rtol * TS:
            break
        w = np.where(live, np.exp(b + h), 0.0)
        denom = float(w @ w) if mode == "minnorm" else float(w.sum())
        if not np.isfinite(denom) or denom <= 0.0:
            break
        step = -residual / denom
        h = h + (step * w if mode == "minnorm" else np.where(live, step, 0.0))
    return h


def deliver(effects, b_ctrl, gate, mass_fix=None):
    """Render every row and read it back the way the scorer does."""
    out = np.empty_like(effects)
    for i in range(effects.shape[0]):
        p = profile_from_effect(effects[i], b_ctrl, gate, mass_fix)
        out[i] = np.log1p(TS * p) - b_ctrl
    return out


def render_diagnostics(effects, b_ctrl, gate):
    """Mass ratio and render residual, the two numbers 156 is about.

    Reported so a fold can be compared against the measured challenge-frame values
    (mass ratio 1.00134, residual 2.83-3.35% of row norm on A/B/C) before any of its
    geometry deltas are believed.
    """
    G = np.where(gate[None, :], effects, 0.0)
    lin = np.clip(np.expm1(b_ctrl[None, :] + G), 0.0, None)
    ratio = lin.sum(1) / TS
    delivered = np.log1p(TS * lin / lin.sum(1, keepdims=True)) - b_ctrl[None, :]
    resid = np.median(np.linalg.norm(delivered - G, axis=1) / np.linalg.norm(G, axis=1))
    return dict(mass_ratio=float(ratio.mean()), mass_ratio_sd=float(ratio.std()),
                render_residual_pct=float(100.0 * resid))


def challenge_bctrl(context: str) -> np.ndarray:
    """The builder's own control profile for a challenge context, on the official axis."""
    import h5py
    cache = WORK / f"bctrl-challenge-{context}.npy"
    if cache.exists():
        return np.load(cache)
    path = ROOT / f"data/challenge/context_{context}.h5ad"
    with h5py.File(path, "r") as h:
        idx = h["X/indices"][:]
        dat = h["X/data"][:].astype(np.float64)
        node = h["var/_index"]
        raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
    assert tuple(genes) == OFFICIAL, "challenge gene axis differs from gene_names.csv"
    totals = np.bincount(idx, weights=dat, minlength=len(OFFICIAL))
    b = np.log1p(TS * totals / totals.sum())
    np.save(cache, b)
    return b


# -------------------------------------------------------------------- endpoints
def _unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def endpoints(P, T):
    """pds (exact rank statistic, ties at half), macro cosine, hubness, sse*, sign credit."""
    A, B = _unit(P), _unit(T)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    pds = float(np.mean((beaten + 0.5 * ties) / (n - 1)))
    G = A @ A.T
    rr = float(G[~np.eye(n, dtype=bool)].mean())
    a = float((P * P).sum()); x = float((P * T).sum()); d = float((T * T).sum())
    sse = (1.0 - x * x / (a * d)) if a > 0 else 1.0
    live = T != 0
    agree = (np.sign(P) == np.sign(T)) & live
    sgn = float(agree.sum() / max(live.sum(), 1))
    wt = np.abs(T) * live
    sgn_w = float((np.where(agree, wt, 0.0)).sum() / max(wt.sum(), 1e-300))
    return dict(pds=pds, cos=float(own.mean()), rr=rr, sse=sse,
                gstar=(x / a) if a > 0 else 0.0, sgn=sgn, sgn_w=sgn_w)


# ------------------------------------------------------------------ fold build
def _b_ctrl_on_axis(bank):
    """Map a donor's control profile onto the official axis and re-impose sum expm1 == TS.

    The donor's own axis is a superset/subset of the official one, so the mapped linear
    profile no longer sums to TS.  `mass_consistent_effect` and the renderer's own
    normalisation both assume it does (156 §1), so the profile is renormalised here and
    the identity is asserted.
    """
    amap = E062.map_bank_to_axis(bank, OFFICIAL)
    raw = np.load(bank.b_ctrl_path, mmap_mode="r")
    b = np.zeros(len(OFFICIAL), dtype=np.float64)
    b[amap.unique_target_positions] = np.asarray(
        raw[amap.unique_source_columns], dtype=np.float64)
    lin = np.clip(np.expm1(b), 0.0, None)
    lin = lin / lin.sum() * TS
    b = np.log1p(lin)
    assert abs(np.expm1(b).sum() - TS) < 1e-6 * TS
    return b, amap.unique_target_positions


def build_fold(held: str, limit_gib: float = 40.0):
    """Everything that does not depend on an arm's knobs."""
    cache = WORK / f"fold-{held}.npz"
    spec = FOLDS[held]
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    held_bank = E062.Bank.load(held, DONOR_PATHS[held])

    # targets: the panel genes this library can adjudicate
    targets = tuple(g for g in PANEL if g in held_bank.pert_pos)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)

    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if tuple(str(t) for t in z["targets"]) == targets:
            print(f"[{held}] cache hit", flush=True)
            return dict(held=held, targets=targets, pre=z["pre"], truth=z["truth"],
                        b_ctrl=z["b_ctrl"], gate=z["gate"], scored=z["scored"],
                        wcols=z["wcols"], op_path=spec["op"], sources=spec["sources"],
                        panel_positions=panel_positions, audit=json.loads(str(z["audit"])))

    prepared = {n: E062.prepare_main_source(b, targets, OFFICIAL, panel_positions,
                                            chunk_rows=256, max_rss_gib=limit_gib)
                for n, b in banks.items()}
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prepared = E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions,
                                           chunk_rows=256, max_rss_gib=limit_gib)
    shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET) \
        if "k562_gwps" in banks else None

    source_prepared = []
    pm_n = 0
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            pm_n = int(both.sum())
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            source_prepared.append(E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(prepared[n].target_available, bool).copy(),
                mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
                raw_median_norm=float(prepared[n].raw_median_norm),
                mapped_gene_count=int(prepared[n].mapped_gene_count),
                duplicate_symbols=tuple(prepared[n].duplicate_symbols)))
        else:
            source_prepared.append(prepared[n])
    pre, votes = E062.build_baseline(source_prepared)

    held_prep = E062.prepare_main_source(held_bank, targets, OFFICIAL, panel_positions,
                                         chunk_rows=256, max_rss_gib=limit_gib)
    truth = E062.raw_truth(held_prep)
    assert held_prep.target_available.all(), "a held target has no truth row"

    b_ctrl, held_cols = _b_ctrl_on_axis(held_bank)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > MIN_CPM) & (b_ctrl > 0)

    # scored columns: what the held library can actually adjudicate, minus the panel.
    # `pds_cosine` drops all 300 panel coordinates; both sides are already zero there.
    scored_mask = np.zeros(len(OFFICIAL), dtype=bool)
    scored_mask[held_cols] = True
    scored_mask[panel_positions] = False
    scored = np.flatnonzero(scored_mask)

    op_path = Path(spec["op"])
    if spec["op_kind"] == "production":
        wcols = E064.operator_columns([banks[n] for n in spec["sources"]],
                                      OFFICIAL, PANEL_SET)
    else:
        _, wcols = E064.load_frozen_operator(
            op_path, [banks[n] for n in spec["sources"]], OFFICIAL, PANEL_SET,
            expected_sha256=E062.sha256_file(op_path))

    audit = dict(
        held=held, sources=list(spec["sources"]), n_targets=len(targets),
        panel_coverage=f"{len(targets)}/300",
        op=str(op_path.relative_to(ROOT)), op_sha256=E062.sha256_file(op_path)[:8],
        op_cols=int(len(wcols)), scored_genes=int(len(scored)),
        gated_genes=int(gate.sum()), pm_supported_targets=pm_n,
        vote_histogram={str(v): int((votes == v).sum()) for v in np.unique(votes)},
        truth_source_overlap={
            n: int(sum(t in banks[n].pert_pos for t in targets)) for n in spec["sources"]},
    )
    np.savez_compressed(cache, pre=pre.astype(np.float32), truth=truth.astype(np.float32),
                        b_ctrl=b_ctrl, gate=gate, scored=scored, wcols=wcols,
                        targets=np.asarray(targets), audit=json.dumps(audit))
    print(f"[{held}] built: {json.dumps(audit)}", flush=True)
    return dict(held=held, targets=targets, pre=pre, truth=truth, b_ctrl=b_ctrl,
                gate=gate, scored=scored, wcols=wcols, op_path=op_path,
                sources=spec["sources"], panel_positions=panel_positions, audit=audit)


class Operator:
    """(C + lam I)^(-gamma/2) with production's candidate-specific on-axis preserve."""

    def __init__(self, path: Path):
        with np.load(path) as z:
            self.evals = np.asarray(z["evals"], dtype=np.float64)
            self.evecs = np.asarray(z["evecs"], dtype=np.float64)

    def apply(self, M, wcols, gamma, lam):
        out = np.asarray(M, dtype=np.float64).copy()
        before = float(np.median(np.linalg.norm(out[:, wcols], axis=1))) or 1.0
        w = (self.evals + lam) ** (-gamma / 2.0)
        t = ((out[:, wcols] @ self.evecs) * w) @ self.evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
        out[:, wcols] = t * (before / after)
        return out


def cd4_for(targets):
    banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, OFFICIAL, PANEL_SET, "gene_norm",
                                 chunk_rows=256, max_rss_gib=40.0)
    complete = np.asarray(prep.target_condition_count) == len(CD4_NAMES)
    mean_effect = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    mean_effect[~complete] = 0.0
    votes[~complete] = 0
    return E062.CD4Prepared(mode=prep.mode, mean_effect=mean_effect, gene_votes=votes,
                            target_condition_count=np.asarray(
                                prep.target_condition_count, dtype=np.uint8).copy(),
                            condition_details=prep.condition_details)


def solve_fraction(whitened, cd4, target_f):
    """Production `solve_cd4_fraction`: f is a ratio of median row norms, linear in fraction."""
    added = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0)
    rn = float(np.median(np.linalg.norm(added, axis=1)))
    bn = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    per_unit = rn / bn if bn else 0.0
    return (target_f / per_unit if per_unit > 0 else 0.0), per_unit

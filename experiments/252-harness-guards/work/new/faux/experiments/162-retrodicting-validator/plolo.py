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

GUARDS (experiment 252, closing 242's findings #4, #5, #6, #13, #14 and the silent zero)
----------------------------------------------------------------------------------------
Every arithmetic path an existing experiment used is bit-identical; what changed is that
the harness now REFUSES inputs it used to accept silently, and names the escape hatches:

  * `Operator(path, name=None)` records the file's sha256 and resolves it against
    `OPERATOR_REGISTRY` (five operators, pinned by sha256 + metadata sha256 + axis digest +
    max eigenvalue).  A file that is neither registered nor self-described by its own
    metadata json (`sha256` + `axis_sha256`/`wcols_sha256`, the 233/235 convention) raises
    `OperatorGuardError`; `name=` asserts a specific operator, so a swapped file raises.
    `apply()` then asserts the axis it is applied on (ordered symbols, 235's `axis_sha256`;
    or the fit-time index digest) matches the operator's record -- the Oct-22 trap fires
    here -- and raises instead of substituting 1.0 when the on-axis median row norm is 0
    (242 #13; the message carries the coverage).  `apply(..., check_axis=False)` is the
    named bypass for a DELIBERATE trap measurement (233's `trapprice233.py`).
  * `build_fold` keys its cache on a digest of every input (`fold_inputs`: held, sources,
    operator path/kind/metadata sha, PM_BLEND, MIN_CPM, donor paths, panel and official-axis
    digests, targets).  A pre-252 cache carries no key and is accepted ONLY after its
    targets, wcols, gate, scored columns, sources and operator are recomputed from the live
    inputs and found identical (`fold_cache_status`); the returned audit records which.
  * `endpoints()` returns NaN for `sgn`/`sgn_w` (and `sse`/`gstar` where undefined) on a
    truth block with no live coordinate instead of 0.0 / ZeroDivisionError;
    `mean_endpoints()` is the aggregation that skips those explicitly and counts them.
  * `solve_fraction(..., sizing=None)`: the default (`SIZING_DEFAULT = "all_rows"`,
    production's A/B/C rule) RAISES `SilentZeroError` on a residual whose all-rows median
    is zero -- a source below 50% coverage used to solve to fraction 0.0 and read as
    identically absent (223 §7.1).  `sizing="carrying"` is 218 §7.1's rule (the builder's
    `solve_carrying_row_fraction`, October's D/E/F sizing); `sizing="legacy_all_rows"`
    reproduces the silent zero by name, loudly.  A reproduction driver for an archive script
    that read the zero as a diagnostic may set `plolo.SIZING_DEFAULT = "legacy_all_rows"`.
"""
from __future__ import annotations

import argparse
import hashlib
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


# ------------------------------------------------------------------- guards (252)
class HarnessGuardError(ValueError):
    """A guard in this harness refused an input it used to accept silently."""


class OperatorGuardError(HarnessGuardError):
    """Operator identity / axis / on-axis-preserve guard."""


class SilentZeroError(HarnessGuardError):
    """`solve_fraction` would have returned fraction 0.0 for a residual that is present."""


def axis_sha256(symbols) -> str:
    """sha256 over an ORDERED symbol axis: newline-joined UTF-8, no trailing newline.

    Byte-identical to `scripts/build_union_whitener.py: axis_sha256` and the builder's
    `operator_axis_sha256` (experiment 235).  On the three 17,942-column folds the digest of
    `OFFICIAL[wcols]` is `edf6ba88…`, the value the builder pins for both production
    operators in `LEGACY_OPERATOR_AXIS_SHA256` -- an independent cross-implementation check.
    """
    return hashlib.sha256("\n".join(map(str, symbols)).encode("utf-8")).hexdigest()


def index_sha256(wcols) -> str:
    """sha256 over the int64 column indices (062's `wcols_sha256`; 233's `axis_sha256`)."""
    return hashlib.sha256(np.ascontiguousarray(wcols, dtype=np.int64).tobytes()).hexdigest()


# The five operators any PLOLO caller may load, pinned by NAME.  `sha256` is the npz, checked
# against `E064.KNOWN_OPERATOR_SHA256` at import; `metadata_sha256` its sidecar json;
# `axis_sha256` the ordered-symbol digest of the axis it was fitted for (`OFFICIAL[wcols]`);
# `wcols_sha256` the fit-time index digest 062 wrote into the json (None where the fitter did
# not write one); `max_eigenvalue` the human-readable discriminator rule 10 cites.
OPERATOR_REGISTRY: dict[str, dict] = {
    "held_k562_gwps": dict(
        path="experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz",
        sha256="13d988fcc3b129d25d254524db73e8f8253143fa493ee64d9654f0f61bc3ac55",
        metadata_sha256="9df1f55ab2672a327228e9bbfc602107897e6cec15e4289c3e92f33746274e83",
        # measured by 252 over OFFICIAL[wcols] of the fold whose index digest equals 062's
        # fit-time `wcols_sha256` below (the independent record)
        axis_sha256="6eaba848c81b77dd1d4b83c0d61b7a8689c1e39762a34da7dcfe26ea1b5d765f",
        wcols_sha256="a000353b3dbeb8a55bdf7ae3bcc60abd13bbc52f714270413c10f7b1c23e1a02",
        n_cols=17806, max_eigenvalue=8.425881805318935,
        fitted_for="held k562_gwps: sources orion_hct116 + orion_hek293t, today's panel (062)"),
    "held_orion_hct116": dict(
        path="experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz",
        sha256="7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc",
        metadata_sha256="c4e00ad1fa5fbba7b9b0f93e4dae91b76da73638ef9431316923bb6ca4b585bc",
        axis_sha256="edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f",
        wcols_sha256="b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d",
        n_cols=17942, max_eigenvalue=8.036339907795812,
        fitted_for="held orion_hct116: sources k562_gwps + orion_hek293t, today's panel (062)"),
    "held_orion_hek293t": dict(
        path="experiments/062-cd4-post-whiten/cloud-output/whiteners/held_orion_hek293t.npz",
        sha256="97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1",
        metadata_sha256="cb22406a7a529380c4348bf636f2532c29a43ab0fc4b00d6f6d98fb5e02360c6",
        axis_sha256="edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f",
        wcols_sha256="b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d",
        n_cols=17942, max_eigenvalue=3.6958415582273543,
        fitted_for="held orion_hek293t: sources k562_gwps + orion_hct116, today's panel (062)"),
    "union_pruned": dict(
        path="data/donors/whitener_union_pruned.npz",
        sha256="e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
        metadata_sha256="838f263c6a4e1bed069d5faeb3bc3bb2b9234c5a5a933219a94241037a801f4d",
        axis_sha256="edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f",
        wcols_sha256=None,
        n_cols=17942, max_eigenvalue=5.894570063345345,
        fitted_for="PRODUCTION A/B/C (a15 -> h1union): all three donors, essential-pruned basis"),
    "union": dict(
        path="data/donors/whitener_union.npz",
        sha256="5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3",
        metadata_sha256="a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6",
        axis_sha256="edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f",
        wcols_sha256=None,
        n_cols=17942, max_eigenvalue=6.256702039557637,
        fitted_for="v11 production: all three donors, full union basis (162's `a2-prune` control)"),
}
_REGISTRY_BY_SHA = {e["sha256"]: name for name, e in OPERATOR_REGISTRY.items()}
assert {n: e["sha256"] for n, e in OPERATOR_REGISTRY.items()} == E064.KNOWN_OPERATOR_SHA256, \
    "plolo.OPERATOR_REGISTRY and E064.KNOWN_OPERATOR_SHA256 disagree"


def registered_operator_name(path) -> str | None:
    """Registry name for an operator PATH (no hashing), or None."""
    p = Path(path).resolve()
    for name, e in OPERATOR_REGISTRY.items():
        if (ROOT / e["path"]).resolve() == p:
            return name
    return None


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
    """pds (exact rank statistic, ties at half), macro cosine, hubness, sse*, sign credit.

    252 (242 #14): a truth block with NO live coordinate cannot carry a sign accuracy.  It
    used to read `sgn = sgn_w = 0.0` (indistinguishable from "every sign wrong") when the
    prediction was also zero, and to raise ZeroDivisionError in `sse` otherwise.  It now
    returns NaN for `sgn`, `sgn_w`, and for `sse`/`gstar` where they are undefined; every
    value on a block with live truth is bit-identical to before.  Aggregate such rows with
    `mean_endpoints`, which skips NaN explicitly and counts what it skipped.
    """
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
    nan = float("nan")
    if a > 0:
        sse = (1.0 - x * x / (a * d)) if d > 0 else nan
    else:
        sse = 1.0
    live = T != 0
    n_live = int(live.sum())
    agree = (np.sign(P) == np.sign(T)) & live
    wt = np.abs(T) * live
    if n_live > 0:
        sgn = float(agree.sum() / max(live.sum(), 1))
        sgn_w = float((np.where(agree, wt, 0.0)).sum() / max(wt.sum(), 1e-300))
    else:
        sgn, sgn_w = nan, nan
    return dict(pds=pds, cos=float(own.mean()), rr=rr, sse=sse,
                gstar=(x / a) if a > 0 else 0.0, sgn=sgn, sgn_w=sgn_w)


def mean_endpoints(rows, keys=None) -> dict:
    """Mean of endpoint dicts with EXPLICIT NaN skipping (252, 242 #14).

    A block whose `sgn`/`sgn_w` is NaN (no live truth) contributes nothing to that key's
    mean; `_skipped[key]` counts how many rows were dropped and `_n` how many were offered,
    so an aggregate can never quietly average a "0.0 = every sign wrong" that was really "no
    truth to grade".  Underscore keys of the input rows are ignored.
    """
    rows = list(rows)
    if not rows:
        raise ValueError("mean_endpoints: no rows")
    keys = list(keys) if keys is not None else [k for k in rows[0] if not str(k).startswith("_")]
    out, skipped = {}, {}
    for k in keys:
        vals = np.asarray([float(r[k]) for r in rows], dtype=np.float64)
        ok = ~np.isnan(vals)
        skipped[k] = int((~ok).sum())
        out[k] = float(vals[ok].mean()) if ok.any() else float("nan")
    out["_skipped"] = skipped
    out["_n"] = len(rows)
    return out


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


def _relpath(p) -> str:
    p = Path(p).resolve()
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def fold_inputs(held: str) -> dict:
    """Every input the cached fold arrays depend on, as a json-able dict (252, 242 #6).

    The old cache test compared the target list alone, so a changed `FOLDS[held]["sources"]`,
    `["op"]`, `PM_BLEND` or `MIN_CPM` printed `cache hit` and restored the stale audit.
    Loading the banks here is metadata-only (perts/genes lists and a memmap handle).
    """
    spec = FOLDS[held]
    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in spec["sources"]}
    held_bank = E062.Bank.load(held, DONOR_PATHS[held])
    targets = tuple(g for g in PANEL if g in held_bank.pert_pos)
    op_path = Path(spec["op"])
    meta_path = op_path.with_suffix(".json")
    payload = dict(
        schema=1, held=held, sources=list(spec["sources"]),
        op=_relpath(op_path), op_kind=spec["op_kind"],
        op_metadata_sha256=E062.sha256_file(meta_path) if meta_path.exists() else None,
        op_registered_name=registered_operator_name(op_path),
        pm_blend=float(PM_BLEND), min_cpm=float(MIN_CPM), pm_path=_relpath(PM_PATH),
        donor_paths={n: _relpath(DONOR_PATHS[n]) for n in (held, *spec["sources"])},
        panel_sha256=axis_sha256(PANEL), official_sha256=axis_sha256(OFFICIAL),
        targets=list(targets),
    )
    return dict(payload=payload, banks=banks, held_bank=held_bank, targets=targets,
                op_path=op_path, spec=spec)


def fold_cache_key(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def fold_cache_status(held: str, inputs: dict | None = None) -> dict:
    """Decide whether `WORK/fold-{held}.npz` may be used for the LIVE inputs, and why.

    status  'hit'      a keyed cache whose key equals the live key
            'legacy'   a pre-252 cache (no key) whose targets, wcols, gate, scored columns,
                       sources and operator all recompute identically from the live inputs
                       (and PM_BLEND is the 0.5 every pre-252 cache was built under)
            'miss'     a cache exists but disagrees; `reason` names the first disagreement
            'absent'   no cache file
    """
    inputs = inputs or fold_inputs(held)
    cache = WORK / f"fold-{held}.npz"
    key = fold_cache_key(inputs["payload"])
    out = dict(held=held, file=_relpath(cache), key=key, status="absent", reason=None,
               validated=[])
    if not cache.exists():
        return out
    z = np.load(cache, allow_pickle=True)
    targets = inputs["targets"]
    if tuple(str(t) for t in z["targets"]) != targets:
        out.update(status="miss", reason="targets differ (held bank / PANEL)")
        return out
    if "cache_key" in z.files:
        if str(z["cache_key"]) == key:
            out["status"] = "hit"
            return out
        stored = json.loads(str(z["inputs"])) if "inputs" in z.files else {}
        live = inputs["payload"]
        diff = [k for k in sorted(set(stored) | set(live)) if stored.get(k) != live.get(k)]
        out.update(status="miss", reason=f"cache key differs; inputs that differ: {diff}")
        return out
    # ---- a pre-252 cache: no key.  Recompute everything cheap and require identity.
    spec, banks, held_bank = inputs["spec"], inputs["banks"], inputs["held_bank"]
    audit = json.loads(str(z["audit"]))
    checks = []
    if audit.get("sources") != list(spec["sources"]):
        out.update(status="miss", reason=f"legacy cache sources {audit.get('sources')} != live {list(spec['sources'])}")
        return out
    checks.append("audit.sources")
    if audit.get("op") != _relpath(inputs["op_path"]):
        out.update(status="miss", reason=f"legacy cache op {audit.get('op')} != live {_relpath(inputs['op_path'])}")
        return out
    checks.append("audit.op")
    name = inputs["payload"]["op_registered_name"]
    if name is not None:
        if audit.get("op_sha256") != OPERATOR_REGISTRY[name]["sha256"][:8]:
            out.update(status="miss", reason="legacy cache op_sha256 prefix != registry pin")
            return out
        checks.append("audit.op_sha256==registry")
    if float(PM_BLEND) != 0.5:
        out.update(status="miss", reason=f"legacy cache was built at PM_BLEND 0.5; live PM_BLEND {PM_BLEND}")
        return out
    checks.append("PM_BLEND==0.5")
    wcols = E064.operator_columns([banks[n] for n in spec["sources"]], OFFICIAL, PANEL_SET)
    if not np.array_equal(np.asarray(z["wcols"]), wcols):
        out.update(status="miss", reason="legacy cache wcols != recomputed from live sources/panel")
        return out
    checks.append("wcols")
    b_ctrl = np.asarray(z["b_ctrl"], dtype=np.float64)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > MIN_CPM) & (b_ctrl > 0)
    if not np.array_equal(np.asarray(z["gate"], dtype=bool), gate):
        out.update(status="miss", reason=f"legacy cache gate != recomputed at MIN_CPM {MIN_CPM}")
        return out
    checks.append("gate@MIN_CPM")
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    _, held_cols = _b_ctrl_on_axis(held_bank)
    scored_mask = np.zeros(len(OFFICIAL), dtype=bool)
    scored_mask[held_cols] = True
    scored_mask[panel_positions] = False
    if not np.array_equal(np.asarray(z["scored"]), np.flatnonzero(scored_mask)):
        out.update(status="miss", reason="legacy cache scored columns != recomputed from held bank/panel")
        return out
    checks.append("scored")
    out.update(status="legacy", validated=checks)
    return out


def build_fold(held: str, limit_gib: float = 40.0):
    """Everything that does not depend on an arm's knobs."""
    cache = WORK / f"fold-{held}.npz"
    inputs = fold_inputs(held)
    spec, banks, held_bank, targets = inputs["spec"], inputs["banks"], inputs["held_bank"], inputs["targets"]
    key = fold_cache_key(inputs["payload"])
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    op_path = inputs["op_path"]
    op_name = inputs["payload"]["op_registered_name"]

    status = fold_cache_status(held, inputs)
    if status["status"] in ("hit", "legacy"):
        z = np.load(cache, allow_pickle=True)
        audit = json.loads(str(z["audit"]))
        audit["cache"] = dict(status=status["status"], key=key, file=status["file"],
                              validated=status["validated"],
                              pm_blend=float(PM_BLEND), min_cpm=float(MIN_CPM))
        print(f"[{held}] cache hit ({status['status']}"
              f"{': ' + ','.join(status['validated']) if status['validated'] else ''})", flush=True)
        return dict(held=held, targets=targets, pre=z["pre"], truth=z["truth"],
                    b_ctrl=z["b_ctrl"], gate=z["gate"], scored=z["scored"],
                    wcols=z["wcols"], op_path=spec["op"], sources=spec["sources"],
                    panel_positions=panel_positions, audit=audit,
                    op_name=op_name, cache_key=key)
    if status["status"] == "miss":
        print(f"[{held}] cache MISS: {status['reason']} -- rebuilding", flush=True)

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

    # The operator's file hash is compared against a PINNED expectation (252, 242 #4): the
    # registry entry for this path, or the operator's own metadata `sha256` (233/235's
    # self-describing fits).  Never the hash of the file being checked.
    meta_path = op_path.with_suffix(".json")
    if op_name is not None:
        expected_sha = OPERATOR_REGISTRY[op_name]["sha256"]
    elif meta_path.exists() and isinstance(json.loads(meta_path.read_text()).get("sha256"), str):
        expected_sha = json.loads(meta_path.read_text())["sha256"]
    else:
        raise OperatorGuardError(
            f"{_relpath(op_path)}: not in plolo.OPERATOR_REGISTRY and its metadata json declares "
            "no 'sha256'; there is no pinned expectation to check the file against")
    if spec["op_kind"] == "production":
        wcols = E064.operator_columns([banks[n] for n in spec["sources"]],
                                      OFFICIAL, PANEL_SET)
        op_sha = E062.sha256_file(op_path)
        if op_sha != expected_sha:
            raise OperatorGuardError(f"{_relpath(op_path)}: sha256 {op_sha} != pinned {expected_sha}")
    else:
        _, wcols = E064.load_frozen_operator(
            op_path, [banks[n] for n in spec["sources"]], OFFICIAL, PANEL_SET,
            expected_sha256=expected_sha)
        op_sha = expected_sha
    if op_name is not None:
        axis = axis_sha256(OFFICIAL[i] for i in wcols)
        if axis != OPERATOR_REGISTRY[op_name]["axis_sha256"]:
            raise OperatorGuardError(
                f"{op_name}: this fold's axis {axis} != the operator's registered axis "
                f"{OPERATOR_REGISTRY[op_name]['axis_sha256']} (fitted for a different panel)")

    audit = dict(
        held=held, sources=list(spec["sources"]), n_targets=len(targets),
        panel_coverage=f"{len(targets)}/300",
        op=str(op_path.relative_to(ROOT)), op_sha256=op_sha[:8],
        op_cols=int(len(wcols)), scored_genes=int(len(scored)),
        gated_genes=int(gate.sum()), pm_supported_targets=pm_n,
        vote_histogram={str(v): int((votes == v).sum()) for v in np.unique(votes)},
        truth_source_overlap={
            n: int(sum(t in banks[n].pert_pos for t in targets)) for n in spec["sources"]},
    )
    np.savez_compressed(cache, pre=pre.astype(np.float32), truth=truth.astype(np.float32),
                        b_ctrl=b_ctrl, gate=gate, scored=scored, wcols=wcols,
                        targets=np.asarray(targets), audit=json.dumps(audit),
                        cache_key=key, inputs=json.dumps(inputs["payload"], sort_keys=True))
    audit["cache"] = dict(status="built", key=key, file=_relpath(cache), validated=[],
                          pm_blend=float(PM_BLEND), min_cpm=float(MIN_CPM))
    print(f"[{held}] built: {json.dumps(audit)}", flush=True)
    return dict(held=held, targets=targets, pre=pre, truth=truth, b_ctrl=b_ctrl,
                gate=gate, scored=scored, wcols=wcols, op_path=op_path,
                sources=spec["sources"], panel_positions=panel_positions, audit=audit,
                op_name=op_name, cache_key=key)


class Operator:
    """(C + lam I)^(-gamma/2) with production's candidate-specific on-axis preserve.

    252: the file is identified by sha256 against `OPERATOR_REGISTRY` (or its own metadata
    json), `name=` asserts a specific operator, and `apply` asserts the axis and refuses a
    zero on-axis median instead of substituting 1.0.  See the module docstring.
    """

    def __init__(self, path: Path, *, name: str | None = None):
        path = Path(path).resolve()
        self.path = path
        self.sha256 = E062.sha256_file(path)
        meta_path = path.with_suffix(".json")
        self.metadata = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        entry_name = _REGISTRY_BY_SHA.get(self.sha256)
        if name is not None:
            if name not in OPERATOR_REGISTRY:
                raise OperatorGuardError(
                    f"Operator(name={name!r}): not a registered operator name; "
                    f"registered: {sorted(OPERATOR_REGISTRY)}")
            if entry_name != name:
                e = OPERATOR_REGISTRY[name]
                raise OperatorGuardError(
                    f"SWAPPED OPERATOR: {_relpath(path)} hashes to {self.sha256}"
                    f"{' (= registered ' + repr(entry_name) + ')' if entry_name else ' (unregistered)'}, "
                    f"but name={name!r} pins {e['sha256']} ({e['path']}, max eigenvalue "
                    f"{e['max_eigenvalue']:.6f})")
        if entry_name is not None:
            e = OPERATOR_REGISTRY[entry_name]
            if meta_path.exists():
                meta_sha = E062.sha256_file(meta_path)
                if meta_sha != e["metadata_sha256"]:
                    raise OperatorGuardError(
                        f"{entry_name}: metadata json {_relpath(meta_path)} hashes to {meta_sha}, "
                        f"registry pins {e['metadata_sha256']}")
            self.name = entry_name
            self.axis_sha256 = e["axis_sha256"]
            self.wcols_sha256 = e["wcols_sha256"]
            self.provenance = "OPERATOR_REGISTRY"
        else:
            declared = self.metadata.get("sha256")
            axis_rec = self.metadata.get("axis_sha256")
            wcols_rec = self.metadata.get("wcols_sha256")
            if declared == self.sha256 and (axis_rec or wcols_rec):
                self.name = f"self-described:{path.stem}"
                self.axis_sha256 = axis_rec
                self.wcols_sha256 = wcols_rec
                self.provenance = f"metadata:{meta_path.name}"
            else:
                raise OperatorGuardError(
                    f"UNREGISTERED OPERATOR: {_relpath(path)} (sha256 {self.sha256}) is not in "
                    f"plolo.OPERATOR_REGISTRY {sorted(OPERATOR_REGISTRY)} and its metadata json "
                    f"{'declares sha256 ' + repr(declared) if declared else 'declares no sha256'}"
                    f"{'' if (axis_rec or wcols_rec) else ' and carries no axis record'}. "
                    "Register it by name (sha256 + metadata sha256 + axis digest), or re-fit it "
                    "with scripts/build_union_whitener.py, which writes 'sha256' and 'axis_sha256'.")
        with np.load(path) as z:
            self.evals = np.asarray(z["evals"], dtype=np.float64)
            self.evecs = np.asarray(z["evecs"], dtype=np.float64)
        self.n_cols = int(self.evals.shape[0])
        if entry_name is not None:
            e = OPERATOR_REGISTRY[entry_name]
            if self.n_cols != e["n_cols"] or abs(float(self.evals.max()) - e["max_eigenvalue"]) > 1e-6:
                raise OperatorGuardError(
                    f"{entry_name}: n_cols {self.n_cols} / max eigenvalue {float(self.evals.max()):.9f} "
                    f"!= registry {e['n_cols']} / {e['max_eigenvalue']:.9f}")
        self._axis_ok: dict[str, str] = {}
        self.axis_checked: str | None = None
        print(f"[operator] {self.name} sha256 {self.sha256[:12]} n_cols {self.n_cols} "
              f"max eig {float(self.evals.max()):.6f} ({self.provenance})", flush=True)

    def _check_axis(self, wcols) -> None:
        axis_rec = getattr(self, "axis_sha256", None)
        wcols_rec = getattr(self, "wcols_sha256", None)
        if not (axis_rec or wcols_rec):
            return                       # a bare object (216's RefitOperator) carries no record
        w = np.asarray(wcols, dtype=np.int64)
        n_cols = int(getattr(self, "n_cols", self.evals.shape[0]))
        if w.shape[0] != n_cols:
            raise OperatorGuardError(
                f"{self.name}: applied on {w.shape[0]} columns, fitted on {n_cols}")
        idx = index_sha256(w)
        if idx in self._axis_ok:
            self.axis_checked = self._axis_ok[idx]
            return
        sym = axis_sha256(OFFICIAL[i] for i in w)
        if sym == axis_rec:
            how = "symbols:axis_sha256"
        elif idx == wcols_rec:
            how = "indices:wcols_sha256"
        elif idx == axis_rec:
            how = "indices:axis_sha256"
        else:
            raise OperatorGuardError(
                f"AXIS MISMATCH: {self.name} ({_relpath(self.path)}) was fitted for a DIFFERENT "
                f"axis than the one it is being applied on.\n"
                f"  applied axis: symbols {sym}, indices {idx} ({n_cols} columns)\n"
                f"  fit record:   axis_sha256 {axis_rec}, wcols_sha256 {wcols_rec}\n"
                "  The operator's axis is union(donor coverage) - panel, a function of the "
                "panel; n_cols is 17,942 on every panel measured, so a column count cannot see "
                "this (AGENTS.md 8.3, 235).  Pass check_axis=False ONLY to reproduce a "
                "deliberate trap measurement (233's trapprice233.py).")
        self._axis_ok[idx] = how
        self.axis_checked = how

    def apply(self, M, wcols, gamma, lam, *, check_axis: bool = True):
        out = np.asarray(M, dtype=np.float64).copy()
        if check_axis:
            self._check_axis(wcols)
        norms_in = np.linalg.norm(out[:, wcols], axis=1)
        before = float(np.median(norms_in))
        if not before > 0.0:
            carrying = int((norms_in > 0).sum())
            raise OperatorGuardError(
                f"ON-AXIS PRESERVE UNDEFINED: the median row norm on the operator axis is "
                f"{before} -- only {carrying}/{norms_in.size} rows "
                f"({100.0 * carrying / max(norms_in.size, 1):.1f}% coverage) carry anything on it. "
                "The pre-252 harness substituted 1.0 here and silently returned the whitened rows "
                "un-preserved (252 measured x5.69 at 40% coverage).  Whiten the carrying rows "
                "only, or size the residual on its carrying rows (solve_fraction sizing='carrying').")
        w = (self.evals + lam) ** (-gamma / 2.0)
        t = ((out[:, wcols] @ self.evecs) * w) @ self.evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1)))
        if not after > 0.0:
            raise OperatorGuardError(
                f"ON-AXIS PRESERVE UNDEFINED: the whitened median row norm is {after}")
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


SIZING_RULES = ("all_rows", "carrying", "legacy_all_rows")
SIZING_DEFAULT = "all_rows"


def solve_fraction(whitened, cd4, target_f, *, sizing: str | None = None):
    """Production `solve_cd4_fraction`: f is a ratio of median row norms, linear in fraction.

    sizing (None -> `plolo.SIZING_DEFAULT`, which is "all_rows"):
      "all_rows"         production's A/B/C rule, the median over EVERY panel row.  Bit-identical
                         to the pre-252 harness wherever it returned a positive fraction; where
                         the median is zero (coverage <= 50%) it RAISES `SilentZeroError` instead
                         of returning fraction 0.0 (223 §7.1: the residual read as absent).
      "carrying"         218 §7.1: the median over the rows the residual carries -- the builder's
                         `solve_carrying_row_fraction`, October's D/E/F sizing (247 §3).
      "legacy_all_rows"  the pre-252 behaviour by name: returns (0.0, per_unit) on a zero median,
                         printing a LEGACY SIZING line.  For reproducing archive rows only.
    Returns (fraction, f_per_unit) exactly as before.
    """
    rule = SIZING_DEFAULT if sizing is None else sizing
    if rule not in SIZING_RULES:
        raise ValueError(f"solve_fraction: sizing must be one of {SIZING_RULES}, got {rule!r}")
    added = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0)
    norms = np.linalg.norm(added, axis=1)
    n_rows = int(norms.size)
    live = norms[norms > 0]
    carrying = int(live.size)
    rn_all = float(np.median(norms))
    rn_carry = float(np.median(live)) if live.size else 0.0
    rn = rn_carry if rule == "carrying" else rn_all
    bn = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    per_unit = rn / bn if bn else 0.0
    if per_unit > 0:
        return target_f / per_unit, per_unit
    coverage = carrying / max(n_rows, 1)
    if rule == "legacy_all_rows":
        print(f"[plolo] LEGACY SIZING (sizing='legacy_all_rows'): all-rows median row norm 0.0 with "
              f"{carrying}/{n_rows} carrying rows ({100 * coverage:.1f}% coverage); returning "
              f"fraction 0.0 -- the residual is NOT applied (223's silent zero, by name)", flush=True)
        return 0.0, per_unit
    if not bn > 0:
        raise SilentZeroError(
            f"solve_fraction: the whitened base has median row norm {bn}; nothing to size against")
    if carrying == 0:
        raise SilentZeroError(
            f"solve_fraction(sizing={rule!r}): the residual carries 0/{n_rows} rows -- it is empty, "
            "and a fraction for it measures nothing")
    carry_pu = rn_carry / bn
    raise SilentZeroError(
        f"SILENT ZERO REFUSED: solve_fraction(sizing={rule!r}) would return fraction 0.0 for a "
        f"residual that IS present -- {carrying}/{n_rows} carrying rows ({100 * coverage:.1f}% "
        f"coverage), all-rows median row norm {rn_all} (zero because coverage <= 50%), "
        f"carrying-row median {rn_carry:.6f}, base median {bn:.6f}.  Under the carrying-row rule "
        f"(218 §7.1) f_per_unit = {carry_pu:.6f} and target_f {target_f} solves to "
        f"{target_f / carry_pu:.6f}: pass sizing='carrying'.  To reproduce the pre-252 zero by "
        "name pass sizing='legacy_all_rows'.")

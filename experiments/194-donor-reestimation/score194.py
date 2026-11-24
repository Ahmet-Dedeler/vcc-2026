"""194 — score every donor-table variant on PLOLO (162's instrument).

DESIGN, and why it is the only fair one available
-------------------------------------------------
An estimator change moves BOTH sides of a leave-one-donor-out comparison: the sources and
the truth are built by the same estimator.  If both move, a variant can score better
merely by making the two sides agree on a shared artefact.  So the PRIMARY arm varies the
SOURCE tables only and holds the TRUTH at the incumbent shipped table, which is the
closest stand-in we have for "what the scorer will actually measure".  A secondary arm
varies both and is reported separately.

Everything else is held byte-fixed against 162: the same folds, the same frozen held-out
whitening operators (fitted WITHOUT the scored line), the same PerturbMatch blend, the
same CD4 residual, the same production renderer, the same challenge frames A/B/C.

The baseline arm is `v0` built by THIS pipeline on THIS restricted row set -- not the
shipped table -- because `prepare_main_source` centres each bank on its own mean and
rescales by its own median centred norm, so a table with a different row set is not
comparable to one with the shipped row set.  All arms therefore share a row set and
differ only in the estimator.  `v0` reproduces the shipped table to 4.8e-7 on the rows
they share (known-answer row, RESULTS.md section 0).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, E064, FOLDS, OFFICIAL, PANEL, PANEL_SET, PM_PATH,  # noqa: E402
                   PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints,
                   solve_fraction)
from tables import SumBank  # noqa: E402

# a15/a23 production geometry.  `mass` is a23's exact-mass repair; it is applied to the
# render, not to the effect, and it is held identical across arms.
PROD = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, mass="minnorm")

DONOR_MODE = {"k562_gwps": "restricted", "orion_hct116": "restricted",
              "orion_hek293t": "restricted"}
VARIANT_ROOT = WORK / "tables"


def variant_path(donor: str, variant: str) -> Path:
    """Materialise `donor`'s table under `variant`, cached on disk.

    `full:<v>` selects a FULL-row table shipped into `data/donors/194-<v>-<donor>/`,
    which is directly comparable to the shipped incumbent (the `shipped` arm) because it
    has the incumbent's own row set.  Everything else is a restricted-row screen table,
    comparable only to the restricted `v0`.
    """
    if variant.startswith("full:"):
        p = ROOT / "data/donors" / f"194-{variant.split(':', 1)[1]}-{donor}"
        if not (p / "delta.npy").exists():
            raise SystemExit(f"missing full table {p}")
        return p
    out = VARIANT_ROOT / donor / variant.replace("@", "_at_")
    if not (out / "delta.npy").exists():
        SumBank(donor, DONOR_MODE[donor]).write(variant, out)
    return out


def donor_paths(variant: str, swap: tuple[str, ...]) -> dict[str, Path]:
    p = dict(plolo.DONOR_PATHS)
    for d in swap:
        p[d] = variant_path(d, variant)
    return p


MIN_CELLS = [0]
_CELLS: dict[str, dict] = {}


def _cell_counts(donor: str) -> dict[str, int] | None:
    """Cells per target as the accumulation pass counted them."""
    if donor not in DONOR_MODE:
        return None
    if donor not in _CELLS:
        d = WORK / f"sums-{donor}-{DONOR_MODE[donor]}"
        meta = json.loads((d / "meta.json").read_text())
        c = np.load(d / "_cells.npy")
        _CELLS[donor] = {n: int(c[i]) for i, n in enumerate(meta["names"])}
    return _CELLS[donor]


_PREP: dict[tuple, object] = {}


def _evict_variants() -> None:
    """Keep only the incumbent/PerturbMatch preparations memoised.

    Each `PreparedSource` is `len(targets) x 18,533` float64 (~43 MB), so holding one per
    (variant, donor, fold) reaches 8 GB across a full screen and trips 062's own guard.
    The incumbent side is what is worth caching -- it is identical for every arm.
    """
    for k in [k for k in _PREP if str(VARIANT_ROOT) in k[1]]:
        del _PREP[k]


def _prepare(name: str, path: Path, held: str, fold: dict):
    """`prepare_main_source`, memoised: the incumbent side is identical across arms."""
    key = (name, str(path), held)
    if key not in _PREP:
        bank = (E064.load_perturbmatch(path) if name == "perturbmatch"
                else E062.Bank.load(name, path))
        _PREP[key] = (bank, E062.prepare_main_source(
            bank, fold["targets"], OFFICIAL, fold["panel_positions"],
            chunk_rows=256, max_rss_gib=12.0))
    return _PREP[key]


def build_pre(held: str, fold: dict, variant: str,
              swap: tuple[str, ...]) -> tuple[np.ndarray, dict]:
    """The production pre-whitening effect matrix, from `variant` source tables."""
    spec = FOLDS[held]
    targets = fold["targets"]
    panel_positions = fold["panel_positions"]
    paths = donor_paths(variant, swap)
    pairs = {n: _prepare(n, paths[n], held, fold) for n in spec["sources"]}
    banks = {n: v[0] for n, v in pairs.items()}
    prepared = {n: v[1] for n, v in pairs.items()}
    source_prepared = []
    pm_n = 0
    if "k562_gwps" in banks:
        pm_bank, pm_prepared = _prepare("perturbmatch", PM_PATH, held, fold)
        shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL,
                                             PANEL_SET)
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            pm_n = int(both.sum())
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - plolo.PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + plolo.PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            source_prepared.append(E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(prepared[n].target_available, bool).copy(),
                mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
                raw_median_norm=float(prepared[n].raw_median_norm),
                mapped_gene_count=int(prepared[n].mapped_gene_count),
                duplicate_symbols=tuple(prepared[n].duplicate_symbols)))
        else:
            source_prepared.append(prepared[n])
    if MIN_CELLS[0]:
        # candidate 5's `min_cells`: a source stops voting for a target it measured in
        # fewer than N cells, so the consensus falls back to the other donors.
        # `target_available` is replaced, never mutated -- the PreparedSource objects are
        # memoised across arms and an in-place edit would leak into every later arm.
        rebuilt = []
        for sp in source_prepared:
            counts = _cell_counts(sp.name.replace("_pm50", ""))
            if counts is None:
                rebuilt.append(sp)
                continue
            avail = np.asarray(sp.target_available, bool).copy()
            avail &= np.array([counts.get(t, 0) >= MIN_CELLS[0] for t in targets])
            rebuilt.append(E062.PreparedSource(
                name=sp.name, effects=sp.effects, target_available=avail,
                mean_on_axis=sp.mean_on_axis, raw_median_norm=sp.raw_median_norm,
                mapped_gene_count=sp.mapped_gene_count,
                duplicate_symbols=sp.duplicate_symbols))
        source_prepared = rebuilt
    pre, votes = E062.build_baseline(source_prepared)
    aud = dict(pm_targets=pm_n, min_cells=MIN_CELLS[0],
               coverage={n: int(prepared[n].target_available.sum()) for n in banks},
               raw_median={n: float(prepared[n].raw_median_norm) for n in banks},
               n_rows={n: int(banks[n].shape[0]) for n in banks})
    return pre, aud


def variant_truth(held: str, fold: dict, variant: str) -> np.ndarray:
    """Truth rebuilt under `variant` (secondary arm only)."""
    bank = E062.Bank.load(held, variant_path(held, variant))
    prep = E062.prepare_main_source(bank, fold["targets"], OFFICIAL,
                                    fold["panel_positions"], chunk_rows=256,
                                    max_rss_gib=12.0)
    return E062.raw_truth(prep)


def pds_per_target(P: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Per-target rank statistic, the quantity `endpoints`' `pds` averages.

    Kept separately so the paired target-cluster bootstrap can put an interval on an
    arm-to-arm difference instead of leaning on the between-fold scatter of four folds.
    """
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return (beaten + 0.5 * ties) / (n - 1)


def score(fold: dict, pre: np.ndarray, op: Operator, cd4, truth_cols: np.ndarray,
          truth: np.ndarray) -> dict:
    """Endpoints in the challenge frames, under both shipped render laws.

    `a23` is the shipped one (exact-mass repair).  `a15` (no repair) is carried too
    because the mass repair adds a large, near-common per-gene term that compresses
    `sgn`/`sgn_w` from ~0.50 to ~0.33 on every arm (162's own table), so it is the wrong
    lens through which to compare estimators on the SIGN channel even though it is the
    right one for `pds`.
    """
    eff = pre * PROD["alpha"]
    whitened = op.apply(eff, fold["wcols"], PROD["gamma"], PROD["lam"])
    frac, _ = solve_fraction(whitened, cd4, PROD["cd4_f"])
    final, _ = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    out: dict[str, float] = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc, per = [], []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            acc.append(endpoints(dl[:, truth_cols], truth))
            per.append(pds_per_target(dl[:, truth_cols], truth))
        for k in acc[0]:
            out[f"{tag}_{k}"] = float(np.mean([a[k] for a in acc]))
        if tag == "a23":
            out["_pds_per_target"] = np.mean(per, axis=0).tolist()
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "sgn", "sgn_w", "rr", "sse")})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=["v0"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--truth", default="incumbent", choices=["incumbent", "variant"])
    ap.add_argument("--swap", nargs="*", default=list(DONOR_MODE),
                    help="which source donors are replaced by the variant table")
    ap.add_argument("--tag", default="")
    ap.add_argument("--min-cells", type=int, default=0)
    ap.add_argument("--out", default=str(WORK / "scores.json"))
    args = ap.parse_args()
    swap = tuple(args.swap)
    MIN_CELLS[0] = args.min_cells
    tag = args.tag or "+".join(s.split("_")[0] for s in swap)

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["variant"], r["truth"], r.get("swap_tag")) for r in rows}

    for held in args.folds:
        _PREP.clear()
        fold = build_fold(held)
        cd4 = cd4_for(fold["targets"])
        op = Operator(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                      else PROD_OPS / "whitener_union_pruned.npz")
        base_truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        for v in args.variants:
            if (held, v, args.truth, tag) in done:
                continue
            t0 = time.time()
            _evict_variants()
            if held in DONOR_MODE and args.truth == "variant" and v != "shipped":
                truth = variant_truth(held, fold, v)[:, fold["scored"]]
            else:
                truth = base_truth
            if v == "shipped":
                # the incumbent full tables, exactly as `build_fold` assembled them.
                # Reproducing 162's a23 endpoints from this row is known-answer row 3.
                pre, aud = np.asarray(fold["pre"], dtype=np.float64), {"source": "162"}
            else:
                pre, aud = build_pre(held, fold, v, swap)
            e = score(fold, pre, op, cd4, fold["scored"], truth)
            rows.append(dict(fold=held, variant=v, truth=args.truth, swap_tag=tag,
                             swap=list(swap), n_targets=len(fold["targets"]),
                             **e, audit=aud, seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{v}/{tag}] pds {e['pds']:.5f} cos {e['cos']:.5f} "
                  f"q {e['a15_sgn']:.5f} q_w {e['a15_sgn_w']:.5f} "
                  f"sse {e['sse']:.6f} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

"""209 Phase 2 — score SUBSPACE DEFLATION of the two Orion donors on PLOLO.

162's instrument, byte-unchanged from 204: four folds, challenge frames A/B/C, shipped
geometry (alpha 0.982, gamma 0.65, lam 0, cd4_f 0.10), a23 render for `pds`/`cos`, a15
render for `q`/`q_w`.  `debench.Surface.pds()` is never called (188 Phase 3).

THE ONLY THING THAT VARIES ACROSS ARMS is a projection applied to the two Orion donors'
prepared effect matrices before `build_baseline` averages them:

    effects_i  <-  ( effects_i - lam * (effects_i V^T) V ) * (norm_before / norm_after)

with `V` an orthonormal basis (rows) on the official axis and `lam` the deflation strength.
`lam = 0` is the identity and reproduces 162's cached `pre` exactly -- the known-answer row.

The trailing rescale restores the donor's own median row norm over its available targets.
That is what production would do by itself: `prepare_main_source` divides by the median
centred row norm of the table it is given, and a projection commutes with the gene-mean
subtraction, so deflating the TABLE and letting production renormalise is exactly the
norm-restored form.  The `_raw` arms drop the rescale and are therefore amplitude-changing
(rule 13's +0.0089 gate rather than +0.0019); both are measured.

K562 IS NEVER TOUCHED (mandate falsifier 4).  The audit records, per fold and per arm, the
max abs deviation of the K562 source's effects from the un-deflated build, and it must be
exactly 0.0.

The whitening operator is the fold's frozen one on every arm and is NEVER refitted -- as in
204, so an arm is one variable rather than two.  Refitting it on deflated donors is a
different arm and is not run here.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, E064, FOLDS, OFFICIAL, PANEL_SET, PM_PATH, PROD_OPS,  # noqa: E402
                   Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction)

PROD = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10)
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ORION = ("orion_hct116", "orion_hek293t")
SEED = 20260903


# --------------------------------------------------------------------------- the bases
def _load_basis(name: str, r: int | None = None) -> np.ndarray:
    """An orthonormal basis (rows) on the 18,533-gene official axis."""
    if name == "blob2":
        z = np.load(WORK / "blob-axis.npz")
        return np.asarray(z["basis"], dtype=np.float64)
    z = np.load(WORK / f"subspace-{name}.npz")
    V = np.asarray(z["basis"], dtype=np.float64)
    cols = np.asarray(z["cols"])
    if r is not None:
        V = V[:r]
    full = np.zeros((V.shape[0], len(OFFICIAL)))
    full[:, cols] = V
    return full


def _random_basis(dim: int, source: str, seed: int) -> np.ndarray:
    """A random orthonormal basis inside the SAME column span as the real subspace."""
    cols = np.asarray(np.load(WORK / f"subspace-{source}.npz")["cols"])
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((len(cols), dim))
    Q, _ = np.linalg.qr(G)
    full = np.zeros((dim, len(OFFICIAL)))
    full[:, cols] = Q.T
    return full


def build_arms() -> dict:
    arms: dict[str, dict] = {"equal": dict(kind="identity", lam=0.0,
                                           note="incumbent, known-answer row")}
    for lam in (0.25, 0.5, 0.75, 1.0):
        arms[f"subC_l{lam:g}".replace(".", "")] = dict(
            kind="basis", basis=("C", 2), lam=lam, donors=ORION,
            note=f"excess subspace (C, r=2), strength {lam}")
    for lam in (0.25, 0.5, 0.75, 1.0):
        arms[f"subW_l{lam:g}".replace(".", "")] = dict(
            kind="basis", basis=("WIDE", 2), lam=lam, donors=ORION,
            note=f"excess subspace (WIDE, r=2), strength {lam}")
    for lam in (0.25, 0.5, 0.75):
        arms[f"blob2_l{lam:g}".replace(".", "")] = dict(
            kind="basis", basis=("blob2", None), lam=lam, donors=ORION,
            note=f"the long-gene axis alone, strength {lam}")
    # dimension sensitivity.  `subCx_r2` re-scores the SAME rank-2 arm off the extended
    # basis and must reproduce `subC_l1` exactly -- the built-in check that `subspace-Cext`
    # is the same estimator.
    for r in (1, 2, 5, 10, 20):
        arms[f"subCx_r{r}"] = dict(kind="basis", basis=("Cext", r), lam=1.0, donors=ORION,
                                   note=f"excess subspace (C), dimension sweep r={r}")
    arms["subC_l1_raw"] = dict(kind="basis", basis=("C", 2), lam=1.0, donors=ORION,
                               restore_norm=False,
                               note="excess subspace (C, r=2), NO norm restoration "
                                    "(amplitude-changing)")
    arms["blob2"] = dict(kind="basis", basis=("blob2", None), lam=1.0, donors=ORION,
                         note="the long-gene axis alone: both Orion top PCs, rank 2")
    arms["blob2_then_subC"] = dict(kind="basis_stack",
                                   stack=[("blob2", None), ("CblobDeflated", 3)],
                                   lam=1.0, donors=ORION,
                                   note="blob first, then the excess re-estimated after it")
    arms["subCblobDefl"] = dict(kind="basis", basis=("CblobDeflated", 3), lam=1.0,
                                donors=ORION,
                                note="ONLY the blob-orthogonal excess (mandate falsifier 2)")
    for s in (1, 2, 3):
        arms[f"permC_s{s}"] = dict(kind="random", dim=2, source="C", seed=SEED + s,
                                   lam=1.0, donors=ORION,
                                   note=f"falsifier: a random 2-D subspace, seed {SEED + s}")
        arms[f"permC_s{s}_l05"] = dict(kind="random", dim=2, source="C", seed=SEED + s,
                                       lam=0.5, donors=ORION,
                                       note="falsifier: random 2-D in C, strength 0.5")
        arms[f"permW_s{s}"] = dict(kind="random", dim=2, source="WIDE", seed=SEED + 10 + s,
                                   lam=1.0, donors=ORION,
                                   note="falsifier: random 2-D in the WIDE span, the "
                                        "matched control for `blob2`")
        arms[f"permW_s{s}_l05"] = dict(kind="random", dim=2, source="WIDE",
                                       seed=SEED + 10 + s, lam=0.5, donors=ORION,
                                       note="falsifier: random 2-D in WIDE, strength 0.5 "
                                            "-- the matched control for `blob2_l05`")
    arms["rank1_175"] = dict(kind="own_pc", k=1, lam=1.0, donors=DONORS,
                             note="175's per-donor rank-1 deflation, all three donors")
    arms["rank1_orion"] = dict(kind="own_pc", k=1, lam=1.0, donors=ORION,
                               note="175's rank-1, Orion only (matched to the subspace arms)")
    return arms


# ------------------------------------------------------------------- source preparation
_PREP: dict[tuple, object] = {}


def _prepare(name: str, held: str, fold: dict):
    key = (name, held)
    if key not in _PREP:
        path = PM_PATH if name == "perturbmatch" else plolo.DONOR_PATHS[name]
        bank = (E064.load_perturbmatch(path) if name == "perturbmatch"
                else E062.Bank.load(name, path))
        _PREP[key] = (bank, E062.prepare_main_source(
            bank, fold["targets"], OFFICIAL, fold["panel_positions"],
            chunk_rows=256, max_rss_gib=9.0))
    return _PREP[key]


def fold_sources(held: str, fold: dict):
    """162's sources with the PerturbMatch blend, verbatim from 204's `score204.py`."""
    spec = FOLDS[held]
    pairs = {n: _prepare(n, held, fold) for n in spec["sources"]}
    banks = {n: v[0] for n, v in pairs.items()}
    prepared = {n: v[1] for n, v in pairs.items()}
    out = []
    for n in spec["sources"]:
        if n == "k562_gwps":
            pm_bank, pm_prepared = _prepare("perturbmatch", held, fold)
            shared = E064.shared_k562_pm_columns(banks[n], pm_bank, OFFICIAL, PANEL_SET)
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - plolo.PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + plolo.PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            out.append((n, E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(prepared[n].target_available, bool).copy(),
                mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
                raw_median_norm=float(prepared[n].raw_median_norm),
                mapped_gene_count=int(prepared[n].mapped_gene_count),
                duplicate_symbols=tuple(prepared[n].duplicate_symbols))))
        else:
            out.append((n, prepared[n]))
    return out


# --------------------------------------------------------------------------- deflation
def own_top_pcs(effects, avail, k):
    """175's object, computed on the fold's own rows: the donor's own dominant directions."""
    X = effects[avail]
    Xc = X - X.mean(0, keepdims=True)
    _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Vt[:k]


def deflate(effects, avail, V, lam, restore_norm=True):
    out = np.asarray(effects, dtype=np.float64).copy()
    rows = np.flatnonzero(avail)
    if len(rows) == 0 or V is None or V.shape[0] == 0:
        return out, dict(before=0.0, after=0.0, scale=1.0)
    before = float(np.median(np.linalg.norm(out[rows], axis=1)))
    out[rows] -= lam * ((out[rows] @ V.T) @ V)
    after = float(np.median(np.linalg.norm(out[rows], axis=1)))
    scale = (before / after) if (restore_norm and after > 0) else 1.0
    if scale != 1.0:
        out[rows] *= scale
    return out, dict(before=before, after=after, scale=scale)


def deflated_sources(srcs, arm: str, spec: dict):
    """Return (sources, audit).  Sources not named in the arm are passed through untouched."""
    lam = float(spec.get("lam", 0.0))
    donors = set(spec.get("donors", ()))
    restore = bool(spec.get("restore_norm", True))
    audit: dict = {"lam": lam, "donors": sorted(donors), "restore_norm": restore,
                   "per_donor": {}}
    if spec["kind"] == "identity" or lam == 0.0:
        return srcs, audit

    V_static = None
    if spec["kind"] == "basis":
        V_static = _load_basis(*spec["basis"])
    elif spec["kind"] == "basis_stack":
        V_static = np.vstack([_load_basis(*b) for b in spec["stack"]])
        Q, _ = np.linalg.qr(V_static.T)
        V_static = Q.T
    elif spec["kind"] == "random":
        V_static = _random_basis(spec["dim"], spec["source"], spec["seed"])
    if V_static is not None:
        audit["dim"] = int(V_static.shape[0])
        audit["orthonormality_err"] = float(
            np.abs(V_static @ V_static.T - np.eye(V_static.shape[0])).max())

    out = []
    for name, sp in srcs:
        base = name.split("_pm")[0] if name.startswith("k562_gwps") else name
        if base not in donors:
            out.append((name, sp))
            continue
        avail = np.asarray(sp.target_available, bool)
        V = (own_top_pcs(sp.effects, avail, spec["k"]) if spec["kind"] == "own_pc"
             else V_static)
        eff, a = deflate(sp.effects, avail, V, lam, restore)
        a["dim"] = int(V.shape[0])
        a["removed_energy_frac"] = float(
            1.0 - (eff[avail] ** 2).sum() / max((np.asarray(sp.effects)[avail] ** 2).sum(),
                                                1e-300))
        audit["per_donor"][base] = a
        out.append((name, E062.PreparedSource(
            name=sp.name, effects=eff,
            target_available=avail.copy(),
            mean_on_axis=np.asarray(sp.mean_on_axis, float).copy(),
            raw_median_norm=float(sp.raw_median_norm),
            mapped_gene_count=int(sp.mapped_gene_count),
            duplicate_symbols=tuple(sp.duplicate_symbols))))
    return out, audit


def build_baseline_equal(sources):
    """`build_baseline` at equal weights, expressed exactly as 204's `weighted_baseline`."""
    n_t, n_g = sources[0][1].effects.shape
    num = np.zeros((n_t, n_g))
    den = np.zeros(n_t)
    for _, sp in sources:
        rows = np.asarray(sp.target_available, bool)
        num[rows] += sp.effects[rows]
        den += rows.astype(float)
    covered = den > 0
    num[covered] /= den[covered, None]
    if (~covered).any():
        mean_response = sum(sp.mean_on_axis for _, sp in sources) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        num[~covered] = plolo.EFFECT_SCALE * mean_response / mean_norm
    return num, dict(n_covered=int(covered.sum()))


def _pds_per_target(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def score(fold, pre, op, cd4, truth_cols, truth) -> dict:
    eff = pre * PROD["alpha"]
    whitened = op.apply(eff, fold["wcols"], PROD["gamma"], PROD["lam"])
    frac, _ = solve_fraction(whitened, cd4, PROD["cd4_f"])
    final, _ = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    out: dict[str, float] = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc, per, perc = [], [], []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            acc.append(endpoints(dl[:, truth_cols], truth))
            p, c = _pds_per_target(dl[:, truth_cols], truth)
            per.append(p)
            perc.append(c)
        for k in acc[0]:
            out[f"{tag}_{k}"] = float(np.mean([a[k] for a in acc]))
        out[f"_pds_per_target_{tag}"] = np.mean(per, axis=0).tolist()
        out[f"_cos_per_target_{tag}"] = np.mean(perc, axis=0).tolist()
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "sgn", "sgn_w", "rr", "sse")})
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    out["whitened_median_row_norm"] = float(np.median(np.linalg.norm(whitened, axis=1)))
    out["cd4_fraction"] = float(frac)
    return out


def main() -> None:
    ARMS = build_arms()
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=list(ARMS))
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "scores209.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        _PREP.clear()
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        base_pre = np.asarray(fold["pre"], dtype=np.float64)
        srcs = fold_sources(held, fold)          # stream FIRST (062's guard reads peak RSS)
        k562_ref = {n: np.asarray(sp.effects).copy() for n, sp in srcs
                    if n.startswith("k562_gwps")}
        cd4 = cd4_for(fold["targets"])
        op = Operator(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                      else PROD_OPS / "whitener_union_pruned.npz")
        for arm in args.arms:
            if (held, arm) in done:
                continue
            t0 = time.time()
            spec = ARMS[arm]
            dsrcs, aud = deflated_sources(srcs, arm, spec)
            # mandate falsifier 4: K562's own rows must be bit-identical on every arm
            k562_dev = 0.0
            for n, sp in dsrcs:
                if n in k562_ref:
                    k562_dev = max(k562_dev, float(
                        np.abs(np.asarray(sp.effects) - k562_ref[n]).max()))
            aud["k562_max_abs_dev"] = k562_dev
            # mandate falsifier 4 applies to every Orion-only arm.  `rank1_175` is 175's
            # own arm and deflates all three donors by construction, so it is exempt and
            # is the only arm allowed a non-zero deviation.
            if "k562_gwps" not in set(spec.get("donors", ())):
                assert k562_dev == 0.0, f"{arm}: K562 was modified ({k562_dev})"
            pre, b = build_baseline_equal(dsrcs)
            aud.update(b)
            aud["max_abs_dev_from_162_pre"] = float(np.max(np.abs(pre - base_pre)))
            aud["pre_median_row_norm"] = float(np.median(np.linalg.norm(pre, axis=1)))
            e = score(fold, pre, op, cd4, fold["scored"], truth)
            rows.append(dict(fold=held, arm=arm, n_targets=len(fold["targets"]),
                             note=spec.get("note", ""), **e, audit=aud,
                             seconds=round(time.time() - t0, 1),
                             peak_rss_gib=resource.getrusage(
                                 resource.RUSAGE_SELF).ru_maxrss / (1 << 30)))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] pds {e['pds']:.6f} a15pds {e['a15_pds']:.6f} "
                  f"cos {e['cos']:.6f} q {e['a15_sgn']:.6f} q_w {e['a15_sgn_w']:.6f} "
                  f"rownorm {e['median_row_norm_final']:.4f} "
                  f"dev {aud['max_abs_dev_from_162_pre']:.2e} "
                  f"k562dev {k562_dev:.1e} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

"""204 Phase 2 — score DONOR WEIGHTINGS on PLOLO (162's instrument, unchanged).

The only thing that varies between arms is the vector of donor weights inside
`build_baseline`.  Everything else -- the shipped tables, the PerturbMatch blend, the
frozen held-out whitening operators, alpha/gamma/lam, the CD4 residual, the renderer, the
challenge frames A/B/C, the read-back -- is byte-identical to 162 and to 194.

The incumbent is `w = (1, 1, 1)`.  `build_baseline` sums the sources that carry a target
and divides by the NUMBER of them; the weighted form sums `w_i * e_i` and divides by
`sum of w_i over the sources that carry the target`.  At `w = 1` the two are the same
expression, so the `equal` arm is a known-answer row against 162's own cached `pre`.

Two properties of the incumbent are deliberately preserved, because changing them would
make the arm two variables sold as one (rule 10):

  * 039's implicit coverage shrinkage.  A donor that does not measure a gene contributes
    a zero on it and still counts in the denominator, so a gene seen by one donor is
    emitted at a fraction of full magnitude.  039 measured that "fixing" this costs
    0.019 macro PDS.  Re-weighting necessarily moves that fraction -- that is a real
    consequence of re-weighting, not a bug -- and §4 of RESULTS.md prices it separately.
  * the uncovered-target fallback (`EFFECT_SCALE * mean_response / |mean_response|`),
    which is left at the incumbent's unweighted donor mean.

    uv run python score204.py --arms equal a_self deatten ...
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

from plolo import (E062, E064, FOLDS, OFFICIAL, PANEL_SET, PM_PATH, PROD_OPS,  # noqa: E402
                   Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction)

PROD = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10)
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
SPAN_PDS = 0.452          # 162's board reference span for `pds`
CALIB = 1.011             # 162's fitted PLOLO -> board scale

# 194 §1 / Phase 0 of this experiment: the three donors' achieved repression.
RESID = dict(k562_gwps=0.155, orion_hct116=0.300, orion_hek293t=0.510)
A_SELF = {d: 1.0 - RESID[d] for d in DONORS}            # (0.845, 0.700, 0.490)

ARMS: dict[str, dict] = {
    # ---- identity ---------------------------------------------------------------------
    "equal": dict(w=(1.0, 1.0, 1.0), note="incumbent, known-answer row"),
    # ---- derived from measured repression -----------------------------------------------
    "a_self": dict(w=tuple(A_SELF[d] for d in DONORS),
                   note="w ∝ a_i, a from the self-gene residual (194 / Phase 0)"),
    "deatten": dict(w=tuple(1.0 / A_SELF[d] for d in DONORS),
                    note="de-attenuate then average: divide each donor by its own a_i"),
    # ---- falsifier 1: permuted attenuation ---------------------------------------------
    "perm_a_1": dict(w=(0.700, 0.490, 0.845), note="falsifier: a cycled"),
    "perm_a_2": dict(w=(0.490, 0.845, 0.700), note="falsifier: a cycled the other way"),
    "perm_a_3": dict(w=(0.700, 0.845, 0.490), note="falsifier: swap the two Orion-adjacent"),
    "perm_a_4": dict(w=(0.490, 0.700, 0.845), note="falsifier: a reversed"),
    "perm_a_5": dict(w=(0.845, 0.490, 0.700), note="falsifier: swap the two Orion lines"),
    # ---- 019 / 010's tuned shares, re-measured on the modern instrument ------------------
    "share_4_1_1": dict(w=(4.0, 1.0, 1.0), note="010's dominant-primary optimum"),
    "share_2_1_1": dict(w=(2.0, 1.0, 1.0), note="019's 2:1"),
    "share_1_2_2": dict(w=(1.0, 2.0, 2.0), note="Orion-primary, the mirror image"),
}

# A 1-D sweep through the incumbent along each of the two axes a weighting can move on,
# so the derived vectors can be located ON the measured surface rather than only compared
# to one point.  `k<x>` moves K562's share against a fixed Orion pair; `h<x>` moves
# HEK293T -- the most attenuated donor -- against a fixed K562/HCT116 pair.
for _x in (0.5, 0.7, 0.85, 1.15, 1.3, 1.5):
    ARMS[f"k{_x:g}".replace(".", "")] = dict(
        w=(_x, 1.0, 1.0), note=f"sweep: K562 share x{_x}")
    ARMS[f"h{_x:g}".replace(".", "")] = dict(
        w=(1.0, 1.0, _x), note=f"sweep: HEK293T share x{_x}")


def set_measured_arms(path: Path = WORK / "atten2_204.json") -> None:
    """Arms whose weights are MEASURED in Phase 1 are filled in from its output."""
    if not path.exists():
        return
    a = json.loads(path.read_text())
    tri = a["routeB_triple_gene_controlled"]["common_genes"]["shared"]
    ARMS["a_slope"] = dict(w=tuple(tri["a_rel"][d] for d in DONORS),
                           note="w ∝ a_i, a from the 3-donor shared-target Gram, "
                                "gene-coverage controlled")
    ARMS["minvar"] = dict(w=tuple(tri["w_minvar_rel"][d] for d in DONORS),
                          note="w ∝ a_i/σ_i², minimum-variance, from the same Gram")
    ARMS["deatten_slope"] = dict(
        w=tuple(1.0 / tri["a_rel"][d] for d in DONORS),
        note="de-attenuate then average, a from the Gram route")
    if "h1_2025" in a.get("routeB_instrument_gene_controlled", {}):
        h = a["routeB_instrument_gene_controlled"]["h1_2025"]["a_rel"]
        ARMS["a_inst"] = dict(w=tuple(h[d] for d in DONORS),
                              note="w ∝ a_i from the h1_2025 fourth-line instrument")
    # minimum-variance built on ROUTE A's attenuation and the Gram's noise energy:
    # the mandate's own a_i/σ_i² with a measured the way the mandate measures it.
    ne = tri["noise_energy"]
    base = A_SELF[DONORS[0]] / ne[DONORS[0]]
    ARMS["minvar_self"] = dict(
        w=tuple((A_SELF[d] / ne[d]) / base for d in DONORS),
        note="w ∝ a_i/σ_i², a from the self-gene residual, σ² from the Gram")


# --------------------------------------------------------------------- per-target weights
_RESID_CACHE: dict[str, tuple[dict[str, float], float]] = {}


def per_target_resid(donor: str) -> dict[str, float]:
    if donor not in _RESID_CACHE:
        r = np.load(WORK / f"resid-{donor}.npy")
        names = (WORK / f"resid-{donor}.txt").read_text().splitlines()
        _RESID_CACHE[donor] = ({n: float(v) for n, v in zip(names, r)}, RESID[donor])
    return _RESID_CACHE[donor][0]


def target_weights(arm: str, donor: str, targets, w_scalar: float) -> np.ndarray:
    """Per-target weight vector for one donor.  Scalar arms broadcast `w_scalar`."""
    n = len(targets)
    if not arm.endswith("_pt"):
        return np.full(n, w_scalar)
    r = per_target_resid(donor)
    med = RESID[donor]
    a = np.array([1.0 - r.get(t, med) for t in targets])
    a = np.clip(a, 0.25, 1.0)                      # 194's floor on the same quantity
    return a if arm.startswith("a_self") else 1.0 / a


PT_ARMS = {
    "a_self_pt": "per-target w ∝ (1 − r_it), floor 0.25",
    "deatten_pt": "per-target de-attenuation w ∝ 1/(1 − r_it), floor 0.25 (194's kdrescale, "
                  "expressed in the mixing frame)",
}

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


def fold_sources(held: str, fold: dict) -> list[tuple[str, object]]:
    """The fold's prepared sources, with the PerturbMatch blend applied exactly as 162."""
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


def weighted_baseline(sources, targets, arm: str, w: dict[str, float]):
    """`build_baseline` with a per-(donor, target) weight.  w == 1 reproduces it exactly."""
    n_t, n_g = sources[0][1].effects.shape
    num = np.zeros((n_t, n_g))
    den = np.zeros(n_t)
    for donor, sp in sources:
        rows = np.asarray(sp.target_available, bool)
        wt = target_weights(arm, donor, targets, w[donor])
        num[rows] += wt[rows, None] * sp.effects[rows]
        den += np.where(rows, wt, 0.0)
    covered = den > 0
    num[covered] /= den[covered, None]
    if (~covered).any():                                  # incumbent fallback, unweighted
        mean_response = sum(sp.mean_on_axis for _, sp in sources) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        num[~covered] = plolo.EFFECT_SCALE * mean_response / mean_norm
    return num, dict(n_covered=int(covered.sum()),
                     mean_denominator=float(den[covered].mean()))


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
    return out


def _pds_per_target(P, T):
    """Per-target rank statistic AND per-target cosine, the two quantities `endpoints`
    macro-averages, kept separately so a paired target bootstrap can put an interval on
    an arm-to-arm difference rather than on the scatter of four folds."""
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def main() -> None:
    set_measured_arms()
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["equal"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "scores204.json"))
    args = ap.parse_args()

    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        _PREP.clear()
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        base_pre = np.asarray(fold["pre"], dtype=np.float64)
        # ORDER MATTERS: 062's guard reads PEAK RSS, so the 2.6 GiB operator and the CD4
        # banks must be allocated AFTER the donor streaming passes, not before.
        srcs = fold_sources(held, fold)
        cd4 = cd4_for(fold["targets"])
        op = Operator(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                      else PROD_OPS / "whitener_union_pruned.npz")
        for arm in args.arms:
            if (held, arm) in done:
                continue
            t0 = time.time()
            spec = ARMS.get(arm)
            w = ({d: 1.0 for d in DONORS} if arm in PT_ARMS
                 else {d: float(spec["w"][i]) for i, d in enumerate(DONORS)})
            pre, aud = weighted_baseline(srcs, fold["targets"], arm, w)
            aud["max_abs_dev_from_162_pre"] = float(np.max(np.abs(pre - base_pre)))
            aud["weights"] = w
            aud["sources"] = [n for n, _ in srcs]
            e = score(fold, pre, op, cd4, fold["scored"], truth)
            rows.append(dict(fold=held, arm=arm, n_targets=len(fold["targets"]),
                             note=(spec or PT_ARMS).get("note", PT_ARMS.get(arm, "")),
                             **e, audit=aud, seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] pds {e['pds']:.6f} a15pds {e['a15_pds']:.6f} "
                  f"cos {e['cos']:.6f} q {e['a15_sgn']:.6f} q_w {e['a15_sgn_w']:.6f} "
                  f"sse {e['sse']:.6f} dev {aud['max_abs_dev_from_162_pre']:.2e} "
                  f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

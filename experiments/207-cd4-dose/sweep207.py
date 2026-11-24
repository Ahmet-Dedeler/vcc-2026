"""207 Phase 0 + Phase 1 — sweep the CD4 dose against COSINE, on 162's PLOLO instrument.

The instrument is imported verbatim from `experiments/162-retrodicting-validator/plolo.py`
(cached folds copied bit-identical into `work/plolo/`).  Nothing about the donor tables, the
whitening operators, the PerturbMatch blend, alpha/gamma/lam, the renderer or the read-back
is touched.  The ONLY thing that varies between arms here is the CD4 residual.

DOSE FRAMES -- there are three and conflating them is the documented trap
------------------------------------------------------------------------
  * `fraction`  -- 062's knob.  `add_post_whiten_residual` adds `fraction * cd4.mean_effect`.
                   062 selected `fraction = 0.05`.
  * `f`         -- 099's knob (`VCC_CD4_TARGET_F`).  The ratio of the added residual's MEDIAN
                   ROW NORM to the whitened base's median row norm.  Production ships `f = 0.10`,
                   which the builder solves to `fraction = 0.0947` (meta.json, abc-knee).
  * `effective_f_at_cd4_fraction` = 0.0528 -- a DIAGNOSTIC printed by the builder: the `f` you
                   would get at the LEGACY default fraction 0.05.  It is not shipped and never
                   was.  a14 shipped `fraction = 0.05` (= f 0.0528); a15 shipped `f = 0.10`
                   (= fraction 0.0947).  So the a14 -> a15 board arm, which 197 measured at
                   +6.8% cosine, is a DOSE INCREASE of 1.89x, not the introduction of CD4.

Endpoints, per arm, per fold:
  `pds`      exact rank statistic, 162's `endpoints`
  `cos`      macro mean of per-target cosines (162's `cos`)
  `cos_agg`  sqrt(1 - sse) = the ENERGY-weighted aggregate cosine, which is the same functional
             197 inverts off the board (`X / sqrt(A * SigmaD)`).  THIS IS THE PRIMARY ENDPOINT.
  `sgn_w`    |truth|-weighted sign accuracy = 187's `q_w`.

    uv run python sweep207.py --arms all
"""
from __future__ import annotations

import argparse
import hashlib
import json
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

from plolo import (E062, FOLDS, OFFICIAL, PANEL_SET, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, deliver, endpoints, solve_fraction)

ALPHA = 0.982
GAMMA = 0.65
LAM = 0.0
SPAN_PDS = 0.452           # 162's board reference span for `pds`
SHIPPED_F = 0.10           # production `VCC_CD4_TARGET_F`
LEGACY_FRACTION = 0.05     # 062's selected `fraction`, shipped in a14

# ------------------------------------------------------------------------------ dose grid
# Log-spaced from the shipped point in both directions, plus 0 and the two historical points.
DOSE_F = [0.0, 0.01, 0.02, 0.035, 0.05, 0.07, 0.10, 0.14, 0.17, 0.20, 0.24, 0.28,
          0.34, 0.40, 0.47, 0.55, 0.68, 0.80, 1.0, 1.3, 1.7, 2.2, 3.0, 4.5, 7.0]


def sha8(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


# --------------------------------------------------------------------------- residual objects
def shuffled_cd4(cd4, seed: int = 20726):
    """Falsifier 1: permute WHOLE ROWS (effect + support together) across targets.

    The geometry, the per-row norms, the gene support pattern and the number of carrying rows
    are all preserved exactly; only the target->row assignment is destroyed.  If this buys the
    same cosine, the gain is not biology.
    """
    rng = np.random.default_rng(seed)
    n = cd4.mean_effect.shape[0]
    carry = np.flatnonzero(np.asarray(cd4.target_condition_count) == 3)
    perm = np.arange(n)
    shuffled = carry.copy()
    while True:                                    # derangement on the carrying rows
        rng.shuffle(shuffled)
        if len(carry) < 2 or not np.any(shuffled == carry):
            break
    perm[carry] = shuffled
    return E062.CD4Prepared(
        mode=cd4.mode + "+shuffled",
        mean_effect=np.asarray(cd4.mean_effect)[perm].copy(),
        gene_votes=np.asarray(cd4.gene_votes)[perm].copy(),
        target_condition_count=np.asarray(cd4.target_condition_count)[perm].copy(),
        condition_details=cd4.condition_details)


def self_residual(cd4, pre):
    """Falsifier: 099's own control -- add back the BASE MIXTURE's un-whitened effect.

    Same support mask, same rows carrying, same `f` solve.  099 measured that this reproduces
    most of CD4's direction curve (x1.122/x1.123 against CD4's x1.134/x1.077), i.e. that most
    of the residual's value is whitening RELAXATION, not the CD4 data.
    """
    eff = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(pre, dtype=np.float64), 0.0)
    return E062.CD4Prepared(
        mode="self_unwhitened", mean_effect=eff,
        gene_votes=np.asarray(cd4.gene_votes).copy(),
        target_condition_count=np.asarray(cd4.target_condition_count).copy(),
        condition_details=[{"condition": "base_mixture_unwhitened"}])


_ALT_CACHE: dict[tuple[str, tuple], np.ndarray] = {}


def alt_line_residual(cd4, name: str, targets):
    """Falsifier 2 (lineage): the same construction with an UNMATCHED extra line.

    KOLF2 is an iPSC line covering 273/300 of the panel and contributed no row to any
    production operator.  If a non-T-lineage fourth line buys the same cosine at the same `f`,
    the T-lineage story is wrong.  Never run on the fold where that line is the truth.
    """
    key = (name, tuple(targets))
    if key not in _ALT_CACHE:
        bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
        official_pos = {g: i for i, g in enumerate(OFFICIAL)}
        panel_positions = np.asarray([official_pos[g] for g in PANEL_SET if g in official_pos],
                                     dtype=np.int64)
        prep = E062.prepare_main_source(bank, targets, OFFICIAL, panel_positions,
                                        chunk_rows=256, max_rss_gib=9.0)
        _ALT_CACHE[key] = (np.asarray(prep.effects, dtype=np.float64),
                           np.asarray(prep.target_available, bool))
    eff, avail = _ALT_CACHE[key]
    votes = np.asarray(cd4.gene_votes).copy()
    votes[~avail] = 0
    out = np.where(votes > 0, eff, 0.0)
    counts = np.asarray(cd4.target_condition_count).copy()
    counts[~avail] = 0
    return E062.CD4Prepared(mode=f"alt_{name}", mean_effect=out, gene_votes=votes,
                            target_condition_count=counts,
                            condition_details=[{"condition": name}])


# ------------------------------------------------------------------------------------ score
def _per_target(P, T):
    """162's `endpoints` internals kept per target, so a PAIRED TARGET bootstrap can put an
    interval on an arm-to-arm difference rather than on the scatter of four folds (204)."""
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def score_effects(final, fold, truth, per_context=True) -> dict:
    out: dict[str, object] = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc, per, perc = {}, [], []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            e = endpoints(dl[:, fold["scored"]], truth)
            e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
            acc[ctx] = e
            p, c = _per_target(dl[:, fold["scored"]], truth)
            per.append(p)
            perc.append(c)
            if per_context:
                for k, v in e.items():
                    out[f"{tag}_{ctx}_{k}"] = float(v)
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in "ABC"]))
        out[f"_pds_per_target_{tag}"] = np.mean(per, axis=0).tolist()
        out[f"_cos_per_target_{tag}"] = np.mean(perc, axis=0).tolist()
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    out["median_row_norm_wcols"] = float(
        np.median(np.linalg.norm(final[:, fold["wcols"]], axis=1)))
    return out


def build_arms(dose_f) -> dict[str, dict]:
    arms: dict[str, dict] = {
        # --- identity / known-answer -------------------------------------------------------
        "a14": dict(kind="cd4", fraction=LEGACY_FRACTION,
                    note="a14: 062's `fraction`=0.05 (f=0.0528). Board raw mse 0.882775, "
                         "cos 0.07418"),
        "a15": dict(kind="cd4", f=SHIPPED_F,
                    note="a15/rc_hybrid SHIPPED: 099's f=0.10 -> fraction 0.0947. Board raw "
                         "mse 0.873805, cos 0.079249"),
    }
    for v in dose_f:
        arms[f"f{v:g}"] = dict(kind="cd4", f=(v if v > 0 else None),
                               note=f"dose sweep, 099 frame f={v:g}")
    arms["f0"]["note"] = "dose sweep, CD4 OFF (f=0)"
    # --- falsifiers, at every dose ---------------------------------------------------------
    for v in dose_f:
        if v <= 0:
            continue
        arms[f"shuf{v:g}"] = dict(kind="shuffled", f=v, note=f"target-shuffled CD4, f={v:g}")
    for v in (0.05, 0.10, 0.20, 0.40):
        arms[f"self{v:g}"] = dict(kind="self", f=v, note=f"099's self-residual control, f={v:g}")
        arms[f"kolf{v:g}"] = dict(kind="alt", line="kolf2", f=v,
                                  note=f"unmatched-lineage line as the residual, f={v:g}")
    return arms


def residual_for(kind: str, cd4, pre, targets, line: str | None):
    if kind == "cd4":
        return cd4
    if kind == "shuffled":
        return shuffled_cd4(cd4)
    if kind == "self":
        return self_residual(cd4, pre)
    if kind == "alt":
        return alt_line_residual(cd4, line, targets)
    raise ValueError(kind)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["all"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "sweep207.json"))
    args = ap.parse_args()

    arms = build_arms(DOSE_F)
    names = list(arms) if args.arms == ["all"] else args.arms
    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(fold["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        op_meta = dict(path=str(op_path.relative_to(ROOT)), sha8=sha8(op_path),
                       shape=[int(op.evecs.shape[0]), int(op.evecs.shape[1])],
                       n_wcols=int(len(fold["wcols"])))
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

        for arm in names:
            if (held, arm) in done:
                continue
            spec = arms[arm]
            if spec.get("kind") == "alt" and spec["line"] == held:
                continue                                    # never leak the truth
            t0 = time.time()
            resid = residual_for(spec["kind"], cd4, pre, fold["targets"], spec.get("line"))
            per_unit = None
            if spec.get("f") is not None:
                frac, per_unit = solve_fraction(whitened, resid, spec["f"])
            else:
                frac = float(spec.get("fraction", 0.0))
            if frac > 0.0:
                final, preserve = E062.add_post_whiten_residual(
                    whitened, resid, frac, fold["wcols"])
            else:
                final, preserve = whitened.copy(), {"post_add_global_preserve_scale": 1.0}
            final[:, fold["panel_positions"]] = 0.0
            e = score_effects(final, fold, truth)
            aud = dict(kind=spec["kind"], fraction=float(frac),
                       f_per_unit=(float(per_unit) if per_unit else None),
                       effective_f=(float(frac * per_unit) if per_unit else None),
                       preserve_scale=float(preserve["post_add_global_preserve_scale"]),
                       carrying_rows=int(np.sum(np.asarray(resid.target_condition_count) > 0)),
                       operator=op_meta)
            rows.append(dict(fold=held, arm=arm, n_targets=len(fold["targets"]),
                             note=spec["note"], **e, audit=aud,
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=2))
            print(f"[{held}/{arm}] frac {frac:.5f} pds {e['pds']:.6f} "
                  f"cos_agg {e['cos_agg']:.6f} cos {e['cos']:.6f} q_w {e['a15_sgn_w']:.6f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        del op, whitened
        _ALT_CACHE.clear()


if __name__ == "__main__":
    main()

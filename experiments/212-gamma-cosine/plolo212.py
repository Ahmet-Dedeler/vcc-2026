"""212 Phase 2 — sweep `WHITEN_GAMMA` on PLOLO with the COSINE as the primary endpoint.

WHY THIS INSTRUMENT AND NOT 170's
---------------------------------
170 swept gamma and priced it on `mse` (through 096's `claim` route) and `pds`. Both are
amplitude-sensitive in ways the *direction* is not:

  * `mse`'s numerator is `s^2 A + Kc - min(claim, s^2 B) - 2 s X`. 170 held `alpha = 0.982`
    while lowering gamma raised `A` by 14%, so what it measured is a gamma move CONFOUNDED
    with a 14% amplitude move that the metric refuses to refund above the rho cliff.
  * `pds` is a rank statistic and 197 measured that we extract 26% MORE retrieval per unit
    of bulk cosine than the median shrunk team -- which is exactly what whitening buys and
    exactly why the `pds` optimum need not be the cosine optimum.

The transfer cosine is **exactly invariant** to any global rescale of the prediction, so it
isolates the one channel AGENTS.md 8.2 says is the entire remaining gap. 204 measured this
instrument's resolution on the same folds: +-0.000044 absolute per arm (4,000-draw paired
target bootstrap), i.e. +-0.20% relative.

WHAT VARIES
-----------
`gamma` alone, `lam = 0`, arm `equal` (the incumbent's donor weights). Everything else --
the four folds, the frozen per-fold operators, PerturbMatch, the CD4 residual re-solved at
each gamma exactly as production does, the renderer, the three challenge frames, the
read-back -- is 162's/204's byte-unchanged code, imported not copied.

KNOWN-ANSWER ROW (rule 10's corollary). `gamma = 0.65` must reproduce 204's own published
`equal` row on all four folds, to 1e-12:

    k562_gwps      a23_pds 0.7891213079583496   a15_pds 0.7957840906989219
                   a23_cos 0.03146796694667673  q 0.5016431146906033
                                                q_w 0.5048799395717629
    orion_hct116   a23_pds 0.7918227737342401   a23_cos 0.026040412606551395
    orion_hek293t  a23_pds 0.7773753414904073   a23_cos 0.02247089282968584
    kolf2          a23_pds 0.6101280255692021   a23_cos 0.007559021017607683

If any of those moves, the pipeline is not 204's and nothing below may be read.

    uv run python plolo212.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))

import plolo  # noqa: E402

# 162's fold caches live in its own work dir and 204 copied them into `204/work/plolo`.
# They are CLONED into this experiment's own work dir (APFS `cp -c`, byte-identical) so that
# nothing this script does can write outside `experiments/212-gamma-cosine/`.
plolo.WORK = WORK / "plolo"
assert (plolo.WORK / "fold-kolf2.npz").exists(), plolo.WORK

import score204  # noqa: E402  (imports plolo AFTER the WORK redirect)

score204.WORK = plolo.WORK.parent            # only read: resid-*.npy, atten2_204.json

from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver  # noqa: E402
from plolo import endpoints, solve_fraction, E062  # noqa: E402
from score204 import DONORS, PROD, _pds_per_target, fold_sources, weighted_baseline  # noqa: E402

GRID = [0.00, 0.30, 0.35, 0.40, 0.45, 0.50, 0.525, 0.55, 0.575, 0.60, 0.625,
        0.65, 0.70, 0.75, 0.80]

# rule 11: path + sha256 for every operator this script touches, asserted before use.
OP_SHA = {
    "held_k562_gwps.npz": None,          # filled at run time and recorded, see below
    "held_orion_hct116.npz": None,
    "held_orion_hek293t.npz": None,
    "whitener_union_pruned.npz":
        "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
}

# 204's own `max_abs_dev_from_162_pre` for the `equal` arm on each fold: the float32
# round-trip of 162's cached `pre`, reproduced exactly here as a second known-answer row.
DEV_KA = {
    "k562_gwps": 4.132209086016303e-08,
    "orion_hct116": 5.510534606401052e-08,
    "orion_hek293t": 5.822390525267451e-08,
    "kolf2": 5.6732583963992056e-08,
}

KA = {
    "k562_gwps": dict(a23_pds=0.7891213079583496, a15_pds=0.7957840906989219,
                      a23_cos=0.03146796694667673, a15_sgn=0.5016431146906033,
                      a15_sgn_w=0.5048799395717629),
    "orion_hct116": dict(a23_pds=0.7918227737342401, a23_cos=0.026040412606551395),
    "orion_hek293t": dict(a23_pds=0.7773753414904073, a23_cos=0.02247089282968584),
    "kolf2": dict(a23_pds=0.6101280255692021, a23_cos=0.007559021017607683),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


def score_at(fold, pre, op, cd4, truth_cols, truth, gamma, lam=0.0) -> dict:
    """204's `score()` with gamma as an argument instead of a module constant.

    Byte-for-byte the same operations in the same order; `PROD['gamma']` replaced by
    `gamma`. The CD4 fraction is RE-SOLVED at each gamma, which is what production does
    (`solve_cd4_fraction` is a ratio of median row norms and the whitened norm moves).
    """
    eff = pre * PROD["alpha"]
    whitened = op.apply(eff, fold["wcols"], gamma, lam)
    frac, per_unit = solve_fraction(whitened, cd4, PROD["cd4_f"])
    final, _ = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    out: dict = {}
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
    out["cd4_solved_fraction"] = float(frac)
    out["cd4_f_per_unit"] = float(per_unit)
    return out


def main() -> None:
    out_path = WORK / "plolo212.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["gamma"]) for r in rows}
    ops_seen = {}

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        srcs = fold_sources(held, fold)
        cd4 = cd4_for(fold["targets"])
        op_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                       else PROD_OPS / "whitener_union_pruned.npz")
        sha = sha256_file(op_path)
        want = OP_SHA.get(op_path.name)
        if want is not None and sha != want:
            raise SystemExit(f"operator sha drift on {op_path}: {sha} != {want}")
        op = Operator(op_path)
        ops_seen[held] = dict(path=str(op_path.relative_to(ROOT)), sha256=sha,
                              max_eig=float(op.evals.max()),
                              min_eig=float(op.evals.min()), n=int(op.evals.size))
        print(f"[{held}] operator {op_path.name} sha {sha[:16]} "
              f"eig [{op.evals.min():.3e}, {op.evals.max():.4f}] n={op.evals.size} "
              f"targets={len(fold['targets'])}", flush=True)
        # `pre` does not depend on gamma -- built once per fold, exactly as 204 does.
        pre, aud = weighted_baseline(srcs, fold["targets"], "equal",
                                     {d: 1.0 for d in DONORS})
        aud["max_abs_dev_from_162_pre"] = float(
            np.max(np.abs(pre - np.asarray(fold["pre"], dtype=np.float64))))
        # 204's own recorded value for this arm/fold. It is the float32 round-trip of the
        # cached `pre`, not a difference in the arithmetic, and it must reproduce EXACTLY.
        assert aud["max_abs_dev_from_162_pre"] == DEV_KA[held], (
            f"`equal` no longer reproduces 204's pre on {held}: "
            f"{aud['max_abs_dev_from_162_pre']!r} != {DEV_KA[held]!r}")
        for gamma in GRID:
            if (held, gamma) in done:
                continue
            t0 = time.time()
            e = score_at(fold, pre, op, cd4, fold["scored"], truth, gamma)
            rows.append(dict(fold=held, gamma=gamma, arm="equal",
                             n_targets=len(fold["targets"]), **e,
                             operator=ops_seen[held],
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=1))
            print(f"[{held}/g{gamma:.3f}] cos {e['a23_cos']:.7f} "
                  f"a15cos {e['a15_cos']:.7f} pds {e['a23_pds']:.6f} "
                  f"a15pds {e['a15_pds']:.6f} q_w {e['a15_sgn_w']:.6f} "
                  f"rr {e['rr']:.6f} f {e['cd4_solved_fraction']:.6f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        del op, srcs, cd4, pre
        import gc
        gc.collect()

    # ---------------------------------------------------------- known-answer adjudication
    print("\n===== KNOWN-ANSWER ROW: gamma 0.65 == 204's published `equal` row =====")
    ok = True
    for held, want in KA.items():
        got = [r for r in rows if r["fold"] == held and r["gamma"] == 0.65][0]
        for k, v in want.items():
            rel = abs(got[k] / v - 1.0)
            good = rel <= 1e-12
            ok &= good
            print(f"  {held:<15} {k:<10} {'PASS' if good else 'FAIL'}  "
                  f"{got[k]!r}  rel {rel:.2e}")
    with open(WORK / "plolo212_ka.json", "w") as fh:
        json.dump({"passed": bool(ok), "operators": ops_seen, "grid": GRID}, fh, indent=1)
    if not ok:
        raise SystemExit("PLOLO KNOWN-ANSWER ROW FAILED — not 204's instrument. STOP.")
    print("  ==> PASSED on all four folds.\n", flush=True)


if __name__ == "__main__":
    main()

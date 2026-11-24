"""216 Phase 0 — known-answer rows.  Abort on failure.

Rebuilds 162's PLOLO folds FROM SCRATCH in this experiment's own work dir (no cache is
copied in) and checks every number the mandate names:

  * PLOLO a15 `pds` 0.795784 and a23 0.789121 on the `k562_gwps` fold (204 §0 gates 2/2b)
  * the a15 -> a23 mass-fix price +0.004014 scaled, 4-fold pre-committed unweighted mean
  * 207's dose sweep at its two shipped points (f = 0.10, f = 0.34) on BOTH cosine
    functionals, because 207 reported `cos_agg` and 212 reported the macro `cos` and the
    mandate quotes one number from each
  * 212's forward +13.83% at f = 0.34 (macro `cos`, 4-fold, gamma 0.65)
  * 197's one-sided cosine floor identity on our own board rows
  * rule 11: path + sha256 for every operator this experiment can reach

    uv run python gate216.py
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

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, OFFICIAL, PROD_OPS, Operator, build_fold,  # noqa: E402
                   cd4_for, deliver, endpoints, solve_fraction)

ALPHA = 0.982
GAMMA = 0.65
LAM = 0.0
SPAN_PDS = 0.452

# ---------------------------------------------------------------- expected values
EXPECT = {
    "a15_pds_k562": 0.79578409,
    "a23_pds_k562": 0.78912131,
    "a23_q_k562": 0.50164311,
    "a23_qw_k562": 0.50487994,
    "massfix_scaled": 0.004014,
    "cos_agg_f0": 0.018091,
    "cos_agg_f010": 0.020086,
    "cos_agg_f034": 0.023698,
    "pds_f0": 0.735107,
    "pds_f010": 0.742112,
    "pds_f034": 0.740275,
    "cos_f010_212": 0.0218418,
    "cos_f034_212": 0.0248622,
    "fwd_212_pct": 13.83,
    "op_prod_sha": "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
    "op_prod_basis": [29394, 17942],
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(16 << 20), b""):
            h.update(block)
    return h.hexdigest()


def score_effects(final, fold, truth):
    """207's `score_effects`, kept verbatim in behaviour: a23 = mass fix, a15 = none."""
    out = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc = {}
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            e = endpoints(dl[:, fold["scored"]], truth)
            e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
            acc[ctx] = e
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in "ABC"]))
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    return out


def main() -> None:
    t0 = time.time()
    report = {"gates": [], "operators": [], "folds": {}}

    # ---------------------------------------------------------------- rule 11 census
    op_paths = sorted(set(
        [PROD_OPS / "whitener_union_pruned.npz"]
        + [Path(FOLDS[f]["op"]) for f in FOLDS]))
    for p in op_paths:
        with np.load(p) as z:
            evals = np.asarray(z["evals"])
            evecs = np.asarray(z["evecs"])
        report["operators"].append(dict(
            path=str(p.relative_to(ROOT)), sha256=sha256_file(p),
            n=int(evecs.shape[0]), max_eig=float(evals.max()), min_eig=float(evals.min())))
        print(f"[op] {p.relative_to(ROOT)}  sha {report['operators'][-1]['sha256'][:8]}  "
              f"n={evecs.shape[0]}  maxeig={evals.max():.6f}", flush=True)

    prod = [o for o in report["operators"] if o["path"].endswith("whitener_union_pruned.npz")][0]
    meta = json.loads((PROD_OPS / "whitener_union_pruned.json").read_text())
    ok_sha = prod["sha256"] == EXPECT["op_prod_sha"]
    ok_basis = list(meta["basis"]) == EXPECT["op_prod_basis"]
    report["gates"].append(dict(name="KA-0 production operator sha256+basis",
                                expected=f"{EXPECT['op_prod_sha'][:8]} / {EXPECT['op_prod_basis']}",
                                measured=f"{prod['sha256'][:8]} / {meta['basis']}",
                                passed=bool(ok_sha and ok_basis)))
    assert ok_sha and ok_basis, "production operator identity FAILED"

    # ------------------------------------------------------------------- fold build
    folds, truths, cd4s, whits = {}, {}, {}, {}
    for held in FOLDS:
        f = build_fold(held)
        folds[held] = f
        truths[held] = np.asarray(f["truth"], dtype=np.float64)[:, f["scored"]]
        cd4s[held] = cd4_for(f["targets"])
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        whits[held] = op.apply(np.asarray(f["pre"], dtype=np.float64) * ALPHA,
                               f["wcols"], GAMMA, LAM)
        del op
        report["folds"][held] = dict(
            n_targets=len(f["targets"]), audit=f["audit"],
            cache_sha256=sha256_file(plolo.WORK / f"fold-{held}.npz"))
        print(f"[fold] {held} built ({time.time()-t0:.0f}s)", flush=True)

    # compare the freshly built caches against 207's, array by array (zip mtime differs)
    ref_dir = ROOT / "experiments/207-cd4-dose/work/plolo"
    for held in FOLDS:
        mine = np.load(plolo.WORK / f"fold-{held}.npz", allow_pickle=True)
        ref = np.load(ref_dir / f"fold-{held}.npz", allow_pickle=True)
        worst = 0.0
        for k in ("pre", "truth", "b_ctrl", "gate", "scored", "wcols"):
            a, b = np.asarray(mine[k]), np.asarray(ref[k])
            assert a.shape == b.shape, (held, k, a.shape, b.shape)
            if a.dtype.kind in "fc":
                worst = max(worst, float(np.max(np.abs(a.astype(np.float64)
                                                       - b.astype(np.float64)))))
            else:
                assert np.array_equal(a, b), (held, k)
        report["gates"].append(dict(name=f"KA-1 fold {held} rebuilt == 207's cache",
                                    expected="max|diff| = 0", measured=f"{worst:.3e}",
                                    passed=bool(worst == 0.0)))
        print(f"[KA-1] fold {held}: max|diff| vs 207 = {worst:.3e}", flush=True)

    # ------------------------------------------------------------- the arms we need
    def run(held, f_val=None, fraction=None):
        fold, cd4, whit = folds[held], cd4s[held], whits[held]
        if f_val is not None:
            frac, per_unit = solve_fraction(whit, cd4, f_val)
        else:
            frac, per_unit = float(fraction or 0.0), None
        if frac > 0.0:
            final, _ = E062.add_post_whiten_residual(whit, cd4, frac, fold["wcols"])
        else:
            final = whit.copy()
        final[:, fold["panel_positions"]] = 0.0
        e = score_effects(final, fold, truths[held])
        e["fraction"] = float(frac)
        e["f_per_unit"] = float(per_unit) if per_unit else None
        return e

    arms = {}
    for held in FOLDS:
        for name, kw in (("f0", dict()), ("f010", dict(f_val=0.10)),
                         ("f034", dict(f_val=0.34)), ("a14", dict(fraction=0.05))):
            arms[(held, name)] = run(held, **kw)
            print(f"[arm] {held}/{name} pds {arms[(held,name)]['pds']:.6f} "
                  f"cos_agg {arms[(held,name)]['cos_agg']:.6f} "
                  f"cos {arms[(held,name)]['cos']:.7f} ({time.time()-t0:.0f}s)", flush=True)

    def gate(name, expected, measured, tol):
        ok = abs(measured - expected) <= tol
        report["gates"].append(dict(name=name, expected=expected, measured=measured,
                                    tol=tol, passed=bool(ok)))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: expected {expected} measured {measured}",
              flush=True)
        return ok

    k = arms[("k562_gwps", "f010")]
    allok = True
    allok &= gate("KA-2 a15 pds, k562_gwps fold", EXPECT["a15_pds_k562"],
                  round(k["a15_pds"], 8), 5e-6)
    allok &= gate("KA-3 a23 pds, k562_gwps fold", EXPECT["a23_pds_k562"],
                  round(k["a23_pds"], 8), 5e-6)
    allok &= gate("KA-3b a23 q_w, k562_gwps fold", EXPECT["a23_qw_k562"],
                  round(k["a23_sgn_w"], 8), 5e-6)

    massfix = float(np.mean([(arms[(h, "f010")]["a23_pds"] - arms[(h, "f010")]["a15_pds"])
                             for h in FOLDS])) / SPAN_PDS
    allok &= gate("KA-4 a15->a23 mass-fix price, scaled, 4-fold",
                  EXPECT["massfix_scaled"], round(massfix, 6), 5e-6)

    for tag, key, exp in (("f0", "cos_agg", "cos_agg_f0"), ("f010", "cos_agg", "cos_agg_f010"),
                          ("f034", "cos_agg", "cos_agg_f034"), ("f0", "pds", "pds_f0"),
                          ("f010", "pds", "pds_f010"), ("f034", "pds", "pds_f034")):
        v = float(np.mean([arms[(h, tag)][key] for h in FOLDS]))
        allok &= gate(f"KA-5 207 sweep {tag} {key}", EXPECT[exp], round(v, 6), 1.5e-5)

    cos010 = float(np.mean([arms[(h, "f010")]["cos"] for h in FOLDS]))
    cos034 = float(np.mean([arms[(h, "f034")]["cos"] for h in FOLDS]))
    allok &= gate("KA-6 212 macro cos at f=0.10", EXPECT["cos_f010_212"],
                  round(cos010, 7), 2e-6)
    allok &= gate("KA-6b 212 macro cos at f=0.34", EXPECT["cos_f034_212"],
                  round(cos034, 7), 2e-6)
    allok &= gate("KA-7 212 forward gain %", EXPECT["fwd_212_pct"],
                  round(100.0 * (cos034 / cos010 - 1.0), 2), 0.02)

    # 197's one-sided floor identity, on our own board rows
    floors = {}
    for name, raw_mse in (("a8", 1.3627), ("a14", 0.882775), ("a15/rc_hybrid", 0.8738045),
                          ("a23/rc_callset", 0.8738045)):
        floors[name] = (1.023611 - raw_mse) / 2.074495
    report["cosine_floors_197"] = floors
    allok &= gate("KA-8 197 floor at a15's raw mse", 0.07221,
                  round(floors["a15/rc_hybrid"], 5), 1e-5)

    report["arms"] = {f"{h}/{n}": v for (h, n), v in arms.items()}
    report["four_fold"] = {
        tag: {k2: float(np.mean([arms[(h, tag)][k2] for h in FOLDS]))
              for k2 in ("pds", "cos", "cos_agg", "a15_pds", "a23_pds", "sgn_w")}
        for tag in ("f0", "f010", "f034", "a14")}
    report["all_gates_passed"] = bool(allok and all(g["passed"] for g in report["gates"]))
    report["seconds"] = round(time.time() - t0, 1)
    (WORK / "gate216.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({"all_gates_passed": report["all_gates_passed"],
                      "seconds": report["seconds"]}, indent=2))
    if not report["all_gates_passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()

"""190 Phase 4 (cont.) — price the regulon DIRECTION correction on PLOLO `pds`.

162 is the only instrument in this repo validated against a board `pds` point (a23: sign
correct, magnitude over-predicted 5.5x in the challenge frame).  It cannot produce `reach`,
`fid` or `jac` at all (162 HANDOFF trap 1), so it is used here for `pds` only; the other
members are priced on h1-test in `comp190.py`.

Arm: a15's exact shipped config, then
        final' = final + c * ||final||_2(row) * unit(CollecTRI correction)
with the correction built exactly as `nets190.ct_prop` builds it, on the fold's 300-panel
axis, from the fold's OWN whitened effect (so the propagated TF activity is the arm's own
prediction, not h1-test's).  c is swept; tydymy published c = 0.25.

KNOWN-ANSWER ROW: the c = 0 arm must reproduce 162's recorded a15 `pds` for this fold, and
`identity/intended` vs `library/delivered` must keep 156's sign.  Aborts otherwise.

Memory: one 17,942^2 float64 operator (2.6 GB) plus the cached fold.  Guarded at 9 GB.
"""
from __future__ import annotations

import json
import os
import resource
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(ROOT, "experiments", "162-retrodicting-validator"))
sys.path.insert(0, HERE)

HELD = "orion_hct116"
CS = (0.0, 0.05, 0.10, 0.25, 0.50)


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 3)


def guard(tag: str) -> None:
    r = rss_gib()
    print(f"   [rss {r:.2f} GiB after {tag}]", flush=True)
    if r > 9.0:
        raise SystemExit("ABORT: peak RSS over 9 GiB, the machine is shared")


def main() -> None:
    import plolo
    from plolo import Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction
    from run_arms import effects_for, frames_for
    import prior190

    cfg = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, mass=None, op="union_pruned")
    fold = build_fold(HELD)
    guard("build_fold")
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    cd4 = cd4_for(fold["targets"])
    guard("cd4")
    op = Operator(plolo.FOLDS[HELD]["op"])
    guard("operator")
    print(f"   operator {plolo.FOLDS[HELD]['op'].name} maxeig {op.evals.max():.4f}")

    final, caud = effects_for(fold, cfg, op, cd4)
    del op
    targets = np.asarray([str(t) for t in fold["targets"]])
    genes = np.asarray(plolo.OFFICIAL)

    # --- the correction, built on this fold's own effect ------------------------------
    d = prior190.load_collectri(weighted=True)
    gi = {g: j for j, g in enumerate(genes)}
    d = d[d["gene"].isin(gi)]                    # only edges landing on the official axis
    tfs = sorted(set(d["tf"]) & set(gi))
    fpos = {f: i for i, f in enumerate(tfs)}
    sub = d[d["tf"].isin(fpos)]
    R = np.zeros((len(tfs), len(genes)), dtype=np.float64)
    R[[fpos[f] for f in sub["tf"]], [gi[g] for g in sub["gene"]]] = (
        sub["mor"].to_numpy() * sub["w"].to_numpy())
    tf_col = np.asarray([gi[f] for f in tfs])
    direct = [t for t in targets if t in set(d["tf"])]
    print(f"   CollecTRI on the 300-panel PLOLO fold: {len(tfs)} TFs measurable on the "
          f"official axis, {len(direct)}/{len(targets)} fold targets are themselves TFs")

    def correction(topk):
        A = final[:, tf_col].copy()
        sd = A.std(1, keepdims=True)
        A = A / np.where(sd == 0, 1.0, sd)
        if topk is not None:
            keep = np.zeros_like(A, dtype=bool)
            np.put_along_axis(keep, np.argsort(-np.abs(A), axis=1)[:, :topk], True, axis=1)
            A = np.where(keep, A, 0.0)
        N = A @ R
        N[:, fold["panel_positions"]] = 0.0     # production zeroes the panel rows
        return N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)

    def direct_correction():
        N = np.zeros_like(final)
        ti = {t: i for i, t in enumerate(targets)}
        hit = d[d["tf"].isin(ti) & d["gene"].isin(gi)]
        N[[ti[t] for t in hit["tf"]], [gi[g] for g in hit["gene"]]] = (
            -hit["mor"].to_numpy() * hit["w"].to_numpy())
        N[:, fold["panel_positions"]] = 0.0
        return N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)

    frames = frames_for(fold)
    n0 = np.linalg.norm(final, axis=1, keepdims=True)
    rows = []
    base_int = None
    print(f"\n   {'correction':16} {'c':>6} {'intended pds':>13} {'chalA pds':>11} "
          f"{'library pds':>12} {'sse':>9}")
    for cname, U in (("ct-direct", direct_correction()), ("ct-prop-k25", correction(25)),
                     ("ct-prop-all", correction(None))):
        for c in CS:
            if c == 0.0 and base_int is not None:
                continue
            P = final + c * n0 * U
            e = endpoints(P[:, fold["scored"]], truth)
            row = dict(corr=cname, c=c, space="intended", **e)
            chal = lib = None
            for fname, b, gate in frames:
                if fname not in ("library", "chalA"):
                    continue
                mat = deliver(P, b, gate, cfg["mass"])
                ed = endpoints(mat[:, fold["scored"]], truth)
                rows.append(dict(corr=cname, c=c, space=f"{fname}/delivered", **ed))
                if fname == "chalA":
                    chal = ed["pds"]
                else:
                    lib = ed["pds"]
            rows.append(row)
            if c == 0.0:
                base_int = (e["pds"], chal, lib)
            print(f"   {cname if c else 'BASE (a15, c=0)':16} {c:6.2f} {e['pds']:13.6f} "
                  f"{chal:11.6f} {lib:12.6f} {e['sse']:9.6f}")
    guard("scoring")

    print(f"\n   KNOWN-ANSWER ROW (162 HANDOFF §2): library/delivered a15 on {HELD} = 0.789536")
    lib0 = [r for r in rows if r["corr"] == "ct-direct" and r["c"] == 0.0
            and r["space"] == "library/delivered"][0]["pds"]
    ok = abs(lib0 - 0.789536) < 5e-5
    print(f"   {'OK ' if ok else 'FAIL'} measured {lib0:.6f}")
    if not ok:
        raise SystemExit("known-answer row failed — the PLOLO numbers are not trustworthy")

    print(f"\n   deltas vs a15 (162 gate: predicted overall = PLOLO d_pds / 6 * 1.011)")
    print(f"   {'correction':16} {'c':>6} {'d pds intended':>15} {'d pds chalA':>13} "
          f"{'-> overall':>12}")
    for cname in ("ct-direct", "ct-prop-k25", "ct-prop-all"):
        for c in CS[1:]:
            i = [r for r in rows if r["corr"] == cname and r["c"] == c
                 and r["space"] == "intended"][0]["pds"]
            a = [r for r in rows if r["corr"] == cname and r["c"] == c
                 and r["space"] == "chalA/delivered"][0]["pds"]
            print(f"   {cname:16} {c:6.2f} {i - base_int[0]:+15.6f} {a - base_int[1]:+13.6f} "
                  f"{(a - base_int[1]) / 6 * 1.011:+12.6f}")

    with open(os.path.join(WORK, "plolo190.json"), "w") as f:
        json.dump({"fold": HELD, "base": base_int, "rows": rows,
                   "n_targets": len(targets), "n_tf_targets": len(direct),
                   "peak_rss_gib": rss_gib()}, f, indent=2, default=float)
    print(f"\n   peak RSS {rss_gib():.2f} GiB")
    print(f"wrote {WORK}/plolo190.json")


if __name__ == "__main__":
    main()

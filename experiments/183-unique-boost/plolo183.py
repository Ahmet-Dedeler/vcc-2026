"""183 — price the unique-component boost's `pds` on PLOLO, in the challenge frame.

Rule 13 names PLOLO (experiment 162) as the only locally priceable `pds` instrument:
four donor folds over the real 300-target panel, rendered into the CHALLENGE control
frames A/B/C, scored against the held library's own deep pseudobulk.  It retrodicts board
`pds` 6/6 on sign with a leave-one-arm-out RMS of 0.001494 scaled.

162's and 169's files are NOT modified; their cached folds are reused.

THREE PLACEMENTS, because 156 §6 is about exactly this
------------------------------------------------------
`deliver()` = `mass_consistent_effect` -> `profile_from_effect` -> the scorer's read-back.
Under a23/production's `minnorm` mass fix that composition is an identity to 1e-12, so a
geometry transform can sit on either side of it.  156 measured that a PRE-render centroid
move is silently undone by the renderer when the mass fix is OFF; the honest thing is to
measure all three placements rather than assume the mass fix rescues it.

  `post`  T(deliver(F))                — the transform in the space the scorer reads.
                                         Pure geometry; NOT shippable on its own, because
                                         T(D) no longer satisfies the mass constraint.
  `re`    deliver(T(deliver(F)))       — the same, pushed back through the renderer.
  `pre`   deliver(T(F))                — what production would actually do: transform the
                                         effect matrix, render once.  THE SHIPPABLE ARM.

If `pre` and `post` disagree, the renderer is eating the transform and the lever is dead
for the same reason Duy's cross-target centring was.

KNOWN-ANSWER ROWS (rule 10's corollary)
---------------------------------------
  1. alpha = 0.00 must reproduce 162's own a23 rows fold-by-fold from `work/arms.json`,
     in all three placements.  ABORT otherwise.
  2. a23 - a15 in the challenge frame must reproduce 162's published mass-fix price
     (+0.004014 scaled = +0.001815 raw pds).
  3. ||h'|| == ||h|| out of the transform, to 1e-12, at every alpha.

No cell-eval2, no build, no submission, no write outside this directory.
"""
from __future__ import annotations

import json
import os
import resource
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(REPO, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
sys.path.insert(0, HERE)

import plolo                                                          # noqa: E402
from plolo import (FOLDS, PROD_OPS, Operator, build_fold, cd4_for,     # noqa: E402
                   deliver, endpoints)
from run_arms import CONFIGS, effects_for, frames_for                  # noqa: E402
from uboost import norm_check, shared_energy, unique_boost             # noqa: E402

ALPHAS = (0.00, 0.25, 0.50, 0.75, 1.00)
PLACEMENTS = ("post", "re", "pre")
BASE = {"a23": CONFIGS["a23"], "a15": CONFIGS["a15"]}
PDS_SPAN = 0.95228 - 0.50020          # debench ANCHOR, the board-fitted pds anchors


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 3)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out_path = os.path.join(WORK, "plolo183.json")
    rows = json.loads(open(out_path).read()) if os.path.exists(out_path) else []
    done = {(r["fold"], r["arm"], r["alpha"], r["space"]) for r in rows}
    worst_norm = 0.0

    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = [f for f in frames_for(fold) if f[0] != "library"]   # challenge frames only
        op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                      if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        print(f"[{held}] operator maxeig {op.evals.max():.4f}, "
              f"{len(fold['targets'])} targets, {len(fold['scored'])} scored genes, "
              f"RSS {rss_gib():.2f} GiB", flush=True)

        for arm, cfg in BASE.items():
            final, _ = effects_for(fold, cfg, op, cd4)
            for fname, b, gate in frames:
                D0 = deliver(final, b, gate, cfg["mass"])             # the delivered object
                for alpha in ALPHAS:
                    if arm == "a15" and alpha != 0.00:
                        continue                  # a15 is only the mass-fix control row
                    variants = {}
                    Dpost = unique_boost(D0, alpha)
                    worst_norm = max(worst_norm, norm_check(D0, Dpost))
                    variants["post"] = Dpost
                    variants["re"] = deliver(Dpost, b, gate, cfg["mass"])
                    Fp = unique_boost(final, alpha)
                    worst_norm = max(worst_norm, norm_check(final, Fp))
                    variants["pre"] = deliver(Fp, b, gate, cfg["mass"])
                    for place, M in variants.items():
                        sp = f"{fname}/{place}"
                        if (held, arm, alpha, sp) in done:
                            continue
                        e = endpoints(M[:, fold["scored"]], truth)
                        g = shared_energy(M[:, fold["scored"]])
                        # how much of the transform the renderer kept
                        rowcos = float(np.mean(
                            (M * Dpost).sum(1)
                            / np.maximum(np.linalg.norm(M, axis=1)
                                         * np.linalg.norm(Dpost, axis=1), 1e-300)))
                        rows.append(dict(fold=held, arm=arm, alpha=alpha, space=sp,
                                         cos_to_post=rowcos, **e, **g))
                        print(f"  {held:<14} {arm} a={alpha:.2f} {sp:<12} "
                              f"pds {e['pds']:.6f}  cos {e['cos']:.5f}  "
                              f"rr {g['mean_pairwise_cos']:+.5f}  "
                              f"cE {g['centroid_energy_frac']:.4f}  "
                              f"|h| {g['median_row_norm']:.5f}  sse {e['sse']:.6f}",
                              flush=True)
                    del variants
                open(out_path, "w").write(json.dumps(rows, indent=1))
                del D0
        del op

    print(f"\nworst |‖h'‖/‖h‖ − 1| over every transform call: {worst_norm:.3e}", flush=True)
    if worst_norm > 1e-12:
        sys.exit("ABORT: the transform is not norm-preserving")
    print(f"peak RSS {rss_gib():.2f} GiB")
    report(rows, worst_norm)


# ------------------------------------------------------------------------------ report
def report(rows, worst_norm) -> None:
    R = {(r["fold"], r["arm"], r["alpha"], r["space"]): r for r in rows}
    folds = list(FOLDS)

    print(f"\n{'=' * 110}\nKNOWN-ANSWER GATE 1: alpha=0 must reproduce 162's own a23 rows")
    ref = json.loads(open(os.path.join(P162, "work/arms.json")).read())
    A162 = {(r["fold"], r["arm"], r["space"]): r for r in ref}
    worst = 0.0
    for f in folds:
        for c in "ABC":
            a = A162.get((f, "a23", f"chal{c}/delivered"))
            for place in PLACEMENTS:
                h = R.get((f, "a23", 0.00, f"chal{c}/{place}"))
                if a is None or h is None:
                    continue
                d = abs(a["pds"] - h["pds"])
                worst = max(worst, d)
    print(f"  worst |d pds| across 4 folds x 3 frames x 3 placements = {worst:.3e}")
    if worst > 1e-9:
        print("  ABORT" if worst > 1e-2 else "  (inside 162 §6's mass-fix fixed-point floor)")
    if worst > 1e-2:
        sys.exit("ABORT: alpha=0 does not reproduce the incumbent")

    def chal(f, arm, alpha, place):
        return float(np.mean([R[(f, arm, alpha, f"chal{c}/{place}")]["pds"] for c in "ABC"]))

    print(f"\n{'=' * 110}\nKNOWN-ANSWER GATE 2: the mass fix's own PLOLO price "
          f"(162 published +0.004014 scaled)")
    d = [chal(f, "a23", 0.00, "pre") - chal(f, "a15", 0.00, "pre") for f in folds]
    print("  per fold raw: " + "  ".join(f"{f}={x:+.6f}" for f, x in zip(folds, d)))
    print(f"  4-fold mean raw {np.mean(d):+.6f}  ->  scaled {np.mean(d) / PDS_SPAN:+.6f}")

    print(f"\n{'=' * 110}\nKNOWN-ANSWER GATE 3: norm preservation  worst |ratio-1| "
          f"= {worst_norm:.3e}  (want < 1e-12)")

    for place in PLACEMENTS:
        print(f"\n{'=' * 110}\nUNIQUE-COMPONENT BOOST on PLOLO, challenge frame, "
              f"placement = {place}")
        print(f"  {'alpha':>6} " + " ".join(f"{f:>14}" for f in folds) +
              f" {'MEAN raw':>10} {'sd':>9} {'MEAN scaled':>12}")
        base = {f: chal(f, "a23", 0.00, place) for f in folds}
        for alpha in ALPHAS:
            v = [chal(f, "a23", alpha, place) - base[f] for f in folds]
            print(f"  {alpha:>6.2f} " + " ".join(f"{x:+14.6f}" for x in v) +
                  f" {np.mean(v):+10.6f} {np.std(v):9.6f} {np.mean(v) / PDS_SPAN:+12.6f}")

    print(f"\n{'=' * 110}\nGEOMETRY, 4-fold x 3-frame mean (placement = pre, the shippable one)")
    print(f"  {'alpha':>6} {'pds':>10} {'d pds':>10} {'cos':>10} {'d cos %':>9} "
          f"{'pair cos':>10} {'centroid E':>11} {'med |h|':>9} {'sse*':>10} "
          f"{'cos to post':>12}")
    for alpha in ALPHAS:
        def m(key, place="pre"):
            return float(np.mean([np.mean([R[(f, "a23", alpha, f"chal{c}/{place}")][key]
                                           for c in "ABC"]) for f in folds]))

        def m0(key, place="pre"):
            return float(np.mean([np.mean([R[(f, "a23", 0.00, f"chal{c}/{place}")][key]
                                           for c in "ABC"]) for f in folds]))
        dc = 100.0 * (m("cos") / m0("cos") - 1.0)
        print(f"  {alpha:>6.2f} {m('pds'):10.6f} {m('pds') - m0('pds'):+10.6f} "
              f"{m('cos'):10.5f} {dc:+9.2f} {m('mean_pairwise_cos'):+10.5f} "
              f"{m('centroid_energy_frac'):11.5f} {m('median_row_norm'):9.5f} "
              f"{m('sse'):10.6f} {m('cos_to_post'):12.6f}")

    print(f"\n{'=' * 110}\nPLACEMENT DISAGREEMENT — does the renderer eat the transform?")
    print(f"  {'alpha':>6} {'post':>10} {'re':>10} {'pre':>10} "
          f"{'|re-post|':>10} {'|pre-post|':>11}")
    for alpha in ALPHAS:
        vals = {}
        for place in PLACEMENTS:
            base = np.mean([chal(f, "a23", 0.00, place) for f in folds])
            vals[place] = np.mean([chal(f, "a23", alpha, place) for f in folds]) - base
        print(f"  {alpha:>6.2f} {vals['post']:+10.6f} {vals['re']:+10.6f} "
              f"{vals['pre']:+10.6f} {abs(vals['re'] - vals['post']):10.6f} "
              f"{abs(vals['pre'] - vals['post']):11.6f}")


if __name__ == "__main__":
    sys.exit(main())

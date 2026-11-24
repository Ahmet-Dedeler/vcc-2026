"""188 — predict `nmae` FROM THE IDENTITY, on every surface, and decompose the disagreement.

169 §5's identity, gated there against `debench.Surface.nmae()` to 0.0e+00 / 1.1e-16:

    nmae_p = 1 + ( U_p - 2 M_p ) / D_p
    U_p = sum_g |h_g|                            the L1 BUDGET we spend on the pool
    M_p = sum_g min(|h_g|,|r_g|) * 1[sign ok]    the refund, capped at |r| per gene
    D_p = sum_g |r_g|                            reference-only, fixed

so a shape change is priced by exactly two numbers, `dU` and `dM`, and the marginal refund rate
`e = dM/dU`:  **nmae improves iff `e > 1/2`.**  That threshold is not fitted (169 §5).

This script computes `U`, `M`, `D` per perturbation on THREE surfaces -- h1-test's emitted DE
tables, the four PLOLO folds, and (for the pool composition only) the challenge panel -- and
reports `e` for each arm.  If the surfaces disagree on `nmae`, the disagreement must live in
`e`, and `e` factorises into sign accuracy and saturation.  That is the comparison AGENTS.md
rule 12 and this experiment's mandate ask for.

KNOWN-ANSWER ROW: (I) must reproduce `debench.Surface.nmae()` on every h1-test table to <1e-12.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
sys.path.insert(0, HERE)
sys.path.insert(0, E169)
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface                                    # noqa: E402
from budget188 import FOLDS, N_DE, deliver, load, nmae_star     # noqa: E402
from mapfit import to_lfc                                       # noqa: E402
from sharp import sharpen, sharpen_lfc                           # noqa: E402

PART = os.environ.get("PART", "all")


def umd(P, R, pool):
    """`U`, `M`, `D` per row over `pool`, and the identity's own nmae."""
    U, M, D = [], [], []
    for i in range(P.shape[0]):
        s = np.flatnonzero(pool[i])
        if s.size < 10:
            continue
        h, r = np.nan_to_num(P[i, s]), R[i, s]
        d = np.abs(r).sum()
        if not np.isfinite(d) or d == 0:
            continue
        m = np.where(np.sign(h) == np.sign(r), np.minimum(np.abs(h), np.abs(r)), 0.0)
        U.append(np.abs(h).sum()); M.append(m.sum()); D.append(d)
    U, M, D = np.array(U), np.array(M), np.array(D)
    return U, M, D, float(np.mean(1 + (U - 2 * M) / D))


def report(tag, rows, base):
    print(f"\n  {tag}")
    print(f"    {'arm':<14}{'nmae (I)':>10}{'d nmae':>10}{'k = U/D':>9}{'dU/D':>9}"
          f"{'dM/D':>9}{'e = dM/dU':>11}{'L1 over k=1':>13}")
    for name, r in rows.items():
        b = rows[base]
        du = float(np.mean((r["U"] - b["U"]) / r["D"]))
        dm = float(np.mean((r["M"] - b["M"]) / r["D"]))
        e = dm / du if abs(du) > 1e-12 else float("nan")
        print(f"    {name:<14}{r['nmae']:>10.5f}{r['nmae']-b['nmae']:>+10.5f}"
              f"{float(np.mean(r['U']/r['D'])):>9.4f}{du:>+9.5f}{dm:>+9.5f}{e:>11.4f}"
              f"{r['over']:>13.4f}")


def main() -> None:
    out = {}
    if PART in ("all", "h1"):
        print("=" * 106)
        print("h1-test — emitted DE tables, the reference's own BH<0.05 pool (the surface "
              "validated for `reach`)")
        tabs = {}
        for p in sorted(glob.glob(os.path.join(E169, "de169.*.npz"))):
            tabs[os.path.basename(p)[len("de169."):-4]] = p
        for p in sorted(glob.glob(os.path.join(HERE, "de188.*.npz"))):
            tabs[os.path.basename(p)[len("de188."):-4]] = p
        want = [a for a in ("q1.00_mf", "q1.15_mf", "q1.15L_mf", "q1.25_mf", "q1.25L_mf",
                            "KA_q1.15_mf", "g10_q1.00_s1", "g10_q1.15_s1", "g20_q1.15_s1",
                            "map_med_s1", "map_bud_s1") if a in tabs]
        rows, gate_ok = {}, True
        for a in want:
            s = Surface(tabs[a])
            pool = s.real_sig & np.isfinite(s.real_lfc)
            U, M, D, nm = umd(s.pred_lfc, s.real_lfc, pool)
            ref = s.nmae()
            if abs(nm - ref) > 1e-12:
                gate_ok = False
            h = np.nan_to_num(s.pred_lfc)
            ov = float(np.abs(h)[pool & (np.abs(h) > np.abs(s.real_lfc))].sum() /
                       max(np.abs(h)[pool].sum(), 1e-12))
            rows[a] = dict(U=U, M=M, D=D, nmae=nm, over=ov, debench=ref,
                           npool=float(np.median(pool.sum(1))))
            print(f"    {a:<14} identity {nm:.9f}  debench {ref:.9f}  "
                  f"|d| {abs(nm-ref):.1e}  median pool {rows[a]['npool']:.0f}")
        print(f"  KNOWN-ANSWER GATE (I) vs debench: {'PASS' if gate_ok else 'FAIL'}")
        if not gate_ok:
            sys.exit("ABORT: identity (I) does not reproduce debench.Surface.nmae()")
        report("h-space and lfc-space arms, baseline `q1.00_mf`", rows, "q1.00_mf")
        if "g10_q1.00_s1" in rows:
            report("GATE-WIDTH arms, baseline `g10_q1.00_s1` (identity transform, same seed)",
                   {k: v for k, v in rows.items()
                    if k in ("g10_q1.00_s1", "g10_q1.15_s1", "g20_q1.15_s1", "q1.15_mf",
                             "map_med_s1", "map_bud_s1")}, "g10_q1.00_s1")
        out["h1"] = {k: dict(nmae=v["nmae"], over=v["over"], npool=v["npool"],
                             kUD=float(np.mean(v["U"] / v["D"])))
                     for k, v in rows.items()}

    if PART in ("all", "plolo"):
        print("\n" + "=" * 106)
        print("PLOLO — four donor folds, truth = the held library's own deep pseudobulk, "
              "pool = its top-2000 |lfc|")
        z = np.load(os.path.join(WORK, "chalframes.npz"))
        cg = z["gA"] & z["gB"] & z["gC"]
        bA, gA = z["bA"], z["gA"]
        acc = {}
        for held in FOLDS:
            F = load(held)
            mk = np.zeros(cg.size, bool)
            mk[F["scored"]] = True
            c = np.flatnonzero(mk & cg)
            LR = to_lfc(F["truth"], F["b"], gA)[:, c]
            pool = np.zeros(LR.shape, bool)
            k = min(N_DE, c.size)
            np.put_along_axis(pool, np.argpartition(-np.abs(LR), k - 1, axis=1)[:, :k],
                              True, 1)
            P0 = F["pred"]
            for name, Q in (("q1.00", np.where(gA[None, :], P0, 0.0)),
                            ("q1.15", sharpen(P0, 1.15, gA, mode="row")),
                            ("q1.15L", sharpen_lfc(P0, 1.15, bA, gA, mode="row")),
                            ("q1.25", sharpen(P0, 1.25, gA, mode="row"))):
                LP = to_lfc(deliver(Q, bA, gA), bA, gA)[:, c]
                U, M, D, nm = umd(LP, LR, pool)
                ov = float(np.abs(LP)[pool & (np.abs(LP) > np.abs(LR))].sum() /
                           max(np.abs(LP)[pool].sum(), 1e-12))
                acc.setdefault(name, []).append(dict(U=U, M=M, D=D, nmae=nm, over=ov))
            print(f"    [{held}] done", flush=True)
        rows = {n: dict(U=np.concatenate([a["U"] for a in v]),
                        M=np.concatenate([a["M"] for a in v]),
                        D=np.concatenate([a["D"] for a in v]),
                        nmae=float(np.mean([a["nmae"] for a in v])),
                        over=float(np.mean([a["over"] for a in v]))) for n, v in acc.items()}
        report("PLOLO folds pooled, baseline `q1.00`", rows, "q1.00")
        out["plolo"] = {k: dict(nmae=v["nmae"], over=v["over"],
                                kUD=float(np.mean(v["U"] / v["D"]))) for k, v in rows.items()}

    tag = "" if PART == "all" else f"-{PART}"
    with open(os.path.join(WORK, f"identity188{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/identity188{tag}.json")


if __name__ == "__main__":
    sys.exit(main())

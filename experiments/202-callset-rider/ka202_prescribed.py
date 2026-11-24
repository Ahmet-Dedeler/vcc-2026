"""202 known-answer row #2 — reproduce 192 §5.3's h1-test row for THIS arm, verbatim.

192's protected design was priced by *prescribing* a call set on 185's frozen h1-test DE
tables (no emission).  The generating code was run inline as a heredoc and is not on disk;
it was recovered verbatim from `experiments/192-predicted-callset/claude-run.jsonl` line 802
and is reproduced below with the loop bodies unchanged.  If a single number of the (K, M)
grid moves, my reading of 192's arm is wrong and the build must not proceed.

The mandate's row is  fid +0.0010  jac +0.0129  reach +0.0112 raw, marginal to `rc_hybrid`,
i.e. `work/protect192.log`'s f = 0.50 block at K = 2,600 / M = 400:

    baseline (rc_hybrid @ n_pred 4,474)  fid 0.50771  jac 0.04058  reach 0.26126
    K=2600 M=400                         fid 0.50868  jac 0.05346  reach 0.27249

This script also reports the arm's CALL-SET COMPOSITION on h1-test (falsifier 3's h1-test
twin), which 192 never printed: the size of each half of the union, how many genes they
share, and the overlap with `rc_hybrid`'s own call set.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

R = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(R, "experiments/202-callset-rider")
sys.path.insert(0, f"{R}/experiments/153-de-bench")
sys.path.insert(0, f"{R}/experiments/192-predicted-callset")
from debench import Surface                                              # noqa: E402
from frontier192 import A15, INV, PowerSurface, order_key, topk_call     # noqa: E402

# The arm, exactly as preregistered.
K_ARM, M_ARM = 2600, 400
# 192's own numbers for the row the mandate names, and the baseline it is marginal to.
EXPECT = dict(fid=0.50868, jac=0.05346, reach=0.27249)
EXPECT_BASE = dict(fid=0.50771, jac=0.04058, reach=0.26126)
EXPECT_DELTA = dict(fid=+0.00097, jac=+0.01288, reach=+0.01123)


def main() -> None:
    E = np.load(f"{R}/experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(A15)
    z = np.load(A15, allow_pickle=True)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)
    absH = np.abs(E["H"][np.ix_(ti, gi)].astype(np.float64))
    T = absH.shape[0]
    inc = -np.log(np.maximum(z["pred_p"], 1e-300))
    ckey = np.broadcast_to(cpm[None, :], absH.shape).copy()
    print(f"h1-test surface: {T} perturbations x {absH.shape[1]} gated genes")

    out = {}
    for f in (0.35, 0.4, 0.5, 0.65):
        sf = PowerSurface(A15, f)
        b = topk_call(inc, np.full(T, 4474))
        a0 = dict(fid=sf.fid(call=b), jac=sf.jac(call=b),
                  reach=sf.reach(key=order_key(b, absH)))
        print(f"\nf = {f}   baseline (rc_hybrid @ n_pred 4,474): fid {a0['fid']:.5f} "
              f"jac {a0['jac']:.5f} reach {a0['reach']:.5f}")
        print(f"  {'K':>6}{'M (|H| protected)':>19}{'fid':>9}{'jac':>9}{'reach':>9}"
              f"{'SUM @1.290':>12}{'SUM @0.4277':>13}")
        for K in (2200, 2600, 3400, 4474):
            for M in (0, 400, 800, 1600):
                if M >= K:
                    continue
                head = topk_call(absH, np.full(T, M)) if M else np.zeros(absH.shape, bool)
                rest = topk_call(np.where(head, -np.inf, ckey), np.full(T, K - M))
                c = head | rest
                v = dict(fid=sf.fid(call=c), jac=sf.jac(call=c),
                         reach=sf.reach(key=order_key(c, absH)))
                dfid = (v["fid"] - a0["fid"]) * INV["fid"] / 6
                djac = (v["jac"] - a0["jac"]) * INV["jac"] / 6
                dr = (v["reach"] - a0["reach"]) * INV["reach"] / 6
                out[f"{f}|{K}|{M}"] = dict(raw=v, base=a0, sum_129=dfid + djac + dr * 1.290,
                                           sum_043=dfid + djac + dr * 0.4277)
                print(f"  {K:>6}{M:>19}{v['fid']:>9.5f}{v['jac']:>9.5f}{v['reach']:>9.5f}"
                      f"{dfid + djac + dr * 1.290:>+12.5f}"
                      f"{dfid + djac + dr * 0.4277:>+13.5f}")

    # ---------------------------------------------------------------- the known-answer gate
    row = out[f"0.5|{K_ARM}|{M_ARM}"]
    print("\n" + "=" * 78)
    print(f"KNOWN-ANSWER ROW #2 — 192 §5.3, f = 0.50, K = {K_ARM}, M = {M_ARM}")
    ok = True
    for m in ("fid", "jac", "reach"):
        got, want = row["raw"][m], EXPECT[m]
        gotb, wantb = row["base"][m], EXPECT_BASE[m]
        d, wd = got - gotb, EXPECT_DELTA[m]
        good = abs(got - want) < 5e-6 and abs(gotb - wantb) < 5e-6
        ok &= good
        print(f"  {m:<6} base {gotb:.5f} (192: {wantb:.5f})   arm {got:.5f} "
              f"(192: {want:.5f})   delta {d:+.5f} (192: {wd:+.5f})   "
              f"{'PASS' if good else 'FAIL'}")
    print(f"  KNOWN-ANSWER ROW #2: {'PASS' if ok else 'FAIL'}")

    # ---------------------------------------------------------------- composition, h1-test
    head = topk_call(absH, np.full(T, M_ARM))
    rest = topk_call(np.where(head, -np.inf, ckey), np.full(T, K_ARM - M_ARM))
    arm = head | rest
    # the unrestricted halves, i.e. the two sets the arm's logic says do different jobs
    raw_cpm = topk_call(ckey, np.full(T, K_ARM - M_ARM))
    inc_call = topk_call(inc, np.full(T, 4474))          # `rc_hybrid`'s own call set
    print("\nCALL-SET COMPOSITION on h1-test (mean over the 85 perturbations)")
    print(f"  |arm|                              {arm.sum(1).mean():8.1f}")
    print(f"  |top {M_ARM} by |H||                    {head.sum(1).mean():8.1f}")
    print(f"  |top {K_ARM - M_ARM} by cpm| (unrestricted)     {raw_cpm.sum(1).mean():8.1f}")
    print(f"  |halves shared| (|H| head n cpm)   {(head & raw_cpm).sum(1).mean():8.1f}")
    print(f"  |arm n rc_hybrid call set|         {(arm & inc_call).sum(1).mean():8.1f}")
    print(f"  jaccard(arm, rc_hybrid call set)   "
          f"{((arm & inc_call).sum(1) / (arm | inc_call).sum(1)).mean():8.4f}")

    with open(os.path.join(HERE, "work/ka202_prescribed.json"), "w") as fh:
        json.dump(dict(grid=out, known_answer_pass=bool(ok),
                       composition=dict(
                           arm=float(arm.sum(1).mean()),
                           head=float(head.sum(1).mean()),
                           cpm_raw=float(raw_cpm.sum(1).mean()),
                           shared=float((head & raw_cpm).sum(1).mean()),
                           overlap_rc=float((arm & inc_call).sum(1).mean()),
                       )), fh, indent=1, default=float)
    if not ok:
        raise SystemExit("known-answer row #2 FAILED — abort")


if __name__ == "__main__":
    main()

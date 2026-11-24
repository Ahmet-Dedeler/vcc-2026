"""224 — the decisive disambiguation: is the PANEL negative substance, or power?

`ladder224.py` on P4/P5/P6 fits the per-gene gain from 73-120 panel targets and reads
-0.04% / +1.77% / +0.45%, against +5.16% on the 968-target broad surface.  Two readings:

  (a) SUBSTANCE — the 300-gene panel deliberately excludes common-essential genes (099 §1.1),
      so the per-gene cross-line reliability ranking learned on the broad surface simply does
      not hold for panel responses;
  (b) POWER — 7,312 per-gene gains cannot be estimated from ~150 rows, and the panel surfaces
      are too thin to fit anything.

This script separates them by fitting the gain on **thousands of NON-panel targets** and
evaluating it on **panel targets of a held line**.  Three of the fit-pool lines
(`rpe1`, `nadig_jurkat`, `nadig_hepg2`) carry **0 of the 300 panel targets**, so that half of
the fit is provably panel-leak-free.

If the cross-fitted gain still reads ~0 on the panel, the answer is (a) and the finding does
not reach production.  If it recovers +5%, the answer is (b) and it does.

    uv run python crosspanel224.py
"""
from __future__ import annotations

import itertools
import json

import numpy as np

from core224 import (DIAG_SHRINK, WORK, decoy_cos, diag_fit, donor_mean,  # noqa: E402
                     inner_pairs, macro_cos, per_target_cos, pds, unit)
from substrate224 import LINES, PANEL_SET, ROOT, load_effects

SEED = 20260904
# lines that carry enough of the 300-gene panel to be a held line or a donor there
PANEL_LINES = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
CD4 = dict(path="data/external/cd4t/cloud-compact/cd4_compact/cd4_rest", group="CD4T")
# lines used ONLY to fit the gain.  rpe1 / nadig_jurkat / nadig_hepg2 carry 0/300 panel
# targets, so they cannot leak the panel into the gain at all.
FIT_LINES = ["k562_gwps", "orion_hct116", "orion_hek293t", "rpe1", "nadig_jurkat",
             "nadig_hepg2"]
# the strictest possible fit pool: every line here carries 0 of the 300 panel targets, so
# the gain cannot have seen a panel response at all, in any line, at any point.
FIT_LINES_ZERO_PANEL = ["rpe1", "nadig_jurkat", "nadig_hepg2"]


def meta(n):
    return CD4 if n == "cd4_rest" else LINES[n]


def perts(n):
    return set(x.strip() for x in
               (ROOT / meta(n)["path"] / "perts.txt").read_text().splitlines() if x.strip())


def genes(n):
    return set(x.strip() for x in
               (ROOT / meta(n)["path"] / "genes.txt").read_text().splitlines() if x.strip())


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--zero-panel-fit", action="store_true",
                    help="fit the gain ONLY on lines carrying 0/300 panel targets")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    global FIT_LINES
    if args.zero_panel_fit:
        FIT_LINES = FIT_LINES_ZERO_PANEL
    rng = np.random.default_rng(args.seed)
    panel_lines = PANEL_LINES + ["cd4_rest"]
    allnames = sorted(set(panel_lines) | set(FIT_LINES))
    gaxis = sorted(set.intersection(*[genes(n) for n in allnames]) - PANEL_SET)
    fit_targets = sorted(set.intersection(*[perts(n) for n in FIT_LINES]) - PANEL_SET)
    eval_targets = sorted(set.intersection(*[perts(n) for n in panel_lines]) & PANEL_SET)
    print(f"gene axis {len(gaxis)} (panel coords dropped)\n"
          f"FIT targets  : {len(fit_targets)} NON-panel, shared by {FIT_LINES}\n"
          f"EVAL targets : {len(eval_targets)} PANEL, shared by {panel_lines}\n"
          f"overlap fit/eval: {len(set(fit_targets) & set(eval_targets))} (must be 0)")
    assert not (set(fit_targets) & set(eval_targets))

    E = {n: load_effects(meta(n)["path"], fit_targets, gaxis) for n in FIT_LINES}
    Ev = {n: load_effects(meta(n)["path"], eval_targets, gaxis) for n in panel_lines}
    # mean-remove each table on its OWN rows of the relevant block
    Cf = {n: v - v.mean(0, keepdims=True) for n, v in E.items()}
    Ce = {n: v - v.mean(0, keepdims=True) for n, v in Ev.items()}
    G = len(gaxis)
    gene_bins = np.arange(G)

    q = rng.permutation(len(fit_targets))
    tsel = np.sort(q[: int(0.6 * len(q))])
    tval = np.sort(q[int(0.6 * len(q)):])

    rows = []
    print(f"\n{'held':16s}{'k':>3}{'donors':<46}{'R0':>9}{'diag':>9}{'rel%':>8}"
          f"{'dpds':>9}{'shrink':>8}")
    for A in panel_lines:
        pool = [n for n in panel_lines if n != A]
        F = [n for n in FIT_LINES if n != A]
        Cfa = np.stack([Cf[n] for n in F])
        # fit the gain by inner leave-one-training-LINE-out on the NON-panel targets
        Ps, Ys, Fs = inner_pairs(Cfa, list(range(len(F))), tsel)
        Pv, Yv, Fv = inner_pairs(Cfa, list(range(len(F))), tval)
        hold = int(rng.integers(len(F)))
        ms, mv = Fs != hold, Fv == hold
        sh = max(DIAG_SHRINK, key=lambda s: macro_cos(
            Pv[mv] * diag_fit(Ps[ms], Ys[ms], gene_bins, G, s), Yv[mv]))
        Pa, Ya, _ = inner_pairs(Cfa, list(range(len(F))), np.arange(len(fit_targets)))
        g = diag_fit(Pa, Ya, gene_bins, G, sh)
        gshuf = g[rng.permutation(G)]                        # 004's gene-label null
        gmm = Ya.std(0) / np.maximum(Pa.std(0), 1e-300)      # moment-match, no target pairing
        gperm = diag_fit(Pa, Ya[rng.permutation(len(Ya))], gene_bins, G, sh)  # label derange

        Cea = np.stack([Ce[n] for n in panel_lines])
        ai = panel_lines.index(A)
        T = Cea[ai]
        for k in range(1, len(pool) + 1):
            for D in itertools.combinations([panel_lines.index(n) for n in pool], k):
                P = donor_mean(Cea, list(D), np.arange(len(eval_targets)))
                b, m = macro_cos(P, T), macro_cos(P * g, T)
                rows.append(dict(held=A, k=k, donors=[panel_lines[i] for i in D],
                                 R0=b, diag=m, rel=100 * (m - b) / b if b else 0.0,
                                 null_geneshuffle=macro_cos(P * gshuf, T),
                                 null_moment_match=macro_cos(P * gmm, T),
                                 null_label_perm=macro_cos(P * gperm, T),
                                 pds0=pds(P, T), pds1=pds(P * g, T),
                                 decoy0=decoy_cos(P, T), decoy1=decoy_cos(P * g, T),
                                 shrink=float(sh), n_fit_lines=len(F),
                                 pt0=per_target_cos(P, T).tolist(),
                                 pt1=per_target_cos(P * g, T).tolist()))
                print(f"{A:16s}{k:>3}{'+'.join(panel_lines[i] for i in D):<46}"
                      f"{b:9.5f}{m:9.5f}{100 * (m - b) / b if b else 0:+8.2f}"
                      f"{rows[-1]['pds1'] - rows[-1]['pds0']:+9.4f}{sh:8.1f}")

    print("\nSUMMARY by k (macro over held lines and donor subsets)")
    print(f"{'k':>3}{'n':>4}{'R0':>10}{'diag':>10}{'rel%':>9}{'dpds':>9}"
          f"{'geneshuf%':>11}{'moment%':>9}{'labperm%':>10}")
    out = {}
    for k in sorted({r["k"] for r in rows}):
        v = [r for r in rows if r["k"] == k]
        b = np.mean([r["R0"] for r in v])
        m = np.mean([r["diag"] for r in v])
        dp = np.mean([r["pds1"] - r["pds0"] for r in v])
        nulls = {n: 100 * (np.mean([r[f"null_{n}"] for r in v]) - b) / b
                 for n in ("geneshuffle", "moment_match", "label_perm")}
        out[k] = dict(n=len(v), R0=b, diag=m, rel=100 * (m - b) / b, dpds=dp, nulls=nulls)
        print(f"{k:>3}{len(v):4d}{b:10.5f}{m:10.5f}{100 * (m - b) / b:+9.2f}{dp:+9.4f}"
              f"{nulls['geneshuffle']:+11.2f}{nulls['moment_match']:+9.2f}"
              f"{nulls['label_perm']:+10.2f}")

    print("\nPER-HELD-LINE (macro over donor subsets at each line's own kmax)")
    print(f"{'held':16s}{'k':>3}{'R0':>10}{'diag':>10}{'rel%':>9}{'dpds':>9}")
    for A in panel_lines:
        v = [r for r in rows if r["held"] == A and r["k"] == max(
            x["k"] for x in rows if x["held"] == A)]
        b = np.mean([r["R0"] for r in v])
        m = np.mean([r["diag"] for r in v])
        print(f"{A:16s}{v[0]['k']:>3}{b:10.5f}{m:10.5f}{100 * (m - b) / b:+9.2f}"
              f"{np.mean([r['pds1'] - r['pds0'] for r in v]):+9.4f}")
    v = [r for r in rows if r["held"] != "kolf2"]
    b, m = np.mean([r["R0"] for r in v]), np.mean([r["diag"] for r in v])
    print(f"{'ALL minus kolf2':16s}{'':>3}{b:10.5f}{m:10.5f}{100 * (m - b) / b:+9.2f}"
          f"{np.mean([r['pds1'] - r['pds0'] for r in v]):+9.4f}   "
          f"(159 §4: kolf2 is anti-aligned and carries its own blob)")

    # paired target bootstrap at k = kmax, targets as the resampling unit
    kmax = max(r["k"] for r in rows)
    A0 = np.array([r["pt0"] for r in rows if r["k"] == kmax])
    A1 = np.array([r["pt1"] for r in rows if r["k"] == kmax])
    idx = rng.integers(0, A0.shape[1], size=(4000, A0.shape[1]))
    rel = 100 * (A1[:, idx].mean(axis=(0, 2)) - A0[:, idx].mean(axis=(0, 2))) \
        / A0[:, idx].mean(axis=(0, 2))
    obs = 100 * (A1.mean() - A0.mean()) / A0.mean()
    print(f"\npaired target bootstrap at k={kmax}: {obs:+.2f}%  95% CI "
          f"[{np.percentile(rel, 2.5):+.2f}, {np.percentile(rel, 97.5):+.2f}]  "
          f"P(>0) {(rel > 0).mean():.3f}  P(>gate 2.9%) {(rel > 2.9).mean():.3f}")

    for r in rows:
        r.pop("pt0")
        r.pop("pt1")
    (WORK / (f"crosspanel224{'-zero' if args.zero_panel_fit else ''}-{args.seed}.json")).write_text(json.dumps(dict(
        gene_axis=len(gaxis), n_fit_targets=len(fit_targets),
        n_eval_targets=len(eval_targets), fit_lines=FIT_LINES, panel_lines=panel_lines,
        rows=rows, summary={str(k): v for k, v in out.items()},
        bootstrap=dict(rel=float(obs), lo=float(np.percentile(rel, 2.5)),
                       hi=float(np.percentile(rel, 97.5)),
                       p_pos=float((rel > 0).mean()),
                       p_gate=float((rel > 2.9).mean()))), indent=1))


if __name__ == "__main__":
    main()

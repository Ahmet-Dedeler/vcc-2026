"""811 diagnostic (input to PROPOSAL.md; prices nothing): how well does each public line, and each consensus of lines, point
along a held cancer line's own response, in raw delta space, per panel (A/B/C + WORKING L1-L6; never L7/L8)?

For each panel P and line group g (base811.GROUPS): line row L_g(t) = normalise(mean of its banks' unit rows), non-panel genes.
  1. pairwise alignment: median over targets carried by both of cos(L_a(t), L_b(t)) on the genes both measure;
  2. consensus -> held line h (k562_gwps / orion_hct116 / orion_hek293t; kolf2 on A/B/C): cos(U(t), Y_h(t)) on h's measured
     genes, for U = sources only (the other cancer lines), all other lines, and all other lines minus one line (ablation);
  3. the across-target shared share n|mean|^2 / sum|row|^2 of Y_h, U_src, U_all (706's statistic; the pds-coherence risk).
Run on EC2 (reads the composition donors' delta rows): uv run python experiments/811-lolo-model/align811.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
OUT = Path(os.environ.get("ALIGN811_OUT", str(HERE / "work" / "align811.json")))
GENES = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
GPOS = {g: i for i, g in enumerate(GENES)}
PANELS = {"ABC": ROOT / "data/challenge/pert_counts.csv",
          **{f"L{k}": ROOT / f"experiments/801-lab-folds/work/L{k}/bundle/pert_counts.csv" for k in range(1, 7)}}
COMPOSITION = ("k562_gwps", "orion_hct116", "orion_hek293t")
# the same groups as base811.GROUPS (copied, so this script needs no engine import)
GROUPS = {
    "K562": ("k562_gwps", "viperturb_k562", "jiang_K562"),
    "HCT116": ("orion_hct116",), "HEK293T": ("orion_hek293t",),
    "CD4": ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"),
    "H1": ("h1_pan2025", "h1_vcc2025"), "KOLF2": ("kolf2",), "TeloHAEC": ("telohaec_pilot",), "DLD1": ("dld1_lowmoi",),
    "A549": ("altair_a549", "jiang_A549"), "RPE1": ("rpe1",), "Jurkat": ("nadig_jurkat",),
    "BXPC3": ("jiang_BXPC3",), "HAP1": ("jiang_HAP1",), "HT29": ("jiang_HT29",), "MCF7": ("jiang_MCF7",),
}
HELD = {"K562": "k562_gwps", "HCT116": "orion_hct116", "HEK293T": "orion_hek293t", "KOLF2": "kolf2"}
CANCER = ("K562", "HCT116", "HEK293T")


def panel(p):
    with open(p) as f:
        return [r["target_gene"] for r in csv.DictReader(f)]


_CACHE = {}


def bank(b):
    """(targets -> raw row on the official axis, measured mask)."""
    if b in _CACHE:
        return _CACHE[b]
    if b in COMPOSITION:
        D = ROOT / "data/donors" / b
        dg = (D / "genes.txt").read_text().split()
        dp = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        cols = np.array([GPOS.get(g, -1) for g in dg])
        ok = cols >= 0
        meas = np.zeros(len(GENES), bool)
        meas[cols[ok]] = True
        want = set().union(*[set(panel(p)) for p in PANELS.values()])
        rows = {}
        for i, p in enumerate(dp):
            if p in want and p not in rows:
                r = np.zeros(len(GENES), np.float32)
                r[cols[ok]] = np.asarray(delta[i], np.float32)[ok]
                rows[p] = r
        _CACHE[b] = (rows, meas)
    else:
        z = np.load(HERE / "work" / "banks" / f"{b}.npz")
        rows = np.asarray(z["rows"])       # read ONCE: z["rows"] inside the loop re-reads the whole array per target (OOM, run 075142Z)
        _CACHE[b] = ({str(t): rows[i] for i, t in enumerate(z["targets"])}, np.asarray(z["measured"], bool))
    return _CACHE[b]


def unitv(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def line_rows(g, targets, npn):
    """{t: unit line row} and the group's measured mask (union of its banks')."""
    meas = np.zeros(len(GENES), bool)
    acc = {}
    for b in GROUPS[g]:
        rows, m = bank(b)
        meas |= m
        for t in targets:
            if t in rows:
                r = np.asarray(rows[t], np.float64) * npn
                if np.any(r):
                    acc.setdefault(t, []).append(unitv(r))
    return {t: unitv(np.mean(v, 0)) for t, v in acc.items()}, meas & (npn > 0)


def cos_on(a, b, mask):
    a, b = a[mask], b[mask]
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan


def consensus(lines, rows_by_line, t):
    vs = [rows_by_line[g][t] for g in lines if t in rows_by_line[g]]
    return unitv(np.mean(vs, 0)) if vs else None


def shared(M):
    M = np.asarray(M)
    mu = M.mean(0)
    return float(len(M) * (mu @ mu) / (M * M).sum())


def main():
    t0 = time.time()
    res = {}
    for P, path in PANELS.items():
        targets = panel(path)
        npn = np.ones(len(GENES))
        for t in targets:
            if t in GPOS:
                npn[GPOS[t]] = 0.0
        R, MEAS = {}, {}
        for g in GROUPS:
            R[g], MEAS[g] = line_rows(g, targets, npn)
        pair = {}
        names = [g for g in GROUPS if R[g]]
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                m = MEAS[a] & MEAS[b]
                c = [cos_on(R[a][t], R[b][t], m) for t in targets if t in R[a] and t in R[b]]
                c = [x for x in c if np.isfinite(x)]
                if len(c) >= 5:
                    pair[f"{a}|{b}"] = dict(n=len(c), median=float(np.median(c)), mean=float(np.mean(c)))
        held_res = {}
        for h in (CANCER + ("KOLF2",) if P == "ABC" else CANCER):
            Y = R[h]
            others = [g for g in names if g != h]
            src = [g for g in CANCER if g != h]
            sets = {"src": src, "all": others, **{f"all-{g}": [x for x in others if x != g] for g in others}}
            # also: cancer-only extras, non-cancer only
            sets["src+CD4"] = src + ["CD4"]
            sets["nonsrc"] = [g for g in others if g not in src]
            out = {}
            for nm, S in sets.items():
                cs, Us = [], []
                for t in targets:
                    if t not in Y:
                        continue
                    U = consensus(S, R, t)
                    if U is None:
                        continue
                    cs.append(cos_on(U, Y[t], MEAS[h]))
                    Us.append(U * MEAS[h])
                cs = [x for x in cs if np.isfinite(x)]
                out[nm] = dict(n=len(cs), median=float(np.median(cs)) if cs else None, mean=float(np.mean(cs)) if cs else None,
                               shared_U=shared(Us) if len(Us) > 2 else None)
            out["shared_Y"] = shared([Y[t] * MEAS[h] for t in targets if t in Y])
            held_res[h] = out
        res[P] = dict(pairs=pair, held=held_res, n_lines={g: len(R[g]) for g in GROUPS})
        print(f"[{P}] " + " | ".join(f"{h}: src {v['src']['mean']:.4f} all {v['all']['mean']:.4f} nonsrc "
                                      f"{(v['nonsrc']['mean'] or float('nan')):.4f} (sharedU src {v['src']['shared_U']:.3f} "
                                      f"all {v['all']['shared_U']:.3f} Y {v['shared_Y']:.3f})" for h, v in held_res.items())
              + f" ({time.time() - t0:.0f}s)", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()

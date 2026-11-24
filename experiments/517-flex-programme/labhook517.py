#!/usr/bin/env python3
"""517 — a cross-line, same-lab analogue of "Arc's programme as the hook source on A/B/C", on public cancer banks.

On A/B/C the design comes from other labs' donor rows, and the hook would take a SAME-LAB, OTHER-LINE programme's sign on
uncalled genes. Xaira's lines (HCT116, HEK293T) and Weissman's (K562-GWPS, RPE1, Jurkat, HepG2) give that setting: for a
held-out line L and its targets t (non-essential ones where the screen has them, the panel's regime),
  design_t   = the OTHER lab's donor rows for t, averaged (target-specific, like our donor-transfer design)
  calls      = head (top 4 % of genes by |design_t|) + block (top 22 % by L's control expression): the production
               proportions (400 + 2,200 of ~10,000 gated genes); complement = the rest (the hook's territory)
  truth      = sign of L's own delta on (t, g) pairs with |delta| >= 0.2 (a significance proxy)
and compares, on the complement and on the block: the design's sign; the same-lab other line's programme (P_down over its
rows); the other lab's programme consensus (CLEAN4's role); all other lines. Coverage: every covered gene, or only the
source's most confident 63 % / 30 % by |z| (Arc's tau .05 / .1 coverage on H1).
On H1, Arc (same line AND lab) beats the design on uncalled gate genes by +0.15 (0.685 vs 0.532 at tau .05); this
measures what survives a line change within a lab. Banks read by row chunks (no mmap). Targets capped at 3,000 per line.
usage: uv run python experiments/517-flex-programme/labhook517.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from xline517 import BANKS, ESS, THR, NpyRows  # noqa: E402

PAIRS = {"HCT116": "HEK293T", "HEK293T": "HCT116", "K562": "RPE1", "RPE1": "K562", "Jurkat": "K562", "HepG2": "K562"}
COVS = (1.0, 0.63, 0.30)


def main():
    gl = {k: [l.strip() for l in open(d / "genes.txt")] for k, (_, d) in BANKS.items()}
    common = sorted(set.intersection(*map(set, gl.values())))
    G = len(common)
    P, bc, R = {}, {}, {}
    for k, (_, d) in BANKS.items():
        pos = {g: i for i, g in enumerate(gl[k])}
        col = np.array([pos[g] for g in common])
        bc[k] = np.load(d / "b_ctrl.npy")[col]
        dl = NpyRows(d / "delta.npy")
        perts = [l.strip() for l in open(d / "perts.txt")]
        dn = np.zeros(G); nz = np.zeros(G)
        for s in range(0, len(perts), 1000):
            X = np.asarray(dl.rows(list(range(s, min(s + 1000, len(perts))))), dtype=np.float32)[:, col]
            dn += (X < 0).sum(0); nz += (X != 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            P[k] = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
        R[k] = (dl, {p: i for i, p in enumerate(perts)}, col)
        print(f"{k}: {len(perts)} rows", flush=True)

    def zc(v):
        return (0.5 - v) / np.nanstd(v)

    def topcov(z, frac):
        """z where |z| is in the top `frac` of covered genes, else NaN (frac 1 = every covered gene)."""
        if frac >= 1.0:
            return z
        a = np.abs(z)
        thr = np.nanquantile(a, 1.0 - frac)
        return np.where(a >= thr, z, np.nan)
    n_head, n_block = int(round(0.04 * G)), int(round(0.22 * G))
    rng = np.random.default_rng(517)
    out = {}
    for L, partner in PAIRS.items():
        lab = BANKS[L][0]
        other_lab = [k for k in BANKS if BANKS[k][0] != lab]
        others = [k for k in BANKS if k != L]
        dl, pidx, col = R[L]
        genome_wide = L in ("K562", "HCT116", "HEK293T")
        tg = [t for t in pidx if (t not in ESS if genome_wide else True) and any(t in R[k][1] for k in other_lab)]
        if len(tg) > 3000:
            tg = list(rng.choice(tg, 3000, replace=False))
        tg = sorted(tg, key=lambda t: pidx[t])
        src0 = {f"same-lab line ({partner})": zc(P[partner]),
                "other-lab consensus": np.nanmean(np.vstack([zc(P[k]) for k in other_lab]), 0),
                "all-other-lines consensus": np.nanmean(np.vstack([zc(P[k]) for k in others]), 0)}
        srcs = {f"{k} @{c:.0%}": topcov(v, c) for k, v in src0.items() for c in COVS}
        block = np.zeros(G, bool); block[np.argsort(-bc[L])[:n_block]] = True
        acc = {f"{reg} | {k}": [0, 0] for reg in ("complement", "block") for k in ["design", *srcs]}
        for s in range(0, len(tg), 500):
            chunk = tg[s:s + 500]
            X = np.asarray(dl.rows([pidx[t] for t in chunk]), dtype=np.float32)[:, col]
            D = {}
            for k in other_lab:
                have = [t for t in chunk if t in R[k][1]]
                if have:
                    Y = np.asarray(R[k][0].rows(sorted(R[k][1][t] for t in have)), dtype=np.float32)[:, R[k][2]]
                    srt = sorted(have, key=lambda t: R[k][1][t])
                    for i_, t in enumerate(srt):
                        D.setdefault(t, []).append(Y[i_])
            for i_, t in enumerate(chunk):
                x = X[i_]
                dsg = np.mean(D[t], 0)
                head = np.zeros(G, bool); head[np.argsort(-np.abs(dsg))[:n_head]] = True
                gate = np.abs(x) >= THR
                for reg, m0 in (("complement", ~head & ~block & gate), ("block", block & ~head & gate)):
                    m = m0 & (dsg != 0)
                    acc[f"{reg} | design"][0] += int((np.sign(dsg[m]) == np.sign(x[m])).sum()); acc[f"{reg} | design"][1] += int(m.sum())
                    for k, v in srcs.items():
                        mk = m0 & np.isfinite(v)
                        acc[f"{reg} | {k}"][0] += int((np.sign(v[mk]) == np.sign(x[mk])).sum()); acc[f"{reg} | {k}"][1] += int(mk.sum())
        res = {k: {"acc": round(a / max(n, 1), 4), "n": n} for k, (a, n) in acc.items()}
        out[L] = {"partner": partner, "targets": len(tg), **res}
        print(f"\n{L} (same-lab partner {partner}; {len(tg)} targets{' non-essential' if genome_wide else ''})", flush=True)
        for k, v in res.items():
            print(f"   {k:55s} {v['acc']:.4f} (n {v['n']})", flush=True)
    (HERE / "work/labhook517.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""517 — how should a lab programme be ESTIMATED so that it carries to another line of the same lab? And how much does
Arc's small perturbation count (275 rows after the panel is excluded) cost? labhook517's instrument, Xaira's two genome-wide
lines (the cleanest same-lab pair): hold out L, design = the other lab's (Weissman) donor rows, truth = L's delta signs on
|delta| >= 0.2 pairs of its non-essential targets (2,000 sampled), scored on the complement (the hook's territory) at the
programme's most confident 63 % of genes (Arc's tau .05 coverage on H1), and on the block.
Programme variants, each from the PARTNER line's rows:
  pdown_all      P_down over all rows (the production rule, 396)
  pdown_N        P_down over N random rows, N = 275 (Arc's count) / 1,000; 5 draws each
  pdown_ess275   275 rows with Arc's essential share (36 %), 5 draws
  pdown_strong   rows in the top half by |delta| norm
  pdown_weak     rows in the bottom half
  meanz          sign of the mean of each row's z-scored delta (magnitude-weighted instead of a sign count)
usage: uv run python experiments/517-flex-programme/progest517.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from xline517 import BANKS, ESS, THR, NpyRows  # noqa: E402

PAIRS = {"HCT116": "HEK293T", "HEK293T": "HCT116"}
COV = 0.63


def main():
    gl = {k: [l.strip() for l in open(d / "genes.txt")] for k, (_, d) in BANKS.items()}
    common = sorted(set.intersection(*map(set, gl.values())))
    G = len(common)
    R, bc = {}, {}
    for k, (_, d) in BANKS.items():
        pos = {g: i for i, g in enumerate(gl[k])}
        col = np.array([pos[g] for g in common])
        bc[k] = np.load(d / "b_ctrl.npy")[col]
        R[k] = (NpyRows(d / "delta.npy"), [l.strip() for l in open(d / "perts.txt")], col)
    rng = np.random.default_rng(5175)

    def read(k, idx):
        dl, perts, col = R[k]
        out = np.empty((len(idx), G), np.float32)
        order = np.argsort(idx)
        srt = [idx[i] for i in order]
        for s in range(0, len(srt), 1000):
            out[order[s:s + 1000]] = np.asarray(dl.rows(srt[s:s + 1000]), dtype=np.float32)[:, col]
        return out

    def pdown(X):
        nz = (X != 0).sum(0); dn = (X < 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)

    def topcov(z, frac=COV):
        a = np.abs(z)
        return np.where(a >= np.nanquantile(a, 1.0 - frac), z, np.nan)
    weissman = [k for k in BANKS if BANKS[k][0] == "Weissman"]
    out = {}
    for L, partner in PAIRS.items():
        # the partner's programme variants
        pp = R[partner][1]
        Xp = read(partner, list(range(len(pp))))                       # all partner rows on the common axis
        ess = np.array([p in ESS for p in pp])
        norms = np.linalg.norm(Xp, axis=1)
        V = {"pdown_all": [pdown(Xp)]}
        for N in (275, 1000):
            V[f"pdown_{N}"] = [pdown(Xp[rng.choice(len(pp), N, replace=False)]) for _ in range(5)]
        ie, ine = np.flatnonzero(ess), np.flatnonzero(~ess)
        V["pdown_ess275"] = [pdown(Xp[np.concatenate([rng.choice(ie, 99, replace=False), rng.choice(ine, 176, replace=False)])])
                             for _ in range(5)]
        hi = norms >= np.median(norms)
        V["pdown_strong"] = [pdown(Xp[hi])]
        V["pdown_weak"] = [pdown(Xp[~hi])]
        sd = Xp.std(0); sd[sd == 0] = 1
        V["meanz"] = [(Xp / sd).mean(0)]            # mean z-scored delta: + = up (used as z directly below)
        # Arc-like small programmes (275 rows, 36 % essential = 99 + 176) and their estimator variants, 5 draws each
        draws = [np.concatenate([rng.choice(ie, 99, replace=False), rng.choice(ine, 176, replace=False)]) for _ in range(5)]
        V["pdown_ne176 (Arc's non-essential rows only)"] = [pdown(Xp[d[99:]]) for d in draws]
        V["pdown_ess99 (Arc's essential rows only)"] = [pdown(Xp[d[:99]]) for d in draws]
        V["meanz_ess275 (Arc-like mix, mean z)"] = [(Xp[d] / sd).mean(0) for d in draws]
        V["meanz_ne176"] = [(Xp[d[99:]] / sd).mean(0) for d in draws]
        V["pdown_strong_of_ess275 (top half by norm)"] = [pdown(Xp[d[np.argsort(-norms[d])[:137]]]) for d in draws]
        V["pdown_strong_of_ne176"] = [pdown(Xp[d[99:][np.argsort(-norms[d[99:]])[:88]]]) for d in draws]
        del Xp
        # L's truth and the other-lab design
        lp = R[L][1]
        tg = [i for i, p in enumerate(lp) if p not in ESS and any(p in R[k][1] for k in weissman)]
        tg = sorted(rng.choice(tg, min(2000, len(tg)), replace=False).tolist())
        XL = read(L, tg)
        n_head, n_block = int(round(0.04 * G)), int(round(0.22 * G))
        block = np.zeros(G, bool); block[np.argsort(-bc[L])[:n_block]] = True
        D = np.zeros((len(tg), G), np.float32); cnt = np.zeros(len(tg))
        for k in weissman:
            have = [(j, R[k][1].index(lp[t])) for j, t in enumerate(tg) if lp[t] in R[k][1]]
            if have:
                Y = read(k, [h for _, h in have])
                for (j, _), y in zip(have, Y):
                    D[j] += y; cnt[j] += 1
        D /= np.maximum(cnt, 1)[:, None]
        res = {}
        base = {"complement": [0, 0], "block": [0, 0]}
        masks = []
        for j in range(len(tg)):
            head = np.zeros(G, bool); head[np.argsort(-np.abs(D[j]))[:n_head]] = True
            gate = np.abs(XL[j]) >= THR
            masks.append((~head & ~block & gate, block & ~head & gate))
            for reg, m in zip(("complement", "block"), masks[-1]):
                mm = m & (D[j] != 0)
                base[reg][0] += int((np.sign(D[j][mm]) == np.sign(XL[j][mm])).sum()); base[reg][1] += int(mm.sum())
        res["design"] = {reg: round(a / max(n, 1), 4) for reg, (a, n) in base.items()}
        for name, draws in V.items():
            accs = {"complement": [], "block": []}
            for v in draws:
                z = v if name.startswith("meanz") else (0.5 - v)
                zc = topcov(z)
                for reg_i, reg in enumerate(("complement", "block")):
                    a = n = 0
                    for j in range(len(tg)):
                        m = masks[j][reg_i] & np.isfinite(zc)
                        a += int((np.sign(zc[m]) == np.sign(XL[j][m])).sum()); n += int(m.sum())
                    accs[reg].append(a / max(n, 1))
            res[name] = {reg: [round(float(np.mean(x)), 4), round(float(np.std(x)), 4)] for reg, x in accs.items()}
        out[L] = {"partner": partner, "targets": len(tg), **res}
        print(f"\n{L} (programme from {partner}; {len(tg)} non-essential targets; design = Weissman rows)", flush=True)
        for k, v in res.items():
            print(f"   {k:45s} complement {v['complement'] if k == 'design' else v['complement'][0]:.4f}"
                  + ("" if k == "design" else f" (sd {v['complement'][1]:.4f})")
                  + f" | block {v['block'] if k == 'design' else v['block'][0]:.4f}", flush=True)
    (HERE / "work/progest517.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""517 §2 — does the recipient-side IDENTITY axis ("genes a line expresses more than other lines go down under perturbation";
441 ps6, the xinyizhanglab preprint) carry across CANCER lines, and does it add to a cross-line consensus there? H1 says
yes (+0.011 gate-block over CLEAN4 in c4id2, not control noise: halves517b). H1 is an ESC; this checks public cancer lines.

Lines (bank rows = pseudobulk deltas vs each screen's own controls): Weissman K562 (GWPS), RPE1, Jurkat, HepG2; Xaira HCT116,
HEK293T. Common gene axis = every bank's genes. For each held-out line L:
  truth    sign of L's delta on (target, gene) pairs with |delta| >= 0.2 (a significance proxy), targets = L's rows
           (all; and the non-essential ones where the screen has them: GWPS, HCT116, HEK293T)
  consensus(other lines) z-mean of the OTHER lines' (0.5 - P_down), same lab + other lab ("CLEAN4-like"), or other lab only
  identity  -I_L = -(log frac_L - mean over the other lines of log frac), frac = control pooled fraction renormalised on the
           common genes; within-lab version against the same lab's other lines only
  combos    z(consensus) + lambda z(-I_L), lambda 0.5 / 1 / 2
Reported: pooled sign accuracy; the target-bootstrap is replaced by a per-target mean (targets vary a lot in size).
Caveat: the truth and I_L share L's control sample (no cells here for 462's halves); on |delta| >= 0.2 pairs H1 showed no
inflation (halves517b), but cancer banks have smaller controls (10-12k for RPE1 / Jurkat / HepG2).
usage: uv run python experiments/517-flex-programme/xline517.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BANKS = {"K562": ("Weissman", ROOT / "data/donors/k562_gwps"), "RPE1": ("Weissman", ROOT / "data/donors/rpe1"),
         "Jurkat": ("Weissman", ROOT / "data/donors/nadig_jurkat"),
         "HepG2": ("Weissman", ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2"),
         "HCT116": ("Xaira", ROOT / "data/donors/orion_hct116"), "HEK293T": ("Xaira", ROOT / "data/donors/orion_hek293t")}
ESS = ({l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
       | {l.strip() for l in open(ROOT / "data/donors/rpe1/perts.txt") if l.strip()})
THR = 0.2


class NpyRows:
    """Row chunks of a 2-D .npy WITHOUT mmap (file-backed mmap pages counted toward RSS: 6.6 GB for the Orion banks)."""
    def __init__(self, path):
        self.f = open(path, "rb")
        ver = np.lib.format.read_magic(self.f)
        rd = np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0
        self.shape, fortran, self.dtype = rd(self.f)
        assert not fortran and len(self.shape) == 2
        self.off = self.f.tell()
        self.rb = self.shape[1] * self.dtype.itemsize

    def rows(self, idx):
        """Rows at the sorted index list idx (contiguous runs read in one call)."""
        out = np.empty((len(idx), self.shape[1]), self.dtype)
        i = 0
        while i < len(idx):
            j = i
            while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
                j += 1
            self.f.seek(self.off + idx[i] * self.rb)
            out[i:j + 1] = np.fromfile(self.f, self.dtype, count=(j - i + 1) * self.shape[1]).reshape(j - i + 1, -1)
            i = j + 1
        return out


def main():
    gl = {k: [l.strip() for l in open(d / "genes.txt")] for k, (_, d) in BANKS.items()}
    common = sorted(set.intersection(*map(set, gl.values())))
    G = len(common)
    print(f"common genes {G}", flush=True)
    P, frac, rows = {}, {}, {}
    for k, (_, d) in BANKS.items():
        pos = {g: i for i, g in enumerate(gl[k])}
        col = np.array([pos[g] for g in common])
        bc = np.load(d / "b_ctrl.npy")[col]
        e = np.expm1(bc); frac[k] = e / e.sum()
        delta = NpyRows(d / "delta.npy")
        perts = [l.strip() for l in open(d / "perts.txt")]
        dn = np.zeros(G); nz = np.zeros(G)
        keep_rows = []
        for s in range(0, len(perts), 1000):
            X = np.asarray(delta.rows(list(range(s, min(s + 1000, len(perts))))), dtype=np.float32)[:, col]
            dn += (X < 0).sum(0); nz += (X != 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            P[k] = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
        rows[k] = (delta, perts, col)
        print(f"{k}: {len(perts)} rows, frac down {np.nanmean(P[k] > 0.5):.3f}", flush=True)
    logf = {k: np.log(np.maximum(v, 1e-9)) for k, v in frac.items()}

    def zc(v):
        return (0.5 - v) / np.nanstd(v)
    out = {}
    for L, (lab, _) in BANKS.items():
        others = [k for k in BANKS if k != L]
        same = [k for k in others if BANKS[k][0] == lab]
        other_lab = [k for k in others if BANKS[k][0] != lab]
        cons_all = np.nanmean(np.vstack([zc(P[k]) for k in others]), 0)
        cons_ol = np.nanmean(np.vstack([zc(P[k]) for k in other_lab]), 0)
        I_all = logf[L] - np.mean(np.vstack([logf[k] for k in others]), 0)
        I_same = (logf[L] - np.mean(np.vstack([logf[k] for k in same]), 0)) if same else None
        preds = {"consensus (all other lines)": cons_all, "consensus (other lab)": cons_ol, "-I (vs all other lines)": -I_all}
        if I_same is not None:
            preds["-I (vs same-lab lines)"] = -I_same
        for k in others:
            preds[f"single line {k} ({BANKS[k][0]})"] = zc(P[k])
        zI = -I_all / np.std(I_all)
        for lam in (0.5, 1.0, 2.0):
            preds[f"consensus(all) + {lam} (-I)"] = cons_all / np.nanstd(cons_all) + lam * zI
            preds[f"consensus(other lab) + {lam} (-I)"] = cons_ol / np.nanstd(cons_ol) + lam * zI
        delta, perts, col = rows[L]
        res = {}
        for sub in ("all", "non-essential"):
            ridx = [i for i, p in enumerate(perts) if sub == "all" or p not in ESS]
            if len(ridx) < 50:
                continue
            hit = {k: 0 for k in preds}; n = 0
            ptarget = {k: [] for k in preds}
            for s in range(0, len(ridx), 1000):
                X = np.asarray(delta.rows(ridx[s:s + 1000]), dtype=np.float32)[:, col]
                m = np.abs(X) >= THR
                n += int(m.sum())
                for k, v in preds.items():
                    sv = np.sign(v)
                    ok = m & np.isfinite(v)[None, :]
                    h = (np.sign(X) == sv[None, :]) & ok
                    hit[k] += int(h.sum())
                    cnt = ok.sum(1)
                    ptarget[k].extend(list((h.sum(1) / np.maximum(cnt, 1))[cnt >= 10]))
            res[sub] = {"pairs": n, "targets": len(ridx),
                        **{k: {"pooled": round(hit[k] / max(n, 1), 4), "per_target_mean": round(float(np.mean(ptarget[k])), 4)}
                           for k in preds}}
            print(f"\n{L} ({lab}) {sub}: {len(ridx)} targets, {n} pairs with |delta| >= {THR}", flush=True)
            for k in preds:
                print(f"   {k:40s} pooled {res[sub][k]['pooled']:.4f}  per-target {res[sub][k]['per_target_mean']:.4f}", flush=True)
        out[L] = res
    (HERE / "work/xline517.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

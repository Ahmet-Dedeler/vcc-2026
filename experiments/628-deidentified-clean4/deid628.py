#!/usr/bin/env python3
"""628: does CLEAN4 carry its BANK lines' own identity erosion, and does removing it transfer better to new lines?

Mechanism (639 / 622 / 627): under perturbation, a line's identity genes (expressed above other lines) tend to go DOWN.
CLEAN4 = the mean P_down of four CRISPRi banks: K562 GWPS, K562 essential, Nadig Jurkat, Nadig HepG2. Each bank's P_down
should then carry +beta_b·z(I_b), with I_b the BANK line's own identity. K562 enters twice. For a new line that
component is noise at best: K562's identity genes are not the new line's identity genes. The new line's own identity is
what 639's rider adds.

Per bank b (3' references only, since the banks are 3'; 622 found 3'/Flex mixing hurts identity):
  I_b(g) = b_ctrl_b(g) - mean over the reference profiles of b_ctrl(g), in log1p CP50k (donor b_ctrl.npy).
  References exclude b's own line and, per evaluation surface s, s's own line and same-study siblings. The surface's
  expression therefore never enters its own prior (639's rider is the only place it may enter).
  beta_b = OLS slope of (P_b - 1/2) on z(I_b) over the genes where both exist.
  P_b'   = P_b - beta_b · z(I_b) where I_b exists, else P_b.
CLEAN4' = mean over the allowed banks (LOLO, 636's BANK_FORBIDDEN) of P_b'.
Scored on 636's four gate surfaces (480's gates: h1 Flex, Orion HCT116 / HEK293T 3', K562 GWPS 3'), paired vs CLEAN4-LOLO
with 636's 2,000-draw target bootstrap:
  (1) CLEAN4' vs CLEAN4 (sign of P - 1/2);
  (2) the 639 rider on each, z(C - 1/2) + 0.5 z(I_s), with I_s = log1p ctrl CPM_s - mean over the other profiles
      (636 surfaces + A/B/C Flex, 639's construction); the rider's I_s is the same in both arms.
Analysis only, < 2 GB.
usage: uv run python experiments/628-deidentified-clean4/deid628.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import BANK_FORBIDDEN, G480, SURF, Inst  # noqa: E402

WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
CLEAN4_BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")
# 3' control profiles: name -> (dir, line)
REF3 = {
    "k562_gwps": (ROOT / "data/donors/k562_gwps", "K562"),
    "k562_essential": (ROOT / "data/donors/k562_essential", "K562"),
    "nadig_jurkat": (ROOT / "data/donors/nadig_jurkat", "Jurkat"),
    "nadig_hepg2": (ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2", "HepG2"),
    "rpe1": (ROOT / "data/donors/rpe1", "RPE1"),
    "orion_hct116": (ROOT / "data/donors/orion_hct116", "HCT116"),
    "orion_hek293t": (ROOT / "data/donors/orion_hek293t", "HEK293T"),
    "telohaec": (ROOT / "data/external/igvf-telohaec", "TeloHAEC"),
    "frangieh_melanoma": (ROOT / "experiments/166-donor-acquisition/donors/frangieh_melanoma", "melanoma"),
}
SURF_LINES = {"h1": {"H1"}, "hct116": {"HCT116", "HEK293T"}, "hek293t": {"HCT116", "HEK293T"}, "k562": {"K562"}}
MIN_REFS = 3


def load_profile(d: Path) -> dict:
    b = np.load(d / "b_ctrl.npy").astype(np.float64)
    genes = [ln.strip() for ln in open(d / "genes.txt") if ln.strip()]
    assert len(genes) == b.size, (d, len(genes), b.size)
    return dict(zip(genes, b))


def bank_identity(bank: str, surface: str, prof: dict, genes: np.ndarray) -> np.ndarray:
    line = REF3[bank][1]
    refs = [n for n, (_, ln) in REF3.items() if ln != line and ln not in SURF_LINES[surface]
            # one profile per line among the references (the two K562 banks are excluded with K562 anyway)
            and not (n == "k562_essential")]
    own = np.array([prof[bank].get(g, np.nan) for g in genes])
    R = np.array([[prof[n].get(g, np.nan) for g in genes] for n in refs])
    R = np.where(R > 0, R, np.nan)                     # a 0 in a padded axis = not measured
    cnt = np.isfinite(R).sum(0)
    mean = np.where(cnt >= MIN_REFS, np.nanmean(np.where(np.isfinite(R), R, np.nan), 0), np.nan)
    return np.where(np.isfinite(own) & (own > 0), own - mean, np.nan), refs


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def surface_identity(I: Inst):
    """639's I_s: log1p ctrl CPM_s - mean over the OTHER profiles (636 surfaces + A/B/C Flex), on genes all measure."""
    prof = {}
    for s, f in SURF.items():
        zz = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
    fr = np.load(ROOT / "experiments/604-identity-prior/work/flex_refs604.npz", allow_pickle=True)
    for c in "ABC":
        prof[f"ctx{c}"] = dict(zip(fr["genes"].astype(str), fr[f"cpm_{c}"].astype(float)))
    common = sorted(set.intersection(*(set(p) for p in prof.values())))
    L = {k: np.log1p(np.array([prof[k][g] for g in common])) for k in prof}
    out = {}
    for s in I.S:
        ident = L[s] - np.mean([L[k] for k in prof if k != s], 0)
        out[s] = I.gene_vec(s, common, ident)
    return out, len(common)


def pseudoP(score):
    """Sign-carrying pseudo P_down for 636's scorer (P > 1/2 -> predict down); NaN stays NaN."""
    return np.where(np.isfinite(score), 0.5 + 0.01 * np.tanh(score), np.nan)


def main() -> None:
    I = Inst()
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    cgenes = cp["genes"].astype(str)
    Pb = {b: cp["P_banks"].astype(float)[names.index(b)] for b in CLEAN4_BANKS}
    prof = {n: load_profile(d) for n, (d, _) in REF3.items()}
    Is, n_common = surface_identity(I)
    out = {"betas": {}, "rows": []}
    rows = []
    for s in I.S:
        allowed = [b for b in CLEAN4_BANKS if s not in BANK_FORBIDDEN[b]]
        Praw, Pdeid = [], []
        for b in allowed:
            Ib, refs = bank_identity(b, s, prof, cgenes)
            zI = z(Ib)
            ok = np.isfinite(Pb[b]) & np.isfinite(zI)
            beta, alpha = np.polyfit(zI[ok], Pb[b][ok] - 0.5, 1)
            r = float(np.corrcoef(zI[ok], Pb[b][ok])[0, 1])
            out["betas"][f"{s}/{b}"] = dict(beta=float(beta), alpha=float(alpha), r=r, n=int(ok.sum()), refs=refs)
            Praw.append(Pb[b])
            Pdeid.append(np.where(np.isfinite(zI), Pb[b] - beta * zI, Pb[b]))
        C4 = I.gene_vec(s, cgenes, np.nanmean(Praw, 0))
        C4d = I.gene_vec(s, cgenes, np.nanmean(Pdeid, 0))
        zi = z(Is[s])
        rider = lambda C: pseudoP(np.where(np.isfinite(C), (C - 0.5) / np.nanstd(C - 0.5) + 0.5 * np.nan_to_num(zi, nan=0.0), np.nan)
                                  + np.where(np.isfinite(zi) | ~np.isfinite(C), 0.0, 0.0))
        both_rider = np.isfinite(zi)
        r_c4 = np.where(both_rider, rider(C4), np.nan)
        r_c4d = np.where(both_rider, rider(C4d), np.nan)
        for label, P, ref in (("CLEAN4' (de-identified) vs CLEAN4", C4d, C4),
                              ("rider: z(CLEAN4') + 0.5 z(I_s) vs z(CLEAN4) + 0.5 z(I_s)", r_c4d, r_c4),
                              ("rider z(CLEAN4) + 0.5 z(I_s) vs CLEAN4 (639 check)", r_c4, np.where(both_rider, C4, np.nan))):
            for row in I.score(label, P_down={s: P}, ref={s: ref}):
                if row["targets"] == "all":
                    rows.append(row)
                    print(f"{s:8s} {label:62s} acc {row['acc']:.4f} d {row.get('d_vs_ref', np.nan):+.4f} "
                          f"[{row['d_ci'][0]:+.4f}, {row['d_ci'][1]:+.4f}] n {row['n_shared']}", flush=True)
    out["rows"] = rows
    print("\nbank identity slopes (P_down - 1/2 on z(I_b); > 0 = the bank's own identity genes go down in the bank):")
    for k, v in out["betas"].items():
        print(f"  {k:28s} beta {v['beta']:+.4f}  r {v['r']:+.3f}  n {v['n']}")
    out["n_common_surface_identity"] = n_common
    (WORK / "deid628.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""629: is Arc's H1 programme anti-informative on other lines BECAUSE it carries H1's identity erosion?

613 / 458d: Arc's own H1 programme (458d's artifact, P_down per gene from Arc's 2025 H1 screen; same lab and chemistry as
A-F) is +0.087 gate sign on H1 itself but −0.17 / −0.19 / −0.09 below CLEAN4 on HCT116 / HEK293T / K562. 458d hurt on
the board.
Hypothesis (622 / 627 / 628): a line's programme = a shared core + its own identity erosion (identity genes go down).
H1 is an ESC line against cancer lines, so its identity term is large and points the wrong way for other lines. If
removing it recovers a transferable core, Arc's own lab programme becomes a sign source for D/E/F.

Sources (DOWN where score > 0; z = scaled by sd, centred at 1/2 for P):
  CLEAN4          z(P_C4 - 1/2), LOLO (636's forbidden banks dropped)
  H1 raw          zP = z(P_H1 - 1/2)
  H1 de-identified  zP - beta·zI_H1, with beta = OLS slope of zP on zI_H1 over H1's genes.
                  I_H1 = log1p CPM_H1 - mean log1p CPM of A/B/C (flex_refs604, Flex).
  H1 expression-detrended  zP - beta_E·zE_H1, E = log1p CPM_H1 (609's expression-only control)
  H1 low-|I|      zP only where |zI_H1| < 0.5 (genes H1 expresses like the cancer lines)
  combos          zC4 + 0.5·zP_deid;  zC4 + 0.5·zI_s (639's rider);  zC4 + 0.5·zP_deid + 0.5·zI_s
Targets:
  (a) 636's gate surfaces HCT116, HEK293T, K562 (q < 0.05); I_s as in 639 (the other 636 surfaces + A/B/C Flex);
  (b) the Broad-16's 12 working lines (authors' fdr < 0.05; the 4 reserved lines never parsed); I_s same-lab (622g).
Read: paired gate-sign Δ vs CLEAN4 on genes both cover, 2,000-draw target bootstrap. Plus the raw programme's accuracy
by I_H1 tercile (where the anti-information sits).
Analysis only, < 4 GB.
usage: uv run python experiments/629-deidentified-h1-programme/deid629.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402
import broad622 as B  # noqa: E402

WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
RNG = np.random.default_rng(6290)
A458 = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def zc_half(P):
    return (P - 0.5) / np.nanstd(P - 0.5)


def slope(y, x):
    ok = np.isfinite(y) & np.isfinite(x)
    return float(np.polyfit(x[ok], y[ok], 1)[0])


def paired(score, ref, truth_dn, gate, n_boot=2000):
    """Δ sign accuracy of `score` vs `ref` (both DOWN where > 0) on gate cells where both are defined; target bootstrap."""
    both = np.isfinite(score) & np.isfinite(ref) & (score != 0) & (ref != 0)
    m = gate & both[None, :]
    a = (((score > 0)[None, :] == truth_dn) & m).sum(1).astype(float)
    b = (((ref > 0)[None, :] == truth_dn) & m).sum(1).astype(float)
    n = m.sum(1).astype(float)
    keep = n > 0
    a, b, n = a[keep], b[keep], n[keep]
    if n.sum() == 0:
        return dict(d=float("nan"), ci=[float("nan")] * 2, n=0, acc=float("nan"))
    idx = RNG.integers(0, len(n), size=(n_boot, len(n)))
    d = (a[idx].sum(1) - b[idx].sum(1)) / n[idx].sum(1)
    return dict(d=float((a.sum() - b.sum()) / n.sum()), ci=[float(np.quantile(d, .025)), float(np.quantile(d, .975))],
                n=int(n.sum()), acc=float(a.sum() / n.sum()), acc_ref=float(b.sum() / n.sum()))


def h1_sources():
    p = np.load(A458, allow_pickle=True)
    g = p["genes"].astype(str)
    P = np.asarray(p["pdown_A"], float)
    fr = np.load(ROOT / "experiments/604-identity-prior/work/flex_refs604.npz", allow_pickle=True)
    assert list(fr["genes"].astype(str)) == list(g)
    h1 = fr["cpm_h1"].astype(float)
    I = np.where(h1 > 5, np.log1p(h1) - np.mean([np.log1p(fr[f"cpm_{c}"].astype(float)) for c in "ABC"], 0), np.nan)
    E = np.where(h1 > 5, np.log1p(h1), np.nan)
    zP, zI, zE = zc_half(P), z(I), z(E)
    bI, bE = slope(zP, zI), slope(zP, zE)
    src = {"H1 raw": zP, "H1 de-identified": np.where(np.isfinite(zI), zP - bI * zI, np.nan),
           "H1 expression-detrended": np.where(np.isfinite(zE), zP - bE * zE, np.nan),
           "H1 low-|I|": np.where(np.abs(zI) < 0.5, zP, np.nan)}
    return g, src, zI, dict(beta_identity=bI, beta_expression=bE, corr_P_I=float(np.corrcoef(zP[np.isfinite(zP) & np.isfinite(zI)],
                                                                                              zI[np.isfinite(zP) & np.isfinite(zI)])[0, 1]))


def main() -> None:
    g, src, zI_h1, info = h1_sources()
    print(f"H1 programme vs its own identity: corr(zP, zI_H1) {info['corr_P_I']:+.3f}, beta {info['beta_identity']:+.3f}; "
          f"vs expression beta {info['beta_expression']:+.3f}", flush=True)
    out = {"info": info, "inst636": {}, "broad16": {}}
    # (a) 636 surfaces
    I = Inst()
    prof = {}
    for s, f in SURF.items():
        zz = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(zz["genes"].astype(str), zz["ctrl_cpm"].astype(float)))
    fr = np.load(ROOT / "experiments/604-identity-prior/work/flex_refs604.npz", allow_pickle=True)
    for c in "ABC":
        prof[f"ctx{c}"] = dict(zip(fr["genes"].astype(str), fr[f"cpm_{c}"].astype(float)))
    common = sorted(set.intersection(*(set(p) for p in prof.values())))
    Lg = {k: np.log1p(np.array([prof[k][x] for x in common])) for k in prof}
    for s in ("hct116", "hek293t", "k562"):
        D = I.S[s]
        truth_dn = D["lfc"] < 0
        gate = D["gate"]
        C4, _ = clean4_lolo(I, s)
        zc4 = zc_half(C4)
        Is = I.gene_vec(s, common, Lg[s] - np.mean([Lg[k] for k in prof if k != s], 0))
        zIs = z(Is)
        S = {k: I.gene_vec(s, g, v) for k, v in src.items()}
        zIh = I.gene_vec(s, g, zI_h1)
        comb = {"C4 + 0.5 H1deid": zc4 + 0.5 * np.nan_to_num(S["H1 de-identified"]),
                "C4 + 0.5 I_s (639 rider)": zc4 + 0.5 * np.nan_to_num(zIs),
                "C4 + 0.5 H1deid + 0.5 I_s": zc4 + 0.5 * np.nan_to_num(S["H1 de-identified"]) + 0.5 * np.nan_to_num(zIs)}
        res = {}
        for k, v in {**S, **comb}.items():
            res[k] = paired(v, zc4, truth_dn, gate)
        res["C4 + 0.5 H1deid + 0.5 I_s vs rider"] = paired(comb["C4 + 0.5 H1deid + 0.5 I_s"], comb["C4 + 0.5 I_s (639 rider)"], truth_dn, gate)
        # where the raw programme's anti-information sits: accuracy by I_H1 tercile (gate cells, both defined)
        terc = {}
        q = np.nanquantile(zIh, [1 / 3, 2 / 3])
        for lab, mk in (("low I_H1", zIh <= q[0]), ("mid", (zIh > q[0]) & (zIh <= q[1])), ("high I_H1", zIh > q[1])):
            sub = np.where(mk, S["H1 raw"], np.nan)
            terc[lab] = paired(sub, np.where(mk, zc4, np.nan), truth_dn, gate, n_boot=500)
        res["_raw_by_IH1_tercile"] = terc
        out["inst636"][s] = res
        print(f"\n== {s} (gate genes; Δ vs CLEAN4 [95% CI], n)")
        for k, r in res.items():
            if k.startswith("_"):
                continue
            print(f"   {k:38s} {r['d']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}]  acc {r['acc']:.3f} vs {r.get('acc_ref', float('nan')):.3f}  n {r['n']}")
        print("   raw H1 programme by I_H1 tercile: " + "; ".join(f"{k} {v['d']:+.4f} (n {v['n']})" for k, v in terc.items()))
    # (b) Broad-16
    gl, LFC, cl = B.load_matrix(B.P / "lfc_matrix.csv", 1)
    gf, FDR, cf = B.load_matrix(B.P / "fdr_matrix.csv", 1)
    gp, PB, cp = B.load_matrix(B.P / "pseudobulk_sum_by_ko.csv", 0)
    lines = sorted({ln for ln, _ in cl}); assert len(lines) == 12 and not (set(lines) & B.RESERVED)
    all_t = set(k for _, k in cl)
    gpos = {x: i for i, x in enumerate(gp)}; rowmap = np.array([gpos.get(x, -1) for x in gl])
    cpm = {}
    for ln in lines:
        cc = [j for j, (l2, k) in enumerate(cp) if l2 == ln and k not in all_t]
        s_ = np.nansum(PB[:, cc], 1); v = 1e6 * s_ / np.nansum(s_)
        cpm[ln] = np.where(rowmap >= 0, v[np.maximum(rowmap, 0)], np.nan)
    c4 = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    nm = list(c4["bank_names"].astype(str))
    Pc = c4["P_banks"].astype(float)[[nm.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    C4b = np.array([dict(zip(c4["genes"].astype(str), Pc)).get(x, np.nan) for x in gl])
    zc4b = zc_half(C4b)
    gmap = {x: i for i, x in enumerate(g)}
    to_b = lambda v: np.array([v[gmap[x]] if x in gmap else np.nan for x in gl])
    Sb = {k: to_b(v) for k, v in src.items()}
    colline = np.array([ln for ln, _ in cl]); coltgt = np.array([k for _, k in cl])
    rows = {}
    for ln in lines:
        jl = np.flatnonzero(colline == ln)
        L_, F_ = LFC[:, jl].T, FDR[:, jl].T
        own = coltgt[jl][:, None] == gl[None, :]
        gate = np.isfinite(L_) & (L_ != 0) & (F_ < 0.05) & ~own
        tdn = L_ < 0
        others = [o for o in lines if o != ln]
        Is = np.where(cpm[ln] > 5, np.log1p(cpm[ln]) - np.nanmean([np.log1p(cpm[o]) for o in others], 0), np.nan)
        zIs = z(Is)
        comb = {"C4 + 0.5 H1deid": zc4b + 0.5 * np.nan_to_num(Sb["H1 de-identified"]),
                "C4 + 0.5 I_s (rider)": zc4b + 0.5 * np.nan_to_num(zIs),
                "C4 + 0.5 H1deid + 0.5 I_s": zc4b + 0.5 * np.nan_to_num(Sb["H1 de-identified"]) + 0.5 * np.nan_to_num(zIs)}
        r = {k: paired(v, zc4b, tdn, gate, n_boot=500) for k, v in {**Sb, **comb}.items()}
        r["C4 + 0.5 H1deid + 0.5 I_s vs rider"] = paired(comb["C4 + 0.5 H1deid + 0.5 I_s"], comb["C4 + 0.5 I_s (rider)"], tdn, gate, n_boot=500)
        rows[ln] = r
    out["broad16"] = rows
    keys = list(rows[lines[0]].keys())
    print("\n== Broad-16, 12 working lines (Δ vs CLEAN4; KO gates lean UP, so read relative orderings)")
    summ = {}
    for k in keys:
        ds = np.array([rows[ln][k]["d"] for ln in lines])
        summ[k] = dict(mean=float(np.nanmean(ds)), lines_pos=int(np.nansum(ds > 0)), min=float(np.nanmin(ds)), max=float(np.nanmax(ds)))
        print(f"   {k:38s} mean {summ[k]['mean']:+.4f}  lines > 0: {summ[k]['lines_pos']}/12  range [{summ[k]['min']:+.4f}, {summ[k]['max']:+.4f}]")
    out["broad16_summary"] = summ
    (WORK / "deid629.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

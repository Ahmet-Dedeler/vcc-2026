"""516 step 3: the cancer folds (480's gates: Orion HCT116 / HEK293T thinned harnesses, K562-GWPS panel cells; 150
cells per target, relative direction only). No emitter output exists there, so production is emulated:

  base_mean     sign of the mean of the fold's available production donors (log2 ratio, pc 1 CPM). 430: on cross-study
                cancer folds the raw donor mean is within +0.0012 [-0.0007, +0.0031] of the whitened production
                composition, so this is the production carrier's analog there.
  base_mcd4     base_mean + 0.34 x the CD4 three-condition mean log2FC (production's CD4 partial rows).
  +c4q80        either, with the CLEAN4-minus-L sign on the fold's CPM block (top 2,200 by the fold's control CPM) except
                the strongest 20 % |analog| per target (430's PPprior emulation of the delivered clean4q80 structure).
Donor sets: CROSS-STUDY (the board's situation: no A/B/C donor shares their study) and LITERAL (adds the same-study Orion
sibling). The held-out line's own bank is never used; VIPerturb K562 is dropped on the K562 fold (same line).
Candidates: Stouffer of standardised uncentred scores over the fold's sources + CD4 x3 (+ VIPerturb), applied where
|Z| >= thr on non-block genes (on H1, block ~ called and rest ~ uncalled) or on all genes; the analog elsewhere.
Priced as pooled pool-sign accuracy by band with a target-bootstrap CI of the difference, target-mean accuracy,
|l|-weighted accuracy and an nmae analog at the shipped MDE law's magnitudes (480's mde_cons).
    PYTHONDONTWRITEBYTECODE=1 uv run python fold516.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

from common516 import ROOT, WORK, boot_mean_delta, sign_table
from common480 import gate_from, nmae_rows

FOLDS = {
    "orion_hct116": dict(own="orion_hct116", sib="orion_hek293t", prod=["k562_gwps", "orion_hek293t"],
                         extra=["viperturb_k562"], c4_drop=[]),
    "orion_hek293t": dict(own="orion_hek293t", sib="orion_hct116", prod=["k562_gwps", "orion_hct116"],
                          extra=["viperturb_k562"], c4_drop=[]),
    "k562_gwps": dict(own="k562_gwps", sib=None, prod=["orion_hct116", "orion_hek293t"], extra=[],
                      c4_drop=["k562_gwps", "k562_essential"]),
}
CD4 = ["cd4_Rest", "cd4_Stim8hr", "cd4_Stim48hr"]
CLEAN4 = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]


def mde_add(cnt, c=0.1, cap=0.1, mn=0.3):
    return np.where(cnt >= mn, np.minimum(c / np.sqrt(np.maximum(cnt, 1e-3)), cap), 0.0)


def load_fold(name):
    g = np.load(ROOT / f"experiments/480-nmae-instrument/work/gate-{name}.npz", allow_pickle=True)
    z = np.load(WORK / f"src-{name}.npz", allow_pickle=True)
    T, G = g["targets"].astype(str), g["genes"].astype(str)
    assert (z["targets"].astype(str) == T).all() and (z["genes"].astype(str) == G).all()
    q, l = g["q"].astype(float), g["lfc"].astype(float)
    own = G[None, :] == T[:, None]
    rs = (q < 0.05) & ~own
    pool = rs & np.isfinite(l) & (l != 0)
    n_real = rs.sum(1)
    band = np.select([n_real < 20, n_real < 100], [0, 1], 2)
    cnt = g["ctrl_cpm"].astype(float) * float(g["lib_med"]) / 1e6
    return dict(name=name, T=T, G=G, q=q, l=l, pool=pool, band=band, truth=np.sign(np.nan_to_num(l)),
                gateT=gate_from(q, l, G, T), ctrl=g["ctrl_cpm"].astype(float), cnt=cnt, z=z)


def raw(z, name):
    return (z[f"{name}__log_fc"] if name.startswith("cd4_") else z[f"{name}__x"]).astype(float)


def std(z, name, cap=3.0, centred=False):
    x = raw(z, name)
    if centred:
        x = x - z[f"{name}__med"][None, :]
    return np.clip(x / (z[f"{name}__mad"][None, :] * z[f"{name}__tau"][:, None]), -cap, cap)


def stouffer(z, names, cap=3.0, weights=None, centred=False):
    Ss = [std(z, n, cap, centred) for n in names]
    w = np.ones(len(names)) if weights is None else np.asarray(weights, float)
    num = sum(wi * np.nan_to_num(s) for wi, s in zip(w, Ss))
    den = np.sqrt(sum((wi ** 2) * np.isfinite(s) for wi, s in zip(w, Ss)))
    return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)


def c4_sign(G, drop):
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    bn = list(cp["bank_names"].astype(str))
    banks = [b for b in CLEAN4 if b not in drop]
    P = cp["P_banks"][[bn.index(b) for b in banks]].mean(0)
    gi = {g: i for i, g in enumerate(cp["genes"].astype(str))}
    return np.array([(-1.0 if P[gi[g]] > 0.5 else 1.0) if g in gi else -1.0 for g in G])   # builder: -1 uncovered


def analogs(F, spec, donors):
    z, G, T = F["z"], F["G"], F["T"]
    pan = np.stack([raw(z, d) for d in donors])
    mean = np.where(np.isfinite(pan).any(0), np.nanmean(pan, 0), np.nan)
    cd4m = np.nanmean(np.stack([raw(z, c) for c in CD4]), 0)
    out = {"base_mean": np.nan_to_num(mean), "base_mcd4": np.nan_to_num(mean) + 0.34 * np.nan_to_num(cd4m)}
    blk = np.zeros(len(G), bool); blk[np.argsort(-F["ctrl"])[:2200]] = True
    c4 = c4_sign(G, spec["c4_drop"])
    for k in list(out):
        R = out[k]; E = np.sign(R)
        bc = np.flatnonzero(blk)
        for i in range(len(T)):
            v = np.abs(R[i, bc]); thr = np.quantile(v, 0.8)
            flip = bc[v <= thr]
            E[i, flip] = c4[flip]
        out[k + "+c4q80"] = E
        out[k] = np.sign(R)
    return out, blk


def price_fold(F, base, cand, label, rows):
    """cand: +1/-1 where it has an opinion, 0 elsewhere -> hybrid with base."""
    H = np.where(cand != 0, cand, base)
    tab = sign_table(F["truth"], F["pool"], F["band"], base, H)
    tabw = sign_table(F["truth"], F["pool"], F["band"], base, H, weight=np.abs(np.nan_to_num(F["l"])))
    m = F["pool"] & (base != 0); n = m.sum(1); ok = n >= 10
    ha = ((base == F["truth"]) & m).sum(1) / np.maximum(n, 1); hb = ((H == F["truth"]) & m).sum(1) / np.maximum(n, 1)
    tm = boot_mean_delta(np.where(ok, ha, np.nan), np.where(ok, hb, np.nan))
    add = mde_add(F["cnt"])[None, :]
    na, nb = nmae_rows(base * add, F["gateT"]), nmae_rows(H * add, F["gateT"])
    nm = boot_mean_delta(na, nb)
    a = tab["all"]
    line = (f"{F['name'][:13]:13s} {label[:58]:58s} pool {a[0]:.4f}->{a[1]:.4f} {a[2]:+.4f} [{a[3]:+.4f},{a[4]:+.4f}] "
            f"b<20 {tab['<20'][2]:+.3f} b20 {tab['20-100'][2]:+.3f} b100 {tab['>=100'][2]:+.4f} | w {tabw['all'][2]:+.4f} "
            f"| tmean {tm[0]:+.4f} [{tm[1]:+.4f},{tm[2]:+.4f}] | nmae(mde) {np.nanmean(na):.4f} {nm[0]:+.5f} "
            f"[{nm[1]:+.5f},{nm[2]:+.5f}]")
    print(line, flush=True)
    rows.append(dict(fold=F["name"], label=label, acc={k: list(map(float, v)) for k, v in tab.items()},
                     acc_w={k: list(map(float, v)) for k, v in tabw.items()}, tmean=tm, nmae_base=float(np.nanmean(na)),
                     nmae_delta=nm))
    return H


def main():
    rows = []
    for name, spec in FOLDS.items():
        F = load_fold(name)
        print(f"\n### {name}: {len(F['T'])} targets, pool cells {int(F['pool'].sum())}, bands (targets) "
              f"{[int((F['band'] == b).sum()) for b in range(3)]}, cells {[int(F['pool'][F['band'] == b].sum()) for b in range(3)]}")
        for view in ("xstudy", "literal"):
            if view == "literal" and spec["sib"] is None:
                continue
            donors = [d for d in spec["prod"] if view == "literal" or d != spec["sib"]]
            A, blk = analogs(F, spec, donors)
            srcs = donors + CD4 + spec["extra"]
            print(f"-- view {view}: production donors {donors}; Stouffer sources {srcs}")
            zero = np.zeros_like(F["truth"])
            for bn, base in A.items():
                price_fold(F, base, zero, f"[{view}] KA {bn} vs itself", rows)
            if view == "literal" or spec["sib"] is None:
                # known answer: 480's mde_cons on this fold = sign(mean of the literal donors) x the MDE law
                want = json.load(open(ROOT / "experiments/480-nmae-instrument/work/fold480.json"))[name]["nmae"]["mde_cons"]
                got = rows[-len(A)]["nmae_base"]
                print(f"KNOWN ANSWER {'PASS' if abs(got - want) < 2e-4 else 'FAIL'}: {name} mde_cons nmae {got:.6f} vs "
                      f"480's {want:.6f}", flush=True)
                if abs(got - want) >= 2e-4:
                    sys.exit(1)
            Zs = {"stz(src)": stouffer(F["z"], srcs), "stz(donors+cd4x3)": stouffer(F["z"], donors + CD4),
                  "stz(donors)": stouffer(F["z"], donors), "stz(cd4x3)": stouffer(F["z"], CD4)}
            for bn in ("base_mean", "base_mcd4+c4q80", "base_mean+c4q80"):
                base = A[bn]
                for zn, Z in Zs.items():
                    for thr in (0.0, 0.5, 1.0):
                        for region, rm in (("rest", ~blk), ("all", np.ones_like(blk))):
                            C = np.where(rm[None, :] & (np.abs(np.nan_to_num(Z)) >= thr), np.sign(np.nan_to_num(Z)), 0.0)
                            price_fold(F, base, C, f"[{view}] {zn} |Z|>={thr} {region} vs {bn}", rows)
    json.dump(rows, open(WORK / "fold516.json", "w"), indent=1)


if __name__ == "__main__":
    main()

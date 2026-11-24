"""159 SS4 -- the line-pair transfer matrix, so the answer is not one accident of tissue distance.

Every donor table is the SAME object in the SAME space as the harness truth
(`scripts/build_donor_table.py`: delta = b_p - b_ctrl, b = log1p(TS*P/sum P), TS = 5e4), so
these numbers are directly comparable to SS3's.  What the donor tables do NOT have is
cell-level data, so no split-half is possible and no pair here can be disattenuated from its
own data.

That is why the primary column is the RATIO same/diff, not the excess: SS3 shows empirically
(a 10x depth sweep on the Orion pair) that the ratio is invariant to attenuation while the
excess is not.  The ratio is therefore the one quantity that can be compared across pairs
whose depths differ by an order of magnitude.

`k562_gwps x k562_essential` is the POSITIVE CONTROL: both are Replogle K562, i.e. the SAME
cell line measured in two different libraries.  Whatever a same-line pair reads is the level a
cross-line pair would have to reach for "the response is a property of the target gene" to be
true.  `k562_gwps x perturbmatch` (if present) is a second, stronger control -- a reanalysis of
the identical cells.
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

SEED = 20260902
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
MAX_PERT = 3000
MIN_PERT = 30

DONORS = {
    "k562_gwps": ("data/donors/k562_gwps", "K562 (Replogle gwps)"),
    "k562_essential": ("data/donors/k562_essential", "K562 (Replogle essential)"),
    "rpe1": ("data/donors/rpe1", "RPE1"),
    "nadig_jurkat": ("data/donors/nadig_jurkat", "Jurkat"),
    "orion_hct116": ("data/donors/orion_hct116", "HCT116 (Orion)"),
    "orion_hek293t": ("data/donors/orion_hek293t", "HEK293T (Orion)"),
    "h1_2025": ("data/donors/h1_2025", "H1 hESC"),
    "kolf2": ("data/external/kolf2", "KOLF2.1J iPSC"),
    "cd4_rest": ("data/external/cd4t/cloud-compact/cd4_compact/cd4_rest", "CD4 T rest"),
    "cd4_stim48hr": ("data/external/cd4t/cloud-compact/cd4_compact/cd4_stim48hr",
                     "CD4 T stim 48h"),
}
SAME_LINE = {("k562_gwps", "k562_essential")}


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def load_donor(rel):
    p = os.path.join(REPO, rel)
    if not os.path.exists(os.path.join(p, "delta.npy")):
        return None
    delta = np.load(os.path.join(p, "delta.npy"))
    perts = [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()]
    genes = [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()]
    cells = None
    mp = os.path.join(p, "meta.json")
    if os.path.exists(mp):
        m = json.load(open(mp))
        cells = (m.get("cells_per_pert") or {}).get("median")
    assert delta.shape == (len(perts), len(genes)), f"{rel}: shape mismatch"
    return dict(delta=np.nan_to_num(delta.astype(np.float32)), perts=perts, genes=genes,
                cells=cells)


def pair(dA, dB, rng):
    pa = {p: i for i, p in enumerate(dA["perts"])}
    pb = {p: i for i, p in enumerate(dB["perts"])}
    sh = sorted(set(pa) & set(pb))
    if len(sh) < MIN_PERT:
        return None
    if len(sh) > MAX_PERT:
        sh = [sh[i] for i in sorted(rng.choice(len(sh), MAX_PERT, replace=False))]
    ga = {g: i for i, g in enumerate(dA["genes"])}
    gb = {g: i for i, g in enumerate(dB["genes"])}
    shg = sorted(set(ga) & set(gb))
    drop = set(sh)
    shg = [g for g in shg if g not in drop]        # pds_cosine drops panel target coordinates
    if len(shg) < 500:
        return None
    A = dA["delta"][np.asarray([pa[p] for p in sh])][:, np.asarray([ga[g] for g in shg])]
    B = dB["delta"][np.asarray([pb[p] for p in sh])][:, np.asarray([gb[g] for g in shg])]
    A, B = A.astype(np.float64), B.astype(np.float64)
    out = {"n_perts": len(sh), "n_genes": len(shg)}
    for space in ("full", "resid"):
        X = A - A.mean(0, keepdims=True) if space == "resid" else A
        Y = B - B.mean(0, keepdims=True) if space == "resid" else B
        a, b = unit(X), unit(Y)
        S = a @ b.T
        same = float(np.diag(S).mean())
        diff = offdiag_mean(S)
        out[space] = {"same": same, "same_median": float(np.median(np.diag(S))),
                      "diff": diff, "excess": same - diff,
                      "ratio": same / diff if abs(diff) > 1e-6 else float("nan"),
                      "within_A": offdiag_mean(a @ a.T), "within_B": offdiag_mean(b @ b.T)}
        # a permutation null on the SAME statistic: destroys the target correspondence only
        if space == "resid":
            perm = rng.permutation(len(sh))
            out[space]["same_permuted"] = float((a * b[perm]).sum(1).mean())
    return out


def main() -> None:
    rng = np.random.default_rng(SEED)
    loaded = {}
    for k, (rel, lbl) in DONORS.items():
        d = load_donor(rel)
        if d is None:
            print(f"  [skip] {k}: no delta.npy at {rel}")
            continue
        d["label"] = lbl
        loaded[k] = d
        print(f"  {k:16} {len(d['perts']):6d} perts x {len(d['genes']):6d} genes  "
              f"median cells/pert {d['cells']}")

    rows = []
    for a, b in itertools.combinations(sorted(loaded), 2):
        r = pair(loaded[a], loaded[b], rng)
        if r is None:
            continue
        r["A"], r["B"] = a, b
        r["same_line"] = (a, b) in SAME_LINE or (b, a) in SAME_LINE
        rows.append(r)

    rows.sort(key=lambda r: -r["full"]["ratio"])
    print(f"\n{'=' * 122}\nLINE-PAIR TRANSFER MATRIX  "
          f"(full-vector space; RATIO is the attenuation-invariant column)")
    print(f"{'line A':<16}{'line B':<16}{'n_pert':>7}{'n_gene':>8}"
          f"{'same':>9}{'diff':>9}{'excess':>10}{'RATIO':>8}"
          f"{'wA':>8}{'wB':>8}{'resid same':>12}{'resid null':>11}")
    for r in rows:
        f, s = r["full"], r["resid"]
        tag = " *SAME LINE*" if r["same_line"] else ""
        print(f"{r['A']:<16}{r['B']:<16}{r['n_perts']:7d}{r['n_genes']:8d}"
              f"{f['same']:+9.4f}{f['diff']:+9.4f}{f['excess']:+10.5f}{f['ratio']:8.4f}"
              f"{f['within_A']:+8.4f}{f['within_B']:+8.4f}"
              f"{s['same']:+12.4f}{s['same_permuted']:+11.4f}{tag}")

    with open(os.path.join(WORK, "donorpairs.json"), "w") as fh:
        json.dump({"rows": rows, "max_pert": MAX_PERT, "seed": SEED}, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'donorpairs.json')}")


if __name__ == "__main__":
    main()

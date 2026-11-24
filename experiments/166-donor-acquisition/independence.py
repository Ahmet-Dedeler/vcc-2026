"""How independent are the two new lines -- in basal space and in effect space?

Independence is the whole point (brief rule 3): a second screen of a line we already hold buys
far less than a new tissue. Two measurements, in the units this repo already uses.

**A. Basal profile.** Cosine of `b_ctrl` between every pair of lines, and between each line and
the three challenge contexts (whose control cells give `b_ctrl` directly). Reported twice: raw,
and after removing the across-line mean profile. The raw number is ~0.9 for any two human cell
types because the housekeeping core dominates; the mean-removed number is the one that
discriminates, exactly as 159 mean-removes before reading any effect statistic.

**B. Effect transfer.** 159's line-pair statistic, same code shape as
`experiments/159-transfer-ceiling/donorpairs.py`: same-target cross-line cosine, decoy
(different-target) cosine, their RATIO -- the attenuation-invariant column -- and a
permutation null on the residual space.

KNOWN-ANSWER ROW (AGENTS.md 11.1 rule 10). This script recomputes three pairs whose values 159
already published and aborts if they do not reproduce:

    k562_gwps x orion_hct116     same +0.0321  diff +0.0020  resid same +0.0321
    k562_gwps x orion_hek293t    same +0.0289  diff +0.0016  resid same +0.0286
    orion_hct116 x orion_hek293t same +0.0306  diff +0.0016  resid same +0.0300

If those three do not come back, nothing else in this file is readable.

Writes `work/independence.json`.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import h5py
import numpy as np

SEED = 20260902
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TS = 5e4
MAX_PERT = 3000
MIN_PERT = 30

DONORS = {
    "k562_gwps": ("data/donors/k562_gwps", "K562 — CML, blood"),
    "k562_essential": ("data/donors/k562_essential", "K562 — CML, blood"),
    "rpe1": ("data/donors/rpe1", "RPE1 — retinal pigment epithelium"),
    "nadig_jurkat": ("data/donors/nadig_jurkat", "Jurkat — T-ALL, blood"),
    "orion_hct116": ("data/donors/orion_hct116", "HCT116 — colorectal"),
    "orion_hek293t": ("data/donors/orion_hek293t", "HEK293T — embryonic kidney"),
    "h1_2025": ("data/donors/h1_2025", "H1 — hESC"),
    "kolf2": ("data/external/kolf2", "KOLF2.1J — iPSC"),
    # the two built here
    "nadig_hepg2": ("experiments/166-donor-acquisition/donors/nadig_hepg2",
                    "HepG2 — hepatoblastoma, LIVER  [NEW]"),
    "frangieh_melanoma": ("experiments/166-donor-acquisition/donors/frangieh_melanoma",
                          "patient melanoma — MELANOCYTE  [NEW]"),
}
NEW = ("nadig_hepg2", "frangieh_melanoma")

KNOWN = {  # 159 RESULTS.md section 4, to 4 dp
    ("k562_gwps", "orion_hct116"): (0.0321, 0.0020, 0.0321),
    ("k562_gwps", "orion_hek293t"): (0.0289, 0.0016, 0.0286),
    ("orion_hct116", "orion_hek293t"): (0.0306, 0.0016, 0.0300),
}
# 159's estimator caps a pair at MAX_PERT targets by a random draw, and the SAME draw also
# sets the gene axis (the sampled targets' own coordinates are dropped). So each published row
# is one realisation, not a point estimate, and a fixed absolute tolerance is the wrong gate:
# measured over 12 draws, `orion_hct116 x orion_hek293t` gives same 0.0285 +- 0.0011, and 159's
# published 0.0306 is the top of that range. The gate therefore re-draws and asks whether 159's
# value lies within 3 sd of the mean -- which still catches a real bug (those miss by ~10x)
# without pretending a 2-sd draw is a reproduction.
KNOWN_DRAWS = 5
KNOWN_TOL_FLOOR = 0.0015
KNOWN_TOL_SD = 3.0


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def load_donor(rel):
    p = REPO / rel
    if not (p / "delta.npy").exists():
        return None
    delta = np.load(p / "delta.npy")
    perts = [l.strip() for l in (p / "perts.txt").read_text().split("\n") if l.strip()]
    genes = [l.strip() for l in (p / "genes.txt").read_text().split("\n") if l.strip()]
    b = np.load(p / "b_ctrl.npy") if (p / "b_ctrl.npy").exists() else None
    assert delta.shape == (len(perts), len(genes)), f"{rel}: shape mismatch"
    return dict(delta=np.nan_to_num(delta.astype(np.float32)), perts=perts, genes=genes,
                b_ctrl=None if b is None else np.asarray(b, dtype=np.float64))


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
    shg = [g for g in sorted(set(ga) & set(gb)) if g not in set(sh)]
    if len(shg) < 500:
        return None
    A = dA["delta"][[pa[p] for p in sh]][:, [ga[g] for g in shg]].astype(np.float64)
    B = dB["delta"][[pb[p] for p in sh]][:, [gb[g] for g in shg]].astype(np.float64)
    out = {"n_perts": len(sh), "n_genes": len(shg)}
    for space in ("full", "resid"):
        X = A - A.mean(0, keepdims=True) if space == "resid" else A
        Y = B - B.mean(0, keepdims=True) if space == "resid" else B
        a, b = unit(X), unit(Y)
        S = a @ b.T
        same, diff = float(np.diag(S).mean()), offdiag_mean(S)
        out[space] = {"same": same, "diff": diff, "excess": same - diff,
                      "ratio": same / diff if abs(diff) > 1e-6 else float("nan"),
                      "within_A": offdiag_mean(a @ a.T), "within_B": offdiag_mean(b @ b.T)}
        if space == "resid":
            out[space]["same_permuted"] = float(
                (a * b[rng.permutation(len(sh))]).sum(1).mean())
    return out


def context_b_ctrl(path: Path):
    """log1p(TS * p) of a challenge context's 18,400 control cells, streamed."""
    with h5py.File(path, "r") as h:
        X = h["X"]
        g = h["var"]
        key = g.attrs.get("_index", "_index")
        key = key.decode() if isinstance(key, bytes) else key
        node = g[key]
        # the challenge contexts store `.var` index as a nullable-string group
        # ({'mask','values'}), not a plain dataset
        raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
        genes = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in raw])
        n = len(genes)
        tot = np.zeros(n, dtype=np.float64)
        d, idx = X["data"], X["indices"]
        step = 20_000_000
        for lo in range(0, d.shape[0], step):
            hi = min(lo + step, d.shape[0])
            tot += np.bincount(idx[lo:hi], weights=d[lo:hi].astype(np.float64), minlength=n)
    return genes, np.log1p(TS * tot / tot.sum())


def basal_block(profiles: dict[str, tuple[np.ndarray, np.ndarray]]):
    """Raw and mean-removed basal cosines on the pairwise-shared gene axis."""
    names = list(profiles)
    raw = {}
    cen = {}
    # a common axis for the centring step, so "the mean profile" is one object
    common = set(profiles[names[0]][0])
    for k in names[1:]:
        common &= set(profiles[k][0])
    common = sorted(common)
    M = np.vstack([profiles[k][1][[dict((g, i) for i, g in enumerate(profiles[k][0]))[g]
                                   for g in common]] for k in names])
    Mc = M - M.mean(0, keepdims=True)
    for i, j in itertools.combinations(range(len(names)), 2):
        a, b = M[i], M[j]
        raw[f"{names[i]}|{names[j]}"] = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
        a, b = Mc[i], Mc[j]
        cen[f"{names[i]}|{names[j]}"] = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
    return common, raw, cen


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
        print(f"  {k:<18} {len(d['perts']):>6} perts x {len(d['genes']):>6} genes   {lbl}")

    # ---------- B. effect transfer, with the known-answer gate first ----------
    rows = {}
    for a, b in itertools.combinations(sorted(loaded), 2):
        r = pair(loaded[a], loaded[b], rng)
        if r is not None:
            rows[(a, b)] = r

    print(f"\n{'='*94}\nKNOWN-ANSWER ROW (159 section 4), {KNOWN_DRAWS} independent target draws each")
    bad, known_report = [], {}
    for key, expect in KNOWN.items():
        a, b = key
        draws = [pair(loaded[a], loaded[b], np.random.default_rng(7000 + i))
                 for i in range(KNOWN_DRAWS)]
        got = np.array([[d["full"]["same"], d["full"]["diff"], d["resid"]["same"]]
                        for d in draws])
        mu, sd = got.mean(0), got.std(0, ddof=1)
        tol = np.maximum(KNOWN_TOL_SD * sd, KNOWN_TOL_FLOOR)
        ok = bool(np.all(np.abs(mu - np.array(expect)) <= tol))
        known_report[f"{a}|{b}"] = dict(published=list(expect), mean=mu.tolist(),
                                        sd=sd.tolist(), tol=tol.tolist(), ok=ok)
        print(f"  {a:<14} x {b:<14}")
        for nm, i in (("same", 0), ("diff", 1), ("resid same", 2)):
            print(f"      {nm:<11} 159 {expect[i]:+.4f}   here {mu[i]:+.4f} +- {sd[i]:.4f}"
                  f"   |d| {abs(mu[i]-expect[i]):.4f} vs tol {tol[i]:.4f}")
        print(f"      -> {'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(key)
    if bad:
        raise SystemExit(f"\nABORT: known-answer row did not reproduce: {bad}")
    print("  known-answer row reproduced; the rest of this file is readable.\n")

    print(f"{'='*118}\nEFFECT TRANSFER -- pairs involving the two NEW lines")
    print(f"{'line A':<18}{'line B':<18}{'n_pert':>7}{'n_gene':>8}{'same':>9}{'diff':>9}"
          f"{'RATIO':>8}{'resid same':>12}{'resid null':>12}")
    new_rows = []
    for (a, b), r in sorted(rows.items(), key=lambda kv: -kv[1]["full"]["ratio"]):
        if a not in NEW and b not in NEW:
            continue
        f, s = r["full"], r["resid"]
        print(f"{a:<18}{b:<18}{r['n_perts']:7d}{r['n_genes']:8d}{f['same']:+9.4f}"
              f"{f['diff']:+9.4f}{f['ratio']:8.2f}{s['same']:+12.4f}"
              f"{s['same_permuted']:+12.4f}")
        new_rows.append(dict(A=a, B=b, **r))

    # ---------- A. basal profiles ----------
    profiles = {k: (np.asarray(v["genes"]), v["b_ctrl"])
                for k, v in loaded.items() if v["b_ctrl"] is not None}
    for c in "ABC":
        g, b = context_b_ctrl(REPO / f"data/challenge/context_{c}.h5ad")
        profiles[f"ctx_{c}"] = (g, b)
        print(f"\n  context_{c}: {len(g)} genes, b_ctrl loaded")
    common, raw, cen = basal_block(profiles)
    print(f"\n{'='*78}\nBASAL PROFILE COSINE on {len(common)} shared genes "
          f"(mean-removed across the {len(profiles)} lines)")
    order = [k for k in profiles if k not in ("ctx_A", "ctx_B", "ctx_C")]
    print(f"{'':<18}" + "".join(f"{c:>10}" for c in ("ctx_A", "ctx_B", "ctx_C")) + "   |  vs NEW")
    for k in order:
        cells = []
        for c in ("ctx_A", "ctx_B", "ctx_C"):
            v = cen.get(f"{k}|{c}", cen.get(f"{c}|{k}"))
            cells.append(f"{v:+10.4f}")
        extra = []
        for n in NEW:
            if n == k or n not in profiles:
                continue
            v = cen.get(f"{k}|{n}", cen.get(f"{n}|{k}"))
            extra.append(f"{n.split('_')[0][:6]}={v:+.4f}")
        print(f"{k:<18}" + "".join(cells) + "   |  " + "  ".join(extra))

    (HERE / "work/independence.json").write_text(json.dumps(dict(
        seed=SEED, known_answer=known_report,
        n_common_basal_genes=len(common),
        basal_cosine_raw=raw, basal_cosine_meanremoved=cen,
        effect_pairs=[{"A": a, "B": b, **r} for (a, b), r in rows.items()],
    ), indent=2))
    print("\nwrote work/independence.json")


if __name__ == "__main__":
    main()

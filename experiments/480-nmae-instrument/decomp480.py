"""480 step 2: where does our remaining nmae gap live? Oracle decomposition by gene class.

nmae = mean_p sum_g |lhat - l| / D_p with D_p = sum_g |l| over target p's gate, so it is additive over any partition
of the gate cells within a target: nmae = sum_C mean_p sum_{g in C} e_pg / D_p. Each class's contribution below is that
term; the contributions add to the surface's nmae exactly (checked).

Per class C (cells of the scored targets' gates):
  mass       share of the gate's |l| mass: mean_p sum_{g in C} |l| / D_p (the zero prediction's contribution)
  ours       our contribution to nmae
  acc        |l|-weighted gate sign accuracy of our sign (a zero prediction scores 1/2)
  rec_mag    recoverable with ORACLE MAGNITUDES at our signs (lhat -> s*|l|): ours - contribution of s*|l|
  rec_sign   recoverable with ORACLE SIGNS at our magnitudes (lhat -> sign(l)*|lhat|)
  rec_abst   recoverable if we could ABSTAIN exactly where our sign is wrong (lhat -> 0 there)
  rec_all    recoverable with both oracles (= ours: the contribution goes to 0)
  best_c     in-sample best class-constant MDE-law add on the base (base + s*c*min(1/sqrt(cnt),1)); its gain vs the
             shipped add (c 0.1) is the most a class-based magnitude rule of that shape can buy (upper bound, in-sample)
Our sign s = sign(ours); where ours is 0, the donor-consensus sign, then +1 (434's convention).

Classes: called by us (our q < 0.05) or not; truth direction; control-expression tertile (count per cell, pooled over
the gate); cis position (TSS distance to the target's TSS, GENCODE v44 from 425); target strength (467's oracle-free
design strength n_gt_0.2, and the oracle gate size n_real); our rank (by q then |lfc|) head/body/uncalled.

    uv run python decomp480.py      # -> work/decomp480.json, work/decomp480-cells.npz
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common480 import HERE, ROOT, load_donor_rows, load_h1, nmae

import fold480 as F

TSS = ROOT / "experiments/425-crispri-cis-neighbours/work/tss_v44.tsv"


def tss_classes(targets, genes):
    """Per (target, gene) cis class: 0 <10 kb, 1 10-100 kb, 2 100 kb-1 Mb, 3 same chrom > 1 Mb, 4 other chrom,
    5 unknown TSS."""
    t = pd.read_csv(TSS, sep="\t").set_index("gene_name")
    chrom = t["chrom"].to_dict()
    pos = t["tss"].to_dict()
    gc = np.array([chrom.get(g, "") for g in genes])
    gp = np.array([pos.get(g, -1) for g in genes], dtype=np.int64)
    out = np.full((len(targets), len(genes)), 5, dtype=np.int8)
    for k, tg in enumerate(targets):
        if tg not in chrom:
            continue
        same = (gc == chrom[tg]) & (gp >= 0)
        d = np.abs(gp - pos[tg])
        cls = np.where(d < 10_000, 0, np.where(d < 100_000, 1, np.where(d < 1_000_000, 2, 3)))
        out[k] = np.where(same, cls, np.where(gp >= 0, 4, 5))
    return out


CIS_NAMES = ["<10 kb", "10-100 kb", "0.1-1 Mb", "same chrom >1 Mb", "other chrom", "no TSS"]


def contributions(E, D, G, cls, n_cls):
    """Per-class contribution: mean over scored targets of sum_{g in C} E / D_p. E, cls are targets x genes."""
    scored = G.any(1)
    out = np.zeros(n_cls)
    for c in range(n_cls):
        m = G & (cls == c)
        out[c] = np.mean((np.where(m, E, 0.0).sum(1) / D)[scored])
    return out


def decompose(T, pred, s, base, cnt, classes, title):
    """Returns {class_var: [row dicts]} with the oracle decomposition."""
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    D = np.abs(Tz).sum(1)
    D[D == 0] = np.nan
    E0 = np.abs(pred - Tz)
    Emag = np.abs(s * np.abs(Tz) - Tz)
    Esign = np.abs(np.sign(Tz) * np.abs(pred) - Tz)
    right = np.sign(pred) == np.sign(Tz)
    Eabst = np.where(right, E0, np.abs(Tz))
    hit = np.where(pred == 0, 0.5, right.astype(float))
    add1 = np.minimum(1 / np.sqrt(np.maximum(cnt, 1e-3)), 1.0)[None, :] * (cnt >= 0.3)[None, :]
    total = float(np.nanmean(np.where(G, E0, 0).sum(1)[G.any(1)] / D[G.any(1)]))
    res = {"_total": total}
    grid = np.linspace(0, 0.6, 61)
    for var, (cls, names) in classes.items():
        n = len(names)
        rows = []
        mass = contributions(np.abs(Tz), D, G, cls, n)
        ours = contributions(E0, D, G, cls, n)
        cm = contributions(Emag, D, G, cls, n)
        cs = contributions(Esign, D, G, cls, n)
        ca = contributions(Eabst, D, G, cls, n)
        # class-constant MDE-law add on the base: scan c per class (separable across classes)
        best = np.zeros(n)
        bestc = np.zeros(n)
        shipped = contributions(np.abs(base + s * 0.1 * add1 - Tz), D, G, cls, n)
        scan = np.stack([contributions(np.abs(base + s * c * add1 - Tz), D, G, cls, n) for c in grid])
        for c in range(n):
            j = int(np.argmin(scan[:, c]))
            best[c], bestc[c] = scan[j, c], grid[j]
        assert abs(ours.sum() - total) < 1e-9, (ours.sum(), total)
        for c in range(n):
            m = G & (cls == c)
            if not m.any():
                continue
            w = np.abs(Tz[m])
            rows.append(dict(cls=names[c], cells=int(m.sum()), mass=float(mass[c]), ours=float(ours[c]),
                             acc=float((hit[m] * w).sum() / w.sum()), acc_unw=float(hit[m].mean()),
                             med_l=float(np.median(w)), med_lhat=float(np.median(np.abs(pred[m]))),
                             rec_mag=float(ours[c] - cm[c]), rec_sign=float(ours[c] - cs[c]),
                             rec_abst=float(ours[c] - ca[c]), rec_all=float(ours[c]),
                             shipped=float(shipped[c]), best_c=float(bestc[c]), best_gain=float(best[c] - shipped[c])))
        res[var] = rows
    print(f"\n### {title}: nmae {total:.4f}")
    for var, rows in res.items():
        if var.startswith("_"):
            continue
        print(f"\n{var}")
        print(f"{'class':18s} {'cells':>7s} {'mass':>6s} {'ours':>7s} {'acc':>6s} {'med|l|':>6s} {'med|l^|':>7s} "
              f"{'recMag':>7s} {'recSign':>7s} {'recAbst':>7s} {'best_c':>6s} {'bestGain':>8s}")
        for r in rows:
            print(f"{r['cls']:18s} {r['cells']:7d} {r['mass']:6.3f} {r['ours']:7.4f} {r['acc']:6.3f} {r['med_l']:6.3f} "
                  f"{r['med_lhat']:7.3f} {r['rec_mag']:7.4f} {r['rec_sign']:7.4f} {r['rec_abst']:7.4f} "
                  f"{r['best_c']:6.2f} {r['best_gain']:+8.4f}")
    return res


def h1_inputs():
    h = load_h1()
    T, cnt, genes, targets = h["T"], h["cnt"], h["genes"], h["targets"]
    st = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
    Lc = np.nan_to_num(np.asarray(st["lfc"], float))
    qc = np.asarray(st["q"], float)
    L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"],
                                  float))
    cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0] for d in F.DON]), 0))
    return h, T, cnt, genes, targets, Lc, qc, L0, cons


def main():
    h, T, cnt, genes, targets, Lc, qc, L0, cons = h1_inputs()
    G = np.isfinite(T)
    s = np.sign(Lc)
    s = np.where(s == 0, np.sign(cons), s)
    s = np.where(s == 0, 1.0, s)
    print(f"KA: ours coarse nmae {nmae(Lc, T):.5f} (expect 0.94270)")
    assert abs(nmae(Lc, T) - 0.94270) < 5e-5
    nT = len(targets)
    # classes
    called = np.where(qc < 0.05, 0, 1).astype(np.int8)
    updown = np.where(np.nan_to_num(T) > 0, 0, 1).astype(np.int8)
    q1, q2 = np.quantile(np.broadcast_to(cnt, T.shape)[G], [1 / 3, 2 / 3])
    expr = np.broadcast_to(np.where(cnt < q1, 0, np.where(cnt < q2, 1, 2)).astype(np.int8), T.shape)
    cis = tss_classes(targets, genes)
    stj = json.load(open(ROOT / "experiments/467-strength-gated-programme/work/strength467.json"))
    order = stj["targets_in_order"]
    assert list(order) == list(targets), "467 target order"
    strength = np.array([stj["per_target"][str(i)]["n_gt_0.2"] for i in range(nT)])
    s_cls = np.where(strength <= 1, 0, np.where(strength <= 5, 1, 2)).astype(np.int8)
    nreal = G.sum(1)
    n_cls = np.digitize(nreal, [100, 500, 2000]).astype(np.int8)
    # our rank by q then |lfc| (434e)
    rank = np.zeros(T.shape, dtype=np.int32)
    for k in range(nT):
        o = np.lexsort((-np.abs(Lc[k]), qc[k]))
        rank[k, o] = np.arange(len(o))
    rk = np.where(rank < 200, 0, np.where(rank < 1000, 1, np.where(rank < 3000, 2, 3))).astype(np.int8)
    classes = {
        "called by us (q<0.05)": (called, ["called", "uncalled"]),
        "truth direction": (updown, ["up", "down"]),
        "control expression tertile": (expr, [f"low (<{q1:.2f}/cell)", f"mid", f"high (>={q2:.2f})"]),
        "cis position (TSS distance)": (cis, CIS_NAMES),
        "design strength (467 n_gt_0.2)": (np.broadcast_to(s_cls[:, None], T.shape), ["0-1", "2-5", ">=6"]),
        "gate size n_real (oracle)": (np.broadcast_to(n_cls[:, None], T.shape), ["10-99", "100-499", "500-1999",
                                                                                 ">=2000"]),
        "our rank (q, |lfc|)": (rk, ["0-199", "200-999", "1000-2999", ">=3000"]),
        "our predicted direction": (np.where(s > 0, 0, 1).astype(np.int8), ["we say up", "we say down"]),
        "donor consensus agrees with our sign": (np.where(cons == 0, 2, np.where(np.sign(cons) == s, 0, 1))
                                                 .astype(np.int8), ["agree", "disagree", "no donor value"]),
        "donor |cons| tertile (gate-pooled)": (np.digitize(np.abs(cons), np.quantile(np.abs(cons)[G], [1 / 3, 2 / 3]))
                                               .astype(np.int8), ["small", "mid", "large"]),
    }
    res = {"h1_ours": decompose(T, Lc, s, L0, cnt, classes, "H1, ours (shipped coarse)")}
    # cross-tab: called x expression, called x direction
    cross = {}
    for nm, (a, an), (b, bn) in [("called x expression", classes["called by us (q<0.05)"],
                                  classes["control expression tertile"]),
                                 ("called x direction", classes["called by us (q<0.05)"], classes["truth direction"]),
                                 ("called x our direction", classes["called by us (q<0.05)"],
                                  classes["our predicted direction"]),
                                 ("our direction x donor agreement", classes["our predicted direction"],
                                  classes["donor consensus agrees with our sign"]),
                                 ("our direction x expression", classes["our predicted direction"],
                                  classes["control expression tertile"])]:
        cc = (a.astype(int) * len(bn) + b).astype(np.int8)
        names = [f"{x} / {y}" for x in an for y in bn]
        cross[nm] = (cc, names)
    res["h1_ours_cross"] = decompose(T, Lc, s, L0, cnt, cross, "H1, ours, cross-tabs")
    json.dump(res, open(HERE / "work/decomp480.json", "w"), indent=1)
    np.savez_compressed(HERE / "work/decomp480-cells.npz", called=called, expr=expr, cis=cis, strength=strength,
                        nreal=nreal, rank=rank)


if __name__ == "__main__":
    main()

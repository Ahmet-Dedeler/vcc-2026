"""285 -- the SECOND readout: the scorer's own per-target `pds` (218's `_per_target`, verbatim; 273 KA-2:
== cell-eval2's `pds_cosine` to 5.6e-17) and the aggregate cosine, for the headline contrasts, on the two
sources whose cohorts are large enough for a rank statistic to mean anything.

View A is the prediction, view B is the truth, every target of that cell is a distractor, midrank ties; the
statistic is taken on the UN-CENTRED effect (profile minus the view's own control), which is what the scorer
reads.  Reported as the mean of pds(A->B) and pds(B->A).

    PYTHONDONTWRITEBYTECODE=1 uv run python pds285.py
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

P = L.Logger(L.WORK / "pds285.log")
NDRAW = 4


def summarise(v, tag, n_t, out, cells):
    out[tag] = dict(n=int(n_t), cells_per_side=cells, pds_mean=float(np.nanmean(v["pds"])),
                    pds_median=float(np.nanmedian(v["pds"])), cos_agg=float(v["cos_agg"]),
                    cos_own_median=float(np.nanmedian(v["cos_own"])))
    P(f"    {tag:38s} n_t {n_t:4d}  {cells:5.0f} cells/side  mean pds {out[tag]['pds_mean']:.4f}  "
      f"cos_agg {out[tag]['cos_agg']:.4f}  own-cosine median {out[tag]['cos_own_median']:+.4f}")


def main() -> None:
    t0 = time.time()
    P(f"285 pds285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    out: dict = {}

    # ---------------------------------------------------------------- K562-GWPS panel
    m = np.load(L.WORK / "k562_gwps_meta.npz", allow_pickle=True)
    genes = m["genes"].astype(str); pan_t = m["panel_targets"].astype(str)
    pert = m["pan_pert"].astype(str); batch = m["pan_batch"].astype(str); guide = m["pan_guide"].astype(str)
    ctrl = np.asarray(m["ctrl"], dtype=np.float64)
    gpos = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(len(genes))) - {gpos[t] for t in pan_t if t in gpos}), dtype=np.int64)
    idx_t = {t: np.flatnonzero(pert == t) for t in pan_t}
    batches = sorted(set(batch.tolist())); bpar = np.asarray([batches.index(b) % 2 for b in batch])
    Xd = np.asarray(np.load(L.WORK / "k562_gwps_panelX.npy", mmap_mode="r"), dtype=np.float32)

    def sums(sel):
        return Xd[np.sort(sel)].sum(0, dtype=np.float64)

    P("\nK562-GWPS, 272 panel targets:")
    out["k562"] = {}
    for n in (25, 50):
        acc = defaultdict(list)
        for d in range(NDRAW):
            r = np.random.default_rng(L.SEED + 11 * d + n)
            keep, SA, SB, GA, GB, k2 = [], [], [], [], [], []
            for j, t in enumerate(pan_t):
                ix = idx_t[t]
                if ix.size >= 2 * n:
                    p = r.permutation(ix)
                    keep.append(t); SA.append(sums(p[:n])); SB.append(sums(p[n:2 * n]))
                a, b = ix[bpar[ix] == 0], ix[bpar[ix] == 1]
                if a.size >= n and b.size >= n:
                    k2.append(t); GA.append(sums(r.choice(a, n, replace=False)))
                    GB.append(sums(r.choice(b, n, replace=False)))
            acc["cells"].append(L.views(np.stack(SA), np.stack(SB), ctrl[0], ctrl[0], scored))
            acc["gemgroup"].append(L.views(np.stack(GA), np.stack(GB), ctrl[1], ctrl[2], scored))
            acc["n_cells"].append(len(keep)); acc["n_gem"].append(len(k2))
        for tag in ("cells", "gemgroup"):
            v = dict(pds=np.nanmean([x["pds"] for x in acc[tag]], 0),
                     cos_own=np.nanmean([x["cos_own"] for x in acc[tag]], 0),
                     cos_agg=float(np.mean([x["cos_agg"] for x in acc[tag]])))
            summarise(v, f"n={n} {tag}", acc[f"n_{'cells' if tag == 'cells' else 'gem'}"][0], out["k562"], n)

    # the across-study rows (blocks are on disk)
    P("\nK562 GWPS vs the essential screen, shared perturbations:")
    ma, mb = m, np.load(L.WORK / "k562_essential_meta.npz", allow_pickle=True)
    ga_, gb_ = ma["genes"].astype(str), mb["genes"].astype(str)
    shared = ma["shared"].astype(str)
    common_g = sorted(set(ga_.tolist()) & set(gb_.tolist()))
    pa = {g: i for i, g in enumerate(ga_)}; pb = {g: i for i, g in enumerate(gb_)}
    ia = np.asarray([pa[g] for g in common_g]); ib = np.asarray([pb[g] for g in common_g])
    tset = set(shared.tolist())
    keep_col = np.flatnonzero([g not in tset for g in common_g])
    BA = np.load(L.WORK / "k562_gwps_study_blocks.npy", mmap_mode="r")
    BB = np.load(L.WORK / "k562_essential_study_blocks.npy", mmap_mode="r")
    nA, nB = np.asarray(ma["ncell_shared"]), np.asarray(mb["ncell_shared"])
    cA = np.asarray(ma["ctrl"], dtype=np.float64)[0][ia]; cB = np.asarray(mb["ctrl"], dtype=np.float64)[0][ib]
    EDG = (0, 25, 50, 100)
    out["k562_study"] = {}
    for n in (25, 50):
        cohort = np.flatnonzero((nA >= 2 * n) & (nB >= 2 * n))
        acc = defaultdict(list)
        for d in range(BA.shape[0]):
            def side(Bk, ii, lo, hi):
                S = np.zeros((cohort.size, len(common_g)))
                for b in range(len(EDG) - 1):
                    if EDG[b] >= lo and EDG[b + 1] <= hi:
                        S += np.asarray(Bk[d, cohort * (len(EDG) - 1) + b][:, ii], dtype=np.float64)
                return S
            A1, A2 = side(BA, ia, 0, n), side(BA, ia, n, 2 * n)
            B1, B2 = side(BB, ib, 0, n), side(BB, ib, n, 2 * n)
            acc["within_gwps"].append(L.views(A1, A2, cA, cA, keep_col))
            acc["across_study"].append(L.views(A1, B1, cA, cB, keep_col))
        for tag in ("within_gwps", "across_study"):
            v = dict(pds=np.nanmean([x["pds"] for x in acc[tag]], 0),
                     cos_own=np.nanmean([x["cos_own"] for x in acc[tag]], 0),
                     cos_agg=float(np.mean([x["cos_agg"] for x in acc[tag]])))
            summarise(v, f"n={n} {tag}", cohort.size, out["k562_study"], n)
    del Xd

    # ---------------------------------------------------------------- Nadig Jurkat
    P("\nNadig Jurkat (challenge line A's tissue), 534 targets used:")
    import h5py

    def col(h, frame, key):
        g = h[frame]
        if "__categories" in g and key in g["__categories"]:
            cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["__categories"][key][:]])
            return cats[g[key][:]]
        v = g[key][:]
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])

    z = np.load(L.WORK / "nadig_jurkat.npz", allow_pickle=True)
    use = list(z["targets"].astype(str)); multi = list(z["multi"].astype(str))
    PATH = L.ROOT / "data/external/nadig2025/jurkat_raw_singlecell.h5ad"
    with h5py.File(PATH, "r") as h:
        gene = col(h, "obs", "gene"); sg = col(h, "obs", "sgID_AB"); gg = col(h, "obs", "gem_group")
        gnames = col(h, "var", "gene_name")
        n_cells, n_genes = h["X"].shape
        ctrl_lab = [t for t in set(gene.tolist()) if t.lower().startswith("non-target")]
        keep = np.isin(gene, use) | np.isin(gene, ctrl_lab)
        blocks = []
        for r0 in range(0, n_cells, 16384):
            r1 = min(r0 + 16384, n_cells)
            mm = keep[r0:r1]
            if mm.any():
                blocks.append(np.asarray(h["X"][r0:r1][mm], dtype=np.float32))
        Xn = np.concatenate(blocks); del blocks
    idx = np.flatnonzero(keep)
    gene, sg, gg = gene[idx], sg[idx], gg[idx]
    is_ctrl = np.isin(gene, ctrl_lab)
    cfull = Xn[is_ctrl].sum(0, dtype=np.float64)
    gp = {g: i for i, g in enumerate(gnames)}
    sc = np.asarray(sorted(set(range(n_genes)) - {gp[t] for t in use if t in gp}), dtype=np.int64)
    it = {t: np.flatnonzero(gene == t) for t in use}
    gems = sorted(set(gg.tolist())); gpar = np.asarray([gems.index(x) % 2 for x in gg])
    out["nadig_jurkat"] = {}
    n = 25
    acc = defaultdict(list)
    for d in range(NDRAW):
        r = np.random.default_rng(L.SEED + 11 * d + n)
        SA, SB, GA, GB, HA, HB, CA_, CB_ = [], [], [], [], [], [], [], []
        nk, ng, nh = 0, 0, 0
        for t in use:
            ix = it[t]
            if ix.size >= 2 * n:
                p = r.permutation(ix); SA.append(Xn[np.sort(p[:n])].sum(0, dtype=np.float64))
                SB.append(Xn[np.sort(p[n:2 * n])].sum(0, dtype=np.float64)); nk += 1
            a, b = ix[gpar[ix] == 0], ix[gpar[ix] == 1]
            if a.size >= n and b.size >= n:
                GA.append(Xn[np.sort(r.choice(a, n, replace=False))].sum(0, dtype=np.float64))
                GB.append(Xn[np.sort(r.choice(b, n, replace=False))].sum(0, dtype=np.float64)); ng += 1
            if t in multi:
                cs = sorted([ix[sg[ix] == g] for g in sorted(set(sg[ix].tolist()))], key=lambda a_: -a_.size)
                if cs[0].size >= 2 * n and cs[1].size >= n:
                    pb = r.permutation(cs[0])
                    HA.append(Xn[np.sort(pb[:n])].sum(0, dtype=np.float64))
                    HB.append(Xn[np.sort(r.choice(cs[1], n, replace=False))].sum(0, dtype=np.float64))
                    CA_.append(Xn[np.sort(pb[:n])].sum(0, dtype=np.float64))
                    CB_.append(Xn[np.sort(pb[n:2 * n])].sum(0, dtype=np.float64)); nh += 1
        acc["cells"].append(L.views(np.stack(SA), np.stack(SB), cfull, cfull, sc))
        acc["gemgroup"].append(L.views(np.stack(GA), np.stack(GB), cfull, cfull, sc))
        acc["guide"].append(L.views(np.stack(HA), np.stack(HB), cfull, cfull, sc))
        acc["cells_one_construct"].append(L.views(np.stack(CA_), np.stack(CB_), cfull, cfull, sc))
        acc["n"] = (nk, ng, nh)
    for tag, nt in (("cells", acc["n"][0]), ("gemgroup", acc["n"][1]), ("guide", acc["n"][2]),
                    ("cells_one_construct", acc["n"][2])):
        v = dict(pds=np.nanmean([x["pds"] for x in acc[tag]], 0),
                 cos_own=np.nanmean([x["cos_own"] for x in acc[tag]], 0),
                 cos_agg=float(np.mean([x["cos_agg"] for x in acc[tag]])))
        summarise(v, f"n={n} {tag}", nt, out["nadig_jurkat"], n)

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    L.Logger.dump(P, out, L.WORK / "pds285.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()

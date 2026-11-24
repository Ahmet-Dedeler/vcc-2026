"""424 validation of data/donors/viperturb_k562 (c): panel coverage, knockdown, and cross-screen effect agreement.

KNOWN ANSWER FIRST: 159's `donorpairs.pair()` on (k562_essential, k562_gwps) must reproduce 159's same-line row
(full-space same +0.2712, diff +0.0329; resid same +0.2857) before any VIPerturb number is read.

  1. coverage   panel ∩ perts (>= 20 cells) == 284 (413's count); the 30 panel targets absent from K562-GWPS.
  2. knockdown  own-gene delta (log-pseudobulk, the table's own space) per target; per-guide knockdown from the
                per-cell own-target counts (cells.tsv.gz); the same own-gene statistic for K562-GWPS as the reference.
  3. agreement  159's line-pair statistic (same / diff / ratio / resid same / permutation null) for VIPerturb vs
                K562-GWPS, K562-essential, Orion HCT116 / HEK293T, kolf2, with 159's ladder rows re-measured alongside;
                plus the full per-target cosine distribution over ALL shared targets and the 254 shared PANEL targets.

    uv run python validate424.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import gzip  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/159-transfer-ceiling"))
import donorpairs  # noqa: E402

PANEL = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
DON = {
    "viperturb_k562": "data/donors/viperturb_k562",
    "k562_gwps": "data/donors/k562_gwps",
    "k562_essential": "data/donors/k562_essential",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
    "kolf2": "data/external/kolf2",
}


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30


class Lazy:
    """A donor table memory-mapped; `sub(perts, genes)` materialises only the needed block."""

    def __init__(self, name):
        p = ROOT / DON[name]
        self.name = name
        self.delta = np.load(p / "delta.npy", mmap_mode="r")
        self.perts = [l.strip() for l in (p / "perts.txt").read_text().splitlines() if l.strip()]
        self.genes = [l.strip() for l in (p / "genes.txt").read_text().splitlines() if l.strip()]
        self.pi = {g: i for i, g in enumerate(self.perts)}
        self.gi = {g: i for i, g in enumerate(self.genes)}
        assert self.delta.shape == (len(self.perts), len(self.genes))
        m = p / "meta.json"
        self.cells = json.loads(m.read_text()).get("cells_per_pert", {}).get("median") if m.exists() else None

    def sub(self, perts, genes):
        r = np.asarray([self.pi[p] for p in perts])
        c = np.asarray([self.gi[g] for g in genes])
        order = np.argsort(r)
        blk = np.asarray(self.delta[np.sort(r)], dtype=np.float32)[:, c]
        out = np.empty_like(blk)
        out[order] = blk
        return np.nan_to_num(out)


def as_pairdict(A: Lazy, B: Lazy, rng, cap=donorpairs.MAX_PERT):
    """The (perts, genes, delta) dicts donorpairs.pair() takes, restricted to what it will read: pair() draws its
    3,000-target subsample from the SAME rng as here, so the draw is reproduced before materialising."""
    sh = sorted(set(A.perts) & set(B.perts))
    shg = sorted(set(A.genes) & set(B.genes))
    state = rng.bit_generator.state
    if len(sh) > cap:
        sh_used = [sh[i] for i in sorted(rng.choice(len(sh), cap, replace=False))]
    else:
        sh_used = sh
    rng.bit_generator.state = state              # pair() will redraw the identical subsample
    genes = [g for g in shg if g not in set(sh_used)]
    return (dict(delta=_pad(A, sh, sh_used, genes), perts=sh, genes=genes, cells=A.cells),
            dict(delta=_pad(B, sh, sh_used, genes), perts=sh, genes=genes, cells=B.cells))


def _pad(D, sh, used, genes):
    """Rows for every shared target (pair() indexes by the full shared list); only the used ones are real."""
    M = np.zeros((len(sh), len(genes)), dtype=np.float32)
    idx = {p: i for i, p in enumerate(sh)}
    M[[idx[p] for p in used]] = D.sub(used, genes)
    return M


def percos(A: Lazy, B: Lazy, targets, centre_on=None):
    """Per-target cosine over shared genes minus the target set's coordinates (159's convention), raw and centred.
    Centring subtracts each table's mean row over `centre_on` (default: the targets themselves)."""
    genes = [g for g in sorted(set(A.genes) & set(B.genes)) if g not in set(targets)]
    Xa, Xb = A.sub(targets, genes).astype(np.float64), B.sub(targets, genes).astype(np.float64)
    u = lambda M: M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-300)
    raw = (u(Xa) * u(Xb)).sum(1)
    Ca, Cb = Xa - Xa.mean(0), Xb - Xb.mean(0)
    cen = (u(Ca) * u(Cb)).sum(1)
    S = u(Ca) @ u(Cb).T
    diff = float(S[~np.eye(len(targets), dtype=bool)].mean())
    pds = float(np.mean([(S[i] < S[i, i]).sum() / (len(targets) - 1) for i in range(len(targets))]))
    q = lambda v: dict(mean=float(v.mean()), median=float(np.median(v)), q10=float(np.quantile(v, .1)),
                       q90=float(np.quantile(v, .9)), frac_pos=float((v > 0).mean()))
    return dict(n=len(targets), n_genes=len(genes), raw=q(raw), centred=q(cen), centred_diff=diff,
                centred_pds_retrieval=pds, per_target=dict(zip(targets, cen.round(5).tolist())))


def knockdown(D: Lazy, perts):
    own = [(p, float(D.delta[D.pi[p], D.gi[p]])) for p in perts if p in D.gi]
    v = np.asarray([x for _, x in own])
    return dict(n_with_own_gene=len(own), n_without=len(perts) - len(own), frac_down=float((v < 0).mean()),
                median=float(np.median(v)), q10=float(np.quantile(v, .1)), q90=float(np.quantile(v, .9)),
                frac_below_m0p1=float((v < -0.1).mean())), dict(own)


def main():
    out = {}
    L = []
    say = lambda *a: (print(*a, flush=True), L.append(" ".join(map(str, a))))
    D = {n: Lazy(n) for n in DON}
    rng = np.random.default_rng(donorpairs.SEED)

    # ---------------------------------------------------------------- KA: 159's same-line row
    a, b = as_pairdict(D["k562_essential"], D["k562_gwps"], rng)
    r = donorpairs.pair(a, b, rng)
    say(f"KA-159 k562_essential x k562_gwps: n {r['n_perts']} x {r['n_genes']}; full same {r['full']['same']:+.4f} "
        f"diff {r['full']['diff']:+.4f} ratio {r['full']['ratio']:.2f}; resid same {r['resid']['same']:+.4f} "
        f"(159: 2003 x 6402, +0.2712 / +0.0329 / 8.26; resid +0.2857)")
    assert (r["n_perts"], r["n_genes"]) == (2003, 6402) and abs(r["full"]["same"] - 0.2712) < 5e-5 \
        and abs(r["full"]["diff"] - 0.0329) < 5e-5 and abs(r["resid"]["same"] - 0.2857) < 5e-5, "KA-159 FAILED"
    out["ka159"] = r

    # ---------------------------------------------------------------- 1. coverage
    V = D["viperturb_k562"]
    k562 = set(D["k562_gwps"].perts)
    cov = sorted(set(PANEL) & set(V.perts))
    fill = sorted(set(cov) - k562)
    both = sorted(set(cov) & k562)
    prov = json.loads((ROOT / DON["viperturb_k562"] / "provenance.json").read_text())
    cells = prov["cells_per_target"]
    out["coverage"] = dict(panel_ge20=len(cov), fill=fill, n_fill=len(fill), n_both=len(both),
                           panel_missing=sorted(set(PANEL) - set(V.perts)),
                           panel_cells=dict(median=float(np.median([cells[t] for t in cov])),
                                            min=int(min(cells[t] for t in cov)), max=int(max(cells[t] for t in cov))),
                           fill_cells={t: cells[t] for t in fill})
    say(f"1. coverage: panel at >= 20 cells {len(cov)} (413: 284); shared with K562-GWPS {len(both)}; fill {len(fill)}: "
        f"{fill}; missing {out['coverage']['panel_missing']}; panel cells median {out['coverage']['panel_cells']['median']}")
    assert len(cov) == 284 and len(fill) == 30, "coverage does not reproduce 413"

    # ---------------------------------------------------------------- 2. knockdown
    kv, ownv = knockdown(V, V.perts)
    kvp, ownvp = knockdown(V, cov)
    kg, owng = knockdown(D["k562_gwps"], D["k562_gwps"].perts)
    kgp, _ = knockdown(D["k562_gwps"], both)
    kvb, _ = knockdown(V, both)
    out["knockdown"] = dict(vip_all=kv, vip_panel=kvp, gwps_all=kg, gwps_panel_shared=kgp, vip_panel_shared=kvb)
    for lab, k in (("VIPerturb all", kv), ("VIPerturb panel (284)", kvp), ("VIPerturb panel shared (254)", kvb),
                   ("K562-GWPS all", kg), ("K562-GWPS panel shared (254)", kgp)):
        say(f"2. knockdown {lab:30s}: own gene on axis {k['n_with_own_gene']:5d} (not {k['n_without']}); "
            f"down {100 * k['frac_down']:.1f}%, < -0.1 {100 * k['frac_below_m0p1']:.1f}%; delta median {k['median']:+.3f} "
            f"[q10 {k['q10']:+.3f}, q90 {k['q90']:+.3f}]")
    # per-guide knockdown from per-cell own-target counts
    g_rows = []
    ctrl = None
    for bn in ("binA", "binB", "binC"):
        c = pd.read_csv(WORK / f"conv/{bn}/cells.tsv.gz", sep="\t", usecols=["gene", "guide", "libsize", "target_count"])
        feats = (WORK / f"conv/{bn}/features.txt").read_text().split("\n")
        if ctrl is None:
            S = np.load(WORK / f"conv/{bn}/sums.npy", mmap_mode="r")
            grp = (WORK / f"conv/{bn}/groups.txt").read_text().split("\n")
            nt = np.asarray(S[grp.index("NO-TARGET")])
            ctrl = dict(zip(feats, nt / nt.sum()))
        c = c[(c.gene != "NO-TARGET") & c.target_count.notna()]
        agg = c.groupby(["gene", "guide"]).agg(n=("libsize", "size"), lib=("libsize", "sum"), t=("target_count", "sum"))
        g_rows.append(agg.reset_index())
    G = pd.concat(g_rows)
    G["frac"] = G.t / G.lib
    G["ctrl"] = G.gene.map(ctrl)
    G = G[G.ctrl > 0]
    G["log2kd"] = np.log2((G.frac + 0.5 / G.lib) / G.ctrl)
    Gp = G[G.gene.isin(cov)]
    out["guides"] = dict(n_guides=int(len(G)), median_log2=float(G.log2kd.median()),
                         frac_guides_log2_gt_m0p5=float((G.log2kd > -0.5).mean()),
                         panel_n_guides=int(len(Gp)), panel_median_log2=float(Gp.log2kd.median()),
                         panel_frac_guides_log2_gt_m0p5=float((Gp.log2kd > -0.5).mean()),
                         panel_targets_all_guides_weak=sorted(Gp.groupby("gene").log2kd.min().loc[lambda s: s > -0.5].index))
    say(f"2b. per-guide knockdown (log2 own-target CPM vs NO-TARGET): {len(G)} guides, median {G.log2kd.median():+.2f}, "
        f"{100 * (G.log2kd > -0.5).mean():.1f}% weaker than -0.5; panel {len(Gp)} guides median {Gp.log2kd.median():+.2f}, "
        f"{100 * (Gp.log2kd > -0.5).mean():.1f}% weak; panel targets whose best guide is weaker than -0.5: "
        f"{out['guides']['panel_targets_all_guides_weak']}")
    kd_tab = pd.DataFrame(dict(target=cov, cells=[cells[t] for t in cov], vip_own_delta=[ownvp.get(t, np.nan) for t in cov],
                               gwps_own_delta=[owng.get(t, np.nan) for t in cov], in_k562_gwps=[t in k562 for t in cov]))
    best = Gp.groupby("gene").log2kd.min()
    kd_tab["best_guide_log2kd"] = kd_tab.target.map(best)
    kd_tab.to_csv(HERE / "knockdown.csv", index=False)
    say(f"   rss {rss():.2f} GiB")

    # ---------------------------------------------------------------- 3. agreement (159's statistic)
    rows = []
    pairs = [("viperturb_k562", "k562_gwps"), ("viperturb_k562", "k562_essential"), ("viperturb_k562", "orion_hct116"),
             ("viperturb_k562", "orion_hek293t"), ("viperturb_k562", "kolf2"), ("k562_gwps", "orion_hct116"),
             ("k562_gwps", "orion_hek293t"), ("k562_essential", "orion_hct116")]
    for x, y in pairs:
        a, b = as_pairdict(D[x], D[y], rng)
        r = donorpairs.pair(a, b, rng)
        del a, b
        if r is None:
            continue
        r["A"], r["B"] = x, y
        rows.append(r)
        say(f"3. {x:15s} x {y:15s} n {r['n_perts']:5d} x {r['n_genes']:6d}: full same {r['full']['same']:+.4f} diff "
            f"{r['full']['diff']:+.4f} ratio {r['full']['ratio']:6.2f} | resid same {r['resid']['same']:+.4f} null "
            f"{r['resid']['same_permuted']:+.4f} | within A {r['full']['within_A']:+.4f} B {r['full']['within_B']:+.4f}  "
            f"(rss {rss():.2f} GiB)")
    out["pairs"] = rows

    # ---------------------------------------------------------------- 3b. full distributions
    sh_all = sorted(set(V.perts) & k562)
    out["dist"] = {}
    for lab, tg, other in (("vip_x_gwps_all_shared", sh_all, "k562_gwps"), ("vip_x_gwps_panel254", both, "k562_gwps"),
                           ("ess_x_gwps_all_shared", sorted(set(D["k562_essential"].perts) & k562), "k562_gwps"),
                           ("vip_x_hct116_panel", sorted(set(cov) & set(D["orion_hct116"].perts)), "orion_hct116"),
                           ("gwps_x_hct116_panel", sorted(set(both) & set(D["orion_hct116"].perts)), "orion_hct116")):
        A = D["k562_essential"] if lab.startswith("ess") else (D["k562_gwps"] if lab.startswith("gwps") else V)
        r = percos(A, D[other], tg)
        out["dist"][lab] = r
        say(f"3b. {lab:24s} n {r['n']:5d} genes {r['n_genes']:6d}: raw cos median {r['raw']['median']:+.4f} "
            f"(>0 {100 * r['raw']['frac_pos']:.0f}%) | centred mean {r['centred']['mean']:+.4f} median "
            f"{r['centred']['median']:+.4f} [q10 {r['centred']['q10']:+.4f}, q90 {r['centred']['q90']:+.4f}] >0 "
            f"{100 * r['centred']['frac_pos']:.0f}% | diff {r['centred_diff']:+.4f} | retrieval pds "
            f"{r['centred_pds_retrieval']:.3f}  (rss {rss():.2f} GiB)")
    kd_tab["cos_vs_gwps_centred"] = kd_tab.target.map(out["dist"]["vip_x_gwps_panel254"]["per_target"])
    kd_tab.to_csv(HERE / "knockdown.csv", index=False)
    out["peak_rss_gib"] = rss()
    (WORK / "validate424.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "validate424.log").write_text("\n".join(L) + "\n")
    say(f"done; peak RSS {rss():.2f} GiB")


if __name__ == "__main__":
    main()

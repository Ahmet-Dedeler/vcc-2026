"""284 — KA rows (PREREG §2), the leak quantification (§3) and the split construction (§4).

  PYTHONDONTWRITEBYTECODE=1 uv run python ka284.py

Outputs: work/ka284.json, work/ka284.log, work/splits284.npz
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

import lib284 as L

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    prereg = L.HERE / "PREREGISTRATION.md"
    h = hashlib.sha256(prereg.read_bytes()).hexdigest()
    frozen = (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]
    P(f"PREREGISTRATION.md sha256 {h}")
    P(f"frozen digest            {frozen} -> {'MATCH' if h == frozen else 'MISMATCH'}")
    if h != frozen:
        sys.exit("preregistration digest mismatch — aborting")

    out: dict = {"seed": L.SEED, "prereg_sha256": h, "ka": {}, "leak": {}, "splits": {}}
    held = L.load_held()
    n = len(held)
    P(f"held genes {n}")

    # ================= KA-1: 281's primary cell, on 281's own split ============================
    te281, take281, block = L.draw_281(held)
    p281 = pd.read_csv(L.PRED281, usecols=["gene", "set", "split_family", "split_random"])
    p281 = p281[p281["set"] == "held"].reset_index(drop=True)
    same_order = bool((p281["gene"].to_numpy() == held["gene"].to_numpy()).all())
    ref_te = (p281["split_family"].to_numpy() == "held")
    match = int((te281 == ref_te).sum())
    P(f"KA-1a split: same gene order {same_order}; my family-held test == pred281 split_family on "
      f"{match}/{n} genes; n_test {int(te281.sum())} in {len(set(block[te281]))} blocks "
      f"({len(take281)} taken)")

    y = held["y_raw"].to_numpy(float)
    X = held[L.FEATS].to_numpy(float)
    XO = held[L.RELIABILITY].to_numpy(float)
    tr, te = ~te281, te281
    pr, mr = L.fit_ridge(X[tr], y[tr], X[te])
    pg, mg, _ = L.fit_gbm(X[tr], y[tr], X[te])
    po, mo = L.fit_ridge(XO[tr], y[tr], XO[te])
    got = dict(ridge_rho=L.spearman(y[te], pr), ridge_r2=L.r2(y[te], pr),
               gbm_rho=L.spearman(y[te], pg), rel_rho=L.spearman(y[te], po),
               ridge_alpha=mr["alpha"], gbm_params=mg["params"], rel_alpha=mo["alpha"])
    req = dict(ridge_rho=0.320420, ridge_r2=0.111332, gbm_rho=0.315974, rel_rho=0.351664)
    fit281 = json.loads((L.ROOT / "experiments/281-what-predicts-conservation/work/fit281.json"
                         ).read_text())["results"]
    exact = dict(ridge_rho=fit281["family_raw"]["ridge"]["spearman"],
                 ridge_r2=fit281["family_raw"]["ridge"]["r2"],
                 gbm_rho=fit281["family_raw"]["gbm"]["spearman"],
                 rel_rho=fit281["family_raw_rel_only"]["ridge"]["spearman"])
    devs = {k: abs(got[k] - exact[k]) for k in req}
    ka1 = dict(split_match=match, n_test=int(te281.sum()), n_blocks=len(set(block[te281])),
               got=got, required=req, fit281_exact=exact, max_abs_dev=float(max(devs.values())),
               deviations=devs,
               pass_=bool(match == n and int(te281.sum()) == 3422
                          and len(set(block[te281])) == 481 and max(devs.values()) == 0.0))
    for k in req:
        P(f"  KA-1 {k:12s} mine {got[k]:+.6f}  281 {exact[k]:+.6f}  |Δ| {devs[k]:.1e}")
    P(f"  KA-1 -> {'PASS' if ka1['pass_'] else 'FAIL'}")
    out["ka"]["KA-1"] = ka1
    if not ka1["pass_"]:
        sys.exit("KA-1 failed — aborting")

    # ================= KA-2: c_held re-derived over all 364 pairs ==============================
    z = np.load(L.HELD279, allow_pickle=True)
    Ccross = np.asarray(z["Ccross"])
    c_held = np.asarray(z["c_held"], dtype=np.float64)
    with np.errstate(invalid="ignore"):
        mine = L.cmap(Ccross, np.arange(Ccross.shape[0]))          # the half-map code, all pairs
        mine64 = np.nanmean(Ccross.astype(np.float64), axis=0)     # float64, for scale only
    d = float(np.nanmax(np.abs(mine - c_held)))
    d64 = float(np.nanmax(np.abs(mine64 - c_held)))
    genes279 = np.asarray(z["genes"]).astype(str)
    same_genes = bool((genes279 == held["gene"].to_numpy().astype(str)).all())
    ka2 = dict(max_abs_diff=d, max_abs_diff_float64=d64, n=int(np.isfinite(mine).sum()),
               median=float(np.nanmedian(mine)), same_gene_order_as_feat281=same_genes,
               pass_=bool(d == 0.0 and int(np.isfinite(mine).sum()) == 17106 and same_genes))
    P(f"KA-2 c_held over all 364 pairs: max|Δ| {d:.1e} (float64 recompute {d64:.1e} — 279 stored "
      f"c_held as float32), n {ka2['n']}, median {ka2['median']:+.6f}, "
      f"gene order matches feat281 {same_genes} -> {'PASS' if ka2['pass_'] else 'FAIL'}")
    out["ka"]["KA-2"] = ka2
    if not ka2["pass_"]:
        sys.exit("KA-2 failed — aborting")

    # ================= KA-3: the balanced core off znf279.csv ==================================
    c3 = held["c3_balanced_core"].to_numpy(float)
    ok3 = np.isfinite(c3)
    isz = held["znf_class"].to_numpy() == "ZNF"
    isn = held["znf_class"].to_numpy() == "NONZNF"
    ka3 = dict(n=int(ok3.sum()), znf_median=float(np.median(c3[isz & ok3])),
               nonznf_median=float(np.median(c3[isn & ok3])))
    ka3["d_median"] = ka3["znf_median"] - ka3["nonznf_median"]
    ka3["pass_"] = bool(ka3["n"] == 8541 and abs(ka3["znf_median"] - 0.011813) < 5e-6
                        and abs(ka3["nonznf_median"] - 0.029663) < 5e-6
                        and abs(ka3["d_median"] + 0.017850) < 5e-6)
    P(f"KA-3 balanced core: n {ka3['n']}, ZNF {ka3['znf_median']:+.6f}, non-ZNF "
      f"{ka3['nonznf_median']:+.6f}, Δ {ka3['d_median']:+.6f} -> "
      f"{'PASS' if ka3['pass_'] else 'FAIL'}")
    out["ka"]["KA-3"] = ka3
    if not ka3["pass_"]:
        sys.exit("KA-3 failed — aborting")

    # ================= §3 the leak, quantified ==================================================
    mem = L.memberships(held)
    nmem = np.array([len(m) for m in mem])
    out["leak"]["membership_counts"] = {str(k): int(v) for k, v in
                                        sorted(Counter(nmem.tolist()).items())}
    out["leak"]["n_multi_family"] = int((nmem >= 2).sum())
    P(f"\n§3 LEAK. Held genes with >= 2 HGNC groups: {int((nmem >= 2).sum())} "
      f"({(nmem >= 2).mean():.1%}); 1 group {int((nmem == 1).sum())}; none {int((nmem == 0).sum())}")

    take_fams = {t for t in take281 if not t.startswith("NOFAM:")}
    train_fams: set[str] = set()
    for i in np.flatnonzero(~te281):
        train_fams.update(mem[i])
    # (b) TEST genes sharing ANY family with a TRAIN gene, under 281's split
    b_cnt = int(sum(1 for i in np.flatnonzero(te281) if set(mem[i]) & train_fams))
    # (c) TRAIN genes carrying a HELD block as a non-first group
    c_cnt = int(sum(1 for i in np.flatnonzero(~te281)
                    if set(mem[i][1:]) & take_fams))
    # how many distinct held families leak into training rows
    leaked_fams = sorted({f for i in np.flatnonzero(~te281) for f in mem[i][1:] if f in take_fams})
    out["leak"]["b_test_genes_sharing_family_with_train"] = b_cnt
    out["leak"]["b_frac_of_test"] = float(b_cnt / int(te281.sum()))
    out["leak"]["c_train_genes_carrying_held_family_non_first"] = c_cnt
    out["leak"]["n_held_families_leaking"] = len(leaked_fams)
    out["leak"]["n_held_families"] = len(take_fams)
    out["leak"]["leaked_families_examples"] = leaked_fams[:25]
    P(f"  (b) under 281's split, TEST genes sharing any HGNC family with a TRAIN gene: {b_cnt} "
      f"of {int(te281.sum())} ({b_cnt / int(te281.sum()):.1%})")
    P(f"  (c) TRAIN genes carrying one of the {len(take_fams)} held families as a NON-FIRST group: "
      f"{c_cnt}")
    P(f"      distinct held families with at least one training row: {len(leaked_fams)} of "
      f"{len(take_fams)}")

    # ================= §4 the splits ============================================================
    teS1, takeS1, _ = L.draw_s1(held)
    comp = L.components(held)
    ncomp = int(comp.max() + 1)
    sizes = np.bincount(comp)
    teS2 = L.draw_components(comp, n)
    teRND = L.draw_random(n)
    fams_S1 = {t for t in takeS1 if not t.startswith("NOFAM:")}

    def leak_check(te: np.ndarray, held_fams: set[str] | None = None) -> dict:
        trf: set[str] = set()
        for i in np.flatnonzero(~te):
            trf.update(mem[i])
        shared = int(sum(1 for i in np.flatnonzero(te) if set(mem[i]) & trf))
        d = dict(n_test=int(te.sum()), frac=float(te.mean()),
                 test_genes_sharing_family_with_train=shared)
        if held_fams is not None:
            d["train_genes_carrying_a_held_family"] = int(
                sum(1 for i in np.flatnonzero(~te) if set(mem[i]) & held_fams))
        # the components the test genes fall into, and the giant one's share of them
        ct = comp[te]
        big = int(np.argmax(np.bincount(comp)))
        d["n_components_in_test"] = int(len(set(ct.tolist())))
        d["giant_component_share_of_test"] = float((ct == big).mean())
        return d

    lk281 = leak_check(te281, take_fams)
    lkS1 = leak_check(teS1, fams_S1)
    lkS2 = leak_check(teS2)
    lkR = leak_check(teRND)
    P(f"\n§4 SPLITS")
    for nm, lk in (("281 (first group)", lk281), ("S1 (all membership)", lkS1),
                   ("S2 (component)", lkS2), ("random (20260907)", lkR)):
        P(f"  {nm:20s} test {lk['n_test']} ({lk['frac']:.3%}); test genes sharing ANY family with "
          f"a train gene {lk['test_genes_sharing_family_with_train']}; train genes carrying a HELD "
          f"family {lk.get('train_genes_carrying_a_held_family', 'n/a')}; test spans "
          f"{lk['n_components_in_test']} components, giant "
          f"{lk['giant_component_share_of_test']:.1%} of them")
    P(f"  components: {ncomp}; largest {int(sizes.max())} genes ({sizes.max() / n:.2%} of held); "
      f"size distribution 1:{int((sizes == 1).sum())} 2-4:{int(((sizes >= 2) & (sizes <= 4)).sum())} "
      f"5-19:{int(((sizes >= 5) & (sizes <= 19)).sum())} 20-99:"
      f"{int(((sizes >= 20) & (sizes <= 99)).sum())} >=100:{int((sizes >= 100).sum())}")
    top = np.argsort(-sizes)[:10]
    P(f"  10 largest components (size, an example gene, n distinct families): " + "; ".join(
        f"{int(sizes[c])}/{held['gene'].to_numpy()[np.flatnonzero(comp == c)[0]]}/"
        f"{len({f for i in np.flatnonzero(comp == c) for f in mem[i]})}" for c in top))

    # components-from-HGNC-only (descriptive: no cytoband edges)
    held_nb = held.copy()
    held_nb["cytoband"] = ""
    comp_hgnc = L.components(held_nb)
    s2 = np.bincount(comp_hgnc)
    P(f"  (descriptive) HGNC-only components: {len(s2)}; largest {int(s2.max())} "
      f"({s2.max() / n:.2%}); singletons {int((s2 == 1).sum())}")

    out["splits"] = dict(
        n_genes=n,
        s281=dict(**lk281, n_blocks_taken=len(take281),
                  n_blocks_in_test=len(set(block[te281]))),
        s1=dict(**lkS1, n_blocks_taken=len(takeS1), n_families_held=len(fams_S1)),
        s2=dict(**lkS2, n_components=ncomp),
        random=lkR,
        components=dict(n=ncomp, largest=int(sizes.max()),
                        largest_share=float(sizes.max() / n),
                        n_singleton=int((sizes == 1).sum()),
                        n_2_4=int(((sizes >= 2) & (sizes <= 4)).sum()),
                        n_5_19=int(((sizes >= 5) & (sizes <= 19)).sum()),
                        n_20_99=int(((sizes >= 20) & (sizes <= 99)).sum()),
                        n_ge100=int((sizes >= 100).sum()),
                        median_size=float(np.median(sizes)),
                        hgnc_only_n=int(len(s2)), hgnc_only_largest=int(s2.max()),
                        hgnc_only_largest_share=float(s2.max() / n)),
    )
    # overlap of the three test sets
    out["splits"]["overlap_s1_281"] = dict(
        s1_superset_of_281=bool(np.all(te281 <= teS1)),
        n_only_in_s1=int((teS1 & ~te281).sum()), n_only_in_281=int((te281 & ~teS1).sum()),
        n_both=int((teS1 & te281).sum()))
    P(f"  S1 vs 281 test sets: both {int((teS1 & te281).sum())}, only S1 "
      f"{int((teS1 & ~te281).sum())}, only 281 {int((te281 & ~teS1).sum())}")

    np.savez_compressed(L.WORK / "splits284.npz", te281=te281, teS1=teS1, teS2=teS2,
                        teRND=teRND, comp=comp, block=np.asarray(block, dtype=object),
                        take281=np.asarray(take281, dtype=object),
                        takeS1=np.asarray(takeS1, dtype=object),
                        gene=held["gene"].to_numpy().astype(str))
    out["log"] = LOG
    (L.WORK / "ka284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "ka284.log").write_text("\n".join(LOG))
    P("wrote work/ka284.json, work/splits284.npz")


if __name__ == "__main__":
    main()

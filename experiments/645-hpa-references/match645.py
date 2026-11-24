"""645: HPA's bulk cell-line atlas as non-CCLE references for day-of identity (642's second stage).

Three reads, all through 590/642's matcher (CCLE cache + axis, 8 methods, Borda, stage-1 weights frozen):
 1. CALIBRATION: HPA's own profiles of lines that are ALSO in CCLE, scored CCLE-only. If bulk HPA is a CCLE-compatible
    platform, they are called as themselves at high scores (that licenses an HPA row as a reference beside CCLE).
 2. HEK293 from HPA: one HEK293_REF row built from HPA's bulk "HEK293" profile, scored on the 426 queries (590's + 597's
    real-Flex pilot). ALL THREE single-cell HEK293 profiles (Orion HEK293T FiCS, XuCao 3', Satija Flex) are now held out.
    Compared with 642's shipped non-Flex row (XuCao + Orion).
 3. VETO PANEL: every HPA non-cancer / immortalized line (no CCLE entry) added at once as NONCCLE:<line> rows. A call to
    one of them (8/8, >= 0.75) would only VETO a same-line switch-on (no bank). Counted: vetoes of queries whose true line
    IS in CCLE (false vetoes, must be 0) and calls on absent-line queries (RPE1, HEK293, iPSC, ...).
    uv run python experiments/645-hpa-references/match645.py   -> work/match645.json, work/match645.md
"""
from __future__ import annotations

import json
import re
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "experiments/590-identity-validation"))
sys.path.insert(0, str(ROOT / "experiments/642-hek293-contingency"))
import def_identity as DI  # noqa: E402
import match642 as MM  # noqa: E402

HPA = ROOT / "data/external/hpa/hpa_cache645.npz"
FLOOR2 = 0.75
NONCCLE = ["HEK293", "HEK TE", "hTERT-RPE1", "HaCaT", "hTERT-HME1", "RPTEC/TERT1", "HBEC3-KT", "TIME", "HUVEC/TERT2",
           "BJ [Human fibroblast]", "BJ1-hTERT", "BJ hTERT+ SV40 Large T+", "fHDF/TERT166", "PODO/TERT256", "PODO/SVTERT152",
           "hTEC/SVTERT24-B", "HBF/TERT88"]


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30


def norm(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", s.upper())


def main():
    t0 = time.time()
    Xfull, ccle_genes, samples, clinical = DI.load_cache()
    official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    g2c = {g: i for i, g in enumerate(ccle_genes)}
    common = [g for g in official if g in g2c]
    X = Xfull[:, np.array([g2c[g] for g in common], dtype=np.int64)]
    genes = np.asarray(common)
    del Xfull
    med = np.median(X, axis=0)
    variance, mean = np.var(X, axis=0), np.mean(X, axis=0)
    non_hk = np.array([not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
    fsets = {"all": np.arange(len(genes))}
    for k in (2000, 5000, 10000):
        el = np.flatnonzero(non_hk & (mean > 0.25))
        fsets[f"hvg{k}"] = el[np.argsort(-variance[el])[:k]]
    z = np.load(HPA, allow_pickle=False)
    H, hlines, hgenes = z["X"].astype(np.float64), z["lines"].astype(str), z["genes"].astype(str)
    hpos = {g: i for i, g in enumerate(hgenes)}
    hcol = np.array([hpos.get(g, -1) for g in genes])

    def hpa_row(line: str, fill=True) -> np.ndarray:
        v = H[list(hlines).index(line)]
        r = np.where(hcol >= 0, v[np.maximum(hcol, 0)], np.nan)
        return np.where(np.isfinite(r), r, med) if fill else r

    out = {"hpa_genes_on_axis": int((hcol >= 0).sum()), "n_axis": int(len(genes))}
    # ---------------------------------------------------------------- 1. calibration: HPA profiles of CCLE lines
    ccle_by = {}
    for s in samples:
        ccle_by.setdefault(norm(str(s).split("_")[0]), []).append(str(s))
    cal_q, cal_true = {}, {}
    for l in hlines:
        hits = ccle_by.get(norm(l), [])
        if len(hits) == 1:
            cal_q[f"hpa:{l}"] = hpa_row(l, fill=False)
            cal_true[f"hpa:{l}"] = hits[0]
    base_q, info = MM.build_queries(genes, official)
    _, qual0 = MM.score_all(X, samples, {k: base_q[k] for k in list(DI.ANCHOR_NAMES)}, fsets)   # weights from the anchors
    rows_cal, _ = MM.score_all(X, samples, {**cal_q, **{k: base_q[k] for k in DI.ANCHOR_NAMES}}, fsets, quality_fixed=qual0)
    cal = []
    for q, true in cal_true.items():
        r = rows_cal[q]
        cal.append(dict(query=q, true=true, call=r["call"], correct=r["call"] == true, unanimous=r["unanimous"],
                        n_agree=r["n_agree"], min_s1=r["min_s1"]))
    acc = float(np.mean([c["correct"] for c in cal]))
    s_ok = [c["min_s1"] for c in cal if c["correct"]]
    out["calibration"] = dict(n=len(cal), top1_acc=acc, unanimous_frac=float(np.mean([c["unanimous"] for c in cal])),
                              min_s1_correct_q=[float(np.quantile(s_ok, q)) for q in (0.01, 0.05, 0.5)] if s_ok else None,
                              wrong=[c for c in cal if not c["correct"]][:40], rows=cal)
    print(f"[645] calibration: {len(cal)} HPA profiles of CCLE lines: top-1 {acc:.3f}, unanimous {out['calibration']['unanimous_frac']:.3f}, "
          f"min_s1 q01/q05/q50 {out['calibration']['min_s1_correct_q']}  rss {rss():.2f} GiB  {time.time() - t0:.0f}s", flush=True)
    # ---------------------------------------------------------------- 2. HEK293 from HPA vs 642's non-Flex row
    queries = base_q
    _, qual = MM.score_all(X, samples, {k: queries[k] for k in DI.ANCHOR_NAMES}, fsets)
    assert qual == qual0
    nf = np.load(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", allow_pickle=False)
    nfpos = {g: i for i, g in enumerate(nf["genes"].astype(str))}
    nfrow = np.array([nf["values"][nfpos[g]] if g in nfpos and np.isfinite(nf["values"][nfpos[g]]) else med[j]
                      for j, g in enumerate(genes)], dtype=np.float64)
    variants = {"hpa_hek293": hpa_row("HEK293"), "nonflex642": nfrow}
    hekq = ["donor:orion_hek293t", "xu:HEK293", "pilot:HEK293"]
    res2 = {}
    for name, row in variants.items():
        Xa = np.vstack([X, row[None, :].astype(X.dtype)])
        sa = np.append(samples.astype(str), "HEK293_REF")
        rows, _ = MM.score_all(Xa, sa, queries, fsets, quality_fixed=qual0)
        acc_q = {q: (rows[q]["call"] == "HEK293_REF" and rows[q]["unanimous"] and rows[q]["min_s1"] >= FLOOR2) for q in queries}
        pos = {q: dict(call=rows[q]["call"], n_agree=rows[q]["n_agree"], min_s1=round(rows[q]["min_s1"], 4), accepted=acc_q[q]) for q in hekq}
        false = [q for q in queries if q not in hekq and acc_q[q]]
        near = sorted(((q, sum(v == "HEK293_REF" for v in rows[q]["top1s"].values()), round(rows[q]["min_s1"], 3))
                       for q in queries if q not in hekq and any(v == "HEK293_REF" for v in rows[q]["top1s"].values())),
                      key=lambda x: (-x[1], -x[2]))
        res2[name] = dict(positives=pos, false_accepted=false, near_misses=near[:15])
        print(f"[645] {name}: positives {pos}; false accepted {len(false)}; top near-misses {near[:5]}  rss {rss():.2f} GiB", flush=True)
    out["hek293"] = res2
    # ---------------------------------------------------------------- 3. veto panel (all non-CCLE HPA lines at once)
    present = [l for l in NONCCLE if l in set(hlines)]
    assert not any(ccle_by.get(norm(l)) for l in present), [l for l in present if ccle_by.get(norm(l))]
    rowsV = np.vstack([hpa_row(l) for l in present]).astype(X.dtype)
    Xa = np.vstack([X, rowsV])
    sa = np.append(samples.astype(str), [f"NONCCLE:{l}" for l in present])
    rows, _ = MM.score_all(Xa, sa, queries, fsets, quality_fixed=qual0)
    sample_set = set(samples.astype(str))
    veto = {}
    for q in queries:
        r = rows[q]
        if r["call"].startswith("NONCCLE:") and r["unanimous"] and r["min_s1"] >= FLOOR2:
            veto[q] = dict(call=r["call"], min_s1=round(r["min_s1"], 4), true=info[q].get("true"),
                           true_in_ccle=info[q].get("true") in sample_set)
    false_veto = {q: v for q, v in veto.items() if v["true_in_ccle"]}
    absent_calls = {q: v for q, v in veto.items() if not v["true_in_ccle"]}
    # anything (non-unanimous) where a NONCCLE row wins the consensus: the soft view
    soft = sorted(((q, r["call"], r["n_agree"], round(r["min_s1"], 3)) for q, r in rows.items() if r["call"].startswith("NONCCLE:")),
                  key=lambda x: -x[2])
    out["veto_panel"] = dict(lines=present, n_vetoes=len(veto), false_vetoes=false_veto, absent_line_calls=absent_calls,
                             consensus_wins=soft[:40])
    print(f"[645] veto panel ({len(present)} lines): accepted calls {len(veto)}; FALSE vetoes (true line in CCLE) {len(false_veto)}; "
          f"absent-line calls {absent_calls}; consensus wins {soft[:10]}  rss {rss():.2f} GiB", flush=True)
    (HERE / "work/match645.json").write_text(json.dumps(out, indent=1, default=str))
    print(f"[645] done {time.time() - t0:.0f}s, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    main()

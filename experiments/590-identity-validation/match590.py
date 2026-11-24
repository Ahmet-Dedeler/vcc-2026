"""590 — score every labelled profile with the day-of identity matcher, exactly as Oct 22 will run it.

The arithmetic is imported from `scripts/def_identity.py` (041's code, unchanged): the same CCLE cache, the same
official ∩ CCLE gene axis, the same 4 feature sets × {pearson, spearman}, the same three calibration anchors (so the
method weights are the Oct 22 weights), the same anchor-weighted Borda consensus and the same unanimity flag that
`def_sameline.py` gates on. Only the query set is larger.

Known answer (must pass before anything is read): the A/B/C rows reproduce a fresh `def_identity.match()` on
`data/challenge` exactly (method_quality, consensus rows, call, unanimity).

Optional `--flexbias <npz>` adds a per-gene log2 bias (Flex − 3') to every non-Flex profile before scoring (the
preregistration's Flex-like read); rows used to fit that bias are flagged by the bias file.

Writes `work/match590[<tag>].json` and `.csv`. Usage: `uv run python match590.py [--flexbias work/flexbias590.npz]`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
import def_identity as DI  # noqa: E402
import def_sameline as DS  # noqa: E402

REGISTRY = set(DS.REGISTRY)
# b_ctrl-based queries beyond the three anchors (path, true CCLE sample or None, chemistry)
BCTRL_QUERIES = {
    "donor:k562_essential": (ROOT / "data/donors/k562_essential", "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "10x 3' v3"),
    "donor:viperturb_k562": (ROOT / "data/donors/viperturb_k562", "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "VIPerturb-seq"),
    "donor:orion_hek293t": (ROOT / "data/donors/orion_hek293t", None, "Orion FiCS"),
    "donor:rpe1": (ROOT / "data/donors/rpe1", None, "10x 3' v3"),
    "donor:nadig_hepg2": (ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2", "HEPG2_LIVER", "10x 3'"),
    "donor:h1_2025": (ROOT / "data/donors/h1_2025", None, "10x Flex (Arc)"),
    "donor:h1_2025_full": (ROOT / "data/donors/h1_2025_full", None, "10x Flex (Arc)"),
    "donor:kolf2": (ROOT / "data/external/kolf2", None, "KOLF2.1J atlas"),
}
ANCHOR_TRUE = {"K562": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HCT116": "HCT116_LARGE_INTESTINE",
               "NADIG_JURKAT": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"}
ABC_TRUE = dict(DI.CALL_041)


def bctrl_query(donor_dir: Path, genes: np.ndarray) -> np.ndarray:
    """def_identity's anchor construction, verbatim: log2 CPM of expm1(b_ctrl) on the donor's genes, NaN elsewhere."""
    g = np.asarray((donor_dir / "genes.txt").read_text().splitlines())
    pos = {x: i for i, x in enumerate(g)}
    q = DI.log2_cpm_from_bctrl(np.load(donor_dir / "b_ctrl.npy"))
    overlap = np.array([gname in pos for gname in genes])
    full = np.full(len(genes), np.nan)
    full[overlap] = q[[pos[gname] for gname in genes[overlap]]]
    return full


def counts_query(counts_official: np.ndarray, official: list[str], genes: np.ndarray) -> np.ndarray:
    """A D/E/F-style query from official-axis counts: log2 CPM over the PRESENT official genes, mapped to `genes`."""
    present = np.isfinite(counts_official)
    tot = counts_official[present].sum()
    l2 = np.full(len(official), np.nan)
    l2[present] = np.log2(1.0 + 1e6 * counts_official[present] / tot)
    pos = {g: i for i, g in enumerate(official)}
    return l2[[pos[g] for g in genes]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flexbias", type=Path, default=None)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    t0 = time.time()

    Xfull, ccle_genes, samples, clinical = DI.load_cache()
    official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    gene_to_ccle = {g: i for i, g in enumerate(ccle_genes)}
    common = [g for g in official if g in gene_to_ccle]
    X = Xfull[:, np.array([gene_to_ccle[g] for g in common], dtype=np.int64)]
    genes = np.asarray(common)
    del Xfull
    sample_set = set(samples.tolist())
    lineage = clinical["ONCOTREE_LINEAGE"].to_dict()

    queries: dict[str, np.ndarray] = {}
    info: dict[str, dict] = {}
    # 1) A/B/C exactly as def_identity.match reads a context
    for ctx in ("A", "B", "C"):
        h = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad", backed="r")
        pos = {g: i for i, g in enumerate(h.var_names.astype(str))}
        q = DI.log2_cpm_from_counts(h.X)
        queries[ctx] = q[[pos[g] for g in genes]]
        info[ctx] = {"dataset": "challenge", "chemistry": "10x Flex (Arc)", "true": ABC_TRUE[ctx], "n_cells": int(h.shape[0])}
        h.file.close()
    # 2) the three anchors, verbatim
    for name, donor in DI.QUERY_DONORS.items():
        queries[name] = bctrl_query(ROOT / f"data/donors/{donor}", genes)
        info[name] = {"dataset": "anchor", "chemistry": "anchor", "true": ANCHOR_TRUE[name], "n_cells": None}
    # 3) other b_ctrl profiles
    for name, (d, true, chem) in BCTRL_QUERIES.items():
        if (d / "b_ctrl.npy").exists():
            queries[name] = bctrl_query(d, genes)
            info[name] = {"dataset": "donor", "chemistry": chem, "true": true, "n_cells": None}
    # 4) single-cell control pseudobulks built by profiles590.py
    for f in sorted((HERE / "work/profiles").glob("*.npz")):
        z = np.load(f)
        meta = json.loads(str(z["meta"]))
        for r, c in zip(meta, z["counts"]):
            queries[r["query"]] = counts_query(c, official, genes)
            info[r["query"]] = {k: r[k] for k in ("dataset", "chemistry", "true", "n_cells")}
    # Flex-like read: add the per-gene bias to non-Flex profiles (never to A/B/C, anchors or Flex rows)
    bias_note = None
    if a.flexbias:
        zb = np.load(a.flexbias, allow_pickle=False)
        bpos = {g: i for i, g in enumerate(zb["genes"].astype(str))}
        b = np.array([zb["bias"][bpos[g]] if g in bpos else 0.0 for g in genes])
        b = np.where(np.isfinite(b), b, 0.0)
        fitted_on = set(zb["fitted_on"].astype(str).tolist())
        for k in queries:
            if info[k]["dataset"] in ("challenge", "anchor") or "Flex" in str(info[k]["chemistry"]):
                continue
            # log2(1+cpm) + bias ~ log2(cpm * 2^bias); re-normalise to a CPM total over present genes
            q = queries[k]
            m = np.isfinite(q)
            cpm = np.expm1(q[m] * np.log(2)) * np.exp2(b[m])
            cpm = 1e6 * cpm / cpm.sum()
            q2 = q.copy()
            q2[m] = np.log2(1.0 + cpm)
            queries[k] = q2
            info[k]["flexbias_fitted_on_this_row"] = k in fitted_on
        bias_note = {"file": str(a.flexbias), "sha256": hashlib.sha256(a.flexbias.read_bytes()).hexdigest(),
                     "fitted_on": sorted(fitted_on)}
    print(f"[590] {len(queries)} queries; matrix {X.shape}; {time.time() - t0:.0f}s", flush=True)

    # feature sets, verbatim
    variance = np.var(X, axis=0)
    mean = np.mean(X, axis=0)
    import re
    non_hk = np.array([not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
    feature_sets = {"all": np.arange(len(genes))}
    for k in (2000, 5000, 10000):
        eligible = np.flatnonzero(non_hk & (mean > 0.25))
        feature_sets[f"hvg{k}"] = eligible[np.argsort(-variance[eligible])[:k]]

    results: dict[str, dict] = {}
    extra: dict[str, dict] = {}
    method_quality: dict[str, float] = {}
    for fname, fidx in feature_sets.items():
        for transform in ("pearson", "spearman"):
            key = f"{fname}:{transform}"
            Xm = X[:, fidx]
            if transform == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            # group queries by their valid-gene mask so only one Xq copy is alive at a time (memory)
            groups: dict[bytes, list[str]] = {}
            for qname, qfull in queries.items():
                valid = np.isfinite(qfull[fidx])
                groups.setdefault(hashlib.sha1(np.packbits(valid).tobytes()).digest(), []).append(qname)
            result, ex, rr = {}, {}, []
            for sig, members in groups.items():
                valid = np.isfinite(queries[members[0]][fidx])
                Xq = Xm[:, valid]
                for qname in members:
                    qfull = queries[qname]
                    q = qfull[fidx][valid]
                    if transform == "spearman":
                        q = rankdata(q)
                    scores = DI.row_pearson(Xq, q)
                    result[qname] = DI.top_rows(scores, samples, clinical)
                    order = np.argsort(-scores)
                    true = info[qname]["true"]
                    e = {"top1": str(samples[order[0]]), "gap12": float(scores[order[0]] - scores[order[1]]),
                         "s1": float(scores[order[0]]), "n_valid": int(valid.sum()),
                         # how far the best match stands above the bulk of CCLE (robust z of s1 vs all 1,600 lines)
                         "z1": float((scores[order[0]] - np.median(scores)) / (1.4826 * np.median(np.abs(scores - np.median(scores))) + 1e-12))}
                    if true is not None and true in sample_set:
                        where = np.flatnonzero(samples[order] == true)
                        e["true_rank"] = int(where[0] + 1)
                        e["true_score"] = float(scores[np.flatnonzero(samples == true)[0]])
                    ex[qname] = e
                del Xq
            for qname in DI.ANCHOR_NAMES:  # the anchors' reciprocal ranks, verbatim (iteration order of ANCHOR_NAMES)
                if qname in queries:
                    target = DI.ANCHOR_NAMES[qname]
                    where = np.flatnonzero(np.array([r["sample"] for r in result[qname]]) == target)
                    # top_rows keeps only 25 rows: recompute the full rank exactly as def_identity does
                    valid = np.isfinite(queries[qname][fidx])
                    q = queries[qname][fidx][valid]
                    if transform == "spearman":
                        q = rankdata(q)
                    sc = DI.row_pearson(Xm[:, valid], q)
                    w2 = np.flatnonzero(samples[np.argsort(-sc)] == target)
                    rank = int(w2[0] + 1) if w2.size else len(samples) + 1
                    rr.append(1.0 / rank)
                    result[qname + "_target_rank"] = rank
            quality = float(np.mean(rr)) if rr else 0.0
            method_quality[key] = quality
            results[key] = result
            extra[key] = ex
            print(f"{key} anchor ranks {[result.get(n + '_target_rank') for n in DI.ANCHOR_NAMES]} quality {quality:.4f}"
                  f"  ({len(groups)} gene masks)", flush=True)

    # consensus, verbatim (Borda with weights quality^2 over each method's top 25; top 40 kept)
    rows = []
    consensus = {}
    for qname in queries:
        acc: dict[str, float] = {}
        total_w = 0.0
        for key, result in results.items():
            w = method_quality[key] ** 2
            total_w += w
            for row in result[qname]:
                acc[row["sample"]] = acc.get(row["sample"], 0.0) + w * (26 - row["rank"])
        ordered = sorted(acc, key=acc.get, reverse=True)[:40]
        cscores = [acc[s] / max(total_w, 1e-12) for s in ordered]
        consensus[qname] = [{"rank": i + 1, "sample": s, "score": float(v)} for i, (s, v) in enumerate(zip(ordered, cscores))]
        call = ordered[0]
        top1s = {k: results[k][qname][0]["sample"] for k in results}
        unanimous = len(set(top1s.values())) == 1
        true = info[qname]["true"]
        true_in = true is not None and true in sample_set
        true_crank = next((i + 1 for i, s in enumerate(ordered) if s == true), None) if true_in else None
        gaps = [extra[k][qname]["gap12"] for k in results]
        rows.append({
            "query": qname, **{k: info[qname].get(k) for k in ("dataset", "chemistry", "n_cells")},
            "true": true, "true_in_ccle": true_in, "true_lineage": lineage.get(true) if true_in else None,
            "call": call, "call_lineage": lineage.get(call), "unanimous": unanimous,
            "n_methods_agree_call": sum(v == call for v in top1s.values()),
            "correct": bool(true_in and call == true),
            "call_in_registry": call in REGISTRY,
            "switch_on": bool(unanimous and call in REGISTRY),
            "c1": cscores[0], "c2": cscores[1] if len(cscores) > 1 else 0.0,
            "margin": cscores[0] - (cscores[1] if len(cscores) > 1 else 0.0),
            "true_consensus_rank": true_crank,
            "min_gap12": float(np.min(gaps)), "median_gap12": float(np.median(gaps)),
            "min_s1": float(np.min([extra[k][qname]["s1"] for k in results])),
            "median_s1": float(np.median([extra[k][qname]["s1"] for k in results])),
            "min_z1": float(np.min([extra[k][qname]["z1"] for k in results])),
            "n_valid_all": extra["all:pearson"][qname]["n_valid"],
            "flexbias_fitted_on_this_row": info[qname].get("flexbias_fitted_on_this_row"),
        })
    out = {"method_quality": method_quality, "n_queries": len(queries), "n_common_genes": int(len(genes)),
           "flexbias": bias_note, "rows": rows, "consensus": consensus,
           "methods_top25": {k: {q: v for q, v in r.items() if q in ("A", "B", "C")} for k, r in results.items()},
           "per_method": extra}

    # known answer: A/B/C vs a fresh def_identity.match on data/challenge
    ref = DI.match(ROOT / "data/challenge", ["A", "B", "C"])
    def top_set(rows_, n):
        return sorted((r["sample"], round(float(r["score"]), 12)) for r in rows_[:n])
    ka = {"method_quality": ref["method_quality"] == method_quality,
          "consensus": all(top_set(ref["consensus"][c], 25) == top_set(consensus[c], len(ref["consensus"][c][:25]))
                           for c in "ABC"),
          "call": {c: ref["call"][c] == consensus[c][0]["sample"] for c in "ABC"},
          "unanimous": {c: ref["call_unanimous_across_methods"][c] == next(r["unanimous"] for r in rows if r["query"] == c)
                        for c in "ABC"}}
    ka["PASS"] = bool(ka["method_quality"] and ka["consensus"] and all(ka["call"].values()) and all(ka["unanimous"].values()))
    out["known_answer"] = ka
    print("[590] KNOWN ANSWER", json.dumps(ka), flush=True)
    tag = a.tag or ("_flexbias" if a.flexbias else "")
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / f"work/match590{tag}.json").write_text(json.dumps(out, indent=1) + "\n")
    pd.DataFrame(rows).to_csv(HERE / f"work/match590{tag}.csv", index=False)
    print(f"[590] wrote work/match590{tag}.json/.csv in {time.time() - t0:.0f}s", flush=True)
    if not ka["PASS"]:
        raise SystemExit("known answer FAILED: do not read the rows")


if __name__ == "__main__":
    main()

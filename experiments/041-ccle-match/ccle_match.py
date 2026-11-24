"""Match VCC control pseudobulks to CCLE while validating the cross-platform metric.

The challenge controls are UMI pseudobulks; CCLE is bulk RNA-seq TPM. Raw
cross-gene correlation can be dominated by platform-specific gene effects, so
we do not trust a matching rule unless it recovers the two lines present in
both collections: K562 and HCT116.

Outputs rankings for several feature sets and similarity transforms. The
challenge-context consensus is weighted only by each method's anchor recovery;
the A/B/C identities never participate in method selection.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[2]
CCLE = ROOT / "data/external/ccle_broad_2025/data_mrna_seq_tpm.txt"
CLIN = ROOT / "data/external/ccle_broad_2025/data_clinical_sample.txt"
OUT = ROOT / "experiments/041-ccle-match"
ANCHOR_NAMES = {
    "K562": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
    "HCT116": "HCT116_LARGE_INTESTINE",
    "NADIG_JURKAT": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
}
QUERY_DONORS = {
    "K562": "k562_gwps",
    "HCT116": "orion_hct116",
    "NADIG_JURKAT": "nadig_jurkat",
}


def log2_cpm_from_counts(x) -> np.ndarray:
    if hasattr(x, "to_memory"):
        x = x.to_memory()
    total = np.asarray(x.sum(axis=0)).ravel().astype(np.float64)
    return np.log2(1.0 + 1e6 * total / total.sum())


def log2_cpm_from_bctrl(b: np.ndarray) -> np.ndarray:
    counts = np.expm1(np.asarray(b, dtype=np.float64))
    return np.log2(1.0 + 1e6 * counts / counts.sum())


def row_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean()
    return (a @ b) / (
        np.linalg.norm(a, axis=1) * max(float(np.linalg.norm(b)), 1e-12) + 1e-12
    )


def top_rows(scores: np.ndarray, samples: np.ndarray, clinical: pd.DataFrame, n=25):
    order = np.argsort(-scores)[:n]
    meta = clinical.reindex(samples[order])
    rows = []
    for rank, (i, sid) in enumerate(zip(order, samples[order]), 1):
        r = meta.loc[sid] if sid in meta.index else pd.Series(dtype=object)
        rows.append({
            "rank": rank,
            "sample": str(sid),
            "score": float(scores[i]),
            "lineage": None if pd.isna(r.get("ONCOTREE_LINEAGE")) else str(r.get("ONCOTREE_LINEAGE")),
            "disease": None if pd.isna(r.get("ONCOTREE_PRIMARY_DISEASE")) else str(r.get("ONCOTREE_PRIMARY_DISEASE")),
            "site": None if pd.isna(r.get("SAMPLE_COLLECTION_SITE")) else str(r.get("SAMPLE_COLLECTION_SITE")),
        })
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("loading CCLE expression", flush=True)
    ccle = pd.read_csv(CCLE, sep="\t")
    genes = ccle.pop("Hugo_Symbol").astype(str).to_numpy()
    ccle.pop("Entrez_Gene_Id")
    samples = ccle.columns.astype(str).to_numpy()
    X = ccle.to_numpy(dtype=np.float32).T
    del ccle
    clinical = pd.read_csv(CLIN, sep="\t", comment="#").set_index("SAMPLE_ID")

    vcc_genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
    gene_to_ccle = {g: i for i, g in enumerate(genes)}
    common = [g for g in vcc_genes if g in gene_to_ccle]
    cidx = np.array([gene_to_ccle[g] for g in common], dtype=np.int64)
    X = X[:, cidx]
    genes = np.asarray(common)
    print(f"matrix {X.shape}; common genes {len(genes)}", flush=True)

    queries: dict[str, np.ndarray] = {}
    for ctx in "ABC":
        a = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad", backed="r")
        pos = {g: i for i, g in enumerate(a.var_names.astype(str))}
        q = log2_cpm_from_counts(a.X)
        queries[ctx] = q[[pos[g] for g in genes]]
    for name, donor in QUERY_DONORS.items():
        donor_dir = ROOT / f"data/donors/{donor}"
        if not (donor_dir / "delta.npy").exists():
            continue
        g = np.asarray((donor_dir / "genes.txt").read_text().splitlines())
        pos = {x: i for i, x in enumerate(g)}
        b = np.load(donor_dir / "b_ctrl.npy")
        q = log2_cpm_from_bctrl(b)
        overlap = np.array([gname in pos for gname in genes])
        full = np.full(len(genes), np.nan)
        full[overlap] = q[[pos[gname] for gname in genes[overlap]]]
        queries[name] = full

    # Exclude genes whose expression is nearly constant across CCLE. The
    # feature sets are selected without looking at any query.
    variance = np.var(X, axis=0)
    mean = np.mean(X, axis=0)
    non_housekeeping = np.array([
        not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes
    ])
    feature_sets = {"all": np.arange(len(genes))}
    for k in (2000, 5000, 10000):
        eligible = np.flatnonzero(non_housekeeping & (mean > 0.25))
        feature_sets[f"hvg{k}"] = eligible[np.argsort(-variance[eligible])[:k]]

    results: dict[str, dict] = {}
    method_quality: dict[str, float] = {}
    for fname, fidx in feature_sets.items():
        for transform in ("pearson", "spearman"):
            key = f"{fname}:{transform}"
            Xm = X[:, fidx]
            if transform == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            result = {}
            reciprocal_ranks = []
            for qname, qfull in queries.items():
                valid = np.isfinite(qfull[fidx])
                Xq = Xm[:, valid]
                q = qfull[fidx][valid]
                if transform == "spearman":
                    q = rankdata(q)
                scores = row_pearson(Xq, q)
                rows = top_rows(scores, samples, clinical)
                result[qname] = rows
                if qname in ANCHOR_NAMES:
                    target = ANCHOR_NAMES[qname]
                    where = np.flatnonzero(samples[np.argsort(-scores)] == target)
                    rank = int(where[0] + 1) if where.size else len(samples) + 1
                    reciprocal_ranks.append(1.0 / rank)
                    result[qname + "_target_rank"] = rank
                    result[qname + "_target_score"] = float(
                        scores[np.flatnonzero(samples == target)[0]]
                    )
            quality = float(np.mean(reciprocal_ranks))
            method_quality[key] = quality
            result["anchor_mean_reciprocal_rank"] = quality
            results[key] = result
            print(
                key,
                "anchor ranks",
                [result.get(name + "_target_rank") for name in ANCHOR_NAMES
                 if name in queries],
                "quality",
                f"{quality:.4f}",
                flush=True,
            )

    # Weighted Borda consensus over methods. Weak anchor-recovery methods get
    # negligible weight; no A/B/C result changes the method weights.
    consensus = {}
    for qname in "ABC":
        acc = {s: 0.0 for s in samples}
        total_w = 0.0
        for key, result in results.items():
            w = method_quality[key] ** 2
            total_w += w
            for row in result[qname]:
                acc[row["sample"]] += w * (26 - row["rank"])
        ordered = sorted(acc, key=acc.get, reverse=True)[:40]
        scores = np.array([acc[s] / max(total_w, 1e-12) for s in ordered])
        consensus[qname] = top_rows(
            scores, np.asarray(ordered), clinical, n=min(25, len(ordered))
        )

    out = {
        "source": "cBioPortal CCLE Broad 2025 mirror (DepMap Public 25Q3)",
        "n_samples": int(len(samples)),
        "n_common_genes": int(len(genes)),
        "method_quality": method_quality,
        "methods": results,
        "consensus": consensus,
    }
    (OUT / "results.json").write_text(json.dumps(out, indent=2))
    for qname in "ABC":
        print(f"\n{qname} consensus")
        for row in consensus[qname][:10]:
            print(
                f"  {row['rank']:2d} {row['sample']:<45} "
                f"{row['lineage']!s:<20} {row['score']:.3f}"
            )


if __name__ == "__main__":
    main()

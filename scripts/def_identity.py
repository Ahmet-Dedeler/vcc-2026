"""Identify the round's cell lines from their control profiles alone, in minutes.

Pre-staged, cache-backed form of `experiments/041-ccle-match/ccle_match.py`, written for
**Oct 22** (experiment 199, transition item 2b). 041 identified the validation contexts as
**A = JURKAT, B = HELA, C = CAL33** with the same top match under all eight feature/rank
variants, after calibrating the cross-platform metric on three known controls. That call is
what motivated the CD4 residual now shipping in a15, so line identity is not trivia: it
selects which auxiliary donors are matched.

The matching arithmetic is 041's, unchanged and copied verbatim, so `--verify` reproduces its
published call exactly. The only thing that changed is the ingest: 041 spent most of its wall
clock in `pd.read_csv` on a 313 MB TSV. That parse is done once, ahead of the round, into
`ccle_cache.npz`; the matcher then loads a float32 array and runs in seconds.

The cache is keyed on CCLE's OWN gene axis, not the challenge's, so it survives a D/E/F
`gene_names.csv` that differs from A/B/C's. The bundle's axis is read from the file at match
time and never assumed.

    # one-off, already done (writes data/external/ccle_broad_2025/ccle_cache.npz)
    uv run python scripts/def_identity.py --build-cache

    # the known-answer row: must print A=JURKAT B=HELA C=CAL33 and PASS
    uv run python scripts/def_identity.py --verify

    # the round
    uv run python scripts/def_identity.py --bundle data/challenge --contexts D,E,F \
        --out experiments/199-def-runbook/work/def-identity.json

**Boundary (041, AGENTS.md §7).** Public data may train or condition a model. It may not be
used to hand-revise an output. A line identity selects donors and adapters; it never licenses
a literature lookup of a measured effect.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os  # 642: the HEK293 contingency (DEF_DAY_ALLOW_NONCCLE_REF)
import re
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
CCLE_TPM = ROOT / "data/external/ccle_broad_2025/data_mrna_seq_tpm.txt"
CCLE_CLIN = ROOT / "data/external/ccle_broad_2025/data_clinical_sample.txt"
CACHE = ROOT / "data/external/ccle_broad_2025/ccle_cache.npz"
REF_041 = ROOT / "experiments/041-ccle-match/results.json"

# The cross-platform calibration. These three are present in BOTH collections, so their rank
# is a measurement of the metric, not of the challenge. No A/B/C (or D/E/F) result is allowed
# to influence which methods are trusted.
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
# 041's published call, used as this script's known-answer row.
CALL_041 = {"A": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
            "B": "HELA_CERVIX",
            "C": "CAL33_UPPER_AERODIGESTIVE_TRACT"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- 041, verbatim

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


def top_rows(scores, samples, clinical, n=25):
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


# --------------------------------------------------------------------------------- the cache

def build_cache() -> dict:
    """Parse the 313 MB CCLE TSV once, ahead of the round, into a float32 .npz.

    Stores CCLE's own gene axis with duplicates in file order, so the `{gene: index}` map the
    matcher rebuilds is the same dict 041 built (last occurrence wins) for ANY bundle axis.
    """
    t0 = time.time()
    for path in (CCLE_TPM, CCLE_CLIN):
        if not path.exists():
            raise SystemExit(
                f"missing {path}. The mirror is cBioPortal's CCLE Broad 2025 study "
                "(DepMap Public 25Q3): https://cbioportal-datahub.s3.amazonaws.com/"
                "ccle_broad_2025.tar.gz . Pin it BEFORE Oct 22 — do not depend on a "
                "download working on the day."
            )
    print(f"[cache] parsing {CCLE_TPM.name} ({CCLE_TPM.stat().st_size / 1e6:.0f} MB)", flush=True)
    ccle = pd.read_csv(CCLE_TPM, sep="\t")
    genes = ccle.pop("Hugo_Symbol").astype(str).to_numpy()
    ccle.pop("Entrez_Gene_Id")
    samples = ccle.columns.astype(str).to_numpy()
    X = ccle.to_numpy(dtype=np.float32).T
    del ccle
    clinical = pd.read_csv(CCLE_CLIN, sep="\t", comment="#").set_index("SAMPLE_ID")
    keep = ["ONCOTREE_LINEAGE", "ONCOTREE_PRIMARY_DISEASE", "SAMPLE_COLLECTION_SITE"]
    clin = clinical.reindex(columns=keep).astype(object).where(pd.notna(clinical.reindex(columns=keep)), "")
    np.savez(
        CACHE,
        X=X,
        genes=genes.astype("U"),
        samples=samples.astype("U"),
        clin_index=np.asarray(clinical.index.astype(str), dtype="U"),
        clin_values=clin.to_numpy().astype("U"),
        clin_columns=np.asarray(keep, dtype="U"),
    )
    meta = {
        "source": "cBioPortal CCLE Broad 2025 mirror (DepMap Public 25Q3)",
        "built_by": "scripts/def_identity.py --build-cache",
        "tpm_txt": {"path": str(CCLE_TPM.relative_to(ROOT)),
                    "bytes": CCLE_TPM.stat().st_size, "sha256": sha256(CCLE_TPM)},
        "clinical_txt": {"path": str(CCLE_CLIN.relative_to(ROOT)),
                         "bytes": CCLE_CLIN.stat().st_size, "sha256": sha256(CCLE_CLIN)},
        "cache": {"path": str(CACHE.relative_to(ROOT)), "bytes": CACHE.stat().st_size,
                  "sha256": sha256(CACHE)},
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "n_samples": int(X.shape[0]),
        "n_ccle_genes": int(X.shape[1]),
        "build_seconds": round(time.time() - t0, 1),
    }
    CACHE.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[cache] {X.shape} -> {CACHE} in {meta['build_seconds']}s "
          f"({CACHE.stat().st_size / 1e6:.0f} MB)", flush=True)
    return meta


def load_cache():
    if not CACHE.exists():
        raise SystemExit(f"no cache at {CACHE}; run: uv run python {Path(__file__).name} --build-cache")
    z = np.load(CACHE, allow_pickle=False)
    clinical = pd.DataFrame(z["clin_values"], index=pd.Index(z["clin_index"]),
                            columns=list(z["clin_columns"]))
    clinical = clinical.replace("", np.nan)
    return z["X"], z["genes"], z["samples"], clinical


# ------------------------------------------------------------------------------- the matcher

# 642: the HEK293 contingency. A line ABSENT from CCLE (Arc's 2026-10-01 post: the six lines span "disease status"; CCLE
# has no HEK293, RPE1, MCF10A, BJ or HaCaT) can only be named by a reference row outside CCLE. The second stage below adds
# ONE pinned row and re-scores the contexts with stage 1's method weights frozen; it never changes a stage-1 key. Off
# unless DEF_DAY_ALLOW_NONCCLE_REF names the reference. Validated in experiments/642-hek293-contingency (real-Flex HEK293
# 8/8 at 0.820; 0 false calls on 423 labelled profiles).
NONCCLE_REFS = {
    "hek293": ("data/identity_refs/hek293_nonflex_ref642.npz", "1b829de027a3fa8ec0c7615b2388f75fef6dc9da46520d2c06bd272500d60704"),
}
NONCCLE_FLOOR = 0.75   # the second stage's own floor (stage 1 keeps 590's 0.65): the true positive is 0.820, the
                       # nearest non-HEK293 profiles with the row top-1 in >= 5/8 methods sit at 0.46-0.69


def load_nonccle_ref(key: str):
    rel, want = NONCCLE_REFS[key]
    got = sha256(ROOT / rel)
    if got != want:
        raise SystemExit(f"642: {rel} sha256 {got[:16]}... is not the pinned {want[:16]}...: refusing an unpinned reference")
    z = np.load(ROOT / rel, allow_pickle=False)
    prov = json.loads(str(z["provenance"]))
    flex = [s.get("name") for s in prov.get("sources", []) if "flex" in str(s.get("chemistry", "")).lower()]
    if flex or not prov.get("sources"):
        # 642: a Flex-built HEK293 row pulled Satija's Flex HAP1 to a unanimous 8/8 call at 0.818 (same-platform references
        # attract same-platform queries; the contexts are Flex). Only non-Flex sources may build a reference.
        raise SystemExit(f"642: reference {key} lists Flex (or no) sources {flex}: refusing (never build one from Flex)")
    return z["genes"].astype(str), z["values"].astype(np.float64), prov


def nonccle_stage(key: str, X, genes, samples, feature_sets: dict, method_quality: dict, queries: dict, names) -> dict:
    """Score `names` (keys of `queries`) against CCLE + the pinned `key` row; method weights = stage 1's (frozen). A gene
    the reference lacks takes CCLE's per-gene median. `call[name]` is the line only when the augmented consensus is the
    row, unanimous over the methods and min_top1_score >= NONCCLE_FLOOR."""
    rgenes, rvals, prov = load_nonccle_ref(key)
    pos = {g: i for i, g in enumerate(rgenes)}
    med = np.median(X, axis=0)
    row = np.array([rvals[pos[g]] if g in pos and np.isfinite(rvals[pos[g]]) else med[j] for j, g in enumerate(genes)])
    rname = f"{key.upper()}_REF"
    Xa = np.vstack([X, row[None, :].astype(X.dtype)])
    sa = np.append(np.asarray(samples).astype(str), rname)
    top = {n: {} for n in names}
    for fname, fidx in feature_sets.items():
        for transform in ("pearson", "spearman"):
            k = f"{fname}:{transform}"
            Xm = Xa[:, fidx]
            if transform == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            for n in names:
                qf = queries[n]
                valid = np.isfinite(qf[fidx])
                q = qf[fidx][valid]
                if transform == "spearman":
                    q = rankdata(q)
                sc = row_pearson(Xm[:, valid], q)
                order = np.argsort(-sc)[:25]
                top[n][k] = [(str(sa[i]), r + 1, float(sc[i])) for r, i in enumerate(order)]
            del Xm
    per = {}
    for n in names:
        acc, tw = {}, 0.0
        for k, rows in top[n].items():
            w = method_quality[k] ** 2
            tw += w
            for s, rank, _ in rows:
                acc[s] = acc.get(s, 0.0) + w * (26 - rank)
        ordered = sorted(acc, key=acc.get, reverse=True)[:5]
        top1 = {k: rows[0][0] for k, rows in top[n].items()}
        s1 = float(min(rows[0][2] for rows in top[n].values()))
        unanimous = len(set(top1.values())) == 1
        per[n] = {"call": ordered[0], "unanimous": unanimous, "n_methods_agree": sum(v == ordered[0] for v in top1.values()),
                  "min_top1_score": s1, "top5": [(s, acc[s] / max(tw, 1e-12)) for s in ordered],
                  "row_top1_methods": sum(v == rname for v in top1.values()),
                  "accepted": bool(ordered[0] == rname and unanimous and s1 >= NONCCLE_FLOOR)}
    return {"reference": key, "row_name": rname, "file": NONCCLE_REFS[key][0], "sha256": NONCCLE_REFS[key][1],
            "provenance": prov, "floor": NONCCLE_FLOOR, "weights": "stage 1 method_quality, frozen", "per_context": per,
            "call": {n: (key.upper() if per[n]["accepted"] else None) for n in names}}


def match(bundle: Path, contexts: list[str], genes_csv: Path | None = None) -> dict:
    t0 = time.time()
    Xfull, ccle_genes, samples, clinical = load_cache()
    t_cache = time.time() - t0

    # The gene axis is READ FROM THE FILE, never assumed: a D/E/F bundle is allowed to ship a
    # different axis and this must not silently mis-align. Both challenge CSVs carry a header
    # row even though `vcc prep --genes` documents "headerless".
    genes_csv = genes_csv or (bundle / "gene_names.csv")
    vcc_genes = pd.read_csv(genes_csv)["gene_name"].astype(str)
    gene_to_ccle = {g: i for i, g in enumerate(ccle_genes)}
    common = [g for g in vcc_genes if g in gene_to_ccle]
    cidx = np.array([gene_to_ccle[g] for g in common], dtype=np.int64)
    X = Xfull[:, cidx]
    genes = np.asarray(common)
    del Xfull
    print(f"[match] cache {time.time() - t0:.1f}s; matrix {X.shape}; "
          f"common genes {len(genes)}/{len(vcc_genes)}", flush=True)

    t1 = time.time()
    queries: dict[str, np.ndarray] = {}
    for ctx in contexts:
        path = bundle / f"context_{ctx}.h5ad"
        a = ad.read_h5ad(path, backed="r")
        pos = {g: i for i, g in enumerate(a.var_names.astype(str))}
        missing = [g for g in genes if g not in pos]
        if missing:
            raise SystemExit(f"{path.name}: {len(missing)} matched genes absent from its var axis")
        q = log2_cpm_from_counts(a.X)
        queries[ctx] = q[[pos[g] for g in genes]]
        print(f"[match] {path.name}: {a.shape[0]} control cells", flush=True)
    for name, donor in QUERY_DONORS.items():
        donor_dir = ROOT / f"data/donors/{donor}"
        if not (donor_dir / "delta.npy").exists():
            print(f"[match] WARNING calibration anchor {name} missing ({donor_dir})", flush=True)
            continue
        g = np.asarray((donor_dir / "genes.txt").read_text().splitlines())
        pos = {x: i for i, x in enumerate(g)}
        b = np.load(donor_dir / "b_ctrl.npy")
        q = log2_cpm_from_bctrl(b)
        overlap = np.array([gname in pos for gname in genes])
        full = np.full(len(genes), np.nan)
        full[overlap] = q[[pos[gname] for gname in genes[overlap]]]
        queries[name] = full
    t_queries = time.time() - t1

    # Feature sets are selected from CCLE alone, without looking at any query.
    t2 = time.time()
    variance = np.var(X, axis=0)
    mean = np.mean(X, axis=0)
    non_housekeeping = np.array([not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
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
                result[qname] = top_rows(scores, samples, clinical)
                if qname in ANCHOR_NAMES:
                    target = ANCHOR_NAMES[qname]
                    where = np.flatnonzero(samples[np.argsort(-scores)] == target)
                    rank = int(where[0] + 1) if where.size else len(samples) + 1
                    reciprocal_ranks.append(1.0 / rank)
                    result[qname + "_target_rank"] = rank
                    result[qname + "_target_score"] = float(scores[np.flatnonzero(samples == target)[0]])
            quality = float(np.mean(reciprocal_ranks)) if reciprocal_ranks else 0.0
            method_quality[key] = quality
            result["anchor_mean_reciprocal_rank"] = quality
            results[key] = result
            print(key, "anchor ranks",
                  [result.get(n + "_target_rank") for n in ANCHOR_NAMES if n in queries],
                  "quality", f"{quality:.4f}", flush=True)

    # Weighted Borda consensus. Weak anchor-recovery methods get negligible weight; no
    # challenge-context result changes a method weight.
    consensus = {}
    for qname in contexts:
        acc = {s: 0.0 for s in samples}
        total_w = 0.0
        for key, result in results.items():
            w = method_quality[key] ** 2
            total_w += w
            for row in result[qname]:
                acc[row["sample"]] += w * (26 - row["rank"])
        ordered = sorted(acc, key=acc.get, reverse=True)[:40]
        scores = np.array([acc[s] / max(total_w, 1e-12) for s in ordered])
        consensus[qname] = top_rows(scores, np.asarray(ordered), clinical, n=min(25, len(ordered)))
    t_score = time.time() - t2

    out = {
        "source": "cBioPortal CCLE Broad 2025 mirror (DepMap Public 25Q3)",
        "bundle": str(bundle),
        "contexts": contexts,
        "gene_axis_file": str(genes_csv),
        "n_samples": int(len(samples)),
        "n_common_genes": int(len(genes)),
        "method_quality": method_quality,
        "methods": results,
        "consensus": consensus,
        "timing_seconds": {"cache_load": round(t_cache, 2), "queries": round(t_queries, 2),
                           "score": round(t_score, 2), "total": round(time.time() - t0, 2)},
        "top1_per_method": {
            ctx: {k: results[k][ctx][0]["sample"] for k in results} for ctx in contexts
        },
        "call": {ctx: consensus[ctx][0]["sample"] for ctx in contexts},
        "call_unanimous_across_methods": {
            ctx: len({results[k][ctx][0]["sample"] for k in results}) == 1 for ctx in contexts
        },
        # 590: the lowest top-1 correlation over the 8 methods. A line absent from CCLE is still called
        # unanimously (HEK293 -> HAP1, 8/8) at 0.51-0.61; correct calls sit >= 0.68 (q01 over 390 labelled
        # profiles), A/B/C at 0.719-0.772. def_sameline.py --min-top1 gates the same-line mode on it.
        "min_top1_score": {ctx: min(results[k][ctx][0]["score"] for k in results) for ctx in contexts},
    }
    for ctx in contexts:
        print(f"\n{ctx} consensus")
        for row in consensus[ctx][:10]:
            print(f"  {row['rank']:2d} {row['sample']:<45} {str(row['lineage']):<20} {row['score']:.3f}")
        print(f"  lowest top-1 correlation over the 8 methods: {out['min_top1_score'][ctx]:.3f} "
              f"(590: same-line rows need >= 0.65; A/B/C 0.719-0.772)")
    # 642: the HEK293 contingency's second stage (off unless DEF_DAY_ALLOW_NONCCLE_REF names a reference). Every key
    # above is final; this only adds `nonccle_ref`.
    _nref = os.environ.get("DEF_DAY_ALLOW_NONCCLE_REF", "")
    if _nref:
        if _nref not in NONCCLE_REFS:
            raise SystemExit(f"642: DEF_DAY_ALLOW_NONCCLE_REF={_nref!r} is not one of {sorted(NONCCLE_REFS)}")
        out["nonccle_ref"] = nonccle_stage(_nref, X, genes, samples, feature_sets, method_quality, queries, list(contexts))
        for ctx in contexts:
            r = out["nonccle_ref"]["per_context"][ctx]
            print(f"[642] {ctx}: with the {_nref} reference row: {r['call']} {r['n_methods_agree']}/8, min top-1 "
                  f"{r['min_top1_score']:.3f}, row top-1 in {r['row_top1_methods']}/8 methods -> "
                  f"{'CALLED ' + _nref.upper() if r['accepted'] else 'not called'}", flush=True)
    return out


def verify(out: dict) -> bool:
    """Known-answer row: reproduce 041's published call and its anchor calibration exactly."""
    ref = json.loads(REF_041.read_text())
    ok = True
    print("\n--- known-answer row: experiments/041-ccle-match/results.json ---")
    for field in ("n_samples", "n_common_genes"):
        good = out[field] == ref[field]
        ok &= good
        print(f"  {field:<16} {out[field]} vs {ref[field]}   {'OK' if good else 'MISMATCH'}")
    for key in sorted(ref["method_quality"]):
        for anchor in ANCHOR_NAMES:
            a, b = out["methods"][key].get(anchor + "_target_rank"), ref["methods"][key].get(anchor + "_target_rank")
            good = a == b
            ok &= good
            if not good:
                print(f"  MISMATCH anchor rank {key} {anchor}: {a} vs {b}")
    print(f"  anchor ranks     8 methods x 3 anchors   "
          f"{'all reproduce' if ok else 'MISMATCH above'}")
    for ctx, expect in CALL_041.items():
        got = out["call"][ctx]
        ref_top = ref["consensus"][ctx][0]["sample"]
        good = got == expect == ref_top
        ok &= good
        print(f"  consensus {ctx}      {got:<45} expected {expect}   {'OK' if good else 'MISMATCH'}")
        unan = out["call_unanimous_across_methods"][ctx]
        print(f"    top-1 identical under all 8 feature/rank variants: {unan}")
        ok &= bool(unan)
    print(f"\nVERIFY {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--build-cache", action="store_true", help="parse the CCLE TSV into ccle_cache.npz and exit")
    p.add_argument("--verify", action="store_true", help="run on A/B/C and reproduce 041's call")
    p.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    p.add_argument("--contexts", default="", help="comma-separated, e.g. D,E,F; default = the bundle manifest's")
    p.add_argument("--genes", type=Path, default=None, help="gene_names.csv; default <bundle>/gene_names.csv")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    if args.build_cache:
        build_cache()
        return

    bundle = args.bundle if args.bundle.is_absolute() else ROOT / args.bundle
    if args.verify:
        contexts = ["A", "B", "C"]
    elif args.contexts:
        contexts = [c.strip() for c in args.contexts.split(",") if c.strip()]
    else:
        contexts = list(json.loads((bundle / "manifest.json").read_text())["contexts"])
    out = match(bundle, contexts, args.genes)
    out["verified_against_041"] = verify(out) if args.verify else None
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.out}")
    if args.verify and not out["verified_against_041"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

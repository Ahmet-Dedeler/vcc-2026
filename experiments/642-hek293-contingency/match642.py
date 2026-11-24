"""642 step 2: does a HEK293 reference row, added to the CCLE matrix, recognise HEK293 without creating a false call?

Same matcher as 590 / def_identity (CCLE cache, official ∩ CCLE axis, 4 feature sets from CCLE's own variance x
{pearson, spearman}, anchor-weighted Borda over each method's top 25, unanimity = all 8 top-1 equal, min_s1 = 590's
min_top1). Only the candidate set grows: one extra row "HEK293_REF:<variant>" in CCLE units (log2(1 + CPM) over the
profile's present genes; a gene absent from the profile takes CCLE's per-gene median). Feature sets stay CCLE-only.

Queries: every 590 query (A/B/C, the three anchors, the donor b_ctrl profiles, 590's 28 single-cell profile files incl.
200 MIX-seq lines) + the three real-Flex pilot profiles (pilot642.py): 597's K562 / HAP1 / HEK293.

Reference variants (all HEK293 / 293T non-targeting controls on disk):
  xu      XuCao2023 HEK293, 10x 3' (590's profile xu.npz)
  orion   X-Atlas/Orion HEK293T b_ctrl, Orion FiCS
  flex    Satija VIPerturb pilot HEK293, 10x Flex v1 (486 NT cells)
  nonflex mean of xu and orion in log2 space (the day-of candidate: no Flex in the reference, Flex only in the query)
Positives are CROSS-SOURCE only (a reference is never scored against the profile it was built from).

Known answer: with no extra row the rows reproduce 590's match590.json (call, unanimous, min_s1) for every common query.
    uv run python experiments/642-hek293-contingency/match642.py      -> work/match642.json, work/match642.md
"""
from __future__ import annotations

import hashlib
import json
import re
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "experiments/590-identity-validation"))
import def_identity as DI  # noqa: E402
import match590 as M  # noqa: E402

REF_PREFIX = "HEK293_REF"
FLOOR = 0.65


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30


def build_queries(genes, official):
    queries, info = {}, {}
    for ctx in ("A", "B", "C"):
        h = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad", backed="r")
        pos = {g: i for i, g in enumerate(h.var_names.astype(str))}
        queries[ctx] = DI.log2_cpm_from_counts(h.X)[[pos[g] for g in genes]]
        info[ctx] = {"dataset": "challenge", "chemistry": "10x Flex (Arc)", "true": M.ABC_TRUE[ctx]}
        h.file.close()
    for name, donor in DI.QUERY_DONORS.items():
        queries[name] = M.bctrl_query(ROOT / f"data/donors/{donor}", genes)
        info[name] = {"dataset": "anchor", "chemistry": "anchor", "true": M.ANCHOR_TRUE[name]}
    for name, (d, true, chem) in M.BCTRL_QUERIES.items():
        if (d / "b_ctrl.npy").exists():
            queries[name] = M.bctrl_query(d, genes)
            info[name] = {"dataset": "donor", "chemistry": chem, "true": true}
    for f in sorted((ROOT / "experiments/590-identity-validation/work/profiles").glob("*.npz")):
        z = np.load(f)
        meta = json.loads(str(z["meta"]))
        for r, c in zip(meta, z["counts"]):
            queries[r["query"]] = M.counts_query(c, official, genes)
            info[r["query"]] = {k: r[k] for k in ("dataset", "chemistry", "true")}
    z = np.load(HERE / "work/flexpilot642.npz")
    for k, c in zip(z["lines"].astype(str), z["counts"]):
        queries[f"pilot:{k}"] = M.counts_query(c, official, genes)
        info[f"pilot:{k}"] = {"dataset": "satija_pilot", "chemistry": "10x Flex v1",
                              "true": {"K562": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HAP1": "HAP1_ENGINEERED"}.get(k)}
    return queries, info


def hek_source(qname: str, info: dict) -> str | None:
    """Which HEK293 source a query comes from (None: not HEK293)."""
    if qname == "donor:orion_hek293t":
        return "orion"
    if qname == "pilot:HEK293":
        return "flex"
    if info.get("dataset", "").lower().startswith("xu") or qname.lower().startswith("xu"):
        return "xu"
    return None


def score_all(X, samples, queries, feature_sets, quality_fixed=None):
    """quality_fixed: use these method weights (590's, from the CCLE-only run) instead of re-deriving them from the
    anchors' ranks over the augmented candidate set, so adding a row cannot re-weight the methods."""
    results, s1s, quality = {}, {}, {}
    for fname, fidx in feature_sets.items():
        for transform in ("pearson", "spearman"):
            key = f"{fname}:{transform}"
            Xm = X[:, fidx]
            if transform == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            groups = {}
            for qname, qfull in queries.items():
                valid = np.isfinite(qfull[fidx])
                groups.setdefault(hashlib.sha1(np.packbits(valid).tobytes()).digest(), []).append(qname)
            res, s1 = {}, {}
            for _, members in groups.items():
                valid = np.isfinite(queries[members[0]][fidx])
                Xq = Xm[:, valid]
                for qname in members:
                    q = queries[qname][fidx][valid]
                    if transform == "spearman":
                        q = rankdata(q)
                    sc = DI.row_pearson(Xq, q)
                    order = np.argsort(-sc)[:25]
                    res[qname] = [(str(samples[i]), r + 1) for r, i in enumerate(order)]
                    s1[qname] = float(sc[order[0]])
                del Xq
            rr = []
            for an, target in (DI.ANCHOR_NAMES.items() if quality_fixed is None else ()):
                valid = np.isfinite(queries[an][fidx])
                q = queries[an][fidx][valid]
                if transform == "spearman":
                    q = rankdata(q)
                sc = DI.row_pearson(Xm[:, valid], q)
                w = np.flatnonzero(samples[np.argsort(-sc)] == target)
                rr.append(1.0 / (int(w[0] + 1) if w.size else len(samples) + 1))
            quality[key] = float(np.mean(rr)) if quality_fixed is None else quality_fixed[key]
            results[key], s1s[key] = res, s1
            del Xm
    rows = {}
    for qname in queries:
        acc, tw = {}, 0.0
        for key, res in results.items():
            w = quality[key] ** 2
            tw += w
            for s, rank in res[qname]:
                acc[s] = acc.get(s, 0.0) + w * (26 - rank)
        ordered = sorted(acc, key=acc.get, reverse=True)[:5]
        top1s = {k: results[k][qname][0][0] for k in results}
        call = ordered[0]
        rows[qname] = dict(call=call, unanimous=len(set(top1s.values())) == 1,
                           n_agree=sum(v == call for v in top1s.values()),
                           min_s1=float(min(s1s[k][qname] for k in results)), top5=ordered, top1s=top1s,
                           c1=acc[ordered[0]] / tw, c2=(acc[ordered[1]] / tw) if len(ordered) > 1 else 0.0)
    return rows, quality


def main():
    t0 = time.time()
    Xfull, ccle_genes, samples, clinical = DI.load_cache()
    official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    g2c = {g: i for i, g in enumerate(ccle_genes)}
    common = [g for g in official if g in g2c]
    X = Xfull[:, np.array([g2c[g] for g in common], dtype=np.int64)]
    genes = np.asarray(common)
    del Xfull
    queries, info = build_queries(genes, official)
    print(f"[642] {len(queries)} queries, CCLE {X.shape}, rss {rss():.2f} GiB, {time.time() - t0:.0f}s", flush=True)
    variance, mean = np.var(X, axis=0), np.mean(X, axis=0)
    non_hk = np.array([not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
    fsets = {"all": np.arange(len(genes))}
    for k in (2000, 5000, 10000):
        el = np.flatnonzero(non_hk & (mean > 0.25))
        fsets[f"hvg{k}"] = el[np.argsort(-variance[el])[:k]]
    med = np.median(X, axis=0)

    def ref_row(v):
        r = np.where(np.isfinite(v), v, med)
        return r.astype(np.float32), int((~np.isfinite(v)).sum())

    xu_q = [q for q in queries if hek_source(q, info.get(q, {})) == "xu"]
    assert len(xu_q) == 1, xu_q
    src = {"xu": queries[xu_q[0]], "orion": queries["donor:orion_hek293t"], "flex": queries["pilot:HEK293"]}
    both = np.where(np.isfinite(src["xu"]) & np.isfinite(src["orion"]), (src["xu"] + src["orion"]) / 2,
                    np.where(np.isfinite(src["xu"]), src["xu"], src["orion"]))
    variants = {"xu": (src["xu"], {"xu"}), "orion": (src["orion"], {"orion"}), "flex": (src["flex"], {"flex"}),
                "nonflex": (both, {"xu", "orion"})}

    out = {"n_queries": len(queries), "queries_hek": {q: hek_source(q, info[q]) for q in queries if hek_source(q, info[q])}}
    base, qual0 = score_all(X, samples, queries, fsets)
    ref590 = {r["query"]: r for r in json.loads((ROOT / "experiments/590-identity-validation/work/match590.json").read_text())["rows"]}
    mism = [q for q in queries if q in ref590 and (ref590[q]["call"] != base[q]["call"] or ref590[q]["unanimous"] != base[q]["unanimous"]
                                                   or abs(ref590[q]["min_s1"] - base[q]["min_s1"]) > 1e-9)]
    out["ka"] = dict(n_common=sum(q in ref590 for q in queries), mismatches=mism, PASS=not mism)
    print(f"[642] KA vs 590: {out['ka']['n_common']} common queries, {len(mism)} mismatches -> {'PASS' if not mism else 'FAIL'}"
          f"  rss {rss():.2f} GiB", flush=True)
    if mism:
        raise SystemExit(f"known answer FAILED: {mism[:5]}")
    out["base"] = {q: base[q] for q in queries}
    lines = ["| variant | NaN-filled genes | method quality = 590's | positives (cross-source): call, n/8, min_s1 | "
             "false HEK293 calls | other calls changed | switch-eligible false calls |", "|---|---|---|---|---|---|---|"]
    for vname, (v, built_from) in variants.items():
        r, nnan = ref_row(v)
        Xa = np.vstack([X, r[None, :]])
        sa = np.append(samples, f"{REF_PREFIX}:{vname}")
        rows, qual = score_all(Xa, sa, queries, fsets, quality_fixed=qual0)
        hekname = f"{REF_PREFIX}:{vname}"
        pos, false, changed = {}, [], []
        for q in queries:
            s = hek_source(q, info[q])
            if s is not None and s not in built_from:
                pos[q] = (rows[q]["call"] == hekname, rows[q]["n_agree"], round(rows[q]["min_s1"], 4), rows[q]["unanimous"])
            elif s is None:
                if rows[q]["call"] == hekname:
                    false.append((q, rows[q]["n_agree"], round(rows[q]["min_s1"], 4), rows[q]["unanimous"]))
                elif rows[q]["call"] != base[q]["call"]:
                    changed.append((q, base[q]["call"], rows[q]["call"]))   # one-stage only; the two-stage rule keeps base
        elig = [f for f in false if f[3] and f[2] >= FLOOR]
        out[vname] = dict(nan_filled=nnan, method_quality_same=qual == qual0, quality=qual, positives=pos, false_calls=false,
                          changed=changed, switch_eligible_false=elig,
                          hek_rank_in_top5={q: (rows[q]["top5"].index(hekname) + 1 if hekname in rows[q]["top5"] else None)
                                            for q in queries},
                          rows={q: rows[q] for q in queries})
        lines.append(f"| {vname} | {nnan} | {qual == qual0} | "
                     + "; ".join(f"{q}: {'HEK' if p[0] else 'other'} {p[1]}/8 {p[2]}{' U' if p[3] else ''}" for q, p in pos.items())
                     + f" | {len(false)} | {len(changed)} | {len(elig)} |")
        print(lines[-1], f"rss {rss():.2f} GiB", flush=True)
        for f in false[:10]:
            print("    false:", f, flush=True)
        for c in changed[:10]:
            print("    changed:", c, flush=True)
    (HERE / "work/match642.json").write_text(json.dumps(out, indent=1, default=str))
    (HERE / "work/match642.md").write_text("\n".join(lines) + "\n")
    print(f"[642] done {time.time() - t0:.0f}s, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    main()

"""597 — real-Flex HEK293 / HAP1 / K562 non-targeting pseudobulks from Satija's VIPerturb pilot, scored with 590's matcher.

Preregistered in PREREGISTRATION.md (e2b4e80e) before the download. The RDS is read with 424's streaming parser
(`rds_stream.RDSStream`: large vectors spill to disk), never `readRDS`, so the process stays light.

  uv run python experiments/597-flex-hek293-floor/score597.py inspect <rds> <spill dir>
      structure only: assays, layers, meta.data columns, value counts of the small categorical columns -> work/inspect597.json
  uv run python experiments/597-flex-hek293-floor/score597.py score <rds> <spill dir> <line col> <guide col> <nt regex>
      per line (values of <line col> naming K562 / HAP1 / HEK293), the cells whose <guide col> matches <nt regex>; sums of the
      RNA assay's raw `counts`; 590's matcher unchanged (same CCLE cache, axis, 8 methods, 3 anchors, Borda); known answer =
      590's method_quality and A/B/C rows -> work/score597.json
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
sys.path.insert(0, str(ROOT / "experiments/424-viperturb-k562"))
sys.path.insert(0, str(ROOT / "experiments/590-identity-validation"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE))
from rds_stream import RList, S4, Spilled  # noqa: E402
from rds_stream597 import RDSStream597  # noqa: E402  (424 reader + R byte code)

OFFICIAL = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()


def log(m: str) -> None:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    print(f"{time.strftime('%H:%M:%S')} [rss {rss:.2f} GiB] {m}", flush=True)


def names_of(x):
    n = getattr(x, "attrs", {}).get("names")
    return list(n) if n is not None else None


def named(x, key):
    return x[names_of(x).index(key)]


def dimnames(x):
    dn = x.attrs.get("dimnames")
    return (list(dn[0]) if dn is not None and dn[0] is not None else None,
            list(dn[1]) if dn is not None and len(dn) > 1 and dn[1] is not None else None)


def column(df, col):
    v = named(df, col)
    lev = getattr(v, "attrs", {}).get("levels")
    if lev is not None:
        lev = list(lev)
        return np.asarray([lev[i - 1] if i > 0 else None for i in np.asarray(v)], dtype=object)
    return np.asarray(v, dtype=object) if isinstance(v, RList) else np.asarray(v)


def parse(rds: Path, spill: Path):
    spill.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    reader = RDSStream597(str(rds), str(spill), log=lambda m: None)
    hdr, obj, eof = reader.read()
    log(f"parsed {rds.name} in {time.time() - t0:.0f}s; header {hdr}; clean EOF {eof}; byte-code items skipped "
        f"{reader.n_bytecode}; class {list(obj.attrs.get('class', []))}")
    if not eof:
        raise SystemExit("the stream did not end where the object did: a misparse")
    return obj


def inspect(rds: Path, spill: Path) -> None:
    obj = parse(rds, spill)
    out = {"class": list(obj.attrs.get("class", [])), "assays": {}}
    assays = obj.attrs["assays"]
    for an in names_of(assays):
        a = named(assays, an)
        rec = {"class": list(a.attrs.get("class", []))}
        if "layers" in a.attrs:
            rec["layers"] = names_of(a.attrs["layers"])
            feats = dimnames(a.attrs["features"])[0] if "features" in a.attrs else None
            rec["n_features"] = len(feats) if feats else None
            rec["features_head"] = feats[:8] if feats else None
        out["assays"][an] = rec
    meta = obj.attrs["meta.data"]
    cols = names_of(meta)
    out["meta_columns"] = cols
    out["n_cells"] = len(list(meta.attrs["row.names"]))
    small = {}
    for c in cols:
        v = column(meta, c)
        u = pd.Series(v.astype(str)).value_counts()
        if len(u) <= 60:
            small[c] = {k: int(n) for k, n in u.items()}
        else:
            small[c] = {"n_unique": int(len(u)), "head": {k: int(n) for k, n in u.head(8).items()}}
    out["meta_values"] = small
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/inspect597.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    log(f"assays {list(out['assays'])}; RNA layers {out['assays'].get('RNA', {}).get('layers')}; "
        f"{out['n_cells']} cells; meta columns {cols}")
    for c, v in small.items():
        log(f"  {c}: {v if len(str(v)) < 400 else str(v)[:400]}")


LINES = {"K562": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "HAP1": "HAP1_ENGINEERED", "HEK293": None}


def score(rds: Path, spill: Path, line_col: str, guide_col: str, nt_regex: str) -> None:
    import def_identity as DI
    import match590 as M
    from scipy.stats import rankdata
    import hashlib
    obj = parse(rds, spill)
    rna = named(obj.attrs["assays"], "RNA")
    cnt = named(rna.attrs["layers"], "counts")
    assert list(cnt.attrs["class"]) == ["dgCMatrix"], cnt.attrs["class"]
    features = dimnames(rna.attrs["features"])[0]
    meta = obj.attrs["meta.data"]
    cells = list(meta.attrs["row.names"])
    dim = [int(v) for v in cnt.attrs["Dim"]]
    assert dim == [len(features), len(cells)], (dim, len(features), len(cells))
    p = np.asarray(cnt.attrs["p"], dtype=np.int64)
    ispill, xspill = cnt.attrs["i"], cnt.attrs["x"]
    line = column(meta, line_col).astype(str)
    guide = column(meta, guide_col).astype(str)
    nt = np.asarray([bool(re.search(nt_regex, g)) for g in guide])
    # line labels -> our three names (case-insensitive substring)
    lab = np.full(len(cells), "", dtype=object)
    for key in LINES:
        lab[[key.lower() in str(x).lower().replace("-", "") for x in line]] = key
    # HEK293 must not swallow other labels; HAP1/K562 are disjoint strings
    groups = [k for k in LINES if ((lab == k) & nt).sum() > 0]
    log(f"line col {line_col!r}, guide col {guide_col!r}, NT regex {nt_regex!r}: "
        + ", ".join(f"{k} {int(((lab == k) & nt).sum())} NT of {int((lab == k).sum())}" for k in LINES))
    sel = {k: np.flatnonzero((lab == k) & nt) for k in groups}
    code = np.full(len(cells), -1, np.int64)
    for gi, k in enumerate(groups):
        code[sel[k]] = gi
    G, F = len(groups), len(features)
    acc = np.zeros((G, F), dtype=np.float64)
    integer = True
    with open(ispill.path, "rb") as fi, open(xspill.path, "rb") as fx:
        for c0 in range(0, len(cells), 4000):
            c1 = min(c0 + 4000, len(cells))
            a, z = int(p[c0]), int(p[c1])
            fi.seek(4 * a); fx.seek(4 * a)
            ii = np.fromfile(fi, dtype=np.int32, count=z - a)
            xx = np.fromfile(fx, dtype=np.float32, count=z - a).astype(np.float64)
            integer &= bool(np.all(xx == np.round(xx)))
            local = np.repeat(np.arange(c1 - c0), np.diff(p[c0:c1 + 1]))
            g = code[c0:c1][local]
            m = g >= 0
            np.add.at(acc, (g[m], ii[m]), xx[m])
    log(f"pseudobulks built; integer counts: {integer}; totals " + ", ".join(f"{k} {acc[i].sum():.0f}" for i, k in enumerate(groups)))
    if not integer:
        raise SystemExit("the counts layer is not integer: stop (preregistration step 3)")
    # onto the official axis, NaN where the probe set lacks a gene
    pos = {g: j for j, g in enumerate(features)}
    official_counts = {}
    for i, k in enumerate(groups):
        v = np.full(len(OFFICIAL), np.nan)
        for j, gname in enumerate(OFFICIAL):
            q = pos.get(gname)
            if q is not None:
                v[j] = acc[i, q]
        official_counts[k] = v
    # ---- 590's matcher, unchanged arithmetic (M.counts_query, DI primitives, the same anchors) ----
    Xfull, ccle_genes, samples, clinical = DI.load_cache()
    g2c = {g: i for i, g in enumerate(ccle_genes)}
    common = [g for g in OFFICIAL if g in g2c]
    X = Xfull[:, np.array([g2c[g] for g in common], dtype=np.int64)]
    genes = np.asarray(common)
    del Xfull
    queries = {}
    import anndata as ad
    for ctx in "ABC":
        h = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad", backed="r")
        posc = {g: i for i, g in enumerate(h.var_names.astype(str))}
        queries[ctx] = DI.log2_cpm_from_counts(h.X)[[posc[g] for g in genes]]
        h.file.close()
    for name, donor in DI.QUERY_DONORS.items():
        queries[name] = M.bctrl_query(ROOT / f"data/donors/{donor}", genes)
    for k in groups:
        queries[f"pilot:{k}"] = M.counts_query(official_counts[k], OFFICIAL, genes)
    import re as _re
    variance, mean = np.var(X, axis=0), np.mean(X, axis=0)
    non_hk = np.array([not _re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
    fsets = {"all": np.arange(len(genes))}
    for kk in (2000, 5000, 10000):
        el = np.flatnonzero(non_hk & (mean > 0.25))
        fsets[f"hvg{kk}"] = el[np.argsort(-variance[el])[:kk]]
    results, quality = {}, {}
    for fname, fidx in fsets.items():
        for tr in ("pearson", "spearman"):
            key = f"{fname}:{tr}"
            Xm = X[:, fidx]
            if tr == "spearman":
                Xm = rankdata(Xm, axis=1).astype(np.float32)
            res, rr = {}, []
            for qn, qf in queries.items():
                valid = np.isfinite(qf[fidx])
                q = qf[fidx][valid]
                if tr == "spearman":
                    q = rankdata(q)
                sc = DI.row_pearson(Xm[:, valid], q)
                res[qn] = DI.top_rows(sc, samples, clinical)
                if qn in DI.ANCHOR_NAMES:
                    w = np.flatnonzero(samples[np.argsort(-sc)] == DI.ANCHOR_NAMES[qn])
                    rr.append(1.0 / (int(w[0] + 1) if w.size else len(samples) + 1))
            quality[key] = float(np.mean(rr))
            results[key] = res
    ref = json.loads((ROOT / "experiments/590-identity-validation/work/match590.json").read_text())
    ka = {"method_quality_equals_590": quality == ref["method_quality"]}
    out_rows = {}
    for qn in queries:
        acc2, tw = {}, 0.0
        for key, res in results.items():
            w = quality[key] ** 2
            tw += w
            for row in res[qn]:
                acc2[row["sample"]] = acc2.get(row["sample"], 0.0) + w * (26 - row["rank"])
        ordered = sorted(acc2, key=acc2.get, reverse=True)[:40]
        top1s = {k: results[k][qn][0]["sample"] for k in results}
        lin = clinical["ONCOTREE_LINEAGE"].to_dict()
        out_rows[qn] = {"call": ordered[0], "unanimous": len(set(top1s.values())) == 1,
                        "n_methods_agree": sum(v == ordered[0] for v in top1s.values()),
                        "min_top1_score": float(min(results[k][qn][0]["score"] for k in results)),
                        "top5": [(s, round(acc2[s] / tw, 3), lin.get(s)) for s in ordered[:5]]}
    want = {"A": ("JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", 0.772), "B": ("HELA_CERVIX", 0.719),
            "C": ("CAL33_UPPER_AERODIGESTIVE_TRACT", 0.733)}
    ka["abc_reproduce"] = all(out_rows[c]["call"] == want[c][0] and out_rows[c]["unanimous"]
                              and abs(out_rows[c]["min_top1_score"] - want[c][1]) < 5e-4 for c in "ABC")
    ka["PASS"] = ka["method_quality_equals_590"] and ka["abc_reproduce"]
    log(f"KNOWN ANSWER {ka}")
    h = hashlib.sha256()
    with open(rds, "rb") as fh:  # streamed: reading the 1.58 GB file whole put the first run over 4 GB
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    del X
    res = {"known_answer": ka, "rds_sha256": h.hexdigest(),
           "line_col": line_col, "guide_col": guide_col, "nt_regex": nt_regex,
           "cells": {k: int(len(sel[k])) for k in groups}, "rows": out_rows}
    (HERE / "work/score597.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    if not ka["PASS"]:
        raise SystemExit("known answer FAILED: do not read the pilot rows")
    for k in groups:
        r = out_rows[f"pilot:{k}"]
        log(f"{k}: {len(sel[k])} NT cells -> {r['call']} unanimous={r['unanimous']} ({r['n_methods_agree']}/8) "
            f"min_top1={r['min_top1_score']:.4f} top5={[s for s, _, _ in r['top5']]}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "inspect":
        inspect(Path(sys.argv[2]), Path(sys.argv[3]))
    elif cmd == "score":
        score(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])

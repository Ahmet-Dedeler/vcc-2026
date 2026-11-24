"""285 item 1 -- INVENTORY FROM METADATA ONLY.  No expression value is read here: obs / var / X layout only.

    PYTHONDONTWRITEBYTECODE=1 uv run python inv285.py

Writes work/inv285.json + work/inv285.log.
"""
from __future__ import annotations

import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/285-reliability-benchmark"
WORK = HERE / "work"
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    if v.dtype.kind in "iuf":
        return np.asarray([str(x) for x in v])
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def col(h, frame: str, key: str):
    """obs/var column, resolving the legacy `__categories` group used by the Nadig files."""
    g = h[frame]
    if key in g and not (isinstance(g[key], h5py.Dataset) and "__categories" in g and key in g["__categories"]):
        pass
    if "__categories" in g and key in g["__categories"]:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in g["__categories"][key][:]])
        return cats[g[key][:]]
    return labels_of(h, f"{frame}/{key}")


def index_of(h, frame: str):
    return col(h, frame, str(h[frame].attrs["_index"]))


def xlayout(h):
    X = h["X"]
    if isinstance(X, h5py.Group):
        return dict(kind=str(X.attrs.get("encoding-type", "?")), shape=[int(x) for x in X.attrs["shape"]],
                    data=str(X["data"].dtype), indices=str(X["indices"].dtype), nnz=int(X["data"].shape[0]),
                    compression=str(X["data"].compression), chunks=str(X["data"].chunks))
    return dict(kind="dense", shape=[int(x) for x in X.shape], data=str(X.dtype),
                compression=str(X.compression), chunks=str(X.chunks))


def panel() -> list[str]:
    import csv
    with open(ROOT / "data/challenge/pert_counts.csv") as fh:
        return [r["target_gene"] for r in csv.DictReader(fh)]


def summarise_counts(c) -> dict:
    c = np.asarray(list(c))
    if c.size == 0:
        return {}
    return dict(n=int(c.size), min=int(c.min()), q25=float(np.percentile(c, 25)), median=float(np.median(c)),
                q75=float(np.percentile(c, 75)), max=int(c.max()),
                ge50=int((c >= 50).sum()), ge100=int((c >= 100).sum()), ge200=int((c >= 200).sum()),
                ge400=int((c >= 400).sum()), ge800=int((c >= 800).sum()))


def main() -> None:
    t0 = time.time()
    PAN = set(panel())
    out: dict = {"panel_size": len(PAN), "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    P(f"285 inv285.py -- {out['utc']}; panel {len(PAN)} targets")

    # ---------------------------------------------------------------- h1-test
    p = ROOT / "data/harness/h1-test.h5ad"
    with h5py.File(p, "r") as h:
        lab = labels_of(h, "obs/target_gene")
        bc = index_of(h, "obs")
        lay, ok = xlayout(h), list(h["obs"].keys())
    suffix = np.asarray([b.split("-", 1)[1] for b in bc])
    c = Counter(lab.tolist())
    tg = sorted(set(lab.tolist()) - {"non-targeting"})
    out["h1_test"] = dict(path=str(p.relative_to(ROOT)), obs_keys=ok, X=lay, n_targets=len(tg),
                          lanes=len(set(suffix.tolist())), cells_per_target=summarise_counts([c[t] for t in tg]),
                          panel_overlap=len(set(tg) & PAN), controls=c["non-targeting"])
    P(f"h1-test: {lay['shape']} {lay['kind']}, {len(tg)} targets, {len(set(suffix.tolist()))} lanes, 400 cells each; "
      f"panel overlap {out['h1_test']['panel_overlap']}")

    # ---------------------------------------------------------------- Replogle K562 (two experiments)
    kset = {}
    for name, fn in (("k562_gwps", "ReplogleWeissman2022_K562_gwps.h5ad"),
                     ("k562_essential", "ReplogleWeissman2022_K562_essential.h5ad")):
        p = ROOT / "data/external/scperturb" / fn
        with h5py.File(p, "r") as h:
            ok, lay = list(h["obs"].keys()), xlayout(h)
            pert = col(h, "obs", "perturbation")
            gene = col(h, "obs", "gene")
            batch = col(h, "obs", "batch")
            guide = col(h, "obs", "guide_id")
            umi = np.asarray(h["obs/UMI_count"][:], dtype=np.float64) if "UMI_count" in h["obs"] else None
            gnames = index_of(h, "var")
        cp = Counter(pert.tolist())
        tg = sorted(set(pert.tolist()) - {"control"})
        pan_t = sorted(set(tg) & PAN)
        # guides per perturbation, and cells per (pert, guide)
        pg = defaultdict(Counter)
        for t, g in zip(pert, guide):
            pg[t][g] += 1
        multi = {t: len(pg[t]) for t in tg if len(pg[t]) > 1}
        kset[name] = set(tg)
        d = dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, obs_keys=ok, X=lay,
                 n_cells=int(pert.size), n_perts=len(tg), n_control=cp["control"],
                 n_genes=len(gnames), cells_per_pert=summarise_counts([cp[t] for t in tg]),
                 n_panel=len(pan_t), cells_per_panel_target=summarise_counts([cp[t] for t in pan_t]),
                 n_batches=len(set(batch.tolist())), guide_key="guide_id",
                 cells_per_batch=summarise_counts(Counter(batch.tolist()).values()),
                 n_control_guides=len(pg["control"]),
                 guides_per_pert_hist=dict(Counter(len(pg[t]) for t in tg)), n_multi_guide=len(multi),
                 multi_guide_panel=sorted(set(multi) & PAN),
                 multi_guide_panel_cells={t: dict(pg[t]) for t in sorted(set(multi) & PAN)},
                 umi_per_cell=dict(median=float(np.median(umi)), q25=float(np.percentile(umi, 25)),
                                   q75=float(np.percentile(umi, 75))) if umi is not None else None)
        # panel targets that reach the ladder rungs on a per-guide basis
        d["panel_cells_ge"] = {n: int(sum(cp[t] >= n for t in pan_t)) for n in (50, 100, 200, 400, 800, 1600)}
        out[name] = d
        P(f"{name}: {lay['shape']} {lay['kind']} {lay['data']} ({lay.get('compression')}), {d['n_perts']} perts "
          f"({d['n_panel']} panel), {d['n_control']} controls in {d['n_control_guides']} guides, {d['n_batches']} gemgroups; "
          f"panel cells/target median {d['cells_per_panel_target'].get('median')}, ≥400 {d['panel_cells_ge'][400]}, "
          f"≥800 {d['panel_cells_ge'][800]}; two-construct panel targets {len(d['multi_guide_panel'])}; "
          f"median UMI/cell {d['umi_per_cell']['median'] if d['umi_per_cell'] else '?'}")
    ga, gb = kset["k562_gwps"], kset["k562_essential"]
    out["k562_study_overlap"] = dict(gwps=len(ga), essential=len(gb), shared=len(ga & gb),
                                     shared_on_panel=sorted(ga & gb & PAN), essential_on_panel=sorted(gb & PAN),
                                     gwps_on_panel=len(ga & PAN))
    P(f"K562 GWPS vs essential: shared perturbations {len(ga & gb)}; shared ∩ panel {len(ga & gb & PAN)} "
      f"{sorted(ga & gb & PAN)[:10]}; essential ∩ panel {len(gb & PAN)}")

    # ---------------------------------------------------------------- Nadig 2025
    for name, fn in (("nadig_jurkat", "jurkat_raw_singlecell.h5ad"), ("nadig_hepg2", "hepg2_raw_singlecell.h5ad")):
        p = ROOT / "data/external/nadig2025" / fn
        with h5py.File(p, "r") as h:
            ok, lay = list(h["obs"].keys()), xlayout(h)
            gene = col(h, "obs", "gene")
            gg = col(h, "obs", "gem_group")
            sg = col(h, "obs", "sgID_AB")
            umi = np.asarray(h["obs/UMI_count"][:], dtype=np.float64) if "UMI_count" in h["obs"] else None
            gnames = col(h, "var", "gene_name")
        cg = Counter(gene.tolist())
        ctrl = [t for t in cg if t.lower().startswith("non-target") or t.lower() in ("control", "ntc", "safe-harbor", "negative_control")]
        real = sorted(set(cg) - set(ctrl))
        pan_t = sorted(set(real) & PAN)
        pg = defaultdict(Counter)
        for t, g in zip(gene, sg):
            pg[t][g] += 1
        d = dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, obs_keys=ok, X=lay, n_cells=int(gene.size),
                 control_labels=ctrl, n_control=int(sum(cg[c] for c in ctrl)), n_perts=len(real),
                 n_genes=len(gnames), cells_per_pert=summarise_counts([cg[t] for t in real]),
                 n_panel=len(pan_t), cells_per_panel_target=summarise_counts([cg[t] for t in pan_t]),
                 n_gem_groups=len(set(gg.tolist())), cells_per_gem_group=summarise_counts(Counter(gg.tolist()).values()),
                 guides_per_pert_hist=dict(Counter(len(pg[t]) for t in real)),
                 panel_multi_guide=sorted([t for t in pan_t if len(pg[t]) > 1]),
                 panel_guide_cells={t: dict(pg[t]) for t in pan_t},
                 panel_cells_ge={n: int(sum(cg[t] >= n for t in pan_t)) for n in (50, 100, 200, 400, 800)},
                 umi_per_cell=dict(median=float(np.median(umi)), q25=float(np.percentile(umi, 25)),
                                   q75=float(np.percentile(umi, 75))) if umi is not None else None)
        out[name] = d
        P(f"{name}: {lay['shape']} {lay['kind']} {lay['data']}, {d['n_perts']} perts ({d['n_panel']} panel), "
          f"controls {ctrl} n={d['n_control']}, {d['n_gem_groups']} gem groups, guides/pert {d['guides_per_pert_hist']}; "
          f"panel cells/target median {d['cells_per_panel_target'].get('median')} (≥200 {d['panel_cells_ge'][200]}, "
          f"≥400 {d['panel_cells_ge'][400]}); multi-guide panel targets {len(d['panel_multi_guide'])}; "
          f"median UMI {d['umi_per_cell']['median'] if d['umi_per_cell'] else '?'}")

    # ---------------------------------------------------------------- 268 Broad-16 halves
    hv = ROOT / "experiments/268-broad16-context-residuals/work/halves"
    lines = {}
    for f in sorted(hv.glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        meta = json.loads(f.with_suffix(".json").read_text()) if f.with_suffix(".json").exists() else {}
        lines[f.stem] = dict(keys=list(z.keys()), shapes={k: list(np.asarray(z[k]).shape) for k in z.keys()},
                             meta={k: v for k, v in meta.items() if len(str(v)) < 500})
    out["broad16_halves"] = lines
    k0 = list(lines)[0]
    P(f"268 halves: {len(lines)} lines {sorted(lines)}; keys {lines[k0]['keys']}; shapes {lines[k0]['shapes']}")

    # ---------------------------------------------------------------- 280 CD4 bank
    z = np.load(ROOT / "experiments/280-gwcd4i-donor-replicates/work/obs280.npz", allow_pickle=True)
    keys = list(z.keys())
    rows = np.load(ROOT / "experiments/280-gwcd4i-donor-replicates/work/bank280_rows.npy")
    b = dict(obs_keys=keys, bank_rows=int(rows.size),
             bank_shape=list(np.load(ROOT / "experiments/280-gwcd4i-donor-replicates/work/bank280_X.npy", mmap_mode="r").shape))
    for k in keys:
        v = np.asarray(z[k])
        if v.ndim == 1 and v.size > 1000:
            u = np.unique(v)
            b[f"{k}_n_unique"] = int(u.size)
            if u.size <= 12:
                b[f"{k}_uniques"] = [str(x) for x in u]
    nc = np.asarray(z["n_cells"])[rows]
    b["bank_n_cells"] = summarise_counts(nc)
    tot = np.asarray(z["total"])[rows].astype(np.float64)
    b["umi_per_cell_bank"] = dict(median=float(np.median(tot / np.maximum(nc, 1))))
    out["cd4_bank280"] = b
    P(f"280 bank: {b['bank_rows']} rows {b['bank_shape']}; obs keys {keys}; n_cells median "
      f"{b['bank_n_cells']['median']} (≥50 {b['bank_n_cells']['ge50']}, ≥200 {b['bank_n_cells']['ge200']}); "
      f"median UMI/cell {b['umi_per_cell_bank']['median']:.0f}")

    # ---------------------------------------------------------------- 267 cached half banks
    h267 = ROOT / "experiments/267-replicate-validated-oracle/work/half"
    c267 = {}
    for line in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        d = h267 / line
        if d.exists():
            c267[line] = {f.name: f.stat().st_size for f in sorted(d.iterdir()) if f.is_file()}
        t = h267 / f"truth-{line}.npz"
        if t.exists():
            zz = np.load(t, allow_pickle=True)
            c267.setdefault(line, {})["truth_keys"] = list(zz.keys())
    out["e267_cache"] = c267
    P("267 cache: " + json.dumps({k: [x for x in v if not x.endswith('.npy')][:12] for k, v in c267.items()}))

    out["seconds"] = round(time.time() - t0, 1)
    (WORK / "inv285.json").write_text(json.dumps(out, indent=1, default=str))
    (WORK / "inv285.log").write_text("\n".join(LOG) + "\n")
    P(f"\ndone in {out['seconds']}s")


if __name__ == "__main__":
    main()

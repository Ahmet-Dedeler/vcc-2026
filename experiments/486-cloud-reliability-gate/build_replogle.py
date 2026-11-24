"""Build per-line log2-CPM fold-change matrices from Replogle 2022 raw pseudobulks.

K562 essential and K562 genome-wide are the same cell line. They are aggregated to
one K562 profile (inverse-variance weighted on shared entries) so a leave-one-line-out
fold cannot see the line it is predicting.

Point estimate: cell-count-weighted mean of guide-level per-cell UMI profiles, then
log2((CPM + 1) / (NTC CPM + 1)).

Sampling variance of that log2 ratio: Poisson delta method on the per-cell means,
plus the across-guide variance of the log2 ratio divided by the number of guides.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np

LN2 = float(np.log(2.0))
MU_FLOOR = 1e-3


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def symbol_from_obs(name: str) -> str:
    """Replogle bulk obs names are '{id}_{symbol}_{promoter or transcript}_{ensembl}'."""
    bits = str(name).split("_")
    if len(bits) < 2:
        return str(name)
    return bits[1]


def collapse_duplicate_symbols(X: np.ndarray, names: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sum columns that share a gene symbol (split Ensembl ids, e.g. HSPA14)."""
    names = np.asarray(names).astype(str)
    uniq, inv = np.unique(names, return_inverse=True)
    if len(uniq) == len(names):
        return np.asarray(X, dtype=np.float64), uniq
    acc = np.zeros((len(uniq), X.shape[0]), dtype=np.float64)
    np.add.at(acc, inv, np.asarray(X, dtype=np.float64).T)
    return acc.T, uniq


def _cpm(mu: np.ndarray) -> np.ndarray:
    total = float(np.sum(mu))
    if total <= 0:
        return np.zeros_like(mu)
    return mu / total * 1e6


def poisson_var_log2(mu_p: np.ndarray, n_p: float, mu_c: np.ndarray, n_c: float) -> np.ndarray:
    mp = np.maximum(mu_p, MU_FLOOR)
    mc = np.maximum(mu_c, MU_FLOOR)
    return (1.0 / (LN2 * LN2)) * (1.0 / (n_p * mp) + 1.0 / (n_c * mc))


def screen_to_arrays(path: Path) -> dict:
    a = ad.read_h5ad(path)
    names = a.var["gene_name"].astype(str).to_numpy()
    X, genes = collapse_duplicate_symbols(np.asarray(a.X), names)
    symbols = np.array([symbol_from_obs(x) for x in a.obs_names], dtype=object)
    nc = a.obs["num_cells_filtered"].to_numpy(dtype=np.float64)
    nu = a.obs["num_cells_unfiltered"].to_numpy(dtype=np.float64)
    n_cells_row = np.where(np.isfinite(nc) & (nc > 0), nc, nu)
    n_cells_row = np.where(np.isfinite(n_cells_row) & (n_cells_row > 0), n_cells_row, 1.0)

    is_ctrl = symbols == "non-targeting"
    if int(is_ctrl.sum()) < 5:
        raise RuntimeError(f"{path.name}: found {int(is_ctrl.sum())} non-targeting rows")
    ctrl_w = n_cells_row[is_ctrl]
    ctrl_mean = (ctrl_w[:, None] * X[is_ctrl]).sum(axis=0) / ctrl_w.sum()
    n_ctrl = float(ctrl_w.sum())
    cpm_c = _cpm(ctrl_mean)
    log2_c = np.log2(cpm_c + 1.0)
    ctrl_log = np.log(cpm_c + 1.0)  # natural log, the feature's units

    # Group guide rows by target symbol.
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(symbols.tolist()):
        if s == "non-targeting":
            continue
        groups.setdefault(s, []).append(i)

    perts = sorted(groups)
    G = len(genes)
    P = len(perts)
    lfc = np.empty((P, G), dtype=np.float64)
    var = np.empty((P, G), dtype=np.float64)
    n_out = np.empty(P, dtype=np.float64)
    n_guides = np.empty(P, dtype=np.int16)
    for pi, pert in enumerate(perts):
        idxs = groups[pert]
        w = n_cells_row[idxs]
        rows = X[idxs]
        mu = (w[:, None] * rows).sum(axis=0) / w.sum()
        n_p = float(w.sum())
        lfc[pi] = np.log2(_cpm(mu) + 1.0) - log2_c
        v = poisson_var_log2(mu, n_p, ctrl_mean, n_ctrl)
        if len(idxs) >= 2:
            lib = np.maximum(rows.sum(axis=1, keepdims=True), 1e-12)
            guide_lfc = np.log2(rows / lib * 1e6 + 1.0) - log2_c
            between = guide_lfc.var(axis=0, ddof=1) / len(idxs)
            v = v + between
        var[pi] = np.maximum(v, 1e-8)
        n_out[pi] = n_p
        n_guides[pi] = len(idxs)

    # Biological known-answer: the targeted gene's own log2FC is a knockdown.
    gindex = {g: i for i, g in enumerate(genes.tolist())}
    own = []
    for pi, pert in enumerate(perts):
        j = gindex.get(pert)
        if j is not None:
            own.append(lfc[pi, j])
    own = np.asarray(own)
    median_own = float(np.median(own)) if len(own) else float("nan")
    if not np.isfinite(median_own) or median_own > -0.3:
        raise RuntimeError(f"{path.name}: median own-gene log2FC is {median_own}, expected a knockdown")

    return {
        "genes": genes.astype(str),
        "perts": np.array(perts, dtype=str),
        "lfc": lfc,
        "var": var,
        "n_cells": n_out,
        "n_guides": n_guides,
        "ctrl_log": ctrl_log,
        "n_ctrl": n_ctrl,
        "median_own_gene_lfc": median_own,
        "n_ctrl_rows": int(is_ctrl.sum()),
        "source": str(path),
    }


def merge_screens(a: dict, b: dict) -> dict:
    """Inverse-variance merge of two screens of the SAME line. One profile per target."""
    genes = sorted(set(a["genes"]).union(b["genes"]))
    perts = sorted(set(a["perts"]).union(b["perts"]))
    ga = {g: i for i, g in enumerate(a["genes"].tolist())}
    gb = {g: i for i, g in enumerate(b["genes"].tolist())}
    pa = {p: i for i, p in enumerate(a["perts"].tolist())}
    pb = {p: i for i, p in enumerate(b["perts"].tolist())}

    def gather(screen, pindex, gindex):
        lfc = np.full((len(perts), len(genes)), np.nan)
        var = np.full_like(lfc, np.nan)
        p_src = np.array([pindex.get(p, -1) for p in perts])
        g_src = np.array([gindex.get(g, -1) for g in genes])
        p_ok = p_src >= 0
        g_ok = g_src >= 0
        lfc[np.ix_(p_ok, g_ok)] = screen["lfc"][np.ix_(p_src[p_ok], g_src[g_ok])]
        var[np.ix_(p_ok, g_ok)] = screen["var"][np.ix_(p_src[p_ok], g_src[g_ok])]
        return lfc, var

    l1, v1 = gather(a, pa, ga)
    l2, v2 = gather(b, pb, gb)
    both = np.isfinite(l1) & np.isfinite(l2) & np.isfinite(v1) & np.isfinite(v2) & (v1 > 0) & (v2 > 0)
    w1 = np.zeros_like(l1)
    w2 = np.zeros_like(l2)
    w1[both] = 1.0 / v1[both]
    w2[both] = 1.0 / v2[both]
    denom = w1 + w2
    lfc = np.where(np.isfinite(l1), l1, l2)
    var = np.where(np.isfinite(v1), v1, v2)
    lfc[both] = (w1[both] * l1[both] + w2[both] * l2[both]) / denom[both]
    var[both] = 1.0 / denom[both]
    var = np.maximum(var, 1e-8)

    def pert_vec(screen, pindex, key):
        out = np.full(len(perts), np.nan)
        for i, p in enumerate(perts):
            pi = pindex.get(p)
            if pi is not None:
                out[i] = screen[key][pi]
        return out

    n1 = pert_vec(a, pa, "n_cells")
    n2 = pert_vec(b, pb, "n_cells")
    g1 = pert_vec(a, pa, "n_guides")
    g2 = pert_vec(b, pb, "n_guides")
    n_cells = np.nansum(np.vstack([n1, n2]), axis=0)
    n_guides = np.nansum(np.vstack([g1, g2]), axis=0)

    def ctrl_on(screen, gindex):
        out = np.full(len(genes), np.nan)
        for j, g in enumerate(genes):
            gi = gindex.get(g)
            if gi is not None:
                out[j] = screen["ctrl_log"][gi]
        return out

    c1 = ctrl_on(a, ga)
    c2 = ctrl_on(b, gb)
    # Weight control logs by control cell counts.
    both_c = np.isfinite(c1) & np.isfinite(c2)
    wa, wb = float(a["n_ctrl"]), float(b["n_ctrl"])
    ctrl = np.where(both_c, (wa * c1 + wb * c2) / (wa + wb), np.where(np.isfinite(c1), c1, c2))
    return {
        "genes": np.array(genes, dtype=str),
        "perts": np.array(perts, dtype=str),
        "lfc": lfc,
        "var": var,
        "n_cells": n_cells,
        "n_guides": n_guides.astype(np.int32),
        "ctrl_log": ctrl,
        "n_ctrl": wa + wb,
    }


def save_line(path: Path, rec: dict, meta: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        genes=np.asarray(rec["genes"]).astype(str),
        perts=np.asarray(rec["perts"]).astype(str),
        lfc=np.asarray(rec["lfc"], dtype=np.float32),
        var=np.asarray(rec["var"], dtype=np.float32),
        n_cells=np.asarray(rec["n_cells"], dtype=np.float32),
        n_guides=np.asarray(rec["n_guides"], dtype=np.int32),
        ctrl_log=np.asarray(rec["ctrl_log"], dtype=np.float32),
        n_ctrl=np.float64(rec["n_ctrl"]),
    )
    meta_path = path.with_suffix(".json")
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def summarise(rec: dict) -> dict:
    own = []
    gindex = {g: i for i, g in enumerate(rec["genes"].tolist())}
    for pi, pert in enumerate(rec["perts"].tolist()):
        j = gindex.get(pert)
        if j is not None and np.isfinite(rec["lfc"][pi, j]):
            own.append(float(rec["lfc"][pi, j]))
    return {
        "n_perts": int(len(rec["perts"])),
        "n_genes": int(len(rec["genes"])),
        "n_ctrl_cells": float(rec["n_ctrl"]),
        "median_n_cells": float(np.median(rec["n_cells"])),
        "median_own_gene_lfc": float(np.median(own)) if own else None,
        "n_perts_with_own_gene": int(len(own)),
        "median_abs_lfc": float(np.median(np.abs(rec["lfc"][np.isfinite(rec["lfc"])]))),
        "frac_lfc_negative": float(np.mean(rec["lfc"][np.isfinite(rec["lfc"])] < 0)),
    }


def main() -> None:
    root = Path(__file__).resolve().parent
    data = root / "work" / "data"
    out = root / "work" / "responses"
    files = {
        "k562_essential": data / "K562_essential_raw_bulk_01.h5ad",
        "k562_gwps": data / "K562_gwps_raw_bulk_01.h5ad",
        "rpe1": data / "rpe1_raw_bulk_01.h5ad",
    }
    built = {}
    for name, path in files.items():
        print(f"building {name} from {path.name}", flush=True)
        digest = sha256_file(path)
        rec = screen_to_arrays(path)
        built[name] = rec
        meta = {
            "line_or_screen": name,
            "url_note": "see RESULTS.md for the figshare URL",
            "sha256": digest,
            "nbytes": path.stat().st_size,
            **summarise(rec),
            "median_own_gene_lfc_screen": rec["median_own_gene_lfc"],
            "n_ctrl_rows": rec["n_ctrl_rows"],
        }
        if name != "k562_essential":
            # gwps and rpe1 are saved; essential is an input to the merged K562 line.
            save_line(out / f"{name}.npz", rec, meta)
        else:
            save_line(out / f"{name}.npz", rec, meta)
        print(json.dumps(meta), flush=True)

    print("merging K562 screens", flush=True)
    merged = merge_screens(built["k562_gwps"], built["k562_essential"])
    meta = {
        "line_or_screen": "k562",
        "merge": "inverse-variance weighted mean of k562_gwps and k562_essential; same cell line",
        "inputs": ["k562_gwps.npz", "k562_essential.npz"],
        **summarise(merged),
    }
    # Knockdown should survive the merge.
    if meta["median_own_gene_lfc"] is None or meta["median_own_gene_lfc"] > -0.3:
        raise RuntimeError(f"merged K562 own-gene median {meta['median_own_gene_lfc']}")
    save_line(out / "k562.npz", merged, meta)
    print(json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()

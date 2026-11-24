"""804 — 702's arm (PIE bulk blend, w 0.25, PIE-seen targets, norm-preserving) re-expressed on a K562 PLOLO fold.

On the board (702), PIE's predicted rows for contexts A/B/C were blended into the cap-6 + P2g bulk dump on the targets in
PIE's training split ("seen", 154/300), u = base/|base|, v = pie/|pie|, new = |base| * normalise((1-w) u + w v) over the
non-panel genes (experiments/700-pie-transfer/make_pie_dump700.py, `pie_effect` and `blend`, copied below verbatim).

On a K562 fold the held-out line is K562 and PIE's K562 rows (replogle_xdataset checkpoint: no Replogle rows in training or
evidence) exist only on K562's measured genes. The fold's truth lives on the same genes, so the fold can see only that part
of the change. The blend is therefore applied on the K562 axis's non-panel genes, norm-preserving there, and the row is
untouched elsewhere ("visible-part" emulation; 702 rule 3 labels the price PARTIAL: on A/B/C, 41 % of 702's change energy
lay off the K562 axis). The effect space uses K562's own control pseudobulk (data/donors/k562_gwps/b_ctrl.npy), as
make_pie_dump700 used each context's own controls.

    rows = pie_rows(pred_parquet, targets, official_genes)            # (n, G) effect rows on the official axis, 0 off-axis
    seen = seen_mask(targets)                                          # PIE train-split perturbations (700's pie_train_perts.json)
    new  = blend_axis(base_rows, rows, 0.25, axis_nonpanel_mask, seen)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TS = 5e4
MIN_CPM = 5.0
K562 = ROOT / "data/donors/k562_gwps"
SEEN = ROOT / "experiments/700-pie-transfer/work/pie_train_perts.json"


def pie_effect(lfc: np.ndarray, bc: np.ndarray) -> np.ndarray:
    """make_pie_dump700.pie_effect, verbatim."""
    p = np.expm1(bc)
    p = p / p.sum()
    lfc = np.where(1e6 * p >= MIN_CPM, lfc, 0.0)
    q = p * np.exp2(lfc)
    q = q / q.sum()
    return np.log1p(TS * q) - np.log1p(TS * p)


def k562_axis(official: list[str]):
    genes = (K562 / "genes.txt").read_text().split()
    bc = np.load(K562 / "b_ctrl.npy").astype(np.float64)
    pos = {g: i for i, g in enumerate(official)}
    on = np.array([g in pos for g in genes])
    return [g for g, o in zip(genes, on) if o], bc[on], np.array([pos[g] for g, o in zip(genes, on) if o])


def pie_rows(pred_parquet: Path, targets: list[str], official: list[str], dataset="replogle", context="k562") -> np.ndarray:
    import pyarrow.parquet as pq
    t = pq.read_table(pred_parquet)
    pg = json.loads(t.schema.metadata[b"pie"])["genes"][dataset]
    df = t.to_pandas()
    df = df[(df["dataset"] == dataset) & (df["context"] == context)].set_index("perturbation")
    kg, kbc, kpos = k562_axis(official)
    gpos = {g: i for i, g in enumerate(pg)}
    have = np.array([g in gpos for g in kg])
    gix = np.array([gpos.get(g, 0) for g in kg])
    out = np.zeros((len(targets), len(official)))
    for i, tg in enumerate(targets):
        if tg not in df.index:
            continue
        lf = np.where(have, np.asarray(df.loc[tg, "lfc_pred"], np.float64)[gix], 0.0)
        out[i, kpos] = pie_effect(lf, kbc)
    return out


def seen_mask(targets: list[str]) -> np.ndarray:
    seen = set(json.loads(SEEN.read_text()))
    return np.array([t in seen for t in targets])


def blend_axis(base: np.ndarray, pie: np.ndarray, w: float, cols: np.ndarray, rows: np.ndarray) -> np.ndarray:
    """make_pie_dump700.blend on the columns `cols` (bool) and the rows `rows` (bool) only; other entries untouched."""
    out = base.copy()
    for i in np.flatnonzero(rows):
        b, v = base[i, cols], pie[i, cols]
        nb, nv = np.linalg.norm(b), np.linalg.norm(v)
        if nb == 0 or nv == 0:
            continue
        m = (1.0 - w) * b / nb + w * v / nv
        nm = np.linalg.norm(m)
        if nm > 0:
            out[i, cols] = nb * m / nm
    return out


def change_share(base: np.ndarray, new: np.ndarray, cols: np.ndarray) -> float:
    """share of the change energy on `cols` (702 rule 3's statistic, the complement is 'outside the axis')."""
    d = new - base
    tot = float((d * d).sum())
    return float((d[:, cols] ** 2).sum()) / tot if tot > 0 else float("nan")

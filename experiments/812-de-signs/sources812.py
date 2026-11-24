"""812 — build skill812 input files (one npz per source and surface: targets, genes, pred) for the step-3 sources.

(a) the per-target mean of the non-H1 public lines' rows (bank812's work/h1banks; banks of one line averaged first, then
    lines averaged with equal weight; a gene counts where >= 1 line measures it):
    pubmean_raw   rows as they are (each bank's own units)
    pubmean_rms   each bank's rows divided by the bank's RMS over all its rows (units equalised)
    pubmean_ctr   as pubmean_rms after subtracting each bank's mean row over ALL its perturbations (per-line centring:
                  the line's shared programme removed, the target-specific part kept)
(b) bulk_sign     the production bulk delta's own sign on H1: run496's `ps_bulk` is the cap15 composition times a positive
                  per-target scalar, and sharpening sign(h)|h|^q keeps signs, so its sign is sign(cap15 effect) (323's /
                  618's effects.npz, h1_2025 never a donor), ungated (5-10 CPM genes included).
Also writes the KA inputs: clean4 (CLEAN4 P_down sign, same everywhere) and emitted_zs (the zs arm's realized lfc).

usage: uv run python experiments/812-de-signs/sources812.py   -> work/sources/<name>_<surface>.npz + sources812.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

OUT = SF.HERE / "work/sources"
OUT.mkdir(parents=True, exist_ok=True)
man = json.loads((SF.HERE / "work/h1banks/manifest.json").read_text())
banks = {b: np.load(SF.HERE / f"work/h1banks/{b}.npz", allow_pickle=False) for b in man}
summary = {}
for surf in ("h1test", "h1tv"):
    S = SF.load(surf, 0)
    targets, genes = S["targets"], S["genes"]
    T, G = len(targets), len(genes)
    ti = {t: i for i, t in enumerate(targets)}
    lines = sorted({man[b]["line"] for b in man})
    acc = {v: {ln: [np.zeros((T, G)), np.zeros((T, G))] for ln in lines} for v in ("raw", "rms", "ctr")}
    for b, z in banks.items():
        ln = man[b]["line"]
        rows, rms, mean = z["rows"].astype(np.float64), float(z["rms"]), z["mean"]
        for k, t in enumerate(z["targets"].astype(str)):
            if t not in ti:
                continue
            r = rows[k]
            m = np.isfinite(r)
            for v, x in (("raw", r), ("rms", r / rms), ("ctr", (r - mean) / rms)):
                acc[v][ln][0][ti[t], m] += x[m]
                acc[v][ln][1][ti[t], m] += 1
    cov = {}
    for v in ("raw", "rms", "ctr"):
        num = np.zeros((T, G))
        cnt = np.zeros((T, G))
        for ln in lines:
            s_, n_ = acc[v][ln]
            has = n_ > 0
            num[has] += s_[has] / n_[has]            # the line's mean over its banks
            cnt[has] += 1
        pred = np.where(cnt > 0, num / np.maximum(cnt, 1), np.nan)
        np.savez_compressed(OUT / f"pubmean_{v}_{surf}.npz", targets=targets, genes=genes, pred=pred.astype(np.float32))
        cov[v] = dict(targets_with_any_line=int((cnt.max(1) > 0).sum()), mean_lines_per_signed_cell=round(float(cnt[cnt > 0].mean()), 2),
                      signed_share=round(float((cnt > 0).mean()), 4))
    # lines per target
    lt = {ln: int((acc["raw"][ln][1].max(1) > 0).sum()) for ln in lines}
    # (b) bulk sign: ungated cap15 effects on the kept axis
    e = np.load(SF.SURF[surf]["eff"])
    assert list(e["targets"].astype(str)) == list(targets)
    pos = {g: i for i, g in enumerate(e["genes"].astype(str))}
    col = np.array([pos.get(g, -1) for g in genes])
    eff = np.full((T, G), np.nan)
    eff[:, col >= 0] = np.asarray(e["cap15"], float)[:, col[col >= 0]]
    np.savez_compressed(OUT / f"bulk_sign_{surf}.npz", targets=targets, genes=genes, pred=eff.astype(np.float32))
    c4 = np.where(np.isfinite(S["c4"]), np.where(S["c4"] > 0.5, -1.0, 1.0), np.nan)
    np.savez_compressed(OUT / f"clean4_{surf}.npz", targets=targets, genes=genes, pred=np.repeat(c4[None, :], T, 0).astype(np.float32))
    np.savez_compressed(OUT / f"emitted_zs_{surf}.npz", targets=targets, genes=genes, pred=S["ol"].astype(np.float32))
    summary[surf] = dict(coverage=cov, targets_per_line=lt, bulk_signed_share=round(float((np.nan_to_num(eff) != 0).mean()), 4))
    print(surf, json.dumps(summary[surf]), flush=True)
(SF.HERE / "work/sources812.json").write_text(json.dumps(summary, indent=1))

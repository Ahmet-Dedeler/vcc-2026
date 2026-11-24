"""812 — extract the H1 surfaces' targets (h1-test 85 + h1-trainval 126) from every NON-H1 public bank on disk (donor format:
delta.npy rows = perturbation means minus control means at target_sum 50k, log1p, as 811's ext811 reads them), onto H1's
5-CPM gene axis (10,779 genes). Read-only. VCC25 H1 (h1_2025, h1_vcc2025) and Pan's H1 (h1_pan2025) are never read.

Per bank b -> work/h1banks/<b>.npz:
  targets (n,)          the H1 targets the bank carries
  rows    (n, G) f32    their rows on the H1 axis, NaN where the bank does not measure the gene
  mean    (G,) f64      the bank's mean row over ALL its perturbations (its line programme; for per-line centring)
  rms     ()            the bank's RMS over all measured entries of all its rows (for scaling)
  line                  the cell line (banks of one line are averaged before lines are averaged)
manifest: work/h1banks/manifest.json (source shas, n targets, n genes measured).

usage: uv run python experiments/812-de-signs/bank812.py
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "work/h1banks"
ST = HERE / "work/stage"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
T223 = ROOT / "experiments/223-silent-zero-audit/work/tables"
BANKS = {  # bank -> (dir, line)
    "k562_gwps": (ROOT / "data/donors/k562_gwps", "K562"),
    "k562_essential": (ROOT / "data/donors/k562_essential", "K562"),
    "viperturb_k562": (ROOT / "data/donors/viperturb_k562", "K562"),
    "jiang_K562": (T225 / "jiang_K562", "K562"),
    "orion_hct116": (ROOT / "data/donors/orion_hct116", "HCT116"),
    "orion_hek293t": (ROOT / "data/donors/orion_hek293t", "HEK293T"),
    "nadig_jurkat": (ROOT / "data/donors/nadig_jurkat", "Jurkat"),
    "rpe1": (ROOT / "data/donors/rpe1", "RPE1"),
    "dld1_lowmoi": (T223 / "dld1_lowmoi", "DLD1"),
    "altair_a549": (T223 / "altair_a549", "A549"),
    "jiang_A549": (T225 / "jiang_A549", "A549"),
    "jiang_BXPC3": (T225 / "jiang_BXPC3", "BXPC3"),
    "jiang_HAP1": (T225 / "jiang_HAP1", "HAP1"),
    "jiang_HT29": (T225 / "jiang_HT29", "HT29"),
    "jiang_MCF7": (T225 / "jiang_MCF7", "MCF7"),
    "telohaec_pilot": (ROOT / "data/external/igvf-telohaec", "TeloHAEC"),
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


genes = np.loadtxt(ST / "genes_keep_h1test.txt", dtype=str)
gpos = {g: i for i, g in enumerate(genes)}
want = set(np.loadtxt(ST / "targets_h1test.txt", dtype=str)) | set(np.loadtxt(ST / "targets_h1tv.txt", dtype=str))
OUT.mkdir(parents=True, exist_ok=True)
man = {}
t0 = time.time()
for b, (D, line) in BANKS.items():
    if not (D / "delta.npy").exists():
        print(f"{b}: MISSING {D}", flush=True)
        continue
    bg = (D / "genes.txt").read_text().split()
    perts = (D / "perts.txt").read_text().split()
    delta = np.load(D / "delta.npy", mmap_mode="r")
    assert delta.shape == (len(perts), len(bg)), f"{b}: {delta.shape} vs {len(perts)} x {len(bg)}"
    cols = np.array([gpos.get(g, -1) for g in bg])
    ok = cols >= 0
    seen = set()
    for j in np.flatnonzero(ok):          # duplicated symbols: first occurrence wins
        if cols[j] in seen:
            ok[j] = False
        seen.add(cols[j])
    jj = np.flatnonzero(ok)
    # bank mean and RMS over all rows, chunked
    s = np.zeros(len(bg))
    ss = 0.0
    nn = 0
    for a in range(0, len(perts), 2000):
        blk = np.asarray(delta[a:a + 2000], np.float64)
        s += blk.sum(0)
        ss += float((blk[:, jj] ** 2).sum())
        nn += blk.shape[0] * jj.size
    mean = np.full(len(genes), np.nan)
    mean[cols[jj]] = (s / len(perts))[jj]
    pidx = {}
    for i, p in enumerate(perts):
        pidx.setdefault(p, i)
    T = sorted(want & set(pidx))
    rows = np.full((len(T), len(genes)), np.nan, np.float32)
    for k, t in enumerate(T):
        rows[k, cols[jj]] = np.asarray(delta[pidx[t]], np.float32)[jj]
    dst = OUT / f"{b}.npz"
    np.savez_compressed(dst, targets=np.asarray(T), rows=rows, mean=mean, rms=np.sqrt(ss / max(nn, 1)), line=line)
    man[b] = dict(dir=str(D.relative_to(ROOT)), line=line, n_perts=len(perts), n_h1_targets=len(T), n_genes_on_h1_axis=int(jj.size),
                  rms=float(np.sqrt(ss / max(nn, 1))), sources={f: sha(D / f) for f in ("delta.npy", "genes.txt", "perts.txt")},
                  out_sha256=sha(dst))
    print(f"{b:16s} {line:9s} perts {len(perts):6d}  H1 targets {len(T):4d}  genes {jj.size:6d}  rms {man[b]['rms']:.4f} ({time.time()-t0:.0f}s)", flush=True)
(OUT / "manifest.json").write_text(json.dumps(man, indent=1))

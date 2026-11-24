"""811 -- extract the EXTRA public banks' rows for the fold panels' targets (A/B/C + WORKING L1-L6; never L7/L8) onto the
official 18,533-gene axis, so the EC2 pricing job needs ~0.3 GB of compact inputs instead of whole banks.

The three composition donors (k562_gwps, orion_hct116, orion_hek293t) are NOT extracted: the EC2 job reads them from the
lab's data mirror exactly as 804's `donor_unit_mean` does.

Per bank -> work/banks/<bank>.npz: targets (n,), rows float32 (n x 18533; 0 where the bank does not measure the gene),
measured bool (18533,). Rows are raw delta values (no panel zeroing, no normalisation: the job does both per fold).
work/banks/manifest.json: per bank the source files' sha256, shape, n targets per panel, and the output sha256.

    bash experiments/602-council-gate/gate602.sh ext811 3 -- uv run python experiments/811-lolo-model/ext811.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work" / "banks"
GENES = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
GPOS = {g: i for i, g in enumerate(GENES)}
PANELS = {"ABC": ROOT / "data/challenge/pert_counts.csv",
          **{f"L{k}": ROOT / f"experiments/801-lab-folds/work/L{k}/bundle/pert_counts.csv" for k in range(1, 7)}}
T223 = ROOT / "experiments/223-silent-zero-audit/work/tables"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
SIMPLE = {   # bank -> donor-format dir (delta.npy rows indexed by perts.txt, columns by genes.txt)
    "viperturb_k562": ROOT / "data/donors/viperturb_k562",
    "h1_pan2025": ROOT / "experiments/115-gse295214-h1/work/compact-donor",
    "h1_vcc2025": ROOT / "data/donors/h1_2025_full",
    "kolf2": ROOT / "data/external/kolf2",
    "telohaec_pilot": ROOT / "data/external/igvf-telohaec",
    "dld1_lowmoi": T223 / "dld1_lowmoi",
    "altair_a549": T223 / "altair_a549",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    **{f"jiang_{ln}": T225 / f"jiang_{ln}" for ln in ("K562", "A549", "BXPC3", "HAP1", "HT29", "MCF7")},
}
CD4_DIRS = [ROOT / "data/external/cd4t/cloud-compact/cd4_compact"] + \
           [ROOT / f"experiments/801-lab-folds/work/L{k}/cd4/cd4_compact" for k in range(1, 7)]
CD4_CONDS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def panel(p: Path) -> list[str]:
    with open(p) as f:
        return [r["target_gene"] for r in csv.DictReader(f)]


def read_bank(dirs: list[Path], want: set[str]):
    """rows for the wanted targets from one or more donor-format dirs (first dir carrying a target wins; dirs that carry the
    same target must agree exactly -- asserted)."""
    rows, measured, src = {}, np.zeros(len(GENES), bool), {}
    for D in dirs:
        genes = (D / "genes.txt").read_text().split()
        perts = (D / "perts.txt").read_text().split()
        delta = np.load(D / "delta.npy", mmap_mode="r")
        assert delta.shape == (len(perts), len(genes)), f"{D}: delta {delta.shape} vs perts {len(perts)} genes {len(genes)}"
        cols = np.array([GPOS.get(g, -1) for g in genes])
        ok = cols >= 0
        # a gene symbol appearing twice in the bank axis: keep the first occurrence (deterministic)
        seen = set()
        for j in np.flatnonzero(ok):
            if cols[j] in seen:
                ok[j] = False
            seen.add(cols[j])
        measured[cols[ok]] = True
        pidx = {}
        for i, p in enumerate(perts):
            pidx.setdefault(p, i)
        for t in sorted(want & set(pidx)):
            r = np.zeros(len(GENES), np.float32)
            r[cols[ok]] = np.asarray(delta[pidx[t]], np.float32)[ok]
            if t in rows:
                assert np.array_equal(rows[t], r), f"{D}: {t} differs from an earlier dir"
                continue
            rows[t] = r
        src[str(D.relative_to(ROOT))] = {f: sha(D / f) for f in ("delta.npy", "genes.txt", "perts.txt")}
    return rows, measured, src


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pans = {k: panel(p) for k, p in PANELS.items()}
    want = set().union(*map(set, pans.values()))
    jobs = {b: [D] for b, D in SIMPLE.items()}
    for c in CD4_CONDS:
        jobs[c] = [D / c for D in CD4_DIRS]
    man = dict(panels={k: str(p.relative_to(ROOT)) for k, p in PANELS.items()}, n_wanted=len(want), banks={})
    for b, dirs in jobs.items():
        rows, measured, src = read_bank(dirs, want)
        T = sorted(rows)
        R = np.vstack([rows[t] for t in T]) if T else np.zeros((0, len(GENES)), np.float32)
        nz = np.abs(R).sum(1) > 0
        dst = OUT / f"{b}.npz"
        np.savez_compressed(dst, targets=np.asarray(T), rows=R, measured=measured)
        man["banks"][b] = dict(sources=src, n_targets=len(T), n_zero_rows=int((~nz).sum()), n_measured=int(measured.sum()),
                               coverage={k: len(set(v) & set(T)) for k, v in pans.items()}, out=str(dst.relative_to(ROOT)),
                               out_sha256=sha(dst))
        print(f"{b:16s} {len(T):5d} targets, {int(measured.sum()):5d} genes measured, zero rows {int((~nz).sum())}, "
              + " ".join(f"{k}:{v}" for k, v in man["banks"][b]["coverage"].items()) + f"  ({time.time() - t0:.0f}s)", flush=True)
    man["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30   # macOS: bytes
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1))
    print(f"done in {time.time() - t0:.0f}s, peak RSS {man['peak_rss_gib']:.2f} GiB", flush=True)


if __name__ == "__main__":
    sys.exit(main())

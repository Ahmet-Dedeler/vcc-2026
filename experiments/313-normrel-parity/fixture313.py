"""313 -- ONE shared synthetic fixture for the normrel parity test.  Deterministic, tiny, no bank is read.

The fixture is a three-donor screen on a 60-gene official axis with 20 panel targets (target i's own gene is
official column i, so columns 0..19 are the panel and 20..59 are the operator columns), built so that every
case the mandate names is present:

  (a) PerturbMatch present for most targets and ABSENT for `T003` and `T019`;
  (b) gene masks: `k562_gwps` covers official 0..49, `orion_hct116` all 60, `orion_hek293t` 10..59 --
      so no donor covers the whole axis the same way and PM covers 5..54;
  (c) `T014`/`T015` are missing from `orion_hek293t` (two carrying donors), `T016`/`T017` from BOTH Orion
      donors (one carrying donor), `T018` from `k562_gwps` (two donors, and the PM branch cannot fire);
  (d) `T005` is a NEAR-ZERO `orion_hct116` row (its delta is set to the mean of that donor's other rows, so the
      centred row is ~1e-17) and `T006` is a small-but-ordinary `orion_hek293t` row (0.05x its own median);
  (e) `T007` is a `k562_gwps` row at ~8x its own panel median -- the case the -1.5 exponent is meant to punish.

Everything is written to `work/fixture313.npz` together with the seed.  Nothing outside this experiment is touched.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

SEED = 3130907
N_GENES = 60
N_TARGETS = 20                      # official columns 0..19 are the panel (target i's own gene is column i)
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")

# donor gene coverage on the official axis (case b)
GENE_COVER = {
    "k562_gwps": (0, 50),
    "orion_hct116": (0, 60),
    "orion_hek293t": (10, 60),
}
# targets each donor does NOT carry (case c)
PERT_MISSING = {
    "k562_gwps": ("T018",),
    "orion_hct116": ("T016", "T017"),
    "orion_hek293t": ("T014", "T015", "T016", "T017"),
}
PM_GENES = (5, 55)                  # PerturbMatch observes official columns 5..54
PM_MISSING = ("T003", "T019")       # targets PerturbMatch does not carry (case a)

NEAR_ZERO = ("orion_hct116", "T005")     # exactly-degenerate row (case d, pathological)
SMALL_ROW = ("orion_hek293t", "T006")    # 0.05x median row (case d, ordinary)
BIG_ROW = ("k562_gwps", "T007")          # 8x median row (case e)


def official_genes(n_genes: int = N_GENES) -> list[str]:
    return [f"G{i:03d}" for i in range(n_genes)]


def targets() -> list[str]:
    # target i's own gene is official column i, so the panel is columns 0..19
    return [f"T{i:03d}" for i in range(N_TARGETS)]


def build() -> dict:
    rng = np.random.default_rng(SEED)
    genes = official_genes()
    tg = targets()
    # a shared cross-donor response component (so the donor rows are correlated, as real ones are) plus
    # donor-specific and target-specific structure
    shared = rng.normal(0.0, 0.35, size=(N_TARGETS, N_GENES))
    out: dict[str, np.ndarray] = {}
    for name in DONORS:
        lo, hi = GENE_COVER[name]
        cover = np.arange(lo, hi)
        keep = [t for t in tg if t not in PERT_MISSING[name]]
        delta = np.zeros((len(keep), len(cover)))
        for r, t in enumerate(keep):
            ti = tg.index(t)
            row = 0.55 * shared[ti, cover] + rng.normal(0.0, 0.30, size=cover.size)
            # the donor's own measured knockdown of the target's own gene lives on the PANEL column ti
            if lo <= ti < hi:
                row[ti - lo] = -2.4 + rng.normal(0.0, 0.4)
            delta[r] = row
        if name == BIG_ROW[0]:
            delta[keep.index(BIG_ROW[1])] *= 8.0
        if name == SMALL_ROW[0]:
            delta[keep.index(SMALL_ROW[1])] *= 0.05
        if name == NEAR_ZERO[0]:
            # the fixed point of `s = delta - delta.mean(0)`: setting a row to the mean of the OTHER rows makes
            # that row's centred value exactly zero in exact arithmetic (~1e-17 in float64).
            r = keep.index(NEAR_ZERO[1])
            other = np.delete(np.arange(len(keep)), r)
            delta[r] = delta[other].mean(0)
        out[f"delta__{name}"] = delta
        out[f"genes__{name}"] = np.asarray([genes[i] for i in cover])
        out[f"perts__{name}"] = np.asarray(keep)
        out[f"b_ctrl__{name}"] = np.abs(rng.normal(2.0, 0.5, size=cover.size))

    # PerturbMatch: its own delta over its own gene window, carrying every target but PM_MISSING
    pm_cover = np.arange(*PM_GENES)
    pm_keep = [t for t in tg if t not in PM_MISSING]
    pm_delta = np.zeros((len(pm_keep), pm_cover.size))
    for r, t in enumerate(pm_keep):
        ti = tg.index(t)
        pm_delta[r] = 0.45 * shared[ti, pm_cover] + rng.normal(0.0, 0.45, size=pm_cover.size)
    out["pm_delta"] = pm_delta
    out["pm_genes"] = np.asarray([genes[i] for i in pm_cover])
    out["pm_perts"] = np.asarray(pm_keep)

    out["official"] = np.asarray(genes)
    out["targets"] = np.asarray(tg)
    out["donors"] = np.asarray(DONORS)
    out["panel_positions"] = np.arange(N_TARGETS, dtype=np.int64)
    out["seed"] = np.asarray([SEED], dtype=np.int64)
    return out


def write() -> Path:
    WORK.mkdir(parents=True, exist_ok=True)
    fx = build()
    path = WORK / "fixture313.npz"
    np.savez(path, **fx)
    return path


def sha256_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


if __name__ == "__main__":
    p = write()
    print(f"wrote {p} sha256 {sha256_of(p)}")

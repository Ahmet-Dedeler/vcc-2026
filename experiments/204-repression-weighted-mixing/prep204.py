"""204 — shared preparation: production-frame donor effects, cached.

Everything downstream (attenuation estimation, weighted mixing, scoring) reads the SAME
object the production consensus averages: `E062.prepare_main_source`'s output, i.e. the
donor's delta rows mapped onto the 18,533-gene official axis, centred on that bank's own
gene mean over ALL its rows, rescaled by `EFFECT_SCALE / raw_median_norm`, with the 300
panel columns zeroed.

Two target sets are cached:
  `fold:<held>`  the PLOLO fold's own target list, for scoring (162's `build_fold`)
  `shared`       panel targets carried by all three donors + a fixed draw of shared
                 non-panel targets, for estimating attenuation with precision

`raw_median_norm` is carried alongside, because the production frame divides it out and
whether that division is itself a de-attenuation is one of the questions.
"""
from __future__ import annotations

import json
import os
import resource
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402
from plolo import DONOR_PATHS, E062, OFFICIAL, PANEL, PANEL_SET  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
SEED = 20260903
N_SHARED_EXTRA = 2000


def rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


_BANKS: dict[str, object] = {}


def bank(name: str):
    if name not in _BANKS:
        _BANKS[name] = E062.Bank.load(name, DONOR_PATHS[name])
    return _BANKS[name]


def shared_targets() -> tuple[str, ...]:
    """Panel targets in all three donors, then a fixed draw of shared non-panel ones."""
    cache = WORK / "shared-targets.json"
    if cache.exists():
        return tuple(json.loads(cache.read_text()))
    sets = [set(bank(d).pert_pos) for d in DONORS]
    common = set.intersection(*sets)
    panel = [g for g in PANEL if g in common]
    extra = sorted(common - PANEL_SET)
    rng = np.random.default_rng(SEED)
    draw = list(rng.choice(np.asarray(extra), size=min(N_SHARED_EXTRA, len(extra)),
                           replace=False)) if extra else []
    out = tuple(panel) + tuple(sorted(str(x) for x in draw))
    cache.write_text(json.dumps(list(out)))
    print(f"shared: {len(panel)} panel + {len(draw)} extra = {len(out)}", flush=True)
    return out


def prepared(name: str, targets: tuple[str, ...], tag: str) -> dict:
    """`prepare_main_source` on the shipped table, cached to disk as float32 + scalars."""
    cache = WORK / f"prep-{tag}-{name}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if tuple(str(t) for t in z["targets"]) == targets:
            return dict(effects=np.asarray(z["effects"], dtype=np.float64),
                        target_available=np.asarray(z["avail"], dtype=bool),
                        mean_on_axis=np.asarray(z["mean_on_axis"], dtype=np.float64),
                        raw_median_norm=float(z["raw_median_norm"]),
                        mapped_gene_count=int(z["mapped_gene_count"]))
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)
    p = E062.prepare_main_source(bank(name), targets, OFFICIAL, panel_positions,
                                 chunk_rows=256, max_rss_gib=9.0)
    np.savez_compressed(cache, effects=np.asarray(p.effects, dtype=np.float32),
                        avail=np.asarray(p.target_available, dtype=bool),
                        mean_on_axis=np.asarray(p.mean_on_axis, dtype=np.float64),
                        raw_median_norm=float(p.raw_median_norm),
                        mapped_gene_count=int(p.mapped_gene_count),
                        targets=np.asarray(targets))
    print(f"[prep/{tag}/{name}] {p.effects.shape} raw_median {p.raw_median_norm:.4f} "
          f"mapped {p.mapped_gene_count} rss {rss_gib():.1f}G", flush=True)
    return dict(effects=np.asarray(p.effects, dtype=np.float64),
                target_available=np.asarray(p.target_available, dtype=bool),
                mean_on_axis=np.asarray(p.mean_on_axis, dtype=np.float64),
                raw_median_norm=float(p.raw_median_norm),
                mapped_gene_count=int(p.mapped_gene_count))


if __name__ == "__main__":
    tg = shared_targets()
    for d in DONORS:
        p = prepared(d, tg, "shared")
        print(d, p["effects"].shape, p["raw_median_norm"],
              int(p["target_available"].sum()), f"rss {rss_gib():.1f}G", flush=True)

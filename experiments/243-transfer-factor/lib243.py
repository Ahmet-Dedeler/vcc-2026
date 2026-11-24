"""243 — shared constants and artifact loaders.  Every number here is READ from an artifact
on disk and the path is recorded next to it; nothing is re-derived that already exists.

Read-only everywhere except `experiments/243-transfer-factor/work/`.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
CTX = ("A", "B", "C")
TS = 50_000.0
N_GENES = 18533

# --- 197's inversion constants ------------------------------------------------------------
# experiments/197-cosine-inversion/work/invert197.json  (sumD from gate197.json G6, K from
# experiments/096-mse-floor/work/board_solve.npz)
P197 = ROOT / "experiments/197-cosine-inversion/work"
INV197 = json.loads((P197 / "invert197.json").read_text())
SUMD = float(INV197["sumD"])            # 0.22914103411323902
K_CAP = float(INV197["K"])              # 0.24652868897523508
FLOOR0 = float(INV197["floor0"])        # 1.023611085832704  = 1 + mean_c Kc_c / sumD
COEFF = float(INV197["coeff"])          # 2.0744945316971    = 2 sqrt(K / sumD)
KC = {"A": 0.005083962657824914, "B": 0.005227704077951962, "C": 0.005919139136949927}  # 096 (e)
# 180's 46-construct NTC jackknife C_p, the cap per row (experiments/211-amplitude-reopened/
# work/bx211.json `cp46`, same values as 191's claim191.json `cap`)
CP46 = {"A": 14.439751080335085, "B": 14.805734135005, "C": 16.439650143279888}

# 215's one fitted constant: d(raw pds)/d(cos) on a14 -> a15 (chan215.py line 50)
PDS_PER_COS_215 = (0.835987 - 0.832702) / (
    INV197["ours"]["a15"]["cos"] - INV197["ours"]["a14"]["cos"])

# --- the board rows, exact, from the leaderboard snapshots -------------------------------
# data/leaderboard/<snapshot>.json, entry id -> our arm.  Raw members to full precision.
BOARD_IDS = {
    "a14": ("0NvkHnMPzLFaJl8phAxS", "lb_20260830T1555.json"),
    "rc_callset": ("31oB9yWy8CSNBDEkY4iW", "2026-09-04T0013Z.json"),
    "cd4f034": ("QSe3JlNUYZ2BzvlcsL14", "2026-09-04T0313Z.json"),
    "h1union": ("iUziqE1gewfcPKzYUDaB", "2026-09-04T0913Z.json"),
    "rc_hybrid": ("nQUJZhZb3MW2uFJHBX5e", "2026-09-03T1610Z.json"),
}
# a15 (08-30, 0.1740865) is not in any snapshot on disk (the board shows only a team's latest
# entry and a23 overwrote it before the next snapshot).  docs/BOARD-LEDGER.md records its raw
# `mse` 0.873805 / `pds` 0.835987, and rc_hybrid / rc_callset reproduce both to 16 s.f.
# (197 §5 falsifier 4, 202 BOARD-RESULT).  rc_callset's row therefore IS a15's profile row.


def board_row(arm: str) -> dict:
    eid, snap = BOARD_IDS[arm]
    d = json.loads((ROOT / "data/leaderboard" / snap).read_text())
    rows = d["entries"] if isinstance(d, dict) else d
    r = next(x for x in rows if x["id"] == eid)
    return dict(
        id=eid, snapshot=snap, overall=r["score_avg"],
        raw=dict(pds=r["pds_cosine"], mse=r["expr_mse_unbiased_capped_norm"],
                 nmae=r["de_wilcoxon_lfc_nmae"], fid=r["de_wilcoxon_direction_fidelity_yield_raw"],
                 reach=r["de_wilcoxon_direction_reach_raw"], jac=r["de_wilcoxon_sig_jaccard"]),
        scaled=dict(pds=r["score_pds"], mse=r["score_mse"], nmae=r["score_nmae"],
                    fid=r["score_fid"], reach=r["score_reach"], jac=r["score_jac"]),
    )


def floor_cos(raw_mse: float) -> float:
    """197 §2.1: cos >= (1 + Kc/sumD - raw) / (2 sqrt(K/sumD)), the one-sided floor."""
    return (FLOOR0 - raw_mse) / COEFF


# --- PLOLO artifacts ----------------------------------------------------------------------
def sweep216():
    rows = json.loads((ROOT / "experiments/216-cd4-hard/work/sweep216.json").read_text())
    return {(r["fold"], r["arm"]): r for r in rows}


def sweep207():
    rows = json.loads((ROOT / "experiments/207-cd4-dose/work/sweep207.json").read_text())
    return {(r["fold"], r["arm"]): r for r in rows}


def compose218():
    rows = json.loads(
        (ROOT / "experiments/218-independent-donors/work/compose218.json").read_text())
    return {(r["fold"], r["arm"]): r for r in rows}


def repro222():
    return json.loads((ROOT / "experiments/222-h1-union-build/work/repro222.json").read_text())


# The functional keys as the PLOLO instrument writes them (162 plolo.py `endpoints`,
# 216 sweep216.py / 218 resid218.py `score_effects`):
#   cos_agg      = sqrt(1 - sse) on the a23 (mass-fixed) tag  == 197's energy-weighted aggregate
#   a15_cos_agg  = the same on the a15 (no mass fix) tag
#   a23_cos      = macro (per-target mean) cosine, a23 tag    == 222's "macro"
#   a15_cos      = macro cosine, a15 tag                      == 212's functional (216 §4.1)
FUNCTIONALS = {
    ("cos_agg", "a23"): "cos_agg",
    ("cos_agg", "a15"): "a15_cos_agg",
    ("macro", "a23"): "a23_cos",
    ("macro", "a15"): "a15_cos",
}
PDS_KEYS = {"a23": "pds", "a15": "a15_pds"}


def rel_gain(before: dict, after: dict, key: str, folds, how: str) -> float:
    """Relative gain of a functional across a fold set.

    how = "pooled"  : ratio of fold-means      (216 §4.4 / board216.py: m(f0.34)/m(f0.1) - 1)
    how = "perfold" : mean of per-fold ratios  (222 repro222.py: mean_f (x_f / base_f - 1))
    """
    b = np.array([before[f][key] for f in folds], dtype=float)
    a = np.array([after[f][key] for f in folds], dtype=float)
    if how == "pooled":
        return float(a.mean() / b.mean() - 1.0)
    if how == "perfold":
        return float(np.mean(a / b - 1.0))
    raise ValueError(how)


def d_pds(before: dict, after: dict, key: str, folds) -> float:
    return float(np.mean([after[f][key] - before[f][key] for f in folds]))


def rss_gb() -> float:
    import resource, sys
    r = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return r / 1024**3 if sys.platform == "darwin" else r / 1024**2

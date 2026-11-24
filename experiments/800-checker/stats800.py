"""800 — the checker's statistics: per-target arrays in, board-unit member deltas with paired target-bootstrap intervals out.

Engines (801's bulk folds, 803's DE surfaces) write one npz per (fold, arm[, seed]) in the STANDARD FORMAT below; this module
never imports an engine. Everything here is light (MBs), so it runs on the Mac.

STANDARD FORMAT (npz):
  targets            (n,) str
  bulk surfaces      x (F, n), a (F, n): <P, T> and <P, P> per frame on the scored genes; d (n,): <T, T>;
                     pds (F, n): literal per-target PDS; pdsl (F, n): board-emulated PDS (lambda 0.5); K (F, n) optional
  DE surfaces        <member>_<regime> (n,) for member in fid, jac, reach, nmae and regime in f1, f05, f035, cmp
                     (NaN where the member is undefined for that target, e.g. nmae's gate < 10); K (n,) optional

UNITS (DESIGN.md section 5). Board raw delta -> scaled: x INV; overall = sum / 6.
  pds  : fold delta of the literal per-target PDS, read as a board raw delta (calibrated later by the known answers)
  mse  : 518's identity, +1 % relative aggregate cosine X  ->  -0.0041 board raw mse (k 1)
  DE   : the surface's raw member delta at the chosen regime (619: f 0.35 is board-scale for nmae), read as board raw
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

INV = dict(pds=2.161362660432289, mse=1.046637855109273,          # 702 vs ALLSW ledger pair (effective board 1/span)
           fid=3.32482, jac=2.72081, reach=1.12927, nmae=1.64787)  # 192's board fit (lib521.INV)
SIGN = dict(pds=+1, mse=-1, fid=+1, jac=+1, reach=+1, nmae=-1)     # +1: higher raw is better
MSE_PER_PCT_X = -0.0041                                             # 518's identity (522 / 706 units, k 1)
DE_MEMBERS = ("fid", "jac", "reach", "nmae")
REGIMES = ("f1", "f05", "f035", "cmp")


def overall_of(raw: dict) -> float:
    """board-overall delta of raw member deltas (missing members count 0)."""
    return float(sum(SIGN[m] * raw[m] * INV[m] for m in raw if m in INV and np.isfinite(raw[m])) / 6.0)


# ------------------------------------------------------------------------------------------------- per-fold statistics
def _bulk_X(z, idx):
    """frame-averaged aggregate cosine over the targets idx (row_score's identity; 522's boot522.cos_agg)."""
    x, a, d = z["x"][:, idx], z["a"][:, idx], z["d"][idx]
    D = d.sum()
    return float(np.mean([x[j].sum() / math.sqrt(a[j].sum() * D) for j in range(x.shape[0])]))


def bulk_fold_raw(base, cand, idx=None) -> dict:
    """raw member deltas on one bulk fold (cand vs base, the same targets): r (fraction), pds, pdsl, and mse via 518."""
    n = base["d"].size
    idx = np.arange(n) if idx is None else idx
    r = _bulk_X(cand, idx) / _bulk_X(base, idx) - 1.0
    dp = float(np.mean(cand["pds"][:, idx] - base["pds"][:, idx]))
    dpl = float(np.mean(cand["pdsl"][:, idx] - base["pdsl"][:, idx])) if "pdsl" in cand and "pdsl" in base else float("nan")
    return dict(r=r, pds=dp, pdsl=dpl, mse=MSE_PER_PCT_X * 100.0 * r)


def de_fold_raw(base_seeds: list, cand_seeds: list, regime: str, idx=None) -> dict:
    """raw DE member deltas on one DE surface: per target, cand - base averaged over the paired seeds, then the member
    mean over the targets where both are defined (NaN-skip), as the scorer averages each member over its own cohort."""
    out = {}
    for m in DE_MEMBERS:
        k = f"{m}_{regime}"
        if k not in base_seeds[0]:
            continue
        B = np.nanmean(np.vstack([np.asarray(s[k], float) for s in base_seeds]), 0)
        C = np.nanmean(np.vstack([np.asarray(s[k], float) for s in cand_seeds]), 0)
        if idx is not None:
            B, C = B[idx], C[idx]
        ok = np.isfinite(B) & np.isfinite(C)
        out[m] = float(np.mean(C[ok]) - np.mean(B[ok])) if ok.any() else float("nan")
    return out


# ------------------------------------------------------------------------------------------------- pooled + bootstrap
@dataclass
class Fold:
    fid: str                      # e.g. "L3/orion_hct116" or "h1tv"
    kind: str                     # "bulk" | "de"
    base: object                  # bulk: dict of arrays; de: list of per-seed dicts
    cand: object
    regime: str = "f035"
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int((self.base["d"] if self.kind == "bulk" else self.base[0][f"nmae_{self.regime}"]).size)

    def raw(self, idx=None) -> dict:
        if self.kind == "bulk":
            return bulk_fold_raw(self.base, self.cand, idx)
        return de_fold_raw(self.base, self.cand, self.regime, idx)


def _fold_overall(raw: dict, kind: str, members: tuple) -> float:
    if kind == "bulk":
        return overall_of({m: raw[m] for m in ("pds", "mse") if m in members})
    return overall_of({m: raw[m] for m in DE_MEMBERS if m in members and m in raw})


def price(folds: list[Fold], members=("pds", "mse", "fid", "jac", "nmae"), n_boot=4000, seed=800,
          top_out=3) -> dict:
    """Pooled price over folds (equal weight per fold) with a paired target bootstrap (every fold resampled independently
    in each draw). `members` chooses what enters `overall` (reach is reported but, by 619's finding, not priced for
    delivery levers unless the caller includes it)."""
    rng = np.random.default_rng(seed)
    pts = []
    for f in folds:
        raw = f.raw()
        pts.append(dict(fold=f.fid, kind=f.kind, n=f.n, raw=raw, overall=_fold_overall(raw, f.kind, members)))
    point = float(np.mean([p["overall"] for p in pts]))
    draws = np.empty(n_boot)
    per_member = {m: np.empty(n_boot) for m in ("pds", "mse", "r", *DE_MEMBERS)}
    for b in range(n_boot):
        ov, mm = [], {m: [] for m in per_member}
        for f in folds:
            idx = rng.integers(0, f.n, f.n)
            raw = f.raw(idx)
            ov.append(_fold_overall(raw, f.kind, members))
            for m in per_member:
                if m in raw:
                    mm[m].append(raw[m])
        draws[b] = np.mean(ov)
        for m in per_member:
            per_member[m][b] = np.mean(mm[m]) if mm[m] else np.nan
    lo, hi = np.quantile(draws, [0.025, 0.975])
    # leave-top-k-targets-out per fold: drop the k targets with the largest positive contribution to the fold's overall
    lto = []
    for f in folds:
        n = f.n
        if n > 600:
            continue
        # the fold's overall without target i; the targets whose removal lowers it most are its top contributors
        loo = np.array([_fold_overall(f.raw(np.delete(np.arange(n), i)), f.kind, members) for i in range(n)])
        top = np.argsort(loo)[:top_out]
        lto.append(_fold_overall(f.raw(np.delete(np.arange(n), top)), f.kind, members))
    mem = {}
    for m, v in per_member.items():
        if np.all(np.isnan(v)):
            continue
        pt = float(np.nanmean([p["raw"][m] for p in pts if m in p["raw"]]))
        q = np.nanquantile(v, [0.025, 0.975])
        mem[m] = dict(raw=pt, ci95=[float(q[0]), float(q[1])],
                      overall=(SIGN[m] * pt * INV[m] / 6.0) if m in INV else None)
    return dict(point=point, ci95=[float(lo), float(hi)], sd=float(np.std(draws)), p_le0=float(np.mean(draws <= 0)),
                draws=draws,
                folds_pos=int(sum(p["overall"] > 0 for p in pts)), n_folds=len(pts),
                leave_top_out_mean=float(np.mean(lto)) if lto else None, per_fold=pts, members=mem,
                n_boot=n_boot, seed=seed, overall_members=list(members))


def combine(blocks: list[dict]) -> dict:
    """the sum of independent blocks (different surfaces: bulk folds + DE surfaces): points add, bootstrap draws add."""
    D = np.sum([b["draws"] for b in blocks], axis=0)
    lo, hi = np.quantile(D, [0.025, 0.975])
    return dict(point=float(sum(b["point"] for b in blocks)), ci95=[float(lo), float(hi)], sd=float(np.std(D)),
                p_le0=float(np.mean(D <= 0)), draws=D, parts=len(blocks))


# ------------------------------------------------------------------------------------------------- best-of-N, calibration
def e_max_normal(n: int) -> float:
    """E[max of n iid standard normals] (exact for n <= 5, Blom's approximation above)."""
    table = {1: 0.0, 2: 0.5641895835, 3: 0.8462843753, 4: 1.0293753730, 5: 1.1629644736}
    if n in table:
        return table[n]
    from scipy.stats import norm
    return float(norm.ppf((n - 0.375) / (n + 0.25)))


def best_of_n(res: dict, n_tries: int) -> dict:
    """AGENTS 11.1c.4: shift a family's winner down by e(N) x its bootstrap sd."""
    s = e_max_normal(max(1, int(n_tries))) * res["sd"]
    return dict(n_tries=int(n_tries), shift=s, point=res["point"] - s, ci95=[res["ci95"][0] - s, res["ci95"][1] - s])


def calibrate(rows: list[dict]) -> dict:
    """Over-read ratio per class from the known answers: rho = checker / board (rows with |board| >= 3e-4 and the same
    sign only enter the ratio; sign failures are counted separately). Leave-one-out spread of the median ratio."""
    out = {}
    for cls in sorted({r["class"] for r in rows}):
        R = [r for r in rows if r["class"] == cls and abs(r["board"]) >= 3e-4]
        rat = [r["checker"] / r["board"] for r in R if np.sign(r["checker"]) == np.sign(r["board"])]
        sign_fail = sum(np.sign(r["checker"]) != np.sign(r["board"]) for r in R)
        if not rat:
            out[cls] = dict(n=len(R), rho=None, sign_fail=int(sign_fail))
            continue
        loo = [float(np.median(rat[:i] + rat[i + 1:])) for i in range(len(rat))] if len(rat) > 1 else [rat[0]]
        out[cls] = dict(n=len(R), rho=float(np.median(rat)), rho_range=[float(min(rat)), float(max(rat))],
                        rho_loo=[float(min(loo)), float(max(loo))], sign_fail=int(sign_fail))
    return out


def project(res: dict, cal: dict | None) -> dict:
    """board projection = checker / rho, with rho's full observed range folded into the interval (DESIGN 7)."""
    if not cal or cal.get("rho") is None:
        return dict(point=None, ci95=None, note="no calibration for this class: the range must include 0 (702 rule 4)")
    lo, hi = res["ci95"]
    rmin, rmax = cal["rho_range"]
    cands = [lo / rmin, lo / rmax, hi / rmin, hi / rmax]
    return dict(point=res["point"] / cal["rho"], ci95=[min(cands), max(cands)], rho=cal["rho"])


def decide(lo: float, point: float, hi: float, harm: float = -0.0005) -> str:
    """the ship rule: price from the lower end (702 rule 1)."""
    if lo > 0:
        return "SHIP"
    if hi < 0 or point <= harm:
        return "HARM"
    return "NULL" if point <= 0 else "WEAK"

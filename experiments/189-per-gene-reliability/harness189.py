"""189 — the h1-test instrument, extended with the per-gene prior.

Wraps `experiments/172-sign-confidence-key/harness172.py` UNCHANGED so that the four
known-answer rows (production 0.18283, `|H|` 0.23533, sig-first 0.21914, oracle 0.64553) and
INBOX-4 §1's two decile rows gate every screen in this directory too.  Adds:

  * `theta`, `raw`, `eb` -- the per-gene prior mapped onto h1-test's 10,779-gene DE axis
    (10,725 carried; the 54 uncovered genes take the prior mean, theta = 0).
  * `reach_vec(key)` -- debench's `reach` decomposed per target, so a PAIRED bootstrap is
    available.  Its mean is asserted equal to `debench.Surface.reach` to 1e-12.
  * `qw(weight)` -- the |p|-weighted sign accuracy over the reference DE pool, the quantity
    187 §5 identifies as the only channel through which `nmae` moves.

DONOR AUDIT (mandate falsifier 3).  The prior is fitted on HCT116 / HEK293T / K562 / RPE1 /
JURKAT.  `data/donors/h1_2025` is BANNED in `common189.py` and is never opened.  h1-test's
recipient line H1 is not a fitting line, so the deployment fold is leave-H1-out by
construction; `audit()` prints the donor-to-surface TARGET overlap for the record.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "experiments", "172-sign-confidence-key"))

from harness172 import H1Test as _H1Test, check as _check  # noqa: E402

from common189 import WORK, load_table  # noqa: E402

H_ORACLE = 0.61862      # 172 §4.2, measured on this instrument (INBOX quotes 0.62119)


class Surface189(_H1Test):
    def __init__(self, prior_file: str = "prior189.npz") -> None:
        super().__init__()
        z = np.load(os.path.join(WORK, prior_file), allow_pickle=True)
        pg = {str(g): i for i, g in enumerate(z["genes"])}
        idx = np.array([pg.get(g, -1) for g in self.genes])
        hit = idx >= 0
        self.prior_covered = hit
        self.theta = np.zeros(len(self.genes))
        self.theta[hit] = z["theta"][idx[hit]]
        self.raw = np.full(len(self.genes), np.nan)
        self.raw[hit] = z["raw"][idx[hit]]
        self.eb = np.full(len(self.genes), np.nan)
        self.eb[hit] = z["eb"][idx[hit]]
        self.n_prior = np.zeros(len(self.genes))
        self.n_prior[hit] = z["n_all"][idx[hit]]
        # everything the reach walk actually reads
        self.n_real = self.S.real_sig.sum(1)
        self.scored = np.array([self.S.real_sig[p].sum() > 0
                                for p in range(len(self.perts))])

    # ---- endpoints -------------------------------------------------------------------

    def reach_vec(self, key=None, call=None) -> np.ndarray:
        """Per-target R_p, in debench's own order, for the targets debench scores."""
        S = self.S
        out = []
        for p in range(len(self.perts)):
            pool = np.flatnonzero(S.real_sig[p])
            if pool.size == 0:
                continue
            pool = pool[S._order(p, pool, key, call)]
            pool = pool[S.adj[p, pool]]
            if pool.size == 0:
                out.append(0.0)
                continue
            ok = np.sign(S.pred_lfc[p, pool]) == np.sign(S.real_lfc[p, pool])
            purity = np.cumsum(ok) / np.arange(1, pool.size + 1)
            deep = np.flatnonzero(purity >= 0.9)
            out.append((int(deep[-1] + 1) if deep.size else 0) / int(S.real_sig[p].sum()))
        return np.asarray(out, dtype=np.float64)

    def qw(self, w=None, pool_only: bool = True) -> float:
        """|p|-weighted sign accuracy over the reference DE pool (187 §5's `q_w`).

        `w` is an optional per-gene multiplier on the |p| weights -- the counterfactual
        "what would the credit-weighted accuracy be if 188 re-shaped magnitudes by w".
        """
        S = self.S
        num = den = 0.0
        for p in range(len(self.perts)):
            idx = np.flatnonzero(self.pool[p] if pool_only else S.adj[p])
            if idx.size == 0:
                continue
            wt = np.abs(S.pred_lfc[p, idx])
            if w is not None:
                wt = wt * np.asarray(w)[idx]
            ok = (np.sign(S.pred_lfc[p, idx]) == np.sign(S.real_lfc[p, idx])).astype(float)
            num += float(np.sum(wt * ok))
            den += float(np.sum(wt))
        return num / max(den, 1e-300)

    def map_gene(self, genes, values, fill: float = 0.0) -> np.ndarray:
        """Put a per-gene vector defined on `genes` onto h1-test's DE gene axis."""
        pg = {str(g): i for i, g in enumerate(genes)}
        idx = np.array([pg.get(g, -1) for g in self.genes])
        out = np.full(len(self.genes), fill, dtype=np.float64)
        hit = idx >= 0
        out[hit] = np.asarray(values)[idx[hit]]
        return out

    def mass_ratio(self, w=None) -> float:
        """K = sum|p| / sum|t| over the reference DE pool (187's `K`), under weights `w`."""
        S = self.S
        num = den = 0.0
        for p in range(len(self.perts)):
            idx = np.flatnonzero(self.pool[p])
            if idx.size == 0:
                continue
            wt = np.abs(S.pred_lfc[p, idx])
            if w is not None:
                wt = wt * np.asarray(w)[idx]
            num += float(wt.sum())
            den += float(np.abs(S.real_lfc[p, idx]).sum())
        return num / max(den, 1e-300)


def paired_bootstrap(a: np.ndarray, b: np.ndarray, n: int = 20000, seed: int = 20260903):
    """172 §4.1's instrument: resample TARGETS, keep both keys on the same draw."""
    rng = np.random.default_rng(seed)
    d = a - b
    idx = rng.integers(0, len(d), size=(n, len(d)))
    boot = d[idx].mean(1)
    return (float(d.mean()), float(boot.std(ddof=1)),
            float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)),
            float((boot > 0).mean()))


def audit(T: Surface189) -> dict:
    """Donor-to-surface target overlap, printed before every screen (mandate falsifier 3)."""
    from common189 import BANNED, TABLES
    tg = set(T.perts)
    rows = []
    print(f"{'donor':<16}{'LINE':<9}{'targets':>9}{'h1-test targets covered':>26}   status")
    for tid in TABLES:
        if tid in BANNED:
            print(f"{tid:<16}{TABLES[tid][1]:<9}{'--':>9}{'85/85 (100.0%)':>26}   "
                  f"BANNED, never opened (AGENTS.md trap #1)")
            rows.append(dict(donor=tid, covered=85, banned=True))
            continue
        tab = load_table(tid)
        if tab is None:
            continue
        c = len(set(tab["perts"]) & tg)
        used = tid in ("orion_hct116", "orion_hek293t", "k562_gwps", "rpe1", "nadig_jurkat")
        print(f"{tid:<16}{tab['line']:<9}{len(tab['perts']):9d}"
              f"{f'{c}/85 ({100*c/85:.1f}%)':>26}   "
              f"{'FITTING LINE -- gene-level only, no target identity' if used else 'unused'}")
        rows.append(dict(donor=tid, covered=c, banned=False, used=used))
    return dict(rows=rows)


def check(T: Surface189) -> None:
    _check(T)
    r = T.reach_vec()
    assert abs(float(r.mean()) - T.reach()) < 1e-12, "reach_vec != debench reach"
    assert len(r) == int(T.scored.sum()), (len(r), int(T.scored.sum()))
    print(f"  OK  reach_vec                {r.mean():.5f}   matches debench over "
          f"{len(r)} scored targets")


if __name__ == "__main__":
    T = Surface189()
    check(T)
    print()
    audit(T)
    print(f"\nprior coverage on the h1-test DE axis: {int(T.prior_covered.sum())}/"
          f"{len(T.genes)} genes;  median n_g {np.median(T.n_prior[T.prior_covered]):.0f}")
    print(f"q_w (|p|-weighted sign accuracy over the pool) {T.qw():.4f};  "
          f"count accuracy {float(T.correct[T.pool].mean()):.4f};  K = {T.mass_ratio():.4f}")

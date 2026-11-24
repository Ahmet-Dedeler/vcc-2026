"""172 — the h1-test evaluation harness for a per-(target, gene) sign-confidence key.

WHAT THIS FILE IS
-----------------
The *instrument*, and nothing else.  It loads the one surface validated against the board for
`reach` on emission-law changes (h1-test, 157 §1 / INBOX-4 §6), aligns a15's DE table with the
design-space effect `H` recovered by 169, and exposes exactly two endpoints:

    reach(key)        -- experiments/153-de-bench/debench.py:Surface.reach, unmodified
    topx(key)         -- sign accuracy over the reference's DE pool in the top X% of the key

Every other script in this experiment imports from here so that no measurement can silently
drift onto a different surface, a different pool, or a different tie convention.

KNOWN-ANSWER ROWS (AGENTS.md §11.1 rule 10) — `--selftest` aborts unless all of these
reproduce:

    reach(production order)                         0.18283   (5 dp)
    reach(|H|)                                      0.23533   (5 dp)
    reach(sig-first, then |H| inside the block)     0.21914   (5 dp)
    reach(oracle: our own signs, perfectly ordered) 0.64553   (5 dp)
    INBOX-4 §1's decile sign-accuracy table for `|H|` and for expression  (3 dp, +/- 0.002)

A screen with no known-answer row cannot detect the class of bug where the wrong numbers look
exactly as reasonable as the right ones (rule 10's corollary; two such bugs were caught in this
project in two days).

THE TOP-X% CONVENTION, AND A CORRECTION TO THE BRIEF I WAS HANDED
-----------------------------------------------------------------
The brief quotes `|H|`'s top-X% sign accuracy as 1% .779 / 2% .729 / 5% .724 / 10% .691 /
25% .646 / 50% .582 / all .536.  That curve does not reproduce here under any convention I could
construct, and its shape is internally inconsistent with INBOX-4 §1's own decile table, which
puts `|H|`'s *first decile* (= top 10%) at **0.621**, not 0.691.

What does reproduce, exactly, is INBOX-4 §1: with the pool taken as the reference's
BH-significant AND adjudicable genes, ordered per target, and accuracy POOLED over
(target, gene) pairs, `|H|` reads 0.621 0.552 0.544 0.537 0.531 0.522 0.521 0.516 0.501 0.499
against the coordinator's 0.621 0.550 0.543 0.536 0.531 0.521 0.521 0.515 0.502 0.499, and
expression reads 0.540 0.551 0.539 0.533 0.535 0.532 0.519 0.537 0.529 0.529 against their
0.540 0.550 0.539 0.532 0.534 0.533 0.519 0.536 0.528 0.528.  Two independent rows, ten cells
each, agreeing to <=0.0021 (the selftest tolerance is 0.003: the coordinator's decile
 edges can differ from mine by one gene in a 89-gene pool, which is worth ~0.001 a cell).  (It also identifies INBOX-4's "`|predicted effect|`" row as `|H|`
and not `|pred_lfc|`, which reads 0.589 in its first decile.)

So this file defines the curve pooled, reports the macro (mean of per-target accuracies) variant
alongside every time, and treats the brief's curve as unreproducible rather than as a target.
The macro variant is the one whose first decile is 0.689 ~ the brief's 0.691, so the brief's
curve is most likely a macro one whose `all` row was taken from the pooled statistic; either way
**every model in this experiment is compared against `|H|` measured the same way, on the same
lines of code**, which is what the comparison needs.

TIES (097 §10).  `topx` never uses ordinal ranks.  A key with a large tied block -- and every
donor-agreement key has one, because the agreement vector is exactly zero wherever a donor is
silent -- would otherwise be handed an arbitrary ramp ordered by position in the gene axis.
`topx` computes the EXPECTED accuracy under a uniformly random tiebreak: the tied block that
straddles the top-k boundary contributes its own mean correctness, weighted by the fraction of
it that fits.  `reach` itself is debench's, which breaks ties by gene name -- deterministic and
identical for every key, so it cannot favour one key over another; `topx` is reported both ways
in the selftest to show the two agree.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, os.path.join(ROOT, "experiments", "153-de-bench"))
from debench import Surface  # noqa: E402

DE_PATH = os.path.join(ROOT, "experiments", "164-callset-artefact", "de164.a15.npz")
EFF_PATH = os.path.join(ROOT, "experiments", "169-sharpening-six-member", "work",
                        "h1test-effect.npz")
TOPX = (0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.00)


class H1Test:
    """a15's h1-test DE table, plus the design-space effect on the same axes."""

    def __init__(self) -> None:
        self.S = Surface(DE_PATH)
        E = np.load(EFF_PATH, allow_pickle=True)
        eg = {str(g): i for i, g in enumerate(E["genes"])}
        et = {str(t): i for i, t in enumerate(E["targets"])}
        missing_g = [g for g in self.S.genes if g not in eg]
        missing_t = [p for p in self.S.perts if p not in et]
        assert not missing_g and not missing_t, (len(missing_g), len(missing_t))
        self.gi = np.array([eg[g] for g in self.S.genes])
        self.pi = np.array([et[p] for p in self.S.perts])
        self.genes = self.S.genes
        self.perts = self.S.perts
        self.H = np.asarray(E["H"])[self.pi][:, self.gi]          # (85, 10779)
        self.b_ctrl = np.asarray(E["b"])[self.gi]                 # (10779,)
        self.cpm = np.asarray(E["cpm"])[self.gi]
        self.gate = np.asarray(E["gate"])[self.gi]
        self.ctrl_totals = np.asarray(E["ctrl_totals"])[self.gi]
        # The population every endpoint reads: the reference's DE pool, adjudicable only.
        self.pool = self.S.real_sig & self.S.adj
        self.correct = (np.sign(self.S.pred_lfc) == np.sign(self.S.real_lfc))

    # ---- endpoints --------------------------------------------------------------------

    def reach(self, key=None) -> float:
        return self.S.reach(key=key)

    def deciles(self, key, macro: bool = False) -> list[float]:
        """Sign accuracy by decile of the ordering key over the reference DE pool.

        d1 is the front of the reach walk.  This is INBOX-4 §1's statistic and is the
        known-answer row that validates the pool, the ordering and the tie convention.
        """
        num = np.zeros(10)
        den = np.zeros(10)
        accs: list[list[float]] = [[] for _ in range(10)]
        for p in range(len(self.perts)):
            idx = np.flatnonzero(self.pool[p])
            if idx.size == 0:
                continue
            k = np.asarray(key[p], dtype=np.float64)[idx]
            ok = self.correct[p, idx].astype(np.float64)
            ok = ok[np.lexsort((self.genes[idx], -k))]
            e = (np.arange(11) * idx.size) // 10
            for d in range(10):
                a, b = e[d], e[d + 1]
                if b <= a:
                    continue
                num[d] += ok[a:b].sum()
                den[d] += b - a
                accs[d].append(float(ok[a:b].mean()))
        return [float(np.mean(accs[d])) if macro else float(num[d] / max(den[d], 1))
                for d in range(10)]

    def topx(self, key, fracs=TOPX, tiebreak: str = "expected") -> list[float]:
        """Sign accuracy over the reference DE pool, pooled over (target, gene) pairs.

        `tiebreak='expected'` averages over random tiebreaks (see module docstring);
        `tiebreak='name'` reproduces debench's own gene-name lexsort.
        """
        num = np.zeros(len(fracs))
        den = np.zeros(len(fracs))
        for p in range(len(self.perts)):
            idx = np.flatnonzero(self.pool[p])
            if idx.size == 0:
                continue
            k = np.asarray(key[p], dtype=np.float64)[idx]
            ok = self.correct[p, idx].astype(np.float64)
            if tiebreak == "name":
                order = np.lexsort((self.genes[idx], -k))
                cum = np.cumsum(ok[order])
                for j, f in enumerate(fracs):
                    n = max(1, int(np.ceil(f * idx.size)))
                    num[j] += cum[n - 1]
                    den[j] += n
                continue
            order = np.argsort(-k, kind="stable")
            ks, oks = k[order], ok[order]
            # group boundaries of equal key values
            newgrp = np.ones(ks.size, bool)
            newgrp[1:] = ks[1:] != ks[:-1]
            gid = np.cumsum(newgrp) - 1
            gsz = np.bincount(gid)
            gsum = np.bincount(gid, weights=oks)
            gstart = np.concatenate(([0], np.cumsum(gsz)[:-1]))
            cum_before = np.concatenate(([0.0], np.cumsum(gsum)[:-1]))
            for j, f in enumerate(fracs):
                n = max(1, int(np.ceil(f * idx.size)))
                g = int(np.searchsorted(gstart, n - 1, side="right") - 1)
                take = n - gstart[g]
                num[j] += cum_before[g] + gsum[g] * (take / gsz[g])
                den[j] += n
        return list(num / np.maximum(den, 1))

    def topx_macro(self, key, fracs=TOPX) -> list[float]:
        """Mean over targets of the per-target top-X% accuracy (the macro variant)."""
        accs: list[list[float]] = [[] for _ in fracs]
        for p in range(len(self.perts)):
            idx = np.flatnonzero(self.pool[p])
            if idx.size == 0:
                continue
            k = np.asarray(key[p], dtype=np.float64)[idx]
            ok = self.correct[p, idx].astype(np.float64)
            cum = np.cumsum(ok[np.lexsort((self.genes[idx], -k))])
            for j, f in enumerate(fracs):
                n = max(1, int(np.ceil(f * idx.size)))
                accs[j].append(float(cum[n - 1] / n))
        return [float(np.mean(a)) for a in accs]

    # ---- reference keys ---------------------------------------------------------------

    def key_absH(self) -> np.ndarray:
        return np.abs(self.H)

    def key_sigfirst_absH(self) -> np.ndarray:
        return np.where(self.S.pred_sig, 1e6 + np.abs(self.H), np.abs(self.H))

    def key_oracle(self) -> np.ndarray:
        return self.correct.astype(np.float64)

    def key_random(self, seed: int = 0) -> np.ndarray:
        return np.random.default_rng(seed).random(self.H.shape)


EXPECTED = {
    "reach_production": 0.18283,
    "reach_absH": 0.23533,
    "reach_sigfirst_absH": 0.21914,
    "reach_oracle": 0.64553,
}
# INBOX-4 §1, the coordinator's published decile table (their "|predicted effect|" row is |H|).
EXPECTED_DECILES = {
    "absH": [0.621, 0.550, 0.543, 0.536, 0.531, 0.521, 0.521, 0.515, 0.502, 0.499],
    "expression": [0.540, 0.550, 0.539, 0.532, 0.534, 0.533, 0.519, 0.536, 0.528, 0.528],
}


def check(T: "H1Test", verbose: bool = True) -> None:
    """Every script in this experiment calls this before reporting anything (rule 10)."""
    bad = []
    got = {
        "reach_production": T.reach(),
        "reach_absH": T.reach(T.key_absH()),
        "reach_sigfirst_absH": T.reach(T.key_sigfirst_absH()),
        "reach_oracle": T.reach(T.key_oracle()),
    }
    for k, v in got.items():
        e = EXPECTED[k]
        ok = abs(v - e) < 5e-6
        bad += [] if ok else [k]
        if verbose:
            print(f"  {'OK ' if ok else 'FAIL'} {k:<24} {v:.5f}   expected {e:.5f}")
    keys = {"absH": T.key_absH(),
            "expression": np.broadcast_to(T.b_ctrl, T.H.shape)}
    for name, k in keys.items():
        cur = T.deciles(k)
        exp = EXPECTED_DECILES[name]
        ok = all(abs(a - b) <= 0.003 for a, b in zip(cur, exp))
        bad += [] if ok else [f"deciles_{name}"]
        if verbose:
            print(f"  {'OK ' if ok else 'FAIL'} deciles_{name:<16} " +
                  " ".join(f"{x:.3f}" for x in cur))
            print(f"      INBOX-4 §1              " + " ".join(f"{x:.3f}" for x in exp))
    if bad:
        raise SystemExit(f"KNOWN-ANSWER ROWS FAILED: {bad} -- aborting (rule 10)")


def selftest() -> None:
    T = H1Test()
    n_real = T.S.real_sig.sum(1)
    n_pool = T.pool.sum(1)
    print(f"surface: {len(T.perts)} targets x {len(T.genes)} genes | "
          f"n_real median {np.median(n_real):.0f} (min {n_real.min()} max {n_real.max()}) | "
          f"adjudicable pool median {np.median(n_pool):.0f}")
    check(T)
    for mode in ("expected", "name"):
        print(f"  --  topx_absH[{mode:<8}]  " + " ".join(f"{x:.3f}" for x in
                                                        T.topx(T.key_absH(), tiebreak=mode)))
    print("  --  topx_absH[macro   ]  " + " ".join(
        f"{x:.3f}" for x in T.topx_macro(T.key_absH())))
    print("  --  deciles_absH[macro]  " + " ".join(
        f"{x:.3f}" for x in T.deciles(T.key_absH(), macro=True)))
    print(f"  --  reach(random seed 0)     {T.reach(T.key_random()):.5f}   "
          f"(INBOX-4 §2 reads 0.16282 for its own random key)")
    print("all known-answer rows reproduce; the instrument is trustworthy on this surface")


if __name__ == "__main__":
    selftest()

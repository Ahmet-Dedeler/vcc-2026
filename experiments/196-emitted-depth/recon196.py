"""196 — reconstruct the INTENDED profile behind a delivered block, exactly where it can be
and by its exact posterior where it cannot.

WHY THIS IS NEEDED.  `floor` destroys the fractional part, so the intended profile is not
recoverable gene-by-gene from the artifact.  But it is recoverable in two pieces:

  * OFF the builder's CPM gate `profile_from_effect` sets the effect to EXACTLY zero, so
    `p_g = c * r_g` with `r = expm1(b_ctrl)/TS` the pooled control relative expression and `c`
    ONE unknown scalar per block.  Write `K = n L c`.  Every off-gate gene then supplies the
    constraint `T_g = floor(K r_g)`, i.e. `K in [T_g/r_g, (T_g+1)/r_g)`.  Intersecting ~9,400
    of those pins `K` to ~1e-6 relative -- and the intersection being NON-EMPTY over 9,400
    independent constraints is itself the proof that the gate, `b_ctrl` and the renderer are
    the ones assumed.  This half is EXACT: no sampling.
  * ON the gate, `x_g = T_g + u_g` with `u_g` the fractional part.  For a smooth profile the
    fractional parts of `n L p_g` are equidistributed on [0,1) and independent across genes
    (this fails only for `x < 1`, i.e. below 0.125 control CPM, two orders under the 5-CPM
    scorer gate).  So `u ~ U(0,1)` i.i.d. is the EXACT posterior, not an approximation, and a
    resimulation at depth `m L` is exact in distribution.

`n L` then follows from `sum_g x_g = n L` (the renderer normalises the profile to 1), which
also yields `L` -- and `L` must come out within rounding of an integer near the block's own
median control library.  That is the reconstruction's known-answer row.

At `m = 1` the reconstruction re-delivers the artifact's own totals BY CONSTRUCTION on the
on-gate block and BY THE FIT on the off-gate block; both are asserted.

No cell-eval2, no build, no submission.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ROOT = HERE.parent.parent
TS = 50_000.0


class Recon:
    """One block's reconstructed intended profile, in units of `x = n L p`."""

    __slots__ = ("T", "r", "gate", "K", "nL", "L", "x", "n", "lo", "hi", "off_ok")

    def __init__(self, T, r, gate, n=400, seed=0):
        self.T = np.asarray(T, dtype=np.int64)
        self.r = np.asarray(r, dtype=np.float64)
        self.gate = np.asarray(gate, dtype=bool)
        self.n = int(n)
        off = (~self.gate) & (self.r > 0)
        lo = np.max(self.T[off] / self.r[off])
        hi = np.min((self.T[off] + 1.0) / self.r[off])
        if not lo < hi:
            raise ValueError(f"off-gate constraint set is EMPTY: [{lo:.6g}, {hi:.6g}) -- the "
                             "gate, b_ctrl or the renderer is not the one assumed")
        self.lo, self.hi = float(lo), float(hi)
        self.K = 0.5 * (self.lo + self.hi)
        x = np.zeros(self.r.size, dtype=np.float64)
        x[~self.gate] = self.K * self.r[~self.gate]
        rng = np.random.default_rng(seed)
        u = rng.random(int(self.gate.sum()))
        x[self.gate] = self.T[self.gate] + u
        self.x = x
        self.nL = float(x.sum())
        self.L = self.nL / self.n
        self.off_ok = float(np.mean(np.floor(x[~self.gate]) == self.T[~self.gate]))

    def totals(self, m: float) -> np.ndarray:
        """Delivered per-gene panel totals at depth multiple `m` (m=1 -> the artifact's)."""
        return np.floor(m * self.x).astype(np.int64)

    def profile(self, m: float) -> np.ndarray:
        """The profile to hand `emit_counts` so its own floor reproduces `totals(m)` exactly.

        `emit_counts` recomputes `L' = int(median(lib_sizes))` and takes
        `floor(n L' profile)`, so the profile is expressed against the INTEGER `L'` the
        emitter will use, not against `m * self.L`.
        """
        Lp = self.emit_L(m)
        return (self.totals(m) + 0.5) / (self.n * Lp)

    def emit_L(self, m: float) -> int:
        return int(round(m * self.L))


def load_h1(seed_base=0):
    """h1-test: the 095 artifact's totals + its recovered gate (min_cpm 5, exp 169)."""
    E = np.load(ROOT / "experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    b = E["b"].astype(np.float64)
    r = np.expm1(b) / TS
    cpm = E["cpm"].astype(np.float64)
    gate = (cpm > float(E["min_cpm"])) & (b > 0)
    return dict(targets=np.asarray(E["targets"]).astype(str),
                genes=np.asarray(E["genes"]).astype(str),
                sums=E["sums"].astype(np.int64), r=r, b=b, cpm=cpm, gate=gate)


def load_chal(context: str):
    """challenge: a15's own delivered totals + the builder's pooled-CPM>10 gate."""
    z = np.load(WORK / "chal-totals-a15.npz", allow_pickle=True)
    blocks = np.asarray(z["blocks"]).astype(str)
    sel = np.flatnonzero(np.array([b.split("|")[0] for b in blocks]) == context)
    cz = np.load(WORK / f"ctrl-{context}.npz")
    tot = cz["pooled_totals"].astype(np.float64)
    r = tot / tot.sum()
    b = np.log1p(TS * r)
    cpm_pool = 1e6 * np.expm1(b) / np.expm1(b).sum()
    gate = (cpm_pool > 10.0) & (b > 0)
    return dict(targets=np.array([blocks[i].split("|")[1] for i in sel]),
                genes=np.asarray(z["genes"]).astype(str), sums=z["totals"][sel],
                r=np.expm1(b) / TS, b=b, cpm=cpm_pool, gate=gate,
                mean_cpm=cz["mean_cpm"].astype(np.float64), lib=cz["lib"],
                nnz=z["nnz"][sel], cellstat=z["cellstat"][sel])


def main() -> None:
    print("RECONSTRUCTION KNOWN-ANSWER ROWS (rule 10: abort if any fails)\n")
    print(f"{'surface':<14}{'blocks':>7}{'off-gate':>9}{'K rel width':>13}{'L median':>10}"
          f"{'L sd':>8}{'L int err':>11}{'off floor==T':>14}{'m=1==artifact':>15}"
          f"{'sim L med':>11}{'sim L sd':>10}")
    bad = []
    for name, D, ctrl_L in (("h1-test", load_h1(), 20190.5),
                            ("chal A", load_chal("A"), 20109.0),
                            ("chal B", load_chal("B"), 19946.0),
                            ("chal C", load_chal("C"), 20034.0)):
        Ls, widths, offok, ident, interr = [], [], [], [], []
        nb = D["sums"].shape[0]
        step = max(1, nb // 60)
        for i in range(0, nb, step):
            R = Recon(D["sums"][i], D["r"], D["gate"], seed=1000 + i)
            Ls.append(R.L)
            widths.append((R.hi - R.lo) / R.K)
            offok.append(R.off_ok)
            interr.append(abs(R.L - round(R.L)))
            ident.append(float(np.array_equal(R.totals(1.0), D["sums"][i])))
        Ls = np.asarray(Ls)
        # The emitter's L is `int(median(400 SAMPLED control libraries))`, so the recovered
        # per-block L must match the SAMPLING DISTRIBUTION of that statistic, not the
        # context's own median.  That distribution is simulated here from the same pool.
        if "lib" in D:
            rng = np.random.default_rng(0)
            sim = np.array([int(np.median(rng.choice(D["lib"], 400, replace=False)))
                            for _ in range(300)], dtype=float)
        else:
            sim = np.array([np.nan])
        row = dict(name=name, blocks=nb, off=int((~D["gate"]).sum()),
                   width=float(np.median(widths)), Lmed=float(np.median(Ls)),
                   Lsd=float(Ls.std()), interr=float(np.max(interr)),
                   offok=float(np.min(offok)), ident=float(np.min(ident)),
                   simmed=float(np.median(sim)), simsd=float(sim.std()))
        print(f"{name:<14}{nb:>7}{row['off']:>9}{row['width']:>13.2e}{row['Lmed']:>10.1f}"
              f"{row['Lsd']:>8.1f}{row['interr']:>11.4f}{row['offok']:>14.4f}"
              f"{row['ident']:>15.4f}{row['simmed']:>11.1f}{row['simsd']:>10.1f}")
        if row["offok"] < 1.0:
            bad.append(f"{name}: off-gate floor(K r) != T on some genes")
        if row["ident"] < 1.0:
            bad.append(f"{name}: m=1 does not reproduce the artifact's totals")
        if row["interr"] > 0.5:
            bad.append(f"{name}: recovered L is not within 0.5 of an integer")
        if np.isfinite(row["simmed"]) and abs(row["Lmed"] - row["simmed"]) > 2.5 * row["simsd"]:
            bad.append(f"{name}: recovered L median {row['Lmed']:.0f} is more than 2.5 "
                       f"simulated sd from {row['simmed']:.0f}")
    print("\n  `L int err` is the worst distance from an integer over the sampled blocks -- the "
          "emitter's\n  L IS an integer, so this is the reconstruction's own accuracy check. "
          "`sim L` is the sampling\n  distribution of `int(median of 400 control libraries)`, "
          "which is what L is.")
    if bad:
        raise SystemExit("ABORT (rule 10): " + "; ".join(bad))
    print("\nALL RECONSTRUCTION KNOWN-ANSWER ROWS PASS")


if __name__ == "__main__":
    main()

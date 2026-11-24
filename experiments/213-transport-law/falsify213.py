"""213 Phase 4 — the falsifiers that do not fall out of the sweep itself.

F2  THE PERMUTATION CONTROL.  A transport law is, by construction, a map that reads the
    RECIPIENT's own basal at each gene.  So derange it: apply a fixed random permutation to
    the recipient's control profile before rendering, keeping the CPM gate on the TRUE basal
    so that exactly the same genes carry effect.  Every law is then transporting a donor
    effect onto a cell that has the right marginal expression distribution and the wrong
    gene-to-expression pairing.

    The prediction is sharp and it is different for the incumbent than for anything else.
    Under `log(p + c0)` the delivered effect is `log1p(TS p') - b = g` up to the mass
    renormalisation, INDEPENDENT of `b`, so the incumbent should barely move.  Under any
    other `c` the recipient's basal enters the delivered effect materially, so a deranged
    recipient should cost that law real cosine.  **Any law whose gain over the incumbent
    SURVIVES the derangement is not measuring transport** (209's random control was the
    maximiser over all 36 of its arms; that is the outcome to watch for).

F2b THE DONOR-SIDE DERANGEMENT.  Same permutation applied to the DONOR's basal inside
    `to_law`, which is the other half of the law.

F5  THE TWO FAMILIES' ADDITIVE LIMITS.  `rawadd` is written directly as
    `p' = p_rec + (p^pert_don - p^ctrl_don)`; `bc1` reaches the same map through the Box-Cox
    parametrisation.  They must agree BIT-FOR-BIT, and they are computed by different code.

F6  THE PERTURBMATCH IMPURITY.  PerturbMatch's table is a limma-voom log-fold-change and has
    no law-space form, so it is left untransformed in every arm (which is also what keeps the
    identity row exact).  `plolo213.py --no-pm` removes it entirely; the verdict must not
    depend on it.

    uv run python falsify213.py --laws logc20 logc40 logc5 mult rawadd
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
import score204  # noqa: E402
score204.WORK = plolo.WORK.parent

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   endpoints, solve_fraction)
from score204 import PROD, _pds_per_target  # noqa: E402
import law213 as L  # noqa: E402
import prep213 as P  # noqa: E402
import plolo213 as S  # noqa: E402


class Deranged:
    """A law whose DONOR-side basal is gene-permuted.  Wraps any law, changes nothing else."""

    def __init__(self, law, perm, seed=0):
        self.law = law
        self.seed = seed
        self._perm: dict[int, np.ndarray] = {}
        self.key = f"{law.key}_dondrg{seed}"
        self.kind = law.kind
        self.param = law.param
        self.note = law.note + " [donor basal deranged]"

    @property
    def is_incumbent(self):
        return False           # force the general path: a deranged basal is not the shipped map

    def asdict(self):
        return dict(key=self.key, kind=self.kind, param=float(self.param), note=self.note,
                    is_incumbent=False)

    def _p(self, n):
        # the donor's basal is deranged on whatever column axis it is handed (the full
        # mapped axis for the streamed statistics, the unique axis for the target rows),
        # so the permutation is cached per size from one seed
        if n not in self._perm:
            self._perm[n] = np.random.default_rng(1000 + self.seed + n).permutation(n)
        return self._perm[n]

    def to_law(self, h, p_ctrl):
        p = np.asarray(p_ctrl)
        return self.law.to_law(h, p[self._p(p.size)])

    def render(self, g, p_ctrl, b_ctrl):
        return self.law.render(g, p_ctrl, b_ctrl)

    def dp_dg(self, p_pert, p_ctrl, g):
        return self.law.dp_dg(p_pert, p_ctrl, g)


def score_recipient_deranged(fold, pre, op, cd4, truth_cols, truth, law, perm_seed):
    """`plolo213.score_law`, a15 geometry, with the RECIPIENT basal gene-permuted.

    The gate stays on the TRUE basal so the arm and its control carry effect on exactly the
    same genes; only the basal the LAW reads is deranged.
    """
    eff = pre * PROD["alpha"]
    whitened = op.apply(eff, fold["wcols"], PROD["gamma"], PROD["lam"])
    frac, _ = solve_fraction(whitened, cd4, PROD["cd4_f"])
    final, _ = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    acc, perc = [], []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        rng = np.random.default_rng(perm_seed) if perm_seed is not None else None
        bb = b if rng is None else b[rng.permutation(b.size)]
        # a deranged basal is no longer the shipped map, so force the general render path
        lw = law if perm_seed is None else _ForceGeneral(law)
        dl = L.deliver_law(final, bb, gate, lw, None)
        acc.append(endpoints(dl[:, truth_cols], truth))
        _, c = _pds_per_target(dl[:, truth_cols], truth)
        perc.append(c)
    out = {f"a15_{k}": float(np.mean([a[k] for a in acc])) for k in acc[0]}
    out["_cos_per_target_a15"] = np.mean(perc, axis=0).tolist()
    return out


class _ForceGeneral:
    """The same law with `is_incumbent` False, so `deliver_law` uses the general path."""

    def __init__(self, law):
        self._l = law
        for a in ("key", "kind", "param", "note"):
            setattr(self, a, getattr(law, a))

    is_incumbent = False

    def to_law(self, h, p):
        return self._l.to_law(h, p)

    def render(self, g, p, b):
        c = getattr(self._l, "param", None)
        if self._l.kind == "logc":
            return (p + c) * np.exp(g) - c
        if self._l.kind == "boxcox" and c == 0.0:
            return (p + L.C0) * np.exp(g) - L.C0
        return self._l.render(g, p, b)

    def dp_dg(self, pp, p, g):
        return self._l.dp_dg(pp, p, g)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*",
                    default=["logc20", "logc40", "logc5", "logc100", "mult", "rawadd"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--seeds", nargs="*", type=int, default=[11, 22, 33])
    args = ap.parse_args()

    laws = L.all_laws()
    out_path = WORK / "falsify213.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["law"], r["control"], r["seed"]) for r in rows}

    for held in args.folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        op_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                       else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        gene_n = len(plolo.OFFICIAL)
        for key in args.laws:
            base = laws[key]
            srcs = S.law_sources(held, fold, base)
            pre, _ = score204.weighted_baseline(srcs, fold["targets"], "equal",
                                                {n: 1.0 for n, _ in srcs})
            # (a) recipient-basal derangement, plus the un-deranged reference through the
            #     SAME general render path so the two differ only in the permutation
            for seed in [None] + list(args.seeds):
                tagc = "recipient_true" if seed is None else "recipient_deranged"
                if (held, key, tagc, seed) in done:
                    continue
                t0 = time.time()
                e = score_recipient_deranged(fold, pre, op, cd4, fold["scored"], truth,
                                             base, seed)
                rows.append(dict(fold=held, law=key, control=tagc, seed=seed, **e,
                                 seconds=round(time.time() - t0, 1)))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/{key}/{tagc}/{seed}] a15cos {e['a15_cos']:.7f} "
                      f"pds {e['a15_pds']:.6f}", flush=True)
            del srcs, pre
            gc.collect()
            # (b) donor-basal derangement
            for seed in args.seeds[:1]:
                if (held, key, "donor_deranged", seed) in done:
                    continue
                t0 = time.time()
                # The donor derangement acts on the DONOR's own column axis inside
                # `to_law`, so the sources are rebuilt here rather than reused.  The
                # PerturbMatch slot is left alone (it carries no law at all), so this
                # control is conservative: it deranges only the half that the law reads.
                srcs = []
                for n in FOLDS[held]["sources"]:
                    blk = S.meta(n)
                    dlaw = Deranged(base, None, seed)
                    sp = P.prepared(n, dlaw, fold["targets"], fold["panel_positions"], blk)
                    srcs.append((n, sp))
                pre2, _ = score204.weighted_baseline(srcs, fold["targets"], "equal",
                                                     {n: 1.0 for n, _ in srcs})
                e = score_recipient_deranged(fold, pre2, op, cd4, fold["scored"], truth,
                                             base, None)
                rows.append(dict(fold=held, law=key, control="donor_deranged", seed=seed,
                                 **e, seconds=round(time.time() - t0, 1)))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/{key}/donor_deranged/{seed}] a15cos {e['a15_cos']:.7f} "
                      f"pds {e['a15_pds']:.6f}", flush=True)
                del srcs, pre2
                gc.collect()
        del op, cd4
        gc.collect()


if __name__ == "__main__":
    main()

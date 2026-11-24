"""153 — DE-bench: deterministic, seed-free evaluation of the three DE members.

WHY THIS EXISTS
---------------
`reach`, `fid` and `jac` have only ever been measured in this repo by emitting 400 count
vectors and running the scorer.  That path carries seed noise an order of magnitude larger
than the effects being chased: 095 §12 measured two gamma-Poisson draws at reach 0.1675 vs
0.1535, and `realize_confirm.log` has a same-policy re-seed pair 0.01222 apart.  Effects of
+0.014 are not resolvable through it.

But the metric definitions do not require the emitter:

  * `nmae`  gate is reference-only and reads our lfc -> a deterministic functional of the
            mean profile.
  * `reach` is a pure ordering statistic over the reference's DE pool.  Given an ordering
            key it is exact.
  * `fid`   and `jac` are exact given a call set.

So all three can be evaluated as closed-form functions of (predicted lfc, call set, ordering
key) against a cached reference DE table, in seconds, with ZERO seed noise.

Definitions follow docs/reference/vcc2026-metrics-BRIEF.md sections 4-7 literally.

VALIDATE THE VALIDATOR (repo rule 3): `python debench.py --selftest` reproduces the official
cell-eval2 numbers for the production prediction on orion-hct116-thinned, which 095 §11
records as reach 0.0770 / fid 0.5408 / jac 0.0120.
"""
from __future__ import annotations
import numpy as np

ALPHA = 0.05


class Surface:
    """A cached reference+prediction DE table for one harness."""

    def __init__(self, path: str, n_real_target: int | None = None):
        """n_real_target -- POWER MATCHING (peer correction, 094 A8-BOARD-RESULT.md:70).

        `fid`, `jac` and `reach`'s normalisation all carry n_real in a denominator or a
        union, and n_real is 68-322 on every harness we hold against ~2,000 inferred on the
        challenge.  So every call-set cut looks free here and is charged there; that artefact
        has already cost three board predictions (v14 fid, v15 fid, a8 jac).  Removing the
        emitter kills the SEED noise but not this REGIME error -- a deterministic fid would
        reproduce the artefact without error bars, which is worse because it looks
        authoritative.

        The challenge reference has more power (deeper data, more cells), not different
        biology, so we emulate it by taking the reference's top-N genes by its own p-value
        instead of its BH-significant set.  Everything downstream is unchanged.
        """
        z = np.load(path, allow_pickle=True)
        self._n_real_target = n_real_target
        self.genes = np.asarray([str(g) for g in z["genes"]])
        self.perts = np.asarray([str(p) for p in z["perts"]])
        self.real_padj, self.real_lfc = z["real_padj"], z["real_lfc"]
        # A reference-only table (de_tables.py with no --pred) has no prediction side; callers
        # that supply their own pred_lfc do not need one.
        nan = np.full(self.real_lfc.shape, np.nan)
        self.pred_padj = z["pred_padj"] if "pred_padj" in z.files else nan.copy()
        self.pred_p = z["pred_p"] if "pred_p" in z.files else nan.copy()
        self.pred_lfc = z["pred_lfc"] if "pred_lfc" in z.files else np.zeros_like(self.real_lfc)
        # Each row's own target is removed from every member (BRIEF 4-7).
        gpos = {g: i for i, g in enumerate(self.genes)}
        self.self_idx = np.array([gpos.get(p, -1) for p in self.perts])
        # "Adjudicable": the reference assigns a defined, non-zero log2 fold change.
        self.adj = np.isfinite(self.real_lfc) & (self.real_lfc != 0.0)
        if n_real_target is None:
            self.real_sig = (self.real_padj < ALPHA) & np.isfinite(self.real_padj)
        else:
            rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
            self.real_sig = np.zeros(rp.shape, bool)
            k = min(n_real_target, rp.shape[1])
            for i in range(rp.shape[0]):
                self.real_sig[i, np.argpartition(rp[i], k - 1)[:k]] = True
            self.real_sig &= np.isfinite(z["real_p"])
        self.pred_sig = (self.pred_padj < ALPHA) & np.isfinite(self.pred_padj)
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.adj[p, j] = self.real_sig[p, j] = self.pred_sig[p, j] = False

    # ---- the three members ------------------------------------------------------------

    def reach(self, key=None, call=None) -> float:
        """BRIEF 5.  Pool = reference-significant genes.  Walk in our confidence order,
        counting only adjudicable genes; k* is the DEEPEST prefix with >=90% sign purity
        ("deepest" is literal -- purity is not monotone, so a prefix that dips and recovers
        counts at its deeper crossing).  R_p = k*/n_real.

        key  : higher = more confident.  None -> the production order, i.e. our significant
               calls first, then ascending pred p_adj, then ascending pred p, then
               descending |lfc|, then gene name.
        call : our significant set.  None -> the prediction's own.
        """
        out = []
        for p in range(len(self.perts)):
            pool = np.flatnonzero(self.real_sig[p])
            if pool.size == 0:
                continue
            order = self._order(p, pool, key, call)
            pool = pool[order]
            pool = pool[self.adj[p, pool]]          # count only adjudicable
            if pool.size == 0:
                out.append(0.0)
                continue
            ok = np.sign(self.pred_lfc[p, pool]) == np.sign(self.real_lfc[p, pool])
            purity = np.cumsum(ok) / np.arange(1, pool.size + 1)
            deep = np.flatnonzero(purity >= 0.9)
            k = int(deep[-1] + 1) if deep.size else 0
            out.append(k / int(self.real_sig[p].sum()))
        return float(np.mean(out))

    def fid(self, call=None) -> float:
        """BRIEF 4.  F_p = k / max(n_pred, n_real); the call set is intersected with the
        adjudicable genes, and a called gene with no direction counts as a miss."""
        out = []
        for p in range(len(self.perts)):
            n_real = int(self.real_sig[p].sum())
            sel = np.flatnonzero((self.pred_sig[p] if call is None else call[p]) & self.adj[p])
            n_pred = sel.size
            if n_pred == 0 and n_real == 0:
                continue
            k = int((np.sign(self.pred_lfc[p, sel]) == np.sign(self.real_lfc[p, sel])).sum())
            out.append(k / max(n_pred, n_real) if max(n_pred, n_real) else 0.0)
        return float(np.mean(out))

    def jac(self, call=None) -> float:
        """BRIEF 6.  Jaccard of the two significant sets; an empty union is defined as 1."""
        out = []
        for p in range(len(self.perts)):
            a = self.real_sig[p]
            b = self.pred_sig[p] if call is None else call[p]
            u = int((a | b).sum())
            out.append(1.0 if u == 0 else int((a & b).sum()) / u)
        return float(np.mean(out))

    def nmae(self, pred_lfc=None) -> float:
        """BRIEF 7.  Gate is the reference's own significant set with finite reference lfc.
        The context value is the MEAN of the per-perturbation ratios (not the ratio of sums
        -- that is `mse`'s rule, and using it here reads 0.0078 high).  A perturbation is
        scored only if its gate holds at least 10 genes; a gene the submission does not
        carry, or carries non-finite, is read as a predicted fold change of zero.

        Deterministic in the mean profile -- no emitter, no seed.  This is the member with
        the second-largest gap to the leader (+0.0896 scaled) and it has never been used as
        a screening endpoint in this repo."""
        P = self.pred_lfc if pred_lfc is None else pred_lfc
        out = []
        for p in range(len(self.perts)):
            s = np.flatnonzero(self.real_sig[p] & np.isfinite(self.real_lfc[p]))
            if s.size < 10:
                continue
            den = np.abs(self.real_lfc[p, s]).sum()
            if not np.isfinite(den) or den == 0:
                continue
            h = np.where(np.isfinite(P[p, s]), P[p, s], 0.0)
            out.append(np.abs(h - self.real_lfc[p, s]).sum() / den)
        return float(np.mean(out))


    def pds(self, pred_lfc=None, exclude_panel=True) -> float:
        """`pds_cosine` (BRIEF 1): for each perturbation, the rank of its OWN true effect among
        all perturbations by cosine distance, mapped to [0,1].  The reference's real_lfc IS the
        truth, so this is exact on this surface -- which lets the pds cost and the reach gain of
        a transform be priced on the SAME data.  088 measured them on different surfaces at
        different alphas; that is the gap this closes.

        Scale-invariant by construction, so alpha never confounds it.
        """
        P = self.pred_lfc if pred_lfc is None else pred_lfc
        g = np.flatnonzero(np.isfinite(self.real_lfc).all(0) & np.isfinite(P).all(0))
        if exclude_panel:                      # pds_cosine removes ALL panel targets
            sidx = set(int(i) for i in self.self_idx if i >= 0)
            g = np.array([i for i in g if i not in sidx])
        A, B = P[:, g].astype(np.float64), self.real_lfc[:, g].astype(np.float64)
        A /= np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
        B /= np.maximum(np.linalg.norm(B, axis=1, keepdims=True), 1e-12)
        S = A @ B.T                            # cosine of every prediction to every truth
        n = S.shape[0]
        own = np.diag(S)
        # fraction of decoys strictly farther than the true match, ties at half
        beaten = (S < own[:, None]).sum(1) - 0
        ties = (S == own[:, None]).sum(1) - 1
        return float(np.mean((beaten + 0.5 * ties) / (n - 1)))

    # ---- ordering ---------------------------------------------------------------------

    def _order(self, p, pool, key, call):
        if key is None:
            sig = (self.pred_sig[p] if call is None else call[p])[pool]
            return np.lexsort((self.genes[pool], -np.abs(self.pred_lfc[p, pool]),
                               self.pred_p[p, pool], self.pred_padj[p, pool], ~sig))
        k = key[p][pool] if key[p].ndim else key[p]
        return np.lexsort((self.genes[pool], -np.asarray(k, dtype=float)))

    # ---- diagnostics ------------------------------------------------------------------

    def sign_accuracy(self, mask=None) -> float:
        """Fraction of correct signs over reference-significant, adjudicable genes."""
        num = den = 0
        for p in range(len(self.perts)):
            s = self.real_sig[p] & self.adj[p]
            if mask is not None:
                s = s & mask[p]
            i = np.flatnonzero(s)
            num += int((np.sign(self.pred_lfc[p, i]) == np.sign(self.real_lfc[p, i])).sum())
            den += i.size
        return num / max(den, 1)

    def counts(self):
        nr = np.array([int(self.real_sig[p].sum()) for p in range(len(self.perts))])
        npd = np.array([int((self.pred_sig[p] & self.adj[p]).sum()) for p in range(len(self.perts))])
        return nr, npd


# ---- scaled-member conversion (anchors fitted on 780 board entries, 094 §1) -------------
ANCHOR = {"pds": (0.50020, 0.95228), "fid": (0.51239, 0.80964),
          "reach": (0.07902, 0.96414), "jac": (0.03038, 0.39727)}


def scaled(member: str, raw: float) -> float:
    b, r = ANCHOR[member]
    return (raw - b) / (r - b)


def overall_delta(d_fid=0.0, d_reach=0.0, d_jac=0.0) -> float:
    """Raw deltas -> overall.  Safe for fid/reach/jac/pds only: 094 warns that score_mse
    and score_nmae are per-perturbation normalised and are NOT functions of the reported
    aggregate raw values, so never do this arithmetic on those two."""
    return (d_fid / (ANCHOR["fid"][1] - ANCHOR["fid"][0])
            + d_reach / (ANCHOR["reach"][1] - ANCHOR["reach"][0])
            + d_jac / (ANCHOR["jac"][1] - ANCHOR["jac"][0])) / 6.0


if __name__ == "__main__":
    import sys
    s = Surface(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-")
                else "../095-clean-harness/de_thinned.npz")
    nr, npd = s.counts()
    print(f"surface: {s.perts.size} perturbations x {s.genes.size} gated genes")
    print(f"n_real  median {np.median(nr):8.1f}   mean {nr.mean():8.1f}")
    print(f"n_pred  median {np.median(npd):8.1f}   mean {npd.mean():8.1f}")
    print()
    print(f"{'member':8} {'DE-bench':>10} {'cell-eval2 (095 §11)':>22}")
    print(f"{'reach':8} {s.reach():10.5f} {0.0770:22.4f}")
    print(f"{'fid':8} {s.fid():10.5f} {0.5408:22.4f}")
    print(f"{'jac':8} {s.jac():10.5f} {0.0120:22.4f}")
    print(f"{'nmae':8} {s.nmae():10.5f} {0.9309:22.4f}")
    print()
    print(f"sign accuracy on the reference's DE set: {s.sign_accuracy():.4f}")

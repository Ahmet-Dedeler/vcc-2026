"""Shrink toward the PANEL MEAN or toward ZERO? `mech_miller.py` only measured the first.

`mech_miller.py` measured the Miller/James-Stein family

    dh' = a (*) dh + (1 - a) (*) dbar          dbar = panel mean of the predicted deltas

and reported `sqrt(W/W_glob)`, `W = X'^2/A'`, which is what the master formula
`MSE_min = 1 + (Kc - 2*sqrt(claim*W))/sumD` scales by. It never measured the other obvious
target, `dh' = a (*) dh`, which shrinks toward the CONTROL. That matters because the two
anchors are priced differently on the board: team "Kun Lab"'s explicit zero-effect
submission scored raw 1.0424, while the scored baseline -- the context's mean perturbation
response assigned identically to every perturbation -- sits at raw ~0.9879. Predicting the
control is measurably WORSE than predicting the mean response, so a shrinkage that walks
low-confidence genes toward zero should lose to one that walks them toward `dbar`.

This script measures both, on both harnesses, at three gammas, with TWO agreement keys:

  donor_agreement  `MultiDonorModel.donor_agreement`, the key the builder actually ships
                   under `--keep-key agreement` / `--shrink-key agreement`. Built here from
                   the two donors that are NOT the harness's own line, so no leak.
  raw-delta        `mech_miller.py`'s own `donor_key`, kept so its published 1.381 / 1.172
                   can be reproduced. It differs from the shipped key in one way that
                   matters: `donor_agreement` scores the gene-specific DEVIATIONS
                   `s = delta - mean_p(delta)`, not the raw `delta`.

Alpha comes from the SHIPPED `agreement_shrink_vector`, so a null here is a null on the code
that would be built, not on a re-derivation of it (rule 3).

`claim`, `Kc` and `sumD` are held at their unshrunk values, exactly as `mech_miller.py` did:
shrinking moves `C_hat_p` and therefore `claim` a little, and neither script models that.

RESULT (2026-08-28), RE-MEASURED after the blank-row bug was fixed, and it goes AGAINST the
panel-mean prior. Shrinking toward ZERO beats shrinking toward the panel MEAN on both
harnesses, at every gamma, with both keys -- 12 of 12 comparisons -- once `B` is recomputed
instead of assumed. At gamma=1.0 with the shipped key, d(mse) against the unshrunk prediction
rescaled onto its own optimum:

    h1-test        toward zero  -0.00856      toward mean  +0.00336   (mean is WORSE than
    orion-thinned  toward zero  -0.05625      toward mean  -0.01181     not shrinking at all)

THE BUG THIS TABLE WAS RE-RUN FOR, AND WHY THIS TABLE NEVER CAUGHT IT. `agreement_shrink_vector`
used to rank descending with ties broken by GENE INDEX, so a row whose key is all zero -- fewer
than two donors carry the target -- got `exp(-gamma * g / 2000)` indexed by position in the gene
axis, and was then scored as though that had helped. On the real 300-target challenge panel that
was 111 of 900 rows and it inflated the family's estimated gain from +0.0055 to +0.0102. On
THESE harnesses it touched 0 of 132 orion rows and 4 of 85 h1-test rows, because donor coverage
here is 132/132 and 81/85. So this comparison was very nearly immune to the bug and could not
have detected it -- the same rule-10 point that applies to the `claim` omission below. The
re-measure moves every number by less than 0.001 and changes no ordering.

What DID change visibly is the alpha distribution, which is the fingerprint of the second half
of the fix: ties now share a midrank, so the ~16.5k genes whose agreement is exactly 0 all
receive one alpha instead of a ramp across them. At gamma=0.5 the quantiles went from
9.7e-03/3.1e-02/9.9e-02/3.1e-01/1.0 to 5.9e-02/6.1e-02/6.1e-02/6.3e-02/1.0. `W` barely moved
because the ramp's average was already near the midrank value; the point of the fix is that the
old spread was a function of array position, not of the key.

Two things drive the zero-over-mean result, and both are visible in the table.

1. Our own panel mean is not the metric's anchor. The metric's 0-anchor is the REFERENCE's
   mean perturbation response; `cos(our dbar, the reference's mean effect)` is +0.0162 on
   h1-test and -0.0664 on orion, so the capture fraction `cos^2` is 0.0% and 0.4%, and the
   two surfaces do not even share a sign. Shrinking toward our own mean buys none of the
   baseline's advantage. (Independently reproduced by experiment 097's `common_mode.py`.)

2. The master-formula shortcut flatters the mean target specifically. `(A - B)/A` is +0.014 /
   +0.012 unshrunk and stays at or BELOW that under the zero target (+0.0123 / +0.0092, flat
   across gamma), but the mean target drives it to +0.057 / +0.053 at gamma=1 and +0.071 /
   +0.072 at gamma=2 -- because the panel mean is the common mode and `_across_pert_budget` is
   a CENTRED sum of squares, so shrinking toward the mean collapses `B` while leaving `A`.
   `A - B` is exactly what `MSE_min = 1 + (Kc - 2*sqrt(claim*W))/sumD` drops. So `sqrt(W/Wg)`
   still rises for the mean target -- 1.137 / 1.706 at gamma=1 -- while the scorer's own
   identity says it lost.

   THAT FLATNESS IS THE OPERATIONALLY IMPORTANT PART. `eps = (A-B)/B` is priced by `claim/B`,
   which is 0.183 here and 0.369 on orion but 0.879 on the challenge at the #247 cap
   (experiment 097's `exact.py`), so an eps-inflating family is charged 2.4-4.8x more on the
   board than these surfaces report. The zero target does not inflate eps at any gamma tested,
   so that leverage is an argument against the MEAN target, not against shrinking, and not a
   reason to prefer the a11 base over a12.

CONSEQUENCE FOR THE RECORD: MSE-MECHANISM.md intervention 5 quotes `1.17-1.38x` on `sqrt(W)`
and `+0.0044 to +0.0097 overall` for the panel-mean family. That number is computed through
the shortcut, on the one family whose `A - B` the shortcut cannot ignore, and it is therefore
an overstatement. The buildable residue of intervention 5 is real, but it lives in the ZERO
target, which that document never measured.

Reproduced here: the raw-delta key at gamma=1.0 toward mean returns sqrt(W/Wg) = 1.175 on
orion and 1.160 on h1-test, against mech_miller.py's published 1.172 / 1.381. Orion matches;
the h1-test gap is the blank-row handling (4 of 85 targets are carried by a single donor and
get alpha == 1 here rather than an all-tied rank ramp).
"""
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/scripts")

import importlib.util

from cell_eval2.metrics.delta import _across_pert_budget          # noqa: E402
from decompose import terms                                        # noqa: E402
from msekit import group_stats                                     # noqa: E402
from vcc2026.donor_transfer import Donor                           # noqa: E402
from vcc2026.multi_donor import MultiDonorModel                    # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "builder", f"{ROOT}/scripts/make_submission_pm_cd4.py"
)
BUILDER = importlib.util.module_from_spec(_spec)
sys.modules["builder"] = BUILDER
_spec.loader.exec_module(BUILDER)

# The two doc anchor points of intervention 1, from MSE-MECHANISM.md's own table:
# multiplier 1.000 -> raw 0.8785, multiplier 1.172 -> raw 0.8535. `MSE_min` is linear in
# the multiplier, so this reproduces the doc's 1.381 -> 0.8232 row to 1e-4.
RAW_AT_1 = 0.8785
RAW_SLOPE = (0.8535 - 0.8785) / (1.172 - 1.000)
# section (d), fitted on six well-behaved board teams
SCORE_A, SCORE_B = 0.987862879, 0.950162566
KEYS = ("k562_gwps", "orion_hek293t")   # neither is hct116, so neither harness leaks


def board_units(mult: float) -> tuple[float, float, float]:
    raw = RAW_AT_1 + RAW_SLOPE * (mult - 1.0)
    score = float(np.clip((SCORE_A - raw) / SCORE_B, 0.0, 1.0))
    return raw, score, score / 6.0


def load(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    dh, d, w = t["dh"], t["d"], t["w"]
    for k, c in t["cols"].items():
        dh[k, c] = 0.0
        d[k, c] = 0.0
    return t, dh, d, w, rg, [str(x) for x in pp[t["keep"]]]


def W_of(dh, d, w):
    n = float((np.einsum("ij,ij->i", dh, d) * w).sum())
    q = float((np.einsum("ij,ij->i", dh, dh) * w).sum())
    return (n * n / q if q > 0 else 0.0)


def raw_delta_key(genes, targets):
    """`mech_miller.py`'s donor_key, verbatim: raw delta, no donor renormalisation."""
    DON = {}
    for name in KEYS:
        root = f"{ROOT}/data/donors/{name}"
        dg = np.array(open(f"{root}/genes.txt").read().split())
        dp = np.array(open(f"{root}/perts.txt").read().split())
        gi = {g: i for i, g in enumerate(dg)}
        DON[name] = (np.load(f"{root}/delta.npy", mmap_mode="r"),
                     np.array([gi.get(g, -1) for g in genes]),
                     {p: i for i, p in enumerate(dp)})
    out = np.zeros((len(targets), genes.size))
    for r, t in enumerate(targets):
        vs = []
        for n in DON:
            delta, cols, pi = DON[n]
            if t not in pi:
                vs = []
                break
            row = np.asarray(delta[pi[t]], float)
            v = np.zeros(genes.size)
            ok = cols >= 0
            v[ok] = row[cols[ok]]
            vs.append(v)
        if len(vs) == 2:
            a, b = vs
            out[r] = np.sqrt(np.abs(a) * np.abs(b)) * (
                (np.sign(a) == np.sign(b)) & (a != 0) & (b != 0)
            )
    return out


def shipped_key(genes, targets):
    """`MultiDonorModel.donor_agreement`, i.e. the key the builder actually ships."""
    print("  [key] loading donors for MultiDonorModel...", flush=True)
    donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in KEYS}
    model = MultiDonorModel(donors, np.asarray(genes), gene_norm=False, eb=False)
    out = np.vstack([model.donor_agreement(str(t)) for t in targets])
    del donors, model
    return out


WORK = f"{ROOT}/experiments/096-mse-floor/work"


def cached_key(tag, fn, genes, targets):
    path = f"{WORK}/shrink_key_{tag}.npy"
    if os.path.exists(path):
        got = np.load(path)
        if got.shape == (len(targets), genes.size):
            return got
    got = fn(genes, targets)
    np.save(path, got)
    return got


def best_mse(dh_sh, d, w, cols, claim, Kc, sumD):
    """min over the global scale s of the SCORER'S OWN identity, B recomputed not assumed.

        MSE(s) = 1 + (s^2 A' - 2 s X' + Kc - min(1, s^2 B'/claim) * claim) / sumD

    `mech_miller.py` and MSE-MECHANISM.md section (c) both report the master-formula
    shortcut `MSE_min = 1 + (Kc - 2*sqrt(claim*W))/sumD`, which is the s-optimum ONLY when
    A' == B', i.e. when the prediction carries no common mode. Unshrunk that holds to ~1%
    (||dbar|| is 11% of mean ||dh_p||). Shrinking toward the panel mean PRESERVES the common
    mode while collapsing the deviations, so A' - B' stops being negligible exactly in the
    family the shortcut is being used to score -- and `_across_pert_budget` is a CENTRED sum
    of squares, so it is the very quantity the mean target shrinks. Hence the exact route.

    `claim`, `Kc`, `sumD` are held at their unshrunk values, as `mech_miller.py` did.
    """
    A = float((np.einsum("ij,ij->i", dh_sh, dh_sh) * w).sum())
    X = float((np.einsum("ij,ij->i", dh_sh, d) * w).sum())
    B = float(_across_pert_budget(dh_sh, cols, w))
    if A <= 0 or B <= 0:
        return float("nan"), float("nan"), A, B, X
    s_c = np.sqrt(claim / B)                       # the rho cliff

    def mse(s):
        return 1.0 + (s * s * A - 2 * s * X + Kc
                      - min(1.0, s * s * B / claim) * claim) / sumD

    cands = [s_c]
    if A > B:
        cands.append(min(X / (A - B), s_c))        # interior optimum below the cliff
    cands.append(max(X / A, s_c))                  # interior optimum above the cliff
    cands += list(np.geomspace(max(s_c * 1e-3, 1e-9), s_c * 1e3, 4001))
    cands = [c for c in cands if c > 0]
    vals = [mse(c) for c in cands]
    j = int(np.argmin(vals))
    return float(vals[j]), float(cands[j]), A, B, X


HARNESSES = [
    (f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
     f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "h1-test"),
    (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/"
     f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "orion-thinned"),
]

GAMMAS = (0.5, 1.0, 2.0)

for real, pred, name in HARNESSES:
    t, dh, d, w, genes, targets = load(real, pred)
    base = W_of(dh, d, w)
    dbar = dh.mean(0)
    claim, Kc, sumD, cols = t["claim"], t["Kc"], t["sumD"], t["cols"]
    ref_mse, ref_s, refA, refB, refX = best_mse(dh, d, w, cols, claim, Kc, sumD)
    print(f"\n===== {name}: {dh.shape[0]} perts, {dh.shape[1]} genes. "
          f"W_global = {base:.6g} =====")
    print(f"  ||dbar|| = {np.linalg.norm(dbar):.6f}, "
          f"mean ||dh_p|| = {np.linalg.norm(dh, axis=1).mean():.6f}, "
          f"cos(dbar, mean d) = "
          f"{float(dbar @ d.mean(0) / (np.linalg.norm(dbar) * np.linalg.norm(d.mean(0)))):.4f}")
    print(f"  unshrunk, rescaled onto its own optimum: mse={ref_mse:.6f} at s*={ref_s:.4f}; "
          f"A={refA:.4f} B={refB:.4f} (A-B)/A={1 - refB / refA:+.4f} X={refX:.4f}; "
          f"claim={claim:.5f} Kc={Kc:.5f} sumD={sumD:.5f}")
    print(f"  harness mse is UNREADABLE IN LEVEL (rule 1); only the column deltas transfer.")
    for key_name, key_fn in (("donor_agreement (shipped)", shipped_key),
                             ("raw-delta (mech_miller)", raw_delta_key)):
        dk = cached_key(f"{name}_{'shipped' if key_fn is shipped_key else 'rawdelta'}",
                        key_fn, genes, targets)
        cov = int((dk.sum(1) > 0).sum())
        print(f"\n  --- key: {key_name}; covers {cov}/{len(targets)} perturbations ---")
        print(f"  {'family':<30}{'sqrt(W/Wg)':>11}{'A-B)/A':>9}{'s*':>9}"
              f"{'harness mse':>13}{'d(mse)':>10}{'doc raw':>9}{'doc ovr':>9}")
        alphas = None
        for g in GAMMAS:
            # No blank-row special case any more: `agreement_shrink_vector` itself returns
            # all-ones for a key with fewer than two DISTINCT values, and ties share a
            # midrank, so a row that carries no ordering is a no-op and the ~16.5k
            # zero-agreement genes in an ordinary row all get the same alpha instead of a
            # ramp indexed by gene position. Both were bugs; this table is the re-measure.
            a = np.vstack([BUILDER.agreement_shrink_vector(row, g) for row in dk])
            if alphas is None:
                alphas = a
            for tgt, shrunk in (("mean", a * dh + (1 - a) * dbar[None, :]),
                                ("zero", a * dh)):
                m = np.sqrt(W_of(shrunk, d, w) / base)
                raw, sc, ov = board_units(m)
                e_mse, e_s, eA, eB, eX = best_mse(shrunk, d, w, cols, claim, Kc, sumD)
                print(f"  {f'gamma={g}, toward {tgt}':<30}{m:>11.3f}"
                      f"{1 - eB / eA:>+9.4f}{e_s:>9.3f}{e_mse:>13.6f}"
                      f"{e_mse - ref_mse:>+10.6f}{raw:>9.4f}{ov:>+9.4f}")
        # rule 3: the intervention actually happened, as a distribution not a summary
        q = np.quantile(alphas, [0.0, 0.25, 0.5, 0.75, 1.0])
        print(f"  alpha at gamma={GAMMAS[0]}: min/q25/med/q75/max = "
              f"{q[0]:.3e}/{q[1]:.3e}/{q[2]:.3e}/{q[3]:.3e}/{q[4]:.3e}; "
              f"{float((alphas > 0.5).sum(1).mean()):.0f} genes/pert keep alpha>0.5")

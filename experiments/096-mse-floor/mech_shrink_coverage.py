"""Does the donor-agreement key cover the challenge's 300 targets, and what does alpha do?

MSE-MECHANISM.md intervention 5 left this open: "the donor key covered 81/85 and 132/132
perturbations on the harnesses and its coverage of the challenge's 300 targets must be
checked". This checks it, on the REAL production three-donor model and the REAL challenge
gene axis, using the shipped `MultiDonorModel.donor_agreement` and the shipped
`agreement_shrink_vector`. Nothing is emitted; this only reads the donors and the bundle.

Rule 3: it prints the DISTRIBUTION of every intervening quantity, so a null on the arm can
be told apart from the shrinkage never having happened.
"""
import importlib.util
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")

from vcc2026.donor_transfer import Donor                           # noqa: E402
from vcc2026.multi_donor import MultiDonorModel                    # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "builder", f"{ROOT}/scripts/make_submission_pm_cd4.py"
)
B = importlib.util.module_from_spec(_spec)
sys.modules["builder"] = B
_spec.loader.exec_module(B)

spec = B.ChallengeSpec.load(B.Path(f"{ROOT}/data/challenge"))
genes = np.asarray(list(map(str, spec.genes)))
perts = [str(p) for p in spec.perts]
print(f"challenge panel: {len(perts)} targets x {genes.size} genes")

donors = {n: Donor.load(f"{ROOT}/data/donors/{n}") for n in B.PRODUCTION_DONORS}
model = MultiDonorModel(donors, genes, gene_norm=False, eb=False)
carried = np.asarray([
    sum(p in model.pert_index[n] for n in model.names) for p in perts
])
agreement = np.vstack([model.donor_agreement(p) for p in perts])
blank = ~agreement.any(axis=1)
print(f"donors carrying each target: "
      f"{np.bincount(carried, minlength=4)[:4].tolist()} for 0/1/2/3 donors")
print(f"agreement key non-blank on {int((~blank).sum())}/{len(perts)} targets; "
      f"{int(blank.sum())} blank (fewer than two donors carry them)")
nz = (agreement > 0).sum(axis=1)
q = np.quantile(nz[~blank], [0.0, 0.25, 0.5, 0.75, 1.0]) if (~blank).any() else np.zeros(5)
print(f"genes with a NON-ZERO agreement score per carried target, "
      f"min/q25/med/q75/max = {q[0]:.0f}/{q[1]:.0f}/{q[2]:.0f}/{q[3]:.0f}/{q[4]:.0f} "
      f"of {genes.size}")

# Where does each target's OWN gene sit in its agreement order? The on-target knockdown is our
# most confident prediction, and a shrinkage that eats it would look fine in every number we
# check. VERIFIED against cell_eval2: the scored consequence is NIL, because all six members
# exclude the coordinate -- `run._ONTARGET_EXCLUDING_FUNCS` covers mse_unbiased_capped,
# de_sig_jaccard and de_lfc_nmae; `metrics.direction._exclude_on_target` covers fidelity and
# reach ("the scored direction metrics remove each RESOLVED target's own gene from N_conf");
# and `metrics.discrimination` takes `exclude_target_gene=True` at `exclusion_scope="panel"`,
# which drops all 300 panel genes from pds. So this is a BUILD GUARD, not a metric risk: a
# submission whose knockdowns silently vanish on some targets is worth catching anyway.
gpos = {g: i for i, g in enumerate(genes)}
own_rank = []
for r, p in enumerate(perts):
    if blank[r] or p not in gpos:
        continue
    scores = agreement[r]
    own_rank.append(int((scores > scores[gpos[p]]).sum()))
own_rank = np.asarray(own_rank)
if own_rank.size:
    print(f"on-target gene's own agreement rank over {own_rank.size} carried targets: "
          f"median {np.median(own_rank):.0f}, q90 {np.quantile(own_rank, 0.9):.0f}, "
          f"max {own_rank.max()}; {int((own_rank < 2000).sum())}/{own_rank.size} inside "
          f"rank 2000 (alpha > 0.368 at gamma=1)")

# ---- unanimity vs pairwise agreement, split by donor count (097's chal_alpha.py claim)
pair_agree = np.vstack([B.donor_agreement_pairwise(model, p) for p in perts])
for label, mask in (("3-donor", carried == 3), ("2-donor", carried == 2)):
    if not mask.any():
        continue
    for name, mat in (("unanimity", agreement), ("pairwise", pair_agree)):
        frac = (mat[mask] > 0).sum(axis=1) / genes.size
        print(f"graded fraction, {label} targets (n={int(mask.sum())}), {name}: "
              f"median {np.median(frac):.3f}  p5 {np.quantile(frac, 0.05):.3f}  "
              f"p95 {np.quantile(frac, 0.95):.3f}")

# ---- the PER-PERTURBATION half (experiment 097's donor-cosine key), cross-checked against
# 097's own independently measured numbers: carriage 269/290/294, key computable on 294/300,
# challenge median cosine +0.0284, mean(z^2) = 0.327 at gamma 2.
cos_score, cos_carried = B.donor_cosine_key(model, perts)
print(f"\ndonor-cosine key computable on {int(cos_carried.sum())}/{len(perts)} targets "
      f"(097 measured 294/300)")
if cos_carried.any():
    cs = cos_score[cos_carried]
    print(f"  key min/med/max = {cs.min():+.4f}/{np.median(cs):+.4f}/{cs.max():+.4f} "
          f"(097 measured median +0.0284, range [-0.285, +0.439])")
# Is the `damp` divergence LIVE or LATENT on the real panel? 097's gate forces an uncarried
# target's WEIGHT to 0.0 before flooring; this builder scores it 0.0 and ranks it with the
# panel. They differ only if rank01(0.0) < FLOOR is false, i.e. only below the floor where a
# 0.0 cosine actually ranks. Measured, not simulated.
if cos_carried.any():
    cs = cos_score[cos_carried]
    frac_neg = float((cs < 0).mean())
    # rank01 a 0.0 score earns among the carried, this builder's own midrank convention
    allv = np.where(cos_carried, cos_score, 0.0)
    d = B._descending_midrank(allv); asc = (allv.size - 1) - d
    r0 = float(((asc + 0.5) / allv.size)[~cos_carried][0])
    print(f"\nfraction of carried targets with donor_cos < 0: {frac_neg:.3f}")
    print(f"rank01 a 0.0 cosine actually earns on the REAL panel: {r0:.4f}")
    print(f"quartiles: {np.quantile(cs,[0,.25,.5,.75,1.0]).round(4).tolist()}")
    for fl in (0.0, 0.10, 0.18, 0.30, 0.50):
        mine = B.pert_shrink_vector(cos_score, cos_carried, 1.0, floor=fl, missing="damp")
        theirs = np.maximum(np.where(cos_carried, (asc + 0.5) / allv.size, 0.0), fl) ** 1.0
        print(f"  FLOOR={fl:.2f}: rows differing {int((np.abs(mine-theirs)>1e-12).sum()):>2}"
              f"   uncarried z mine {mine[~cos_carried][0]:.4f} / theirs {theirs[~cos_carried][0]:.4f}")

for gamma in (1.0, 2.0, 3.0):
    z = B.pert_shrink_vector(cos_score, cos_carried, gamma)
    zq = np.quantile(z, [0.0, 0.25, 0.5, 0.75, 1.0])
    print(f"  gamma={gamma}: z min/q25/med/q75/max = {zq[0]:.4f}/{zq[1]:.4f}/{zq[2]:.4f}/"
          f"{zq[3]:.4f}/{zq[4]:.4f}; mean z^2 = {float((z**2).mean()):.4f}; "
          f"amplitude compensation 1/sqrt(mean z^2) = {1/np.sqrt((z**2).mean()):.3f}x; "
          f"{int((z < 0.05).sum())} targets below z=0.05")
    zf = B.pert_shrink_vector(cos_score, cos_carried, gamma, floor=0.2)
    print(f"           floored at 0.2: mean z^2 = {float((zf**2).mean()):.4f}; "
          f"compensation {1/np.sqrt((zf**2).mean()):.3f}x")

for gamma in (0.05, 0.25, 0.5, 1.0, 2.0):
    alpha = np.vstack([B.agreement_shrink_vector(row, gamma) for row in agreement])
    kept = alpha[~blank]
    aq = np.quantile(kept, [0.0, 0.25, 0.5, 0.75, 1.0])
    print(f"gamma={gamma}: alpha min/q25/med/q75/max = "
          f"{aq[0]:.3e}/{aq[1]:.3e}/{aq[2]:.3e}/{aq[3]:.3e}/{aq[4]:.3e}; "
          f"{float((kept > 0.5).sum(1).mean()):.0f} genes/target keep alpha>0.5; "
          f"mean alpha {kept.mean():.4f}")

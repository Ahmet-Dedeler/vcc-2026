# INBOX for 178 — a fourth rider, measured after your launch (coordinator, 2026-09-02 19:45)

`experiments/175-donor-pathway-match/deflate175.py`. Per-donor principal directions fitted on a
4,000-perturbation subsample of ALL the donor's rows (not the 85 h1 targets), on the shared DE-universe
gene axis; sign accuracy on the h1-test reference DE pool where all three donors carry the pair
(n = 39,945). Known-answer rows: K562 0.5670 / HCT116 0.5375 / HEK293T 0.5308 / equal mean 0.5697 at
k = 0 — 175's numbers, reproduced.

## The Orion pair shares one dominant direction that K562 does not have

| donor | variance on its own top-1 PC | top-3 |
|---|---:|---:|
| K562 | 5.6% | 10.3% |
| HCT116 | 8.6% | 12.7% |
| **HEK293T** | **19.7%** | 22.3% |

Cosine between donors' top-1 directions: **HCT116 × HEK293T 0.872**; K562 × HCT116 0.307; K562 ×
HEK293T 0.234. Two different tissues from the same atlas share their dominant axis almost exactly,
and the independent screen does not — that is a **platform/protocol axis of the X-Atlas screens**,
not biology. It is the donor-table face of the 26.4% / 44.1% target-independent blob 157/159 found
in the Orion harnesses (which mean-centring removes only if its amplitude is constant across
perturbations — and it evidently is not).

## Removing exactly that one direction per donor, before the consensus

| k deflated | K562 | HCT116 | HEK293T | equal mean | (2,1,1) |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.5670 | 0.5375 | 0.5308 | 0.5697 | 0.5740 |
| **1** | 0.5680 | **0.5454** | **0.5373** | **0.5770** | **0.5782** |
| 2 | 0.5560 | 0.5448 | 0.5339 | 0.5683 | 0.5689 |
| 3 | 0.5539 | 0.5453 | 0.5325 | 0.5650 | 0.5655 |
| 5 | 0.5486 | 0.5383 | 0.5276 | 0.5567 | 0.5564 |

**k = 1 helps every donor and the consensus (+0.0073 on the equal mean, +0.0042 on (2,1,1)); k ≥ 2
removes signal.** The interior optimum at exactly one component, with the direction shared across
the two Orion lines and absent from K562, is the signature of a single artefact rather than of
"low-rank biology" (004 showed the transferable signal is full-rank — this is not a contradiction,
it is one artefactual direction on top of a full-rank signal).

## Why this is a rider for YOUR arm and not a separate experiment

Production whitens on the union covariance, which *down-weights* this direction along with every
other high-variance direction but does not remove it — and 097/172/173 showed whitening's own cost is
paid in sign accuracy. Per-donor rank-1 deflation **before** whitening removes an artefact the
whitener is otherwise spending its budget on. It composes naturally with your three switches and it
is a different kind of move: it removes noise rather than trading `pds` for accuracy, so its `pds`
cost may be ~0 or even negative. 156 closed rank-k deflation **of the prediction** (post-consensus,
post-whitening, did not survive rendering); this is pre-consensus, per-donor, and it targets a
measured artefact. Not the same closure.

**Add it as a fourth env-gated switch** (e.g. `VCC_DONOR_DEFLATE=1`: subtract each donor's own top-1
principal direction, fitted on all of that donor's perturbations on the production gene axis, from
its rows before the consensus; default off = bit-identical). Price it alone and in the composite,
same six members. The PCs must be fitted on the donor's full row set, never on the panel or the
held line — record the fit's seed and subsample so it is reproducible.

Caveat to test rather than assume: the direction may carry some real shared toxicity response that
the recipient also has. Sign accuracy says removing it helps on H1; `pds` and `mse` on the challenge
frame will say whether that holds for the profile members.

## Addendum (19:55) — the axis is identified, and it is absent from H1's truth

Against the h1-test truth (85 targets, centred; its own top PC holds 15.7% of variance):

| donor top-1 PC | cos with h1-truth top-1 PC | cos with h1-truth mean direction |
|---|---:|---:|
| K562 | 0.173 | 0.102 |
| HCT116 | 0.111 | 0.087 |
| **HEK293T** | **0.128** | **0.052** |

The Orion axis is essentially orthogonal to everything in H1's response. **Deflating it removes no
H1 signal.** The caveat in the previous section is resolved on this surface.

**What the axis is.** Top-30 loading genes on HEK293T's top PC, all negative: AUTS2, GRIP1, GPC6,
RAD51B, RORA, EXT1, TTC28, ARL15, NRXN3, EFNA5, PRKG1, LARGE1, NLGN1, KLF12, ERC2, CACNA2D3, IGF1R,
ST6GALNAC3, SLIT2, ADGRL3, AGAP1, WWOX, FAT3, TMEM178B, PARD3B, PTPRG, ROR1, MAGI1, ZFAND3, SFMBT2.
**Every one is a very long, intron-rich gene (0.3–1.7 Mb).** That is the well-known
nuclear-RNA-fraction / long-gene technical axis of single-cell data — cells (or perturbations that
change cell state) with more unspliced nuclear RNA show coordinated long-gene signal. It is a
property of the X-Atlas libraries, shared by both Orion lines (cos 0.872), not a transcriptional
response to any target. Production's whitener is currently spending part of its budget flattening a
sequencing artefact.

Two consequences for your switch:
1. Fit the direction per donor on the full row set, but you may also **check it against the
   long-gene signature directly** (gene length from the GTF / Ensembl, or simply the
   gene-length-vs-loading correlation) — if the correlation is strong, a length-regression
   residualisation is an alternative, more principled form of the same removal and it needs no PCA.
2. **Check the challenge controls.** If the A/B/C control cells show the same long-gene axis of
   variation across cells (they are 10x Flex too), the reference's own DE calls may carry it as a
   target-independent component — which would be a reason NOT to remove it from the donors after
   all, since the scorer would reward matching it. Measure the axis's variance share in the
   challenge control cells before deciding the default.

## Addendum 2 (20:05) — the challenge-controls check, done

Cell-to-cell PC1 of the control cells vs the Orion axis: context_A **0.331** (PC1 = 5.9% of variance),
context_B **0.315** (8.1%), h1-test controls **0.425** (2.9%). So the long-gene axis is present as
*within-population* cell-state variation in Arc's 10x Flex libraries too — it is a component of the
dominant cell-state PC everywhere. **But H1's perturbation EFFECTS do not load on it (cos ≤ 0.13,
previous addendum), while Orion's perturbation deltas do (19.7% of HEK293T's between-perturbation
variance).** In Arc's data perturbations do not systematically shift cells along this axis; in the
X-Atlas screens they do — which reads as batch structure across perturbations rather than biology.

**Conclusion for the switch: removal from the donors is supported.** The scorer's reference (Arc data)
will not reward matching a between-perturbation long-gene shift because its own effects do not carry
one. Default the switch ON in your composite arm and price it; keep the OFF arm as the control.

## Addendum 3 (20:20) — CORRECTION: whitening already removes the axis from the shipped effect

Measured on the 095 h1-test artifact's recovered effect: **energy on the long-gene axis is 9.76% in
the un-whitened consensus and 0.11% in the shipped whitened `H`.** Projecting the axis out of `H`
leaves its sign accuracy unchanged (0.5373 → 0.5373); projecting it out of `cons` lifts `cons` from
0.5449 → **0.5570**. So:

- **The fourth switch will do ≈ nothing to the production profile directly** — the whitener at
  gamma 0.65 already annihilates this direction. Do not expect a `pds`/`mse` gain from it alone. The
  only channels through which it can act are (a) the whitener's *refit* on a covariance without the
  artefact (budget re-allocated to real directions — measure, may be ~0) and (b) **the sign-repair
  rider**, which takes signs from `cons`: use the DEFLATED `cons` there, which is +0.012 more accurate
  and should make 173's +0.002 rider larger. Wire the deflation upstream of the sign source.
- Reprioritise: the composite's value is still in gamma + sign repair + weights; deflation is a
  cheap improvement to the sign source, not a fifth lever.

For the record (mechanism, not a rider): the Orion HARNESS blob 157/159 measured is this axis —
HEK293T truth: blob direction × axis cos **0.832**, top PC × axis **0.877**, **36.3% of truth energy**
on it; HCT116: 0.642 / 0.755 / 13.8%; H1: <2%. That is why the Orion surfaces overstate everything
and can invert sign (168d): their "truth" rewards predicting a sequencing artefact.

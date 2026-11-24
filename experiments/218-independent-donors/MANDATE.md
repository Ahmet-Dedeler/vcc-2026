# 218 — CD4's value is INDEPENDENCE, not lineage. Find the next maximally-independent source.

You are a research/acquisition subagent. Work ONLY in `experiments/218-independent-donors/`
(downloads to `data/external/` under a new subdirectory you name — that is the ONE exception).
No `vcc`, no `cell-eval2`, no submission build, no production changes.
Peak RSS ≤ 12 GB. Download budget ≤ 80 GB; report the inventory before exceeding 15 GB.
`uv run python`. **Post nothing anywhere, contact nobody, create no accounts** (`AGENTS.md` §14).

## READ FIRST
1. `AGENTS.md` **§8.4** — why direction arms were under-priced 4.5× and what that reopens.
2. `experiments/207-cd4-dose/RESULTS.md`, especially its self-residual table — **the measurement
   this mandate turns on.**
3. `experiments/210-cd4-build/BOARD-RESULT.md` — the board proof that this channel converts.
4. `experiments/166-donor-acquisition/RESULTS.md` — the previous sweep, **and the criterion it
   used**, which is the thing this mandate changes.
5. `experiments/159-transfer-ceiling/RESULTS.md` (κ², the `1/k` law), `204` (the Orion↔Orion excess),
   `209` (it is the long-gene axis at rank 2), `205` (X-Cell/Pisces do not exist).

## The argument

`cd4f034` moved the board **+0.0131** and lifted the transfer cosine **×1.249** — the first arm ever
to move our direction channel. So the obvious next question is *what else does that*.

**My first hypothesis was that the mechanism is "un-whitened", not "CD4" — and 207 already refuted
it.** At matched `f = 0.10`, with the same support mask and the same solve:

| arm | Δcos_agg | Δ`pds` |
|---|---:|---:|
| CD4 residual | +0.001966 | **positive** |
| 099's self-residual (our own un-whitened effect) | **+0.000873 (44% of CD4's)** | **−0.007788** |
| target-shuffled CD4 | −0.000149 | −0.004641 |

**56% of the cosine and 100% of the `pds` is CD4-specific**, and the self-residual has the *opposite
sign* on `pds`. Whitening-relaxation is real but dominated. **CD4 is carrying genuinely new
information, and the question is what property makes it do that.**

**It is not lineage.** 207 measured the gain at **+7.6% to +15.2% on four unrelated non-T lines** —
lineage-generic. So the T-cell match to context A (Jurkat) is not the mechanism.

**The hypothesis this mandate tests: the property is INDEPENDENCE.** Our three production donors are
**two platforms** — Replogle K562, and Xaira Orion ×2 — and 204/209 measured that the Orion pair
shares an excess subspace the third donor does not (the long-gene technical axis at rank 2). CD4 is
a **third platform** *and* **primary human cells** rather than an immortalised line. Under 159's
`1/k` law what counts is *independent* lines; a source that shares neither platform nor
transformation state with the existing three is worth more per unit of coverage than a fourth
cancer line would be.

**166 searched for the wrong thing.** Its criterion was *genome-wide human Perturb-seq covering the
panel*, and it correctly concluded none remained. **This mandate's criterion is different:
maximally independent from Replogle/Orion, with usable panel coverage** — which admits primary
cells, iPSC-derived and differentiated systems, organoids, other chemistries (CRISPR KO, CRISPRa,
base editing), other readouts, and other labs' platforms. CD4 itself would have been missed by
166's criterion and it is our single best component.

## Work

### Phase 0 — quantify what "independent" is worth, on data we already hold (do this FIRST)
Before searching, **measure the thing you are searching for.** For each existing source
(`k562_gwps`, `orion_hct116`, `orion_hek293t`, CD4, and the training-only tables `nadig_jurkat`,
`rpe1`, `k562_essential`, `nadig_hepg2`, `frangieh_melanoma`, `h1_2025`), compute a **pairwise
independence matrix** on shared targets: effect-vector cosine, the shared-subspace energy 209
measured, and platform/lineage/transformation metadata. Then regress **measured marginal cosine
contribution** against independence, using the points we have:
- CD4 → +11.0% at `f = 0.10` (207);
- the self-residual → +4.8% (207);
- the Orion pair's mutual excess (204: z +5.1 to +6.8; 209: it is the long-gene axis).

**If independence does not predict contribution on the sources we already hold, this mandate's
premise is dead and you should say so and stop.** That is the honest known-answer row and it costs
one afternoon rather than 80 GB.

### Phase 1 — the search, under the new criterion
Sweep for human perturbation-expression resources that are **maximally independent** of
Replogle/Orion and carry usable panel coverage. Cover at minimum: HuggingFace, GEO, ArrayExpress /
BioStudies, Zenodo, figshare, the Arc Virtual Cell Atlas, IGVF, CZ CELLxGENE, scPerturb, and
Europe PMC / bioRxiv since 2025-06. **Explicitly include** primary-cell screens (the CD4 atlas is
the existence proof), iPSC-derived and differentiated systems, organoids, CRISPRa and base-editing
screens, and non-10x chemistries.

For each hit report: organism, cell system, transformation state, platform/chemistry, perturbation
count, **overlap with `data/challenge/pert_counts.csv` (300)**, cells per target, licence, size,
URL, and a one-line independence argument. **Run 166's `census.py` for the coverage number** — do
not eyeball it. 166's own correction is instructive: it dismissed X-Atlas/Pisces' Jurkat arm as a
duplicate of `nadig_jurkat`, but that table is a common-essential library at **0/300**, so it was
never a duplicate.

### Phase 2 — build and measure the best candidate, if there is one
Build a donor table with `scripts/build_donor_table.py` (do **not** modify anything under
`data/donors/`; write to your own directory), then measure it **the way CD4 is measured**: as an
**un-whitened post-whitening residual** at a swept fraction, priced on **Δcos on PLOLO** (204:
±0.20% relative) with `pds` secondary, plus `q_w`.

**Price it under §8.4's corrected rule**, not the old one: `board_cosine ≈ 1.80 × PLOLO_cosine`, and
a direction arm's overall ≈ **2.51 ×** its `mse`-only price because `nmae`, `fid`, `reach` and `jac`
all move with sign accuracy. **Never carry those four at zero for a direction arm** — that
convention is exactly what cost this project eleven closures.

Falsifiers: a **target-shuffled** control at every fraction (207's CD4 control was strictly negative
at all 12 doses with no permutation tie); the **leak-free `kolf2` fold**; and a **non-vacuity
assertion**, because the median-based sizing rule can zero a residual silently when coverage is thin
(200 found it, 207 verified it, and it is the failure mode most likely to waste this run).

## Deliverable
`RESULTS.md`: the Phase-0 independence-vs-contribution regression (**the load-bearing result**), the
Phase-1 inventory with coverage numbers and licences, and Phase 2's measurement if any candidate
earned it. **"Independence does not predict contribution" or "no independent source with usable
coverage exists" are full results** — the first would refute the mechanism, the second would close
donor acquisition properly under the right criterion instead of the wrong one.

## Territory
216 owns the CD4 dose and form and the published-atlas question. 217 is building `f = 0.80`.
213 owns the transport law. **You own donor independence and acquisition.** Do not sweep the CD4
dose — that is 216's, and duplicating it wastes the machine.

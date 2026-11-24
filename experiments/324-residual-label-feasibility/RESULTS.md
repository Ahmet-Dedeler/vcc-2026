# 324 — Multiple studies and lines are not crossed replication

Completed September 9, 2026 by Codex alone. **Metadata-only**: no expression,
DE table, pseudobulk, reserved outcome, model fitting or new data download.
Protocol/code/input hashes were committed at `1cbda12b` before measurement.
Seven known-answer metadata checks and seven synthetic/unit tests pass.

## Decision

**PROVENANCE_UNRESOLVED_NO_RESPONSE_LABEL_PILOT.**

The proposed independently repeatable context-residual label is not justified
by the current metadata inventory. GPSAdb's historical **50 panel targets
observed in multiple accessions and multiple cell lines** do not form even one
target-specific **two-context × the same two-accession** rectangle.
Only **six** panel targets have same-context, same-method multi-accession
observations, and neither sample/control overlap nor construct/publication
independence is documented in these metadata fields.

This does **not** close all response prediction, unbalanced multi-study
models, single-study guide-controlled development, or all external datasets.
It rejects interpreting nominal study/line counts as evidence for this
particular matched independent-residual design. Metadata alone cannot
establish repeatability or statistical power at any sample count.

Primary record: `work/results.json`. Exact supporting groups are in
`work/panel__id_fallback__same_method__pairs.csv`; four complete panel
pair/rectangle tables preserve the identity and method sensitivities.

## 1. What was counted

270's retained browse table reproduces **7,665 dataset groups**, **4,623
accessions**, **329 panel groups**, **96/300 panel targets**, and **50 panel
targets with at least two accessions**. All 50 also have two raw line labels.
Multiple treat/control comparisons within one accession are deduplicated
before counting study edges.

Cellosaurus identifiers are used when present; exact line labels are retained
as explicit fallbacks. An identifier match is not evidence of matching
culture conditions, treatment duration, construct, or control population.
Accessions are **accessions**, not automatically independent publications.

| population / identity / method restriction | targets repeated in the same context | repeated target-method-contexts | accession pairs | crossed rectangles |
|---|---:|---:|---:|---:|
| panel, raw line names, methods mixed | 9 | 15 | 19 | **0** |
| panel, IDs with name fallback, methods mixed | 9 | 15 | 19 | **0** |
| panel, IDs with name fallback, same method | **6** | **10** | **10** | **0** |
| panel, stable IDs only, same method | **6** | **9** | **9** | **0** |
| all metadata, IDs with name fallback, same method | 169 | 209 | 334 | **2** |
| all metadata, stable IDs only, same method | 153 | 190 | 311 | **1** |

With methods mixed, the full metadata has five rectangles, still zero on the
panel. Across all metadata the two same-method rectangles are **MESP2**
(knockout, exact-name iPSCs/PSM) and **RUVBL1** (shRNA, two stable line IDs).
Their accession pairs are respectively GSE116930/GSE116934 and
GSE224646/GSE224654. Their independence is **unresolved**; different accession
numbers alone do not establish it. Neither target is in the challenge panel.

## 2. The six panel targets, without promoting their independence

| target | method | contexts with accession repeats |
|---|---|---|
| EZH2 | knockout; shRNA separately | K-562, HAP1; MCF7, NHDF (name only) |
| KDM5A | shRNA | HEK293T |
| KHSRP | shRNA | K-562, Hep G2 |
| MTA1 | knockout | HCT 116 |
| SRSF5 | shRNA | Hep G2 |
| TAF15 | shRNA | K-562 |

Even KHSRP's two contexts use different ENCODE accessions on one side:
that is not a literal two-accession rectangle. It **could** become a
publication-level crossing if source-level provenance linked those
experiments appropriately; this inventory has not checked that. Therefore
the zero-rectangle count is not a proof that a broader matched design cannot
be recovered. This is precisely the distinction the next source check must
resolve.

Missing metadata fields: publication ID, sample IDs, control sample IDs,
construct ID, guide sequence and treatment duration. Dataset-group
`nSample` counts do not supply these identities. Same method is necessary
for the registered sensitivity, not sufficient to identify equivalent
perturbations.

The panel's 329 groups include only **seven CRISPRi** groups; the rest are
136 shRNA, 120 knockout, 60 siRNA and six inhibition groups. The same-method
repeat table contains no CRISPRi pair. A bulk RNAi/KO repeat would not
automatically qualify as a CRISPRi/Flex validation surface.

## 3. What the existing single-cell inventory supplies

These are **inherited metadata facts**, not new response measurements:

- Replogle K562 GWPS/essential: **2,053** shared targets, **zero** shared
  challenge-panel targets, two experiments in **one publication**.
- Nadig Jurkat/HepG2: two lines with **zero panel overlap** each.
- Broad: **12 working lines**, one publication, guide/cell splits already
  analyzed by 268/285. Reserved lines remain excluded. The known line/batch
  nesting is not resolved by adding random cell halves.
- CD4: four biological donors, three activation states, one publication.
  These are valuable units, but not independent study replicates.

Sources: 285's frozen `work/inv285.json`, 268's `work/units268.md`,
and their previously read `RESULTS.md`. No inventory claim is an exhaustive
September 9 release search.

## 4. Avoid manufacturing a biological residual

The deterministic fixture has zero true context residual and orthogonal
response errors. Subtracting one shared noisy estimate gives a residual
cross-product of **4**; independent estimate views give **0**, the correct
value. This demonstrates a design failure mode, not a finding that 268's
measurements contain it.

For predicting the error of a *fixed baseline*, its shared bias can itself
be a real, repeatable modeling target. The problem is calling that error
**biological context specificity**. A later pilot must declare its estimand
and use independent control/shared-response views or a validated correction
when making the biological claim.

Novelty stays scoped to this crossed-metadata inventory and its design gate.
132/161/224/268/269/283/285/270 were checked before registration; no scalar
sign-stability feature, basal ridge, or successful 224 diagonal is relabeled
as a new method. No new result enters the manuscript.

## 5. Concrete next step

**One metadata-only provenance pass on the ten same-method panel accession
pairs**, prioritizing KHSRP's two contexts and the nearby-accession MTA1
pair. Use official GEO/ENCODE study and sample metadata only: map accession
to publication/experiment, identify shared controls/samples, and identify
constructs and treatment timing where explicitly recorded. A missing field
stays unknown. Check whether accession grouping conceals a genuine
publication-level crossing rather than treating accession counts as final.

This is a separately registered next action, **not performed here**. No DE
or count download, reserve reopening, broad retraining, purchase, or
submission is warranted by the present result. Only verified independent
units could earn a subsequently registered development repeatability pilot.

## 6. Reproduction and state

```sh
uv run --frozen python -m unittest discover \
  -s experiments/324-residual-label-feasibility -p 'test*.py' -v
uv run --frozen python experiments/324-residual-label-feasibility/audit324.py run
```

Input manifest SHA-256:
`5894e2306dd45491de9b3518e60c0a6abc0dbee60445944f6bf33c73be3cb3eb`.
Retained result/CSV output: approximately **40 KiB**. Stdlib only.
Repeat run: **0.12 s**, peak process RSS **40.7 MiB**, results JSON
SHA-256 unchanged (`e8974ede554d38f8b71d03704c631bf902946b6ea1da59d338ee5e52b8a31376`).
CSV line endings are normalized in Git; the frozen computation is unchanged.
All six competition metrics are **unmeasured**; no prediction was changed.
No production edit, submission, reserve access, paid resources or push.

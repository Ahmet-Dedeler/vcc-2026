# 321 — RESULTS: the manuscript corrected where 320 shows it overclaims, and Round 7 folded in

**Reading and editing only.** No `vcc`, no `cell-eval2`, no build, nothing written under `submissions/` or
`data/`. `uv run python` was used only to read board-entry JSON records and do arithmetic on them. Local peak
RSS negligible. **Nothing was published** — the rebuilt page sits at
`experiments/321-writeup-correction/work/manuscript-v10.html` and the coordinator owns the republish.

---

## THE LINE

**Every claim 320 named as too strong is now weaker or gone.** F16 is presented as a **re-measurement**
(NEW 0 / EXTENDED 1 / KNOWN 4 / CONTRADICTED 0) with its provenance table in the body; the partner-construct
control is credited to Jackson 2003 and is **no longer claimed as ours**; the headline p is stated as a
**gene-sampling** p in four places; the **$3–6 k sequence-mutant recommendation is withdrawn** because
Rohatgi 2024 already demonstrates the mechanism at the level of dCas9 binding in this library's parent
collection; six references are added with DOIs; and §3(h) carries Round 7's five scored entries, including the
new best at **0.2186597**. **23 passages changed across the three manuscript files, plus the build script's
masthead. Every number in them was read from the file it came from and logged in `work/verify321.log`. Not one
number was added that strengthens a claim.**

---

## 1. Files changed

| file | passages | authority |
|---|---|---|
| `experiments/294-writeup/PAPER.md` | 16 | mandate §"you may edit" |
| `experiments/294-writeup/CLAIMS.md` | 6 (C-55, C-56, C-57 and the preamble) | mandate §"you may edit" |
| `experiments/294-writeup/REFERENCES.md` | 1 new section, 5 new rows + 1 augmented + 1 withdrawal note | mandate §"you may edit" |
| `experiments/294-writeup/build_html.py` | masthead only (2 strings) | mandate §Rules: *"update its masthead the way earlier rounds did"* — **flagged here because this is a fourth file, outside the three named; no logic was touched** |

Deliverables written inside `experiments/321-writeup-correction/`: `work/verify321.log`,
`work/manuscript-v10.html` (1,627 KB, 26 TOC entries, 8 figure captions), this `RESULTS.md`.

---

## 2. Every sentence changed — old → new

### Correction 1 — the partner-construct control is not ours (320 §3(a), §3(b))

**2.1 `PAPER.md` §3(g), the partner paragraph lead.** 320's wording used, not mine.

> **OLD:** *The partner-construct test is the design's load-bearing element.* Each construct's candidate genes
> were re-scored in the cells of its **partner construct** …

> **NEW:** *The partner-construct test is the design's load-bearing element — and it is not ours.* Reading the
> same candidate genes in the cells of a different guide against the same target is the standard specificity
> control of the seed off-target field: it is [JACKSON03]'s founding observation that different siRNAs against
> one transcript give *siRNA-specific rather than target-specific* signatures, and it is used on CRISPRi
> Perturb-seq for this exact purpose by [HARTMAN26] and as the discovery device by [ROHATGI24]. What is ours is
> only its systematic paired application to all 34 constructs. Each construct's candidate genes were re-scored …

320 §3(b) is followed exactly: "load-bearing element" is **kept** (it describes the role in *this* design and is
accurate); only the authorship is corrected.

**2.2 `PAPER.md` §3(g), the F16 marker.** The sentence claiming the control is replaced by 320 §3(g)'s
replacement text verbatim.

> **OLD:** … Provenance: this result postdates both the 287 audit and its 297 extension and therefore carries no
> audit verdict; under §1.4's rule we say so rather than claiming novelty. What is not new is the phenomenon —
> seed-mediated off-target CRISPRi activity is characterised and quantified in the literature ([SEED24],
> [HARTMAN26]). What is ours is a preregistered aggregate test with the *construct* as the experimental unit on
> this particular library, a selection-matched joint null with an attainable floor three orders of magnitude
> below α, **and a partner-construct control.**

> **NEW:** … Audited by 320: **KNOWN 4, EXTENDED 1, NEW 0.** The phenomenon ([ROHATGI24]; [HARTMAN26]), the
> seed-plus-promoter rule and its measured predictive precision, the partner-guide specificity control
> ([JACKSON03]; [HARTMAN26]) and the aggregate-over-per-item statistical logic are all published. What is ours
> is the design — construct as the unit, with a joint null over candidate sets matched on a measured
> library-wide off-target propensity — and a re-measurement of the phenomenon on the challenge's own library,
> whose specificity its producers did not assess.

**2.3 `PAPER.md` §3(g), the headline clause** (320 §3(c) — the "protospacer not target" observation is the
field's premise, not our result).

> **OLD:** … in aggregate, **and the association follows the protospacer rather than the intended target.**

> **NEW:** … in aggregate, **and, as the seed off-target literature has held since Jackson *et al.* 2003
> [JACKSON03] and as Rohatgi *et al.* 2024 report for CRISPRi [ROHATGI24], the effect follows the protospacer
> rather than the intended target.**

**2.4 `PAPER.md` §4.1** — the second sentence that claimed the control.

> **OLD:** **One result postdates both audits and therefore carries no provenance verdict: F16 (§3(g)).** We
> state its status rather than assert its novelty. … what we have not found done before on this library is the
> construct-as-unit aggregate test with a selection-matched joint null **and a partner-construct specificity
> control**, and until that has been audited … that is all we claim.

> **NEW:** **The one result that postdated both audits has since been audited, and nothing in it is new: F16
> (§3(g)).** 320 graded it **KNOWN 4, EXTENDED 1, NEW 0.** [320 §3(g)'s replacement text, as in 2.2] … It
> therefore belongs in §4.2 with the other re-measurements and not in this section; an earlier revision of this
> paragraph reserved judgement on its novelty, and the reservation is now spent.

The honesty-condition sentence (C8 aborting itself, the repair moving the estimate down 9.5 %) is kept verbatim.

**2.5 `PAPER.md` abstract.** "while the identical genes read at background level in the cells of the *partner*
construct …" now ends "…, **which is the field's standard specificity control and not ours**."

**2.6 `PAPER.md` §4.4.** "the partner-construct control attaches it to the protospacer sequence" → "the
partner-construct control **is a standard one borrowed from [JACKSON03] and [HARTMAN26]**, attaches the
association to the protospacer sequence".

### Correction 2 — F16 as a re-measurement, with its provenance table

**2.7 `PAPER.md` §3(g)** — a new paragraph and a new five-row table immediately after the marker:

> *The provenance, item by item (320).* Nothing biological in §3(g) survives as new; the one **EXTENDED**
> column is the design, and it is a recombination of two published methods rather than a new one.

with rows: the phenomenon (**KNOWN**, [ROHATGI24] with ChIP-seq occupancy, 98–1,141 peaks per guide, 10–33 %
promoter-proximal, in the parent library two years earlier; [HARTMAN26]); the 12-nt seed + NGG promoter rule and
its **published low precision** (**KNOWN**; 21 of ≈ 562 at one promoter; 47.2–59.8 % of guides matching against
0.83 % with a nominated event); the construct-as-unit matched-set design (**EXTENDED**, a recombination of
[SEEDMATCHR24] and [MATCHRANGES23], new only in the *measured* 209,070-spacer covariate, **resting on a negative
literature search**); the partner-construct test (**KNOWN**, [JACKSON03]); aggregate-beats-per-item
(**KNOWN**, [GOEMAN07]).

**2.8 `PAPER.md` §4.2** — "Four items are publishable only as re-measurement" → "**Five** items … (287 §2,
Tier 2, **plus F16 from 320**)", with a new lead bullet crediting Rohatgi, Hartman, Jackson, SeedMatchR and
matchRanges and stating that what is ours is the design and the sample.

**2.9 `PAPER.md` §1.4** — the provenance-discipline sentence now names the second audit:

> **OLD:** … the audit graded its 24 items **NEW 4, EXTENDED 11, KNOWN 7, CONTRADICTED 2**, with one
> contradiction subsequently reconciled — and where a claim postdates the audit and has no verdict, we say so …

> **NEW:** … with one contradiction subsequently reconciled, **and a sixteenth finding that postdated it (F16,
> §3(g)) was audited separately at NEW 0, EXTENDED 1, KNOWN 4, CONTRADICTED 0** — and where a claim postdates an
> audit and has no verdict, we say so …

**2.10 `PAPER.md` abstract** — "we measure in the benchmark's own 2025 CRISPRi library" → "we **re-measure a
published phenomenon on a new sample** — the benchmark's own 2025 CRISPRi library, whose specificity its
producers did not assess"; and "it does not establish off-target binding" → "it does not establish off-target
binding **— already shown by dCas9 ChIP-seq in this library's parent collection**".

**2.11 `PAPER.md` §4.5 T7(iii)** — "…passing at S = +0.713 … after a calibration repair that moved the estimate
down 9.5 %." → "… **— a re-measurement of a published phenomenon (320: NEW 0), not a discovery, and with a
headline p that resamples genes rather than constructs.**"

**2.12 `PAPER.md` plain-language summary** — the last F16 sentence is replaced (see 2.16, where the withdrawal
lives too), adding: "**and it is not a new discovery. The effect, the sequence rule and the two-guide control
are all in the published literature; one of those papers watched guides physically stick to the wrong promoters
in the very library ours descends from, two years earlier. What is new here is only that it has now been
measured on this particular library, which nobody had checked.**"

### Correction 3 — the headline p is a gene-sampling p (320 §3(d))

**2.13 `PAPER.md` §3(g)** — a new paragraph inserted after the C9 result, using 320's sentence:

> *What that p resamples, and what it does not (320).* That p is a **gene-sampling** p — the resampled unit is
> the candidate gene, not the construct — and its three-orders-below-α resolution is a property of the 100,000
> draws, not of the 34 units. In [GOEMAN07]'s terms it is a *competitive, gene-sampling* null, and such nulls
> are anti-conservative under residual dependence between the resampled items; the covariate matching removes
> the confounders it measures, but it does not make genes exchangeable. The construct-level readings are the
> construct bootstrap interval **[+0.32419, +1.07351]**, the paired partner-construct Wilcoxon **p = 0.00459**,
> and, in C8 only, the wrong-recipient construct permutation at its floor **p = 2.0 × 10⁻⁴**; under C9's
> per-construct centring that permutation is the same object as KA-4R and is no longer independent, so C9's
> construct-level inference rests on the bootstrap interval and the partner Wilcoxon. Those, not 2.0 × 10⁻⁵, are
> the numbers a construct-level null supports. This weakens nothing about the conclusion — the interval excludes
> zero and the paired test is significant — but it corrects which number carries the weight, and the phrase "the
> construct as the experimental unit" describes the *estimand and the interval*, not the headline p.

**The construct-unit design is not deleted**, as the mandate requires: the estimand, the bootstrap and the n = 34
interval-width argument in §4.6 are untouched; only what the *p-value* resamples is now stated.

**2.14 `PAPER.md` §4.4** — the F16 bullet gains: "The headline p = 2.0 × 10⁻⁵ is a **gene-sampling** p, not a
construct-level one; the construct-level readings are the bootstrap interval and the paired partner Wilcoxon
(§3(g))."

**2.15 `PAPER.md` §3(g), the null paragraph** — the comparability note 320 §3(f) asks for, plus one referent fix
("That is the design difference" → "**That floor** is the design difference", because a sentence was inserted
between):

> **ADDED:** This rate is not comparable with published seed-match rates: our per-spacer L = 12 promoter-hit
> probability is **0.132353** for the 68 real 2025 spacers and **0.160879** for the pool, against the 47.2–59.8 %
> of guides [HARTMAN26] report with a ≥ 12 bp match within 1 kb of a TSS, because our rule is confined to the
> CRISPRi active window (TSS offset −50 … +300) and applies a 10 kb distal filter. Ours is the stricter rule, by
> design.

### Correction 4 — the wet-lab recommendation withdrawn

**2.16 `PAPER.md` §3(g)**, the whole closing paragraph.

> **OLD:** *What would resolve it, and what it costs.* The decisive experiment is a **sequence-mutant control**:
> … About 24 constructs × ~500 cells, one 10x lane, order **$3–6 k** plus cloning. A **dCas9 occupancy assay**
> …, order **$5–10 k**, answers "is it bound" rather than "is it repressed". … None of these is more analysis of
> these cells, and no amount of that would settle it.

> **NEW:** *What a new experiment would add, and the earlier recommendation withdrawn.* An earlier version of
> this section proposed a **sequence-mutant control** … as *the decisive next experiment*, and a dCas9 occupancy
> assay at $5–10 k beside it. **That recommendation is withdrawn as written.** The mechanism it was meant to
> establish is not open: [ROHATGI24] already demonstrates seed-mediated off-target activity at the level of
> **binding**, with dCas9 ChIP-seq on four seed-off-target guides (98–1,141 reproducible peaks each, 10–33 %
> promoter-proximal, 12 bp of seed complementarity sometimes crucial for strong occupancy), in this library
> family, two years before this measurement. An occupancy assay here would add construct-specific confirmation,
> not the mechanism. The sequence-mutant control remains the only design that separates seed-mediated binding
> from a construct artefact travelling with spacer identity **for these particular constructs**, so what it buys
> is local attribution on a 34-construct sample of a phenomenon already characterised — **which is little, and
> we say so rather than price it.** The cheap route, if it ever exists, is an independent screen carrying the
> same library … §3(h) reports that none does. What is certain is that no further analysis of these cells
> settles it.

The mandate asked for "an honest statement of what, if anything, a new experiment would add … and if the answer
is 'little', **say that**." **The answer is little, and the paper now says so in those words.**

**2.17 `PAPER.md` §4.5 T7**, the closing sentence.

> **OLD:** The decisive next step is a sequence-mutant control, not further analysis of these cells (§3(g)).

> **NEW:** The earlier recommendation of a sequence-mutant control as *the decisive next experiment* is
> withdrawn: [ROHATGI24] already demonstrates the mechanism at the level of dCas9 binding, with occupancy data,
> in this library's own parent collection, so a new experiment here would buy local attribution on 34 constructs
> rather than a mechanism the literature lacks (§3(g)). No further analysis of these cells settles it either.

**2.18 `PAPER.md` plain-language summary** — "settling which of the two it is needs a new experiment — re-making
a few guides with the matching part silently changed — rather than more analysis of these cells" → "**An earlier
draft of this paper called for a new laboratory experiment to settle the mechanism; that call is withdrawn,
because the mechanism was already settled elsewhere.**"

**2.19 `CLAIMS.md` C-57**, the "what would move it" cell — the $3–6 k / $5–10 k prescription is replaced by the
withdrawal, keeping only the free route (an independent screen carrying this library, which does not exist).

### Correction 5 — Round 7's board results in §3(h)

**2.20 `PAPER.md` §3(h)** — a new paragraph, a four-row table and three readings appended after the existing
leaderboard sentence (which is **left exactly as it was**; see §4 for why). Every score is that entry's own
`score_avg` from its captured record:

| entry | form | `score_avg` | Δ vs 0.2160944 |
|---|---|---:|---:|
| `Q519I2pLWRCW7PAEx5ul` | `exp` = −1.5, as shipped | 0.2177189 | +0.0016245 |
| `1n5LFk6PSEHq15UHENIE` | `exp` = −3.0 | 0.2164511 | +0.0003567 |
| `2K5ppxkowMSj2ugoZqYb` | parity-correct, exactly shape-only where the norm match is applied | 0.2169733 | +0.0008789 |
| **`xuDelEik7zGWMXdNu0R4`** | **suppression only, `r ← max(1, r)`** | **0.2186597** | **+0.0025653** |

and the three readings: the dose curve **turns over at −1.5** (−3.0 keeps about a fifth of the gain); the parity
arm **decomposes** the shipped gain into weighting **+0.0008789** and an accidental per-target magnitude
redistribution **+0.0007456**, so the correct form scores below the incorrect one; and the cap form is
**+0.0009409** over the plain form and **1.58×** the previous best arm, paying `pds` (−0.00208) to buy `nmae`
(+0.00475) and `reach` (+0.00193). Then the three axes closed on the real six-member scorer: **emission**
+0.000131 = 7 % of the +0.0019 shape-only gate; **amplitude / operating point** +0.002374 = 27 % of the +0.0089
gate *even counting `mse` free*, at 28.6–40.8 : 1 against; **per-target magnitude** +0.000093 = 1.0 % of the gate
at z = −0.06 against three entropy-matched shuffles of its own multipliers, the best shuffle (+0.000543) beating
the real profile. The paragraph opens "recorded for completeness and not because they are science" and closes
"None of this bears on §3(a)–(g), and none of §3(a)–(g) was built into any of these entries."

**2.21 `PAPER.md` header** — a Round 7 revision block (Draft v10) naming the four corrections, in the style of
the Round 5 and Round 6 blocks above it.

### `CLAIMS.md` — the C8 / C9 provenance rows

- **preamble:** "Rows C-55, C-56 and C-57 postdate both the audit and its extension and therefore carry **no
  provenance verdict at all**" → "Rows C-55 and C-56 … were **audited separately by 320** … **KNOWN 4,
  EXTENDED 1, NEW 0, CONTRADICTED 0** — F16 is a re-measurement of a known phenomenon on a new sample … C-57
  claims nothing and needs no verdict."
- **C-55 provenance:** "postdates the 287 audit … no audit verdict … The phenomenon … is published: [SEED24],
  [HARTMAN26]. The construct-as-unit aggregate test on this library, with a selection-matched joint null **and a
  partner-construct control**, is what is not" → the full 320 grading, item by item, with the partner control
  moved into the **KNOWN** list and the design named as the only EXTENDED column, "resting on a negative
  literature search".
- **C-55 unit · n:** now ends "**The unit of the estimand and of the bootstrap is the construct; the joint null
  resamples GENES**".
- **C-56 provenance:** "as C-55 — **no audit verdict**" → "as C-55 — **320: KNOWN 4, EXTENDED 1, NEW 0.**
  Nothing biological survives as new; the design is the one EXTENDED column".
- **C-56 interval:** gains the gene-sampling statement, [GOEMAN07]'s anti-conservativity caveat, the two
  construct-level readings, and the note that C8's wrong-recipient permutation is not independent under the
  repair.
- **C-57:** see 2.19.

### `REFERENCES.md` — six references added

New section **"Added from the F16 novelty audit (320, 2026-09-08)"**, every row transcribed from
`experiments/320-f16-novelty/work/refs320.md` with its provenance code; **nothing was fetched from the web**.

| key | citation | DOI |
|---|---|---|
| **ROHATGI24** | Rohatgi N, Fortin J-P, Lau T, Ying Y, Zhang Y, Lee BL, Costa MR, Reja R. *Seed sequences mediate off-target activity in the CRISPR-interference system.* Cell Genomics **4**(11):100693 (2024) | 10.1016/j.xgen.2024.100693 (PMID 39510079) |
| **JACKSON03** | Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS. *Expression profiling reveals off-target gene regulation by RNAi.* Nat Biotechnol **21**(6):635–637 (2003) | 10.1038/nbt831 (PMID 12754523) |
| **SEEDMATCHR24** | Cazares T, Higgs RE, Wang J, Ozer HG. *SeedMatchR…* Bioinformatics **40**(1):btae011 (2024) | 10.1093/bioinformatics/btae011 |
| **MATCHRANGES23** | Davis ES, Mu W, Lee S, Dozmorov MG, Love MI, Phanstiel DH. *matchRanges…* Bioinformatics **39**(5):btad197 (2023) | 10.1093/bioinformatics/btad197 (PMID 37084270) |
| **GOEMAN07** | Goeman JJ, Bühlmann P. *Analyzing gene expression data in terms of gene sets: methodological issues.* Bioinformatics **23**(8):980–987 (2007) | **none recorded in `refs320.md` — none written; a DOI was not invented** |
| **HARTMAN26** | already present in §A; the audit's readings (the ≥ 12 bp threshold, 47.2–59.8 % vs 0.83 %, the partner-guide control, unit = guide–gene pair) are added beside it, with its [S] author-list caveat repeated | bioRxiv 10.64898/2026.03.27.714658v2 |

Plus a **withdrawal note**: `[SEED24]` was an [S] (search-index) row with no author list and no DOI, and is the
same paper as `[ROHATGI24]`, which 320 read in full on the preprint. `PAPER.md` and `CLAIMS.md` now cite
`ROHATGI24`; `SEED24` is retained in `REFERENCES.md` only so earlier revisions resolve. **`[SEED24]` no longer
appears in `PAPER.md` or `CLAIMS.md`.**

---

## 3. Judged too weakly supported to keep, or deliberately not written

1. **"…and a partner-construct control" as a claim of ours** — removed from `PAPER.md` §3(g) and §4.1 and from
   `CLAIMS.md` C-55. It is 23 years old ([JACKSON03]) and already in use for this exact purpose on CRISPRi
   Perturb-seq ([HARTMAN26]).
2. **"this result … carries no audit verdict"** — removed everywhere (§3(g), §4.1, `CLAIMS.md` preamble, C-55,
   C-56). It now has one, and it is NEW 0.
3. **"The decisive experiment is a sequence-mutant control … order $3–6 k"** — removed as a recommendation from
   §3(g), §4.5 and `CLAIMS.md` C-57, and from the plain-language summary. Kept only as the description of what
   such an experiment *would* buy (local attribution), with the explicit judgement that it is little.
4. **The bare "p = 2.0 × 10⁻⁵ … three orders of magnitude below α" as the construct-level result** — the
   sentence is literally true about Monte-Carlo resolution and is kept, but it no longer stands alone: four
   passages now say what it resamples.
5. **313's fixture numbers** (delivered row-norm ratio 0.806–2.155, whole-matrix +6.0 %) — **not written**.
   They are Correction 5 in `docs/LOG.md`: measured on a synthetic 20-target × 60-gene fixture, while the real
   panel reads 0.9321–1.2349 and 1.0085 (318). Neither pair appears in the manuscript; §3(h) names the *stage*
   at which the parity arm is shape-only instead ("exactly shape-only **where the norm match is applied**").
6. **A live rank or gap for the new best entry** — **not written.** The newest leaderboard snapshot on disk
   (`data/leaderboard/2026-09-08T0353Z.json`, 787 teams, leader 0.2466093) predates the 05:08:50Z submission of
   `xuDelEik7zGWMXdNu0R4`, and the entry's own record carries a rank but no team count and no leader score.
   Mixing the two epochs would produce a number that is true of neither, so §3(h)'s existing 2026-09-07
   15:48:55Z board sentence is left untouched and the new paragraph quotes only epoch-free per-entry scores.
7. **320 §2's "one place F16 is too hard on itself"** (the non-monotonic seed-length ordering is *consistent*
   with the published account, which hedges the 12-nt threshold) — **not applied.** It would strengthen the
   draft, and this edit only removes claims. The existing text already hedges correctly: it says the ordering
   sits awkwardly with a "*simple* seed-binding reading", which is what 320 says is true.
8. **317's cost and worker count, and 318's six-times-dose −0.001101** — read and verified, not written;
   engineering housekeeping that §3(h) does not need.

---

## 4. Limits of this correction

- **One file outside the three named was touched:** `build_html.py`'s masthead (eyebrow date; "Sources:
  FINDINGS v3.1 · … audits 287, 297, 304, 320"; the Round 7 line and "Draft v10"). The mandate's own rebuild
  instruction directs this ("update its masthead the way earlier rounds did") and the coordinator republishes
  from that script, so leaving it at v9 would have shipped a stale masthead. **No logic, CSS or figure handling
  was changed.**
- **The audit's own limits are inherited.** 320 read Rohatgi 2024 on the **bioRxiv preprint** (both journal
  routes returned HTTP 403), Hartman 2026 is a **preprint** whose author list is search-index-derived, and
  item 3's EXTENDED rests on a **negative literature search**. `REFERENCES.md` and the §3(g) table both say so;
  the paper does not claim the design is unpublished, only that the audit did not find it published.
- **No re-analysis was done.** Every F16 number in the paper is the number that was already there; 321 changed
  what is claimed *about* those numbers, and added only the four board scores, their deltas and 317/318's
  closure figures.
- **Nothing was published.** The page is built at `work/manuscript-v10.html` and not deployed.

---

## 5. Files

| file | what |
|---|---|
| `RESULTS.md` | this: every sentence changed old → new, every reference added, everything dropped |
| `work/verify321.log` | every number written, with the file and line it was read from (§A board, §B literature, §C DOIs, §E what was deliberately not written) |
| `work/manuscript-v10.html` | the rebuilt page, 1,627 KB, 26 TOC entries, 8 figures. **Not published.** |

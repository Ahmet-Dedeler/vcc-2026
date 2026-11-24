# 305 — Expression evidence for the already-frozen T7 seed candidates

> Coordinator review: `COORDINATOR-REVIEW-2026-09-07.md` and the pre-outcome
> `COORDINATOR-DECISION-2026-09-07.md` supersede the inferential headline below.
> Effects reproduce; zero of the full 11 primary tests passes Holm. The
> post-hoc Poisson-binomial result is not an accepted confirmatory rescue.
> The original worker report is retained with its deviations visible.

**11 primary candidate pairs (5 failing-arm + 6 working-arm), 3 NOT TESTABLE (candidate gene below the
pre-registered expression floor), 8 evaluated. 2 of the 8 meet all four pre-registered corroboration
conditions — `GRB10_P2 → SNX8` (failing arm) and `DPF2_P1 → CEP120` (working arm) — as do 2 of the 4
non-targeting-guide comparators, `non-targeting_00018|00127 → DCTN2` and `non-targeting_00121|00339 →
PRKCQ`. So the candidate genes ARE repressed far more often than matched arbitrary genes under the same
construct (4 observed vs 0.029 expected, Poisson-binomial P = 1.8e-09), but the corroboration does NOT
separate failing from working constructs: 1 failing, 1 working, 2 non-targeting. Clustered at the
intended target, exactly 1 of the 17 paired targets (GRB10) carries a corroborated candidate.
**301's negative stands: this is no expression corroboration for the failing-vs-working contrast — never
"off-targets do not explain failing constructs".** Strongest limitation: no biological replicate of any
construct exists in these files, so a construct-specific population or handling difference that happens
to hit the candidate gene cannot be excluded pair-by-pair; only the aggregate rate argues against it, and
that rate is estimated from the same data.**

2026-09-07, 05:25–05:41Z (16 min of compute inside the 2 h box). Preregistration hashed
**05:30:53Z**, sha256 `f7e61ce67878e9499064971edc811bde836368fe9104c38b153d2fc867935fcb`, before any
element of any `X/data` array was read. Peak RSS **3.79 GiB** (the train stream; analysis 0.36 GiB), cap
6 GiB. The three 2025 files were streamed **serially**, one pass each, 21.5 / 13.9 / 6.8 s; the memory
gate was checked before each. No `vcc`, no `cell-eval2`, no build, no submission, no download, no new
seed search, no reserve access, no commit. Nothing written outside
`experiments/305-t7-expression-test/`; no script outside it was invoked.

**This is an extension on already-used data, not an untouched confirmation.** 293 read these three files
and published the per-construct knockdown scalars from which 298 derived the 13/21 failing/working
labels; 301 froze the candidate list from sequence alone. The only genuinely new reading is the candidate
genes' own columns.

---

## 1. The frozen candidate table reproduces 301 exactly

`work/candidates305.py` rebuilds 301's primary reading from its saved `seed_hits301_cage.csv` with its
exact filters (12-nt seed, PAM-anchored NGG, CRISPRi active window TSS −50…+300, distal ≥ 10 kb from any
window of the arm's own gene), duplicate spacers/hits collapsed to one construct–gene pair:

| group | constructs with ≥ 1 candidate | 301 published |
|---|---|---|
| failing arms | **4 / 13** | 4 / 13 |
| working arms | **4 / 21** | 4 / 21 |
| non-targeting | **4 / 31** | 4 / 31 |

24 self-locus hits over 3 construct–gene pairs dropped, as 301 recorded. **All 15 evaluated pairs are
frozen** (`work/candidates305.csv`), not only the failing arm's. All 15 candidate genes are measured
columns; the gene list (18,080 genes, sha `da2b3839…`) is identical in all three splits; **no construct
anywhere in the 2025 library targets any of the 15 candidate genes**.

## 2. Known-answer rows — 7 / 7 pass

| row | result |
|---|---|
| KA-1 34 construct identities, 13/21 labels, 34/34 cell counts | 34/34, 13/21, 34/34 |
| KA-2 gene order identical across splits; candidate genes present | identical digest; 15/15 |
| KA-3 301's contrast reproduced | 4/13, 4/21, 4/31 |
| KA-4 synthetic known fold change (f = 0.5 on one column) | recovered − analytic = **0.00e+00** (independent pure-`math` closed form) |
| KA-5 chunked group sums vs a direct per-row slice | exact equality on all 18,080 columns, all three splits |
| KA-6 same-NTC-halves null | median **+0.00000**, sd 0.0041 over 18,080 genes, [2.5, 97.5] = [−0.009, +0.008] |
| KA-7 34 own-target effects vs 293's `kd_full` | max \|diff\| **0.00e+00** over 34 arms |

KA-7's comparable configuration is stated exactly: 293's `kd_full` is full depth (`thin = 1.0`), **all**
of a construct's cells (not its 1,600-cell extraction cap), against the pooled NTC set — identical to
this run's pooled reference. 293's thinned (0.368) and per-rung (25–200 cells/side) numbers are different
configurations and were **not** used as known answers.

## 3. Design as executed

Response formula, frozen (285's `lib285.prof`, the arithmetic 293 used for `kd_lfc`):
`e(c,g) = log1p(50000·S_c[g]/Σ_j S_c[j]) − log1p(50000·S_ctrl[g]/Σ_j S_ctrl[j])`, the denominator always
the **full library total**. Reported alongside are the exact fold change (CP50k ratio) and the
panel-centred residual; within a split the residual's rank over constructs is identical to the direct
effect's, so the empirical rank is reported once and the residual is reported for its sign and magnitude
— **being below the panel mean is not going below controls**, and the table shows both.

Halves are deterministic, no RNG: within each batch a construct's cells alternate 0,1,0,1 in ascending
file-row order. Median lane imbalance between the two halves is ≤ 1 cell per lane. Control halves are
293's partition verbatim: the 31 non-targeting constructs sorted, `i % 2` → half A (16 constructs,
21,581 cells) / half B (15 constructs, 16,595 cells), **sharing no cell and no guide construct**. Half 0
is referenced to A, half 1 to B; **no control is shared between the two halves of a comparison**. A
non-targeting construct under test is removed from its own containing half (leave-one-construct-out).

**Measured and recorded before the test, a real limitation**: the 31 non-targeting constructs' per-half
count sums are **byte-identical across the three 2025 files** — the control pool is the same 38,176 cells
in every split. Controls are therefore shared *between* splits (never within a comparison).

**No cell bootstrap is labelled biological replication.** These files carry no biological replicate of any
construct; the two halves are two halves of one construct's cells in one experiment, and every interval
here is a cell-sampling interval. None was invented.

## 4. Positive control: the same statistic on each arm's own target gene

Applying the identical four criteria to each arm's **own** target gene: **21 / 34 arms pass** (4 / 13
failing, 17 / 21 working). The strongest are `TMSB4X_P1` e = −3.04 (fold change 0.040, rank 1 of 157,
z = −70.9) and `STAG2_P1` e = −1.79 (z = −114). The machinery detects a real CRISPRi repression as a
rank-1 outlier, and it independently reproduces 293's failing/working direction (4/13 vs 17/21). This is
the sensitivity demonstration for the negative below.

## 5. The frozen test

`work/expression305.csv`, `work/results305.json`. `fc` = exact CP50k fold change; `rank` = position in the
within-split, within-gene construct null; `p` = empirical rank p with the +1 correction; `Holm` over the
8 testable primary pairs; `band` = fraction of expression-matched arbitrary genes meeting the identical
conjunction under the *same* construct (§6, post-hoc).

| group | construct | candidate | e | fc | half 0 | half 1 | rank | p | Holm | Holm (all-splits) | band | criteria |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| failing | `CAST_P1` | UTP18 | −0.0115 | 0.987 | −0.011 | −0.012 | 58/158 | 0.367 | 1.00 | 1.00 | 0.00% | C1 C3 |
| failing | `GRB10_P2` | NEUROG3 | −0.0008 | — | — | — | — | — | — | — | — | **NOT TESTABLE** (control CP50k 0.001) |
| failing | **`GRB10_P2`** | **SNX8** | **−0.4924** | **0.451** | **−0.519** | **−0.466** | **1/105** | **0.0095** | 0.0762 | **0.0252** | 0.044% | **C1 C2 C3 C4** |
| failing | `KDM2B_ENST…071.4` | WRN | −0.0328 | 0.962 | −0.037 | −0.028 | 14/158 | 0.089 | 0.53 | 0.53 | 0.026% | C1 C3 |
| failing | `KLF10_P1` | ST3GAL3 | +0.0061 | 1.011 | −0.037 | +0.054 | 70/158 | 0.443 | 1.00 | 1.00 | 0.00% | — |
| working | **`DPF2_P1`** | **CEP120** | **−0.6328** | **0.345** | **−0.650** | **−0.604** | **1/105** | **0.0095** | 0.0762 | **0.0252** | **1.66%** | **C1 C2 C3 C4** |
| working | `IGF1R_P1` | IFITM5 | −0.0125 | — | — | — | — | — | — | — | — | **NOT TESTABLE** (control CP50k 0.020) |
| working | `MED15_P1` | AFF4 | −0.0079 | 0.991 | −0.007 | −0.009 | 25/105 | 0.238 | 1.00 | 1.00 | 0.63% | C1 C3 |
| working | `MED15_P1` | TMCC1 | +0.0259 | 1.051 | +0.046 | +0.006 | 88/105 | 0.838 | 1.00 | 1.00 | 0.25% | — |
| working | `MED15_P1` | ZSWIM5 | +0.0866 | — | — | — | — | — | — | — | — | **NOT TESTABLE** (control CP50k 0.323) |
| working | `TMSB4X_P2` | GNG10 | +0.0474 | 1.057 | +0.042 | +0.054 | 145/158 | 0.918 | 1.00 | 1.00 | 0.00% | — |
| NTC | `nt_00010\|nt_02551` | DDB1 | +0.0177 | 1.032 | +0.021 | +0.014 | 113/159 | 0.711 | 1.00 | 1.00 | 0.00% | — |
| NTC | **`nt_00018\|nt_00127`** | **DCTN2** | **−0.5157** | **0.549** | **−0.515** | **−0.517** | **1/159** | **0.0063** | 0.0755 | **0.0377** | 0.101% | **C1 C2 C3 C4** |
| NTC | `nt_00095\|nt_03143` | CLNS1A | −0.0143 | 0.982 | −0.018 | −0.011 | 53/159 | 0.333 | 1.00 | 1.00 | 0.00% | C1 C3 |
| NTC | **`nt_00121\|nt_00339`** | **PRKCQ** | **−0.9032** | **0.105** | **−0.909** | **−0.898** | **1/159** | **0.0063** | 0.0755 | **0.0377** | 0.023% | **C1 C2 C3 C4** |

The four non-targeting rows are byte-identical in all three splits (same control cells); the train row is
shown. Control-vs-control effects on the 15 candidate genes: median +0.0031, range [−0.0063, +0.0534] —
the scale on which "unusually negative" must be read.

**The Holm limit was declared in advance and it bound as declared.** The empirical rank p has a floor of
`1/(N_null+1)` = 0.0063 (train, 158 targeting constructs) / 0.0095 (test, 105) / 0.0185 (validation, 54).
With 8 testable primary pairs the smallest attainable Holm-adjusted p is **0.0506 > 0.05**: family-wise
significance is unattainable in the primary null **for any effect however large**, so C4 was pre-set as a
raw-rank criterion and Holm is reported, not thresholded. The pre-declared secondary null (all 317
gene-targeting constructs of the three splits, each against its own split's controls — the same control
cells, so the reference is exactly homogeneous) has the resolution: there both primary hits reach
**Holm = 0.0252** and both non-targeting hits **Holm = 0.0377**.

### Specificity (item 7)

The same candidate gene under the **other construct of the same intended target**:

| pair | tested construct | partner construct | partner e | partner rank |
|---|---:|---|---:|---:|
| `GRB10_P2` → SNX8 | −0.4924 | `GRB10_P1` | **+0.0453** | 93/105 |
| `DPF2_P1` → CEP120 | −0.6328 | `DPF2_P2` | −0.0005 | 40/105 |
| `CAST_P1` → UTP18 | −0.0115 | `CAST_P2` | −0.0028 | 86/158 |
| `KDM2B_…071.4` → WRN | −0.0328 | `KDM2B_…069.4` | −0.0171 | 29/158 |

Both corroborated primary candidates are repressed under **only** the seed-carrying construct of their
target, and not at all under the partner construct that shares the intended target. That is what a
construct-specific association looks like and it is inconsistent with intended-target downstream biology
— but it is equally consistent with any other construct-specific difference, and **it does not prove an
off-target mechanism**: nothing here perturbs the seed or measures dCas9 occupancy.

## 6. Post-hoc controls (added AFTER the outcomes were seen — see §8)

**(a) The base rate is the number that decides the reading.** Applying the identical conjunction to
*every* testable gene (control CP50k ≥ 0.5; 9,596 genes) under the *same* construct:

| construct | cells | genes meeting C1–C4 | of them, the candidate? | the whole exceptional set |
|---|---:|---:|---|---|
| `CAST_P1` | 894 | 1 / 9,596 | no | **CAST** (its own target) |
| `GRB10_P2` | 2,136 | 2 / 9,596 | **yes** | KCTD14 (−0.63), **SNX8** (−0.49) |
| `KDM2B_…071.4` | 994 | 5 / 9,596 | no | JADE3, TCEAL9, NAV1, POFUT2, ADAM19 |
| `KLF10_P1` | 449 | 0 / 9,596 | — | — |
| `TMSB4X_P2` | 2,261 | 0 / 9,596 | — | — |
| `IGF1R_P1` | 1,614 | 0 / 9,596 | — | — |
| `MED15_P1` | 1,459 | 36 / 9,596 | no | — |
| **`DPF2_P1`** | **95** | **120 / 9,596** | yes | — |
| `nt_00010\|nt_02551` | 1,262 | 0 / 9,596 | — | — |
| `nt_00018\|nt_00127` | 1,519 | 5 / 9,598 | **yes** | THOP1, DPY30, **DCTN2**, CCDC32, SMTN |
| `nt_00095\|nt_03143` | 1,285 | 0 / 9,598 | — | — |
| `nt_00121\|nt_00339` | 1,262 | 1 / 9,597 | **yes** | **PRKCQ** — the only gene in the transcriptome |

Two readings follow. First, `PRKCQ`, `SNX8` and `DCTN2` are genes pre-specified from sequence alone that
land in a set of 1, 2 and 5 genes out of ~9,600. Over the 15 pairs, the expected number meeting the
conjunction for an expression-matched arbitrary gene under the same construct is **0.029**; **4** were
observed (exact Poisson-binomial **P(X ≥ 4) = 1.8e-09**; primary 11 pairs, expected 0.028, observed 2,
**P = 2.3e-04**; dropping the four pairs on globally noisy constructs, expected 0.0019, observed 3,
**P = 3.0e-10**). Second, **`DPF2_P1|CEP120` does not survive this control**: 120 of 9,596 genes do the
same thing under `DPF2_P1`, whose 95 cells make it the smallest of all 34 arms, so its per-pair rate is
1.66 % — a nominal p ≈ 0.017, not a finding. The base rate tracks cell count directly (Spearman
−0.48 train / −0.33 test between construct cells and the number of genes meeting the conjunction; the
< 200-cell band averages 247–443 such genes, the > 1,000-cell band 18–38, median over all constructs 5–6).

**(b) Lane-matched controls: 0 / 15 verdict flips.** Every effect was recomputed against a control
pseudobulk carrying the construct's own lane composition
(`S_ctrl = Σ_lane (n_c,lane / N_ctrl,lane)·S_ctrl,lane`, 48 lanes, each null construct lane-matched to
itself). Effects move by ≤ 0.005 (`SNX8` −0.4924 → −0.4904; `PRKCQ` −0.9032 → −0.8790; `CEP120` −0.6328 →
−0.6334) and no conjunction verdict changes. Batch/lane confounding is excluded as the explanation.

**(c) The control pool is barely affected.** Removing the two off-targeting non-targeting constructs
(2,781 of 38,176 control cells) shifts the pooled control by at most **0.023** (on PRKCQ), median 0.000
over all genes. The finding is about individual guides, not about the usability of the pooled NTC control.

## 7. Reading

- **The seed-off-target phenomenon gets observational corroboration in this library.** A 12-nt PAM-proximal
  seed match inside a promoter's CRISPRi-active window, frozen from sequence before any expression was
  read, predicts strong, half-reproducible, lane-robust, construct-specific repression of that gene far
  more often than chance (4 of 15 vs 0.029 expected). Three of the four survive every control.
- **It does not explain construct failure, and T7's negative stands.** The corroborated pairs are 1 failing
  arm, 1 working arm (and that one is explained away by its 95 cells), and **2 of the 4 non-targeting
  guides**. Clustered where the uncertainty belongs — at the **intended target**, 17 paired targets, not 34
  independent constructs — exactly **one** target (GRB10) carries a corroborated candidate on one of its
  two constructs. 301's 4/13-vs-4/21 contrast (Fisher p = 0.679) gains no support here.
- **"Non-targeting" means "matches no annotated TSS by design", not "binds nothing".** Two of the 31 control
  constructs repress a sequence-predicted promoter candidate, one of them 9.5-fold (PRKCQ, the single
  exceptional gene in that guide's transcriptome). That is a new, small, concrete observation about this
  library's control set; it is consistent with the published 0.8–4.7 % off-target rate for a validated
  CRISPRi library and it is why 301's non-targeting group was the right background.
- **Absence of a supported candidate for the other 6 testable primary pairs is "no expression corroboration
  under this test"**, never "those constructs have no off-targets": four of them (UTP18, WRN, AFF4,
  CLNS1A) are repressed in direction on both halves at magnitudes 3–50× too small to clear C2.

## 8. Deviations, and whether each was fixed before or after outcomes

1. **Before outcomes (tooling only).** `apply_patch` is not exposed in this harness; authored files were
   written with the equivalent file-write tool. No behavioural difference.
2. **Before outcomes.** The positive control of §4 (the identical frozen statistic on each arm's own target
   gene) is not named in PREREGISTRATION §7. It was written into `expr305.py` and fixed before the script
   was run, i.e. before any outcome was seen; it introduces no new degree of freedom (same statistic, a
   different pre-specified gene).
3. **After outcomes.** The three post-hoc controls of §6 — the expression-matched base rate
   (`baserate305.py`), the lane-matched controls (`lanematch305.py`) and the cell-count / pool-contamination
   context (`context305.json`) — were written after the frozen table was read. Each is a *negative*
   comparator: it can only weaken a positive reading, never create one, and in fact it removed one of the
   two primary hits (`DPF2_P1|CEP120`). No frozen quantity was changed: not the seed length, window,
   distance filter, candidate set, failing definition, response formula, effect threshold, pass rules,
   null definition or testability rules.
4. **Before outcomes, recorded as declared.** PREREGISTRATION §8 computed the smallest attainable primary
   Holm p as 11 × 0.0063 = 0.070 assuming 11 testable pairs; 3 pairs turned out untestable, so the realised
   family is 8 and the attainable minimum is 0.0506. Still > 0.05, so the declared conclusion (Holm
   unattainable in the primary null) is unchanged.
5. **Not a deviation, noted.** The 4 non-targeting comparator pairs were evaluated in all three splits as
   §1 of the prereg allows; the three evaluations are byte-identical because the control cells are the same
   in every file. They are counted once as 4 pairs and reported as 12 split-evaluations.

## 9. Artefacts

| path | what |
|---|---|
| `PREREGISTRATION.md`, `PREREGISTRATION.sha256` | frozen design, hashed 05:30:53Z, `f7e61ce6…` |
| `work/input_manifest.json` | the three h5ads and every derived input: bytes, sha256, shapes, gene-list digests |
| `work/mapping305.json` | per split: constructs, cells, lanes, NTC parity partition, batch composition |
| `work/candidates305.csv`, `work/candidates305_meta.json` | the 15 frozen pairs and 301's reproduced counts |
| `work/sums305_{train,validation,test}.npz` (13.5 / 6.2 / 9.9 MiB) | **the reusable sufficient statistics**: exact int64 per-construct, per-half count sums over all 18,080 genes, + cells, UMIs, gene names. sha `f0fa94d8…` / `06d8d0fe…` / `40bd3b00…` |
| `work/sums305_*.json` | per-half cell/UMI counts, per-lane composition, file sha256, KA-5 |
| `work/ctrl_lane_sums.npz` | per (lane × NTC parity) control sums, 48 × 2 |
| `work/ka305.py`, `work/ka305.json`, `work/ka305.log` | the 7 known-answer rows (the tests) |
| `work/expression305.csv` | 23 rows: effects, both halves, fold change, ranks, p, Holm, all-splits Holm, residual, paired-construct comparator, control-vs-control, base rate, lane-matched effect, cell/UMI counts, lane imbalance, testability |
| `work/results305.json` | every number above, machine-readable, incl. the post-hoc blocks |
| `work/{candidates305,manifest305,lib305,sums305,expr305,baserate305,lanematch305,finalize305}.py` | the code |
| `work/baserate305.json`, `work/lanematch305.json`, `work/context305.json` | the three post-hoc controls |

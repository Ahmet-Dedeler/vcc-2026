# 353 — Recipient-conditioning scout: three basal-state experiments that 325–341 did not run

**Status: COMPLETE** (literature + repo synthesis; proposal memo, no fitted arms).

## Known-answer reproduction (own code, `ka353.py`)

| quantity | reproduced | recorded | match |
|---|---|---|---|
| mean NMAE, cap15, 58 eligible targets | 0.9508364630325662 | 323 `metrics-cap15.json` raw nmae 0.9508364630325662 | **exact (16 digits)** |
| macro \|lfc\|-weighted delivered sign accuracy | 0.5879679616799467 | 343 `signs-cap15.json` `w_all` 0.5879679616799467 | **exact** |
| eligible targets / mean gate | 58 / 1119.5 | 343: 58 / 1120 | match |
| micro-averaged (pooled) weighted sign | 0.5619 | not recorded (new) | descriptive |

Note for consumers: 343's 0.588 is a **macro** mean over targets; the pooled
micro average is 0.562. Both are development statistics on the reused H1 panel,
not A/B/C forecasts. Gate: reference `p_adj<0.05`, finite nonzero lfc, gene ≠
target, ≥10 genes/target. `reference.npz` opened mmap'd (only `genes`/`keep`
touched).

## What was measured

Read-only synthesis. Local: RESEARCH-PROGRAMME-2026-09, SPRINT-HANDOFF-2026-09-09
(325–342 ledger), 343, 344, 246 §0–§1, 315 §0, 314 §4–§5, `core341.py`,
334 `selection.json`. Web (all accessed 2026-09-10): STATE, TxPert, GEARS, CPA,
scLAMBDA, scGPT/scFoundation benchmark, Tahoe-100M/Tx1/TX1-CD, X-Cell, scDEBART,
AROMA, TwinCell, LGR, IRIS, PerturbNet, GPerturb, PS-score, CRISPR-GEM. No
process exceeded the tiny-check budget; no bank was fitted.

## The gap this memo serves

343: delivered sign accuracy 0.588 macro (0.562 micro) on the reference's
significant genes; reference-**up** genes are a coin flip (0.516); amplitude is
exhausted (best scalar → 0.928; 227/317). 341 proved sign accuracy *can* move
(+0.033) but only with the recipient's own perturbation **labels** — a regime we
will never have in October. The untested middle: the recipient's **label-free
basal state** (18,400 control cells per A/B/C context, currently used only for
means, CPM gate, noise proxy) as a conditioner of *which* genes respond and in
which direction.

## Literature survey: what each model consumes about the target line

| model (source) | recipient information consumed | transfer justification / validation |
|---|---|---|
| STATE/ST+SE (bioRxiv 2025.06.26.661135; Cell 2026; https://github.com/ArcInstitute/state) | the full **set of control cells** from the target context + SE foundation embedding (167M observational cells) | 5 held-out Tahoe lines vs context-mean and perturbation-mean baselines; SE explicitly credited for no-perturbation-data contexts. MMD loss over cell sets, self-attention over the control set |
| TxPert (arXiv 2505.14919; https://github.com/valence-labs/TxPert) | batch-matched control profile via MLP basal encoder (or raw controls, or scVI embeddings) + GNN perturbation embedding (STRING/GO/PxMap/TxMap) | 4-line leave-one-line-out beats general baseline all 4 lines. **Key ablation**: for unseen lines, *no* basal encoder beats learned MLP ("by far"); MLP only helps on scVI embeddings. Lesson: with ≤4-line diversity, learned recipient encoders fail |
| GEARS (Nat Biotech 2023, s41587-023-01905-6) | **none from a new line** — trained per dataset (K562, RPE1 separately); prediction added onto a sampled control cell; co-expression graph from *training* line | unseen-*gene* transfer within line only. No cross-line mechanism exists in GEARS |
| CPA (Mol Sys Biol 2023, msb.202211517) | the cell's own expression, adversarially disentangled into basal state; counterfactual by swapping perturbation embedding | combos/doses within observed lines. Nonlinear decoder "enables cell-type-specific effects" — claimed, not held-line validated |
| scLAMBDA (bioRxiv 2024.12.04.626878) | basal latent z (VAE) disentangled from GenePT-embedding perturbation code via MINE | unseen genes/combos within line. No cross-line |
| scGPT/scFoundation (BMC Genomics 2025, s12864-025-11600-2) | unperturbed cells + perturbation token | **Negative**: train-mean baseline beat both; GO-feature random forest beat scGPT "by a large margin". Perturbation representation > basal encoding; benchmarks low-variance |
| Tahoe-100M/Tx1/TX1-CD (bioRxiv 2025.02.20.639398; 2025.10.23.683759; https://blog.tahoebio.ai/p/mapping-drug-response-cascades-with) | DMSO-cell basal decoder (frozen) + TX1 embeddings that "smooth over assay/batch noise"; mutation status (BRAF/KRAS-stratified response shown) | Tx1+ST matches direct expression, beats SE-600M on held-out contexts; few-shot + zero-shot context tasks. Modality is chemical (see rejection 4) |
| X-Cell/Pisces (bioRxiv 2026.03.18.712807) | control→perturbed transition refined by cross-attention to multimodal priors (text, protein LM, networks, DepMap, morphology) | 5× PearsonΔ; zero-shot to melanocyte progenitors + primary CD4s; power-law scaling; test-time adaptation helps |
| scDEBART (ICML 2026 poster 61521) | **basal expression as the conditioning input**, pretrained to predict logFC (scVI-denoised) at scale: 6.28M change profiles / 66.6M cells | EF 11.96, 4–7× scGPT/GEARS; dose-dependent cross-modal (drug) transfer. Direct precedent that basal-conditioned change pretraining beats static co-expression pretraining |
| AROMA (arXiv 2604.20263) | per-(perturbation, target, **cell line**) triple classification (up/down/non-DE) with retrieval-grounded multimodal evidence | zero-shot unseen line 0.77→0.73 macro-F1, still beating baselines. Closest published endpoint to our sign problem |
| DYNAMO (ICML 2026 poster 70802) | GRN-topology operator only; **no recipient data at all** | K562→RPE1 ρΔ 0.576 (51% retained vs 30% unstructured). Structured propagation transfers; GRN topology +19% over text/co-expression/FM embeddings |
| TwinCell (bioRxiv 2026.01.29.702072) | foundation-model embeddings as context for interactome-constrained path decomposition (inverted task: states→regulators) | zero-shot lines+perturbations winner on TwinBench; reports STATE "on par with random HVG" on its metric — cautionary about benchmark-dependence |
| LGR-Tahoe (arXiv 2608.01734) | cell-line-restricted retrieval pool (observed deltas in the *target* line) selected by LLM | sign-accuracy gains >20pp over cell mean; "retrieval quality rather than predictor complexity is key" |
| IRIS (Nat Methods 2026, s41592-026-03213-8) | control expression projected into conserved signaling-response signature space (CVAE trained on hESC perturbation atlas) | cross-cell-type, cross-species signature transfer; endoderm→mesoderm generalisation. Supports a conserved-component view |
| PS-score (Nat Cell Biol 2025, s41556-025-01626-9) | **within-line** cell state (TF activities, state, cofactor context) explains heterogeneous response to the *same* perturbation | response strength varies with cell state inside one line — the fine-grained analogue of our cross-line question |
| CRISPR-GEM (PMC11244939) | basal expression explicitly: CRISPRa modeled inversely, CRISPRi positively correlated with initial expression | weak-mechanism precedent that basal level gates response magnitude/direction |
| GPerturb (Nat Commun 2025, s41467-025-61165-7) | basal component (GP over cell parameters) + sparse perturbation component, additively separated | directionality disagrees across methods/preprocessing — cautionary on sign fragility |

Two structural lessons span the table. (1) **Nobody with ≤6 lines has made a
*learned* recipient encoder work** — TxPert's ablation is explicit, and our 332
is its in-house replication (TabICL+recipient retrieval .676 vs source-mean
.733). Any proposal that fits recipient parameters on 6 lines is pre-refuted.
(2) What *does* transfer across lines is **structure with the recipient entering
parameter-free**: DYNAMO's frozen GRN operator (51%), IRIS's conserved
signatures, LGR's target-line retrieval pool, scDEBART's basal-conditioned
(rather than basal-fitted) change model.

## Proposed experiments (ranked)

### 1. Expression-gated network propagation — zero fitted recipient parameters

**Mechanism.** A knockdown's effect propagates through the regulatory network,
but an edge can only carry response in a line where its regulator/pathway is
active — readable label-free from the recipient's basal expression (e.g. a
repressor that is unexpressed in Jurkat cannot mediate down-regulation there,
which is exactly the kind of asymmetry behind 343's 0.516 up-gene coin flip).
Donor lines fit *how strongly* each edge/pathway propagates; the held-out line
contributes only a binary/continuous activity gate over the same edges. The
recipient thus conditions *which paths exist* without contributing a single
fitted degree of freedom, dodging the TxPert/332 failure mode by construction.

**Estimator.** LOLO over the six-CRISPRi-line frame (`experiments/334-text-conditioned-transfer/selection.json`:
lines k562_gwps / orion_hct116 / orion_hek293t / rpe1 / nadig_jurkat /
nadig_hepg2; 1,024 shared targets; 512 coordinates). (a) Fit per-edge (or
per-pathway) propagation strengths by ridge on 5 donors' responses
(`experiments/341-known-line-calibration/core341.py` `fit`/`predict` as the
regression harness). Graph: STRING + CollecTRI/DoRothEA prior edges (the same
priors TX1-CD and 334's STRING control draw on). (b) Gate the fitted graph by
the held-out line's basal expression (basal pseudobulks already on disk per
donor bank; challenge-side analogue: control pseudobulks / CCLE profiles) and
predict response signs by propagation from the knocked-down target. Compare
against the ungated (donor-mean-topology) propagation and the per-gene diagonal.

**Primary endpoint.** \|lfc\|-weighted sign accuracy on the held-out line's
reference-significant genes (macro over targets, 343's `w_all` functional);
retrieval reported as guard only, never selected on.

**Cheap falsifier (<1 h).** Wrong-recipient swap: gate with each of the other
5 lines' basal profiles. **Kill if** own-basal-gated minus best-wrong-basal
Δq_w ≤ +0.005, or if gated ≤ ungated diagonal. (Scale: 232's one-sided rate
+0.01 q_w ≈ −0.0026 raw nmae; +0.005 is half a meaningful nmae step and 1/6 of
341's labeled +0.033.)

**Why 325–341 do not answer it.** 340 fitted *response-gene* couplings
(undirected regression on responses; retrieval −0.0052 vs diagonal, all 6
lines) — no basal gating, no directed propagation, no wrong-recipient control.
334 used STRING only as a control for a *text* model. 339's gene-wise marginal
maps had no network structure at all (+0.0029, below floor). GEARS itself never
tested cross-line. Nothing has put the recipient's basal profile inside the
topology.

### 2. Supervised per-(target, gene) sign classification with recipient-basal + prior features

**Mechanism.** Sign is a discrete property and may be predictable without
modeling magnitudes at all — AROMA's triple-classification result (0.73
zero-shot macro-F1) is the existence proof. Features per (t, g, L): donor
consensus sign and agreement, basal expression of t and g in L, network
distance t→g, shared-pathway indicators, orthologous DE-direction priors. A
shallow classifier (logistic/GBM) learns how these combine; magnitudes never
enter, so 228's magnitude-`nmae` coupling cannot fire.

**Estimator.** Same six-line LOLO frame and `core341.py` interval machinery;
labels are donor lines' observed response signs (no recipient labels ever —
unlike 341). Fit on 5 lines' (t, g) pairs, predict held-out line signs.
Banks: the six CRISPRi donor tables behind `selection.json`; STRING/GO priors
as in TxPert/AROMA; basal features from each line's control pseudobulk.

**Primary endpoint.** Same \|lfc\|-weighted sign accuracy + retrieval guard.

**Cheap falsifier (<1 h).** Two controls: shuffled-recipient basal features,
and the diagonal-sign baseline. **Kill if** cross-line Δq_w vs diagonal ≤
+0.005, or if shuffled-recipient is within noise of own-recipient (recipient
features carry nothing beyond donor consensus).

**Why 325–341 do not answer it.** 341 is the lookalike and the difference is
the whole point: it used the recipient's own perturbation *labels*
(known-line regime) and still failed retrieval guards — it never tested
label-free basal features cross-line. 269 routed among frozen *actions* per
target (0/3 settings), not per-(t,g) sign classification — a different object
with different features. 332's TabICL predicted response *values* with a
generic tabular prior and was judged on retrieval, not sign accuracy.

### 3. Genotype-stratified response priors (DepMap/CCLE mutation status)

**Mechanism.** Tahoe-100M's BRAF/KRAS stratification shows driver mutations
rewire response direction along pathways; the recipient's mutation profile
(which we hold for A/B/C-grade lines via CCLE/DepMap, both on disk:
`data/external/ccle`, `data/external/depmap26q1`) selects *which donor's sign
to trust* per pathway — a discrete, parameter-free conditioner like Exp. 1 but
over genotypes instead of expression.

**Estimator.** Per-pathway donor sign consensus stratified by mutation status
of upstream drivers (KRAS/BRAF/TP53/CDKN2A panel from the Tahoe precedent);
predict held-out line signs with the mutation-matched stratum. Same LOLO
frame; mutation labels from CCLE/DepMap metadata (no fitting on responses).

**Primary endpoint.** Same weighted sign accuracy + retrieval guard.

**Cheap falsifier (<1 h).** Permuted mutation labels across lines. **Kill if**
mutation-matched ≤ unstratified consensus, or if permutation matches the
real stratification (labels carry no information at 6 lines — the likely
outcome, stated in advance).

**Why 325–341 do not answer it.** No ledger experiment used genotype metadata
at all; the closest (289) used lineage/proliferation/p53 *expression*
programmes (0/25 cells transcriptome-wide) — mutation status as a discrete
stratum selector is untested. Ranked third because at 6 lines the strata are
thin (mostly 1–2 lines each), so the falsifier is the expected outcome; it is
cheap enough to earn its place.

## Considered and rejected (with reason)

- **Learned basal encoder (MLP/surrogate on control profiles):** pre-refuted —
  TxPert's own ablation (no-encoder best at 4-line diversity) + our 332
  (recipient-conditioned TabICL retrieval −0.057 vs source mean).
- **STATE checkpoint fine-tuning:** closed by 182 (cannot be conditioned on our
  targets); and 314/315 show the available sequence/chromatin priors are
  cell-type-blind (PromoterAI) or undefined for held-out lines (ChromBPNet).
- **Text/semantic features (GenePT, TF-IDF, PubMedBERT):** closed by 333
  (vectors quarantined: exact collisions on 34 challenge targets), 334, 336
  (both lose retrieval incl. vs shuffled text).
- **Gene-wise marginal/basal→response maps:** closed by 339 (+0.0029, below
  floor, diagonal unbeaten) — Exp. 1–2 differ by adding network structure
  (Exp. 1) or supervised sign combination (Exp. 2).
- **Response-neighborhood coupling:** closed by 340 (retrieval-negative all 6
  lines) — undirected response regression, not basal-gated propagation.
- **Own-line perturbation labels:** 341's regime — informative (+0.033 sign)
  but unbuildable for unseen contexts; Exps. 1–3 are strictly label-free.
- **Uniform additive down-shift:** closed by 155 (reach +0.0047 for pds
  −0.0103, 4:1 against); Exps. 1–3 predict per-gene signs, not a global shift.
- **CPM-gate moves / basal gating beyond 10 CPM:** closed by 246-H3
  (088/164/213/006); Exps. 1 conditions topology, not call-set membership.
- **Sign-consistent consensus / majority vote:** closed by 246-H4
  (232 `maj_strict` q_w −0.18; 234 pool-blind flips all-negative).
- **Empirical-null / control-DE keys:** dead by 344 (NTC constructs carry no
  DE) — Exps. 1–2 read control *levels*, not control *DE structure*.
- **Per-CPM-stratum amplitude:** oracle ceiling +0.009 raw (227) — amplitude
  is not the gap; all three experiments target signs only.
- **Tahoe drug-response transfer:** closed by 184 (cross-line adaptation) +
  chemical-vs-CRISPRi modality gap; only Tahoe's *genotype-stratification
  idea* is reused (Exp. 3), not its response values.
- **Full basal-conditioned pretraining (scDEBART-style):** infeasible in this
  programme's budget and an engineering project, not a 1-hour-falsifiable
  experiment; Exps. 1–2 are its cheap linear/supervised shadows.
- **Control-cell covariance / co-expression subspaces:** closed by 246-H6b
  (155 null split-half r=0.028; 160; 174) as operator/subspace objects; Exp. 1
  uses basal *levels* as gates, a different functional of the controls.
- **scGPT/scFoundation-style FM fine-tuning:** benchmark-negative (train-mean
  wins; GO-RF wins big) and perturbation-blind pretraining; no mechanism for
  our sign gap.

## Verdict against the PLAN.md kill criteria

- Known-answer gate: **PASS** (both numbers bit-exact, §Known-answer).
- Literature kill check ("recipient conditioning needs ≥10 lines, untestable
  at 6"): **not triggered** — all three proposals are fittable/testable in the
  six-line LOLO frame with parameter-free recipient entry, and each carries a
  sub-hour falsifier with a numeric kill threshold.
- CLOSED-list collision check: **PASS** — every rejected idea cites its
  closing experiment above; the three proposals each name the lookalike
  experiment and the exact difference (340 vs gated-propagation; 341-labels vs
  341-basal-features; 289-expression-programmes vs mutation strata).

## Limits

H1-panel development reasoning only; no A/B/C forecast is made. No retrieval,
cosine or proxy gain is claimed as a score gain. LOLO sign-accuracy deltas are
conditional on the six-line frame's correlated targets and study/line
confounding (per handoff §ledger). The ≤15-word web-quote rule was observed by
paraphrase; URLs + access dates are in the survey table. Genotype proposal
(Exp. 3) is a priori unlikely at 6 lines — included only because its falsifier
costs minutes.

## Exact commands to reproduce

```bash
# known-answer (writes nothing; reads 323 work files only)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/353-recipient-conditioning-scout/ka353.py
# expected: n_targets 58, mean_gate ~1119.5,
#   weighted_sign_all(micro) 0.5619153131793206,
#   macro 0.5879679616799467 == 343 signs-cap15.json w_all,
#   mean_nmae 0.9508364630325662 == 323 metrics-cap15.json raw nmae
```

## Peak memory and wall time

`ka353.py`: <60 s wall; inputs 34 MB (`stats-cap15.npz`) + 18 MB
(`reference_de.parquet`), `reference.npz` mmap'd (small-key touch only) —
well under 1 GiB RSS, under the 2 GiB / 5 min tiny-check budget. No other runs.
Memory gate at start: 75% free.

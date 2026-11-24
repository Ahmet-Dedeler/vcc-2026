# 813 PROPOSAL — TeloHAEC genome-wide CRISPRi as a fourth held-out line (BULK folds) and an 811 consensus line

Written 2026-10-16 ~09:05Z by the 813 worker, for the judge (`experiments/800-checker/JUDGE.md`). The inventory behind it
is in `RESULTS.md`. No data file has been downloaded and no cloud has been spent ($0). Every number below that is not
cited is a prediction, written before the data exist.

## 0. The ask in four lines
- **Data:** IGVF's TeloHAEC genome-wide CRISPRi Perturb-seq (Engreitz lab, released 2026-10-12).
  - 19,276 genes and 114,885 guides.
  - The A/B/C panel is 300/300 in the library, and L1–L6 are 298–300/300.
  - Resting rep 1 alone has ≥ 100 cells for 295–299 targets per panel (`work/meta/telohaec_rep1_coverage_v2.json`).
- **Use (a):** a fourth held-out line on every BULK panel (L1–L6 → 24 W folds) and a fifth BULK-ABC fold. It is built with
  801's conventions and carries no refit, because the line is in no operator, key or composition fit.
- **Use (b):** an extra near-full line for 811's many-line consensus, delivered as a bank. 811 owns the try.
- **Cost:** stage 0 (the kill gate) ≤ $3; everything ≤ $15; cap requested **$20**. Zero configurations priced on BULK-W by
  this lane; the retrodiction rows are known answers, not tries.

## 1. Mechanism (falsifiable)
Every held-out read in the lab rests on three cancer lines: K562, HCT116 and HEK293T. The last two come from one lab
(Orion), and 811 measured that pair as the strongest mutual partners (+0.011 … +0.015, 811 RESULTS §3). So the W-folds
partly re-measure one lab's biology, and a change that helps "lines like HCT116/HEK293T" reads as transfer.

D/E/F are three unseen lines. A fold line from a different lineage (endothelium: not haematopoietic, not epithelial
cancer, not pluripotent) with near-full panel coverage tests transfer to an unseen line more strictly than the current
folds.
- **(a)** The checker reports TeloHAEC folds beside the cancer folds. A candidate whose gain also holds on TeloHAEC is less
  likely to be a cancer-panel or Orion-lab artefact.
- **(b)** In 811's consensus, a line covering ~all of any panel adds an independent vote on every target. The pilot that
  811 has (2,167 targets, 65/300 panel, median 67 cells) cannot do that.

**Falsifiable:** if the TeloHAEC-held folds give the wrong sign on ≥ 2 of the 3 priceable bulk known answers (§6 P2), the
line is not a board-like instrument. It then stays a diagnostic and never enters a ship decision.

## 2. Class
**Transfer** (AGENTS 11.1c.1): (a) is an instrument for pricing transfer; (b) is a public data source for a transfer
model. Nothing is fitted to A/B/C. The data are public, fixed, and available on the day.

## 3. The four problems the lab lead named, head-on

### 3.1 Multi-guide cells (mean 4.70 guides per cell; median 4; 10 % carry 0; rep 1)
A cell that carries target t also carries ~3.7 other random perturbations. A naive "carriers vs zero-guide cells" row
therefore estimates δ_t + (k̄−1)·μ̄ + technical(k), where μ̄ is the mean single-perturbation effect and technical(k) is
anything that tracks guide load. Guide load can track capture, cell size, depth or stress; GSE337988's paper reports that
guide overload induces stress and slows the cell cycle. Every target would then share that offset: it is a blob by
construction. Two estimators, the choice between them fixed by simulation only:

- **A (primary): stratum-matched background.**
  - Strata s = (n_guides_total, pool, rep).
  - For target t, carriers C_t are the cells whose assignment includes any of t's guides. The background is all non-carriers
    of t, re-weighted to the carriers' stratum mixture: P_bg(t) = Σ_s (|C_t∩s| / |s∖C_t|) · S(s∖C_t), with S the summed
    counts.
  - δ̂_t = log1p(TS·P_t/ΣP_t) − log1p(TS·P_bg/ΣP_bg), with TS 5e4: the donor-table space (`scripts/build_donor_table.py`)
    and the scorer's pds space.
  - Carriers and background then carry the same number of guides and the same technical stratum. The other guides' mean
    effect and guide-load effects cancel to first order. What is left is δ_t − μ̄: one random perturbation's mean effect,
    the same for every target. Most genes are null, so μ̄ should be small. The blob test (§3.3) decides whether it matters.
- **B (check): sparse deconvolution.** A linear model on per-cell normalised expression with per-stratum intercepts and one
  coefficient per target, using the cells × targets assignment (nnz 6.28 M in rep 1). It is solved as XᵀX β = Xᵀy by
  conjugate gradient over gene blocks. XᵀX is 19 k × 19 k, with ~3.5e7 nonzeros (Σ_c k_c²). B removes the
  co-occurrence term exactly rather than on average. Its coefficients are mapped back to the TS-5e4 log-ratio through
  each target's own carriers.
- **Known answer KA-SIM (decides A vs B, before any real row is read):**
  - Take rep 1's real assignment matrix and real library sizes.
  - Give each targeted gene a true row: K562-GWPS donor rows for the genes it carries (~9.8 k, realistic structure and
    magnitude), zero for the rest.
  - Base profile: TeloHAEC's background pseudobulk. Cell rates are L_c · softmax(log p + Σ_{t∈c} δ_t). Group sums are
    drawn Poisson from the summed rates on a 3,000-gene axis.
  - Simulate the same design twice: once with the real multi-guide sets, once with each cell keeping one random guide
    (same incidences per target, i.e. the single-guide twin).
  - **Pass:**
    1. A's median cosine(δ̂_t, δ_t) over the A/B/C panel targets ≥ 0.9 × the single-guide twin's;
    2. ‖mean_t(δ̂_t − δ_t)‖ ≤ 0.10 × median‖δ_t‖, i.e. no constructed blob;
    3. the code recovers δ exactly (cos > 0.9999) on a noise-free run.
  - B replaces A as primary only if it beats A's median cosine by ≥ +0.02 on the simulation. Nothing on the real folds
    chooses the estimator.
- **Known answer KA-SG (real data, single-guide subset):**
  - Rep 1's cells with exactly one assigned guide. Their count is unknown; the 25th percentile of n_guides is 2, so at most
    ~15 %.
  - The single-guide estimate (vs single-guide-stratum background) is noisy per target, so it is read only in aggregate:
    the regression slope of A's rows on the single-guide rows over the 100 stage-0 targets is in [0.8, 1.25], and the
    cross-target mean of (A − single) is within its bootstrap 95 % interval of 0.
  - A failure means the multi-guide rows are biased, not just noisy. STOP.

### 3.2 Shallow depth (median 5,452 UMI per cell, 2,590 genes; rep 1, all 62,757 features)
- **Per target, the pseudobulk depth is comparable to the folds in use, because each cell counts toward ~4.7 targets:**

  | source | cells or incidences per target | UMI per cell | UMIs per target |
  |---|---|---|---|
  | TeloHAEC resting, both reps (EST from file size) | ~550 incidences | 5.5 k | **~3.0 M** |
  | K562-GWPS | 175 | ~11 k | ~2.0 M |
  | Orion HCT116 | 154 | ~15 k | ~2.3 M |
  | Arc | 400 | ~20 k | ~8 M |

  The price: the other perturbations add variance (not bias, after 3.1) to each carrier pseudobulk.
- **Measured, not assumed:** stage 0 reports the rows' split-half residual cosine. The definition is 271's (two disjoint
  random carrier halves, background fixed, the cross-target mean removed), so it compares directly with the folds in use:
  271 measured +0.047 / +0.028 / +0.028 for K562 / HCT116 / HEK293T at 69–96 cells per half. It is reported by CPM bin of
  the gene, where shallow depth shows first.

### 3.3 615's blob test, before the line may serve as truth
- The 615 statistic, unchanged (`experiments/615-orion-hek-harness/PREREGISTRATION-blob.md`): split-half share
  n·⟨μ_A, μ_B⟩ / Σ_t⟨d_tA, d_tB⟩ in TS-5e4 log-pseudobulk space, on the gene axis minus panel genes.
- Data: 100 seeded A/B/C panel targets (seed 813, sorted list, ≥ 100 incidences), each capped at a seeded 250 carrier cells.
- **Decides:** estimator A's control (the stratum-matched background). Bar: **≤ 0.10**, as 615 required. For reference,
  615 measured raw HEK P′ 0.060, h1-test 0.049 and the thinned surfaces 0.82–0.88.
- **Reported, not deciding:** the naive zero-guide control, to show the mechanism if a blob exists (the 615 lesson: the
  blob came from control selection); the NTC-only control (64,464 NTC incidences in rep 1); and the top-PC share and its
  cosine with the long-gene axis, as 615 reported.

### 3.4 The 5-CPM gate on its axis
- **The scored axis** = official-axis genes whose mean per-cell CPM over TeloHAEC's background cells is > 5. This is the
  scorer's `filter_gene_min_cpm_cell` (5 CPM, control cells only; `docs/reference/vcc2026-metrics-BRIEF.md` l.83), applied
  to the line's own control-equivalent. Panel columns are out.
- One rule, fixed now. The 801 folds use the held bank's genes on the axis. At this depth a gene below 5 CPM gives
  ~9 counts per target pseudobulk at most (5 CPM × ~1.9 M UMIs in rep 1), i.e. a log-ratio sd of ≥ ~0.33, so it is
  counting noise.
- Stage 0 reports how many genes pass and the truth's split-half reliability per CPM bin ([5,20), [20,100), ≥ 100).
- `b_ctrl` for the fold is the background pseudobulk of all non-carrier cells, mixed to the panel's carrier strata.

## 4. Stages, gates and costs (all on EC2 through `labcloud.sh`, on-demand, Project=lab; nothing to the Mac but < 200 MB of tables)

| stage | input (fetched EC2-side from IGVF's presigned S3, us-west-2) | work | instance, cap | cost |
|---|---|---|---|---|
| **0: kill gate** | resting rep 1 only, IGVFFI2522NCFR (66.3 GB) | one streamed CSR pass: stratum sums plus carrier sums (halves) for the 100 stage-0 targets and all 300 A/B/C targets; KA-SIM (3,000 genes); KA-SG | r7i.2xlarge (64 GiB), 200 GB gp3, 120 min cap | **≤ $1.5** per launch; ≤ $3 with one retry |
| 1: the bank and the folds | resting rep 2 (44.5 GB) added | estimator A (or B by KA-SIM) for all 19,276 targets on both reps; bank in `build_donor_table` format (`delta.npy`, `b_ctrl.npy`, `perts.txt`, `genes.txt`, `meta.json`) to `s3://…/output/lab/813/`; 801 dumps for L1–L6/telohaec and ABC/telohaec with the production operator (`--fold-ops prod`; nothing to hold out); the engine's KAs; the retrodiction rows | same, ~3 h | ≤ $5 |
| 2 (optional, for 811 only) | oxLDL + IL-1β, 4 files, ~154 GB | the same bank per condition (a state axis for 811) | same | ≤ $5 |

### Stage-0 kill gate (all four must pass; otherwise STOP and record the negative with its mechanism)
1. **G0-KD knockdown:** on the 100 targets, the own-gene residual (carriers / matched background, the target's own CPM ≥ 10)
   has median ≤ 0.60 (615's Orion HEK: 0.48), with ≥ 60 % of targets < 0.70.
2. **G0-POWER:** median split-half residual cosine (271's definition) ≥ **+0.028**, the weakest fold line in use
   (HCT116/HEK293T), at the line's own halves. Share ≤ 0 must be ≤ 34 % (HEK293T's).
3. **G0-BLOB:** split-half share ≤ **0.10** (3.3).
4. **KA-SIM and KA-SG pass** (3.1).

The stage-0 targets and seed are fixed above. Stage 0 reads no fold, no arm and no board number.

## 5. Day-of path (Oct 22 → Nov 5)
- The bank is public and fixed. For D/E/F's own panel (released Oct 22), TeloHAEC rows exist for ~all 300 targets (the
  library is genome-wide), so a **TeloHAEC-held fold on D/E/F's real panel** can be built the same day. That gives 708's
  day-of re-pricing a fourth fold line next to the three cancer lines, at full coverage.
- Nothing in it is fitted to A/B/C. The estimator choice comes from simulation, the gate is the scorer's rule, and the
  blob test is preregistered.
- For (b): 811 reads the D/E/F panel targets' rows from the bank like any other line, centred per line by its own
  factor C.

## 6. Predictions, bins and falsifiers (written before any row exists)
- **P1, stage 0** (the bins are the gate's): G0-KD median residual 0.35–0.60; G0-POWER median split-half residual cosine
  +0.03 … +0.08; G0-BLOB 0.03–0.10; KA-SIM A ≥ 0.9× the single-guide twin. My probability that all four pass: ~0.6. The
  likeliest failure is G0-BLOB through a guide-load technical effect that stratum matching does not fully absorb.
  **Falsifier: any gate fails → the line is not truth.**
- **P2, the retrodiction rows** (stage 1, the engine's own named arms; the signs come from 801's known-answer list).
  ABC/telohaec reads:
  - **495 cap 6** (vs cap 3) **positive**;
  - **498p2g** (vs cap 6) **positive**;
  - **702** (PIE w 0.25 vs ALLSW) **≤ 0**. PIE never saw TeloHAEC, so this read is fair for it (702 rule 2).

  The same three signs on ≥ 4/6 of the L-panel TeloHAEC folds (for 702: where its rows exist).
  - Bins: **TRUSTED** = 3/3 on ABC/telohaec and ≥ 4/6 on L; **DIAGNOSTIC** = 2/3; **REJECTED** = ≤ 1/3.
  - 479sp2 is not priceable in bulk (PHASEB), so it is not counted.
  - The over-read ratio ρ for the TeloHAEC folds is reported per class next to the cancer folds'.
- **P3, levels:** the incumbent's aggregate cosine on the TeloHAEC folds is below every cancer fold's on the same panel
  (cross-lineage transfer is harder), at 0.4–0.9× of their mean. If it reads above every cancer fold, suspect a shared
  component (a residual blob that passed 0.10), and re-run the blob test on the fold's scored axis.
- **P4, for 811 (b): one configuration** (M811's chosen model + TeloHAEC rows), counted in 811's family, not here.
  - Prediction: Δoverall on W-conf in [−0.0005, +0.0015] vs the same model without it.
  - Its per-line shared share is below 0.10 after factor C. KOLF2's 0.245 is the cautionary case.
  - Falsifier: a dropping-it ablation that gains on ≥ 2/3 A/B/C cancer folds → excluded from 811's list. That is how KOLF2
    z-scores were caught.

## 7. Novelty against the record (what this is NOT)
- **Not the TeloHAEC pilot as a stack donor.** 223/225 priced `telohaec@0.28` at **−0.00120** (225 RESULTS l.22, l.159):
  65 panel targets, median 67 cells, a 5th donor in the whitened stack. That falls under "sources saturate"
  (218/225/229/248), and this proposal adds no donor to the production stack. Use (a) is a held-out truth line. Use (b) is
  a vote in 811's raw consensus, a construction 811's PROPOSAL §5 already separates from the stack closures.
- **Not 122/135's TeloHAEC panel atlas** (135: terminal NO-GO on 65 targets × 67 cells, an a15 integer-transport
  endpoint). That was a prediction source scored against TeloHAEC truth, with the noisy 2kG pilot as both.
- **Not new truth for DE members.** These are bulk folds only. DE members stay with DE-H1 / DE-CD4 (DESIGN §4). 3′ chemistry
  at 5.5 k UMI with multi-guide cells is not an emitter surface, and I do not propose one.
- **New:** the release is four days old. It is genome-wide (19,276 genes, ~300/300 on every panel, ≥ 183 incidences at
  the 10th percentile in one replicate), not the 2,167-gene pilot. No experiment on disk has read it.

## 8. Cost and tries (best-of-N accounting)
- **Cloud:** stage 0 ≤ $3 (the kill gate); stage 1 ≤ $5; stage 2 ≤ $5. **Cap $20** (lane budget $50; spent so far $0).
  Every launch ≤ $50 and is logged in `docs/LAB-LOG.md` by the launcher.
- **Configurations:**
  - Estimator: 2 (A, B), chosen on simulation only (KA-SIM), never on a fold.
  - Control, gate and stage-0 targets: 1 each, fixed above.
  - Priced arms on BULK-W from this lane: **0**. The retrodiction reads are the engine's named known-answer arms, and they
    choose nothing. If a gate or a KA fails, nothing is re-tuned; the negative is recorded.
- **Wear:** the TeloHAEC folds start unworn. Their first reads are the known answers. Their wear entries go to 800's ledger
  as `surface=telohaec`, so later candidates' best-of-N counts include them.

## 9. DLD-1 (GSE337988): not proposed now
TeloHAEC wins on value per dollar.
- **Coverage:** DLD-1 overlaps each L panel at only 65–82 targets, against ~300; it would need two new DLD-1-specific
  panels (pool 2,421) to be a full fold.
- **Independence:** DLD-1 is MSI colorectal, the same lineage and MSI status as HCT116, so it adds less independence than
  endothelium.
- **Power:** 116 measured a median split-half cosine of 0.070 on the author LFCs, and its per-promoter cells are an
  estimate (~90–150).
- Its one edge is a cleaner design (MOI 0.3, single guide) in a cancer line. **It becomes the fallback** if TeloHAEC fails
  stage 0 on G0-BLOB or KA-SG, the multi-guide failure modes: ~5 GB of processed objects, ≤ $5, the same gates.

## 10. Not done, and not claimed
- No IGVF or GEO data file has been downloaded. The facts above come from catalogue records and byte-capped metadata range
  reads (`RESULTS.md` §0).
- The IGVF data-use summary reads "no certificate". The consortium's external-use terms are re-read at stage 0, before
  any fetch. Pisces is not part of this proposal (still empty at 06:19Z).

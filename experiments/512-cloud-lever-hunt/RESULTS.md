# 512 — fresh-model lever hunt

**Decision.** Against the best entry (446bj2, 0.2723), Illumina is **+0.0850 overall, of which +0.0676 (79.5%) is the four DE members**. No unused property of the scorer, and no public corpus the 09-26 hunt missed, is large enough to close that. The one unread public object that already posts leader-like direction on H1’s own line is the 09-30 slot (`experiments/458-arc-lab-programme`). Two cheap falsifications are worth pricing. Both have priors of a few thousandths of overall, and a large chance of a null. A renamed closed item is not among them.

Tags: **VERIFIED** = computed here from `data/leaderboard/*.json` or read from a cited RESULTS/BRIEF. **REPORTED** = a secondary page or an earlier experiment’s number, not recomputed on this VM.

Script: `board512.py` (stdlib). Output: `work/forensics512.json`. `gaps_vs_our_best` in that file is against the **displayed** row (0.2704006), not 446bj2. The table below recomputes against 446bj2.

---

## 1. Board forensics

**Census (VERIFIED).** 180 snapshots, 2026-08-21T19:06Z → 2026-09-27T07:44Z. 160,841 rows, **11,966 unique entries** (deduped by `id`), **1,255 teams**, 10,711 consecutive steps. `max |score_avg − mean of six scaled members| = 3.33e-16`. Bit-identical = `|Δ| ≤ 1e-9`. A member “moves” above pds 5e-4, mse 1e-3, nmae 1e-3, fid 5e-4, reach 1e-3, jac 3e-4. **298 mixed steps** (some raw members bit-identical and some moved). The API stores only each team’s current entry, so a submission that was replaced before the next snapshot is absent: captured / `submission_count` is Illumina 38/50, cqawesome 50/75, Testh 29/44, us 29/56.

**446bj2 is not in the snapshots (VERIFIED).** `data/board-entries/YvVru3DAF89Jiyr34olR.json`, captured 2026-09-27T00:27Z, submitted 00:01:30Z, `score_avg` 0.2722766434916974, rank 6 at capture. 446br2 replaced it at 00:28:08Z (`the pipette rack learned origami`, 0.270400605862736) before the 01:44Z snapshot. Displayed rank at 07:44Z is **7th**. People-key dedupe does not merge Illumina with anyone else in the top 12. Ranks 11–12 are Ava Brown 0.2535 and HY H 0.2533 (the pack, outside this top 10).

**Gap vs 446bj2 (VERIFIED).** Overall points = Δscaled / 6, using each team’s published scaled members. Positive = they are ahead. Raw Δ is theirs − ours. nmae raw down is a gain.

| team | overall | Δ overall | bulk pts | DE pts | DE share | raw nmae | raw fid | raw reach | raw pds | raw mse | raw jac |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina AI | 0.3573 | +0.0850 | +0.0174 | **+0.0676** | **79.5%** | 0.8138 | 0.5466 | 0.4277 | 0.8877 | 0.5711 | 0.0383 |
| cqawesome | 0.3212 | +0.0489 | +0.0084 | +0.0405 | 82.7% | 0.8273 | 0.5380 | 0.3366 | 0.8731 | 0.5900 | 0.0375 |
| Testh | 0.3061 | +0.0339 | +0.0071 | +0.0268 | 79.1% | 0.8512 | 0.5280 | 0.3302 | 0.8615 | 0.5742 | 0.0364 |
| CellWM | 0.2794 | +0.0071 | −0.0138 | +0.0209 | (bulk is ours) | 0.8676 | 0.5255 | 0.3088 | 0.8694 | 0.7093 | **0.0451** |
| nailong king | 0.2747 | +0.0024 | −0.0138 | +0.0162 | (bulk is ours) | 0.8513 | 0.5216 | 0.2884 | 0.8641 | 0.6976 | 0.0389 |
| chill | 0.2708 | **−0.0014** | −0.0302 | +0.0288 | (bulk is ours) | **0.8351** | 0.5303 | 0.3076 | 0.8436 | 0.7488 | 0.0368 |
| **446bj2** | **0.2723** | 0 | 0 | 0 | | **0.8971** | **0.5144** | **0.2852** | **0.8674** | **0.6264** | **0.0408** |
| displayed (446br2) | 0.2704 | −0.0019 | | | | 0.8956 | 0.5144 | 0.2872 | 0.8672 | 0.6410 | 0.0408 |

Illumina’s DE points: nmae +0.0235, fid +0.0184, reach +0.0270, jac −0.0012. Raw Δ: nmae −0.0833, fid +0.0322, reach +0.1424. cqawesome: nmae +0.0192, fid +0.0133, reach +0.0096 points. Testh’s pds is worse than ours (−0.0021 points); the gap is mse + nmae + fid + reach. CellWM is the only top-10 team ahead of us on jac (+0.0020 points, raw 0.0451 vs 0.0408).

**Dominant partial move in the field (VERIFIED).** Of 298 mixed steps, the most common moved-sets with pds held are `nmae+fid+reach+jac` 62, `fid+reach+jac` 32, `reach` only 31, `nmae+reach` 20. That is the unwelded construction we already ship (389/392): log-pseudobulk sums pinned, DE channel free. `pds+mse` alone is 3 steps.

### 1.1 Steps where some members stayed bit-identical

Bounds use `docs/reference/vcc2026-metrics-BRIEF.md` (read this session).

**cqawesome, three nmae-primary steps, pds, mse, fid, jac bit-identical (VERIFIED).** Whole board, only five steps have pds+mse+fid+jac bit and `|Δnmae| ≥ 0.001`. Three are cqawesome:

| date | models | Δnmae raw | Δreach raw | Δ overall |
|---|---|---:|---:|---:|
| 2026-09-18T06:27 | m0920d → m0921a | **−0.02064** | +0.00041 (under threshold) | +0.00588 |
| 2026-09-22T04:24 | m0923r → m0923s | −0.00158 | +0.00357 | +0.00104 |
| 2026-09-22T17:05 | m0923s → m0924kc | −0.00687 | +0.00207 | +0.00238 |

The other two are PumpkinChicken (v44→v45, nmae −0.0021, reach +0.008; v50→v51, nmae +0.0023, reach bit). Bound: pds and mse bit ⇒ the log-pseudobulk sums and the jackknife credit did not move. fid and jac bit ⇒ the call set and the count of correct call signs did not move (a mean of 300 Jaccards matching to 1e-9). nmae moved ⇒ arithmetic means of reference-DE genes changed, including genes outside the call set. Reach staying near-flat means the 90% purity depth did not move. This is a pure magnitude edit on the gate with call signs held. **Tried:** coarse MDE (433a, nmae −0.0168 raw, pds/mse bit; FINDINGS §0) and zero-gene repair (433h, nmae −0.0247). 480: past the MDE, at our signs, oracle magnitudes recover only −0.0985 on H1 and no class rule beats the shipped MDE; the −0.254 pot is signs. cqawesome’s −0.0206 matches a first MDE-like fix, not a construction we lack. 480 already recorded the −0.0206 and −0.0069 siblings. cqawesome publishes no descriptions (0 of 50 captured entries).

cqawesome also has 12 mixed steps with pds and mse both bit, from 09-07 on (fid/jac only, or nmae+fid+reach+jac). Largest: 09-07 fid +0.0120, jac +0.0010, overall +0.0072; 09-19 fid −0.0212, jac +0.0033, overall −0.0106 (a losing call-set edit); 09-24 nmae −0.0016, fid +0.0074, reach +0.0022, jac −0.0007, overall +0.0048. 09-13 mse −0.0535 with pds flat (not bit) is the bulk-magnitude step FINDINGS already notes. **Tried** as the unwelded DE channel plus call-set edits (192, 385). 385: shrinking K loses −0.033 overall; precision stays 0.507–0.517 at every key and every size (192).

**Testh 09-23 is not a bit-pure DE step (VERIFIED).** 2026-09-23T01:18, code `pmnfrj`, overall +0.0262. Raw: pds **+0.00201**, mse **−0.00199**, nmae **−0.04807**, fid +0.00538, reach **+0.04364**, jac +0.00093. BOARD-INTEL’s “pds/mse flat” is the right scale and the wrong identity test. Bound: jac +0.0009 ⇒ the call set barely moved. fid +0.005 is a handful of extra correct calls if n_pred is ~2,600. nmae −0.048 and reach +0.044 are far larger than a call-sign tweak, so arithmetic means and signs of gate genes moved, with a small pseudobulk leak. That is the 496 channel at ~30× the hook. **Tried as delivery, not as a source:** 496 with CLEAN4 is nmae −0.00146, overall-equiv +0.00018, P(≤0) 0.41 (WEAK). The source that matches leader fid (~0.55) is Arc’s programme (462: fid-like 0.5545 vs CLEAN4 0.5244 vs oracle 0.5544 on H1’s own line), already the 09-30 read.

Testh steps that are bit-pure on pds and mse:

| date | Δnmae | Δfid | Δreach | Δjac | Δ overall | bound |
|---|---:|---:|---:|---:|---:|---|
| 09-20T09:29 | −0.01646 | +0.00464 | −0.00539 | +0.00141 | +0.0069 | MDE / zero-repair shape. 480 recorded −0.0165. **Tried** (433a, 433h). |
| 09-25T06:34 | +0.00001 | +0.00522 | **+0.02792** | +0.00183 | +0.0090 | Means pinned. Ranks and the call set moved. Descriptions already silent (last text 09-15). |
| 09-25T10:32 | −0.00866 | −0.00026 | −0.00545 | ~0 | +0.0012 | Means of the gate moved; calls almost fixed; reach got worse. Magnitude edit. **Tried.** |

Last Testh description (09-15, REPORTED from the stored API text): library-size tilt removed, and a projection on the top 14 directions. The tilt claim is the unwelded channel. The projection is 483, closed at his disclosed k=14, β=0.5 (his own LOCO basis untested; it moves pds, which is not the DE gap). Earlier descriptions (09-12, 09-14) name the reliability weight `0.15+0.85·z²/(z²+8)·n/(n+4)` and an 8-feature learned gate. 473 re-priced that weight on the sharpened bulk and closed it while q ≥ 1.5.

**chill’s 0.835 is the 09-11 entry, and the informative partial step is the next day (VERIFIED).** Literal first submission 2026-09-05 is overall 0.0844, nmae 0.9407, pds 0.6417, mse 5.80, not 0.835. First time nmae hits 0.835: **2026-09-11T18:32**, BioDose-JointDecision, overall 0.2291, nmae **0.8348**, pds **0.8239**, mse **0.8451**. The step into it moves everything (pds +0.0644, nmae −0.0282, fid +0.0104, reach +0.0507). Weak bulk is real: pds 0.824 and mse 0.845 against our then-better bulk. Stored description: family-balanced reference aggregation, signed conditional responses, separate DE membership and direction ranking, jointly constrained counts. No query perturbation outcomes.

The partial step is **2026-09-12T19:19**: pds and mse **bit-identical**, nmae −0.00016 (flat), fid +0.00386, reach **+0.03381**, jac +0.00347, overall **+0.01014**. Description (truncated in the snapshot, the constraint sentence is complete): DE membership trained with **equal total weight per source study family**, separate directional ordering, bounded joint integer refinement, membership rank weight 0.5, **each cell and gene total preserved**, cell-mean LFC drift bounded by 0.01. Bound: fixed sums (pds, mse) and fixed means (nmae) with ranks, the call set, and order free. That is a fixed-mean free-rank edit. Later chill steps improve bulk (09-18 mse −0.042; 09-19 pds +0.009, mse −0.036; 09-21 pds +0.006, mse −0.035) while nmae stays ~0.83–0.835. 09-24 again has pds and mse bit, reach +0.0143, fid +0.0050, jac −0.0043. 434f/480: nmae 0.835 maps to gate sign accuracy ~0.65–0.66 if magnitudes are L1-optimal, and that accuracy is an **upper bound** if magnitudes follow confidence (REPORTED from 480).

**Illumina has no four-member-bit step (VERIFIED).** Zero steps with pds, mse, fid, and jac all bit-identical. Trajectory is full-profile: broken test1 on 09-02 (overall −0.966, mse 425, nmae 6.88) → test2 0.138; PromoterAI entries swinging fid by −0.22 then +0.25; 09-11 mostly a bulk step (mse −0.218, nmae +0.001); 09-18 a −0.043 then a +0.065 recovery (mse −0.191, nmae −0.036, fid +0.022, reach +0.029) with pds not in the moved set; 09-21 reaches 0.3598; later steps trade mse and nmae together. No partial-move construction to copy. The proprietary-atlas reading is REPORTED (FINDINGS §1), not re-verified here. No descriptions.

**Our captured steps match the ledger at the big deltas (step numbers VERIFIED; experiment ids REPORTED from FINDINGS).** 09-13 stapler: pds and mse bit, nmae +0.0041, fid +0.0082, reach −0.0039 (392’s CLEAN4 prior; H1-BOARD-CALIBRATION quotes the same three deltas). 09-22 spreadsheet: overall +0.0142, pds −0.0036, mse −0.079, reach +0.0088 (410b ballast). 09-24 centrifuge: nmae −0.0258, mse −0.026 (433h repair + relevel). 09-25 incubator: pds +0.0051, mse −0.031, nmae −0.024 (q plus the repair, if the previous snapshot was pre-repair).

---

## 2. Structural levers

From the BRIEF. A construction is legitimate when it is a fitted prediction from public data and the recipient’s controls, not a lookup of the recipient’s real effects.

**Already used.**

| property | member | where |
|---|---|---|
| Ballast cells: a few deep cells set the pseudobulk and carry 1/400 of the weight in means and ranks | mse without the nmae cost | 410b |
| Unwelded per-cell libraries | pds/mse bit, arithmetic means free | 389/392 |
| Exact-U mass moves inside control intervals | nmae −0.0018 to −0.0026 on H1; U, ties, calls, signs held | 423/429 |
| Coarse MDE add | nmae −0.0168 raw on the board; pds/mse bit | 426h, 433a |
| Zero-gene repair, floor −5 | nmae −0.0247 raw | 428, 433h |
| Call set = top 400 by \|H\| ∪ top 2,200 by control CPM | fid/jac/reach | 192 |
| q80 protection of the strongest calls | reach | 395 |

**Looked at, not a new lever.**

- **Fixed-mean, free-rank room.** chill’s 09-12 step is the existence proof on the board (+0.010 overall, means and sums pinned). 338 measured the other direction (fixed rank, free mean: median 0.27% symmetric CPM; 20% without the rank constraint) and did not build this dual. 385 already emits an alternative call set with pds, mse, and nmae pinned, so the missing object is the score, not a solver. That score is proposal 1. A new LP under a new name is 338/the rank-moment LP, which is on the closed list (325–342).
- **Reach’s deepest prefix can recover after early errors** (purity ≥ 0.9, not the first break). One wrong call in the first 9 can still recover at k=10. 406: the walk breaks inside the head (mean k* ≈ 13–14), and no order key extends it (neighbourhood keys ≤ +0.0014, CIs cross 0). 348’s best non-oracle reorder is +0.0264 raw reach, concentrated in pools of ≤20 genes, and missed a +0.03 bar. Uncalled-gene signs cannot extend k* until the head stays ≥90% pure.
- **fid’s precision asymptote is not a knob.** At n_pred ≈ 2,600 and median n_real 31–45, most calls sit on noise signs, so precision sticks at 0.507–0.517 (192, 385). Under-calling is free on fid until n_pred hits n_real and fatal on the heavy tail.
- **nmae weights perturbations equally**, so a 10-gene gate counts as a 2,000-gene gate. 467/470/480 already use that frame. Predicting 0 everywhere scores 1. L1-optimal magnitude given sign accuracy a is positive only for a > 0.5.
- **Genes at ≤5 CPM in control are invisible to all four DE members** (same gene set for every submission). Counts spent only there move pds/mse. Not audited in the emitter on this VM (no control matrices).
- **Wilcoxon tie term.** The FULL doc says “normal approximation to the rank sum” and does not state a tie correction. A code search of ArcInstitute/cell-eval2 returned HTTP 429. **UNVERIFIED.** 423 holds ties fixed. Do not manufacture ties without reading the scorer.
- **Jackknife rho.** Real signal has rho = 1, so noise cannot farm mse (FINDINGS). Closed.
- **“Number of predicted cells unconstrained”** in the metric text. The CLI contract is exactly 400. Not available.
- **Panel genes are dropped from pds only.** They still enter DE when they pass the CPM filter and are not the perturbed gene. Whether our emitter zeros them was not checked (no build on this VM).

---

## 3. Literature, 2025–2026 — what is new

Nothing in this pass is a new information source. The new object is a **training target** chill already described in public (proposal 1), which is a label, not a corpus.

| item | what it adds | status |
|---|---|---|
| State (Cell 2026, [S0092-8674(26)00921-9](https://www.cell.com/cell/fulltext/S0092-8674(26)00921-9); bioRxiv [2025.06.26.661135](https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2)) | Context transfer from a large observational embedding | **Closed here** (FINDINGS: pds ≈ 0.50). REPORTED. |
| Ahlmann-Eltze 2025, Nature Methods | Deep models do not beat linear/additive baselines on unseen single genes | Already the reason the stack is a fitted linear transfer. REPORTED. |
| ResearchSquare rs-10434123 (page text read this session) | Unseen single-gene: State +26% vs linear +19% vs GEARS +22% over a cell-mean baseline; on Replogle the three are within noise (0.36 / 0.34 / 0.33 pΔ). Cross-context, every method loses 0.17–0.25 pΔ; State keeps more of its within-context number (0.55→0.38) than linear (0.49→0.24) | **REPORTED** from that PDF, not recomputed. The metric is top-20 DE Pearson, not our nmae/reach. The residual over linear is a few points and does not identify weak-gene signs. Not a lever. |
| TxPert (Nat Biotech 2026, [s41587-026-03113-4](https://www.nature.com/articles/s41587-026-03113-4); arXiv 2505.14919) | STRING+GO graphs; GEARS falls below the non-learned baseline on K562/RPE1/HepG2/Jurkat | Graph signs already lost to \|H\| (190). “No basal encoder at ≤4 lines” is 332/353. REPORTED. |
| VCC 2025 winners ([Arc wrap-up](https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up)) | BM_xTVC: DEG frequency and mean expression as features, because pure deep models lost to statistical baselines on DES and MAE. XLearning: pseudobulk residual FCN. Outlier/TransPert ([arXiv 2511.14954](https://arxiv.org/abs/2511.14954) class): summary-level Wilcoxon plus a global scale. Altos PRiMeFlow ([arXiv 2604.13986](https://arxiv.org/html/2604.13986v1)): flow matching; the winning move was locating held-out perturbations in other public lines | DEG frequency is 192’s control-CPM detectability prior, not a new source. The 2025 scale trick is dead: 2026 pds is scale-invariant. Same-line location is 451, already the October mode. REPORTED (482c plus this session’s wrap-up snippet). |
| AROMA (arXiv 2604.20263) | Three-way up / down / non-DE classification | The label class of proposal 1, not a new corpus. REPORTED from 353. |
| 482c (2026-09-26), re-read this session | No new public genome-scale CRISPRi transcriptome for a cancer line outside banks already held. X-Atlas/Pisces Jurkat+HepG2 is an empty Hugging Face repo (`usedStorage` 0 since 2026-03-25). A549 Cas9 KO CROP-seq GSE345058: median on-target residual 0.58 of NTC CPM. HeLa and CAL-33 have no public knockdown atlas | **REPORTED** from 482c’s RESULTS (that session verified the files). Not re-downloaded here. |
| CD4 Cell 2026 atlas; KOLF2 atlas | Already the production CD4 bank. KOLF2 is NO-GO as a donor (−0.00174) and covers 181–251/300 of synthetic October panels | REPORTED, FINDINGS §0 and 381. |

---

## 4. Proposals

Ranked by expected transfer to a new panel and new lines. Neither is a path to 0.357.

### 1. Family-balanced source-study DE membership as the call key

**Mechanism.** Reach orders significant calls first, then by p, then by |lfc|. fid = precision × min(1, n_pred / n_real). jac is the set overlap. If cell and gene totals stay fixed, pds and mse stay bit-identical. If arithmetic means stay fixed, nmae stays flat. chill’s 09-12 step is that signature on the board: reach +0.0338, fid +0.0039, jac +0.0035, nmae −0.00016, pds and mse bit, overall +0.0101.

**Information.** Public source screens. The label is “this gene was significant in the source study,” with equal total weight per study family, not |H| and not CLEAN4’s shrunk P_down. Recipient controls enter only as the CPM gate the scorer already applies. No lookup of A/B/C effects.

**Members.** Reach up if the new calls are sign-purer in the head. fid up only if precision at matched size rises (192 says it does not, for every key tried). jac either way. nmae, pds, mse pinned if the emitter holds means and sums, which 385 already does for an alternate key.

**Closest closed item, and the difference.** 192’s call set is top 400 by |H| ∪ top 2,200 by control CPM, and predicted |lfc| added nothing to membership over control CPM. This is a different label (source significance, family-balanced). It will still correlate with control CPM, which is why the prior is small. 385 closed the **size** K, not which genes. 406 and 348 reorder a fixed set; this changes which genes are called. 09-13’s “prior signs on non-called genes” and 496 move means; this holds means.

**How to price.** First `experiments/192-predicted-callset/enrich192.py`: overlap of the source-DE key with the reference call set at matched K, against control CPM. **Kill if it does not beat control CPM.** Only then change the key inside `experiments/385-per-target-callset-size/run385.py` (not K) and walk reach with `experiments/348-reach-ordering/reach348.py`.

**Prior.** +0.001 to +0.004 overall if the key beats CPM, and ~0 if it does not. chill’s +0.010 is on their base; the reach piece is +0.0063 points (scaled reach +0.038 / 6). Ordering transfers about 0.4× (H1-BOARD-CALIBRATION), which is ~+0.0025 from reach plus a small fid term. 348’s best reorder on our own set was +0.026 raw reach on H1 and still missed the bar. Do not spend a slot on it before the overlap test.

### 2. Genotype-stratified donor signs

**Mechanism.** nmae is L1 on the reference significant set. 480: each +0.01 of gate-wide sign accuracy is −0.006 to −0.008 raw nmae. A discrete match on TP53 / KRAS / BRAF / CDKN2A chooses which donor’s sign to trust. Parameter-free. The recipient’s mutation status is basal identity, not a perturbed profile.

**Information.** CCLE/DepMap mutation calls for the donor lines and, once a line is identified, for D/E/F. For A/B/C the genotypes are already known.

**Members.** nmae down only if gate sign accuracy rises. fid is nearly blind at our call-set size (+0.08 gate signs → +0.007 fid, 480). Reach moves only if the new signs land in the head.

**Closest closed item, and the difference.** The August closed list includes basal-similarity donor weighting, and 379 killed continuous basal similarity (best +0.004 sign accuracy, below +0.005). 289 used p53 **expression** programmes. Neither is a discrete mutation stratum. 353 named this construction and did not run it; 379 reconfirmed it was still unexecuted.

**How to price.** The six-line sign frame in `experiments/379-similarity-weighted-donors/run379.py` (`experiments/334-text-conditioned-transfer/selection.json`). Replace the similarity weight with a mutation-matched indicator. Falsifier, written in 353 before any fit: permute the mutation labels. **Kill if matched ≤ unstratified consensus, or if the permutation matches.**

**Prior.** 0 to +0.002 overall. 353’s own reason: at six lines the strata are one or two lines, so the falsifier is the expected outcome. Nadig’s essential-gene line correlations do not apply to this panel (FINDINGS: 0/300 essential). Price it only because it is cheap and unnamed in §2. Do not build past the permutation.

### Not proposed

- Another MDE, exact-U pass, or zero repair. cqawesome’s nmae-only steps are that signature, and 433a / 433h already shipped it. 480: magnitude past the MDE at our signs is the smaller pot.
- A second copy of Testh’s 09-23 shape (nmae −0.048, reach +0.044, jac +0.0009, bulk almost flat) with CLEAN4 as the source. That is 496, priced at +0.00018 overall-equiv.
- Arc’s H1 programme as a new prior. It is 458d, the 09-30 slot. 462 already measured the fid-like row at 0.5545 on H1’s own line.
- Testh’s 8-feature reliability score as a DE rank key. Same functional as the order keys 406 closed, and the same features 473 closed on the bulk while q ≥ 1.5. 172’s learned sign-confidence key was inside the noise of |H|.
- Testh’s top-14 projection. 483 closed it at the disclosed k and β. His own basis is untested and would move pds, which is not the gap.
- Expression-gated STRING / CollecTRI propagation (353’s first item). 190: where the network speaks, its sign is less accurate than |H|; the direct-regulon oracle ceiling is reach +0.00078, +0.00015 overall.
- A learned sign or magnitude model. 434: the sign-confidence GBM (AUC 0.589) plus a conditional-median magnitude lost to the MDE on panel-like targets (−0.0066 vs −0.0111). 480 closed the 24-parameter class rules.
- STATE, scGPT, GEARS, flow matching, GenePT. Closed (383; sprint 325–342; FINDINGS pds ≈ 0.50). The ResearchSquare residual over linear is not our metric.
- A new public CRISPRi line. 482c: none. A549 knockout residual 0.58 is the wrong modality. Pisces is empty. CD4 2026 and KOLF2 are already in hand.
- Per-target K or a smaller call set. 385: shrinking K loses −0.033 overall.

---

## Coordinator's disposition (2026-09-27 12:05Z, on master)

- **Reproduced locally:** `board512.py` on 181 snapshots (one more than the agent's) gives 12,020 unique entries and 1,258
  teams, consistent with 11,966 / 1,255 on 180. chill's 2026-09-12T19:19 step reproduces exactly: pds and mse
  bit-identical, reach +0.03381, fid +0.00386, jac +0.00347, overall +0.01014. The committed `work/forensics512.json` is
  the agent's 180-snapshot output.
- **Headline accepted:** against 446bj2, 79.5 % of Illumina's +0.0850 is the four DE members; nothing structural or
  public closes it. It agrees with 480 and 482c.
- **Proposal 1 (family-balanced source-significance call key): not queued.** Its reach mechanism is contradicted on H1 by
  427f: ordering our own call set by donor agreement read 0.2547 against our p-order's 0.2747, and 406/348 closed every order
  key. chill's +0.034 was measured from its own earlier base, and chill's reach today (0.308) is only +0.023 above ours. The
  kill test (source-significance enrichment vs control CPM at matched K) needs per-(target, gene) donor Wilcoxon tables
  that are not on disk. Revisit only if a later read shows our head's sign purity, not membership, is the reach limit.
- **Proposal 2 (genotype-stratified donor signs): not queued.** Its own falsifier (permuted mutation labels) is the
  expected outcome with one or two lines per stratum (353's reason), and 379 closed continuous basal similarity.
- Useful corrections to the record: Testh's 09-23 step was not bit-pure on pds/mse (+0.00201 / −0.00199 raw), and chill's
  literal first entry (09-05) was 0.0844; its nmae first reached 0.835 on 09-11.
- **§2's unaudited item ("whether our emitter zeros the panel genes") is closed: it does not.** In 446bj2's stack dump
  (`417/work/dump_h1x2_cd4x15_lodo05.npz`) the off-diagonal panel coordinates (each row's other 299 panel genes) are
  populated: 89,570 of 89,700 nonzero in both `bulk` and `carriers`, 1.8 % / 2.0 % of the row energy, median |value|
  1.24× the non-panel median (panel genes are highly expressed by selection). The builder's "zero panel support" applies to
  the whitened residual space only; the rows that reach the emitter carry panel predictions. 127's ceiling (2.3 % of DE
  mass) therefore bounds a correction of those predictions, not a missing block.

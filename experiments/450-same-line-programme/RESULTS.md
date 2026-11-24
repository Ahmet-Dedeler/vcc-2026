# 450 — A same-line bank from a DIFFERENT study carries the recipient's response: programme signs, per-target signs, and per-target direction (H1, analytic)

**Date:** 2026-09-25 02:35–02:50Z (coordinator). Surface: h1-test (the only board-validated local surface), truth = the
scorer's own Wilcoxon table (`323-current-family-validation/work/reference.npz`). Same-line bank: Pan 2025's hESC CRISPRi
screen (GSE295214 pseudobulk database, `data/external/gse295214/pseudo_bulk_perturbation_database.h5ad`, dataset
`Pan2025_hESC`, 2,354 perturbations). It is a different lab, library and chemistry from the 2025 VCC H1 Flex data that
h1-test comes from. Scripts: `sameline450.py`, `common450.py`, `signs450.py`, `bulk450.py` (each runs in < 1 min, < 3 GB).

## Why
441 put the DE gap to the leaders in the recipient's own shared response programme (the leave-one-target-out mean
response reaches 0.645 sign accuracy on gate-block genes, CLEAN4 0.580) and found no cross-line route to it. It never
tried a same-line bank from another study. DEF-RUNBOOK §3 row 220 says a D/E/F line with its own public Perturb-seq would
be in-domain training data, needing "its own A/B/C-style validation — which must be finished before Oct 22 or not
attempted". Nobody had done it. H1 is the one surface where we hold such a bank and the truth.

## 1. The programme: one sign per gene, the same for every target (Pan rows for h1-test's 85 targets excluded)
Known answers: CLEAN4 0.580 on gate-block and the oracle 0.645 both reproduce 441 §2/§3.

| predictor | gate-block | gate-nonblock | all gate | block, all adjudicable genes (fid-like) |
|---|---:|---:|---:|---:|
| CLEAN4 (shipped; 4 cancer banks) | 0.580 (cov 0.75) | 0.567 | 0.573 | 0.524 |
| **Pan hESC P_down (same line, other study)** | **0.619** (cov 0.96) | 0.570 | 0.587 | 0.539 |
| Pan mean delta sign | 0.614 | 0.568 | 0.584 | 0.536 |
| z(CLEAN4) + z(Pan) | 0.604 | 0.576 | 0.586 | 0.534 |
| H1 leave-one-target-out mean (oracle) | 0.645 | 0.614 | 0.625 | 0.554 |

On the genes both priors cover (`common450.py`):
- gate-block: CLEAN4 0.5807, Pan 0.6253, oracle 0.6529. Pan closes **62 %** of the gap.
- all gate: 0.5736 / 0.5944 / 0.6378, 32 % of the gap.
- fid-like block: 0.5244 / 0.5390 / 0.5554, 47 % of the gap.
- Per target (49 targets with ≥ 10 common gate-block genes): **Pan − CLEAN4 = +0.057 [95 % bootstrap +0.024, +0.091];
  Pan better on 37 / 49.**
- Mixing CLEAN4 in makes Pan worse on the block. The same-line bank should replace the cross-line prior, not join it.

## 2. Per-target signs: Pan's own row for the target (30 of 85 h1-test targets covered)
On the covered targets' gate genes (`signs450.py`), unweighted / |lfc|-weighted:

| predictor | gate | gate-block |
|---|---:|---:|
| ours (cross-line composition, 323 `cap15`) | 0.551 / 0.572 | 0.551 / 0.577 |
| CLEAN4 programme | 0.578 / 0.573 | 0.586 / 0.582 |
| Pan programme | 0.589 / 0.580 | 0.621 / 0.605 |
| **Pan's own row for the target** | **0.764 / 0.810** | **0.765 / 0.819** |

Per target (26 with ≥ 10 gate genes): Pan's row 0.773 against ours 0.554; Pan better on 26 / 26. The leaders' inferred
gate sign accuracy (434f) is 0.62–0.66.

## 3. Per-target direction for the bulk (`bulk450.py`, the same 30 targets, 10,588 genes ≥ 10 control CPM, panel out)
| predictor | cosine | cosine, target-centred | pds-like (centred, among the 30) |
|---|---:|---:|---:|
| ours (cross-line) | 0.011 | 0.069 | 0.947 |
| Pan's row (same line) | 0.038 | **0.260** | 0.966 |
| unit(ours) + 0.5·unit(Pan) | 0.025 | 0.171 | **0.985** |
| unit(ours) + 1·unit(Pan) | 0.033 | 0.221 | 0.984 |
| unit(ours) + 2·unit(Pan) | 0.037 | 0.253 | 0.978 |

Pan's centred cosine beats ours on 30 / 30 targets (3.8× on average). Blends do best on discrimination (pds-like); pure
Pan does best on cosine.

## 4. Replication in a cancer line (`k562rep450.py`)
Truth: VIPerturb-seq K562 (Satija lab 2026; a different lab and technology; median 48 cells per perturbation). Predictors
are three banks, on 8,055 common targets and 7,539 common genes, each target's own gene out: Replogle K562-GWPS (same
line), Orion HCT116 and Orion HEK293T (the production donors).

| stratum (sample) | predictor | cosine, centred | sign on the truth's top-100 genes |
|---|---|---:|---:|
| K562-essential (1,294) | **K562-GWPS** | **0.062** | **0.565** |
| | HCT116 | 0.026 | 0.526 |
| | HEK293T | 0.022 | 0.523 |
| non-essential (1,500) | K562-GWPS | 0.008 | 0.514 |
| | HCT116 | 0.004 | 0.510 |

- Essential stratum: the same-line bank beats HCT116 on 956 / 1,294 targets, 2.4× in direction.
- Non-essential stratum: the 48-cell truth is noise-dominated for weak responses and cannot rank predictors. On its
  strongest quartile, K562-GWPS reads 0.012 vs 0.005.
- The effect replicates wherever the truth has signal.

## 5. Does an essential-only same-line bank carry the programme? (`esssplit450.py`)
Context A's only same-line bank, Nadig Jurkat, is an essential screen. On H1, Pan's programme recomputed from subsets of
its perturbations (common coverage with CLEAN4 0.5807 gate-block / 0.5244 fid-like):

| Pan perturbations used | n | gate-block | fid-like block |
|---|---:|---:|---:|
| K562-essential-list only | 1,284 | **0.6267** | 0.5388 |
| the rest | 1,040 | 0.6124 | 0.5370 |
| all | 2,324 | 0.6253 | 0.5390 |

The programme is carried by essential knockdowns at least as well. So an essential screen of the right line
(Nadig Jurkat for A, Nadig HepG2 for a HepG2 context) should carry it. One caveat for A: CLEAN4 already gives Nadig
Jurkat a quarter of its weight, so A's headroom is smaller than H1's, where CLEAN4 holds no H1 data.

## 6. When the recipient's line is already a donor (`mixrep450.py`)
K562 and HCT116 are production donors, so an October K562 or HCT116 context already gets its own bank at about a third
of the mix. Setup: VIPerturb K562 truth, essential stratum, 1,294 targets, centred cosine, simple unit-mean mixes (the
production mix is whitened, which may suppress line-specific signal more).

| predictor | centred cosine |
|---|---:|
| HCT116 + HEK293T (no same-line bank) | 0.032 |
| K562-GWPS + HCT116 + HEK293T (same line at 1/3) | 0.057 |
| K562-GWPS alone | 0.062 |
| mix + 1 × K562-GWPS | **0.066** |
| mix + 2 × K562-GWPS | 0.066 |

- The same-line blend adds about +15 % direction when the line is already in the mix.
- It adds the full effect (about 2×, and 3.8× on H1) when the bank is new information: HepG2, DLD-1, A549 and the other
  non-donor lines.
- The programme part (§1, §5) applies either way, because the shipped programme is CLEAN4, not the recipient's.

## 7. The bulk's shared shift is NOT same-line transferable; centre the bank rows (`bulkprog450.py`)
- In h1-test's truth (log pseudobulk CPM vs control, genes ≥ 10 CPM), the mean over the 85 targets carries **97 %** of
  the rows' energy. Our cross-line prediction's mean carries 1 %.
- That shared shift is platform- or batch-specific:
  - cos(truth mean, Pan's mean delta) = **−0.12**;
  - cos(truth mean, our mean) = +0.02.
  - It may be particular to the 2025 H1 data. On A/B/C the baseline (the mean perturbed profile) reads mse 0.986–0.992,
    so no comparable shared shift is visible there in the scorer's own space.
  - Either way, centring the bank rows is the safe choice.
- The same-line bank carries the DE-sign programme (§1) and the target-specific direction (§3). It does not carry the
  bulk shared shift.
- So the builder blends **programme-centred** bank rows: each row minus the bank's own mean over its non-panel
  perturbations. On the 30 covered targets:

| predictor | raw cosine | target-centred cosine |
|---|---:|---:|
| ours | 0.011 | 0.069 |
| Pan row, raw | 0.038 | 0.260 |
| Pan row minus Pan's mean | **0.063** | 0.260 |
| ours + 1 × raw row | 0.033 | 0.218 |
| ours + 1 × centred row | 0.050 | 0.217 |

  Per-target DE signs keep the raw row, because the programme is informative for signs.
  **Superseded by §9's emitter-level prices:** through the real emitter, RAW rows beat centred rows (+0.0806 vs
  +0.0741). Centred rows gain pds (+0.0010) but lose mse (+0.013 raw) and nmae (+0.010). `def_sameline.py` therefore
  stores raw rows by default; `--centre` is optional.

## 8. Estimator and near-line checks (`progest450.py`, `nearline450.py`)
- **Programme estimator.** On the genes CLEAN4 covers, gate-block accuracy by estimator:

  | estimator | gate-block |
  |---|---:|
  | P_down over all rows (396's rule) | 0.625 |
  | P_down over the strong half | 0.629 |
  | median delta | 0.626 |
  | mean delta | 0.623 |
  | mean of unit rows | 0.618 |
  | P_down over the weak half | 0.598 |

  Keep 396's rule. It is simple, known-answer-pinned, and within 0.004 of the best.
- **Near lines (iPSC banks for H1 hESC): INCONCLUSIVE, not negative.**
  - kolf2's programme reads 0.441 (below chance). Its matrix is per-gene z-scores (`zscore_filt_log1p`), which centre out
    exactly the shared shift a programme is.
  - Tian 2019 iPSC reads 0.516, but its own-target deltas are negative in only 61 % of rows (median −0.02), so its
    knockdowns are barely visible.
  - Neither is a fair near-line bank. Whether a same-tissue line carries the programme is still open; the same-line
    mode stays same-line only.

- **Near line, same tissue (`nearline_colon450.py`): no support either.**
  - Setup: Orion HCT116 truth; DLD-1 (colon, MSI-high, another lab) vs K562 and Orion HEK293T; centred cosine on the
    essential stratum.
  - DLD-1 reads 0.067, K562 0.087, HEK293T 0.107.
  - DLD-1's table is a sparse DE-lfc matrix (only significant entries non-zero), which penalises it. The best
    predictor of HCT116 is HEK293T, from the same lab and protocol (Xaira's Orion), so lab and protocol matter for
    cross-line transfer too.
  - Same tissue alone is not shown to help. The mode stays same-line only.

- **LINCS L1000 as a same-line bank (`cmapsl450.py`): NO.**
  - Setup: Nadig HepG2 CRISPRi Perturb-seq truth; 346 targets every CMap core line covers; landmark genes; CGS shRNA
    consensus signatures.
  - CMap HEPG2 (same line) reads centred cosine 0.027, below A549, MCF7, HT29, PC3 and A375 (0.028–0.072).
  - Per target, HEPG2 minus the mean of the others is −0.018 [−0.026, −0.011].
  - L1000's same-line data carries no same-line advantage, so the 17 CMap lines do not join the registry. The same-line
    effect needs a Perturb-seq-grade bank.

## 9. The emitter-level H1 price (the shipped harness, `run450.py` / `run450b.py` / `run450c.py`; paired vs 435's seeds)
Board units use 435's H1→board slopes. The rows blend Pan's own row into our effect for the 30 covered targets (of 85),
norm-preserving. "pt" means the covered targets' block prior is the row's own sign. The default is the prior for every
other gene and target. The known answers (the CLEAN4 prior through each wrapper, and w 0) reproduce 435 bit for bit.

| configuration | seeds | pds | mse | nmae | fid | reach | jac | overall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| programme prior only (all targets) | 3 | 0 | 0 | −0.0086 | +0.0096 | −0.0407 | −0.0001 | +0.0001 |
| rows w 0.5, CLEAN4 default | 1 | +0.0112 | −0.1287 | −0.0480 | +0.0106 | +0.0344 | +0.0017 | +0.0531 |
| rows w 1, CLEAN4 default | 3 | +0.0115 | −0.1744 | −0.0640 | +0.0144 | +0.0493 | +0.0027 | +0.0710 |
| rows w 2, CLEAN4 default | 1 | +0.0104 | −0.1944 | −0.0708 | +0.0164 | +0.0581 | +0.0033 | +0.0790 |
| rows w 4, CLEAN4 default | 1 | +0.0092 | −0.1963 | −0.0714 | +0.0171 | +0.0617 | +0.0035 | +0.0803 |
| rows w 8, CLEAN4 default | 1 | +0.0088 | −0.1937 | −0.0705 | +0.0171 | +0.0539 | +0.0034 | +0.0779 |
| rows w 32 (≈ replacement), CLEAN4 default | 1 | +0.0085 | −0.1899 | −0.0698 | +0.0173 | +0.0552 | +0.0034 | +0.0773 |
| rows w 1, programme default | 3 | +0.0115 | −0.1744 | −0.0706 | +0.0236 | +0.0104 | +0.0027 | +0.0706 |
| rows w 1 + pt, programme default | 3 | +0.0115 | −0.1744 | −0.0772 | +0.0332 | +0.0245 | +0.0027 | +0.0806 |
| rows w 1 + pt, CLEAN4 default | 3 | +0.0115 | −0.1744 | −0.0743 | +0.0267 | +0.0675 | +0.0018 | +0.0837 |
| rows w 2 + pt, programme default | 1 | +0.0104 | −0.1944 | −0.0836 | +0.0361 | +0.0336 | +0.0037 | +0.0892 |
| **rows w 2 + pt, CLEAN4 default** | 3 | +0.0102 | −0.1944 | −0.0806 | +0.0311 | +0.0761 | +0.0036 | **+0.0934 ± 0.0002** |
| rows w 4 + pt, CLEAN4 default | 3 | +0.0092 | −0.1963 | −0.0814 | +0.0301 | +0.0793 | +0.0028 | +0.0933 ± 0.0005 |
| rows w 1 centred + pt, programme default | 3 | +0.0125 | −0.1618 | −0.0670 | +0.0332 | +0.0143 | +0.0031 | +0.0741 |

- **The rows are the lever.** Blending the same-line rows into 35 % of the targets moves every member the right way:
  mse by −0.17 raw, nmae −0.064, reach +0.05, pds +0.011, fid +0.014.
- The per-target prior adds about +0.01 on top.
- The programme default trades fid/nmae for reach on the uncovered targets (the §9-top row, and CNN3/COL6A1-type
  single-gene coin flips). On H1 it is net negative inside the full mode (+0.0806 vs +0.0840).
- The October default is rows + pt, with the CLEAN4 default (`:pt:noprog`). The programme switches on only if 446bja's
  board read is positive.
- **Dose.** Rows alone peak at w 2–4, and near-replacement (w 32) is only 0.003 lower. For covered targets, the
  cross-line composition adds little once a same-line bank exists. The full mode is flat from w 2 to w 4; w 2 keeps more
  pds.
- **The programme default costs 0.003–0.004 at every dose** (reach), inside the full mode.
- **October default: `:2:pt:noprog`, +0.093 on H1 at 35 % coverage.** `scripts/def_day.sh` defaults to it. The
  programme is added only if 446bja reads positive.

## 10. The dose for a WEAKER bank (`drive450f.sh`, 2026-09-25 22:15Z)
Pan's rows are weakened by flipping the sign of a random fraction p of each row's measured genes. Magnitudes and
sparsity are kept; Gaussian noise instead broke the emitter's per-cell totals at w 1, a harness artefact. Full mode
(rows + per-target prior, CLEAN4 default), seed 0, board units:

| bank quality (gate signs) | w 0.5 | w 1 | w 2 |
|---|---:|---:|---:|
| clean Pan (0.764) | +0.0531 (rows only) | +0.0837 | **+0.0934** |
| p 0.2 (≈ 0.66) | **+0.0323** | +0.0316 | +0.0268 |
| p 0.3 (≈ 0.61) | **+0.0166** | +0.0086 | +0.0007 (mse +0.020 worse) |

- **The best w falls with bank quality**: 2 for a Pan-quality bank, 0.5–1 at 0.66, 0.5 at 0.61.
- The blend never hurts at w 0.5 down to a 0.61 bank. At w 2 a 0.61 bank is neutral overall and costs mse.
- So w should be per context. 451's bank audit puts the Jiang lines, VIPerturb and Myllia at 0.62–0.83 own-target
  knockdown, the 'weaker' class.
- Until the builder takes a per-context w, the day-of rule is: w 2 if every matched bank is strong (Replogle, Nadig,
  Orion, Sunshine), else w 0.5.

## Reading
1. **A same-line bank from another lab and chemistry carries the line's biology.** It carries the shared programme (most
   of the oracle's advantage on the block genes where our DE calls live), per-target signs far above anything the field
   shows, and the bulk direction.
2. **For A/B/C this can be tested on context A only, and only for the programme.** Nadig 2025 Jurkat is an essential
   screen and covers 0 of the essential-free panel, so no per-target rows. The builder already has the knob
   (`VCC_DE_SIGN_PRIOR=clean4q80ja`, 396's artifact: Nadig Jurkat's P_down for A, CLEAN4 for B/C), built on 09-16 and
   never scored. 370's within-Jurkat split measured 0.577 vs cross-line ≈ 0.55 on A's block. The arm is 446bja; its
   preregistration is written after the H1 emitter price below.
3. **For October this is the largest lever we have, if a D/E/F line has a bank.** "Six lines from different tissues of
   origin", and blood is taken by Jurkat. The likely tissues have banks: colon (HCT116 Orion genome-wide; DLD-1
   low-MOI), liver (HepG2, Nadig essential: programme only), lung (A549, Myllia/Altair), kidney (HEK293T Orion). A
   same-line mode needs:
   (a) the day's identity call (`def_identity.py`, already a driver stage);
   (b) a per-context programme prior (generalise 396's `psign_A` to any context and bank);
   (c) per-target same-line rows into that context's effect (a per-context union: the H1 union's post-whiten residual
       mechanism, which has a board point), dose priced on H1 by `run450b.py`.
   It has to be built and rehearsed before Oct 22 (DEF-RUNBOOK §3 row 220).
4. Rules: this is public Perturb-seq feeding a fitted transfer model, exactly like the K562/Orion donors we ship. It is not
   a hand-curated lookup and not the challenge's own data. The final write-up should frame the per-context weight as a
   basal-similarity-conditioned donor weight fitted on H1.

## Queued (H1 emitter prices; `drive450.sh`, `drive450b.sh`, serial after 446br2 and 448's lvl3)
- `run450.py`: the shipped harness with the block prior swapped to Pan's programme (seeds 0–2) and Pan-with-CLEAN4-fill.
  Known answer: the CLEAN4 prior through the wrapper equals 435's stats bit for bit.
- `run450b.py`: Pan's rows blended into the effect for the 30 covered targets at w 0.5 / 1 / 2, and w 1 with Pan's
  programme prior (seeds 0–2). Known answer: w 0 equals 435.

## 11. After 446bja's null: the programme belongs to the lab more than the line (2026-09-26 04:00Z)
**450h, `labline450.py`.** Per-gene P_down (def_sameline's rule, A/B/C panel rows out) for eight banks from three labs:
Weissman (K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2, RPE1), Xaira (Orion HCT116, HEK293T), Satija (VIPerturb
K562). Pairwise agreement of (P_down − 0.5):

| pair type | pairs | mean r | sign agreement (both \|P−0.5\| ≥ 0.05) |
|---|---:|---:|---:|
| same line, same lab (GWPS / K562-ess) | 1 | **+0.632** | 0.981 |
| **same line, different lab** (K562 Weissman / VIPerturb) | 2 | **+0.090** | 0.713 |
| **different line, same lab** (Weissman lines; HCT116 / HEK293T) | 10 | **+0.505** | 0.845 |
| different line, different lab | 15 | −0.031 | 0.512 |

- Among cancer and immortalised lines, a bank's programme is mostly its lab's and protocol's signature.
  - Five Weissman lines agree with each other (r 0.43–0.71).
  - Two K562 screens from different labs barely agree (r 0.09).
- This explains 446bja's null without any H1 contradiction. Nadig's Jurkat programme is largely the Weissman
  signature, which CLEAN4 (four Weissman banks) already carries.

**450i, `protocol450.py`: does a protocol-matched bank substitute? No.**
- Setup: the VCC H1 truth is 10x Flex, and VIPerturb is a Flex-style probe-based screen (in K562).
- h1-test gate-block sign accuracy of each programme, with the h1-test targets excluded from the bank rows:
  - VIPerturb K562 **0.539**;
  - Orion HEK293T 0.540, Orion HCT116 0.565;
  - Weissman lines 0.557–0.579;
  - CLEAN4 0.581;
  - **Pan H1 (same line, other lab, 10x 5′) 0.625**.
- On the 6,006 genes all four cover: CLEAN4 0.581 / Pan 0.625 / VIPerturb 0.539 / HCT116 0.561.
- **So the H1 transfer is a same-line effect that survived a lab change.** A strong ESC-specific shared response
  (knockdown-driven exit from pluripotency) plausibly dominates any lab signature there. For cancer lines no such
  line-specific programme stands out of the lab signature. A Flex-matched bank from another line does worse than the
  Weissman consensus.

**Consequences.**
- October keeps `noprog`, now with a mechanism.
- The rows part of the same-line mode is per-target direction, and 450 §4 replicated it across labs on K562 (2.4× in
  direction on essential targets). It is a different object from the programme and stays on.
- Keep CLEAN4 as the default block prior: it beats every single line, and every other lab's programme, on H1.

**450j, `multilab450.py`: a multi-lab consensus does not beat CLEAN4 by a usable margin.** h1-test gate-block accuracy on
CLEAN4's 16,697 gate-block pairs (per-lab z(P_down − 0.5), lab as the unit):

| consensus | accuracy |
|---|---:|
| CLEAN4 (Weissman) | 0.5799 |
| Xaira (HCT116 + HEK293T) | 0.5657 |
| Satija (VIPerturb) | 0.5384 |
| Weissman + Xaira | **0.5856** |
| Weissman + Satija | 0.5800 |
| all three labs | 0.5826 |

- The best is +0.006 over CLEAN4. At 450's H1 conversion (Pan's +0.044 gate-block netted +0.0001 overall), that is
  about 0 overall and inside noise.
- No arm. Averaging labs does not recover a line's programme. Only the same line did, on H1.

**450k, `labline_rows450.py --common`: per-target rows keep a line effect across labs (the programme does not).**
Setup: the same eight banks, on the 1,129 targets all eight screened and 6,083 shared genes (the A/B/C panel excluded).
Rows are centred by each bank's own mean row; the statistic is the mean per-target cosine (strongest quartile in brackets).
The pairwise-target run, without `--common`, is confounded by each pair's own target set: genome-wide pairs carry many
weak targets. The common set removes that.

| pair type | pairs | cos (top quartile) |
|---|---:|---:|
| same line, same lab (GWPS / K562-ess) | 1 | 0.216 (0.356) |
| same line, different lab (VIPerturb / the two Weissman K562 screens) | 2 | **0.057 (0.109)** |
| different line, same lab | 10 | 0.129 (0.189) |
| different line, different lab | 15 | 0.065 (0.120) |

- Within VIPerturb's own row, which holds its noise fixed, the same-line partners read 0.066 / 0.048. Its five
  other-line partners read 0.029–0.040 (mean 0.033). **The line effect survives the lab change at ~1.7×.**
- Within the Weissman lab, same line vs other lines: 0.216 vs 0.129, also ~1.7×. Lab matters too: same-lab other-line
  pairs (0.129) beat cross-lab other-line pairs (0.065).
- Contrast with the programme (450h): there, same line across labs reads r 0.09 against 0.51 for other lines in one lab.
  The line effect vanished.
- For October: the same-line **rows** carry line information across labs (as 450 §4's VIPerturb-truth replicate showed,
  2.4×). The same-line **programme** does not. The mode's design (rows on, `noprog`) matches the data.

**450l, `arcprog450.py`: Arc's own programme.** P_down from Arc's 2025 H1 screen (`h1_2025`: 255 perturbations, Flex, 400
cells each; A/B/C panel rows out) against the other programmes:

| vs | r | sign agreement |
|---|---:|---:|
| Pan H1 (same line, Whitehead) | **+0.261** | 0.755 |
| K562-ess / RPE1 / K562-GWPS / HepG2 / Jurkat (Weissman) | +0.20 / +0.20 / +0.17 / +0.16 / +0.15 | 0.62–0.72 |
| HEK293T / HCT116 (Xaira) | +0.15 / +0.13 | 0.66–0.67 |
| VIPerturb K562 (Satija, Flex-style) | +0.008 | 0.533 |

- Arc's programme shares a modest component with both Weissman and Xaira (r 0.13–0.20), although those two labs do not
  agree with each other (r ≈ 0).
- It shares nothing with the only other Flex-style screen. So "Flex chemistry" is not the shared factor. A common
  CRISPRi construct lineage is a candidate: the 2025 challenge library derives from Horlbeck 2016 (Weissman). So is a
  generic proliferation response.
- **Hypothesis for a later arm (458, not built): Arc's own H1 programme as the block prior for all three contexts.**
  - If programmes are lab signatures, the recipient's lab is the right source, and it would carry to D/E/F.
  - Against it: ESC-specific content (r 0.26 with Pan's H1). And the prior-swap route's cap: 450's full-panel H1 price
    of a much better programme netted ~0 overall because reach pays on small gates.
  - Needs a prior artifact and a builder registry entry like 396's `clean4q80ja`, plus a known answer. Price it before
    spending a slot.

**450h v2, `labline450_v2.py` (07:55Z): four labs, fifteen banks.** It adds Jiang 2025's six lines (K562, A549, MCF7, HT29,
HAP1, BxPC3; 218-target pathway screens) and a second HEK293-family screen (xucao_hek293).

| pair type | pairs | mean r of (P_down − 0.5) |
|---|---:|---:|
| same line, same lab | 1 | +0.632 |
| **same line, different lab** (K562 × 4 labs; HEK293 × 2 labs) | **6** | **−0.020** (range −0.283 … +0.092) |
| different line, same lab | 25 | +0.291 (Weissman +0.506 ×9; Xaira +0.499 ×1; Jiang +0.148 ×15) |
| different line, different lab | 73 | −0.003 |

- Every lab with more than one line shows positive cross-line agreement. The noisy, small Jiang screens still reach +0.15.
- No same-line pair across labs clears +0.1. The two HEK293-family screens from different labs anti-correlate (−0.28).
- The programme is a lab and protocol signature, now across four labs and two lines with cross-lab pairs.

**450h v3, `labline450_v3.py` (09:30Z): CORRECTION. With the perturbation set controlled, the lab dominates but the line is
not nil.** v1/v2 computed each bank's programme over its whole screen. Banks with different perturbation sets were then
compared: essential-only vs genome-wide vs 218-gene pathway screens. v3 computes both programmes of a pair over only the
perturbations both screened (A/B/C panel excluded).

| pair type (pairs sharing ≥ 1,000 perturbations) | pairs | mean r |
|---|---:|---:|
| same line, same lab (GWPS / K562-ess) | 1 | +0.806 |
| same line, different lab (K562-ess / VIPerturb 0.291 on 1,443 essential; GWPS / VIPerturb 0.118 on 8,358 genome-wide) | 2 | +0.205 |
| different line, same lab (Weissman ×9, Xaira ×1) | 10 | **+0.544** |
| different line, different lab | 15 | +0.175 |

- **Lab or protocol is the larger factor.** Across lines, same-lab pairs sit about 0.37 above cross-lab pairs.
  - Weissman essential screens vs Xaira on shared essential perturbations: 0.24–0.40.
  - Pairs involving VIPerturb: −0.14 … +0.12. Its programme is idiosyncratic.
- **The line adds a smaller but real amount.** Inside VIPerturb's own row, with matched set types:
  - essential sets: same-line K562-ess 0.291 against other-line Jurkat / RPE1 / HepG2 at 0.122 / 0.060 / 0.042;
  - genome-wide sets: same-line GWPS 0.118 against HCT116 / HEK293T at −0.100 / −0.141.
  - So the same line adds about +0.2 in r.
- Over all pairs sharing ≥ 100 perturbations (the Jiang pathway screens and the HEK293 pair included; median ~200 shared):
  same line, different lab +0.014 (5 pairs); different line, same lab +0.307 (25); different line, different lab +0.057 (44).
  The small sets are noisy.
- **v2's "same line across labs r −0.02" is withdrawn as a headline.** It mixed unmatched perturbation sets. v3 is the
  controlled statement: the programme is predominantly a lab and protocol signature, with a smaller line component (~+0.2).
- Consistency with 446bja's null: CLEAN4 already averages four Weissman banks, which carry Weissman's shared signature.
  Moving A to Nadig Jurkat alone adds the line component (~+0.2 r) but gives up the averaging. A net of ≈ 0 is what
  the board read.

**450l on shared sets (09:45Z).** Arc-H1's programme against each bank, both computed over their shared perturbations
(57–227):
- Weissman: GWPS 0.150, K562-ess 0.183, Jurkat 0.178, HepG2 0.152, RPE1 0.249.
- Xaira: HCT116 0.157, HEK293T 0.150.
- Satija: VIPerturb 0.103.

This is the generic cross-lab level (~0.18 on matched sets), aligned with no one lab. So Arc behaves like a separate lab.
If its within-lab agreement matches the others' (~0.5), Arc's own programme is the best available prior for Arc's lines.
That is 458's test.

# 622 — The identity sign source in its board / day-of form: Flex-only references, stable, key-gated

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** Here every HCT116 / HEK293T (Orion) row here is VOID. That includes the identity-sign gains on those folds and the hook rows' orion_hct116 / orion_hek293t lines. The claim "public 3′ references alone are harmful" rested on HEK293T only and is WITHDRAWN. The H1, K562 and Broad-16 readings (622m–p, including 622n's low-identity UP tail and 622o's idall3 stress tail) stand.

2026-10-02, council session [7772dc], identity track (handed over by [d4cd5f] after 639). Plumbing (prior files, builder
knob, def_day patch + KA) is [4d0207]'s. This experiment is the science that sets the score. Analysis only: no build, no
`vcc`, no `cell-eval2`. Every script < 3 GB, seconds. Scripts: `stab622.py`, `refs622.py`, `hyb622.py` (run with cwd
`experiments/480-nmae-instrument`). Outputs: `work/stab622.json`, `work/refs622.json`, `work/hyb622.json`, and the A/B/C
control CPM cache `work/abc_cpm622.npz`.

Score (639): S = z(P_C4 − 0.5) + λ·z(I), positive → predict down. Identity: I = log1p(ctrl CPM of the context) − mean
over reference profiles of log1p(ctrl CPM).

## (b) Stability inside A/B/C: disjoint NTC-guide halves of each context's controls (stab622)

| context | genes > 5 CPM | Spearman(I half 1, half 2) | top-decile overlap | S sign agreement λ 0.5 / λ 1 | genes whose sign flips vs CLEAN4, λ 0.5 / λ 1 |
|---|---:|---:|---:|---|---|
| A | 9,929 | 0.9997 | 0.989 | 0.998 / 0.997 | 9.7 % / 17.2 % |
| B | 9,626 | 0.9994 | 0.996 | 0.997 / 0.995 | 8.7 % / 16.4 % |
| C | 10,124 | 0.9997 | 0.992 | 0.997 / 0.994 | 11.3 % / 20.8 % |

I is a line property, not control sampling noise. 639's control-halves check on h1 said the same about the truth
side.

## (c) Which references: gate sign accuracy vs CLEAN4 (LOLO) on 636's gates (refs622)
Δ at λ 0.5 / λ 1, with target-bootstrap 95 % CIs in `work/refs622.json`:

| surface | mixed (639: every other profile) | **Flex only** (A, B, C [+ h1]) | public 3′ only |
|---|---|---|---|
| h1 (Flex) | +0.014 / +0.016 | +0.009 / +0.014 (CIs touch 0; 3 refs) | +0.006 / +0.004 |
| hct116 | +0.039 / +0.015 | +0.031 / +0.015 | +0.004 / −0.008 |
| hek293t | +0.017 / +0.024 | +0.027 / **+0.042** | **−0.051 / −0.085** |
| k562 (low power) | −0.005 / −0.010 | +0.001 / −0.001 | −0.002 / −0.019 |

- **Public 3′ references alone are harmful** (hek293t −0.05 … −0.085). The A/B/C Flex profiles are what make I work.
  Flex-only references match or beat the mixed set on the 3′ surfaces. On h1 they keep most of the gain with only 3
  references.
- Leave-one-reference-out of the mixed set (`refs622.json`, "mixed-minus-*"): every variant stays positive on h1,
  hct116 and hek293t (+0.008 … +0.039), most CIs exclude 0, and none flips sign. The smallest h1 value comes from
  dropping C (CAL33): +0.008 / +0.011.
- Worst case over h1/hct116/hek293t with Flex references: λ 1 +0.014, λ 0.5 +0.009.

## (a) Key-gated HYBRID + identity on 612's folds: nmae, MDE-law magnitude, only signs differ (hyb622)
612's rule: the key sign where |key| is in the target's top half, the programme in the low half. Here the low half
uses sign(S) instead of sign(CLEAN4); identity alone where CLEAN4 is undefined. Δ is mean(ref − X), positive = better,
target bootstrap.

| fold | vs 612's hybrid (CLEAN4 in the low half): Flex λ 0.5 | Flex λ 1 | mixed λ 0.5 | vs key-only: Flex λ 0.5 |
|---|---|---|---|---|
| h1 | **+0.0048 [+0.0016, +0.0086]** | +0.0048 [+0.0014, +0.0087] | +0.0015 [−0.0008, +0.0038] | +0.0057 [+0.0015, +0.0099] |
| orion_hct116 | **+0.0031 [+0.0011, +0.0050]** | +0.0003 [−0.0019, +0.0026] | +0.0050 [+0.0034, +0.0066] | +0.0210 |
| orion_hek293t | **+0.0040 [+0.0019, +0.0063]** | +0.0070 [+0.0044, +0.0097] | +0.0050 [+0.0033, +0.0068] | +0.0203 |
| k562_gwps (low power) | −0.0013 [−0.0031, +0.0006] | −0.0010 [−0.0034, +0.0013] | −0.0018 [−0.0037, +0.0000] | +0.0003 |

- Programme-only signs (S with no key half) LOSE on h1 and k562 (−0.003 … −0.020 vs key-only). **Keep the key half.**
- Flex λ 0.5 is the min-max-regret choice on nmae: positive with CIs excluding 0 on 3/4 folds, and the worst case is
  k562's −0.0013 with a CI spanning 0. Against key-only it is ≥ 0 on all four folds.

## Recommendation (sent to [4d0207] and [d4cd5f])
- **Score:** S = z(P_C4 − 0.5) + λ·z(I), positive → down, with λ ∈ {0.5, 1} as the preregistered 2-point dose. λ 0.5
  is preferred on nmae; λ 1 has the better worst case on gate-sign accuracy. Identity alone where CLEAN4 is undefined.
- **I uses Flex references only, never public 3′ profiles:**
  - board (A/B/C): each context vs the other two A/B/C contexts (+ H1's Flex controls);
  - day-of (D/E/F): each context vs every other Flex context (A, B, C, the other two of D/E/F, + H1).
  - z over genes that all reference profiles measure.
- **Delivery, key-gated as in 612:** the key's sign where |key| is in the target's top half of nonzero |key|, S's sign
  in the low half (and where the key is zero). For the block prior, S replaces CLEAN4 only where today's design does
  NOT protect the gene (q80).
- **Class:** transfer. Nothing is fitted to A/B/C. I is a function of the round's own controls, so it re-fits on the day.
- **Expected size:** on nmae it is the identity's increment over 612's hybrid, +0.003 … +0.005 raw on three of four
  folds (≈ +0.001 overall at H1's magnitude multiplier 1.2–1.8; more if a rank-neutral UP delivery raises the
  magnitudes, as in 632). On gate-sign accuracy it is +0.009 … +0.042. A prior-only board arm reads it exactly
  (shared draws).

## (d) More references: Satija's Flex profiles (`morerefs622.py`, `work/morerefs622.json`)
Adding 597's three Satija Flex pseudobulks (K562, HAP1, HEK293; Flex v1, another lab) to the Flex references is mixed:
h1 ≈ same (+0.0088 vs +0.0089 at λ 0.5), hct116 lower (+0.026 vs +0.031), hek293t higher (+0.036 vs +0.027), and
k562 worse (−0.010 vs +0.001). **A reference of the SAME line as the context removes its identity signal** (sat_K562
on the k562 surface). **Do not add other labs' Flex profiles.** Keep Arc's own Flex contexts, and never include a
reference whose line matches the context.

## (e) Fewer references: the board form has only 2–3 (`fewrefs622.py`, `work/morerefs622b.json`)
Gate-sign Δ vs CLEAN4 at λ 0.5 / λ 1:

| refs | hct116 | hek293t |
|---|---|---|
| {B, C} | +0.033 / +0.020 | +0.024 / +0.030 |
| {B, C, h1} | +0.032 / +0.014 | +0.026 / +0.039 |
| {A, C, h1} | +0.032 / +0.010 | +0.026 / +0.055 |
| {A, B, h1} | +0.034 / +0.017 | +0.025 / +0.030 |
| {A, B, C} | +0.031 / +0.019 | +0.022 / +0.029 |

- Every λ 0.5 variant has a CI excluding 0. k562 is null throughout (low power). With 2–3 Arc Flex references the
  identity works as well as with four, and including H1 makes little difference.
- **The board arm's exact form (each A/B/C context vs the other two, ± H1) is supported.**

## (f) Does |I| predict DE MEMBERSHIP? (`member622.py`, `work/member622.json`)
Pooled AUROC for gate membership, Flex references:

| surface | control log-CPM | \|I\| | rank-mix | \|I\| within expression deciles (mean) |
|---|---:|---:|---:|---:|
| h1 (full-power Flex) | 0.611 | 0.571 | **0.631** | **0.582** |
| hct116 | 0.814 | 0.593 | 0.752 | 0.514 |
| hek293t | 0.795 | 0.581 | 0.737 | 0.529 |
| k562 | 0.716 | 0.570 | 0.681 | 0.540 |

On the board-like surface (h1), identity strength carries membership information beyond expression. On the thinned 3′
surfaces it does not. This is a footnote, not a lever: 621b showed that at K 2,600 every ranking reaches φ ≈ 0.07.

## (g) Independent confirmation on the Broad/DepMap 16-line screen (`broad622.py`, `work/broad622.json`)
One lab, 90 Cas9-KO targets per line, 3′ chemistry. That is a different lab, modality and chemistry from every surface
above, never tuned on by this work.
- Only the 12 working lines are read. The 4 reserved lines (JVE127, MG63, SLR23, UMRC3) are dropped by `usecols` before
  any value is parsed.
- Gates are the authors' fdr < 0.05.
- I_L is taken against the other 11 working lines' OR-control pseudobulks: same-lab references, the analogue of Arc-only
  Flex references.

| source (gate-sign accuracy, Δ vs CLEAN4 on shared genes) | mean Δ over 12 lines | lines CI > 0 | lines CI < 0 |
|---|---:|---:|---:|
| identity alone | −0.146 (acc 0.46–0.61; > 0.5 on 10/12 lines) | 0 | 6 |
| all-down | −0.240 (these KO gates lean UP) | 0 | 12 |
| **z(C4) + 0.5·z(I)** | **+0.009** (−0.012 … +0.036) | 2 | **0** |
| z(C4) + 1·z(I) | +0.003 (−0.029 … +0.057) | 1 | 1 (T3M5) |
| same-lab LOLO programme (reference only; not available for D/E/F) | +0.087 | 11 | 0 |

- The principle generalises weakly: identity alone beats a coin flip on 10/12 lines. As a rider on CLEAN4 it is mildly
  positive at λ 0.5 and never significantly harmful there. At λ 1 one line loses significantly.
- **λ 0.5 is the dose to ship.** λ 1 is a second dose point only. The effect here is smaller than on 636's surfaces
  (+0.009 … +0.04), so expect the board read to sit at the low end.
- The same-lab programme row shows what a lab's own multi-line truth would buy (+0.087). For D/E/F that would be A/B/C's
  truth, which we do not have. This is the lab-signature effect 450h described.

## (h) H-633's decisive test: does a SAME-LAB programme from only 3 lines transfer? (`lab3_622.py`, `work/lab3_622*.json`)
This is the A/B/C → D/E/F analogue, run on the Broad-16's 12 working lines (reserved 4 never read).
- For each held-out line L and 40 random triples of other lines, the programme is P_down(g) = the fraction of the
  triple's gated responses that go down.
- The coarse version is the per-module majority sign over 16 truth-free modules (CLEAN4 octile × expression half),
  which is roughly what a ~16-arm board readout could decode.
- Gate-sign accuracy is read paired against CLEAN4.

| programme | mean Δ vs CLEAN4 (12 lines) | all triple sources | cross-batch triples only |
|---|---|---:|---:|
| full 11-line same-lab | +0.087 (+0.082 cross-batch) | 11/12 lines CI > 0 | |
| **3-line same-lab, gene-level** | **+0.054** | positive on 12/12; 5–95 % of triples > 0 on 8/12 | **+0.054** |
| 3-line, 16-module majority sign | +0.020 | | +0.021 |
| CLEAN4 with modules flipped where the triple disagrees | +0.016 | | +0.017 |

- **A lab's programme learned from three lines transfers to new lines of the same screen at +0.054 gate-sign.** That
  is about 2/3 of the 11-line value, and it is not a batch artifact (the cross-batch triples give the same number).
  The prior H-633 used (466: same-lab cross-line agreement 0.17 → ~+0.005 fid-like) is too pessimistic for one
  screen campaign. 458d's harm (Arc's 2025 H1 programme on 2026 cancer lines) is a different situation: another year,
  another screen, an ESC line.
- **A coarse 16-module readout keeps ~40 % (+0.020).** That is ≈ +0.003 … +0.004 overall on D/E/F by 635's exchange
  rate (+0.045 gate-sign ≈ +0.008 … +0.010), against ~16 prior-only A/B/C slots.
- Rule check (AGENTS 11.1c rule 3): mining A/B/C is barred unless a held-out test shows transfer, and this is that
  test. Whether to spend slots on a module readout, and whether that is within the spirit of the rules, is the
  coordinator's and Ahmet's decision (H-631's governance applies by analogy).

## (i) Is identity redundant once a same-lab programme exists? (`combo622.py`, `work/combo622.json`)
Broad-16, held-out line L, 40 random 3-line programmes per line, gate-sign Δ against the 3-line programme alone:
- + 0.5·z(I): +0.0011 (positive on 7/12 lines);
- + 1·z(I): −0.0037 (3/12);
- + 0.5·z(C4): +0.0034 (9/12);
- + 0.5·z(C4) + 0.5·z(I): +0.0046 (9/12).

**Identity mostly substitutes for a missing lab programme.** That is October's situation, since we cannot see A/B/C's
truth. Once a same-lab programme exists, identity adds ≈ 0. If a board readout of A/B/C's programme is ever done
(622h), re-price identity on top of it rather than assuming its full value.

## (j) Re-check in the DELIVERY-FAITHFUL shape (`hyb622b.py`, `work/hyb622b.json`)
[21dd4c] pointed out that kgh612/hyb622's fold model puts sign × MDE on a ZERO base, so the hook source sets every
uncalled gene's sign. In the real emitter the hook only sets the direction of a small MDE add on top of a key-signed base
lfc. 480's `small_base+mde` shape, P = BETA·cons + s·add (BETA 0.5025), matches the real delivery.

Δ nmae, positive = better, Flex references:

| fold | kgI λ0.5 vs no hook | kgI λ0.5 vs key-gated CLEAN4 | kgI λ0.5 vs plain CLEAN4 hook | plain identity hook (no gate) vs plain CLEAN4 |
|---|---|---|---|---|
| h1 | +0.0057 [+0.0015, +0.0099] | **+0.0058 [+0.0016, +0.0108]** | +0.0064 | +0.0007 (ns) |
| hct116 | +0.0210 | **+0.0068 [+0.0041, +0.0096]** | +0.0061 | +0.0061 |
| hek293t | +0.0203 | **+0.0075 [+0.0046, +0.0104]** | +0.0055 | +0.0032 (ns) |
| k562 | +0.0003 | −0.0010 (ns) | +0.0105 | **−0.0076** |

- **Independent replication.** These numbers equal [21dd4c]'s own `612/kgh612c.py` run (kgID vs no hook +0.0057 /
  +0.0210 / +0.0203 / +0.0003; kgID − kgC4 +0.0058 / +0.0068 / +0.0075 / −0.0010) to the digits printed, so two
  implementations agree.
- **The identity rider survives real-delivery modelling, and its increment grows**: +0.006 … +0.0075 over key-gated
  CLEAN4 on 3/4 folds, against +0.003 … +0.005 in the zero-base model. It dominates the plain CLEAN4 hook everywhere.
- **The key gate is mandatory**: without it the identity hook loses on k562 (−0.0076). λ 1 behaves like λ 0.5 (better on
  hek293t, worse on hct116); λ 0.5 stays the ship point.

## (k) Independent review of [4d0207]'s 604 artifact (`review604.py`, `work/review604.json`)
Each context's identity-augmented prior is recomputed with 622's code path from 604's pinned Flex profiles, using 604's
references (the other two contexts + H1 Flex), λ 0.5, and a quantile map onto CLEAN4's P.

| context | genes (mine / 604) | sign agreement | Spearman | max \|ΔP\| on CLEAN4-covered genes | genes flipped vs CLEAN4 |
|---|---|---:|---:|---:|---:|
| A | 10,962 / 10,962 | 1.00000 | 1.000000 | **0** | 8.8 % / 8.8 % |
| B | 10,962 / 10,962 | 1.00000 | 1.000000 | **0** | 9.4 % / 9.4 % |
| C | 10,962 / 10,962 | 1.00000 | 1.000000 | **0** | 10.2 % / 10.2 % |

- **604's ident604_flex_lam0.5_ABC.npz reproduces 622's science exactly.**
- The only difference is upstream and it is mine. 622's cached `abc_cpm622.npz` accumulated per-cell CPM in float32
  (row sums 1,000,000.06 … 1,000,000.44; relative error 10⁻⁷ … 10⁻⁵ against 604's float64 profiles, which sum to 10⁶
  exactly). It cannot move a rank or a z-score, and no 622 conclusion changes. 604's profiles are the precise ones.
  The "~1e-12" known-answer hint 622 sent to 4d0207 was too strict at the CPM level; the P-level equality above is the
  right KA.

## (l) Where the gain lives (`tercile622.py`, `work/tercile622.json`)
By target gate-size tercile (truth-derived; interpretation only), Δ gate-sign of z(C4) + 0.5·z(I) vs CLEAN4:

| surface | weak T1 | mid T2 | strong T3 |
|---|---:|---:|---:|
| hct116 | +0.084 | +0.040 | +0.019 |
| hek293t | +0.074 | +0.049 | +0.019 |
| k562 | +0.030 | −0.013 | +0.002 |
| h1 | +0.000 (per-target −0.023; the weak targets' few gate genes are noisy) | +0.007 | +0.009 |

On the 3′ surfaces the identity helps weak targets most, where the generic programme dominates their few DE genes. On
h1 the help is on mid and strong targets. This is not decision-changing.

## CORRECTION to (j) (2026-10-02 ~13:00Z, from [21dd4c]'s 612 CORRECTION 2, 22691ce5): the delivery claims are WITHDRAWN
The "delivery-faithful" fold shape (P = BETA·cons + s·add on EVERY gene) is mis-specified for hook decisions.
- Production adds the hook's sign only on UNCALLED genes (rc_emit_k.py:1059, `use_s = ~head & ~blk & (asg != 0)`).
  The 400-gene head and the 2,200-gene CPM block are steered by the key and the block prior, not by the hook.
- In a call-set-aware model (612/ds612c.py), the identity hook vs no hook reads:

  | fold | zero/all-gene model (j) | call-set-aware |
  |---|---:|---:|
  | h1 | +0.0057 | +0.0059 |
  | hct116 | +0.0210 | +0.0013 |
  | hek293t | +0.0203 | +0.0022 |
  | k562 | +0.0003 | −0.0021 |

- Most of (j)'s Orion gain sat on weak-key CPM-block genes, which the hook never touches. The call-set-aware model also
  fails a pre-committed H1 real-emitter retrodiction, so **no fold model is validated for hook-class levers**.
- **Withdrawn:**
  - "the identity rider survives real-delivery modelling and its increment grows";
  - "dominates the plain CLEAN4 hook everywhere";
  - "the key gate is mandatory" (the k562 −0.0076).
  Our independent replication reproduced the same mis-specified model, so it is not evidence.
- **Unaffected:** every gate-gene SIGN-ACCURACY result in (a)–(i) (refs622, broad622, lab3_622, combo622, member622,
  tercile622). These measure signs, not delivered nmae. The identity source's sign quality stands.
- **The delivered value is open:**
  - through the HOOK it is small on H1's real emitter (21dd4c's idkg seeds);
  - through the BLOCK prior (604's identity-augmented prior on the CPM block) it needs a real-emitter / prior-only
    board read.

## CORRECTION 2 to the recommendation (2026-10-02 ~13:12Z, from [d4cd5f]'s 639d route split)
- **Withdrawn: "identity alone where CLEAN4 is undefined".** In CLEAN4's gaps, identity-only signs are right only
  0.27–0.34 on HEK293T (head / block / uncalled) and 0.49 on HCT116's uncalled genes. The score must apply ONLY on
  CLEAN4-covered genes and leave every other gene as today. 639d sent this correction to [4d0207] for 604's spec.
  622's own numbers in (c)–(f) were all read on genes where both sources exist, so they stand.
- **The delivery route is the BLOCK prior.** On the cancer folds the identity gain sits on the CPM block (≥ 0.83 of the
  |l|-weighted net on HCT116 / HEK293T). On H1 (ESC) it sits on uncalled genes. A–F are cancer lines, so 604's block-prior
  arm is the primary delivery and the hook is secondary. This agrees with CORRECTION (j): the hook barely moves it.
- The ship point stays λ 0.5. λ 1 is the dose partner on the block only.

## (m) The block-restricted tail on the Broad-16 (2026-10-02 ~13:30Z; `block622.py`, `work/block622_{gated,raw}.json`)
Asked by [d4cd5f] for 639's block prereg (the hostile-line tail, [21dd4c]'s review item 4), plus the coordinator's lineage
condition. The setup is 622g's: 12 working lines (the 4 reserved lines are never parsed), authors' fdr < 0.05 gates,
same-lab references, λ 0.5, and no identity-only extension (only CLEAN4-covered genes are paired).

**Mask (approximate production mask, labelled as such).** Per (line, target):
- head = 400 genes by |consensus lfc| of that target over the other 11 lines;
- block = top 2,200 by the held-out line's control CPM after removing the head and the target;
- uncalled = the rest.
Two key analogues were read:
- `gated` (primary): the consensus of the other lines' DE calls (lfc · [fdr < 0.05]), the closest analogue of the donor
  key;
- `raw`: the mean raw lfc. It picks noisy low-count genes, so its head barely meets the gate (16–71 gate pairs per line).
Every target had ≥ 2 other lines, so every target has a head.

| region (gated head) | mean Δ vs CLEAN4 | lines negative | min | lines CI < 0 / CI > 0 | pooled share of the net gain |
|---|---:|---:|---:|---:|---:|
| **block** | **−0.0022** | **7/12** | **−0.0265 (T3M5)** | 2 (T3M5, HKGZCC) / 1 | −0.11 |
| head | +0.0133 | 3/12 | −0.0242 (C4I) | 0 / 3 | +0.28 |
| uncalled | **+0.0204** | **0/12** | +0.0000 (HCT15, SNU8) | 0 / 3 | **+0.83** |
| all gate genes (622g) | +0.0081 | 1/12 | −0.0011 | 0 / 1 | 1 |

The raw head gives the same picture: block −0.0024 (7/12 negative, min −0.0240 T3M5), uncalled +0.0221 (1/12
negative), pooled shares block −0.15 / head −0.02 / uncalled +1.16.

**By lineage (DepMap OncotreeLineage; gated head, block region).** The working set has **no blood lines**, so it cannot test
the K562 / Jurkat question.

| group | lines | block mean | block negative | uncalled mean |
|---|---|---:|---:|---:|
| epithelial carcinoma | HCT15 bowel −0.002, KYSE140 −0.006, KYSE450 +0.014 (oesophagus), SKGII −0.014, C4I +0.002 (cervix), T3M5 head & neck −0.027, KMRC20 kidney +0.006, SNU8 ovary −0.020, HKGZCC biliary −0.015 | −0.0069 | 6/9 | +0.0187 (0/9 negative) |
| other | RPMI7951 melanoma +0.019, U343 glioma +0.019, ONE58 mesothelioma −0.002 | +0.0119 | 1/3 | +0.0255 |
| blood | none | | | |

**Reading.**
1. **On these 12 cancer lines the identity gain does NOT sit on the CPM block.** The block is null-to-negative (7/12
   negative). The gain sits on uncalled genes (+0.020, 0/12 negative) and partly on the head.
   This is H1's and K562's pattern in 639d, not HCT116's / HEK293T's.
   The negatives cluster in EPITHELIAL carcinomas (6/9), including the two lines most like A/B/C's HeLa and CAL33 (cervix
   SKGII −0.014, head & neck T3M5 −0.027, the minimum). They are not a blood effect.
2. **CORRECTION 2's route claim ("on cancer lines the gain sits on the block") rests on the two Orion lines.**
   Measured so far: block-positive on 2 lines (HCT116, HEK293T, one lab); block-null or negative on 13 (12 Broad + K562).
   The route is UNRESOLVED. The block-form arm's held-out tail is "mean −0.002, 7/12 negative, min −0.027, 2/12
   significantly negative".
   The real-emitter H1 runs (625d) and a prior-only board read decide the route, not this table.
3. **The gain is identity, not expression level** (partial-out, gated head, all gate genes):
   - z(C4) + 0.5·z(E) with E = log1p CPM: +0.0011 (6/12 negative);
   - identity rider: +0.0081 (1/12 negative);
   - identity with expression regressed out: +0.0052 (1/12 negative).
   Expression alone predicts "down" at 0.49–0.62 on these gates, but CLEAN4 already carries that.
4. Caveats:
   - Cas9 KO, 3′ chemistry; the gates are the authors' fdr, not the scorer's Wilcoxon / BH.
   - The mask is an approximation: production's head comes from the donor key, not same-screen truth from other lines.
   - The control-halves check (`../627-tahoe-identity/halves622g.py`) leaves 622g's all-gate numbers unchanged (rider Δ
     +0.0093 with halves vs +0.0091 full; identity alone 0.5436 vs 0.5440), so none of this is shared-control noise.

## (n) Identity's TAILS on uncalled gate genes: a UP-side signal for 630's pass; no DOWN pass (2026-10-02 ~14:05Z; `tails622.py`, `work/tails622.json`)
Question: does identity concentrate at the extremes of I on UNCALLED gate genes? 609 (4d0207) and 622m put the gain
there. If so, 630's UP pass (UP side) or a top-I DOWN pass could deliver it with real magnitude.

**H1, real production emitter** (625d `none` s0; called = q < 0.05; I against A/B/C's Flex controls):

| uncalled gate genes | f 1 | f 0.35 | today's delivered sign accuracy | mean \|l_real\| |
|---|---:|---:|---:|---:|
| bottom 10 % of I → truth UP | 0.609 (4,661 pairs) | 0.711 (1,251) | 0.52–0.53 | |
| top 10 % of I → truth DOWN | 0.597 (6,966) | 0.578 (2,904) | 0.53–0.55 | 0.53 / 0.76 |

**Broad-16, approximate production mask** (622m's gated head + 2,200 block; the rest is "uncalled"):

| | value (mean of 12 lines) |
|---|---:|
| base rate of DOWN on uncalled gate genes | 0.267 (the KO gates lean UP) |
| bottom 10 % of I → UP | **0.781**, above 0.5 on 12/12 lines; CLEAN4 0.621 on the same pairs |
| top 10 % of I → DOWN | **0.296**, above the base rate on 7/12 lines only; still mostly UP |

**Reading.**
1. **The low-I tail goes UP on both instruments**: H1 0.61–0.71 against today's 0.52, Broad 0.78. Low-I genes are those a
   line expresses below other lines. 630's UP pass is the route that gives this real magnitude. It is the mechanism
   behind 630 × KGID compounding (nmae −0.0105 vs −0.0080).
2. **No identity DOWN pass.** The top-I tail goes DOWN on H1 (0.58–0.60) but NOT on the Broad-16 KO lines (0.30). That
   matches the DOWN harm in 603 / 632. Base rates of DOWN on uncalled genes differ by modality (H1 ~0.5, Broad KO 0.27),
   so a rule keyed to absolute tails would not transfer. The z-score rider's gain (622m, +0.020 on Broad uncalled genes)
   comes from modulating the sign around each surface's base, not from a tail rule.
3. Sent to [d4cd5f] (630) together with 625e's delivery calibration.
4. **Expression-only control (4d0207's 609b, 7b419f72).** On clean H1 the low-I → UP effect is identity, not expression.
   I's bottom decile gains +0.082 over CLEAN4 (n 2,562), while the expression-only bottom decile is flat (0.580, equal
   to CLEAN4).
   - On the Orion folds, low-expression uncalled genes skew UP (0.83 / 0.68), a depth confound. Any UP-sign claim read on
     Orion needs the expression-only control.
   - The Broad-16 support should be read as the INCREMENT over each line's uncalled base rate (UP 0.78 vs base 0.73), with
     12/12 lines above 0.5.
   - No block override: in the block the low-I tail adds only +0.02 on H1, and on Orion it is mostly DOWN.

5. **RE-CENTRING CHECK (10-04, `recentre622.py`, `work/recentre622.{out,json}`), prompted by 657's calib657:** the Broad-16 KO
   data carry a global UP offset that no CRISPRi bank shares, so "> 0.5 on 12/12" is not the right null.
   - Against each line's OWN uncalled UP base rate, the bottom-10 % I → UP tail is above base on **8/12 lines**: mean
     **+0.048**, median +0.060, pair-weighted +0.081. The bottom 5 % is above base on 9/12.
   - The honest Broad support is therefore "8/12, about +0.05 over base", not "12/12, 0.781". It is weaker, but it survives.
   - H1 (Arc chemistry, no offset; 0.61–0.71 vs today's 0.52) and 609b's expression-only control stand unchanged.
   - (o)'s hook rows already include base-rate-balanced versions, so idall3's prereg evidence does not change.

## (o) How often is a new line "K562-like" for identity hooks? (2026-10-02 ~14:20Z; `hook622.py`, `work/hook622_{samelab,donor}.json`; asked by [21dd4c] for idall3's prereg)
Instrument: 613d's add rule (`pred_dnmae`, copied verbatim and checked equal to block613d's on random inputs). It is the
predicted raw Δnmae on UNCALLED gate genes (+ = worse), for 12 Broad working lines.
- MDE add from fold480 at the line's own mean UMIs (11k–76k), with a beside-run at 20k.
- Pools 521 f 1 / 0.5 / 0.35, ordered by the authors' fdr.
- Hooks:
  - c4h = CLEAN4 sign at |P − ½| ≥ 0.05;
  - idall3 = sign of z(C4 − ½) + 0.5 z(I) on every CLEAN4-covered gene, ungated, with same-lab I.
- Decomposition: coverage = c4all − c4h (CLEAN4's sign extended to weak priors); identity = idall3 − c4all.

**Two KEYS (they decide the answer).**
- `samelab`, as specified: head and key from the other 11 lines' same-KO response. That is a near-oracle same-lab
  programme, so every hook loses to it.
  - c4h worse than no hook on 12/12 lines, by +0.0045–0.0057.
  - idall3 worse than c4h on 12/12, by +0.0046–0.0058.
  - This is not a D/E/F-like setting.
- `donor`, production-like: head and key from K562 GWPS's CRISPRi response to the same target (another lab and modality,
  as our donor key is for D/E/F).

| donor key, own depth | f 1 | f 0.5 | f 0.35 |
|---|---|---|---|
| (a) c4h − no hook, lines worse / mean | 6/12, −0.0001 | 7/12, +0.0002 | 7/12, +0.0002 |
| **(b) idall3 − c4h, lines worse / mean** | **8/12, +0.0002** (worst T3M5 +0.0030) | 5/12, +0.0001 | **7/12, +0.0005** (worst ONE58 +0.0027) |
| coverage (c4all − c4h) | 10/12 worse, +0.0012 | 10/12, +0.0010 | 9/12, +0.0014 |
| identity (idall3 − c4all) | 3/12 worse, **−0.0010** | 1/12, −0.0010 | 0/12, −0.0009 |
| c4h + identity TAIL extension (weak-prior genes signed only where \|z(I)\| ≥ 1.5) − c4h | **1/12 worse, −0.00017** | 1/12, −0.00015 | 3/12, −0.00011 |

The base-rate-balanced versions (UP and DOWN truth weighted equally per target) give the same picture: idall3 − c4h is
worse on 4–6/12 lines with a mean of ≈ 0, identity is −0.0005 to −0.0007, and coverage is +0.0004 to +0.0007.

**Reading (for idall3's prereg).**
1. **The K562-like stress tail is common, not rare.** With a production-like key, idall3 is worse than the plain CLEAN4
   hook on 5–8 of 12 new cancer lines (mean +0.0000 … +0.0005 raw). One or two lines per 12 lose as much as the K562
   fold (+0.0027–0.0037 vs K562's +0.0024).
2. **Mechanism: idall3's identity part HELPS on 9–12 of 12 lines (−0.001). Its coverage expansion HURTS (+0.001, 8–10/12).**
   Signing the weak-prior genes (|P − ½| < 0.05), which c4h leaves unsigned, costs more than the identity earns there.
   Inside c4h's own coverage the identity rarely changes a sign (≈ 0).
3. A form that keeps c4h's coverage and extends it only to identity TAILS (|z(I)| ≥ 1.5) beats c4h on 9–11 of 12 lines.
   The effect is tiny (−0.0001 to −0.00017 raw). It is a candidate shape, not a lever.
4. Caveats:
   - KO gates lean UP (the balanced rows address this);
   - 613d's rule is the add-only instrument, with its board retrodiction pending (496a, 10-05);
   - Broad's KO targets differ from the panel;
   - K562 GWPS covers most but not all KO targets; uncovered targets get no head and key 0.

## (p) No head override: on HEAD genes the donor key beats the identity rider everywhere (2026-10-02 ~14:15Z; analysis inline, H1 625d `none` s0)
Head = each target's 400 genes by |key| (lib521 state, the leakage-free H1 design key). Sign accuracy on gate genes:

| | n | key | rider z(C4)+0.5z(I) | CLEAN4 | key right where rider ≠ key |
|---|---:|---:|---:|---:|---:|
| head, f 1 | 2,239 | **0.746** | 0.576 | 0.575 | 0.683 (n 1,039) |
| head, low-\|key\| half, f 1 | 923 | **0.679** | 0.560 | 0.574 | 0.621 |
| head, f 0.35 | 1,102 | **0.838** | 0.592 | 0.569 | 0.776 |

Identity's source gain in the head (639d) does not translate into a lever, because production's head signs come from
the key, which is right 0.68–0.89. Closed: no identity override of head signs, not even on weak-key head genes.

# 528 — The DE call set's composition: a call-propensity block instead of the top-CPM block

2026-09-29, research worker (Claude Opus 5.5), mandate `MANDATE.md`. Analysis, the real production emitter on H1, and
patches that are NOT applied. No build, no `vcc`, no `cell-eval2` CLI, no git; nothing written outside this directory.
Peak RSS 7.0–7.54 GB per real-emitter run, ≤ 3.9 GB elsewhere (§8 for the one rule deviation).

## Answer

**POSITIVE on H1, measured on the real production emitter, IF the swapped-in genes do not take the block's CLEAN4
steering.** Deliverability is not the obstacle; which genes are steered is.

1. **The block key.** key = log10 control CPM + 22 · excess dropout. Excess dropout is 517's Poisson-expected minus
   observed detection rate in the recipient's OWN controls. Beta 22 was chosen by the capture of h1-trainval's 126 other
   perturbations, so no h1-test truth enters it. The key swaps 616 of the 2,200 block genes per target (28 %).
2. **Best form, `lcap-kept`.** The call set is head ∪ the propensity block, and only the retained CPM genes keep the CLEAN4
   steering. On H1 over 3 emission seeds, vs 435 (446bj2z's DE construction):

   | regime | overall | seeds | Δjac | Δfid | Δreach | Δnmae |
   |---|---|---|---:|---:|---:|---:|
   | f 1 | **+0.0051** [+0.0033, +0.0070] | +0.0052 / +0.0049 / +0.0053 | +0.0091 | −0.0013 | +0.0096 | +0.0004 |
   | f 0.5 | +0.0048 | | | | | |
   | f 0.35 | +0.0045 | | | | | |
   | cmp | +0.0053 [+0.0033, +0.0073] | | | | | |

   - The member deltas are raw values at f 1.
   - pds and mse are bit-identical on every run.
   - The gain is broad: 58 of 85 targets gain, and the top 3 carry 18 %.
3. **Deliverability.**
   - 98.1 % of the swapped-in genes are called: ≥ 0.97 in every count stratum, including < 1.5 counts per cell.
   - Today's block is called at 0.984.
   - Calls that fail to materialise cost 0.4 % of the prescribed set's jac (0.08380 realised vs 0.08415 prescribed).
4. **The steering is what matters.** Two other forms were measured:
   - Coupled (`lcap`: the swapped-in genes also take CLEAN4 plus the kf shrink): **+0.0027**. The shrink cuts these
     responsive genes' delivered |lfc| from 0.121 to 0.085 (truth median 0.314), so nmae rises +0.0041 raw.
   - Decoupled (`lcap-dec`: steer the CPM block, call the propensity block): nmae is pinned, but 29 % of the
     steered-but-uncalled genes leak through BH (181 calls per target). It prices +0.0063 at f 1 but only +0.0033 at
     f ≤ 0.5.
5. **Mechanism, falsifier and dose.**
   - At matched expression, dropout-prone genes are cell-state programme genes: G2/M, immediate-early/ISR, interferon and
     hypoxia genes. They move ×1.4 more (mean |lfc|) and are called 1.3–3× as often.
   - An expression-matched shuffle of the dropout feature costs **−0.0058** overall (jac −0.0052 raw, exactly as
     emulated).
   - Dose, kept form, swapping 10 / 25 / 100 % of the block: +0.0025 / +0.0041 / +0.0051.
6. **Interactions.**
   - **On top of 527's hook** (446bj2zmz's construction; KA-hook reproduces 517's run bit for bit) the kept block is worth
     more: **+0.0085** [+0.0064, +0.0108] at f 1 over 3 seeds (+0.0083 / +0.0079 / +0.0093), +0.0074 / +0.0075 /
     +0.0082 at f 0.5 / 0.35 / cmp; fid +0.0032, jac +0.0089 (board 0.72 × fid+jac ≈ +0.0042). The hook's programme
     signs reach the called swapped-in genes, so fid turns positive.
   - The coupled form collapses on the hook (+0.0012 at seed 0, nmae +0.009): it pulls the hook's best genes out of the
     hook's complement.
   - **With 523's K 1,600 it is not additive.** At CallSet(400, 1200) the kept swap makes ~377 unprescribed calls per target
     (base 66; low-to-mid-CPM complement genes, not the swapped ones): −0.0012 / −0.0027 / −0.0027 at f 0.5 / 0.35 / cmp
     (2 seeds), +0.0054 at f 1.
7. **Transfer.**
   - On public cancer truths, excess dropout beats its shuffle on 3 of 3 folds.
   - The block beats the CPM block on 2 of 3: K562 +0.071 and HEK293T +0.039 capture; HCT116 −0.012 (an Orion blob of
     high-CPM, low-dropout genes).
   - A/B/C's dropout scale is 1.7–2.3× H1's, so the key is median-rescaled. That moved every fold the right way and gives
     A/B/C swaps of 546 / 572 / 684.
   - The swapped genes are context-specific (Jaccard 0.06–0.14 between contexts).
8. **Board prediction** for a DE-only arm vs 446bj2z (fid and jac at 0.72×, the 523 convention):
   - **Δ(fid + jac) +0.0020 … +0.0025 overall** across the four regimes.
   - With reach at ×1 / ×0 / ×−0.76: +0.0037 … +0.0044 / +0.0020 … +0.0028 / +0.0007 … +0.0016.

**Recommendation.**
- **Build 446bj2zcp** = 446bj2z + `VCC_RC_CALLSET_BLOCK=prop:<abc artifact>:<sha>` (steer `kept`, swap 1) as a DE-only
  information slot. The draft bins are in `PREREGISTRATION-draft.md`.
- If 446bj2zmz (the hook) becomes the incumbent first, build the same knob on top of it instead: its H1 price is larger.
- Keep K 2,600 for this arm. If 523 adopts K 1,600, this knob needs its own re-design before it rides.
- The patches, tests, known-answer design, A/B/C artifacts and October pass-through are in §5.

---

## 1. Known answers (all pass)
| check | result |
|---|---|
| KA-a: 517's callprop517 (2-fold, block 2,200 among non-panel genes by pooled CPM vs 517's propensity) | capture 0.3385689510035993 / 0.3671678299471995, per-target 0.3896447928 / 0.4462257692: **exactly 517's** (`callprop528.py`) |
| **real emitter: `cpm` mode seed 0 (528's patched emitter, block key unset) == 435's production run** | **KA PASS, pre-zs and zs stats (sums, jk, p, q, lfc) bit-equal** (12:45Z) |
| KA-hook: `cpm+mzg` seed 0 (528's harness + 517's Arc-sourced hook as 527 builds it) == 517's real `c4-mzg0.05` run | **PASS, pre-zs and zs bit-equal** |
| V3CHECK: the FINAL emitter (patch_emitter528.py's output, the one the builder patch pins) re-runs `lcap` (steer key; ran on v1) and `lcap-dec` (steer cpm; ran on v2) at seed 0 | **PASS, both, pre-zs and zs bit-equal** |
| price528: 435's four DE members recomputed from its stats with 521's functions (realized order) | equal to 435's metrics json on seeds 0/1/2 (≤ 1e-12) |
| real-emitter K 1,600 base (`cpm-n1200` vs 435, seeds 0 and 1) vs 521 §5's analytic global-K price | −0.0130 [−0.0217, −0.0050] / −0.0008 / +0.0025 / +0.0051 [+0.0027, +0.0075] at f 1 / 0.5 / 0.35 / cmp (seed spread ≤ 0.0009) vs 521's −0.0155 / −0.0003 / +0.0048 / +0.0061: the first real-emitter read of 523's K 1,600 on H1 agrees with 521, a little less positive at f 0.35 |
| feat528 on h1-test vs 517's `excess_dropout` | x / det / exp_det / pooled CPM max \|diff\| **0.0**; the 5-CPM gate == 323's `keep`; mean CPM == `emit_ref_lin` to 1.7e-15 rel., ranking identical |
| truth528 (h1-trainval's scorer calls vs h1-test's controls): the same code on 6 h1-test targets | q max \|diff\| **0.0**, lfc ≤ 3.2e-15, calls identical; ref_lin bit-equal |
| lib521's production block reconstruction vs `block_by(ref_lin)` | identical (asserted) |
| emitter patch unit tests (`test_emitter528.py`, `test_emitter528v2.py`) | 23/23 and 26/26: unset / key = ref_lin / keep = n_cpm bit-identical to production; the three steer forms; refusals |
| builder patch unit tests (`test_patch528.py` on the committed post-527 builder 39c1dc66) | **22/22 PASS** (after `k_emit528.py` was regenerated to the final form) |
| def_day patch tests (`test_defday528.py` on def_day 0b53ba89, 526-patched) | **29/29 PASS** |
| callblock528 `--ka` / inert528 `--selfcheck` | the H1 key gives the priced `lcap` blocks on 85/85 targets; rc_emit's control-CPM arithmetic reproduces h1-test's `emit_ref_lin` bit for bit |
| a pseudo-arm (435 seed 1 filed as seed 0) through price528 | pure re-roll noise: fid+jac overall +0.0001, reach ±0.005 raw (the reach noise scale) |

## 2. The capture screen on the EXACT production call set (H1; `callprop528.py`, `keys528.py`)
Production construction (lib521's reconstruction of 435: head = top 400 by |key| per target, CPM block = top 2,200 by the
emitter's `ref_lin` among the rest). Every key uses only the recipient's OWN controls as features (`lc` = log10 control CPM,
`x` = excess dropout) and never the held-out target's truth: `tv` keys are fitted on h1-trainval's 126 perturbations
(`truth528.py`: their scorer calls against h1-test's 18,400 controls), `2f` keys are 517's 2-fold on h1-test.
`jac_emul` = the jac member if the prescribed set were delivered exactly.

| block key | capture of gate genes | per-target | Δ jac_emul raw [95 % target CI] | swapped / target |
|---|---:|---:|---|---:|
| control CPM (shipped) | 0.3327 | 0.3641 | 0 (jac_emul 0.07476; delivered 0.07470) | 0 |
| `ltv`: lc + 11.54 x (linear logit fitted on trainval) | 0.3599 | 0.4194 | +0.0068 [+0.0051, +0.0086] | 327 |
| `ptv`: 517's 4-term logistic fitted on trainval | 0.3588 | 0.4173 | +0.0065 [+0.0048, +0.0083] | 310 |
| `p2f`: 517's logistic, 2-fold on h1-test | 0.3602 | 0.4202 | +0.0069 [+0.0052, +0.0087] | 322 |
| **`lcap`: lc + 22 x (beta maximising trainval capture)** | **0.3703** | **0.4369** | **+0.0094 [+0.0069, +0.0121]** | 616 |
| `lcap-k1980` (bottom 10 % of the CPM block swapped) | 0.3523 | 0.4095 | +0.0049 [+0.0037, +0.0062] | 205 |
| `lcap-k1650` (25 %) | 0.3656 | 0.4260 | +0.0082 [+0.0061, +0.0104] | 431 |
| `lcap-k1100` (50 %; ≈ the full swap) | 0.3705 | 0.4371 | +0.0094 | 594 |
| **falsifier `shuf`: lc + 22 x_shuffled within 20 expression bins** | 0.3108 | 0.3405 | **−0.0052 [−0.0067, −0.0038]** | 619 |
| matched random swap from the CPM band just below the block | 0.3300 | 0.3694 | −0.0007 [−0.0011, −0.0003] | 213 |

- The coefficients are stable: the logit's beta is 11.5 (trainval), 12.7 / 10.8 (the two h1-test halves), 12.0 (all 211).
  The capture-optimal beta is ~2× the logit's on BOTH truths (trainval peak 0.3602 at 22, h1-test 0.3703 at 20–25).
- The dose curve saturates by the 50 % partial: every key's top 1,100 CPM genes stay in the block.
- Regime-robust (full-delivery jac_emul for `lcap`: +0.0094 / +0.0092 / +0.0087 / +0.0096 at f 1 / 0.5 / 0.35 / cmp); the
  gain lives in strong targets (n_real ≥ 1,000: 65 %; 100–1,000: 30 %; < 100: 5 %).
- Ceiling: a block from each gene's call rate in the SAME line (h1-trainval's 126 perturbations) captures 0.4238 (jac_emul
  +0.023). Not buildable for A/B/C, and H1's call rate does not transfer across lines: as a block key on the public folds it
  loses 0.28–0.33 capture on every fold (192's finding for Orion, now on three folds; mixing it into the CPM or dropout key
  hurts 2/3 folds). Closed.
- Mechanism (H1, strong targets): at matched expression the most dropout-prone tercile has ~1.4× the mean |lfc|
  (quintile 3: 0.139 vs 0.097) and 1.3–3× the call rate: these genes RESPOND more; it is not only test power. The
  swapped-in genes are cell-state programme genes: G2/M (UBE2C, AURKA/B, CENPA/E, KIF20A, HMMR, CCNA2, NDC80, CDCA*),
  immediate-early / ISR (EGR1, FOS, HES1, ATF5, CHAC1, SLC7A11, GDF15), interferon (ISG15, RSAD2, OAS1), hypoxia (CA9,
  PFKFB4, SLC16A3, PDK1), histones; in A (Jurkat) also T-cell subset genes (GZMA, CD4, CD40LG, PTCRA).

## 3. Transfer checks (`abc528.py`, `fold528.py`, `fold528b.py`)
**A/B/C (no truth; from their own controls).** x agrees across contexts (r A–B 0.58, A–C 0.51, B–C 0.69; with H1 0.34–0.41;
517's numbers reproduced); the key lc + 22 x agrees at r 0.44–0.68. A/B/C's dropout scale is 1.7–2.3× H1's (median x 0.020
/ 0.016 / 0.015 vs 0.0088), so the raw beta-22 key would swap 1,097 / 789 / 977 genes (H1: 626). With x rescaled by
median_H1 / median_context the swap is 546 / 572 / 684 (callblock528's artifact). The swapped-in sets are context-specific:
Jaccard 0.06–0.14 between contexts (CPM blocks: 0.32–0.43).

**Public cancer-line truths (480's folds; other chemistries, 150 cells, ~12k UMI; Orion truths carry a technical blob:
relative direction only).** Capture of the fold's gate genes by a 2,200-gene block (Δ vs its CPM block):

| fold (median x; gate/target) | raw lc + 22 x | norm (x · median_H1/median_fold) | falsifier (x shuffled) | best beta, norm form |
|---|---:|---:|---:|---|
| K562-GWPS, cross-study, no blob (0.0050; 3) | +0.058 | **+0.071** | −0.003 | ≥ 40 (+0.081) |
| Orion HEK293T (0.0124; 141) | +0.034 | **+0.039** | −0.111 | 22 |
| Orion HCT116 (0.0156; 68) | −0.033 | **−0.012** | −0.059 | 5 (+0.009) |

x carries call information on every fold (every beta-22 block beats its expression-matched shuffle by +0.03 … +0.15); the
block beats the CPM block on K562 and HEK293T, not on HCT116, whose truth is dominated by high-expression, LOW-dropout genes
(quintile-5 call rate 0.10 / 0.03 / 0.02 by dropout tercile: the blob). The median rescale moves all three folds the right
way (+0.013 / +0.004 / +0.021 vs raw); on H1 it is 1 by construction, so the H1 price is the price of the deployed key.

## 4. Real production-emitter H1 runs (`run528.py`, `drive528*.sh`, `price528.py`, `summ528.py`)
Harness: 435's (run496, RUN426_EMITTER=prod: CallSet(400, 2200, 1.9), CLEAN4 q80 prior, coarse MDE, + the zs pass = 446bj2z's
DE construction on H1) with 528's patched emitter and the block key per target. Each arm vs 435's run of the same seed (or
528's own `cpm+mzg` / `cpm-n1200` run for the interaction arms). Members with 521's exact functions (realized order) on the
truth pools at f 1 / 0.5 / 0.35 (PowerSurface) and cmp. Units: overall (192's INV / 6); board columns = global521.board_cols
(fid and jac at 0.72×, reach ×1 | ×0 | ×−0.76, nmae ×1). CIs: 4,000 target resamples of the seed-averaged rows. All numbers
are the zs form (the zs pass moves nmae by ≤ 0.00006 and nothing else).

### 4.1 Six members (pds / mse are bit-identical on every run: max |Δ| 0.0, per-gene totals equal)
| arm (vs its base) | seeds | overall f 1.0 [95 % CI] | f 0.5 | f 0.35 | cmp [95 % CI] | Δfid / Δjac / Δreach / Δnmae raw (f 1) | fid+jac overall f 1 / cmp | board f 1 [reach ×1 \| ×0 \| ×−0.76] |
|---|---|---|---|---|---|---|---|---|
| `lcap-kept` (vs base) | 3: +0.0052 / +0.0049 / +0.0053 | **+0.0051** [+0.0033, +0.0070] | +0.0048 | +0.0045 | +0.0053 [+0.0033, +0.0073] | -0.0013 / +0.0091 / +0.0096 / +0.0004 | +0.0034 / +0.0031 | +0.0041 \| +0.0023 \| +0.0010 |
| `lcap` (vs base) | 3: +0.0027 / +0.0024 / +0.0028 | **+0.0027** [-0.0005, +0.0057] | +0.0020 | +0.0024 | +0.0034 [+0.0005, +0.0062] | -0.0016 / +0.0090 / +0.0034 / +0.0041 | +0.0032 / +0.0035 | +0.0018 \| +0.0011 \| +0.0007 |
| `lcap-dec` (vs base) | 3: +0.0064 / +0.0060 / +0.0064 | **+0.0063** [+0.0042, +0.0084] | +0.0033 | +0.0033 | +0.0035 [+0.0017, +0.0055] | +0.0007 / +0.0092 / +0.0090 / +0.0000 | +0.0046 / +0.0019 | +0.0050 \| +0.0033 \| +0.0020 |
| `lcap-k1980` (vs base) | 3: +0.0018 / +0.0002 / +0.0012 | **+0.0011** [-0.0016, +0.0030] | -0.0012 | +0.0011 | +0.0011 [-0.0011, +0.0030] | -0.0004 / +0.0047 / -0.0027 / +0.0011 | +0.0019 / +0.0021 | +0.0005 \| +0.0010 \| +0.0014 |
| `lcap-k1650` (vs base) | 3: +0.0019 / +0.0014 / +0.0008 | **+0.0014** [-0.0020, +0.0041] | +0.0003 | +0.0014 | +0.0026 [-0.0001, +0.0050] | -0.0015 / +0.0079 / -0.0043 / +0.0021 | +0.0027 / +0.0030 | +0.0006 \| +0.0014 \| +0.0020 |
| `lcap-k1100` (vs base) | 1: +0.0023 | **+0.0023** [-0.0013, +0.0054] | +0.0010 | +0.0018 | +0.0031 [+0.0001, +0.0059] | -0.0018 / +0.0090 / +0.0012 / +0.0040 | +0.0031 / +0.0034 | +0.0014 \| +0.0011 \| +0.0010 |
| `lcap-k1980-kept` (vs base) | 1: +0.0025 | **+0.0025** [+0.0012, +0.0039] | +0.0014 | +0.0026 | +0.0023 [+0.0011, +0.0036] | -0.0006 / +0.0049 / +0.0045 / +0.0006 | +0.0019 / +0.0018 | +0.0020 \| +0.0012 \| +0.0005 |
| `lcap-k1650-kept` (vs base) | 1: +0.0041 | **+0.0041** [+0.0021, +0.0060] | +0.0024 | +0.0027 | +0.0045 [+0.0028, +0.0063] | -0.0012 / +0.0080 / +0.0057 / -0.0001 | +0.0030 / +0.0028 | +0.0032 \| +0.0022 \| +0.0013 |
| `ltv` (vs base) | 3: +0.0028 / +0.0014 / +0.0015 | **+0.0019** [-0.0013, +0.0045] | +0.0006 | +0.0016 | +0.0025 [-0.0000, +0.0048] | -0.0008 / +0.0065 / -0.0015 / +0.0012 | +0.0025 / +0.0027 | +0.0012 \| +0.0015 \| +0.0017 |
| `p2f` (vs base) | 1: +0.0027 | **+0.0027** [-0.0008, +0.0055] | +0.0010 | +0.0011 | +0.0028 [+0.0002, +0.0051] | -0.0005 / +0.0065 / +0.0018 / +0.0012 | +0.0027 / +0.0028 | +0.0019 \| +0.0016 \| +0.0013 |
| `shuf` (vs base) | 3: -0.0059 / -0.0063 / -0.0054 | **-0.0058** [-0.0119, -0.0017] | -0.0052 | -0.0029 | -0.0021 [-0.0043, -0.0001] | -0.0001 / -0.0052 / -0.0177 / +0.0003 | -0.0024 / -0.0017 | -0.0052 \| -0.0018 \| +0.0007 |
| `lcap+mzg` (vs cpm+mzg) | 1: +0.0012 | **+0.0012** [-0.0023, +0.0043] | -0.0002 | +0.0002 | +0.0016 [-0.0015, +0.0045] | -0.0018 / +0.0091 / +0.0036 / +0.0093 | +0.0031 / +0.0034 | +0.0003 \| -0.0003 \| -0.0009 |
| `lcap-kept+mzg` (vs cpm+mzg) | 3: +0.0083 / +0.0079 / +0.0093 | **+0.0085** [+0.0064, +0.0108] | +0.0074 | +0.0075 | +0.0082 [+0.0060, +0.0104] | +0.0032 / +0.0089 / +0.0120 / -0.0015 | +0.0058 / +0.0058 | +0.0069 \| +0.0046 \| +0.0029 |
| `lcap-n1200` (vs cpm-n1200) | 1: -0.0012 | **-0.0012** [-0.0059, +0.0030] | +0.0012 | +0.0033 | +0.0018 [-0.0031, +0.0064] | -0.0062 / +0.0085 / -0.0047 / +0.0026 | +0.0004 / +0.0039 | -0.0013 \| -0.0004 \| +0.0002 |
| `lcap-kept-n1200` (vs cpm-n1200) | 2: +0.0052 / +0.0056 | **+0.0054** [-0.0001, +0.0119] | -0.0012 | -0.0027 | -0.0027 [-0.0054, +0.0001] | +0.0048 / +0.0076 / -0.0024 / +0.0011 | +0.0061 / -0.0019 | +0.0037 \| +0.0041 \| +0.0045 |

### 4.2 Deliverability (seed 0; `price528.deliverability`)
| arm | pds / mse Δ (all seeds) | per-gene totals == base | swapped-in / target | swapped-in call rate | missed / target (base) | unsuppressible / target (base) | swapped-out now called | realized vs prescribed jac (f 1) |
|---|---|---|---|---|---|---|---|---|
| `lcap-kept` | max 0.0e+00 | True | 616 | 0.9807 | 34.0 (35.3) | 20.4 (7.0) | 0.0019 | 0.08380 / 0.08415 |
| `lcap` | max 0.0e+00 | True | 616 | 0.9706 | 37.4 (35.3) | 4.8 (7.0) | 0.0004 | 0.08370 / 0.08415 |
| `lcap-dec` | max 0.0e+00 | True | 616 | 0.9898 | 26.0 (35.3) | 188.6 (7.0) | 0.2937 | 0.08388 / 0.08415 |
| `lcap-k1980` | max 0.0e+00 | True | 205 | 0.9730 | 33.6 (35.3) | 6.1 (7.0) | 0.0007 | 0.07937 / 0.07968 |
| `lcap-k1650` | max 0.0e+00 | True | 431 | 0.9692 | 35.6 (35.3) | 5.5 (7.0) | 0.0007 | 0.08257 / 0.08302 |
| `lcap-k1100` | max 0.0e+00 | True | 594 | 0.9701 | 36.7 (35.3) | 4.8 (7.0) | 0.0004 | 0.08372 / 0.08417 |
| `lcap-k1980-kept` | max 0.0e+00 | True | 205 | 0.9954 | 30.9 (35.3) | 9.9 (7.0) | 0.0015 | 0.07955 / 0.07968 |
| `lcap-k1650-kept` | max 0.0e+00 | True | 431 | 0.9857 | 30.5 (35.3) | 15.5 (7.0) | 0.0021 | 0.08269 / 0.08302 |
| `ltv` | max 0.0e+00 | True | 327 | 0.9733 | 34.0 (35.3) | 5.9 (7.0) | 0.0008 | 0.08115 / 0.08156 |
| `p2f` | max 0.0e+00 | True | 322 | 0.9730 | 33.9 (35.3) | 5.8 (7.0) | 0.0008 | 0.08120 / 0.08164 |
| `shuf` | max 0.0e+00 | True | 619 | 0.9608 | 47.6 (35.3) | 7.1 (7.0) | 0.0004 | 0.06949 / 0.06957 |
| `lcap+mzg` | max 0.0e+00 | True | 616 | 0.9687 | 40.1 (38.1) | 6.3 (11.3) | 0.0006 | 0.08366 / 0.08415 |
| `lcap-kept+mzg` | max 0.0e+00 | True | 616 | 0.9593 | 48.6 (38.1) | 30.1 (11.3) | 0.0021 | 0.08349 / 0.08415 |
| `lcap-n1200` | max 0.0e+00 | True | 504 | 0.9999 | 0.0 (0.0) | 36.0 (66.3) | 0.0000 | 0.07664 / 0.07636 |
| `lcap-kept-n1200` | max 0.0e+00 | True | 504 | 1.0000 | 0.0 (0.0) | 377.4 (66.3) | 0.0003 | 0.07587 / 0.07636 |

Swapped-in genes by control count per cell (lcap-kept): 0.985 at < 1.5, 0.973 at 1.5–2, 0.985 at 2–3 (today's block: 0.975 at
2–3, 0.981 at 3–5, 0.991 at 5–10, 0.998 above). The emitter can make these lower-count genes significant; 441 PS3's loss
(the prior-confidence block) came from its genes rarely being true DE genes, not from calls failing to materialise.

### 4.3 Why the steering form decides the price (seed 0, `lcap` vs base)
| gene set (per target) | form | sign accuracy all / gate | median delivered \|lfc\| on gate genes (truth 0.314 / 0.167) | nmae contribution Δ |
|---|---|---|---|---|
| swapped-in (616; 75.6 gate) | base: uncalled, design sign + MDE | 0.508 / 0.533 | 0.121 | — |
| | coupled: called, CLEAN4 + kf shrink | 0.509 / 0.541 | 0.085 | **+0.0035** |
| | kept: called, design sign + MDE | 0.508 / 0.533 | 0.122 | +0.0000 |
| swapped-out (616; 46.9 gate) | base: called, CLEAN4 | 0.516 / 0.583 | 0.078 | — |
| | coupled / kept: uncalled, design + MDE | 0.508 / 0.552 | 0.096 | +0.0004 |
| | decoupled: uncalled but still CLEAN4-steered | — | — | 0 (but 29 % leak through BH: 181 calls / target, uniform across target strength; same count and \|lfc\| as the held ones) |

fid's small loss in every form is the swapped-in genes' lower sign predictability (0.508–0.509 vs the swapped-out
CLEAN4-signed 0.516); jac's gain is the same in all three forms (+0.0090 … +0.0092).

### 4.4 Dose and falsifier
Kept form, swapping 10 / 25 / 100 % of the block (205 / 431 / 616 genes): +0.0025 / +0.0041 / +0.0051 (f 1; fid+jac
+0.0019 / +0.0030 / +0.0034). Coupled: +0.0011 / +0.0014 / +0.0023 (50 %) / +0.0027. The logit-fitted key (`ltv`, beta
11.5, coupled) +0.0019 and 517's own 2-fold form (`p2f`, coupled) +0.0027 sit on the same curve. The shuffle falsifier
loses −0.0058 overall with jac −0.0052 raw on all three seeds, exactly the emulated value.

### 4.5 Interactions (each vs 528's own base run with the same construction)
**On top of 527's hook** (`+mzg` = run517's `c4-mzg0.05`: the MDE add on complement genes takes Arc's meanz_ne programme sign at tau 0.05 on 467's strong targets and CLEAN4's on the weak third = 446bj2zmz's DE construction; KA-hook PASS), each vs `cpm+mzg` of the same seed:

| arm | regime | overall [95 % CI] (seeds); members raw |
|---|---|---|
| `lcap-kept+mzg` | 1.0 | +0.0085 [+0.0064, +0.0108] (seeds +0.0083 / +0.0079 / +0.0093); fid +0.0032, jac +0.0089, reach +0.0120, nmae -0.0015; fid+jac overall +0.0058 |
| `lcap-kept+mzg` | 0.5 | +0.0074 [+0.0042, +0.0105] (seeds +0.0073 / +0.0080 / +0.0071); fid +0.0031, jac +0.0088, reach +0.0068, nmae -0.0016; fid+jac overall +0.0057 |
| `lcap-kept+mzg` | 0.35 | +0.0075 [+0.0050, +0.0104] (seeds +0.0068 / +0.0078 / +0.0078); fid +0.0031, jac +0.0083, reach +0.0088, nmae -0.0014; fid+jac overall +0.0054 |
| `lcap-kept+mzg` | cmp | +0.0082 [+0.0060, +0.0104] (seeds +0.0081 / +0.0084 / +0.0081); fid +0.0031, jac +0.0091, reach +0.0095, nmae -0.0022; fid+jac overall +0.0058 |
| `lcap+mzg (coupled)` | 1.0 | +0.0012 [-0.0023, +0.0043] (seeds +0.0012); fid -0.0018, jac +0.0091, reach +0.0036, nmae +0.0093; fid+jac overall +0.0031 |
| `lcap+mzg (coupled)` | 0.5 | -0.0002 [-0.0058, +0.0045] (seeds -0.0002); fid -0.0017, jac +0.0091, reach -0.0073, nmae +0.0072; fid+jac overall +0.0031 |
| `lcap+mzg (coupled)` | 0.35 | +0.0002 [-0.0043, +0.0039] (seeds +0.0002); fid -0.0016, jac +0.0086, reach -0.0043, nmae +0.0073; fid+jac overall +0.0030 |
| `lcap+mzg (coupled)` | cmp | +0.0016 [-0.0015, +0.0045] (seeds +0.0016); fid -0.0016, jac +0.0095, reach +0.0047, nmae +0.0099; fid+jac overall +0.0034 |

With the hook the called swapped-in genes carry the programme's signs where it is confident: sign accuracy 0.5277 (all) / 0.5946 (gate) against 0.508 / 0.533 with design signs, so fid turns positive; the swapped-in call rate is 0.9593 and unsuppressible calls 30.1 / target (base 11.3). The coupled form moves the swapped-in genes into the steered block, i.e. OUT of the hook's complement, which is where the hook earns its nmae (517 §2.1): nmae +0.009.

**At 523's K 1,600** (`-n1200` = CallSet(400, 1200, 1.9)), each vs `cpm-n1200` of the same seed:

| arm | regime | overall [95 % CI] (seeds); members raw |
|---|---|---|
| `lcap-kept-n1200` | 1.0 | +0.0054 [-0.0001, +0.0119] (seeds +0.0052 / +0.0056); fid +0.0048, jac +0.0076, reach -0.0024, nmae +0.0011; fid+jac overall +0.0061 |
| `lcap-kept-n1200` | 0.5 | -0.0012 [-0.0043, +0.0020] (seeds -0.0015 / -0.0008); fid -0.0067, jac +0.0050, reach +0.0021, nmae +0.0004; fid+jac overall -0.0015 |
| `lcap-kept-n1200` | 0.35 | -0.0027 [-0.0056, +0.0001] (seeds -0.0027 / -0.0028); fid -0.0077, jac +0.0045, reach -0.0019, nmae +0.0006; fid+jac overall -0.0022 |
| `lcap-kept-n1200` | cmp | -0.0027 [-0.0054, +0.0001] (seeds -0.0029 / -0.0024); fid -0.0078, jac +0.0052, reach -0.0046, nmae -0.0006; fid+jac overall -0.0019 |
| `lcap-n1200 (coupled)` | 1.0 | -0.0012 [-0.0059, +0.0030] (seeds -0.0012); fid -0.0062, jac +0.0085, reach -0.0047, nmae +0.0026; fid+jac overall +0.0004 |
| `lcap-n1200 (coupled)` | 0.5 | +0.0012 [-0.0045, +0.0065] (seeds +0.0012); fid -0.0053, jac +0.0101, reach -0.0033, nmae -0.0007; fid+jac overall +0.0016 |
| `lcap-n1200 (coupled)` | 0.35 | +0.0033 [-0.0013, +0.0082] (seeds +0.0033); fid -0.0036, jac +0.0106, reach +0.0040, nmae +0.0011; fid+jac overall +0.0028 |
| `lcap-n1200 (coupled)` | cmp | +0.0018 [-0.0031, +0.0064] (seeds +0.0018); fid -0.0030, jac +0.0122, reach -0.0067, nmae +0.0031; fid+jac overall +0.0039 |

At K 1,600 the key swaps 504 of the 1,200 block genes (top-CPM genes leave) and every prescribed gene is called, but the kept form makes 377 unsuppressible calls / target (base 66): low-to-mid CPM complement genes (median CPM rank ~4,300), spread over targets, 83 % of them not leaking in the base, and NOT the swapped-out genes (call rate 0.0003). The coupled form leaks less (36) but loses fid. The size and the composition are not additive under today's emitter.

## 5. Deliverables (NOT applied)
| file | what | tests |
|---|---|---|
| `patch_emitter528.py` → `k_emit528.py` (final form, byte-identical to `k_emit528v3.py`) | `rc_hybrid_block_k(..., block_key=None, block_keep=0, block_steer="key")`: the CallSet's CPM block by a key (top `block_keep` CPM genes kept, the rest by the key; stable ties like the CPM key); `block_steer` = key (the steered block is the called block) / cpm (steer the CPM block) / kept (steer only the called block's CPM part). Unset = the verbatim code path. Source pinned: production 1658ae2c → final 0f56209d | `test_emitter528v2.py` 26/26 (unset, key = ref_lin and keep = n_cpm all bit-identical to production; the three steer forms; v3 == v1 / v2 on their modes); `test_emitter528.py` 23/23; V3CHECK on real H1 runs (§1) |
| `patch_builder528.py` | `VCC_RC_CALLSET_BLOCK = cpm (default) \| prop:<npz>:<sha256>[:<swap>[:<steer>]]` (steer default `kept`); per context `key_<c>[rc_control.keep]` → the emitter; keywords reach the emitter only when set; emitter pin → the final emitter; config `rc_callset_block` (source, path, sha256, swap, block_keep, steer, per-context audit) | `test_patch528.py` on the committed post-527 builder 39c1dc66: parse (steer forms, 523's size knob composes: swap 0.25 at CPM 1600 → keep 1200), 8 refusals, exactly 6 hunks, compiles, emitter pin, k_emit528 == the patch of production |
| `callblock528.py` | the per-context artifact from each context's OWN controls: key = log10 control CPM + 22 · x · (median_H1 / median_context); beta 22 and median_H1 = 0.00877448551593013 pinned (no H1 file on the day) | `--ka`: the H1 key gives the priced `lcap` blocks on 85/85 targets; deterministic (rerun: same sha) |
| `work/callblock528_abc.npz` | A/B/C artifact, sha256 `04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd`; swaps (no head) A 546, B 572, C 684 | — |
| `inert528.py` → `work/callblock528_inert_abc.npz` | KA-2's inert artifact: key = the context's control CPM in `rc_emit.control_block_from_counts`' own arithmetic, sha256 `e62ee76539f6932728f9c58bd37d72a8e4c8cd74d4657ae7251cecfe62d23b69` | `--selfcheck`: reproduces h1-test's `emit_ref_lin` bit for bit |
| `patch_defday528.py` | October: `DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]]` runs callblock528 on the day's bundle (`$B`, `$CTXS`), exports `VCC_RC_CALLSET_BLOCK` for lvl0 + re-level builds, checks the final meta carries it (sha, swap, steer), records it in def_day.final.env, allows `rc_callset_block` in the rule-10 diff | `test_defday528.py` on def_day 0b53ba89: 7 insertions only, bash -n, validation (dies on the unpatched builder), the compute block on the real bundle reproduces the A/B/C artifact byte for byte, the final-meta check |
| `PREREGISTRATION-draft.md` | draft board bins for 446bj2zcp vs 446bj2z (and the numbers for the variant on top of 446bj2zmz) | — |

**Known-answer design for the coordinator** (the committed builder + `patch_builder528.py`, the emitter +
`patch_emitter528.py`; builds via 527's `run_ka_527.sh` pattern, output sha256 compared):
- **KA-1:** 446bj2z's env with `VCC_RC_CALLSET_BLOCK` unset → output `ec3a5ae6…` byte for byte (the emitter's unset path is
  verbatim and the builder passes no new keyword).
- **KA-2 (inert, exercises the whole keyed path, steer `kept`):** 446bj2z's env +
  `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_inert_abc.npz:e62ee76539f6932728f9c58bd37d72a8e4c8cd74d4657ae7251cecfe62d23b69`
  → output `ec3a5ae6…` byte for byte (the key's order is the CPM order: the called block, the steered block, every draw and
  every count are today's).
- **Arm 446bj2zcp:** 446bj2z's env +
  `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd`.
  Expected config diff vs 446bj2z: exactly `rc_callset_block` (+ `de_sign_emitter_sha256`, 499 F2's no-output-effect pin).
  pds and mse must equal 446bj2z's to ≤ 1e-6 raw.
- On top of 446bj2zmz: the same knob added to 446bj2zmz's env; KA-1 then rebuilds 446bj2zmz's output.

## 6. Recommendation
- Build **446bj2zcp** (steer `kept`, swap 1) as a DE-only information slot vs 446bj2z (bins: `PREREGISTRATION-draft.md`;
  primary Δ(fid + jac) +0.0020 … +0.0025 overall). pds/mse are pinned, so the read is clean.
- If 446bj2zmz becomes the incumbent first, arm the same knob on top of it: on H1 the block's increment over the hook
  (+0.0085) is larger than over 446bj2z (+0.0051) (§4.5), because the hook supplies better signs for exactly the genes the
  new block calls.
- October: if the board read is POSITIVE, apply `patch_defday528.py` and run `DEF_DAY_CALLSET_BLOCK=prop` (the day's keys
  come from D/E/F's own controls; nothing is fitted on the day).
- Do not stack with 523's K 1,600 as is: at CallSet(400, 1200) the kept swap makes ~377 extra calls per target (low-CPM
  complement genes that the base does not call; not the swapped genes). If 523 reads K 1,600 positive, this knob needs its
  own re-design first.
- The coupled form (`steer key`) and the decoupled form (`steer cpm`) stay available in the code for diagnosis only.

## 7. Files (all in this directory; outputs in `work/`)
| script | what | peak RSS |
|---|---|---|
| `feat528.py` | per-gene control features (pooled / mean CPM, detection, Poisson-expected detection, excess dropout) for H1 and A/B/C | 1.77 GB (H1), 0.89 GB (A/B/C) |
| `callprop528.py` | KA-a (517's callprop, exact) and the first capture screen | 0.3 GB |
| `truth528.py` | h1-trainval's 126 perturbations' scorer calls vs h1-test's controls (+ KA on 6 h1-test targets) | **3.88 GB** (§8) |
| `keys528.py` | the per-target block keys (tv / 2f / lcap / shuf / rand) and their capture / jac emulation | 0.34 GB |
| `fold528.py`, `fold528b.py` | the public-fold transfer check (raw, norm, qm, shuffle; the K562 controls streamed from the dense source) | 1.64 GB / 0.16 GB |
| `abc528.py` | A/B/C feature stability and swap overlap | 0.05 GB |
| `callblock528.py`, `inert528.py` | the builder artifacts (and KA / self-check) | 1.32 GB / 0.77 GB |
| `patch_emitter528.py`, `k_emit528.py` (final), `k_emit528v2.py`, `k_emit528v3.py`, `test_emitter528.py`, `test_emitter528v2.py` | the emitter patch, its generated copies (v1 ran the coupled / partial / ltv / p2f / shuf / cpm runs, v2 the -dec runs, v3 the -kept runs; V3CHECK and the unit tests show they agree on each other's modes), unit tests | < 0.3 GB |
| `patch_builder528.py`, `test_patch528.py`, `patch_defday528.py`, `test_defday528.py` | builder and October patches and their tests | < 1.4 GB |
| `run528.py`, `drive528*.sh` | the real-emitter H1 harness (modes: cpm, lcap, ltv, p2f, shuf, rand; -k<keep>, -dec, -kept, -n<cpm>, +mzg) and its gated drivers | 7.0–7.54 GB per run |
| `price528.py`, `summ528.py`, `fill528.py` | pricing (six members, regimes, CIs, board columns, deliverability, sign accuracy), the tables, and the RESULTS fill | 0.62 GB |
| `work/price528.json`, `work/summ528.md`, `work/*.log`, `work/*.json` | every number quoted | — |
| `work/h1_<mode>/` | the runs' stats / metrics / diag (≈ 66 MB per seed; regenerable) | — |

## 8. Deviations
- **Memory.** `truth528.py` peaked at **3.88 GB for ~2 min, twice** (11:12–11:19Z) while 527's build (4.1–4.3 GB) and a 517 H1
  job (~7 GB) were alive and free+inactive was 25–27 GiB. The mandate asks for ≥ 33 GiB free+inactive before a step expected
  above 2 GB; I expected < 2 GB (a gene-chunked design) and did not check. It stayed under the 4 GB cap and nothing else was
  affected. Every other step before 527's build exited was ≤ 1.77 GB.
- **Scheduling.** After 527's build exited (12:39Z), 526's def_day rehearsal started building (12:40Z) and would have held a
  "no builder" gate for ~2 h. I replaced the driver (drive528 → c → d → e, each killed between runs; skip-if-done) with the
  mandate's post-527 rule only: one run of mine at a time, ≤ 10 GB, free+inactive ≥ 23 GiB (25 while a builder is alive),
  never while a 517 H1 job is alive. Every run started with free+inactive 27–34 GiB. Two runs (shuf s0, lcap-k1980 s0)
  finished as orphans of a killed driver; their outputs are complete (checked by the next driver's skip logic).
- **Design iteration on H1.** The `kept` steering was added after the coupled and decoupled forms were read on H1 (a
  principled third form, but a choice among three on the same instrument). The board read is the test; the draft bins
  say so.

# 428 — Our A/B/C packages carry 24 genes at delivered lfc ≈ −35; the zero-gene repair fixes them (`VCC_DE_ZERO_REPAIR`)

Status: 2026-09-23 05:25Z, PI.
- Audit on the restored 426a package, which shares its carriers' per-gene totals with 417a/417e (the A-channel does not
  move totals).
- Repair tested on H1 and in production as `VCC_DE_ZERO_REPAIR` (builder `986e8fad…`, module
  `src/vcc2026/zero_repair.py` `a6208eaa…`; emitter unchanged `72be7e16…`).

## The defect
A gated gene (control mean CPM > 5, target gene excluded) that a block's 400 cells never express has a delivered
arithmetic-mean lfc of log2(1e-9 / ctrl) ≈ −35.
- It happens when the carrier design's knockdown is so deep that pass 1 emits zero counts for the gene.
- Ratio-mode ballast copies the carriers' composition (× bulk/DE ratio), so it gives the gene zero as well.
- If the reference calls the gene, it adds ~30 to that target's nmae numerator. The Wilcoxon call and the sign are
  unaffected (the gene is called down either way).

**`audit428.py` on 426a:**
- 24 (context, target, gene) cells in 24 of 900 blocks: A 6, B 11, C 7. Examples: ANKZF1→ATG9A (the 66-bp cis
  neighbour; H1's reference had it at −4.30), PPIP5K2→GIN1 (11 bp), ZDHHC4→SOD1, WBP4→NDUFA13, ZBTB48→MTFP1, UNK→H3F3B.
- 29 cells sit below lfc −5 in total.
- Only 2 of the 24 are cis-neighbours. The rest look like donor-row knockdowns (off-target-like), which may or may not
  be in A/B/C's gates.
- H1's shipped construction has none: its only −34.7 cell is NGRN in its own knockdown, excluded from every metric.

## The repair (`repair428.py` → `src/vcc2026/zero_repair.py`)
- For every such gene, the carriers get n = max(1, round(ctrl_cpm · 2^floor · Σ carrier libs / 1e6)) single counts
  (floor −5 log2, a 97 % knockdown).
- Each count goes into a distinct random carrier cell (seeded per block). Each cell pays with one count from its own
  largest sub-gate gene (control CPM ≤ 5, never tested).
- Every library and every other tested gene is unchanged (asserted), and so is the ballast.

**H1 (`run426.py --zero-repair -5`):**
- **Knock-out test:** 427's κ 3 carrier raise zeroes 61 genes and costs **nmae +0.0816**. With the repair it is
  **+0.0005**, and 60 genes are repaired (1 cell stays at ≤ −8 because its carriers had no sub-gate count to pay with).
- **Known answer:** base + repair is bit-identical to base (no zeros to repair).

**Board prize (estimate, not measured):**
- On H1 a zeroed gate cell cost ≈ 0.0013 nmae on the 58-target average. On A/B/C, with ~250 scored targets per
  context, a cell costs ≈ 30/Σ|l|_target/250 ≈ 0.001 per context.
- The 24 cells then cost up to ~0.008 raw nmae on the board if they sit in the gates (≈ +0.002 overall), and 0 if they
  do not.
- The repair cannot hurt: a non-gate cell is not scored, and a gate cell's error drops from ~30 to ~1–3. The bulk moves
  by ≤ n counts per gene (pds/mse negligible).

**Arm.** 428a = 426a + `VCC_DE_ZERO_REPAIR=-5`, one variable vs 426a (config diff exactly `de_zero_repair`). The build
runs the audit on its own artifact. It is a 09-25 candidate, a pure-nmae read against 426a.

Files: `audit428.py`, `repair428.py`, `patch_zr428.py`, `test_zr428.py`, `run428_arm_build.sh`.

## Does the zero artifact explain 417c's board nmae loss (γ 0 carriers)? No
- On H1, 416's designs carry non-target zeros: γ 0 has 3, ×1.25 carriers 8, 228's gain 17, the base 0.
- On A/B/C, from restored packages, audited and deleted:
  - **417b** (standard carriers): 24 zeros, 31 cells < −5;
  - **417d** (= 417b + γ 0 carriers): **23 zeros**, 26 cells < −5.
- γ 0 did not add zeros on A/B/C, so 417c's board nmae +0.0089 is a real design effect. The rule that H1 cannot price
  nmae for carrier-design changes stands.
- Neither package has any extreme positive lfc: max +2.73 / +2.86, 0 cells > +3.

## 428a built (06:06Z)
Package `c7bf5d47…`; config diff exactly `de_zero_repair` vs 426a; novelty rc 0.
- **Audit of its own artifact: 0 zeroed gated genes (426a had 24).** lfc < −5: 13 cells (repaired genes land near the
  −5 floor). Max +2.30.
- Repair diagnostics: at most 1 gene repaired and 9 counts moved per block.
- CLIFF 0.9476 / 0.9407 / 0.9447 vs 426a's 0.9477 / 0.9408 / 0.9447. The bulk moves in the 4th decimal, because the
  repaired genes' pseudobulk leaves exact zero. Expect |Δpds| and |Δmse| ≲ 1e-4 against 426a.
- **Board prediction vs 426a:** Δnmae 0 … −0.008 (−0.008 only if all 24 cells sit in reference gates). fid, jac and
  reach are unchanged: the repaired genes stay called down.

## Refinement idea (not built): repair the bulk too, on the ballast side
- The carrier-side repair (428a) runs after `ballast_rows`. It fixes the DE mean (lfc ≈ −5), but the zeroed gene's
  bulk coordinate stays at about −b_ctrl, a complete knockdown in the pseudobulk. For SOD1 in ZDHHC4 that is ≈ −3.3,
  whatever the bulk design wanted.
- A **ballast-side design fallback** would fix both. In `ballast_rows`, genes with t_c = 0 and profile_bulk > 0 would
  take the design count profile_bulk × total instead of the ratio-mode zero.
  - The bulk coordinate returns to the design.
  - The DE mean becomes (K/400)·ballast CPM, which gives lfc ≈ log2(0.025 · bulk FC) ≈ −5.3. That is the same floor
    the carrier repair targets.
- pds: 24 spurious coordinates across 900 blocks are worth ≲ +0.0005 overall on a rough ranking estimate. Worth it
  only as a rider on a later bulk arm.
- Order caveat: running the carrier repair BEFORE ratio-mode ballast would let the ratio (profile_bulk / profile_de,
  huge where the DE design is extreme) amplify the repair counts into the bulk. Keep the carrier repair after the
  ballast, as built.

## Where the zeroed genes come from, and do they replicate?
- The 20 unique zeroed (target, gene) pairs carry strong K562-GWPS knockdowns (δ −0.5 … −1.4) and **kolf2 z ≈ 0**
  (−0.32 … +0.16) on every pair it covers.
- Only 2 targets are in H1-2025 (Arc's guides), and both replicate strongly there: AKT2→CNFN −0.80 (20 row-sd) and
  ANKZF1→ATG9A −0.77. Those two are probably in A/B/C's gates.
- **Across all 237 targets shared by H1-2025 and K562-GWPS, strong K562 knockdowns (δ < −0.5) of non-cis genes
  replicate in H1 only 12 % of the time, and in HCT116 10 %.** Cis genes (≤ 5 kb) replicate 100 % / 88 %. So most
  strong non-cis K562 knockdowns are K562-specific (guide off-targets or K562 biology).
- Do they matter for our bulk? K562-only strong knockdowns (Orion HCT116 and HEK293T |δ| < 0.1) number 320 cells over
  259 targets (1.2 per target) and carry 0.17 % of 417a's non-panel bulk row energy on average (max 6.2 % on one
  target). Too small to be an arm; noted, not built. 325–342's full coordinate median was closed.

## 09-25 option builds (06:07–07:33Z), via `run_arm_generic.sh`
Both built on builder `986e8fad`, emitter `72be7e16`.

| package | = | config diff vs parent | zero audit | CLIFF (A/B/C) | sha |
|---|---|---|---|---|---|
| **428h** | 428a + 417h's re-level | exactly `bulk_ctx_scale` (vs 428a) | 0 zeros; 13 cells < −5; max +2.45 | 0.9638 / 1.0337 / 1.0837 | `c15d3680…` |
| **428g** | 417g + MDE + repair | exactly `de_mde`, `de_zero_repair`, emitter sha (vs 417g) | 0 zeros; 15 cells < −5; max +2.43 | 0.9646 / 1.0345 / 1.0850 | `f59ac498…` |

428k (417k's bulk + MDE + repair) is building next (`chain428k.sh`).

## The bulk side of the zeros is harmless: the ballast-side fallback is closed (07:45Z)
- For the 24 zeroed cells, the pseudobulk coordinate is −b_ctrl = −0.42 … −1.08 (log1p-CP50k). These genes are only
  moderately expressed.
- 426a's loaded bulk design wanted −0.10 … −2.14 on the same genes. The design is itself an extreme knockdown, and
  zero expression only clamps it.
- Excess squared error from the zero, with truth ≈ 0, relative to each row's design energy: median −0.01, max 0.88.
  The max is ANKZF1→ATG9A, a real cis knockdown whose truth is not 0.
- The defect is a DE-mean problem only (lfc −35). The carrier-side repair (428a) is the whole fix.
- **Single-donor outliers in general are not a bulk lever.**
  - Over all panel targets, K562-only knockdowns (K562 δ < −0.3, both Orion lines |δ| < 0.1) number 4,277 cells.
  - They keep a median 0.04× their K562 value in our designs (q90 0.22×), and carry 0.33–0.36 % of the non-panel bulk
    row energy on average (max 4 %), in both 417a's bulk and the stack.
  - The composition already shrinks them, so an outlier-suppression arm would have almost nothing to act on.
| **428k** | 417k + MDE + repair | exactly `de_mde`, `de_zero_repair`, emitter sha (vs 417k) | 0 zeros; 17 cells < −5; 1 cell > +3 (max +3.09) | 0.9600 / 1.0278 / 1.0778 | `338aa226…` |

## Note 2026-09-25 01:30Z (after 433h's board read): a shallower repair floor is a candidate rider
433h's repair (floor −5) gained nmae −0.0247 on the board. The old coordinator session flagged, as it stood down, that the
repaired genes now deliver ≈ −5 (97 % knockdown) and are mostly K562-only knockdowns that replicate rarely. If their A/B/C
truth is moderate (−1 … −2), a floor of −2 would recover the residual |−5 − l|. On H1, the κ 3 zeros repaired to −5 left
+0.0005 of their +0.0816 cost (0.6 %). The board residual is likely similar in share, which suggests ~+0.0002 … +0.001
overall. Low priority: price the floor on H1 (floor −1 / −2 / −3 against −5) before any build.

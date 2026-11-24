# 462 — Programmes and lab signatures on the Broad-16 (one lab, 12 cancer lines, matched targets). Analytic, no build.

2026-09-26 10:20Z, coordinator. Scripts: `lab462.py` (pairwise programme agreement by lab and line on a matched target
set), `halves462.py` (cell halves with disjoint controls). Outputs: `work/lab462.{json,md}`, `work/halves462.json`.

## Why
- 458 (09-30 slot) bets that a bank's per-gene programme is mostly its lab's signature (450h–l), so Arc's own H1
  programme would be the right block prior for Arc's cancer contexts.
- The Broad-16 release (Cas9 knockout, one lab and protocol, 16 cancer lines, the same 100 targets; 268 built its units)
  is the cleanest matched design we hold. 71 of its targets are also screened by K562-GWPS, Orion HCT116/HEK293T and
  VIPerturb K562.
- Its targets are 0/300 of the A/B/C panel, so it is not an October same-line bank. It answers a mechanism question.

## 1. Pairwise programme agreement (P_down − 0.5, Pearson r), 71 matched targets, 7,015 genes
| pair type | labs | pairs | mean r | sign agreement |
|---|---|---:|---:|---:|
| different line, same lab | Broad × Broad | 66 | +0.061 | 0.537 |
| different line, same lab | Xaira × Xaira | 1 | +0.053 | 0.556 |
| same line, different lab | Weissman K562 × Satija K562 | 1 | +0.100 | 0.635 |
| different line, different lab | Broad × Weissman / Xaira / Satija | 12 / 24 / 12 | +0.083 / +0.067 / +0.033 | 0.54–0.57 |

On the smaller set that also includes the essential screens (23 targets): Weissman × Weissman other lines r +0.221 (9
pairs, sign agreement 0.645); Broad × Broad +0.102; Broad × Weissman +0.129.

## 2. The Broad lines' programme reliability is mostly their control's noise
Each (line, target) response is measured against the line's own 10 olfactory-receptor-KO controls (about 500 cells, vs
43–129 per target). Every target in a line shares that control, so control noise enters as a common per-gene offset: it
IS a "programme".

| within-line reliability of P_down (82 targets, 12,199 genes, 12 lines) | mean r |
|---|---:|
| target halves, both against the SAME full control (what a split of targets reads) | **0.801** |
| cell halves, each against its OWN half of the controls (disjoint on both sides) | **0.121** (Spearman–Brown 0.21) |
| cross-line, line i half A vs line j half B (no shared noise) | 0.043 |
| cross-line, full programmes | 0.057 |

- About 85 % of a Broad line's apparent programme reliability is its control's noise.
- Honestly measured, lines share about a third of their programme within one lab: cross-line 0.043 vs within-line 0.121.
- Any programme from a small-control screen is dominated by that control. Target split-halves cannot detect it. Cell
  (or control) halves can.

## Reading
- **For 458.** One lab and protocol gives a moderate common component across its cancer lines (~1/3 of the reliable
  programme here, with Cas9 KO's noisier responses), not a dominant one. The Weissman lines' higher cross-line agreement
  (450h) comes from whole-screen programmes with large controls. 458's prior should be modest: fid near the low end of
  its preregistered range. The prediction and bins stay as written; this is recorded before the entry exists.
- **Arc's H1 bank** has 18,399 control cells against 400 per target: the control-to-target noise ratio is
  sqrt(400/18,399) = 0.15, vs sqrt(50/500) ≈ 0.32 for a Broad line. The first reading here was that its programme
  carries far less control noise. **§3 corrects that:** the programme is a weak common signal, so even that ratio leaves
  about a third of it as control noise. Every single-control programme needs the control-halves test.
- **Method note for any future programme source:** measure its reliability with control halves, never with target
  halves.

## 3. The same test on Arc's H1 bank, the programme 458 ships (`h1halves462.py`, 462c)
`data/harness/h1-full.h5ad` is the `h1_2025` bank's source: 255 perturbations × 400 cells, 18,399 control rows. KA: the
rebuilt full programme equals the bank's P_down at r 0.99999998.

| split-half reliability of P_down | expressed genes (CPM ≥ 5, 10,788) | A / B / C blocks (2,185 each) |
|---|---:|---:|
| target halves vs the SAME full control | 0.793 | 0.884 / 0.862 / 0.875 |
| target halves vs DISJOINT control halves | 0.386 | 0.483 / 0.476 / 0.497 |

- By a variance decomposition, the full-cell bank programme is about 56 % signal, 33 % control-sample noise and 12 %
  target noise on expressed genes.
- Even 18,399 controls leave a third of a weak common programme as noise. The bank's control rows are 15,571 unique
  cells, some resampled, thinned to ~16k UMI.
- **458d** (`458/denoise458d.py`) rebuilds the programme from all public H1 cells: 300 perturbations × ~1,071 cells, 38,176
  controls, raw depth. Disjoint-control reliability **0.687** (shared 0.919). It flips ~21 % of 458's block signs and
  agrees more with CLEAN4 (A 0.650 vs 0.604). It replaces 458 in the 09-30 slot.

## 4. Negative control: 441's oracle gap is NOT control noise (`oracle462.py`, 462d)
441's oracle predicts a held-out h1-test target's gate signs from the leave-one-target-out mean of the other targets'
reference fold changes. All of those share the recipient's 18,400 control cells, so the oracle could be "predicting" the
control sample's per-gene deviation. Test: split the controls into halves, and take the oracle against one half and the
truth against the other.
- KA: the recomputed fold changes equal the scorer's reference table to float32 accumulation precision (median 1.6e-6,
  max 1.2e-5 in log2). The all-controls oracle reproduces 441's 0.625 / 0.632 (here 0.6245 / 0.6324, 58 targets, 64,931
  gate pairs).

| oracle vs truth | gate sign accuracy (unweighted / \|l\|-weighted) |
|---|---|
| both against control half A (shared) | 0.6201 / 0.6283 |
| both against half B (shared) | 0.6167 / 0.6226 |
| oracle from B, truth vs A (disjoint) | 0.6166 / 0.6226 |
| oracle from A, truth vs B (disjoint) | 0.6200 / 0.6283 |
| CLEAN4 (does not see the recipient's controls) | 0.5733 / 0.5713 (49.7 % of gate pairs covered) |

- Shared 0.6184 vs disjoint 0.6183 on average: no inflation.
- On gate genes, which are significant for their target, the target's own effect sets the truth's sign, not the
  control's small offset. So 441's gap (CLEAN4 0.573 → oracle 0.625) is honest recipient signal, and F18's prize stands.
- Control noise distorts programme ESTIMATES over all genes (§2–3), where most entries are weak. It does not distort sign
  agreement on significant genes. For a prior, that is the argument for estimating the programme from as many cells
  and controls as possible (458d), not a correction to the oracle.

## 5. 458d vs 458's programme on H1's own gate (`h1prog462.py`, 462e; same line and lab, an upper reference)
h1-test gate signs (58 targets) from P_down over the OTHER H1 perturbations, with h1-test's 85 targets excluded:

| programme | rows | gate accuracy (unweighted / \|l\|-weighted) |
|---|---:|---|
| 458's bank (`h1_2025`: 400 cells/target, thinned) | 170 | 0.6134 / 0.6073 |
| **458d's (`h1_2025_full`: all cells, 38,176 controls)** | 215 | **0.6400 / 0.6384** |
| LOTO oracle of the other h1-test targets (462d) | 84 | 0.6245 / 0.6324 |
| CLEAN4 | — | 0.5733 / 0.5713 |

- Denoising adds +0.027 on the same line. The clean programme beats even the LOTO oracle, because it averages more
  cells and targets.
- Its share of shared controls cannot explain this: 462d measured no inflation from them.
- For A/B/C the question is cross-line (the board, 09-30). Relative to 458, 458d carries about 1.24× the signal over
  chance ((0.640 − 0.5) / (0.613 − 0.5)).

## 6. The fid-like reading on H1 (`fidlike462.py`, 462f): 458d reaches the oracle
Sign agreement with the reference lfc on ALL block genes (2,200 top control-CPM, own gene out, lfc ≠ 0), 85 targets,
187,000 pairs. This is 441 §3's fid-like column. fid is precision over a large call set, so this is its proxy.

| programme | accuracy |
|---|---:|
| CLEAN4 (KA: 441's 0.5244) | 0.5244 |
| 458's bank programme | 0.5442 |
| **458d's full-cell programme** | **0.5545** |
| LOTO mean of the other h1-test targets | 0.5544 |

- On its own line and lab, the clean Arc programme matches the recipient oracle, +0.030 over CLEAN4. That is the fid
  territory the board leaders occupy (~0.55 raw; ours 0.514).
- For A/B/C the gain scales with the cross-line transfer fraction. At 462's one-lab third, with the prior applied to
  ~80 % of block genes (q80) and the block ~85 % of calls, that is roughly +0.007 raw fid, about +0.004 overall.
- 458d's read on 09-30 is the DE channel's most informative read of the week.

## 7. A CLEAN4 + Arc-full blend, priced on H1 for 458's "partial" bin (`blend462.py`, 462g)
z = w·z(Arc) + (1 − w)·z(CLEAN4) on P_down − 0.5; Arc alone where CLEAN4 is empty. w 0 reproduces CLEAN4 (KA).

| w (Arc) | fid-like (all block genes) | gate |
|---:|---:|---:|
| 0 | 0.5244 | 0.5733 |
| 0.25 | 0.5361 | 0.6103 |
| 0.5 | 0.5456 | 0.6267 |
| 0.75 | 0.5514 | 0.6369 |
| 1 | 0.5545 | 0.6400 |

On H1 (Arc's own line) more Arc is monotonically better. On A/B/C the weight is the board's to set. If 458d reads
partial (+0.0005 … +0.0015 fid), the natural next arm is w 0.5 (`experiments/458-arc-lab-programme/`, a new artifact
in artifact458's format).

## 8. CRISPRi with cells: the Weissman K562 and RPE1 essential screens (`weissman462.py`, 462h; `weissxaira462.py`, 462i)
scPerturb's single-cell Replogle 2022 K562-essential (310k cells, 10,691 controls) and RPE1-essential (248k cells, 11,485
controls): the same lab and library in two lines, 2,055 matched targets, 7,226 genes (all ≥ 5 CPM in both).

| reliability of P_down | K562 | RPE1 |
|---|---:|---:|
| target halves, shared control | 0.977 | 0.995 |
| target halves, disjoint control halves (honest) | **0.741** | **0.962** |

- Cross-line with no shared noise (K562 half vs RPE1 half): **0.453**; full vs full 0.489.
- Honest cross-line / within-line = **0.53**, against 0.35 on the Broad-16 knockouts (§2).
- With ~2,000 targets and ~11k controls, control noise is a smaller share than in the Broad or H1-harness designs, but
  still visible: K562 0.98 → 0.74.

**Is it the lab or essential-gene biology? (462i)** The same comparison on the 1,309 essential targets that every bank
carries (6,581 genes, full-cell programmes; cross-bank r is honest but attenuated):

| pair type | pairs | mean r |
|---|---:|---:|
| same line, same lab (K562-ess vs K562-GWPS) | 1 | 0.702 |
| different line, same lab (Weissman K562/RPE1/GWPS; Xaira HCT116/HEK293T) | 3 | 0.431 |
| same line, different lab (Weissman K562 vs Satija K562) | 2 | 0.386 |
| different line, different lab | 9 | 0.376 |

- On essential targets the shared programme is mostly biology: other labs' lines agree nearly as well as the same lab's
  (0.38 vs 0.43). Knocking down an essential gene has a strong common response.
- 450h's large lab effect (0.54 vs 0.18) came from broad target sets, where most responses are weak and the programme
  is dominated by lab and protocol.
- Arc's panel is non-essential, which is 450h's regime, so 458d's premise stands. It cannot be settled on public data.
  The board decides (09-30).

## 9. A non-essential CLEAN4? No support on H1 (`essnoness462.py`, 462j)
Three of CLEAN4's four banks are essential screens, so it is largely the essential-knockdown programme. Split by the
Replogle essential list, on H1's fid-like and gate rows:

| programme | rows | fid-like | gate |
|---|---:|---:|---:|
| K562-GWPS, all / essential / non-essential | 9,705 / 1,985 / 7,720 | 0.5204 / 0.5219 / 0.5178 | 0.5506 / 0.5631 / 0.5332 |
| Orion HCT116, all / essential / non-essential | 17,683 / 1,783 / 15,900 | 0.5193 / 0.5188 / 0.5181 | 0.5395 / 0.5443 / 0.5328 |

The non-essential programmes are no better, and on the gate they are worse. H1's truth is an ESC, so this is a weak
proxy, but there is nothing here to justify an arm. Closed. CLEAN4 stays.

## 10. Arc's programme from its non-essential rows only? No (`arcness462.py`, 462k)
36 % of Arc's 275 off-panel H1 perturbations are Replogle essentials; the A/B/C panel has none. On h1-test's 64
non-essential targets (fid-like / gate):

| programme | rows | non-essential targets | essential targets (21) |
|---|---:|---|---|
| all rows (458d) | 215 | **0.5438 / 0.6266** | 0.5874 / 0.6489 |
| non-essential rows | 137 | 0.5399 / 0.6217 | 0.5762 / 0.6290 |
| essential rows | 78 | 0.5410 / 0.6060 | 0.5815 / 0.6389 |

More rows beat a better-matched subset, even on non-essential targets. 458d keeps all rows. (215 rows here because
h1-test's own targets are excluded.)

## 11. Where a programme prior's reach cost comes from (450's H1 Pan run, per target)
Pan's same-line programme on H1 (450: fid +0.0097, reach −0.041 raw): per-target reach vs the shipped CLEAN4 run, 82
targets, 3 seeds.
- 12–15 targets change. **The four largest losses carry ~100 % of the mean loss in every seed**: CNN3 −1.0 and
  COL6A1 −1.0 (n_real 2 each), ZNF664 −0.71 (n_real 8), AASS −0.63 … −0.81 (n_real 17).
- The mechanism: in a small-gate target, the prior flips the sign of a top-ranked reference DE gene that the design had
  right. reach needs 90 % purity from the first position, so one flip takes it from 1 toward 0.
- A panel of 300 targets per context has more small-gate targets than H1's 82, so 458d can pay the same cost on the board.
  That is what the reach contingency in `458/PREREGISTRATION.md` (11:31Z) guards. A reach-protecting variant would keep
  the design's sign on each target's top-ranked genes and apply the programme to the rest. Build it only if 458d reads
  fid-positive and reach-negative.

## 12. The H1 six-member profile of Arc's own programme as the block prior (462l, `drive462.sh`)
435's production harness (`rc_emit_k` edb5b8dc, coarse MDE), with only the block prior swapped. KA: clean4 through the
wrapper equals 435 seed 0 bit for bit (PASS). Raw deltas vs the shipped CLEAN4 run, same seeds, pds/mse pinned:

| prior (h1-test targets excluded) | seeds | fid | nmae | reach | jac |
|---|---:|---:|---:|---:|---:|
| **Arc all cells (458d's bank)** | 3 | **+0.0198** | **−0.0148** | −0.0388 | −0.0026 |
| Arc harness bank (458's) | 1 | +0.0128 | −0.0111 | −0.0249 | −0.0031 |
| Pan H1 (450, same line, other lab) | 3 | +0.0096 | −0.0086 | −0.0407 | −0.0001 |

- On its own line, the clean Arc programme doubles Pan's fid and nmae gains. Denoising adds about half again over the
  bank programme.
- **The reach loss is structural, not a quality signal.** The same four small-gate targets carry it for Arc's
  programme and for Pan's, in every seed: CNN3 −1.0 and COL6A1 −1.0 (n_real 2), ZNF664 −0.71 (8), AASS −0.62 … −0.75
  (17).
  - Their few true DE genes move against the line's common programme. Any per-gene prior flips them, and the design's
    weak-but-right sign is unprotected (below q80).
  - The small gains (OAZ1, TOR1A, RFX5 +0.10 … +0.14) do not offset them.
- Scaled roughly on H1's spans (fid 0.30, nmae 0.60, reach 0.90, jac 0.37): fid +0.066, nmae +0.025, reach −0.043, jac
  −0.007. The net is positive, about +0.007 overall-equivalent, even at the oracle's quality.
- For 458d on the board: expect fid up and nmae down (better), with a reach cost of this kind. Judge it on fid and
  overall, as preregistered, not on reach alone.

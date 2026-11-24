# 814 PROPOSAL-2: the lab's panel draw rule. A = the corrected rule (no judge needed); B = A/B/C-fitted weights (judged)

Frozen 2026-10-16 ~11:55Z by the 814 worker. It answers `VERDICT-1.md` (REVISE, text only), items 1–5, and gives two costed
variants. Class: **instrument**: it prices no candidate. Evidence: `RESULTS.md` (cited as R§). The text work used no
compute beyond file reads. One small listing was made from files on disk (§1.2), and its sha is below.

## 0. Summary
- **Part A is a recorded rule that arc800 broke, now corrected in DESIGN §3.1 by the lead's 11:48Z amendment.** It needs no
  validation to be correct. §1.2 has one fix to that amendment: use the screens' full **libraries**, not the filtered
  `data/donors/*/perts.txt`. That changes the pool from 4,355 to 4,348.
- **Part B (propensity-decile weights) is fitted on A/B/C's 300 targets.** A validation PASS on ABCp shows that B emulates
  A/B/C's panel on public lines. It does not show that B transfers to another Arc panel. That test is on D/E/F's panel on
  10-22 (§3.2).
- **Variants** (§4):
  - **A-only:** the lab's default from now on. The lead may draw it without a judge.
  - **A+B:** needs the judge's GO.
  - Each costs ≤ $5.

## 1. Novelty against the record (JUDGE item 4)
### 1.1 Part A is FINDINGS §8.1g's rule
- `docs/LOG.md` l.200–212 (2026-08-21): A/B/C's panel has 0 hits in each of Replogle K562-essential (2,058), RPE1 (2,394),
  Nadig HepG2 (2,394) and Nadig Jurkat (2,394), against about 31–36 each by chance. The verdict there: "deliberately built to
  avoid common-essential genes".
- `docs/FINDINGS.md` §8.1g (l.1647–1697): "a synthetic panel must satisfy ARC'S STATED SELECTION CRITERION, not merely be
  disjoint".
- arc800 applied the K562 list only (`800-checker/DESIGN.md` §3.1, l.52), so it broke the rule. 814 re-found the gap from
  the data (R§2: A/B/C 0 vs 16.9 ± 3.7 per arc800 draw).
- A applies the recorded criterion. It is not fitted on anything and needs no validation to be correct. What A's folds
  will show is how much of the over-read the rule fix removes (§4).

### 1.2 The HepG2 gap: closed from the source files (file reads only, no download)
The four screens share one library. I read the perturbation category lists of the source h5ad files (obs categories only;
X never opened):

| screen (file) | targets (excl. controls) | relation |
|---|---:|---|
| Replogle RPE1 (`data/external/scperturb/ReplogleWeissman2022_rpe1.h5ad`, `perturbation`) | 2,393 | |
| Nadig HepG2 (`data/external/nadig2025/hepg2_raw_singlecell.h5ad`, `gene`) | 2,393 | equals Jurkat except **CDCA8** (HepG2 only) |
| Nadig Jurkat (`data/external/nadig2025/jurkat_raw_singlecell.h5ad`, `gene`) | 2,393 | equals HepG2 except **INCENP** (Jurkat only) |
| union | **2,394** | RPE1's library ∪ {CDCA8, INCENP}, i.e. one essential-gene library |

- Nadig used the same essential-gene library in HepG2 and Jurkat, minus one target each. Replogle's RPE1 used the same
  library. This matches the LOG's identical counts (2,394 = 2,393 targets plus the control label).
- **The fix to DESIGN §3.1's amendment.** `data/donors/rpe1/perts.txt` and `nadig_jurkat/perts.txt` are the donor banks.
  They are filtered to targets with ≥ 20 cells (meta.json: 189 and 123 targets skipped), so their union holds 2,351
  targets and **misses 43 library targets**. Low-cell targets are the most lethal knockdowns, the ones Arc's selection
  removes first.
  - 7 of the 43 are in arc800's null pool: CDK7, CEBPG, CMTR2, POLR3A, POT1, RBM22, ZNF830. Of the 43 missing
    library targets, L4 holds 1 and L5 holds 2.
  - The rule should remove the **library union**: `814/work/flags/esslib.txt`, 2,394 genes, sha256
    `24fffbfee3d3fa3b038d1a7f3f665772ec171a2341ac35c03dafe6eb0a3d55e7`, written from the three source files (sizes
    1,293,034,868 / 5,614,460,941 / 9,366,490,264 bytes).
  - A/B/C ∩ library union = 0. The K562 essential list (`b1e52795…`) stays removed as before.
- I did not edit DESIGN.md (not my file). The lead decides whether the amendment cites `esslib.txt`.

### 1.3 Part B's precedent: 522 §4b, and what prop10 adds
- `522-october-bulk-key/RESULTS.md` §4b (l.243–262) re-drew panels to A/B/C's mix of mean robust strength over the fold's
  own sources, in quintile bins. That is 814's `resp5` construction. 522 also ran an "Arc-like" mode: essential-like
  targets removed, then re-drawn to A/B/C's mix.
- **Result: key tails stayed heavy, cap hi 9.5–10.9 against A/B/C's 6.15–7.03.** 814's preview agrees. The key class stays
  low under every rule: 0.17 → 0.24 (resp5) / 0.26 (prop10); 705 is at the bottom in 3/3 cells (R§3). This is the prior
  for "B does not fix the key class".
- **What prop10 adds beyond 522's strength mix** (`work/flags/prop10.json` `1ec590e2…`; standardised coefficients):
  - CD4 guide-level knockdown, the largest coefficient (+0.87);
  - weak knockdown in K562 / HCT116 (−0.29);
  - three coverage terms: Pan 2025 H1 (+0.20), VIPerturb K562 (+0.12), VCC-2025 H1 control expression (+0.11).
  - Responsiveness (+0.50) is 522's axis.
  - What is new is the KD and coverage axes, the measurable side of Arc's "≥ 80 % knockdown" and "representation in
    public datasets".
  - prop10 also differs from 522 in its essential definition (the screen library, not 522's essential-like flag) and in
    keeping arc800's coverage × CD4-bin strata.
- **Relation to 238's efp / efr** (`238-essential-free-panel/RESULTS.md` §1):
  - efp and efr removed the 2,021-gene K562 list (and A/B/C) from donor-coverage pools and drew uniformly. They had no
    strength strata, no basal floor and no screen-library cut.
  - So they still hold **9/300** RPE1 ∪ HepG2 ∪ Jurkat targets each (238's own table; the lead's 11:41Z count; the
    library-union count here is also 9 / 9).
  - v2-A is efp/efr's correction carried through, with arc800's coverage and strength strata and its floor kept.
  - v2-B adds what neither 238 nor arc800 had: the within-cell knockdown and coverage weights.
- No closed family in FINDINGS §2 covers a panel-draw rule.

## 2. The rule, pinned (JUDGE item 4: pool by sha)
- **Pool.** arc800's null pool, reconstructed by `feat814.py` from 801's `work/draws.json` (sha
  `8c8e67d2…`), minus the library union.
  - 801's record holds the excluded pert_counts files at draw time (n_excluded 4,373), the basal floor
    13.412203349132618 and the stratum counts. Pool0 = 5,125. KA PASS: L1–L6 re-drawn bit-exact.
  - The sealed L7/L8 lists are removed as an exclusion only.
  - Null pool: 4,525 genes, sha256 of the sorted, newline-joined list `5d1ca325…6534`.
  - **v2 pool = null pool − `esslib.txt` (`24fffbfe…`) = 4,348 genes, sha256
    `5eb2568ecb0e0c50c5a8ad9502dccf50a35babb01bb3a72bc077032c02c52604`.**
  - For the record: the donor-list cut that sim814 used gives 4,355 genes (`e3b6b04c…f183d`). That pool is superseded by
    §1.2.
- **KA at draw time: abort on any failure.**
  - The script recomputes the pool from the same inputs. A pool sha other than `5eb2568e…` is a **KA FAIL and aborts**;
    it must never silently become a new rule.
  - It also asserts: 0 K562-essential; 0 library; 0 overlap with any excluded panel (L7/L8 included); every cell's count
    equals its quota.
- **Strata.** arc800's: coverage class × CD4 Rest bin, with A/B/C's quotas after draw801's moves (the L panels' cell
  counts).
- **Part B weights.** Within each cell, w(decile of s) = A/B/C share / pool share, with the decile edges frozen as the pool
  deciles of s over the 4,355 rows used in the fit.
  - s is prop10's frozen score: coefficients, means and sds in `prop10.json` `1ec590e2…`, features from
    `features814.csv.gz` `bdee8c2f…`. The deciles per gene are in `prop10.tsv` `66412560…`.
  - The weights are recomputed on the 4,348 pool by the same formula. The 4,355-pool weights were p0 0.00 … p9 3.46; the
    7 removed genes change them only slightly. Both pools put 20.0 % in deciles 9–10.
- **The pre-recorded V2p panel shas** (`sim2/sim814.json`, `546e658f…` … `5939c5bc…`) were drawn from the superseded
  4,355 pool, so they are no longer the KA. The new KA is a cross-check: `sim814.py --indep6`, with the pool switched to
  `esslib.txt`, must give the same six panel shas as the draw script. Both run before any fold is built. That is a
  seconds-long local or EC2 job, done after GO, not now.

## 3. Predictions and falsifiers
### 3.1 The validation bar (the lead's, fixed in the mandate before any v2 fold exists)
- panelrank800 is run with the six new panels × three held lines in place of L1–L6, on the same seven arm pairs.
- **PASS:** ABCp mean percentile ≥ 0.35 and bottom count ≤ 6/21.
- **PARTIAL:** PASS on ABCp, but LABC < 0.30 or LABC bottom > 9.
- **FAIL:** anything else.
- Each class's percentile (key / direction / shape) and the re-measured over-read factors, with CIs, are reported either
  way.

### 3.2 Two 10-22 falsifiers, read on D/E/F's 300 targets before any v2 read informs an October or 6a decision
| prediction | support | partial | fail | on fail |
|---|---|---|---|---|
| **P-A:** D/E/F targets in the essential library union (`esslib.txt`) | ≤ 3 | 4–7 | ≥ 8 (arc800 expects 16.9 on the donor lists) | Arc changed its criterion. A's premise does not carry to D/E/F, so v2-A factors are not applied to D/E/F decisions, and the lab prices on D/E/F's own panel (708's route) |
| **P-B:** D/E/F targets in prop10's frozen pool deciles 9–10. Scored with prop10.json's coefficients and the frozen edges; D/E/F's features come from feat814's bank code on D/E/F's target list (one EC2 run, ~$0.02) | ≥ 130 | 96–129 | ≤ 95 | **v2-B factors are withdrawn for October decisions; A-only stands** |

How the P-B bins were set:
- A/B/C has 164/300 in deciles 9–10, against 60 for a pool-uniform draw.
- The binomial sd is 8.6 at A/B/C's rate and 6.9 at 0.2. Support allows a second Arc draw to land about 4 sd below
  A/B/C. Fail is within about 5 sd of the pool rate.
- On P-B partial, v2-B reads are reported only beside the A-only read and decide nothing.

P-A carries little information, since it is near-certain by Arc's stated criterion. P-B is the real out-of-sample test of B.

### 3.3 Same-panel branch (JUDGE item 3)
If D/E/F reuses A/B/C's 300 targets (DEF-RUNBOOK's same-panel branch), **both P-A and P-B are void**. They would be
satisfied trivially and say nothing. In that branch ABCp and LABC are D/E/F's own panel on held-out lines, and they are the
deciding instruments. v2 reads are then reported only.

## 4. The two variants, costed
| | **A-only** (lab default; no judge needed) | **A+B** (judge's GO) |
|---|---|---|
| pool | v2 pool, 4,348 (`5eb2568e…`) | the same |
| within-cell draw | uniform (draw801's `rng.choice`) | prop10 decile weights (§2) |
| six panels | **sequential, disjoint from each other**, draw801's deficit moves; seeds 814201–814206. sim814's R1 sequential run on the 4,355 pool needed 0 moves through panel 5 and 25 at panel 6 | **independent** draws, seeds 814101–814106. Disjoint weighted draws exhaust the responsive targets by panel 3 (R§4); expected pairwise overlap ~53/300 |
| overlap with L1–L6 | allowed, about 180–187 per panel. A draw that excludes L1–L6 is infeasible: C3\|Q4 keeps 7 members for a quota of 47 (R§4). This is the one departure from draw801's "minus every earlier lab draw" | the same |
| fitted on A/B/C | nothing | prop10's 6 coefficients and 10 weights; chosen as the best of 12 previews |
| preview (thinned L folds; not a validation) | 0.357 / 5 bottom (LABC 0.303 / 5) | 0.412 / 2 (LABC 0.412 / 2); best-of-12, so it reads high |
| folds | 801's recipe: CD4 banks, day-recipe operator fits, 18 dumps (801's run g did 18 for $1.98) | the same |
| arms | ka804's L set (495, 479p2g, 479sp2, 705, 707), 810 sgu1, 811 c0_l25 (~$1.5) | the same |
| **cost** | **≤ $5** | **≤ $5**; ≤ $9 for both, since the fold code is shared and only the panels differ |
| tries | the rule fix is not a try; its read is try 1 of the A family | try 1 of the B-validation family; the design used 5 simulated rules and 15 + 4 previews (R§5) |
| 10-22 gate | P-A | P-A and P-B |

**Order if both run:** A-only first, as the default working set. A+B second, judged against A-only's read: B is useful only
if it moves ABCp's percentile above A-only's.

## 5. Decision role (JUDGE item 5)
- **A PASS makes the v2 folds a deciding surface beside ABCp and LABC for key and direction candidates. It does not replace
  them.**
  - All three are reported for every such candidate.
  - A candidate needs a non-negative read on each (the checker's lower-bound rule on each surface after best-of-N).
- **Per-class over-read factors are re-measured on v2 and applied** to v2 reads: key, direction and shape, each v2 vs ABCp
  with its CI, replacing the L-fold ×3.8 / ×4.7.
- **PASS on the mean does not certify a class.** A class whose ABCp percentile on v2 stays below 0.25 keeps its measured
  factor > 1 and is labelled so. The key class previews at 0.17–0.26, with 705 at the bottom under every rule, so it is
  expected to keep one.
- **B's PASS clause.** A PASS of A+B validates the **emulation of A/B/C's panel** by a rule fitted on A/B/C. It does not
  validate B's transfer. A+B factors enter October or 6a decisions only after P-B supports on 10-22 (§3.2). Until then
  they are reported beside A-only's.
- **FAIL:**
  - A-only FAIL: A stays the draw rule, since it is the recorded criterion and correct regardless. Key and direction stay
    on ABCp + LABC (PHASEB), and the residual over-read is recorded as a property of the public lines.
  - A+B FAIL: B is closed for this lane; a re-proposal needs a new mechanism.
- **Wear.** v2 folds share about 60 % of their targets with L1–L6, so part of L's wear carries over. No new candidate is
  priced on v2 before its bar is read. The validation arms are known answers already priced on L, and 811's and sgu1's
  L-selected gains bias the read against PASS, not toward it.

## 6. Files this proposal relies on (sha256)
| file | sha256 |
|---|---|
| `work/flags/esslib.txt` (new, §1.2) | `24fffbfee3d3fa3b038d1a7f3f665772ec171a2341ac35c03dafe6eb0a3d55e7` |
| `work/flags/prop10.json` / `prop10.tsv` | `1ec590e2…` / `66412560…` |
| `work/run1/features814.csv.gz` | `bdee8c2fedfab527c774bcc7648ff597706ec3f03af2c68986697fefe2ba8cd5` |
| `experiments/801-lab-folds/work/draws.json` | `8c8e67d2…` |
| v2 pool (sorted, newline-joined list) | `5eb2568ecb0e0c50c5a8ad9502dccf50a35babb01bb3a72bc077032c02c52604` (4,348) |
| superseded donor-list pool | `e3b6b04c110b6c9a893ccf0ef4e31bd07b210abba99de47ebea1a7af803f183d` (4,355) |
| null pool | `5d1ca325ce6f3f4b2516fc1c3fbb81b1c8e3b2817ae41176ed4301ddd28e6534` (4,525) |
| `data/donors/rpe1/perts.txt` / `nadig_jurkat/perts.txt` / `k562_essential/perts.txt` | `4377d2dc…` / `f5a5c1bb…` / `b1e52795…` |

# VERDICT — 814 draw rule v2 ("arc814")

Judge run, written 2026-10-16 11:46Z (by the clock). Judged: the frozen copy
`experiments/819-judge/814-panel-rule-20261016T114114Z/PROPOSAL-frozen.md`, sha256
`1e7e89dfc96e2cacbf23a81c6218160531ec05a4c4cf83367f3e8db7907cf9bf` (re-hashed by the judge; byte-identical to the live
`experiments/814-panel-rule/PROPOSAL.md` at 11:46Z). Rules: `experiments/800-checker/JUDGE.md`. This is a first submission:
there is no earlier `VERDICT-*.md` in `experiments/814-panel-rule/`.

Read: the 814 MANDATE.md and RESULTS.md; PHASEB.md; DESIGN.md §3.1 and §4; `docs/LAB-PANEL-REPORT-2026-10-16.md`;
`docs/LAB-LOG.md` (both known-answer reports and the 10:41Z, 10:46Z, 10:51Z, 11:14Z and 11:41Z notes); `docs/FINDINGS.md` §1–§2,
§8.1g and H.1; `docs/LOG.md` ~l.200; `experiments/522-october-bulk-key/RESULTS.md` §4a–§4b; `experiments/438-…/RESULTS.md`
l.8–10; `experiments/800-checker/work/panelrank800.txt`. The only computation was a set lookup on small text files (below).

## Verdict: **REVISE**

The proposal is cheap (≤ $5), honest about its tries, and has a falsifier. The essential-screen cut (part A) is correct and
should be in any rule the lab draws. But the proposal leaves out the novelty section JUDGE.md requires (item 4). That
section matters here: part A is a rule the record already states and arc800 failed to apply, and part B's closest
construction (re-drawing to A/B/C's source-strength mix) was already run by 522 and left the key tail heavy. Part B is also
fitted to A/B/C's 300 targets, and the validation bar reads A/B/C again, so the validation is in-sample for B and the
proposal does not say so. Every missing item is text; none needs compute. One resubmission.

## The checks (JUDGE.md "The judge's checks")

### 1. Is the mechanism closed already? **Not closed, but its main finding is already in the record, and the proposal does not cite it.**
- **Part A (minus the RPE1 ∪ Jurkat essential-screen targets) is FINDINGS §8.1g's rule.** It is not a new discovery:
  - `docs/LOG.md` l.200–212: A/B/C's panel has 0 hits in Replogle K562-essential (2,058), RPE1 (2,394), Nadig HepG2 (2,394)
    and Nadig Jurkat (2,394), against ~31–36 each by chance. Verdict there: "deliberately built to avoid common-essential
    genes".
  - `docs/FINDINGS.md` l.1689: "**THE RULE: a synthetic panel must satisfy ARC'S STATED SELECTION CRITERION, not merely be
    disjoint.**" l.1695: "Before building on a synthetic panel … grep the log for it first."
  - arc800 applied only the K562 list (`experiments/800-checker/DESIGN.md` l.52). The lead's mandate repeated that narrowing
    ("Essential genes are NOT the difference … 0 of the 2,021 K562-essential perts", `814/MANDATE.md`, "Already known").
    814 corrected it from the data, which is to its credit. But the proposal presents it as a finding needing validation
    (`PROPOSAL-frozen.md` l.27–30). It is a recorded rule that the lab's instrument broke.
  - The judge's lookup confirms the files: RPE1 ∪ Jurkat perts = 2,351 genes, 365 beyond the K562 list; A/B/C ∩ = 0
    (`data/donors/rpe1/perts.txt` `4377d2dc…`, `data/donors/nadig_jurkat/perts.txt` `f5a5c1bb…`,
    `data/donors/k562_essential/perts.txt` `b1e52795…`, `data/challenge/pert_counts.csv`).
- **Part B (responsiveness-weighted draws) has a direct precedent: 522 §4b**, which the proposal does not cite.
  - 522 re-drew panels "to A/B/C's mix" on the mean robust strength over the fold's own sources, in quintile bins
    (`experiments/522-october-bulk-key/RESULTS.md` l.243–257). That is 814's `resp5` construction.
  - 522 also ran it with essential-like targets removed ("Arc-like" row, l.254).
  - Result: those panels "still carry heavier tails than A/B/C: hi 9.5–10.9" (l.260), against A/B/C's 6.15–7.03 (l.212–214).
  - This agrees with 814's own preview, where key stays low after every rule (0.17 → 0.24–0.26; R§3). It is the prior for
    "B does not fix the key class". prop10 differs from 522 by adding CD4 guide-level KD (the largest coefficient, +0.87),
    weak KD and three coverage features. The proposal must say that this, not responsiveness, is what is new.
- 238 / 438's essential-free panels (efp, efr; `438/RESULTS.md` l.8–10) are a further precedent for an essential-free draw.
  The proposal should say how v2 relates to them. As far as the record shows, they used the K562 list and DepMap rather
  than the RPE1 / Jurkat screens.
- No closed family in FINDINGS §2 (l.795–926) covers a panel-draw rule. Not KILL.

### 2. Prediction with falsifier and bins before any number? **Yes for the validation; partial for part B out-of-sample.**
- The bar was written by the lead before any v2 fold exists: mean percentile ≥ 0.35 and bottom ≤ 6/21 on ABCp, with
  LABC < 0.30 or > 9 bottom → PARTIAL (`PROPOSAL-frozen.md` l.89–92). PASS and FAIL consequences are at l.93–96. Good.
- A free 10-22 prediction exists (l.107): ≤ 3 essential-screen targets, and a responsiveness median above the arc800 null's
  95th percentile.
  - Part A's half of it is near-certain by Arc's stated criterion (FINDINGS l.1680), so it carries little information.
  - The responsiveness half is a real test of B, but of one feature only. B's weights come from a six-feature score.
  - The prediction is void if D/E/F reuses A/B/C's panel (DEF-RUNBOOK keeps a same-panel branch), and the proposal does not
    say so.

### 3. Measurable on the checker's surfaces? **Yes.**
It prices no candidate. It re-reads panelrank800 with unchanged scripts and arms (l.84–88), in the same bulk units as
panelrank800 (`800-checker/work/panelrank800.txt` l.1–24). DE members are not involved: the over-read is a bulk key /
direction factor (PHASEB l.24–31).

### 4. Does it use seen data as evidence of gain? **Not in the KILL sense, but part B's validation is in-sample, and the proposal must say so.**
- No board read, training line or sealed reserve enters. L7/L8 and R6/R7 are used only as exclusions (RESULTS §5), as
  draw801 does and as the mandate allows.
- What B is fitted on:
  - prop10's six coefficients and ten decile weights (`work/flags/prop10.json` `1ec590e2…`) are fitted on A/B/C's 300
    against the pool;
  - B was chosen as the best of 12 preview configurations, ranked by ABCp's percentile on the L folds (RESULTS §3, §5);
  - the v2 folds share about 183 of 300 targets with L1–L6 (l.45).
- So passing the bar on ABCp mainly confirms that the thinning preview (0.412 / 2) survives a real redraw with refitted
  operators. It does not show that B transfers to an Arc panel B was not fitted on.
- The only such panel is D/E/F. The judge checked whether VCC-2025 H1 could serve: it cannot. It holds 71 of 255 in the
  essential screens and 57 K562-essential (same files; `data/donors/h1_2025/perts.txt`), so Arc's 2025 rule differed.
- B's PASS consequence (l.93–94) needs that label and an out-of-sample gate before any October / 6a use.

### 5. Is the try budget honest? **Yes.**
- 5 simulated rules plus 15 + 4 preview configurations are counted. The preview is labelled best-of-12 (l.73–76, l.105–106).
- The validation is one preregistered configuration (V2p, panel shas pre-recorded, l.47–54), and a second rule is
  pre-declared as try 2 (l.95–96).
- Note in the proposal's favour: arms 811 c0_l25 (N 11) and sgu1 (N 8) were selected on L folds. Their worn-high L gains
  carry into v2 through the ~60 % overlap, which biases the read against PASS, not toward it.

### 6. Worth its cost? **Yes.**
- ≤ $5, against the instrument that decides every key and direction candidate.
- The over-read it attacks is 3.08× key and ≈ 6.5× direction (`docs/LAB-PANEL-REPORT-2026-10-16.md` l.43–44). That is far
  above any class noise.

### JUDGE.md required contents
| item | status |
|---|---|
| 1 Mechanism | present (l.29–30, l.39–41) |
| 2 Class | "instrument" (l.3, l.103), a fourth class that JUDGE.md does not list; acceptable for a lane that prices nothing |
| 3 Prediction / bins / falsifier | present (l.89–96, l.107) |
| 4 Novelty against the record | **missing**: no section; FINDINGS §8.1g, `docs/LOG.md` l.200, 522 §4b and 238 are not cited |
| 5 Cost and tries | present (l.100–106) |
| 6 Day-of path | present: D/E/F's panel is Arc's; untouched (l.103–104) |

## What is missing (fix all, resubmit once; text only, no compute)
1. **Novelty section** (JUDGE item 4).
   - Cite `docs/LOG.md` l.200–212 and FINDINGS §8.1g (l.1647–1697) for part A. State that A applies a recorded rule arc800
     broke (DESIGN l.52), so A needs no fitting and no validation to be correct.
   - Close the HepG2 gap. The LOG table lists Nadig HepG2 (2,394) as a fourth essential screen, and there is no HepG2 list
     under `data/donors/`. Either show that its targets are identical to Jurkat's (Nadig used one DepMap-essential set for
     both) or name the list used.
   - Cite 522 §4b (l.243–262) as the precedent for B's responsiveness matching, and its result (key tails still heavy,
     hi 9.5–10.9). Say what prop10 adds beyond it (CD4 guide KD, weak KD, coverage).
   - Say how v2 differs from 238's efp / efr essential-free draws.
2. **Label B as A/B/C-fitted, and gate it out-of-sample.**
   - State in the PASS clause that a PASS validates the emulation of A/B/C's panel, not B's transfer.
   - Add a second 10-22 falsifier with bins, read on D/E/F's 300 targets before any v2 read is used for a 6a decision:
     D/E/F's distribution over prop10's frozen pool deciles. A/B/C has 164 / 300 in deciles 9–10 against 60 expected; give
     the support / fail thresholds now.
   - State what happens on a fail (v2 factors withdrawn for October decisions; A-only stands).
3. **Same-panel branch.** State that both 10-22 predictions are void if D/E/F reuses A/B/C's panel.
4. **Pin the pool.** l.25 says "every pert_counts.csv on disk … as recorded at draw time", but the KA (l.47–54) requires
   sim814's exact pool.
   - Name the exclusion set by sha (801's `draws.json` `8c8e67d2…` plus the recorded list shas sim814 used, giving pool
     4,355).
   - Say that a different pool at draw time is a KA FAIL and aborts. It must not become a silent new rule.
5. **Decision role and the key class.** PHASEB l.24–31 (lead, 11:14Z) decides key and direction candidates on ABCp + LABC.
   - Say whether a PASS makes v2 a deciding surface or one reported beside ABCp / LABC.
   - The bar is a mean over 21 cells, and the key class previews at 0.26 (with 705 at the bottom in 3/3 cells under every
     rule; RESULTS §3). State that per-class over-read factors are re-measured on v2 and applied. Also state that a class
     whose percentile stays below ~0.25 keeps a factor > 1: PASS on the mean does not certify the key class.

## Note to the lab lead (not a condition)
Part A needs no judge to be correct: it is FINDINGS §8.1g's rule. It can enter DESIGN §3.1 by dated amendment now,
independent of B's fate. The sealed LZ folds (L7, L8) were drawn under arc800 and so carry about 17 essential-screen targets
each. As RESULTS §7 says, count them at the one access before the price is read.

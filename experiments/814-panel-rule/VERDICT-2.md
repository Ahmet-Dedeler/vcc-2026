# VERDICT: 814 PROPOSAL-2 (panel draw rule: A = the essential-library cut, B = prop10 decile weights)

Judge run, written 2026-10-16 ~11:58Z (from the clock: `date -u` 11:57Z). The judged file is the frozen copy
`experiments/819-judge/814-panel-rule-20261016T115334Z/PROPOSAL-frozen.md`, sha256
`06d085a2b6edb0319ba335a4733c90a074d3e5026f1946dc62c5c54fe87721c1`. The judge re-hashed it, and it is byte-identical to the live
`experiments/814-panel-rule/PROPOSAL-2.md`. Rules: `experiments/800-checker/JUDGE.md`. This is the resubmission after
`experiments/814-panel-rule/VERDICT-1.md` (REVISE, 11:46Z).

What I read:
- 814's MANDATE.md, VERDICT-1.md and PROPOSAL.md (v1);
- PHASEB.md and DESIGN.md §3.1, including the lead's 11:48Z amendment;
- `docs/LAB-LOG.md`: both known-answer reports (l.136–203), the lab queue (l.27–35), and the notes at 11:14Z, 11:41Z, 11:48Z and
  11:49Z (l.424–471);
- the charter's interface and timeline (`docs/LAB-2026-10-16.md` l.53–64);
- `docs/FINDINGS.md` §2 (l.799–930) and §8.1g (l.1647–1697);
- `522-october-bulk-key/RESULTS.md` §4b (l.243–262);
- the RESULTS.md hits for "propensity" / "panel draw" / "draw rule" (322 l.103 and 481 l.150 are panel-mix studies, not
  closures).

Lookups only, nothing heavier. Every sha the proposal cites in §6 matches the file on disk:
- `esslib.txt` `24fffbfe…` has 2,394 lines;
- `prop10.json` is `1ec590e2…`, `prop10.tsv` is `66412560…`, `features814.csv.gz` is `bdee8c2f…`, and 801's `draws.json` is `8c8e67d2…`;
- the donor perts files are `4377d2dc…` / `f5a5c1bb…` / `b1e52795…`.

Set checks:
- A/B/C's panel ∩ `esslib.txt` = **0**;
- `rpe1 ∪ nadig_jurkat` perts = 2,351 genes, all inside `esslib.txt`, which adds **43** more;
- 238's `ef-bundle` and `ef-bundle-r` panels each hold **9/300** library targets;
- `prop10.json` `abc_decile_share` puts 60 + 104 = **164** of A/B/C's targets in deciles 9–10.

## Verdict: **REVISE** (part B). Part A is not a judge matter and gets no separate 814 budget

PROPOSAL-2 fixes all five items of VERDICT-1 (table below). The problems are new ones, all in part B, plus a collision with
work the lead preregistered six minutes before this freeze:
1. **B's own claim has no falsifier.** The only lab test B gets is the shared absolute bar. A-only already previews at that bar
   (0.357 / 5 against ≥ 0.35 / ≤ 6). So an A+B PASS cannot tell B from A, and the B-over-A comparison that would
   ("B is useful only if it moves ABCp's percentile above A-only's", l.153–154) has no margin. It is also confounded by draw
   mode (B independent, 814's A-only sequential, l.143).
2. **No decision B could change has been named.** By its own §5 (l.166–168), B decides nothing until P-B is read on 10-22. From
   10-22 14:00Z the lab builds D/E/F-panel folds on Arc's actual panel (LAB-LOG l.32–33), and those should not carry the
   over-read (l.427). An emulator of Arc's panel is then superseded by the panel itself. As written, the ≤ $5 fold build buys
   a number that cannot enter any decision.
3. **814's A-only variant duplicates 801 follow-up 5** (LAB-LOG l.468–471, preregistered 11:49Z). That draw uses a different
   pool, a different draw mode and different seeds, and the proposal does not cite it.

JUDGE.md allows one resubmission, and this was it. Whether a third text round is worth having is the lead's call; my
recommendation is in the last section. The P-B prediction (l.123–129) is already frozen with bins here, and it stands as the
lab's free 10-22 test of B whatever happens to the fold build.

## VERDICT-1 items: all fixed

| VERDICT-1 item | where PROPOSAL-2 fixes it | status |
|---|---|---|
| 1a cite LOG l.200–212 + FINDINGS §8.1g; A is a recorded rule needing no validation | l.20–29 | fixed |
| 1b close the HepG2 gap | l.31–54: one library; HepG2 = Jurkat + CDCA8 (− INCENP); union 2,394; judge's lookup confirms the union ⊇ both donor lists, and A/B/C ∩ = 0 | fixed, **and it adds a correction** to DESIGN's amendment: the donor perts lists miss 43 library targets, 7 of them in the null pool (l.44–53) |
| 1c cite 522 §4b, its result, and what prop10 adds | l.56–70: CD4 guide KD (+0.87), weak KD (−0.29), three coverage terms; matches `prop10.json` | fixed |
| 1d how v2 differs from 238 efp / efr | l.72–78: K562 list only, uniform, no strata; 9/300 library targets each (judge's lookup: 9 / 9) | fixed |
| 2 label B A/B/C-fitted; second 10-22 falsifier with bins; fail consequence | l.11–13, l.166–168 (label); l.123–129 (P-B: support ≥ 130, partial 96–129, fail ≤ 95; fail → v2-B withdrawn, A-only stands) | fixed |
| 3 same-panel branch | l.133–136: P-A and P-B are void, and v2 reads are reported only | fixed |
| 4 pin the pool; a different pool is a KA FAIL | l.81–96: v2 pool 4,348 `5eb2568e…`; recomputed at draw time; any other sha aborts | fixed (the pool changed from 4,355 to 4,348 for a stated reason, §1.2) |
| 5 decision role; per-class factors; class < 0.25 keeps a factor > 1 | l.156–175 | fixed |

## The JUDGE.md checks

### 1. Is the mechanism closed already? **No.**
- No family in FINDINGS §2 (l.799–930) covers a panel-draw rule. The only panel mentions there (l.850 veto panel, l.885–887
  Testh basis, l.915 238's "gain is common-essential") are about other constructions.
- **Part A** is FINDINGS §8.1g's rule (l.1647–1697), now in DESIGN §3.1 by the lead's 11:48Z amendment. It is a rule
  restored, not an idea under judgement (proposal l.8–10, l.28).
- **Part B**'s nearest precedent is 522 §4b (l.243–262): re-drawing to A/B/C's strength mix, essential-like removed, left the
  key tails heavy (hi 9.5–10.9 against A/B/C's 6.15–7.03). The proposal cites it and names what is new: the KD and coverage
  axes (l.63–69). That differs enough not to be a re-named closure.
- The prior from 522, and from 814's own preview (the key class stays at 0.17–0.26 under every rule, l.164–165), is that B
  does not fix the key class. The proposal says so itself.

### 2. Does the prediction have a falsifier and bins written before any number? **For A and for P-B, yes. For B's in-lab claim, no.**
- The validation bar (l.111–117) was written by the lead before any v2 fold existed (MANDATE l.47–48), and so was P-B
  (l.123–129). Both have bins.
- **The gap is B's increment over A.** The absolute bar does not discriminate:
  - A-only previews at 0.357 / 5 bottom (l.146), which already meets ≥ 0.35 / ≤ 6.
  - So a PASS of A+B is expected whether or not B adds anything.
- The only B-specific statement is l.153–154: "B is useful only if it moves ABCp's percentile above A-only's". It has no
  margin and no bins. It has to clear three things:
  - B's own selection bias: best of 12 previews, and the preview gap 0.412 − 0.357 = 0.055 is a best-of-12 figure;
  - the panel-sampling noise of a 21-cell mean percentile, where each cell ranks ABCp among 7;
  - a comparator confound. 814's A-only uses sequential disjoint draws and A+B uses independent draws (l.143), so the
    comparison changes two things (AGENTS 11.1 rule 5).
- Under JUDGE.md this is a REVISE: "Does the prediction have a falsifier and bins written before any number? If not, REVISE."

### 3. Measurable on the checker's surfaces? **Yes.**
It prices no candidate. It re-runs panelrank800 with the same arms in bulk key / direction / shape units (l.112, l.148). DE
members are not involved: the over-read is a bulk factor (PHASEB l.24–31; VERDICT-1 §3).

**P-B's mechanics are not fully pinned yet:**
- `prop10.json` holds coefficients, means, sds and `abc_decile_share`, but **no decile edges**, and `prop10.tsv` holds only
  decile labels. D/E/F targets outside the 4,355 fit rows need numeric edges (l.99–101, l.123).
- The edges can be recomputed from sha-pinned inputs, but they should be recorded as numbers before 10-22.
- The proposal also does not say how a D/E/F target with missing bank features is counted toward the 300.

### 4. Does it use seen data as evidence of gain? **No (not KILL). B is fitted in-sample, and the proposal now labels that.**
- prop10's 6 coefficients and 10 weights are fitted on A/B/C's 300, and B was picked as the best of 12 previews on ABCp
  (l.145, l.150).
- The proposal states that an A+B PASS validates emulation, not transfer, and that only P-B on D/E/F's unseen panel can
  validate transfer (l.11–13, l.131, l.166–168).
- No board read and no sealed reserve enter. L7/L8 are an exclusion only (l.86).

### 5. Is the try budget honest? **Yes for B, no for the A family.**
- B: 5 simulated rules plus 15 + 4 previews are counted, and the validation is one configuration (l.150).
- **The A family:** the lead preregistered "A-only panels by 801, follow-up 5" at 11:49Z (LAB-LOG l.468–471). It specifies six
  **independent** draws, seeds 801301–801306, the donor-perts essential cut (DESIGN amendment), and the same bar.
- PROPOSAL-2, frozen at 11:53Z, proposes a second A-only read:
  - sequential and disjoint, seeds 814201–814206 (l.143);
  - on a pool that differs by the 7 library-only genes (l.48–53);
  - without citing follow-up 5.
- "The rule fix is not a try" (l.150) holds for one A read, not two. Two A-family reads with different pools and draw modes
  are two tries, and the better one would be read with best-of-2 bias.

### 6. Is the expected gain worth its cost? **The dollars, yes (≤ $5). The decision value has not been shown.**
- The over-read under attack is far above any class noise: key 3.08×, direction ≈ 6.5× (PHASEB l.26–27). For part A that
  makes the instrument fix worth doing, and 801 follow-up 5 is doing it.
- **For B**, the proposal itself makes v2-B factors non-decisive until P-B supports on 10-22 (l.166–168). From 10-22 14:00Z:
  - the lab's queue builds **D/E/F-panel folds**, Arc's real panel on the public lines, "the most October-like instrument …
    the 6a decision would rest on it" (LAB-LOG l.32–33);
  - the 11:14Z note expects those folds not to carry the over-read (l.427).
- So B's factors have no window in which they decide anything. The proposal does not name one: neither a pre-10-22 decision
  (where B is report-only) nor a post-10-22 one that D/E/F-panel folds cannot serve.
- The ordering also wastes money. P-B costs ~$0.02 and needs no folds (l.123). The ≤ $5 fold build runs before it, so it is
  spent even in the branch where P-B withdraws B.

### JUDGE.md required contents
| item | status |
|---|---|
| 1 Mechanism | present: A = Arc's stated criterion (l.20–29); B = KD / coverage propensity (l.63–69) |
| 2 Class | "instrument" (l.4), accepted as in VERDICT-1 for a lane that prices nothing |
| 3 Prediction / bins / falsifier | A and P-B present (l.111–129). **B over A missing** (check 2) |
| 4 Novelty | present (l.19–79) |
| 5 Cost and tries | cost present (l.149). Tries collide with 801 follow-up 5 (check 5) |
| 6 Day-of path | P-B computed from D/E/F's target list and public banks only (l.123); B's weights are A/B/C-fitted and labelled so |

## What is missing (part B; text only, no compute)
1. **A falsifier for B over A, with bins, before any number.**
   - Name the comparator. It must be an A-only read with the **same draw mode and pool**, so that the weights are the only
     difference. In practice that means 801 follow-up 5's independent draws (re-pooled on `esslib.txt` if the lead adopts it),
     not 814's sequential A-only.
   - Give a margin on ABCp's mean percentile, and on the key-class percentile, that clears best-of-12 selection and the
     21-cell panel noise.
   - State what happens below it: B closed for this lane, and A-only factors stand.
2. **The decision B serves.** Name one decision before 11-04 12:00Z that v2-B factors would change and that the D/E/F-panel
   folds (LAB-LOG l.32–33) cannot serve, or withdraw the fold build.
   - If such a decision exists, order the work **P-B first** (10-22, ~$0.02) and build B's folds only on P-B support.
3. **Reconcile with 801 follow-up 5.**
   - Withdraw 814's A-only variant (l.138–154), because 801 owns the A read.
   - Send the `esslib.txt` finding (the library union, 43 targets beyond the donor perts lists, 7 of them in the null pool)
     to the lead as an input to DESIGN §3.1 **before follow-up 5 draws**. Changing follow-up 5's pool is the lead's
     decision, made before any of its numbers exist.
4. **Pin P-B's scoring.**
   - Record prop10's ten frozen decile edges as numbers, recomputed from `features814.csv.gz` `bdee8c2f…` over the 4,355
     fit rows `e3b6b04c…`.
   - State that a D/E/F target with missing features counts as outside deciles 9–10, with the denominator held at 300 (or
     give another rule, fixed now).
5. **A real known answer for the draw code.** The six new panel shas are not pre-recorded (l.105–108), and two scripts agreeing
   with each other is not a known answer.
   - Before drawing on the 4,348 pool, the draw script must reproduce the pre-recorded V2p shas on the superseded 4,355 pool
     (`sim2/sim814.json`, `546e658f…` … `5939c5bc…`).
   - Record the six new shas in a dated note before any fold is built.

## Note to the lab lead (not a condition)
- **Part A is done as far as this judge is concerned.** It is FINDINGS §8.1g's rule, already in DESIGN §3.1, and its read is
  801 follow-up 5's. Only one open point matters: §1.2's evidence that the screens' full library (`esslib.txt`, 2,394,
  `24fffbfe…`) is the right cut, not the ≥ 20-cell donor banks (2,351). The judge's set lookup confirms the numbers. The
  difference inside the pool is 7 genes (CDK7, CEBPG, CMTR2, POLR3A, POT1, RBM22, ZNF830). Settle it before follow-up 5
  draws, because afterwards it becomes a second try.
- **Recommendation on B:** no third text round now. P-B is frozen in this proposal with bins (l.123–129) and costs ~$0.02 on
  10-22. Read it then, beside P-A, as the lab's record of whether Arc's selection follows prop10's knockdown and coverage
  axes. That is a scientific result in its own right. B's folds come back only as a new proposal, and only if P-B supports
  and a decision is named that D/E/F-panel folds cannot serve.
- The sealed LZ folds (L7, L8) were drawn under arc800. VERDICT-1's note stands: count their library targets (`esslib.txt`,
  not only the donor lists) at the one access, before the price is read.

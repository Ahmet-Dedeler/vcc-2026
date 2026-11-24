# DRAFT §3(m) for Manuscript v17: what October's transfer science adds (science lane, 700–708)

Status: draft by the science session [d80a10], 2026-10-15 (F25 added 10-16), on a copy (not merged into PAPER.md). Claims are numbered
provisionally C-99 … C-108; their provenance rows are at the end. Every number below has a file behind it. The 702 board
read is final; 705 / 707 / 708 are preregistered and unread at the time of writing, and are marked [[PENDING]].

---

### 3(m) What an external pretrained model, a held-out over-read and a panel's draw teach about zero-shot transfer

**F21: A pretrained perturbation model is a better direction than a weak baseline and a worse one than a donor consensus.
It is not a free upgrade. Its one held-out edge is ranking which genes respond, not their direction or size.
Provenance: candidate NEW (an external, preregistered evaluation of Arc's PIE on data it never saw, ended by a board test).
Confirmation: held-out reads on K562 GWPS (an unseen context) and h1-test (a seen one); one board entry (702, HARM).**

*Setup (700; C-99).* PIE (Verma et al. 2026) outputs, per (context, perturbation, gene), a DE probability, a log2 fold change
of mean CP10K and a δ-expression. Its fold-change head is supervised only on tested genes (≥ 5 CPM).
- We ran its cross-dataset checkpoint with its own evidence cache (built on a cloud worker; ~68 GB RSS) and context texts
  for every CCLE line. The texts were embedded with the model's own text encoder, reproduced to cos 0.99997.
- An audit fixed what PIE had seen. K562 GWPS is an unseen context; h1-test's context (H1) is seen but its responses are
  not; Orion and Jiang lines are contaminated.
- Below 5 CPM the fold-change head is unconstrained: |lfc| 8–26, A/B/C rows cos 0.02 to anything. Those genes were zeroed
  before any use.

*What it predicts better and worse than the donor stack (700; C-100).*
- **Signs, worse.** On h1-test's called non-DE genes (n 195,858) the sign accuracy is: PIE 0.5115; the 4-line consensus
  CLEAN4 0.5171; our emitted signs 0.5154. A cross-fitted logistic over every source, PIE included, ties CLEAN4 (0.5170).
  Only the line's own true mean response (an oracle) reaches 0.552.
- **Target strength, worse.** Σ p_de predicts h1-test's number of significant genes at Spearman 0.42, against 0.55–0.60 for
  the donors' noise-normalised strength key (521). Per-target sizing stays closed.
- **Which genes respond, better.** AUPRC of the truth's significant genes is 0.175 for PIE vs 0.148 for our design key
  (random 0.074); 0.151 vs 0.137 on unseen targets. This is measured only on a seen context.

*The board test (702; C-101).* PIE's per-target direction, blended into our bulk at weight 0.25 on PIE-seen targets
(norm-preserving), raised the held-out energy-weighted cosine: +3.7 % on K562 (all targets), +8.8 % on h1-test. It was
priced at +0.0008 … +0.0017. The board read −0.0025: pds −0.0042 raw and mse +0.0049 raw, both worse, DE flat. A
competitor whose base was weaker reported a +0.047 gain from the same model. So a pretrained direction helps in proportion
to how far its user's base sits below it, and a three-line whitened donor composition already sits above it on these lines.

**F22: Held-out instruments over-read a new source's direction for four separable reasons, and each gives a pricing rule.
Provenance: NEW (a post-mortem with bootstrap intervals and an axis decomposition). Confirmation: the board's sign
reversal on 702; the rules have not yet been re-tested on a second family.**

*Decomposition (702 RESULTS; C-102).*
1. **The unseen-context read was weak, and it was priced at its point estimate.** K562's gain was +3.7 % [+0.45 %, +6.8 %]
   (4,000-draw target bootstrap). At its lower bound the same translation gives primary ≈ −0.0003 (NULL), not +0.0008.
2. **The seen context supplied the upper end.** h1-test's gain was 6.6× K562's on the same construction, and positive in
   every expression bin, while K562's low-expression bin was null. It predicted the pds cost correctly (−0.0014 vs board
   −0.0015) but not the gain.
3. **The instrument saw only part of the change.** 41 % of the change's energy on A/B/C lies on genes outside K562's truth
   axis. PIE also moved 5–20 CPM genes at twice the base's energy share.
4. **The translation had no point of its own family to retrodict.**

The rules (adopted in the project's record): price from the target-bootstrap lower end; let a seen-context surface bound a
cost but never supply a gain; report the share of a change outside the instrument's gene axis (above ~25 %, the price is
partial); give a family without a board point a range that includes zero. 707's preregistration applied all four before its
read [[PENDING: 707 board read]].

**F23: In the strong tier, the remaining gap is direction, not amplitude. A leader whose advantage is "better amplitude on
strong targets" cannot be matched by re-weighting public features.
Provenance: NEW (13 held-out folds on exact fold surfaces with their known answers). Confirmation: analytic; the
amplitude rider's board read is a safety read only (705) [[PENDING]].**

*Budget (706; C-103).* The scorer's mse is Σ N_p / Σ D_p, so energetic targets carry it.
- On 13 held-out folds (522's surfaces: 4 on A/B/C's panel, 9 on October-like draws) we decomposed the aggregate cosine X
  of the production bulk:
  - the weakest third by donor strength is over-weighted on every fold (optimal relative scale 0.18–0.85);
  - the strongest third is slightly under-weighted on 10 of 13.
- A perfect three-tier amplitude rule is worth +0.2 … +2.3 % X on A/B/C-like folds (≈ +0.0001 … +0.0016 overall). A
  per-target oracle is worth +20 … +84 %, but it needs knowledge that no public feature carries: truth energy follows the
  key at Spearman 0.0–0.2 (522).
- So amplitude is near its predictable optimum, and the strong targets' error is direction.
- The board's leader posts mse 0.79 scaled at a field-level pds. Its gain sits on energetic targets, which a re-weighting
  of the same public evidence cannot reach.

*A direction lead (706; C-104).*
- Blending the source donors' raw (unwhitened) unit-mean into the top donor-strength tertile (β 0.25) raised X on 9/9
  October-like folds (+0.003 … +0.004 in board units), and survived centring and the removal of essential-like targets.
- On A/B/C's panel it gained only +0.0010 (cancer folds) and lost on the stem-cell fold (−0.0019).
- About fifty configurations were tried, and the best-of-N-corrected A/B/C prediction is ≈ 0. The board read decides
  [[PENDING: 707]].

**F24: Panel-dependent switches should be decided on the panel that will be scored, from held-out public lines, after its
release.
Provenance: NEW as a design rule (built on 438 §4.4's rules-vs-panel decomposition and 706). Confirmation: a rehearsal that
must reproduce 438's A/B/C-panel choice exactly [[PENDING: 708 rehearsal]]; the day-of read itself is unscored.**

*Why (C-105).*
- Three switches price differently on random October-like draws than on the one panel Arc built:
  - the key's cap (522: the draws prefer uncapped, A/B/C does not);
  - CD4's dose (438: ×1 on the draws, ×1.5 on A/B/C);
  - 706's direction swap.
- 438 §4.4 showed that A/B/C's panel under the day-of composition rules prices within 0.00014 of A/B/C's own. The
  difference is the panel's draw, not the pipeline.
- Arc selects its panels the same way each time, so the test panel's own folds are the right instrument. They need no
  truth from the test lines.

*How (708; C-106).*
- After the decided recipe is submitted, three held-out fold operators are fitted on the test panel. 438's re-pricer
  re-decides the stack weights under its fixed rule: move only by ≥ +0.0010, on 3/3 folds, and positive after
  leave-top-3-out. The same folds price the late switches.
- A move ships only as a later final entry.
- On A/B/C's own panel the rule moves to H1 ×2.5 / CD4 ×1.5 (+0.0022 over the October default); on the draws it stays.
- Which regime the test panel is in is unknown until 10-22 [[PENDING]].
- Every 708 action must also pass on the targets the held line expresses (AMENDMENT 4, below). This can only block a move.

**F25 (provisional, 10-16): A fold scores only what its held-out line can show, and a panel-dependence result must name
its instrument.
Provenance: NEW (an expression-floor test on 438's fold surfaces with score708's exact arithmetic; a tail diagnostic; an
operator audit). Confirmation: the floor result is measured on three panels; the operator question is open [[PENDING: the
lab's refit of A/B/C's fold operators]].**

*The floor (708 AMENDMENT 4; C-107).*
- In a PLOLO fold the held-out public line is the context being predicted. A target it barely expresses cannot be knocked
  down there.
- Below Arc's floor (the panel's 5th-percentile control CPM, 13.41 on A/B/C), no target reaches an own-gene log-effect of
  −0.5. Above it, 24–76 % do.
- Such targets make up 13–24 % of the Orion folds on Arc-like panels and about half on one older October-like draw. They
  are also less predictable from the donors.
- Restricting the folds to expressed targets changes the direction swap's price by −35 % to +16 % and never flips a
  decision. The floor is a correctness condition for the instrument, not the source of the panel effect.

*The tail, and an instrument confound (712; C-108).*
- The LODO key's tail differs sharply: the fold cap is 5.2–7.0 on A/B/C against 7.8–19.9 on every draw.
- Raw target strength in the donors has the same tail on Arc's 2026 panel, Arc's 2025 H1 panel and all draws (P96.6 /
  median 1.9–2.4). It does not track the cap (Spearman +0.08).
- A/B/C's fold operators were fitted with the full donor basis. The draws' operators, production's and the final day's
  drop the 2,021 K562-essential perturbations and cap each donor at 12k rows. So "A/B/C's panel is different" is partly
  untested: the same contrast also compares two operator recipes. Resolved the same day: rebuilt with the day recipe
  (LABC), A/B/C's panel still ranks last on its held line in 14 of 21 cells. The tail is the panel's, on the source side,
  and Arc-drawn panels should not carry the fresh-panel over-read (key ×3.08, direction ≥ ×2, shape ×1.2; the lab's 800).

---

## Provenance rows (to append to CLAIMS.md)
| claim | statement | source | verdict |
|---|---|---|---|
| C-99 | PIE setup, audit (unseen K562 GWPS, seen H1), < 5 CPM lfc unconstrained | 700 RESULTS (setup, audit700.json, cloud run 20261014T225548Z) | measured |
| C-100 | signs 0.5115 vs CLEAN4 0.5171; Σp_de Spearman 0.42 vs 0.55–0.60; AUPRC 0.175 vs 0.148 | 700 work/signsplit700.json, a0sources700, eval700_h1_v2.json | measured (h1-test = seen context) |
| C-101 | blend w 0.25: held-out +3.7 % / +8.8 %; board −0.0025464 (pds −0.0042, mse +0.0049 raw) | 702 RESULTS, read702.py, entry U7BlIZ53wTlLbaZTTET2 | board-confirmed HARM |
| C-102 | over-read decomposition: K562 CI [+0.45, +6.8] %, H1 6.6×, 41 % off-axis, K562 low-CPM null | 702 overread702.py → work/overread702.json | measured |
| C-103 | amplitude budget: 3-tier oracle +0.2 … +2.3 % (A/B/C-like), per-target oracle +20 … +84 % | 706 budget706.json | measured (13 folds) |
| C-104 | direction swap: 9/9 October-like +0.003 … +0.004; A/B/C +0.0010, kolf2 −0.0019; ~50 configs | 706 dirswap706*.json, h1swap706.json | measured; board read pending (707) |
| C-105 | panel dependence of key cap / CD4 dose / swap; abcdef ≈ abc within 0.00014 | 522 RESULTS, 438 §4.4–4.5, 706 | measured |
| C-106 | day-of re-pricing design; A/B/C-panel choice H1 ×2.5 / CD4 ×1.5 (+0.00217) | 708 PREREGISTRATION, work/expected_abcdef.txt | design + re-print of stored outputs |
| C-107 | held-line floor: 0 % of below-floor targets reach own-gene lfc < −0.5 vs 24–76 % above; 13–24 % of Arc-like Orion folds; swap −35 % … +16 %, no verdict flips | 708 floor708.json, floorprice708.json, RESULTS (c450c150) | measured (3 panels) |
| C-108 | fold cap 5.2–7.0 (A/B/C) vs 7.8–19.9 (draws); strength tails equal incl. Arc H1 2025 (Spearman +0.08 to the cap); A/B/C fold operators unpruned vs pruned elsewhere | 712 RESULTS, tail712.json; 062 / 233 / 238 operator jsons | measured; operator question resolved: the panel (LAB-PANEL-REPORT, 569978b9) |

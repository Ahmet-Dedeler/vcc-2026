"""533: Draft v15 text blocks, part D (§4.1, §4.4, §4.5, §4.6, the new §4.7, §5.6, §5.8, the figure index)."""

S41_F18_ADD = """
  - **[[PENDING 09-30]]** The direct test of the laboratory component on the challenge's own lines is the producer's H1
    programme as every context's block prior (458d) and its gate (458g), §3(j)."""

S41_ADD = """

**Later candidates, not yet audited (2026-09-30), stated with their caveats.**
- **The uncalled-gene hook's bar (C-76, C-77).** A same-lab programme can pay on uncalled genes where it cannot pay as a
  block prior, because there it has to beat the design's own signs (0.53 on H1) rather than a cross-line consensus.
  - Caveats: H1 only, the programme's own line, so an upper reference. The estimator was chosen on the same instruments.
    The attenuation is an emulation (random sign flips, not a real change of line). The cross-line analogue rests on two
    laboratories. The board read is pending (527).
- **Reach is sign-limited, and per-target call-set sizing does not transfer across lines (F20, C-86, C-87).**
  - Caveats: an analytic walk over one pipeline's order on H1, with exact known answers; the cancer folds' truths are
    150-cell pseudobulks. The global size's direction on the board is unread (523).
- **At its cliff, the profile-error member reads the aggregate cosine (C-84).**
  - Caveats: an identity under this pipeline's re-levelling and ballast construction. The rivals' implied cosines assume
    their mean-row and control terms equal ours."""

S44_ADD = """
- **No claim that public data cannot supply fold-change signs in general.** 516 closes the constructions it tested:
  target-specific public rows; their Stouffer, majority and stacked combinations over up to seven sources; and denoised and
  related-target-borrowed rows. It tested them on H1 and on three cross-study cancer folds, whose truths are 150-cell
  pseudobulks (two carry a technical component). 517 closes three target-independent sources at their constructions: the
  identity axis, library-size composition and dropout.
- **No claim about any competitor's method** beyond what its opened code does (525). The six-team `cellpred` cluster is
  private and its fold-change edge unexplained. §3(i)'s mapping of the leaders' nmae to sign accuracy is an upper bound on
  what they need, not a measurement of their signs (C-59).
- **Per-target call-set sizing is not claimed as a lever**, and its H1 predictability (r 0.55–0.60) is not claimed as a
  property of cancer lines (r 0.19–0.43; C-87). The global size's sign on the board is unknown until 523 reads.
- **H1 prices of DE sign sources are not board claims.** 517's hook (nmae −0.0183), 517's block prior (+0.0169
  overall-equivalent, one seed, estimator chosen after seeing reach) and 516's +0.016 are H1 numbers. Sign priors have
  transferred 0.07–0.67× of their H1 fid to the board, and their reach and nmae effects are not priceable on H1 (T9).
- **The pipeline refit's equality (C-82) is a statement about A/B/C relabelled**, not about a new panel. The test panel's
  composition (CD4 coverage, the bulk key's tails; 438, 522) is what the day-of rules price, on local instruments only."""

S45_ADD = """
- **T9 — local prices against the board, 2026-09-27 … 09-29.** The six reads since the H1-to-board calibration table
  (`docs/H1-BOARD-CALIBRATION.md`, 2026-09-26) was compiled; raw members, arm minus base, predictions as preregistered:

  | read | class | local price | board | verdict |
  |---|---|---|---|---|
  | 446bj2z, exact-rank pass | DE magnitude | nmae −0.00183 (H1, 3 seeds) | −0.002179 | transferred, 1.19× |
  | 446br2, repair floor −2 | DE delivery defect | nmae −0.002 … −0.008 | −0.001413 | partial, below the range |
  | 446bj2, LODO bulk key | bulk magnitude key | mse −0.004 … −0.012; pds +0.001 … +0.004 (cancer folds) | mse −0.014473; pds +0.000208 | mse past the range; pds below it |
  | 453dr, refit path | bulk composition (CD4 ×1) | overall −0.0024 (band −0.0045 … −0.0005) | −0.0036816 | inside the band |
  | 457dr, CD4 ×1.5 | bulk composition | overall +0.0017 (k 1) … +0.0022 (k 1.3) | +0.0031484 | transferred, 1.83× the k 1 price |
  | 446bd2, fifth line | bulk composition | mse −0.001 … −0.004; pds −0.0005 … +0.002 | mse −0.000471; pds +0.002148 | pds as priced, mse did not transfer |

  Pure magnitude changes keep transferring at or above their local prices (with 433a's 1.76× and 410's 1.6×). A delivery
  defect scales with its prevalence. A composition step adding a fifth line moved pds as priced and mse not at all. Sign
  sources remain what local surfaces cannot price: 446bja's programme carried 0.07× of its H1 fid (C-63).
- **T10 — within the record: the recipient's identity axis.** On H1 a line's own control identity axis ("genes higher in
  this line go down") signs its responses better than the consensus, and control halves show the gain is honest there (517).
  Across six public lines it holds for HEK293T and weakly HCT116, and reverses for K562, RPE1 and HepG2 (C-80). Resolved as
  a property of some lines, unusable where the line's kind is unknown.
- **T11 — 516's H1 gain and the board.** Most of it (about +0.017) is the unwhitened-carrier effect that inverted nmae on
  the board (417c: −0.0068 predicted on H1, +0.0089 raw on the board). This is not a contradiction inside 516: its
  cross-study cancer folds read the same ≈ 0 that the board implies."""

S46_OLD = ("The programme test (§3(j)) and the two-functional decomposition (§3(k)) are, as scored board entries on those three "
           "contexts.")
S46_NEW = ("The programme tests (§3(j)), the two-functional decomposition (§3(k)), the pipeline refit (§3(h)) and the second "
           "fold-change delivery (§3(i)) are, as scored board entries on those three contexts. §3(l) and the H1 prices in "
           "§3(i)–(j) are measured on H1 and public folds only.")

S47 = """

### 4.7 Open problems

What would move the frontier, stated only from the record.

1. **Signs on weak and uncalled responses.** On H1, oracle signs at this pipeline's magnitudes recover −0.254 raw nmae,
   2.6× the oracle-magnitude pot (−0.0985). 65 % of that sign pot sits on genes the pipeline does not call, where its sign
   accuracy is 0.54 (480). Two public-data leaders post board nmae 0.827–0.835 (2026-09-29, against 0.895 here); under
   §3(i)'s iid map that is a gate-wide sign accuracy of about 0.65–0.66, against 0.56 for this pipeline on H1. That is +0.06
   for the group of entries at 0.87–0.88 and +0.08 … +0.10 for those two teams (434f; 480's rate gives ~+0.08 for one of them). The board's
   leader, which very likely trains on proprietary many-line data, posted 0.804 on 09-29 (532's board read; about 0.69 on
   the same map, by interpolation). These are upper bounds on what the entries need (C-59). Each +0.01 of gate-wide sign accuracy is worth
   −0.006 … −0.008 raw nmae (480).
2. **Public information shown insufficient.**
   - Target-specific public rows and their combinations: +0.016 on H1 and ≈ 0 on three cancer folds; an oracle choice of
     source per target caps the H1 gain at +0.030 (516).
   - Other laboratories' programmes: among cancer lines a programme is mostly its laboratory's signature on non-essential
     responses (450h, 462, 466).
   - An independent screen of the same line: Jurkat, null on the board (446bja, C-63).
   - The recipient's own control cells: co-expression carries only the top two cell-state axes, not regulation (174), and the
     identity axis reverses across lines (517).
   - Many-line drug atlases: a basal-conditioned transfer operator exists across 50 lines, but it is an amplitude law and is
     absent from the CRISPRi donors (184).
   - Public emission designs (525) and technical features of the chemistry (517).
3. **What kind of data would supply it.** The one source that closed most of the gap was a screen of the same line (Pan's H1
   screen: gate signs 0.764 against 0.551 on covered targets; +0.093 through the scorer at 35 % coverage; H1 only). None
   exists for the final lines (520). The route still open is the producer's own screen of another line.
   **[[PENDING 09-30: if 458d carries, a screen's laboratory and chemistry are worth more than its line, and a benchmark
   producer's own screens are the most valuable public resource for its benchmark; if it reads null, a line's programme is
   available only from a screen of that line.]]** Commercial many-line atlases (a leader appears to train on one) are the
   other route, and they are not public.
4. **The magnitude pot.** At our directions, per-target magnitudes are worth +0.018 … +0.023 overall (518), more than the
   profile-error gap to two of the leading teams (+0.0095 … +0.0098 overall). Half of what remains after the current keys
   sits in 11–19 targets per line that no available feature predicts, and donor-side features, dependency included, do not
   supply it (514, 518). Being line-idiosyncratic, it would need a measurement in the line itself.
5. **The call-set regime.** The optimal global call-set size depends on the reference's distribution of response sizes,
   which the board does not publish; the preregistered dose pair (523) will measure it. For benchmark designers, fixing the
   call set relative to the reference's counts, or reporting reach by response-size band, would remove a calibration term
   from two members (§3(l)).
6. **The construct.** §3(h)'s ceiling still stands far above every entry: a gene-level cosine of 0.21–0.57, against the
   median per-target transfer cosine of ≈ 0.04 measured there (the current pipeline's aggregate cosine on H1 is 0.11; 518).
   Knockdown-validated constructs, and a second construct per target in any new data, remain the cheapest improvement to the
   references themselves."""

S56_ADD = (" Two lessons from the engineering track's known-answer rows (2026-09-27 … 09-29) generalise. First, a known-answer "
           "row must check the driver's exit status as well as the output hash: one run matched its expected hash while the "
           "driver had raised after writing it, and every later build would have died (479; found by a review from a "
           "different model, 499). Second, a gate must compare recorded values, not the presence of a record. Two reviews of "
           "the final-phase path by a different model found a post-build check that, once the final panel's download replaced "
           "the reference file, compared the new target list with itself and would have failed every correct final build (519), "
           "and a gate that passed artifacts missing the recipe's same-line record or carrying the wrong parameter values (526). "
           "Each defect was reproduced by a script that exits 0 while the defect is present, and fixed before any final-phase "
           "build.")

S58_ADD = (" Round 14's derived objects are `experiments/513-reach-walk-anatomy/work/{walk513,demote513}.json`, "
           "`experiments/516-de-sign-sources/work/{best516.json,fold516.log,power516.log,stack516.log}`, "
           "`experiments/517-flex-programme/work/{real517,grid517,hookgrid517,mzgrid517,labhook517,progest517,xline517,tech517,"
           "dropout517,callprop517}.json`, `experiments/518-bulk-mse-gap/work/{board518,folds518,h1_518,remain518}.json`, "
           "`experiments/521-strength-nreal/work/{pred521,purity521,global521}.log`, and the captured board entries in "
           "`data/board-entries/`; `experiments/533-science-lead/work/verify_board.py` recomputes every board delta this "
           "revision states.")

FIG_OLD_INTRO = ("All eight are regenerated from the saved tables by `make_figs.py`; full captions with units, n and interval "
                 "source are in `figs/MANIFEST.md`.")
FIG_NEW_INTRO = ("All thirteen are regenerated from the saved tables by `make_figs.py`; full captions with units, n and interval "
                 "source are in `figs/MANIFEST.md`.")
FIG10_OLD = "| `E450/work/labline450.json` |"
FIG10_NEW = "| `E450/work/labline450_v3.json` |"
FIG_ROWS_ADD = """
| **F-12** | `figs/F-12.{png,svg}` | a zero-shot pipeline re-fitted from the panel alone reproduces its tuned self; six single-channel entries read by the two channels | `data/board-entries/*.json` (seven entries) |
| **F-13** | `figs/F-13.{png,svg}` | the uncalled-gene hook's bar is lower than the block prior's, as the programme's information is attenuated (H1) | `E517/work/{grid517,hookgrid517,mzgrid517,real517}.json` |"""

JIANG_OLD = "four labs: five Weissman-lab lines, two Xaira lines, six Jiang 2025 lines, one Satija K562 screen in Flex-style probe chemistry,"
JIANG_NEW = "four labs: five Weissman-lab lines, two Xaira lines, six Jiang 2025 lines [JIANG25], one Satija K562 screen in Flex-style probe chemistry,"

LEADERS_OLD = "- The public board's leaders sit at 0.82–0.84, the 0.65 point of that curve (0.56 here)."
LEADERS_NEW = ("- The public board's leaders sat at 0.82–0.84 when this section was written (2026-09-24), the 0.65 point of that curve (0.56 "
               "here); on 09-29 the "
               "leader posted 0.804 (§4.7).")

ANCHORS_D = [
    (LEADERS_OLD, "replace", LEADERS_NEW),
    (JIANG_OLD, "replace", JIANG_NEW),
    ("  - The closest prior work names the component (Systema); what is candidate-new is where it lives.", "after",
     S41_F18_ADD),
    ("  - Why it matters: any \"universal response\" or programme prior estimated from a single screen needs control halves\n"
     "    to be judged.", "after", S41_ADD),
    ("- **Claim C6** (overexpression reverses knockdown in hESC) is VOID, not decided: the reserve was spent under a rule its "
     "own known-answer row voided; what is reported is a bounded null.", "after", S44_ADD),
    ("it says nothing against reliability-based filtering of training data. Naming the operation is the whole fix.",
     "after", S45_ADD),
    (S46_OLD, "replace", S46_NEW),
    ("this programme's own record contains two claims that had to be retracted for exactly that kind of extrapolation.",
     "after", S47),
    ("and a pass rule on a multi-line or multi-study set must be depth-conditioned before it is frozen (this is why C2a's "
     "miss was predictable).", "after", S56_ADD),
    ("with no expression file reopened.", "after", S58_ADD),
    (FIG_OLD_INTRO, "replace", FIG_NEW_INTRO),
    (FIG10_OLD, "replace", FIG10_NEW),
    ("| `E462/work/*.json`, `E458/work/denoise458d.json` |", "after", FIG_ROWS_ADD),
]

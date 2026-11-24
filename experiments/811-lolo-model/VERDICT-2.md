# VERDICT: 811 PROPOSAL-2 (M811, frozen copy) — **GO, with five binding conditions (G1–G5)**

Judge run 819. Written 2026-10-16 08:29–08:40Z by the clock. I judged the FROZEN copy
`experiments/819-judge/811-lolo-model-20261016T082515Z/PROPOSAL-frozen.md`, sha256
`8336f260d0cb7fc78e3eaf94f4dca2477b856c700e34c0405ac0a8029e3b50c7`. That is the sha the mandate names, and I re-hashed it.
The earlier verdict is `experiments/811-lolo-model/VERDICT-1.md` (REVISE). This is the one resubmission JUDGE.md allows.

I computed only lookups:
- sha256 of every file in the proposal's §0. All ten match: PROPOSAL.md `cc5b8a8f…`, RESULTS.md `f15f0298…`,
  PREREGISTRATION.md `88440507…`, table-base811.md `ffeee455…`, run-base811/SHA256SUMS.txt `05643d07…`, align811.json
  `28b88bd2…`, base811.py `c29f9f97…`, **m811.py `4818048e…`**, agg811.py `d34d1e26…`, banks/manifest.json `37e34335…`;
- reading `m811.py`, `base811.py` and `engine801.py`;
- sub-second reads of fields in `work/align/align811.json` (102 kB).

**Budget granted:** the proposal's two stages, ≲ $1 together (`PROPOSAL-frozen.md:163-165`), and ≤ $3 for the lane
including one re-run. Stage 2's ABCp half waits until 801's ABCp dumps are in its registry. At 08:29Z,
`experiments/801-lab-folds/work/dumps/abcp/` does not exist locally; 801's `lab801-abcp-a` run terminated at 08:28Z.

**Why GO.**
- Every item of VERDICT-1 (A–H) is fixed or answered; see the table below.
- The construction is not closed (check 1).
- The decision rule is written before the data (`:143-161`) and has a real falsifier for the October claim (F2, `:141`).
- The run costs about $1.
- Its outcome is useful either way. It is the first price of the many-line consensus on the D/E/F day's own bank list,
  read on day-recipe A/B/C operators (ABCp). A NULL closes the family for October with a mechanism.

**Why not a clean GO.** Five defects bind how the run is done and read. None changes the idea or the grid of 4 arms.
- **G1.** The frozen `m811.py` cannot produce the decision set. On ABCp folds it voids every fold and uses the wrong
  incumbent.
- **G2.** Falsifier F1 has no power against the rival it is meant to detect.
- **G3.** The selection rule can pick an arm that the proposal's own rule 6 bars from a dossier.
- **G4.** The BULK-ABC try count is one short in the ABCp-after-062 branch. One mechanism sentence is factually wrong.
- **G5.** Wear rows must be recorded.

The proposal itself pre-commits that a bug found after GO means a new sha, reported before any number is read
(`:21`). G1–G4 go into one frozen addendum before the first launch.

---

## VERDICT-1 items: did the resubmission fix them?

| item | fixed? | where |
|---|---|---|
| A. One bank list on every surface; KOLF2 later as a counted try; D/E/F identity rule | **Yes** | `DAY_GROUPS = B.GROUPS minus KOLF2` (`m811.py:41`) is used on every fold (`m811.py:112`) and on the day path (`m811.py:86-93`). "With KOLF2" is deferred as a separately counted try (`PROPOSAL-frozen.md:24-26`). The identity rule (`:70-76`; `m811.py:44-56`) holds out the matched line's group, so it is price-faithful: every priced fold held out its own line's group |
| B. Mechanism paragraph corrected | **Yes, with one false sentence** (G4) | Per-target raw cosine stated as contradicting the hypothesis (`:38-41`). Shared-energy ratios 3.8–7.9× (A/B/C) and 1.3–4.7× (L) match align811.json unrounded: A/B/C 3.79 / 7.90 / 4.50; L min 1.28, max 4.73. 706 cited (`:47-49`). P1 revised by an explicit bound (`:50-51`; asin√0.025 = 9.1°, λ·0.158 ≈ 0.04: correct). **Wrong:** `:45-46` says U's shared share on A/B/C without KOLF2 is "unmeasured". align811.json has it: `ABC/held/{K562,HCT116,HEK293T}/all-KOLF2/shared_U` = 0.0203 / 0.0219 / 0.0207, i.e. 1.8× / 4.0× / 2.0× the truth's |
| C. Drop factor W | **Yes** | `:29`; arms are C × λ only (`:64`, `m811.py:43`) |
| D. Predictions computable; r-gain baseline defined; pdsl beside | **Yes on computability; F1's reference is wrong** (G2) | P1 and F1 are on W-tune, where all 4 arms are priced (`:123-124, 139-140`). Every prediction gives literal pds and pdsl side by side (`:118-131`). The r gain is defined against `inc` (`:139`). That definition makes F1 powerless; see check 2 |
| E. Honest tries | **Yes, one inconsistency** (G4) | N per surface set with e(N) (`:110-116`). e(10) 1.547, e(7) 1.364 and e(9) 1.494 match Blom's formula (DESIGN §8 values). W-tune's N = 10 is conservative: c0 equals the baseline's all_lX on L folds, so it is not a new configuration there. `:92` promises "+1 read on the same truths, counted" if ABCp follows a 062 read, but the table's 9 (`:115`) does not include it |
| F. Decision rule before data, closure, W-conf miss, board projection | **Yes** | `:143-161`. Board projection: "a range that includes 0" (`:153-154`). Closure outcomes and the W-conf ±30 % miss: `:155-161` |
| G. DE-H1 hygiene | **Yes** | The whole H1 group is out, with bank coverage on h1-test / h1-trainval listed (`:100-104`). The limitation that 803 has no hook for external rows is stated, and it blocks any dossier when \|dK_rel\| > 5 % (`:106-108, 152`) |
| H. Provenance and instrument | **Yes** | §0 shas all re-hashed and matching (above). The arm module is named and hashed (`:17`). The surfaces are named: ABCp is the decision set and 062 is reported beside it (`:86-94`). The file was frozen before the judge read it |

---

## The JUDGE.md checks

### 1. Is the mechanism closed already? **No, not at this construction. Factor C sits next to several nulls, and the proposal now predicts a null for it.**
- **The equal-weight many-line consensus blended at λ** is not closed:
  - `docs/FINDINGS.md:885` ("Direction via a fifth donor line: sources saturate") covers donors inside the whitened stack.
  - `docs/FINDINGS.md:870` (343–382 "shared-across-targets priors … coherence is the price") covers shared priors. U(t) is
    per target.
  - 706/707 use the sources-only mean on the top tier (`experiments/706-strong-tier-budget/RESULTS.md:34-71`).
  - 811's own baseline is the first price of this construction: all vs src +0.0018 / +0.0050 / +0.0105 on BULK-W, 17/18
    folds (`experiments/811-lolo-model/RESULTS.md:114-120`).
  - The agreement weighting that VERDICT-1 flagged (016 / 142 / 271) has been dropped.
- **Across-target centring has a long null record.** The proposal cites only 706:
  - 706: "centring U changes nothing" (`706-strong-tier-budget/RESULTS.md:40-41`), cited at `:47`;
  - production's donor tables are already centred across each donor's perturbations, and adding their mean back does
    nothing (`experiments/153-de-bench/RESULTS.md:957-961`; `experiments/059-feng-ipsc/RESULTS.md:255`);
  - adding the mean response changes the centred cosine by +0.0002 (`experiments/004-donor-transfer-h1/RESULTS.md:25-26`);
  - rescaling the bulk's shared mean-across-targets component gives ≤ +0.0015 (393; `docs/FINDINGS.md:879`).
- **What differs here:** the mean is taken per LINE over the PANEL's carried targets, before a many-line average whose
  non-source lines carry the excess shared energy (`nonsrc/shared_U` ≈ 0.097 on A/B/C with KOLF2 in align811.json).
- The proposal does not claim a gain from C. It predicts \|ΔC\| ≤ 0.0010 / 0.0020 (`:123-124`) and uses C as a mechanism
  read (`:52-53`). That is legitimate, but only if the read can discriminate; check 2 and G2 cover this. Two configurations
  of a mechanism read do not re-name a closed construction. Not a KILL.

### 2. Falsifier and bins written before any number? **Yes for the October claim (F2 and the decision rule). F1 as written cannot fire under its own rival (G2).**
- Bins per arm and surface: `:118-131`. Expected bin on BULK-ABC is NULL or WEAK, "no dossier from this run" (`:131-132`).
  F2: `:141`. The decision rule with closure: `:143-161`.
- **F1 is mis-referenced.** F1 says the many-line gain is a shared component if "c1's r gain over inc < 30 % of c0's r gain
  over inc" (`:139-140`). But c0's gain over inc is mostly the sources' blend, not the many-line increment:
  - at λ 0.10: all vs inc r +5.86 % (`811/RESULTS.md:103`), src vs inc +3.89 % (`:110`), all vs src +1.92 % (`:114`);
  - at λ 0.25: +14.02 % (`:106`), +8.67 % (`:112`), +5.04 % (`:116`).
- Suppose centring removed the ENTIRE many-line increment and nothing else, which is the rival reading `:42-43, 52` names.
  Then c1 would keep about 3.89/5.86 = 66 % (λ 0.10) or 8.67/14.02 = 62 % (λ 0.25) of c0's gain over inc. That is far
  above F1's 30 % line, so F1 stays silent.
- For F1 to fire, centring would also have to destroy most of the sources' own gain. 706 already showed it does not
  (`706/RESULTS.md:40-41`).
- So F1 as written cannot fire under the rival it is meant to detect. The class label "transfer if F1 does not fire" (`:64`)
  would be granted by default. The fix (G2) uses arrays that already exist: base811's `src_l10` / `src_l25` on every
  W-tune fold (`work/run-base811/SHA256SUMS.txt`, 66 `src_l*` rows).
- P1's "C on keeps ≥ 70 % of C off's r gain over inc" (`:124`) is the sensitive threshold here: 62–66 % would breach it.
  It should be reported that way.

### 3. Measurable on the checker's surfaces? **Yes. One code defect makes the decision surface unreachable (G1).**
- mse is read through r, with literal pds and pdsl side by side (PHASEB `:8-11, 42-53`).
- The DE members are pinned by ballast K, and the \|dK_rel\| > 5 % trigger is stated (`:105`; DESIGN §4,
  `experiments/800-checker/DESIGN.md:70`).
- The PARTIAL label is in the rule (`:151`; DESIGN §7 `:91-92`).
- The DE-H1 limitation is honest: no 803 hook, so no dossier (`:106-108`).
- **ABCp path defect (decisive for the decision set):**
  - `engine801.AbcFold` keeps `self.panel = "ABC"` for ABCp folds. Only `fold_id` says `ABCp/<held>`
    (`experiments/801-lab-folds/engine801.py:178-179`).
  - In `m811.py:132-134`, `abc = F.panel == "ABC"` is therefore True on ABCp. So `known` is True, and `inc` is built as
    `cap6+P2g(abc)`. §3.2 says the incumbent there is "the October incumbent pct+P2g" (`PROPOSAL-frozen.md:89-90`).
  - Worse, `inc` and `707` are then compared with 804 references at `work/ref804/ABCp__<held>/` (`m811.py:159, 168-169`).
    Those do not exist: `work/ref804/` holds only `ABC__*` and `L*__*`.
  - `B.compare` raises on the missing file. The except branch marks it FAILED (`m811.py:176-178`), and the assert voids the
    fold (`m811.py:187-188`).
  - Every ABCp fold would exit non-zero. The proposal's rule 1 makes ABCp mandatory and says "062 alone never suffices"
    (`:145-146`). So, as frozen, the run can never reach a dossier.
  - Cosmetic: the ABCp `ka811` note is added after the meta file is written (`m811.py:184-186`), so it is never saved.
- **Selection vs rule 6:**
  - Selection takes the highest W-tune bound (`:96-98`), and the proposal expects a λ 0.25 winner (`:125`).
  - The baseline at λ 0.25 moved K by up to +5.9 % (`811/RESULTS.md:136`).
  - Rule 6 with no 803 hook means "no dossier" (`:152`).
  - So the preregistered selection will probably pick an arm the decision rule bars before BULK-ABC is read. The family
    would then "close" on the selection's dose bias, not on the mechanism (G3).

### 4. Does it use seen data as evidence of gain? **No. Not a KILL. Wear is charged.**
- U never holds the held line's group (`m811.py:112`; `base811.py:50, 160`).
- KOLF2 is never in U, so ABC/kolf2 is a clean extra held line (`:94`).
- H1 is out on DE-H1 (`:100-101`).
- No A/B/C board read and no sealed reserve (L7/L8, LZ, R6/R7) is used.
- Nadig Jurkat sits in the bank list, and A is Jurkat. The BULK-ABC truths are the held public lines (K562/HCT116/HEK293T/
  KOLF2), not A/B/C's contexts, so nothing seen enters a gain read. On the D/E/F day the identity rule covers a matched line.
- Wear:
  - the grid (λ 0.10/0.25, C) was formed after reading all 22 baseline folds and align811's A/B/C truths
    (`811/RESULTS.md:166-171`);
  - the proposal charges this: W-conf N = 7, BULK-ABC N = 9 including align811 (`:113-116`);
  - if the W-tune winner is a c0 arm, its W-conf read is the baseline's already-priced all_lX (bit-identical on L, `:59-60`),
    so it adds no new information and must be reported as such.

### 5. Is the try budget honest? **Yes. N ≤ 10 on every surface set, with nested selection (L1–L3 → L4–L6, ABC). One branch is short by 1 (G4).**
- W-tune N = 10, conservatively including c0 = all_lX (`:113`). W-conf N = 7 (`:114`). BULK-ABC N = 9 (`:115`).
- In the branch where ABCp is read after a 062 read, `:92` adds "+1 read on the same truths, counted". N becomes 10, and
  e(10) = 1.547 must be applied there.

### 6. Is the expected gain worth its cost? **Yes for $1.**
- The predicted BULK-ABC gain is +0.0003 … +0.0015 at λ 0.10 and −0.002 … +0.001 at λ 0.25 (`:130-132`). The upper half
  is above bulk direction's ~0.0003 readability (JUDGE.md). The proposal expects no dossier.
- The value of the run:
  - the first day-recipe (ABCp) price of the consensus on the day bank list;
  - a discriminating C mechanism read once G2 is applied;
  - an October closure with its mechanism if NULL.
- Calibration: bulk-external has one board point, 702, negative (DESIGN §9; LAB-LOG Known-answer report part 1: ρ 0.93 on
  BULK-ABC, 0.19 on BULK-W). The proposal claims no board projection beyond a range that includes 0 (`:153-154`): correct.

### Day-of path
Acceptable. Rows come from public banks on disk (`DAY_GROUPS`). Centring runs over D/E/F's own 300 targets. λ is fixed on
L1–L3. Nothing is fitted to A/B/C (`:66-68`; `m811.py:86-93`). The identity rule is price-faithful (`:70-76`).
Class: transfer, provided restated F1 (G2) does not fire and the held-out gates (`:145-149`) pass; otherwise it is scorer
mechanics or closed, as `:64, 158` say.

---

## Binding conditions (one frozen addendum, `experiments/811-lolo-model/PROPOSAL-2-ADDENDUM.md` with its sha256, written before the first M811 launch; the lab lead checks it)

**G1. Fix the ABCp path in `m811.py`; new sha before any number.**
- Detect ABCp by `F.fold_id.startswith("ABCp/")` (or `getattr(F, "variant", None) == "ABCp"`), never by `F.panel`.
- On ABCp, build `inc` as §3.2 says (`pct+P2g`). If `pct+P2g` cannot be built on an ABCp fold, say so in the addendum and
  name the incumbent used, with the reason, before any number.
- Skip the 804-reference KA there, so that only the engine's own fold KAs stand (`:91`). Keep "a missing reference = a
  failed KA" on `ABC/*` and `L*/*` folds.
- Write the ABCp note before saving the meta file.
- Record the new sha256. The KA `c0_lXX == base811 all_lXX` on L folds (`:82`) must still pass at \|d\| ≤ 1e-12 on 804's
  instance type and thread count.

**G2. Restate F1 on the many-line increment.**
- Copy base811's `src_l10` / `src_l25` arrays (sha-checked against `work/run-base811/SHA256SUMS.txt`) into the W-tune run
  dir. Aggregate `c0_lX:src_lX` and `c1_lX:src_lX`; inc is identical by KA, so the pairing is valid.
- **F1 fires** if, at either λ, c1's r gain over src_lX is < 30 % of c0's r gain over src_lX.
- Report the inc-based ratio beside it, labelled as P1's read: 62–66 % is what full removal of the increment would give.
- Caveat to state: centring also touches the source lines inside U. 706 found that changes nothing on the sources-only mean.
- The class label follows restated F1.

**G3. Selection eligibility.**
- Before selection, an arm whose W-tune max over folds of \|dK_rel\| exceeds 5 % is ineligible as the winner, since rule 6
  with no 803 hook bars it from a dossier.
- It is still reported on W-tune.
- If no arm is eligible, report "no DE-eligible arm". Price only the λ 0.10 C-pair on BULK-ABC as the mechanism read,
  within the same N.
- A closure written in RESULTS names the λ and C that were actually read on the decision set. It closes no other
  configuration.

**G4. Counts and facts.**
- BULK-ABC N = 10 (e 1.547) in the branch where ABCp is read after a 062 read (`:92`).
- Correct `:45-46` in the addendum: A/B/C's U without KOLF2 has shared share 0.020–0.022, which is 1.8–4.0× the truth, in
  align811.json `ABC/held/*/all-KOLF2/shared_U`. So C's expected effect on A/B/C is the same small size as on L. This
  supports P1.

**G5. Wear.**
- Append this run's would-be wear rows to `work/wear811.jsonl` (family `811-consensus`, N per surface set as above) for the
  lab lead to add to 800's ledger.
- `check800.py aggregate` stays unused, as in the baseline (`811/RESULTS.md:6-8`).

Then stop. Leave no instance running (`labcloud.sh status`).

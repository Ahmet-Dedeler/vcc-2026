# VERDICT: 810 PROPOSAL-2 (frozen copy) — **GO, with four binding conditions (G1–G4)**

Judge run 819, written 2026-10-16 08:16–08:25Z by the clock. I judged the FROZEN copy
`experiments/819-judge/810-strong-direction-20261016T081233Z/PROPOSAL-frozen.md`, sha256
`3153d6c8b37e2eef900e187694e15a15571f2dfcc914a3ba5591a914f12e51de`. The live `experiments/810-strong-direction/PROPOSAL-2.md`
has the same sha256 at 08:14Z, and `experiments/810-strong-direction/work/SHA256SUMS.txt:1` records it. I computed only
lookups:
- sha256 of the cited files;
- reading the code;
- two sub-second reads of stored dumps, checking that the panel-gene columns of `fin` are zero (see check 3).

**Budget granted:** the proposal's single EC2 run (22 folds × 5 arms, r7i.8xlarge, ~30 min cap, ≈ $1.5; ≤ $5 with one
re-run; `PROPOSAL-frozen.md:138-142`). The ABCp re-read in G1 is a second launch of ≤ $5, and is bought only if the 062
run reaches DOSSIER or the HOLD-707 → DOSSIER branch.

**Why GO.** Every item of VERDICT-1 (A–E) is fixed; see the table below. The construction is not closed. The run is
cheap, and its STOP comparator is worth having even if the candidate dies: it tells the lab whether the bulk checker
sees the γ knife-edge that closed global un-whitening (212 / 097). That holds whatever sgu1 reads.

**Why not a clean GO.** Four issues bind how the result is read:
- The deciding surface whitens the base with a different operator from the one the candidate would act on in production
  (G1). This matters most for a mechanism that is "undo the whitener".
- The HOLD-707 branch can inflate a near-zero checker read through a ρ < 1 taken from one uncalibrated board point (G2).
- The 804-array known answer silently skips when its file is absent (G3).
- The closure list omits 063 (G4).

None of these needs a change to the pricing code (`price810.py` da001560… stays as frozen).

---

## VERDICT-1 items: did the resubmission fix them?

| item | fixed? | where |
|---|---|---|
| A. Same-run closed global comparator, its expected answer, STOP rule, falsifier 2 rewritten, ship condition vs sgu1all | **Yes** | arm `sgu1all` = the same swap on every target with a donor row: `PROPOSAL-frozen.md:53`, `price810.py:93` (`sgu_rows(P, fin, have, 1.0, npn)`). Expected answer from 212 / 097: `:69-83`. STOP: `:79-81` and code `decide810.py:112-114`. C2 against sgu1all in checker units: `:101-103`, `decide810.py:121-122`. C3 (sgu1 − sgu1all ≤ 0 closes): `:104`, `decide810.py:123-124` |
| B. Mechanism numbers corrected; gate justified only by pds exposure | **Yes** | `:9-31`. The all-tier gains +0.0055 / +0.0074 / +0.0203 (W) and +0.0066 / +0.0090 / +0.0092 (ABC) match `810/RESULTS.md:23-25` and 99-100. The pre-only-swap numbers are re-attributed (`:21-23`). The claim is stated as an energy-weighted vs count-weighted exchange (`:24-32`) |
| C. Honest tries; r legs labelled not blind | **Yes** | One priced arm, N = 6, e(6) = 1.28 sd (`:58-61`; `stats800.e_max_normal(6)` = Blom ppf(0.9) = 1.2816; `decide810.py:32, 89-90`). β 0.5 and P2 are withdrawn (`:60`). "Not blind" is stated on the mse legs (`:76-77, 86-87`) |
| D. Decision rule before data; literal pds decides; pdsl beside; 707-wait outcome | **Yes** | `:108-136`, implemented as `decide810.py:108-136`. Literal pds decides, and pdsl > 0 is required for DOSSIER (`:110-111`, `decide810.py:131`). HOLD-707 resolution: `:128-136` (see G2) |
| E. Provenance; frozen file | **Yes** | `810/RESULTS.md:192-215`. I re-hashed price810.py `da001560…`, decide810.py `724e8059…`, decomp810.py `afe28579…`, analyse810.py `0b510e66…`, foldlevel810.py `65b4d723…`, job810p.sh `ea9e743f…`, jobs-p.txt `ee4fe683…`, sync810p.txt `c09135b5…`, work/foldlevel810.json `4c8d7f90…`, work/table810.md `0305e1ba…`, work/s3parts810.tsv `211b8689…` (S3 prefix), work/analyse810.json `1d506cb2…`: all match. The proposal is frozen, and the judge read a copy |

---

## The JUDGE.md checks

### 1. Is the mechanism closed already? **Not at this construction. Its global form is closed, and the run carries that form as a known answer.**
- **Closed, the global form:**
  - `docs/FINDINGS.md:577-578`: "the operator geometry is closed on all three axes (193, 212, 226)".
  - `docs/FINDINGS.md:785-787`: "residual de-whitening for mse … Operator geometry: γ/λ (193), exponent (212)".
  - `experiments/212-gamma-cosine/RESULTS.md:26-30`: d mse/d cos +0.3655 vs d pds/d cos −0.3517, cancelling to 3.8 %.
  - `:46-48`: "do not reopen it again".
  - `experiments/097-direction/RESULTS.md:1272-1277`: gain/spend 0.77 / 0.37.
  - `experiments/170-gamma-mse/RESULTS.md:1-17`: V at γ 0.65, through the ρ-cliff amplitude mechanism.
  - `experiments/063-pds-mse-pareto/RESULTS.md:7-12, 70`: "residual de-whitening is closed for this production geometry".
    That is a global β blend of the raw consensus, with pds monotone down at every β on all 3 folds.
- **Why this construction is not those:**
  - It swaps only the base part, on a fixed top-z_src tertile. It keeps cd4 / h1 / part and |fin_i| (`price810.py:53-64`).
  - So the key, P2g's terciles and the energy A do not move. 170's mechanism (the amplitude A pushing past the cliff)
    does not apply by construction.
  - Its case rests on beating the closed global form on the same folds, in the same run (C2, C3), which is the right
    test.
- **Searched, not closed:**
  - `grep -riE "per-target (γ|gamma)|tier(-specific)? (γ|gamma|whiten)|gated (un-?whiten|gamma)|(γ|gamma) per (target|tier)"`
    over `experiments/*/RESULTS.md` and FINDINGS. The only hits are 468 / 467 (a DE programme prior, an unrelated class).
  - 246's row-subset un-whitened control `C2_partial_ROWSELF` (`experiments/246-direction-hypotheses/RESULTS.md:272-276`:
    +1.36 % cos, pds −0.00012, 2/3) is the nearest prior. It selected rows by CD4-partial availability, not strength.
  - 246 §3.2 (`:36-41`) found that un-whitening's pds cost sits in target-shared directions. That is a live risk for
    strong targets that share programmes, and it is exactly what C2 tests.
- **Learned routing (066 / 263 / 269) is closed, and this is not routing.** The gate is a fixed public statistic
  (`PROPOSAL-frozen.md:148`, 522's def_strength_key522).

### 2. Falsifier and bins written before any number? **Yes.**
- Five falsifiers in checker units on BULK-ABC cancer plus BULK-W: `PROPOSAL-frozen.md:99-106`. Any one closes the family
  for October.
- Expected points and ranges per arm and surface: `:69-97`. The non-blind parts are labelled.
- A STOP rule for an instrument that does not see the closure: `:79-82`.
- The verdict ladder is code written before any price: `decide810.py:108-136`. Its plumbing was checked on 804's existing
  arrays without new pricing (`810/RESULTS.md:156-159`).
- Minor gaps, not blocking:
  - The bins are given as point ranges plus the expected ladder outcome ("most likely HOLD-707", `:96-97`), not labelled
    SHIP / WEAK / NULL / HARM per surface. Read them as: BULK-ABC WEAK expected (lower bound ≤ 0 after e(6)); BULK-W SHIP
    expected.
  - The "net-positive but lower bound ≤ 0 → flag as unexpected" clause (`:82`) is not in `decide810.py`. The worker must
    write that flag by hand in RESULTS.

### 3. Measurable on the checker's surfaces? **mse yes. pds is the weak leg, guarded by the STOP comparator. DE members are pinned. The BULK-ABC operator is a confound for THIS mechanism (G1).**
- **mse:** read through r. BULK-ABC's mse tracks the board with ρ_mse 0.89–1.65 (`experiments/800-checker/PHASEB.md:9-11`;
  `docs/LAB-LOG.md` part 1).
- **pds:** literal fold pds does not track the board for direction changes (`PHASEB.md:42-53`, over-read ~5× for
  sharpening). Here pds is the cost. If the checker does not see un-whitening's pds price at all, `sgu1all` reads
  clearly positive and the rule STOPs (`:79-81`). This is the correct guard.
- **DE members and ballast:** fid / jac / reach / nmae are pinned by unchanged carriers. In engine801, the arm changes
  `fin` only and `P_DE` comes from `fPP` (`experiments/801-lab-folds/engine801.py:130, 142-143`). |dK_rel| > 5 % routes the
  arm to DE-H1 first (`:119`, `decide810.py:115-116`).
- **Row norms (my lookup):** `sgu_rows` rescales the non-panel part to the FULL row norm |fin_i| (`price810.py:63`). That
  preserves |fin_i| only if fin's panel-gene columns are zero. I checked one dump of each format:
  - `experiments/801-lab-folds/work/dumps/dump708-efp-k562_gwps.npz` `fin`: 0 / 146 rows with panel mass;
  - `experiments/518-bulk-mse-gap/work/dump518-k562_gwps.npz` `fST`: 0 / 269 rows.
  - So the norm claim holds, and the sgu0 == inc known answer should pass. If any fold has panel mass, that KA aborts the
    run (it fails closed).
- **The operator confound, logged 08:05Z, six minutes before PROPOSAL-2 was written (not addressed in it):**
  - `docs/LAB-LOG.md` (08:05Z caveat) and `experiments/712-*/RESULTS.md:20-31`: BULK-ABC's fold operators are 062's cloud
    whiteners, whose basis is every donor row, K562 essentials in, with no 12k cap.
  - The L folds, production e2f65ecf and the 10-22 day path use `--basis-prune essential` with the 12k cap.
  - 810's ABC parts are built on 518's fST, i.e. on the 062 operator (`810/RESULTS.md:57`).
  - 810's measured quantity is literally `cos(pre,T) − cos(base,T)`, the whitener's effect. So the deciding surface
    measures the swap against a whitener that neither the board entries nor the day path use.
  - The STOP and C2 / C3 comparisons are internally consistent: sgu1 and sgu1all share the same 062 base. The absolute
    ABC read and any DOSSIER are not.
  - 801 is building ABCp folds with the day recipe (`engine801.py` AbcFold variant, `dump518p.py`). None is registered
    yet (`registry.json` has no ABCp key at 08:15Z). → **G1.**

### 4. Does it use seen data as evidence of gain? **No.**
- The arms use only `F.sources` and the fold's own objects (`price810.py:70-93`). BULK-ABC's truths are public held-out
  lines on A/B/C's panel (the checker's design, `PHASEB.md:7-8`), not A/B/C truth.
- No board read is evidence of gain:
  - 707's board read enters only to close (HARM → CLOSE) or as a class ρ for the projection (`:128-136`).
  - Calibrating ρ on a board point is the lab's protocol (LAB-LOG part 2, "divides by its class ρ"). The way this
    proposal uses it can amplify a read (G2).
- No sealed reserve is used. BULK-W's L folds are the lab's working set. The sealed L7 / L8 are not named.

### 5. Is the try budget honest? **Yes.**
- One priced configuration, charged N = 6 for the six composition-level forms read on these 22 folds
  (`PROPOSAL-frozen.md:58-61`; `810/RESULTS.md:162-167`).
- The comparators (inc, 707, sgu1all) are not tries.
- That is well under JUDGE.md's > 10 threshold, so no nested-fold scheme is needed.
- check800's ledger will log N = 1 against decide810's N = 6 (`:113-114`). → **G4** asks that N = 6 be written into
  candidate.json so the ledger is not understated.

### 6. Is the expected gain worth its cost? **Yes for the falsifiers and the STOP read. Marginal as a ship candidate.**
- **Expected:** BULK-ABC net +0.0015 (0 … +0.0035), bootstrap sd ≈ 0.0009 (`:91, 96`).
  - After e(6) the point is ≈ +0.0004. That is just above the bulk-direction readability of ~0.0003 (JUDGE.md).
  - The lower bound is expected ≤ 0. The proposal itself predicts HOLD-707 (`:96-97`).
- **Cost:** ≈ $1.5 of a $50 lane with $0.73 spent (`:138-142`).
- **What it buys:**
  - a well-posed close / park decision for this family;
  - a direct test of whether BULK-ABC sees the closed γ knife-edge, which every future direction proposal leans on.

### Class and day-of path. **Accepted.**
- Class is transfer (`:39-41`). Day-of: parts recomposed from the day's bundle, 3-donor operator and CD4 banks, with a
  sum-to-dump KA ≤ 1e-12 (measured 3.6e-14). z_src comes from 522's def_strength_key522 (`:144-149`). Nothing is fitted
  to A/B/C.
- One requirement for any later DOSSIER: the implementation must show the builder's DE carriers are unchanged by the
  swap (a KA on the carrier rows). The bulk engine pins them by construction (`fPP`), but the builder's derivation of
  carriers from the stack must be checked, not assumed.

---

## Binding conditions on this GO (no change to the frozen pricing code)

**G1. No DOSSIER on 062-operator BULK-ABC alone.**
- If decide810 returns DOSSIER, or the HOLD-707 branch reaches DOSSIER, re-price the same five arms on 801's ABCp folds
  (day-recipe held-out operators) before any dossier is written:
  - same code, same N = 6, parts recomposed by decomp810 on ABCp;
  - one launch, ≤ $5.
- Require D1–D3 and C1–C4 to hold there too.
- CLOSE, STOP, DE-H1 FIRST and park outcomes read as written on the 062 run. A false close only forgoes a small gain,
  and the asymmetry is acceptable.
- Report BULK-W's own E(sgu1) vs E(sgu1all) and its paired sgu1 − sgu1all beside the ABC ones. decide810 already
  computes them for every set (`decide810.py:145`). BULK-W is the surface on the day-recipe operator.

**G2. No amplification by a single board point.**
- In the HOLD-707 projection (`PROPOSAL-frozen.md:135-136`), use projected = (sgu1's BULK-ABC point after e(6)) / max(ρ, 1).
- Reason: c707's own checker interval is [−0.00092, +0.00178] (`810/RESULTS.md:48`, 804 prospective). It spans zero, so
  ρ = c707 / b707 is not identified below 1. A board read larger than the checker's must not turn a near-zero checker
  point into a DOSSIER.
- The CLOSE-on-HARM and park clauses stand as written.

**G3. The known-answer gate must be complete.**
- `price810.py:107-111` skips the inc / 707 == 804 KA when 804's array is missing on the worker. It does not fail.
- Before reading any number, confirm that every fold's `meta/<fold>.json` `ka` dict carries `inc == 804 {x,a,pds,pdsl}`
  and `707 == 804 {x,a,pds,pdsl}`: 22 folds × 8 keys, all < 1e-12, plus `parts sum == F.fin` and `sgu0 == inc` ×5.
- A missing key counts as a failed KA (rule 0, `:117`). The worker reports the count in RESULTS.

**G4. Record keeping.**
- Write N = 6 into candidate.json (or its notes) for the check800 aggregate, so the wear ledger is not understated.
- Add 063 (`experiments/063-pds-mse-pareto/RESULTS.md:70`; `docs/FINDINGS.md:785` "residual de-whitening for mse") to
  RESULTS' closed-family list next to 097 / 099 §5.3 / 170 / 212. `sgu1all` stands in for it.
- If sgu1all reads net-positive with lower bound ≤ 0, write the "unexpected" flag (`:82`) by hand.
- Fix the cosmetic "D1-D4" in decide810's DOSSIER message, or note that it means D1–D3.

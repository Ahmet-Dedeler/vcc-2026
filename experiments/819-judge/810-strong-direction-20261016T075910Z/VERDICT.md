# VERDICT: experiments/810-strong-direction/PROPOSAL.md — **REVISE**

Judge run 819 (written 2026-10-16 ~08:10Z). I computed nothing beyond lookups in files already on disk.

**Which file I judged.** The mandate names PROPOSAL.md sha256 `d38ed7e32296d325…`. When I read it, the file was sha256
`a421765b87d936640b9daed92c1f03228d2b13373b8a7e1f88dd59e6792ae341`, mtime 07:59Z. That is the minute the judge launched,
and the 810 worker (pid 6481) is still running. This verdict covers the `a421765b…` text. The header says "~08:15Z",
but the file existed at 07:59Z. The resubmission must be a frozen file that nobody edits while the judge reads it.

**Summary.** I do not KILL it. The closed construction is GLOBAL whitening relaxation. A strength-gated swap of only the
base part, which keeps H1 and CD4 and preserves each row's norm, is not closed anywhere in the record. The 810 MANDATE
lists "tier-specific whitening γ" as an admissible idea, and the run is cheap. It does not get GO yet, for four reasons:
- its central falsifier cannot be computed from the run it plans;
- the screen has no arm that reproduces the closure its novelty rests on;
- the mechanism paragraph attributes a predecessor construction's numbers to the candidate;
- the try count leaves out the composition-level reads of these same arms on these same 22 folds.

---

## The JUDGE.md checks

### 1. Is the mechanism closed already? **Not at this construction. The global form it borrows from is closed.**
- **Closed:**
  - FINDINGS §2 `docs/FINDINGS.md:787` ("Operator geometry: γ/λ (193), exponent (212)…") and §1 `docs/FINDINGS.md:578`
    ("the operator geometry is closed on all three axes").
  - 212 `experiments/212-gamma-cosine/RESULTS.md:26-30`: a global γ change cancels to 3.8 % (mse +0.3655 / pds −0.3517
    per unit cos).
  - 212 `:46-48`: "do not reopen it again".
  - 097 §9h `experiments/097-direction/RESULTS.md:1272-1277`: γ 0 nets −0.0033 (h1-test) and −0.0080 (orion). Gain/spend
    ratio 0.77 / 0.37.
  - Residual de-whitening and per-target arm routing (066): `docs/FINDINGS.md:785-786`.
- **Related, not closures of this construction:**
  - 263 and 269 are per-target donor selection and a learned router; they are closed
    (`experiments/269-utility-router-feasibility/RESULTS.md:1`, `…/263-donor-selection-arbiters/RESULTS.md:1`).
  - 246 §3.2: the pds cost of un-whitened rows sits in their top-k target-shared directions
    (`experiments/246-direction-hypotheses/RESULTS.md:36-41`).
  - 246 §3.7: the row-subset un-whitened control `C2_partial_ROWSELF` read +1.36 % cos at pds −0.00012
    (`:272`).
  - 810 is a fixed tertile gate, not a learned one, and it swaps only the base part. The proposal's §5
    (`PROPOSAL.md:93-108`) names the differences correctly.
- **Why it is close to the closure:** the swap 810 proposes is the whitener's input against its output, i.e. γ → 0 on
  the base part. Its own global form has the larger composition gain:
  - "pre1 all" +41.87 % on BULK-W vs "pre1 strong" +33.34 %;
  - on BULK-ABC cancer, +18.37 % vs +8.60 % (`experiments/810-strong-direction/work/foldlevel810.txt:1-2`).

  So the whole case is the claim that the gate improves the mse/pds exchange relative to the closed global form. That
  claim has to be measured on the same folds in the same run (item A below).

### 2. Falsifier and bins written before any number? **Partly. One falsifier cannot be computed, and the r leg is not blind.**
- **Bins:** per-arm bins exist for BULK-ABC and BULK-W (`PROPOSAL.md:72-76`). Four falsifiers are listed (`:84-91`).
- **Falsifier 2 cannot be computed** (`PROPOSAL.md:87-89`, "against 212's global η"):
  - The run prices inc, 707, sgu0 (the KA arm), sgu05, sgu1 and sga1 (`experiments/810-strong-direction/price810.py:118-123`).
    There is no global-unwhitening arm, so "the panel-average pds cost of global unwhitening" does not exist on these
    folds.
  - 212's η is overall pds per unit of aggregate cos, on 162's 4-fold PLOLO and a September stack. The falsifier compares
    it with raw Δpds per moved target. The units differ and so do the instruments.
- **The r projections are not predictions before data:**
  - The projected realised r (`PROPOSAL.md:62-66`) is the composition-level read of the same three arms on the same 22
    folds (`foldlevel810.txt:1-3`), multiplied by 707's realised/composition ratio. 707's realised values are confirmed
    in `experiments/804-known-answers/work/prospective.json`: rows/707 r_pct 7.967 (W) and 1.98 (ABC).
  - Only the pds leg and the render/ballast realisation are blind. The proposal should say this.
- **No ship rule:**
  - On BULK-ABC, the board's panel, the proposal predicts WEAK for all three arms (`:74-76`).
  - It says BULK-W "are not board sizes" (`:79-80`).
  - It says the direction class has no ρ until 707 reads (`:81-82`).
  - It never says which result produces a DOSSIER, and which result means "hold for 707's board read".

### 3. Measurable on the checker's surfaces? **Yes for mse. pds is the deciding leg and is the weak instrument. The DE members are pinned acceptably.**
- mse is read through r with its calibration. BULK-ABC's mse side tracks the board (`experiments/800-checker/PHASEB.md:9-11`).
- The net rests on the literal fold pds leg, which PHASEB says does not track the board for direction changes
  (`PHASEB.md:42-53`: "price direction on mse (r)… treat a pds-only gain as unvalidated").
  - Here pds is a cost. A fold over-read of pds (~5× for sharpening) would bias the read against the candidate, so it is
    conservative.
  - That only holds if the checker sees the pds mechanism of un-whitening at all. Item A's global comparator is the check.
- fid, jac, reach and nmae are pinned by the unchanged carriers. Ballast K is reported, and |dK| > 5 % sends the arm to
  DE-H1 (`PROPOSAL.md:33-34, 91`; DESIGN §4 `experiments/800-checker/DESIGN.md` "Engines" table). That is acceptable.
- It is not the unpriceable 479sp2 class: the key w and per-target q are kept and rows keep their norm
  (`price810.py:4-5, 102-103`). The class does not appear in `PHASEB.md:14` or in the LAB-LOG 07:55Z note.

### 4. Does it use seen data as evidence of gain? **No (not a KILL). It does select on the pricing surface.**
- It uses no A/B/C board reads, no sealed reserve and no held line among the sources.
- `price810.py:51-76` and the fold objects use `F.sources` only.
- The agreement gate (P2) and the β grid {0.5, 1} were chosen after reading composition-level outcomes of these
  constructions on all 22 folds. That is wear on BULK-W and BULK-ABC, so item 5 charges it.

### 5. Is the try budget honest? **No, it is under-counted.**
- The proposal counts 3 configurations, e(3) = 0.85 (`PROPOSAL.md:43-53`).
- Before the proposal, the lane read on the same 22 pricing folds:
  - the six forms in `work/foldlevel810.txt:1-3` / `work/foldlevel810.json`;
  - the 16-column × 3-tier ablation grid in `work/table810.md:14-16, 38-40`.
- The agreement cut was made after reading the strong & agree split.
- At composition level P2 is below P1: "pre1 strong&agree" reads 29.58 % (W) and 7.55 % (ABC), against 33.34 % and 8.60 %
  for "pre1 strong" (`foldlevel810.txt:1-2`). So P2 exists for its pds exposure, not its cosine, and the proposal does not
  say so.
- None of this is over 10 *priced* configurations. But AGENTS 11.1c.4 counts tries on every local instrument.

### 6. Is the expected gain worth its cost? **Yes on cost, marginal on readability.**
- Cost: ≈ $1.5, ≤ $5 with one re-run, on r7i.8xlarge for ~30 min (`PROPOSAL.md:110-113`). That is far under the $50
  per-launch cap and the lane's $50.
- Predicted gain on BULK-ABC is +0.0009 … +0.0015, with the CI spanning zero.
  - After e(3) = 0.85 sd (sd ≈ 0.00065 from the stated ±0.0013), that is about +0.0004 … +0.0010.
  - That is above the bulk-direction readability of ~0.0003, but every arm's interval includes zero. With no class ρ,
    702 rule 4 keeps any board projection spanning zero.
- The read is worth buying for its falsifiers, not as a ship candidate, provided item A makes the decisive comparison
  computable.

### Day-of path
Acceptable as written (`PROPOSAL.md:115-123`):
- parts recomposed from the day's bundle, the day's 3-donor operator and the CD4 banks, with a sum-to-dump KA ≤ 1e-12
  (measured 3.6e-14);
- z_src computed by 522's def_strength_key522;
- agreement taken from the day's raw donor rows;
- nothing fitted to A/B/C. Class: transfer (`:39-41`). Agreed.

---

## What is missing (fix all of these, then resubmit once)

**A. Add the closed global construction as a same-run comparator and known-answer row.**
- Arm `sgu1all` is base→pre at β 1 on every target with a donor row, on all 22 folds, through the same engine. It is not
  a try of the family.
- State its expected answer from the closures:
  - on BULK-ABC cancer, net ≤ ≈ 0, with the pds leg eating ≥ ~95 % of the mse leg (212 `RESULTS.md:26-30`; 097
    `:1272-1277`).
  - If the checker instead reads `sgu1all` clearly net-positive (lower bound > 0), the checker does not see the
    mechanism that closed γ on this class. Then STOP: the gated read is unpriceable.
- Rewrite falsifier 2 against that arm, in the checker's units. For example: "sgu1's (overall pds leg / overall mse leg)
  on BULK-ABC is not smaller than sgu1all's, paired bootstrap point" → the gate does not concentrate cosine relative to
  pds → close.
- Add a ship condition next to it. The best gated arm must beat `sgu1all` paired on BULK-ABC, in addition to falsifier 3
  against 707.

**B. Correct the mechanism paragraph** (`PROPOSAL.md:8-18`).
- For the candidate (base→pre, β 1), the composition gain is positive on every tier. It is not ≈ 0 on weak and mid:
  - BULK-W weak / mid: +0.0055 (14/18) and +0.0074 (15/18);
  - BULK-ABC weak / mid: +0.0066 (2/3) and +0.0090 (3/3) (`work/table810.md:14-16, 38-40`).
- The per-target whitening loss is +0.0096 / +0.0103 / +0.0193 by tier (`work/explore810-L.txt:1, 6, 11`).
- The "≈ 0 on weak/mid" and "+0.0055 vs +0.0005 agree/disagree" numbers belong to the pre-only swap, a different
  construction. Say so.
- Then state the actual claim: the gate is justified only by pds exposure (energy-weighted mse vs count-weighted pds).

**C. Honest try count.**
- Count the composition-level reads of the six `foldlevel810` forms (and the selection of β and the agreement gate on
  the same 22 folds), and apply e(N) with N ≥ 6 (1.27 sd).
- Or reduce the priced family to ONE preregistered arm (sgu1) and keep sgu05 / sga1 out.
- Say plainly that the r projections are composition reads of these arms × 707's ratio, i.e. not blind.

**D. Decision rule.**
- State, before the run, which outcome produces a DOSSIER and which closes the family:
  - suggested DOSSIER condition: BULK-ABC lower bound > 0 after e(N), AND paired > `sgu1all`, AND paired > 707;
  - say which pds unit decides (literal, with pdsl beside, per `PHASEB.md:53`).
- State which outcome waits for 707's board read (10-17/18). The board projection stays a range including zero until then.

**E. Provenance.**
- The proposal cites "RESULTS.md §1–§2" and "price810.py's sha256 in RESULTS.md" (`PROPOSAL.md:4-5`), but no RESULTS.md
  exists in `experiments/810-strong-direction/`.
- Record the sha256 of the following in the proposal (or in a RESULTS.md written before resubmission):
  - `price810.py`, `decomp810.py`, `analyse810.py`, `foldlevel810.py`;
  - `work/analyse810.json`, `work/foldlevel810.json`, `work/table810.md`;
  - the S3 parts prefix.
- Freeze the file while the judge reads it.

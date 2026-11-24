# 625 — PREREGISTRATION: an A/B/C module readout of Arc's 2026 lab programme, carried to D/E/F

2026-10-02 ~12:20Z, council session [7772dc]. Council portfolio item #5 (filler; lowest value per slot; runs LAST). Accepted
by [d4cd5f] under four conditions and flagged to [4d0207] for Ahmet. **Written before any readout arm exists and before
any read.**

## 1. Why (held-out evidence; AGENTS 11.1c rule 3's condition)
- The Broad-16 (one lab, 12 working cancer lines; its 4 reserved lines never read) shows that a same-lab response
  programme learned from 3 lines transfers to held-out lines (622h):
  - gene level: +0.054 gate-sign accuracy vs CLEAN4, positive on 12/12 lines, unchanged with cross-batch triples;
  - 16-module majority: +0.020.
- Simulated with this exact decode-and-flip rule (`sim625.py`, `work/sim625.json`): at bar 0.10 the held-out
  gate-sign gain is **+0.018** (|l|-weighted +0.023), positive on 12/12 lines, 1 % of triples negative, ~4 of 16 modules
  flipped. At bar 0 it is +0.017 (worst line −0.0005).
- A/B/C and D/E/F are one Arc 2026 campaign (same protocol, cancer lines). 458d's harm came from Arc's 2025 H1 ESC
  programme, a different year, screen and line type.
- Value on D/E/F: ≈ +0.003 … +0.004 overall (635's rate: +0.045 gate-sign ≈ +0.008 … +0.010).
- Class: transfer, a lab-level quantity. On the A/B/C board this is model selection with held-out support. It is not
  D/E/F decoding (H-631 untouched).

## 2. Fixed objects (sha-pinned; `modules625.py`)
- `work/modules625.npz`, sha256 `3bce71790fc1fd675990118d6490b87b71d84010c2d37edfeb74e3f73e091576`.
  - module(g) = 2·octile(production CLEAN4 P_down; 391's `clean4_prior_challenge.npz`) + half(A/B/C-pooled control
    log-CPM; 622's `abc_cpm622.npz`), over the 6,346 genes CLEAN4 defines; −1 elsewhere. Sizes 279 … 515.
  - These are FIXED GENE SETS. The same sets are flipped on D/E/F.
- Design: the 16×16 Sylvester-Hadamard H. Arm 1 flips all 16 modules. Arm k (k = 2 … 16) flips the 8 modules with
  H[k−1, M] = −1. Two ADDITIVITY check arms flip module 0 alone and module 15 alone. **18 arms in total.**
- "Flip module M" = P_C4(g) → 1 − P_C4(g) for g ∈ M, in whatever CLEAN4 component the base's DE prior uses: the
  block prior, plus the hook source if the base has the hook (then inside S = z(P_C4 − ½) + 0.5·z(I) if the identity
  rider is on). Nothing else changes.

## 3. The base and the arms
- **Base B** is a scored production entry whose DE prior uses CLEAN4. Every readout arm is a PRIOR-ONLY variant of B:
  same builder and emitter sha, same seeds, so every emission draw is shared (468). The arm's diff against B's recorded
  config must be exactly the flip mask (AGENTS 11.1 rule 10).
- **Never mix the UP pass.** Either B and all 18 arms carry `VCC_DE_UPRN`, or none does (d4cd5f: delivered magnitude
  differs for UP vs DOWN).
- KA-0, before any slot: arm "no flips" rebuilds B's package byte for byte.

## 4. Decode (written now)
For each member m ∈ {nmae (PRIMARY), fid (cross-check), reach (report only)}, Δ_k = m(arm k) − m(B), raw, board-pooled
over A/B/C. With x_M = the effect of flipping module M alone, assuming additivity (disjoint genes):
- arm 1: S = Σ_M x_M = Δ_1;
- arm k ≥ 2: Σ_M H[k−1, M]·x_M = S − 2·Δ_k;
- so y = (S, S − 2Δ_2, …, S − 2Δ_16) and **x = Hᵀ y / 16**.

Reads are exact for this realisation (shared draws), so there is no seed noise in Δ. What remains is A/B/C's
finite-sample programme, which the transfer bar handles.

## 5. Additivity known answer (gates the transfer)
The two single-module arms must reproduce the decoded x_0 and x_15:
|Δ_check − x_M| ≤ 0.10·|x_M| + 1e-5 on nmae and fid. If either fails, delivered calls or magnitudes are not additive
across modules. Then use only modules whose nmae and fid decoded SIGNS agree AND whose |x_M^nmae| ranks in the top
half, and report the failure.

## 6. Transfer rule (D/E/F; preregistered)
Flip module M on D/E/F iff ALL of:
1. x_M^nmae < 0: flipping M lowers nmae on A/B/C, i.e. A/B/C's gate genes disagree with CLEAN4's sign in M;
2. x_M^fid > 0: flipping M raises fid (the second, independent member agrees);
3. |x_M^nmae| ≥ 0.25 × median_M |x_M^nmae|: no near-zero modules, the analogue of sim625's bar.

October path, on the day: P_C4 → 1 − P_C4 for genes in the flipped modules, then the day-of prior
S = z(P_C4' − ½) + 0.5·z(I) and the 612 key gate exactly as without the readout. The flip mask is a fixed 16-bit
vector recorded in `docs/october-decisions.json` with this file's sha.

## 7. Slots and order (d4cd5f's conditions)
- Only A/B/C slots the coordinator's arm pipeline leaves empty between 10-08 and 10-21. Never displaces 630 / 612 /
  identity / K. Runs last.
- If fewer than 18 slots are available, run arm 1 plus as many Hadamard rows as fit, and decode the partial system by
  least squares, reporting the rank. Not reaching 16 + 2 means no transfer: an incomplete decode stays report-only.
- Probe entries on the public board (~18): Ahmet was told (via [4d0207] and directly) and can veto.

## 8. Readings, bins written before the reads
- **GO:** the additivity KA passes and ≥ 1 module meets §6 → the mask enters the October recipe as a transfer-class
  switch, with this file + 622h + sim625 as its held-out evidence.
- **NULL:** no module meets §6, so A/B/C's programme agrees with CLEAN4 module by module → nothing changes. That is
  informative too: CLEAN4 is already Arc-like at this resolution.
- **ABORT:** KA-0 or the additivity KA fails twice → no transfer, and the failure is recorded.

## 9. What would make this wrong
- The Broad-16 lab-programme transfer may be weaker across Arc's lines. Arc picked D/E/F for diversity, and the
  Broad-16's 12 lines are diverse too, so this is partly covered.
- The fid / nmae reads on A/B/C pool three contexts. A programme that differs by context averages out; a per-context
  decode is report-only, if per-context members are available.
- Magnitude asymmetry (UP vs DOWN delivery) can bias the nmae read. §5's check is there for that, along with
  never mixing UPRN.

## Amendment 1 (2026-10-02 ~12:35Z, BEFORE any readout arm or read; from [21dd4c]'s review)
1. **Capacity weighting.** x_M is the DELIVERED change. The hook moves signs only through a small MDE add (diag612:
   2–6 % of signs flip, mostly by shrinkage), and block genes carry more delivered magnitude than uncalled genes. So
   |x_M| measures "A/B/C disagree with CLEAN4 in M" × "how strongly M's genes are delivered", and §6's relative bar
   favours high-capacity modules. On D/E/F the flip applies to every gene of a passing module, whatever its delivered
   capacity there. That is conservative for modules that read low on A/B/C because of capacity, not agreement. No
   change to the rule. This is a reading caveat.
2. **Pooled read only.** Captured board records carry members pooled over A/B/C, with no per-context values
   (`data/board-entries/*.json` checked), so a per-context sign-consistency requirement is impossible. A module flip that
   helps one context and hurts two can pass on the pool, and a line-specific programme is exactly what must not
   transfer. Mitigation already in the evidence: sim625 decodes from POOLED triples and is scored on held-out lines
   (+0.018, 12/12, 1 % of triples negative), so the pooling risk is priced in at the Broad-16's level of line diversity.
3. **Value per slot** ≈ +0.0002, 5–10× below the hook / UP / K arms. It stays LAST, and a partial decode (< 16 + 2) is
   report-only (§7, restated).

## Amendment 2 (2026-10-02 ~12:55Z, BEFORE any readout arm or read): 8 CLEAN4 octiles, 10 entries, the October decode rehearsed exactly
Context: Ahmet approved decoding D/E/F's own module programme through final entries (DEF-RUNBOOK §6 item 8, coordinator
9ce09a9f). The coordinator asked that the A/B/C readout rehearse exactly the October design.
- **Modules:** 8 CLEAN4 octiles. They use the same CLEAN4 file and octile edges as above; the expression half is
  dropped. File `work/modules625_v8.npz`, sha256 `c24e3304ffd63d77d76b4a74341e09a2ee8396d7e372955e9c6feee167133dc1`
  (`modules625_v8.py`).
- **Hadamard column map** (column j ↔ octile [3, 0, 1, 2, 4, 5, 6, 7][j]). Column 0, the noisiest decoded column
  (1.20σ vs 0.66σ), holds octile 3, the weakest prior (median |P − ½| 0.005).
- **Arms:** 8 Sylvester rows (row 0 flips all 8; rows 1–7 flip 4 each) + 2 single-octile additivity arms (octiles 0 and
  7, columns 1 and 7) = **10 entries**. On D/E/F's fused design, an optional 11th entry repeats row 0 with a different
  bulk switch to cut column 0's noise. On A/B/C a repeat is uninformative, because prior-only reads are exact.
- **Why:** Broad-16 transfer simulation: 8 octiles +0.0165 (worst line +0.0025, 12/12) vs the 16-module design's
  +0.0172. Same value, about half the entries.
- **Decode:** `decode625.py --design v8` (default). `--selftest` passes on both designs. It also checks the
  coordinator's noise arithmetic under iid per-entry noise: v8 column 0 sd 1.217σ (theory 1.199), others 0.665σ
  (0.661).
- **Unchanged:** §4's algebra, §5's additivity rule, §6's transfer rule (with 8 modules), §7/§8, and Amendment 1.

## Amendment 3 (2026-10-02 ~14:05Z, BEFORE any A/B/C readout arm or read; from 625d / 625e on H1's real emitter)
1. **Design: SINGLE-OCTILE prior-only probes replace the Hadamard rows.**
   - A prior-only arm's board read is exact for its realisation: shared emission draws (468). 625d confirms the bulk
     sums are bit-identical to the base in every flip run.
   - With exact reads, Hadamard averaging buys nothing, and single flips need no additivity assumption.
   - Arms: octile k flipped alone, k = 0 … 7 (8 entries), then ONE combination entry carrying the decided mask. The
     combination reads additivity of the chosen set directly.
   - That is 9 entries, against 10 under Amendment 2. §5's check arms are dropped.
   - `readdecode625.py --design single` (selftest PASS, 8cddcca2) does the read. Artifacts:
     `experiments/604-identity-prior/build604.py --lam <base λ> --flip <mask>` (4d0207, 0b87b2ba; flip semantics
     independently verified by 7772dc: the octile-2 mask changes exactly its 793 covered genes per context, pdown = 1 − P,
     psign negated).
   - Board arms: `VCC_SAMELINE=<npz>:<sha>:0`, with the hook read from the same artifact. The pre-625 KAs on the
     coordinator's list run before the first arm (NEXT-SLOT-DECISION 13:40Z).
2. **Decision member: nmae ALONE; fid becomes a no-harm check; reach is excluded.**
   - Flipping DOWN-prior octiles to UP EXPANDS the realised call set: octile 5 on H1 gave +281 calls per target, almost
     all off-octile. fid therefore reads call-set mechanics, not sign quality.
   - H1's Δreach comes from a few tiny-n_real targets' discrete jumps.
   - New §6: flip octile k iff
     - x_k^nmae < 0, AND
     - |x_k^nmae| ≥ max(0.25 × median_k |x_k^nmae|, 2·σ_k), with σ_k = H1's 625d seed sd of that octile's Δnmae ×
       √(85/900) (board pooling: 3 contexts × 300 targets vs H1's 85; the number is recorded below when the seeds land),
       AND
     - x_k^fid ≥ −0.002 raw (no material fid harm).
   - Each probe's realised calls per target vs the base (rc diagnostics) is recorded next to its read.
3. **Value, revised by 625e's delivery calibration.**
   - +0.01 of SOURCE gate-sign is worth ≈ +0.0011 overall through today's routes (block prior + hook), not the
     +0.0017–0.0022 assumed.
   - Broad-16's 8-octile transfer (+0.0165 gate sign) is therefore worth ≈ **+0.002 overall on D/E/F if A/B/C reads GO**,
     ~0 if NULL (H1 reads NULL: CLEAN4 agrees on every octile).
   - It stays LAST, in empty A/B/C slots only.
4. **D/E/F use (final phase).**
   - If A/B/C reads GO, the October recipe carries the mask (DEF_DAY_DE_PRIOR_FLIP, 604f, 10-08+), and ONE D/E/F
     prior-only entry (R0 with vs without the mask) verifies it in-sample. The full 8-probe decode on D/E/F is NOT planned:
     at ≈ +0.002 its value per final-phase entry is below FINAL-PHASE-FEEDBACK-PLAN §2's rows.
   - If A/B/C reads NULL, nothing enters the recipe and no D/E/F entry is spent.
5. Unchanged: modules (modules625_v8.npz, octile field), the base rule (§3: one frozen scored base by sha; never mix UPRN),
   KA-0, slots/order (§7), Amendment 1's caveats.

### Amendment 3, numbers recorded (2026-10-02 ~16:05Z, from H1 625d before any A/B/C read)
- σ_k (board-scale realisation sd of Δnmae) = H1 seed sd × √(85/900): octile 0 ≈ 1.0e-5, octile 5 ≈ 1.2e-5. Use **σ_k =
  1.2e-5 for every octile**: the larger measured value, applied uniformly. The 2σ bar (2.4e-5) is far below the
  0.25 × median bar in practice.
- Value if GO, revised by 625f's improvement-direction exchange: ≈ +0.003 overall on D/E/F, previously +0.002.
  Everything else in Amendment 3 is unchanged.

## Amendment 4 (2026-10-03 21:41Z; written BEFORE any 625 A/B/C entry exists)
**What changed.**
- [21dd4c]'s 661 (16570222) priced all 8 single-octile prior-only flips on 618 (fresh h1-trainval surface, faithful key,
  zero repair, 3 seeds).
  - Every octile is measurable.
  - At f 0.35 none passes §6. Octiles 0–3 and 6–7 degrade. 4 and 5 lower nmae but fail fid's no-harm check.
  - Octiles 4–7 expand the realised call set by +156 / +353 / +1,186 / +1,238 calls per target; octiles 0–3 by ≤ 5. That is
    the same mechanism as 625d's octile 5 (+281 on h1-test).

**The coordinator's decision (NEXT-SLOT-DECISION, 10-03):**
1. 625 shrinks to 5 slots: octiles 0–3 + the combination.
2. A passing octile also needs |Δcalls| ≤ 10 per target in its board rc diagnostics.
3. It stays LAST, in empty slots only.

**7772dc's evidence that octiles 0–3 cannot pay** (`octflip625b.py`, `work/octflip625b.{out,json}`). Per octile: how often
flipping would help a held-out line, i.e. CLEAN4's sign agreement with the line's truth on that octile's gate cells < 0.5.

| octile | Broad-16 lines where flipping helps (count / \|lfc\|-w) | mean agreement | H1 agreement |
|---|---|---:|---:|
| 0 | 0/12 / 0/12 | 0.847 | 0.623 |
| 1 | 0/12 / 0/12 | 0.768 | 0.584 |
| 2 | 0/12 / 0/12 | 0.743 | 0.564 |
| 3 | 0/12 / 0/12 | 0.583 | 0.515 |
| 4 | **12/12 / 12/12** | **0.399** | 0.502 |
| 5 | 5/12 / 6/12 | 0.483 | 0.545 |
| 6 | 3/12 / 4/12 | 0.597 | 0.594 |
| 7 | 0/12 / 0/12 | 0.801 | 0.633 |

- On octiles 0–3, CLEAN4 is right on every held-out line: 12 Broad KO lines (no line below 0.55), H1 and 618.
- An A/B/C read of 0–3 is therefore expected NULL on every octile. P(GO) is near 0 and the 5 entries buy almost no
  information.
- The transferable sign signal sits in octile 4 (and partly 5–6). Those are exactly the octiles whose block-prior flip
  expands calls, and there the fid cost outweighs the nmae gain at f 0.35 (661).
- **Recommendation: DROP 625's A/B/C readout.** The octile-4 finding moves to a separate, call-neutral proposal: deliver
  octile 4's UP sign on uncalled genes only (hook or 630-style UP pass), priced on 618. It is not part of 625. The
  coordinator decides.

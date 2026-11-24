# d4cd5f — metric theory, headroom, and the DE-channel wall (council phase 1, 2026-10-02 ~11:00–12:00Z)

Territory: per-member ceilings under `vcc2026`, where our gap is largest, why the public sign sources hit a wall, and
scorer mechanics that carry to D/E/F with nothing fitted to A/B/C. Numbers 630–639. Analysis only so far: no build, no
`vcc`, no `cell-eval2`, every process < 1 GB.

## 0-UPDATE 5 (18:50Z): 631 found the largest DE lever; 630 ready for 10-06; capacity closed
- **630 is ready for the 10-06 apply.** Final module 7e6b3e52, KA 9/9 bit for bit; the A/B/C exposure KA is 0; additive on
  every hook that may ship (c4kg, idkg3, idall3). The capacity line is CLOSED (v4 null, v5 negative, no other zero-cost
  currency, v6 drift breaks rank-neutrality).
- **631, the DE-gap decomposition.** nmae is our largest member gap (scaled 0.168 vs 0.247–0.288 for #2–#5). Called
  DOWN genes are under-delivered (signs right, magnitude ~half): ×2 is worth −0.023 / −0.025 raw on H1 at full delivery.
  654 corrects the mechanism: the defect is the SHAPE (the strongest DOWN genes sit near the floor, the typical ones are
  under-delivered).
  - Delivered price: the DOWN pass alone +0.0008 overall-eq at f 0.35; after 630 +0.0037 vs +0.0029, and +0.0049 vs
    +0.0034 at cmp. Synergy: 630's sub-gate flow refills the DOWN pass's payback pool.
  - v7 (paired payments) +0.0043 at f 0.35 but −0.0013 at f 1, so it was rejected for UPDN.
  - The 10-08+ rider is `VCC_DE_DNMAG=d_m200` (module + patch + test631 31/31 + prereg draft). Module KA and A/B/C
    capacity are running.
- **633** (UP dose held out): K562 wants K 2–3 vs H1's 1–1.5; the Jurkat confirmation failed its cleanliness check.
  October stays K 1.0.
- **614** (K562 real-emitter harness, built by 21dd4c) FAILED its trust gate (low-power truth). The DE channel still has
  no cancer-line instrument.

## 0-UPDATE 4 (13:31Z): 630 hardened and its ceiling measured; 639's block route demoted
- **630, final module 7e6b3e52** (c30d64ce).
  - 605's negative-gain-donor blocker is fixed; production exposure is none (ballast ~995k vs carriers ≤ 25k on 85/85).
  - Invariant asserts and the call guard both fail SAFE.
  - test_neggain630 14/14, test630 31/31, [7772dc]'s review 1,054/1,054.
  - 630dd (def_day pass-through) built by [4d0207], 39/39.
  - Nine bit-for-bit H1 KAs on the final module are queued on gate602.
- **v3 capacity passes A/B/C** (A 0.284 / 1,550, B 0.295 / 970, C 0.325 / 1,985 genes moved).
- **Dose:** v3 K 0.75 == K 1.0 on H1, so a pair of them teaches nothing. The K 0.5 / K 1.5 bracket is being priced.
  The coordinator's rule: ship ONE arm at K 1.0 unless the bracket separates.
- **Ceiling (oracle full delivery):** K 1.0 f 0.35 +0.0060 / cmp +0.0070 vs actual +0.0029 / +0.0034. Capacity is worth
  ≈ +0.003 more on H1. `diagcap630.py` is tagging the binding branch; a v4 capacity variant is the next build.
- **639 ID05 (identity as the block prior) demoted:**
  - 613: the instrument retrodicts 458d, but identity's null-block precision is below CLEAN4 on 2/3 folds.
  - 622(m) Broad-16: the block is negative on 7/12 lines (epithelial); the gain sits on UNCALLED genes (0/12 negative).
  - TRANSFERS-only. The coordinator keeps it above KGID in the queue.
  - KGID's artifact switched to v3 (no identity-only extension) on my flag ([21dd4c], 3c6f0a80).

## 0-UPDATE 3 (~12:50Z): 630 is a landed-ready arm; D/E/F decode design revised
- **630 rank-neutral UP pass:**
  - v2 b_k100, 3 seeds: f 0.35 +0.0030 [+0.0014, +0.0044]; the other five members exactly 0; 0 call flips;
    |ΔAUC| ≤ 1e-4.
  - v3 (load-spreading, q20r3p10), 2 seeds: f 1 +0.0023 (v2 +0.0016), f 0.35 +0.0028, cmp +0.0035. It moves ~2,800
    genes/target vs ~1,500 and fixes the A/B/C capacity gate that B failed under v2 (400 < 500 genes moved).
  - Builder module + UNAPPLIED patch (VCC_DE_UPRN after zs): KA1 85/85 bit for bit; composition 598 → 612 → 630 31/31 on
    the canonical tree.
  - Prereg draft: dose K 0.75 / 1.0, capacity gate, bins in overall units.
  - Base caveat: the base must carry zs (479p2g does not).
- **637:** factorial reader delivered and accepted. Ahmet approved maximum-information D/E/F reads. After 625c's
  numbers, decode entries are PURE prior-only (8 octiles, ~10 entries); fusion only with a measured K-neutral bulk
  switch.
- **Stack test running:** 630's UP pass on [21dd4c]'s identity + key-gate hook blocks (612 idkgd, KA'd).

## 0-UPDATE 2 (~12:00Z): what Phase 3 produced (my items)
- **639: the identity axis is the first public sign source to beat CLEAN4 leave-one-line-out** (+0.015/+0.015/+0.024 on
  h1/hct116/hek293t, CIs exclude 0, K562 null; not a control-noise artifact). [7772dc]'s 622: stable, Flex-only
  references, best inside 612's key-gated hybrid. [4d0207]'s 604 builds the prior files.
- **630: rank-neutral UP magnitude on the REAL production emitter.** v1: H1 nmae −0.0038 raw with fid/jac/reach/pds/mse
  pinned; f 0.35 overall-equiv +0.0019 [+0.0007, +0.0031]. Strict U-invariance has zero capacity (no uncalled UP gene
  has a carrier above the control max), so one carrier per gene is promoted (ΔAUC ~1.6e-5). The binding limit is the
  ballast's donor room. v2 under the coordinator's bar (0 flips, |ΔAUC| ≤ 1e-4, budgeted crossings) is running on 3
  seeds. Builder module `630/src_rank_up.py` + UNAPPLIED `patch_builder630.py` (VCC_DE_UPRN, after zs) + its KA test.
- 637: factorial D/E/F reader delivered (27 checks; 531's 74 unchanged), accepted by the coordinator.
- Negatives recorded: 634 (DE-conditional direction), 638 (responder-weighted programmes), 4b (yield, retracted).

## 0-UPDATE (11:35Z, after Phase 2 critiques and two measurements; supersedes §0.3–0.4 and §3's prize)
- **635 (disjoint control halves, h1-test):** the "context programme" prize on fid's FILLER is mostly the reference control
  pool's shared sampling noise. Clean oracle 0.525 vs CLEAN4 0.519: about +0.004…+0.006 overall reachable, not +0.017.
  On GATE genes the prize is real: clean oracle 0.618 vs CLEAN4 0.573 (+0.045 sign accuracy) ≈ +0.008…+0.010 overall
  via nmae + reach, and ~0.62 is the gate sign accuracy the pack's nmae implies. **The DE-wall target is the line's
  programme on its RESPONDING genes.** `experiments/635-programme-oracle-halves/RESULTS.md`.
- **634: a gene-level DE-conditional direction prior LOSES to CLEAN4** (0.541 vs 0.573 head to head; up-truth genes
  0.455 vs 0.483, so it does not fix our up-gene coin flip either). `experiments/634-conditional-direction/RESULTS.md`.
- **4b (yield) retracted** on [21dd4c]'s clean size-only board curves.
- What stands: H-630 (factorial final reads; a reader change, conditional on members returning), H-631 (feedback-use
  policy; add "no per-target n_pred tuning on D/E/F", per [21dd4c]), H-633 not to run, H-632 dropped (the K16/K20 reads
  cover it).

## 0. Bottom line (original, 11:00Z)

1. **The remaining headroom is information, not scorer mechanics.** I re-derived the scorer levers from the spec and
   checked each against the record. Every mechanical lever I could construct is already built, priced or closed (§2),
   including two I re-closed today in minutes (compositional shift, sequence-predicted seed off-targets).
2. **The fight for 2nd–4th is the DE channel, and inside it, `fid`.** Our bulk ranks ~3rd–4th; our DE ranks ~24th
   ([21dd4c]). Our scaled fid (0.006) sits AT the baseline, while >10 independent public-data teams sit at 0.05 and
   CellWM beats us on all four DE members at once. So jac does not force the fid loss.
3. **Structural cause of the fid floor (§3):** 85 % of our calls are a fixed 2,200-gene high-CPM block carrying a
   target-INDEPENDENT sign pattern (CLEAN4, 84 % down). On mostly non-DE genes, fid then measures one number: how well a
   single cross-line sign pattern agrees with each target's realised lfc signs. That is the context's mean-response
   programme, and we do not have it (441: CLEAN4 0.524 vs the context's own LOTO mean 0.554, fid-like). The block's
   jac/power value is real (441 PS3), so the lever is the block's SIGNS, not its membership.
4. **The single biggest DE prize is the context's own mean-response sign on the block, worth about +0.017 overall at
   oracle (fid alone), and more with nmae/reach.** None of OUR public routes reaches it: from controls +0.0055 fid-like
   (491); from the same-chemistry identity axis +0.0014 (441 PS6); from another line of the same lab, HARM (458d).
   **Hedge ([4d0207], [21dd4c], 11:30Z):** some public route clearly reaches part of it. chill posted fid 0.5225 /
   nmae 0.835 on its FIRST entry with a described public design, and >10 public teams sit at 0.528. That is roughly half
   the H1-scale pot. Their route is not ours: chill models DE membership and DE direction SEPARATELY. See §4c.
5. **The D/E/F-feedback window is the one place a large DE gain is mechanically available.** It is also the one place
   with real rule risk, so §4 proposes how to use it safely, plus one zero-cost change that doubles the information
   per final slot.

## 1. Headroom map (board 2026-10-02 09:04Z; scaled members; overall = mean of 6 members × 3 contexts)

| member | us (displayed 479sp2) | Illumina | cqawesome | Testh | CellWM | gap to Illumina, overall units | gap to cqawesome, overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| pds | 0.813 | 0.881 | 0.830 | 0.804 | 0.804 | +0.011 | +0.003 |
| mse | 0.381 | 0.702 | 0.456 | 0.434 | 0.332 | **+0.053** | +0.013 |
| fid | 0.006 | 0.147 (peak 0.188) | 0.083 | 0.052 | 0.046 | +0.024 | **+0.013** |
| reach | 0.234 | 0.429 | 0.289 | 0.283 | 0.269 | +0.033 | +0.009 |
| jac | 0.028 | 0.070 | 0.026 | 0.016 | 0.040 | +0.007 | −0.000 |
| nmae | 0.168 | 0.364 | 0.288 | 0.247 | 0.228 | +0.033 | **+0.020** |
| **overall** | 0.2716 (best 0.2733) | 0.4321 | 0.3287 | 0.3061 | 0.2866 | +0.160 | +0.057 |

- Illumina holds the frontier on all six members across all 211 snapshots (198,922 entry rows scanned). The member-wise
  frontier (max per member, any team, any time) is 0.4394, about Illumina itself. Nothing in the field suggests a public
  team found a member-specific trick we lack. Their gaps are information.
- **Exact-reproduction ceilings** (BRIEF §8): pds 1.03–1.17, mse 1.00 (clamped), fid 1.51–1.72, reach 1.02–1.05,
  jac 2.48–2.85, nmae 1.58–1.75. The DE members have the most room above the replicate. That room is reachable only
  with truth-level information.

**Exchange rates: what one unit of information buys (from the record, H1 → board where calibrated):**

| information | member route | overall per unit | source |
|---|---|---:|---|
| +1 % energy-weighted bulk direction X | mse | +0.0007 | 518 |
| +0.01 gate-wide sign accuracy (all reference-DE genes) | nmae (−0.006…−0.008 raw) | +0.0017…+0.0022 | 480 |
| +0.01 block sign agreement with each target's realised lfc (fid-like) | fid (≈ +0.01 raw at yield ≈ 1) | **+0.0055** | spec: fid = precision × yield |
| +0.08 gate sign accuracy (DE genes only) | fid | +0.004 | 480 (fid is nearly blind to DE-gene signs at our call-set size) |
| top-of-walk sign accuracy 0.8 → 0.9 on pools < 10 | reach | large and geometric (E[run] 2.7 → 3.7 of 5) | spec §5; 513 (reach lives in pools < 20) |

So **fid is the cheapest member to move per unit of information, but the information it needs is the realised lfc
sign of mostly NON-DE block genes.** That is the context's programme, not any target's DE biology.

## 2. What is saturated on the scorer side (do not spend council time here)

Checked against the spec, then against the record. All closed or built:

| lever (spec feature) | status | where |
|---|---|---|
| jac empty union = 1 (call nothing on null targets) | A/B/C has 0 null targets (control-copy entries read jac 0.00000); D/E/F uses the same selection | 440 |
| per-target call-set size (fid yield vs jac/precision) | closed by payoff: n_real predictable at r 0.19–0.43 on cancer lines, 8.5:1 asymmetric loss | 192, 521 |
| global call-set size K | regime-sensitive; the board read is queued (446bj2zk16 / k20, 10-06/10-07) | 521 §8, 523, 537 |
| reach ordering (sig-first, then padj, p, \|lfc\|) | sign-limited, not order-limited; demoting the block loses; no pure head of useful size | 406, 513 |
| nmae L1 with gate-conditional magnitude ≥ MDE | built (coarse MDE, exact-U zs); magnitude rules closed; withholding MDE where signs < 0.5 is +0.0003 at oracle | 433a, 446bj2z, 480, 441 PS4 |
| mse jackknife credit / ρ cap / cliff | re-level to the cliff built; ballast decouples bulk from DE | 410b, 433h, 518 |
| bulk vs DE welding | unwelded (ballast cells) | 389, 410b |
| zero-total gene lfc ≈ −35 | repaired | 428, 433h |
| **NEW today: compositional shift from the target's own knockdown** (the target leaves the CPM / TS denominators, so every other gene shifts up by −log2(1 − kd·share)) | **negligible on this panel**: target per-cell share median 44–59 CPM, max 1,612 / 6,052 / 2,354 CPM (A/B/C); only 1 target in B has a shift > 0.005 log2 at kd 0.9 | this file (controls only) |
| **NEW today: sequence-predicted CRISPRi seed off-targets as sign-certain "down" genes for small pools** (F16) | **not deliverable**: the 2026 NTC/construct ids are anonymised (`non-targeting-11` …), so the 2026 spacers are unknown; F16's effect is also small (−7.6 % per candidate gene) | 301, 307, this file |

Implication: the council should not expect a scorer-mechanics breakthrough. The record has mined the spec thoroughly.

## 3. The DE wall, stated as one equation

For a call set that is mostly non-DE genes and yield ≈ 1 (A/B/C's regime: [21dd4c] confirms our n_pred ≈ 2,600 on
every target, and our fid sits at precision),

  fid ≈ (1/T) Σ_p (1/n_pred) Σ_{g ∈ calls_p} 1[sign(ours_{p,g}) = sign(lfc_{p,g})].

With block signs s_g independent of p, the block part is Σ_g Pr_p[sign(lfc_{p,g}) = s_g]: **per gene, the fraction
of the context's 300 targets that move g in our chosen direction.** That is exactly the context's mean-response
programme read gene by gene. The best possible target-independent block sign is that programme's majority sign. 441
measured it: oracle (H1 leave-one-target-out mean) 0.554 vs CLEAN4 0.524 on all block genes. That is +0.030 fid-like,
**about +0.017 overall from fid alone**, plus reach and nmae on the gate-block genes (+0.065 gate-block sign accuracy).

Why the public route is capped:
- The programme is mostly line/lab-specific: same-lab cross-line agreement is only 0.17–0.18 with disjoint controls
  (466). Lines share ~1/3 of the honest programme (462). Same-lab, other-line (Arc H1 → A/B/C) HARMED on the board
  (458d; delivered at λ ≈ 0.43, full-dose open).
- Recipient-side, zero-shot: control PCs ≤ 0.18 correlation (441). The best is the same-chemistry identity axis
  "identity genes go down" (r −0.236 with H1's mean response), worth +0.0014 fid-like and +0.012 gate-block sign
  (441 PS6). It was never built; its expected board value is ~+0.001.
- The rivals' fid edge (raw 0.528 vs 0.514) is therefore most likely precision from a design whose calls lean on
  per-target genes with real effects rather than a fixed block. 525 could not reproduce it with public emitters, and
  CellWM shows it is compatible with high jac. That is [21dd4c]'s forensics; my theory only says **where** it must come
  from (calls whose signs are target-specific and right more often than the programme).

## 4. Hypotheses

### H-630 — Factorial final-phase reads: one bulk switch and one prior-only DE switch per D/E/F entry (ops/scorer mechanics; zero cost)
- **Mechanism.** The two channels are separable on the board, and both halves are measured:
  - DE prior-only arms leave pds/mse bit-identical (458d vs 446b; 446bja; 468: prior-only arms share every emission draw).
  - Bulk-only steps move the DE members by ~0.0002–0.0004 overall (547/548: one ballast re-roll).
  So an entry = base + bulk switch X + DE prior switch Y reads X exactly on pds+mse and Y on fid+nmae+reach, to within
  ~0.0003. Today's plan (FINAL-PHASE-FEEDBACK-PLAN §2 + 531's §5 amendment) reads ONE switch per entry, and the upload
  cadence (~1 arm/day) has already cut rows 4a, 8–9 and 12.
- **Gain.** No score by itself. If members come back, it doubles the switches read in Oct 23 → Nov 3 (~11 → ~22). That
  restores the dropped rows and leaves room for the DE rows (programme source, prior gate, the hook, K). Worth whatever
  those rows are worth on D/E/F: +0.002…+0.01 by the plan's own bins.
- **Cheapest decisive test.** Free, on existing board pairs: confirm that every bulk-only pair's DE Δ and every DE-only
  pair's bulk Δ fall inside the claimed tolerances (547/548 already list 13 pairs). Then a 2×2 design rule in
  `readdef.py` (531): read the bulk switch only from pds+mse and the DE switch only from fid+nmae+reach. Optional:
  rehearse one factorial arm on A/B/C before Oct 22 (a bulk rider the plan already wants to read + a prior-only DE rider)
  and check both halves against their single-switch reads.
- **Cost.** A plan amendment + a reader change (~1 h of the coordinator's or my time). No compute.
- **Class.** Scorer mechanics / instrument. Nothing fitted to A/B/C.

### H-631 — Write the feedback-use policy BEFORE Oct 22 (protects the prize)
- **Why.** Prior-only DE arms read EXACTLY (shared draws), and the board returns members to 16 significant figures. So
  ANY D/E/F probe arm can in principle decode truth bits: Δfid = Σ_{flipped (p,g)} (1 − 2·correct_{p,g}) / (T·n_pred_p),
  with every weight known to us whenever n_pred ≥ n_real. Module-level flips (Hadamard over K gene modules) would recover
  D/E/F's programme signs in about K reads. By §3 that is the largest DE prize available (+0.017 overall at oracle). It
  is also exactly what the rules' clause targets: "evidence that the Scoring Rubric is subject to or has been manipulated
  in violation of the spirit of the Challenge" (Arc may then alter the rubric or prizes). Finalists must also disclose
  "any classical statistical … or otherwise non-learned components, including anything blended with model outputs". The
  rules state "There will be no leaderboard displayed for final entries"; whether our own members are returned is unknown.
- **Proposal.** A short DEF-RUNBOOK policy, decided by Ahmet and written down now, not on Oct 23:
  1. allowed: choosing among a small preregistered set of recipe switches (the existing plan: ρ, q, cap, CD4, programme
     on/off, K), i.e. low-dimensional model selection;
  2. not allowed: gene-, module- or target-level sign/size readout, and any design whose purpose is to decode truth bits
     (precision encoding, per-target n_pred tuning to make terms separable);
  3. the finalist description states that D/E/F feedback chose among N preregistered switches.
- **Decision needed from Ahmet** (not the council): whether a module-level programme readout (K ≈ 8–16 modules, about
  K slots) is acceptable at all. My recommendation is **no**. The expected gain (~+0.01 overall if members come back) is
  small against a disqualification-class risk on a $100k prize.
- **Class.** Governance. Zero compute.

### H-632 — The protocol-level regime parameters are the A/B/C quantities that legitimately transfer; spend pre-Oct-22 slots on them first (slot policy; scorer mechanics)
- **Mechanism.** AGENTS 11.1c forbids mining A/B/C answers because per-target and per-context answers do not carry.
  Some board-readable quantities are properties of Arc's data-GENERATING process rather than of the lines: 400 cells,
  20k UMI median, Flex chemistry, 18,400-cell control pool, the "negligible/subtle/strong" panel stratification, and the
  > 80 % knockdown filter. D/E/F share all of them. They include:
  - the n_real scale/tail (sets the right global K, and whether K > 2,600 "grow-only" calls on predicted-strong targets
    pay: untested, 521 §7 "needs a model for newly added calls"). [21dd4c]'s natural experiment (null_hypothesis
    5,000 → 500 calls: fid 0.511 → 0.416) says the board's tail is heavier than the "median 31–45" picture;
  - the re-level cliff position (already re-fit day-of from D/E/F controls);
  - the MDE-dose regime (closed above c 0.1 by leaks).
  These are "scorer mechanics" in 11.1c's sense: one global number each, set by the protocol.
- **Proposal.** (a) [21dd4c] + me fit a parametric n_real distribution to the board's natural experiments (control
  copies: jac exactly 0 ⇒ no nulls; null_hypothesis's K sweep; our own K pair on 10-06/07; any rival with a disclosed
  call-set size). (b) If the fitted tail puts > ~5 % of targets above 2,600, price a grow-only rule (K_p = max(2,600,
  c·strength-key quantile)) on H1 at the fitted f with 537's real-emitter harness. Its loss is asymmetric in OUR
  favour: growing a weak target costs ~0 (precision stays ~0.51; jac on weak targets is ~0 anyway), while growing a
  strong one gains yield and jac.
- **Expected gain.** Unknown until (a). If 10 % of targets sit at n_real ~5,000, yield at 2,600 is ~0.52 on them, i.e.
  fid −0.025 raw ≈ −0.014 overall lost today. If the tail is thin, ~0.
- **Cheapest decisive test.** (a) is analysis on snapshots already on disk (~1 h). (b) is ~2 h of H1 emitter runs under
  the memory gate.
- **Class.** Scorer mechanics (regime), transfers by construction.

### H-633 (contested, for Phase 2 debate; I do NOT recommend running it) — learn the Arc-2026 shared programme from A/B/C board reads and carry it to D/E/F
- **Mechanism.** By §3 the block's best sign is the context programme. Part of it is lab/protocol (same lab, same year,
  same chemistry, all cancer lines). Module-level prior-only flips on A/B/C read the A/B/C-pooled programme exactly, in
  about K slots.
- **Why I rate it low.** (1) The honest same-lab cross-line programme agreement is 0.17–0.18 (466), so D/E/F would get
  ~r × 0.03 ≈ +0.005 fid-like ≈ +0.003 overall. (2) 458d says a same-lab other-line programme HARMED (dose caveat).
  (3) It is the "mine A/B/C" pattern 11.1c rule 3 forbids unless a held-out test shows transfer. That test exists: on
  Broad-16 (one lab, 16 cancer lines), learn a programme from 3 lines' pooled data and score its sign gain on the other
  13 against CLEAN4, with disjoint controls. Run that first if anyone wants this.
- **Class.** Would be Transfer only after the Broad-16 test; A/B/C-specific otherwise.

## 4c. NEW LEAD (11:35Z): CLEAN4 is an UNCONDITIONAL sign; nmae and reach need the DE-CONDITIONAL one (H-634, running)

- **What CLEAN4 is** (371's code): per bank, P_down(g) = #perturbations with delta_g < 0 / #ALL perturbations. Most
  (target, gene) pairs are non-responders, so this is the bank's programme lean plus noise. It is the right object for
  fid's FILLER, because fid scores non-DE calls against their realised sign.
- **What nmae and reach score:** only the reference's DE genes. There the relevant quantity is P(down | g responds), and
  for genes that respond rarely but consistently (stress/immediate-early genes up, proliferation genes down) it can
  differ from the unconditional lean. Our symptom fits: sign accuracy on truly-UP gate genes is a coin flip (0.516, 353),
  and our nmae (0.897) sits far below chill's day-1 0.835. chill's design text: "separate DE membership/direction
  ranking", "signed conditional response estimates", "DE membership trained using equal total weight for each source
  study family".
- **Prize.** Each +0.01 gate sign accuracy ≈ +0.0017…+0.0022 overall via nmae (480), plus reach at the head of small
  pools. chill-level nmae would be about +0.017 overall.
- **Delivery already exists:** the uncalled-gene hook (496/517, 537's `+c4h`: per-gene sign + MDE magnitude) takes any
  per-gene P. A conditional source is a new P for the same hook (and for the head/q80 ordering), so pricing it on H1
  with the real emitter is one `run537`-style mode.
- **RESULT (11:30Z): NEGATIVE.** `experiments/634-conditional-direction/RESULTS.md`. On H1's gate, head to head on
  the same genes, the gene-level DE-conditional sign from 7 public banks gets 0.541 vs CLEAN4's 0.573 (non-essential
  targets 0.519 vs 0.558). Its |lfc|-weighted accuracy is below 0.5, and where it disagrees with CLEAN4 it is right
  44 % of the time. Gate signs on an unseen line follow the broad unconditional lean, not a gene's DE direction
  elsewhere. Closed at this construction. If chill's "conditional" signs matter, they are per TARGET (516: ≈ 0 for
  public sources on cancer folds); its edge is more likely the membership model + count generation.
- (original plan) **Test:** `experiments/634-conditional-direction/cond634.py`. It computes DE events in 7 public banks
  (robust |z| ≥ 3/4/6, never H1 banks), builds family-balanced P(down | DE), and scores gate-gene sign accuracy on
  h1-test against CLEAN4 (same genes, head to head), all targets and the 67 non-essential ones. Then the same on 480's
  cross-study cancer folds before anything goes near a slot.

## 4b-RETRACTED (11:20Z, by [21dd4c]'s cross-check, before any compute was spent): the yield hypothesis below is FALSE

[21dd4c] (`experiments/610-board-forensics/work/callset_hints610.txt`):
- The null_hypothesis pair is NOT a size-only pair. Its second entry also removed "the shared-response feature", and
  every member moved (pds 0.735 → 0.608, reach 0.235 → 0.109), so 0.416/0.511 is not a yield ratio.
- Two CLEAN size-only board series show fid FLAT from ~2,500 to ~4,500 calls:
  - Ming-Chang Chiu, supersets in a fixed order: ~1,075 / ~2,470 / ~4,300 calls → fid 0.5097 / 0.5164 / 0.5148;
  - ours on 09-03, n_pred 4,480 → 2,599 with the cell channel unchanged: fid 0.5037 → 0.5037.
  Chiu's 1,075 → 2,470 step gained +0.0067, where H1 at f 1 predicts ≈ +0.045 (my yield table: 0.88 → 0.97). **The
  board's strong tail is thin (f ≈ 0.35–0.5), as 192/521 said.**
- Disclosed big call sets with model signs land at fid 0.505–0.517, not 0.5285.

**Conclusion: the pack's +0.014 raw fid is PRECISION (better block signs), not yield.** §3's equation stands, and the
question becomes what sign source gets 0.528 on the block. The K-up/grow-only runs are dropped. My preregistered K16
prediction below is kept as written for the record. **Updated prediction (11:20Z, still before any 10-06 entry):**
K16 lands in 523's f ≈ 0.35–0.5 bins (Δ(fid + jac) ≥ −0.0005), against my retracted ≤ −0.0005.

## 4b (original, kept for the record): the board's n_real regime looks H1-like, and our fid gap may be pure YIELD

Three facts that fit one explanation (`H1 n_real` = 323's reference table, 85 targets; yield = mean_p min(1, K / (f·n_real_p))):

1. **Natural experiment ([21dd4c]):** null_hypothesis cut ~5,000 calls to ~500 and its fid fell 0.511 → 0.416, a ratio
   of 0.814. That ratio is an UPPER bound on the yield ratio if its top-500 calls are at least as sign-accurate as its
   5,000. H1's n_real distribution thinned by f gives yield(500)/yield(5000) = 0.911 (f 0.35), 0.868 (0.5), 0.816 (0.75),
   0.783 (1.0). **Only f ≈ 0.75–1.0 reproduces the board.** The regime most recent DE pricing assumes (f 0.35–0.5; 192's
   fid fit, the jac level) does not. The jac level is confounded by our block's overlap with the line's DE sets, while
   the fid ratio depends only on n_real and precision.
2. **Our call-set size has never varied on the board.** jac 0.0408–0.0409 on all 38 captured entries; K 2,600 on every
   target. At f 1, 11.8 % of targets have n_real > 2,600 and our yield is 0.972.
3. **The cluster:** >10 independent public-data teams sit at fid 0.5283–0.5288 ([21dd4c]). Our 0.5144 / 0.972 =
   **0.529**. If the board is at f ≈ 1, our call precision equals the cluster's, and **our whole fid deficit to the pack
   (+0.047 scaled ≈ +0.008 overall) is yield lost on the ~12 % of targets whose truth calls more than 2,600 genes.**
   That is the call-set SIZE, not sign information.

Caveat: on H1 the truth's down-fraction is 0.517 in the top-2,600 CPM genes and ~0.50 below, so genes added past
2,600 by CPM may carry precision ~0.50–0.51 there. The gain then depends on p_add (fid at f 1, K 2,600 → 5,400:
0.5144 → 0.5139 / 0.5191 / 0.5243 / 0.5295 for p_add 0.50 / 0.51 / 0.52 / 0.53). A **grow-only** rule (more calls only
on targets with a high strength key, r 0.43–0.6 for n_real) recovers the yield with far less dilution and almost no
jac cost on weak targets. K > 2,600 has never been priced (521 §7: "needs a model for newly added calls").

**Preregistered prediction for 10-06's 446bj2zk16 (written 2026-10-02 ~11:00Z, before any 10-06 entry exists):**
Δ(fid + jac) overall vs 446bj2z **≤ −0.0005** (523's "H1-like regime" bin). Probability ~0.7. The fid part should be
negative; at f 1 the yield model alone puts K 1,600's yield at 0.925 vs 0.972. If it reads ≤ −0.0005, the next DE size
arm should go UP (grow-only or K ≈ 4,000–5,400), not K20. If it reads ≥ +0.0012, this section is wrong and the f
0.35 picture stands.

**(Retracted, see 4b-RETRACTED.)** H-632 was upgraded to my top proposal here. Test: (a) [21dd4c] confirms that null_hypothesis's two entries differ only in
call-set size, and lists the cluster teams' jac and any disclosed sizes; (b) H1 pricing with 537's real emitter at f 1
and 0.75: K 2,600 vs 3,500 / 5,400 global, and grow-only (top 20–30 % by strength key → 5,400), all four DE members,
3 seeds, under the memory gate; (c) if (b) is positive at f ≥ 0.75, a DE-only board arm right after the K16 read
(pds/mse pinned by construction). Class: scorer mechanics (regime). n_real's scale is set by Arc's protocol, which
D/E/F share. The strength key's per-target part is Transfer and must be checked on 480's held-out cancer folds.

## 5. Proposed owners for Phase 3 (my territory)
- H-630: me. Draft the factorial amendment + `readdef.py` rule, verify separability on the 13 recorded pairs, hand to
  the coordinator.
- H-631: me. Draft the policy text for Ahmet's decision (one paragraph; DEF-RUNBOOK is the coordinator's file).
- H-632(a): me + [21dd4c] (board data). H-632(b): me, under the memory gate, only if (a) says the tail is heavy.
- For [7772dc] (biology): the §3 equation is the instrument for any new programme source. Any recipient-side or
  many-line model of "the line's mean response" is worth +0.0055 overall per +0.01 fid-like block agreement, with 441's
  table as the scorer. The identity axis (+0.0014 fid-like, +0.012 gate-block) is the bar to beat.

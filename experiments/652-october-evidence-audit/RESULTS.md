# 652 — October recipe evidence audit (AGENTS 11.1c): no DECIDED switch rests on a VOID read; eight flags, four need action before the freezes

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** In this
> audit only the 604 identity-block-prior row cites Orion fold numbers (613's HCT116 −0.014 / HEK293T −0.008). Those are VOID.
> The row's NULL expectation rests on h1 +0.006 and Broad-16 −0.001. The validation status of 613's instrument is being
> re-ruled by the coordinator. F1–F8 do not rest on Orion DE truth. F7 concerns Orion's bank rows as same-line
> inputs, not its DE truth.

2026-10-02 16:55–17:04Z (commit time), council session [7772dc], at the coordinator's request. Read-only: nothing was built, scored or
changed outside this directory.

**Sources:**
- `docs/october-decisions.json` (HEAD 60de7801);
- 530's and 538's switch additions (`patch_bundle530.py`, `patch_extras538.py`), which land on 10-05;
- `experiments/613-block-prior-validation/INSTRUMENTS.md`;
- 534's §2–§5;
- `docs/BOARD-LEDGER.md`;
- `docs/NEXT-SLOT-DECISION.md` through 16:40Z;
- `docs/FINAL-PHASE-FEEDBACK-PLAN.md`;
- each rider's own RESULTS: 598, 604, 612, 630, 642, 644, 537 §3/§5, 493, 522, 250/251/256;
- `scripts/def_pipeline.py` SHIPPED_ENV / DEF_ENV.

## 0. Verdict

1. **No DECIDED switch rests on a VOID read.** Each of the ten DECIDED switches is in one of three states:
   - it has an A/B/C board read AND a held-out read from an instrument with a retrodiction record;
   - it is operational;
   - it was decided OFF on a board HARM that a validated instrument retrodicts (613 → 458d).
2. **The VOID instruments touched only things that are not in the October path.**
   - The VOID instruments: the call-set-blind fold nmae models (BETA·consensus: kgh612b/c/d, hyb622b, ds612c), the
     zero-base sign × MDE model (kgh612, dosecurve610), 631f's fold model (failed its retrodiction) and 625c's module-net
     fold model (superseded by 625 FINAL's real-emitter rate).
   - What they touched: 612's KEYGATE expected-value case (withdrawn; the knob is inert and KGID retired), 631's ×2
     called-DOWN dose (not a switch), and 625's first sizing (replaced by the real-emitter rate).
   - The coordinator's 612 CORRECTION 2 consequences (NEXT-SLOT-DECISION, 10-02) record that the plain CLEAN4 hook's
     October case (`mde_sign`, 496a on 10-05) never rested on the folds. Its held-out read is H1's real emitter.
3. **Eight flags. Four need action before the 10-05 / 10-06 freezes:**
   - **F1:** Arc hook / `hook_estimator`, MZ, CPMZ: no valid held-out read. Drop it.
   - **F2:** `callset_block` / CPH. Its price comes from the provisional instrument fitted on its own failed predecessor
     (cp).
   - **F3:** `rc_callset_cpm`. H1 flips sign across regimes, so it is board-only.
   - **F4:** tercilec (598, decided at 10-06 step 9 on KAs alone). Its DE side is unpriced and can re-tier up to half of
     an October panel.

   The other four are awareness items or record fixes: F5 (held-out-only switches), F6 (A/B/C-calibrated constants),
   F7 (642's same-line value), F8 (612 KEYGATE's record).

## 1. Instruments and their status
| instrument | prices | retrodiction record | status |
|---|---|---|---|
| 438 levelled route, k = 1, on real / efp / efr (3 cancer folds each; 522 reuses it) | bulk shape, key and stack steps | 433k (q 1.5) predicted +0.0044, board +0.0046; 446b (q 1.75) +0.0021 vs +0.0023; 446bj2 (LODO) +0.0027 vs +0.0022. It over-admits: in 5 of 534's 12 rejections a local held-out read predicted a gain and the board said no | **validated** for bulk steps. Its panels are tuned-on (11.1c.5) |
| kolf2 PLOLO fold | bulk shape | 417a: pds +0.0063 predicted vs +0.0064; the 2-donor cancer folds over-read 2.5× | validated (shape) |
| Sealed R6 (534 Part B; C10a run 09-29, C10b due 10-04) | the October bulk recipe | sealed, read once | **the confirmation instrument** (11.1c.5) |
| H1 real emitter at f 0.35 (611) | DE delivery levers (calls unchanged) | zs nmae 1.18×, coarse MDE 1.59× | **validated** for delivery. Reach provisional |
| H1 real emitter, call-set changes (611) | fid × ~2, jac × ~0.5 | one pair (cp). That pair was a NULL read against a +0.0020 … +0.0025 prediction | **provisional**. K16 / K20 are its first out-of-sample tests |
| 613 restricted-LOLO block | block-prior sources | 458d HARM; 392 fid gain | **validated** (2 board points) |
| 613d add-rule hook | hook sources on uncalled genes | 2 out-of-sample passes on H1's real emitter (idkg3, idall3, within 0.0003) | validated on H1. **Cross-surface exam = 496a on 10-05** (registered +0.00017) |
| Fold nmae models; zero-base sign × MDE; 631f; 625c | — | wrong-signed or failed retrodictions | **VOID** |
| H1 for H1's own programme | — | circular: H1 is the source's own line | **not a held-out read** |

## 2. Every switch and rider

Tries = arms on the board / local variants for this switch (11.1c.4); "≥" where the record does not give an exact count.
Ship = what 11.1c.2 needs beyond what exists.

### 2a. `docs/october-decisions.json`, DECIDED
| switch = value | class | A/B/C board read | held-out read | instrument status | tries | ship | flag |
|---|---|---|---|---|---|---|---|
| `q` = 1.75 | transfer | 417a +0.0070 (q 1.25); 433k +0.0046 (q 1.5); 446b +0.0023 (q 1.75) | 438 §3.4: q 1.5 step +0.0041 … +0.0049, q 1.75 +0.0020 … +0.0030, 3/3 folds per panel; R6 retention 0.53 / 0.61 | validated + sealed | 3 / 438's dose curve | shipped | — |
| `lodo` = 1 | transfer | 446bj2 +0.0022; 460dr (October path) +0.0029 | 438 §3.5 +0.0031 / +0.0032 / +0.0017; R6 +0.00099 [+0.00033, +0.00181], retention 0.30 (the weakest), 7/12; permuted key −0.0102 | validated + sealed | 2 / ≥ 5 (438; 541/542/544 aggregates, closed) | shipped | — |
| `sharpen_pt` = tercile:2.5,2.0,1.75 | transfer | 479p2g +0.0010 TRANSFERS; 479sp2 −0.0007 FAILS | PLOLO on the LODO stack P2g +0.0035, 3/3; kolf2 +0.0020; 9/9 October-like folds (493, 540c) | validated (bulk). The same instrument ranked S_P2 above P2g (+0.0038), and S_P2 failed on a DE side-effect the bulk instrument cannot see (low-tier ballast nmae +0.0018 raw) | 2 / ≥ 6 rules (493) | shipped | F4 (upgrade) |
| `cd4_stack` = ×1 | ×1.5 is A/B/C-specific | 457dr: ×1.5 +0.0031 | 438 −0.0011 / −0.0024 / −0.0016; R6 −0.00087 | validated + sealed | 2 (417e, 457dr) / 438 | shipped ×1 | — (11.1c's model case) |
| `de_mag` = dn_k250_zs | scorer mechanics | 446bj2z nmae −0.0022 raw; other five members bit-identical | H1 real emitter: −0.0013 … −0.0029 predicted | validated (one of its two retrodictions) | 1 / 456 (zs on 3 bases) | shipped | — (A/B/C capacity 0.70× H1's; D/E/F capacity unmeasured; 554's dose guard catches a collapse) |
| `programme_source` = off | (transfer, rejected) | 458d −0.0085, HARM | 613: the H1 programme is below CLEAN4 on all 3 cancer folds; 629: 0/12 Broad lines in any form, de-identified included | validated (613 retrodicts 458d) | 3 (446bja, 458d, 458g) | decided off | — |
| `prior_gate` = off | — | 458g +0.0036 vs 458d (still −0.0049 vs 446b) | as above | — | 1 | moot | — |
| `sameline_w` = 0.5, `sameline_flags` = noprog | transfer (same-line rows) | none possible: no A/B/C line has a covering bank (Nadig Jurkat is essentials-only). The programme form failed for A (446bja −0.0011) | 450: an other-lab same-line bank on the H1 truth, +0.053 … +0.079 at 35 % coverage, plus a K562 replication; w chosen in 487 | h1-test is board-validated; 450 is one surface | 1 (446bja) / 487's w grid | held-out-only by necessity, written in 487 | F7 (HEK293) |
| `max_levels` = 4; `cd4_thresholds` (day-of) | operational | — | 448 / 452 rehearsals; 481 F1 | — | — | — | — |

### 2b. PENDING (2 on HEAD, plus the six that 530 / 538 add on 10-05)
| switch | class | A/B/C board read | held-out read | instrument status | tries | ship | flag |
|---|---|---|---|---|---|---|---|
| `lodo_hi_pct` | transfer | 495 on 10-03 (on A/B/C the pct cap equals cap 6; predicted +0.0019 … +0.0029 vs cap 3) | 522: vs cap 3 +19.4 / +11.6 / +14.1 % direction; vs cap 6 +5.95 / +6.16 / +5.41 %; R6 +0.00694 [+0.00517, +0.00828], retention 3.31× | validated route + sealed | 0 / ≥ 3 cap forms (522, 540) | 10-03 bins (a flat read fails the instrument's own prediction) | — |
| `mde_sign` (clean4:0.05) | transfer (cross-line prior) | 496a on 10-05 | H1 real emitter, 3 seeds: nmae −0.00146 raw; 613d | validated (H1 delivery). 613d's cross-surface exam is this read | 0 / ≥ 6 hook forms (496, 612, 613d) | 10-05 bins. A form change to idall3 is a second read (row 604-hook) | — |
| `rho_scale` | scorer mechanics with an A/B/C-calibrated target | none possible | none possible: it needs D/E/F truth | — | 0 | FINAL-PHASE row 2 only | F6 |
| `rc_callset_cpm` | scorer mechanics (call-set size) | K16 on 10-06, K20 on 10-07 | H1 at f 0.35: K 1,600 +0.0048 [+0.0009, +0.0100] (521); f 1: −0.015 | **provisional**, and the sign flips across regimes | 0 / ≥ 4 K values (521, 537 P2/P3, 533 D10) | board read + a written reason | **F3** |
| `hook_estimator` (meanz_ne, the Arc hook; arms MZ / CPMZ) | transfer (same-lab programme) | MZ / CPMZ sit in 537's 10-07 / 10-08 branches | H1 only (circular). Its A/B/C attenuation r 0.37 is FITTED to 458d's board read. 517 §2.5's same-lab partner analogues (3/3) are unvalidated. P12: the floor without sign information is −0.003 vs CPH. 629: Arc's programme is uninformative off-line in every form | **none valid** | 0 / ≥ 4 (517, 527, 537 P4 / P12) | — | **F1** |
| `lodo_strength` (θ 0.25) | transfer (off-board value) | ~0 on A/B/C by construction | 522: +0.7 … +1.9 % direction on 9/9 October-like folds | validated route, but tuned-on panels (11.1c.5) | 0 / 522's priced arm | sealed R6 (C10b) + a written reason | F5 |
| `lodo_hi` (fixed cap) | transfer | none | none (531 #12: it must clear a held-out line) | — | 0 / 0 | cannot ship before Oct 22; FINAL-PHASE row 4a only | F5 |
| `callset_block` (528's call-propensity block; CPH = block + CLEAN4 hook) | scorer mechanics (call composition) | 446bj2zcp NULL: −0.0008 overall, Δ(fid + jac) −0.00018 vs predicted +0.0020 … +0.0025. CPH on 10-07 | H1 only (528, 537 §3.2) | **provisional**, fitted on cp, the arm that failed | 1 / ≥ 3 (528: lcap, lcap-dec, kept; 549's CPH) | a written reason, or 614's K562 price after its retrodiction | **F2** |

### 2c. The 10-06 riders (and 604f, 10-08+)
| rider | class | A/B/C board read | held-out read | instrument status | tries | ship | flag |
|---|---|---|---|---|---|---|---|
| 598 tercilec (`sharpen_pt` upgrade; decided at 10-06 step 9 "on KA evidence") | transfer-rule fix | none: it re-tiers 4 of 300 A/B/C targets (E2F4, SLC25A13 HIGH → MID; EPHB2, TMEM209 LOW → MID) | 493h: efp +0.84 % rel cos (3/3) vs plain +0.48 % (2/3). real / efr have 0 non-carried targets, so no read there | validated route, ONE panel, best of 4 tie rules | 0 / 4 | see F4 | **F4** |
| 612 / 612dd `VCC_DE_MDE_SIGN_KEYGATE` | robustness (EV case withdrawn) | none | the fold evidence is VOID (612 CORRECTION 2). H1 real emitter: KGID ≈ the plain hook. KGID retired 15:22Z | VOID / retired | 0 / ≥ 5 (qscan, KGID v1 / v3, ds612 / ds612b) | never set | F8 |
| 630 / 630dd `VCC_DE_UPRN` = b_k100_q20r3p10 | scorer mechanics (delivery) | the 10-08 pair | H1 real emitter v3, 3 seeds: f 1 +0.0023, f 0.35 +0.0029 [+0.0015, +0.0044], cmp +0.0034 overall-eq; 0 call flips; stacks on every hook. Dose K 1.0 chosen on H1 (633: K562 descriptive 1.5–3; 633j unclean) | validated for delivery (nmae; reach provisional) | 0 / ≥ 8 (v1 K 1.0 / 1.25, v2, v3, v4 NULL, v5 negative, v6 breaks neutrality, K 0.75) | the 10-08 read | — (D/E/F capacity unmeasured; the A/B/C gate is ≥ 50 % of H1's) |
| 642 HEK293 contingency (default OFF) | day-of instrument | — | the identity call: real-Flex HEK293 8/8 at 0.820, 0 false calls of 423 (nonflex reference) | validated (KA + cross-source positives) | 4 reference variants | a written decision on the day | F7 |
| 604 identity hook, idall3 (`mde_sign` = hookbank:identity:0) | transfer | slot 4, read vs 496a | H1 real emitter, 3 seeds, vs c4h: +0.0014 / +0.0011 / +0.0006 overall (f 1 / 0.5 / 0.35). 613d's registered prediction held within 0.0003. Broad-12 (622o): worse than c4h on 5–8 of 12 lines with a production-like key | 613d validated on H1; cross-surface pending 496a | 0 / 613d best-of-5 (stated in its prereg) | board WIN; if 496a fails 613d, a written reason | — (tails stated in its prereg) |
| 604 identity block prior (ID10 first, then ID05) | transfer | empty slots only | 613 block instrument: null-block precision below CLEAN4 (HCT116 −0.014, HEK293T −0.008); block reads h1 +0.006, Broad-16 −0.001; NULL expected | validated | 0 / 639's λ × 4, 604 v1–v3 | a TRANSFERS read only | — |
| 604f flip mask (625 readout; 10-08+) | transfer (same-lab) | 625's 8 + 1 A/B/C probes | Broad-16 8 octiles +0.0165 gate sign (held out); exchange rate from H1's real emitter (625 FINAL) | validated (real emitter). Broad-16 is a KO modality (627) | 0 / preregistered (transfer3 rule) | a GO read + ONE D/E/F verification | — |
| 644 v3 `DEF_DAY_SAMELINE_SKIP` | operations | — | 4/4 tests | — | — | — | — |

The 10-05 bundle's guards and instruments (538's pass-throughs, 554's dose guard, 590's identity floor, 641's depth row)
change no prediction. Each carries KAs, so no evidence row is owed beyond the tests in the v2 checklist.

### 2d. Base-recipe constants in the October path (not switches; listed so the record is complete)
| constant | origin | off-board evidence | note |
|---|---|---|---|
| ρ* = 0.9793 (re-level target, NTC-proxy margin) | A/B/C: the EV margin below the proxy cliff, with the proxy error calibrated on A/B/C (ctxlevel417 / 187) | none possible before Oct 22 | F6. FINAL-PHASE row 2 re-reads it if members return |
| `VCC_ALPHA` 0.982 | the a11 → a15 board series (A/B/C), closed by 211 | none | the re-level rescales the bulk effect row, so its bulk role is overridden. Whether it reaches the DE channel's delivered magnitudes is not documented (F6, low) |
| `VCC_CD4_DEFLATE_BASIS` = V5_abc (A/B/C's pinned basis) | A/B/C residual SVD | 250: real pds +0.0153, efr +0.0107 at −0.1 … −2.0 % cosine | OK |
| normrel cap form (319), CD4 partial rows (319p) | board-derived / A/B/C | 534 C: both ≈ 0 off A/B/C (3/9, 4/9 folds) | kept with a written reason (lower max regret, below dayof438's 0.0010 threshold) |

## 3. Flags and recommendations
**F1 — the Arc hook (`hook_estimator` = meanz_ne; MZ / CPMZ in 537's 10-07 / 10-08 branches). Recommend: DROP from the
October path.**
- It has no valid held-out read. H1 is the programme's own line.
- The A/B/C value rests on an attenuation (r 0.37) fitted to 458d's board read, so an A/B/C "win" partly re-reads its
  own calibration.
- 537's P12 puts its D/E/F floor at −0.003 vs CPH when sign information is absent.
- 629 (today) closes the remaining route: Arc's H1 programme, raw, de-identified, detrended or restricted, is 0.06–0.17
  gate sign below CLEAN4 on 3 folds and 0/12 Broad lines.
- If MZ / CPMZ keep a slot as pure information, preregister that no A/B/C read switches October to the Arc hook, so it
  cannot become a board-only switch by default (11.1c.2–3).

**F2 — `callset_block` / CPH (10-07). Recommend: read CPH as information; October adoption needs a written reason or 614.**
- The only instrument that priced it is H1's call-set mode, fitted on one pair (cp). That pair is this lever's own
  failed predecessor: NULL against a +0.0020 … +0.0025 prediction, with jac over-read 2.3× and fid missed.
- A CPH win would therefore be board-only under 11.1c.2.
- The cancer-line DE harness 614 (K562 real emitter; trust verdict due 10-04 05:00Z) is the first instrument that could
  give it a held-out read. Price CPH there after 614 passes its 611-style retrodiction, before 10-07.

**F3 — `rc_callset_cpm` (K16 / K20). Recommend: board-only with a written reason, plus a day-of gate.**
- H1 changes the sign across regimes (K 1,600: f 0.35 +0.0048; f 1 −0.015). The board's regime is the open question
  these reads answer (PREDICTIONS-611-K16-K20.md grades the instrument).
- A K picked on A/B/C's n_real regime is depth-tuned. 641's depth row must stay a hard stop for any non-default K on
  D/E/F.
- Recommended written rule: if D/E/F's control depth or the predicted n_real distribution leaves A/B/C's range, K
  reverts to 2,200.

**F4 — 598 tercilec at 10-06 step 9. Recommend: price its DE side before deciding it; add it to C10b.**
- Step 9 decides `sharpen_pt` tercilec "with --evidence <KA598b/c>". KAs prove the code, not the effect.
- Its effect evidence is ONE panel (efp, 3 folds) and the best of 4 tie rules.
- On A/B/C it re-tiers 4 targets, so no board read can show it.
- On an October panel with thin donor coverage it puts every non-carried target in MID. On efp ~160 of 300 targets are
  non-carried, so up to half the panel can change tier.
- A tier change moves ballast. That is exactly the DE side-effect that failed S_P2 on the board (low-tier ballast nmae
  +0.0018 raw), and the bulk instrument cannot see it.
- Cheapest fix: one H1 real-emitter run of tercilec vs tercile with a coverage-thinned key (so non-carried targets
  exist), priced on nmae at f 0.35, plus tercilec in R6's C10b.
- Otherwise decide it with a written reason. The tie at w = 1 is an artefact, and MID is the neutral tier.

**F5 — held-out-only switches (`lodo_strength`, `lodo_hi`).** An A/B/C read cannot show either one.
- `lodo_strength`'s 9/9 folds are on panels that have already chosen several settings (11.1c.5). It needs R6 (C10b)
  before it is decided, plus a written reason.
- `lodo_hi` has no held-out read. It stays a FINAL-PHASE row 4a candidate only, and must not be decided before Oct 22.

**F6 — A/B/C-calibrated constants (ρ* 0.9793; α 0.982).**
- ρ* is already named "the single most D/E/F-specific number in the recipe" (FINAL-PHASE row 2). Accept it with that
  written reason: the margin protects an asymmetric cliff, and the proxy is re-fitted from D/E/F's controls.
- α: one grep or a KA to confirm it does not reach the DE channel's delivered magnitudes. Low priority.

**F7 — 642's same-line value is unpriced.**
- The HEK293 identity call is validated. The value of switching same-line ON with Orion HEK293T rows for a Flex
  HEK293 line has never been measured: no Flex HEK293 perturbation truth exists, and Orion's rows carry its technical
  axis.
- If it fires on Oct 22, the coordinator's written decision should cite 450's mechanism and this gap, and consider
  w < 0.5.

**F8 — 612 KEYGATE's record.** 644's PENDING section still lists 612 as awaiting a read. KGID is retired, and the
fold evidence is VOID. Mark it RETIRED (never set on the day), so nobody switches it on from the runbook.

## 4. Answers to the brief, in one line each
- Any October decision silently on a void read? **No.**
- A/B/C-specific in the October path? Three items, all already disclosed with written reasons:
  - CD4 ×1.5 is excluded;
  - the cap form and partial rows are kept on a regret argument;
  - ρ* is A/B/C-calibrated.
  α is a fourth, smaller item, to confirm.
- Switches that would ship board-only if they win: CPH (F2), K16 / K20 (F3), MZ / CPMZ (F1). Write the reason before the
  read, or drop.
- A decision on KAs alone: tercilec (F4).

## 5. Addendum (17:25Z): F6's α check, answered from the code and the record
- **α reaches the DE channel, but only as the carriers' design amplitude.** In `scripts/make_submission_pm_cd4.py`, `VCC_ALPHA`
  scales the effect matrix. That matrix becomes the DE carriers' design (`context_rows` → `profile` → `rc_emit_k.rc_hybrid_block_k`,
  ~5987–6030) and enters the ballast count (`ballast_count(profile, profile_bulk, …)`). The re-level (`VCC_BULK_CTX_SCALE`)
  scales only `bulk_rows`, i.e. the ballast render. So the bulk members follow the re-levelled bulk, and α keeps only the
  carriers' amplitude.
- **The DE channel's sensitivity to that amplitude is measured (410 §7, H1, bulk pinned by ballast):**
  - carriers ×1.25 → nmae +0.0025, fid +0.0001, reach −0.0105 raw;
  - ×1.5 → nmae +0.0123, reach −0.0137;
  - ×2 → nmae +0.0287, reach −0.0107 ("DE amplitude is not the lever").
- **No A/B/C-fitted DE choice rides on α.** Its board fit (a11 → a15, 211) chose the WELDED bulk's mse cliff, a role the
  re-level now takes from each D/E/F context's own controls. Its DE role was never tuned on a board read. **No action for
  October.**
- **Open, unpriced:** nothing below ×1 has been measured. 410 §7's trend is monotone-worse above ×1. A carriers ×0.8 vs ×1.25
  pair on today's construction (coarse MDE, zs, q80; the H1 real emitter, which is validated for delivery levers) would show
  whether lower amplitude helps. If so, that is a scorer-mechanics knob. Offered to the coordinator; it needs an experiment
  number.

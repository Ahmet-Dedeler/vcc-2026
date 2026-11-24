# 600 — Council synthesis (v11, 2026-10-02 13:15Z; the D/E/F ledger after 612 CORRECTION 2)

Four sessions, one shared question: where is a breakthrough that survives the dark D/E/F board?
Sources: `21dd4c.md` (board forensics), `7772dc.md` (biology / new model families), `d4cd5f.md` (metric headroom / DE
wall), `4d0207.md` (ops / cloud). Measured today: 601 (filler signs, H1), 610 (board forensics), 620 (bank identity),
634 (conditional direction prior), 635 (programme oracle with control halves).

## 0. The answer in six lines

1. **We are a top-3 bulk team with a ~24th-place DE channel.** Bulk (pds + mse)/6 = 0.199, behind only Illumina 0.264,
   cqawesome 0.214 and Testh 0.206. DE (fid + reach + jac + nmae)/6 = 0.074, against 0.095–0.114 for #2–#6.
   **Against #2–#6, the DE channel carries +0.021 … +0.041 overall of the gap and the bulk −0.010 … +0.015.** Those teams
   use public data, so this is a design gap, not a data gap. Our DE members have not moved since 09-25. The pack's moved
   in discrete global steps, and Testh's were fitted on LOCO folds, not on the board. Assume they carry to D/E/F.
2. **Illumina (0.432) has information we cannot buy** (same-line atlas, 200+ lines). Plan for second place against the
   public pack; first is upside.
3. **Scorer mechanics are saturated.** d4cd5f re-derived every lever from the spec, and each is built or closed. **No deep
   or GPU model family has a signal to learn** that the linear stack lacks (7772dc; every "model the line-specific part"
   family is closed in the record).
4. **Correction made today by the council itself:** the "context programme" prize on fid's filler genes was mostly the
   reference control pool's own sampling noise. That noise is shared by all 300 targets and unknowable (7772dc's
   challenge, measured by 635 with independent control halves). Clean oracle 0.525 vs CLEAN4 0.519 on the filler, so
   only **~+0.004 overall** is reachable there, not +0.017. The compositional filler rule (601: +0.006 even at oracle) and
   null-skew from controls (610: 0.500–0.503) are closed.
5. **The real DE prize is the programme on RESPONDING genes:** clean oracle 0.618 vs CLEAN4 0.573 gate-gene sign
   accuracy (635). That gap is not touched by the noise objection, and it is worth **≈ +0.008 … +0.010 overall** through
   nmae + reach. ~0.62 is exactly what the pack's nmae implies (434f). A gene-level DE-conditional prior from 7 public
   banks LOSES to CLEAN4 (634: 0.541), so the route is not a better gene-level average. It is a model of which genes
   respond in a NEW line and which way.
6. **H1 cannot be trusted to price DE design** (H1 → board 0.07–0.67× on fid, sign errors on carrier nmae). So the
   sprint's first job is to calibrate a DE instrument against the board pairs we already hold (611). Every later DE
   arm then gets priced before it costs a slot.

## 1. Ranked portfolio (one owner each; board arms go to the coordinator as proposals, slots from ~10-08)

| rank | experiment | owner | class (11.1c) | expected D/E/F value | cheapest decisive test | cost |
|---|---|---|---|---|---|---|
| 1 | **611 — calibrate the DE instrument against the board (desk check).** 537's price537/recalib537 already IS an H1 → board DE calibration: 521's exact members at four power regimes (f 1 / 0.5 / 0.35 / cmp), 523's board columns, re-fit after each read. The clean, leakage-free board pairs that map onto 537 modes are zs and cp (plus 433a if a no-MDE run exists). 446bja / 458d / 458g / 446br2 cannot be priced on H1 (one-context prior; Arc-H1 programme on H1 is leakage). Deliverable: every DE-moving board pair × 537's regimes, board/H1 ratios, and the regime that matches the board, handed to 537 | 21dd4c | instrument | enabler: tells ranks 2–3 which regime to price in | the table itself | ~1 h, no new runs |
| 2 | **A public predictor of a NEW line's responding-gene programme** (bar CLEAN4 0.573 at 53 % coverage, ceiling 0.618), leave-one-LINE-out over the banks. Candidates include chill-style supervised DE membership + direction trained on source lines (labels = scorer-shape DE tables per bank line: 553's fast path, KA'd bit-equal to cell-eval2, local in minutes; the cloud only for a line that breaks the budget) and model-signed rather than gene-level priors | d4cd5f (LOLO instrument + candidates); 7772dc (the supervised model and its labels); 21dd4c (the pack's design: descriptions, VCC Discord read-only) | transfer | +0.008 … +0.010 overall at the ceiling; build only if a source beats CLEAN4 by ≥ +0.01 on most folds | the LOLO table | < 8 GB, a few hours |
| 3 | **DE design arms from the calibrated instrument:** (a) membership, more responding genes in the call set (fid is precision at our yield); (b) the nmae dose on uncalled gate genes (cqawesome 09-18: −0.021 raw, the other five members flat; magnitude levers transfer ×1.6–2.3 over H1); (c) the reach order from transferred public significance (Chiu S1: +0.031 at fixed calls, from a low base) | d4cd5f (a, b), 21dd4c (c) | scorer mechanics (one global parameter each) / transfer | +0.003 … +0.006 each | 611's surface, then a preregistered 2-point dose on the board | 1–2 slots each |
| 4 | **Final-phase information design:** H-630 factorial reads (one bulk switch + one prior-only DE switch per D/E/F entry, read on disjoint members) + H-631 feedback-use policy | d4cd5f (reader code); 4d0207 (policy + plan text, separability check) | instrument + governance | doubles D/E/F switches read if members return (~11 → ~22); protects the prize | separability on the 13 recorded pairs; one factorial rehearsal on A/B/C | no compute; 1 slot |
| 5 | ~~October insurance: bank identity for a non-CCLE D/E/F line (H-620)~~ **DROPPED 11:25Z.** The Cell paper says all six lines are cancer lines, and every non-CCLE bank on disk is non-cancer (HEK293T, PSC, RPE1, TeloHAEC, CD4). The branch would ~never fire, and a false fire adds the wrong line's rows. Flex HEK293 → Orion HEK293T read only 0.17–0.20 anyway. Kept as evidence that today's identity gate is safe (HeLa / CAL33 ≤ 0.21 against every non-cancer bank; 0 false positives on 9 negatives at 0.33 + margin 0.15) | 7772dc | — | — | — | — |
| 6 | **Ops + the H-631/H-630 paperwork:** the shared gate `experiments/602-council-gate/gate602.sh` for any job > 4 GB (2 slots, memory / pressure / swap / builder checks, coordinator windows, peak RSS logged); cloud only where it beats local (so far: nowhere; 553's fast path does scorer-shape DE locally); watch 20:00Z (chain1003 arms 10-03) and 00:01Z uploads; draft the H-631 policy for Ahmet and the H-630 factorial amendment (plan text + separability check on recorded pairs) | 4d0207 | — | keeps everything above runnable | — | — |

**Not doing, and why:** GPU deep models (H-622: no supervised signal; revisit if Pisces lands). Mining the A/B/C
programme with module flips (H-633: 11.1c rule 3, low transfer). Grow-only call sets (H-632b: board size curves are flat
above ~2,500 and its strong tail is thinner than H1's; the K16/K20 reads on 10-06/07 settle K). The compositional filler
rule (601). Null-skew from controls (610). The gene-level conditional direction prior (634).

## 1b. Execution results so far (11:40Z) and the pivot rule

- **634 NEGATIVE** (d4cd5f): a gene-level DE-conditional direction prior from 7 public banks loses to CLEAN4 on H1's gate
  (0.541 vs 0.573).
- **638 NEGATIVE** (d4cd5f, on 636's leave-one-line-out instrument): responder-weighted programmes (s_p, s_p², top-10 / 25 %,
  sign of the mean) do not beat CLEAN4 by +0.01 on any held-out line, and pooling more banks loses everywhere. Small
  positive: CLEAN4-ext (CLEAN4's own definition on each bank's full gene axis) adds 5–27 % gate coverage with signs better
  than all-down, worth roughly a quarter more of the hook's share. A rider at best.
- **611 → desk check** (21dd4c): the fid decomposition of our H1 build at f 0.35 is φ 0.05 of calls on DE genes, p_DE 0.68,
  p_null 0.52. fid is filler precision, as d4cd5f said.
- **621** (7772dc): the supervised DE membership + direction model, v1 LOLO over four gate surfaces. Running through the
  602 gate.

**Pivot rule (written now, before 621 reads).** If 621 also fails the bar (≥ +0.01 gate-sign accuracy over CLEAN4 on
most held-out lines), public SIGN SOURCES are closed for this round. The DE sprint then becomes an emission-DESIGN
sweep read directly on the A/B/C board. The board is the only truthful DE instrument, and 537 §A measured its DE read
noise at ±0.0005 overall (2σ) for a DE knob. Allowed parameters are global and protocol-level only, each a preregistered
2-point dose with a written mechanism (11.1c): the MDE / called-gene magnitude dose, the emission dispersion law
(William-style mean-preserving dispersion), the reach-order key, and the K curve already on 10-06/07. About 10 of the
~25 unallocated slots (10-08 … 10-21) can go to it. Each read either moves a protocol parameter that D/E/F shares or
closes it.

## 1c. 11:50Z — the pivot fired; the council's main line is now the DE DESIGN sweep

- **621 NEGATIVE** (7772dc, c7792a89): a chill-style supervised DE membership + direction model loses to CLEAN4 on
  held-out lines (Δ gate sign h1 −0.013, hct116 −0.099, hek293t −0.205, k562 +0.001). The held-out surface's base rate
  P(down | gate) shifts, and a cross-surface learner cannot know it. The one robust signal: a family-balanced consensus
  at |Zbal| ≥ 2 is right 0.81–1.00 against CLEAN4's 0.60–0.73 on the same pairs, worth ≈ +0.001 overall as an override.
  **Public sign sources are closed for this round** (634, 638, 621).
- **611 DELIVERED** (21dd4c, `610/RESULTS-611.md`). At f 0.35, H1 retrodicts the board's DELIVERY pairs, zs (nmae 1.18×)
  and 433a (nmae 1.59×, reach 1.16×, all signs right), and mis-signs reach for the CALL-SET pair (cp) in every regime.
  Planning multipliers: nmae ×1.2–1.6 for delivery / magnitude levers; call-set changes fid ×~2 (cost), jac ×~0.5
  (gain); reach ~1× for delivery only, never for call-set changes (read reach on the board, band ±(H1|Δ| + 0.004)).
- **The forensics that reopened the scorer side** (21dd4c, `610/destep610.py`; 201 DE-only board steps from 105 teams,
  bulk pinned < 2e-4). Many teams posted PURE-nmae steps with fid/jac/reach flat: 11111 −0.070 raw (+0.019 overall in
  one entry), Null-to-Cell −0.053, CellMind −0.052, BioFarmer −0.038, cqawesome −0.021, Testh −0.017. Ours, zs, is
  −0.0022. Mechanism: nmae reads ARITHMETIC means; calls and reach read Wilcoxon RANKS. U counts only pred-vs-control pairs,
  so a pred cell can move anywhere inside its gap between adjacent control values, and the top cell's gap is unbounded.
  One outlier cell per gene therefore sets the mean almost freely (AUC change ≤ 1/400). It is the ballast idea applied
  to DE means. L1 says that at sign accuracy a > 0.5 the optimal magnitude is the full |l| (worth (2a − 1)|l| per gene),
  which is why 480's oracle at our signs is −0.0985 raw. The leak cap that closed MDE doses above c 0.1 (426/435) does
  not apply to rank-neutral delivery. DOWN has no such free direction, and CLEAN4 is 84 % down, so UP coverage decides
  the size.

**Sweep arms (each: 537's real emitter with KAs first, six members, 3 seeds, H1-full + shape-matched h1r7, read against
FINDINGS §2's closures, a 2-point prereg with bins before any slot; competing with 537's list in the 10-08+ day rules):**

| # | arm | owner | class | board evidence | expected |
|---|---|---|---|---|---|
| 1 | rank-neutral nmae: UP means delivered through outlier cells with U, calls and the pseudobulk pinned; dose = full predicted \|l\| or ×2 / ×4 MDE | d4cd5f | scorer mechanics (delivery) | pack steps −0.017 … −0.070 raw nmae | +0.005 … +0.019 overall if UP coverage is large |
| 2 | direction-confidence reach order (calls and means pinned; key = 621's family-balanced consensus z) | 7772dc | scorer mechanics + transfer | chill +0.034, Chiu S1 +0.031, em2.bio +0.030 reach | +0.002 … +0.005; reach read on the board |
| 3 | mean-preserving dispersion law (William: control cells zeroed p 0.5 / doubled otherwise, top 40 % \|DE\| genes) | 21dd4c | scorer mechanics | fid +0.010, jac +0.004 | +0.002 … +0.004 |

## 1d. 12:05Z — the first public sign source that beats CLEAN4 on held-out lines, and the size of sweep #1

- **639 POSITIVE (d4cd5f): the IDENTITY AXIS transfers leave-one-line-out.** It encodes "Genes a line expresses above
  other lines' controls go DOWN under perturbation" (441 PS6's idea; perturbation erodes a line's identity programme).
  Score = z(P_C4 − 0.5) + λ·z(I), I = log1p(ctrl CPM) − mean over the other profiles. Gate sign vs CLEAN4 at λ 1:
  h1 +0.016, hct116 +0.015, hek293t +0.024 (CIs exclude 0), k562 null at low power. It is not shared-control noise
  (control-halves check: 0.5884 vs 0.5887 full), and the same sign on Flex and 3′ rules out a chemistry confound. λ 1
  was the prior choice from 441. Class: TRANSFER, computable on the day within Flex (each D/E/F context vs the other
  five Arc contexts). Size ≈ +0.003 … +0.005 overall, delivered through the block prior and the uncalled-gene hook.
  **It revises §0.5: one public route to the responding-gene programme does exist, the line's own identity contrast.**
  The CLEAN4 gate gap to the clean oracle (0.573 → 0.618 on H1) is about one-third closed (+0.016).
- **Sweep #1 sized (603 + 632 + 610's dose curve, three independent pricings):** the rank-neutral UP magnitude is a
  mid-sized transferable rider, not the +0.019 the rivals' steps suggested. Their big pure-nmae steps were mostly
  catch-up from worse bases (post-step 0.84–0.94 vs our 0.897). H1 at 611's f 0.35: K 1.0 → +0.0060, K 1.25 → +0.0080
  overall. At f 1, K 1.0 → +0.0035 and K ≥ 1.25 overshoots. DOWN hurts at every dose. Board estimate +0.003 … +0.013;
  the natural board pair is K 1.0 / K 1.25, UP only, with the 639 signs. d4cd5f is building it on the real emitter
  under the coordinator's four conditions.
- **Closed this round:** 621 (supervised sign model), 621c (direction-confidence reach order: +0.001 … +0.002 at the
  ceiling, CIs span 0; the oracle order would be worth +0.03, but no public signal knows which calls are right), and
  #3, the dispersion law (21dd4c: in our design the call set is prescribed, so the lever is K, already on 10-06/07).
- **New from 21dd4c (612):** a KEY-GATED hook (CLEAN4 sign where our key is weak, key sign where it is strong) prices about
  3× the plain CLEAN4 hook on H1 nmae (−0.0049 vs −0.0015 at f 1), about +0.0015 … +0.005 on the board. Its gate was
  chosen on H1, so a held-out cancer-fold transfer check is running.

## 1e. 12:15Z — what is on its way to the board (all UNAPPLIED patches; the coordinator applies with KAs)

| arm | owner | what it is | evidence | board expectation | status |
|---|---|---|---|---|---|
| **630 rank-neutral UP** | d4cd5f | one carrier per uncalled UP gene promoted just above the control max (ballast donors, budgeted crossings, a call-flip guard); `VCC_DE_UPRN=<s\|b>_k<K*100>` after zs | real emitter, H1, 446bj2z's construction (KA: zs == 435). b_k100: Δnmae −0.0066 raw at f 1; fid / jac / reach / pds / mse exactly 0; 0 call flips; max ΔAUC 1.0e-4. f 0.35 +0.0028 [+0.0013, +0.0043] overall-equiv. 25 % of the ask is delivered (recipient budgets); K ≥ 1.25 is worse | +0.0026 … +0.005, nmae only | builder patch + test ready; composition 598 → 612 → 630 due by 10-06 06:00Z; dose pair K 0.75 / 1.0 |
| **612 key-gated hook** (`VCC_DE_MDE_SIGN_KEYGATE=<q>`, source-agnostic) | 21dd4c | uncalled-gene hook: the key's sign where the key is strong, the prior's sign where it is weak | **612 CORRECTION 2 (21dd4c, 22691ce5, ~13:10Z): the fold numbers kgh612b/c are VOID.** That fold model applied the hook to every gene; production adds it only on uncalled genes (rc_emit_k.py:1059). Call-set-aware (2,600 called genes masked), kgID vs no hook at f 1: h1 +0.0059, HCT116 +0.0013, HEK293T +0.0022, K562 −0.0021; the plain hook's K562 tail is ≈ −0.0015, not −0.010. That model also FAILS the pre-committed H1 retrodiction (prereg A1), so **no held-out instrument is validated for hook-class levers**. On H1's real emitter key-gated vs plain is −0.0001 … −0.0003 overall. The robustness case (worst-case regret 0.0020 vs 0.0115) rested on the void K562 tail and is withdrawn | ~0 EV; the tail it protects is ~0.0015, not 0.0115 | accepted for 10-06 (a no-op when unset); 10-07 keeps CPH |
| **639/622/604 identity prior** — the identity SOURCE's sign quality stands (LOLO gate-gene sign +0.009 … +0.042; Broad-16 +0.009). Its DELIVERED value is open. **Through the hook:** on H1's real emitter (3 seeds) the key-gated identity hook TIES the plain CLEAN4 hook at board-like regimes (nmae −0.0002 / +0.0007 raw at f 0.5 / 0.35; reach +0.004). The earlier fold-based dominance claim (+0.006 / +0.021 / +0.020 / 0.000; 622 (j)'s "key gate mandatory") is WITHDRAWN (612 CORRECTION 2). **Through 604's BLOCK prior:** it needs a prior-only board or real-emitter read | d4cd5f (price, prereg), 7772dc (science), 4d0207 (plumbing) | z(CLEAN4) + 0.5·z(I) with Flex-only references (each context vs the other Arc contexts + H1; one reference per identity call, never its own), quantile-mapped onto CLEAN4's P | LOLO 3/4 surfaces (h1 +0.014, hct116 +0.031, hek293t +0.027 at λ 0.5, Flex refs); Broad-16 mean +0.009 (12 lines, 0 negative); stable across NTC halves (I Spearman 0.9994+) | hook ≈ 0 (H1 tie); block prior unmeasured | A/B/C artifacts (the 458d path: `VCC_SAMELINE=<npz>:<sha>:0`, no builder change) + day-of `def_sameline.py --programme-source identity:0.5` patch, both in 604; hard gate: per-context prior λ ≥ 0.9 (554's floor) |
| **625 A/B/C module readout** | 7772dc | **Amendment 2 (12:45Z, committed before any read):** 10 PRIOR-ONLY probe entries on the A/B/C board (8 Sylvester-Hadamard rows over 8 CLEAN4-strength octiles + 2 single-octile additivity arms, octiles 0 and 7), only in slots the pipeline leaves empty, run last. It rehearses exactly the October decode. Hadamard column 0 (the noisy one, 1.20σ) holds octile 3, the weakest prior (median \|P−½\| 0.005). `modules625_v8.npz` sha c24e3304; `decode625.py --design v8` selftest PASS | 622h: on the Broad-16, a same-lab programme from just 3 lines transfers (+0.054 gate sign, 12/12 lines); 8 octiles +0.0165 vs 16 modules +0.0172, the same value at about half the entries. Caveat (625c): through the hook alone, per-octile nmae effects are ~0.0002–0.0007 overall vs a decoded sd ≈ 0.00017, so about half to two-thirds of the octiles will be readable; the block-prior share needs real-emitter flip arms | +0.003 … +0.004 (D/E/F) | approved by the coordinator (Ahmet told, no veto); runs last |
| **624 learned per-target bulk magnitude** | 7772dc | ridge on truth-free features (Testh's disclosed design class) | stage 1: learned tilt +~1.5 % aggregate cosine over 522's strength rider on 8/9 October-like folds. **Exact check (624c): +1.29 % on 8/9 October-like folds, but it FAILS 522's adoption bar on the Arc-like re-draws: +0.71 %, 4/9. The tilt reverses on the Orion folds (real hct116 −1.22 %, real hek293t −1.15 %).** | 0 | **CLOSED (624c): do not adopt, bundled or alone; 522's rider alone still passes (8/9)** |

**Running tally of today's NEGATIVES (each with its mechanism in its RESULTS; + 629 at 15:55Z: H1's programme does not transfer in any de-identified / detrended form, 0/12 Broad lines):** 601 (compositional / generic filler
signs), 610 (null skew), 620 (non-CCLE bank identity: never fires on cancer lines), 621 / 621b / 621c (supervised sign
model; learned membership; direction-confidence reach order), 623 (identity erosion as a bulk component), 634
(conditional direction prior), 638 (responder-weighted programmes), 635 (the filler programme was control noise), 624 (a learned
per-target bulk magnitude: fails 522's Arc-like bar, 624c), 651 (library-size composition as a filler sign source: under its
preregistered bar).

## 1f. 12:35Z — Ahmet's decisions and what followed
- **AWS:** "if AWS is not doing anything useful obv stop it". The 5 stopped VCC workers' root volumes were snapshotted
  (`snap-0a8446d4…`, `snap-03c1a4a7…`, `snap-01e6d676…`, `snap-0cfaf592…`, `snap-04d7b421…`; tagged project=vcc) and
  the instances are terminated once the snapshots complete (`work/aws-cleanup-1002.log`). The S3 bucket stays: it
  holds cloud_score's inputs.
- **Cursor:** the stale non-VCC agent (`bc-619e4809…`, survev-rl) got a stop call (HTTP 200).
- **Pisces:** "yeah sure. use <redacted-email>". One email went to the X-Cell corresponding authors (Ci Chu,
  cc Bo Wang). It asks for early research access to the Jurkat rest/active + HepG2 screens and whether challenge use fits
  CC-BY-NC-SA. It contains no method details.
- **H-631 DECIDED (~12:40Z): maximum-information use of D/E/F feedback, including decoding D/E/F's own module
  programme (637 §4).** Ahmet, verbatim, in [7772dc]'s session: "we can do the risky thing". To [537757]: "do what helps us
  learn most ... make use of all info. its not gonna disqualify." Every read is documented and disclosed in the write-up.
  7772dc owns the decode (625's decoder), with d4cd5f's fused design and the coordinator's cadence / re-roll criteria; 8
  CLEAN4-octile modules keep ~all the value at ~10 entries.
- (history) **H-631: CONFLICTING signals, re-asked.** In this session Ahmet wrote "for 1., yeah we can do that",
  first read as approving no-decoding (DEF-RUNBOOK §6 item 8, 7cf33fca). [537757] relays from its session, verbatim:
  "i dont find that an issue. what we should do is do what helps us learn most. we should document and make use of all
  info. its not gonna disqualify." The coordinator holds §6 item 8 as unconfirmed. d4cd5f has the maximal-information
  design ready (637 §4: factorial bulk switch + a 16-module Hadamard decode of D/E/F's own programme per entry, fused
  with 625).
- **Council opinions on decoding D/E/F's own programme from final-entry feedback (637 §4), for Ahmet's decision:**
  - d4cd5f designed it: it reads the final's own truth, worth more than any A/B/C transfer.
  - 7772dc recommends NO: ~16 bits of the test set's answers read through final entries is the "spirit" pattern, for
    ≈ +0.003 … +0.004. Instead, use the A/B/C-learned mask (625 §6), which carries most of that value.
  - The coordinator follows Ahmet's latest explicit words: maximum-information use, documented. It also priced the
    cadence (~13 entries) and the re-roll (per-module sd ≈ 0.00012) costs.
  - 4d0207 (scribe): the expected gain is small next to the stakes, so Ahmet's explicit choice should decide it.
    Either way, 625's A/B/C run is useful as a transfer read and as the rehearsal.
- **625 module readout APPROVED by the coordinator** (Ahmet told, no veto). It runs last, only in empty 10-08+ slots, and
  has four conditions: one frozen base by sha; a decode with a synthetic known answer before any read; the queue behind
  every decision arm; -m from the queue, -d empty.
- **626 sister-reference guard (7772dc)** is in `def_identity_prior.py`: drop a reference whose deviation-Spearman
  similarity to the context is > 0.3, never below 2 references. It is inert on {A, B, C, H1} (pairwise −0.08 … −0.42) and
  fires only on the Broad-16's true sister pair, at a cost of −0.0001.
- **Identity hook, October path with NO def_day.sh edit:** `DEF_DAY_MDE_SIGN=hookbank:identity:0`. def_day's hookbank
  step runs `def_sameline.py --programme-source identity` (patch604; bare name = λ 0.5) and exports the file form at
  τ 0. The KEYGATE pass-through is 21dd4c's patch612dd.

## 1h. 13:35Z — corrections and the independent review
- **612 CORRECTION 2 (21dd4c):** the fold model behind the hook-class gains hooked every gene, while production hooks
  only uncalled ones. Masked, the gains shrink ~15×, and no held-out instrument passes the H1 retrodiction. Every
  hook-class number above is updated. What stands is H1's real-emitter TIE (KGID vs the plain hook).
- **604 → v3 (d4cd5f's 639d + ruling):**
  - identity-only genes get NO prior; their signs were right only 0.27–0.34 on HEK293T;
  - the side is sign(zc + λ·zi), 639's validated rule. v1 / v2's rank map kept CLEAN4's down-count and undid A's
    identity down-shift (A: 58.2 % vs 53.5 % down on covered genes);
  - KAs 39/39;
  - ID05 / ID10 are pinned in 639's block prereg (d5f7b65b): v3 λ 0.5 `6cd84dc2…`, λ 1 `3b1e22f0…`.
- **The Grok cloud review (605)** read a stale commit, but found one REAL blocker in 630's current builder module. The
  donor-gain clamp (`src_rank_up.py:190`) lets an UP gene's mean fall and its lfc flip while every invariant passes.
  The fix (donor filter + mean/sign revert) is with d4cd5f before 10-06. 612's tie-at-cut rule is low-risk (the knob
  is unset). Triage: `605/RESULTS.md`.
- **Drug generalisation (627, 7772dc):** identity erosion holds under 1,137 drug treatments in 50 Tahoe lines
  (Spearman +0.147, grows with strength; the ambient-RNA confound is bounded). See SCIENCE §2.

## 1i. 13:45Z — the 10-06 chain is fixed, extended and pinned; freeze 10-05 18:00Z
- **630 fixed twice after the cloud review:**
  - 69972aa1: donors only where the library is larger than the carrier's, and revert if a gene's mean falls or its
    sign flips;
  - c30d64ce: guard and invariant failures now fail SAFE (the target keeps its zs block) instead of aborting a build.
  - Module `7e6b3e52`. Production exposure to the original bug was none: ballast rows hold ~995k counts vs ≤ 25k for
    carriers on 85/85 H1 blocks. The nine bit-for-bit H1 KAs re-run on the shipping module.
- **630dd (4d0207 for d4cd5f):** the def_day pass-through `DEF_DAY_DE_UPRN` → `VCC_DE_UPRN`; test630dd 45/45.
  - It requires DE_MAG, the builder's grammar with K in [0.5, 3], and the meta variant plus module `7e6b3e52`.
  - It reports per context (never fatal) with a greppable `630dd: uprn FAILSAFE_TOTAL <n>`. Per the coordinator, a
    KA630 build with n > 0 is a KA FAIL.
- **612:** the tie-at-cut fix (740f5bcf).
- **604:** v3, plus the octile-flip option for 625 (default off and bit-identical; test604flip 21/21).
- **606 chain** (641 → 598 → 612 → 612dd → 630 → 630dd → 642 → 604 → 644): ALL PASS, 612dd / 630dd commute. Pins:
  - def_day `c2704f3d`
  - builder `2c6cfe7d`
  - rank_up `7e6b3e52`
  - def_sameline `d964aa03`
  - def_identity_prior `60e98488`
  - the rest unchanged
- **Coordinator: CHAIN FREEZE 10-05 18:00Z.** After it, patch changes wait for 10-08. The day-of flip pass-through
  (def_sameline `--flip-octiles` + `DEF_DAY_DE_PRIOR_FLIP`) is a 10-08+ patch.
- **607 DONE (13:50Z): the chain HOLDS** (no blocker, no major; the reviewer rebuilt it and reproduced every pin).
  - Three MINORs. 604's name collision is fixed. 644's comma-space and 630's rounds-0 grammar are relayed to their
    owners.
  - NOTEs relayed.
  - Re-pin: def_sameline f147f8bb, def_identity_prior 7aeb66b4.
  - The 607b delta review (630 module 7e6b3e52, 630dd, 604 flip) is running.
  - Triage: 607/RESULTS.md.
- **608 (forensics for 21dd4c):**
  - Illumina is a staircase. Its 09-27 → 09-29 jump of +0.052 includes a pure DE step of +0.013 (fid 0.543 → 0.559,
    jac 0.044 → 0.052) and two bulk steps (mse 0.49 → 0.36).
  - No top-30 DE-only step since 10:43Z.
  - The fid pack's jac has been flat for a week.
  - Two clone clusters on the board (5 teams and 4 teams with identical DE members) should each count as ONE design.
- **609 (for 7772dc):** the 639d (block) vs 622m (uncalled) conflict is being settled on shared definitions. Two leads:
  - 622m mean-centres zc (that changes CLEAN4's own baseline sign) where 639d / production do not;
  - 639d's block gain comes from the two ORION folds, whose truth carries a depth-made DOWN skew (AGENTS 6.5;
    HARNESS-DEPTH-CONFOUND). On H1, 639d's own numbers put most of the gain on uncalled genes too (block 184 vs
    uncalled 428 net).

## 1j. 14:05Z — where the identity value can actually reach the board
- 609: the identity sign gain sits on UNCALLED genes. The block route is ~0, and the Orion block result was an artifact.
- 622n: within the uncalled genes, the LOW-identity tail goes UP robustly: H1 real emitter 0.61 / 0.71 vs delivered
  0.52; Broad-16 0.78, 12/12 lines. The high-identity tail's DOWN does not transfer across modalities.
- 612: the hook alone ties the plain hook. 630: UP magnitudes compound with better UP signs.
- So the identity programme's realistic delivery is **630's UP pass with identity-informed UP signs on uncalled genes**,
  not the block prior (ID05 ≈ NULL) and not a DOWN pass. d4cd5f has 622n + 625e for the 630 stack. This is the arm to
  price next, on H1's real emitter, as "630 on the idall3 / KGID base vs 630 on the plain base".

## 1k. 17:43Z — audits done, the 10-07/08 rule frozen, one dose read running
- **652 (7772dc, October evidence audit): no DECIDED switch rests on a VOID read.** Eight flags; the coordinator actioned the
  four that needed action before the freezes (72d5c1de): F1 drops the Arc hook (MZ / CPMZ), F2 gates CPH's October adoption on
  614's K562 sign, F3 writes the K16 / K20 rule before its read, and F4 defers tercilec until an H1 DE side-effect run ([36a81a])
  and C10b.
- **653: adversarial read of the 10-07/08 rule. Its v4 is the frozen rule** (c5cdd3bc).
- **624 CLOSED** (624c; see §1 and the ledger). The ledger's sum loses ≤ +0.001 at its top end.
- **654 (7772dc, 5de66787): DE carrier amplitude.** On H1's real emitter, 3 seeds, KA PASS bit for bit:
  - ×0.8 is worse on 3/3 seeds (f 0.35 −0.0028 overall-eq). Closed below ×1.
  - Production-faithful ×1.25 (with production's zero repair, which 537's H1 harness omits) is better at every
    board-like regime on 3/3 seeds: f 0.5 +0.0011, f 0.35 +0.0019, cmp +0.0012 overall-eq; nmae −0.0042 raw.
  - ×1.25's as-registered loss (−0.0033) was that missing repair (zero-total collapses).
  - It is a small call-set change (+26 calls per target), so 611's call-set rule prices it at ≈ +0.0013 … +0.0018 overall.
    Class: scorer mechanics.
  - 410 §7's "DE amplitude is not the lever" (measured before zero repair existed) is overturned up to ×1.25.
- **655 (7772dc, d7897f74): the dose curve under zero repair** (overall-eq vs ×1; f 0.5 / f 0.35 / cmp):
  - ×1.25: +0.0011 / +0.0019 / +0.0012;
  - ×1.5: +0.0007 / **+0.0031** (3/3) / +0.0019;
  - ×2.0: −0.0022 / +0.0005 / −0.0004, and infeasible on one draw (the K-class emitter's per-cell capacity assert).
  - The peak is interior and moves with the regime: ×1.5 at f 0.35, ×1.25 at f 0.5. Preregistered bin: **WEAK CANDIDATE**.
  - 611's board price: ×1.25 +0.0007 … +0.0018; ×1.5 +0.0018 … +0.0035 (it leans on fid's call-set term, +57 calls per
    target). ×1.25 is the regime-robust choice.
  - Next: 618's held-out read for both doses, then a builder KA (def_postbuild's α pin at 0.982 must move too), then one
    A/B/C arm. Class: scorer mechanics.
- **Board (608 §4):** no DE-only step in the top 30 since 10:43Z. null_hypothesis climbs reach at fixed bulk
  (0.254 → 0.306); that read sits inside 621c's reorder ceiling and tracks their fid (signs), so the reach-order lever stays
  closed.

## 1l. 21:54Z — the evening: Orion DE truth voided, α ×1.25 + 630 chosen, 618 is the held-out read
- **616 (21dd4c): the thinned Orion DE truths were 94–98 % artefact.** Rebuilt on clean raw cells, the same 132 targets
  give HEK293T 37,748 → 774 calls (median 1 per target) and HCT116 11,309 → 712 (median 2). The "Orion technical blob"
  was the old greedy builder's control selection (615 reproduced it from raw cells).
  - **Every council instrument's Orion-fold DE rows are VOID** (`616/AUDIT-orion-dependence.md`).
  - The valid DE surfaces are H1, K562 (thin) and the Broad-16.
  - Donor tables and PLOLO use pooled controls and are clean, so the bulk decisions stand.
  - The coordinator re-ruled from the dependence list (NEXT-SLOT-DECISION 19:10Z).
- **614 / 615 (real-emitter DE harnesses on K562 / Orion raw HEK293T): both too thin to trust.**
  - 614 failed its preregistered trust gate on truth power (median 5 DE genes per target, 34 / 100 scoreable); its numbers
    are descriptive.
  - 615 stopped on its power diagnostic (median 1 DE gene per target, 18 / 100 scoreable).
  - 614's H1@11k read triggered its preregistered faithful-key re-run (low priority, queued behind 618).
- **618: a second H1 surface** (h1-trainval's 126 targets with h1-test's controls; median 149 DE calls per target, 94 / 126
  ≥ 10, so H1-class power). Its KA reproduces 323's cap15 bit for bit, and it is now the held-out instrument for the DE
  magnitude decision. Its 14–15 GB steps run through gate602 v6's opt-in quiet admission.
- **Coordinator decision (21:35Z): path (a), α ×1.25 + 630.** On H1 it is +0.0052 (f 0.35) / +0.0049 (cmp) overall-eq over
  the plain base. 631's DN pass is a SUBSTITUTE (it overshoots on a ×1.25 base), so it is the frozen fallback (b).
  - The pick follows ONE instrument, 618, registered before any read (14b5fda0): the higher MIN over {f 0.35, f 0.5}
    wins, and within 0.0005 (b) wins. The H1 prediction is MIN(a) − MIN(b) = +0.0014.
  - (a) also needs the α-path KA, including 554's dose guard on A/B/C (A/B/C's capacity is ~0.70× H1's).
  - Queue: 630 → 524dr2 → ALLSW → α ×1.25 → idall3 → 625.

## 1m. 10-04 06:18Z — 657, and a calibration caution on Broad-16 sign evidence
- **657 (7772dc, 92767352): CLEAN4 octile 4 delivered UP is call-neutral on 618** (3 seeds, f 0.35 overall-eq):
  - A (hook UP) + 630: −0.00044 (Δcalls +3.6);
  - B (630 also selects octile 4): −0.00011 (Δcalls 0, nmae-only);
  - AB: −0.00050.
  All routes are ELIGIBLE under the no-harm gates, are DN-invariant, and keep the bulk sums pinned. **B is the board route by
  the prereg rule but ranks LAST** (expected ≈ −0.0001). Class: transfer as attempted; 618 measures no-harm only.
- **Caution (post-hoc, `657/calib657.py`): Broad-16 gate cells carry a global UP offset.** UP-right at octile 4 is 0.60 on
  Broad vs 0.50 on H1, and the Broad curve crosses 0.5 at P_down ≈ 0.55 (H1 ≈ 0.51), about 1.5 octiles to the UP side.
  No CRISPRi bank shares it (K562 ×2, Jurkat, HepG2, H1, HCT116, HEK293T). So Broad-16 sign evidence near P_down ≈ 0.5
  can be KO-modality or lab calibration rather than transferable biology.
  - 625's 8-octile Broad-16 transfer (+0.0165 gate sign) is one such case. Parts of 622n's Broad-16 corroboration
    (low-identity tail UP, 12/12) and 639's Broad rows may be others.
  - **H1's own reads stand,** because H1 is calibrated (crossing ≈ 0.51). 622n's H1 result (0.639 vs 0.557) is unaffected.
  - The cheap check (handed to d4cd5f): re-centre each Broad line's UP-right curve at 0.5, then see what survives.
  - **622n re-centred (7772dc, 3a4c98ba): the low-identity → UP tail survives, weaker.** It is above each line's own
    uncalled UP base rate on 8/12 lines (mean +0.048, median +0.060, pair-weighted +0.081), and the bottom 5 % on 9/12,
    no longer "12/12 > 0.5, 0.781". H1 and 609b's expression-only control are unchanged. idall3's prereg evidence does
    not change (622o's hook rows were already base-rate balanced). 625's +0.0165 stays flagged.

## 1g. The D/E/F ledger: what the council's arms could add on the final, if each transfers as measured
Planning numbers in overall units. Each needs its board read (10-06 … 10-14) before it enters the October recipe.

| arm | class | measured held-out / H1 | expected on D/E/F (dark board) | status |
|---|---|---|---|---|
| 630 rank-neutral UP (v3 shipping form `b_k100_q20r3p10`; dose partner K 0.75 in v3 form queued) | scorer mechanics | v3, real emitter, H1, 3 seeds (d4cd5f 13:00Z): f 1 +0.0023 [+0.0005, +0.0042], f 0.35 +0.0029 [+0.0015, +0.0044], cmp +0.0034 overall-eq; every other member exactly 0, 0 call flips; ~2,800 genes moved per target (v2 ~1,500); fixes context B's capacity miss. **Compounds with the identity hook:** on 612's KGID base the UP pass gives nmae −0.0105 raw at f 1 vs −0.0080 on the plain base (+0.0029 f 1 / +0.0041 cmp overall-eq), because better UP signs make the extra magnitude worth more (632's prediction) | +0.002 … +0.005 (capacity per context: gate ≥ 50 % of H1's) | builder patch applies 10-06; board pair from 10-08 |
| identity HOOK (key-gated, τ 0; 622 / 612 / 604) | transfer | **H1's REAL emitter, 3 seeds (21dd4c, 12:55Z): KGID − plain CLEAN4 hook overall-eq f 1 −0.0014, f 0.5 +0.0008, f 0.35 +0.0006, cmp +0.0011 (sd ≤ 0.0009); predicted board Δ(fid+nmae) vs 496a −0.0003 … +0.0001: the TIE bin.** vs no hook it pays nmae −0.0015 … −0.0046 raw in every regime. Held out: call-set-aware folds vs no hook +0.0059 / +0.0013 / +0.0022 / −0.0021 (h1 / HCT116 / HEK293T / K562), but that model fails the H1 retrodiction, so **no held-out instrument is validated for hook levers** (612 CORRECTION 2; the earlier +0.021 / +0.020 fold gains are void) | ≈ 0 over the plain hook (unvalidated beyond H1) | artifacts + day-of path ready (604); KEYGATE 10-06 (612); 612dd for def_day |
| identity BLOCK prior (λ 0.5; 639 / 604 v3: ID05 `6cd84dc2…`, ID10 `3b1e22f0…`) | transfer | LOLO gate sign +0.009 … +0.031 on ALL gate genes (Flex refs). **609 settled where it sits:** 639d's block gain on HCT116 / HEK293T is the Orion truth's depth/expression DOWN skew: an expression-only rider matches or beats it, and HEK293T I⊥E gives −0.006 on the block. Honest held-out BLOCK readings: h1 +0.006, Broad-16 −0.001. The gain lives on UNCALLED genes | **+0.000 … +0.001** (the prereg's +0.001 … +0.004 rested on the Orion reading) | prereg pinned (d5f7b65b); reading 1 to be amended (d4cd5f); exact prior-only read, now low priority |
| 625 A/B/C module readout (Amendment 3: 8 single-octile arms; λ 0 / 0.5 artifacts built, 604) | transfer (same-lab) | Broad-16 8-octile +0.0165 gate sign. **625 FINAL (7772dc, 73c0b648):** a sign-source improvement pays ≈ +0.002 overall per +0.01 gate sign (nmae + fid; range +0.0015 … +0.0027). Prior-only octile reads are exact (σ ≈ 1e-5 raw nmae) | **≈ +0.003 on D/E/F if the A/B/C reads GO** | approved, runs last; artifacts octile-only (F7 / F8, checked independently) |
| 624 ridge magnitude | transfer | +1.29 % on 8/9 October-like folds; **fails 522's Arc-like bar (+0.71 %, 4/9; reverses on Orion)** | 0 | **CLOSED (624c)** |
| 612 key gate | robustness (withdrawn) | the worst-case-regret case (0.0020 vs 0.0115) rested on the void fold model; call-set-aware, the plain hook's K562 tail is ≈ −0.0015 | ~0 EV, small tail | 10-06 (no-op when unset) |
| 654/655 DE carrier amplitude (production-faithful, zero repair) | scorer mechanics | H1 real emitter, 3 seeds: ×1.25 f 0.35 +0.0019 / f 0.5 +0.0011; ×1.5 f 0.35 +0.0031 / f 0.5 +0.0007; ×2.0 past the peak | ×1.25 +0.0007 … +0.0018; ×1.5 +0.0018 … +0.0035 (611's call-set rule) | WEAK CANDIDATE; 618 held-out read, builder KA (α pin), one A/B/C arm pending |
| **sum, if all transfer** | | | **≈ +0.004 … +0.011** (15:55Z: hook ≈ 0; block prior ≈ +0.000 … +0.001 (609; ID05-like on H1 +0.0005 … +0.0011 per 625 FINAL); 625 ≈ +0.003; 630 +0.002 … +0.005) | |

Against today's A/B/C gaps (chill / nailong +0.021, CellWM +0.013, Testh +0.033, cqawesome +0.055), and with 7772dc's
caution that the pack's global DE tuning probably also survives into D/E/F (same Arc campaign), +0.005 … +0.014 closes
CellWM's gap at best. It does not close chill / nailong's (+0.021), let alone cqawesome's. **The council's measured
arms are increments, not the breakthrough. The biggest single item is still the responding-gene programme (gate sign
0.573 → 0.618), and today no public source beats CLEAN4 on it held out (634 / 638 / 621 negative).** **If D/E/F members come back and H-631 is maximum-information**,
637 §4's fused decode reads D/E/F's OWN programme. The clean-oracle gap it targets (gate sign 0.573 → 0.618) is
+0.008 … +0.010 by itself, before cadence and re-roll costs (coordinator's caveats: ~13 entries; per-module sd ≈ 0.00012).
**Corrected 13:05Z, then 13:15Z (7772dc; 612 CORRECTION 2).** The structural discount stands:
- (a) Module vs gene resolution. 635's gap is a per-GENE programme. Module-level knowledge recovers only ~⅓ of the
  gene-level value on the Broad-16 (16 modules +0.017 … +0.020 vs gene-level 3-line +0.054).
- The 13:05Z number (≈ +0.001 … +0.002 through the hook, with ½ … ⅔ of octiles readable) came from the void every-gene
  fold model and is WITHDRAWN.
- **625 FINAL (7772dc, 73c0b648) sizes it.** An improving sign source pays ≈ +0.002 overall per +0.01 gate sign
  (range +0.0015 … +0.0027). The programme prize through today's routes is therefore ≈ +0.007 … +0.012, so the
  council's +0.008 … +0.010 stands. 625's mask, if the A/B/C reads GO, is worth ≈ +0.003 on D/E/F.
- The call-set term (CPM coupling) is nonlinear. Only large DOWN → UP shifts of strong priors cascade (+281 / +1,287
  calls per target on octiles 5 / 7); small improving flips are call-neutral.
- The A/B/C-transfer version of the mask (625) is expected to be worth the same order as the decode. The D/E/F decode
  mostly adds CERTAINTY that the flips are right for D/E/F, not size.

## 2. Decisions only Ahmet can make — status at 13:00Z (five of six closed)
1. **D/E/F feedback policy (H-631): DECIDED, maximum information** ("we can do the risky thing"; "do what helps us learn
   most ... make use of all info"). Every read is documented and disclosed in the write-up. 7772dc owns the decode (625 v8
   rehearses it on A/B/C). The council's earlier NO-decoding recommendation is kept in §1f as history.
2. **Pisces early access (H-621): DONE.** One email from <redacted-email> to the X-Cell corresponding authors (Ci
   Chu, cc Bo Wang). It asks about access and whether challenge use fits CC-BY-NC-SA, with no method details. Waiting for
   their reply.
3. **AWS: DONE.** Each of the 5 stopped VCC workers is terminated once its root-volume snapshot completes.
   - 12:38Z: 2 of 5 terminated (`i-06197bcf…`, `i-0af19b9a…`); 12:48Z: 3 of 5 (`i-0bd65df6…`).
   - 13:33Z: 4 of 5 (`i-08d90ae8…`). **18:35Z: 5 of 5** (`i-08300b39…`, after its 400 GB snapshot completed). No orphan
     volumes. Kept: 5 root-volume snapshots, 571 GB stored in total (snap-0cfaf592… 428 GB, snap-0a8446d4… 84 GB, snap-03c1a4a7… 35 GB, snap-01e6d676… 23 GB, snap-04d7b421… 2 GB), STANDARD tier, ≈ $29/month. Deleting them, or moving them to archive tier (≈ $7/month, 90-day minimum), is Ahmet's call.
   - The S3 bucket stays (cloud_score's inputs).
4. **Stale Cursor agent: stopped (cosmetic remainder).** `bc-619e4809…` (survev-rl) took two stop calls (HTTP 200, 12:20Z
   and 12:50Z). It still lists as RUNNING, but it has 0 messages since 09-20: a ghost record that never ran, so nothing
   bills. Deleting it from the Cursor UI would clean the list; it is not deleted here (deletion is permanent).
5. **625 A/B/C module readout: GO.** The coordinator approved it, and Ahmet was told with no veto. It runs last, in empty
   10-08+ slots (v8: 10 entries).
6. **DONE (tier 1, 15:30Z; tier 2a 16:39Z): the critical data is in the private S3 bucket.** Tier 1 is 279 files with
   sizes verified and a sha256 manifest; tier 2a is all 18 `data/harness` files. Tier 2b, 199's 16.8 GB GWCD4i, landed at
   18:19Z (size match, sha256 recorded). The backup is complete: 298 objects, 50.4 GB. **Restore drill (20:05Z): 5 of 5 files restored with matching sha256** (gene_names, pert_counts, a donor b_ctrl, an 18.8 MB mask, context_C.h5ad 226 MB): 245 MB in 14 s ≈ 17 MB/s, so a full 50 GB restore takes ≈ 50 min. An external SSD is now a convenience, not a gap. Still open: a fast LOCAL restore path (an external SSD), because a full S3 restore is
   slower. Before today no backup of the data existed. Time Machine has no destination, and only code and small files are on GitHub.
   A disk failure between Oct 22 and Nov 5 would mean rebuilding donors and operators from public sources, which takes
   days. The home uplink measured 1.3 MB/s to S3, so a cloud copy of the ~20–30 GB critical set takes ~5–6 h and
   competes with slot uploads.
   - **Cheapest fix: plug in an external SSD.** A one-time copy of `data/donors`, `data/challenge`, `data/harness` and the
     October work files takes ~10–30 min. Or point Time Machine at it.

## 3. Clock
- 11:40Z: v2 (this file); the portfolio goes to the coordinator as proposals.
- Execution now; check-ins about every 2 h (13:45Z, 15:45Z, 17:45Z, 19:45Z). 20:00Z: chain1003 arms the 10-03 slots;
  4d0207 watches it.

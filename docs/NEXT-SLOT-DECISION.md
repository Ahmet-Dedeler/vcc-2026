# Slot decision tree — written 2026-09-09 22:42Z, before either result exists

Supersedes the 2026-09-08 tree (316/317), which is fully resolved: parity lost (−0.00075), the operating-point
trade does not exist (317), and the cap arm won (+0.00094, entry `xuDelEik7zGWMXdNu0R4`, **0.21865974**).

## Tonight, 2026-09-09 UTC (reset 00:00Z)

**Slot 1 is spent**: `319d3-cap30.vcc` (cap −3.0) uploaded 22:36Z as "a moth in the lint trap".
Bins are frozen in `experiments/319-cap-variant/PREREGISTRATION-dose2.md`; Δ is against 0.21865974.
Capture the record the moment it publishes (`scripts/capture_board_entry.py <id> --package ... --sha256 6ddbf76d…`).

**RESOLVED 23:02Z: 319d3 scored 0.2183911, Δ −0.00027 = bin 2 (flat). Family closed; the −4.0 build was stopped, nothing packaged.**

**Slot 2 RE-DECIDED 23:10Z (Ahmet asked why it should expire with 54 min left).** An inventory of `submissions/`
against `SUBMITTED_LEDGER.json` found four never-scored packages. Three are dead by the record (229-stack4 priced
−0.006 "do not spend the slot"; 255 k=10 "dead on the realised factor"; 229-ka is a known-answer rebuild). The
fourth, **254 `cd4svd5+partial`**, is one variable on the scored `cd4svd5` entry, fully verified, drafted at
+0.0009…+0.0018, and the 09-05 tree said it "needs its own transfer reading" — never bought because the base moved.
Bins frozen in `experiments/254-partial-rows/PREREGISTRATION.md` §7 (commit 0a5eb220) and uploaded at 23:09:41Z as
"the intern watered the plastic plant". It is a mechanism reading on its own base, not a candidate: the displayed
row will fall below 319 and 319d3 and is not to be corrected. The −4.0 tree below is kept for the record.

**Slot 2 was NOT committed.** `319d45-cap45.vcc` (cap −4.5) is building speculatively
(`experiments/319-cap-variant/run319d45.sh`, ~60 min). The rule is in `PREREGISTRATION-dose3.md` §2:

| dose-2 result (S3 − S1) | slot 2 |
|---|---|
| **≥ +0.0005** (bin 1: cap axis open) | **upload −4.5** via `scripts/arm_submit_after.sh <d3-entry> 202609092300 submissions/319d45-cap45.vcc <sha> "the toaster filed a complaint" <log>`; bins in `PREREGISTRATION-dose3.md` §3. If prep overruns 00:00Z it becomes tomorrow's slot 1 — still correct. |
| **−0.0005 … +0.0005** (flat) | slot expires. Family closes at +0.0026 (dose-2 bin 2). Keep the −4.5 package, unsubmitted. |
| **≤ −0.0005 but > −0.00127** (softened turnover) | slot expires. Family closes at −1.5. |
| **≤ −0.00127** (cap falls as hard as plain) | slot expires. The up-weighting story for 319's gain is wrong; reopen the mechanism before any further arm. |

An expiring slot is not a reason to override a frozen reading (`docs/READ-FIRST.md` §1.0 correction).

## 2026-09-10 UTC — both slots spent (rewritten 04:20Z)
Slot 1: 319p = **0.2198190, +0.00116 vs 319, bin 1 — new best.** Slot 2: 319pA (Jurkat-calibrated context-A rows,
β=0.35) = 0.2185396, −0.00128 vs 319p, bin 3 — pds pays; closed at this construction.

## 2026-09-11 UTC — slot 1 spent, slot 2 expires (rewritten 20:10Z)
**Slot 1: 319p2 = 0.2201159 (`lqWNXX5YaYArdglu6Stn`), +0.000297 vs 319p, bin 2 (flat) — best-scored package, rank 15/893;
the partial-rows dose family is closed.** The 00:00Z waiter (pid 55426) died in the 2026-09-10 19:56Z reboot fifteen
minutes after arming; re-armed manually at 19:10Z, published 19:41Z. Record: `experiments/319-cap-variant/BOARD-RESULT-p2.md`.

**Slot 2: nothing to build that teaches.** Candidates examined and rejected today (all on disk):
per-target amplitude allocation for `mse` (H1 oracle ceiling −0.012 raw; `scratchpad/mse_alloc.py`), a lower whitening
gamma (170: V-vertex at 0.65), a smaller call set (192 §5.3 priced K at four n_real levels; the yield term makes K < n_real
catastrophic on the heavy tail), the H1-union dose (258), CD4 k=10 (255), and the down-first `reach` order
(`experiments/319-cap-variant/ORDER-ARM-NOTE.md`: −0.0135 raw reach on H1 for the emitter's own order; knob
`VCC_RC_ORDER` left inert). A STATE-blend arm needs STATE predictions for A/B/C that do not exist yet (383 is producing them).

## 2026-09-12 UTC (reset 00:00Z)
**Slot 1 candidate: a STATE-blend on context A (or all contexts) via the existing `VCC_CTX_A_PRIOR` `blend_rows` path**, if
383 delivers identity-checked STATE deltas and the H1 blend cosine is not negative; preregister bins vs 319p2. Otherwise
the best of 384 (PromoterAI-derived per-gene prior as a reach key) / 385 (per-target call-set size), each only if it
clears its own H1 gate. **Slot 2: dose 2 of whichever slot-1 arm lands in bin 1**, else expire.

## Standing constraints, do not re-derive
- Slots before Oct 22 are prize-irrelevant and expire nightly. Spend them on a preregistered arm that buys
  information. Never on rank; never re-submit a scored package.
- The board shows the LATEST entry, not the best. A losing arm drops the displayed row. Costless; not corrected.
- Capture the record immediately after every entry publishes; the API resolves only the current entry.

## 2026-09-12 UTC — both slots UNSPENT (nobody at the machine; the 2026-09-12 00:22Z reboot killed 383's STATE training)

## 2026-09-13 UTC (reset 00:00Z) — written 06:25Z, before any result exists
**Board 05:44Z (920 entries): 319p2 is 19th at 0.2201159.** Illumina's displayed row fell to 0.2629 (PerturbationAI;
their best 0.2872 stands); cqawesome 0.2618 (#2, mse raw 0.711 → 0.658 with pds flat); chill 0.2393 (#3, nmae 0.835
= board best with pds 0.824); a pack of eight at 0.231–0.235. Every top-20 team beats us on nmae; the pack also on
mse/fid/reach; we lead on pds/jac. Gap to the pack ≈ +0.011, all DE-side. Snapshot `data/leaderboard/2026-09-13T0544Z.json`.

**Slot 1: 387 sign-repair arm (`VCC_SIGN_REPAIR=cons` on 319p2), building since 06:17Z (~28 min + 12–17 min upload).**
Preregistration `experiments/387-sign-repair-arm/PREREGISTRATION.md`: the reading is the MEMBER SPLIT — the DE gain of
consensus signs (Δnmae + Δfid + Δreach, H1 says ≈ +0.0015) against the pds cost (H1 −0.0016, PLOLO up to −0.010);
the displayed row is expected to FALL and is not corrected. Why it teaches: today's lab (docs/lab-2026-09-13/) builds a
DE channel decoupled from the bulk (unwelded per-cell libraries, 389) and needs the A/B/C value of the cheapest sign
estimator (388 prices estimators on H1). Submit via `scripts/arm_submit.sh 202609130700 submissions/387-cap15-partial50-
signcons.vcc <sha> "the lint roller has opinions" experiments/387-sign-repair-arm/work/submit387.log` once the build's
one-key config diff passes; capture the prior entry first (`scripts/capture_board_entry.py lqWNXX5YaYArdglu6Stn`).
**Slot 2: expires** unless 389 delivers a KA-passing K=3 block AND a six-member H1 gain by ~20:00Z (not expected today).

**Tomorrow (09-14):** if 388/389 clear, the first unwelded arm = 319p2 + K=3 libraries + LFC channel at the 388-chosen
estimator/amplitude, one variable, preregistered vs 319p2 with the member split as the primary reading.

### 2026-09-13 slots, as executed (appended ~07:25Z)
**Slot 1: 387** (`VCC_SIGN_REPAIR=cons` on 319p2) = **0.2186049, Δ −0.00151 vs 319p2, lower bin** (pds −0.0038, DE +0.0014 inside
noise; `experiments/387-sign-repair-arm/BOARD-RESULT.md`). **Slot 2: 392** (`VCC_DE_SIGN_PRIOR=clean4` on 319p2; CPM-block signs from the 4-bank
cross-line prior through the unwelded DE channel; H1 end to end fid +0.0122 resolved, pds/mse exactly pinned; prereg
commit 175d8454) building since 07:14Z; `experiments/392-de-sign-prior-arm/arm392.sh` waits for the build gates and 387's
entry id, then `arm_submit_after.sh` submits once 387 is terminal.

## 2026-09-14 UTC — tree written before 392 scores (Δ vs 319p2 unless noted)
First check 392's falsifiers: score_pds and score_mse equal 319p2's (≤ 1e-9), jac within ±0.0005 raw. A failure voids
the reading; slot 1 then rebuilds 392 after fixing the bulk leak.
| 392 Δscore_fid | slot 1 | slot 2 |
|---|---|---|
| ≥ +0.020 (transfers) | the q80-protected variant as a one-variable step on 392 — PRICED end to end (394: H1-local +0.0071 [+0.0009,+0.0139] vs cap15, not separated from 392's construction); needs the `protect_q` option ported into `src/vcc2026/rc_emit_k.py` + knob value `clean4q80` after 392's build releases the builder, then a build | expire — the prior-agreement reach order (390 §8) and the mse shared-component scale (393) are both closed |
| −0.010 … +0.020 (no transfer) | read per-context members if exposed; otherwise try K562-essential alone (the best single donor on 4/7 held-out lines in 371) | expire |
| ≤ −0.010 (anti-informative) | close CLEAN4 on A/B/C; no preregistered arm is ready (393 closed: the shared component holds < 0.9 % of predicted energy) — expire unless a new candidate clears H1 | expire |

## 2026-09-14 / 09-15 UTC — all four slots UNSPENT (nobody at the machine after 392 published 09-13 08:10Z)

## 2026-09-16 UTC (reset 00:00Z) — written 05:45Z, before any result exists
**392 read (09-13 slot 2): 0.2226221 = +0.00251 vs 319p2, NEW BEST and incumbent.** Falsifiers pass (pds/mse exact, jac
−0.000006); Δscore_fid +0.0262 → top bin (the prior's signs transfer); nmae +0.0041 raw and reach −0.0039 raw where H1 said
+0.0024 / +0.0175. Record `experiments/392-de-sign-prior-arm/BOARD-RESULT.md`. Board 03:14Z: 13th of 974; Illumina 0.2979,
cqawesome 0.2741, CellWM 0.2473, nailong king 0.2459, chill 0.2426.

**Slot 1: 395 = `VCC_DE_SIGN_PRIOR=clean4q80` on 392** (394's construction: strongest 20 % of block calls keep the delivered
sign), preregistered `experiments/395-protected-signs-arm/PREREGISTRATION.md` (commit d86fd5b3) before the build; build
05:18Z (`run395.sh`), `arm395.sh` submits the moment the build passes its gates (no entry in flight). Predicted +0.0011;
primary reading = whether nmae and reach recover (the strong-gene flips) — bins in the registration.

**Slot 2 (rewritten 05:45Z, before 395 reads): 399 = the per-target neighbourhood tilt on 395** (`VCC_DE_SIGN_PRIOR=knn25o400q80`):
for each panel target its 25 nearest K562-GWPS perturbations vote on every gene through the six-bank family (CLEAN4 + raw Orion
HCT116/HEK293T), shrunk toward the family marginal (weight 400). Evidence, all pre-entry: 397 H1 analytic (called-block precision
0.534 vs CLEAN4 0.522, fid +0.0207 vs E0), 398/398c/398d leave-one-bank-out (over the CLEAN4 marginal: Jurkat +0.031, HepG2 +0.023,
RPE1 +0.032; random neighbourhoods lose everywhere). Registered vs 395's entry whatever its bin (`experiments/399-knn-tilt-arm/
PREREGISTRATION.md`); predicted +0.007 overall; the reading is Δfid. Build starts the moment run395.sh restores the builder
(`patch_builder_prior2.py --write`, unit tests, commit, `run399.sh`); `arm399.sh` submits after 395 is terminal via
arm_submit_after.sh. The Jurkat-only tilt for A (396, artifact + knob ready) moves to tomorrow, as does the 397c variant
(neighbourhood signs on uncalled genes: resolved nmae −0.0015 on H1, needs an emitter extension).

**October note (D/E/F):** the 399 artifact is panel-specific (rows are indexed by target); rebuild it on Oct 22 with
`uv run python experiments/399-knn-tilt-arm/build_prior399.py 25 400 orion` after the new `pert_counts.csv` is in place, re-pin
its sha in the builder (`DE_SIGN_PRIOR_SHA256` for the knn family), and re-run the builder unit tests. The bank hashes it records
must match the JSON sidecar of the A/B/C artifact.

### 2026-09-16 as executed (appended 06:25Z)
**Slot 1: 395 = 0.2239363 (`RzcGA26xZEigSuXgUp5t`), +0.00131 vs 392 — NEW BEST, incumbent.** pds/mse exact; reach +0.0100 raw
carried it, nmae −0.0017, fid −0.0019 (`experiments/395-protected-signs-arm/BOARD-RESULT.md`). The q50 dose is authorised by the
registration (reach carried the gain). **Slot 2: 399** (per-target neighbourhood tilt, six-bank family, shrink 400, on 395)
rebuilding after the emitter guard (06:01Z); `arm399.sh` submits when the build passes its gates (395 is terminal). Read its
`prior_lambda` distribution first (dose), then Δfid vs 395.

## 2026-09-17 UTC (reset 00:00Z) — tree written before 399 reads (Δ vs the incumbent = the better of 395 / 399)
| 399 result | slot 1 | slot 2 |
|---|---|---|
| fid ≥ +0.010 scaled and overall ≥ +0.0005 (transfers; 399 incumbent) | **402**: context A from Jurkat-only neighbourhood votes on 399 (`knn25o400jaq80`; registered `experiments/402-ctxA-jurkat-knn/PREREGISTRATION.md`; builder patch `patch_builder_402.py`, scripts ready) | q50 protection on the better of 399/402 (dose step authorised by 395), or 404's better neighbourhood if it clears its gate |
| flat / partial (399 not incumbent) | read `prior_lambda`: if under-dosed (many targets λ ≤ 0.5) → 399b with a less down-heavy shrinkage target (CLEAN4-family marginal as the default and shrink target, six-bank votes) on 395; else 402 on 395's construction is NOT available (per-target rows need the knn knob) → q50 on 395 | 402 only if slot 1 was a 399 variant that read positive |
| harm (≤ −0.0005) | q50 on 395 (protection dose) | expire unless 404 clears |
Standing: never on rank; capture every record; the board shows the LATEST entry.

## 2026-09-17 … 09-21 UTC — all ten slots UNSPENT (no session ran; 399 = 0.2223667 was read on 09-18, see memory)
399 read (09-16 slot 2): −0.00157 vs 395 (harm; reach gave back q80's gain). 402 was built on 09-16 but never submitted and its
package was swept on 09-19. 395 = 0.2239363 stays the incumbent.

## 2026-09-22 UTC (reset 00:00Z) — written 18:25Z (new session; Ahmet: "take over, go wild")
**The new line: per-target BULK magnitude + BALLAST unwelding (experiments/410-pertscale-revisit).**
- Key: key_p = mean_{i<j} sym<W d_i, d_j> / |W mean d|^2 over the three production donors (raw rows, W = union_pruned) —
  the cross-fitted plug-in of the mse-optimal per-target scale. PLOLO deep-truth folds: rho(key, s*) +0.58/+0.60/+0.63,
  aggregate cosine +27–30 % at clip(0.3, 3), +34–39 % uncapped, CIs above zero on every fold. pds invariant.
- Welded (H1): mse −0.049 raw but nmae +0.024 (the median-target shrink; a shuffled control pays the same) → net ≈ 0.
- Ballast (H1): carriers = 395's DE design on 400−K cells + K ≈995k-count cells putting per-gene totals on the scaled
  profile (ratio mode) + pair-wise anti-correlated ballast noise that restores the jackknife (credit) → pds +0.0056,
  mse −0.0496, nmae/fid/jac pinned, board-slope overall **+0.0097**. A/B/C jackknife margin C_hat/C_ctrl400 1.58–1.86.
**Slot 1: 410 = VCC_PERT_SCALE=wshare_lo0.3_hi3 on 395 (WELDED)**, submitted 18:10Z as "the houseplant filed its taxes".
Primary reading = raw mse vs 0.809898 (≤ −0.030 transfers; > −0.010 no transfer); DE members expected to pay.
**Slot 2: 410b = 410 + VCC_PERT_SCALE_BALLAST=1 (UNWELDED)**, building since 18:16Z (builder 40772dd1); `arm410b.sh`
submits after 410 is captured AND work/GO-410b exists (written 18:20Z). Readings: DE members back to 395's (±noise);
pds/mse at 410's.

## 2026-09-23 UTC — tree written before 410/410b read
| 410 Δ raw mse | 410b DE members vs 395 | slot 1 | slot 2 |
|---|---|---|---|
| ≤ −0.030 (key transfers) | pinned (unwelding works) | 410c = ballast + UNCAPPED key (H1 welded mse −0.089) | 410d = 410b/410c + 228's per-gene gain in the BULK (228: pds +0.0018, mse +0.0023 overall; its nmae −0.0104 was a weld cost) |
| ≤ −0.030 | leaks | fix the leak (K, ballast composition, compensator cells) on H1 first | expire or α step |
| −0.030 … −0.010 | pinned | production-LODO key (P^{−j} projected on d_j, includes CD4/H1 components) | 228 gain in bulk |
| > −0.010 | any | per-target key closed on A/B/C; ballast still usable for 228's per-gene gain | expire |
Standing: never on rank; the board shows the latest entry; capture every record.
**Added 18:40Z (H1, before 410 reads):** the UNCAPPED key through ballast (`bal_lo03_inf_ratio_jk`) reads mse −0.0913 raw,
pds +0.0007, nmae/fid/jac pinned, C_hat/C min 1.64 → board-slope +0.0159 (vs +0.0097 for clip 0.3/3). 410c's build
script is ready (`run410c_build.sh`, pert_scale hi3 → inf vs 410b). The key does NOT predict per-target sign accuracy on
the reference DE genes (H1 Spearman −0.03) and scaling the CARRIERS by it costs nmae +0.013: it is a bulk key only.
**DE-channel queue (ballast makes each context's carriers designable without touching pds/mse):**
(1) context-specific carrier designs — A = Jurkat: CD4-heavier signs; validate off-panel first (predict Nadig Jurkat's
essential-gene responses from the 3 donors ± CD4 at several weights, score |t|-weighted sign accuracy on Jurkat's own DE
genes); B/C epithelial: HCT116/HEK293T-heavier. (2) a DE-side per-target reliability signal (donor sign agreement on each
target's head genes) for nmae magnitudes. (3) 228's per-gene gain goes in the BULK only (pds +0.0018, mse +0.0023).

### 2026-09-22 as executed (appended 19:45Z)
**Slot 1: 410 (welded) = 0.2195152** (EETok0iqIbzQFg4CTYMQ): raw mse −0.0797 (the key transfers), weld cost nmae +0.054.
**Slot 2: 410b (ballast) = 0.2365722, NEW BEST, +0.0126 vs 395** (4WgVi6FIE4iVw9PmRUP3, rank 15): DE members within 0.0007
raw of 395, mse −0.0789, pds −0.0036. Predicted +0.0125.

## 2026-09-23 UTC — armed tonight
**Slot 1: 410c = 410b with the UNCAPPED key** (package built + prepped 19:42Z, `arm410c.sh` submits at 00:01Z; GO-410c
written under the preregistered rule). Bins vs 410b: ≥ +0.003 uncapped becomes the base; ±0.003 flat → keep capped;
≤ −0.003 → the tail over-weights.
**Slot 2: 410d = 410c + 228's per-gene gain on the bulk** (building; `arm410d.sh` waits for 410c's capture and GO-410d).
GO-410d rule (fixed now, before 410c reads): write GO-410d unless 410c shows a DE-member leak (|Δnmae| > 0.010 or
|Δreach| > 0.012 raw vs 410b) or a pds collapse (< −0.010 raw vs 410b). If 410c reads worse than 410b on overall by more
than 0.003, 410d still goes (its reading is the gene gain's increment over 410c, one variable) — the base choice for the
following days is made from both readings together.
Workers running (Opus 5.5 xhigh, Ahmet's pick): 413 source scan, 414 bulk re-pricing under ballast, 415 LODO key.

### 2026-09-23 plan REVISED 20:35Z (PI, before any 09-23 read)
**410d withdrawn** (NOGO-410d; 414 `URGENT-410d.md`): its effect-space energy match raises render-space energy +14..20 %,
past context A's cliff; predicted −0.003 to −0.021 overall. The render-gated fix is in the builder (`g228reg_t025rg`,
410d2) but is now fourth in the queue (+0.0025..+0.0046).
**Slot 1 → 417a = 410b + bulk-only sharpening q 1.25** (`experiments/417-bulk-sharpen/`; builder 400ef6e3, commit
907d6eb0; `arm417a.sh` submits at 00:01Z after its build gates, else execs the untouched `arm410c.sh`). 414 prices it
at +0.006..+0.016 overall vs 410c ≈ 410b, with 4/4 PLOLO folds positive and both controls negative. 414 prices 410c
itself at ≈ flat, so 410c no longer earns the first slot.
Bins and predictions: `experiments/417-bulk-sharpen/PREREGISTRATION.md`.
**Slot 2, rule fixed now:**
- 417a ≥ +0.006 vs 410b → **417b = q 1.5** (build queued after 417a).
- 0..+0.006 → **410c**.
- < 0 → **410c**, and the sharpening class closes.
- A DE member outside the noise band → no bulk arm until the ballast diagnostics are read. 410c goes.
**Slot-2 rule AMENDED 20:57Z (before any 417a read; 417 PREREGISTRATION addendum):**
- 417a < +0.006 → **417c = 417a + `VCC_CARRIER_GAMMA=0`** instead of 410c. 416's H1 gamma dose response on the carriers:
  γ0 nmae −0.0069 / fid +0.0025; γ0.9 nmae +0.0077.
- ≥ +0.006 → 417b.
- DE leak → 410c.
- 410c stays the fallback for any failed build gate.
- Waiter: `arm_slot2_0923b.sh`.
**Slot-2 rule AMENDED AGAIN 22:21Z (before any read; 417 PREREGISTRATION addendum 2, from 420's T-read table):** 417b
only if 417a − 410b ≥ **+0.009**, i.e. strong sharpening transfer. Otherwise 417c; DE leak → 410c. Waiter:
`arm_slot2_0923c.sh`.
Next-day ladder (420): stack (H1 ×2 + CD4 ×1.5 via LOAD) → + q 1.5 → + re-level, cumulative +0.015 / +0.018.

### 2026-09-23 as executed
- **Slot 1: 417a = 0.2435334, NEW BEST** (+0.0070 vs 410b; entry T88GsCjp53Bmk8xQcJMS). Bulk-only sharpening q 1.25:
  pds +0.0064, mse −0.0249 raw, DE pinned. Implied transfer k_sh 0.42 (pds) / 0.65 (mse).
- **Slot 2: 417c = 0.2418746** (−0.0017 vs 417a; entry 0xe3j7aqD2bypXKkCTBA). γ 0 carriers: fid +0.0012, reach +0.0024,
  nmae +0.0089 raw. H1 had predicted nmae −0.0068, so it is sign-inverted on nmae again.

## 2026-09-24 UTC plan (fixed 01:00Z 09-23, before any 09-24 read)
Built and gated, all bulk-only (the DE channel is pinned by ballast):
- 417e = 417a + stack LOAD (H1 ×2 + CD4 ×1.5), cliff 0.948/0.941/0.945;
- 417g = 417a + re-level (A 1.002 / B 1.041 / C 1.065), building;
- 417h = 417e + re-level (A 1.0088 / B 1.0500 / C 1.0736), queued after 417g.

**Slot 1 → 417e**, one variable vs 417a.
- 420 predicts +0.0037 (k = 1) / +0.0056 (k = 1.407); 418's history calibration gives +0.0044.
- Sharpening transferred at about half the PLOLO rate, so expect +0.002 … +0.004.

**Slot 2 rule:**
- 417e − 417a ≥ +0.0015 → **417h** (417e + re-level, one variable vs 417e).
- Otherwise → **417g** (417a + re-level, one variable vs 417a).
- A DE member of 417e outside the band (|Δnmae| > 0.010, |Δreach| > 0.012, |Δfid| > 0.003, |Δjac| > 0.001) → no automatic
  slot 2.
- Re-level predicted +0.0019 (417g) / +0.0024 (417h). The cliff model is validated by 417a's mse read (k_sh 0.65, consistent
  with pds).

Not sent: 417b, 417d (the q 1.5 dose is only worth it re-levelled at k_sh ≈ 0.5; later) and 417f.

### 2026-09-25 candidates (drafted 01:35Z 09-23; the rule is fixed only after 09-24 reads)
All bulk-only (DE pinned by ballast) and built overnight:
- **417k** = 417f (417b + stack LOAD, q 1.5) + its own re-level (A 1.0464 / B 1.0920 / C 1.1140). Against 417h this is
  the LEVELLED q step; 420 puts it at about +0.004 at k_sh 0.55.
- **417j** = 417h with the stack dumped under 420's LODO key at clip(0.5, 3) (one variable vs 417h: the dump). 420:
  about +0.003.
- **417i** = 417a + that key (one variable vs 417a).

Draft rule:
- If 417h is best after 09-24: slot 1 = 417k, slot 2 = 417j. Or 417j first, if 417e's reading makes the stack look
  shaky.
- If 417g is best: re-derive (417k/417j carry the stack), and 417i becomes the key test.

### 2026-09-24 slot 2 — REVISED 04:30Z 09-23 (before any read): 426a replaces the 417h/417g rule
**Why.** 425/426 found a DE-channel lever with a larger and newer reading than the re-level.
- nmae is scored only on the reference's significant genes, so its L1 optimum sits at the minimum detectable effect,
  ~1/sqrt(count).
- An MDE-structured add on the delivered arithmetic-mean LFC moves only the A-channel, so it is a pure-nmae arm.
- H1, 3 seeds: nmae −0.0096 ± 0.0002, fid/jac −0.0003/−0.0001, reach −0.003 ± 0.002, pds/mse 0: **+0.0019 overall**
  (the board has read pure magnitude at up to 2.3× H1).
- The re-level (417h/417g, predicted +0.0019/+0.0024) runs on a validated cliff model, so it loses nothing by waiting
  a day. The MDE is untested on the board and opens a dose family.

**Rule** (`experiments/426-mde-magnitude/arm0924_s2b.sh`; the old waiter pid 27781 was killed while still waiting for
417e's capture):
- 417e DE member outside the band vs 417a → no automatic slot 2.
- Otherwise → **426a** (417e + `VCC_DE_MDE=0.1,0.1,0.3`), one variable vs 417e whatever 417e's overall.
- Build gates: config diff exactly `de_mde` + emitter sha, CLIFF OK, novelty rc 0, package sha.
- Name "a western blot of my feelings", description empty.
- Bins in `experiments/426-mde-magnitude/PREREGISTRATION.md`: pds identical; nmae −0.010 … −0.022; overall +0.0018 …
  +0.0055.

The displayed row may fall on 09-24: 426a sits on 417e, whose stack is untested.

**Waiter fix (06:10Z 09-23).**
- `arm0924_s1.sh`'s clock loop gave up after 1,440 × 20 s = 8 h. It was armed at 00:57Z, so it would have exited
  around 09:00Z 09-23 without submitting, and slot 2 would then have waited in vain.
- Replaced by `experiments/417-bulk-sharpen/arm0924_s1b.sh`: same gates and log file, a 30 h clock cap, and a
  `NOT SUBMITTING` fallback that s2b recognises. The old pid 27780 was killed; the new pid is 83676.
- **Rule for every future waiter:** the clock-loop cap must cover the full time to the reset plus margin. Compute it
  from the target timestamp, never a fixed iteration count.

### 2026-09-25 (redrafted 06:55Z 09-23, before the 09-24 reads)
**Paired arms.** A DE-channel change (nmae/fid/reach/jac) and a bulk change (pds/mse) in the same arm give two
separable reads, because ballast pins each channel against the other. This is how the next steps get two reads per
slot. Built or building:

| package | = | vs its parent, DE part | vs its parent, bulk part |
|---|---|---|---|
| **428a** `c7bf5d47` | 426a + zero repair | nmae read of the repair | — |
| **428h** `c15d3680` | 428a + 417h's re-level (A 1.0088 / B 1.05 / C 1.0736) | vs 426a: the repair | vs 426a: the re-level (mse) |
| **428g** `f59ac498` | 417g + MDE + repair | vs 417g: MDE + repair | — |
| **428k** `338aa226` | 417k (stack, q 1.5, own re-level) + MDE + repair | same DE as 428h | vs 428h: q 1.25 → 1.5 (+ 417k's re-level) |
| 429a (building) | 428h + `VCC_DE_MAG=dn_k250_zs` | vs 428h: the zs pass (nmae only) | — |

**Draft rule, to be fixed after the 09-24 reads:**
- **417e ≥ 417a + 0.0015 (stack good) and 426a's nmae in the transfer or partial bin:**
  - slot 1 = **428h**: the repair (nmae) and the re-level (mse) on 426a;
  - slot 2 = **428k**: q 1.5 levelled vs 428h (pds/mse).
- **Stack bad:**
  - slot 1 = **428g**, on 417a's bulk: MDE + repair + re-level;
  - slot 2 = the best 417a-bulk step available.
- **426a's nmae in the no-transfer bin:** rebuild without the MDE (repair only) before choosing.

### 2026-09-25 (earlier draft, 05:30Z 09-23)
- Candidates:
  - **428a** = 426a + `VCC_DE_ZERO_REPAIR=-5` (building): a pure-nmae read of the zero-gene defect, which puts 24 gated
    genes at lfc ≈ −35;
  - 417h or 417g (re-level; the re-level can ride any DE arm, because it moves mse only);
  - 417k, 417j, 417i (bulk shape and key).
- Draft rule: whatever 426a reads, 428a goes (the repair cannot hurt). Its partner slot is the best bulk step, carried
  with the MDE (if 426a transfers) and the repair.
- The rule is fixed after the 09-24 reads.

### 2026-09-24 slot 2 — REVISED AGAIN 10:05Z 09-23 (before any read): 433a replaces 426a
**Why.** 431's per-block dump (429a's build) shows the plain MDE delivers **no dose in context A** (`mde_lambda` 0 in
92 % of A blocks; B 0.94, C 0.86). 426a's read would have been a B+C read of a delivery we would not ship.
- 426h's `coarse` delivery (production `5568ee84`) keeps H1's nmae (−0.00953 vs −0.0096) with 17× less class pressure.
- On the real data it restores the dose: 433a (= 417e + `VCC_DE_MDE=0.1,0.1,0.3,coarse`) delivers A 0.954 / B 0.975 /
  C 0.991 (means), ~1.65× 426a's MDE.

**Rule** (`experiments/433-coarse-mde-arm/arm0924_s2d.sh`, pid 65018; s2c pid 83836 was killed while still waiting for
417e's capture):
- 417e DE member outside the band vs 417a → no automatic slot 2.
- Otherwise → **433a**, one variable vs 417e: config diff exactly `de_mde` + emitter sha; CLIFF equal to 417e's.
- Package `d879cf9f…`. Name "the centrifuge wrote a memoir", description empty.
- Bins in `experiments/433-coarse-mde-arm/PREREGISTRATION.md`:
  - pds identical;
  - nmae −0.009 … −0.022 (central −0.012);
  - overall +0.0016 … +0.0055.
- 426a stays built and unsent. It is the plain-delivery twin: a later 426a-vs-433a read would price context A alone,
  and that is not worth a slot.

### 2026-09-25 (redrafted 10:05Z 09-23, before the 09-24 reads): the paired arms rebuilt with `coarse`
`chain433hk.sh` builds these after 433a (433h started 09:54Z):

| package | = | vs its parent, DE part | vs its parent, bulk part |
|---|---|---|---|
| **433h** | 428h with coarse: 417e stack, coarse MDE, zero repair −5, 417h's re-level | vs 433a: the repair (nmae) | vs 433a: the re-level (mse) |
| **433k** | 428k with coarse: 417k (q 1.5, its own re-level), coarse MDE, repair | same DE as 433h | vs 433h: q 1.25 → 1.5 + 417k's re-level |

**Draft rule, to be fixed after the 09-24 reads:**
- **Stack good (417e ≥ 417a + 0.0015) and 433a's nmae in the transfer or partial bin:** slot 1 = 433h, slot 2 = 433k.
- **Stack bad:** rebuild on 417a's bulk (417g's re-level + coarse MDE + repair) before choosing.
- **433a's nmae in the no-transfer bin:** the MDE family closes. The repair (from 428a) and the bulk steps go without it.
- zs (`VCC_DE_MAG`) rides only after 435 prices it on top of `coarse` (H1, running) and the MDE has transferred.
- **11:21Z 09-23, both 09-25 candidates built** (433k rebuilt after the 10:39Z reboot):
  - 433h `e7495644…`: config diff vs 428h exactly `de_mde` + emitter; 0 zero genes; CLIFF 0.9638 / 1.0337 / 1.0837;
  - 433k `f1dc369a…`: config diff vs 428k exactly `de_mde` + emitter; 0 zero genes; CLIFF 0.96 / 1.0278 / 1.0778.
  - Both are novelty rc 0.
- The fallback and the next rider are queued behind worker 432 (`chain433gz.sh`):
  - 433g: 428g with coarse, for the "stack bad" branch;
  - 433z: 433h + zs, for 09-26.

### 2026-09-24 00:49Z — 417e read: 0.2485853, NEW BEST (+0.0050519 vs 417a); the stack is good
- mse −0.0262 raw, pds +0.0014 raw; DE members pinned (|Δ| ≤ 0.0005). Record: `417-bulk-sharpen/BOARD-RESULT-417e.md`.
- 09-25 branch: **"stack good" (≥ +0.0015)**. The 09-25 pair is decided by 433a's nmae bin:
  - transfer or partial → slot 1 = 433h, slot 2 = 433k;
  - no transfer → rebuild without the MDE before choosing.
- 433a (coarse MDE on 417e) was submitted at 00:49:17Z by `arm0924_s2d.sh`.

### 2026-09-24 01:30Z — 433a read: 0.2536523, NEW BEST (+0.0050669 vs 417e). The MDE transfers (nmae −0.0168, 1.76× H1)
- pds and mse bit-identical to 417e; reach +0.0024. Record: `433-coarse-mde-arm/BOARD-RESULT-433a.md`.
- **09-25 rule FIXED ("stack good" + "MDE transfers"):**
  - slot 1 = **433h** (433a + zero repair + 417h's re-level) at 00:01Z. Waiter `433-coarse-mde-arm/arm0925_s1.sh`,
    pid 927; name "the flow cytometer took up jazz".
  - slot 2 = **433k** (433h with q 1.5 + 417k's re-level) after 433h's capture, unless 433h reads nmae > +0.001 or
    |Δpds| > 0.003 vs 433a. Waiter `arm0925_s2.sh`, pid 928; name "the incubator adopted a goldfish".
  - Bins: `PREREGISTRATION-0925.md`.
- **09-26 candidates** (bins in `PREREGISTRATION-0926.md`): 433z (+zs, pure nmae) and 433j (+LODO key). If 433k wins,
  both need rebuilding on 433k's bulk; decide after the 09-25 reads.
- **03:01Z 09-24: the 09-26 candidates are built for both 09-25 outcomes.** These are the first arms built through the
  new driver `5463fe28` and builder `63dfd7d0`. All are 0 zero genes and novelty rc 0.

  | if the 09-25 winner is | zs rider (pure nmae) | LODO key rider (bulk) |
  |---|---|---|
  | 433h | 433z `b4598c84`: vs 433h exactly `de_mag` | 433j `e1a34076`: vs 433h exactly `bulk_load` |
  | 433k | 433kz `e0c34441`: vs 433k exactly `de_mag` | 433kj `d34e68ff`: vs 433k exactly `bulk_load`; CLIFF 0.962 / 1.032 / 1.080 |

  Bins: `PREREGISTRATION-0926.md`. Read the 433k rows against 433k instead of 433h.

### 2026-09-24 20:20Z — new coordinator session (takeover after a 16 h gap); state and queue before the 09-25 reads
- Waiters alive: `arm0925_s1.sh` pid 927 (433h at 00:01Z), `arm0925_s2.sh` pid 928 (433k after 433h's capture);
  caffeinate pid 1032 until ~04:00Z; board_snapshot / pisces_watch alive. No reboot since 09-23 10:39Z.
- **Closed tonight (analytic / PLOLO, no slot):**
  - 440: the empty-union Jaccard lever (J(∅,∅) = 1). ~50 board control-copy entries read jac = 0.00000, so A/B/C has no
    target with an empty reference set;
  - 441: the DE sign gap is the context's own mean response on UNCALLED genes (oracle 0.62 vs ours 0.54); no donor,
    control-cell or panel-restricted route reaches it. Prior-confidence block selection trades jac for fid, net negative;
  - 443: 228's per-gene gain on today's sharpened stack (420 harness, KA |d| = 0): ~+1 % cos, kolf2 pds −0.005, closed;
  - 444: soft/hard thresholds and hybrids lose to the power law; the next shape step is q 1.75;
  - 445: DLD-1 merged into the H1 union dilutes it. As its own residual (builder `VCC_UNION2_*`, 248) it gives +0.8 % cos
    on the cancer folds, ~+0.001 overall. Low-priority rider.
- **Running:** 446u, a MEASUREMENT build of 433k with q 1.75 and no ctx scale. Its cliffcheck gives the q 1.75 re-level
  factors (the 417b → 417k route). The package is deleted on completion.
- **Plan for the reads:** read 433h / 433k against `PREREGISTRATION-0925.md`. The 09-26 slots follow the prebuilt table
  above: zs rider (pure nmae) and LODO-key rider (bulk), each one variable on the 09-25 winner's bulk. If 433k's q step
  transfers (433k ≥ 433h + 0.003, or pds Δ ≥ +0.004 raw), build q 1.75 levelled on the winner for 09-27.

### 2026-09-25 01:05Z — 09-25 read: 433h 0.2631120, then 433k 0.2677594, two NEW BESTS; 09-26 rule FIXED (branch A)
- 433h vs 433a: +0.0094597.
  - The zero repair gave nmae −0.0247 raw, 3× its top bin.
  - The re-level gave mse −0.0134 raw (transfers).
- 433k vs 433h: +0.0046475.
  - pds +0.0046 and mse −0.0179 raw; the DE members are pinned.
  - That is the calibrated price of 438's levelled q step.
- **09-26, armed (branch A):**
  - **slot 1 = 446b** (433k with q 1.75 at 433k's matched ρ, sha `e171cb39…`, prereg
    `446-q175-level/PREREGISTRATION.md`). `arm0926_s1.sh`, pid 66841, "the autoclave learned to knit".
  - **slot 2 = the LODO-key rider** after 446b's capture: 446bj2 (446b + LODO key re-levelled to 446b's ρ, `c7206c46…`)
    if 446b − 433k ≥ +0.001 with pds ≥ −0.0005 and DE members pinned, else 433kj (`d34e68ff…`). `arm0926_s2.sh`, pid
    66853, "the petri dish hired an agent". Key-rider bins: `433-coarse-mde-arm/PREREGISTRATION-0926.md`.
  - caffeinate pid 62208 until 04:09Z 09-26.
- **Built and waiting for 09-27+** (all on 446b's bulk, each at matched ρ where the bulk moves):
  - 446bz (zs, `8020a8d5…`, bulk identical);
  - 446bd2 (DLD-1 union 2, building 00:21Z; prereg `445-dld1-bulk-union/PREREGISTRATION.md`).
- 446bj and 446bd (unmatched ρ) were deleted unsubmitted.

### 2026-09-25 03:45Z — 09-26 slot 2 RE-ARMED to 446bja (before any 09-26 entry); 09-27 draft
- **Why.** 450 (H1) found that a same-line bank from another lab and chemistry carries the recipient's programme:
  gate-block sign 0.625 vs CLEAN4 0.581. 451 is building an October same-line mode on that. 446bja (446b + `clean4q80ja`:
  A's block prior from Nadig Jurkat) is the only board test of its programme part.
  - The H1 emitter price (Pan programme vs CLEAN4, 3 seeds) is fid +0.0096 and nmae −0.0086, both significant, and
    reach −0.041. The reach loss is four small-gate coin flips. Net +0.0001; +0.0048 without the two single-gene targets.
  - Prediction for A only: overall +0.0005 … +0.0030 (`450-same-line-programme/PREREGISTRATION-446bja.md`).
- **09-26 slot 1:** 446b, unchanged (`arm0926_s1.sh`, pid 66841).
- **09-26 slot 2:** 446bja vs 446b, whatever 446b scores (`arm0926_s2b.sh`, pid 86417, "the vortex mixer wrote a
  sonnet"). If 446bja's build is missing or fails a gate, the waiter falls back to the old rule (446bj2 / 433kj). The
  old s2 waiter (pid 66853) was stopped.
- **Built and waiting for 09-27** (one variable each vs 446b):
  - 446bj2 (LODO key, prereg `433/PREREGISTRATION-0926.md`);
  - 446br2 (repair floor −2, `9bbc8d2a…`, prereg `449/PREREGISTRATION.md`);
  - 446bz (zs);
  - 446bd2 (DLD-1 union 2).
- **09-27 draft rule:**
  - slot 1 = 446bj2 (the LODO rider moved from 09-26);
  - slot 2 = 446br2.
  - If 446bja reads overall ≥ +0.0010, it becomes the base and the 09-27 riders are rebuilt on it on 09-26 (≈ 40 min
    each).
  - If 446b reads flat or negative (q 1.75 closes), the riders move to 433k's bulk (433kj exists; 433kr2 would be
    built).

### 2026-09-25 18:40Z — the 09-27 rule, FIXED before any 09-26 entry exists
Inputs: 446b's entry (09-26 s1) vs 433k (T76JWD3tOajk5IYvzwUW), and 446bja's entry (09-26 s2) vs 446b.
- **Branch 1 — q 1.75 transfers** (446b − 433k ≥ +0.001 with pds ≥ −0.0005; 446's bins):
  - **1a — 446bja − 446b ≥ +0.0010 overall.** The Jurkat programme becomes part of the base.
    - s1 = 446bja + LODO key ("446bjaj2", built on 09-26 after the read);
    - s2 = 446bja + repair floor −2 ("446bjar2", built on 09-26).
  - **1b — otherwise.** The base is 446b.
    - s1 = 446bj2 (built, `c7206c46…`, 433's key-rider bins);
    - s2 = 446br2 (built, `9bbc8d2a…`, 449's bins).
- **Branch 2 — q 1.75 flat or negative.** The base is 433k (q 1.5).
  - s1 = 433kj (built, `d34e68ff…`);
  - s2 = 433kr2 (433k + floor −2, built on 09-26).
  - 446bja's programme reading still stands on its own (vs 446b). Carry it into the 433k base on 09-28 if 1a's
    condition held.
- In every branch the SHIPPED q for October follows 446b's read: 1.75 in Branch 1, 1.5 in Branch 2. It goes into
  DEF-RUNBOOK "THE DAY IN ONE PAGE" and `def_day.sh`'s argument, and `def_pipeline.py`'s SHIPPED_ENV follows with 437's
  row-1 known answer.

### 2026-09-26 01:10Z — 09-26 as executed; 09-27 armed (Branch 1b); 09-28 draft
- **s1 446b = 0.2700428, NEW BEST** (+0.0023 vs 433k). pds +0.0032 and mse −0.0076 raw are in their transfer bins.
  nmae +0.0013 raw is the larger ballast (K 4 → 6), −0.0004 overall. **q 1.75 is decided for October**
  (`446-q175-level/BOARD-RESULT-446b.md`).
- **s2 446bja = 0.2689767** (−0.0011 vs 446b). The same-line programme does not transfer on A: fid +0.0002 against
  +0.0015 … +0.0040 predicted; reach −0.0079. **October keeps `noprog`** (`450-same-line-programme/BOARD-RESULT-446bja.md`).
- **09-27 (the 18:40Z rule, Branch 1b), armed by `chain0926.sh` at 01:04Z:**
  - s1 = **446bj2** (446b + LODO key at matched ρ; `c7206c46…`; bins `433-coarse-mde-arm/PREREGISTRATION-0926.md`;
    "the petri dish hired an agent"; pid 50989);
  - s2 = **446br2** (446b + repair floor −2; `9bbc8d2a…`; bins `449-repair-floor/PREREGISTRATION.md`; "the pipette rack
    learned origami"; pid 50990).
  - Reader: `446-q175-level/read0927.py key <s1> UDgjOKL9gI1udAexq6ue floor <s2> UDgjOKL9gI1udAexq6ue`.
  - caffeinate pid 14235 runs to about 04:50Z 09-27. Re-arm before the 09-28 slots.
- **09-28 draft (the rule is fixed after the 09-27 reads):**
  - s1 = **453**, the October path run blind on A/B/C (`def_day.sh` on 180's relabelled bundle, relabelled back; q 1.75,
    `pt:noprog`). Building since 01:05Z (`453-abc-dress-rehearsal/chain453.sh`). Prereg `453/PREREGISTRATION.md`,
    counterpart 446b: predicted −0.0024 overall, almost all bulk (CD4 ×1.5 → ×1), DE pinned. It is the only board
    known-answer for the final entry's code path. It goes whatever the riders read, because its information does not
    depend on them.
  - s2 = **446bz** (446b + zs pass, pure nmae; `8020a8d5…`; bins in `433/PREREGISTRATION-0926.md`) if the base stays
    446b. If a 09-27 rider becomes the base, zs is rebuilt on it (~40 min).

### 2026-09-26 01:40Z — the 09-28 rule, FIXED before any 09-27 entry exists
Inputs are the 09-27 entries, both against 446b (`UDgjOKL9gI1udAexq6ue`):
- **key_ok:** 446bj2 − 446b ≥ +0.001 overall with pds Δ ≥ 0 (433's PREREGISTRATION-0926 bins).
- **floor_ok:** 446br2 nmae raw Δ ≤ −0.002 (449's bins).
- A missing or unscored 09-27 entry counts as not ok.

**Slot 1 = 453** (the October path run blind on A/B/C; `453/PREREGISTRATION.md`; "the lab notebook ran a marathon"),
whatever the riders read.
- Condition: its gate log `453/work/gates453.log` carries PASS, CLIFF OK (every final ρ ≤ 1), novelty rc 0 and the
  package sha. The coordinator writes that log only after checking the chain's audits.
- If 453 is not ready, slot 1 takes the zs arm below. Slot 2 is then 446bd2 (DLD-1 union 2; "the microscope joined a
  book club") when the base is 446b, else empty.

**Slot 2 = the zs pass (`VCC_DE_MAG=dn_k250_zs`, pure nmae) on the best 09-27 base**, one variable vs a scored entry:
- neither ok → **446bz** vs 446b ("the pipette tip went on strike");
- key only → **446bj2z** vs 446bj2 ("the fume hood took up yoga");
- floor only → **446br2z** vs 446br2 ("the gel box wrote a haiku");
- both → the one on the larger overall gain.
- If the chosen arm failed its build (456), use 446bz.

zs bins (433's PREREGISTRATION-0926): Δnmae ≤ −0.001 means zs transfers and rides the October stack; −0.001 < Δ ≤ 0 means
below readability (keep it only if ≤ 0); a Δ above 0 closes it.

Automation: `456-zs-riders/chain0928.sh` reads the 09-27 entries, applies this rule, and arms `456-zs-riders/arm0928.sh`
s1/s2. The 09-28 slots fall after this session's window (it ends about 19:25Z 09-27), so the chain and a caffeinate must
be live before then.
- Armed 01:24Z 09-26:
  - `456-zs-riders/chain0928.sh` (pid 53630) applies this rule after the 09-27 second slot. Its rule code was dry-run
    on real and synthetic entries; every branch resolves as written.
  - `456/chain456.sh` (pid 53074) builds 446bj2z and 446br2z after 455, which runs after 453.
  - caffeinate pid 53641 runs to 07:24Z 09-28.
  - bash here is 3.2, which has no associative arrays: `chain0928.sh` uses case lookups.
- 02:35Z: **453 is verified and ready for 09-28 slot 1.**
  - Package `a76283bd…`. `gates453.log` is written: every gate passes, and the config diff vs 446b is exactly the
    preregistered set.
  - It converged in 2 levelled builds at ρ 0.9793 / 0.9794 / 0.9795.
  - The day-of refit (operator, key, CD4 banks, census) reproduced the pinned A/B/C artifacts byte for byte.
  - 455 (the driver's SHIPPED q → 1.75, row-1 known answer = 446b) started 02:31Z. Then 456 (446bj2z / 446br2z), then
    452c.

### 2026-09-26 03:15Z — 09-29 DRAFT (fixed after the 09-27 reads, before any 09-28 entry)
- **s1 = 457**: the CD4 ×1.5 twin of 453 under October rules, one variable vs 453 (the dump's CD4 dose). Prereg
  `457-cd4-half-def/PREREGISTRATION.md`: +0.0017 … +0.0022 predicted, bins ±0.0010. It calibrates 438's instrument on the
  axis October's stack decision rests on. It builds after 452c (about 06:15–08:00Z today).
- **s2 = 446bd2**: DLD-1 union 2, one variable vs 446b (a scored entry), whatever the base; bulk riders are near
  additive. Prereg `445-dld1-bulk-union/PREREGISTRATION.md`. DLD-1's coverage grows on October's panel (67 → 73), which
  is why it is worth a slot.
- If the LODO key transfers on 09-27, its October day-of fit is written on 09-27. 438's code already re-fits 415's
  LODO-V1 key on any panel; the known answer must reproduce the A/B/C artifact `8c82bef3`.

### 2026-09-26 04:30Z — the 09-29 rule, FIXED before any 09-27 entry exists
- **s1 = 457**: the CD4 ×1.5 twin of 453, vs 453, bins in `457-cd4-half-def/PREREGISTRATION.md`; "the thermocycler
  joined a band". It goes only if `457/work/gates457.log` passes. The coordinator writes that log after checking the
  build with `gate457.py`.
- **s2 = 446bd2**: DLD-1 union 2, vs 446b, bins in `445-dld1-bulk-union/PREREGISTRATION.md`; "the microscope joined a
  book club". It goes only if it was built and never submitted: 09-28's rule can use it as a fallback.
- **If 457 is not ready,** s1 = 446bd2 and s2 is empty.
- Neither arm depends on a 09-27 or 09-28 read, so `457/chain0929.sh` arms `457/arm0929.sh` as soon as the 09-28
  waiters exit, whatever they read.
- caffeinate pid 73126 runs to 2026-09-30T05:54:18Z, covering the 09-29 slots.

### 2026-09-26 05:10Z — the 09-30 rule, FIXED before any 09-27 entry exists
- **s1 = 458:** Arc's own H1 programme as every context's block prior, vs 446b. Bins in
  `458-arc-lab-programme/PREREGISTRATION.md`: the lab-signature hypothesis, and an October lever for all contexts if
  fid ≥ +0.0015. "the autoclave read the manual".
- **s2 = 459:** the H1 ×3 twin of 453 under October rules, vs 453. Bins in `459-h1-axis-def/PREREGISTRATION.md`: 438's H1
  axis, predicted +0.0004. "the spectrometer took a nap". It goes only if `459/work/gates459.log` passes; the coordinator
  writes that log after `gate459.py`.
- **If 458 is not ready,** s1 = 459 and s2 is empty.
- Independent of every 09-27…09-29 read. `459/chain0930.sh` arms `459/arm0930.sh` once the 09-29 waiters exit.
- The caffeinate (pid 73126) runs to 05:54Z 09-30.
- Build queue from 05:00Z 09-26: 452c → 457 → 458 → 459, strictly serial. Each chain waits for the previous chain's own
  END and process exit.

### HANDOFF CHECKLIST (drafted 2026-09-26 05:30Z; refresh before the session ends ~19:25Z 09-27)
**Armed, no session needed:**
- 09-27: `446-q175-level/arm0927.sh` s1 446bj2 (pid 50989), s2 446br2 (pid 50990).
- 09-28: `456-zs-riders/chain0928.sh` (pid 53630) → `arm0928.sh` s1 453, s2 zs on the best base.
- 09-29: `457-cd4-half-def/chain0929.sh` (pid 73089) → `arm0929.sh` s1 457 (needs `gates457.log`), s2 446bd2.
- 09-30: `459-h1-axis-def/chain0930.sh` (pid 76522) → `arm0930.sh` s1 458, s2 459 (needs `gates459.log`).
- Builds, strictly serial: 452c → `457/chain457.sh` → `458/chain458.sh` → `459/chain459.sh`.
- caffeinate pid 73126 runs to 05:54Z 09-30. Watchers: board_snapshot (2910), pisces_watch (2911).

**Manual, by whoever is on:**
1. After chain457 ENDs: run `uv run python experiments/457-cd4-half-def/gate457.py`, and if every gate passes add
   `--write`. Without it, 09-29 s1 falls back to 446bd2.
2. After chain459 ENDs: the same with `experiments/459-h1-axis-def/gate459.py`.
3. After each day's entries: read them with the readers written before the entries existed:
   - `446-q175-level/read0927.py`;
   - `456-zs-riders/read0928.py` (dress 453 vs 446b; zs vs its base);
   - 457 vs 453, 446bd2 vs 446b, 458 vs 446b, 459 vs 453, each against its PREREGISTRATION bins.
   Then write BOARD-RESULT files, the BOARD-LEDGER top section and FINDINGS §0.
4. Follow-ups the reads trigger:
   - The LODO key transfers → write its October day-of fit. 438's code re-fits LODO-V1 on any panel; the known answer
     is the A/B/C artifact `8c82bef3`.
     Recipe (read 06:40Z 09-26):
     - `experiments/438-october-stack-repricing/full438.py` computes the full-panel key: `kv1 = L.lodo_key(FC.final("h1c1"), ...)`,
       then `L.weights(kv1, carried, 0.5, 3.0)`. On `abc` it asserts kv1 == 415pp's key to 1e-12.
     - Steps: (1) register the day's bundle as a `lib438.bind` spec (bundle, DEF rules, day operator, `<work>/cd4_compact`);
       (2) write kv1's (0.5, 3) weights in 410's artifact format (targets, w, key, carried, n_donors);
       (3) pin it under a new id (e.g. `lodo415_lo0.5_hi3_def`) with `def_pin_panel_artifacts.py`;
       (4) export that id as VCC_PERT_SCALE for the day's dump stage in `def_day.sh`.
     - Known answer before use: steps 1–2 on `data/challenge` reproduce `pert_scale_415pp_lo0.5_hi3.npz` (8c82bef3) array for array.
   - Floor −2 or zs transfers → move the SHIPPED_ENV default in `def_pipeline.py` and re-prove row 1 (a 455-style chain;
     row 1 = the adopted arm's own artifact sha).
   - 458 positive (fid ≥ +0.0015 vs 446b) → switch October to Arc's programme for every context. Path (~1 h plus a
     rehearsal):
     1. Add `--programme-source <bank>` to `scripts/def_sameline.py`, default off. For EVERY bundle context it writes
        `programme_<c>` / `pdown_<c>` from that bank (396's rule, the day's panel excluded) with `programme_ok` true.
        Unmatched contexts get 0 rows and w 0; matched contexts keep their same-line rows.
        Known answer: with the flag off, arrays equal 452c's `sameline.npz` (compare arrays; zip timestamps differ).
     2. Pass it through `def_day.sh` (a new env, e.g. `DEF_DAY_PROGRAMME_SOURCE=h1_2025`), and run the day with flags
        `pt` instead of `pt:noprog`.
     3. Rehearse on `ef-bundle-r` (a 452d); the programme changes nothing in the bulk, so the level is unchanged.
     4. `artifact458.py` is the A/B/C template: same rule, same bank.
   - **458b DROPPED (10:23Z 09-26, `458/PREREGISTRATION.md` addendum):** it flips only 8 % of 458's block signs, so
     it cannot read distinguishably from 458. The original note follows for the record.
   - 458 null or harm → a direct test of its stated risk, the ESC-specific content (458b):
     - Take the z-scored Arc-H1 P_down and regress out Pan-H1's (same line, Whitehead). That removes the component
       both H1 screens share (pluripotency exit) and keeps Arc's lab-specific part.
     - Build it with `artifact458.py`'s recipe on the residual's sign, and read it against 446b. Only worth a slot if
       458's fid read is not clearly harmful.
5. 10-01 onward: nothing is armed. Fix the rule before the previous day's entries exist.

### 2026-09-26 10:32Z — the 09-30 rule REVISED (before any 09-27 entry exists): 458d replaces 458 in slot 1
- **s1 = 458d** ("the pipette learned to meditate"). It is Arc's H1 programme from ALL public H1 cells: 300
  perturbations × ~1,071 cells against 38,176 controls. Vs 446b, with 458's bins unchanged.
  - Why: 462c measured that 458's harness-bank programme is about a third control-sample noise (disjoint-control
    split-half 0.39 vs shared 0.79). 458d's is 0.69.
  - A 458 null would be ambiguous and a 458 positive would still need the clean programme, so 458d is the single test.
  - Fallback: 458 if 458d is not built and gated. `458/PREREGISTRATION.md` addendum 10:30Z.
- **s2 = 459dr**, unchanged.
- Automation:
  - `459/chain0930b.sh` (pid 36293) replaces `chain0930.sh`, which was killed. Same logic, with 458d added first.
  - `458/chain458d.sh` (pid 36294) builds 458d via `run_arm_generic.sh` once 460's rehearsal ends and the builder is
    restored.
- October, if it carries: `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` with flags `pt` (DEF-RUNBOOK).
- Build queue from 10:05Z 09-26: 460 (LODO October-path rehearsal) → 458d → 461 (October path + Arc's programme, the
  twin of 458d; chain written after 460 releases `def_day.sh`).

### HANDOFF CHECKLIST — refreshed 2026-09-26 10:42Z (supersedes the 05:30Z list where they differ)
**Armed, no session needed:**
- 09-27: `446-q175-level/arm0927.sh` s1 446bj2 (pid 50989), s2 446br2 (pid 50990).
- 09-28: `456-zs-riders/chain0928.sh` (pid 53630) → `arm0928.sh` s1 453dr, s2 zs on the best base.
- 09-29: `457-cd4-half-def/chain0929.sh` (pid 73089) → s1 457dr (needs `gates457.log`, written), s2 446bd2.
- 09-30: **`459-h1-axis-def/chain0930b.sh` (pid 36293)** → s1 **458d** (else 458), s2 459dr (`gates459.log`, written).
  `chain0930.sh` was killed at 10:31Z.
- Builds, strictly serial:
  - `460/chain460.sh` (running): 460dr, the LODO October path.
  - `458/chain458d.sh` (pid 36294): 458d via run_arm_generic.
  - `461/chain461.sh` (pid 37105): 461dr, the Arc-programme October path. It needs the `def_day.sh` pass-through
    below.
- caffeinate pid 73126 runs to 05:54Z 09-30.

**Manual, by whoever is on:**
1. After chain460 ENDs:
   - `uv run python experiments/460-lodo-dayof/gate460.py`, then `--write` if it passes.
   - Add `def_day.sh`'s `DEF_DAY_PROGRAMME_SOURCE` pass-through (`${DEF_DAY_PROGRAMME_SOURCE:+--programme-source …}`
     on the `def_sameline.py` call) and commit it.
   - Do both before chain458d ENDs, or chain461 stops fail-closed.
2. After chain461 ENDs: `gate461.py`, then `--write`.
3. 09-27 reads (`446-q175-level/read0927.py key <s1> UDgjOKL9gI1udAexq6ue floor <s2> UDgjOKL9gI1udAexq6ue`):
   - key_ok → DEF-RUNBOOK: `DEF_DAY_LODO=1`; 460dr goes to a later slot as the October-path confirmation;
   - floor_ok → `bash experiments/463-floor-adopt/chain463.sh "<delta text>"` (SHIPPED_ENV −5 → −2, row 1 vs 446br2);
   - write the BOARD-RESULT files, the BOARD-LEDGER top section and FINDINGS §0.
4. 10-01 onward: nothing is armed. Fix the rule before the previous day's entries exist.
   - Candidates: 460dr if key_ok; 461dr if 458d reads positive; the final October-configuration rehearsal once
     09-27…09-30 are read.
   - Queued idea: swap `h1_2025` in the H1 union for `h1_2025_full` (25 panel rows vs 20, ~2.7× the cells). Needs a
     union-registry entry and a PLOLO/board price.
- 10:49Z addition: `458/chain458e.sh` (pid 38635) builds 458e after 461. 458e is the 0.5 CLEAN4 + Arc-full blend, the
  dose-curve partner of 458d. It is a 10-01 candidate only if 458d reads positive or partial (`458/PREREGISTRATION.md`
  addendum 10:49Z).
  - Build queue: 460 → 458d → 461 → 458e, about 17:00Z.
  - 10-01 draft: s1 461dr if 458d is positive (else 460dr if key_ok); s2 458e if 458d is positive or partial.
  - Fix the rule before 09-30's entries exist.
- 10:50Z: readers for the 458 family are ready: `458/read458.py prog <458d> <446b>`, `dose <458e> <446b> <458d>`,
  `path <461dr> <453dr> <458d> <446b>`. Its known answer: `prog` on 446bja vs 446b reproduces that read (fid +0.0002,
  null, bulk pinned).

### 2026-09-26 12:05Z — the 10-01 rule, FIXED before any 09-27 entry exists
Inputs:
- **key_ok**: 09-27 s1 (446bj2) − 446b ≥ +0.001 overall with pds ≥ 0.
- **P**: the 09-30 s1 entry (458d, or the fallback 458) vs 446b by raw fid:
  - positive ≥ +0.0015; it is **positive_reach** if also reach ≤ −0.005 raw and overall ≤ 0 (458's addendum 11:31Z);
  - partial +0.0005 … +0.0015;
  - otherwise null, harm or missing.

Candidates, in order; the first two that are ready (gate, package, novelty) take s1 and s2:

| P | candidates |
|---|---|
| positive | 461dr ("the bunsen burner took a gap year", vs 453dr) → 458e ("the lab coat filed a complaint", vs 446b) → 460dr if key_ok |
| positive_reach, partial | 458e → 460dr if key_ok |
| null, harm, missing | 460dr if key_ok ("the fridge magnet learned calculus", vs 453dr) |

With no eligible arm nothing is armed, and the next session decides.
- Automation: `experiments/464-slot-1001/chain1001.sh` (pid 53393) waits for the 09-30 slot 2, applies the rule and
  arms `arm1001.sh`. Its rule code was dry-run on real and synthetic entries; every branch resolves as written.
- caffeinate pid 53403 runs to 12:00Z 10-01.
- Readers: `458/read458.py` (path, dose), `460/PREREGISTRATION.md` bins for 460dr.
- **Manual before then:** `461/gate461.py --write` once 461 ENDs (today). Without `gates461.log`, 461dr is never eligible.
- 12:10Z build queue (all serial, automated):
  - 461dr (running; `gate461.py --write` after it, which 10-01 needs);
  - 458e (`458/chain458e.sh`, pid 38635);
  - 462l, the H1 six-member profile of Arc's programme (`462/drive462.sh`, pid 47555);
  - 465, a rehearsal with BOTH new switches on the efr panel (`465-both-switches-efr/chain465.sh`, pid 54337; not
    submittable).
- 13:22Z: **461dr PASS, `gates461.log` written** (the manual step for 10-01 is done). The config diff vs 453dr is
  exactly `sameline`, and its re-level retraced 453's exactly. 458e is building (started 13:21Z).
- 13:31Z: `466-lab-signature-cells/chain466.sh` (pid 66129) queued after 465. It measures where a programme's lab
  signature lives, with single cells: Xaira HCT116/HEK293T vs Weissman K562, non-essential vs essential targets,
  disjoint controls. Analysis only.
- 15:56Z: 465 PASS (both switches on efr; see DEF-RUNBOOK). `467-strength-gated-programme/chain467.sh` (pid 80009)
  queued after 466. It is the H1 price of Arc's programme applied only to design-strong targets: 462l's reach cost is
  four small-gate targets, and reach lives in weak targets. If it holds, it is a candidate construction for a 458d
  follow-up (and a builder knob: the prior per target by predicted strength).

### 2026-09-26 17:04Z — the 09-30 rule REVISED again (before any 09-27 entry exists): s2 = 458g, else 459dr
- **s1 = 458d**, unchanged ("the pipette learned to meditate"; fallback 458).
- **s2 = 458g** ("the cryostat started a podcast"; renamed 17:15Z from "the petri dish joined a book club") when s1 is 458d and 458g is built and gated. Otherwise s2 is
  459dr as before.
  - 458g = 458d + `VCC_DE_PRIOR_GATE=q33`. Arc's programme goes only to each context's design-strong two thirds, and
    the weakest third keeps CLEAN4.
  - Read against 458d, the same day. Bins in `468-prior-gate/PREREGISTRATION.md`: falsifiers pds/mse ≤ 1e-6 and jac
    ±0.0010. Primary Δoverall: ≥ +0.0015 transfers; +0.0005 … +0.0015 weak positive, keep; ±0.0005 null, ship
    ungated; ≤ −0.0005 close.
  - Why: 467 on H1 (3 seeds): gated33 +0.0110 vs ungated +0.0067 overall-equiv; paired +0.0043, seed sd 0.0002. The
    reach cost (−0.0388 → +0.0015) is exactly the four small-gate targets, all with strength 1. A sign prior's reach is
    not priceable on H1 (H1-BOARD-CALIBRATION), so the board is the test. The pair isolates the gate: same seeds,
    same design, same day.
  - 459dr (expected +0.0004, a calibration of 438's H1 axis) moves to a later day; its package and gate log stand.
- Automation: `459/chain0930c.sh` (pid 88484 after the 17:15Z rename; first armed as 86637) replaces `chain0930b.sh` (killed 17:02Z). Every branch was dry-run:
  - now: s1 458d, s2 459dr (458g not yet built);
  - with 458g ready: s1 458d, s2 458g;
  - with 458d missing: s1 458, s2 459dr.
- 458g is built by `468-prior-gate/chain468.sh` (pid 85577) after its two known answers:
  - KA-1: all knobs unset rebuilds 446b byte for byte;
  - KA-2: VCC_SAMELINE without the gate rebuilds 458d byte for byte.
- **10-01 follow-up (fix before 09-30's entries exist):** the 10-01 table does not know 458g yet.
  - If `DEF_DAY_PRIOR_GATE` is built and rehearsed by then, a gated October path (461g) joins the candidates after
    461dr.
  - Its eligibility would depend on 458g's read, which is only known after 09-30 s2. So the 10-01 rule can use it only
    through the chain's run-time read.

### 2026-09-26 17:16Z — the 10-01 rule REVISED (before any 09-27 entry exists): 458g's read and 469dr join it
Inputs as at 12:05Z (key_ok from 09-27 s1; P from 458d vs 446b), plus **G** = 458g − 458d from 09-30 s2, by 468's bins:
- **bug** if pds or mse moved more than 1e-6 raw, or |jac| > 0.0010;
- **pos** at ≥ +0.0005 overall; **null** within ±0.0005; **harm** at ≤ −0.0005;
- **none** without a 458g entry.

**carried** = 458g − 446b fid ≥ +0.0015 raw and overall ≥ 0 (468's table, row 3).

| P | G | candidates (first two ready) |
|---|---|---|
| positive / positive_reach / partial | pos | **469dr** ("the sonicator adopted a cat", vs 453dr and 461dr) → 460dr if key_ok |
| positive | not pos | 461dr → 458e → 460dr if key_ok (the 12:05Z rule) |
| positive_reach / partial | not pos | 458e → 460dr if key_ok (the 12:05Z rule) |
| null / none | pos and carried | 469dr → 460dr if key_ok |
| otherwise | | 460dr if key_ok |

- 469dr is the October path with the programme and the gate on the relabelled A/B/C bundle: `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=q33 def_day.sh`, flags `pt`.
  - `469-gated-october-path/chain469.sh` builds it after chain468.
  - `gate469.py`: its config must differ from 461dr's only in `de_prior_gate` (the same-line artifact identical
    apart from its path).
  - Without `gates469.log` it is never eligible.
- `def_day.sh` gains `DEF_DAY_PRIOR_GATE=q<NN>` (468). It is refused without `DEF_DAY_PROGRAMME_SOURCE`, and exported
  only after the stack dump. Unset, the script is unchanged: every added line is inside a guard on
  `DEF_DAY_PRIOR_GATE`, or on `VCC_DE_PRIOR_GATE`, which the script's VCC_* sweep removes. 469 is its end-to-end test.
- Automation: `464-slot-1001/chain1001b.sh` + `rule1001b.py` replace `chain1001.sh` (killed 17:16Z). A synthetic-entry
  dry run passed all nine branches (the eight above plus no entries at all), with `BOARD_DIR` pointing at synthetic
  entries in the scratchpad.

### HANDOFF CHECKLIST — refreshed 2026-09-26 17:17Z (supersedes the 10:42Z list where they differ)
**Armed, no session needed:**
- 09-27: `446-q175-level/arm0927.sh` s1 446bj2 (pid 50989), s2 446br2 (pid 50990).
- 09-28: `456-zs-riders/chain0928.sh` (pid 53630) → s1 453dr, s2 zs on the best base.
- 09-29: `457-cd4-half-def/chain0929.sh` (pid 73089) → s1 457dr, s2 446bd2.
- 09-30: **`459-h1-axis-def/chain0930c.sh` (pid 88484)** → s1 458d, s2 **458g** (else 459dr).
- 10-01: **`464-slot-1001/chain1001b.sh` (pid 88678)** with `rule1001b.py` (reads 458g; 469dr is a candidate).
- Builds and runs, strictly serial:
  1. `468-prior-gate/chain468.sh` (pid 85577): tests (passed 16:57Z), KA-1 (unset == 446b), KA-2 (VCC_SAMELINE, no gate
     == 458d), then the 458g build. Any failure restores `4feae27d` from `468/work/builder-4feae27d.py`, and 09-30 s2
     falls back to 459dr by itself.
  2. `469-gated-october-path/chain469.sh` (pid 87285): 469dr, the gated October path. **It waits until HEAD's builder
     is `693fde15`** (deadline 00:30Z 09-27).
  3. `470-weak-target-prior/chain470.sh` (pid 88244): H1, which prior the design-weak third wants (preregistered).
- caffeinate runs to 12:00Z 10-01 (pid 53403).

**Manual, by whoever is on:**
1. **The builder `scripts/make_submission_pm_cd4.py` is PATCHED AND UNCOMMITTED:** `693fde15` = `4feae27d` +
   `468/patch_builder468.py`.
   - After chain468 logs KA-1 PASS and KA-2 PASS and ends (builder restored, not pinned), commit it, staging the path
     explicitly, with both KA output shas in the message. That releases chain469.
   - If either KA failed, the chain already restored `4feae27d`: do not commit, and diagnose.
2. `scripts/def_day.sh` is also uncommitted (`5f9bfa26`, `DEF_DAY_PRIOR_GATE`). After chain469 ENDs:
   - run `uv run python experiments/469-gated-october-path/gate469.py`, then `--write` if it passes;
   - commit `def_day.sh` with 469's result.
3. After chain470 ENDs: write 470's RESULTS.md against its preregistration.
4. The 09-27 reads, as in the 10:42Z list, item 3 (`read0927.py key … floor …`, BOARD-RESULT files, BOARD-LEDGER,
   FINDINGS §0; key_ok → `DEF_DAY_LODO=1`; floor_ok → chain463).
5. Queued, not built: the h1_2025_full union swap. It needs `b_ctrl.npy` (the control pseudobulk; `denoise458d.py`'s
   loop), a registry entry and a dump. It adds 3 panel targets (HSBP1, SHPRH, ZFP62) and cleaner rows for 22.
- 18:36Z: **458g built and gated.**
  - Config diff vs 458d exactly `de_prior_gate`; package `eb502099…`.
  - Gate threshold 2 in A/B/C: 86 / 83 / 83 weak targets keep CLEAN4.
  - `chain0930c`'s readiness check passes for 458d and 458g, so 09-30 = s1 458d, s2 458g.
  - Builder `693fde15` committed (dc05673e) with KA-1 (446b) and KA-2 (458d) byte for byte. chain469 started 18:35Z.
- 18:41Z: readiness dry run of the armed chains, using their own `pkg`/`rl`/`ready` functions. Everything they
  may arm passes (package present, PASS, CLIFF OK, novelty rc 0):
  - 09-28: 453dr, 446bz, 446bj2z, 446br2z, 446bd2;
  - 09-29: 457dr, 446bd2;
  - 09-30: 458d, 458g.
  - 10-01: 469dr is still building. 461dr, 458e and 460dr are built.
- 458g's pair has no emission re-roll noise. `compare458g.py` found an exact split: every 458g block equals 458d's
  (strong) or 446b's (weak). So 458g − 446b is the programme on strong targets, and 458d − 458g on weak ones.
- 19:56Z: **469dr PASS, `gates469.log` written.**
  - Config diff vs 461dr is exactly `de_prior_gate`; re-level retraced 461dr; package `917c217b…`.
  - `def_day.sh` `5f9bfa26` committed (8b8ef484).
  - 469dr passes chain1001b's readiness check, so every 10-01 candidate is now built: 469dr, 461dr, 458e, 460dr.
  - chain470 started 19:53Z (H1 weak-target priors); chain471 follows it.
- 20:00Z: the queue is 470 (running; KA PASS; the first all-down seed ties CLEAN4) → 471 (the H1 bank swap in
  the stack union: b_ctrl, registry patch, KA-1, dump KA, arm) → 470b (H1: the strong side's programme from Arc's
  non-essential vs essential perturbations; `chain470b.sh` pid 7834). All preregistered. After 471 ENDs, commit its
  builder with its KA-1 and dump-KA results.
- 20:16Z: **peer review** (`docs/REVIEW-2026-09-26-second-opinion.md`, from a second session Ahmet started). Acted on:
  - 1. the October gate should be absolute: `le1`. 472 is queued after 470b; KA-le1 must equal 458g byte for byte.
  - 2. target-bootstrap intervals are now in 467, 468's addendum, FINDINGS and C-71. The gate's own H1 benefit is not
    established (P(≤ 0) ≈ 0.15), so the 09-30 pair is the test.
  - 3. same-line rows and ballast compliance: Ahmet ruled both allowed.
  - The 41 local commits were pushed.
  - Queue: 470 → 471 → 470b → 472.
- 20:41Z: **board intelligence** from the peer session (`docs/BOARD-INTEL-2026-09-26.md`: 176 snapshots,
  11,731 entries, 3,934 descriptions). Coordinator's judgment:
  - **Lead 1 (re-price 410d2), queued as 473 for 09-27 daytime.** The per-(target, gene) reliability multiplier (228),
    bounded to m ∈ [1/2, 2] and amplitude-matched, was priced +0.0025 … +0.0046 on 410c's bulk (414). It was shelved
    for queue position on 09-23, never closed, never built. Testh's board history shows a per-(target, gene)
    reliability weight took pds +0.015 (09-12), and a learned version with a per-target weight took mse −0.11
    (09-14). The bulk has since moved (stack, q 1.75, re-level), so it needs a re-price on 446b's bulk first, with
    414's PLOLO tooling (`emul410d2.py`).
  - **Lead 2 (a learned DE-membership model, chill/Testh), not queued.** It overlaps closed items:
    - 406 (order keys);
    - 434 (a learned L1 magnitude with panel donors: −0.0066 vs the MDE's −0.0111);
    - 430 (sign sources at the public-data limit);
    - 358 (label-free sign classification: dead).
    Cross-line DE membership replicates poorly: strong non-cis K562 knockdowns are DE in H1 12 % of the time. That
    bounds any membership model built from public cross-line data. Worth a note, not a slot.
  - **chill's nmae (0.835 from its first entry) is UNEXPLAINED** (corrected after the peer's caveat):
    - its bulk was weaker than ours then (pds 0.824, mse 0.845 on 09-11), so better data should have shown there too;
    - it describes itself as source-trained with no query outcomes;
    - 434f's accuracy map assumes iid signs with constant magnitudes. Magnitudes that follow confidence can reach 0.835
      at a lower average accuracy.
    434's conditional-median L1 model adapts magnitude to confidence, so it is the closest test we ran. It lost to the
    MDE with panel donors.
- 21:03Z: the peer session launched three analysis-only workers at Ahmet's request: 480 nmae instrument, 481
  October constants audit, 482 literature/data hunt (`cursor-agent`, isolated worktree). Numbers 480–482 are theirs;
  474–479 stay ours.
  - Memory is the shared resource. free+inactive was 23.2 GiB at 21:02Z (Chrome 11.7 GB, compressor 11.8 GiB), right at
    `def_pipeline`'s 23 GiB refusal.
  - 472's KA runner now waits for memory, like the dump and arm runners.
  - The peer was asked to keep the workers light or gated on vm_stat until ~23:45Z.
- 21:34Z: the peer session's cloud agents (Grok 4.7 xhigh on Cursor Cloud VMs, no Mac load, own `cloud/*`
  branches):
  - 484 builder bug hunt;
  - 485 October day-of path bug hunt;
  - 486 public-data LOCO test of a learned per-(target, gene) reliability gate (the companion to our 473).
  - 482c (literature/data) finished on `cloud/482c-lit-data`. All negative:
    - no new genome-scale cancer-line CRISPRi is public;
    - X-Atlas/Pisces is still empty (our pisces_watch covers it);
    - HeLa and CAL-33 have no public atlas;
    - Arc names no October lines.
    One new A549 Cas9-KO CROP-seq (GSE345058, 1,000 genes) is a weak donor. Pull the branch when convenient.
  - Numbers 480–486 are the peer's; 473–479 are ours.
- 21:59Z: **cloud 484 (builder bug hunt)**, brought into master from `cloud/484-builder-review`:
  - **CONFIRMED**: q80 protection does not protect down genes' magnitudes. `side` is computed after `tgt = E0`, so protected
    down genes take the mass-balance kf (median 0.60 on 432 b2), while protected up genes keep E0. It has been in the
    scored path since 394/395. Only nmae can move.
  - Priced as **475** on H1 (a patched emitter copy; production untouched) after 474. Its rule is preregistered.
  - Smaller items: a q-form threshold ≤ 0 is a silent no-op (472's le1 avoids it; a q guard can ride with the next
    builder patch); the zero-repair residual for highly expressed knockouts (−7.4, not −35); two CPM>5 gates; a ballast
    remainder crash on adversarial input; the target taking a head slot.
  - 482c's data hunt was also brought in from `cloud/482c-lit-data` (all negative).
- **Queue (22:0xZ):** 470b (running) → 472 (le gate KAs) → 471b (re-levelled 471 arm) → 474 (481's patches + three
  def_day rehearsals) → 475 (q80 fix on H1). Builder commits come after 472's and 474's KAs.
- 22:05Z: **cloud 485 (October day-of review)**, brought in from `cloud/485-october-review`. The peer verified
  the core in code:
  - F1 critical: `build()` discards the driver's rc and only tests the h5ad's existence, and the driver exits when the
    h5ad already exists. So a rerun in the same work dir can ship the previous build (old q/factors). The dump's rc is
    ignored the same way.
  - F2 high: an out-of-bound secant step dies before the safe-build fallback, and lvl0 is never linked as final.
  - F3–F8 medium/small.
  - Decision: **not folded into tonight's 474**. Its `def_pipeline.py` patch changes the driver sha that
    `run_arm_generic.sh` pins (471b, next in line, would stop), and its def_day.sh edits overlap 481's in the same hunks.
  - **476 (09-27 daytime):** merge 485 onto the 481-patched files (474), run repro_fixes.py, update run_arm_generic's
    driver pin list, then rehearsals: a same-dir rerun (F1), a forced out-of-bound secant (F2), and the 465 retrace plus
    efp relaxed again.
  - **Interim rule for the day (DEF-RUNBOOK):** always rerun `def_day.sh` in a FRESH work dir. Never "lower q and
    rerun" in place.
- 22:14Z: **research queue after tonight** (our numbers 476–479):
  - **476**: merge 485's def_day fixes onto 474's, with rehearsals (see above).
  - **477**: 480's rider (+0.2 × the panel donors' same-target consensus log2 ratio on the delivered lfc; H1 −0.0030 raw).
    Price it beside 475 as a separate variable.
  - **478, the big one**: a delivery route for SIGNS on uncalled genes.
    - 480 puts 65 % of the nmae sign pot there (our accuracy 0.54).
    - 426 showed our class structure cannot carry sign flips on uncalled genes: prior-sign MDE over-fills the classes,
      and confident-only flips lose to the plain MDE.
    - Candidates: the exact-U mass-holder pass (`VCC_DE_MAG`, 423/429) extended to sign flips, or carrier-total moves
      paid by the ballast.
    - Source: the line's programme if 09-30 says it carries, else CLEAN4 (the analytic prize with CLEAN4 was −0.0128 vs
      −0.0098).
  - **473**: re-price 410d2 on 446b's bulk.
- 22:16Z: **experiment-number allocation across the three sessions Ahmet runs tonight:**
  - this coordinator: 468–479 (473, 476–479 planned);
  - the peer ("Rehearsal script bugs and re-level fixes"): 480–490, then 500+ for new work;
  - a third session ("Parallel agent work exploration"): 491–499. 491 = the recipient's own control cells as the
    DE block-programme source (fluctuation-response), analysis only. If positive, a preregistered per-context prior
    npz in 458d's format is handed to the coordinator.

  All three keep Mac processes ≤ 4 GB (or gated on vm_stat free+inactive ≥ 33 GiB) until 04:30Z. The coordinator's
  next block after 479 is to be agreed before use.
- 22:24Z: cloud 488 (scorer conformance, `cloud/488-scorer-conformance`) ran synthetic matrices through
  cell-eval2 0.16.0's own functions. The emitter's scorer model CONFORMS on every member: p/q/lfc to 1e-23, the control
  gate, midrank Wilcoxon, per-pert BH, reach ordering, fid n_pred, nmae fill, empty-set jac, pds. Three small
  differences, none worth an arm:
  - pooled vs mean per-cell CPM in the builder's render gates: 5–30 of ~10k genes, all at CPM 5.0–5.4;
  - the exact-U pass moves p through the tie term (known; 0 call flips);
  - a NaN jackknife when one cell holds every count (never on the scored path).
  Bring the branch in with the others.
- 22:26Z: cloud 489 (`rank_mag` review, `cloud/489-rankmag-review`, 42/42 checks): **no stop for the 09-28 zs
  riders.**
  - `rank_mag.apply` keeps every library, gene total, ballast row and full midrank R1 exactly. The module sha equals the
    builder pin `ed672832`.
  - The one real finding (429's "a repaired 1-count is frozen" needs a control at or below that cell's CPM) is measured
    not to touch A/B/C: all 24 zero-repaired cells have 9,112–15,635 zero-count control cells.
  - Its patch would change the module sha, so it is not for 09-28. Record it for a later build.
- 22:31Z: **473 (the 410d2 re-price) delegated** to the third session, analysis only, in
  `experiments/473-410d2-reprice/`. Method constraints were sent:
  - known answer: 414's emul410d2 rows and its PLOLO price, reproduced first;
  - base: 446b's actual bulk (dump + q 1.75 + ctx scales), priced as the increment over 446b;
  - variant: m clipped [1/2, 2], amplitude-matched, with 414's GO rule (energy ratio 0.97–1.03, max per-target ratio
    ≤ 1.5);
  - PLOLO folds with a target bootstrap; preregistered before the 446b-base number.
  The coordinator keeps the decision and any build. The third session's 491 (the line's own-control programme): fid-like
  +0.0055 [+0.001, +0.010] on H1, real but small, no arm.
- 22:56Z: **473 NO-GO** (third session, 2aecd167): the per-(target, gene) reliability arm on 446b's sharpened
  bulk prices −0.0008 re-levelled (q 1.75 already reweights the same genes). Closed in FINDINGS §2. Registry defects for
  October (Jiang/DLD1 log2FC space; panel-column holes in the dld1/myllia rows) went into DEF-RUNBOOK §1.
- 23:03Z: **492 (third session): per-target sharpening by the 410 key tercile, GO for a build (479).**
  - P2 (q 2.5 / 2.0 / 1.75 on the low / mid / high key tercile) prices +0.0034 overall re-levelled on PLOLO: 3/3 cancer
    folds positive with target CIs above 0, kolf2 +0.0011, pooled rel cos +0.95 % [+0.67, +1.27].
  - The permuted control loses direction. A/B/C render energy ×0.99, max per-target ×1.00. The ballast nmae cost is
    −0.0002 overall.
  - Board-calibrated expectation **+0.0013 … +0.0032**. All four preregistered GO items pass.
  - **479 (09-27 daytime):**
    1. builder knob `VCC_BULK_SHARPEN_Q=tercile:<low>,<mid>,<high>`, keyed on the pinned per-target key (`VCC_PERT_SCALE`)
       and its terciles over the panel;
    2. KA: `tercile:1.75,1.75,1.75` rebuilds 446b byte for byte, and unset is unchanged;
    3. arm 446b + P2, re-levelled (ntclevel); `def_day.sh` pass-through, since the day computes its own key so the
       tercile ports;
    4. a slot from 10-02.
  - 479's knob definition (from 492/RESULTS.md "Knob definition", e1b6c478):
    - key = the pinned `VCC_PERT_SCALE` artifact's clipped per-target w (non-carried targets 1.0);
    - edges = `np.quantile(w, [1/3, 2/3])`, tercile = `np.digitize(w, edges)`; on A/B/C 100/100/100 (q1 0.5234, q2 1.4605);
    - October guard: `low = (w < q1) | (w <= 0.3 + 1e-12)`, so a floor-heavy panel cannot empty LOW (identical on A/B/C);
    - rows keep their own gated non-panel norm after sharpening, and panel columns are untouched, as the global knob does.
- 23:06Z: **477 NOT GO** (third session, 1b2c9f71; analytic via 490's reconstruction). The rider's sign flips make
  reach seed-volatile. The magnitude-only piece is worth ≤ +0.001, so it is not a lever on its own. Closed in FINDINGS §2.
  - 492b (after the preregistered read, so information, not a rule): L30 = `tercile:3.0,1.75,1.75` prices +0.0023 (3/3
    folds, kolf2 +0.0017). Energy is exactly ×1.000 and the ballast unchanged; its permuted control loses outright.
    - P2 stays the build.
    - Consider P2 + L30 as a same-day dose pair: the same knob, one variable apart. It isolates the middle tercile and
      keeps a zero-energy-risk option for October.
    - Decide in 479's preregistration, before either is built.
  - 492c (after the read; information): raising the 410 key floor 0.3 → 0.5 stacks almost additively with per-target
    sharpening on PLOLO: F5P2 +0.0051 (kolf2 +0.0012), floor 0.5 alone +0.0015.
    - Caveat (420): the floor's pds value is entirely β-emulated (literal pds +0.0000). What is measured is mse
      (−0.002 raw).
    - Keep it a separate, later variable after P2/L30 read.
    - Link: 09-27 s1's LODO key rider (446bj2, `lodo415_lo0.5_hi3`) already carries a 0.5 floor with a different
      fit. Read its mse against 492c's floor-only +0.0015 when interpreting it.
  - 492d (after the read; information): the same q values by response-strength tercile (the unkeyed row norm, rank corr
    with the key 0.15–0.23), S_P2, prices +0.0039 (kolf2 +0.0032, cos +1.76 % [+1.14, +2.43]; A/B/C energy ×0.93, which
    the re-level refunds). Sharpening the STRONGEST third instead costs cos −10 % on every fold.
    - Mechanism: denoise the weak or unreliable rows.
    - The third session stopped here to avoid a forking-path hunt. **479 builds the preregistered P2**, optionally with
      L30 as its dose partner.
    - A strength-based or combined assignment needs its own later preregistration, after P2 reads.
- 23:30Z: **493 (third session): per-target sharpening on October-like panels** (438's h2c1 stack, 3 folds each on
  efp/efr/real; KA q 1.5 == 438's stored arm, |d| 0 on 12/12).
  - The mechanism ports: direction positive 9/9, board-emulated pds positive 9/9, controls lose.
  - **The floor guard is mandatory**: the K3 floor holds 18–41 % of targets there. Without the guard LOW empties (L30 is a
    no-op on efr). P2g (with the guard) ≥ P2 on every October panel, identical on A/B/C.
  - The strength-tercile assignment S_P2 (unkeyed stack row norm; continuous, no ties) is the strongest and most robust:
    direction +1.85 / +2.01 / +1.06 %, pds +0.005 … +0.007 on efp/efr/real, 9/9.
  - **479's plan, revised:**
    - one builder patch with both assignment modes: `VCC_BULK_SHARPEN_Q=tercile:...` (key, floor guard) and
      `stercile:...` (strength = unkeyed stack row norm);
    - KA: all-1.75 in both modes rebuilds 446b byte for byte; unset unchanged;
    - arms P2g and S_P2 (2.5/2.0/1.75), each re-levelled;
    - a same-day board pair, preregistered before either entry exists: one variable apart (the assignment rule), both
      read vs 446b;
    - October's day-of default follows the pair.
  - 479's `stercile:` definition (492/RESULTS.md "Knob definition, strength assignment", befe8fbd):
    - strength_t = the L2 norm of target t's UNKEYED stack row over ALL non-panel columns, computed once on the shared
      matrix (on A/B/C: `np.linalg.norm(dump["carriers"][:, nonpanel], axis=1)`);
    - edges = `np.quantile(strength, [1/3, 2/3])` with digitize: low (weakest) 2.5, mid 2.0, high 1.75. Continuous,
      no ties, no guard needed;
    - rows keep their gated non-panel norm after sharpening, and panel columns are untouched.
    The third session is calibrating its instrument on the board's q 1.5 → 1.75 step (433k → 446b, +0.0023) for 479's
    prediction.
  - 479's board-calibrated predictions (492/RESULTS.md "Board-calibrated predictions", d5debfe6; fitted on the board's two
    levelled q steps: k_pds 0.36, stable; k_X 0.61, pessimistic 0.12):
    - P2g vs 446b: pds +0.0027 raw, mse −0.0022; overall net of the ballast cost ~+0.0012 / +0.0009;
    - S_P2: pds +0.0026, mse −0.0042; net ~+0.0015 / +0.0009.
    Suggested bins, with pds the clean read (bulk members reproduce to 16 figures):
    - transfers at Δpds ≥ +0.0015 raw;
    - partial at 0 … +0.0015;
    - fails below 0;
    - nmae should move only by the ballast share (+0.0005 … +0.001 raw).
    Use these in 479's preregistration.
- 23:41Z: **478 analytic (third session, b4d0fe36): prior signs on uncalled genes are deliverable through the
  coarse MDE** (H1 +0.0005 … +0.0011 overall-equiv with CLEAN4). Next, after 09-30's read:
  - an emitter hook (the sign source for the coarse path's non-block genes, with a confidence threshold), default off,
    with a KA;
  - an H1 harness run for leaked calls;
  - the source is Arc's programme if 458d carries, else CLEAN4 confident-only.
  - The coordinator needs a new experiment-number block after 479. Agree it with the peer and the third session.
- 23:42Z: numbering update:
  - the coordinator also holds **495–499** (given by the third session, which used 494 and moves its new work to 540+);
  - 510–529 is requested from the peer, not yet confirmed (500-cloud-bulk-review exists in its range).
  - 478 addendum (1fbddeaf): the recipient-side programme (491c) as the uncalled-sign source gains nmae (−0.0047 at
    |P − 0.5| ≥ 0.2) but costs reach on every seed. Net ~+0.001, no better than CLEAN4-confident.
  - So the hook's source is CLEAN4 (|P − 0.5| ≥ 0.05) or Arc's programme, decided by 09-30.
- 23:46Z: third session:
  - 494: per-target q80 protection by design-strength tercile is null on H1 (≤ +0.0003; the reverse control −0.0005).
    Closed: q80 is near-optimal per target.
  - 492f: under the strength assignment, the strong third stays at q 1.75. Lowering it to 1.5/1.25 loses direction on the
    cancer folds. S_P2 = (2.5, 2.0, 1.75) stands for 479.
- 23:52Z: 493c (third session, exploratory): on 438's October panels, floor 0.5 stacks with per-target q:
  - F5+P2g direction +0.74 / +2.06 / +1.15 %, F5+S_P2 +2.21 / +3.14 / +1.03 % (efp / efr / real), 8/9 panel-folds
    positive;
  - floor alone +0.3 / +1.1 / +0.2 %, with literal pds ~0, as on A/B/C.
  - **October bulk sequence:** the P2g/S_P2 pair reads (479), then the winning assignment, then the floor as its own
    preregistered step.
- 23:55Z: **471b fails its ρ check** (A 1.005, B 0.996, C 1.026 vs 0.9793: past the cliff in A and C; never
  ship `37d1919d`). The energy-preserving key had already kept the total energy. **471c**, re-levelled from ntclevel's
  EV-scale (A 1.0993, B 1.1557, C 1.1587), is queued after 475 (chain471c pid 67526).
  - **474:** F2 KA PASS (A/B/C Cp46 lines identical, --ka PASS); F7 KA PASS (465's sameline.npz byte for byte); F1/F5
    applied; rehearsal (a) PASS (efp stops at setup with the explicit CD4-gate message). (b) the 465 retrace is running,
    then (c) efp relaxed.
  - 09-27 s1 uploads at 00:01Z. Watchers are on both captures.

### 2026-09-27 01:08Z — the 09-27 reads; 09-28 armed; overnight queue
- **s1 446bj2 = 0.2722766, NEW BEST** (+0.0022338 vs 446b). The LODO key transfers: mse −0.0145 raw, pds +0.0002, DE
  pinned. key_ok, so October runs **`DEF_DAY_LODO=1`**.
- **s2 446br2 = 0.2704006** (+0.0004). The repair floor −2 is partial (nmae −0.0014 raw < the −0.002 bin), so **−5
  stays** and floor_ok is False. 463 is NOT launched.
- **09-28 armed** by chain0928's rule (key_ok, not floor_ok): s1 **453dr** (pid 72218; the October path rehearsal on
  A/B/C, vs 446b), s2 **446bj2z** (pid 72220; the zs rider on the LODO base, vs 446bj2).
- 09-29 (chain0929: 457dr, 446bd2), 09-30 (chain0930c: 458d, 458g) and 10-01 (chain1001b; key_ok, so 460dr is appended)
  stay armed.
- **Overnight queue:**
  - 474 (481's October patches): F2, F7 and (a) PASS; (b) the 465 retrace is running; (c) efp relaxed next;
  - then 475 (the q80 fix on H1; 490 predicts a fail);
  - then 471c (the re-levelled union-swap arm).
  - Builder commits come after each KA set. 474's three patched files are committed after (c) passes.
- **Next builds (09-27 daytime):**
  - 479, per-target sharpening P2g + S_P2 (definitions and predictions in 492/RESULTS.md). The incumbent is now
    446bj2 (the LODO key), so the tercile key is the LODO key. Price the stack on it before building (the third session
    can do it).
  - 476, 485's def_day fixes merged onto 474, with rehearsals.
  - The uncalled-sign emitter hook (478's route) after 09-30.
- 02:38Z: **474 ALL PASS and committed (847cced5)**: F2/F7 KAs; efp strict stops with the CD4 message; 465 retraced
  byte for byte with postbuild PASS; efp relaxed DONE with ρ ≤ 1 and postbuild PASS. DEF-RUNBOOK's one-page day documents
  `DEF_DAY_CD4_THRESHOLDS=relaxed`. 475 is next, then 471c.
- 04:34Z: **492g (third session, 6f7cf8d5): 479 is GO on 446bj2 (the LODO stack).** Calibrated vs 446bj2 (pds 0.34 ×
  literal; mse k_X 0.61 / 0.12; net of the ballast nmae):
  - **P2g +0.0012 / +0.0008**: Δpds +0.0026, Δnmae +0.0010 raw; PLOLO k 1 +0.0035, 3/3 folds, CI [+0.0023, +0.0049];
  - **S_P2 +0.0013 / +0.0008**: energy ×0.94, refunded by the re-level; PLOLO +0.0038, kolf2 +0.0032;
  - CTRL nets +0.0000; the floor and the per-target q do not overlap.
  - Suggested bins per arm vs 446bj2:
    - transfers at overall ≥ +0.0008 with Δpds ≥ +0.0015;
    - partial at +0.0003 … +0.0008, or pds +0.0005 … +0.0015;
    - fails below that;
    - nmae +0.0006 … +0.0010 raw expected; above +0.0015 means the ballast is under-priced.
  - **479's build plan (09-27):**
    - the knob reads everything from the loaded dump. `tercile:` uses r_t = ‖bulk_t‖/‖carriers_t‖ over non-panel
      columns (∝ the key's w for either key), with the floor guard = targets at min r → LOW. `stercile:` uses
      ‖carriers_t‖;
    - KA: tercile/stercile at all-1.75 rebuild 446bj2 byte for byte; unset unchanged;
    - arms P2g and S_P2 on 446bj2, each built, ntclevel'd and re-levelled to 0.9793;
    - a same-day pair preregistered vs 446bj2.
- 04:40Z: third session:
  - **493h (the tie rule):**
    - 479 on A/B/C keeps the plain rule as priced. The patch's tie rule is a no-op there (`tie_to_mid` 0; the 6
      non-carried targets are MID by value).
    - October: stercile is the default candidate (9/9 folds positive, tie-free). If the key tercile is ever used there,
      implement carried-only cut points with non-carried → MID, reading `carried` from the dump's key artifact. Do not
      infer ties from r: a cap block (w = 3) also exists.
  - **540 (preregistered): the LODO key's cap of 3 binds and costs.**
    - Cap 6 prices **+0.0025 overall net vs 446bj2** (k_X 1.31, the key class the board just measured; +0.0019 at k 1).
    - Direction +2.6 % (3/3 cancer folds), kolf2 +8.3 %. All four preregistered rules pass. A/B/C energy ≤ 1.6 %, max
      row ×2.04.
    - Orthogonal to 479: capped targets are all in the HIGH tercile.
    - Coordinator build **495**:
      1. the LODO key artifact at (0.5, 6) (`def_lodo_key.py` with hi 6), registered in the builder's pert-scale
         registry with a KA;
      2. a new stack dump;
      3. the arm 446bj2 + cap 6, re-levelled;
      4. October: `def_lodo_key.py` needs the hi parameter plumbed through `DEF_DAY_LODO`.
- 04:48Z: **540o + 540c (third session, da71f6e4)**:
  - **October:** cap 6 helps more on 438's October-like panels than on A/B/C (direction +12.7 / +5.1 / +8.2 % on
    efp/efr/real, 9/9). Higher caps win only there (A/B/C loses from 8 up), so 6 is the robust setting. After 495's board
    read: `DEF_DAY_LODO` passes `--hi 6`.
  - **Stacking (A/B/C folds):** 479 × 495 stack cleanly.
    - P2g's increment on cap 6 is unchanged (+0.0012 / +0.0009); S_P2's shrinks to +0.0010.
    - Totals vs 446bj2: cap 6 + P2g **+0.0031 … +0.0037**, cap 6 + S_P2 +0.0029 … +0.0035.
    - On October panels the parts multiply cleanly.
  - **Board plan from 10-02:**
    1. the 479 pair (P2g + S_P2, same day, vs 446bj2);
    2. 495 (cap 6) with 471c;
    3. the combined arm (cap 6 + the winning per-target assignment).

    Then the October configuration rehearsal carries every adopted switch. Each day's rule is fixed before the previous
    day's entries exist.
- 05:02Z: **543 (third session, c8ecd6f2): the October form of the cap is a percentile.**
  - hi = percentile(key/median over carried targets, 96.647). On A/B/C that is exactly cap 6 (reproduces 495's
    `dc374281` to 4e-15), so 495's read calibrates it.
  - On 438's heavier-tailed October-like panels it caps at 8–16 and adds +5 … +6 direction points over fixed cap 6
    (9/9 folds). p99 starts losing A/B/C's k562 fold.
  - **If 495 reads "transfers":** `def_lodo_key.py` gets `--hi-pct 96.647` (KA: A/B/C == `dc374281`), and
    `DEF_DAY_LODO` uses it. If 495 fails, neither form goes to October.
  - Caveat: the October-like tails may partly be construction. The rule keeps A/B/C fixed.
  - Closed (third session): 541, robust LODO aggregates, do not beat the mean; 542, kolf2 as an extra truth proxy in the
    key, hurts (−1.1 … −2.4 %).
- 05:13Z: **546 (third session): the October bulk candidate as a whole** is positive on all 9 October-like folds.
  - The LODO key with the p96.647 percentile cap: +19.4 / +11.6 / +14.1 % direction vs today's October default (cap 3).
  - The per-target arm on top still adds on every fold (S_P2 +0.7 / +2.6 / +0.2 %).
  - Totals with S_P2: +20.3 / +14.4 / +14.3 %.
  - **Candidate October bulk configuration, gated on the A/B/C board reads of 479's pair and 495:** `DEF_DAY_LODO=1`
    with `def_lodo_key.py`'s cap at p96.647 (= cap 6 on A/B/C), plus 479's winning per-target arm (S_P2 by default,
    stercile, tie-free).
  - Implementation after the reads:
    - `def_lodo_key.py --hi-pct`;
    - a `def_day.sh` pass-through for the per-target sharpening (`VCC_BULK_SHARPEN_Q=stercile:...` into the lvl builds);
    - a full rehearsal on the relabelled A/B/C bundle and on efr.

### 2026-09-27 05:28Z — the 10-02 rule, FIXED before any 09-28 entry exists
- **s1 479p2g (P2g), s2 479sp2 (S_P2)**: the per-target sharpening pair vs 446bj2, bins and the October pick rule in
  `experiments/479-pertarget-sharpen/PREREGISTRATION.md` (b692c2c5). Builds: chain479 (KAs) → chain479b (the two arms,
  each built at 446bj2's factors, ntclevel'd and rebuilt at the EV-scale factors).
- Fallbacks in order if an arm is not READY: **495** (cap-6 key, chain495 after 479b; own PREREGISTRATION), then **471c**.
- READY = package present, config diff PASS, CLIFF OK, novelty rc 0, novelty re-checked, every context's ntclevel ρ within
  0.5 % of 0.9793. `chain1002.sh` (pid 98741) applies it at 10-01 20:00Z and arms `arm1002.sh`
  s1/s2. The rule reads no 10-01 entry. Dry run now: 471c READY (ρ 0.9812/0.9804/0.9745), the others not built yet.
- Board plan after that (fixed later, each day before the previous day's entries): 10-03 495 + 471c; 10-04 the combined
  arm (cap 6 + the winning per-target assignment).

### 2026-09-27 06:07Z — the 10-03 rule, FIXED before any 09-28 entry exists
- **s1 495 (LODO key at cap 6), s2 471c (h1_2025_full union swap, re-levelled)**, both vs 446bj2, bins in their own
  preregistrations (495: 706ad5f0; 471: the 471b/471c addenda). Both are bulk questions that do not depend on 10-02's
  reads.
- Fallback if either is not READY: **497dr** (the October path + percentile cap + S_P2, vs 460dr; `gate497.py --write`
  at arming time). A package already scored (e.g. 471c if 10-02 had to fall back to it) fails the novelty re-check and is
  skipped.
- `chain1003.sh` (pid 5183) applies it at 10-02 20:00Z and arms `arm1003.sh`. Dry run now: 471c READY;
  495 and 497dr not built yet.
- 10-04 is NOT fixed yet: it waits for 09-28's read of 453dr (the October path itself) and is fixed before any 10-03
  entry exists. Candidates: 498p2g / 498sp2 (cap 6 + the pair-winning per-target form; chain498) and 497dr.
- Build queue as armed (strictly serial): chain479 (KAs) → chain479b (the 10-02 pair) → chain495 → chain476 (485's
  fixes merged onto 474, four rehearsals) → chain497 (October pass-throughs, inert KA, 497dr) → chain498 (combined arms).
  The coordinator commits the builder after chain479 and after chain495, and 476's five files after chain476; each
  later chain waits for its predecessor's process to exit and for HEAD to match.
- 06:17Z: **the build queue no longer needs a live session.** `scripts/commit_after_chain.sh` committers are
  armed behind chain479 (pid 6766), chain495 (6767), chain476 (6768) and chain497 (6769). Each waits for its chain's
  process to exit, refuses on RESTORED or any missing known-answer line, checks exact file shas (d09a986c; b4975024;
  476's e3247feb / 346495d7 / c46beb82 / bd10f8f2; 497's fb8f4264 / 60b2a777), then commits only its paths with the
  KA lines in the message. Logs: `<chain>.commit.log` beside each chain log. Do not commit these files by hand.
- 06:23Z: **-m names now come from `docs/MODEL-NAMES.md`'s queue** (AGENTS §14; I had baked invented
  lab-object names into chain1002/chain1003 this morning). Claimed with `scripts/next_model_name.py --claim`: 479p2g #1,
  479sp2 #2, 495 #3, 471c #4, 497dr #5 (`docs/model-names-claimed.tsv`). chain1002 (pid 7821) and chain1003 (pid 7822)
  were killed and re-armed with the same logic and the queue names. The build chains' h5ad names are local file names
  and stay as they are.

### 2026-09-27 06:27Z — the 10-04 rule, FIXED before any 09-28 entry exists (`498-combined-bulk/rule1004.py`)
Inputs, read at 10-03 20:00Z by `chain1004.sh` (pid 8621): 453dr (09-28), the 479 pair (10-02),
495 (10-03). All three are entries that do not exist yet.
- P = 479's pair pick (read1002's rule). C = 495's bin vs 446bj2. O = the October path's status from 453dr vs 446b
  (BAD if bulk < −0.0060, |DE| > 0.009, or overall < −0.02, per 453's bins).
- **s1 = 498<P>** (cap-6 dump + the picked per-target form, re-levelled) if P ≠ NONE and C ≠ FAILS.
- **s2 = 497dr** (October path + percentile cap + S_P2, vs 460dr) if O ≠ BAD, C ≠ FAILS and P = S_P2.
- Empty slots take the first READY of: 498<other form> (if C ≠ FAILS and its own 479 bin was ≥ PARTIAL), 460dr, 461dr,
  469dr, 458e. Anything already scored fails the novelty re-check.
- -m: 498p2g #6, 498sp2 #7, 497dr #5 from the queue; the four fallbacks keep the names chain1001b baked on 09-26.
- Synthetic-entry test of the rule (5 scenarios) and a dry run of the chain on today's state (→ 460dr, 461dr) pass.

### HANDOFF CHECKLIST — refreshed 2026-09-27 06:31Z (supersedes 09-26 17:17Z where they differ)
**Armed, no session needed** (pids as of this entry; after a reboot re-arm with the same scripts, AGENTS §6.2):
- Slots: 09-28 `456/arm0928.sh` s1 453dr (72218), s2 446bj2z (72220) · 09-29 `457/chain0929.sh` (73089) ·
  09-30 `459/chain0930c.sh` (88484) · 10-01 `464/chain1001b.sh` (88678) · 10-02 `479/chain1002.sh` (7821) ·
  10-03 `495/chain1003.sh` (7822) · 10-04 `498/chain1004.sh` (8621, rule1004.py).
- Serial build queue, each link waiting for its predecessor's PROCESS to exit:
  chain479 KAs (91855; KA-S, KA-T PASS) → chain479b (97871; 10-02's pair) → chain495 (98390) → chain476 (99563) →
  chain497 (4693) → chain498 (4922) → chain496p (9132; the hook into production + 496a) → chain497e (9293; all October
  switches on efr).
- Committers (`scripts/commit_after_chain.sh`) behind chain479 (6766), chain495 (6767), chain476 (6768), chain497
  (6769), chain496p (9148): exact-sha and known-answer gated. **Do not commit those files by hand.**
- Watchers: board_snapshot (2910), pisces_watch (2911); caffeinate to 10-05 00:00Z (9424).
- Cloud: 499 (Cursor agent `bc-0c5cd3a8`, grok-4.7-xhigh) reviews 476/497/496 on branch `cloud/499-october-code-review`.
  Read its RESULTS.md when it finishes; nothing it finds is applied without a local reproduction.

**Manual, by whoever is on:**
1. Board reads as they land (each reader exists): 09-28 `456/read0928.py`; 09-30 `468/read468.py`; 10-02
   `479/read1002.py`; 10-03 `495/read1003.py`. Write BOARD-RESULT files, BOARD-LEDGER, FINDINGS §0, and the
   DEF-RUNBOOK decisions (LODO cap, per-target sharpening, programme source, gate).
2. Fix the 10-05 rule before any 10-04 entry exists (candidates: 496a, the DE-only hook arm on 446bj2, bins in
   496/PREREGISTRATION.md; whatever of 10-04's lists is unused).
3. After 497e: record the all-switches efr rehearsal in 497/RESULTS.md.
4. The manuscript artifact (v29) is due a republish with C-73 (`build_html.py`, same URL).
- 06:57Z: **499 (cloud review) caught a blocker; the queue is re-ordered.** The driver's arm assert did
  `float()` on 479's per-target form, so every per-target build (the 10-02 pair, 498, 497dr, 497e) would have exited rc 1
  before its postbuild (KA-T/KA-ST show it; their hashes matched only because the KA runner ignored the rc). Driver fixed
  (44714b9a, cb9e432a; run_arm_generic pins it). Details and the other fixes: `experiments/499-october-code-review/RESULTS.md`.
  - Queue now: chain479 DONE (KA-S/T/ST PASS; builder d09a986c committed 24dd26ae by its committer) → **chain495
    (running; builder b4975024 patched)** → **chain479b (pid 14854; 10-02's pair)** → chain476 → chain497 → chain498 →
    chain496p → chain497e.
  - Committers (spec files; argv no longer matches any gate): 495 (10763), 476 (11292), 497 (11634), 496p (11635).
  - Cloud: 501 (Cursor `bc-89f0f1e2`) is reviewing the slot chains, gates, committer and 479's tercile code.
- 06:48Z: 496z, the hook stacks exactly on the zs rider on H1 (nmae −0.00146 on either base; reach identical; 3/3 seeds
  deterministic).

### 2026-09-27 06:59Z — the 10-05 rule, FIXED before any 09-28 entry exists (`496-uncalled-sign-hook/rule1005.py`)
- **s1 = 496a** (446bj2 + the uncalled-gene sign hook, DE-only; board bins in `496/PREREGISTRATION-board.md`; -m #8).
- **s2** = the first READY, not-yet-scored of: 497dr (only if C ≠ FAILS and P = S_P2 and O ≠ BAD, rule1004's P/C/O),
  498sp2 / 498p2g (only if C ≠ FAILS and that form's own 479 bin was ≥ PARTIAL), 469dr, 461dr, 460dr, 458e. If 496a is not
  READY, s1 takes the first of that list too.
- `chain1005.sh` (pid 15700) applies it at 10-04 20:00Z. Dry run on today's state: 496a not built yet →
  s1 469dr, s2 461dr.
- 07:42Z: **experiment numbers (confirmed with the peer session "Rehearsal script bugs and re-level fixes").**
  The coordinator's range is 510-529, then 530+. The peer's is 500-509, the third session's 491-499. Already taken in those
  ranges by the coordinator before the confirmation, and live: 495-499, 501, 502 (the peer was told, and asked to relay to
  the third session). Next coordinator experiment: 510.
- 07:46Z: **487 (the peer's Opus review of the October recipe) processed, and two October decisions changed
  before any 09-28 entry exists:**
  1. **Same-line rows:** `def_day.sh` args 5–6 = `0.5` / `noprog`, or `0.5` / `""` with the programme source (w 0.5, no pt).
     The likely October matches are banks already in our mix (lift ×1.16–1.29 vs Pan-on-H1's ×4.37). DEF-RUNBOOK updated.
  2. **Prior gate:** a NULL 458g read keeps `le1`; only HARM closes it (468 addendum; `read468.py`). `rule1001b.py`
     is amended to match, so P ∈ {positive, positive_reach, partial} with G ∈ {pos, null} picks 469dr, the gated October
     path, for 10-01. Tested on synthetic reads: null → 469dr, pos → 469dr, harm → 461dr, 458e.
  - Rehearsal 510 goes into the build queue after 497e.
- 07:54Z: **process map after the 499/501/487 re-arms** (this supersedes the pids above):
  - build queue: chain495 98390 (running) → chain479b 25674 →
    chain476 25695 → chain497 25715 →
    chain498 25739 → chain496p 25759 →
    chain497e 25781 → chain510 25804 →
    chain511 26434.
  - committers (spec files): 495, 476, 497, 496p, 511. Slot chains: 0929 73089, 0930c 88484,
    1001b 88678, 1002 7821, 1003 7822,
    1004 18096, 1005 18097. 09-28 waiters 72218, 72220 (unchanged).
- 08:45Z: **495 round 2 left C at ρ 0.9627 (−1.7 %).** A levelling round 3 (`495/relevel3/chain495.sh`, pid
  37277) is inserted after chain479b and before chain476, at A 1.0680 B 1.1215 C 1.1235. It rebuilds under the same tag, name and
  stem, so chain1003, chain498 and read1003 are unchanged. 495's builder (b4975024) was committed at 08:44Z; chain479b has been
  building since 08:44Z (a log line at 08:44:45Z saying it was killed is corrected in the next line: it was not).
- 11:40Z: **495 round 3 READY** (package `ce0a41af`, ρ A 0.9805 B 0.9812 C 0.9788, all within 0.5 % of 0.9793; config
  diff PASS, CLIFF OK, novelty rc 0). chain1003's readiness test dry-run by hand: PASS, so 10-03 s1 = 495. **chain476
  started 11:26:56Z**: the merged day-of files are applied (def_day e3247feb, helpers 346495d7, driver 9b2c7eca, relevel
  bd10f8f2; repro_fixes PASS) and (a) efp strict PASSED (stops at setup with the CD4-gate message).
  - Disk: the arm runner keeps each ~12 GB h5ad and gates at >= 20 GiB free (then gives up after 2 h). Free was 67 GiB
    against ~80 GB queued, so the packaged dr h5ads and their run-dir copies were swept (now 143.6 GiB free;
    `docs/STORAGE-SWEEP-2026-09-27.md`). **Keep 497dr's h5ad through 10-05**: gate497.py re-reads it at each readiness check.
  - `468/read468.py`'s October line now follows the 09-27 addendum (NULL gate keeps le1), tested on synthetic entries.
    `rule1001b.py` already did. 10-02/10-03 readers checked against their preregistrations: consistent.
  - Process map (pids): build queue chain476 25695 (running) → chain497 25715 → chain498 25739 → chain496p 25759 →
    chain497e 25781 → chain510 25804 → chain511 26434; committers 476 25862, 497 25863, 496p 25864, 511 26450.
    Slots: 09-28 72218/72220, 0929 73089, 0930c 88484, 1001b 88678, 1002 7821, 1003 7822, 1004 18096, 1005 18097.
    Watchers 2910/2911; caffeinate 9424 (to 10-05) and others.
  - Candidate follow-up, not built: if 09-30 reads the programme positive AND 10-05 reads the hook positive, the October
    hook should probably take its sign from the active block prior (Arc's programme) instead of CLEAN4
    (`patch_builder496.py` hard-wires `clean4`). Only the board can price that (H1 is Arc's own line).
- 11:58Z: **queue of the next three (AGENTS §11.1b), after the automatic build queue drains (~02:30Z 09-28):**
  1. **09-28 reads (~00:45Z):** `456/read0928.py dress <453dr> <446b> zs <446bj2z> <446bj2>`; then BOARD-RESULT files,
     BOARD-LEDGER, FINDINGS §0, and `scripts/october_recipe.py decide de_mag <dn_k250_zs|off> --evidence ...`. If zs
     TRANSFERS, apply 502's `DEF_DAY_DE_MAG` pass-through with its inert known answer (a def_day retrace) before 10-04.
  2. **The recipe dr for 10-06** (design after 09-30's read, rule fixed before 10-05 00:00Z): the October path on 180's
     relabelled A/B/C bundle with exactly `october_recipe.py command --allow-pending`'s env as of 10-05 02:00Z, relabelled
     and gated like 497dr. It is the one end-to-end A/B/C proof of the exact final configuration (the DE switches and the
     bulk switches have only been read separately).
  3. **Hook source follows the programme** (only if 09-30 AND 10-05 both read positive; low value, ±0.0005): not built.
  - 512 (cloud lever hunt) is done and recorded: no structural or public lever closes the DE gap; its two proposals are
    not queued (reasons in its RESULTS). Panel-gene coordinates are populated in production rows (audited).
  - `docs/october-decisions.json` + `scripts/october_recipe.py`: record every read's decision there as it lands.

### 2026-09-28 12:20Z — reboot recovery (new coordinator session; the three 09-27 sessions are gone)
- **The machine rebooted ~20:23Z 09-27** (uptime 15:45 at 12:08Z). Every waiter, chain, committer and watcher died,
  including both 09-28 slot waiters (72218/72220) before the 00:00Z reset. **The 09-28 slots were unspent at 12:08Z.**
  Builder, def_day, helpers, driver and build_submission all equalled HEAD (no pinned builder left behind).
- **Re-armed (AGENTS §6.2 order):**
  - 09-28 slots with chain0928's exact arguments: s1 453dr (pid 55341; submitting since 12:11:29Z), s2 446bj2z
    (pid 55362). `arm0928.sh`'s slot-2 cap moved 12:00Z → 23:30Z (slot 1's loop was already satisfied).
  - Slot chains, same scripts: 0929 55686, 0930c 55698, 1001b 55709, 1002 55722, 1003 55731,
    1004 55738, 1005 55746. Every candidate package for 09-29 … 10-03 is on disk.
  - Build queue: chain498 56031 (building 498p2gr1 since 12:13Z; the interrupted dir is
    `build_498p2gr1.aborted-20260927T2021Z`) → chain496p 56131 → chain497e 56152 → chain510 56166 → chain511 56184;
    committers 496 and 511 re-armed from their spec files.
  - Watchers: board_snapshot, pisces_watch; caffeinate to 11-06 00:00Z.
- **497dr is READY** (package `6ed67b2a`, gates497 ALL PASS, config diff vs 460dr exactly bulk_ctx_scale,
  bulk_load_sha256, bulk_sharpen_pt, bulk_sharpen_q). Its 20:20:39Z STOP was NOT the reboot: `def_postbuild.py`'s two
  "axis moved off A/B/C's" rows fail BY CONSTRUCTION on 180's relabelled bundle (A/B/C's own 300 targets give A/B/C's
  axis). 453dr/460dr/461dr/469dr predate the postbuild step that 476 merged (481 F5), so 497dr was the first dr through it.
  - Fix: the rows key on the TARGET set; new contexts on A/B/C's targets assert the axis IS A/B/C's. This also removes a
    latent October-day failure if D/E/F ever reused A/B/C's targets (the site says D/E/F have "their own perturbation list").
  - Proof: the unpatched script, re-run against the reconstructed pinned builder (pin dry-run diff on b4975024 →
    1f23b265 = the run's recorded pin), reproduced the run's 43 verdicts exactly; patched: 42/42. Old-vs-new on
    446bj2/data/challenge and 474efp/ef-bundle: identical verdict lists. Relabel, prep, novelty and gate497 then ran as
    chain497's (b) would (`497/finish497dr.sh`).
- **The official pages (read 12:30Z):** Oct 22 releases only D/E/F's controls plus "their own perturbation list"; no
  A/B/C truth release is mentioned anywhere (rules, About the Data, FAQ). No free same-lab bank arrives in October.
- 13:05Z: **09-28 reads.** 453dr 0.2663613: the October path is VALIDATED (bulk −0.0031 in band, DE pinned; O = OK for
  rule1004/1005). **446bj2z 0.2728260 NEW BEST**: zs TRANSFERS (nmae −0.00218 raw, the other five bit-identical);
  `october_recipe.py decide de_mag dn_k250_zs` done. Records: 453/BOARD-RESULT-453dr.md, 456/BOARD-RESULT-446bj2z.md.
- 13:09Z: **chain502 armed** (pid 68231): after chain511 exits and the 496/511 committers land, it applies 502's
  DEF_DAY_DE_MAG / DEF_DAY_MDE_SIGN pass-throughs (dry-applied today == repro/*.new; test_riders502 9/9), runs the inert
  KA (efp relaxed == ac49d319) and a rider-on rehearsal; committer from `502/commit502.spec`. The 10-06 recipe dr needs it.
- 13:20Z: **parallel research launched (Ahmet: local Opus 5.5 max + cloud Grok 4.7 xhigh).** Local workers (MANDATE.md in
  each dir; memory ≤ 4 GB/process until 02:00Z 09-29): 516 target-specific DE sign sources (pid 68842), 517 lab/Flex programme
  + arms for every 09-30 outcome, deadline 09-30 00:00Z (pid 68874), 518 the bulk mse gap (pid 68908). Cloud: 519
  October code review of the 09-28 changes (bc-3ade4f0d, cloud/519-october-review), 520 D/E/F line + same-line data hunt
  (bc-aa71e11f, cloud/520-def-lines). Nothing a worker returns is applied without a local reproduction.
- 14:35Z: **worker/agent results and re-arms since 13:20Z.**
  - 520 (cloud) done: no public source names D/E/F; high-coverage banks are the least likely matches (FINDINGS §0).
  - 519 (cloud) done: five defects, all reproduced locally, all applied with tests (commit 2d1fc214): the pinned A/B/C target
    digest in `def_postbuild.py` (the morning's fix would have failed every Oct 22 build), fail-closed chain502, gate497
    reads `postbuild.json`, chains 1003–1005 reject a failed gate, chain1001b date floor, the resolver prints refused commands
    commented out. Re-armed: chain502 81230 (+ committer), chain1001b 81703, chain1003 81553,
    chain1004 81568, chain1005 81583.
  - 516 (Opus) done, NEGATIVE: no target-specific public sign source transfers (oracle per-target ceiling +0.030).
  - 518 (Opus) done: the bulk mse identity and a +0.018 … +0.023 magnitude pot; its strength key is a 495 contingency.
  - **521 (Opus, pid 89379) launched**: 518's strength key predicts H1 n_real at r +0.60 vs 192's +0.24, so
    per-target call-set sizing (oracle +0.036) is re-priced on all four DE members at the board's n_real regime.
  - Still running: 517 (programme arms for 09-30, deadline 09-30 00:00Z), the build queue (498sp2 now).
  - 09-29 reads: `457-cd4-half-def/read0929.py cd4 <457dr> dld1 <446bd2>` (written and self-tested today).
- 15:10Z: **queue for 10-06 onward (rules fixed later, each before the previous day's entries exist):**
  1. **10-06 s1 = the October recipe dr** (previous coordinator's plan): `def_day.sh` on 180's relabelled bundle with exactly
     `october_recipe.py command --allow-pending`'s env as of 10-05 02:00Z (needs chain502's pass-through), 453's day operator
     cloned, relabelled and gated like 497dr; board-read vs 453dr (the sum of the adopted switches on the October path).
  2. **The production "all switches" A/B/C arm** (446bj2z + the 10-02/10-03 bulk picks + the 10-05 hook if adopted + the
     09-30 programme if adopted): the likely new best; with (1) on the same day it measures the October path's cost at the
     full recipe (453dr measured it at the plain recipe: −0.0037).
  3. **Information slot: the global DE call-set size at K 1,600** (521 §8.1): DE-only, pds/mse pinned; reads the board's
     n_real regime through fid's coverage term (192's fit says ~0, the jac level says +0.0013 … +0.0049); informs D/E/F sizing.
     Needs bins drawn from 521's tables before arming.
  4. Contingencies: 518's strength-key arm S1 only if 495 reads FAILS on 10-03; 517's programme arms per the 09-30 read.
  - 498sp2 round 3 (relevel3/chain498.sh, C 1.1149) is armed after chain502 so the 10-04 rule can use 498sp2.
- 19:45Z: **evening state.**
  - Both October rehearsals PASS: 497e (every candidate switch on efr) and 510 (the decided recipe, forced in-mix same-line
    call, w 0.5 / programme on / gate le1 / pct cap / S_P2): converged, ρ ≤ 1, postbuild PASS, audit 0 (RESULTS in each).
  - 496p landed (8960c102: emitter 1658ae2c + builder 24b95dec; KA-H1, KA-1 PASS); 496a READY for 10-05.
  - Workers: 518, 521, 522 done (FINDINGS §0). 517 running; its pre-read note became 458's jac-leak addendum (49d63258);
    its own H1 confirmation runs (`517/drive517.sh`) start after 02:00Z 09-29 under the memory gate.
  - Queue: chain511 (running) → chain502 → relevel_sp2 (498sp2 C round 3) → chain523 (call-set knob + K pair).
  - **October pass-through bundle (after the 10-05 reads, ONE def_day patch + ONE inert KA + one rider-on rehearsal):**
    the final-phase feedback knobs (ρ* scale, CD4 dose; FINAL-PHASE-FEEDBACK-PLAN rows 2 and 4), the 522 strength rider
    (if 495 transfers), and `DEF_DAY_RC_CALLSET_CPM` (if 523's K pair reads positive). Each default off, bit-identical.

### 2026-09-29 01:10Z — the 10-06 rule, FIXED before any 09-29 entry exists (`524-october-recipe-dr/chain1006.sh`)
- **s1 = 524dr**, the October recipe dress rehearsal: `recipe_dr524.sh` (pid 40679) waits for ≥ 10-05 02:00Z, no other
  build, def_day/helpers/builder/emitter == HEAD and ≥ 30 GiB free+inactive; resolves `october_recipe.py command
  --allow-pending --json` on 180's def-bundle (switches still PENDING then run OFF and are listed), clones 453's operator,
  runs def_day.sh with exactly the resolved env/args, relabels, preps, novelty, `gate524.py --write` (diff vs 453dr ⊆ the
  keys the ON switches explain; every ON switch visible in the meta; smoke-tested on 497dr's run: PASS with the matching
  recipe, FAIL on a recipe claiming zs, FAIL on one omitting S_P2). Read vs 453dr with 524/PREREGISTRATION.md's bins
  (Σ of each ON switch's own A/B/C board read ± 0.0020; Σ = +0.0027832 with the two switches decided today).
- **s2 = 446bj2zk16** (523's K 1,600 vs 446bj2z). If 524dr is not READY by the rule's run, s1 = 446bj2zk16 and
  s2 = 446bj2zk20: the whole preregistered K pair on one day.
- `chain1006.sh` (pid 40686) runs at 10-05 20:00Z after recipe_dr524's process exits (23:30Z at the latest) and arms
  `arm1006.sh` (arm1005.sh with the dates moved). Names claimed from the queue (#19–#21). 10-07 (not yet fixed): the other
  K arm or K20, and the production all-switches arm built after the 10-05 read.
- **Model-name queue:** after these three claims only #22–#23 remain; extend it from Ahmet's Notion quotes (same screen)
  before the next chain is written.

### 2026-09-29 10:46Z — 09-29 state
- 09-29 s1 **457dr = 0.2695096** (+0.0031 vs 453dr; bin ≥ +0.0010; October keeps CD4 ×1). s2 446bd2 uploaded 03:08Z →
  ~10:40Z (network crawl), now scoring; read with `read0929.py dld1 <entry>`.
- 525 (cloud) negative: public DE designs lose to ours on all four DE members; `cellpred` is private.
- 526 (cloud, bc-55f10455) reviews the October path changes since 2d1fc214.
- **build527.sh** (pid 90765) + committer: patch_builder517 → suite → 517's tests → KA-1 (446bj2z) → KA-496 (496a) → arm
  **446bj2zmz** (Arc-sourced hook, 527/PREREGISTRATION.md). Candidate for 10-07 with the production all-switches arm and K20;
  the 10-07 rule is fixed after the 09-30 read and before any 10-06 entry.
- 10:54Z: 446bd2 = 0.2704629 (+0.0004 vs 446b): FLAT, DLD-1 stays out. chain0930c armed 09-30: s1 458d, s2 458g.
- **10:56Z: the 10-01 rule AMENDED before any 09-30 entry exists** (`464/rule1001b.py`): when 458d reads null (fid − 446b
  > −0.0010, i.e. not HARM) or is missing, 446bj2zmz (527) is appended as the last candidate. In that branch the list held
  only 460dr, so slot 2 was empty; 517 §6 makes the Arc-sourced hook arm 1 there. Tested on fake boards (null/none →
  `460dr 446bj2zmz`; harm → `460dr`; positive unchanged). chain1001b killed and re-armed (pid 93212) with the tag's
  package and queue name (#22); logic otherwise unchanged.
- 11:14Z: **526 (cloud review, round 2) applied.** Four defects, all reproduced by its repro scripts and all fixed:
  (1) `gate524.py` passed artifacts that did not carry the recipe (missing same-line record, wrong tau / le, unnamed
  programme bank): now checks the recorded values; (2) a control file without `ntc_id` passed preflight and silently
  switched the builder to the random C_p proxy: `def_preflight.py` now fails closed (`DEF_DAY_ALLOW_NO_NTC_ID=1` overrides);
  (3) a reboot mid-pin left the builder pinned with no record, and a fresh `def_day` then treated the pin as the original:
  `def_pin_panel_artifacts.py` writes record + before-image before replacing, `def_day.sh` refuses a non-HEAD builder;
  (4) `recipe_dr524.sh` missed python builders in its wait and exited 0 on a failed gate after prep had overwritten the
  package: fixed and re-armed (pid 99304). (1)/(4) committed (150b334a); (2)/(3) are production files, committed by
  `commit526.spec` after `fix526.sh` (pid 99565) runs the inert KA (efp relaxed == ac49d319) once build527 has landed.
  Note for 529: gate524's MDE-sign check parses `<source>:<tau>` only; the hookbank form needs its own row.
- 13:38Z: 526's production fixes committed (600adea8) after KA526 PASS (efp relaxed through the patched def_day.sh on builder
  39c1dc66 == ac49d319: also proves 523's and 527's builder knobs inert on the October path). 527 committed the builder
  (67df0421); 446bj2zmz READY (ed90e3b3).
- **Roster at 13:50Z.** Coordinator (this session). Local headless Opus 5.5 workers: 528 (call-set composition),
  533 (science lead: audited manuscript update, works on a copy in `533-science-lead/paper/`). Cloud Grok 4.7 xhigh: 530
  (October bundle patches), 531 (final-phase tooling + tabletop), 532 (literature/data refresh). **Peer interactive session
  "Public data scope" owns 534 (generalization audit, at Ahmet's request): AGENTS.md §11.1c + CONFIRMATION-LEDGER + its own
  dir; Part B (sealed mock final, whitener refits ~13 GiB) only in free windows (not 09-30 00:30–06:00Z, not 10-05).**
- 14:20Z: **534 (peer session "Public data scope")**: Part A integrated (held_out lines in october-decisions.json, b8681568;
  its §4 corrected on lodo_hi_pct after my note: the 10-03 read is cap 6 vs cap 3 and a flat read counts against the switch).
  Part B (sealed mock final R6, whitener fits) runs in today's free window. **Part C** (preregistered ad88f51f) re-prices the
  three board-chosen base steps (normrel −1.5, its cap form, CD4 partial rows 0.5) on real/efp/efr before **10-03 00:00Z**.
  A step it marks FAILS gets a default-off def_day knob in the post-10-05 pass-through bundle (529/530), with an inert KA and
  a rehearsal, so it can come out of October.
- 16:30Z: **528 done (positive on H1, see FINDINGS §0); 446bj2zcp preregistered (535) and build535.sh armed (pid 73717,
  starts ≥ 09-30 06:00Z: emitter 0f56209d + builder 2973c48b, 528's tests, KA unset / inert artifact / 446bj2zmz, then the
  arm) with its committer.** 533's Draft v15 republished (artifact version 32). The 10-07 rule (fixed after the 09-30 read,
  before any 10-06 entry) chooses among 446bj2zmz (if not on 10-01), 446bj2zcp, K20 and the all-switches production arm.

### 2026-09-30 06:38Z — the 09-30 reads and three rule amendments (before any 10-01 entry exists)
- 458d 0.2615288 (HARM: Arc's H1 programme as the DE block prior, fid −0.0055 raw), 458g 0.2650949 (gate +0.0036 vs 458d).
  Decided: `programme_source off`, `prior_gate off`. Record: `458/BOARD-RESULT-458d-458g.md`.
- 10-01 (chain1001b, unchanged logic): s1 = 460dr (pid 94238, the 536-hardened waiter); the rule gave no s2 (HARM → no
  446bj2zmz, per 517's HARM branch it waits for 496a's 10-05 read).
- **Amendment 1:** 10-01 s2 = **446bj2zcp** (535, preregistered 09-29 16:45Z; bins unchanged) via `535/arm1001cp.sh` (pid in its
  work dir): after build535 exits, only if its arm passed (PASS, CLIFF OK, novelty, package sha), it arms 464's arm1001.sh as
  slot 2 with the queue name #23. If build535 is not done by 10-01 00:00Z, no s2.
- **Amendments 2–3:** `498/rule1004.py` and `496/rule1005.py` fallbacks drop 461dr / 469dr / 458e (they carry the programme
  that read HARM); 460dr stays.
- Slot 1's first 09-30 upload died when the Mac slept on battery (02:39Z); 536 hardened arm1001–arm1006 (retry/resume).

### 2026-09-30 10:55Z — 446bj2zcp rebuilt for 10-01 slot 2 (no rule change)
- build535 passed (tests, KA535a/b/c byte for byte, arm gates PASS / audit 0 / CLIFF OK), but `run_arm_generic.sh`'s prep
  gate never opened: the data volume was 99 % full (17.7 GiB free, gate needs ≥ 20). build535's unconditional last line then
  deleted the only copy of the prediction, so `arm1001cp.sh` armed nothing (10:38Z).
- Freed 52 GiB (AGENTS §10: `git clean -xdf` on 417's 09-23 build/dump trees, 428a/g/h, 399's failed build). `rebuild535.sh`
  (pid in `535/work/rebuild535.pid`) rebuilds the arm with build535's env on the committed builder 4daae237 / emitter
  0f56209d and must reproduce the first artifact `0b3733d5` byte for byte, else it writes `464/work/NOGO-1001-s2`.
  `arm1001cp.sh` re-armed on the rebuild's pid (pid 66296): same gates, same slot, same name #23. Record: `535/RESULTS.md`.
- Commit `2e448683` (the 535 production change) was amended before its push to name the production builder 4daae237
  (the message had test_patch528's test copy, 2973c48b): `06120021`.
- 11:12Z: **roster.** Local headless Opus 5.5 worker **537** (`experiments/537-arm-pipeline/`, mandate there): inventory,
  price and specify the candidate arms for 10-07 → 10-21 (decision tree on the 10-01 … 10-06 reads, env diffs vs 446bj2z,
  bins); first ranked list due 10-02 12:00Z. No other worker runs; 528/533/534 and cloud 530–532 are merged.
- 11:33Z: **10-01 slot 2 armed: 446bj2zcp** (pid 74867; package `fa556a05…`). The rebuild reproduced the first artifact
  `0b3733d5` byte for byte; gates PASS / audit 0 / CLIFF OK / novelty 0. Reader `464/read1001.py` (lodo + cp), written 11:11Z.
- 11:40Z: **coordinator queue to 10-07** (AGENTS §11.1b; heavy local jobs serialise through the memory gate):
  1. **10-01 reads** (`464/read1001.py lodo …` and `cp …`), BOARD-RESULT files, October switch updates. 460dr decides
     `DEF_DAY_LODO` confirmation; 446bj2zcp decides the call-set block's October rider (and opens 446bj2zmzcp).
  2. **10-05 (crowded):** 02:00–~04:00Z `recipe_dr524.sh` (automatic) → apply the October bundle (530 + 528's
     `patch_defday528.py` + cloud 538's extras, via 538's applier: `--check`, apply, tests, inert KA efp relaxed ==
     `ac49d319`, committer) → **R6 claim C10b** (534's protocol, the final bulk configuration, once, after the 10-05 read) →
     one rider-on rehearsal. 537's H1 jobs yield through the memory gate.
  3. **The 10-07 rule** before any 10-06 entry (10-06 00:01Z), from 537's ranked list (due 10-02 12:00Z) and the 10-01 …
     10-05 reads. 547 (and 537's final version) decide whether DE bins can use σ ≈ 0.0003 instead of 0.0045.
  Workers: 537 (local Opus 5.5, resumed 11:37Z with amendment 1: memory gate 24 GiB + no build); cloud 538 (Grok 4.7 xhigh,
  October bundle extras, branch `cloud/538-bundle-extras`). 539 (rival trajectories) and 547 (DE-noise check) are recorded.
- 11:44Z: **cloud 538 delivered and is merged** (patch scripts only, NOT applied; local 13/13 tests OK, `--check` clean, live
  files unchanged). On 10-05, after `recipe_dr524.sh`'s process exits, the bundle goes in with ONE applier that runs 530, then
  528's call-block (re-anchored: 528's own patch misses its anchor after 530), then the extras (`DEF_DAY_LODO_HI`,
  `DEF_DAY_CENSUS_MAX_ZERO`, `def_postbuild --require-config-diff`, REPO from the script path). Commands:
  `experiments/538-cloud-bundle-extras/RESULTS.md` "Commands the coordinator runs on 10-05". New switches `lodo_hi` and
  `callset_block` land PENDING. Do not also run `patch_bundle530.py`.
- 11:46Z: **cloud 548** (Grok 4.7 xhigh, branch `cloud/548-de-noise-audit`): adversarial audit of AGENTS rule 13's DE noise (0.0045) against 547 and 537's claim; verdict + proposed wording for rule 13 and FINAL-PHASE-FEEDBACK-PLAN §1.4. Nothing changes in the rules until 537's final version and 548 agree.
- 15:13Z: **to do before any 10-06 entry (deadline 10-06 00:01Z): amend 523's K-pair dose check.** 523 says K 2,000 must sit
  between 446bj2z and K 1,600 or the pair is "noise-dominated: no adoption". On the real production emitter (537 v3 §3.2 P2/P3,
  H1, 2–3 seeds) K 2,000's Δ(fid + jac) is ≥ K 1,600's in every regime (+0.0011 vs +0.0002 at f 0.5; tied ~+0.003 at f 0.35
  and cmp), and the DE read noise for this law is ~0.0003 (547/548): a non-monotone pair is H1's prediction, not noise.
  Draft: adopt the K with the larger Δ(fid + jac) if ≥ +0.0012; K 2,000 wins ties within ±0.0005; opposite signs beyond
  ±0.0012 → diagnose. **Done 09-30 21:10Z** (523/PREREGISTRATION.md, "Amendment").
- 15:13Z: **549** (`build549.sh`, 446bj2zcph = block + CLEAN4 hook) is written and committed, not launched: launch by hand only
  if 446bj2zcp reads ≥ PARTIAL on 10-01 (549/PLAN.md). It pauses 537's driver between jobs and resumes it on exit.
- 21:13Z: **peer session "Explore contribution areas"** (started by Ahmet) joins under the parallel-session protocol: numbers
  550–559 (the coordinator continues at 560+), commits only its own dirs, cloud Grok agents for code/public data, no heavy local
  jobs during builds or on 10-05 00:00–18:00Z, disk footprint ≤ 15 GiB. Territory: (b) a new sealed reserve **R7** (Arc-like
  panels disjoint from every tuned panel and R6, plus a DE-capable surface if one is untouched), declared before any draw,
  not accessed; (c) a cloud adversarial review of the 10-05 procedure (538 applier, inert KA, 549's path) by 10-04 12:00Z;
  (d) **R6's claim C10b** (534's protocol), moved by a ledger amendment to run after the 10-04 read, when the bulk
  configuration is final (10-05's read is DE-only); the coordinator sends the exact bulk env after the 10-04 read. (a), a fid
  decomposition, was declined: 361, 528, 525 and 480 cover it.
- 21:44Z: **machine windows agreed with the peer.** 10-04 05:00–12:00Z (after the coordinator's config message ~04:30Z; must
  end by 16:00Z): the peer's R6 C10b (experiments/551-r6-c10b; ~2.5 h of operator re-fits that count as builders). 10-04
  12:00–24:00Z and 10-05 00:00–18:00Z: coordinator builds only. The peer's R7: R7a (4 fresh Arc-like panels, declared
  b4dc26f9, drawn 85afb153) and R7b (Sunshine 2023 Calu-3 CRISPRi, a DE surface in a SARS-CoV-2-infected state: restricted
  to mock/uninfected cells if the data allow, fixed before building; asymmetric rule: it can veto a DE switch only on a
  preregistered HARM beyond its own noise band). R7 access only after the DE switches settle (not before the 10-07 reads).
- 21:58Z: **550 (peer's cloud review of the 10-05 procedure, 09bae587) accepted.** On 10-05: follow
  `experiments/550-cloud-1005-review/RESULTS.md` "10-05 command sequence" exactly; the inert KA by hand in a fresh dir with
  `env -u` on every DEF_DAY_* knob (NEVER re-run chain502/fix526/chain476: chain502 copies pre-502 files over live scripts
  before its wait); only 538's applier (530 + 528 re-anchored + extras); only after recipe_dr524.sh's process exits. Every
  DE-only arm gets `550/repro/bulk_identity550.py` (per-(context, target) summed counts, h5ad vs h5ad) and
  `config_keys550.py` before it is slotted: run_arm_generic's CLIFF band (0.85–1.15) lets a 10 % bulk move pass. The
  hook-conditional resolver patch is decided at the 10-05 `decide` step.
- 22:05Z: **446bj2zcph READY** (549): package `de1e4ab9…`, gates PASS / audit 0 / CLIFF OK / novelty 0; config diff vs 446bj2zcp
  exactly `de_mde_sign`; bulk vs 446bj2zcp identical in 895/900 groups, 5 groups swap one count between two low-count genes
  (the hook's whole-count moves; expected pds/mse ≪ 1e-6). Packaged by hand: `run_arm_generic.sh`'s prep gate
  (free+inactive ≥ 23 GiB) cannot open while ~18.5 GiB of other apps' pages sit in the compressor (memory_pressure ~60 % free,
  no swap). **To do (daylight, with a test):** give that gate a pressure-based alternative (memory_pressure ≥ 55 % and swap
  free ≥ 400 MiB) so the next builds do not stall. The peer's gates approved: ≤ 8 GB jobs at ≥ 20 GiB + pressure ≥ 55 %; C10b
  fits at ≥ 22 GiB + pressure ≥ 55 % + swap free ≥ 400 MiB.
- 10-01 01:04Z: **reads done.** 460dr 0.2692306 (+0.0028694 vs 453dr; LODO key transfers on the October path; `DEF_DAY_LODO=1`
  stands). 446bj2zcp 0.2720318 (−0.0007942 vs 446bj2z; NULL on Δ(fid + jac); jac +0.0040, fid −0.0036 raw): the block alone
  is not an October rider; 537's rescue condition is met, so **446bj2zcph (READY, `de1e4ab9…`) is the 10-07 s1 candidate**
  (its board price to be re-derived from this read before the 10-07 rule). Displayed row now 0.2720318 (rank 9).
- 10-01 01:04Z: **disk plan for 10-05** (48 GiB free after sweeping callprop535.h5ad): sweep 497dr's source h5ad after its read
  (10-03/04, +11 GiB) and `git clean -xdf` the 428 build trees from 09-23…09-28 as they pass 7 days (~20–28 GiB) on 10-04,
  for ~79 GiB before 10-05's three builds. 383-state (~109 GiB) still awaits Ahmet's OK.
- 10-01 08:38Z: **October robustness (peer worker 553):** the K-class emitter's steering dose (CLEAN4 block prior, MDE, hook)
  collapses SILENTLY to λ 0 when a context's control-library spread q75/q25 exceeds ~3.5–5 (A/B/C 1.57–1.95; placement
  crashes by 7.24), and the MDE/hook dose erodes with low sub-gate mass per cell. Asked the peer for an unapplied guard patch
  (def_preflight: q75/q25 ≤ 3.0 or DEF_DAY_ALLOW_WIDE_LIBS=1; def_postbuild: per-context median prior_lambda / mde_lambda ≥ 0.9,
  call_set_missed ≤ 100 or DEF_DAY_ALLOW_LOW_DOSE=1), calibrated on every meta we have, by 10-08 12:00Z; applied after the
  10-05 bundle with an inert KA. R7b (Calu-3, ≥ 5,000-count cells, 106 targets) is ready for access after the 10-07 reads.
- 10-01 08:59Z: **554 dose guard accepted** (peer, 138572b4; UNAPPLIED). Coordinator re-ran it locally: test_guard554 5/5 OK;
  identity554 (real data): preflight text + rc identical on 6/6 bundles (spreads A/B/C and every relabelled clone 1.57–1.95);
  postbuild identical on 11/12 metas, the 12th (474's efr rehearsal, programme prior pooled λ 0.278) gains exactly the guard's
  FAIL row; R7b-as-built (7.24) and 429a (A mde 0, C 0.856) trip as designed. **10-05 sequence:** 538's applier → 554's
  `patch_guard554.py` (--check, apply) → 538's 13 tests + 554's 5 + identity554 → `bash -n` → the inert KA by hand (550 step 4)
  with the audit dump now always on: must still give `ac49d319`, and postbuild.json must carry dose_guard per context →
  committer. Add 524dr's meta to calib554/identity554 once it exists.
- 10-01 12:59Z: **second peer session "Explore contribution areas [57577c]"** (started by Ahmet; numbers 590–599; own dirs only, no
  builds / vcc / cell-eval2 / shared-doc edits, ≤ 4 GB jobs outside the reserved windows): **590** validates
  `scripts/def_identity.py` (the day-of CCLE call that switches same-line rows on) on MIX-seq's 209 labelled lines and other
  panels, including the false-positive rate of a unanimous registry match and a Flex-like read; **591** the
  `run_arm_generic.sh` pressure-based gate as an UNAPPLIED patch + stubbed tests; optional later: manuscript v16 on a copy.
  537 v5 (01:20Z) fires the rescue branch: 10-07 s1 = 446bj2zcph; s2 = K20 if k16 reads ≥ MID on 10-06; it predicts the block
  does not reach its October shipping bar even with the hook (CPH vs 446bj2z Δ(fid + jac) −0.0001 … +0.0009 < +0.0012).
- 10-01 13:14Z: **disk: 99.6 GiB free.** Ahmet OK'd the 383 cleanup conditional on an independent review; the review
  (read-only Opus agent) agreed and corrected the gain: the two support-set folders were APFS clones of each other, so deleting
  them plus `runs/` freed 55.4 GiB (not ~109). Zip kept (intact, hash-pinned). The 10-04 sweeps (497dr's h5ad, aged build
  trees) are still planned but no longer critical. GitHub caught up: `0e4a6e82` (2,283 small outputs since 09-28).
- 10-01 13:25Z: **590 accepted** (second peer, c3f4c5d5; UNAPPLIED). `def_identity.py` validated on 402 labelled public profiles:
  top-1 0.985 when the line is in CCLE (0.995 when unanimous; every miss a Cellosaurus duplicate or same-patient pair); 0
  false switch-ons on 146 unregistered sister lines; THE failure is a line absent from CCLE (HEK293 → HAP1_ENGINEERED, 8/8
  unanimous, a registry line). Gate: registered AND unanimous AND min_top1_score ≥ 0.65 removes 2/2 false switch-ons (3' and
  Flex-simulated reads), costs 2/19 correct ones (none D/E/F-shaped); A/B/C clear it by ≥ 0.069 (0.772 / 0.719 / 0.733).
  Coordinator re-run: test_identity590 ALL PASS; `--check` clean on the live and post-530 def_sameline; live files unchanged.
  **10-05 sequence gains a step:** 538's applier → 554 → `patch_identity590.py --apply` → all tests → inert KA (the gate is
  inert there: only A reaches it, 0.772). Day-of: log each context's min_top1_score; eyeball the consensus top-5 if 0.60–0.70.
- 10-01 13:32Z: **591 applied** (run_arm_generic.sh c58b9553 → e023aebc): prep gate has a pressure branch (≥ 55 % + swap free ≥ 400 MiB); build gate unchanged (def_pipeline refuses rather than waits). test_gate591 ALL PASS.

### 10-05 CHECKLIST (consolidated 2026-10-01 14:03Z) — SUPERSEDED by "10-05 CHECKLIST v2" at the bottom (650 tabletop, 10-02 16:40Z)
0. `recipe_dr524.sh` builds 524dr from 02:00Z: wait until its PROCESS has exited (`experiments/550-cloud-1005-review/repro/procs_clear550.sh`
   prints CLEAR); no run_arm_generic.sh running. Confirm live == HEAD on the 11 paths (`git diff --quiet -- <paths>`; NOT
   594's frozen table: `october-decisions.json` moves with every `decide`). HEAD must be the pre-bundle commit until step 9.
1. BEFORE applying: `uv run python experiments/590-identity-validation/test_identity590.py` (590's suite is pre-application only).
2. 538: `apply_october538.py --root . --check`, then `--root .`; `test_extras538.py` (13/13). Never `patch_bundle530.py` separately.
3. 554: `patch_guard554.py --check`, then apply; `test_guard554.py` (5/5), `identity554.py`. (554 refuses before 538.)
4. 590: `patch_identity590.py --check`, then `--apply`; `scripts/def_identity.py --verify` (prints 0.772 / 0.719 / 0.733).
5. `bash -n scripts/def_day.sh`; `uv run python scripts/october_recipe.py check --allow-pending`.
5b. (not a gate) `scripts/test_october_recipe.py`: after 538, five tests fail because 538's new PENDING switches (lodo_hi,
   callset_block) give rc 2 where they expect 0/3 (598's flag): add the new switches to their decide lists, keep 6/6 green.
6. `uv run python experiments/594-def-bundle-composition/post594.py` → ALL PASS (b425c8bc: it rebuilds the expected 11 files
   from `git show HEAD:<file>` + 538 → 554 → 590 and compares live to that; 594's table for reference: def_day.sh 5851f84b,
   def_day_helpers e3ff45bc, relevel_iter e360a1d9, def_sameline 22c6f0f0, october_recipe 503c45a0, october-decisions 8ee9e74d,
   gate524 7a43d2f5, def_postbuild 557c2b52, def_pipeline 9a718b9e, def_preflight 2d181194, def_identity c98a608f).
6b. 641 (depth row, [537757], c12cce73; applies AFTER post594 because it changes def_preflight.py + def_day.sh):
   `uv run python experiments/641-depth-row/patch_depth641.py --root . --check`, then without --check;
   `uv run python experiments/641-depth-row/test_depth641.py` (4/4). A/B/C ratios ~1.00: silent in the KA.
7. Inert KA BY HAND (550 RESULTS "10-05 command sequence" step 4): efp relaxed, fresh work dir, `env -u` on every DEF_DAY_*
   knob, audit dump on → output must be `ac49d319…`; postbuild.json carries `dose_guard` per context. NEVER re-run
   chain502 / fix526 / chain476.
8. `run_arm_generic.sh`: add the post-538 driver `9a718b9eaa24f10cfbf255ce830ea4a748f61252ec503bbb6db69d3ddf5ca9fd` to its
   driver case list (line ~58; only while no run_arm_generic.sh runs); `bash -n`; `bash experiments/591-prep-gate/test_gate591.sh`.
9. Commit the 11 files + run_arm_generic.sh with the KA hash in the message (committer spec or by hand, explicit paths).
10. `october_recipe.py decide` for the switches the 10-01 … 10-05 reads settled (sharpen_pt 10-02, lodo_hi_pct 10-03, mde_sign
    10-05; lodo_hi and callset_block stay PENDING unless read), then ONE rider-on rehearsal on efr with the resolved command.
11. 20:00Z: chain1006 arms 10-06 (524dr + K16) automatically. Before 10-06 00:01Z: write the 10-07 rule (537's tree, v5+).
- 10-01 14:05Z: **second peer's queue:** 593 (cloud Grok, adversarial audit of manuscript v16's new claims, 503-style); 595 (foreign-control day-of rehearsal on scratch copies with 538+554+590 applied: mock D/E/F from Orion HCT116 NTC, Arc H1 2025 controls, Nadig Jurkat NTC on the current panel; preflight spread guard, identity incl. min_top1 + top-5, sameline, census; no builds, ≤ 4 GB); 596 (cloud inventory of public 10x Flex datasets of named cancer lines, to re-check 590's 0.65 floor on true Flex profiles).
- 10-01 14:19Z: **595 (second peer, c692306a): the day-of control stages behave on FOREIGN controls** (mock D/E/F bundles in exact
  format from six other lines, 538+554+590 applied on scratch copies; no builds). Registry lines switch same-line ON, unanimous
  above the floor (HCT116 0.794, A549 0.773, Calu-3 0.765, THP-1 0.665); non-CCLE lines stay off (H1 split call 0.653;
  **HEK293T → HAP1 unanimous 0.616 blocked only by 590's floor**: the hazard confirmed in situ); library spread 1.57–2.10 on five
  lines / four chemistries (H1 Flex 1.575), Calu-3 7.53 fails preflight and the override works; census clean. **Oct 22 rule
  (add to DEF-RUNBOOK when 590 is applied on 10-05):** log every context's min_top1_score and identity top-5; for a value in
  0.60–0.70, read the top-5: a correct registry call can sit just above the floor (THP-1 0.665), and a correct one that falls
  just below it is blocked safely (no same-line rows): force it only if the top-5 is unanimous and lineage-consistent, and
  record why. Peer's memory deviation (bundle build 15.9 GB RSS, outside the ≤ 4 GB budget; no harm) noted in 595 RESULTS.
- 10-01 14:26Z: **596 (cloud Flex inventory, 0ec58ac6):** 8 public Flex datasets with named human lines (Jurkat, K562, A375, MCF7, PC-9, REH, K562-CRISPRi NTC, VIPerturb); none names a likely D/E/F line (no HeLa / CAL33). 590 addendum: donor viperturb_k562 is real GEM-X Flex (true min_top1 0.765, unanimous K562). **597 planned** (pending Ahmet's OK for one 1.58 GB download, zenodo 18460279, CC BY 4.0): the VIPerturb pilot's HEK293 + HAP1 + K562 Flex NTC cells through 590's matcher, preregistered: HEK293 < 0.65 confirms the floor on true Flex.
- 10-01 14:38Z: **597 (second peer, 34a00b6b): 590's floor CONFIRMED on real Flex** (VIPerturb pilot, Flex v1 NTC cells; download
  OK'd by Ahmet; known answers PASS): K562 0.693 and HAP1 0.747 unanimous and correct; HEK293 (absent from CCLE) 0.6472 and SPLIT
  (SNU-398 5/8, HAP1 2nd): blocked by the floor AND by unanimity. Real Flex raises an absent line's score by ~0.1 over 3′ and the
  590 simulated bias (unsafe direction); on Flex, correct calls 0.693–0.772, absent lines 0.647–0.665, all split.
  **Decision: keep the preregistered 0.65** (no post-hoc retune on two points; unanimity blocks both Flex absent-line cases;
  the Oct 22 reading rule covers 0.60–0.70). Oct 22 note: on Flex, a unanimous call at 0.65–0.70 gets its top-5 read before
  the day proceeds.
- 10-02 00:59Z: **reads done** (479p2g 0.2733131 NEW BEST, P2g TRANSFERS; 479sp2 FAILS); `sharpen_pt` DECIDED tercile:2.5,2.0,1.75 (00:58Z; the 10-05 dress rehearsal now resolves with it). **Queued: 598** (second peer, unapplied patch + tests by 10-07 12:00Z): builder mode `tercilec` (493h's carried-only cut points, non-carried → MID), def_day validation on the post-594 file, KA plan (unset == ec3a5ae6; tercile == 479p2g's artifact; tercilec moves 4 A/B/C targets); coordinator applies ~10-08 with the KAs, then `decide sharpen_pt tercilec:2.5,2.0,1.75`. 10-03 slots arm at 20:00Z (495 / 471c). The best entry has no zs rider: P2g + zs is a candidate for the all-switches arm.
- 10-02 01:10Z: **553 done (first peer): R7b's sealed read is ready** (one command, ~25 min, < 4 GB, truth unopened); access after the 10-07 reads with two preregistration rules (noise band includes control re-draws, band_total; retrodict against the SHAPE-MATCHED H1 rehearsal). **Pricing insight:** H1 cut to R7b's cell/control counts reproduces the board's ratios for the call-propensity block far better than H1-full (jac 0.45× / fid 1.4× / reach −0.05× vs board 0.48× / 1.95× / −0.38×): consider shape-matched H1 for pricing CPH and the K arms before the 10-07 rule (asked the peer for the recipe). Expected at R7b's shape: CP +0.0014, CPH +0.0015, K 2,000 −0.0016, hook ~0.
- 10-02 01:16Z: **shape-matched DE pricing (553, 777f907a):** recipe in 553 RESULTS (last section); prices on disk (draws553.json): CP +0.0014, CPH +0.0015, K2000 −0.0016, CPHK20 −0.0001, K-on-CPH −0.0017, hook ~0, K1600 −0.0030 (1 draw). It is a lower-power regime instrument: it matches the block's jac but under-reads fid's cost and misses reach's sign (CP read −0.0008 on the board): use it only as far as it retrodicts the 10-01 / 10-05 / 10-06 reads. 553's 31 GiB work/ pruned by the peer (keep h1r7 surfaces, R7b access inputs, KA references).
- 10-02 01:31Z: **598 accepted** (second peer, 9c608e0f; UNAPPLIED; coordinator re-run 43/43 PASS, live untouched): builder mode
  `tercilec` (carried-only cut points, non-carried → MID; `carried` from the DUMP's key artifact, sha + target order checked, and a
  loud refusal if ||bulk_t||/||carriers_t|| ≠ c·w_t to 1e-6). On A/B/C it moves exactly 4 targets (E2F4, SLC25A13 HIGH→MID;
  EPHB2, TMEM209 LOW→MID); on day-of dumps 0–4 move (value on thin-coverage panels). Touches 6 files (builder, def_day.sh
  anchored on post-bundle 5851f84b, helpers meta_matches, def_pipeline september_rows, october_recipe, gate524).
  **Apply on 10-06 (after the 10-05 bundle):** `patch598.py` → `post598.py` → KA598a (unset == ec3a5ae6), KA598b (tercile ==
  479p2g's dd82740d: the first rebuild of 479p2g on builder 4daae237), KA598c (tercilec: `ka598.py <meta>` PASS) via
  run_ka_527.sh → commit → `october_recipe.py decide sharpen_pt tercilec:2.5,2.0,1.75`. run_arm_generic's driver list needs
  def_pipeline's post-598 sha too: take the sha post598.py prints at apply time ([1] live == canonical), not 598's
  table value (76314af6… was computed from HEAD 54b8886d + 538 and moves if def_pipeline changes first).
- 10-02 01:41Z: **R7 machinery ready (first peer, 21381a77 / 9725b988; no sealed data touched):** R7a driver `551/r7a551.py` (plan / KA-CD4 pass; fresh fits on T1–T4; ~3–4 h), R7b one-command read (~25 min), claim draft `552/CLAIM-R7-DRAFT.md` with rules fixed now (R7a PRIMARY = F vs M5, 12 folds, paired bootstrap seed 552; R7b board-likeness first, HARM only at ≤ −2.5 × band_total). C10b's P2g = 598's tercilec rules (same cut points, floor, ties). Cloud review 555 running. **Tentative window: 10-07 06:00–12:00Z**, confirmed after the 10-07 read with the final F and DE arms.
- 10-02 02:23Z: **555 (cloud review of R7) merged and fixed** (first peer, d9cbc6f0 / 90ec5d3d; RESPONSE.md): 3 blockers closed (R7a code pin now covers the full import closure; sealed access arms only if the claim sha is in the ledger's R7 section; un-emitted switch pairs are listed in plan.json with band_total None, never dropped). Draft claim rules: K switches REPORT-ONLY on R7b; block / hook / CPH veto at ≤ −2.5 × band_total (false-veto 0.6 % at the band, 5.5 % for CP at the rehearsal spread); board-likeness = raw jac > 0 and raw fid < +0.001. **Coordinator's 10-07 ~04:30Z message must give F AND the list of DECIDED DE switches** (only block / hook switches can veto).
- 10-02 02:27Z: **10-07 ~04:30Z message to the first peer (it pastes it verbatim into CLAIM-R7 §0; fill_claim552.py, e465c987):**
  exactly (1) F = the full `october_recipe.py command` line after the 10-07 decides, (2) block on/off, (3) hook on/off,
  (4) K default / 2000 / 1600. R7b pairs: block + hook judged at the default size (B+H final: veto pairs CPH vs hook and CPH vs
  CP); K 2,000 on top is report-only (CPHK20 vs CPH). If only K or nothing is decided, **R7b is NOT accessed** (stays sealed).
  The peer has a wake task at 10-07 04:20Z.
- 10-02 10:34Z: **peer addresses changed after Ahmet's Claude Code update:** first peer (550–559; C10b 10-04, R7 10-07) = "Explore contribution areas [10eca4]" (socket /tmp/cc-socks/28465.sock); second peer (590–599) = [36a81a]. Addresses can change at every restart: before sending the 10-04 C10b config (~04:30Z) and the 10-07 R7 message (~04:30Z), run ListAgents and send to the row that is the 550–559 peer (reply to its latest from= address).
- 10-02 10:48Z: **council session ([4d0207] ops; [21dd4c] board forensics; [7772dc] biology / new model families; [d4cd5f] metric headroom / DE wall; numbers 600–639, own dirs only)** found recipe_dr524's free+inactive ≥ 30 GiB floor could keep 524dr from building (29 GiB now). Fixed: floor 25 GiB (above def_pipeline's 23 GiB build-stage refusal), re-armed pid 33446. AWS: 5 stopped VCC instances bill ~$96/month EBS (council asks Ahmet before terminating).
- 10-02 10:54Z: **session [537757] "Virtual Cell Challenge 2026"** (Ahmet; Arc's 10-01 "Behind the Data 2026" post): panel = targets with ≥ 80 % median knockdown (median 89 %), 400 of a median 500 cells; guides at week −2, collection week 0 (~14 d; K562 GWPS day 8); low MOI ~10 %; same sgRNAs in all six lines; full depth 29,482 UMI/cell, A/B/C controls downsampled to 20.0k (expect D/E/F ~20k); 100 NTCs / ~50k controls per line (we get 46 × 400); an RPE1 / hTERT figure is a weak hint at a non-cancer D/E/F line. 610 (to be renumbered 640: 600–639 belong to the council): October-like draws carry ~2× more weak-knockdown targets than Arc's panel; 438's CD4 ×1.5 loss survives knockdown-matched subsets → October CD4 ×1 stands. Asked for: a day-of depth row (unapplied, after 554) and an RPE1 check (identity + same-line coverage).
- 10-02 10:57Z: **640 accepted** ([537757], 88d9f8d3): October CD4 ×1 survives Arc's selection rules (knockdown-matched and strength-reweighted subsets all keep 438's ×1.5 loss; abc's CD4 gain sits in the weak energy tercile: a target-identity effect). RPE1: not in our CCLE reference (590: abstains at 0.466); Replogle RPE1 covers 0/300 abc, 8–9/300 October-like: no action. **Contingency queued (low priority):** an UNAPPLIED HEK293 identity reference + Orion-HEK293T registry row (294/300 coverage), OFF by default, switched on only by a unanimous above-floor HEK293 call AND the coordinator's explicit decision.
- 10-02 11:17Z: **council portfolio v2** (ac737573) recorded in FINDINGS §0: the gap is the DE channel; 635's correction (responding-gene programme is the prize, +0.008 … +0.010); 634/601/610 closed. 611 runs local real-emitter jobs (~7.5 GB each, 2 in parallel under the gate) today, outside my windows; DE slot proposals from ~10-08. HEK293: a Flex HEK293 profile matches Orion's HEK293T at only 0.17–0.20 (council 620): [537757]'s 642 contingency is likely moot.
- 10-02 11:20Z: **637 (council [d4cd5f], 3fc82ed9) accepted in principle:** a factorial final phase (each D/E/F entry after R0 carries one bulk switch, read from pds + mse exactly, plus one prior-only DE switch, read from the DE members minus a −0.0003 ± 0.0002 bulk-leak allowance; separability measured on 8 board pairs) doubles the switches read Oct 23 → Nov 3 (~11 → ~22) IF members come back. **Its A/B/C validation slot = the all-switches arm** (446bj2z + P2g as bulk X + the adopted DE switch with a single 446bj2z read as Y), preregistered with 637's pass criteria in the 10-08 … 10-10 rule. readdef.py --factorial patch (unapplied) requested by 10-09. H-631 feedback-use policy (preregistered switch selection only; no gene / target probes, no per-target n_pred tuning, no precision encoding) is with Ahmet.
- 10-02 11:23Z: **642 (HEK293 contingency) goes ahead** ([537757]; correction to my earlier pointer: 620's 0.17–0.20 was its own Spearman route, which 620 itself calls unreliable): 590's matcher + one non-Flex HEK293_REF row (XuCao 3′ + Orion HEK293T) calls the held-out real-Flex HEK293 8/8 at 0.820 with 0/423 false calls; a Flex-built reference is UNSAFE (attracts Flex HAP1). Two-stage, default-OFF patch (DEF_DAY_ALLOW_NONCCLE_REF=hek293), anchored after the post-bundle files; open question: whether Orion HEK293T rows transfer to an Arc-Flex HEK293 context; switch-on is the coordinator's explicit decision.
- 10-02 11:25Z: **637 reader accepted** (council, c958b9c0; local 27/27 + 531's 74/74): `readdef_factorial637.py` allows bulk rows 2/3/4 × prior-only DE rows 6/7, refuses the rest, and requires every verdict to hold across the leak band (else "inconclusive_leak_band", nothing adopted); `--validate` implements the preregistered pass rule. `readdef637.diff` (5-line dispatch) stays UNAPPLIED until the final-phase tooling freeze. Council SYNTHESIS v4 pivot rule: if 621's supervised DE model misses its bar, a preregistered global DE dose sweep on the board from 10-08, conditions: priced first on 611's instrument (retrodicting 446bja / 458d/g / cp / zs), global protocol-level only, two-point doses with fixed bins.
- 10-02 11:28Z: **council DE candidate (610 forensics, [21dd4c]): a RANK-NEUTRAL mean shift moves nmae alone.** 201 DE-only board steps from 105 teams show many pure-nmae steps with fid/jac/reach flat (up to −0.070 raw nmae, +0.019 overall in one entry; ours zs −0.0022). nmae reads arithmetic means, calls and reach read Wilcoxon ranks: inflating cells already above every control value raises the mean without moving a rank (UP unbounded, DOWN bounded). If real: scorer mechanics, +0.005 … +0.019 overall; risk: CLEAN4 is ~84 % down. Being priced by [d4cd5f] on 537's real emitter under conditions (KAs, six members, H1 + h1r7, per-gene U and calls unchanged, bulk ≤ 1e-6, escape argument vs 426/435/426b/434/480, 2-point preregistration). Also queued: a confidence-ordered reach (chill / Chiu S1: reach +0.03) and a mean-preserving dispersion law (William: fid +0.010 / jac +0.004).
- 10-02 11:29Z: **611 delivered** (H1 f 0.35 retrodicts delivery pairs; not call-set reach) and **621 negative** (supervised sign model; public sign sources closed): recorded in FINDINGS §0. Council design sweep, ranked: #1 rank-neutral nmae ([d4cd5f], under the 4 conditions; 611 meets condition 1), #2 calls-unchanged confidence-ordered reach ([7772dc], keyed on 621's consensus z), #3 mean-preserving dispersion law ([21dd4c]). Slot proposals with preregistrations from ~10-08.
- 10-02 11:34Z: **639 identity-axis sign prior: POSITIVE held-out (3/4 lines)** → candidate October DE switch (transfer class). Needs: a delivery spec (per-context prior files A/B/C from each context's controls vs the others; builder support for per-context priors beyond 446bja's context-A knob; the hook's sign source), an emitter price, a preregistered board arm for the 10-08+ rules, and the October path (def_day computes I from D/E/F's own controls on the day).
- 10-02 11:35Z: **612 (council [21dd4c]): the KEY-GATED uncalled-gene hook** (CLEAN4 sign where our key is in the target's weak half, the key's sign elsewhere) is min-max-regret on 480's held-out cancer folds: ≥ 0 vs key-like production on all 4 surfaces (+0.016 … +0.018 nmae on HCT116 / HEK293T, CIs excl. 0); vs the plain CLEAN4 hook +0.012 … +0.020 on h1 / K562, −0.003 … +0.001 on Orion; worst regret 0.0032. **Decision: 10-07 stays as 537's tree** (CPH tests the block rescue; swapping the hook form would confound it, and the key-gated form needs an emitter mode with KAs). It becomes the leading October-hook candidate for a 10-08+ arm, consolidated with 639's identity score (one board arm), read against the plain hook's own board effect (496a, 10-05).
- 10-02 11:36Z: **639 delivery plan** ([4d0207], no builder change): a same-line-format artifact gives every context the identity-augmented prior as its programme default (VCC_SAMELINE=<npz>:<sha>:0, w 0, no rows), P quantile-mapped onto CLEAN4's own P distribution (q80 / shrink floor / τ unchanged); KA λ 0 == CLEAN4's P exactly. **Hard gate before any slot: per-context median prior_lambda ≥ 0.9 on the A/B/C build (554's audit dump)**, so 458d's 0.43 cannot repeat. Roles: 639 plumbing [4d0207]; price / prereg + rank-neutral #1 [d4cd5f]; key-gated hook mode + consolidation [21dd4c]; 639 stress tests [7772dc].
- 10-02 11:36Z: **612 implementation plan** ([21dd4c]): builder-only form `VCC_DE_MDE_SIGN=clean4kg:<tau>` (per-target mask: CLEAN4 sign where |design row| is below the target's median among nonzero entries, else keep the delivered key sign; emitter unchanged 0f56209d; unset and clean4:<tau> byte-identical). H1 analytic: KG hook −0.0049 nmae at f 1 vs the plain hook's −0.0015 → board ≈ +0.0011 … +0.0015 overall over the plain hook (611 multipliers). UNAPPLIED patch anchored on the POST-598 builder (composition test 598 → 612), KA-1 / KA-2 + H1 price + the 639 consolidation, due 10-06 12:00Z; applied 10-06 after 598.
- 10-02 11:42Z: **642 accepted** ([537757], eea35ae2; UNAPPLIED, default OFF): two-stage identity (stage-2 HEK293_REF from NON-Flex sources only; vetoes the CCLE call when it fires) + def_sameline orion_hek293t row (294/300 at w 2) behind DEF_DAY_ALLOW_NONCCLE_REF=hek293 AND a recorded DEF_DAY_NONCCLE_DECISION. 6/6 tests: 590 reproduced 423/423, 0/423 false HEK293 calls, held-out Flex HEK293 8/8 at 0.820, Flex-built references refused. **Change requested: stage-2 floor 0.75** (new, unused parameter; near-misses at 0.67–0.69 on shallow 3′ A549 / MCF7). Apply 10-06 after 598 / 612. Open question: Orion HEK293T rows' value for an Arc-Flex HEK293 context is unmeasured.
- 10-02 11:45Z: 642 updated (478d1f8d): stage-2 NONCCLE_FLOOR 0.75 (stage 1 keeps 0.65); patch sha256 f62b631d…, reference 1b829de0… unchanged; 6/6 tests (every near-miss 0.46–0.69 now below the floor). Apply 10-06: `patch_hek293_642.py --root . --check`, then without.
- 10-02 11:48Z: **rank-neutral nmae (council #1) bar fixed:** strict U-invariance has zero capacity on the real emitter (no uncalled-UP carrier above the control max; top carrier a median 116 controls below it), so one carrier per gene is promoted above the max (ΔAUC ~1.6e-5). Condition (3) = calls exactly unchanged (0 flips), |ΔAUC| ≤ 1e-4 per gene, pds / mse ≤ 1e-6 raw, all six members measured, reach read on the board with 611's band. 604 identity-prior artifacts built (λ 0 == CLEAN4 exactly; λ 1 flips 1,606–1,843 genes per context; block down-fraction 0.66–0.68 → 0.76–0.84, favouring full dose); being rebuilt on 622's Flex-only references.
- 10-02 11:56Z: **612 redesigned as an orthogonal knob** ([21dd4c]): `VCC_DE_MDE_SIGN_KEYGATE=<q>` masks the hook's per-target sign, whatever its source, to each target's weak-key genes (strong half → keep delivered); composes with 604's identity prior (file:/programme: + KEYGATE = the identity hook) and with clean4:0.05 (the CLEAN4 hook). test612 17/17 (incl. 598 → 612 composition; unset returns the same object for every source; H1 real-emitter KA bit for bit). KA plan: KA612a unset == ec3a5ae6; KA612b clean4:0.05 == e21ff752; KA612c clean4 + KEYGATE 0.5 with ka612.py (config diff vs 496a exactly de_mde_sign_keygate). **10-06 apply day:** 598 (KA598a/b/c) → 612 (KA612a/b/c) → 642 (default off) ≈ 6 KA builds (~4 h), then decides; the ONE-form hook recommendation + 10-08+ prereg arrive from the council first.
- 10-02 12:05Z: **council:** (1) 622g independent confirmation of the identity rider (Broad-16, other lab, Cas9 KO, 3′): z(C4) + 0.5·z(I) +0.009 mean over 12 lines, 0 significant losses; λ 1 has one significant loss: **ship point λ 0.5**, λ 1 the second dose; guard: drop any reference whose identity call matches the context (same-line reference kills the signal). (2) Identity-type sign sources are REPORT-ONLY on R7b (chemistry-mixed I; asked [10eca4] to add it to the claim). (3) **630 rank-neutral UP magnitude ([d4cd5f]) PASSES the bar on the real emitter:** b_k100 Δnmae −0.0066 raw at f 1, fid / jac / reach / pds / mse EXACTLY 0 on both seeds, 0 call flips, max ΔAUC 1.0e-4; f 0.35 +0.0028 [+0.0013, +0.0043] overall-equiv; board estimate +0.0026 … +0.005, nmae only. Builder path UNAPPLIED (src/vcc2026/rank_up.py + patch_builder630.py, VCC_DE_UPRN=<s|b>_k<K×100>, after zs); dose pair K 0.75 / 1.0. Applied 10-06 if its 598 → 612 → 630 composition test lands by 10-06 06:00Z, else 10-08.
- 10-02 12:10Z: **644 queued** ([537757]): DEF-RUNBOOK amendment draft (diff, coordinator merges) for 554 / 590 (+ the reading rule) / 641 / 642 / 637, placeholders for 612 / 630 / 639, each override's legitimate use and what to record, every command verified in --help / --check / dry-run form, plus a one-page Oct 22 decision card; due 10-08 (after the 10-05 bundle).
- 10-02 12:11Z: **612 CORRECTED by its author ([21dd4c]; boxed in 612 RESULTS):** the earlier "3× the plain hook / +0.016–0.018" came from a zero-base model; in rc_emit_k the hook only sets the direction of a small MDE add on a key-signed base (it flips 2–6 % of signs). Real emitter, H1 s0: key-gated vs plain hook −0.0001 (f 1) / −0.0003 (f 0.35) overall. Delivery-faithful folds: the plain hook pays where CLEAN4 fits the line (HCT116 / HEK293T +0.015) and hurts where it does not (K562 −0.010); the gate keeps most of the gain and removes the downside (worst-case regret 0.0020 vs 0.0115). **So the key gate is a ROBUSTNESS lever for unseen D/E/F lines (~0 EV on CLEAN4-friendly lines), not a large gain.** The October hook form (plain vs gated) is decided after 496a's 10-05 read and the council's seeds 1–2 / idkg runs. Lesson again: price on the real emitter before believing an analytic gain.

### 10-06 APPLY DAY (written 10-02 12:12Z) — SUPERSEDED by "10-06 APPLY DAY v2" at the bottom (650 tabletop, 10-02 16:40Z)
1. 598 `experiments/598-tercilec-knob/patch598.py` (--check, apply) → test598 → post598 (note its def_pipeline sha).
2. 612 `experiments/612-key-gated-hook/patch612.py` (--check, apply) → test612.
2b. 612dd `experiments/612-key-gated-hook/patch612dd.py` (--check, apply; def_day.sh pass-through DEF_DAY_MDE_SIGN_KEYGATE →
   VCC_DE_MDE_SIGN_KEYGATE, insert-only) → test612 (22/22 incl. T8–T11).
3. 630 `experiments/630-rank-neutral-up/patch_builder630.py --root .` (--check, apply; installs src/vcc2026/rank_up.py) → test630 (27/27
   verified locally 10-02) + test_rank_up630.
3b. 630dd `patch630dd` (sha 4db2a0c5, commit 626949ef; def_day.sh pass-through DEF_DAY_DE_UPRN → VCC_DE_UPRN; needs DEF_DAY_DE_MAG;
   unset-inert; the final meta must carry `de_uprn.variant` + module d8c39c3f) → test630dd (39/39). 612dd and 630dd commute.
4. 642 `experiments/642-hek293-contingency/patch_hek293_642.py --root .` (--check, apply; default OFF; stage-2 floor 0.75) → its tests.
4a. 604 `experiments/604-identity-prior/patch604.py` (--check, apply; installs scripts/def_identity_prior.py + the identity path in
   the POST-BUNDLE def_sameline.py 22c6f0f0; unset-inert) → test604 (27/27 since a43eae97: + T6b forced same-line match,
   T6c 626 sister guard) + test604_guard.py (5/5). The October identity hook needs 604 + 612 + 612dd.
4b. 644's `experiments/644-runbook-amendments/patch_skip644.py` **v3** (sha256 65fab4e9…, d4a3d275; v2 4c0bb12e and v1 11716c46
   superseded; 607's whitespace + log fixes): `DEF_DAY_SAMELINE_SKIP=<ctx>`
   → def_sameline `--skip-contexts` (a skipped context takes the exact unmatched path; with a programme source it keeps the
   source's programme at w 0). Touches def_day.sh AND def_sameline.py; anchors on 642's lines, so it applies last, after 604.
   test_skip644 v2 4/4 on 606's full-chain trees. Final pins are in the 606 block below (def_sameline f147f8bb after 607's fixes).
   **Composition (606, [4d0207], pushed):** 641 (10-05) → 598 → 612 → 612dd → 630 → 642 → 604 applies clean on test598's canonical
   post-bundle tree; 604 and 642 commute, as do 644 and 604. Re-run with 4b (skip644 v2) as the last link: ALL PASS (1b177a5c).
   **Pins: the authoritative table is the LAST "PINNED" block in `experiments/606-applyday-composition/RESULTS.md`** ([4d0207]
   maintains it). The apply checks every touched file against the block current at the freeze. As of f152e4c6:
   - def_day.sh 7075d83d, builder fe4678c3, def_sameline f147f8bb, def_identity_prior 7aeb66b4, src/vcc2026/rank_up.py 7e6b3e52;
   - def_identity 7073406f, helpers fb715e38, def_pipeline 76314af6 (= post-598), october_recipe 9dbcd914, def_preflight
     17f0836f, def_postbuild 557c2b52, gate524 b16e60c0.
   **Chain freeze: 10-05 18:00Z.** After that, patch changes wait for 10-08; the apply uses the pin current at the freeze.
   Re-run compose606.py after any further patch change, and re-pin.
5. Builder suite (`experiments/410-pertscale-revisit/run_builder_tests.py`) + every patch's suite green.
6. KA builds (run_ka_527.sh style, serial, memory gate): ONE all-unset KA on the composed builder == `ec3a5ae6` (each patch is unset-inert, so
   one proves all); KA598b tercile == 479p2g's `dd82740d`; KA598c tercilec (`ka598.py` PASS); KA612b clean4:0.05 == `e21ff752`; KA612c
   clean4:0.05 + KEYGATE 0.5 (`ka612.py` PASS); KA630 K 1.0 (config diff exactly `de_uprn`, AND 630dd's audit shows 0 fail-safe
   triggers: the line `630dd: uprn FAILSAFE_TOTAL <n>` must be PRESENT and n == 0. A missing line, n > 0, or n = −1 (missing or
   unreadable dump; 607b hardening, bd76246b) = KA FAIL). ≈ 6 builds, ~4 h.
   607's "not checkable here" items:
   - (1) the all-unset KA == ec3a5ae6 above is exactly its "same recipe, same hash vs the pre-chain builder".
   - (2) Before the apply: H1 per-target min(ballast lib) > max(carrier lib) on the seed-0 dumps ([d4cd5f]), i.e. 630's bug
     condition never occurs.
   - (3) Count the priced 612 targets that tie at the key-gate cut ([21dd4c]).
7. run_arm_generic.sh driver list += def_pipeline's post-598 sha (from post598). bash -n; test_gate591.
8. Commit every touched file with the KA hashes in the message; tell [537757] the final shas (644's card).
9. `october_recipe.py decide sharpen_pt tercilec:2.5,2.0,1.75`.
- 10-02 12:13Z: **[VOID at 12:58Z: call-set-blind fold model, 612 CORRECTION 2]** **October hook form, decision table ([21dd4c], 612/kgh612c, delivery-faithful folds; Δnmae vs no hook):** plain CLEAN4 −0.0006 / +0.0149 / +0.0148 / −0.0102 (h1 / HCT116 / HEK293T / K562); key-gated CLEAN4 −0.0001 / +0.0142 / +0.0128 / +0.0013; **key-gated IDENTITY (622's z(P_C4 − 0.5) + 0.5 z(I), Flex refs, weak half) +0.0057 [+0.0015, +0.0099] / +0.0210 / +0.0203 / +0.0003: dominates the plain hook on every surface**, ≈ +0.0017 … +0.003 overall over it (×1.2–1.6 per 611). Leading candidate = 604's identity artifact as the hook source + VCC_DE_MDE_SIGN_KEYGATE=0.5 (612 + 604 plumbing, no further builder change). **Gate before any slot: the H1 REAL-EMITTER idkg seeds must confirm it** (the analytic 612 numbers did not survive the emitter); then a 10-08+ preregistered arm read against 496a.
- 10-02 12:34Z: **idkg seed 0 on H1's REAL emitter does not show dominance** ([21dd4c], diag612b). idkg − plain hook, overall-eq: f 1 −0.0015,
  f 0.5 +0.0007, f 0.35 +0.0005, cmp +0.0008 (one seed). The 12:13Z "dominates on every surface" was the fold reading. Diagnosis: the
  identity source helps where CLEAN4 is blind (sign accuracy 0.507 → 0.541 there); the key gate gives up the plain hook's
  shrink-where-sources-disagree magnitude benefit. [VOID at 12:58Z] The folds still favour KGID: K562 tail −0.010 for the plain hook vs −0.0013 (ns).
  PREREGISTRATION-draft-KGID.md (55451545) now calls it an EV-vs-tail choice, with its TIE bin resting on min-max regret.
  Seeds 1–2 are due ~13:00–13:15Z; the prereg freezes after them. 10-06 APPLY DAY gains 612dd and 604 (steps 2b/4a) and, conditionally,
  644's SAMELINE_SKIP (4b), in 606's tested order.
- 10-02 12:43Z: **KGID 3-seed H1 real-emitter table** ([21dd4c]). idkg − c4h overall-eq (sd over seeds): f 1 −0.0014 (0.0009),
  f 0.5 +0.0008 (0.0001), f 0.35 +0.0006 (0.00005), cmp +0.0011 (0.0004). Raw nmae −0.0032 / −0.0002 / +0.0007 / −0.0037 and reach
  −0.0121 / +0.0037 / +0.0042 / +0.0006. fid and jac within ±0.0001 everywhere (delivery-class).
  - Decomposition: the key gate costs nmae on H1 (+0.0003 … +0.0013) and the identity source pays it back (−0.0003 … −0.0033).
    KGID vs no hook: nmae −0.0015 … −0.0046 in every regime.
  - Predicted board primary Δ(fid+nmae) KGID − 496a: −0.0003 … +0.0001 overall at f 0.35–0.5, so the TIE bin is most likely.
  - Disagreement shrink is CLOSED on the fold model (FINDINGS §2).
  - **Open validity issue:** the fold model gets KGID − plain's f 0.35 sign wrong on h1 (fold −0.004 for plain; real emitter
    +0.0007 for plain). The KGID TIE rule's min-max-regret basis (the K562 tail) comes from that model. [21dd4c] is building a
    call-set-aware fold model (ds612c). Its retrodiction criterion and the fallback are fixed BEFORE its fold rows are read
    (see the 12:44Z message to [21dd4c]).
- 10-02 12:53Z: **KGID prereg amendment A1** (3efebfd7, committed 12:52:10Z; verified: no ds612c file existed yet).
  - Retrodiction criterion (a) f 1 kgID − plain in [−0.0063, −0.0016]; (b) f 0.35 > 0; (c) f 0.5 ≤ +0.0002 and monotone in f.
  - Mechanism check (d): kgC4_05 − plain > 0 at f 0.5 and f 0.35 (the gate costs nmae on the real emitter).
  - Known answer: an empty call set must reproduce ds612b exactly.
  - Fallback on FAIL: TIE and BELOW → plain hook (A/B/C and October); KGID needs WIN.
  - On PASS: TIE → the smaller max regret over 4 folds × {f 0.5, 0.35} (equal → plain). BELOW → plain on A/B/C; October takes KGID only
    if it leads on all 4 folds at both regimes and the board gap is ≤ 0.0010. kgh612c's board expectation is SUSPENDED until A1 passes.
- 10-02 12:53Z: 645 (HPA references, [537757], 191ade02): 642 unchanged; the broad non-CCLE veto panel is rejected (FINDINGS §2).
- 10-02 12:58Z: **CORRECTION: the hook-class fold evidence is VOID; KGID's A1 FAIL branch is in force** ([21dd4c], 612 CORRECTION 2, 22691ce5).
  - The call-set-aware fold model ds612c fails A1's retrodiction:
    - f 0.35: −0.0027 vs the real +0.0007 (wrong sign);
    - f 0.5: −0.0019 vs −0.0002 (not monotone);
    - f 1 passes; (d) passes.
  - Root cause: kgh612b/c/d, qscan612 and hyb622b hooked EVERY gene, while production hooks only uncalled genes
    (rc_emit_k.py:1059 `use_s = ~head & ~blk`). Masking the 2,600 called genes shrinks kgID's fold gains:
    - HCT116 +0.021 → +0.0013; HEK293T +0.020 → +0.0022;
    - the plain hook's K562 "tail" becomes −0.0013 … −0.0019.
  - **Consequences:**
    - (1) KGID needs a board WIN over 496a (Δ(fid+nmae) ≥ +0.0002). H1 predicts −0.0003 … +0.0001, so P(WIN) is low; the
      KGID read ranks BELOW every arm with a validated mechanism.
    - (2) The plain CLEAN4 hook's October case (mde_sign, decided by 496a vs 446bj2 on 10-05) never rested on the folds:
      its held-out read is H1's real emitter (3 seeds, nmae −0.00146; october-decisions.json). It stands as preregistered.
    - (3) The identity axis keeps its sign-accuracy result (639, 622g). Its delivery through the hook is measured (≈ the plain
      hook); delivery through the BLOCK prior is unpriced (611) and is now the identity axis's open question.
    - (4) patch612 / 612dd stay applied on 10-06: tested, inert-when-unset engineering.

- 10-02 13:08Z: **Queue after the 10-07 / 10-08 tree (coordinator).** 537 §5's tree keeps 10-07 and the 10-08 branch (CPH, K20 / MZ,
  459dr, ALLSW). After that, in this order:
  1. 630's dose pair (K 0.75 / 1.0; exact members, nmae only).
  2. 639 **ID05**, the identity BLOCK prior, once its prereg freezes on the four review items; **ID10** only if ID05 ≥ WEAK.
  3. The factorial-validation / all-switches arm.
  4. 524dr2.
  5. The KGID read (A1 FAIL branch: needs a WIN).
  6. 625's module readout (8 modules + 2 additivity arms, last; it doubles as the final-phase decode rehearsal).
  Each arm is armed with its prereg's base = the best scored entry at arming time.
- 10-02 13:11Z: **639 ID05/ID10 artifacts pinned** (d5f7b65b; shas re-verified by the coordinator), all in
  `experiments/604-identity-prior/work/`:
  - ID05 `ident604v3_flex_lam0.5_ABC.npz` 6cd84dc258cb609d…
  - ID10 `ident604v3_flex_lam1_ABC.npz` 3b1e22f037aff26b…
  - KA (λ 0 == CLEAN4) `ident604v3_flex_lam0_ABC.npz` f65d12fc3a7d…, byte-identical to v2's.

  The prereg freezes when the numbers for its null-block, 458d-retrodiction and Broad-16-by-lineage readings exist (rules
  fixed at d5f7b65b).
- 10-02 13:13Z: **630 has a real bug** (605 Grok review, triaged by [4d0207], 539c5dfa; `src_rank_up.py:190`).
  - When a carrier's library is ≥ the ballast libraries, the donor-gain clamp leaves the ask unbounded.
  - `fill` also takes from negative-gain donors, and nothing reverts a gene whose mean fell. An UP lfc can turn DOWN with every
    invariant green.
  - [d4cd5f] holds the fix spec (donor filter + mean/sign revert). 630 joins the 10-06 apply only if the fix, test630, its KA1,
    the H1 re-price (both seeds, reported against the pre-fix numbers) and 606's re-run all land by 10-06 06:00Z. Otherwise it
    applies 10-08.
  - 604 v3 replaces v1 (identity-only genes get no prior; side = sign(zc + λ·zi)); def_identity_prior a2d48da0.
  - chain1003 pre-check ([4d0207]): 495 and 471c pass every ready() condition, so s1 495 / s2 471c at 20:00Z.
- 10-02 13:15Z: **Queue re-ordered after 613.** 613 shows ID05 ships only on TRANSFERS, and pre-read NULL/HARM is about as likely as
  TRANSFERS. The factorial validation has a deadline: it must read before Oct 22 for the final-phase reads. New order after the
  10-07 / 10-08 tree:
  1. 630's dose pair (after its fix);
  2. the factorial-validation / all-switches arm;
  3. 639 ID05 (ID10 only after a TRANSFERS read);
  4. 524dr2;
  5. KGID;
  6. 625.
  - 612's tie-at-cut fix ([21dd4c]; no change on tie-free keys) changes patch612's sha, so 606 re-pins before the 10-06 apply.

- 10-02 13:21Z: **Queue re-ordered again ([d4cd5f]'s 630 / 639 status).**
  - 630's fix is in (69972aa1, module d8c39c3f; test630 31/31, test_neggain630 11/11). The old code could not flip a gene on H1:
    ballast ~995k vs carriers ≤ 25k, so every donor gain was positive. A bit-for-bit re-run of the 9 priced H1 runs is in progress.
  - Under v3, H1 cannot separate K 0.75 from K 1.0 (+0.0030 / +0.0029 at f 0.35). K 0.5 and K 1.5 report ~14:00Z, then one arm or a
    pair H1 actually separates.
  - 639 ID05 is DEMOTED: route unresolved; only TRANSFERS ships. Evidence:
    - 458d retrodiction passes;
    - null-block precision below CLEAN4 on 2 of 3 folds;
    - Broad-16 block mean −0.0022, 7 of 12 negative, epithelial like B/C;
    - block route seen on 2 lines vs 13 elsewhere.
  - KGID drops to empty-slot-only. Its read is delivery-class, so 611's validated H1 multipliers already predict it (TIE/BELOW →
    plain); the slot would buy almost no information. ID05 stays above it: a prior-type change H1 cannot price, so the board is
    the only way to learn it.
  - **Order after the 10-07 / 10-08 tree:**
    1. 630 (one arm, or a separated pair);
    2. 524dr2 (the applied October recipe end to end; early, so a defect leaves time to fix);
    3. the factorial validation (deadline Oct 22);
    4. 625 (8 + 2 arms);
    5. ID05;
    6. KGID (empty slots only).
- 10-02 13:29Z: **10-05's 496a read also discriminates two hook instruments** ([21dd4c], KGID A2 3c6f0a80 / A3 2c24e322, committed
  13:25Z / 13:28Z; idkg3's run started 13:24Z with no stats yet).
  - **H1 real emitter** predicts a clear improvement for the CLEAN4 hook: Δnmae −0.0025 … −0.0034 raw.
  - **613d** (426's add rule on uncalled gate genes; known-answer rows on H1 within 0.0004) predicts ≈ 0 on cancer lines:
    board Δnmae +0.00017 raw.
  - Rule, fixed now: if 496a's Δnmae < −0.0005 raw, 613d is closed for hooks and A3's PASS branch is void. mde_sign is still
    decided by 496's own preregistered bins; this reading is extra and must be written in 496a's BOARD-RESULT.
  - A3 test: real idkg3 − c4h within ±0.0006 of 613d's −0.00202 / −0.00034 / +0.00006 (f 1 / 0.5 / 0.35). On PASS, 613d's cancer
    folds become KGID's held-out read; on FAIL, A1's FAIL branch stays. KGID now names 604 v3 (6cd84dc2); idkg3 seeds ETA ~13:45Z.
- 10-02 13:37Z: **Instrument card for the K16 (10-06) and K20 (10-07 s2) reads:** `experiments/610-board-forensics/PREDICTIONS-611-K16-K20.md`
  (a995f0e2, committed 13:36Z before either entry exists). 611's call-set multipliers on 537's H1 rows predict:
  - K16: Δfid +0.0069 raw (+0.0034 … +0.0097), Δjac +0.0010, Δ(fid+jac) +0.0043 overall.
  - K20: Δfid +0.0070, Δjac +0.0012, Δ(fid+jac) +0.0044 overall.
  - The readers grade it alongside 523/537's slot bins, which it does not change. A board Δfid < 0 means every f 0.35 call-set
    price is wrong-signed; board/H1 fid ≈ 2 confirms ×2 beyond cp, ≈ 1 means size changes do not inherit it.
- 10-02 13:40Z: **Pre-625 KAs** ([4d0207]). These run before 625's first A/B/C arm, NOT on 10-06:
  - (1) λ 0 + no flip + `VCC_DE_MDE_SIGN=programme:τ` reading the artifact's pdown == `clean4:τ`, byte for byte;
  - (2) the same with 612's KEYGATE set.
  The flip arms take the hook from the artifact, so both must pass. def_identity_prior is re-pinned to 60e98488; test604 on it is
  queued, and [4d0207] confirms.

- 10-02 13:43Z: **630 is ONE arm: b_k100_q20r3p10** (module 7e6b3e52; prereg amended at 27621cc0 before any entry). Board Δnmae raw
  −0.004 … −0.015, center ≈ −0.009; bins unchanged.
  - v3 dose bracket, 3 seeds, f 0.35 overall-eq: K 0.5 +0.0020, K 0.75 +0.0030, K 1.0 +0.0029, K 1.5 +0.0025. It is flat above 0.75
    because carriers saturate: they hold 0.55 % of their library in sub-gate genes, the pass spends 93 %, and 76 % end fully drained.
  - The freed slot goes to 524dr2.
  - The full-delivery oracle is +0.0060 vs +0.0029, so a new compensation currency would be worth ≈ +0.003 on H1. [d4cd5f] is
    measuring one (v4: DOWN-gene counts parked on zs's mass holders, paid carrier → ballast). It enters 10-06 only if pinned
    before the 10-05 18:00Z freeze, otherwise 10-08, with a prereg amendment before its arm.
  - Nine bit-for-bit KAs on 7e6b3e52 are queued on gate602.
- 10-02 13:46Z: **Candidate hook form idall3: UNGATED v3 identity hook** (`VCC_DE_MDE_SIGN=file:<ident604v3 λ0.5>:<sha>:0`, no KEYGATE;
  [21dd4c] 613d4/613d5). Prediction registered at 633f3a39 before its real-emitter run: H1 idall3 − c4h = −0.00160 / −0.00232 /
  −0.00225 ± 0.0006 raw nmae (f 1 / 0.5 / 0.35); 3 seeds in chain612d.
  - **PASS** → idall3 REPLACES KGID in the queue; KGID retires, since idall3 beats it everywhere but K562. Its board prereg vs 496a
    must be written and committed BEFORE 10-05's 496a read, with the K562 stress tail (+0.0005 … +0.0024 raw) stated.
  - **FAIL** → the hook family stays as is.
  - Caveats:
    - Its fold rows come from 613d, which is still unvalidated out of sample. If 496a closes 613d for hooks, idall3's October case
      rests on H1 + board only and ships only with a written reason.
    - Hook forms tried so far are ≥ 8 (c4h, c4kg, idkg, idkg3, idall3 and the shrink variants), so its H1 edge is best-of-N biased
      (11.1c.4). The board read is the unbiased one.
- 10-02 13:48Z: **ID05 → empty-slot only** (609: the Orion block route is a depth/expression artifact; expected ≈ +0.000 … +0.001,
  below its TRANSFERS bar). **Queue after the 10-07 / 10-08 tree:**
  1. 630 (one arm);
  2. 524dr2;
  3. factorial validation;
  4. idall3 (only if its registered H1 prediction passes; its prereg is committed before 10-05's read);
  5. 625.
  Empty slots only: ID05, then KGID (retired if idall3 passes).
- 10-02 13:54Z: **607 (Grok review pinned 7dbb1007): the 10-06 chain HOLDS**, no blocker and no major (2ece7b66). The reviewer rebuilt the
  chain and reproduced all twelve pins. 604's three items are fixed, with byte-identical artifacts.
  - Re-pin: def_sameline f147f8bb, def_identity_prior 7aeb66b4.
  - 644's two MINOR items (strip() on DEF_DAY_SAMELINE_SKIP; log "unmatched (644)") were relayed to [537757].
  - 604f, the day-of flip pass-through (`--flip-octiles` / DEF_DAY_DE_PRIOR_FLIP), is for 10-08+ only, outside the frozen chain.
- 10-02 13:57Z: **625 Amendment 3** ([7772dc], fe2064b7, before any A/B/C read): 9 prior-only entries (8 single-octile probes masks
  10000000 … 00000001, then 1 combination with the decided mask), still LAST and empty-slot-only from 10-08.
  - Reader: `readdecode625.py --design single --rule transfer3` (selftest PASS, 899d1534). nmae alone decides; fid is a no-harm
    check (DOWN→UP flips grow the realised call set); reach is excluded.
  - The pre-625 KAs (λ 0 + no flip == clean4:τ, ± KEYGATE) gate its first arm.
  - D/E/F: no decode; one verification entry only if GO (FINAL-PHASE-FEEDBACK-PLAN §6 update).

- 10-02 14:26Z: **630 READY for the 10-06 apply** ([d4cd5f]; rank_up 7e6b3e52).
  - KA 9/9 bit-for-bit vs every priced H1 run, so all prices stand; the pre-fix code never flipped or lowered a gene.
  - A/B/C exposure KA on 446bj2z's real blocks: negative-donor rows, mean reverts and guard failures are all 0; capacity is unchanged.
  - Pins as in 606's latest PINNED block. The arm is ONE entry, VCC_DE_UPRN=b_k100_q20r3p10 on a zs-carrying base (prereg + 27621cc0).
- 10-02 14:28Z: **633 (held-out UP-dose reading; rules fixed before the run, e17f0f60): October's UP dose is REOPENED, conditional on a
  second clean cancer surface.**
  - Retrodiction: the fold model reproduces 630's real-emitter oracle in shape (signs 8/8, same best K at f 0.35, same f 1
    overshoot); its size runs ~2× at f 1.
  - K562 at A/B/C depth: L1-optimal K is 2 (f 1) and 3 (f 0.5 / 0.35), against H1's 1–1.5. This matches the sign accuracy of
    uncalled-UP genes: K562 0.73 / 0.89 vs H1 0.54 / 0.58.
  - 10-06 is unchanged: the single arm stays b_k100_q20r3p10, because v3 saturates above K 0.75 and cannot deliver K 2.
  - It raises the value of capacity work: at full delivery K562 gains ~2.5× K 1.0's Δ. [d4cd5f]'s v5 (DOWN-gene currency) is its top
    item; 10-06 only if pinned before the 10-05 18:00Z freeze.
- 10-02 14:50Z: **10-08+ apply list (outside the frozen chain):** 604f `patch604f.py` (8465eb14; test604f 24/24). It adds def_sameline
  `--flip-octiles` and DEF_DAY_DE_PRIOR_FLIP, applies after the full 10-06 chain (needs 630dd + 644 v3; composes ALL PASS), and
  is unset-inert. Apply it ONLY if 625 reads GO on A/B/C: its mask then enters R0 through this knob. test604 is 35/35 on the
  pinned 7aeb66b4.
- 10-02 14:53Z: **KGID A3 PASSES** (c964229b): 613d's registered H1 predictions held within 0.0003 in every regime.
  - KGID's A1 PASS-branch rules now apply, with 613d's cancer folds as its held-out read:
    - TIE → KGID (max regret 0.00006 vs the plain hook's 0.00204 on K562 at f 0.5);
    - BELOW → plain;
    - WIN → KGID.
    All of this is void if 496a reads < −0.0005 on 10-05.
  - idkg3 on H1 vs c4h, overall: +0.0009 / +0.0006 / −0.0001. Dropping v1's identity-only extension recovered its f 1 reach loss
    (+0.0138).
  - **Queue change:** the hook-form arm is back in the main queue. A board read is now needed to switch October's hook form
    (11.1c.2), and a TIE would move October to KGID. Slot 4 = idall3 if it passes its registered H1 test, else KGID.
  - Order after the 10-07 / 10-08 tree:
    1. 630;
    2. 524dr2;
    3. factorial validation;
    4. the hook-form arm;
    5. 625.
    Empty slots only: ID05.
- 10-02 14:58Z: **October UP dose stays K 1.0.** 633j (Jurkat, Nadig 2025) was meant as the second clean cancer surface, but it fails its
  preregistered cleanliness check (block gate genes 70 % DOWN; limit 0.62). So there is no confirmation. Descriptively it agrees with
  K562 (d20 optimum K 1.5–3), and it is not re-read post hoc. 630 v4 is closed, NULL (FINDINGS §2). v6, the last capacity probe, is
  queued on H1; if it fails, 630 ships as is.
- 10-02 15:22Z: **idall3 PASSES its registered H1 test and takes queue slot 4; KGID is RETIRED.** Its board prereg is FINAL
  (612/PREREGISTRATION-draft-IDALL3.md, cf3f8004), committed before 10-05.
  - H1 real emitter, 3 seeds, vs c4h, overall: +0.0014 / +0.0011 / +0.0006 (f 1 / 0.5 / 0.35). reach +0.0045 / +0.0029 / +0.0002; fid and
    jac flat. vs KGID v3: +0.0005 … +0.0007.
  - Arm: best base + `VCC_DE_MDE_SIGN=file:experiments/604-identity-prior/work/ident604v3_flex_lam0.5_ABC.npz:6cd84dc2…:0`, no KEYGATE.
    The config diff must be exactly `de_mde_sign`. October form: `DEF_DAY_MDE_SIGN=hookbank:identity:0` (604, applied 10-06).
  - Predicted Δ(fid+nmae) vs Δ_496a ≈ +0.0007 … +0.0010 overall; reach band ±0.0042 raw.
  - Bins:
    - WIN (≥ +0.0002) → idall3. If 496a closes 613d, October needs a written reason (11.1c.2).
    - TIE / BELOW → plain.
    - HARM → no identity hook.
  - Stated risks: the K562 stress tail (+0.0005 … +0.0024 raw); Broad-12 KO evidence is mixed (5–8 of 12 lines worse, ~10–15 %
    K562-like); 613d's best-of-5 selection.

- 10-02 15:44Z: **630 stacks additively on every hook that may ship** ([d4cd5f], H1 seed 0, 0 calls changed). Overall-eq at
  f 1 / f 0.35 / cmp:
  - on the plain hook: +0.0022 / +0.0027 / +0.0034;
  - on c4kg: +0.0024 / +0.0030 / +0.0037;
  - on idkg3: +0.0028 / +0.0029 / +0.0038.
  Better signs compound mildly. 630's board read stays clean on any hook base, since the knob is the only config diff. idall3's stack
  (3 seeds) is queued.
- 10-02 16:00Z: **630's capacity line is CLOSED; 630 ships as v3** (FINDINGS §2). 639 is re-sized with 625 FINAL's 0.35 factor:
  ID05 +0.0005 / +0.0011, ID10 +0.0009 / +0.0017 (H1 f 1 / f 0.35); NULL stays modal. **Coordinator amendment, before any read:**
  if an empty slot opens for the identity block prior, run ID10 FIRST. It has the larger dose and more power to show whether the
  block route exists on Arc lines at all. ID05 runs only if ID10 reads TRANSFERS, to pick the dose. [d4cd5f] amends the prereg.

- 10-02 16:22Z: **631f (held-out dose for the called-DOWN lever) FAILS its retrodiction.** The fold model's H1 row puts m* at 1–1.5; the
  real emitter says 2. Its base cannot represent the emitter's DOWN under-delivery.
  - Descriptively, the cancer folds want ~2× H1's called-DOWN magnitude (m 2.5–3), the same direction as 633.
  - Mechanism: the mass-balance kf (~0.60) deliberately keeps DOWN genes less down, because freed mass would inflate other
    genes' CPM and make false calls (610). Magnitude can only grow where freed mass is invisible: sub-gate pools, or v7's
    pairing with UP receipts.
  - **614 assigned (offered to [21dd4c]; [d4cd5f] after v7 if declined): a REAL-EMITTER DE harness on K562.** It is the
    cancer-line held-out DE surface 534 flagged as missing. It must pass a 611-style retrodiction on the board's known DE pairs
    (zs, coarse MDE, cp) before it prices 630's dose, 631/v7 and the hook forms. Target for the trust verdict: 10-04 05:00Z.

### 10-05 CHECKLIST v2 (10-02 16:40Z; 650's tabletop fixes B1, L1–L7 folded in; run in this order, every command from the repo root)
0. Wait for `recipe_dr524.sh` (builds 524dr from 02:00Z) to EXIT: `bash experiments/550-cloud-1005-review/repro/procs_clear550.sh`
   prints CLEAR, and no run_arm_generic.sh is running. Live == HEAD on the 11 paths (`git diff --quiet -- <paths>`). HEAD stays the
   pre-bundle commit until step 9.
1. **BEFORE ANY APPLY, run every bundle suite.** They build temp trees from the LIVE files and re-apply the patches themselves, so
   after the apply they die "already patched" (B1):
   - `uv run python experiments/590-identity-validation/test_identity590.py` → ALL PASS 8/8 (0.772 / 0.719 / 0.733);
   - `uv run python experiments/538-cloud-bundle-extras/test_extras538.py` → 13/13 (the L3 fix is in; verified 10-02);
   - `uv run python experiments/554-dose-guard/test_guard554.py` → 5/5;
   - `uv run python experiments/554-dose-guard/identity554.py` → prints `all identical: False`, which is CORRECT (preflight 6/6
     identical; postbuild 11/12, 465efr adds its dose-guard FAIL row at prior 0.278; `--rc-audit` on 433a adds no row, 429a adds
     FAIL then ALLOWED). It rewrites the tracked `experiments/554-dose-guard/work/identity554.json`:
     `git checkout -- experiments/554-dose-guard/work/identity554.json` afterwards;
   - `uv run python experiments/641-depth-row/test_depth641.py` → 4/4.
2. 538: `uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root . --check`, then `--root .` (def_pipeline →
   9a718b9e). Never run `patch_bundle530.py` separately.
3. 554: `uv run python experiments/554-dose-guard/patch_guard554.py --root . --check`, then `--root .` (def_postbuild 557c2b52,
   def_day 5851f84b). `--root` is required; 554 refuses before 538.
4. 590: `uv run python experiments/590-identity-validation/patch_identity590.py --check`, then `--apply`; then
   `uv run python scripts/def_identity.py --verify` → VERIFY PASS 0.772 / 0.719 / 0.733.
5. `bash -n scripts/def_day.sh`; `uv run python scripts/october_recipe.py check --allow-pending` → rc 0.
5b. (not a gate) `scripts/test_october_recipe.py`: after 538, 4 tests fail (all_switches_on, riders_need_502, dependency_violations,
   empty_flags) because 538 adds six PENDING switches. Add `("rho_scale", "off"), ("rc_callset_cpm", "off"), ("hook_estimator",
   "off"), ("lodo_strength", "off"), ("lodo_hi", "off"), ("callset_block", "off")` to those four tests' `self.decide([...])` → 6/6
   (650-verified).
6. `uv run python experiments/594-def-bundle-composition/post594.py` → ALL PASS.
6b. 641: `uv run python experiments/641-depth-row/patch_depth641.py --root . --check`, then without `--check` (def_preflight 17f0836f,
   def_day 738cb0ae). Its suite already ran in step 1.
7. Inert KA by hand (550 RESULTS "10-05 command sequence" step 4): efp relaxed, fresh work dir, audit dump on.
   - Clear EVERY DEF_DAY_* knob, including the five that 550's list misses: `DEF_DAY_ALLOW_CELLS_PER_PERT`, `DEF_DAY_ALLOW_DEPTH`,
     `DEF_DAY_ALLOW_LOW_DOSE`, `DEF_DAY_ALLOW_WIDE_LIBS`, `DEF_DAY_INCUMBENT_META` (`env -u` each, or `env -i HOME=$HOME PATH=$PATH …`
     plus the relaxed setting).
   - The output must be `ac49d319…` (526/efp_inert/474efp_final.h5ad.meta.json); postbuild.json carries `dose_guard` per context;
     the log has `DONE. Next`.
   - NEVER re-run chain502 / fix526 / chain476.
8. `experiments/428-zero-gene-audit/run_arm_generic.sh`: add `9a718b9eaa24f10cfbf255ce830ea4a748f61252ec503bbb6db69d3ddf5ca9fd` to
   its driver case list (line ~58; only while no run_arm_generic.sh runs); `bash -n`;
   `bash experiments/591-prep-gate/test_gate591.sh` → ALL PASS.
9. Commit with a pathspec (AGENTS §10): the 11 files + run_arm_generic.sh + scripts/test_october_recipe.py, with the KA hash in the
   message.
10. `october_recipe.py decide` the switches the reads settled (lodo_hi_pct 10-03, mde_sign 10-05). `command` still refuses on six
    PENDING switches (rho_scale, rc_callset_cpm, hook_estimator, lodo_strength, lodo_hi, callset_block), so the ONE rider-on
    rehearsal on efr uses `command --allow-pending`, which prints pending switches as OFF.
11. 20:00Z: chain1006 arms 10-06 (524dr + K16) automatically. Before 10-06 00:01Z: write the 10-07 rule (537's tree, v5+).

### 10-06 APPLY DAY v2 (10-02 16:40Z; 650's B2–B4 and L8–L11 folded in)
**Preconditions, both before the 10-05 18:00Z chain freeze:**
- **(P1) DONE** (26c8eddb; proof 980c0de2 on a simulated post-10-05 HEAD: test604 35/35, test604f 24/24, test_ka630_audit 11/11).
  Every chain-tree builder skips a patch whose `--check` says it is already in the tree ("partially" still fails); verify644 has
  the same rule (3ea6debe).
- **(P2) DONE:** `uv run python experiments/630-rank-neutral-up/ka630_audit.py <build>.rc_audit.jsonl` (26c8eddb, 9e239a41). It runs
  630dd's own snippet and exits 0 only if `630dd: uprn FAILSAFE_TOTAL 0` AND every context clears 630's capacity floor (mean
  ur_frac_got ≥ 0.127 and ur_n_moved ≥ 500 per target).

**Steps** (pins: compare every touched file against the LAST PINNED block in `experiments/606-applyday-composition/RESULTS.md`, as of
the freeze; all 12 must match):
1. 598: `uv run python experiments/598-tercilec-knob/patch598.py --root . --check` → `test598.py` (43/43; it must run BEFORE the apply,
   because T0 compares the live builder) → apply `--root .` → `post598.py` (ALL PASS; def_pipeline 76314af6).
2. 612: `patch612.py --root .` (--check, apply) → `test612.py` (23/23).
2b. 612dd: `patch612dd.py --root .` (--check, apply) → test612 (23/23).
3. 630: `uv run python experiments/630-rank-neutral-up/patch_builder630.py --root .` (--check, apply; installs src/vcc2026/rank_up.py
   7e6b3e52) → `test630.py` (35/35) + `test_rank_up630.py 0 b_k100_q20r3p10 85` (PASS bit for bit, ~5 min, ~4 GB).
3b. 630dd: `uv run python experiments/630-rank-neutral-up/patch630dd.py --root .` (--check, apply; sha per the PINNED block) →
    `test630dd.py` (50/50; needs P1).
4. 642: `patch_hek293_642.py --root .` (--check, apply; default OFF) → `test_hek293_642.py` (needs P1).
4a. 604: `patch604.py` (--check, apply). It prints `WARNING def_sameline.py is ff898053…, not the post-bundle 22c6f0f0`, which is
    BENIGN in this order (642 ran first; 604 and 642 commute per 606). → `test604.py` (35/35) + `test604_guard.py` (5/5).
4b. 644 v3: `patch_skip644.py --root .` (--check, apply) → `test_skip644.py` (4/4; needs P1).
5. `uv run python experiments/410-pertscale-revisit/run_builder_tests.py` → 56 passed / 0 failed / 1 skipped; plus every suite
   above green.
6. KA builds (run_ka_527.sh style, serial, memory gate):
   - all-unset == `ec3a5ae6` (446bj2z meta; this is also 607 item (1));
   - KA598b tercile == `dd82740d` (479p2g);
   - KA598c tercilec (`ka598.py` PASS);
   - KA612b clean4:0.05 == `e21ff752` (496a);
   - KA612c clean4:0.05 + KEYGATE 0.5 (`ka612.py` PASS);
   - **KA630:** `VCC_DE_UPRN=b_k100_q20r3p10` on the all-unset recipe (446bj2z, zs on), with `VCC_RC_AUDIT_DUMP=<fresh path>.jsonl`
     (the builder APPENDS, so delete a stale file first). Its config diff vs the all-unset KA's meta must be exactly `de_uprn`, and
     `ka630_audit.py <dump>` must exit 0 (FAILSAFE_TOTAL 0 and every context above the capacity floor). Anything else is a KA FAIL.
   - 607 item (2) is DONE (650: on 85/85 H1 blocks the ballast library exceeds the carrier library by ≥ 968,667).
   - 607 item (3) is DONE: 0 of 85 priced H1 targets tie at the KEYGATE q 0.5 cut, and the f5ecaf2b fill changes no weak set
     (`experiments/612-key-gated-hook/ties612.py`, d1e8cb1f). Ties now break deterministically; KA612c covers the rest.
7. `experiments/428-zero-gene-audit/run_arm_generic.sh` driver list += the full post-598 def_pipeline sha (from post598); `bash -n`;
   `bash experiments/591-prep-gate/test_gate591.sh`.
8. Commit every touched file with a pathspec, KA hashes in the message; send [537757] the final shas.
9. **DEFERRED (652 F4):** do NOT decide `sharpen_pt tercilec` on KAs alone. It is decided only after an H1 real-emitter
   tercilec-vs-tercile run with a coverage-thinned key ([36a81a]; the DE side-effect of tier moves, which the bulk instrument
   cannot see) AND C10b's read (10-04). C10b's claim row IS tercilec's held-out bulk read: designs551's P2g = the carried-only
   form = 598's tercilec, priced vs uniform q 1.75 on R6's 12 sealed folds ([10eca4]; written reason in 551 RESULTS, f96f8a4e).
   No extra C10b row. Off-reserve form comparison: 493h efp, carried +0.84 % rel cos (3/3 folds) vs plain +0.48 % (2/3).
   **C10b DONE (10-04 07:24Z): P2g@M6 TRANSFERS on sealed R6** (+0.00208 overall-eq, 11 / 12 folds; board projection +0.00093). With 599
   clean, step 9 is UNBLOCKED: decide `sharpen_pt tercilec:2.5,2.0,1.75` on 10-06 once KA598b / KA598c pass, citing C10b
   (f6aa2ed9) + 599 (7db7938c).
   **599 DONE (20:24Z, 7db7938c): R1, the DE side is clean.** On H1's real emitter, 5 seeds, coverage-thinned (43/85 non-carried; plain
   tercile reproduces efp's empty MID): tercilec − tercile DE overall-eq at f 0.35 = +0.00002 (SE 0.00006), nmae +0.00002.
   Bulk: pds +0.0029 (5/5), mse +0.0006 (slightly worse), ≈ +0.0010 overall-eq. Step 9 now waits ONLY on C10b's 10-04 read.
- 10-02 17:05Z: **652 (October evidence audit, [7772dc], c1301c3c): no DECIDED switch rests on a VOID read.** Decisions on its flags:
  - **F1, DROP the Arc hook** (hook_estimator meanz_ne; MZ / CPMZ). Its held-out read is circular, its A/B/C value is fitted to
    458d's read, P12's floor is −0.003, and 629 finds Arc's programme uninformative off-line. The frozen 10-07 rule drops MZ / CPMZ
    from 537's tree, and at the 10-05 decide step hook_estimator is set to its non-Arc value, citing 652 F1 and 629.
  - **F2, CPH (10-07 s1):** the read stands. A CPH WIN enters October only if 614's K562 price agrees in sign (if 614 passes its
    trust gate); otherwise it needs a written reason that names the call-set mechanism.
  - **F3, K16 / K20 (rc_callset_cpm),** written BEFORE the 10-06 K16 read: a board-only switch (H1 flips sign across regimes),
    graded with the 610 instrument card. If adopted, October reverts to the default K (2,200) whenever D/E/F's depth (641's row) or
    predicted n_real falls outside A/B/C's range.
  - **F4, tercilec:** 10-06 step 9 is DEFERRED (APPLY DAY v2 step 9 amended).
  - **F5 / F6 / F7** noted:
    - lodo_strength needs C10b; lodo_hi is feedback-only.
    - ρ* / α are A/B/C-calibrated; one check that α does not reach DE magnitudes.
    - 642's Orion HEK293T rows on a Flex line are unpriced (contingency only).
  - **F8:** the runbook's PENDING 612 is marked RETIRED (idall3 replaces KGID).


### 10-02 17:14Z — the 10-07 / 10-08 rule, FIXED before any 10-06 entry exists (653's v4 of 537's draft)
- **Rule file:** `experiments/653-rule1007-review/rule1007_v4_653.py`, sha256 `e00af6c7d126d4768dfc9f8c5d1541f44d80d98bba40f10fefe0383fd18128ce`. v3 (537's
  `rule1007_draft537.py`, 8367c89d) is superseded. v4 = v3 plus R1 (the 10-08 list excludes the arms actually scored,
  `--scored=`), R2 (a missing capture reads as UNKNOWN), R3 (LEAK bin if pds/mse move > 1e-6 raw vs the base), R4, R8 and R11
  (anchor_version equality). Selftest PASS. With `--arc-hook` and `--scored` = the actual picks, it equals v3 in 384/384 states
  (equiv653, re-run by the coordinator).
- **Invocation: WITHOUT `--arc-hook`** (652 F1: the Arc hook is dropped, so MZ / CPMZ never enter a list). On 10-08, always pass
  `--scored=<the tags actually scored on 10-07>`.
- **R4 amendment to CPH's SHIP bin, fixed now, before any CPH read:** SHIP needs Δ(fid + jac) ≥ +0.0012 vs 446bj2z (535) AND
  overall Δ score_avg vs 446bj2z > 0. Even then, October adoption needs 614's K562 sign or a written reason (652 F2).
- **R6 hand step:** read 524dr (10-06 s1) with 524's reader BEFORE the 10-07 waiter fires (10-06 ~20:00Z). On LAG or a falsifier,
  10-07 s1 becomes the diagnosis arm or is held, and the rule's 10-07 s1 is not submitted.
- **Expected 10-07 picks**, as 653 verified for every reachable input: s1 = CPH; s2 = K20 if k16 ≥ MID, else 459dr. All are READY
  packages.
- **10-08 and after (R7):** the rule's 10-08 list, then the post-tree queue, each taken only when READY and not yet scored:
  630 → 524dr2 → factorial validation (ALLSW) → idall3 → 625. ALLSW and S1 are not built yet: a 10-08 s2 with nothing READY takes
  the next READY queue arm. Capture k16 immediately after it publishes (R9). 459dr is a pre-LODO package, so a ≥ +0.0015 H1-dose
  result must be confirmed inside 524dr2 (R10).
- **Chain:** chain1007 (to write by 10-05) runs the rule at 10-06 ~20:00Z on the captured 10-01 … 10-06 entries and arms 10-07's
  two submits through `scripts/arm_submit_after.sh` (FLOOR 202610070000; names claimed from MODEL-NAMES; `-d` empty).
- 10-02 17:19Z: **631's DOWN pass, H1 delivered level** ([d4cd5f], c5d6eefb; predictions registered first; 3 seeds).
  - DOWN alone delivers 0.14 of the ×2 need, overall-eq +0.0002 / +0.0008 / +0.0018 (f 1 / f 0.35 / cmp).
  - With 630 v3 it delivers 0.28 (630's sub-gate flow refills the ballast pool): +0.0023 / +0.0037 / +0.0049, vs 630 alone
    +0.0023 / +0.0029 / +0.0034. So DN adds +0.0008 … +0.0015 at f 0.35 / cmp and ~0 at f 1.
  - fid and jac are exactly 0; reach moves ≤ 0.0014 raw.
  - **Decision:** an emitter-level fix only if v7 (paired payments) lifts delivery well beyond 0.28. Otherwise DN becomes a 10-08+
    rider on 630, read against 630's entry after 614's cancer-line read and its own prereg. In the queue it goes after idall3 and
    before 625.
- 10-02 17:47Z: **614 trust gate FAILS** (4db74d02): zs PASS (1.9×), coarse MDE FAIL (1/22), cp FAIL (1/13). The truth's power is
  the cause; the plumbing is fine.
  - **Coordinator decision on the preregistered faithful-key re-run:** conditional WAIVER. Run it only if the H1@11k control
    reproduces the board's MDE size (≥ 1/3 of −0.0168), which would point at the K562 line or key rather than depth or power.
    Otherwise waive it, because a better key cannot add gated genes. Either way, record it in 614's TRUSTGATE.
  - **Consequences (no K562 held-out read exists for delivery or call-set levers):**
    - CPH's October adoption (652 F2) needs a WRITTEN REASON;
    - DN's rider case rests on H1 + its board read + a written reason;
    - 630's October dose stays K 1.0;
    - 614's hook-form rows are descriptive.

- 10-02 17:51Z: **615 assigned ([21dd4c]): the 614 harness on Orion RAW HEK293T** (200 cells/target, 20.7k UMI, 219k controls). It is the
  only public cancer-like surface with ≥ 100 cells per target AND ~A/B/C depth; K562 GWPS, K562 essential, Nadig, RPE1, Frangieh
  and Myllia all miss on cells or depth (614 RESULTS). Order:
  - (0) A preregistered technical-blob test: batch-matched vs pooled controls. Stop if the blob stays near the thinned 26–44 %.
  - (1) Key from K562 + Nadig donors only.
  - (2) 614's gate with bands committed first, plus the low-power diagnostic up front.
  - (3) HCT116 raw only if HEK passes.
  Heavy stages run outside 10-04 05:00Z – 10-05 18:00Z. If it passes, it becomes the held-out read for CPH, DN, K and 630's dose.
- 10-02 18:13Z: **631 v7 (paired payments) is NOT adopted; there is no emitter-level work** (DOWN delivery not lifted well beyond 0.28).
  - Overall-eq at f 1 / 0.5 / 0.35 / cmp, 3 seeds, registered first:
    - v7: +0.0010 / +0.0038 / +0.0043 / +0.0043;
    - 630 v3: +0.0023 / +0.0030 / +0.0029 / +0.0034;
    - UPDN (630 + DN): +0.0023 / +0.0030 / +0.0037 / +0.0049.
  - v7 wins at f 0.35 but loses at f 1, where its payments overshoot weak called-DOWN genes.
  - **The 10-08+ rider is UPDN:** 630 frozen plus a separate `VCC_DE_DNMAG` knob (module src_dn_mag631). It is ≥ v3 in every
    regime and simpler. It switches to v7+DN only if that beats UPDN at f 0.35 without the f 1 loss (its run is registered and in
    progress).
  - The read is against 630's own entry, with a raw-unit prereg before any build.
- 10-02 18:38Z: **616 approved ([21dd4c]): re-check every council instrument built on the THINNED Orion harnesses on CLEAN raw folds.**
  615 traced the Orion blob to the thinned builder, so its DOWN bias could have favoured DOWN-leaning priors. Instruments affected:
  480's gates, 636, 613 / 613c / 613d, 633's Orion d20 rows and 609.
  - Every row reports Δ(clean − thinned).
  - A clean-fold version of 613d's 496a prediction is registered BEFORE 10-04 05:00Z, next to the thinned one already registered.
  - 613's 458d retrodiction is re-checked on the clean folds.
  - Any sign flip on an open call (idall3's held-out read, KGID's regret table, 613c, 630's dose) goes to its owner and to the
    coordinator.
- 10-02 18:50Z: **DNMAG rider (631) base and day fixed by the coordinator.**
  - Base = 630's scored entry (`VCC_DE_UPRN=b_k100_q20r3p10` on the best base at 630's arming). The arm = that base +
    `VCC_DE_DNMAG`; its config diff must be exactly `de_dnmag`.
  - Day = the first free slot after 630's read and idall3's, per the queue (before 625).
  - Apply: `patch_builder631.py` joins the 10-08+ apply list after 604f's slot. It needs a 606-style composition re-run on the
    applied 10-06 tree, plus its KA0 (unset byte-identical).
  - Prereg: 631/PREREGISTRATION-DNMAG-DRAFT.md (70eb605a). Bins on Δnmae raw vs 630's entry: ≤ −0.0025 TRANSFERS, to −0.0008
    PARTIAL, NULL, > +0.0008 HARM.
  - H1 definitive: UPDN − v3 = +0.0008 / +0.0015 overall-eq at f 0.35 / cmp, ~0 at f 1. With A/B/C capacity ~½ H1's, PARTIAL is the
    likely board bin.
- 10-02 18:50Z: **DNMAG prereg FROZEN:** `experiments/631-nmae-headroom/PREREGISTRATION-DNMAG.md` (bb9b9189), sha256
  2019e24e9941e5a0f334d68c8be7507a3fce9c8cf297d83d528a99e891197336, re-verified by the coordinator. The expected bin, written
  before any read, is PARTIAL (centre ≈ −0.002 raw). Interaction rule with 654: if DN adds < ½ its standalone gain on a ×1.25
  base, the predicted range halves before the read; the bins stay unchanged.
- 10-02 19:04Z: **615 gate STOPPED by the coordinator before any gate row.** Its preregistered pre-gate power diagnostic reads worse
  than K562, whose gate failed on power: median 1 DE gene per target and 18 of 100 targets nmae-scored, against K562's 5 and 34.
  The wrapper KA is kept, since it validates run615.
  - **Bigger finding:** the THINNED Orion HEK gate table had a median of 141.5 calls per target, against 1 on raw cells with more
    cells and more depth. The thinned Orion DE truths were almost entirely the builder's selection artifact.
  - **Pending 616's clean gate tables (≤ 1 h):** if they are near-empty, Orion fold rows carry no usable DE truth. [21dd4c] then
    lists every instrument reading and decision that used them (613's 458d retrodiction, 613d's cancer folds behind idall3's
    held-out read, 633's dose rows, 609 / 639d), with what remains for each. The coordinator re-rules the open decisions from
    that list.
- 10-02 19:10Z: **Coordinator re-ruling after 616's Orion-dependence audit** (`experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`,
  7dac7a98; nothing recomputed).
  - **Untouched:** every BULK decision (q, lodo, lodo_hi_pct, sharpen_pt, cd4_stack), production's Orion donors and PLOLO's Orion
    folds (the donor tables pool all controls, 615 config P, clean), h1-test, K562 gates, Broad lines, de_mag and K16/K20.
  - **613 block instrument** is downgraded to "consistent on K562 only": 2 of its 3 retrodiction folds were thinned Orion.
    ID05/ID10 stay empty-slot only. Their held-out reading is K562 (null block +0.0001) and Broad-16 (block −0.0022, 7/12 negative).
  - **mde_sign / 496a:** decided by 496's own bins, unchanged. The 613d-vs-H1 discrimination rule ("496a < −0.0005 closes 613d")
    is WITHDRAWN: 613d's registered 496a number averages two void Orion rows with K562 (+0.0013). The number stays on record, but no
    conclusion about 613d is drawn from it. 613d rests on its two out-of-sample H1 passes.
  - **idall3** keeps queue slot 4 and its board bins. Its held-out case is H1 (registered PASS), Broad-12 neutral (10–15 % K562-like
    tail) and K562 against. A board WIN enters October only with a written reason; the strongest one available is that A/B/C's
    read is an out-of-sample test on three Arc Flex lines.
  - **programme_source = off** stays (now K562-only evidence; off is the conservative side). **604 v3's identity-only exclusion**
    stays (conservative; its HEK evidence is void). **639** stays positive on H1 +0.016 and Broad-16 +0.009. **630** is unchanged
    (H1 + K562).
  - **Pre-council DE conclusions** from thinned harnesses (~30 experiments, 031–194 and 273) are void for the DE members.
    [21dd4c] checks which FINDINGS §2 closures rest ONLY on them; those become re-openable.

- 10-02 19:16Z: **618 approved, faithful key ([21dd4c]): a SECOND H1 real-emitter surface on unworn targets** (h1-trainval, 126 targets ×
  400 cells, same chemistry). h1-test is now the only valid real-emitter DE surface and the most tuned-on (11.1c.5).
  - Conditions:
    - a KA reproducing 323's cap15 on h1-test's targets bit for bit before trainval's key is built;
    - shared controls are declared;
    - call-set / lcap rows are marked WORN (528's β was selected on trainval);
    - 614's trust gate, with bands committed first;
    - a fresh-draw re-check of 613d's registered H1 passes (idkg3, idall3);
    - heavy stages outside 10-04 05:00Z – 10-05 18:00Z.
  - It buys a panel-draw noise estimate for every H1 price and an unworn surface for October confirmations.
- 10-02 19:49Z: **idall3's K562 line is INCONCLUSIVE, not "against."** On 614's real K562 emitter (ec622f62; descriptive, since 614's gate
  failed on truth power), idall3 − c4h = −0.00026 nmae (seeds −0.0002 / −0.0002 / −0.0004), so idall3 ≈ no hook; c4h − no hook =
  +0.00027. 613d's analytic K562 row said idall3 was worse by +0.0005 … +0.0024. Both read a low-power truth (median 5 DE genes per
  target), and the hook acts on weak genes, the lever class 614's gate failed on.
  - idall3's October written reason states both rows and calls K562 inconclusive.
  - Its held-out case is therefore H1 (registered PASS) + Broad-12 (neutral) + K562 (inconclusive).
- 10-02 20:03Z: **chain1007 ARMED** (pid 39901; `experiments/660-slot-1007/`).
  - It waits for 10-06 20:00Z, checks the rule sha and runs the frozen v4 rule without --arc-hook.
  - **R6:** `work/HOLD-1007-s1` holds slot 1 (write it if 524dr reads LAG or hits a falsifier).
  - Dry run today: s1 446bj2zcph, s2 459dr (k16 unread).
- 10-02 20:16Z: **Due tonight (heartbeat):**
  - 656's verdict on the same-line w re-price with zero repair (registered 28d99181; ~22:30Z+; `work/read656.out`). Rule: the October
    w moves only if some w beats 0.5 at BOTH p by > 2 per-seed sd. If [7772dc] has not committed it by 10-03 02:00Z, the coordinator
    reads and commits it.
  - 655's carrier-amplitude dose curve (×1.5 / ×2 with zero repair).
  - 618's KA and ref_tv (waiting on memory).
- 10-02 20:18Z: **614's conditional waiver did NOT apply. The faithful-key re-run runs.** The H1@11k control gives G2 coarse MDE −0.00648
  (≤ the preset −0.0056 trigger; H1@20k −0.00966; G1 zs −0.00337 vs −0.00183 at 20k), so depth alone does not kill the lever. The
  rule was set before the number; no post-hoc waiver.
  - Scope (written before the run): it tests only the KEY. Cells per target and line biology remain as explanations.
  - Priority: behind 618, 655 and 656, with no deadline. FAIL closes 614; PASS gives a validated K562 harness.
- 10-02 20:27Z: **DNMAG gates pass** ([d4cd5f]): the UPDN KA1 is bit for bit on 85/85 targets × 3 seeds. A/B/C capacity after 630 on
  446bj2z's real blocks: dn_frac A 0.196 / B 0.220 / C 0.185, ~1,500 genes per target, fail-safe 0, floor met everywhere. v7+DN
  overshoots (f 1 −0.0008), so UPDN stays the rider.
  - **630 float-noise note:** carriers are promoted for genes already at target (dA ≈ +1e-16). This is inside every priced and KA'd
    630 number, so 10-06 is unaffected. Coordinator rule: NO 630 module revision before 630's board read (keep that read clean); a
    dA > 1e-9·A guard is reconsidered afterwards, only if its re-price moves a number.
- 10-02 21:09Z: **655 carrier-amplitude dose curve (H1 real emitter, zero repair; d7897f74).** Overall-eq vs ×1 at f 0.5 / f 0.35 /
  cmp:
  - ×1.25: +0.0011 / +0.0019 / +0.0012 (+26 calls per target);
  - ×1.5: +0.0007 / +0.0031 / +0.0019 (+57);
  - ×2.0: −0.0022 / +0.0005 / −0.0004, and infeasible on one draw (the emitter's per-cell totals assert fired).
  Preregistered bin: s* = ×1.5, WEAK CANDIDATE (min over regimes +0.0007).
  - **Coordinator:** carry ×1.25 AND ×1.5 to 618's held-out read with zero repair, and pick by the regime-robust max–min rule.
    ×1.25 is the default lean: it has half the call-set change, and ×1.5's board price leans on 611's provisional call-set term.
  - Before any arm: a builder KA for the α path, and def_postbuild's alpha pin (0.982, def_postbuild.py:116) moved for any October
    value (10-08+ apply list).
  - Queue: the α arm goes AFTER ALLSW and before idall3: 630 → 524dr2 → ALLSW → α → idall3 → DN → 625.
  - Tries: 5 amplitudes on H1, 0 on the board.
- 10-02 21:11Z: **α arm plan (655 RESULTS, f15380e2; unapplied).**
  - **The A/B/C arm needs NO code change:** `VCC_ALPHA=1.2275` (×1.25), EXPECT config diff exactly `alpha`, plus
    `--allow-arm-drift VCC_ALPHA` if def_postbuild runs. The 0.982 pin stays until October adopts.
  - **KA before the arm:**
    - (a) staging: DE rows at 1.2275 = 1.25 × rows at 0.982 (≤ 1e-12 relative), with the loaded bulk IDENTICAL;
    - (b) a full build with the rc audit: config diff exactly {alpha}, block-total composition pinned;
    - (c) 554's dose guard PASSES in every context. A/B/C's capacity is ~0.70× H1's, so this is the check most likely to bite.
  - **October adoption is one commit:** def_postbuild SHIPPED_ARMS["alpha"], def_pipeline SHIPPED_ENV["VCC_ALPHA"], a fresh
    --incumbent meta built with the new α, and ARM_BOARD_DELTA's note. On D/E/F a failed dose guard falls back to 0.982 by written
    rule.
- 10-02 21:35Z: **DNMAG and α ×1.25 are SUBSTITUTES** ([d4cd5f], prediction registered first). 630 on ×1.25 adds (+0.0033 at
  f 0.35), but DN's increment on ×1.25 is −0.0005 / −0.0003 / +0.0009 (f 1 / 0.35 / cmp): it overshoots. The frozen interaction
  rule triggers, so DNMAG is never read or adopted on a ×1.25 base.
  - **Coordinator decision: path (a), α ×1.25 + 630.** It gives +0.0052 (f 0.35) / +0.0049 (cmp) over the plain base and is larger
    in the board-like regimes, while (b) ×1 + 630 + DN gives +0.0037 / +0.0049. It is ONE setting with no code change, and its
    call-set term looks benign on H1 (fid / jac ~0; reach untrusted per 611).
  - **Conditions:** 618's held-out read confirms ×1.25, and the α-path KA passes, INCLUDING 554's dose guard on A/B/C.
  - **Fallback (b):** DN on a ×1 base is read ONLY if α fails those conditions or reads NULL / HARM on the board. Only one of the two
    is ever read.
  - **Queue:** 630 → 524dr2 → ALLSW → α (×1.25) → idall3 → 625. DN is fallback only.
- 10-02 21:35Z: **Amendment to the α / DN pick, made before any 618 read:** both stacks are priced on 618's fresh-H1 surface
  (zero repair, 3 seeds, predictions registered first), and the pick follows ONE instrument:
  - (a) ×1.25 + 630, vs (b) ×1 + 630 + DN;
  - rule: the higher MIN over {f 0.35, f 0.5} wins;
  - within 0.0005 → (b), which carries no capacity risk;
  - (a) still needs 554's dose guard on A/B/C.
  Why: at ×1.25 the call-set label is nominal (fid −0.00002, jac −0.00006 raw), so its gain is nmae, priced like delivery. Its real
  risks are carrier capacity (the dose guard) and unpriced reach at +26 calls per target; DN has no capacity risk. The queue slot
  stays one arm: whichever stack wins.
- 10-02 23:20Z: **656 verdict: same-line w = 0.5 STANDS** (zero repair, 3 seeds; `work/read656.out`). In board units vs 435's base:
  - p 0.2: w 1 −0.00071 (sd 0.00016), w 2 −0.00486;
  - p 0.3: w 1 −0.00799, w 2 −0.01386;
  - zero non-own zero totals on every arm.
  450's choice holds without the artifact, and the October w is unchanged. [7772dc] commits 656's outputs; the coordinator does it if
  they are still uncommitted at 10-03 02:00Z.
- 10-02 23:34Z: **618 TRUST GATE PASS** on all three pairs, delivery-only too (e9f0342f).
  - G1 zs −0.00163 (board / 618 1.34); G2 coarse MDE −0.01473 (1.14); G3 cp fid −0.00408 vs the board's −0.00358.
  - Zero repair leaves 0 zero totals.
  - 618 is now the held-out instrument for the α / DN stack pick, the 613d re-check and the hook forms (all running in that order).

- 10-03 10:09Z: **Overnight incident.** The Mac ran on battery with the lid closed and slept ~01:45Z–09:59Z 10-03.
  - 10-03 s1 (495) finished on wake through the 536 resume and is a NEW BEST at 0.2739107; lodo_hi_pct is DECIDED. s2 (471c) has
    been uploading since 10:00:45Z (cap 12:00Z).
  - Battery was at 82 % and discharging at 10:09Z; Ahmet was pushed twice to plug in.
  - **Worker 537 is DEAD** (pid 75622 gone; log silent since ~01:30Z; no result event). It is NOT relaunched: its ranked list (v5),
    the frozen 10-07 rule and the council's pricing cover its remaining scope. The coordinator takes over its last open item,
    537 §2.5's preregistered noise check on the 10-01 … 10-04 reads and the rule 13 rewrite, after the 10-04 reads.
- 10-03 10:54Z: **10-03 s2 (471c) upload PAUSED by the coordinator at 10:53Z** (entry yzNoMVih0pKqoDjdqrXm is in 10-03's allowance; the
  pending upload is kept). The Mac is on an iPhone hotspot (gateway 172.20.10.1) at ~0.28 MB/s with ~1.6 GB left. The pause was
  made on a stale battery reading: AC was back at 10:48Z. It now stands only to protect Ahmet's cellular data.
  - **To finish it:** `bash experiments/495-lodo-cap6/resume1003s2.sh`. It refuses on battery or on the hotspot unless FORCE=1; on
    rc 0 it records and captures the entry. Then read it with `read1003.py none <id>`.
  - **Data warning for the next nights:** each day's two uploads are ~5.5 GB. Ahmet decides whether they run over the hotspot.
- 10-03 12:07Z: **618 hook rows (held-out H1, faithful key; 0cd71e16).** Real emitter, raw nmae at f 1 / 0.5 / 0.35:
  - c4h − no hook: −0.0023 / −0.0007 / −0.0014;
  - **idall3 − c4h: −0.0024 / −0.0006 / −0.0007**, about +0.0002 overall at f 0.35 under 619's ×1.0–1.2. This is a held-out H1
    confirmation of idall3's sign in every regime.
  - idkg3 − c4h: −0.0012 / −0.0008 / −0.0007.
  - The registered 613d re-check FAILS for idall3 at f 1 / f 0.5 (613d is conservative on a new panel), so by its branch hook prices
    come from the real emitter on each surface. c4h and idkg3 pass.
  - Unpriced reach differs between the two identity forms: idall3 −0.0048 at f 0.35, idkg3 +0.0040 at f 1. Since 619, reach is not
    a price on this surface, so idall3's 10-05+ board read settles it, and its written October reason cites this row.

- 10-03 12:09Z: **Ahmet: hotspot data is fine.** 471c's upload was RESUMED at 12:08:34Z (`FORCE=1 resume1003s2.sh`) and should
  finish around ~13:45Z at ~0.28 MB/s. Read it with `read1003.py none <yzNoMVih0pKqoDjdqrXm>` once CAPTURED. No future upload is held
  for cellular data. AC + an open lid remain the rule.
- 10-03 17:05Z: **10-04 slot 2 FILL (coordinator).** Preview of the frozen rule1004 on the captured inputs (453dr q3VuTI…, 479p2g ioNMsG…,
  479sp2 NPiJoP…, 495 9Hxef2…): `498p2g 460dr|460dr|P=P2g,C=TRANSFERS,O=OK`.
  - s1 = 498p2g is READY: gates PASS, CLIFF OK, novelty OK, ρ A/B/C 0.97930 / 0.97930 / 0.98012.
  - The s2 fill list is exhausted: its only entry, 460dr, was scored on 10-01, and 497dr needs P = S_P2. So chain1004 arms s1 alone.
  - **Fill: 446bj2zk20 (K20, READY; gates PASS, CLIFF OK, novelty OK) under its own prereg** (523's bins + the 610 instrument card,
    a995f0e2). It moves K20 up from 10-07 s2, where the v4 rule then takes 459dr (the k20 id is passed). No rule's decision changes.
  - **Due at the 20:17Z heartbeat**, after chain1004.log shows "armed 10-04: s1 498p2g":
    `nohup bash experiments/498-combined-bulk/arm1004.sh 2 446bj2zk20 experiments/428-zero-gene-audit/work/build_446bj2zk20/run446bj2zk20.log submissions/446bj2zk20-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cpm1600.vcc "You're behind—very behind." "446bj2zk20 slot 2 on 10-04 (coordinator fill)" > experiments/498-combined-bulk/work/arm1004_s2.nohup 2>&1 &`
- 10-03 17:04Z: **α / DN pick (618, mechanical, 08ebeed1): (b) 630 + DN WINS.** MIN over {f 0.5, f 0.35}, overall-eq vs plain:
  (a) ×1.25 + 630 +0.00360 vs (b) 630 + DN +0.00404. The difference −0.00044 is inside the 0.0005 band, so (b) wins either way; the
  nmae-only reading agrees.
  - H1's registered "(a) by +0.0014" did NOT transfer: ×1.25 alone is only +0.0003 … +0.0009 on 618.
  - **×1.25 is NOT bulk-neutral in 618's harness** (block totals move in ~1.35M cells), so it was priced without price630's
    sums-pinned assert, and its pds / mse effect is outside these numbers.
  - **Queue:** the α arm is DROPPED. The magnitude arm after 630 is DNMAG on 630's entry (prereg bb9b9189, frozen):
    630 → 524dr2 → ALLSW → DNMAG → idall3 → 625.
  - 654 / 655's α results stay recorded with this caveat.
- 10-03 19:03Z: **661 board predictions registered from 618** (00ad7165; 619's multipliers; reach unpriced). Each later read
  also tests 618's transfer.
  - **630:** +0.0030 overall (+0.0027 … +0.0032), nmae only, rank-neutral.
  - **DNMAG on 630's entry:** +0.0015 (+0.0013 … +0.0016). Note 618 puts it in the TRANSFERS bin, where the prereg expected PARTIAL;
    it is regime-dependent (+0.0009 if the board reads like f 0.5).
  - **idall3:** 618 predicts a **TIE** vs Δ_496a (+0.00007, band −0.00002 … +0.00015), where h1-test predicted a WIN. Its read now
    mainly discriminates the two surfaces; TIE → the plain hook stays.
  - **625:** octile scan on 618 running (~22:00Z). Predicted: {0, 1, 2, 6, 7} degrade; only {5} (maybe 3 / 4) can pass 625's rule.
- 10-03 20:20Z: **10-04 ARMED:** s1 498p2g (pid 43945, chain1004, 20:01Z); s2 446bj2zk20 (pid 43983), the coordinator's fill
  through fill1004s2.sh (20:01:17Z). Both submit at 00:01Z. Over the hotspot (~0.3–0.5 MB/s) each 2.77 GB upload takes ~1.5–2.5 h,
  so expect s1 ~02:00Z and s2 ~04:00–05:00Z (cap 12:00Z). Readers: 498's `read1004` path in rule1004 / 498 PREREG for 498p2g
  (vs 495 = the per-target increment on the cap-6 base); 523's bins + the 610 card for K20 (vs 446bj2z).
- 10-03 21:38Z: **625 octiles on 618 (661 item 4, 16570222).** All 8 are measurable, but **none passes 625's rule at f 0.35.**
  - Octiles 0–3 and 6–7 degrade nmae.
  - Octiles 4 and 5 lower nmae but fail fid no-harm; they "pass" at f 1 only through call-set expansion (+156 … +1,238 calls per
    target for octiles 4–7). Octiles 0–3 move calls by ≤ 5.
  - **Coordinator decision:** 625 shrinks to **5 slots: octiles 0–3 + the combination.** Its passing rule gains a call-set guard
    (|Δcalls| ≤ 10 per target in the board rc diagnostics), so a D/E/F mask never carries call expansion. It stays LAST and
    empty-slot only. [7772dc] amends 625's prereg before any read.
- 10-03 21:42Z: **625's A/B/C readout is DROPPED** (owner's Amendment 4 answer, ff3786b1; octflip625b.py).
  - Flipping octiles 0–3 helps on 0 / 12 Broad-16 lines (agreement 0.58–0.85, H1 agrees), so an A/B/C read is expected NULL. The 5
    slots are freed.
  - **The transferable signal is octile 4,** CLEAN4's weakest DOWN prior: wrong on 12 / 12 Broad lines (agreement 0.40), a coin flip
    on H1. Its block-prior flip expands calls, though.
  - **New item 657 ([7772dc], with [d4cd5f]):** deliver octile 4's UP sign CALL-NEUTRALLY, as the hook sign on uncalled genes and/or
    UP signs inside 630's rank-neutral pass. Price it on 618 against the current hook / 630 base with a registered prediction and
    |Δcalls| ≤ 10 per target. Class: transfer (12 independent KO lines + neutral H1). An A/B/C arm only after a 618 PASS.
  - **Queue:** 630 → 524dr2 → ALLSW → DNMAG → idall3 → (657 if it passes 618). 625 is removed.
- 10-04 05:25Z: **498p2g's reading rule, FIXED NOW before its entry is scored** (it had no preregistered bins; C10b's F depends on it).
  Increment Δ = 498p2g − 495 (entry 9Hxef2F6Clf4G2hIkNPq; the per-target sharpening on the cap-6 base):
  - **P2g HOLDS** (F = M6+P2g, sharpen_pt stays DECIDED) unless both hold: Δ overall ≤ −0.0005 AND Δ bulk (pds + mse scaled,
    overall units) < 0. The DE members carry only the re-roll (sd 0.00025, rule 13), so the bulk decides.
  - **P2g FAILS** (F = M6; sharpen_pt is re-opened and October runs it off pending a re-decide) only in that case.
  - A decided switch is not reversed on noise. 479p2g's own increment on cap 3 was +0.0010 overall (pds +0.0016, mse −0.0037 raw).
- 10-04 06:06Z: **Build plan for the 10-08+ arms (the coordinator owns it; worker 537 is dead).** All builds use the POST-10-06 builder
  (production path), are serial under the memory gate, and stay clear of C10b (10-04) and R7 (10-07 06:00–12:00Z).
  - **10-07, after the 10-06 apply + KAs:**
    1. **630 arm** = 446bj2z's recipe + `VCC_DE_UPRN=b_k100_q20r3p10`. 446bj2z is the scored zs-carrying base; config diff exactly
       `de_uprn`; 630dd audit FAILSAFE 0 + capacity floor.
    2. **ALLSW** = 498p2g's env (chain498.sh) + `VCC_DE_MAG=dn_k250_zs` + the DE switches decided by then (mde_sign per 496a on
       10-05; tercilec per step 9). REF_META 498p2g's. It doubles as 637's factorial validation (bulk X = cap6+P2g, read; DE Y =
       zs, read).
    3. **524dr2** = def_day on efr / efp with the resolved post-10-06 recipe (the October dress rehearsal on the applied tree).
  - **10-08 → 10-09:** DNMAG (after 630's read, patch_builder631 applied with its KA0), then idall3 (base = the best entry at slot
    time).
  - **chain1008** (10-08 slots) is written by 10-06 evening: v4's 10-08 list with `--scored`, then the READY-gated queue tail.
- 10-04 06:18Z: **657 CLOSED** (every route is slightly negative on 618; the octile-4 Broad signal looks like a KO-chemistry offset). It
  gets no slot, not even an empty one, since its expected sign is negative. The queue stays 630 → 524dr2 → ALLSW → DNMAG → idall3.
  **658 assigned to [7772dc]:** preregistrations for the two 10-07 builds, due 10-06 18:00Z:
  - (a) ALLSW (498p2g env + zs + the DE switches decided by then), with 637's `--validate` additivity reading and raw-unit bins;
  - (b) 524dr2 (the post-10-06 October dress rehearsal), with 524's LAG / falsifier bins carried over and amended for the new
    switches.

- 10-04 06:44Z: **658 preregs accepted** ([7772dc], cdc1445c).
  - **ALLSW** = 498p2g's round-2 build + zs (+ the clean4:0.05 hook only if 496a TRANSFERS). Prediction 0.2754683 (+0.00055 vs 498p2g,
    nmae only).
    - **Deviation ACCEPTED: tercilec stays OUT of ALLSW.** It re-tiers only 4 A/B/C targets, and including it would re-roll the DE
      and void the exact interaction read; 524dr2 carries it.
    - Reads: F1–F3; R1 = 637 `--validate`; R2 = the exact interaction (ADDITIVE ≤ 0.00025). Pre-slot gate: `bulkid658.py` must print
      BULK IDENTICAL vs 498p2g.
  - **524dr2:** 524's bins vs 453dr, plus a primary read vs 524dr by branch N0–N3 (N0 = nothing new → no slot) and a culprit map.
  - **R6 hand step on 10-06 (the 524dr read):**
    `uv run python experiments/658-prereg-1007-builds/read658.py dr <524dr id> --recipe experiments/524-october-recipe-dr/work/recipe524.json [--hook <496a id>]`.
    Exit 0 = read, 2 = refusal / hidden members, 3 = falsifier. On LAG or exit 3, write `experiments/660-slot-1007/work/HOLD-1007-s1`.
  - 637's leak band is off-centre against 9 board pairs (mean −0.00017, sd 0.00021). [d4cd5f] is asked to re-centre it to
    −0.0002 ± 0.0004 before Oct 22.
- 10-04 06:46Z: 637's leak band is re-centred and versioned (0b484be2; readdef_factorial637 01eb7374; test 32/32): v2 is the default (centre
  −0.0002); v1 stays reproducible with `--leak-band v1`. ALLSW's R1 moves to v2 by a pre-read amendment ([7772dc]). [d4cd5f] writes
  build_630arm.sh by 10-06 18:00Z.
- 10-04 06:56Z: **10-07 build drivers:**
  - `experiments/630-rank-neutral-up/build_630arm.sh` ([d4cd5f], committed, UNARMED; DRY=1 → NOT READY until the 10-06 apply, as
    expected). Launch after APPLY DAY v2 step 8; it needs step 7's post-598 driver sha in run_arm_generic.sh. Terminal line:
    `READY 630arm <sha>` / `NOT READY <reason>`.
  - ALLSW and 524dr2 drivers ([7772dc], 658's dir) are in progress; ALLSW's R1 is amended to the v2 band (9e3b0c70).
- 10-04 07:02Z: **All three 10-07 build drivers committed, unarmed** (e01f5662; test658 39/39). Launch order on 10-07, after APPLY DAY step 8:
  1. `bash experiments/630-rank-neutral-up/build_630arm.sh`;
  2. `ALLSW_BUILDER_SHA=<post-10-06 builder> ALLSW_EMITTER_SHA=<emitter> bash experiments/658-prereg-1007-builds/build_allsw658.sh`;
  3. `DR2_BUILDER_SHA=… DR2_EMITTER_SHA=… bash experiments/658-prereg-1007-builds/recipe_dr524dr2.sh`.

  Serial, one at a time. Each ends with one `READY <tag> <sha>` / `NOT READY <reason>` line. 524dr2 classifies its recipe vs 524dr's
  (N0 → no slot). chain1008 waits on these lines.
- 10-04 09:10Z: **chain1008 + arm1008 drafted** ([7772dc], 66cf2c5a; test659 13/13; UNARMED until the coordinator arms it by 10-07 evening).
  **Order decision, before any 10-07 entry exists: QUEUEFIRST-1008 set.** 10-08 takes the rule's DE follow-ups first (CPHK20 / K20 …),
  then the post-tree queue (630arm → 524dr2 → allsw658), then the rule's generic fallbacks (459dr, ALLSW).
  - 630 reads early because DNMAG's base is 630's entry.
  - 524dr2 reads early so a defect in the applied October recipe leaves time to fix.
  - ALLSW still reads well before Oct 22.
  This amends my 10-02 freeze note ("the rule's 10-08 list, then the post-tree queue") for the generic fallbacks only. The rule's
  own decisions are untouched.
- 10-04 13:58Z: **10-04 s2 (K20) is crawling.** It sits at 82.2 % committed (2.29 / 2.79 GB). On the iPhone hotspot the session's offset query keeps
  failing ([Errno 8] inside httpx), and the link runs ~0.1 MB/s with an RTT up to 770 ms. resume1004s2 retries until 23:00Z.
  - **`vcc cancel` exists:** cancelling a stuck submission does NOT count against the daily limit; only scored entries do.
  - **guard1005.sh is ARMED:** at 23:30Z, if K20 is still pending, it cancels the entry so 10-05's 00:01Z submit (496a, the mde_sign
    read) is not blocked by "one in flight". K20's package stays novel and can go on a later day.
- 10-04 17:51Z: **10-05 slot 2 fill = idall3** (663, 9421c53d). rule1005 gives `496a 498p2g 460dr|498p2g 460dr`, and both s2 candidates are
  scored, so s2 would be empty.
  - Build: idall3 = 498p2g (the best entry) + the ungated v3 identity hook, the prereg's arm.
  - It is gated on memory: free+inactive was 20.3 GiB against the 23 needed, held by apps and the compressor.
  - **cutoff_idall3.sh** kills it if it has not started by 22:30Z, so recipe_dr524 (02:00Z) is never delayed.
  - **fill1005s2.sh** arms arm1005 slot 2 with idall3 once chain1005 exits with s2 empty and the build prints READY (cap 10-05 04:00Z).
  - **Risk to watch:** recipe_dr524 needs free+inactive ≥ 25 GiB at 02:00Z. With the same apps open it will wait, and 10-06's s1
    (524dr) depends on it.
- 10-04 20:04Z: **10-05 ARMED: s1 496a** (chain1005, pid 35703, 20:00:39Z; rule `496a 498p2g 460dr|498p2g 460dr`; both s2 candidates
  "not ready" because they are scored). **s2 = idall3** through fill1005s2 (pid 91656), once its build prints READY.
  - idall3 built and passed its gates at 19:23Z (config diff exactly `de_mde_sign` + `de_sign_emitter_sha256`; CLIFF OK, equal to
    498p2g's). Its PREP gate has stayed closed since 19:23:58Z. The iOS Simulator holds ~17.5 GiB in the compressor:
    free+inactive 18.6–21.8 GiB against fi23's 23, and swap free 386 MiB against pressure55's 400. The gate expires at ~21:24Z.
  - **prep_fallback663.sh** (pid 34190; b66accb0) acts only on an expired prep gate. It takes over fill1005s2's role, waits for the
    same gate until 01:30Z, then preps on memory_pressure free ≥ 50 % as a logged deviation (25 s, 549 precedent). It gives up at
    01:45Z, so it never meets recipe_dr524's 02:00Z build. Ahmet was asked by push to quit the Simulator.
  - **10-05 reads:** `uv run python experiments/663-idall3-arm/read1005.py --a496 <496a id> [--idall3 <idall3 id>]`. It encodes 496's
    board bins, PRED-662-496a's instrument reading and CPH's T, and 612's IDALL3 bins against Δ_496a with 618's and I-662's
    predictions. KA: base-vs-itself reads TIE with 0 deltas.
  - **mde_sign decide waits for recipe_dr524's resolve.** 658's 524dr2 design assumes 524dr ran mde_sign OFF (PREREGISTRATION-524dr2
    line 23), and CHECKLIST v2 step 10 decides after dr524. Read 496a and write its BOARD-RESULT as soon as it is captured. Run
    `october_recipe.py decide mde_sign …` only after `recipe_dr524.log` has its `recipe:` line, or the CHECKLIST reaches step 10.
- 10-04 20:19Z: **idall3 READY, 10-05 s2 ARMED.** The prep gate opened on the normal fi23 branch at 20:17:51Z (free+inactive
  23.3 GiB, memory_pressure free 66 %), so no deviation was needed. Package `dad0e10eb85edd71a2ce8e7d539974651c096f1464abcf52e2254aa6ab4ee4af`
  (2,773,841,920 B; meta beside it 227dab65; novelty rc 0; prep peak footprint 24.1 GB). fill1005s2 armed arm1005 slot 2 at
  20:19:09Z (pid 41181); it waits for 496a's capture. prep_fallback663 read READY and exited ("nothing to do").
  - Both 10-05 slots: s1 496a (pid 35703) → s2 idall3 (pid 41181). Readers: `read1005.py --a496 <id> --idall3 <id>`.
  - recipe_dr524 (02:00Z) needs free+inactive ≥ 25 GiB. At 20:19Z it is 22.6 GiB, with the Simulator at 33.6 GiB RSS. It waits
    rather than fails (cap 10-05 20:00Z).
- 10-04 20:23Z: **10-05 CHECKLIST v2 preflight on tonight's live tree (read-only; nothing applied).** Step 1's five suites all pass:
  590 ALL PASS (min_top1 0.7716 / 0.7191 / 0.7329); 538 13/13; 554 5/5; identity554 "all identical: False" (correct; its json
  differed only in the runtime field and was restored from a byte copy); 641 4/4. Peak RSS ≤ 2.74 GiB. Step 2's
  `apply_october538.py --root . --check` reports "patches apply cleanly (not written)", def_pipeline 9b2c7eca → 9a718b9e as the
  checklist expects. Untracked `scripts/def_pipeline.py.orig` (09-27 08:36, sha 44714b9a, a pre-485 patch backup) is harmless;
  left in place. Tomorrow's run still starts at step 0, after recipe_dr524 exits.
- 10-04 20:25Z: **chain1008 reviewed and ARMED** (pid 48648; waits for 10-07 20:00Z; a wake after 10-08 06:00Z arms nothing).
  test659 13/13 on the live tree. A live DRY (`work/dry-chain1008.log`) now sees K20's capture (kJfe2dP9…). It runs the
  sha-pinned rule (`10-08: 459dr ALLSW` at today's inputs, H/K/CPH UNKNOWN until 10-05/06/07). With QUEUEFIRST-1008 present, its
  candidate order is 630arm → 524dr2 → allsw658 → 459dr. Today it would arm s1 459dr and s2 none, because no tail arm is built yet.
  The three 10-07 builds must print READY before 10-07 20:00Z to be picked. `tail1008.tsv` is read at that time, so it can still be
  edited. To change chain1008.sh itself, kill pid 48648 and re-arm it; never edit the live file.
- 10-05 02:02Z: **10-05 reads done** (BOARD-RESULT-496a / -idall3; eb329569, e3f0b16b). 496a TRANSFERS (+0.0005552 vs 446bj2;
  Δnmae −0.00158 raw) → mde_sign = clean4:0.05 at CHECKLIST step 10. idall3 BELOW the plain hook (−0.0002458 vs 498p2g) → no
  identity hook. recipe_dr524 released at 02:00:41Z; its recipe resolved with mde_sign PENDING → OFF, as 658 assumes.
  - chain1007 DRY with tonight's inputs: rule `CPH 459dr ALLSW|459dr ALLSW|10-08: ALLSW` (H=TRANSFERS, K=UNKNOWN until K16).
    It would arm s1 CPH and s2 459dr. If K16 reads ≥ MID, K20 heads s2. K20 was scored on 10-04, so the novelty gate skips it
    and 459dr follows (rule lines 153–157 always append the tail), so 10-07 = CPH + 459dr either way.
- 10-05 02:05Z: **recipe_dr524 FAILED on memory at 02:01:12Z.** def_pipeline refused at free+inactive 22.5 < 23 GiB, 31 s after
  the 25 GiB release; two GlowUp iOS Simulators hold the memory. The builder was restored == HEAD. **664's rerun_dr524.sh is armed**
  (pid 18948): it re-runs the unchanged script once free+inactive ≥ 27 GiB twice, 60 s apart, retries only on that refusal, cap
  19:00Z. mde_sign stays PENDING until it exits. CHECKLIST v2 waits too. Ahmet was pushed to quit the Simulator.
- 10-05 05:20Z: **524dr re-run started** (664 attempt 1, 05:18:20Z, at free+inactive 27.0 → 28.9 GiB). The recipe resolved as at 02:00Z
  (mde_sign OFF). Setup passed this time (repin-panel-support 15.0 GiB peak); the LODO key is on, the dump is f39b7f77, and the level
  builds are running. The CHECKLIST v2 and the mde_sign decide follow its exit.
- 10-05 05:21Z: **Network block: Cold Turkey Blocker (Ahmet's app) refuses connections to virtualcellchallenge.org and github.com**
  ("Connection refused" after ~21 ms from 192.168.0.171; google.com and storage.googleapis.com work). It began between 10-04 23:20Z
  (the site answered then) and 10-05 05:20Z. git pushes fail; commits stay local until it lifts. vcc needs the site for whoami,
  submit, status and capture, so **the 10-06 00:01Z submit fails if the block is on then** (arm1006 retries until its 12:00Z cap).
  Ahmet was asked by push to allow both domains. Never change the blocker ourselves (his app and his settings). Before each slot
  window, check `curl -s -o /dev/null -w '%{http_code}' https://virtualcellchallenge.org/api/leaderboard` (200 = open).
- 10-05 07:13Z: **524dr READY** (re-run 07:05Z; gate524 PASS, CLIFF OK 0.9735 / 0.9794 / 0.9884, postbuild 42/42, package cbeadd65;
  8ac61842). **CHECKLIST v2 running (665):** steps 0–6b and 8 PASS (538, 554, 590 and 641 applied; every hash as expected; recipe
  tests 6/6; post594 ALL PASS; gate591 ALL PASS). Step 7's inert KA waits for free+inactive ≥ 27 GiB (pid 83715, cap 16:00Z); then
  steps 9 and 10. The working tree carries the applied bundle, uncommitted until step 9. idall3's build h5ad was swept (scored) for
  disk: 44 GiB free.
- 10-05 09:04Z: **Storage.** Peer [ae1165] (Ahmet's cleanup thread) prepared a sweep that Ahmet runs himself: 85.1 GB of byte-identical
  copies (aa7ebefb). It also backfilled 2,009 small outputs (ef754c86). The applied bundle files are untouched (all 13 verified
  against byte copies). The coordinator swept 497dr's source h5ad: never submitted, no rule after 10-05 lists it (line 1259 kept it
  through 10-05 for gate497's readiness checks), and its .vcc sha 6ed67b2a equals gates497.log. cph549 (10-07 s1) and 524dr
  (10-06 s1) keep their h5ads until scored.
- 10-05 10:36Z: **10-05 CHECKLIST v2 COMPLETE.** Step 7's inert KA PASS (09:39–10:35Z after both GlowUp simulators were shut down by
  their sessions): output ac49d319 == 526/efp_inert, dose_guard per context. Step 9 committed 633d0ab6 (538 → 554 → 590 → 641 +
  run_arm_generic 9a718b9e + test_october_recipe OFF538). Step 10: **`mde_sign` DECIDED clean4:0.05** (496a TRANSFERS) and
  **`hook_estimator` DECIDED off** (652 F1). The recipe now resolves with `DEF_DAY_MDE_SIGN=clean4:0.05`. Still PENDING: rho_scale,
  rc_callset_cpm (K16 10-06), lodo_strength, lodo_hi, callset_block (CPH 10-07). Step 11: chain1006 arms 10-06 at 20:00Z
  (s1 524dr READY, s2 K16). Next: **10-06 APPLY DAY v2**.
- 10-05 10:37Z: **PLAN FOR TONIGHT AND 10-06 (coordinator; read this first after any restart or compaction).**
  1. **≥ 18:00Z (the 606 freeze; last PINNED block 10-02 14:56Z, 85d4be68): start "10-06 APPLY DAY v2" early.** Steps 1–4b are
     patches plus suites. Step 5 is run_builder_tests. Step 6 is six KAs, serial, under the memory gate (each a full build, ~20 GiB,
     30–75 min): all-unset ec3a5ae6, KA598b dd82740d, KA598c, KA612b e21ff752, KA612c, KA630 + ka630_audit. Step 7 adds the post-598
     def_pipeline sha to run_arm_generic. Step 8 commits and sends [537757] the shas. Step 9 decides sharpen_pt tercilec after
     KA598b/c (C10b f6aa2ed9 + 599 7db7938c). arm1006 checks only package shas, and chain1006 reads only gate logs, so applying
     tonight cannot touch 10-06's slots.
  2. **20:00Z:** chain1006 (pid 40686) arms 10-06: s1 524dr (READY, cbeadd65), s2 K16.
  3. **00:01Z 10-06:** both submit. **Read 524dr** with `uv run python experiments/658-prereg-1007-builds/read658.py dr <524dr id>
     --recipe experiments/524-october-recipe-dr/work/recipe524.json --hook HoQY4xrzYer5y232oh85`, plus gate524 / 524's own prereg.
     On LAG or exit 3, write experiments/660-slot-1007/work/HOLD-1007-s1 BEFORE 10-06 20:00Z (R6). **Read K16** with 523's bins,
     PRED-662-K16 (I-662 −0.0067 overall vs f1 −0.0030) and the 610 card; October K stays 2,200 unless the read says otherwise.
  4. **10-06, once APPLY DAY step 8 is committed:** launch the three 10-07 builds in order (build_630arm.sh → build_allsw658.sh with
     the post-apply builder/emitter shas → recipe_dr524dr2.sh). Each must print READY before 10-07 20:00Z (chain1008 pid 48648).
  5. **10-06 20:00Z:** chain1007 (pid 39901) arms 10-07: CPH + 459dr. **10-07 ~04:30Z:** the R7 message to [10eca4].
  6. Before each slot window, curl the site (Cold Turkey). Re-create the heartbeat cron by 10-09.
- 10-05 10:38Z: **APPLY DAY v2 preflight on the real post-10-05 HEAD (10156278), read-only:** 598 `--check` clean and test598 43/43
  (T0 against the live builder). compose606 ALL PASS: 641 skipped (on HEAD), then 598 → 612 → 612dd → 630 → 630dd → 642 → 604 → 644,
  each with --check plus apply in a temp tree; `bash -n` OK. Final shas match the PINNED 14:56Z block (def_day 54fcf040, def_sameline
  5cd3ab07, def_identity 7073406f, def_day_helpers fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914). The apply starts at
  18:00Z as planned.
- 10-05 14:29Z: **Disk plan through 10-07** (41 GiB free now, sliding ~1 GiB/h from other sessions). Tonight's APPLY DAY KAs each write
  one ~12 GB h5ad and delete it after hashing (run_ka_527.sh line 34), so 41 GiB is enough. The three 10-07 builds need ~45–60 GB:
  630arm and ALLSW keep 12 GB h5ad + 2.8 GB vcc each; 524dr2 is a def_day with ~36 GB of level h5ads at peak. Before launching
  them on 10-06, sweep (all .vcc + sha logged):
  1. the inert KA's 474efp_final.h5ad (11.9 GB, one hard-linked inode with lvl1; hash ac49d319 recorded in 665, as 526 kept only
     metas);
  2. submissions/524dr-october-recipe-abc.h5ad once 524dr is captured (12 GB);
  3. cph549.h5ad once CPH is captured on 10-07.
  [ae1165]'s 85 GB sweep (Ahmet runs it) would make this moot.
- 10-05 18:19Z: **APPLY DAY v2 started at the freeze (666).** Steps 1–5 and 7 PASS: 598, 612, 612dd, 630, 630dd, 642, 604 and 644 v3
  applied; every suite green; 12/12 shas == the PINNED block; builder tests 56/0/1. test_hek293_642 ran 6/6 in a pre-apply
  worktree: its setUp re-applies 642, so it must run before the apply (checklist defect, recorded in 666). Step 6's six KAs are
  running (pid 26113, ~3.5 h). Step 8 (commit + shas to [537757]) and step 9 (sharpen_pt tercilec) follow KAS ALL PASS. The live
  tree carries the patches uncommitted until then (byte copy in the scratchpad).
- 10-05 20:29Z: **10-06 ARMED** by chain1006 at 20:01:28Z: s1 524dr (pid 71273), s2 K16 446bj2zk16 (pid 71275). Both submit at
  00:01Z; s2 waits for s1's capture. APPLY DAY KAs: KAu PASS (== ec3a5ae6, all knobs inert), KAb PASS (== dd82740d), KAc PASS
  (ka598.py). KAh is building, then KAg and KAr.
- 10-05 21:56Z: **APPLY DAY KAs: KAu, KAb, KAc, KAh and KAg PASS. KAr failed on a FULL DISK, not on 630.** Its build crashed with
  OSError Errno 28 "No space left on device" while writing the h5ad at 21:46:50Z, then SIGSEGV at shutdown, rc 139. Within ~20 min
  the volume went from 33 to 4.2 GiB:
  - another session downloaded 11 GB of ACE-Step checkpoints into ~/Code/_vendor (done 21:40Z);
  - a 1.2 GB git pack appeared in this repo at 21:29Z;
  - the six KA dirs held six real 1.3 GB copies of whitener e2f65ecf.
  Fixes, no information lost:
  - the six whiteners → one file plus five APFS clones (all e2f65ecf), +6.4 GB;
  - the inert KA's verified 474efp_final.h5ad removed (hard-linked with lvl1; ac49d319 recorded in 665 and 633d0ab6), +11.9 GB.
  Disk is now 33 GiB. **KAr re-run as KA1006r2** (run_kar2.sh, 21:55:48Z, same env and checks).
  Lesson: run_ka_527 leaves a 1.3 GB whitener copy per KA dir. Clone them after a KA batch, and check disk before each heavy step
  when other sessions are active.
- 10-05 23:21Z: **APPLY DAY v2 COMPLETE (666).** KAS ALL PASS (KAr via KA1006r2 after the full-disk failure). Step 8 committed as
  76d0b2ed (pins 12/12); step 9 decided **sharpen_pt tercilec:2.5,2.0,1.75** (5ee58ffc). [537757] has the shas. The October recipe
  now resolves to: DEF_DAY_LODO=1, LODO_HI_PCT 96.647, SHARPEN_PT tercilec:2.5,2.0,1.75, DE_MAG dn_k250_zs, MDE_SIGN clean4:0.05
  (pending: rho_scale, rc_callset_cpm, lodo_strength, lodo_hi, callset_block).
  - **10-07 builds started early:** build_630arm.sh launched at 23:20:50Z (pid 29138; DRY gates PASS). Disk was freed for it:
    cph549.h5ad removed (CPH's .vcc sha de1e4ab9 matches its run log; restorable), giving 44 GiB.
  - **Next:** after 630arm is READY, build_allsw658.sh with ALLSW_BUILDER_SHA=fe4678c31b3c1f282976ad4fbbc23b44ffe1efc43f151fd832e0b139d1a44a89
    ALLSW_EMITTER_SHA=0f56209d57e8385abd08a7a15692ae3241b9373b436b4921546e546d571e5876, then recipe_dr524dr2.sh with the same shas
    (DR2_*). Before each, keep ≥ 35 GiB free (sweep 524dr's h5ad after its capture). Ahmet's 85 GB sweep (aa7ebefb) would end
    the squeeze.
  - 00:01Z: the 10-06 slots (524dr, K16) submit; read them per the PLAN entry.
- 10-05 23:25Z: **644 re-verify on the applied tree done early by [537757]** (3c063fac): verify644 37/37, including new rows K1–K5.
  K1: the recipe emits tercilec / LODO pct / zs / clean4:0.05. K2–K5: def_day, def_sameline and the builder carry 598 / 612dd /
  630dd / 644 / 604 / KEYGATE / UPRN / tercilec. The 641 / 642 / 644 pre-apply suites now SKIP with a reason on an applied tree.
  The coordinator re-ran verify644: all_ok True; its rewritten run records were restored to 3c063fac's. Re-run verify644 (K1)
  whenever KEYGATE / UPRN / an identity prior value is decided.
- 10-06 01:59Z: **630arm READY** (00:42:28Z): package 579132e7665ba64640091418b287c6bc535adbb35aef34c82801407973eb34fd; CLIFF OK, prep rc 0,
  novelty rc 0, audit per build_630arm. **ALLSW building** since 01:59:10Z (build_allsw658.sh with ALLSW_FLOOR = now, overriding
  the script's 10-07 00:00Z default floor: ALLSW needs only decided switches; the hook is ON because mde_sign is DECIDED;
  EXPECT = de_mag, de_mde_sign, de_sign_emitter_sha256). 524dr (10-06 s1) was submitted at 00:01:42Z, uploaded, and is still
  scoring at 02:00Z (496a took ~62 min). K16 (s2) waits behind it. Disk 55 GiB.
- 10-06 03:33Z: **ALLSW READY by coordinator override** (package c86515cc; experiments/658-prereg-1007-builds/DIAGNOSIS-ALLSW-bulkid.md).
  bulkid658 found 28 paired ±1 count swaps (net 0 per block, 11 of 900 blocks): rounding ties under emitter 0f56209d vs 498p2g's
  edb5b8dc, a move EXPECT allows. F1 is a board rule (1e-6 raw), predicted PASS (idall3 precedent: pds identical, mse +9.9e-8), and
  read658 still checks it on the read. Every other gate PASS. **524dr2 building** since 03:33Z: branch N3 (tercilec + mde_sign new
  vs 524dr), DR2_FLOOR = now. Disk: 630arm's and ALLSW's h5ads swept (both .vcc shas == logged), so 65 GiB. 524dr (10-06 s1) is
  still in flight at 03:3xZ (server shows in_flight 1, spent 0). K16 waits behind it.
- 10-06 03:34Z: **Cold Turkey is blocking virtualcellchallenge.org and github.com again** (refused after ~23 ms; the site answered at
  02:29Z), the same morning window as 10-05 (blocked by 05:20Z, open by ~07:40Z). 524dr's `vcc submit --wait` survives it:
  connection errors are a retryable ApiError, and only 400/401/403/404 are fatal. It resumes when the site opens. K16 (s2) cannot
  submit until then (arm1006 cap 12:00Z). Pushes queue locally (e0e47ff0 and later). Ahmet was pushed to allow both domains. If
  the capture step fails while blocked, the arm logs the entry id; re-run `scripts/capture_board_entry.py <id>` once open.
- 10-06 04:30Z: **524dr looks STUCK server-side.** It was submitted at 00:01:42Z and is still in flight at 04:30Z (limits: in_flight 1,
  spent 0), while other teams' entries submitted up to 03:36Z have been scored. Normal is 30–62 min (496a 62, idall3 33). The site
  answers again (Cold Turkey's block on it lifted by ~04:29Z; GitHub is still blocked). The entry id exists only in the waiting
  `vcc submit --wait` process (pid 35534). The web app (Home, Leaderboard, Team, Credentials) has no submissions list.
  **Decision point ~06:30Z:** if still in flight, kill pid 35534 so arm1006 s1 logs vcc's output, which carries the entry id; then
  `vcc status <id>`. If it is stuck, `vcc cancel <id>` (frees the slot; only scored entries count) and resubmit 524dr's package
  by hand with the same -m name. It was never scored or ledgered, so novelty passes. K16 (s2) "goes alone" once s1 ends, retries
  until 12:00Z and probably takes the freed slot first; 524dr then follows as the day's second entry. Its R6 read must land
  before 20:00Z.
- 10-06 05:32Z: **10-06 slot recovery (667).** 524dr's entry **VqBgraiBxlaK0WP5efOF** was stuck in "launching" for 5.5 h: its Cloud
  Batch scoring job (us-east1) never got a machine, while other teams' later entries scored.
  - The waiting vcc was SIGTERMed. SIGINT is ignored by nohup'd background jobs; click.echo had already flushed "entry: <id>"
    into the arm's log.
  - `vcc cancel --yes` returned "stopped the queued scoring job before it started running". in_flight 0, spent 0, so no slot was
    used.
  - **K16 took the freed slot** (submitting 05:31:20Z, arm pid 71275).
  - **524dr is resubmitted after K16** by experiments/667-slot1006-recovery/resubmit524dr.sh (pid 81769). It runs chain1006's
    own slot-1 command once K16's arm exits; cap 11:30Z.
  - Reads: 524dr with read658 dr (R6, before 20:00Z); K16 with 523's bins + PRED-662-K16.
  - **All three 10-07 builds READY:** 630arm 579132e7; ALLSW c86515cc (override); 524dr2 93e49b81 (05:15Z, gate524 rc 0).
  - Lesson: a stuck "launching" entry shows as in_flight 1 / spent 0 in /api/cli/submissions/limits. Its id is only in the waiting
    vcc's output. SIGTERM the wait, read "entry:" from the arm's log, `vcc status`, then `vcc cancel --yes` while it is still queued.
- 10-06 06:30Z: **K16 also in flight 1 h** (submitted 05:31:20Z; limits in_flight 1, spent 0). Other teams' entries submitted at
  05:26–06:03Z are already scored, so our scoring jobs, not the scorer, are failing to get machines. Hypothesis, not verified: our
  ~2.8 GB packages need a large machine type that is short in us-east1 tonight (496a and idall3 scored in 33–62 min on 10-05).
  Plan: wait. If K16 is still "launching" by ~08:30Z, get its id the same way (SIGTERM the wait, read "entry:" from the arm's log),
  run `vcc status`, and tell Ahmet. Contacting support (help@virtualcellchallenge.org / Discord) is his decision. 524dr's
  resubmission (667) still waits behind K16.
- 10-06 10:19Z: **K16 queued 4 h 47 min** (submitted 05:31:20Z; in_flight 1, spent 0). The 08:30Z check was missed: the coordinator
  session restarted and was idle 07:29–10:18Z. Both 10-06 entries got stuck in the scoring queue, while other teams score normally.
  - Likely mechanism, from the client: `launch_submission` sends the package's nnz (~2e9 for ours), so the server presumably sizes
    the Cloud Batch machine by it. A large-memory type that is short in us-east1 would queue only big entries. On 10-04/05 our
    jobs got machines within 33–62 min.
  - **Escalated to Ahmet** (contacting VCC support is his call). The push failed (Remote Control inactive), so it is in the chat.
  - Plan: keep K16 queued, so it keeps its position. 524dr resubmits only if K16 finishes by 11:30Z (667 cap). Otherwise 524dr is
    not read today, and 524dr2 (READY) supersedes it in the 10-08 queue.
  - **Risk for 10-07:** if the queue problem persists, a still-queued K16 blocks 10-07's 00:01Z submits ("one in flight").
    Decide by ~23:00Z whether to cancel K16 (cancel works while queued and spends no slot).
- 10-06 10:20Z: **guard1006 armed** (667, pid 29602). At 23:00Z, if K16's arm is alive, nothing was captured and in_flight is 1, it
  SIGTERMs K16's vcc wait, reads "entry: <id>" from submit1006s2.log, and runs `vcc cancel --yes` only if the status is still
  "launching". This keeps 10-07's 00:01Z submits unblocked. If 524dr is resubmitted (only if K16 finishes by 11:30Z) and then gets
  stuck, cancel it by hand before 23:30Z: the guard handles K16 only.
- 10-06 15:29Z: K16 still queued (10 h). **guard_inflight.sh** (667, generic form of guard1006) armed for **10-07 23:00Z** (pid 80420)
  over submit1007s1/s2.log. If our 10-07 entry is still queued then, it cancels it so 10-08's 00:01Z submits are not blocked. Arm
  the same for 10-08 (chain1008's logs) if the queue problem persists.
- 10-06 18:32Z: **K16 read** (Imj3EXViNco8rT4eLLlE, 0.2707192): K 1,600 hurts, so **rc_callset_cpm DECIDED off** (e2d03ab1). **524dr
  resubmitted with today's second slot** at 18:32:02Z through scripts/arm_submit.sh (log 667/work/submit524dr_late.log). Its board
  entry is the base of 524dr2's preregistered primary read (N3 vs 524dr). **guard_inflight armed for 23:30Z today** (pid 31874):
  if 524dr is still queued then, it is cancelled, so CPH's 00:01Z submit is not blocked.
- 10-06 22:11Z: **REBOOT at 18:53Z killed every waiter.** Recovered per AGENTS:
  - builder, emitter, def_day and def_pipeline == HEAD (nothing pinned);
  - **chain1007** re-armed: it fired at once and armed **10-07: s1 CPH (pid 26353), s2 459dr (pid 26355)** (rule: K=NEG from K16,
    H=TRANSFERS);
  - **chain1008** re-armed (pid 26439; 10-07 20:00Z);
  - guard_inflight for 10-07 23:00Z re-armed (pid 26440);
  - caffeinate, board_snapshot.sh and pisces_watch.sh restarted.
  **524dr's 18:32Z upload was cut by the reboot** (entry KTzIDR0Occ2jYzklkARE, pending upload). resume524dr.sh resumes that same entry
  (pid of vcc 26595), records novelty and captures. guard524dr.sh (pid 26745) cancels it at 23:30Z if it is still uploading or
  launching, so CPH's 00:01Z submit is not blocked. The generic guard for 23:30Z was replaced: a resumed submit logs "resuming
  entry", not "entry:".
- 10-07 04:31Z: **R7 POSTPONED (the reserve stays sealed).** The ledger (CONFIRMATION-LEDGER line 1050) says R7 is accessed only
  after the DE switches settle, not before the 10-07 reads. CPH (10-07 s1), which decides `callset_block`, has been queued on the
  scorer since 00:01Z and is unread. The 06:00–12:00Z window was only "tentatively agreed" (552 §2). Also, the R7 runner session
  ("Explore contribution areas" [10eca4]) is absent after the 10-06 reboot: ListAgents shows new refs, and this coordinator is
  now [50e998]. Inputs decided so far: F = M6+P2g (lodo_hi_pct 10-03T10:02:43Z; sharpen_pt tercilec 10-05T23:19:45Z), H on
  (mde_sign clean4:0.05, 10-05T10:36:41Z), K default (rc_callset_cpm off, 10-06T18:30:57Z), hook_estimator off. B is PENDING until
  CPH reads. **Next:** once CPH is read and callset_block decided, fill R7-MESSAGE-DRAFT.md and send it to the R7 peer (Ahmet must
  reopen that session), or run 552's §6 runbook in the coordinator session, with a new window.
- 10-14 22:23Z: **RESTART after a 7-day idle gap (coordinator [400b6e], one of 5 sessions Ahmet started 10-14 ~22:10Z).** Nothing ran
  10-08 ~16:00Z → 10-14 22:12Z: the last commit was 10-07 04:31Z, the Mac rebooted ~10-11 18:30Z (uptime 3 d), and no waiter, chain,
  caffeinate or snapshot was alive. The 10-07 / 10-08 entries had scored and been captured (data/board-entries/) but were never read:
  now in BOARD-LEDGER (CPH: callset_block does not ship; 630arm PARTIAL; 524dr2 +0.0007; 459dr not October material). Slots 10-09 →
  10-13 (10 slots) were lost.
  - **Roles (agreed 22:20Z; [020377]'s deciding message):** [400b6e] coordinator + slots (670–679); [020377] news + intel (680–689;
    docs/NEWS-2026-10-15.md); [17adab] D/E/F day-of readiness (690–699); [bf003b] science / transfer (700–719; the field gap first).
    Only the coordinator asks Ahmet anything. One heavy job machine-wide; builds are serial (def_day pins the builder).
  - **The field moved:** 61 teams above our best (30 on 10-08), 19 above 0.30, Illumina AI 0.4729. Our displayed row 524dr2 0.2731 is
    rank 66. The pack beats us on every member but jac; the biggest gap is fid (0.04–0.09 scaled vs our 0.005). [020377]: the
    look-alike rows are multi-account operators, not a public recipe; cell-eval2 0.18.0 (10-01) changes no scored number; anchors are
    still r4. Science's first job is the fid/nmae gap.
  - **10-14 slot (expiring at 00:00Z): ALLSW** c86515cc (658's frozen prereg, reader `read658.py allsw`; predicted ≈ 0.2755 + the hook's
    +0.0006; R1 decides whether a D/E/F entry may pair one bulk + one DE switch). Submitted 22:16Z via arm_submit.sh
    (log experiments/670-coord-1015/work/submit1014_allsw.log). It is ALSO the 10-14 board read; nothing else was ready.
  - **10-15 s1: DN (631's FROZEN prereg, base = 630arm's entry, bins on Δnmae raw: ≤ −0.0025 TRANSFERS, −0.0025 … −0.0008 PARTIAL).**
    patch_builder631 applied 22:20Z (builder fe4678c3 → 8f8f6d5c, dn_mag.py 0a8cc560; run_builder_tests 56/0/1; test631 26/31, the 5
    fails are its pre-apply re-apply checks, as with 642). experiments/671-dnmag-arm/dn671.sh builds DN, then KA1015u (all-unset ==
    ec3a5ae6). The patch is committed only after KAu PASS.
  - **10-15 s2: science's nomination** (a field-gap probe, preregistered), else the next queue item.
- 10-15 01:54Z: **ALLSW read** 0.2760851 NEW BEST (R1 PASS, R2 ADDITIVE). **DN (s1) read** 0.2741292, PARTIAL → de_dnmag d_m200 flagged.
  631 committed (17adc547, KA1015u PASS); 690a/c/g applied (da052119); de_uprn ON, mde_sign_keygate off. s2: choose1015s2.sh (PIE 702 if
  READY before 18:00Z, else ku30). Day plan: docs/PLAN-2026-10-15.md.
- 10-15 13:33Z: **SLEEP + SHUTDOWN + REBOOT.** pmset log: at 05:32Z the lid was closed and AC unplugged ("Clamshell Sleep", then "Using
  Batt"). The Mac slept with only maintenance wakes until it shut down at ~12:55Z, then rebooted at 13:03Z, still on battery (70 %, ~1 h
  left at 13:30Z). Board snapshots stop at 04:16Z. Lost: DNv2 (ka693b's second KA, never finished: it slept, the module is not slow),
  R7a (never reached its fits; resumable, still the single access), run690, rabc694 (NOT READY on the shutdown), every waiter, watch691 +
  r0ready692. Recovered: builder restored from build_dn631v2's pin record (== HEAD e99bbd4a); caffeinate, board_snapshot, pisces_watch
  restarted; 10-16 s2 (ku30) and 10-17 s1 (lstr705) waiters re-armed. **No builds until AC is back.** Peer sessions restarted with new
  refs ([dc639e], [642c06], [3c04af]); the coordinator is now [7ab092]. Ahmet notified to plug in and keep the lid open.
  Re-plan once AC is back: DNv2 → R-ABC (10-16 s1) → K-UP (10-16 s2), with R7a interleaving on gates; R-NEW / R-KA move after them.
- 10-15 21:21Z: **Second loss of the day:** every job started at 13:33Z (caffeinate, snapshots, both waiters) was dead by 21:19Z
  without a reboot (lid closed again ~13:40Z; the app was most likely quit, killing its process group). AC is back. Fix: scripts/detach.sh
  (setsid), now in AGENTS §6.2. Relaunched detached: caffeinate, board_snapshot, pisces_watch, ka693b (KAU2 + DNv2, 21:20Z), waiters
  1016s2 (ku30) and 1017s1 (lstr705). **Builder order:** ka693b → R-ABC (10-16 s1) → build701 K-UP (10-16 s2) → R7a → build705 (10-17 s1) →
  run690 (R-NEW + R-KA) → build707. Sessions restarted again: coordinator [2d6f46], science [d80a10], peers [4592ee] / [730b99].
- 10-16 04:22Z: **R-ABC (s1) read:** 0.2743318, +0.0012429 vs 524dr2, COMPOSES (bulk exact, nmae −0.0045 raw). s2 = K-UP (build701) when READY.
- 10-16 04:34Z: **No swap file since the reboot** (vm.swapusage total 0.00M) made every "swap free ≥ 400 MiB" gate unpassable. gate602.sh fixed
  (total 0 = pass; 1e7c26b2). The same patch to gate551.py broke R7a: gate551 is inside R7a's pinned instrument closure (r7a551 asserts byte-identity
  to the efp-KA commit e7be21ca), so R7a's access failed its integrity check at plan stage (04:33Z, no truth opened). gate551.py reverted
  byte-for-byte (3f3dcf58); R7a's route is [d80a10]'s call. **Lesson: never patch a file inside a sealed instrument's pinned closure; check
  r7a551's instrument_files() first.** K-UP (ku30) READY 04:33Z and submitting for 10-16 s2.
- 10-16 05:00Z: **K-UP (s2) read:** 0.2740045 (+0.0011785 vs 446bj2z); primary Δ(fid+jac) +0.000123 = NULL → K stays 2,600. Unregistered: reach +0.0054
  raw (~4σ of the re-roll), asymmetric vs K16/K20. Proposed 10-18 s2: a preregistered pass-1 re-draw replication at CPM 3,000 (primary Δreach, guard Δoverall).
- 10-16 07:31Z: **R-NEW (new-panel rehearsal) FAILED in 6 min:** def_lodo_key.py → sens415.cd4_frozen_and_partial → 062 prepare_cd4_mean(max_rss_gib=8.0):
  "MemoryError: cd4_rest centering: peak RSS 8.530 GiB reached 8.000 GiB". A hard-coded self-guard that a NEW panel exceeds (A/B/C stayed under), so R0 on a
  new D/E/F panel could die the same way. Routed to [4592ee]: make the cap a parameter (default 16 GiB), KA the key byte-identical, grep the day path for
  other 8 GiB caps, re-run R-NEW after R-KA. This is what the rehearsal was for.
- 10-17 00:31Z: **lstr705 (10-17 s1) read: 0.2757977** (−0.0002874 vs ALLSW); primary Δ(pds+mse) −0.000240 > −0.0005 → **CLEARED**: θ 0.25 ships for
  October (DEF_DAY_LODO_STRENGTH=0.25), into R-ABC2 and rider690. Near-free option (Arc-drawn LABC +0.00005; board −0.00029). sw707 (10-17 s2)
  uploading since 00:29:51Z; its read decides 707's hook (TRANSFERS ≥ +0.0010 → R0; not HARM → 708's day-of swap route).
- 10-17 01:33Z: **sw707 (10-17 s2) read: 0.2771181, a NEW BEST** (+0.0010330 vs ALLSW); primary Δ(pds+mse) +0.000842 = **WEAK** → the hook is NOT adopted
  into R0. The swap stays live via 708's day-of rule (not HARM) as a 6a replacement entry after 10-22. **R-ABC2 = R-ABC + θ 0.25 only**:
  pred = 0.2743318 − 0.0002874 = 0.2740444 (±0.0005). Next: 10-18 s1 R-ABC2, 10-18 s2 ku38 (dial test), 10-19 s1 h1s710.
- 10-18 01:45Z: **R-ABC2 (10-18 s1) read: 0.2739849, D −0.0000596 vs pred → COMPOSES.** θ 0.25 stays in R0's recipe; R0 = R-ABC + θ, board-read on
  def_day. 10-18 s2 ku38 (the dial test) has been uploading since 01:40:40Z; read711 decides October K. Then 10-19 s1 h1s710, 10-19 s2 R-ABC3 (only on
  ADOPT K 3,400) or R-ABC2H, 10-20 s1 R-ABC2H (pred 0.2750179).
- 10-18 03:20Z: **ku38 (10-18 s2) read: 0.2736785; read711 → October K stays 2,600** (the dial test refutes f 0.61: Δreach +0.0094 > +0.0019;
  condition (i) at f 0.87 fails). R-ABC3 is not built. The board's reach rises monotonically with K while both H1 dials said otherwise, so H1
  can't price size. **10-19 s2 → R-ABC2H** (moved up from 10-20 s1; pred 0.2750179). 10-19 s1 h1s710 unchanged. 10-20 / 10-21: open.
- 10-19 01:33Z: **h1s710 (10-19 s1) read: 0.2765675 (+0.0004824 vs ALLSW) → POSITIVE, R710 = no block** for 708's day-of H1 move. R0 unchanged.
  10-19 s2 R-ABC2H has been uploading since 01:21:20Z (pred 0.2750179); then 10-20 s1 C1 (pred = R-ABC2H + 0.00089).
- 10-19 01:52Z: **R-ABC2H (10-19 s2) read: 0.2748936, D −0.0001244 vs pred 0.2750179 → COMPOSES.** The 6a candidate (R0 + swap) holds on A/B/C.
  10-20 s1 C1 (armed; floor 00:01Z): pred 0.2757814 (lab 333dddc8); HARM if ≤ 0.2752814 or < 0.2743936 → 811's 6a route closed. 10-20 s2 /
  10-21: unspent unless a decision-relevant arm appears (the October recipe is frozen).
- 10-20 02:40Z: **C1 (10-20 s1) read: 0.2777649, D +0.0019835 vs pred 0.2757814 → ABOVE** (+0.0028713 over R-ABC2H, all bulk/mse). 811's 6a
  route stays open; SIXA-DECISION input 3 = not HARM. X2 still needs the lab's DEFo/DEF dossier + the LZ sign check. 10-20 s2 / 10-21: see below.
- 10-20 02:41Z: **10-20 s2 → C0 (699-c0811)**, armed: arm1015 1020s2 pid 98102 (ready-log = 699's own c0811_699.log, so a NOT READY stops it;
  floor 02:40Z, cap 10-21 18:00Z, so a late READY goes out in 10-21 s1/s2), name_guard 98364, name claimed from the queue, guard_release
  covers 1020s2. Prior = a stub (1020s1 terminal + captured) so an old pid in submit-1020s1.log cannot be recycled into a hold. 699 starts
  only on readiness's GO-c0811 after the lab's READY (hook + KAs + prediction). Reading: SIXA-DECISION Amendment 1 (HARM: D ≤ −0.0005 vs the
  lab's pred, or C0 < 0.2734849 → swap-OFF X2 closed; else the swap-OFF X2 has its A/B/C read). Board row: likely below C1's 0.2778.
- 10-20 02:53Z: **C0 GO.** Lab READY (e603fe73 / 16fb7665 / 017cc763; hook dumphook811c0.sh 297c3832, hook811.py c5a5614a, KA-on 96204332); relayed to
  readiness for GO-c0811. **Pred 0.2763696** (ρ 1.1280, the class median of 707 and C1). HARM if C0 ≤ 0.2758696 or < 0.2734849. Reader `read695.py --c0`.
- 10-20 03:33Z: **10-21 s1 → W1 (673)**, preregistered (5c4f4bc5): R0's exact recipe + DEF_DAY_H1_STACK=x2.5 on def_day / A/B/C, pred 0.2744673; BELOW
  (D ≤ −0.0005) blocks any 708 H1 move above ×2 in the final, otherwise no change (guard only). Asked readiness to build it after C0's READY (one builder
  at a time; no start after 10-21 12:00Z). Reader `read695.py --w1` (dry-run 5 cases). 10-21 s2: open.
- 10-20 03:39Z: **10-21 s1 armed for W1**: arm1015 1021s1 (ready-log = 673's own w1_673.log; floor 10-21 00:01Z, cap 18:00Z), name_guard 84272, queue
  name claimed. Prior = a stub that `prior_relay.sh` (new, 4-case tested) writes only once 1020s2 is terminal AND no `bash arm1015.sh 1020s2` is alive.
  Pids wrapped past 99999 on 10-20, so the first pid in submit-1020s2.log could be recycled into a hold. W1's driver (readiness, pid 83780) waits for C0's
  driver to exit (one builder). If C0 slips past midnight, it takes 10-21 s1 and W1 goes after it. 10-21 s2: open.
- 10-20 04:32Z: **10-21 s2 → W2 (674)**, preregistered (12562add): W1 + DEF_DAY_CD4_STACK=x1.5, read vs W1; pred = W1 + 457dr's +0.0031484; BELOW → X1 may not
  carry CD4 > ×1 (SIXA Amendment 3, applied by the coordinator; decide708 unchanged). Asked readiness to build it after W1 (no start after 10-21 12:00Z).
- 10-20 04:36Z: **10-21 s2 armed for W2**: arm1015 1021s2 (ready-log = 674's w2_674.log; floor 10-21 00:01Z, cap 18:00Z), name_guard 93757, queue name
  claimed, prior = prior_relay stub of 1021s1. Build queue (readiness): C0 699 (94469) → W1 673 (83780) → W2 674 (93562, starts only after "READY w1"),
  ending ~15:30Z on 10-20. Slots: C0 → 10-20 s2, W1 → 10-21 s1, W2 → 10-21 s2. If C0 slips past midnight, everything shifts one slot and W2 is not
  submitted (no allowance; it stops at its cap). No-bare-kill rule for all three drivers (their EXIT traps write into the ready-logs).
- 10-20 07:05Z: **C0 READY** 52410f55 at 07:03:56Z (699: CONVERGED at lvl2, gate699 ALL PASS on the first run, dump == the lab's 96204332). arm1015 1020s2
  submitted at 07:04:22Z: entry **4rq3qHIuv8bHAxpTkUeu** created (name passed the filter), uploading. Read with read695.py --c0 once captured
  (pred 0.2763696; HARM ≤ 0.2758696 or < 0.2734849). W1 (673) def_day started 07:04:20Z → READY ≈ 11:15Z → 10-21 s1.
- 10-20 07:31Z: **C0 (10-20 s2) read: 0.2780624, D +0.0016928 vs pred 0.2763696 → ABOVE** (+0.0040775 over R-ABC2, all mse). Swap-OFF X2 has its A/B/C
  read; both X2 branches now do. C1 − C0 = −0.0002975 (the swap adds nothing once 811 is in, on A/B/C). Team best. Next: 10-21 s1 W1 (building; READY ≈ 11:15Z),
  10-21 s2 W2. prior_relay wrote prior-1021s1.log (entry 4rq3qHIu) at 07:29:55Z, as designed.
- 10-20 07:40Z: SIXA Amendment 4 (before D/E/F data): on a swap-ON day the swap stays in X2 only if C1 − C0 on DEFo > 0 on ≥ 2/3 folds; else X2 = C0 if it
  passes its own bars. Trigger: board bulk C1 − C0 = −0.00031 (readiness). Relayed via OUTBOX (messaging paused).
- 10-20 11:34Z: **W1's gate673 false-FAILed** (11:15Z): the config diff vs R-ABC2 also listed bulk_sharpen_pt, whose DATA-DERIVED fields (edges, key_ratio_c,
  key_ratio_maxdev) move when the H1 dose changes the stack dump; every behavioural field is identical. **Re-gated on the finished build** by
  experiments/670-coord-1015/regate673.py (698 precedent, no rebuild): READY w1 1337cc9b. **10-21 s1 re-armed** (arm1015 on regate673.ready.log /
  regate673.run.log, same claimed name, prior = C0's stub), name_guard + w1route re-armed; the first attempt's log is submit-1021s1.attempt1.log.
  **W2 (674) dropped:** its start rule needed W1's own READY line, it is the lowest-value guard, and its gate has the same spec error. 10-21 s2 is
  unspent. SIXA Amendment 3 therefore never applies (decide6a reads W2 as NOT READ). The same fix went into OCT22-CHECKLIST 6a (a) for an X1 weights move.
- 10-20 15:20Z: SIXA **Amendment 4b** (lab's calibration note; before D/E/F data): on a swap-ON day X2 keeps the swap only if C1:C0 on DEFo has ci95[0] > 0 (both Arc
  panels lean positive while the board reads −0.0003), else C0. decide6a updated (17/17).
- 10-20 19:24Z: guard_release false fire (controls-final listed, not available) cancelled W1's 10-21 s1 waiter at 18:58Z. Guard fixed to require available true;
  NOGOs cleared; **10-21 s1 W1 re-armed** (00:01Z), w1route re-armed.
- 10-21 00:55Z: **W1 (10-21 s1) read: 0.2746084, D +0.0001410 vs pred 0.2744674 → COMPOSES**: no block on a 708 H1 move (w1route left wait708 at R710=POSITIVE).
  10-21 s2: unspent (W2 dropped; nothing decision-relevant can be built before the release). **That was the last A/B/C slot.** Next entry = R0 on D/E/F.
- 10-22 06:09Z: **R0 SUBMITTED as the first final entry** (k2OKr0zPvJq48BCtdUx8, published 06:05Z, members hidden). Final-phase slots: 2/day, only the last entry counts, deadline
  Nov 5 23:59Z. Next: a 6a replacement only if SIXA-DECISION (decide6a.py) says one is eligible: 708's decision (day708 running since 05:21Z) + the lab's DEFo/DEF/LZ dossier.
- 10-22 18:21Z: **X2 SUBMITTED as the replacement final entry** (S5RxOrm4ZbqtnPRIX8NX, published 18:20Z, members hidden). It is now the team's counting final entry. The 6a route
  is complete (one replacement). Until Nov 5 23:59Z: no further final submission unless a NEW decision, fixed in advance and not based on any D/E/F score (none are visible),
  justifies one. The reserved final name stays unused.
- 11-03 01:45Z: **Final-week check (coordinator restart after an 11-day gap; experiments/722-final-week-1103).** X2 S5RxOrm4ZbqtnPRIX8NX is
  published, FINAL, and still the counting entry (package e39d2ba6 == ledger). An independent semantic check of the submitted prediction
  (h5ad 94ae91b6) passes: each context's pseudobulk matches its own controls (r 0.9999 vs 0.74–0.92 cross), depths match, all 300 labels per
  context carry their own knockdown, counts clean and under the cap. No candidate exists (every post-811 lab lane closed or WEAK, LZ spent,
  board dark), so **no further final submission before Nov 5 23:59Z.** Nothing is armed that could submit.

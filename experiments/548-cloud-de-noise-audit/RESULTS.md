# 548 — Is 0.0045 the re-roll of one DE design?

**Verdict.** The claim's first half holds, and it is not a close call. The second half holds only for a narrow class of arms. Relaxing every DE bin by 10–20× would be the wrong repair.

The number in AGENTS rule 13 is the sample standard deviation of six **different emission laws**, not the re-roll of one design. Recomputed from `experiments/230-toll-magnitude/work/fit230.json` (`falsifier2.zero_sd`, line 262): sd of scaled Δ`reach` = **0.027010964677330575**, and 0.027010964677330575 / 6 = **0.004501827**. Those six arms do not share an emission design (`fit230.py:57-68`).

Arms that keep the **current** DE design (same emitter, sign prior, MDE, zero repair, DE magnitude, call set, ballast and seeds) and change only the bulk move the four DE members by **rms 0.000256 overall** across the four scored steps from 433h to 446bd2. That is the ~0.0002–0.0004 in the claim. It is leakage of specific bulk knobs, measured on four steps, not a re-roll sd.

Two things the claim does not survive:

- The same "DE design held fixed" test on the earlier family (222 through 319p2) gives DE-overall rms **0.00321**, almost all of it `nmae`, because `nmae` reads the profile. Ordering members on that family (fid, reach, jac) are rms **0.00034** overall. A single 0.0003 sd does not describe every bulk arm.
- No tracked submission pair differs only in a seed. Every meta that records `emission_seed` has the value 1000003 (99/99). The re-roll of one design is not a board measurement. On H1, the production construction's seed sd is **0.00069 overall-equivalent** over 85 targets (`435`, three seeds, DE members only). Dividing by sqrt(900/85) to get ~0.00021 on the board is an assumption. Reach's seed difference is carried by about five targets.

`experiments/547-de-noise-check/` and `experiments/537-arm-pipeline/` are not in this checkout (`origin/master` at `c7db045c`). The pair list below is reconstructed from tracked metas and `data/board-entries/`, not replayed from 547. Script: `experiments/548-cloud-de-noise-audit/audit548.py`.

---

## 1. What "rms" and "quadrature" mean here

For a set of scored pairs, Δ is arm B minus arm A in **scaled** board members. DE-overall is (Δfid + Δreach + Δjac + Δnmae) / 6, which is the pair's contribution to overall from the four DE members (and equals Δoverall when pds and mse do not move).

- **rms** = sqrt(mean of squares) of that DE-overall delta. Not a sample sd around the mean. These deltas are not centered: several point the same way.
- **quadrature** = sqrt(rms(Δfid)² + rms(Δreach)² + rms(Δjac)² + rms(Δnmae)²) / 6, with each member still in scaled units. It treats the four members as if they were independent. When one member dominates a step, quadrature and the direct rms are close. When they move together, quadrature is larger than the sd of the sum.

Primary pair set: time-consecutive arms inside a clique whose DE fingerprint matches. All-pairs numbers are also reported; they double-count a clique.

The fingerprint (`audit548.py`, `fingerprint`) includes the emitter and `de_sign_emitter_sha256`, the emission law (`emission`, `overdispersion`, `dispersion_argument`, `spread*`, `resample_*`, `tail_spread`), every seed key (`emission_seed`, `library_pick_seed`, `seed`, `rc_jitter_seed_base`), the call set (`rc_callset_*`, `rc_key`, `rc_sign_policy`, `rc_scorer_min_cpm`, `keep_top`, `keep_spread`, `keep_key`), the sign prior (name, sha, protect q, shrink floor, classes, context-A prior, `ctx_a_prior_*`), MDE (`c`, `cap`, `coarse`, `min_count`, `scope`), zero-repair floor, DE-magnitude variant and module sha, ballast (`kmax`, `lb`), `min_cpm`, `gate_unexpressed`, `mass_fix`, `sign_repair`, `carrier_gamma`, and sameline (file sha, `prog`, `pt`, `w`, and per-context `programme_default` / `w` / `pt`) plus the prior gate's `q` / `tau` / `weak` / `strength`. Missing is not equal to an explicit value. Bulk knobs (sharpen q, bulk load, CD4 dose, normrel, whitener, alpha, `pert_scale` weights, union g) are outside the fingerprint on purpose: they are what a "same DE design" pair is allowed to change.

A path inside a strict pair whose name looks DE-related would mean the fingerprint missed a key. There are **0** such paths (`pairs548.json`, `missed_deish_paths_inside_strict_pairs`).

---

## 2. Pairs, and what each questionable inclusion does

34/34 captured board entries join to a meta. Three cliques have size ≥ 2.

| clique | arms, in submission order | DE-overall sd across arms |
|---|---|---:|
| early profile family | 222, 229, 249, 312, 312b, 316, 319, 319d3, 254, 319p, 319p2 | 0.00208 |
| unwelded sharpen / stack | 410b, 417a, 417e | 0.000096 |
| current design | 433h, 433k, 446b, 446bj2, 446bd2 | 0.000400 |

### 2.1 Pooled consecutive pairs — the number you get if you do not split families

16 consecutive pairs. DE-overall rms **0.002542**. Quadrature **0.002627**. All-unordered-pairs (68) rms **0.002652**, quadrature **0.002661**. None of these pairs has |Δpds| and |Δmse| both under 1e-9: a bulk step moves the profile, so "profile pinned" is the wrong filter for this question. pds and mse moving is the bulk arm working.

Leave-one-pair-out on the 16: dropping 222→229 (DE-overall −0.00700, an `nmae` move of −0.0430 scaled) lowers the rms from 0.00254 to **0.00191**. Dropping 229→249 (+0.00661) lowers it to **0.00200**. Dropping any modern step does not move the pooled rms (`pairs548.json`, `leave_one_out_consecutive`). The pooled 0.0025 is the early family's `nmae`, not a re-roll.

### 2.2 The current design, which is the claim

Four steps, `families548.json` `modern_433h_to_446bd2`:

| step | what the config changed | Δ DE-overall | ordering (fid+reach+jac)/6 | nmae/6 |
|---|---|---:|---:|---:|
| 433h → 433k | `bulk_ctx_scale`, `bulk_sharpen_q` | −0.000213 | +0.000028 | −0.000241 |
| 433k → 446b | same two | −0.000241 | +0.000137 | −0.000378 |
| 446b → 446bj2 | relevel + LODO load | −0.000388 | −0.000372 | −0.000016 |
| 446bj2 → 446bd2 | relevel + DLD-1 load | −0.000093 | −0.000075 | −0.000018 |

- DE-overall rms **0.000256**. Quadrature **0.000300**.
- Ordering-only rms **0.000202**. `nmae`/6 rms **0.000224**.
- Member rms, scaled: fid 0.000046, reach **0.001190**, jac 0.000041, nmae **0.001346**.
- Leave-one-step-out: dropping 446b→446bj2 (the reach step, −0.00220 scaled reach) lowers the rms to **0.000193**. Dropping 446bj2→446bd2 raises it to **0.000291**. No single step is the whole result, and none of them is near 0.0045.

433k→446b's `nmae` is the ballast, not a mystery residual. Scaled Δnmae −0.00227 / slope −1.6532 = **+0.00137 raw**. FINDINGS §0 records that step as nmae +0.0013 raw from median K 4→6 (`docs/FINDINGS.md:131`; `docs/BOARD-LEDGER.md:56`). 446b→446bj2 is a reach step (−0.00220 scaled ≈ −0.00195 raw), which the same section called pinned next to an mse gain of −0.0145 raw. Both are real file differences at a fixed seed. They are the right order of magnitude for the claim. They are not Gaussian re-roll draws.

The sharpen clique is smaller still: 410b→417a +0.000180, 417a→417e −0.000031, rms **0.000129**, quadrature **0.000170**.

### 2.3 The early family, which breaks a universal 0.0003

Ten steps, 222→319p2. DE-overall rms **0.003211**. Quadrature **0.003317**. Of that, `nmae`/6 rms is **0.003304** and ordering-only rms is **0.000339**. Scaled `nmae` rms is 0.01982. The two union/CD4 steps are Δnmae −0.0430 and +0.0427 scaled (222→229, 229→249), about **0.007 overall** from `nmae` alone. Those arms share an emission design and still move a DE member, because `nmae` is a magnitude metric on the reference's DE genes (`AGENTS.md` §4). Calling 0.0003 the sd of "the DE members" on every bulk arm would hide a real `nmae` effect.

### 2.4 Sameline, priors, ballast, repair, magnitude — inclusions that are not noise

Dropping one fingerprint component and re-chaining consecutive pairs (`families` via `pairs548.json` `sensitivity_drop_one_component`):

| dropped from the fingerprint | extra consecutive links | what those links are | their Δ DE-overall |
|---|---:|---|---:|
| sameline (and the prior gate) | 5 | 446bd2→458d, 458d→458g, plus three October-path links | **−0.00803**, **+0.00357**, and three links of −0.00021 / +0.00005 / +0.00007 |
| `ctx_a_prior_*` | 2 | 319p→319pA→319p2 | −0.00021, +0.00033 |
| ballast | 1 | 410→410b | **+0.01710** |
| emitter sha | 1 | 395→410 (sha and `pert_scale`) | **−0.01692** |
| zero repair | 3 | 433a→433h; 446bj2→446br2→446bd2 | **+0.00691**; +0.00076; −0.00085 |
| DE magnitude | 2 | 446bj2→446bj2z (pds/mse pinned); then on to 446bd2 | **+0.00055**; −0.00064 |
| MDE, call set, or the seed keys | 0 | already identical inside the cliques | — |

Putting 458d and 458g into the modern chain changes that chain's rms from **0.000256 to 0.00359** (`families548.json` `modern_plus_programme_arms`). 458d is the programme prior (FINDINGS §0: fid −0.0055 raw, overall −0.0085 vs 446b). 458g vs 458d is +0.003566 DE-overall with pds and mse unchanged at the 1e-9 test: the gate's effect, profile-pinned, exactly the class FINDINGS §0 says shares its draws (`docs/FINDINGS.md:418-423`). Treating sameline as a non-DE flag manufactures a "noise" of several thousandths out of the prior.

446b→446bja is not in the strict set. The fingerprint differs in `de_sign_prior` (`clean4q80` vs `clean4q80ja`), `de_sign_prior_context_A` (`nadig_jurkat`), and the prior sha. Board delta: pds **0**, mse **−2.5e-8**, reach **−0.00847** scaled, DE-overall **−0.00107**. That is the Jurkat-prior effect (FINDINGS §0: reach −0.0079 raw; −0.00847 / 1.1292 = −0.00750 raw on this pairing). It is larger than the bulk drift and far smaller than 0.009. A rule that needs >0.009 before a DE delta counts would have called this noise. The read at the time did not, because the prediction was about fid and the reach move had a sign. The arithmetic point stands: a shared-draw prior delta of 0.001 overall is a measurement.

453dr vs 457dr differ in the sameline **file sha** with `prog: false` on both (relabelled context D, per-target rows, no programme). DE-overall **+0.000048**. Splitting them on the sha is conservative and does not move the rms. Dropping sameline entirely is what pulls 458d in. The sha belongs in the fingerprint; the programme flag is not optional.

410→410b at +0.017 DE-overall is the welded `nmae` failure versus the ballast that unwelded it. Ballast is part of the DE design. An audit that classifies ballast as a bulk knob will report a noise rms an order of magnitude too high.

Zero repair and `de_mag` are post-passes. They are effects (433h's repair, 446bj2z's zs at +0.00055 DE-overall with the profile pinned; zs was scored as nmae −0.00218 raw, and −0.00218 × −1.6532 / 6 = +0.00060). They do not belong in a noise pool. Whether they re-draw the multinomial is a separate question: the seed keys do not change, and zs is a post-pass on a fixed draw. Their board delta is the knob, with no second seed beside it.

Seeds: dropping the seed keys joins **no** new pair. See §4.

---

## 3. Where 0.0045 comes from

`fit230.py:247-252` prints the sd of board Δ`reach` over the arms whose profile relL1 is 0. The recorded value is `fit230.json:262`, `zero_sd` **0.027010964677330575**, `zero_n` 6. RESULTS.md of that experiment states the conversion at lines 241-242: "sd `0.027011` scaled = **`0.004502` on overall**". FINDINGS §8.1b repeats it as the ordering channel (`docs/FINDINGS.md:1077`) and the rule-13 amendment repeats the arithmetic (`docs/FINDINGS.md:1744-1745`: "sd **0.027011 scaled = 0.004502 overall**"). AGENTS rule 13 (`AGENTS.md:384-385`) and the final-phase plan (`docs/FINAL-PHASE-FEEDBACK-PLAN.md:23-24`) both describe that figure as the sd "between digest-identical arms" and, in the plan, as "the emission re-roll".

The six inputs, from `fit230.py:57-68`, cross-checked against the board tables:

| arm | vs | scaled Δreach | emission law, in 230's own words | cross-check |
|---|---|---:|---|---|
| a8 | a2-prune | +0.021398 | overdispersion 0 → 4.0 | `docs/BOARD-LEDGER.md:171` minus `:170` = 0.124349 − 0.102951. Also `docs/LOG.md:1068` |
| a10 | a8 | −0.000754 | spread 0.25→0.15, phi 4→16 | `:172` minus `:171` = 0.123595 − 0.124349 |
| a21 | a15 | −0.018482 | resample x = 0.14 | `docs/LOG.md:2802` |
| a22 | a15 | −0.010656 | resample, the complement | `docs/LOG.md:2865` |
| rc_hybrid | a15 | +0.054940 | \|z\| re-assigned by \|H\| | `docs/BOARD-LEDGER.md:253` minus a15 at `:177` = 0.170308 − 0.115368 |
| rc_callset | rc_hybrid | +0.022590 | call set, n_pred 4480→2599 | `:253`, 0.192898 − 0.170308 |

Cross-check absolute errors on the four ledger subtractions are 0 or < 1e-16 (`work/trace230.json`). a21 and a22 are taken from LOG.md because the ledger prints them already as deltas (`docs/BOARD-LEDGER.md:179` gives a21's reach only to four decimals, −0.0185).

Sample sd, divisor n−1, of those six reach deltas is 0.027010964677330575. Dividing by 6 puts one member onto the overall scale. That division is why the figure is "on overall". It is not the sd of overall. The six overall deltas, transcribed from the same rows (a8 +0.004052, a10 −0.000177, a21 −0.021458 from `docs/LOG.md:2867`, a22 −0.002518 from the same line, rc_hybrid +0.0093325 from `docs/BOARD-LEDGER.md:241`, rc_callset +0.0084448 from `:242`), have sample sd **0.01132**. a21's overall move is five times its reach/6 term because nmae and fid moved with the resample (`docs/LOG.md:2803-2804`).

They are not replicates. Each delta is a treatment effect against a different incumbent. rc_callset is nested inside rc_hybrid, which is nested inside a15. a21 and a22 share a base. A sample sd over that list is a description of how far those laws landed from each other, which is what 230 used it for ("unexplained spread", RESULTS.md:241). It is not a variance component of one law redrawn.

Sensitivity of that sd, still on scaled reach / 6: drop rc_hybrid and it is still **0.00310** (five laws). a10 alone is −0.000177 overall and is the only small one, and a10 still changed phi and the spread. There is no subset of this table that is a re-roll.

230 already knew the scorer itself does not jitter. `docs/BOARD-LEDGER.md:156-160`: a byte-identical resubmission of a15 reproduced every member to float representation, and "Re-seeded twins of the same config are a different matter and do differ." The re-seed was named and not measured. The number that entered rule 13 is the other one.

---

## 4. The re-roll of one design

### 4.1 Board seeds

Across 106 metas, the seed keys take one value each (`work/seed_census.json`):

| key | metas that have it | values |
|---|---:|---|
| `emission_seed` | 99 | 1000003 only |
| `library_pick_seed` | 99 | 0 only |
| `seed` | 99 | 0 only |
| `rc_jitter_seed_base` | 90 | 1000 only |

Groups of metas that match on every other config path and differ in a seed key: **0**. There is no board estimate of a same-design re-roll. Direction members being identical to 16 significant figures (`docs/BOARD-LEDGER.md:256-260`, rc_hybrid and rc_callset vs a15) is a statement about the exact-total repair under a **changed** law, not about two seeds of one law.

### 4.2 H1 seed triplets

Tracked `metrics-*.json` files with a `raw` six-vector and a `seed`, under 435, 517, 528, 480, 490, 467 and 429: 155 files, **40** groups with ≥2 seeds. **480 and 490 contribute none.** 429, 435, 467, 517 and 528 do.

Slopes are the board's, `docs/BOARD-LEDGER.md:162-163` and `experiments/429-zs-production/RESULTS.md:46`: pds 2.2138, mse −1.051, nmae −1.6532, fid 3.3645, reach 1.1292, jac 2.7263. Overall-equivalent of a raw vector is the mean of slope × raw. The sd across seeds of that number is the seed sd of one arm. DE-only uses the four DE slopes and still divides by 6, so a pinned pds and mse do not inflate it.

The production construction, 435 coarse MDE, three seeds (`experiments/435-zs-on-coarse/work/h1/metrics-base_s{0,1,2}_bps_prodem_mde0.1cap0.1min0.3coarse.json`):

| | seed 0 | seed 1 | seed 2 | sample sd |
|---|---:|---:|---:|---:|
| reach raw | 0.269257 | 0.274033 | 0.267435 | 0.003407 |
| fid raw | | | | 0.000085 |
| jac raw | | | | 0.000018 |
| nmae raw | | | | 0.000178 |
| pds raw | 0.895238 | 0.894818 | 0.894538 | 0.000352 |
| mse raw | | | | 0.000082 |
| overall-equivalent (6 members) | | | | **0.000732** |
| DE-only / 6 | | | | **0.000693** |

pds is not seed-invariant on this H1 instrument (a few units in the fourth decimal). On the board, emission-law arms that preserve totals held pds and mse to 16 significant figures. The DE-only figure is the one that corresponds to a board re-roll of this emitter. Reach is most of it: 0.003407 × 1.1292 / 6 = 0.000641.

467's own sentence is "seed sd of one arm's overall-equiv: 0.0004–0.0007" (`experiments/467-strength-gated-programme/RESULTS.md:47`), with reach sd 0.002–0.003 raw. The recomputed production-construction sd sits at the top of that band, slightly above it (0.00073 vs their 0.0007 ceiling). Same object.

Across 36 groups with three seeds, the DE-only sd has minimum 0.000304, median 0.000674, maximum **0.001175** (517's `h1_c4-mzg0.05` hook, reach raw sd 0.00593). The re-roll depends on the design. The production construction is 0.00069. A hook design on the same 85 targets reached 0.00118. Nothing in the tracked triplets reaches 0.0045.

Per-target reach on the 435 triplet (82 of 85 targets finite; 3 are NaN in the csv): the top 5 targets carry **0.927, 0.959 and 0.936** of the absolute per-target seed difference on the three seed pairs (`h1_seeds.json` `per_target_435`). The panel-mean sd (0.00341) matches rms(per-target diff) / sqrt(82) (about 0.0038–0.0046), so the mean really is an average of heavy-tailed target noise, not an extra shared shock. A handful of targets dominate each draw. That does not by itself cancel a 1/sqrt(n) scaling: if a similar fraction of targets is jumpy on a larger panel, the mean's sd still falls as 1/sqrt(n). It does mean the scaling factor is only as good as the assumption that the board's targets are as jumpy as H1's. The board's `n_real` mix is not H1's (`docs/FINDINGS.md:1196-1203`). That factor is unmeasured here.

Scaling the production DE-only sd by sqrt(85/900), under that assumption: 0.000693 × sqrt(85/900) = **0.000213**. By sqrt(85/300), if the three contexts do not add independent draws: **0.000369**. Refusing the credit and keeping the H1 panel sd: **0.000693**.

### 4.3 Shared-draw paired effect, recomputed

467 gated33 minus 435's CLEAN4 base, same three seeds, 429 slopes (`families548.json` `paired_h1`). pds and mse deltas are **0** on every seed. Overall-equivalent of the difference: +0.011088, +0.010971, +0.011215. Mean **+0.01109**, seed sd **0.000122**.

467's published paired sd is a different contrast: gated minus ungated, +0.0043 with seed sd **0.0002** (`RESULTS.md:45`), and the target bootstrap on that contrast is [−0.0034, +0.0132] (`RESULTS.md:9-10`). The emission re-roll of a shared-draw effect is ~0.0001–0.0002 on H1. The uncertainty that actually binds is which targets carry it, and that interval is wider than 0.0045. Those are different questions. Rule 13 and §1.4 are about the emission re-roll.

---

## 5. Verdict table

Sd in overall units. "Use" is the number a bin should be written against. The H1→board factor is stated beside it rather than buried.

| case | sd to use | what it is | what it is not |
|---|---:|---|---|
| (i) bulk step, **current** DE design, seeds held | **0.0004** | Clique sd of DE-overall across 433h, 433k, 446b, 446bj2, 446bd2 is 0.000400. Step rms is 0.000256; ordering-only rms is 0.000202. The sharpen clique is 0.00013. | Not a re-roll. Two of the four steps are identified (`nmae` from ballast K; reach from the LODO load). Not the early family: there, ordering-only rms is 0.00034 but `nmae`/6 rms is 0.0033, and one step is 0.007 overall. Do not apply 0.0004 to `nmae` on a profile-changing arm. |
| (ii) prior-only, same draws (468) | **0** on the board difference | 458g vs 458d is +0.00357 DE-overall with pds and mse flat. 446b vs 446bja is −0.00107 with pds at 0. FINDINGS `docs/FINDINGS.md:420-423`: the blocks match byte for byte, so the difference has no re-roll in it. A new draw of that effect, on H1, has seed sd 0.00012 (this audit, gated33 − CLEAN4) to 0.0002 (467's gated − ungated). | Not 0.0045, and not a reason to ignore the target bootstrap. The effect's composition uncertainty on H1 was several times 0.004. |
| (iii) a knob that re-draws the emission | **0.0007** on an 85-target panel; board value **not measured** | 435 production construction, DE-only seed sd 0.000693 (3 seeds). Other tracked 3-seed designs span 0.00030–0.00118. If per-target noise matches and rows are independent, × sqrt(85/900) = 0.00021. | Not 0.0045. The sqrt(n) figure is an assumption about the board's `n_real` mix and about contexts adding independent draws. No meta pair isolates a seed. |

The 0.009 adoption bar is 2 × 0.0045 (`docs/FINAL-PHASE-FEEDBACK-PLAN.md:24`; `AGENTS.md:432` uses the same sd to say the best of 5 identical arms reads about +0.005). Against case (iii)'s H1 sd, 2 sd is about **0.0014**, and the best of 5 is about **0.0008**, not +0.005. Against the sqrt(n) guess, 2 sd is about **0.0004**. The second of those is 20× below 0.009 and it is the one that is not yet a measurement. Shipping it as the bar would be repeating 230's mistake: promoting a spread computed on the wrong sample.

---

## 6. Replacement wording

Rule 13, the ordering sentence (`AGENTS.md:384-385`). One sentence:

> The figure 0.0045 overall is the sample sd of scaled Δreach/6 across the six digest-identical arms in experiment 230 (a8, a10, a21, a22, rc_hybrid, rc_callset), which are six emission laws rather than re-rolls of one law; a bulk step that keeps the current DE design, including ballast and the seeds, has moved the DE members by about 0.0003 overall on the board, a shared-draw prior arm has no re-roll in its board difference, and a re-drawn emission has seed sd about 0.0007 overall-equivalent on H1's 85 targets (board re-roll unmeasured), so an ordering gain is read against that case and an `nmae` move on a profile-changing arm is the arm rather than this noise.

`docs/FINAL-PHASE-FEEDBACK-PLAN.md` §1.4, replacing the noise bullet (lines 22-25):

> **Noise.** pds and mse reproduce to 16 significant figures between re-builds of one configuration, so a bulk arm is read exactly on those two from a single entry. The 0.0045 figure is the spread of six different emission laws in experiment 230, not the re-roll of the DE design we ship. A bulk step that keeps that design (emitter, sign prior, MDE, zero repair, DE magnitude, call set, ballast, seeds) has moved the four DE members by ≤ ~0.0004 overall on A/B/C; do not apply that tolerance to `nmae` when the arm changes the profile. A prior-only arm that shares draws with its base (468) is read with no re-roll penalty. An arm that re-draws the emission needs an effect > ~0.0014 overall (twice the H1 seed sd of the production construction, 0.0007, on 85 targets) or two seeds; the further reduction to ~0.0002 from averaging 900 rows is not granted until a same-design re-roll is measured on the board or on a panel with the board's `n_real` mix.

§11.1c.4 (`AGENTS.md:432`) still says the best of 5 identical arms reads about +0.005 at sd 0.0045. That sentence inherits the mislabel. The same replacement, if anyone edits it: at seed sd 0.0007 the best of 5 reads about +0.0008, and at the unmeasured board scaling about +0.00025.

The current wording should not stay. It names a design-spread as a re-roll, and the adoption threshold built on it is several times any re-roll this checkout can measure. Replacing it with a single 20× smaller constant would overshoot, for the reasons in §5.

---

## 7. What this checkout cannot settle

- 547's pair list, its rms and its quadrature. The directory is not here, so §2 is an independent reconstruction. If 547's pairs are the modern four, their rms should be 0.000256 and their quadrature 0.000300. If they pooled 222 with 446b, they should have got about 0.0025, and that number is `nmae` on profile arms.
- 537's draft. Same absence.
- A board re-roll. It is not hiding in the metas. Every recorded emission seed is 1000003.
- The H1→900 scaling. Reach noise on 435 is a few targets per draw (top 5 carry ~93% of the absolute mass). Whether the board's smaller `n_real` values make each target jumpier is not in the tracked metric files.
- 480 and 490. No tracked per-seed six-member table in either directory.
- Whether a post-pass (zs, zero repair, MDE) consumes extra random numbers. The seed keys do not change, and the board deltas of those knobs are effects. A second seed of zs on the board does not exist.

## Coordinator's local check (2026-09-30T12:08Z)
- Branch `cloud/548-de-noise-audit` (e90bc4fe, cut from c7db045c: the agent's checkout predated 547 and 537's snapshot, so
  §2 is an independent reconstruction). Taken onto master file by file.
- `audit548.py` re-run on this Mac (0.1 s, tracked inputs only): `families548.json`, `h1_seeds.json`, `pairs548.json`,
  `trace230.json` byte-identical to the agent's; modern-design DE-overall rms 0.000256 reproduced.
- `seed_census.json` differs only by coverage: this checkout has 122 metas (16 newer). Every seed key is still
  single-valued (`emission_seed` 1000003 in 115/115). The committed file is this checkout's.
- It agrees with 547's independent pairs (bulk-only on the current design ~0.0002–0.0004 overall; the `nmae` part
  systematic). Rule 13 / FINAL-PHASE §1.4 are rewritten only after 537's final version (due 10-02 12:00Z) is read against
  §5–§6 here.

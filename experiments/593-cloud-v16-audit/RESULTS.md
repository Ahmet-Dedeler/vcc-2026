# 593 — Adversarial audit of manuscript Draft v16

**Verdict: the new board arithmetic is right, and seven statements are not.** Every scaled delta below was recomputed from `data/board-entries/<id>.json` (`score_*` / 6 for a member's share of overall; raw = the `cell-eval2` member). Dose figures come from `output.rc_callset.prior_lambda` in the build metas. The 198-block bound holds. The October path has no DE cost once the key matches; calling that gap "446b's own realisation" is the reading the entries support, with the limit in attack (b). The 0.00025 re-roll figure matches the seven bulk-only pairs. It does not match a √(85/900) rescaling of the stored H1 panel sd, and no builder-seed pair exists on the board.

Re-run: `uv run python experiments/593-cloud-v16-audit/audit593.py` (log `work/audit593.out`). Line numbers are the v16 files on `origin/master` (`experiments/592-manuscript-v16/paper/`). The paper was not edited.

## Statements that are not confirmed

| location | statement | recomputed | primary | verdict | corrected wording |
|---|---|---|---|---|---|
| `PAPER.md:864`; `CLAIMS.md:210` (C-97) | On H1 the same programme priors were delivered at λ 1.0 on 85 of 85 targets (462, 467, 517). | 83/85 in every programme diag those experiments store. Below 1: MED23 at 0.849 (arcbank) or 0.881–0.883 (arcfull, gated), and PHF12 at 0.970–0.996. The 13 diags that are 85/85 are consensus or non-programme runs (`h1_arc-none-q0.7`, 475's fix, 528's block, 416), not the programme. | `experiments/{462,467,517}-*/work/h1_*/diag-*.json` | **WRONG** | "On H1 the same programme priors were delivered at λ 1.0 on 83 of 85 targets (462, 467, 517). The two exceptions are MED23 (0.85–0.88) and PHF12 (0.97–1.00)." |
| `PAPER.md:1306–1308`; `CLAIMS.md:211` (C-98) | Every wrong call on an in-CCLE line names a documented duplicate or a same-patient line (PC-14, KCI-MOH1, KPL-1, GP2d/GP5d). | Plain cohort, 390 in-CCLE lines: 6 wrong top-1 calls. Five are those duplicates. The sixth is `schraivogel8:K562` → `CORL23_LUNG`, non-unanimous, min top-1 correlation 0.334. The Flex-simulated read adds further non-duplicate misses (HCC1806 → KYSE180; BICR6 → TE15 and HSC3). | `experiments/590-identity-validation/work/match590.csv`, fold CSVs, `flexbias590_fold{0,1}.npz` (cohort as in `analyse590.py`) | **WRONG** | "Of the six wrong top-1 calls on an in-CCLE line, five name a documented duplicate or same-patient line: PC-14/PC-9, KCI-MOH1/HPAC, KPL-1/MCF-7, and GP2d/GP5d (two profiles). The sixth is a non-unanimous K562 → COR-L23 call (min top-1 correlation 0.334)." |
| `PAPER.md:1311–1313`; `CLAIMS.md:211` (C-98) | The lowest top-1 correlation is ≥ 0.68 for 99 % of correct calls. | On the plain cohort the 1st percentile is 0.681, and 98.9 % are ≥ 0.68, only among the 365 correct *unanimous* calls. Across all 384 correct calls the 1st percentile is 0.527 and 97.1 % are ≥ 0.68. Flex-simulated correct-unanimous q01 is 0.678. | same CSVs; `summary590.json` records the unanimous quantile | **OVERSTATED** | "It is ≥ 0.68 for 99 % of correct unanimous calls (1st percentile 0.681 on the plain read; 0.678 after the Flex-minus-3′ bias). Across every correct call, including non-unanimous ones, 97 % are ≥ 0.68 and the 1st percentile is 0.53." |
| `PAPER.md:1315–1316`; `CLAIMS.md:211` (C-98) | A floor of 0.65 removes every false switch at the cost of 2 of 19 correct switches (a partial gene axis; sci-RNA-seq3 or ECCITE). | Plain: 19 correct registry switches, 17 still on at 0.65. The two losses are Nadig HepG2 (10x 3′, 0.540) and sci-RNA-seq3 A549 (0.636). Flex-simulated: 18 correct registry switches, 16 still on. The two losses are the same HepG2 profile (0.396) and ECCITE THP-1 (0.636). ECCITE is not one of the plain 19. | same | **OVERSTATED** | "On the plain read the floor removes both false switches and costs 2 of 19 correct registry switches (Nadig HepG2, partial axis, 0.540; sci-RNA-seq3 A549, 0.636). On the Flex-simulated read the cost is 2 of 18 (that HepG2 profile at 0.396, and ECCITE THP-1 at 0.636)." |
| `CLAIMS.md:209` (C-96); the same gloss is looser in `PAPER.md:952` | H1's re-roll scaled by √(85/900) predicts 0.00123 raw for reach. | The stored prediction is 0.001230 and equals √2 × 0.026098 / √900 (`per_target_sd` on 82 finite targets). `h1_panel_sd` × √(85/900) = 0.001047. | `experiments/537-arm-pipeline/work/noise537b.json` key `h1.base_reroll_independent.1.0.reach`; formula in `noise537.py` | **OVERSTATED** | "H1's per-target re-roll, as a paired difference over the board's 900 rows (√2 × 0.02610 / √900), predicts 0.00123 raw for reach." |
| `PAPER.md:949` | A bulk change re-draws pass 1 on 83–84 of 85 targets, with per-target reach uncorrelated with the base. | Pass 1 stays equal on 1, 2 and 1 of 85 targets (seeds 0, 1, 2), so 84/83/84 are re-drawn. Same-seed correlation of per-target reach is −0.211; nmae +0.146; fid +0.421; jac +0.597. | `experiments/537-arm-pipeline/work/redraw537.json`, bulk-only step | **OVERSTATED** | "A bulk change re-draws pass 1 on 83–84 of 85 targets. Per-target reach correlates −0.21 with the base (nmae +0.15). Fid and jac stay correlated at +0.42 and +0.60, because the call set is prescribed." |
| `PAPER.md:1091` | Both entries delivered the programme well under full dose. | 458d's pooled median is 0.4331. 458g's pooled median is 0.9526. The under-half dose on 458g is the strong-block bound (attack a), not the pooled median. | the two build metas | **OVERSTATED** | "458d's pooled median dose is 0.43. 458g's pooled median is 0.95, lifted by its 252 consensus blocks; at least 198 of the 648 programme blocks it shares with 458d are at λ ≤ 0.433 (C-97)." |
| `PAPER.md:863`; `CLAIMS.md:210` (C-97) | Pooled median λ is 1.0 for 446b and every consensus-prior entry. | 54 build metas have no `programme_default`. 53 have median 1.0. The exception is 392: prior `clean4`, emitter `a12cf8f9` (not the current emitter `edb5b8dc`), median 0.901, min 0.087. On the current emitter every non-programme meta has median 1.0. Those medians are not "every block": 446b's own minimum is 0. | `submissions/*.h5ad.meta.json` | **OVERSTATED** | "Pooled median λ is 1.0 for 446b and for every other non-programme build on the current emitter. One earlier consensus build (392, clean4) has median 0.901. A median of 1.0 still leaves individual blocks clipped: 446b's minimum is 0." |

## Confirmed (recomputed)

| location | statement | recomputed value |
|---|---|---|
| `PAPER.md:7–8,467–472`; C-94 | 460dr − 453dr overall +0.0028694; mse −0.013060 raw = 0.902× the key's −0.014473; pds +0.001044 raw. 460dr − 446bj2: DE −0.0000139 overall, bulk −0.0030321. Bands \|Δfid\| 2.7e−6, \|Δjac\| 5.4e−6, \|Δreach\| 0.00122, all inside 0.0001 / 0.00005 / 0.0040. | `LrkhsIXWMjEsvRXSGToa`, `q3VuTImh4FRGKXH2waCp`, `YvVru3DAF89Jiyr34olR`, `UDgjOKL9gI1udAexq6ue` |
| `PAPER.md:457–466`; C-82 | 453dr − 446b overall −0.0036816 (bulk −0.0030866, DE −0.0005950). 457dr − 453dr overall +0.0031484 (bulk +0.0031006; 0.0031484/0.00172 = 1.83). 457dr − 446b overall −0.0005332 (bulk +0.0000141, DE −0.0005473). DE vs 446b: 446bj2 −0.000388, 446bd2 −0.000481, 453dr −0.000595, 457dr −0.000547, and 460dr −0.000402. | same entries plus `R5uE5Rve3R9nWXGbeSLM`, `1FfirZsGsAZvaXYTd39J` |
| `PAPER.md:860–871`; C-89, C-97 | 458d − 446b overall −0.008514 (fid −0.005545 raw). 458g − 458d +0.003566. 458g − 446b −0.004948. Medians: 446b 1.0, 458 0.2427, 458d 0.43310546875 = 887/2048, 458g 0.9526. Strong/weak 214+217+217 = 648/252. Programme down-fraction on 458d context A is 0.568. | metas; `experiments/468-prior-gate/RESULTS.md` |
| `PAPER.md:868–871` | 458g's programme blocks equal 458d's. | 468: CSR data and indices, 0 of 900 blocks equal neither. All 648 strong blocks equal 458d. One strong block in A is also equal to 446b. |
| `PAPER.md:940–955`; C-83, C-96 | Reach rms of the seven pairs 0.001295, χ² 95 % [0.000860, 0.00264]. Fid rms 2.41e−5, jac rms 1.30e−5, DE-overall rms 0.000268, reach's share rms 0.000242. 230's six emission-law pairs: sample sd of scaled reach 0.027010964677331, /6 = 0.004501827446222, matching `fit230.json`. | board JSON; `noise537b.json`; `fit230.json` |
| `PAPER.md:1044–1057`; C-95; T12 | 446bj2zcp − 446bj2z: overall −0.000794; Δ(fid+jac) = −0.0001765 overall; jac +0.003986 raw (0.438× H1's +0.0091); fid −0.003580 (2.75× H1's −0.0013); reach −0.003353; pds 0; mse −2.32e−8. Swapped-in sign accuracy 0.508. A/B/C swap sizes 546/572/684. | `K1Q0y2OFYjTNRqh9jKfk`, `be90smZUDcIV7j5MgQ03`; `528/RESULTS.md` |
| `PAPER.md:1298–1318`; C-98, the counts | 402 profiles, 390 in CCLE, top-1 0.9846, accuracy given unanimity 0.9946, sister rows 146 with 0 false switches, false switches 2 (HEK293 0.554 and 0.606; Flex-simulated 0.508 and 0.561), A/B/C 0.772/0.719/0.733, primary T cells 0.719, official ∩ CCLE axis 17,900 genes. | the 590 CSVs and fold masks |
| `PAPER.md:119,132–133`; `figs/MANIFEST.md:73` | "Under half" is 458d's median 0.433. F-12's −0.0037 / −0.0031 / +0.0031 / −0.0005 are the rounded board deltas above. | metas; board JSON |
| `REFERENCES.md:245–246` | The two new keys are CELLOSAURUS and CCLE-2025. The duplicate IDs the reference names are the five calls the CSV supports. The reference does not mention COR-L23. | the CSV, as above |

## Attack (a) — the bound of 198

The arithmetic holds, and it uses nothing beyond the pooled median plus 468's split.

458d records min 0, median 0.43310546875, max 1, over 900 blocks. That median is exactly 887/2048. A 10-step bisection lives on the grid k/1024, and 887/2048 is the average of 443/1024 and 444/1024. For an even sample the median is the average of the 450th and 451st order statistics, so those two order statistics are 0.432617 and 0.433594. At least 450 blocks therefore have λ ≤ 0.433.

468's gate puts 86 + 83 + 83 = 252 blocks on the weak third and 648 on the programme. If every weak block sits in the lower 450, 450 − 252 = 198 programme blocks remain there. That is a lower bound. The measured count could be anything from 198 up to 450; the metas do not store per-block λ, and the paper says so in C-97's uncertainty column. No further assumption is required.

468 restored the three packages and compared all 900 blocks on CSR data and indices. Every 458g block equals 458d (strong targets) or 446b (weak targets), and none equals neither. The programme blocks are the same counts in both scored files, so the same delivered dose is in both. "Byte-identical" here means that CSR equality, which is what 468 recorded. One strong block in context A is identical in all three packages; it still equals 458d.

458g's own pooled median, 0.953, is a different statistic. It is lifted by the 252 consensus blocks and does not say the programme blocks were delivered near full dose.

## Attack (b) — is the DE gap 446b's draw, or a cost of the path?

DE share of overall, recomputed:

| pair | DE overall | what changed |
|---|---:|---|
| 417a − 410b | +0.000182 | bulk, earlier in the chain |
| 417e − 417a | −0.000031 | bulk |
| 433k − 433h | −0.000184 | bulk |
| 446b − 433k | −0.000233 | bulk; 446b sits *below* its parent |
| 446bj2 − 446b | −0.000388 | LODO key, tuned path |
| 446bd2 − 446b | −0.000481 | fifth donor, tuned path |
| 453dr − 446b | −0.000595 | October path, CD4 ×1, no day-of key |
| 457dr − 453dr | +0.000048 | CD4 ×1.5 on that path |
| 457dr − 446b | −0.000547 | the net the previous draft charged to the path |
| 460dr − 453dr | +0.000193 | day-of key on the October path |
| 460dr − 446bj2 | −0.000014 | October path vs tuned path, key on both sides |
| 460dr − 446b | −0.000402 | October path with the key, against 446b |

A systematic DE cost of the October path would show up in 460dr − 446bj2. That difference is −0.000014 overall, about a twentieth of the 0.00025 re-roll figure, and the whole of the overall gap (−0.003046) is bulk (−0.003032). The entries do not support a path cost.

They do support a gap that is attached to 446b. The three builds that leave 446b by a bulk change only — the LODO key, the fifth donor, and the October path — all land 0.00039 to 0.00060 below it, and the later October-path re-rolls (457dr, 460dr) stay in that band. Earlier bulk steps in the same emission law move DE in both directions, so the one-sidedness is not "any bulk change costs DE".

The limit: this is three bulk departures from one file, not a seed pair, and 457dr and 460dr are re-rolls of the October family rather than new draws against 446b. 446b itself is 0.00023 below 433k, so the draw is high relative to its descendants, not relative to the whole chain. "Most likely 446b's own realisation" is the reading that survives the path-cost alternative. It is not a demonstration that no later code change is shared by those descendants.

## Attack (c) — the 0.00025 figure

The seven bulk-only pairs in 537 §2, recomputed from the board JSON, give reach rms 0.00130 raw with the χ² interval [0.00086, 0.00264], reach's contribution to overall rms 0.00024, and the four DE members together rms 0.00027. "About 0.00025 overall" is that pair of figures. 230's 0.004502 overall is the sample sd of six pairs that each changed the emission law (overdispersion, resampling, the |z| assignment, the call set), divided by 6. The revision note is right to separate them.

Whether those seven pairs are re-rolls: 537 §2.6 and 547's correction say yes, because one emission stream serves the file and a moved ballast count shifts every later block. The H1 test of that mechanism re-draws pass 1 on 83–84 of 85 targets, including targets whose own K did not move. v16 follows that correction. 548's earlier conclusion that the same pairs are knob leakage rather than a re-roll is the sentence 547 withdrew. 548's other caveat still stands and v16 does not break it: the earlier family (222 through 319p2) has DE-overall rms 0.00321, almost all of it nmae, and v16 scopes 0.00025 to the law in use since 410b, with an explicit ballast term in nmae.

Two sentences are stronger than those records.

- C-96's "scaled by √(85/900)" does not produce 0.00123 from the stored H1 panel sd (that product is 0.001047). The stored 0.001230 is √2 × (per-target sd 0.026098) / √900, on 82 finite targets. §3(k)'s "scaled from 85 to 900 target rows, predicts 0.00123" names the right number and the wrong count of finite targets.
- "Uncorrelated" overstates reach ρ = −0.21 and skips fid and jac, which stay correlated because the call set is not drawn. 548's count is also still true: 99 of 99 metas use emission seed 1000003, so no seed-only board pair exists. v16's own falsifier for C-96 asks for that pair. The 0.00025 figure is the bulk-change re-roll 537 measured, not a seed re-roll.

The eighth pair named in §3(k), 460dr − 453dr, has raw deltas fid +2.7e−6, jac +5.4e−6, reach +0.00122, all inside the preregistered 3σ bands. That check is real.

## Unchanged text

No unchanged sentence still prices a re-roll of this emission law at 0.0045. The remaining 0.0045 figures are the preregistered composition band on 453dr (−0.0045 … −0.0005) and the revision note that identifies the old number as 230's law-change spread. C-93 still says the key's DE members sit "inside the ordering noise" without restating 0.0045; the number that now belongs there is 0.00025. The dose caveat was written into the abstract, the plain-language summary, §4.1, §4.7 and F-12. §4.1's "both entries … well under full dose" is the sentence that overshoots §3(j)'s own median for 458g; it is in the table above.

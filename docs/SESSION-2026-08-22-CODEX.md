# VCC research continuation — 2026-08-22

## Starting state

- Live leaderboard: rank 3 of 128, score `0.1377129868`.
- Current board artifact: v8, model name `We regret to inform the genes`.
- Staged candidate: v9, union-axis whitening at gamma `0.65`, packaged and
  previously validated by `vcc prep --dry-run`; it has not been submitted.
- Submission policy: do not submit without Ahmet's explicit authorization.

## Memory-safety policy

The prior local run exhausted system memory. This continuation therefore uses:

- no parallel local `cell-eval2` jobs;
- a 4 GiB process-tree guard for substantial local Python work;
- 1 GiB guards for synthetic tests;
- cloud workers for full eigendecompositions and count-level scoring;
- periodic `memory_pressure -Q` and process inspection.

Observed local free-memory readings stayed between 85% and 89%. No local VCC
process exceeded 1 GiB during this continuation.

## Completed

### Experiment 047 — primary CD4 is a real additional donor

The compact primary-CD4 extraction was recovered from AWS and verified:

- 291/300 panel targets in the condition union;
- 10,113 mapped genes;
- extraction peak RSS 0.34 GiB;
- compact archive SHA-256
  `ed694e2d6c978d88a2a5065bc4ce2ca075604f5421c165e686702133c4354e61`.

The frozen direct-transfer rule selected on HCT116 + HEK293T was a
gene-normalized equal mean of the three CD4 conditions at fraction `0.50`.
It improved the two-fold development mean by `+0.03085`, the worst development
fold by `+0.01867`, and the K562 target-role holdout by `+0.01970`.

This is not a fully independent K562 donor holdout because K562 appears as a
source during development LODO. It is enough to advance CD4 to a
production-union test, not enough to ship it directly.

### Experiment 051 — per-perturbation CD4 breadth weighting fails

The full sweep completed locally under guard:

- 242 shared targets;
- 17,942 non-panel union genes;
- internal peak RSS 2.82 GiB;
- synthetic tests: 4 passed.

The selected adaptive rule helped both development folds but failed on K562:

- adaptive versus matched fixed fraction: `-0.00279` PDS;
- adaptive versus shuffled-breadth mean: `+0.00037`;
- shuffled p95: `0.77483`, above adaptive `0.76993`;
- one-sided empirical permutation `p = 0.50`.

Decision: reject response breadth/QC as a per-target blend-weight signal.
Retain the fixed `0.50` three-condition CD4 blend.

### New data scan

The most actionable newly identified source is the KOLF2.1J genome-scale
CRISPRi atlas:

- 10,167 perturbations x 38,606 genes;
- 273/300 VCC targets;
- public browser matrix supports byte-range extraction of only the needed
  rows, avoiding a full large download;
- experiment 055 is extracting and testing it under a 4 GiB cap.

Other findings:

- PerturbMatch provides a new 748 MB K562-GWPS LFC reanalysis with 251/300
  panel targets and is the next donor-quality experiment.
- VIPerturb-seq provides an independent K562 CRISPRi measurement, but only 124
  panel targets survive its filtered object; it is lower priority than KOLF.
- Broad's 16-line DepMap pilot contains only 90 CRISPR knockout targets and
  currently requires Taiga credentials; it is a metadata probe, not tonight's
  primary acquisition.
- The IGVF production TF atlas spans four differentiated lineages and its
  1,983-target library overlaps 91/300 VCC targets, but the current GCS/Synapse
  artifacts were not anonymously accessible.

## Running

### Experiment 040 — clean HCT gamma sweep

Four production-union candidates (`0.55`, `0.65`, `0.75`, `0.85`) are built
on AWS. `cell-eval2` is scoring them serially. This replaces the quarantined
H1 comparison that excluded the wrong target-coordinate panel.

### Experiment 052 — source subspace versus full union whitening

Complete and negative. The fixed diagnostic `subspace_then_union` composition
was `+0.00503` macro PDS over union-only, but the leak-free nested selector
could not identify it. Under both split seeds, selected macro PDS was
`0.77148`, or `-0.01399` versus union-only, with `0/3` positive outer folds.
The predeclared gate failed, so v9's geometry remains unchanged.

The corrected synthetic suite has 10 passing tests. The complete cloud run
peaked at 3.97 GiB in-process; repeat it only on cloud.

### Experiment 054 — fixed CD4 blend under the clean HCT union whitener

Complete and negative. On the real 132-target HCT panel:

- union-whitened K562 + HEK baseline: `0.80095` PDS;
- frozen gene-normalized CD4 blend: `0.76544` (`-0.03551`);
- raw CD4 blend: `0.73791` (`-0.06303`).

The loss persists on CD4-covered and all-three-condition target subsets.
Production self-renormalization and fixed-baseline-scale diagnostics agree.
CD4 helps before whitening but does not compose with v9's geometry; do not add
it to the submission.

### Experiment 055 — compact KOLF2.1J donor

Complete and weakly positive. Strict range extraction fetched only the 273
needed panel rows (10,539,438 bytes) and produced a 30 MB compact donor over
18,106 mapped VCC genes.

The rule selected on HCT116 + HEK293T only was raw KOLF at fraction `0.10`:

- development mean delta: `+0.00049` PDS;
- worst development fold: `-0.00115`;
- frozen K562 target-role holdout: `+0.00399`;
- three-fold macro delta: `+0.00166`.

The gain passes the predeclared direct-transfer gate but is small and has one
negative development fold. Experiment 057 is testing whether it survives
production-style union whitening on the frozen K562 role. Do not add KOLF to
v9 from the direct result alone.

### Experiment 057 — KOLF fails under union whitening

Complete and negative on the frozen K562 target-role holdout:

- no whitening: `+0.00387` PDS, reproducing experiment 055;
- gamma `0.55`: `-0.00805`;
- fixed production gamma `0.65`: `-0.00687`;
- gamma `0.75`: `-0.00685`;
- gamma `0.85`: `-0.00696`.

The operator used HCT116 + HEK293T rows only, with zero K562/KOLF fit rows and
all 300 panel coordinates removed. Decision: do not add KOLF to the
union-whitened model.

### Experiment 058 — PerturbMatch processing ensemble

Complete and promising. PerturbMatch is a different analysis of the same
Replogle K562 assay, not an independent biological donor. Replacing the
incumbent K562 table loses, but averaging the two K562 estimates 50:50 before
combining with the other source donor improves source-whitened PDS on both
valid Orion holdouts:

- HCT116: `+0.00488`;
- HEK293T: `+0.00546`;
- macro: `+0.00517`;
- paired target-cluster bootstrap interval:
  `[-0.00076, +0.01131]`.

The result uses a complete 238-target retrieval gallery, 6,949 shared
non-panel genes, and a candidate-independent gamma-0.65 operator. It is not
yet production-approved because the confidence interval crosses zero.
Experiment 061 is the full HCT116 count-level confirmation.

### Experiment 059 — Feng 34-line iPSC atlas

The official pooled effect table was found, checksum-verified, and streamed
into a compact local matrix:

- 6,673 source targets x 6,471 biological genes;
- 182/300 VCC targets;
- 6,151/18,533 VCC genes;
- extraction peak RSS `0.14 GiB`.

The best predeclared raw-space candidate, a 2% blend restricted to targets
with significant on-target knockdown, moved three-fold macro PDS by only
`+0.000070`—approximately a handful of rank steps and below practical
resolution. Larger fractions were harmful. Retain the data artifact, but do
not add Feng to v9/v10 or spend cloud scorer time on it.

### Experiment 041 — context identity

Cross-platform CCLE matching, validated by recovering known K562, HCT116, and
Jurkat controls, gives stable top matches under all eight methods:

- context A: Jurkat;
- context B: HeLa;
- context C: CAL33 / head-and-neck squamous.

This is being used only to motivate a predeclared proxy test. Experiment 060
tests whether the fixed primary-CD4 donor helps a held K562 hematopoietic
target under union whitening; no challenge output is changed directly from
the identity match.

### Experiment 060 — context-A CD4 routing rejected

The frozen held-K562 test is complete and strongly negative after source-only
union whitening:

- HCT116 + HEK293T baseline PDS: `0.810157`;
- gene-normalized CD4 mean3 at fraction `0.50`: `0.667072`
  (`-0.143085`);
- raw-CD4 diagnostic: `0.717191` (`-0.092967`).

Before whitening, the primary CD4 arm was positive (`+0.022119`), so this is a
specific non-composition result rather than evidence that the CD4 data are
poor. The production operator amplifies CD4 directions that do not transfer
to held K562. Do not route context A to CD4 and do not add CD4 to v9/v10.

### Experiment 062 — post-whitening CD4 residual works

Changing the order of operations rescued the CD4 signal. Whitening only the
ordinary donor consensus and then adding a small unwhitened CD4 residual gave:

| held line | delta PDS |
|---|---:|
| HCT116 | `+0.004818` |
| HEK293T | `+0.002418` |
| frozen K562 | `+0.004407` |

The fixed rule is CD4 mean3, gene-normalized, fraction `0.05`. It was selected
by worst development-fold gain before opening K562. The K562 paired bootstrap
95% interval was `[+0.001286, +0.007613]`. This is a real effect-level
candidate, but still needs count-level validation.

### Experiment 063 — residual de-whitening cannot recover MSE

A frozen blend from the production union-whitened direction back toward the
raw donor consensus was monotonic in the wrong score tradeoff. Raw share
`0 -> 1` changed macro PDS `0.794219 -> 0.705406`, while the best-case
nonnegative-alpha MSE floor improved only `0.998385 -> 0.995052`. Every point
remains above the official `0.989` MSE baseline and therefore earns zero
scaled MSE. Keep the production direction unchanged.

### Experiment 064 — PerturbMatch and post-whitening CD4 compose

The two frozen changes were tested as a no-retuning 2x2 factorial on held
HCT116 and HEK293T:

| arm | macro PDS | delta |
|---|---:|---:|
| incumbent | `0.804309` | — |
| PerturbMatch only | `0.808896` | `+0.004587` |
| CD4 only | `0.808063` | `+0.003754` |
| both | **`0.812875`** | **`+0.008566`** |

Both folds improved, the interaction was slightly positive (`+0.000226`), and
the paired target-cluster bootstrap interval for both versus incumbent was
`[+0.001532, +0.015724]`. Experiment 065 is prepared to test this exact joint
candidate at HCT116 count level after experiment 061 finishes.

### Experiment 066 — per-target arm routing rejected

A leakage-safe nested selector tested whether inference-time features could
choose among experiment 064's four frozen arms per perturbation. The best gate
only learned when to omit PerturbMatch while retaining CD4:

- fixed combined arm D: `0.812875360` macro PDS;
- nested gated result: `0.814191691`;
- gain: `+0.001316331`;
- paired bootstrap interval: `[-0.001912793, +0.004781983]`.

The gain is below the `+0.003` practical floor, its interval crosses zero, and
the selected feature changed between donor folds. A same-target lookup trained
on the other donor lost `-0.006735911`, confirming weak cross-donor stability.
Keep the fixed combined arm D; reopen routing only with a third independent
somatic holdout.

The authoritative run peaked at `0.213 GiB`. Two preliminary read-only mmap
probes reached about `2.60 GiB` before the reader was replaced with bounded
streaming; they exited and wrote nothing.

### Experiment 068 — seed-off-target row ablation rejected

Public seed-effect tables identified four panel targets with high-confidence
source-guide off-target candidates: `PQBP1`, `SHPRH`, `FZD2`, and `GLYR1`.
Surgically deleting the implicated source donor row was strongly harmful:

- four-target mean PDS delta: `-0.135226360`;
- native experiment-064 case mean: `-0.143344710`;
- improved targets: `0/4`;
- matched-null `P(null >= observed)`: `0.933603`.

`FZD2` improved on HEK293T but reversed strongly on independent KOLF truth.
Whole-row deletion discards useful perturbation identity signal along with the
suspected off-target component. Keep arm D unchanged. The successful run
peaked at `0.608 GiB`; its guarded preliminary attempt stopped below the
`0.75 GiB` cap.

### Experiment 067 — context-control variance transport rejected

The public `sd_target_control / sd_source_control` idea was tested as a
diagonal effect-coordinate transport using disjoint control-cell replicates,
robust SD/MAD scales, clipped ratios, and placement before or after the frozen
source-only union whitener.

The leakage-safe source-only selector was inconsistent:

- K562 selected MAD post-whitening at eta `0.5` and gained about `+0.002`;
- HCT116 selected identity;
- HEK293T selected SD pre-whitening at eta `0.25` and lost about `-0.003`;
- replicate-A macro delta: `-0.000249419`;
- replicate-B macro delta: `-0.000394082`;
- no selected nonidentity fold beat its coordinate-permutation null.

The optimistic fixed-grid winner, MAD post-whitening eta `0.5`, gained only
`+0.001008` on replicate A and `+0.000668` on replicate B, below the
predeclared `+0.002` gate. Reject variance transport.

Successful final runs stayed at or below `0.59 GiB` process-tree RSS. An
abandoned prototype briefly reached `0.81 GiB` before its guard killed it; no
artifact from that attempt was used.

### Experiment 040 — clean HCT gamma-0.65 baseline completed

The production-shaped 132-target HCT116 count harness completed successfully
on 2026-08-23:

- raw PDS: `0.7978834143`;
- official-anchor scaled PDS: `0.6583058879`;
- raw MSE: `2.0537529303`, scaled `0`;
- scaled six-metric average: `0.1216402396`;
- wall time: `3:27:46`;
- peak scorer RSS: about `17.0 GiB`;
- exit code: `0`.

Compact outputs and hashes are preserved in
`experiments/040-union-gamma/cloud-output/g0.65_q1.0_a1.08_spread0.25/`.
This is now the fixed comparator for experiments 053, 061, and 065.

The first experiment-053 watcher attempt failed before heavy work because
systemd did not inherit the path to `uv`. The watcher was patched, restaged,
and relaunched. Its 133-block template transform completed under `0.5 GiB`
RSS and the template scorer started.

Experiment 061 initially received S3 `403` while trying to read the clean
baseline. A temporary read-only permission for exactly `output-clean-hct/*`
was added to the worker role, the baseline was copied, scoring started, and
the permission was immediately removed. No candidate rebuild was needed.

### Live state checkpoint

At 2026-08-23 00:30 UTC the public leaderboard had 130 teams. Ahmet remained
rank 3 at `0.1377129868`; Ibrahim Mansour led at `0.1447796056` and VIVAI was
second at `0.1411686718`.

The staged v9 package remains unsubmitted:

- `submissions/v9.vcc`: 3,666,636,800 bytes;
- SHA-256:
  `3a35d3981360fe2f8d6ab90754722c592af9d67cf5728dfe542ad43c2891efe8`.

## Current decision

Keep v9 staged and unchanged while the paired count-level experiments finish.
The clean HCT comparator is now available. CD4 and subspace geometry are closed
when CD4 is blended before whitening; a small post-whitening CD4 residual is
positive and composes with PerturbMatch at effect level. KOLF, Feng, and
residual de-whitening for MSE are closed. Template emission and PerturbMatch
alone are actively scoring; the combined PerturbMatch + CD4 candidate remains
gated on experiment 061. None of these results is submission authority. No
challenge submission has been made in this continuation.

### 2026-08-23 04:15 UTC continuation

Both outstanding paired HCT116 scores completed successfully:

- Experiment 053 template emission is rejected. PDS was exactly unchanged,
  but fidelity collapsed enough to lower the scaled six-metric average from
  `0.1216402396` to `0.0509509792` (`-0.0706892605`). Keep quasi emission.
- Experiment 061 PerturbMatch is positive but practically flat. Raw PDS moved
  `0.7978834143 -> 0.7981725653` (`+0.0002891511`), and scaled six-metric
  average moved `0.1216402396 -> 0.1217308225` (`+0.0000905829`). It passes
  the predeclared non-negative prerequisite but is not worth shipping alone.

Experiment 065, the frozen PerturbMatch plus post-whitening CD4 candidate, was
launched serially on the restarted experiment-061 worker:

- run id: `20260823T040540Z`;
- worker: `i-0343d04439564b709`;
- S3 output:
  `s3://<your-s3-bucket>/output/exp065-pm-cd4-hct-count/20260823T040540Z/`.

The build completed and its single `cell-eval2` scorer started at approximately
04:11 UTC. No challenge submission was made.

The freed template worker is running experiment 069, a gated effect-level test
of fitting the covariance to the same frozen PerturbMatch K562 ensemble used
by the prediction:

- run id: `20260823T041514Z`;
- worker: `i-08d90ae88c7b210ed`;
- S3 output:
  `s3://<your-s3-bucket>/output/exp069-perturbmatch-aligned-whitener/20260823T041514Z/`.

The first exp069 dispatch failed before tests or scientific computation because
the source archive omitted the repository `README.md` required by
`pyproject.toml`. The launcher was fixed and the clean run above was started.

Experiment 069 then completed successfully and strongly rejected the
PerturbMatch-aligned covariance:

- incumbent arm-D macro PDS: `0.7913002211`;
- aligned-whitener macro PDS: `0.7826016436`;
- delta: `-0.0086985775`;
- HCT116 delta: `-0.0079777009`;
- HEK293T delta: `-0.0094194541`;
- paired target-cluster bootstrap interval:
  `[-0.0133902826, -0.0041270185]`;
- peak process-tree RSS: `8.862 GiB`;
- exit code: `0`.

Do not put PerturbMatch into covariance fitting and do not run a count-level
confirmation of experiment 069.

A separate production adapter now exists at
`scripts/make_submission_pm_cd4.py`. It preserves v9's three-donor weighting,
operator, alpha, CPM gate, quasi emission, and split RNG streams while adding
the frozen PerturbMatch replacement and post-whitening CD4 residual. It writes
atomically and refuses to overwrite v9. Eight synthetic tests pass at
approximately `0.13 GiB`; no full challenge build has started.

### 2026-08-23 05:50 UTC continuation

Experiment 072 completed its full Jiang 2025 pathway-DE archive audit under a
hard 1 GiB local process-tree limit. The first implementation used memory maps
for sequential scans of the 2.6 GB Orion donor matrices; macOS retained mapped
pages and tripped the experiment's 0.90 GiB internal ceiling. Replacing that
path with the existing `os.pread` row reader reduced the completed run to
0.56 GiB peak process-tree RSS.

The archive contains 271 DE tables, 218 targets, six cell lines, and 9,540
mapped VCC genes, but only nine VCC panel targets. The predeclared minimum was
30, so held-panel outcomes were never opened and Jiang is rejected as a
production donor. On 171 non-panel development targets, a raw 0.20 blend was
positive in all three held-donor folds (`+0.006318` macro PDS), but its
target-cluster bootstrap interval crossed zero
(`[-0.002603, +0.015079]`). This is suggestive biology, not submission
evidence.

Experiment 065 remained healthy at 05:41 UTC: the scorer was active at 100%
CPU, 7.53 GiB peak RSS, and 53 GiB available memory on
`i-0343d04439564b709`. No final aggregate existed yet.

The public leaderboard had 138 teams at 05:42 UTC. Ahmet was rank 4 with the
unchanged v8 score `0.1377129868` and remained the raw-PDS leader at
`0.8163285`. The staged v9 package was unchanged and unsubmitted.

### 2026-08-23 06:36 UTC continuation

Experiment 073, post-whitening CSLS-density distillation, completed
successfully and was decisively rejected:

- baseline macro PDS `0.7903534764`;
- nested selected macro PDS `0.7793740584`;
- delta `-0.0109794180`;
- target-cluster bootstrap 95% interval
  `[-0.0193550278, -0.0024592699]`;
- exit code `0`, peak process-tree RSS `2.708 GiB`.

Experiment 075, query-only reciprocal local centering, completed after two
pre-compute integration fixes. The frozen candidate `k=64, beta=0.5`
improved all three outer folds:

- K562-GWPS `+0.0013019664`;
- Orion HCT116 `+0.0040555506`;
- Orion HEK293T `+0.0034719104`;
- macro `+0.0029431425`;
- bootstrap 95% interval `[+0.0010475592, +0.0048985863]`.

It failed only the predeclared practical gate `macro >= +0.003`, by
`0.0000568575`. Do not promote it directly. A frozen independent KOLF2.1J
confirmation is being implemented as experiment 079.

Experiment 076's public-data scout found one immediately usable new context:
the released IGVF WTC11-derived cardiomyocyte CRISPRi filtered global-DE
table. The 1.89 MB artifact contains 87/300 panel targets; its guide manifest
contains 91/300. CellClip and PRiMeFlow add no new qualifying donor, and
X-Atlas/Pisces still has no public payload. A leakage-safe sparse-DE donor
test is being implemented as experiment 081.

Experiment 070's corrected VIPerturb run failed before `readRDS()` because the
source manifest repeats the single `NO-TARGET` control once per source bin.
The parser now removes recognized controls before enforcing uniqueness on the
6,724 real perturbations. Fifteen synthetic tests and three cloud R schema
tests pass. Corrected run `20260823T063300Z` is active on
`i-08d90ae88c7b210ed`.

Experiment 077's low-alpha `0.35` and `0.40` count-level arms are actively
scoring concurrently on isolated Spot worker `i-0e4070204e83188fb`.
Experiment 065 remains active and healthy on `i-0343d04439564b709`.

The production challenge controls and exact v9 union whitener are now staged
under `input/production-challenge/`. Local and staged source hashes are:

- whitener NPZ
  `5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3`;
- whitener JSON
  `a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6`;
- staged v9 package
  `3a35d3981360fe2f8d6ab90754722c592af9d67cf5728dfe542ad43c2891efe8`.

No challenge submission was made. The Mac remained free of heavy VCC jobs and
reported 90% system-wide memory free.

### 2026-08-23 07:05 UTC continuation

Experiment 070's corrected evaluation-only retry completed successfully from
the preserved VIPerturb compact without rereading the RDS. VIPerturb is
decisively rejected:

- Orion HCT116 PDS delta `-0.0055267205`;
- Orion HEK293T PDS delta `-0.0083982122`;
- macro delta `-0.0069624664`;
- paired-target bootstrap 95% interval
  `[-0.0108373282, -0.0034542003]`;
- resetting all modified low-cell targets still left a negative
  `-0.0046015955` delta.

All promotion checks failed. No count-level confirmation or production
integration is allowed. The retry exited `0`, peaked at `4.5166 GiB`, and all
downloaded artifacts verified against their manifest under
`experiments/070-viperturb-k562/cloud-output/20260823T063300Z/evaluation-retry/20260823T070153Z/`.
The now-idle on-demand worker `i-08d90ae88c7b210ed` was stopped after
preservation.

An independent review of experiment 080 found multiple protocol blockers,
including the wrong final retrieval axis, target-coordinate mixing through
the dense operator, an incomplete permutation-selection null, partial
galleries during tune/pretest, and post-gate refitting. Experiment 080 was not
launched; the implementation was returned for correction and expanded
regression tests.

Experiment 065 and both experiment-077 scorers remain healthy on AWS. No
challenge submission was made.

### 2026-08-23 08:17 UTC continuation

Experiment 065 completed successfully and is positive:

- raw PDS `0.7978834143 -> 0.8025098311`
  (`+0.0046264168`);
- official-anchor scaled PDS `+0.0102241256`;
- scaled six-metric average `0.1216402396 -> 0.1236484596`
  (`+0.0020082200`);
- NMAE and Jaccard improved, fidelity and reach regressed slightly, and MSE
  remained clamped at zero;
- build peak process-tree RSS `8.8168 GiB`, scorer peak `15.6500 GiB`, exit
  code `0`.

Compact terminal evidence was mirrored to
`experiments/065-pm-cd4-hct-count/cloud-output/20260823T040540Z/`; all 18
downloaded entries in `artifact-sha256.txt` verified. The 2.28 GB candidate
H5AD remains in S3 with its checksum preserved.

The positive result opened experiment 071's exploratory packaging gate. Its
first local authorization attempt failed closed before cloud dispatch because
the validator expected `output/065-...` while the actual frozen prefix is
`output/exp065-...`. The validator and regression fixture were corrected; all
15 tests then passed under a 1 GiB guard at `0.16 GiB` peak.

Experiment 071 launched at `2026-08-23 08:10:17 UTC`:

- run id `20260823T080948Z`;
- worker `i-0343d04439564b709`;
- output
  `s3://<your-s3-bucket>/output/exp071-production-pm-cd4/20260823T080948Z/`.

At 08:16 UTC context A was complete, the output directory was 8.1 GB, peak
process-tree RSS was 14.18 GiB, and the worker had 59 GiB available. This run
can build and validate a package but cannot submit it.

Experiment 077 remained healthy with both alpha `0.35` and `0.40` scorers
active concurrently on its isolated 128 GiB worker. At 08:16 UTC combined
process-tree RSS was about 58 GiB, peak 59.92 GiB, with 97 GiB available.

The live public leaderboard had 144 entries. Ahmet remained rank 4 at
`0.1377129868` and retained the best raw PDS at `0.8163285024`; Ibrahim
Mansour led at `0.1466648890`.

Experiment 080 is still under implementation repair and has not passed its
post-change tests. Experiment 082's independent audit found unresolved real
staging, score-axis, target-coordinate preprocessing, fold-aggregation, null,
and S3-reuse blockers despite its synthetic suite passing; it remains a
no-launch path. No challenge submission was made.

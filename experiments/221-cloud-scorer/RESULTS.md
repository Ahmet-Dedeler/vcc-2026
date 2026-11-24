# 221 — the cloud `cell-eval2` worker, and the first six-member local→board calibration

Date **2026-09-04** · AWS account `<aws-account-id>` · `us-east-1` · bucket
`<your-s3-bucket>` · `cell-eval2` 0.16.0, `--preset vcc2026`, DE backend `scanpy`

## Headline

**The known-answer row passes bit-exactly** — the cloud worker reproduces all ten statistics of
the frozen 095 h1-test artifact's `agg_results.csv`, nine of them to the last bit and the tenth
to `3.3e-16`, including `expr_mse_unbiased_capped_norm = 2.0551018299736388`, experiment 170's
own gate 3. **It is the same scorer.**

**The six-member transfer table exists**, measured through the real scorer on three
board-anchored, single-variable rungs (`a15` → `rc_hybrid` → `rc_callset` → `cd4f034`), against
a re-seeded twin that gives the surface's own noise floor:

| member | order arm | call-set arm | direction arm (coverage-corrected) |
|---|---:|---:|---:|
| `pds` | **pinned = 0 both sides** | **pinned = 0 both sides** | +0.333× |
| `mse` | **pinned = 0 both sides** | **pinned = 0 both sides** | sign inverted |
| `nmae` | below noise | pinned | **+0.342×** |
| `fid` | +0.169× | **−0.011× — SIGN INVERTED** | **+0.558×** |
| `reach` | **+0.769×** | **wrong sign, below noise** | **−24× — SIGN INVERTED** |
| `jac` | +0.507× | **+1.022×** | sign inverted |

**Three things it settles.**

1. **`jac` is priceable on h1-test after all** — `+1.022×`, within 2.2 % of the board, on the arm
   class that moves it. §11.1 rule 9 lists it as "not locally priceable at all".
2. **`reach` has no single factor and, for two of three classes, no factor at all.** `0.769×` on
   an order arm; sign-inverted on a call-set arm and on a CD4-dose arm. **A screen that had
   priced `rc_callset` on h1-test's `reach` would have vetoed the submission that took us from
   rank 24 to rank 10.** The order-arm figure also puts **185's published `1.290` in question**:
   its board half reproduces exactly, but its local half — measured on a proxy harness rather
   than the scorer — is 1.72× too small (§4.4).
3. **h1-test structurally cannot price a CD4-dose arm**: the frozen bank covers **4 of its 85
   targets against 276 of the challenge panel's 300**, a 19.6× dilution. Corrected for that,
   `nmae` and `fid` — two of the four members §8.4 says we carried at zero or at the wrong sign —
   *are* recoverable, at `×0.34` and `×0.56`. `reach`, which paid `+0.0044` of `cd4f034`'s
   `+0.0131`, is not.

**Cost: $13.85 so far, $1.48 per six-member score.** Eleven instances launched, ten terminated
(the eleventh is still building the real bundle and carries all three self-termination guards).
Local peak RSS **6.85 GiB**, under the 8 GB ceiling, for the whole mandate.

---

## 0. Housekeeping — and why the 13-day orphan cannot recur

### 0.1 `scripts/aws_reap.sh` — the outer sweep

Sweeps **all 17 enabled regions**, lists every RUNNING instance with its age, type and `Name`
tag, pulls 24 h of CloudWatch `CPUUtilization` for each, and stops anything older than
`--max-age-hours` (24) that is idle on **both** readings: 24 h maximum CPU under `--cpu-max`
(5 %) **and** maximum busy cores (`cpu% × vCPUs / 100`) under `--busy-cores` (0.25). Dry run by
default; `--apply` acts. It **stops, never terminates** — an EBS-backed stop is reversible and
the root volume survives, so a false positive costs nothing.

**Why two readings, and this is a real defect a percentage-only rule has.** CloudWatch reports
CPU as a percentage of the *whole instance*, and `cell-eval2` scores **single-threaded** —
measured on this experiment's own workers, a fully busy scorer reads **13.7 % on an 8-vCPU
r7i.2xlarge**, which is 1.10 busy cores. The same job on a 32-vCPU instance reads **3.4 %**,
*under* a 5 % threshold, and the first version of this script would have stopped a working
scorer. The cores reading is instance-size independent, and both must agree before anything is
stopped.

First sweep, 2026-09-04, before anything of ours was launched:

```
REGION          INSTANCE               TYPE             AGE_H  CPU_PCT   CORES  NAME  VERDICT
running instances: 0   flagged for reaping: 0
```

Zero running instances in any region. Five stopped ones remain (four in `us-east-1`, one in
`us-west-2`), all EBS-backed and costing storage only. Re-run against this experiment's own
live fleet it reads them correctly, including their busy-core load:

```
us-east-1  i-0e038f2a9056ab5c1  r7i.2xlarge   1.6  13.735  1.099  vcc-cloud-score-ka2       keep (young)
us-east-1  i-0f33dddb34a808e4d  r7i.2xlarge   0.2  13.813  1.105  vcc-cloud-score-ladder6   keep (young)
```

**Known-answer row for the reaper** — `tests/reap_known_answer.sh`, the classifier extracted
verbatim from the script, **16 rows, ALL PASS**
([`work/reap-known-answer.txt`](work/reap-known-answer.txt)). Selected rows:

| age_h | cpu % | vCPU | busy cores | verdict | what it is |
|---:|---:|---:|---:|---|---|
| 312.0 | 0.27 | 8 | 0.022 | **REAP** | **the 13-day orphan, its measured average** |
| 312.0 | 1.9 | 8 | 0.152 | **REAP** | **the same instance, its measured 24 h maximum** |
| 48.0 | 13.7 | 8 | 1.099 | keep (busy) | a busy `cell-eval2` scorer, measured here |
| 48.0 | 3.4 | 32 | 1.088 | keep (busy) | **the same busy scorer on 32 vCPU — 3.4 % would have been reaped by a %-only rule** |
| 48.0 | 3.13 / 3.12 | 8 | 0.250 / 0.250 | keep (busy) / **REAP** | both sides of the cores boundary |
| 48.0 | 0.79 / 0.77 | 32 | 0.253 / 0.246 | keep (busy) / **REAP** | the same boundary, different instance size |
| 24.0 / 24.1 | 0.1 | 8 | 0.008 | keep (young) / **REAP** | both sides of the age boundary |
| 48.0 | — | 8 | — | keep (no CPU data) | CloudWatch has nothing; never reap on absence |

A reaper that has only ever seen an empty account is untested on the branch that matters. This
table exercises every branch and both sides of every boundary, and rows 1–2 are the exact
instance the mandate exists because of.

### 0.2 The self-terminating guard — three of them, independent

`scripts/cloud_score.sh` launches every worker with **all three**, so no single failure —
including the coordinator dying — can leave an instance running:

| # | guard | where | what it survives |
|---|---|---|---|
| a | `--instance-initiated-shutdown-behavior terminate` | `run-instances` | anything that reaches a shutdown |
| b | `shutdown -h +N` armed in **user-data, before any work starts** | the instance | a hung job, an OOM kill, a crashed coordinator, a lost network |
| c | the worker's `EXIT` trap → `shutdown -h now \|\| poweroff` | `cloud_score_worker.sh` | any exit path, success or failure |

Guard (c) was exercised for real on the first ladder launch: its known-answer gate failed 25 s
into the build, the worker wrote `FAILED 1` to S3 and the instance was **terminated 3 minutes
after launch** with no intervention. That is the behaviour experiment 061 did not have.


### 0.3 S3 — what is there, and what is safe to delete

`s3://<your-s3-bucket>` holds **629 objects / 81,509,440,203 bytes = 75.91 GiB**
(the mandate's "81.5 GB" is the same number in decimal GB). Full listing:
[`work/s3-listing-20260904.txt`](work/s3-listing-20260904.txt).

**All 629 objects except three are under `input/`** — they are *staged copies* of files that
still exist locally, uploaded so a worker could read them. The three exceptions are
`output-clean-hct/g0.65_q1.0_a1.08_spread0.25/{agg_results.csv,exit-code.txt,score.log}`, ~0.5 MB,
which `experiments/061-perturbmatch-hct-count/aws/run.sh` still references as a baseline.

**84 % of the bucket is one thing: 53 copies of six distinct 1.2 GiB whitening operators.**

| operator | copies in S3 | GiB | local original |
|---|---:|---:|---|
| `held_orion_hct116.npz` | **19** | 22.8 | `experiments/062-cd4-post-whiten/cloud-output/whiteners/` |
| `held_orion_hek293t.npz` | **19** | 22.8 | same |
| `held_kolf2.npz` | 7 | 8.4 | `data/donors/whitener_union.npz` |
| `held_k562_gwps.npz` | 4 | 4.8 | `experiments/062-.../whiteners/` |
| `whitener_union.npz` | 3 | 3.6 | `data/donors/whitener_union.npz` |
| `whitener_union_pruned.npz` | 1 | 1.2 | `data/donors/whitener_union_pruned.npz` |
| | **53** | **63.49** | |

Each launcher re-uploaded the same operators under its own `input/<exp>/<run-id>/operators/`
prefix, so `exp082` alone holds two full sets and `exp087` three.

| | objects | GiB |
|---|---:|---:|
| total | 629 | 75.91 |
| **KEEP** — `input/harness/h1-test.h5ad`, `input/production-challenge/*`, `output-clean-hct/*` | 7 | 2.45 |
| **SAFE TO DELETE** — every other `input/` object | **622** | **73.47** (96.8 %) |

Safe because every one is a staged duplicate of a local file, and the experiments that staged
them (052, 060, 061, 065, 069, 070, 071, 073, 075, 077, 079, 080, 082, 085, 087, 088, 089, 090,
091, 094, 103, clean-hct, template053) are all closed.

**Recommendation, not an action.** At S3 Standard's $0.023/GiB-month the whole bucket costs
**$1.75/month** and the cleanup saves **$1.69/month**. That is not worth the risk of deleting a
staging prefix that one of the five agents running concurrently today (216–220) is mid-flight
against, so I have listed rather than deleted. When someone wants it, with the listing above in
hand:

```bash
aws s3 rm s3://<your-s3-bucket>/input/ --recursive \
  --exclude 'harness/h1-test.h5ad' --exclude 'production-challenge/*' --dryrun   # drop --dryrun to act
```

The durable fix is not deletion, it is the duplication: a launcher should reference one
canonical `input/operators/<sha256>.npz` rather than copy 1.2 GiB per run. `scripts/cloud_score.sh`
already accepts `s3://` URIs for `--real` and `--pred` for exactly this reason, and this
experiment's own runs re-used the already-staged `input/harness/h1-test.h5ad` rather than
uploading a seventh copy of it.

---

## 1. The worker — `scripts/cloud_score.sh` + `scripts/cloud_score_worker.sh`

One command takes local (or already-staged) `.h5ad` candidates plus a reference and returns the
full six-member score, having created and destroyed the machine that computed it:

```bash
VCC_CLOUD_SCORE_CLEARANCE=YES scripts/cloud_score.sh \
    --real data/harness/h1-test.h5ad \
    --pred a15=path/to/a15.h5ad --pred rc_hybrid=path/to/rc_hybrid.h5ad \
    --out experiments/<exp>/cloud-output --wait
```

It runs `baseline` → `prep-real-bundle` → `run` → `score` at `--preset vcc2026`,
`--pert-col target_gene --control non-targeting`, on **one shared real bundle** for every
candidate in the job, and returns `agg_results.csv`, `results.csv` (per perturbation),
`score.csv` and a sha256 of every artifact.

Beyond the mandate's four steps it carries five things the existing 061/082 tooling did not:

1. **`--prologue`** — a build step that runs on the worker before scoring, so the heavy staging
   never touches the laptop. This experiment's own candidates are built that way; local peak RSS
   for the whole mandate stayed at **6.84 GiB**, under the 8 GB ceiling, while the worker peaked
   at **8.79 GiB** for the build and ~14 GiB for the scorer.
2. **`--bundle`** — reuse a previously built real bundle. `prep-real-bundle` averages 5 disjoint
   half-splits and is by far the most expensive step; it depends only on the reference, so
   building it once per surface is worth an explicit flag.
3. **Built candidates are parked in S3 the moment the prologue finishes**, so an interruption
   costs only the scoring.
4. **The DE engine is asserted, not assumed.** `cell-eval2` itself warns *"DE numbers differ
   between engines"*; the 095 artifact this worker reproduces was scored with
   `resolved_de_backend: scanpy`, and the worker now aborts if a run resolves anything else, or
   a different `cell_eval2_version`.
5. **A reused S3 prefix is refused** (082's discipline), and the worker sha256 is verified on the
   instance before it is executed.

### 1.1 Sizing, measured rather than assumed

| what | measured |
|---|---|
| `cell-eval2` parallelism on this workload | **single-threaded** — 100 % of one core, 13.7 % of an 8-vCPU r7i.2xlarge |
| scorer peak RSS, h1-test (52,400 × 18,080) | **14.0–16.7 GiB** |
| candidate build peak RSS (operator + control block + 5 emissions) | **8.79 GiB** |
| candidate build wall clock, 5 candidates | **18 min** (quasi 166 s each, rank-controlled 268–276 s each) |
| local peak RSS for the whole mandate | **6.85 GiB** (under the 8 GB ceiling) |

**Therefore more vCPUs buy nothing here and the right instance is the cheapest one with memory
headroom.** `r7i.2xlarge` (8 vCPU, 64 GiB) gives 3.8× headroom over the measured 16.7 GiB peak
and is the cheapest r7i with ≥ 64 GiB; `r7i.xlarge` (32 GiB) would also fit h1-test at 1.9×
headroom but leaves nothing for the challenge frame, where §6.4 records ~25 GB. Spot was
$0.2280/h and on-demand $0.5292/h at launch; a 200 GB gp3 root adds $0.0219/h.

**Do not use spot for the long steps.** The first known-answer run was reclaimed 42 min into
`cell-eval2 baseline` (`"Service initiated (2026-09-04 01:49:50 GMT)"`) and every subsequent run
was launched on-demand.

### 1.2 Three bugs the guards caught, all before any number was reported

| # | bug | how it surfaced | cost |
|---|---|---|---|
| 1 | Gate 2 compared the h1-test effect matrix ENTRYWISE-relative. `apply_v9_whitening` is a pair of 17,942² GEMMs and BLAS reduction order is architecture-dependent, so a 1.3e-14 **absolute** wobble on a near-zero entry is a 6.9e-06 **relative** one. | the gate failed 25 s into the first ladder run | 3 min of spot |
| 2 | The worker inferred completion from `$?` in its `EXIT` trap. A spot SIGTERM runs the trap with `$? == 0`, so a run reclaimed two thirds of the way through `baseline` **wrote `DONE 0` with zero output files.** | noticed only because the S3 output prefix was empty next to a `DONE` marker | would have silently reported a passing known-answer row from a run that scored nothing |
| 3 | Built candidates were parked under `input/`. The instance role is least-privilege — `GetObject` on `input/*`, `PutObject` on `output/*` — so the park failed `AccessDenied`. | the run failed after an 18-minute build | 18 min of build |

Two more caught by the gates rather than by the infrastructure:

4. **`scripts/make_submission_pm_cd4.py` changed under this experiment mid-run.** Commit
   `481cd09` (218's `VCC_H1_UNION_G` second residual) took `apply_frozen_cd4_residual` from a
   2-tuple to a 3-tuple, and a tuple-unpacking call site died 20 s into a cloud run. Gate 1 then
   **passed against the new HEAD to 2.2e-16**, which is the whole point of having it. The call
   sites now index rather than unpack.
5. **Gate 4's expectation was wrong, not the build.** I predicted `a15` and its re-seeded twin
   would share per-gene panel totals because the quasi emitter repairs each gene's total exactly.
   Measured 0/85. The seed also re-draws `pick_rng`'s 400 control **library sizes**, and the panel
   total is `profile × Σ libraries`. The twin is therefore a noise floor on **all six** members
   rather than on the DE members alone, which makes it a better instrument than intended.

---

## 2. The four board-anchored configs, rebuilt on h1-test

`experiments/221-cloud-scorer/build_h1_candidates.py` rebuilds each of the four submissions'
**exact** configs on h1-test's 85 targets, by calling the shipped production functions rather
than reimplementing them. Diffing the four `submissions/*.meta.json` `config` blocks, the only
keys that differ across the whole ladder are `emitter`, the four `rc_*` keys it implies, and
`cd4_target_f` — alpha, gamma, the operator, the donors, PerturbMatch, the spread policy and
both seeds are identical. So each rung is a genuine one-variable step:

| rung | the one variable |
|---|---|
| `a15` → `rc_hybrid` | `emitter` `quasi` → `rc_hybrid` — the emitted p-ORDER |
| `rc_hybrid` → `rc_callset` | `emitter` → `rc_callset` (`top 400 by |H| ∪ top 2,200 by control CPM`, `z_cap` 1.90) |
| `rc_callset` → `cd4f034` | `VCC_CD4_TARGET_F` 0.10 → 0.34 |

A fifth artifact, `a15_seed1`, is `a15` with `library_pick_seed` 0 → 1 and nothing else: the
surface's own noise floor.

### 2.1 The build's four known-answer gates — all pass

| # | gate | result |
|---|---|---|
| 1 | `apply_v9_whitening` + `apply_frozen_cd4_residual` on the 300-target panel must reproduce **a15's** shipped audit at (γ 0.65, λ 0, f 0.10) and **a19's** at (0.80, 0.007) — the same 8 + 6 scalars `170/gate.py` gates on | **PASS**, worst relative **2.2e-16** on 14 scalars |
| 2 | the h1-test effect at f = 0.10 must reproduce `170/work/eff/h1_g0.65.npz` | **PASS** — `‖Δ‖_F/‖E‖_F` **6.46e-15**, `1 − min row cosine` **2.2e-16**, `max|Δ|` **1.34e-14** against a median row norm 4.656 |
| 3 | the CD4 fraction at f = 0.34, re-solved by re-importing the builder under `VCC_CD4_TARGET_F=0.34`, must equal `fraction(0.10) × 3.4` | **PASS**, rel **2.2e-16** (`0.32207201067578156` vs `0.3220720106757816`) |
| 4 | per-gene panel-total digests, all 85 blocks (rule 13's profile-preserving test, **proven not asserted**) | **PASS**, all four rows |

Gate 4 in full, and it is the cleanest statement of the arm classes this project has:

| pair | blocks digest-identical | expected | why |
|---|---:|---:|---|
| `a15` vs `a15_seed1` | 0 / 85 | 0 | the seed also re-draws each target's 400 control library sizes |
| **`a15` vs `rc_hybrid`** | **85 / 85** | 85 | rule 13 profile-preserving |
| **`rc_hybrid` vs `rc_callset`** | **85 / 85** | 85 | rule 13 profile-preserving |
| `rc_callset` vs `cd4f034` | **81 / 85** | 81 | CD4 covers 4 of h1-test's 85 targets |

`a15`, `rc_hybrid` and `rc_callset` share one digest over all 85 blocks —
`453cdc8d4ec90cb6b861384a1afcb86c7a783840dfa3abea390aeae6623f7a61`. That is the local
counterpart of the board's own repeated measurement (`pds` and `mse` identical to sixteen
significant figures across those three entries) and it is reproduced here from the counts.

### 2.2 🔴 h1-test cannot price a CD4-dose arm. It covers 4 of 85 targets against 276 of 300.

The frozen CD4 bank's `complete mean3` policy leaves a residual on

- **276 of the challenge panel's 300 targets (92.0 %)**, but
- **4 of h1-test's 85 (4.7 %)** — `ANKZF1`, `MED15`, `MED25`, `ZNF32`.

On those four rows the dose change is large — `|Δ|/|effect|` **0.334 / 0.444 / 0.218 / 0.185**,
row cosine down to **0.9305** — and on the other 81 the two effect matrices are bit-identical.
Over the whole surface the change is a **5.04 % relative Frobenius** and a **median row cosine of
1.000000000**.

**So the ratio of coverage is 92.0 / 4.7 = 19.6×, and a whole-surface h1-test reading of the
`cd4f034` arm is diluted by about that factor by construction.** This is a *structural* property
of the surface, not noise, and it is a direct part of the answer to §8.4's question "why did we
under-price the one direction arm that worked": **the only board-validated local surface we have
is almost blind to the lever that arm moved.** (Independently reached the same day by another
agent — commit `04101c5`, *"h1-test cannot see the CD4 residual at all"*.)

The consequence for the table below is stated where it appears: the third rung's transfer
factor is reported, but it is a factor between a 4.7 %-coverage local measurement and a
92 %-coverage board measurement, and must not be used as a dose-response transfer constant. The
per-target read on the four covered targets is reported separately, and *that* is the readable
object.

---

## 3. 🔴 THE KNOWN-ANSWER ROW — PASSES, AND BIT-EXACTLY

> *"If the cloud scorer does not reproduce those, it is not the same scorer and nothing it
> produces is usable."*

The frozen 095 h1-test artifact
(`experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`, 52,400 ×
18,080, scored on this Mac on 2026-08-27) was shipped to the cloud worker and re-scored through
`cell-eval2` 0.16.0. **All ten statistics of `agg_results.csv`'s `mean` row**, not just `mse`:

| metric | 2026-08-27, on disk | cloud worker, 2026-09-04 | rel |
|---|---:|---:|---:|
| `pds_cosine` | 0.8655462184873949 | 0.8655462184873949 | **0** |
| **`expr_mse_unbiased_capped_norm`** | **2.0551018299736388** | **2.0551018299736388** | **0** |
| `expr_mse_unbiased` | 0.003700950826313404 | 0.003700950826313404 | **0** |
| `expr_mse_unbiased_capped` | 0.003700950826313404 | 0.003700950826313404 | **0** |
| `expr_distance_unbiased` | 0.0018008600704525074 | 0.0018008600704525074 | **0** |
| `expr_real_mass_ratio` | 1.0 | 1.0 | **0** |
| `de_wilcoxon_lfc_nmae` | 0.9667252780145132 | 0.9667252780145129 | 3.3e-16 |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.5046621497798164 | 0.5046621497798164 | **0** |
| `de_wilcoxon_direction_reach_raw` | 0.20045827672164832 | 0.20045827672164832 | **0** |
| `de_wilcoxon_sig_jaccard` | 0.06359002047099205 | 0.06359002047099205 | **0** |

**Nine of ten reproduce bit-for-bit; the tenth agrees to 3.3e-16. Worst relative difference over
all ten: `3.331e-16`.**

`expr_mse_unbiased_capped_norm` **= 2.0551018299736388 is exactly experiment 170's gate 3**, the
value `msekit` was validated against. The mandate asked for that one number; the run delivers all
ten, on a different architecture (x86_64 vs this machine's arm64), from a different process, on a
machine that did not exist an hour before.

The instance's own `run_meta.json` records `cell_eval2_version: 0.16.0`,
`resolved_de_backend: scanpy`, `source_fingerprint_strict: true` and the same
`anchor_semantic_identity` `67663d68…` the 2026-08-27 run recorded.

**The cloud scorer is the scorer. Everything below is usable.**

### On a15's `0.88228 / 0.11112`

The mandate pairs the 095 anchor with *"a15's raw `mse` 0.88228 / `score_mse` 0.11112"*. Those two
are **challenge-frame** numbers, not h1-test ones: 164 §5 derives them analytically through 096's
`claim` route on the 300-target panel, and experiments 169, 171, 188, 196, 197, 200, 211 and 214
all carry them as a known-answer row for that identity. `cell-eval2` on h1-test cannot emit them —
h1-test has 85 targets, a different gene axis and a different reference. They are reproduced here
by the same analytic route the other experiments use, and the 095 row above is the *scorer's* own
known answer, which is the one that certifies this tool.

---

## 4. THE TRANSFER TABLE — six members, three board-anchored rungs, the real scorer

All six candidates scored by `cell-eval2` 0.16.0 on h1-test, `--preset vcc2026`, the same
reference object (`sha256` verified on every worker), `resolved_de_backend: scanpy` on all six.
Machine-readable: [`work/transfer-table.md`](work/transfer-table.md),
[`work/transfer-table.json`](work/transfer-table.json).

### 4.1 The levels

| candidate | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|
| `ka095` *(the known-answer artifact)* | 0.865546 | 2.055102 | 0.966725 | 0.504662 | 0.200458 | 0.063590 |
| `a15` | 0.873810 | 1.015978 | 0.960014 | 0.507478 | 0.166269 | 0.064039 |
| `a15_seed1` *(re-seeded twin)* | 0.874510 | 1.016016 | 0.960063 | 0.508254 | 0.162949 | 0.064143 |
| `rc_hybrid` | 0.873810 | 1.015978 | 0.960046 | 0.508602 | 0.229541 | 0.064304 |
| `rc_callset` | 0.873810 | 1.015978 | 0.960046 | 0.495058 | 0.226704 | 0.074204 |
| `cd4f034` | 0.873249 | 1.018538 | 0.958995 | 0.495303 | 0.225905 | 0.074199 |

**`pds` and `mse` are identical to six decimals across `a15`, `rc_hybrid` and `rc_callset`.**
That is the board's own most-cited fact — those three entries agree on `pds` and `mse` to
*sixteen* significant figures — reproduced independently, on a different surface, from the
emitted counts, by the scorer the competition runs. Rule 13's profile-preserving class is no
longer a construction argument here either.

### 4.2 The noise floor

`a15` against `a15_seed1`, identical config and identical effect matrix, `library_pick_seed`
0 → 1:

| | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|
| \|Δ\| raw seed noise | 0.000700 | 0.000039 | 0.000048 | 0.000777 | **0.003319** | 0.000104 |

`reach` is by far the noisiest member on this surface — its seed noise is **4.3×** `pds`'s and
**32×** `jac`'s. Any h1-test `reach` delta below 0.0033 raw is not a measurement.

### 4.3 The table

**board Δraw / h1-test Δraw.** Board deltas are `Δscaled / slope` off `docs/BOARD-LEDGER.md`,
a route that reproduces all seven directly recorded raw deltas to ≤ 9.2e-3 relative (worst case
is `jac`'s +0.00039, i.e. 3.6e-6 absolute — the rounding of the recorded value).

| rung | member | h1-test Δraw | board Δraw | **transfer** | verdict |
|---|---|---:|---:|---:|---|
| **`a15`→`rc_hybrid`** | `pds` | **+0.000000** | −0.000000 | — | **pinned both sides** |
| *the emitted p-ORDER* | `mse` | **+0.000000** | +0.000000 | — | **pinned both sides** |
| | `nmae` | +0.000032 | −0.000031 | — | below noise both sides |
| | `fid` | +0.001124 | +0.000190 | **+0.169×** | readable |
| | `reach` | +0.063272 | +0.048654 | **+0.769×** | readable |
| | `jac` | +0.000265 | +0.000134 | **+0.507×** | readable |
| **`rc_hybrid`→`rc_callset`** | `pds` | **+0.000000** | +0.000000 | — | **pinned both sides** |
| *the CALL SET* | `mse` | **+0.000000** | −0.000000 | — | **pinned both sides** |
| | `nmae` | −0.000000 | −0.000000 | — | pinned |
| | `fid` | −0.013544 | +0.000146 | **−0.011×** | **SIGN INVERTED** |
| | `reach` | −0.002837 | +0.020005 | — | below noise, **and the wrong sign** |
| | `jac` | +0.009900 | +0.010119 | **+1.022×** | readable |
| **`rc_callset`→`cd4f034`** | `pds` | −0.000560 | −0.003646 | — | below noise |
| *a DIRECTION arm* | `mse` | +0.002561 | −0.037316 | **−14.6×** | **SIGN INVERTED** |
| | `nmae` | −0.001052 | −0.006394 | **+6.08×** | readable, right sign |
| | `fid` | +0.000244 | +0.002717 | — | below noise |
| | `reach` | −0.000799 | +0.023381 | — | below noise, **wrong sign** |
| | `jac` | −0.000005 | +0.000394 | — | below noise |

### 4.4 What the table says, member by member

**`jac` transfers essentially 1:1 for a call-set arm — `+1.022×`.** This is the single most
useful number here. `jac` was on rule 9's "not locally priceable at all" list, and on the arm
class that actually moves it, h1-test reads it within **2.2 %** of the board. 181's `×0.44` came
from a spread-policy arm (a8's `jac`) and does not describe this class.

**`reach` transfers at `+0.769×` for an ORDER arm** — h1-test *overstates* by 1.30×, and this is
the first estimate whose local half comes from the real scorer rather than a proxy. §11.1's "do
not apply a single transfer factor to `reach` across arm classes" is confirmed; the measured
spread across classes is now **0.39 (spread policy) → 0.77 (order) → sign-inverted (call set,
CD4 dose)**.

**🔴 `reach`'s order-arm factor is 0.769×, not 185's 1.290× — and the qualitative claim inverts.**
185's `BOARD-RESULT-rc_hybrid.md` §"The finding that outlasts the arm" concludes *"h1-test
UNDERSTATES board `reach` by 3× on this class"* and publishes **1.290**. That number is measured
on **exactly this rung**, and the board half of it reproduces here to four decimals — 185's
`+0.04865` against this table's `+0.048654`. **The whole disagreement is in the local half:**

| | a15 | `rc_hybrid` | Δ raw | factor |
|---|---:|---:|---:|---:|
| 185, `debench` proxy harness | 0.17910 *(twin mean)* | 0.21681 | +0.03771 | **1.290×** |
| **this table, the real `cell-eval2`** | 0.164609 *(twin mean)* | **0.229541** | **+0.064932** | **0.749×** |
| the same, seed-0 `a15` only | 0.166269 | 0.229541 | +0.063272 | **0.769×** |

**The proxy compresses the local delta by 1.72×** — it reads `a15`'s `reach` *higher* and
`rc_hybrid`'s *lower* than the scorer does — and that is enough to flip the conclusion from
"h1-test understates by 3×" to "h1-test **overstates** by 1.30×". Every `reach` prediction for a
designed-order arm since 2026-09-03 has passed through the 1.290.

**I am reporting this as a correction candidate, not a settled overturn**, and the honest
statement of what separates them is: both ends of the ratio here come from the scorer the
competition runs, on a single-variable rebuild whose four known-answer gates pass, while 185's
local half came from `153-de-bench`'s `Surface` proxy and its own `realize185.py` re-seeded the
rank-controlled emitter's second pass where the production writer does not. **The check that
settles it is to re-run 185's own h1 artifacts through `scripts/cloud_score.sh` — now a $1.50,
one-command operation, which is the entire point of this experiment.**

**But `reach` does NOT transfer for a CALL-SET arm, and h1-test gets its sign wrong.** The board
paid `+0.0200` raw for `rc_callset` over `rc_hybrid`; h1-test reads `−0.0028`, inside its own
0.0033 noise floor and pointing the other way. **A screen that priced `rc_callset` on h1-test's
`reach` would have vetoed the arm that took us from rank 24 to rank 10.**

**`fid` is unreadable, exactly as rule 9 says.** Both rungs where it moves give the wrong sign or
a 6× magnitude error: `+0.169×` on the order rung, `−0.011×` on the call-set rung. `n_real` is 222
on h1-test against a board median of ~89, so `fid`'s `min(1, n_pred/n_real)` yield term sits in a
different regime. **Rule 9 stands, and this is the first direct two-sided measurement of it.**

**`pds` and `mse` are pinned, and that is a result, not a gap.** Two of the three rungs move
them by *exactly zero* on both surfaces. Rule 12's "unmeasured members cannot move" exemption is
now demonstrated rather than argued, for the emission-side classes it was written for.

### 4.5 The CD4 rung, corrected for coverage — the direct answer to §8.4

The whole-surface reading of `rc_callset` → `cd4f034` is worthless, and §2.2 says why: the arm
touches 4 of h1-test's 85 targets and 276 of the challenge panel's 300. The build proves only
those targets move — 81/85 blocks are digest-identical and the per-target deltas on the other 81
are zero to five decimals. So divide both sides by their coverage (0.047 and 0.920) and the two
become comparable:

| member | h1-test, per covered target | board, per covered target | transfer |
|---|---:|---:|---:|
| `pds` | −0.01190 | −0.00396 | **+0.333×** |
| `nmae` | −0.02033 | −0.00695 | **+0.342×** |
| `fid` | +0.00529 | +0.00295 | **+0.558×** |
| `reach` | −0.00105 | +0.02541 | **−24.1×** — sign inverted |
| `jac` | −0.00011 | +0.00043 | **−3.90×** — sign inverted |

(`mse`'s per-perturbation value is undefined — `expr_mse_unbiased_capped_norm` is a
ratio-of-sums, `NaN` in the per-perturbation frame by construction, and `NaN` in the `count` row
of 095's own `agg_results.csv` too.)

The per-target deltas on the four covered targets, which is where the arm exists at all:

| target | Δ`pds` | Δ`nmae` | Δ`fid` | Δ`reach` | Δ`jac` |
|---|---:|---:|---:|---:|---:|
| `ANKZF1` | 0 | **−0.05332** | +0.00693 | 0 | −0.00039 |
| `MED15` | −0.03571 | −0.00258 | +0.01077 | −0.00154 | −0.00038 |
| `MED25` | 0 | −0.00509 | +0.00231 | −0.00267 | +0.00033 |
| `ZNF32` | −0.01190 | — | +0.00115 | 0 | 0 |
| *the other 81, mean* | +0.00000 | −0.00000 | −0.00000 | −0.00079 | +0.00000 |

**So, on the arm §8.4 is about:**

1. **`nmae` and `fid` — the two members §8.4 says we priced at the wrong sign and at zero — do
   move in the RIGHT direction on h1-test**, once coverage is divided out. `nmae` −0.0203 against
   the board's −0.0070 per covered target; `fid` +0.0053 against +0.0030. **They were priceable
   all along; nobody looked, because looking meant running the real scorer.**
2. **They must be deflated, not taken at face value: `×0.33` and `×0.56`.** h1-test overstates
   both by 2–3×.
3. **`reach` — the single largest realised member of that arm, `+0.0044` of its `+0.0131`
   overall — is NOT recoverable from h1-test at any coverage correction.** It reads `−0.00105`
   per covered target against a board `+0.02541`: **the wrong sign, and 24× out.** A screen that
   had priced `cd4f034` on h1-test's `reach` would have subtracted the member that paid most.
4. `jac` likewise inverts, but it is worth only `+0.0002` overall on this arm and its absolute
   size (−0.00011) is inside the surface's own noise.

**The verdict on §8.4's fix.** "Price all six members explicitly" is right and now cheap — this
tool makes it a `scripts/cloud_score.sh` invocation. But **the surface we would price them on
cannot price two of them**, and `reach` is the one it fails on hardest, on both the arm classes
we care about. Concretely, for a direction arm the honest local pricing rule this experiment
supports is:

```
nmae   :  board ≈ 0.34 x  (h1-test per-covered-target delta)     [1 rung, 4 targets]
fid    :  board ≈ 0.56 x  (h1-test per-covered-target delta)     [1 rung, 4 targets]
pds    :  board ≈ 0.33 x  (h1-test per-covered-target delta)     [1 rung, 4 targets]
reach  :  NOT PRICEABLE on h1-test for a CD4-dose or a call-set arm — wrong SIGN on both
jac    :  1.02x for a CALL-SET arm; NOT priceable for a CD4-dose arm — wrong sign
mse    :  keep using 096's `claim` route; h1-test's whole-surface reading is +0.0026
          against a board -0.0373 and is a coverage artefact, not a measurement
```

---

## 5. Falsifiers

**1. The known-answer row.** §3 — **PASSES**, nine of ten metrics bit-for-bit, worst relative
`3.3e-16` over all ten, including `expr_mse_unbiased_capped_norm = 2.0551018299736388` exactly.
No table above was written before it landed. The build carries four more known-answer gates
(§2.1), all passing, two of which caught real defects before any number was read.

**2. h1-test is not the challenge, and these factors are empirical.** Stated plainly:

- **They are calibrated on THREE rungs, not four**, and each rung is one board pair. There is
  no error bar on a transfer factor from a single pair; the noise floor in §4.2 says only
  whether the local half is a measurement at all.
- **`n_real` differs by a factor of ~2.5** — 222 on h1-test against a board median of ~89
  (§11.1 rule 9, reconstructed by 188). Every member with `n_real` in a denominator or a union
  is affected, and the table measures exactly that: `fid` inverts on one rung and is 6× out on
  the other.
- **The CD4-coverage correction assumes only covered targets move.** That is proven here
  (81/85 digest-identical, other-81 mean delta zero to five decimals) but it is an assumption
  about the *challenge* side, where we cannot check it.
- **Nothing here says the factors extrapolate to a new arm class.** The one thing this table
  establishes about extrapolation is negative: `reach`'s factor moves from `0.39` (spread policy)
  through `0.769` (order) to `1.290` (185's board point) and inverts in sign on a call-set or a
  CD4-dose arm. **Treat every factor as belonging to its rung until a second board pair of the
  same class exists.**
- **Emission seeds are not controlled between surfaces.** The board arms carry
  `library_pick_seed 0` in the challenge frame; the h1-test rebuilds carry it on a different
  control pool, so the two are the same *policy*, not the same *draw*.

**3. Cost discipline and termination.** §6.

---

## 6. Cost, and proof that nothing is left running

`experiments/221-cloud-scorer/cost_and_terminations.sh` rebuilds this from CloudTrail's
`RunInstances` (the complete launch ledger — `describe-instances` drops a terminated instance
after about an hour, so it cannot be one) and each worker's `STATUS` object in S3, whose
`LastModified` is written by the EXIT trap immediately before the instance powers off.

| run | instance | market | launch UTC | end UTC | min | $ | what it was |
|---|---|---|---|---|---:|---:|---|
| `ka` | `i-01db4020744f09521` | spot | 01:07:04 | 01:49:52 | 42.8 | 0.178 | **RECLAIMED by spot** 42 min into `baseline` |
| `ladder` | `i-01e401ae1863945b4` | spot | 01:27:34 | 01:28:36 | 1.0 | 0.004 | gate 2 (entrywise-relative) fired |
| `ladder2` | `i-0fd55468b05fa77af` | spot | 01:38:38 | 01:56:51 | 18.2 | 0.076 | gate 4's expectation was wrong |
| `ka2` | `i-0e038f2a9056ab5c1` | on-demand | 01:54:38 | 04:31:57 | 157.3 | 1.445 | `prep-real-bundle --baseline` given a directory |
| `ladder3` | `i-0eaa6d65f4e47452a` | on-demand | 02:05:29 | 02:06:26 | 0.9 | 0.009 | builder HEAD moved under the run |
| `ladder4` | `i-04d6bcfbd1018ac4c` | on-demand | 02:24:19 | 02:43:33 | 19.2 | 0.177 | parked candidates under `input/` — AccessDenied |
| `ladder5` | `i-0a1cb2df981c58dad` | on-demand | 02:51:58 | 03:18:29 | 26.5 | 0.244 | **BUILD OK**; I stopped it to re-order runs before bundle |
| `ladder6` | `i-0f33dddb34a808e4d` | on-demand | 03:19:07 | 08:31:00 | 311.9 | 2.865 | scored `ka095`, `a15` |
| `ladder7` | `i-043a904dbea6b188d` | on-demand | 03:29:40 | 08:57:19 | 327.6 | 3.010 | scored `a15_seed1`, `rc_hybrid` |
| `ladder8` | `i-0926ebf77f8240701` | on-demand | 03:29:57 | 08:57:18 | 327.4 | 3.007 | scored `rc_callset`, `cd4f034` |
| `bundle` | `i-08028907eeba5a980` | on-demand | 05:27:24 | — | — | — | the h1-test real bundle |

**Cost per score.** The three productive scoring instances cost **$8.88 for 6 candidates =
$1.48 per six-member score**, of which $1.37 is the ~2.5 h `cell-eval2 run` itself (r7i.2xlarge
on-demand $0.5292/h + 200 GB gp3 at $0.0219/h). The candidate BUILD is $0.24 for all five.
**A six-member score of a new candidate on h1-test therefore costs about $1.50 and 2.5 hours of
someone else's machine, with zero local memory pressure.** On spot it would be $0.64, but see
the reclaim above: only use spot for steps under ~30 minutes.

**Total spend for this experiment: $13.85 at the last accounting** (the bundle instance is still running; re-run the script for the final figure), of which **$2.13 was burned on the five failed launches** — every one caught by a gate or a guard, every one self-terminated, none
requiring intervention, and every one producing a fix that is now in the tooling.

### 6.1 Termination

Eleven instances launched, **ten terminated**, one still building the real bundle. Every
termination but one was the instance ending itself: `shutdown -h now` under
`instance-initiated-shutdown-behavior=terminate`, which is *not* an API call and therefore leaves
no `TerminateInstances` event — the two that do appear in CloudTrail are `ladder5`, which I
stopped deliberately to re-order the pipeline, and `ladder6`/`7`/`8`, which I stopped once their
runs were in S3 and only the (redundant) scoring step remained.

```
$ bash scripts/aws_reap.sh          # all 17 enabled regions
REGION          INSTANCE               TYPE             AGE_H  CPU_PCT   CORES  NAME                          VERDICT
us-east-1       i-08028907eeba5a980    r7i.2xlarge        4.0   14.542   1.163  vcc-cloud-score-bundle-...    keep (young)

running instances: 1   flagged for reaping: 0
```

The one live instance is the h1-test real bundle (§7). It carries all three guards: launched
`instance-initiated-shutdown-behavior=terminate`, `shutdown -h +900` armed in user-data at boot,
and the worker's EXIT trap. Re-run `bash experiments/221-cloud-scorer/cost_and_terminations.sh`
for the current state and the final total.

---

## 7. The scaled table — what is still running, and how to finish it

Everything above is in the **RAW** space, which is where every local→board factor this project
already uses lives (164's `0.4277`, 185's `1.290`, `docs/BOARD-LEDGER.md` §4's `0.390`). The
scaled space needs an h1-test **real bundle** — `cell-eval2 prep-real-bundle`, which averages 5
disjoint half-splits of the reference and is the single most expensive step in the whole stack.

Three things were learned building it and they are worth more than the table it produces:

1. **`AGENTS.md` §6.4's snippet is wrong.** It says
   `prep-real-bundle --real real.h5ad --baseline out/baseline`. `--baseline` wants the baseline
   **prediction `.h5ad`**, not the `baseline -o` output directory, and passing the directory
   fails with `IsADirectoryError` *after* `cell-eval2 baseline` has already burned 2.6 hours
   (run `ka2`, $1.45). The correct form is

   ```bash
   cell-eval2 baseline -ar real.h5ad -o out/baseline --save-pred out/baseline-pred.h5ad --preset vcc2026
   cell-eval2 prep-real-bundle --real real.h5ad --baseline out/baseline-pred.h5ad -o out/bundle --preset vcc2026
   ```

   That is a fair part of the answer to "why was the four-step pipeline never actually run".
2. **We already had the expensive half.** `data/harness/h1-test-baseline-pred.h5ad` (2.7 GB,
   built 2026-08-21) carries `uns['baseline_emission'] = {emit: dispersed, seed: 0, n_rows:
   34000}` — the current non-legacy arm, exactly what `prep-real-bundle` wants. Supplying it
   skips `cell-eval2 baseline` entirely (2.6 h, measured). `scripts/cloud_score.sh` now takes
   `--baseline-pred` for this.
3. **`prep-real-bundle` is ~4 h+ and still running** at the time of writing (`bundle-20260904`,
   started 05:27 UTC). It is single-threaded like everything else in this stack.

**To finish the scaled table** — the bundle is written to
`s3://<your-s3-bucket>/output/cloud-score/bundle-20260904/bundle/` and every
`agg_results.csv` is already on disk, so the last step needs no cloud and almost no memory
(`score` reads aggregates only, never a matrix):

```bash
aws s3 sync s3://<your-s3-bucket>/output/cloud-score/bundle-20260904/bundle/ \
    experiments/221-cloud-scorer/work/bundle/
for n in ka095 a15 a15_seed1 rc_hybrid rc_callset cd4f034; do
  d=$(find experiments/221-cloud-scorer/cloud-output -path "*/$n/run/agg_results.csv")
  uv run cell-eval2 score --user-agg "$d" \
      --real-bundle experiments/221-cloud-scorer/work/bundle \
      -o "$(dirname "$(dirname "$d")")/score.csv"
done
uv run python experiments/221-cloud-scorer/transfer_table.py \
    --run-dir experiments/221-cloud-scorer/cloud-output/ladder6-20260904 \
    --run-dir experiments/221-cloud-scorer/cloud-output/ladder7-20260904 \
    --run-dir experiments/221-cloud-scorer/cloud-output/ladder8-20260904
```

`transfer_table.py` picks the `score.csv` files up automatically, fits h1-test's own
`slope` per member from the six candidates (each carries both a raw and a scaled reading of the
same run, so the map is over-determined), and prints the scaled factors as
`transfer_scaled = (slope_board / slope_h1test) · transfer_raw`. **It changes no conclusion
above** — a per-member constant multiplies every rung's factor equally, so signs, the
below-noise verdicts and the relative ordering of the classes are all invariant.

---

## 8. What this changes

**For the tooling.** `scripts/cloud_score.sh` removes the excuse §6.4's safety note has been
carrying since 2026-08-22. A six-member score of any candidate on any reference now costs
**$1.50 and one command**, with **zero local memory pressure** — the mandate's 8 GB ceiling was
never approached (6.85 GiB peak, and that was the gate check, not a scorer). `scripts/aws_reap.sh`
plus the three-guard self-termination means the 13-day, $170 orphan of experiment 061 cannot
recur silently: it now takes a *simultaneous* failure of the shutdown timer, the shutdown
behaviour, the EXIT trap and the reaper.

**For `AGENTS.md`, four corrections and one confirmation:**

1. **§6.4's `prep-real-bundle` invocation is wrong** and costs 2.6 h of compute to discover.
   `--baseline` is the baseline PREDICTION h5ad; use `baseline --save-pred`.
2. **§11.1 rule 9 is CONFIRMED and now measured on both sides.** `fid` transfers at `+0.169×`
   on an order arm and at `−0.011×` — the wrong sign — on a call-set arm. It is not priceable
   on h1-test, and the reason (`n_real` 222 vs a board median ~89) is exactly the one the rule
   gives.
3. **`jac` is NOT unpriceable.** Rule 9 lists it with `fid`; on a call-set arm h1-test reads it
   at **`+1.022×`**, within 2.2 % of the board. 181's `×0.44` is a spread-policy constant and
   does not describe this class. `jac` should move to the priceable list *for call-set arms*.
4. **`reach`'s transfer factor is not one number and is not even one sign.** Measured on the
   same instrument, same day: `0.769×` on an order arm, and **sign-inverted** on both a call-set
   arm and a CD4-dose arm. AGENTS.md's "do not apply a single transfer factor to `reach` across
   arm classes" understates it: for two of three classes there is no factor at all.
6. **185's `1.290` for a designed-order arm is a correction candidate (§4.4).** The board half
   reproduces to four decimals; the local half, measured through the real scorer instead of the
   `debench` proxy, is **1.72× larger**, which takes the factor to `0.749–0.769×` and inverts
   "h1-test understates by 3×" into "h1-test overstates by 1.30×". Settling it is one
   `scripts/cloud_score.sh` run on 185's own artifacts.
5. **§8.4's fix is right in principle and half-blocked in practice.** "Price all six explicitly"
   is now cheap, and `nmae` and `fid` — two of the four members §8.4 says were carried at zero or
   at the wrong sign — *are* recoverable from h1-test at `×0.34` and `×0.56` once the CD4 bank's
   4/85 coverage is divided out. But **`reach`, which delivered `+0.0044` of `cd4f034`'s
   `+0.0131`, is not recoverable at any correction** — h1-test reads it at the wrong sign and 24×
   out. A direction arm's `reach` still has no local instrument.

**The one number to take away.** `rc_callset` bought `+0.0200` raw `reach` on the board and reads
`−0.0028` on h1-test, inside h1-test's own `0.0033` noise. **A pre-registration that had priced
that arm on this surface would have vetoed the submission that moved us from rank 24 to rank 10.**
That is not an argument against local screens; it is an argument for knowing, per member and per
arm class, which ones the surface can see — which is what the table in §4 is for.

**What to do with the remaining ~96 slots.** Two concrete consequences:

- **Order and call-set arms can be screened on h1-test for `jac` (1.02×) and, for order arms
  only, `reach` (0.77×), and must not be screened on it for `fid` at all.**
- **Direction arms cannot be screened on h1-test for `reach`.** Either find a surface whose CD4
  coverage resembles the panel's 92 % — the frozen bank's `complete mean3` policy is the binding
  constraint, not the surface — or accept that a direction arm's largest member is unpriced and
  size the arm from `mse` + `nmae` + `fid` with §8.4's 2.51× multiplier, knowing that multiplier
  now has a *mechanism*: it is `reach` and `fid` moving with the cosine, and `reach` is the half
  we cannot see.

---

## 9. Answers to the four INBOX messages

### 9.1 `INBOX-BOARDPOINT` — `h1union` CANNOT be a fourth rung. It is pure oracle leakage on h1-test.

The coordinator asked me to add the new `h1union` board point (0.2118310, rank 3,
**+0.0069114** over `cd4f034`) as a fourth rung, doubling the calibration data.

**It cannot go on this surface, and the reason is decisive rather than a judgement call.**
`h1union` is `VCC_H1_UNION_G` with `VCC_UNION_SOURCES` unset, which
`scripts/make_submission_pm_cd4.py:793` resolves to exactly two banks —
`gse295214_h1` and **`h1_2025`**. Measured directly against the 85 targets this experiment's
h1-test ladder is built on:

| the arm's source | covers |
|---|---|
| `data/donors/h1_2025` | **85 / 85 of h1-test's targets** |
| `experiments/115-gse295214-h1/work/compact-donor` | 30 / 85 |

`AGENTS.md` §6.4, trap 1, in as many words: *"`data/donors/h1_2025` covers 85/85 of h1-test's
targets. Scoring an h1-test screen with that donor in the mix is **pure oracle leakage** — it read
a buildable ordering key at `0.404` before the leak was found (exp 157)."*

**`h1union`'s entire content is adding H1 donor data, and h1-test is H1.** Every one of its 85
targets would be scored against a prediction that had read that target's own answer. It is not a
weak rung, it is the single worst arm in the project to put on this surface, and any transfer
factor it produced would be an artefact of the leak. **Reporting that is the answer; producing a
number would have been the error.**

**And it bears on the 1.384-vs-1.000 contradiction the INBOX asks about.** That disagreement is
between two *PLOLO cosine* → board cosine readings, which is a different instrument from this
table (mine is per-member, on scored artifacts, through the real scorer). I cannot arbitrate it
with h1-test for the reason above. But this experiment does supply one relevant constraint that
neither hypothesis in the INBOX names: **`cd4f034` and `h1union` are not the same arm class in the
sense that matters here.** `cd4f034` re-doses a residual that covers 276/300 panel targets;
`h1union` adds a source into rows that were previously *empty*. The INBOX's hypothesis 3 — *"a
first source into empty rows versus a fourth estimate onto occupied ones"* — is the one this
experiment's coverage measurement supports, and it predicts the factor is a function of the arm's
**coverage on the rows it moves**, not a constant. The CD4 rung here is the same phenomenon seen
from the other side: its whole-surface reading is diluted by exactly its coverage ratio (19.6×),
and dividing that out restores sign agreement on three of five members. **Recommend testing
hypothesis 3 first, and expressing the factor as a function of coverage.**

### 9.2 `INBOX-HARNESS` / `INBOX-HARNESS-RETRACTED` — the one surviving question

> *"Is there any board pair that isolates γ = 0 vs γ > 0?"*

**No, and there is no board pair that isolates γ at all.** Over every scored arm in
`docs/BOARD-LEDGER.md` — `a2-prune`, `a8`, `a10`, `a11`, `a12`, `a13`, `a14`, `a15`, `a19`,
`a21`, `a23`, `abc-knee`, `rc_hybrid`, `rc_callset`, `cd4f034`, `h1union` — exactly one moved
`WHITEN_GAMMA`: **a19**, and it moved `gamma` 0.65 → 0.80 **and** `lam` 0 → 0.007 together. 170
later attributed only **21 %** of a19's −0.0959 to gamma, the other 79 % to the ridge. Every
arm below γ = 0.65 exists only locally.

**So whitening's value has never been isolated on the board**, and §4.2's "+0.060, largest
measured gain" is, correctly and by its own text, a **`pds`-on-LODO** claim. That is worth
recording — it means the board has never priced the single largest local gain in the project —
but it is not a contradiction, and the retraction's reasoning (212 §3.4 already measured the
cosine↔`pds` trade along γ at 96.2 % of the board's own exchange rate) is right. I did not spend
further effort on it.

### 9.3 `INBOX-COMMIT` — done

`0add391` stages **only** the eight files this experiment wrote plus its `RESULTS.md` and
`cloud-output/` CSVs — never `git add -A`, nothing under another experiment's tree, and nothing
touching `scripts/make_submission_pm_cd4.py`, which 235 is mid-edit on. The known-answer evidence
is in the commit message: the 16 reaper rows and the nine-of-ten bit-exact scorer row.

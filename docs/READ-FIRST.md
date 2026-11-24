# READ-FIRST — current state and reading map

Updated September 9, 2026, 21:10 UTC (September 10, 00:10 Istanbul).
This page supersedes dated narrative in `AGENTS.md`,
`CURRENT_STATE.md`, and earlier versions of this page. It is a router, not a lab
notebook. The previous page is preserved byte-for-byte at
`docs/archive/READ-FIRST-through-2026-09-08.md`; historical line citations to
READ-FIRST resolve there.

**September 9 research window has ended:** Ahmet authorized work through
**20:58 UTC / 23:58 Istanbul, September 9**. The final comparison is
`docs/SPRINT-HANDOFF-2026-09-09.md`; 18 experiments are terminal and none
cleared a new submission candidate. Ahmet assesses and decides submissions. Codex
coordinates; **no built-in Codex subagents**. Scoped **Claude Code** workers
are now allowed (Fable 5.1 xhigh; Opus 5 max fallback after verifying availability).
This supersedes earlier blanket no-subagent language, not the no-submission,
no-paid-cloud, no-reserve-access rules. Read `docs/SPRINT-2026-09-09.md` first
for the current queue, runtime status and cutoff. CLI authorization was pending
at setup; at **17:45 UTC** the existing Max login was confirmed and the first
restricted read-only `claude-fable-5-1` xhigh review started. Codex still owns
implementation, numerical runs and verification. No paid API or overage is
being used in its initial rate-limit record.

## 1. Current state — verified September 9

Authenticated/public status refreshed at **21:10 UTC, September 9 /
00:10 Istanbul, September 10**:
319 remains published, **rank 13 of 839**, with unchanged **0.21865973944346428**.
No new submission. This supersedes the **08:19:43 UTC** rank-6 record in
`experiments/323-current-family-validation/work/final-verification.json`.
Rank movement is not a model change. Current receipt:
`docs/sprint-2026-09-09/live-verification-closure.json`.
No research was launched after the cutoff; only final verification/closure.

Public board fetched at **04:27:12 UTC / 07:27:12 Istanbul**:
**5th of 821**, **0.21865973944346428**, entry
`xuDelEik7zGWMXdNu0R4`, model `the parrot approved my expenses`.
This is the best score in the retained history, not a regression. Leader:
`ViRtUaL cElL builders`, **0.2620924638410216**; gap **0.0434327243975573**.
Both values and the exact six-member decomposition are in
`docs/review-2026-09-09/state.json`. A fresh `vcc --json status` independently
confirmed our entry is published and rank 5.

| UTC date | scored arm | score | record |
|---|---|---:|---|
| Sep 7 | 312, plain exponent −1.5 | 0.2177188848 | `data/board-entries/Q519I2pLWRCW7PAEx5ul.json` |
| Sep 7 | 312b, plain exponent −3.0 | 0.2164511429 | `data/board-entries/1n5LFk6PSEHq15UHENIE.json` |
| Sep 8 | 316, PM parity at −1.5 | 0.2169732875 | `data/board-entries/2K5ppxkowMSj2ugoZqYb.json` |
| Sep 8 | **319, cap at −1.5** | **0.2186597394** | `data/board-entries/xuDelEik7zGWMXdNu0R4.json` |

Both September 8 slots were used. No September 9 submission is recorded.
On September 9 the authenticated CLI reports `can_submit=true` and the installed
client's read-only limits method reports `limit_reached=false`; that endpoint
does **not** return a numeric remaining allowance. Reset is **00:00 UTC =
03:00 Istanbul**, not local midnight.

No VCC builder, scorer, submission waiter, or headless research worker was found
in the takeover process check. An existing `board_snapshot.sh` process remains;
it was not modified. Other interactive app sessions are not evidence of active
VCC research. State can change: recheck immediately before action.

**Initial user direction:** main-agent work with no subagents.
The later sprint authorization above now
allows scoped Claude Code workers only; built-in Codex subagents remain forbidden.
Do not auto-schedule submissions at reset or spend two slots merely
because they become available. Prefer scientific progress and a considered
information-gain decision. No submission or new automation was authorized by
this audit.

## 2. Best artifact, pending artifact, and defaults are different

- Best A/B/C package: `submissions/319-normrel15-cap.vcc`,
  SHA-256 `57585d558f9777f066c02f5b299cacccca0dd4954b4201b8fa653bd0f6816dd2`.
  Its sibling `.h5ad.meta.json` is the build record.
- Its extra settings are `VCC_NORMREL_EXP=-1.5`,
  `VCC_NORMREL_MODE=cap`, and `VCC_NORMREL_PM` **unset**, on the existing
  `cd4svd5` stack. These are not new default settings in the builder.
- **319d3 is built and preregistered, NOT submitted or armed**:
  `submissions/319d3-cap30.vcc`,
  SHA-256 `6ddbf76d2c7b914fae09f606f5da44262f25e415870552f211d606fa08c33d74`;
  only `normrel_exp` changes to −3.0.
  Read `experiments/319-cap-variant/PREREGISTRATION-dose2.md` and the
  September 9 coordinator review before deciding on it. A useful dose contrast
  is not, by itself, a causal test of why smaller donor rows were harmful.
  326 also measures resolved local raw-MSE harm: **+.0151952
  [+.0029296,+.0290382]**, masked at the H1 point score by its clamp.
  This is a concrete downside risk, not a calibrated A/B/C forecast.
- Source prediction `.h5ad` files were deleted in the September 8 sweep.
  Recover content with `scripts/restore_prediction_h5ad.py`, preserving the
  `.vcc` and build metadata. **Recovered content is not a promise of the original
  HDF5 file SHA.** See the storage audit before rerunning an old script.
- The fresh package/config audit is
  `docs/sprint-2026-09-09/package-verification-2005.json`. It explicitly
  resolves three CD4 peak-memory metadata differences and reconstructs
  each historically pinned builder SHA; it does not pretend the input
  metadata or historical/current builder bytes were identical.
- D/E/F driver defaults remain the separately rehearsed configuration:
  pinned `V5_abc` CD4 basis and carrying-row sizing. A/B/C gains do not
  automatically authorize adoption in October. Use `docs/DEF-RUNBOOK.md`.

## 3. Scientific interpretation that must not drift again

- **PDS is a cosine-based retrieval rank, not own-target cosine.** Higher PDS
  does not imply higher biological alignment. A fixed-gallery counterexample is
  reproduced by `scripts/audit_takeover_state.py`.
- Do not infer another team's cosine from our fitted MSE inversion. Positive
  `score_mse`, or similar raw MSE, does not establish matching emission scale,
  jackknife correction, or context weighting. Read public `score_*` columns
  directly for board gaps.
- **316 is norm-matched before whitening, not after emission.** 318 measured
  delivered amplitude redistribution in both 312 and 316. Their score difference
  cannot be called an isolated magnitude effect.
- **317/318 are negative for their tested constructions on the historical
  `rc_callset` H1 harness**, not a proof that every ordering/amplitude model on
  the current 319 base is impossible. Three shuffles provide a descriptive
  comparison, not a well-estimated null distribution.
- **322's frozen KOLF bank is panel-limited, not the upstream source.** The
  upstream atlas has 10,167 perturbations and was re-extracted in 322. The tested
  additions still failed; do not repeat them under a new “breadth” label.
- **F16 is a re-measurement, not a new biological discovery** (320/321).
  C9's small headline p is gene-sampling, not independent construct replication.
  The prior wet-lab-spend recommendation is withdrawn. R5 and the other used
  reserves stay spent; consult `docs/CONFIRMATION-LEDGER.md` before any new access.

Full reasoning, next experiments and stop conditions:
`docs/SCIENTIFIC-REVIEW-2026-09-09.md`.

### September 9 research continuation: experiment 323 complete

`experiments/323-current-family-validation/RESULTS.md` is the new primary record.
The current-family, union-OFF H1 benchmark now uses a correctly refitted operator:
the historical H1 adapter let 80 excluded H1 target coordinates enter its dense
whitener. A separate evaluator repair matches Scanpy's normalization/tie
arithmetic without changing production emission.

Construction gates and the historical six-metric check pass. Six actual-count
arms, H1-local anchors and 10,000 paired target bootstraps are complete.
**4/4 known board point signs match; 0/4 intervals resolve the small changes.**
The target-shuffle negative control passes. This H1 panel is **not a standalone
selector for small norm-weighting gains**; no new submission candidate is cleared.

**324 also complete:** metadata-only residual-label feasibility. GPSAdb's 50
panel targets spanning multiple accessions/lines contain **zero literal
two-accession × two-context rectangles**; only six targets have same-context,
same-method repeats. Publication/sample/construct independence is unresolved.
No response-label pilot launched. Next: a bounded official-metadata provenance
pass on the ten repeated accession pairs, particularly KHSRP, where distinct
ENCODE experiments might group into a genuine publication-level crossing.
Primary record: `experiments/324-residual-label-feasibility/RESULTS.md`.

**Sprint update:** 325 completed that official-metadata pass. KHSRP's apparent
publication crossing reuses ENCODE perturbation material; no independent-label
pilot was earned. SRSF5 retains a within-study distinct-construct lead.
326 completed the current-family cap30 diagnostic: local overall
**+.0001924 [-.0046193,+.0031757]**, unresolved, and raw MSE worsens despite
both arms' point MSE scores being clamped to zero on H1. Read each experiment's
`RESULTS.md`; neither clears an A/B/C candidate. 327's coordinate-median
correction is also terminal: unresolved local gain and no advantage over its
magnitude-only/shuffled controls. **328's source gates pass after a disclosed
numerical continuation**, while its original numerical NO-GO is preserved.
The fit converged at iteration 588, with unchanged controls. **329 then failed
its count-stage advancement gates**: local overall -.0056597 with an interval
spanning zero, resolved NMAE/fidelity harm, and no established overall advantage
over the fitted controls. Source predictive density is not a competition gain.
See 328's `RECOVERY-RESULTS.md`, 329's `RESULTS.md` and the sprint queue.
330's fixed non-panel source gene scales also completed: favorable raw point
estimates on all six members, but overall +.000835 [-.000302,+.002166] and
unresolved required control contrast. **No advancement.** No posterior or
scale construction is promoted; see 330's `RESULTS.md`. The sprint then continued
with distinct-method feasibility checks, not retuning these closed arms.
331's isolated CPU TabICL test is now **engineering-feasible only**: two
synthetic sanity workloads pass at approximately 1 GiB peak RSS. The root
scorer environment is unchanged. This is not biological evidence; read
331's `RESULTS.md` before reusing the model.
**332 is complete, no advancement:** absolute error improves but recipient
advantage is unresolved; cosine and retrieval worsen. All 24 fits and
independent endpoint/hash checks are complete. **333 acquired the GenePT
release, but its vectors are quarantined**: distinct primary descriptions
share identical vectors, including 34 challenge targets. Source mapping
inspection does not resolve assignment. Read both `RESULTS.md` records.
**334 is also complete, no advancement:** a direct-description TF-IDF/SVD
model with an explicitly frozen empty-description fallback loses retrieval
to source-only and the diagonal, and does not beat shuffled text.
49,152 endpoint rows and raw-source training features were independently
reconstructed. Read its `RESULTS.md`; no heavy job is active at the 15:05 UTC
checkpoint, and no new candidate is cleared.
**335/336 are now complete:** local biomedical-text encoding passes its
engineering gates, but the matched response model also loses retrieval,
including versus shuffled semantic features. All twelve source-only
replays are exact; 67,584 endpoint rows were checked. No text model is
promoted. See their result records and the sprint's next-method checkpoint.
**337 also complete:** availability-normalized consensus has H1-local
overall +.00322188 with a zero-spanning interval. Raw MSE improves, but
NMAE/fidelity/reach worsen. No advancement. The unchanged baseline and
independent scorer/dense-budget audit pass; no new candidate is cleared.
**338 also complete:** control-only fixed-rank histogram bounds fail the
useful-freedom gate in three of four targets. All 768 LP extrema and their
independent primal/dual certificates pass, so this is a measured negative
for the frozen relaxation, not a solver failure. No integer emitter or
accuracy test is earned. See its `RESULTS.md`.
**339 also complete:** zero-anchored gene-wise marginal calibration gains
.00293255 source retrieval, below its floor, with no advantage over the
diagonal/Gaussian controls and a cosine cost. No advancement. All twelve
known-answer prediction replays and the independent reconstruction pass.
**340 also complete:** sparse response-neighborhood regression improves
cosine/error but loses retrieval to the diagonal on all six lines, with
a resolved -.00520658 macro contrast. Sixty models, 43,008 endpoints and
independent uncertainty checks pass. No advancement or new candidate.
**341 also complete:** own-line training labels improve weighted sign by
.03254973 versus the cross-line diagonal, but retrieval loses .04829784.
Both retrieval safeguards fail despite passing sign gates. Independent
reconstruction of 36 fits, 36,864 endpoint rows and 37 intervals passes.
This is known-line held-target development, not unseen-line/cross-study
evidence. No advancement. **342 also complete, no advancement:** the
current-family Hamilton remeasurement has local overall -.00238441 with
a zero-spanning interval, no advantage over random allocation and resolved
fidelity harm. Exact floor replay, six-metric/dense-bootstrap checks and
independent reconstruction of all 255 total vectors pass. No candidate.

Ahmet subsequently authorized sole-coordinator continuation, and the existing
`VCC research watch` resumed this task hourly during the bounded sprint.
It is now **PAUSED**, verified in the final closure receipt. Its prompt,
cadence, target and notification preferences were preserved when pausing.
No built-in Codex subagents, automatic submissions, paid cloud/assays or new
reserve access were used. The stale a14/a15 rollback instruction had been
removed at takeover; this was never a reset-time submission queue.

## 4. Read only what the next action needs

| Task | Sources |
|---|---|
| Current scores / provenance | this page; `docs/BOARD-LEDGER.md`; `data/board-entries/`; `submissions/SUBMITTED_LEDGER.json` |
| Storage safety / recovery | `docs/STORAGE-AUDIT-2026-09-09.md`; `docs/storage-sweep-2026-09-08/`; restore helper |
| Plan research | September 9 scientific review; `docs/RESEARCH-PROGRAMME-2026-09.md`; relevant experiment `RESULTS.md` |
| Current-family H1 validation | experiment 323 `RESULTS.md`, registrations, input manifests and `work/results.json` |
| Next response-label feasibility | experiment 324 `RESULTS.md`, registration, and exact accession-pair CSVs |
| Biological claims | `docs/FINDINGS-2026-09.md`; `experiments/294-writeup/CLAIMS.md`; original results and coordinator reviews |
| Held data / confirmation | `docs/CONFIRMATION-LEDGER.md`; the exact preregistration, before accessing values |
| October build | `docs/DEF-RUNBOOK.md`; `scripts/def_pipeline.py`; experiments 235/240/247/251/256/311 |
| Metrics | `docs/reference/vcc2026-metrics-BRIEF.md`, then FULL for derivations |
| Historical reasoning | `docs/LOG.md`, historical `CURRENT_STATE.md`, archived READ-FIRST; never treat dates as live state |

## 5. Binding workflow

Read `AGENTS.md` §§2/4/6/7/11.1/14 for the full rules. In short: freeze the
hypothesis and controls before measurement; known-answer tests first; held
contexts excluded from fitting; one intended variable; inspect the **built
artifact**, not just the command; price or explicitly bound all six members.
Label measured, fitted, synthetic-fixture, historical, and proposed quantities.
Do not turn a failed construction into a biological impossibility claim.

One local heavy job at a time, with the documented memory gate (≥40% free,
free+inactive ≥20 GiB, swap free ≥400 MiB). Preserve other work. No blanket
staging or destructive cleanup. No re-submission to restore rank. Two slots/day
is a cap, not a target. Check live state and novelty before a separately
authorized submission. Public model name: opaque joke; description: empty.

Neither the audit nor experiment 323 changed the production model, uploaded a
submission, reopened reserve data, or published the manuscript. The initial
audit started no automation; the subsequent explicit takeover updated the
existing hourly continuation, which is now paused as recorded above.

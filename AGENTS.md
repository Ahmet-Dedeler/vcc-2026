# AGENTS.md — Virtual Cell Challenge 2026

Working rules for any agent (Claude Code, Codex, OpenCode, muse, Cursor) or human in this repo.
This file is stable rules only. What we have learned lives elsewhere:

- **Current state and findings:** [`docs/FINDINGS.md`](docs/FINDINGS.md) — current best entry,
  current conclusions, the closed-families list, and the full historical record (old §4.1x, §8.x,
  §12 live there under the same numbers).
- **Research programme:** [`docs/RESEARCH-PROGRAMME-2026-09.md`](docs/RESEARCH-PROGRAMME-2026-09.md).
- **Live slot plan and board record:** `docs/NEXT-SLOT-DECISION.md` (newest section at the bottom),
  `docs/BOARD-LEDGER.md` (newest at the top), `submissions/SUBMITTED_LEDGER.json`.
- **October D/E/F procedure:** `docs/DEF-RUNBOOK.md`.
- **Metric spec:** `docs/reference/vcc2026-metrics-BRIEF.md` (read before designing anything) and
  `vcc2026-metrics-FULL.md`.

Before you premise work on a fact about the data or a past result, grep `docs/` and the relevant
`experiments/*/RESULTS.md`. With 400+ experiments, the record knows more than you remember.

---

## 1. The task

Arc Institute ran CRISPRi Perturb-seq (10x Flex) on six human cancer cell lines. Three are the
validation panel (anonymised contexts **A, B, C** = Jurkat, HeLa, CAL33 by our CCLE match), three
are the final test panel (**D, E, F**, released Oct 22). For each context you get only the
non-targeting control cells plus 300 target genes, and you predict the perturbed single-cell
profiles. There is no challenge training set: this is zero-shot transfer of perturbation biology
from public data into a line you have never seen perturbed. The basal profile is the context's
only identity.

## 2. Hard numbers

| Quantity | Value |
|---|---|
| Contexts per round | 3 (A,B,C validation → D,E,F final) |
| Perturbations per context | 300, same set across a round's contexts |
| Cells per perturbation | exactly 400 |
| Rows in a submission | 360,000 (300 × 400 × 3) |
| Genes | 18,533, fixed order from `gene_names.csv` |
| Control cells per context | 18,400 (46 NTC guides × 400) |
| Max stored entries in `.X` | 4,750,000,000 (dense is 1.40× over: always sparse) |
| Max counts per cell | 1,000,000 |
| Submissions | 2 per day per team, reset 00:00 UTC, one in flight at a time |

`context_A.h5ad`: CSR float32 integer counts, median library 20,109 (range 3,275–52,420), 109.9M
nonzeros. Control pseudobulk Pearson A↔B 0.742, A↔C 0.685, B↔C 0.793: different tissues.

### 2.1 How the panel was built

Targets were drawn from an H1 pilot screen, binned by number of DE genes (negligible / subtle /
strong) and sampled for diverse response modes, with clean guide assignment, high basal target
expression, and adequate coverage in public datasets (hence K562-gwps covers 269/300). Two
consequences: `n_real` spans orders of magnitude by design, so sizing must be per perturbation, not
one global `n_pred`; and every target is expressed in the reference lines by selection, so a
missing expression/effect relationship is range restriction, not biology. Full text: FINDINGS H.1.

## 3. Timeline

| Date | Event |
|---|---|
| Aug 20, 2026 | Validation data and submissions opened |
| **Oct 22, 2026** | Final test set (D/E/F) released; leaderboard goes dark |
| **Nov 5, 2026, 23:59 UTC** | Final submission deadline, hard |
| Nov 20, 2026, 00:00 UTC | Winners announced at a virtual ceremony (Arc email 10-15) |

Only the last final entry counts for prizes ($100k / $50k / $25k). Leaderboard entries before
Oct 22 are never prize-eligible.

## 4. Scoring

Scored by `cell-eval2` with the `vcc2026` preset. Each metric is scaled `s = (u − b)/(r − b)`:
`b` = baseline (predict the context-mean perturbed profile for every target, score 0), `r` =
replicate (real data against itself over 5 half-splits, score 1). **Overall = unweighted mean of 6
metrics × 3 contexts.** Above 1 is possible; below 0 is normal.

| Board | `cell-eval2` name | Measures | `b` | `r` |
|---|---|---|---|---|
| `pds` | `pds_cosine` | Is the prediction nearest its own real effect by cosine, among 300 | 0.500 | 0.927–0.984 |
| `mse` | `expr_mse_unbiased_capped_norm` | Jackknife-corrected squared profile error / real effect size | 0.986–0.992 | 0.028–0.045 |
| `fid` | `de_wilcoxon_direction_fidelity_yield_raw` | Sign accuracy of your DE calls × `min(1, n_pred/n_real)` | 0.505–0.522 | 0.795–0.832 |
| `reach` | `de_wilcoxon_direction_reach_raw` | How deep your ranking of the reference's DE genes stays ≥90% sign-pure | 0.047–0.097 | 0.958–0.978 |
| `jac` | `de_wilcoxon_sig_jaccard` | Jaccard of DE gene sets | 0.021–0.037 | 0.375–0.423 |
| `nmae` | `de_wilcoxon_lfc_nmae` | Normalised MAE of log2FC on the reference's DE genes (lower is better) | ~1.001 | 0.369–0.431 |

Clamps: `mse` to [0,1], `nmae` floored at −6.

What falls out of the spec:

1. The target gene's own row is excluded everywhere; `pds` also drops all 300 panel genes. The
   whole signal is the downstream response.
2. Doing nothing (control copy) is the baseline by construction.
3. `pds` is exactly scale-invariant; `mse` is quadratic in magnitude. `mse` also reads the bulk,
   so the bulk (pds, mse: log-pseudobulk of summed counts) and the DE channel (fid, reach, jac,
   nmae: per-cell arithmetic CPM and Wilcoxon ranks) are different functionals of the same counts
   and can be designed separately (see FINDINGS §1).
4. DE is recomputed from your counts: Wilcoxon against the real controls, BH α=0.05, genes > 5 CPM
   in control. Your within-perturbation cell-to-cell variance decides significance, so emission is
   a modelling decision. 400 identical cells make everything significant and destroy `jac`.
5. `fid` penalises under-calling; over-calling costs `jac`. Call-set size is a real tradeoff.
6. `mse`'s jackknife correction cannot be farmed by emitting noise.

## 5. Data

`data/challenge/` (matrices gitignored, the three small CSV/JSON files tracked; refresh with `vcc datasets download controls`): `context_{A,B,C}.h5ad`
(18,400 × 18,533 raw counts; `.obs` has `target_gene`=`non-targeting`, `context`, `ntc_id`),
`gene_names.csv`, `pert_counts.csv`, `manifest.json`. Both CSVs have a header row (`gene_name`,
`target_gene`); read by column name, and always read gene order from the file, never from docs.

External data (links in FINDINGS H.2b): 2025 VCC H1 (same chemistry, with the held-out test
split; `scripts/download_2025_h1.sh`), Replogle K562/RPE1, Nadig HepG2/Jurkat, X-Atlas/Orion
HCT116/HEK293T, Jiang 2025. Direction quality is limited by the number of independent lines.

## 6. Submission mechanics

### 6.1 Format contract (all enforced)

- One `.h5ad` for all three contexts; `.obs['context']` uses the control files' labels.
- `.obs['target_gene']` = gene symbols, exactly the 300 in `pert_counts.csv`, 400 cells each.
- No `non-targeting` rows. Controls are inputs only.
- `.var` index = all 18,533 genes in `gene_names.csv` order.
- `.X` = raw non-negative integer counts (no normalisation, no log), sparse CSR without explicit
  zeros, no cell above 1,000,000 counts.

### 6.2 Commands and slot rules

```bash
vcc prep prediction.h5ad -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv -o submissions/prediction.vcc --dry-run
vcc submit submissions/prediction.vcc -m "Towards a unified theory of meh" --wait
```

Always `--dry-run` first. `vcc status <id>`, `vcc whoami`, `vcc sample` (random valid file for smoke
tests). CLI guide: https://vcc-cli-wiki.virtualcellchallenge.org/ · skill `/vcc`; re-run
`vcc skill install` after CLI upgrades. `-m` and `-d` are public: see §14.

Production builds go through `scripts/make_submission_pm_cd4.py` with env knobs (`VCC_*`); the
current recipe is in the best entry's `experiments/<NNN>/` and its `*.meta.json`. Unset knobs are
bit-identical to the previous default; prove that with a test, never assert it. Arms are armed with
`scripts/arm_submit.sh` / `scripts/arm_submit_after.sh`, which run `check_submission_novelty.py` and
record `submissions/SUBMITTED_LEDGER.json`.

**Slot rules (Ahmet's standing rules):**

- **Never spend a slot on rank, and never re-submit an already-scored package.** If a bad arm drops
  the displayed row, leave it. The board shows a team's latest entry, not its best; only the last
  entry after Oct 22 matters. The goal is first place on Nov 5, not a better daily row.
- **Do not hoard slots either.** Every slot before Oct 22 expires at midnight and is prize-irrelevant.
  An unspent slot teaches nothing. Each day, ask what a slot would teach, write the reading and its
  bins down before the entry exists, then submit the best available arm. Losing arms are fine; tell
  Ahmet in advance when the displayed row will fall. Two slots on one preregistered dose curve beat
  two unrelated guesses.
- Before saying "nothing to submit", list `submissions/*.vcc` against the ledger and read each
  unsubmitted package's registration.
- Uploads take 11–19 min for a 2.5–2.8 GB package. Start by T−25 min before reset.
- Capture every board record immediately (`scripts/capture_board_entry.py`): the API only resolves a
  team's current entry, so the previous one's members vanish when the next publishes.

**Waiters:**

- A waiter for tomorrow's slot gets `FLOOR` = the next `YYYYmmdd0000` UTC, never "now". An allowance
  refusal is terminal in `arm_submit_after.sh`; there is no retry after reset.
- A waiter's wait loops must be capped by the target timestamp plus margin, never by a fixed iteration count. On
  2026-09-23 an 8 h cap (1,440 × 20 s) on a waiter armed 23 h before the reset would have exited silently before its
  slot.
- Sessions die without warning (safeguard blocks, relogins, reboots). The moment a package passes
  its gates, leave a waiter that can finish alone, and write state to the repo before starting the
  next thread. After any gap, check `uptime`, `ps` for the waiter pid, and the tail of its log.
- **A sleeping Mac pauses every waiter and drops uploads** (2026-09-30: on battery the Mac slept ~02:40–05:29Z; the 458d
  upload died with a connection reset and the waiter's retry sat in `sleep 90` for two hours). `caffeinate -s` holds only on
  AC, and a closed lid sleeps anyway: keep the Mac on AC with the lid open through slot windows and Oct 22 – Nov 5. Waiters
  from arm1001 on retry up to 6 times before their cap and resume an interrupted upload with `vcc submit --resume <entry>`
  (`experiments/536-waiter-retry/patch_retry536.py`); copy that block into any new waiter.
- **Launch every long-lived job with `nohup bash scripts/detach.sh <cmd> … > log 2>&1 &`** (setsid: a new session). On 2026-10-15
  every plain `nohup` waiter, caffeinate and snapshot loop started from an agent shell at 13:33Z was dead by 21:19Z with NO reboot
  (uptime 8 h): quitting or restarting the app signals its whole process group, and nohup only ignores SIGHUP.
- **A reboot kills every waiter, chain, worker and watcher** (2026-09-23 10:39Z). Recovery order:
  1. builder and emitter shas against the last commit. A build killed after its pin stage leaves the builder pinned:
     restore with `scripts/def_pin_panel_artifacts.py --restore --record <build>/pin-artifacts.json`;
  2. re-arm the day's waiters (same scripts, same arguments; they re-read their logs);
  3. `scripts/board_snapshot.sh` and `scripts/pisces_watch.sh`, plus `caffeinate -is -t <s>` until the slots are done;
  4. move an interrupted build dir aside (`build_X.aborted-…`) and rebuild;
  5. relaunch workers with a RESUME note in their mandate. Copy `claude-run.jsonl` aside first: the launcher
     truncates it.
  6. 10-21 → 10-23: re-launch `experiments/691-release-watcher/watch691.sh` (AUTO_R0=1) with the same env. It keeps its
     baseline.json and refuses an existing data/challenge-def.
- Waiters run under bash with macOS BSD `/usr/bin/grep`. In the agent's own shell, `grep` may be a GNU-style wrapper, so a
  pattern that matches there can silently fail in the waiter. BSD grep ignores a `$` anchor placed before `\|`: the
  pattern `" END$\|stopping"` never matches " END" (2026-09-25: two waiters would have hung). Use `grep -E " END$|stopping"`,
  and test a new waiter's conditions with `bash -c '…'`.
- `bash` here is macOS `/bin/bash` 3.2, which has no associative arrays. `declare -A` errors, and then `X[key]=v` silently
  writes index 0, so every keyed lookup returns the last value (caught in `456/chain0928.sh`, 2026-09-26). Use `case`
  functions, and zsh's `${=var}` or a bash loop when a tool-shell test needs word splitting.
- A log that quotes another log's tail can contain that log's `... END` line. On 2026-09-26 chain456 waited for `" END$"`
  in chain455.log, and 455 had logged 453's `...02:31:12Z END` on its first line, so 456 started beside 455's build. Its
  run_arm gate parked it at "other builders 3" and it was killed before pinning anything. Wait for a line only the
  watched script writes last, AND for its process to exit.
- A waiter must fire when the watched chain's PROCESS exits, then read the log's tail. Never require a success-only END
  line: on 2026-09-27 chain475 failed and exited on an error path that wrote no bare `END`, so both its waiter and the next
  chain's wait sat idle for 2 h. Every chain's failure path also ends with a terminal line (`<ts> END`).
- Never edit a script file that a live bash waiter is executing (bash reads scripts incrementally).
  Change what it calls later, or kill it via its `.pid` and re-arm with the same arguments. Python
  helpers it calls are read fresh at start and are safe to edit.
- **`pgrep -f <name>` matches any process whose arguments contain the string**, including your own shell command and
  every headless `claude -p` worker whose mandate text mentions it. On 2026-09-23 a worker's build chain deadlocked
  waiting for "no process matching chain433k_resume.sh", because its own prompt named that script. Match the
  interpreter plus the path instead (`ps -axo comm,args | awk '$1 ~ /^(bash|python[0-9.]*)$/ && /<script>/'`). Don't
  put a gate's literal pattern into a mandate. The same goes for any long-lived helper's argv: on 2026-09-27 two
  committers carrying `run_arm_generic.sh` and `def_day.sh` as path arguments would have deadlocked the build queue.
  Helpers read their paths from a spec file (`scripts/commit_after_chain.sh <spec>`). The same holds for a chain's own FILE NAME: a
  chain that waits at the end of a queue must not be named like the builds the earlier gates look for (`chain49[5-8].sh`);
  on 2026-09-28 a waiting `relevel3/chain498.sh` held chain497e, which it was itself waiting behind.
- A known-answer run passes only if the driver exited 0 AND the output hash matches. On 2026-09-27 479's KA-T
  matched 446bj2's hash while the driver's arm assert had raised (rc 1). Every per-target build would have died
  after it; a cloud review (499) caught it.

### 6.3 Use the writer

`scripts/build_submission.py` owns the format contract and streams CSR into HDF5 (a full
submission is ~2e9 nonzeros; never materialise it). `ChallengeSpec.load`, `write_submission`,
`check_submission`. Emission dispersion is a real hyperparameter.

### 6.4 Local and cloud scoring

`cell-eval2` reproduces the scorer offline on any surface with ground truth (in practice
`data/harness/h1-test.h5ad`). The correct four steps:

```bash
uv run cell-eval2 baseline -ar real.h5ad -o out/baseline --save-pred out/baseline-pred.h5ad --preset vcc2026
uv run cell-eval2 prep-real-bundle --real real.h5ad --baseline out/baseline-pred.h5ad -o out/bundle --preset vcc2026
uv run cell-eval2 run -ap pred.h5ad -ar real.h5ad -o out/run --preset vcc2026
uv run cell-eval2 score --user-agg out/run/agg_results.csv --real-bundle out/bundle -o out/score.csv
```

`prep-real-bundle --baseline` takes the baseline **prediction h5ad**, not the `baseline -o`
directory; passing the directory fails with `IsADirectoryError` after hours of `baseline`. For
h1-test, `data/harness/h1-test-baseline-pred.h5ad` already exists and skips that step.

- **Locally:** at most one `cell-eval2` process at a time (~25 GB RSS each). Use
  `scripts/run_with_memory_guard.py`.
- **Cloud:** `VCC_CLOUD_SCORE_CLEARANCE=YES scripts/cloud_score.sh --real … --pred name=path --out … --wait`
  runs the full pipeline on a disposable EC2 worker and terminates it: about $1.50 and 2.5 h per
  six-member score, no local memory pressure, known-answer row reproduces 095's h1-test numbers.
  `cell-eval2` is single-threaded here, so pick the cheapest instance with ≥64 GiB; avoid spot for
  anything over ~30 min (a known-answer run was reclaimed 42 min in). `scripts/aws_reap.sh` sweeps
  strays. "We cannot run the real scorer" is not a reason to price members through proxies.

### 6.5 Harness caveats (details in FINDINGS H.4)

- `h1-test` and, since 2026-10-02, 618's `h1-trainval` real-emitter surface (126 unworn targets; trust gate PASS; call-set
  rows worn) are the surfaces validated against the board. Orion thinned harnesses carry a
  technical blob in their truth: relative direction only, never quote absolutes from them. The blob is
  an artifact of their builder's control selection (615, 2026-10-02). Raw-cell truths with pooled
  controls carry h1-test's level and are usable once they pass their own trust gate.
  `data/harness/orion-hct116.h5ad` is depth-confounded (`docs/HARNESS-DEPTH-CONFOUND.md`).
- Never score h1-test with `data/donors/h1_2025` in the donor mix: it covers all of h1-test's
  targets (oracle leakage). Audit donor↔surface target overlap before every screen.
- To power-match `reach`, scale each perturbation's pool proportionally; a uniform
  `n_real_target` deletes the small-`n_real` stratum where most reach lives.
- When validating whitening on a historical panel, exclude that reference's own targets from the
  operator.
- H1 real-emitter pricing must match production's zero repair (`run496 --zero-repair`; default OFF there). Without it a lever
  that creates zero-total genes reads falsely bad, ~+31 error per collapse vs production's −5 floor (654, 2026-10-02). Count
  cells with lfc < −30 (`experiments/654-carrier-amplitude/zerocount654.py`) on any DOWN-magnitude price.
- Any screen whose arm family already has a board point must retrodict it as its known-answer row
  before its numbers are read.

## 7. Rules that constrain modelling

Full rules: https://virtualcellchallenge.org/app/rules

- Machine-learning predictions only. Lab results or literature values may train or fine-tune a
  model but may not be submitted or used to revise output. A hand-curated DE lookup table is out;
  a model trained on public Perturb-seq is in. Frame the donor-transfer stack as a fitted model in
  the final description.
- Proprietary models and data are fine if we have the right to use them.
- Code submission is not required; finalists publish a high-level method description.
- Arc's State model is allowed but "not able to win"; it is also closed here (pds ≈ 0.50 on A/B/C).
- Team: "Ahmet Dedeler" (solo), account `<redacted-email>`.

## 8. Findings (and old §12)

Moved to [`docs/FINDINGS.md`](docs/FINDINGS.md). Read its §0–§2 (current state, current
conclusions, closed families) at the start of any session that will design or price an arm.

## 9. Environment

Apple M5 Max, 18 cores, 64 GB RAM. `uv sync`; run everything with `uv run python …` (env in `.venv`).
`cell-eval2` is in the venv; `vcc` is a global tool (`~/.local/bin/vcc`, keyring auth). `torch` has
MPS; `cell-eval2[gpu]` does not help on this Mac. Stream large h5ad files with `backed='r'` and
precompute pseudobulk/DE tables instead of reloading cell matrices.

**Memory gate.** Proceed with a heavy job when `memory_pressure` reports ≥ 40% free AND `vm_stat`
free+inactive ≥ 20 GiB AND swap free ≥ 400 MiB (a swap TOTAL of 0, as after a reboot, passes: nothing is swapped). `vm.swapusage` free alone is not a pressure signal
(the swap file grows on demand; a "< 2 GB swap free" rule stalls forever). Run heavy jobs strictly
serially, report peak RSS with `/usr/bin/time -l`, and if you deviate from a mandate's memory rule,
say so in RESULTS.md. A job that runs out of memory here dies without a traceback: re-run it with
headroom before debugging it as code.

**Shell.** The Bash tool's shell is zsh, which does not word-split unquoted variables: `for x in $VAR`
runs once. Glob directly (`for p in experiments/2*/x; do`) or use `${(f)VAR}`.

## 10. Repo layout, git, storage

```
AGENTS.md            this file
docs/                FINDINGS.md, runbooks, ledgers; docs/reference/ = vendored specs
data/challenge/      official controls (h5ad ignored) data/external/  public datasets (gitignored)
data/donors/         donor tables and whitening operators   data/harness/  local scoring surfaces
data/leaderboard/    board snapshots (scripts/board_snapshot.sh)   data/board-entries/  captured records
src/vcc2026/         model code          scripts/  builders, waiters, scoring, downloads
experiments/NNN-name/  one dir per experiment: MANDATE/PLAN, PREREGISTRATION, RESULTS.md, work/
submissions/         .vcc + .h5ad ignored; ledger, *.meta.json, logs tracked
```

Every experiment writes `experiments/<name>/RESULTS.md`: config, scaled scores per member, and one
honest sentence on what changed. Negative results get written down with their mechanism. Nothing is
submitted without a local score first.

**Git.**
- Stage explicit paths (`git add docs/ experiments/<mine>/ AGENTS.md`). Never `git add -A` while
  other agents run: sweeps capture other agents' in-progress edits and mislabel history.
- **All sessions share ONE git index.** `git add X && git commit` also commits whatever another session
  has staged at that moment. On 2026-10-02 two commits carried other sessions' files: 54ea28f3 and 4cd40535.
  Always commit with a pathspec: `git commit -m "…" -- <paths>`. It records only the named paths' working-tree
  content and leaves everyone else's staged changes alone. For a new or ignored file, `git add -f <path>`
  first, then commit with the same pathspec.
- HEAD is not a stable baseline during a session. For a transient edit to a shared file, save a byte
  copy first and restore from it, not from HEAD; compare bit-identity proofs against a pinned rev.
  Sweeps can also clobber a live file, so copy long-lived deliverables to your scratchpad after
  substantive edits and re-check them before declaring done.
- Edit the production builder only after the running build logs "builder after restore". An agent
  that changes production code commits it as part of "done", staging the path explicitly and
  putting its known-answer hashes in the message.
- Cloud agents (Cursor etc.) fork GitHub, so anything not pushed does not exist for them. Commit and push each
  experiment's small outputs (RESULTS, logs, result json/csv/txt, board snapshots, fitted arrays the builder reads) when
  it finishes: any file under ~5 MB, `git add -f` an ignored array under ~50 MB. Bigger data stays local.
- Never `git add` a file over 100 MB. GitHub rejects the push and it can fail silently for weeks.
  After pushing, check `git rev-list --objects origin/master..master | git cat-file --batch-check`.
  Never push `backup/pre-filter-2026-09-08`.

**Storage.** Each arm writes a ~12 GB source h5ad and a ~2.8 GB `.vcc`.
- **After the final (2026-11-06):** `data/external/`, experiment arrays and old `.vcc` packages were deleted (`experiments/723-post-final-cleanup/`). Any rebuild
  or re-run needs the public data re-downloaded first (links in FINDINGS H.2b). The rules below describe the competition period.
- Watch the data volume: `df -h /System/Volumes/Data` (the `/` line is the system volume). It reached 99 % on 2026-09-30.
  Start a build only with ≥ 35 GiB free: the h5ad takes 12 GB and `run_arm_generic.sh`'s prep gate then needs ≥ 20 GiB.
  Keep ≥ 60 GiB free through the final phase.
- Delete a prediction h5ad only after its `.vcc` exists and its sha is logged, and make every chain's cleanup step
  conditional on that. On 2026-09-30 an unconditional `rm -f` deleted the only copy of an arm whose prep gate never opened.
- Once a `.vcc` exists, its source `.h5ad` is redundant: `scripts/restore_prediction_h5ad.py <vcc> <out>`
  restores the content (not necessarily the original byte layout or file sha). Sweep `pred.h5ad`
  once its scores are recorded.
- Once an arm is scored and its numbers are in RESULTS.md, clean its build tree with
  `git clean -xdf experiments/<arm>/work/build*`, never `rm -rf`: tracked provenance files
  (`*.meta.json`, `SHA256SUMS.txt`, `postbuild.json`, `frozen_inputs.json`, `cd4_compact/`) live
  inside. Leave anything touched in the last 7 days alone.
- The anti-re-submission record is `submissions/SUBMITTED_LEDGER.json` (git-tracked, with content
  hashes), not the `.vcc` files; it is not a complete history, so absence from it does not prove a
  package was never submitted. Keep every `*.meta.json` sidecar.
- Keep `199/work/GWCD4i.DE_stats.h5ad`, known-answer rehearsal outputs and `data/external/`.
  Duplicates become APFS clones (`cp -c`), never symlinks. History: `docs/STORAGE-AUDIT-2026-09-09.md`.

## 11. How to work

### 11.0 Roles

Ahmet funds the project; the agent running it is the lead scientist. He is not technical in this
area and has said "the calls and shots are urs." Make the scientific decisions and execute: no
opening clarifying questions, no "shall I?", no menus of options. Any method is allowed (new
architectures, training, papers, compute). The goal is first place, and 0.25+ overall. Beyond
that, he wants scientifically useful results, not only a rank; think long-term. Write findings into
the repo (`docs/FINDINGS.md`, the experiment's RESULTS.md), not only into chat: sessions move
between runtimes and the repo is the only shared memory. He reads summaries, so lead with the
number and the decision.

**Standing permissions: never ask first** (Ahmet, 2026-10-02: "dont ask for approval on everything").
- **Downloads.** Download any data, model, paper or supplement you need, including the D/E/F bundle on
  Oct 22, as long as `df -h /System/Volumes/Data` still shows ≥ 60 GiB free afterwards (§10's
  final-phase floor). Record each download's path, source and size in the experiment's RESULTS.md.
- **Upgrades.** Upgrade any library, CLI or tool as needed: `vcc` (then `vcc skill install`),
  `cell-eval2`, uv packages, Homebrew tools.

Two technical guards, which are not approval steps:
- Prefer official sources, and never execute a downloaded binary from an untrusted one.
- After an upgrade that can touch the builder's or the scorer's numerics (numpy, scipy, anndata,
  h5py, cell-eval2), run the all-unset known-answer rebuild before the next production build.

### 11.1 Method rules

Hold weak priors and strong measurements: any idea may be tested, none believed without a number.
Cheap, "naive" ideas have repeatedly beaten sophisticated ones here. The numbering below is cited
elsewhere in the repo; worked examples are in FINDINGS H.5.

1. Negative results are results. Record what failed and why, or the next agent re-runs it.
2. A metric you cannot move must not drive a shared parameter. Check reachability before optimising.
3. Validate the validator: a harness is trusted only after it retrodicts a known board result.
4. Score the thing you are ranked on. Proxy objectives drift (mean cosine is not `pds`).
5. Change one thing per run.
6. Differences inside seed noise (0.2–4.8% of span) are not findings.
7. Local runs are unlimited; slots are two a day. A slot must buy information (§6.2).
8. Never average across an axis the metric refuses to average across. DE members are scored per
   perturbation; pooled statistics answer a different, easier question.
9. Any metric with `n_real` in its denominator or union (`fid`, `jac`, `reach`'s normalisation) is
   unreadable on our harnesses; only the ordering-driven parts transfer.
10. Diff the finished artifact's recorded config against the incumbent's and require exactly the
    intended number of differences. Assert the operator sha. Every screen carries a configuration
    whose answer you already know and aborts if it does not reproduce.
11. Record a path and a sha, never a description, and never infer an operator from a filename.
    Several whitening operators with near-identical names exist. Un-whiten with the operator that
    whitened.
12. Price every arm on all six members, or state which are unmeasured and the mechanism that pins
    them.
13. Before a slot, price the arm in units of the noise on its channel. Direction members (pds, mse)
    reproduce on the board to 16 significant figures, so small direction gains are readable. DE-member
    noise depends on the arm's class (537 §2.5, CONFIRMED 2026-10-03 on 12 board pairs):
    - an arm that keeps the emission LAW but re-draws pass 1 (any bulk change that moves ballast K,
      any K or seed change) carries a re-roll sd ≈ 0.00025 overall (reach rms 0.0013 raw, fid/jac
      ≲ 1e-4 raw);
    - a prior-only DE arm that keeps pass 1 reads near-exactly (≈ 1e-5 overall);
    - only a change of emission LAW carries 230's ≈ 0.0045 spread.

    A gain must clear ~2 sd of its own class. There is no validated `reach` toll; do not charge one.
    Re-price every rider against the current incumbent.

Also: verify substance, not narration. A subagent's (or your own) timings, effort and "I ran N
draws" are unverified unless a log shows them. Re-run its known-answer rows yourself before acting
on a result. Before reporting a ratio, name the units on both sides: local instruments report raw
statistics, the board reports scaled members. Before calling a property verified, say at which
pipeline stage it was measured. Borrow freely from other fields (retrieval, econometrics, signal
processing); if an analogy suggests an experiment, run it.

### 11.1b Never idle while work runs

An armed monitor is not permission to stop. While background jobs or workers run, keep launching
independent work: more workers on disjoint mandates, data and literature hunting, long cloud jobs,
analysis of what is already on disk. Nights and any time Ahmet is away are for the slow, expensive
things. Keep a queue of the next three experiments. Fill the machine up to the memory gate. Never
report "waiting" as a status; say what you are doing with the interval.

Supervise running workers. After a session restart or relogin, detached workers can stall mid tool
call with a live pid: judge liveness by log mtime (older than ~10 min means stalled), kill and
relaunch, redo any gating checks they owned yourself, and re-arm monitors first.
`scripts/agent_supervisor.sh` relaunches runs that died on transient provider errors;
`scripts/heartbeat_subagents.sh` is the session heartbeat.

### 11.1c The final is a different exam

The prize is scored on D/E/F: three new cell lines, a new 300-target panel, a dark board. An A/B/C board gain says a
change worked on these three lines and this panel. It says the change will work on D/E/F only as far as its mechanism
carries. Ahmet (2026-09-29): gains must come from solving the task, not from fitting this round.

1. **Say why a change helps.** Every PREREGISTRATION and RESULTS puts the change in one of three classes:
   - **Transfer:** a better prediction of the responses (donors, weights, keys, signs, shape). It must also gain on data
     the model never saw, before it enters the October recipe: a held-out line (a PLOLO fold, H1 without H1 donors) or
     an October-like panel. The board alone does not show transfer.
   - **Scorer mechanics:** how the counts are written for this scorer (emission, ballast, re-level, MDE, delivery
     fixes). It carries to D/E/F because the scorer does not change, but only if nothing in it is fitted to A/B/C.
     Re-fit any fitted part from D/E/F's own controls on the day, as the re-level's NTC proxy does. It does not make
     the model better at biology; say so when reporting it.
   - **A/B/C-specific:** tuned per context or per target on board reads, or valid only for A/B/C's lines. It never
     enters the October recipe.
2. **An October switch needs two readings:** its A/B/C board read, and a held-out read (instrument, number, path).
   `docs/october-decisions.json` records both. A board-only switch ships only with a written reason. 438 is the
   model: CD4 ×1.5 won on A/B/C, lost on every October-like panel, and October ships ×1.
3. **Never mine the board for A/B/C's answers.** Per-target or per-context probes, and dose hill-climbs on board reads
   without a mechanism priced first on held-out data, can raise the displayed row but cannot reach D/E/F. D/E/F
   feedback after Oct 22 is different: those reads are the final's own truth (`docs/FINAL-PHASE-FEEDBACK-PLAN.md`).
4. **Count the tries.** Every adopted change states how many arms were tried for its switch, on the board and on each
   local instrument. The best of N noisy reads is biased upward: the best of 5 identical arms reads about +1.2 sd by
   luck alone. That is ≈ +0.0003 overall at the re-roll sd 0.00025 (same emission law) and ≈ +0.005 for changes of
   emission law (sd 0.0045). Rule 13 lists the classes.
5. **Instruments wear out.** A surface that has chosen many settings (h1-test, the PLOLO folds, 438's `real` / `efp` /
   `efr`) is tuned-on, so its later readings are optimistic. Building fresh test sets from public data (held-out lines,
   new Arc-like panel draws) is the right fix. Before the final recipe freezes, confirm it once on a sealed reserve
   that nobody has tuned on (`docs/CONFIRMATION-LEDGER.md`; the current one is 534's R6).
6. **Report both numbers.** In a progress summary for Ahmet, put the held-out or sealed number next to the board
   number, and name the class of change that moved it.

### 11.2 Subagents and workers

- **Ask Ahmet which runner to use before dispatching workers** (native Claude subagents, the `muse`
  CLI, Antigravity `agy`, `cursor-agent`, Codex): he picks by credits. One short question, then go.
- **Always the newest, strongest model at high/xhigh effort.** List the CLI's models first
  (`agy models`, the vendor list) and pick the newest; never pick from memory. **Never Sonnet** for
  any task here. If a task is mechanical, do it yourself instead of spawning a cheap agent.
- `muse` (Muse Spark 1.3 Contributor, `--reasoning-effort xhigh`, which is its ceiling) runs headless
  with tools only with `--yolo`: `muse exec --json --workspace <repo> --model muse-spark-1.3-contributor
  --reasoning-effort xhigh --user-input-auto-resolve --yolo --prompt-file <file>`. Launcher:
  `docs/lab-2026-09-13/launch_muse.sh`.
- Claude Code workers: `bash scripts/claude_subagent.sh <exp-dir> <mandate.md>` (run with `bash`, not
  `source`; zsh has no `BASH_SOURCE`). It logs to `experiments/<exp>/claude-run.jsonl`. Use
  `CLAUDE_SUBAGENT_MODEL=opus` for biology- or dataset-name-heavy mandates: the default model's
  safety classifier false-positives on them (`[bio]`). A printed PID is not proof of startup: within
  the first minute check that the log is growing with nonzero tool calls and no early `result` event.
- Workers write only in their own experiment dir: no git, no `vcc`, no builder edits, no
  `cell-eval2` unless the mandate says so. The coordinator commits.

### 11.3 Fan-out

Fan out for breadth, not repetition. Give every worker territory nobody else has (by dataset, line,
metric, method family or source) and say what the others cover. Two stages work: fan out to find
(concrete inventories: URLs, sizes, lines, licence, why it could help), then hand each find to a
different worker to wire into the harness and return a scored number. Every mandate states the
objective, the harness that produces the number, the incumbent number, its territory, what others
do, and what to return. Adversarial review by a different model is cheap and has paid off.

Where to look: `gh search`, `hf datasets ls --search … --limit 50` (`ls`, not `list`; default
limit 10), arXiv/bioRxiv supplements, GEO, Zenodo, scPerturb, competitors' repos and model strings.

## 13. Sources

- Challenge: https://virtualcellchallenge.org/app · evaluation https://virtualcellchallenge.org/evaluation
  · rules https://virtualcellchallenge.org/app/rules · leaderboard JSON
  `curl https://virtualcellchallenge.org/api/leaderboard` (raw and scaled members, model strings)
- Scorer: https://github.com/ArcInstitute/cell-eval2 · CLI guide https://vcc-cli-wiki.virtualcellchallenge.org/
- Support: help@virtualcellchallenge.org · Discord https://discord.com/invite/f2aWGPXwej
- Full source list: FINDINGS History §13.

## 14. Naming discipline and the description field

Ahmet's standing rule: nothing about our method goes on the leaderboard, ever. Both the model name
and the `-d/--description` field are public and permanent, and we have read competitors' whole
methods out of theirs.

1. **`-m` is a random joke**: a fake paper title or misdirection ("Attention is not all you
   needed", "the State model (lying)", "the houseplant filed its taxes"). No method, donor,
   coefficient, emission, gate, env var, version token, or number that could be α / top-k / CPM /
   γ / φ. Don't reuse a name. **Don't invent one either:** take the next name from the ordered
   queue in [`docs/MODEL-NAMES.md`](docs/MODEL-NAMES.md) with
   `uv run python scripts/next_model_name.py --claim <arm>` when you bake it into a chain.
2. **`-d` is left empty.** Not vague, empty.
3. Never confirm or deny another team's guess, and never publish a number that lets someone invert
   one of ours.
4. The real config goes in `experiments/<name>/RESULTS.md` and the `*.meta.json` sidecar. The method
   is described in the finalists' write-up after the final test set closes.
5. This applies to every surface: Discord, GitHub, X, issues on `ArcInstitute/*`, public repos,
   preprints before Nov 5.

v1's model name and the a21/a22 descriptions already leaked recipes; do not add to it.

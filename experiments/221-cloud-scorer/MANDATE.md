# 221 — Stand up a cloud `cell-eval2` worker. It fixes the root cause of §8.4's 4.5× error.

You are an infrastructure + measurement subagent. Work in `experiments/221-cloud-scorer/` and
`scripts/` (a new reusable launcher is part of the deliverable). AWS is authorised for this mandate.
No `vcc submit`, no challenge submission build. Local peak RSS ≤ 8 GB — **the heavy work goes to
the cloud, which is the point.** `uv run python`.

## READ FIRST
1. `AGENTS.md` **§8.4** — why this matters — plus §6.4 and the **memory-and-validation safety note**,
   which says in as many words: *"never run more than one `cell-eval2` process locally. A scorer can
   consume about 25 GB RSS… Use `scripts/run_with_memory_guard.py`, **or score on an isolated cloud
   worker.**"* That sentence has been in the file since 2026-08-22 and nobody has built the worker.
2. `experiments/061-perturbmatch-hct-count/aws/run.sh` and `experiments/082-biological-target-kernel/aws/`
   — **the existing AWS tooling. Reuse it; do not reinvent it.** 082's run was guarded, hashed and
   recovered cleanly (27 output hashes verified, worker stopped after recovery) — copy that discipline.
3. `experiments/210-cd4-build/BOARD-RESULT.md` — the board point that exposed the problem.
4. `AGENTS.md` §11.1 rule 9 (`fid`, `jac` and `reach`'s normalisation are unreadable on our
   harnesses) and rule 12 (price all six).

## Why this, and why now

§8.4: `cd4f034` predicted +0.00437 and delivered **+0.0130557**, and **60% of the realised gain came
from the four members we carried at zero** — `nmae`, `fid`, `reach`, `jac`. We carry them at zero
because we have never had a way to measure all six members of a candidate build together: the real
scorer needs ~25 GB and the standing rule forbids running it alongside anything else, so in practice
it is never run at all and we price members piecemeal through proxies.

**That is a tooling constraint producing a scientific error, and it has now cost us a measured
4.5× under-pricing on the only channel that matters.** There are ~96 slots left before Oct 22 and
Ahmet has AWS credits sitting idle. This is the right thing to spend them on.

## Work

### Phase 0 — housekeeping first, and make it permanent
An instance (`i-0af19b9a537e08d4b`, `r7i.2xlarge`, us-west-2) was left **running for 13 days** after
experiment 061 finished on 2026-08-22 — CPU max 1.9%, average 0.27%. I have stopped it. **Your first
deliverable is that this cannot recur:**
- a `scripts/aws_reap.sh` that lists every running instance across regions with its age and its
  `Name` tag, and stops anything older than a threshold whose CPU has been under 5% for 24 h;
- every launcher you write must register a **self-terminating** guard — a shutdown timer set at boot
  (`shutdown -h +N`) — so a worker dies on its own even if the coordinator crashes;
- report S3 state too: `s3://<your-s3-bucket>` currently holds **629 objects / 81.5 GB**
  from the 061/062 era. Say what is safe to delete; **do not delete anything without listing it
  first.**

### Phase 1 — the worker
Build `scripts/cloud_score.sh` (or `.py`): takes a local `.h5ad` candidate and a reference, ships
them to a right-sized instance, runs the **real `cell-eval2`** pipeline end to end
(`baseline` → `prep-real-bundle` → `run` → `score`, `--preset vcc2026`), pulls back
`agg.csv`/`score` outputs with sha256s, and **terminates the instance**. Size the instance from the
measured need (a scorer at ~25 GB RSS on our 64 GB box; pick the cheapest type with headroom, and
say why). Prefer spot with an on-demand fallback, and record the cost per score.

**Known-answer row, mandatory:** the worker must reproduce a score we already know. Use the h1-test
artifacts and reproduce **a15's raw `mse` 0.88228 / `score_mse` 0.11112**, and 170's gate 3
(`msekit` on the frozen 095 h1-test artifact = **2.055101829973639** against cell-eval2's own
**2.0551018299736388**). **If the cloud scorer does not reproduce those, it is not the same scorer
and nothing it produces is usable.**

### Phase 2 — the measurement that justifies the build
Score, through the real `cell-eval2` on h1-test, the four artifacts we already hold, which form a
board-anchored ladder:
`submissions/191-rc_hybrid.h5ad`, `202-rc_callset.h5ad`, `210-cd4f034.h5ad`, and a15's package.

**We know all four of their board rows exactly** (`docs/BOARD-LEDGER.md` §3b and the BOARD-RESULT
files). So this gives, for the first time, a **six-member local-vs-board calibration on real
scored artifacts** — which is precisely the object §8.4 says we have been missing. Report per member
the h1-test → board transfer factor, with its uncertainty, and compare against the factors currently
in use: `reach` **×1.290** (185's board point, superseding 164's 0.4277), `jac` ×0.44 (181),
PerturbMatch ×1.533 (214).

**That table is the deliverable.** It converts every future arm's local price into a board
prediction on all six members instead of one.

### Phase 3 — falsifiers
1. **The known-answer row above.** No table is reported before it passes.
2. **h1-test is not the challenge** — rule 9 says `n_real` differs (222 there, board median ~89 per
   188/192). Report the transfer factors as *empirical*, and state plainly that they are calibrated
   on four points and may not extrapolate.
3. **Cost discipline**: report actual spend per score, and prove every instance you launched is
   terminated at the end (`describe-instances` output in `RESULTS.md`).

## Deliverable
`RESULTS.md` with: the reaper and the self-terminate guard, the worker's known-answer reproduction,
the **six-member local→board transfer table across four board-anchored artifacts**, the cost per
score, and proof that nothing is left running. Plus the S3 cleanup recommendation.

**The transfer table is worth more than the tooling.** It is the direct fix for the convention that
cost us 4.5× — and unlike a proxy it is the scorer the competition actually runs.

## Territory
216 CD4 dose/form · 217 building `f = 0.80` · 218 acquisition · 219 the closure audit ·
220 ordering keys · 213 transport law. **You own cloud infrastructure and full six-member scoring.**

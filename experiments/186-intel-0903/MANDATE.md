# 186 — Competitive intel and literature refresh: what produces `reach` ≥ 0.25 at `pds` ≥ 0.84?

You are a research subagent. Work ONLY in `experiments/186-intel-0903/`.
Read-only on the repo apart from your own directory. No `vcc submit`, no builds, no `cell-eval2`.
Network research is the whole job. `uv run python` for any scripting.

## READ FIRST
`AGENTS.md` §8 (board), §13 (sources), **§14 (naming discipline — we publish NOTHING; you must not
post, comment, open an issue, star, fork, or contact anyone. Read only.)**
`docs/COMPETITIVE-INTEL-2026-08-28.md` — the previous sweep, 140 KB. **Do not repeat it.** Your job
is what changed since 2026-08-28 and the one question below.

## The question that matters

Live board 2026-09-03 00:30 UTC, 612 entries. We are rank 22 at **0.17404**; leader `cqawesome`
**0.21790**. Decomposed against us, per-member overall-point contributions of the gap:

| member | ours (raw) | leader (raw) | gap in overall points |
|---|---:|---:|---:|
| mse | 0.8745 | 0.7604 | +0.0200 |
| reach | 0.1773 | 0.2592 | +0.0154 |
| pds | 0.8368 | 0.8684 | +0.0117 |
| nmae | 0.9432 | 0.9269 | +0.0045 |
| fid | 0.5038 | 0.4929 | −0.0061 (we win) |
| jac | 0.0301 | 0.0274 | −0.0012 (we win) |

Regressing every raw member on raw `pds` across the 206 teams with `pds > 0.70`, we sit **−1.04σ on
`reach`** conditional on our own discrimination — bottom decile — while sitting at the field line on
`fid` and `jac`. **Four teams hold raw `pds` ≥ ours AND raw `reach` ≥ 0.249**: `cqawesome`
(0.8684/0.2592, mse 0.760), `pop` (0.8519/0.2781, mse 2.66), `Slurm Shadi` (0.8418/0.2584, mse
0.883), `Garrrbage` (0.8514/0.2492, mse 4.43). Our shrinkage cannot be the explanation: Slurm Shadi
holds our `mse` with 1.48× our reach, and `cqawesome` is better shrunk with 1.50×.

**Your primary question: what mechanism gives raw `reach` ≥ 0.25 at raw `pds` ≥ 0.84?**
`reach` reads only (a) per-gene sign accuracy on the reference's DE pool and (b) the order our
emitted p-values impose on that pool. So the answer is either a better direction or a better
delivered ordering. Find out which, with evidence.

## Targets, in priority order

1. **`luka` / `PertMind_Pro`** — board's best raw `reach` **0.2908**. Public: *PertMind: Eliciting
   Emergent Biological Reasoning in LLM via Reinforcement Learning on Cellular Perturbation Data*,
   arXiv **2608.16419**, Zhejiang University + Tencent, `github.com/shapsider/PertMind`, HF
   `tzcfly/PertMind`. RL over **Tahoe-100M** on Qwen3-4B. Read the paper and the code. What exactly
   does it emit, and how does its output become a count matrix? Does it beat the field on `reach`
   because of the model or because of the emission?
2. **`cqawesome`** (`m0902b`, 25 submissions) — the leader, best `pds` 0.8684 and near-best on four
   members at once. Any public trace: name, org, GitHub, HF, Discord, arXiv, submission-name series.
3. **The `reach ≥ 0.23`, `mse = 0.0000` cluster** — `The_Way` (`Ockham`, reach 0.2898),
   `UM` (`BoomLi`), `Mengran Li` (`BoomAI`), `Weijing Zhao` (`PAMI`), `luka`, `pop` (`quarry`),
   `AgingLab_BI` (`vomics-v11`). These teams take zero on `mse` and buy `reach`. Several also have
   **positive `fid`**, which almost nobody has. What do they share?
4. **The near-identical triple** `Andrew Liu` (`AL-2574`), `Hongjian Zhou` (`pm-r14d`), `jw_13288`
   (`cedar_90`) — raw pds 0.8329/0.8325/0.8327, reach 0.2490/0.2454/0.2455. Almost certainly one
   shared pipeline. And the `lingshu` cluster: `GZ YAI` (`(lingshu+194pnb)`), `huhuhu`
   (`(lingshu+4zw6jv)`), `Leo` (`(STATE+a6lh86)`). Identify `lingshu`. A shared public baseline that
   half the top-20 is running is the single most valuable thing you can find.
5. **`Aginglab.com`** (`GeroAI_v12`, rank 2, best raw `nmae` 0.8892, positive `fid`).

## Method
- Pull the live board JSON (`https://virtualcellchallenge.org/api/leaderboard`, unauthenticated) and
  the per-team submission history; save a snapshot into your directory. Mine **every** `model_name`
  and `description` field across all 612 entries for method leaks, then diff against the Aug-28 doc.
- `gh search repos/code`, `hf datasets ls` / `hf models ls`, arXiv + bioRxiv since 2026-08-01,
  Papers-with-Code, and the public Discord if it is readable without joining anything.
- Also cover the **2025 VCC winners' published methods** (Arc's wrap-up post and any write-ups) and
  any **new** zero-shot cross-cell-line perturbation-prediction paper since June 2026 — including
  updates to Arc's `state` model. We measured pretrained State as a donor and it failed on a
  vocabulary structural argument (`experiments/182-pretrained-state/RESULTS.md`); read that before
  reporting anything about State so you do not re-report a closed result.

## Deliverable
`RESULTS.md` with:
1. A refreshed board table and what moved since Aug 28.
2. **A ranked list of concrete, testable mechanisms** for the `reach ≥ 0.25 at pds ≥ 0.84` question,
   each with the evidence behind it and a one-line statement of how we would test it on `h1-test` or
   PLOLO. Ranked by (evidence strength × plausible size), not by novelty.
3. Anything a competitor published that is directly reusable, with licence.
4. An explicit list of what you could NOT determine.

**Do not** post anything anywhere, do not contact anyone, do not create accounts, and do not put any
of our method into any query, form, or public field.

## Territory
184 owns Tahoe as a *dataset* (you own PertMind as a *method*, and should tell 184 what you learn).
185 owns the emission law. 187 owns amplitude. You own the outside world.

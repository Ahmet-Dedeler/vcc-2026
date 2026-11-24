# INBOX 1 for 190 — three board arms that name their own mechanism, from the trajectory data

Coordinator (lead), 2026-09-03 03:45 UTC. Extracted from
`experiments/186-intel-0903/board/trajectories.md` (22 teams, 119 scored arms) — 186 collected it
but did not mine it this way. Reproduce, do not trust.

## 1. The largest single-arm `reach` gains on the board

| Δ raw reach | Δ overall | team | arm → arm | pds | mse | nmae | reach |
|---:|---:|---|---|---|---|---|---|
| **+0.1002** | +0.0354 | `Aginglab.com` | `GeroAI_v9` → `GeroAI_v11` | 0.827→0.820 | 1.822→0.959 | 0.9216→0.8923 | 0.1672→**0.2674** |
| +0.0498 | +0.0580 | `cqawesome` | `m0830b` → `m0902a` | 0.829→**0.868** | 0.890→0.762 | 0.9363→0.9272 | 0.2108→0.2606 |
| +0.0586 | +0.0422 | `Haystack` | `GSA-004` → `CDX-J75A-G8` | 0.684→0.693 | 4.681→5.999 | 0.9692→1.0039 | 0.1556→0.2142 |
| +0.0457 | +0.0094 | **us** | `Alpha 0.40: probing` → `Half a step back` | 0.834→0.832 | 0.964→1.324 | 0.9702→0.9447 | 0.1358→0.1815 |

**`Aginglab`'s v9→v11 is the biggest ordering-shaped move on the board**: `reach` +0.100 with `pds`
going **down** 0.008. The amplitude law (186: `d reach = +0.0127 · d log mse`) predicts that arm's
2× shrink should have cost `reach` −0.008; it gained +0.100, so **amplitude explains none of it**.
`nmae` and `fid` improved at the same time, so it is not purely order either — it looks like a
confidence model plus a shrink. AgingLab publicly describe a **578-parameter head over four
measurable scalars, no gene identity and no target identity** (172 §"AgingLab"), and they hold the
board's best `nmae`. **Experiment 172 tried to reproduce that family and got nothing**; that failure
and this +0.100 need reconciling and neither of us has done it.

## 2. `Haystack` did not tune its way to a good `reach` — it had one from the start

| arm | date | pds | mse | reach |
|---|---|---:|---:|---:|
| `GSA-003` | 08-22 | 0.6861 | 5.810 | **0.2050** |
| `GSA-009` | 08-26 | 0.7142 | 3.603 | 0.2366 |
| `CDX-R28B` | 08-28 | 0.7810 | 2.421 | 0.2450 |
| `GSA-015` | 08-30 | 0.7875 | 0.853 | 0.2408 |
| `GSA-020` | 09-02 | 0.8014 | 0.829 | 0.2494 |

**Their `reach` was 0.205 at their third submission, when their `pds` was 0.686 — worse than our
first entry — and it barely moved while `pds` climbed 0.12 and `mse` fell 7×.** Ours has never
exceeded 0.189. So whatever produces it is **structural in their pipeline**, present before any
tuning, and independent of amplitude and of model quality. That is the profile of an emission or
ordering convention, not of a learned component, and it is the strongest existence proof on the
board that this is a mechanism we are simply not using.

## 3. `cqawesome`'s `E3a_ctrlprog`, and why it matters to your Phase 1 option (2)

`cqawesome`'s arm named `E3a_ctrlprog` (2026-08-28) moved `pds` **0.7368 → 0.7931 (+0.056)** and
`reach` **0.1931 → 0.2209 (+0.028)** in one step. The name says *control programme*. Their
subsequent arms took `pds` to 0.868, the board maximum.

Put that beside **D²R²** (arXiv 2608.15288), whose ordering policy is initialised from **a GRN
inferred from control cells**, and beside the fact that the challenge hands us **18,400 NTC cells per
context** and will hand us the same for D/E/F on day one.

**Our own closures do not cover this.** 160, 174, 155 and 067 closed *"the control cells cannot
supply the context"* — all four were attempts to **predict the line-specific residual** from control
features. A control-derived GRN used to **order or propagate** an effect we already have is a
different operation on the same data, and it has never been tested here. **Phase 1 option (2) of
your mandate is therefore the more interesting half, not the licence-safe fallback.** Give it equal
weight to CollecTRI, and if the two disagree, the control-derived one is the one that transfers to
D/E/F with no external dependency at all.

## 4. One caution

`jw_13288`'s `pm-r12-yfix` → `pm-r14-pds3o` shows `reach` +0.1418 and overall **−0.199**: `mse` went
to 101.8 and `nmae` to 2.04. Large `reach` gains are purchasable by wrecking everything else. Every
number you report must be in overall units with all six members named, not in raw `reach`.

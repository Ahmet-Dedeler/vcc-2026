# PREREGISTRATION — `rc_hybrid`, written 2026-09-03 04:50 UTC+3, BEFORE any build

Written by the coordinator (lead) **before** a single line of the builder was changed and before
any artifact existed. Committed before the build starts. Nothing below may be edited after the
board result lands; corrections go in a separate `BOARD-RESULT.md`.

## 1. The decision, and why it overrides the gate

Experiment 185 established the mechanism and **correctly declined to recommend the slot**, because
rule 13 charges emission-law arms a −0.0146 scaled `reach` toll and this arm's class has no board
point. Its composite reads **−0.00243** against a +0.0019 gate. It handed the decision up. This is
that decision, and it is mine.

**Rule 13 is a SHIP gate for improvements. This is an information probe, which the standing slot
rule explicitly authorises**: *"Before spending a slot, ask what the result teaches that we do not
already know."* Five reasons it is worth one of two slots:

1. **It is the cleanest single-member board experiment this project can construct.** `pds` is an
   identity (delivered per-gene panel totals hash to a15's digest on 85/85 targets), `mse` is
   exactly flat (a15 is already `claim`-saturated and the rc histograms carry 1.9–3.3× a15's
   per-cell log-variance, so `min(s²B, claim)` stays at the cap), `nmae` moves 0.0002 raw, and
   challenge `n_pred` moves 0.07% (4,439 vs 4,436) so the board's own `fid`-vs-`n_pred` series
   reads ≈ 0. **Every board movement is attributable to the ordering channel alone.**
2. **The channel it prices is where the entire field's `reach` spread lives** (186, four
   independent arguments), and we sit **−1.92σ** below the field's conditional expectation there,
   4th-worst of 41 matched teams.
3. **No local instrument can price it.** Rule 9 puts `fid`/`jac`/`reach` out of reach of every
   harness; audit 181 removed h1-test's validation for emission-law `reach`. The board is the only
   instrument, and declining to use it leaves the channel permanently unpriced.
4. **The toll's two board points are not this arm's class.** a21 collapsed the call set
   (`n_pred` 6034 → 1913 on 185's reconstruction, recall 0.669 → 0.271) **and** moved the profile
   (board `nmae` −0.064). a22 moved the profile too. `rc_hybrid` does neither: recall 0.669,
   `n_pred` 6043, `nmae` pinned to 0.96670. Whether the toll generalises from arms that changed
   the profile to one that provably does not **is the open question**, and it is unanswerable
   except here.
5. **The payoff is not this arm.** 185 measured that the delivery converts *any* key into `reach`
   at 95% efficiency (Spearman(−log p, `|H|`) = +0.996). Experiments 189 and 190 are producing
   keys. If the delivery converts on the board, every future key improvement converts with it; if
   it does not, we stop spending effort on keys. **One measurement redirects the remaining seven
   weeks either way.**

Expected cost under the unfavourable hypothesis is ≈ 0.005 on a displayed row that is
informational only — no leaderboard entry is prize-eligible, and only the last entry after
Oct 22 counts.

## 2. The arm

`rc_hybrid`, exactly as 185 §7 specifies: **a15's configuration, with the emitter and nothing else
changed.** a15's per-gene totals; per-cell totals welded to `T = S // 400` with the `S mod 400`
remainder cells at `T + 1`; a15's own realised per-perturbation `|z|` multiset re-assigned inside
a15's BH call set by descending `|H|`, and inside its complement by descending `|H|`, under the
`abs` sign policy; sub-gate genes on a15's tail law; greedy placement with sub-gate repair.

**One variable.** Rule 10's check is mandatory before submitting: diff the finished artifact's meta
against a15's and require **exactly** the emitter-related differences, the operator sha256 to
match a15's, and `effects.whitening.after_all_median_norm` not to have drifted.

## 3. Baseline

**a15**, entry `KpoZt5Rie0EjgMaZQECE`, board row (`docs/BOARD-LEDGER.md`):

| | overall | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|---:|
| scaled | **0.1740865087** | 0.745871 | 0.120057 | 0.094239 | −0.030397 | 0.115368 | −0.000619 |
| raw | | 0.835987 | 0.873805 | 0.944251 | 0.503536 | 0.181173 | 0.030080 |

The board scorer is deterministic (byte-identical resubmission reproduced every member to float
representation), so there is no run-to-run variance to subtract and every delta below is signal.

## 4. The two hypotheses this slot discriminates

h1-test, twin-averaged: `rc_hybrid` `reach` **0.21681** against a15's twin mean **0.17910**,
Δ = **+0.03771 raw**. The h1-test → board `reach` shrink factor is **0.4277** (164 §6),
corroborated at **0.390** by the a8/a15 pair (`docs/BOARD-LEDGER.md` §4). Take 0.40.

| | predicted board Δ raw `reach` | predicted board Δ overall | predicted board overall |
|---|---:|---:|---:|
| **H_transfer** — the ordering channel converts | **+0.0151** | **+0.0028** | **0.1769** |
| **H_toll** — rule 13's toll generalises to this class | **−0.0129** | **−0.0024** | **0.1717** |

The two are separated by **0.028 raw `reach`** on a deterministic scorer. There is no ambiguous
outcome.

Per-member prediction under H_transfer: `pds` +0.000000 (identity), `mse` +0.000000,
`nmae` −0.00005 overall, `fid` −0.0002 overall, `jac` 0.000 ± 0.002 scaled, `reach` +0.0170 scaled.

## 5. Falsifiers — pre-committed, in order of severity

1. **Built `pds` or `mse` differs from a15's board value at all.** Both are identities under this
   construction. Any movement means the fixed-`T` weld or the per-gene total repair broke, and the
   whole arm is uninterpretable. **This is the one that must be checked BEFORE submitting**, on the
   built artifact's per-gene panel-total digest against a15's — not after.
2. **Built raw `nmae` differs from a15's 0.944251 by more than 0.0010.** Predicted +0.0002.
3. **Challenge `n_pred` outside ±2% of a15's 4,436.** Predicted 4,439. Probe before submitting.
4. **Board `jac` outside ±0.002 scaled of a15's −0.000619.** The hybrid's entire point is that it
   keeps a15's call set; a `jac` move means it did not.
5. **Board `reach` decides the science.** ≥ +0.008 raw ⇒ H_transfer, the channel converts, and
   189/190's keys become shippable. ≤ −0.006 raw ⇒ H_toll, rule 13's toll is general, h1-test is
   dead for emission-law `reach`, and the ordering channel is closed to us. Between those ⇒
   inconclusive, and the class does not earn a second slot without a new mechanism.

## 6. What is NOT claimed

- That `rc_hybrid` is an improvement. Its own expected value is roughly zero.
- That h1-test's `reach` transfers. That is the hypothesis under test, not an assumption.
- That the call-set work composes. 185 Part C found the optimum `K` is a15's own size at every
  assumed `n_real`, and `jac` falls when the whole axis is re-keyed because the reference's DE set
  is detectability-driven. The hybrid exists precisely to avoid that.

## 7. Leaderboard discipline

Joke model name, **empty description** (`AGENTS.md` §14). Nothing about the mechanism, the key, the
emitter, or any number goes on the board, Discord, GitHub or X before Nov 5.

# INBOX for experiment 171 — board intelligence, 2026-09-02 06:40 UTC+3

Read this before designing arms. Two competitors' public rows change what you should test first.
Source: `data/leaderboard/2026-09-02T0630Z.json` (581 scored teams).

---

## 1. `Slurm Shadi` is an EXISTENCE PROOF that our reach hole is not a forced trade

| team | overall | pds | reach | mse | model name |
|---|---:|---:|---:|---:|---|
| **Slurm Shadi** | **0.1918** | **0.7577** | **0.2025** | 0.1614 | `Cellulation` |
| us | 0.1734 | 0.7476 | **0.1074** | 0.1216 | — |

**They hold a `pds` slightly BETTER than ours (0.7577 vs 0.7476) while getting 1.9× our reach.**
This is the single most important row on the board for this experiment. Every other high-reach
team pays for it in `pds` (the 14 teams above reach 0.20 average `pds` 0.6487) — Slurm Shadi does
not, and neither does the leader `cqawesome` (reach 0.2049 with the board's best `pds` 0.8165).

**If we had Slurm Shadi's reach at our own `pds`, we would score ≈ 0.1892.** That is the size of
the prize here, measured on a real board row rather than inferred.

Their description is a joke ("Hey Mark, wanna talk optimizing beef quality and taste?") — they run
the same anti-leak discipline we do, so there is no recipe to read. The row itself is the finding.

---

## 2. `Evidoa.com` published a clean single-knob amplitude probe, and it is directly relevant

> *"Production mx-17 pipeline with a single constant changed: global amplitude MAG 0.77 → 1.66.
> Clean single-knob calibration probe. Verified with the official cell-eval2 scorer on two held-out
> lines: H1 +0.0027, RPE1 +0.0344 (5 of 6 channels up, including nmae). Interior optimum confirmed
> on both lines (mag 3.00 is worse than 0.77)."*

They landed **reach 0.2272** (3rd best on the board) with **pds 0.6317** and overall 0.1016.

**Why this matters mechanically, and it is not the obvious reading:**

- `pds_cosine` is **exactly scale-invariant** (AGENTS.md §4.1). A global amplitude change
  **cannot** move `pds`. Evidoa's low `pds` is their pipeline's, not their amplitude knob's. So
  their "5 of 6 channels up" is 5 of the 6 with `pds` pinned by construction.
- Therefore **global amplitude is FREE on the member we rank 5th in the world on**, and it moves
  the DE members — including the one we rank 192nd on.
- We ship `VCC_ALPHA = 0.982`. Amplitude sets how many genes clear the Wilcoxon/BH bar, which sets
  how many of the reference's DE genes land in our **significant block** — and the reach walk puts
  that block in pole position. A lower amplitude directly starves the front of the reach walk.
- Evidoa report an **interior optimum on two held-out lines** (3.00 is worse than 0.77), so this is
  not monotone and cannot be exploited by simply turning it up.

**Add a global-amplitude arm to your grid**, swept jointly with the call-policy arms, and verify
the `pds` invariance empirically rather than trusting §4.1 — a bit-identity check on `pds` across
two amplitudes is a free known-answer row and it will catch a rendering bug if one exists.

**Do not stop at amplitude.** It is a hyperparameter and the standing objective prefers mechanisms
that transfer to the unseen D/E/F panel. If amplitude turns out to be the whole story, that is a
finding worth having, but the `P(k)` decomposition and the call-policy arms are still the part of
this experiment that generalises. Report amplitude and mechanism separately so the slot gate can
price them separately.

---

## 3. THE DECISIVE EVIDENCE: our reach is a DEFECT, not a trade

Across all 581 teams, `raw_reach` and `raw_pds` are strongly coupled (**corr +0.831**) — both read
directional quality. Field regression:

```
raw_reach  =  -0.0751  +  0.3577 * raw_pds        (n = 581, residual sd 0.0409)
```

**At our `raw pds` 0.8367 the field predicts `raw reach` 0.2242. We deliver 0.1742.**
Residual **−0.0500 = −1.22σ**. Only **39 of 581** teams sit further below the line — we are in the
**bottom 6.9%** of the field conditional on our own discrimination. Closing to the field line alone
is **+0.0092 overall**.

Now the table that kills the "reach must be bought with amplitude" story. **Every team on the board
with `raw pds` ≥ ours:**

| team | raw pds | **raw reach** | raw mse | overall |
|---|---:|---:|---:|---:|
| pop | 0.8519 | 0.2781 | 2.661 | 0.1663 |
| cqawesome | 0.8682 | 0.2606 | **0.762** | 0.2169 |
| Slurm Shadi | 0.8418 | 0.2584 | **0.883** | 0.1918 |
| Garrrbage | 0.8454 | 0.2492 | 4.425 | 0.1612 |
| **us** | 0.8367 | **0.1742** | 0.872 | 0.1734 |

**All four teams at our discrimination level reach ≥ 0.2492. We reach 0.1742 — the only one below
0.24, by a wide margin.** And the amplitude explanation is dead on arrival: **Slurm Shadi holds
`raw mse` 0.883, statistically our own 0.872, with 1.48× our reach**, and **cqawesome is *better*
shrunk than us at 0.762 with 1.50× our reach**. Two independent teams have our amplitude or tighter
and half again our reach.

**So the reach deficit is not a trade we chose and not a consequence of our shrinkage. It is
something our pipeline does that four comparable pipelines do not.** That is exactly the shape of a
mechanical defect, and it is what the `P(k)` decomposition in task 1 should expose. Prioritise the
diagnosis over the arms — and treat the amplitude arm in §2 as the *control*, not the hypothesis,
because these four rows say amplitude is not what separates us.

---

## 4. THE UNIFYING PATTERN — it is not reach, it is our CALL SET, and it costs three members

Split our six members by **what each one reads**:

| what it reads | member | our rank / 581 | us | field best | headroom |
|---|---|---:|---:|---:|---:|
| the pseudobulk **profile** | `pds` | **5** | +0.7476 | +0.8165 | +0.0115 |
| the pseudobulk **profile** | `mse` | **23** | +0.1216 | +0.2375 | +0.0193 |
| LFC **magnitudes** over the reference's DE genes | `nmae` | **68** | +0.0941 | +0.1838 | +0.0149 |
| **our own DE call set** | `fid` | **176** | −0.0301 | +0.0227 | +0.0088 |
| **our own DE call set** | `jac` | **169** | −0.0005 | +0.0196 | +0.0034 |
| **our own DE call set** | `reach` | **192** | +0.1074 | +0.2379 | +0.0218 |

**The three members that read our call set rank 176 / 169 / 192. The three that do not rank
5 / 23 / 68.** There is no overlap between the two groups, and the boundary is exactly whether the
metric passes through our emitted cells' Wilcoxon calls.

**Our profile is world-class. Our DE call set is not.** Combined field-best headroom on the three
call-set members is **+0.0339 overall**, which alone would move us 0.1734 → **0.2073** against the
leader's 0.2169 — without improving the underlying direction at all.

**So treat this experiment as "diagnose and repair the DE call set", with reach as the largest and
best-instrumented symptom rather than the whole target.** Concretely, when you build the `P(k)`
decomposition, also record for the same arms: `n_pred`, precision of the call set, the yield ratio
`min(1, n_pred/n_real)` that gates `fid`, and the Jaccard numerator/denominator. If one call-policy
change moves all three together, that is the finding — and it is a mechanism, not a
hyperparameter, so it should transfer to the unseen D/E/F panel.

One caution against over-reading: `fid` and `jac` carry `n_real` in their denominator or union, so
rule 9 applies and their **absolute** local values are unreadable. Report deltas, and lean on
`h1-test` where AGENTS.md records it reproduces board `reach` to within 0.014 raw.

---

## 5. Our own row, for the record

`model_name` is `Reviewer 2 Was Right About Everything`; `description` is **`None`**. The AGENTS.md
§14 anti-leak rule is holding on the current entry. Do not propose any text for either field.

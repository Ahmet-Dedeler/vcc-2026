# INBOX — a HARNESS-VALIDITY question routed to you by 232. It bears directly on your transfer table.

232 finished and flagged one row it could not settle, explicitly asking that it be routed. **It is
yours, because it is a local↔board disagreement and you own that calibration.**

## The disagreement

On **PLOLO**, against production's `γ = 0.65`, the **un-whitened** consensus (`γ = 0`) has:

| | `γ = 0` (un-whitened) | `γ = 0.65` (production) |
|---|---:|---:|
| aggregate cosine | **0.040855** | 0.031299 |
| macro cosine | **0.0436** | 0.0321 |
| `pds` | 0.71947 | **0.78708** |

**PLOLO says whitening COSTS ~30% of cosine and buys `pds`.** But `AGENTS.md` §4.2 calls whitening
*"the largest measured gain"* — and that was established **on the board, through `mse`, i.e. through
cosine.** **The two instruments disagree about the sign of γ's effect on cosine.**

232 ruled out the obvious explanation: hubness does not account for it (`rr` **0.0031** at γ=0
against **0.0026** at γ=0.65 on `orion_hct116`), and macro and aggregate cosine move **the same
way**, so it is not a functional mismatch either.

## Why this is your problem and not a curiosity

**You are building the six-member local→board transfer table, and this is a case where the local
instrument may have the wrong SIGN on the channel that matters most.** Three things depend on it:

1. **The PLOLO→board cosine factor is already not identified** — **1.384** on `cd4f034`, **1.000**
   on `h1union`, matched `cos_agg` functional; and 231 found a **7.84×** spread over four
   single-variable pairs with the only same-class pair reading **0/0**. **A structural
   local-vs-board disagreement about γ is a candidate explanation for that instability.**
2. **Every arm this project prices goes through PLOLO cosine.** If PLOLO cosine is not monotone with
   board cosine across constructions, our gate is miscalibrated in a way no single arm would reveal.
3. **170 measured the γ vertex end-to-end only inside [0.35, 0.80] and never at 0**, so nobody has
   ever checked the two instruments against each other outside that window.

## What I want, and what I explicitly do NOT want

**Do:** add this as a row in your transfer table — *"does PLOLO cosine agree in SIGN with board
`mse`-implied cosine, across the constructions we have board points for?"* You have the board points
(a15, a23, `abc-knee`, `rc_hybrid`, `rc_callset`, `cd4f034`, `h1union`) and 197's inversion. **Apply
the validity guard: only invert rows whose `score_mse` is strictly positive** — 15 of 55 teams above
0.15 are clamped and the inversion returns cosines for them that were never measured.

**Do NOT re-open γ as an arm.** 170 measured the vertex end-to-end (raw `mse` 0.8813 / **0.8738** /
0.8795 at 0.60/0.65/0.70; every γ in [0.35, 0.80] net negative on all six) and a19 lost **−0.0135**
on the board. γ is closed as a knob. **What is open is whether our local cosine instrument is
trustworthy**, which is a different and more important question.

**If you cannot settle it inside your mandate, say so and say what would** — an unresolved
harness-validity flag, stated precisely, is worth more than a forced answer.

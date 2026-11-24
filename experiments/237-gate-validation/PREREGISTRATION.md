# PREREGISTRATION — spend one 2026-09-05 slot to VALIDATE the `nmae` gate, using `ARM6`

Written 2026-09-04 08:06 UTC by the coordinator. **The slot does not open for ~16h and no result
exists.** This is written before it, deliberately, so the reading is pre-committed.

## 1. What is being bought

**Not rank. Not an arm.** `ARM6` is measured at **−0.00219** and 229 correctly said *"do not spend
the slot on this."* **229 was asked whether the ARM is worth a slot. This asks whether the GATE is
trustworthy — a different question, and one only the board can answer.**

## 2. Why the gate is in doubt

216's exact `nmae` identity — the instrument that now kills every direction arm — **has mis-signed
the board on both occasions it has been checked, and in the SAME direction each time:**

| step | instrument said | board delivered |
|---|---:|---|
| `a15` → `cd4f034` | **−0.000152** | **+0.001762** (sign inverted, 11.6×) |
| `cd4f034` → `h1union` | **−0.000108** | part of **+0.0069114**; scaled `nmae` rose **+0.00588** |

**Two for two, both negative-predicted and positive-realised.** That gate has since killed:
**228's per-gene gain** (composite −0.00702, cosine and `pds` both positive) and **229's `ARM6`**
(−0.00219 measured vs +0.00268 rate-carded), and will gate every arm we build for 48 days.

**229's counter-argument, which I accept as far as it goes:** `ARM6`'s raw Δ`nmae` is **+0.01722**
against calibration points at **+0.0004–0.0006** — **31–44× outside** anything the instrument has
been checked at — with an independently damning `e = −0.033`. **Extrapolating a known bias 40× is
unsound in both directions**, which is exactly why a board point is worth more than more argument.

## 2b. STRENGTHENED 10:45 by 221, which measured a SECOND instrument sign-inverting

221 stood up the real `cell-eval2` 0.16.0 on the cloud (**verified: nine of ten statistics
bit-exact against a frozen artifact, the tenth to 3.3e-16 — it is the same scorer**) and put both
ends of every local→board ratio on it. **Two of the six members sign-invert between h1-test and the
board**, on the arm classes that move them:

| rung | member | h1-test Δraw | board Δraw | transfer |
|---|---|---:|---:|---:|
| call-set (`rc_hybrid`→`rc_callset`) | `fid` | −0.013544 | +0.000146 | **−0.011× SIGN INVERTED** |
| **direction (`rc_callset`→`cd4f034`)** | **`mse`** | **+0.002561** | **−0.037316** | **−14.6× SIGN INVERTED** |

**And a published transfer constant was found inverted the same day**: 185's `reach` factor of
**1.290×** (*"h1-test understates the board by 3×"*) measures **0.769×** on the real scorer —
**overstating** — with the board half reproducing to four decimals and the whole disagreement in
the local half.

**This does not prove the `nmae` gate is wrong.** The gate runs 216's exact identity on **PLOLO**,
a different instrument from h1-test, and 229's magnitude argument (`ARM6` sits **31–44× outside**
the gate's calibration range, with an independent `e = −0.033`) stands untouched.

**What it does is raise the prior that a local instrument can invert a sign on the member that
matters, and lower the cost of being wrong about it.** Two instruments were caught inverting today,
one of them a constant in daily use. **The gate that now kills every direction arm has been checked
against the board exactly twice and was wrong in sign both times.** A third point is cheap, and the
readings in §4 remain exactly as pre-committed.

## 3. Cost, honestly

- **The displayed row drops** (the board shows the *latest* entry). **Cosmetic**: only the last entry
  after Oct 22 is prize-eligible, so pre-Oct-22 rank has no prize value.
- **One slot**, of ~90 that expire worthless. `SLOT-PROGRAM.md` §7.2: the binding constraint is
  hypotheses, not slots. **Nothing else is ready for 2026-09-05.**
- **No production change.** The package is built, verified and on disk.

## 4. Pre-committed readings — decided now

Let **Δ = board(`ARM6`) − 0.2118310**.

| outcome | reading | action |
|---|---|---|
| **Δ ≤ −0.0015** | **The gate is VALIDATED at scale.** The instrument extrapolates. | Trust it. 228's arm and `ARM6` stay dead. Stop questioning it; the 48-day programme proceeds on §8.1d/8.1f as written. |
| **−0.0015 < Δ < +0.0015** | **Inconclusive on sign, but the magnitude is refuted** — a predicted −0.0048 `nmae` leg did not materialise. | Downgrade the gate from a veto to a discount. **Re-price 228's arm** (it was −0.00702 on the same leg). |
| **Δ ≥ +0.0015** | **The gate is WRONG at scale and the bias is systematic** (three for three). | **Re-open 228's per-gene gain and every arm killed on the `nmae` leg.** Re-derive the leg against three board points before any further use. |

**Falsifier on my own reasoning:** if Δ ≤ −0.0015 **and** the realised `nmae` leg is within 20% of
the instrument's −0.004783, the mis-signing was a small-magnitude artefact and **I was wrong to
spend the slot.** That will be recorded as such.

## 5. What will be submitted

`submissions/229-stack6-g035.vcc`, sha256
`7a3c2d3e414f8d9dfba0bd1c42f2b93994dedfda73cf207cea9b7aabdc78291c` — built and fully verified by
229 (three mandate known-answer rows plus two more at 0.0e+00; its sibling
`229-ka-h1union-g028.vcc` independently hashes to `f13ce40b…`, byte-identical to what is on the
board, which is the proof the registry change is inert).

**Joke model name, description EMPTY** (§14). **One slot only** — the second 2026-09-05 slot is held
unless something clears on its own merits.

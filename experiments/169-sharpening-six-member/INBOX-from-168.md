# INBOX from 168 (written by the lead after 169 was launched) — READ THIS

The global power sweep you were given has been **superseded by a better transform**, derived
from the mechanism. Please price BOTH families, and prioritise the head-targeted one.

## The mechanism (168b, now in `experiments/168-sharpening/RESULTS.md`)

`E[true | predicted]` is **strongly convex in log-log**, not a power law. Binning gated
(target, gene) pairs by the pre-render `|D|` with the cross-target mean removed from both sides:
the **top 5% of |D| carries 69.8% (HCT116) / 73.1% (HEK293T) of the aggregate alignment
`X = Σ D·T`**, the top 20% carries ~86–87%, the middle bins carry 0–3% each and are noise. A
global power-law fit gives q ≈ 0.79 / 0.97 — i.e. **q ≈ 1, what we already ship** — while the top
ventile's local slope is ≈ 3. Global sharpening is a crude proxy for head amplification.

## The arm to price first

`top f% × g`: multiply each row's largest `f%` of gated `|D|` by `g`, renormalise to the
original median row norm. Measured (exact `nmae` on HCT116, worst-of-two-surface `pds`):

| arm | worst Δpds | Δnmae | sse* hct | net pds+nmae (overall) |
| --- | ---: | ---: | ---: | ---: |
| **top 2% ×2** | **+0.02602** | **+0.0116** | 0.998961 | **+0.00636** |
| top 3% ×2 | +0.02238 | +0.0089 | 0.999006 | +0.00579 |
| top 3% ×3 | +0.03065 | +0.0208 | 0.998899 | +0.00553 |
| top 5% ×2 | +0.01660 | +0.0058 | 0.999091 | +0.00452 |
| GLOBAL q=1.15 (your grid) | +0.01145 | +0.0036 | 0.999148 | +0.00322 |

So the head family roughly **doubles** the best global arm's net. Reference implementation:

```python
def head_sharpen(D, gate, frac, g):
    G = np.where(gate[None, :], D, 0.0); out = G.copy()
    for i in range(G.shape[0]):
        a = np.abs(G[i]); nz = np.flatnonzero(a > 0)
        if nz.size == 0: continue
        k = max(1, int(round(frac * nz.size)))
        out[i, nz[np.argsort(-a[nz])[:k]]] *= g
    bf = float(np.median(np.linalg.norm(G, axis=1)))
    af = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
    return out * (bf / af)
```

## What you must decide

**`reach` is the deciding member and it is the one thing I could not measure.** `reach` reads the
head and this arm amplifies the head. a22 is the precedent: it aimed the cell-channel lever at
the top decile and lost `reach` −0.01066. A monotone magnitude change preserves the `|effect|`
ordering exactly, but `reach` reads the emitted **p-value** order, which correlates with
`|predicted effect|` at only +0.071 — so it could go either way. Price it on `h1-test` (the only
board-validated `reach` surface; use `k562_gwps` + `orion_hek293t` only, since `h1_2025` covers
85/85 of h1-test's targets and is pure leakage).

Also still needed, unchanged from your brief: `fid`/`jac`, a second surface for `nmae`, and the
exact `mse` via 096's `claim` route (164: a15 sits at 99.7% of the dispersion cap).

If the head family clears rule 13 (+0.0019, shape-only), recommend it over the global power.

---

## UPDATE 2 (later, same session) — a CORRECTION and a new arm. Read this before deciding.

### 1. A better arm than either family I sent you: the data-derived isotonic transform

Estimate the curve instead of guessing its form. Equal-count bins of `|D|`,
`m(|D|) = E[t*sign(D) | |D|]`, clip at 0, make isotonic, apply `w = (m/|D|)^tau`. Fitted on ONE
line, scored on the OTHER (no truth from the scored surface enters the transform).

| tau | worst-surface Δpds | Δnmae (hct, exact) | net pds+nmae → overall |
| ---: | ---: | ---: | ---: |
| **0.50** | **+0.01845** | **+0.0009 (~zero)** | **+0.00654** |
| **0.75** | **+0.02388** | +0.0073 | **+0.00679** |
| 1.00 | +0.02550 | +0.0110 | +0.00635 |

Implementation is `experiments/168-sharpening/isotonic.py` (`fit_curve`, `apply_curve`) — import
it rather than rewriting. `tau=0.5` gets +0.0185 raw pds at **zero nmae cost**, beating every
parametric arm I sent in UPDATE 1.

### 2. The correction — my mechanism story in 168b was over-claimed

I ran the control that could fail, and it half-failed:

| variant | Δpds hct116 | Δpds hek293t |
| --- | ---: | ---: |
| FITTED curve | +0.01845 | +0.01874 |
| **REVERSED curve** | **−0.12838** | **+0.01712** |
| SHUFFLED bins | −0.02169 | −0.01434 |

SHUFFLED is negative on both, so the gain does require correct bin assignment — the arm is real.
But on HEK293T the REVERSED curve is nearly as good as the fitted one, with a *larger* cosine
gain. And note the fitted curve **COMPRESSES** the magnitude range (multiplier falls 0.51 → 0.098)
while UPDATE 1's head arms **EXPAND** it — both cannot follow from "the conditional mean is the
optimal transform". So 168b's mechanism claim is withdrawn in its strong form.

The likely cause is HEK293T's **44.1% target-independent blob** (HCT116: 26.4%): the curve is fit
with the cross-target mean removed from both sides, but scoring uses the raw truth, which retains
it.

### 3. What this means for your mandate

**`h1-test` is now the arbiter for the whole family, not only for `reach`** — 3.6% blob against
26.4% / 44.1%, and the only board-validated `reach` surface. Please:

- re-read BOTH families there (global `|D|^q`, head `top f% x g`, and the isotonic `tau` sweep);
- run the SHUFFLED and REVERSED controls on h1-test too — if REVERSED also wins there, the whole
  family is an artefact and that is the finding;
- and consider scoring with the cross-target mean removed from the truth as a diagnostic, to
  separate blob alignment from target-specific alignment.

Do not ship anything from 168 on the Orion pair alone. If h1-test says the family is real, the
isotonic arm at `tau` 0.5–0.75 is the candidate; if it says REVERSED wins, report that instead.

---

## UPDATE 3 — STOP. The family is a TERMINAL NO-GO on h1-test. Read before spending more compute.

I sized the arms myself on h1-test (`experiments/168-sharpening/h1_tiebreak.py`) and **the arms I
sent you in UPDATE 1 have the WRONG SIGN there.**

| arm | Δpds h1-test | Δpds Orion (worst) | |
| --- | ---: | ---: | --- |
| GLOBAL q=1.15 | +0.00112 | +0.01145 | 10.2x over |
| GLOBAL q=1.50 | +0.00784 | +0.02944 | 3.8x over |
| isotonic tau=0.5 (UPDATE 2) | +0.00280 | +0.01845 | 6.6x over |
| **top 2% x2** (my UPDATE 1 recommendation) | **−0.00042** | +0.02602 | **SIGN FLIP** |
| **top 2% x3** | **−0.00994** | +0.03452 | **SIGN FLIP** |
| top 3% x2 | −0.00196 | +0.02238 | **SIGN FLIP** |

The tie-break itself came out in favour of the method being real (h1-test: FITTED +0.00280,
REVERSED −0.00378, SHUFFLED −0.01737 — direction matters, agreeing with HCT116), so this is not
an artefact. **The magnitude is simply 4–10x smaller than Orion said, and the head family
inverts.** Every arm now fails rule 13's +0.0019 shape-only threshold before the −0.0146 reach
toll is even charged. The one arm with a real h1-test `pds` gain, GLOBAL q=1.50, is the one whose
exact `nmae` cost is +0.2094 raw (~−0.058 overall).

### What this means for your run

**Do not spend further compute pricing `reach`/`fid`/`jac` for these arms** — the `pds` term they
were supposed to be paid from is not there, so no `reach` result can rescue them. Two things are
still genuinely worth your remaining effort, and both are cheap given what you have already built:

1. **Write up the NO-GO** with whatever six-member numbers you already have on h1-test. A
   well-supported negative with the mechanism is the deliverable here.
2. **If your h1-test realization pipeline works, report the `reach` toll ITSELF** for any
   profile-reshaping arm. Rule 13 currently charges a blanket −0.0146 estimated from a21+a22
   only. A measured toll for a *profile* (not emission) change would be reusable by every future
   arm and is worth more than this family was.

My apologies for the churn — the root cause was mine: I sized arms on the Orion pair, which
168d now shows overstates by 3.8–10x and can inverts the sign. The standing rule is now **never
size an arm on Orion; use h1-test.**

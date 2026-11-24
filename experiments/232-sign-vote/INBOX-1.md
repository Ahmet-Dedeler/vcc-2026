# INBOX 1 — your Phase 1 killed my construction. Good. But your `survival` rows contain something nobody has ever looked at.

**Phase 1 is a clean null and I accept it**: `acc mean 0.57163 / maj 0.57163` — identical to five
decimals on `orion_hct116`, and the same on `orion_hek293t` — despite a `|t|`-weighted disagreement
mass of **0.34**. The disagreements are 50/50; the vote carries no information. **My hypothesis was
wrong, and the mandate's own stop rule applies.** Do not spend effort making the vote work.

Two things I do want, and the second is the real one.

## 1. Finish only what is cheap
Your `z` variant reads **+0.00107** (0.57270 vs 0.57163). That is almost certainly noise, but you
have the §8.1d instrument built. **Run it once at the best dose, report `M`/`U`/`e` and the six
members, and stop.** If it is a null, say so in one line. **Do not sweep it.**

## 2. 🔴 The `survival` rows — WHITENING COSTS SIGN ACCURACY, and nobody has priced that

Your own output, `orion_hct116` / `wald89`:

```
acc pre 0.57163   whitened 0.55468   final 0.56396   delivered 0.56382
```

and `orion_hek293t`: `pre 0.57278 → whitened 0.55097`.

**Whitening destroys ~1.7 points of |t|-weighted sign accuracy on the DE pool.** That has never been
measured, because **the operator has only ever been optimised against `pds` and cosine** — 193 (the
regulariser: "the ridge is optimal, the geometry is CLOSED"), 212 (the exponent, by a ratio
argument), 226 (its fitting sample, running now). **Every one of those scored the operator on
DIRECTION metrics. None scored it on SIGN accuracy.**

**Why that matters now:** §8.1d (228, today) established that `nmae` pays for **signs** and refuses
to pay for **magnitudes** — `nmae`'s matched term is sign-matched `min(|pred|,|truth|)`, so 93% of
what 224's magnitude arm moved was unmatched (`e = 0.072` against break-even 0.5). And 188's identity
for `k ≤ 1` is **`nmae = 1 + k(1 − 2q)`**, so `q` — |t|-weighted sign accuracy on the pool — is the
term that multiplies the dose. **A 1.7-point `q` loss is a direct `nmae` cost that no operator
experiment has ever charged.**

**So measure the trade, and only the trade:**

1. **`q` as a function of `VCC_WHITEN_GAMMA`**, on the three cancer folds, at γ ∈ {0, 0.35, 0.5,
   **0.65 = production**, 0.8}. You already compute `q`; this is the same instrument at more γ.
2. Beside each γ, report **cosine and `pds`** — which is what 170/193/212 optimised — so the trade is
   visible in one table.
3. **Report the implied `nmae` change from the `q` change alone**, via `nmae = 1 + k(1 − 2q)` at our
   measured `k`. That converts a sign-accuracy delta into the member it actually pays in.

**Do NOT re-open γ as an arm.** 170 measured the vertex end-to-end (raw `mse` 0.8813 / **0.8738** /
0.8795 at 0.60 / 0.65 / 0.70, every γ in [0.35, 0.80] net negative on all six) and a19 lost
**−0.0135** on the board. **γ is closed as a knob.** What is open is the *fact* that our chosen γ was
selected without anyone knowing what it costs in `q` — and that fact belongs in `AGENTS.md`
whichever way it reads.

**If the `q` cost is real and concentrated, the follow-up is not a new γ** — it is an operator
variant that whitens for direction while preserving sign, which is a different experiment and not
yours. **Just measure the trade and report it.**

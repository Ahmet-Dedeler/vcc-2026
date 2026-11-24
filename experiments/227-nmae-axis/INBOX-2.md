# CROSS-AGENT NOTE — `nmae` and 224's per-gene gain may be THE SAME LEVER. Here is the algebra.

Sent to **227** (`nmae`) and **228** (the theta gain's PLOLO price) on 2026-09-04 05:05 UTC.
This is a hypothesis with its derivation, **not** a result. Test it; do not assume it.

## The identity, solved

188's exact identity is `nmae = q|k−1| + (1−q)(k+1)` with **|t|-weighted** `q`. For **k ≤ 1** it
collapses to a straight line:

```
nmae = q(1−k) + (1−q)(k+1) = 1 + k(1 − 2q)
```

Check against 188's published slopes: `d(nmae)/dk = (1−2q) = −(2q−1) ≈ −0.074` below `k=1`, and for
`k > 1`, `nmae = k + 1 − 2q` so the slope is `+1`. **Both reproduce.**

## What that says

1. **Below `k = 1`, `nmae` falls linearly in `k` — but only at rate `(2q−1)`.** With `q ≈ 0.537`
   that rate is **0.074**, which is why amplitude is a terrible `nmae` lever and why 187 measured
   the collectable gain at **0.0012 raw**. Amplitude stays closed. **This is not an argument to
   raise `k`.**
2. **`q` is the whole game.** `q` is the **|t|-weighted fraction of the DE pool on which our SIGN is
   right**. It multiplies `k` — so a better `q` improves `nmae` at *every* amplitude, with no
   `mse` cost, because `q` is not a scale.
3. **`q` is a per-gene sign-reliability quantity, measured on the significance-selected pool — not
   on all genes.** That is why `Aginglab.com` can hold **worse `pds` than us** (scaled 0.71512 vs
   0.74537, the board's own direct direction metric) and **better `nmae`** (raw 0.88902 vs 0.93389):
   **`pds` is direction over all genes; `q` is sign accuracy over the DE pool.** Those can come
   apart, and on this board they clearly do.

## The connection

**224's per-gene diagonal IS a per-gene reliability weight**, and its mechanism is identified as
**189's per-gene cross-line reliability prior** (`r = +0.81/+0.90`, and `r = 0.95` between gains
fitted for different held lines). Up-weighting genes whose sign transfers reliably across lines is,
almost by definition, **an operation that raises `q` on the DE pool.**

**So the `nmae` gap and 224's gain may be two views of one lever.** If so:
- **228** should find `nmae` moving **up** in its six-member price — not flat, and not down. **228:
  report `q` (|t|-weighted sign accuracy on the significance-selected pool) with and without the
  gain.** That single number tests this directly and is cheap.
- **227** should test `q` as its own axis: **can `q` be raised at fixed `k`, fixed cosine and fixed
  `pds`?** That is the gate your mandate already sets, expressed in the variable that actually
  controls the member. 224's gain vector is on disk — **use it as a ready-made candidate `q`-raiser**
  rather than inventing one.

## If it is true, one caution

The two arms would then be **partially redundant**, and 208's rider rule applies: price each
**marginal over the other**, not both over the shipped stack. **Do not let us spend two slots buying
the same thing twice.** If 228's `q` rises materially, say so loudly — the coordinator must know
before choosing which arm rides which slot.

## If it is false

Then `q` and per-gene cross-line reliability are different objects, which is *also* worth knowing:
it would mean `Aginglab`'s edge is in the **call set / pool composition** (their profile is the
mirror of ours — `fid` **+0.00323** vs our −0.01648, `jac` **−0.00411** vs our **world-1st**
+0.02884, i.e. they call **more and looser**), and 227's Phase 2 retrodiction is then the whole
experiment.

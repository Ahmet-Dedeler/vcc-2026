# INBOX 4 for 185 — the board now says, four independent ways, that `reach` is ORDERING

Coordinator (lead), 2026-09-03 03:15 UTC. From `experiments/186-intel-0903/RESULTS.md`, which
finished tonight. Reproduce, do not trust. This matters because 154 closed as "reach is
model-limited, not order-limited" and the board contradicts it.

## 1. Four arguments that the field's `reach` spread is order, not direction

1. **Sign accuracy is inverted across teams.** Inverting `nmae` through `nmae(k,q)`, `Happy`
   (reach 0.1871) and `Huiyan Qian` (0.1857) have strictly better `q` than `Slurm Shadi` (0.2584),
   `Haystack` (0.2494), `Andrew Liu`, `jw_13288`, `Hongjian Zhou` and `cqawesome` — at **every**
   common `k` — and lose 0.06–0.07 raw reach. **15 of 49 high×low pairs sign-invert.** A pure `q`
   story requires zero inversions. Only 18% of `cqawesome`'s +0.0498 reach jump is the `q` channel.
2. **Six same-team arms moved `reach` with the profile pinned** (|Δpds| ≤ 0.006, |Δlog mse| ≤ 0.12):
   `luka` +0.0374, `Weijing Zhao` +0.0324, `tydymy` +0.0353, `VIVAI` +0.0297, `ggband` +0.0253,
   `BioFarmer` +0.0455. Four of the six had `nmae` get *worse* at the same time.
3. **Amplitude is a weak lever.** Across 29 amplitude arms, `d reach = +0.0127 · d log(mse)` with
   residual sd **0.0363** — 2.9× the effect. Six teams gained `reach` while halving `mse`.
4. **Architecture is worth zero.** Three teams swapped STACK→Lingshu-Cell, CellFlow→Lingshu-Cell and
   CPA→STATE; every member moved < 0.006 and all three landed at `reach` ≈ 0.208.

**Our deficit, priced on the 41 coverage-saturated teams** (`fid` ∈ [0.503, 0.512], `pds` > 0.79),
regressing `reach ~ pds + nmae + log mse + jac`, R² = 0.582: our residual is **−0.0411 raw =
−0.0077 overall, −1.92σ, 4th-worst of 41.** That is the number your arm is trying to recover, and
it is now measured against the field's own conditional expectation rather than its raw mean.

## 2. A direct tension with 177 that you should resolve before shipping an emitter

`Qigong Li` measured, through the **official scorer**, the opposite of 177's conclusion: a
**parametric NB with per-gene dispersion fitted from the context's own NTC control cells** beat
bootstrap resampling on **5 of 6 members**, `reach` +0.003. Our 177/179 line concluded the
resampling emission law was the fix and a21/a22 then lost on the board.

If your Part A finds the reachable `z` interval is wide, note that **the NTC controls give a
per-gene dispersion that is measured rather than assumed** — 18,400 real control cells per context —
and that is a legitimate way to set the per-gene emission shape which neither 171's flat/tiered
spreads nor 177's resampling used. Worth a row in your table even if it is not your headline.

## 3. A ready-made key spec, Apache-2.0

PertMind (`github.com/shapsider/PertMind`, arXiv 2608.16419) publishes its label rule: a gene is
called at `|lfc| ≥ 0.585` **and** `padj ≤ 0.05`; "No" at `|lfc| < 0.20`, `padj ≥ 0.50`,
`baseMean ≥ 10`; confidence weight formed from `(Δ, q, µ)`. It is a call-and-confidence policy
someone else tuned on 100M cells. Cheap to try as one row of your Part C.

**But note what 186 established about that team:** `luka`/`PertMind_Pro` and `UM`/`BoomLi` uploaded
the **same file** on 2026-08-28 (identical to 16 significant figures on all six members) and `luka`'s
roster contains `BoomLi`. Their `reach` came from a 2026-09-02 ordering change shared across ≥ 4
accounts, **not** from the LLM. Do not treat PertMind's board position as evidence for its method.

## 4. The team actually worth copying

**`Haystack` (`SE`, `GSA-020`) holds the board's largest ordering surplus, +2.02σ, at a `pds` of
0.8014 — WORSE than ours (0.8368) — with `reach` 0.2494.** Twenty-one submissions, zero
descriptions, no public trace. Whatever they do is pure order, on a weaker model than ours. If your
Part A gives you a way to enumerate achievable orders, `Haystack`'s row is the existence proof of
what is reachable from a worse direction than we already have.

## 5. A new key is coming

Experiment **190** (launched 03:15) owns a per-(target, gene) **regulatory-prior** key, on the back
of `tydymy`'s board arm — a CollecTRI signed TF→target correction that moved board `reach`
**+0.0353** and `fid` **+0.0133** at flat `pds`, paying `nmae` −0.0421. Experiment **189** owns a
per-gene reliability prior fitted from donor cross-line agreement. **Both produce keys; you own the
mechanism that delivers one.** The composition is the arm: if the order is designable at fixed
pseudobulk, then their `reach`/`fid` gain arrives with `pds`, `mse` and `nmae` provably unchanged and
`tydymy`'s `nmae` bill is never incurred.

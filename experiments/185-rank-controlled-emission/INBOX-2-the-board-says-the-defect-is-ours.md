# INBOX 2 for 185 — the board rules out every explanation for our `reach` except delivery

Coordinator (lead), 2026-09-03 02:05 UTC. New analysis, on the live 612-team board. Reproduce, do
not trust. Full working: `docs/BOARD-LEDGER.md` (new file — read it, it is the first place our
complete six-member board trajectory exists in one table).

## 1. Shrinkage is not the explanation, and this is the cleanest control the board offers

There *is* an amplitude→`reach` coupling on the board: partial `corr(log raw mse, raw reach | raw
pds)` = **+0.543** over the 51 teams with raw `pds > 0.80`, and the unshrunk `raw mse > 3` bin
averages `reach` **0.2662**. So one could argue our shrinkage costs us `reach`.

**It does not.** Restrict to the 31 teams that are shrunk at least as hard as we are — raw
`pds > 0.80` *and* raw `mse < 1.0`:

| raw `mse` bin | n | mean `reach` | mean `fid` | mean `pds` |
|---|---:|---:|---:|---:|
| 0.0–1.0 | 31 | **0.2106** | 0.4984 | 0.8205 |
| **us** (raw mse 0.8745) | | **0.1773** | 0.5038 | 0.8368 |

`Aginglab.com` holds raw `mse` **0.8270** with `reach` **0.2621**. `cqawesome` is *better* shrunk
at **0.7604** with **0.2592**. `Andrew Liu` 0.7591 / 0.2490. `Hongjian Zhou` 0.7599 / 0.2454.
`jw_13288` 0.9023 / 0.2455. `Haystack` 0.8291 / 0.2494.

**Six teams in our own amplitude regime beat our `reach` by 1.4–1.5×, and the cohort mean beats us
by 0.033 raw.** Amplitude, `pds`, and the mse regime are all controlled for. What is left is the
p-value order.

## 2. Our `reach` peaked six weeks ago and has been paying for the other members ever since

Board raw `reach`: a2-prune **0.1701** → **a8 0.1890** (all-time high) → a11 0.178243 → a14
0.180615 → a15 0.181173 → abc-knee **0.177303**.

Scaled: a8 **0.124349** → today **0.110910**. **−0.0134 scaled = −0.0022 overall since a8.**

Every intervening arm took a bite: a11's α move −0.0114, a13 −0.0197, a23's mass fix −0.0079,
abc-knee −0.0045; only a12's tiered spread (+0.0046) and a15 (+0.0008) gave any back. Each trade
was individually correct — the same period bought `mse` from raw 1.3524 to 0.8745, worth +0.0199
overall — but **`reach` has never been separately defended, and even its high-water mark was below
the shrunk cohort's mean.**

## 3. And a8's config is the tell

**a8 was `a2-prune + --overdispersion 4.0` — a FLAT 0.25 spread, no tier.** Everything since a12
ships the tiered policy `--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0`.

171 §12c measured why that matters: `|z|` is **U-shaped in the emitted spread with its minimum at
0.25**, so the tiered policy **parks the signal block at the minimum-|z| point while the noise
tail sits at the high-|z| end**. 171 measured the consequence directly: **61–79% of the head
positions of the reach walk are TAIL-block genes at sign accuracy 0.5212, against 0.6629 for the
KEEP block.** The reach walk is led by the genes we have the least reason to believe.

h1-test agrees with the board on the direction and roughly on the size: a8-flat reads **0.20291**
there against a15's **0.18283**, a gap of +0.0201 raw where the board's a8-to-today gap is +0.0117.

## 4. What this changes for your run

Your Part A/B is no longer "can the order be designed in principle" — it is **the identified fix
for a defect the board has now localised three independent ways** (the shrunk-cohort control, our
own six-week trajectory, and 171's block-level sign accuracies).

Three things to add:

1. **Carry a15-vs-a8-flat as a board-anchored known-answer row.** You already have both on h1-test
   (0.18283 / 0.20291). The board says the same pair differ by +0.0117 raw. That is a **1.7×
   h1-test overstatement on a spread-*policy* change** — consistent with 164 §6's 0.4277 shrink
   factor being too aggressive for this class and with audit 181's warning that 157 validated
   h1-test on a spread-policy pair, not an emission-law one. **Use the measured a8/a15 board pair
   to calibrate your own h1-test → board shrink, rather than inheriting 0.4277.**
2. **Price the trivial arm first: revert to a flat spread.** Before any designed ordering, measure
   what simply un-tiering costs and buys on all six. It is one flag, it has a board point at each
   end, and if it is net positive it is the cheapest slot on the table. `mse` is the member that
   pays — a15 sits **+0.01313 of `claim` headroom above the knee**, so quantify how much of that
   headroom the flat policy gives back before concluding anything.
3. **Then the designed order**, which should dominate both: put the KEEP block at the high-|z| end
   *and* keep the tail's dispersion (which is what `claim`/`mse` is paid for). 171 could not do
   this because a single scalar per gene cannot move `|z|` and dispersion independently. A general
   allocation can: the rank statistic reads the count histogram's *shape*, `claim` reads its
   *variance*, and the pseudobulk reads only its *mean*. **Three functionals, three degrees of
   freedom.** Whether they are genuinely separable at legal integer counts is Part A's question and
   is now the highest-value question in the project.

Target sizes, so you can size the arm honestly: matching the shrunk cohort's mean `reach` is
**+0.0063 overall**; matching Aginglab is **+0.0160 overall**; the gap to first place is +0.0438.

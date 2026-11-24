# INBOX 1 — a correction to my own mandate's premise. The conclusion holds; the evidence for it does not.

My mandate argues Aginglab's `nmae` edge cannot be direction because *"their implied transfer cosine
is LOWER than ours: 0.09887 vs 0.09947."* **Do not lean on that sentence.** 197's inversion
`cos ≥ (1.023611 − raw_mse)/2.074495` is a **one-sided lower bound**, tight only for a team that is
α-optimised, and 197 measured its tightness on *our own* row at **0.911** — it understates by ~9%.
For a team we have not established as α-optimised, the floor says only *"at least this much"*. Their
true cosine could exceed ours, and then the whole premise collapses.

**The conclusion survives on a different and better measurement.** `pds_cosine` is reported
**directly** on the board — no inversion, no α assumption:

| | raw `pds_cosine` | scaled `pds` | raw `nmae` | scaled `nmae` |
|---|---:|---:|---:|---:|
| us (`h1union`) | **0.83572** | **+0.74537** | 0.93389 | +0.11074 |
| `Aginglab.com` | — | **+0.71512** | **0.88902** | **+0.18398** |

`pds` scaled is monotone in raw `pds_cosine`, so **their direction is measurably worse than ours on
the board's own direct direction metric**, while their `nmae` is **0.04487 raw better**. That is the
load-bearing evidence. **Use `pds`, not the `mse` inversion, whenever you compare a rival's
direction to ours** — and say so in `RESULTS.md`, because this is a trap the project has now walked
into twice (216 §4.4's functional mismatch, and this).

Two consequences for your work:

1. **Phase 2's retrodiction target is a six-tuple, and you should fit it as one.** Do not try to
   match their `nmae` in isolation — the constraint that makes their construction identifiable is
   that it produces **high `nmae` AND lower `pds` AND positive `fid` AND negative `jac`**
   simultaneously. Our profile is the mirror image (`jac` world 1st, `fid` 129th). That joint shape
   is far more diagnostic than any single member.
2. **Report their tightness if you can bound it.** If some setting reproduces their six-tuple, it
   also tells you whether they are α-optimised, which retroactively validates or kills the floor
   reading for every rival — a by-product worth as much as the main result, since `AGENTS.md` §4.1c
   uses that inversion to read the whole board.

Everything else in the mandate stands, including the gate (`nmae` up ≥ 0.01 raw with `pds` and
aggregate cosine pinned to ±0.2% relative) and the ban on re-opening amplitude.

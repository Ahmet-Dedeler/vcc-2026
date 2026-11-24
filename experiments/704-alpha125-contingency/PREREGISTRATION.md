# 704 — carrier amplitude ×1.25 (654) as DN's contingency substitute (preregistration)

Written 2026-10-14 ~23:30Z by the science session [bf003b], on the coordinator [400b6e]'s request. At this time DN (671, 10-15
s1) is still building, so its entry does not exist. 704 is built and submitted ONLY IF DN reads NULL or HARM by its frozen
prereg (631's PREREGISTRATION-DNMAG.md). If DN reads TRANSFERS or WEAK-adopted, 704 is not built: the two are substitutes
(631 / 618: both lift called-DOWN magnitudes, and the DN pass on a ×1.25 base overshoots).

## The arm (one variable vs 630arm, entry QJRJK6OC6Zw1kjwnZqid, 0.2737556 = DN's own base)
- **704a125** = 630arm's exact env + `VCC_ALPHA=1.2275`. Under VCC_BULK_LOAD, VCC_ALPHA scales only the DE carriers and the
  ballast count (654 §3), and 0.982 × 1.25 = 1.2275.
- Built by run_arm_generic.sh against 630arm's meta. Expected config diff: `alpha` (+ `de_sign_emitter_sha256` if the pin
  moved). No re-level: 630arm's context factors are kept, because the knob moves the bulk only through ballast-K drift (654 §2:
  mean |Δ log1p CP10k| 2.3e-4).
- Same base as DN, so the two substitutes are read against one entry.

## Class and noise
**Scorer mechanics** (the carriers' design amplitude; nothing fitted to A/B/C). Same emission law; the ballast K moves, so
pass 2 is re-drawn: σ ≈ 0.00025 overall (537 §2.5). It also moves the bulk slightly (618: not bulk-neutral; bulk term
priced at +0.0005 on H1, indicative only).

## Held-out reads already on file (the 11.1c.2 requirement)
- H1 real emitter, production-faithful (654, 3 seeds): ×1.25 − ×1 = +0.0011 / +0.0019 overall-eq at f 0.5 / 0.35; nmae
  −0.0042 raw; head-DOWN delivered |lfc| +22 %.
- 618 (fresh h1-trainval surface, 126 unworn targets): (a) ×1.25 + 630 = +0.00360 at MIN(f 0.5, f 0.35) vs 630 alone
  +0.00269, so the increment over 630 is **+0.0009**. Plus the bulk term (+0.0005, indicative).
- 655's dose curve puts the peak at ×1.25 for f 0.5 and at ×1.5 for f 0.35. ×1.25 is the regime-robust dose; ×2.0 is past
  the peak and can exceed the emitter's per-cell capacity.

## Primary statistic, prediction, bins (board, vs 630arm's own entry)
- **Primary: Δoverall (all six members)**, since the knob moves both channels.
- **Prediction:** +0.0005 … +0.0015 overall. That is 618's +0.0009 DE increment through 611's call-set rule (nmae × 1.2–1.6,
  fid × ~2, jac × ~0.5), plus up to +0.0005 bulk.
- **TRANSFERS, ≥ +0.0008 (≥ 3σ):** adopt α ×1.25 for October in place of DN. The knob reaches def_day as VCC_ALPHA under
  BULK_LOAD; a pass-through and an inert KA are needed ([17adab]).
- **WEAK, +0.0003 … +0.0008:** adopt only if Δnmae ≤ −0.0010 raw (the mechanism's member).
- **NULL, −0.0003 … +0.0003:** not adopted.
- **HARM, ≤ −0.0003:** not adopted.
- **Falsifiers (stop and diagnose):**
  - more than 6 gated zero-total cells in the build's zero-gene audit (654: ×1.25 + zero repair leaves only own-target
    zeros);
  - a config diff beyond `alpha`;
  - the emitter's per-cell-totals assert (the ×2.0 failure mode).

## Tries
- Board on this switch: 1.
- H1: 654 (×0.8 / ×1.25, 3 seeds), 655 (×1.5 / ×2.0), 618 (×1.25 ± 630, 3 seeds).

## AMENDMENT 1 (2026-10-15T05:06Z, [bf003b]; written before any 704 build or entry): head-to-head with DN, a filler only
- DN read **PARTIAL** on 10-15 (0.2741292, +0.0003736 vs 630arm; Δnmae raw −0.001646, the other five members bit-identical),
  and October adopts `de_dnmag d_m200`. So 704's original trigger (DN NULL/HARM) did not fire.
- The coordinator [400b6e] (10-15 ~05:00Z): 704 is a **filler**, built only for an A/B/C slot that would otherwise be empty.
- **Decision variable:** Δ704 − ΔDN, both against 630arm's own entry (QJRJK6OC6Zw1kjwnZqid), in overall units.
  - DN keeps pass 1 and reads near-exactly, while 704 moves ballast K and re-draws pass 2. So the difference carries 704's
    re-roll only: σ ≈ 0.00025 (537 §2.5).
- **Swap rule:** October replaces `de_dnmag d_m200` by α ×1.25 (VCC_ALPHA=1.2275 under BULK_LOAD) iff
  **Δ704 − ΔDN ≥ +0.0005 (2σ)** and 704's own falsifiers are clean. Otherwise DN stays. A swap needs [17adab]'s VCC_ALPHA
  pass-through with an inert known answer and one rehearsal before the freeze. Without them, DN stays.
- Expected (amendment-time): 618's DE increment +0.0009 × A/B/C delivery ≈ 0.55 = +0.0005, plus up to +0.0005 bulk
  (indicative), against DN's +0.00037. So Δ704 − ΔDN ≈ +0.0001 … +0.0006: below or at the swap bar. Most likely outcome:
  DN stays, which is why 704 is only a filler.
- The original bins (vs 630arm) are still reported for the record. They no longer decide anything on their own.

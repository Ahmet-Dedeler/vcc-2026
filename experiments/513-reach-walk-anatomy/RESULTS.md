# 513 — Anatomy of our reach walk on H1: reach lives in tiny pools; no class reorder or donor-agreement subset builds a pure head

2026-09-28 12:40–13:05Z, coordinator. Analysis only (no build, no `vcc`, no `cell-eval2`). Inputs: 435's production-emitter
H1 stats (`stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz`, the shipped DE construction on H1) and 323's reference
table. Scripts: `walk513.py`, `keys513.py`, `block513.py`, `demote513.py`, `purity513.py`. Every process < 0.2 GB.

## Known answers
`walk513.py` recomputes the four DE members from the saved stats with the BRIEF's definitions: fid 0.5089291277653674,
jac 0.07469853420125332, nmae 0.9427040517722616, reach 0.26925697759272765 — all equal 435's metrics json (≤ 2e-16).

## Why this was asked
The board says three public-data teams moved reach fast with means pinned (Testh 09-24/25: +0.030 and +0.028 raw; chill
09-12: +0.034), and our reach (0.285) trails them (0.30–0.34). The question: is our reach limited by ORDER (fixable by a
fixed-mean free-rank edit) or by SIGNS?

## Results (H1, 82 targets with a non-empty pool)
1. **Reach lives in tiny pools.** n_real < 20 (32 targets): reach 0.643, walk sign accuracy 0.81. 20–100 (9): 0.078, 0.51.
   100–1000 (19): 0.027. ≥ 1000 (22): 0.014. The median first wrong sign in the walk is at position 4; median k* is 4.
   A/B/C's median n_real is 31–45 (AGENTS §2), so most A/B/C targets sit in the band where we get ~0.08.
2. **Sign accuracy decays with our rank** (targets with n_real ≥ 100): 0.80 at ranks 0–10, 0.67 at 10–25, 0.645 at 25–50,
   0.61 at 50–100, 0.59 at 100–250, 0.55 beyond. Our order carries confidence; it never reaches the 0.9 the walk needs.
   Ordering by the truth's own |lfc| with our signs is WORSE (0.234): our own order is already the better confidence key.
   Ordering by |our lfc| ignoring the sig-first rule: 0.2794 (+0.010, all from small pools; not deliverable as such).
3. **The high-CPM block crowds the head with worse signs** (block = top 2,200 by control CPM, called for every target with
   the shrunk CLEAN4 sign): in the 20–100 band the top-10 walk slots hold 69 block genes at 0.58 accuracy and 21 per-target
   genes at 0.95; at ranks 10–30 the block is at 0.458.
4. **But demoting the block loses** (analytic reorders keeping means, calls and signs, so every other member is pinned):
   Δreach −0.003 … −0.038 across block sizes 1,000–4,000 and three demotion schemes. The 0.95 per-target genes are only
   ~2 per target; behind them our per-target genes are at chance, so swapping block genes out brings worse ones in.
5. **No identifiable pure-head subset of useful size** (sign accuracy on pool genes; donors k562_gwps, orion_hct116,
   orion_hek293t, rpe1; H1 itself never used):
   - ≥ 3 donors agree and |consensus| ≥ 1 log2: 0.972, n 72 (almost all in pools ≥ 100: ~1.6 genes per target);
     ≥ 4 agree, |cons| ≥ 0.5: 0.908, n 261; ≥ 3 agree, |cons| ≥ 0.5: 0.830, n 576. The 20–100 band has 7 such genes.
   - CLEAN4 by confidence |P_down − 0.5|: 0.56 [0, 0.05) n 24,854; 0.68 [0.05, 0.1) n 8,508; 0.77 [0.1, 0.2) n 1,198.
   - Our emitted sign: 0.561 on all pool genes; 0.621 where ≥ 2 donors agree with it; **0.461 where ≥ 2 donors agree
     against it** (3,988 genes, 6 % of the pool). Flipping those is the consensus rider 477/480 already priced (≈ −0.003
     nmae on H1).

## Reading
Our reach is sign-limited, not order-limited. The fixed-mean free-rank edits that moved rivals' reach need a head of
≥ 0.9-pure genes several deep per target; on H1 our sources (donor rows, CLEAN4, our design) supply that only for the
strongest targets and only ~1–2 genes deep. The medium band (where A/B/C concentrate) has essentially no confident genes.
This agrees with 406/348/480 and closes "reorder the block" at its construction. What would change it is a better sign
source for weak responses (the 09-30 read of Arc's H1 programme is the live candidate).

Negative results recorded so they are not re-run: block demotion (all sizes/schemes), truth-|lfc| order, and the pure-head
subset search by donor agreement and CLEAN4 confidence.

## Addendum 513f (13:10Z) — called-gene magnitudes are at the H1 L1 optimum (`scale513.py`)
Analytic nmae with the delivered lfc of called genes multiplied by s (calls, signs, order unchanged; base 0.94270 = the
known answer): s 0.5 +0.0157, s 1.5 +0.0103, s 2 +0.0356, s 3 +0.1180. Block-only and rest-only scalings behave the same
way. Called gate genes: sign accuracy 0.600 (uncalled 0.537); where our sign is right, |ours|/|truth| median 0.435
(IQR 0.28–0.71). With 40 % wrong signs, s = 1 is the L1 optimum on H1 (the curve bottoms near s ≈ 1.1).

**Transfer caveat, and a candidate slot (not armed).** Magnitude levers have transferred at 1.6–2.3× their H1 price
(`docs/H1-BOARD-CALIBRATION.md`: MDE 1.76×, zero repair 3×). Under the L1 algebra a magnitude gain scales with (2a − 1),
so that ratio says A/B/C's gate sign accuracy is ~0.62 against H1's 0.56–0.60, and A/B/C's optimum s sits above 1. A
same-day pair s = 1.25 / 1.6 on the called genes' delivered lfc (nmae-only; bulk pinned by the ballast; reach may move a
little) would map it on the board. Price at H1 × 2 before arming; candidate for a slot after 10-05.

## Correction (2026-09-29, the 533 claims audit)
D8 (block demotion Δreach is −0.0007 … −0.0383 per work/demote513.json, not −0.003 … −0.038; every scheme still loses) and D9 (the A/B/C median n_real 31–45 is 192's board-inferred estimate, not an AGENTS §2 number). Source: `experiments/533-science-lead/RESULTS.md` §2.3.

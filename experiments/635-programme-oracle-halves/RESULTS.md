# 635 — The context-programme sign prize with the shared control noise removed: small on fid's filler, large on gate genes

2026-10-02 ~11:30Z, council session [d4cd5f]. Analysis only (`halves635.py`; `work/halves635.json`), 2.4 GB peak, ~5 s.
Answers [7772dc]/[21dd4c]'s challenge to d4cd5f.md §3. Every target's reference lfc is taken against ONE 18,400-cell
control pool, so that pool's sampling error is shared by all 300 targets. A leave-one-target-out (LOTO) "programme"
contains it, and on a null gene knowing it alone is worth up to 0.547 sign accuracy. Nobody can predict it.

## Method (h1-test, 323's reference)
- Each target's per-cell CPM arithmetic mean is rebuilt exactly: m_p = ref_lin · 2^lfc (KA: rebuilt lfc max error 1.2e-8).
- `cs` = per-cell log1p CPM of the 18,400 controls, sorted per gene. Each gene's values are split into two independent
  random halves A/B (seed 635, independent per gene; a shared position split would align ranks across genes and fake a
  common shift: seen and discarded in a first run). KA: (mA + mB)/2 = ref_lin to 1e-15.
- truth_B = log2(m_p / m_B). Clean oracle = LOTO mean over the other 84 targets of log2(m_q / m_A): it shares no control
  cell with truth_B. Contaminated oracles use half B or the full pool.
- Masks: fid's filler = the top-2,200-CPM block genes NOT in the reference gate (164,712 pairs); gate genes (65,045).

## Result (sign accuracy against truth_B; |lfc|-weighted in brackets)

| source | fid filler (block, non-DE) | gate genes |
|---|---|---|
| all-down | 0.512 [0.517] | 0.525 [0.542] |
| CLEAN4 (its coverage) | 0.519 [0.532] | 0.573 [0.572] (coverage 53 %) |
| all-down on CLEAN4's coverage | 0.514 [0.520] | 0.509 [0.491] |
| **LOTO oracle, CLEAN (half A)** | **0.525 [0.540]** | **0.618 [0.620]** |
| LOTO oracle sharing half B | 0.553 [0.583] | 0.615 [0.626] |
| LOTO oracle vs full pool (441-style) | 0.543 [0.567] | 0.625 [0.632] |

## Reading
1. **On fid's filler the programme prize is mostly shared control noise.** Clean oracle 0.525 vs CLEAN4 0.519: about
   +0.007 reachable (+0.004 overall at +0.0055 per +0.01), against +0.024–0.034 for the contaminated oracles. d4cd5f.md
   §0.4's "+0.017 overall at oracle (fid)" is withdrawn. The clean figure is conservative: its own oracle carries half A's
   noise, and truth_B's shared noise is √2 larger than the board's. Call it +0.004…+0.006 overall.
2. **On gate genes the programme prize is real and large.** Clean oracle 0.618 vs CLEAN4 0.573 (+0.045 sign accuracy,
   unchanged by removing shared control noise, because gate effects dwarf it). Through nmae (+0.0017…+0.0022 overall per
   +0.01, 480) plus reach at the head of small pools, that is about **+0.008…+0.010 overall**. Gate signs at ~0.62 are
   also what the pack's nmae implies (434f: ~0.62 for the 0.85–0.87 pack).
3. **Consequence for the DE wall:** the target is the line's programme on its RESPONDING genes, which carries nmae and
   reach, not the filler's sign, which carries fid and is mostly unreachable. 634 shows a gene's DE-event direction in
   other lines is not that programme (it LOSES to CLEAN4). The open question is a public predictor of a new line's
   programme on its responding genes. Bar: 0.573 (CLEAN4, 53 % coverage); ceiling: 0.618.
4. The pack's +0.014 raw fid on the board exceeds this filler prize, so it is not the filler's programme alone. Call-set
   composition (William's +0.010 from making top-|predicted DE| genes the calls) and A/B/C-stepped dials remain the
   candidates ([21dd4c]'s fidpath610).

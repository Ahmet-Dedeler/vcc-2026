# 321 — Correct the manuscript where 320's audit shows it overclaims, and fold in Round 7's board results. Rebuild the page. Do NOT publish.

You are a research subagent. Work in `experiments/321-writeup-correction/` **and** you may edit
`experiments/294-writeup/PAPER.md`, `CLAIMS.md` and `REFERENCES.md` — those three only. **Edit nothing else.**
No `vcc`, no `cell-eval2`, no build, nothing under `submissions/` or `data/`, **no publishing** (you have no
publishing tool; the coordinator republishes). `uv run python`. **Time box: 2.5 hours.**

## Why — the draft currently claims more than the evidence supports
**320** (`experiments/320-f16-novelty/RESULTS.md`) audited F16 adversarially and returned
**NEW 0 / EXTENDED 1 / KNOWN 4 / CONTRADICTED 0**. Its own summary: *"Nothing biological survives as new. F16 is a
re-measurement of a known phenomenon on a new sample — the challenge's own 2025 library. That is worth having, and
it is worth saying that it is what it is."* `docs/FINDINGS-2026-09.md` **v3.1** carries the re-grade; the draft
does not.

## The corrections, in priority order
1. **§3(g) and §4.1 claim the partner-construct control as ours. It is not.** 320 §3 quotes both sentences and
   supplies the weaker true version — **use 320's wording, do not invent your own.** The control is
   **Jackson et al., *Nat Biotechnol* 21(6):635–637 (2003)**, the founding observation of the seed off-target
   field, and it is already used for exactly this purpose on CRISPRi Perturb-seq (Hartman 2026) and as the
   discovery device in Rohatgi 2024.
2. **F16 must be presented as a RE-MEASUREMENT, with its provenance table**, not as a finding whose novelty is
   unstated. The phenomenon is **Rohatgi et al., *Cell Genomics* 4(11):100693 (2024)** — screens + RNA-seq **+
   dCas9 ChIP-seq occupancy**, 98–1,141 reproducible peaks per guide, **in our constructs' own parent
   hCRISPRi-v2 library, two years earlier** — and **Hartman et al. 2026**. The seed rule, its window, its 12-nt
   threshold **and its published low precision** are theirs too. Only the **design** is EXTENDED (a recombination
   of Cazares 2024's aggregate seed-shift test and Davis 2023's matched-null method, new only in using a
   *measured* 209,070-spacer propensity as the matching covariate).
3. **The statistical overreach: the headline p is a GENE-SAMPLING p, not a construct-level one** (320 §3(d)).
   Wherever the draft leans on "the construct as the experimental unit" to carry the p-value, say what the p
   actually resamples. Do not delete the construct-unit design — describe it accurately.
4. **Withdraw the wet-lab recommendation as written.** §4/§5 propose a ~$3–6 k sequence-mutant control as the
   decisive next experiment. Rohatgi 2024 already demonstrates the mechanism at the level of **binding**. Replace
   with an honest statement of what, if anything, a new experiment would add that the literature lacks — and if
   the answer is "little", **say that**.
5. **Fold in Round 7's board results** where the draft discusses the competition track (§3(h)): the per-target
   weighting family, its dose curve turning over at −1.5, the parity result, and **the new best entry
   `xuDelEik7zGWMXdNu0R4` at 0.2186597 (+0.0025653 over `cd4svd5`) from the suppression-only `cap` form**, plus the
   three axes 317/318 closed on the real six-member scorer. Keep it brief and factual; this is not the paper's
   subject.
6. `REFERENCES.md`: add Rohatgi 2024, Hartman 2026, Jackson 2003, Cazares 2024, Davis 2023, Goeman & Bühlmann 2007
   with DOIs from `experiments/320-f16-novelty/work/refs320.md`. `CLAIMS.md`: update the C8/C9 rows with the
   provenance verdicts.

## Rules
- **Every number you write must be read from the file it came from**; keep a log `work/verify321.log` naming the
  source file and line for each. **Do not quote a number measured on 313's synthetic fixture as a delivered
  value** (that error is Correction 5 in `docs/LOG.md`).
- **Prefer the weaker true sentence everywhere.** This edit removes claims; it does not add any.
- Rebuild the page to `experiments/321-writeup-correction/work/manuscript-v10.html` with
  `uv run python experiments/294-writeup/build_html.py <out>` (read that script first; update its masthead the way
  earlier rounds did). **Do not publish.**

## Deliverable
The three edited files, `work/verify321.log`, the rebuilt HTML, and `RESULTS.md` listing every sentence changed
(old → new), every reference added, and anything you judged too weakly supported to keep. Then stop.

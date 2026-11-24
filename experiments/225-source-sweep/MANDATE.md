# 225 — The exhaustive public-source sweep: is the leaders' climb just sources added one at a time?

You are a research subagent. Work in `experiments/225-source-sweep/`. No `vcc`, no submission build,
no production change. Peak RSS ≤ 25 GB, and **never while another agent is running `cell-eval2`** —
check first. `uv run python`. Downloads are authorised; put bulk data under `data/` and record
provenance (URL, accession, sha256, size, `LastModified`) for every file you fetch.

## Why this experiment exists — read this carefully, it is the whole idea

215 decomposed 707 board arm-pairs and found the leaders' direction **arrived in many small steps**:
`cqawesome` in **5** separately resolved steps, **largest only 28%** of the total; `Haystack` in 3;
**no leader shows a single step carrying their direction.** I have been reading that as iterative
model retraining (224 is testing that). **But there is a simpler reading that fits the shape better:
they added a public dataset, then another, then another.** Five sources, five small steps.

Two things make this the stronger hypothesis:
1. **A VCC participant argued publicly that data generation, not computation, is the bottleneck** —
   that with no training set it is "luck whether or not you have in-house data that's closer or
   farther from the test cell lines." Whoever believes that spends their effort on **sources**.
2. **We just re-ran the arm and it worked.** 223 priced **DLD-1 at +0.00228 — 1.20× the shape-only
   gate — with discrimination up on all three cancer-line folds, and both of its independent
   libraries separately positive.** That is a live source found on the *first* systematic look.
   223 only tested **9** candidates. **Nobody has ever swept the space.**

## READ FIRST
- `AGENTS.md` **§8.1, §8.2, §8.4** — the corrected **2.51×** multiplier and the **+2.5–2.9%**
  shape-only bar. **Price every candidate against that bar, not against zero.**
- `experiments/223-*/RESULTS.md` — the 9 candidates, the 7 deaths, DLD-1's two-library falsifier,
  and **the exact harness it used. Reuse that harness; do not write a new one.**
- `experiments/208-*/RESULTS.md` — **the rider rule**: a source must be re-priced against the
  *current* shipped stack, not against the stack it was measured on. 222 is shipping an H1 union
  right now, so **your baseline is h1union, not cd4f034.**
- `experiments/162-retrodicting-validator/` — PLOLO, the only pds instrument validated against a
  board point. `debench.Surface.pds()` must **never** size an arm (188 Phase 3).
- The **`solve_cd4_fraction` sizing trap** (`AGENTS.md`, and 218/223): the harness mirror at
  `experiments/162-retrodicting-validator/plolo.py:363-369` takes a median over **every** panel row,
  so a source with ≤50% coverage silently reads **f_per_unit = 0.0** — a real source scores exactly
  zero and looks dead. **Low-coverage sources are exactly what this sweep will find. Use the
  carrying-row median, and assert loudly on zero, or you will discard live sources as dead.**

## Work

### Phase 0 — known-answer rows (mandatory, before any new source)
Reproduce **223's DLD-1 number (+0.00228) and at least two of its seven deaths** through your own
harness. If you cannot reproduce them, stop and report — the harness is wrong and nothing after it
counts.

### Phase 1 — enumerate the space, exhaustively
Build a census of **every public perturbation dataset that shares targets with the 300-gene panel**.
Not the 9 that were convenient — all of them. At minimum sweep: GEO/SRA, ArrayExpress/BioStudies,
scPerturb, Perturb-seq/CRISPRi/CRISPRa depositions, the Tahoe release, Replogle's full set, Norman,
Adamson, Dixit, Datlinger, Papalexi, Frangieh, Jost, Nadig, X-Atlas/Pisces, CZI CELLxGENE, Broad
DepMap perturbation arms, and any 2025–2026 deposition. For each record: cell line, modality
(CRISPRi/a/KO/chemical), panel-target overlap count, cell count, and whether it is fetchable.
**Deliver the census table even for sources you cannot afford to fetch — it is a durable asset.**

### Phase 2 — price the fetchable ones, cheapest-per-overlap first
For each: measure transfer cosine gain on PLOLO against the **h1union** baseline, and price against
§8.4's **+2.5–2.9%** shape-only bar. Report the **six scored members**, not just cosine — 216 priced
one arm at −0.00642 *because nmae turns over*, and a source that lifts cosine can still lose overall.

### Phase 3 — falsifiers, all mandatory
1. **Two-library split** — the check that made DLD-1 credible. If a source has independent libraries
   or donors, both halves must be **separately positive**, and the two-split mean should beat both.
2. **Wrong-recipient null** — pair the source with the wrong held line. 184 died here: arm +0.0665,
   control +0.0667, null +0.0663, **identical to three decimals**.
3. **Row-matched control** — 223 killed a candidate this way *at every dose* because the arm in
   flight had already taken the shared direction.
4. **Marginal, not standalone** — a source that duplicates a direction we already hold is worth
   **zero**. Report each survivor's gain **on top of the full stack**, and the stack's gain with it
   **removed**.
5. **Leak audit** — `data/donors/h1_2025` covers **85/85** of h1-test's targets and is inside 222's
   shipped union. State every donor↔surface overlap for every source you keep.

### Phase 4 — the stack
If ≥2 sources survive, report the **best achievable stack** and whether gains are additive or
saturating. **Saturation is the key number**: if five sources give 5× DLD-1 the leaders' shape is
explained and we know what to do with the remaining 48 days. If they saturate at 1.5×, the
hypothesis is dead and 224's reading wins.

## Deliverable
`RESULTS.md`: the reproduced known-answer rows, the **full census table**, the priced candidates with
all six members, the four falsifiers per survivor, the stack curve, and a one-line verdict naming
which sources (if any) should ride the next slot.

**"The space is swept and DLD-1 is the only live source" is a complete result** — it would mean the
data-bottleneck reading is wrong for us and would hand the remaining time to 224's channel.

## Territory
221 owns cloud scoring, 222 owns the H1 build, **224 owns learned transfer — you own sources.**
224 may ask you for the LOLO substrate; share it. Do not fit models; do not touch the builder.

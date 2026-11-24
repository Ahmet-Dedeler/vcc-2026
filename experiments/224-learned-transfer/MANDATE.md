# 224 — Can ANY trained model beat the donor mean on held-line transfer? A scoped feasibility study.

You are a research subagent. Work in `experiments/224-learned-transfer/`; cloud GPU is authorised
(Ahmet has AWS credits and asked for them to be used). No `vcc`, no submission build, no production
change. Local peak RSS ≤ 10 GB — **the training goes to the cloud.** `uv run python`.

**Every instance you launch MUST carry a boot-time `shutdown -h +N` with
instance-initiated-shutdown-behavior=terminate**, and you must prove in `RESULTS.md` that nothing is
left running. On 2026-09-04 an orphan from experiment 061 was found **13 days** after that
experiment finished. Prefer spot. Report cost.

## READ FIRST — and know exactly what has already failed, because most of this space is closed
1. `AGENTS.md` **§8.1, §8.2, §8.4** (the 2.51× multiplier and the corrected **+11.6–13.6%** /
   **+2.5–2.9%** bars), §11.1, and `docs/CURRENT_STATE.md`'s header.
2. **The closures you must not re-run:** `195-retrieval-transform` (a learned transform, **negative
   in 12/12 cells at the cell-line boundary**: +0.00904 across unseen *perturbations*, **−0.00618
   across an unseen line**); `184-tahoe-transfer` (a basal-conditioned operator on **50** Tahoe
   lines — spectacular there, **absent in CRISPRi**: arm +0.0665 vs control +0.0667 vs
   wrong-recipient null +0.0663, identical to three decimals); `132`, `172`, `203`, `082`, `087`,
   `114`, `125` (small learned models and kernels, all NO-GO); `205` (X-Cell's repo raises
   `NotImplementedError` in both public functions).
3. `experiments/215-leader-decomposition/RESULTS.md` — **why this experiment exists.**
4. `experiments/159-transfer-ceiling/RESULTS.md` — κ² ≈ 0.34, the `1/k` law, and the ceiling any
   model is competing against.

## The question, narrowly

215 decomposed 707 board arm-pairs and found the leaders' direction **arrived in many small steps
inside the readable window** — `cqawesome` in **5** separately resolved steps, largest only 28%;
`Haystack` in 3; **nobody shows a single step carrying their direction.** 216 has now killed the
atlas explanation with byte-level evidence (**zero objects moved on that bucket between 2026-08-26
and 09-03**), and 218 killed the independence explanation.

215 named the remaining reading itself: *"iterative retraining of a learned transfer model, which
two hobbyist-shaped accounts with versioned model names would also produce."* **That is the last
standing hypothesis and this project has never seriously tested it** — every learned attempt here
was a low-capacity layer bolted **on top of** the consensus, never a model trained to replace it.

**So the question is not "build a foundation model". It is:**

> **Does ANY trained model beat the equal-weight donor mean on held-LINE transfer cosine, using the
> data we already hold — and if so, by how much, and does the advantage grow with training lines?**

## Work

### Phase 0 — the baseline and the protocol, before any model (mandatory)
- Assemble the leave-one-LINE-out substrate from the tables on disk: `k562_gwps`, `orion_hct116`,
  `orion_hek293t`, CD4, `h1_2025`, `nadig_jurkat`, `nadig_hepg2`, `rpe1`, `k562_essential`,
  `frangieh_melanoma`, `dld1_lowmoi`, `telohaec`, and anything 218/223 built. Report shared-target
  counts per pair.
- **The baseline is the equal-weight mean of the other lines**, scored as **transfer cosine** on the
  held line, macro over targets. Reproduce 159's κ² and the `1/k` curve as your known-answer rows.
- **Pre-register the gate before you fit anything:** a model must beat the donor mean by
  **≥ +2.9% relative cosine** (§8.4's corrected shape-only bar) **on held LINES**, not held
  perturbations. 195's failure was precisely that it cleared on perturbations and failed on lines.

### Phase 1 — the ladder, cheapest first, and stop at the first that fails to clear
Do **not** start with a transformer. Fit in this order, each leave-one-line-out:
1. **Ridge on target features → response**, features from what is public and already used here:
   ESM-2 embeddings, STRING, GenePT, DepMap co-essentiality (**132 did exactly this and failed —
   reproduce its failure as a known-answer row**, then say what you changed);
2. the same **conditioned on the recipient's basal profile** (184's shape, on CRISPRi);
3. a **small MLP / factorised bilinear** over (target embedding × recipient basal);
4. only if 1–3 show a *positive slope in training lines* — a larger model on GPU.

**The slope in `k` is the decision variable, not the level.** 184's whole result was that Tahoe's
advantage grew with lines and CRISPRi's did not. **If the advantage is flat or negative in `k` at
step 1–3, stop and report — a bigger model will not fix a flat slope, and that is the finding.**

### Phase 2 — falsifiers, all mandatory
1. **Wrong-recipient null**: pair each prediction with the wrong held line's basal. 184's arm,
   control and null read +0.0665 / +0.0667 / +0.0663 — **identical to three decimals**, and that is
   what killed it. Any model whose gain survives this is measuring amplitude, not transfer.
2. **Permutation over target labels.**
3. **Held-perturbation vs held-line, reported separately.** 195 was +0.009 on one and −0.006 on the
   other; a single pooled number hides the wall.
4. **Leak audit**: `data/donors/h1_2025` covers 85/85 of h1-test's targets, and CD4/`h1_2025` are
   inside 222's shipped union. State every donor↔surface overlap. PLOLO's challenge frame is the
   scoring instrument.

## Deliverable
`RESULTS.md`: the LOLO substrate, the reproduced 132 failure, the ladder with **cosine vs number of
training lines** for each rung, all four falsifiers, the cost, and proof nothing is left running.

**"No trained model beats the donor mean on held lines, and the slope in `k` is flat" is a full
result** — and given 195, 184, 132, 172 and 203 it is the likely one. It would close the last
standing explanation for the leaders' climb and mean the remaining 48 days belong to the source
stack (223's DLD-1 line) and to D/E/F robustness.

**Do not let this become a moonshot.** It is a feasibility measurement with a pre-registered gate
and a stop rule at every rung.

## Territory
221 owns cloud scoring infrastructure — **reuse its launcher and its shutdown guard rather than
writing your own.** 222 is building the H1 union. **You own learned transfer.**

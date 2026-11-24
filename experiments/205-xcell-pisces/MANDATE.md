# 205 — X-Cell and X-Atlas/Pisces: a foundation model trained on 16 contexts INCLUDING JURKAT, which is context A

You are a research subagent. Work ONLY in `experiments/205-xcell-pisces/` (downloads to
`data/external/xcell/` — that is the ONE exception).
No `vcc`, no `cell-eval2`, no submission build, no production changes.
Peak RSS ≤ 12 GB. **Download budget ≤ 120 GB; report the file inventory before exceeding 20 GB.**
`uv run python`. **Post nothing anywhere and put none of our method in any query or form**
(`AGENTS.md` §14).

## READ FIRST
`AGENTS.md` §5.3, §7 (**the rules on proprietary models and data — this mandate turns on them**),
§8.2, §11.1 rules 1–13. Then:
- `experiments/182-pretrained-state/RESULTS.md` **in full** — it closed pretrained State on a
  *structural* argument and its harness is directly reusable here. **Its gates are your gates.**
- `experiments/166-donor-acquisition/RESULTS.md` — the acquisition sweep that found Pisces empty on
  2026-09-02, and still empty today (`used_storage: 0`, unchanged since 2026-03-25).
- `experiments/197-cosine-inversion/RESULTS.md` — our cosine is **0.079249**, the leading direction
  is **≥1.66× ours**, worth **+0.00365 overall per +0.01 of cosine**, and the leader's level implies
  **~10.9 effective donor lines** against our 3.09.
- `experiments/184-tahoe-transfer/RESULTS.md` — cross-line adaptation is closed on 50 lines.

## Why this is the highest-value open lead

**Every route to a better direction has closed except one, and the reason is always the same: three
cell lines.** 166 proved no fourth panel-covering CRISPRi screen exists publicly. 184 closed
adaptation. 195 closed learned transforms *at the cell-line boundary* (+0.009 across unseen
perturbations, −0.006 across an unseen line, negative in 12/12 cells). 193 closed the operator.
194 closed the estimators. 197 says the leader sits at ~10.9 effective lines and **we cannot get
there by donor count.**

**Xaira released, in March 2026, exactly the thing that would break that constraint:**

- **X-Atlas/Pisces** — *"25.6 million perturbed single-cell transcriptomes across **16
  biologically diverse contexts**, including widely used cell lines (HCT116, HEK293T, **HepG2**),
  induced pluripotent stem cells (iPSCs), **resting and CD3/CD28 activated Jurkat T lymphoma
  cells**, and multi-lineage differentiating iPSCs"* (dataset card).
  **Jurkat is context A** (experiment 041 identified A = JURKAT, B = HELA, C = CAL33 under all eight
  feature/rank variants). Pisces also covers HepG2, which our essential-gene tables cannot serve.
- **X-Cell** — *"Scaling Causal Perturbation Prediction Across Diverse Cellular Contexts via
  Diffusion Language Models"*, bioRxiv **2026.03.18.712807**, Xaira's virtual cell model trained on
  Pisces. Press coverage states Xaira is *"making a subset of the Pisces dataset and X-Cell model
  available to the scientific community."*

## Work — four gates, in order. Stop at the first that fails and write it up.

### Gate 1 — AVAILABILITY (cheap, do this first)
- Is the **Pisces subset** actually downloadable anywhere: the HF repo (still `used_storage: 0`),
  a second HF repo, Figshare+ (Orion's processed data lives there —
  `plus.figshare.com/articles/dataset/…/29190726`), GEO, Zenodo, or the paper's data-availability
  statement? **Read the bioRxiv paper's data and code availability sections verbatim.**
- Are **X-Cell weights and inference code** available, and where?
- **Report exact URLs, sizes, file lists and shas.** If nothing is available, say so and stop —
  that is a complete result and closes the lead cleanly.

### Gate 2 — LICENCE (this can kill it regardless of the science)
Orion and Pisces are **CC-BY-NC-SA-4.0** (non-commercial, share-alike). `AGENTS.md` §7 says
proprietary models and data are allowed *"if you have the right to use them"*, and 182 flagged the
same question for State: **whether submitting output derived from a non-commercial-licensed model or
dataset to a $100,000-prize competition is inside the grant is a genuine open question, and
share-alike adds a second condition.** Read the actual licence text and X-Cell's own licence.
**Report what the licence says and what the open question is. Do not resolve it yourself and do not
argue it away — flag it for Ahmet as a decision he makes, exactly as 182 did.**

Note the asymmetry that makes the work worth doing anyway: **using it to measure a ceiling is not
the same as shipping its output.** Even if shipping is blocked, a measurement tells us whether 16
contexts closes the 1.66× direction gap, which decides how the remaining time is spent.

### Gate 3 — THE STRUCTURAL TEST THAT KILLED STATE
182's decisive finding was structural, not empirical: every published Replogle State checkpoint
conditions on the perturbation via `pert_rep: onehot` over a frozen **2,024-gene vocabulary** whose
intersection with the VCC 2026 panel is **0 / 300** — and a one-hot has no out-of-vocabulary path.

**Ask the same question of X-Cell first, before any inference run.** How is the perturbation
represented? If it is a one-hot over a fixed vocabulary, what is the vocabulary and what is its
intersection with `data/challenge/pert_counts.csv`? Pisces targets were genome-wide in Orion
(18,903 protein-coding genes), so the intersection may be near-complete — **but confirm it from the
checkpoint's own config, not from the paper's prose.** If there is no out-of-vocabulary path and the
intersection is small, the lead dies here as State's did.

Also check how the **context** is conditioned. Zero-shot into an unseen line requires conditioning on
the basal state; if X-Cell conditions on a context *identifier*, it cannot be pointed at D/E/F.

### Gate 4 — MEASUREMENT, with 182's gates
Reuse 182's harness. Its structure is what makes the result readable and you must keep all of it:
- an **in-distribution row** (a context X-Cell trained on) to prove the inference path works;
- the **held-line row** (h1-test — donors `k562_gwps` + `orion_hek293t` only; `data/donors/h1_2025`
  covers 85/85 of h1-test's targets and is pure leakage);
- **our own consensus on the same pool** as the comparator (182 measured 0.5787 sign accuracy,
  0.9500 discrimination on its pool);
- an **excess-over-context-mean-null** column — 182's whole point was that Arc's published held-line
  `de_direction_match` of 0.67 *ships without a null*, and its own K562 row put that null at 0.5200;
- **182's promotion gates: sign accuracy ≥ 0.56 and pds ≥ 0.55.**
- Report **degeneracy**: 182 checked that State's H1 predictions were *not* degenerate (20 distinct
  effects at full magnitude), so State was confidently wrong rather than silent. Do the same.

Then convert to our units: the **implied transfer cosine** by 197's route, and **|t|-weighted `q_w`**
(187 §5). The bar that matters is §8.2's: **+31% cosine for a single arm to clear rule 13's gate.**

**And if Pisces data (not the model) is available, that is the better prize.** A Jurkat or HepG2
donor table is a genuine fourth and fifth panel-covering line — 197 says each independent line moves
the effective count that everything else is capped by. Run 166's `census.py` on it against
`pert_counts.csv` and report panel coverage per context. Build a donor table only if coverage is
material, and never modify an existing file under `data/donors/`.

## Deliverable
`RESULTS.md`: the availability inventory with URLs and shas, the licence position stated as an open
question for Ahmet, the structural test, the measurement table with 182's null columns and gates,
and the conversion into cosine / `q_w` / overall units. **A clean "nothing is available" or "the
vocabulary does not cover the panel" is a full result.**

## Territory
196 emitted depth · 200 shipped-stack ablation · 202 the call-set build · 204 donor mixing.
**You own X-Cell, Pisces, and any new external model or atlas.** Do not touch production.

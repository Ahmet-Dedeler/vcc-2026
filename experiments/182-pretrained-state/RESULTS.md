# 182 — Arc's published pretrained State Transition checkpoints: NO-GO as a sign donor, and structurally undeployable on this panel

**Date:** 2026-09-03 · CPU, deterministic (seed `20260903`) · **no submission, no `vcc`, no
`cell-eval2`, no build.** Nothing outside `experiments/182-pretrained-state/` and
`data/external/state-pretrained/` was written. Total download **471 MB** (cap 40 GB).

---

## Verdict

**NO-GO. Public pretrained State is closed as a donor**, on two independent legs, the first of
which does not depend on any accuracy measurement:

1. **It cannot be conditioned on a single one of our targets.** Every published Replogle ST
   checkpoint uses `pert_rep: onehot` over a frozen **2,024-perturbation vocabulary**, and the
   intersection of that vocabulary with the VCC 2026 panel is **0 / 300**. A one-hot has no
   representation for an out-of-vocabulary gene, so there is no inference path to a challenge
   target at all — not a weak one, none. (The *competition* State recipe avoids this by
   conditioning on ESM2 features of the perturbed gene; the *published pretrained* weights do not.)
2. **On a line it never saw, its per-gene signs are at chance**, while the same checkpoint on a
   line inside its training set is clearly informative. On the h1-test reference DE pool it scores
   **0.4959 / 0.5135** against our production consensus's **0.5787** on the identical pool, and its
   excess over its own context-mean null is **negative** in both arms.

Both verdict-rule gates fail: sign accuracy ≥ 0.56 → measured **0.51**; `pds`-analogue ≥ 0.55 →
measured **0.500 / 0.571** against our consensus's **0.950** on the same 20 targets.

This is *not* a re-run of 106/113. Those used a small locally-trained `state_sm`. This is Arc's
published multi-line checkpoint, and it is demonstrably a working model — which is what makes the
H1 collapse informative rather than a null result about our own plumbing.

---

## 0. Known-answer rows (AGENTS.md §11.1 rule 10)

Both reproduce, so the numbers below are readable.

| known-answer row | expected | measured |
|---|---|---|
| production consensus, h1-test pool, all-three-carrying | **0.5697** at n=39,945 (exp 175) | **0.5697** at **n=39,945** |
| K562 donor alone, same mask | **0.5670** (exp 175) | **0.5670** |
| ST on a line **inside** its training set beats chance by a lot | ≫ 0.50 | **0.5962**, excess **+0.0762**, discrimination **0.8640** |
| checkpoint loads with a strict state-dict load | no shape errors | passed |

The third row is the one the mandate asked for, and it is the load-bearing one: **the inference
path works.** On K562 — a line in `zeroshot/hepg2`'s training set — the checkpoint recovers
per-gene direction at 0.5962 against its own context-mean null of 0.5200, and ranks each target's
own real effect above 86.4% of the others. Nothing about the pipeline is broken.

---

## 1. Inventory (HuggingFace `arcinstitute/*`, checked 2026-09-03)

| repo | what it is | pert rep | output space | held-out | size | licence |
|---|---|---|---|---|---:|---|
| `ST-HVG-Replogle` | ST on HVG input, `zeroshot`/`fewshot` × 4 lines | onehot 2,024 | 2,000 HVGs | **a cell line** | 57 GB | `MODEL_LICENSE.md` |
| `ST-SE-Replogle` | ST on SE-600M embeddings | onehot 2,024 | gene | **a cell line** | 95 GB | `MODEL_LICENSE.md` |
| `st-x-replogle-full` | per-line data-scaling sweep (`*_0.99`) | onehot 2,024 | 6,546 genes | **only perturbations, not a line** | 39 GB | **none in repo** |
| `st-se-replogle-full` | as above, SE input | onehot 2,024 | gene | as above | 46 GB | none in repo |
| `ST-*-Tahoe`, `ST-*-Parse` | other atlases | — | — | — | up to 1.5 TB | `MODEL_LICENSE.md` |
| `SE-600M`, `SE-100M` | State *Embedding*, not a transition model | n/a | n/a | n/a | 29 GB | `MODEL_LICENSE.md` |

**The mandate preferred `st-x-replogle-full`; I did not use it, for two reasons.** Its `config.yaml`
shows `toml_config_path: /data/replogle_llm/{line}.toml` with `downsample: 0.99` and tag
`replogle_downsampling` — each `*_0.99` model is trained on **one line** and evaluated on held-out
*perturbations* of that same line. It has no held-out **line**, so it cannot answer the question.
It also ships **no licence file at all**. `ST-HVG-Replogle/zeroshot/*` is the real cross-line model,
carries the licence, and Arc publishes its own held-line eval alongside it.

**Chosen checkpoint** (one, as instructed):
`data/external/state-pretrained/ST-HVG-Replogle/zeroshot/hepg2/checkpoints/best.ckpt`
sha256 `fe6d7fc71948509117cf6ed2207a9c4136023b54ad05803dda507b87d7468c42`
(trained on jurkat + k562 + rpe1; hepg2 held out — so **K562 is in-train**, giving the known-answer row).

### Licence

`MODEL_LICENSE.md` is the **Arc Research Institute State Model Non-Commercial License**
(last updated 2025-06-23). §2.1 grants use "solely for Non-Commercial Purposes"; §1.4 defines that
as a purpose "not undertaken for direct or indirect monetary compensation, financial consideration,
or any other commercial undertaking or advantage."

- **This experiment is clearly inside the grant**: evaluation and research by an individual, no
  submission, no compensation.
- **Shipping State-derived output into the challenge is a genuine open question, not a settled
  yes.** The challenge awards **$100k / $50k / $25k**, so a prize-seeking submission is arguably
  "undertaken for direct monetary compensation." Cutting the other way, Arc is the licensor *and*
  the organiser, and its rules explicitly permit State while stating it is **"not able to win"**
  (`AGENTS.md` §7) — which, read plainly, removes the prize and therefore the conflict, but also
  removes the reason to use it.
- **This is moot for us**, because leg 1 above closes deployment regardless. I flag it only so the
  question is not silently inherited if a future experiment fine-tunes State weights.

---

## 2. Arc's own held-line numbers, and why they do not transfer

The repos ship tiny per-perturbation eval CSVs (2.9 MB) for the held-out line. Read before
downloading any weights:

| held-out line | `de_direction_match` | `discrimination_score_cosine` | n perts |
|---|---:|---:|---:|
| hepg2 | 0.6766 | 0.7247 | 1,340 |
| k562 | 0.6688 | 0.7071 | 1,383 |
| jurkat | 0.6926 | 0.6871 | 1,537 |
| rpe1 | 0.6737 | 0.6871 | 1,499 |

At face value 0.67–0.69 direction match on a held-out line looks far above our 0.537. **It does not
survive contact with our setting**, and the reason is the whole finding: those four lines are the
*same platform, same lab, same protocol, same gene panel*. "Held-out line" inside Replogle is a much
smaller shift than Replogle → H1 hESC. Our own measurement below is the like-for-like version, and
it lands at chance. Note also that our LODO consensus scores **pds 0.885** on held-out K562
(161 §0) against ST's 0.707 discrimination on the same held-line task — even on Arc's home turf the
published checkpoint does not beat our incumbent at retrieval.

---

## 3. Method

`prep.py` → `infer.py` → `score.py`, all in this directory.

- **Preprocessing** reproduces `state tx preprocess_train` (`sc.pp.normalize_total` → `log1p`),
  except HVG *selection* is replaced by the checkpoint's own frozen 2,000-gene list in its own
  order; genes absent from a source are zero columns (H1: 1,859/2,000 present; K562: 2,000/2,000).
- **Inference** is the official `state._cli._tx._infer.run_tx_infer`, not a hand-rolled forward.
  One compatibility shim: the installed `transformers` rejects `hidden_size=328` with 12 heads,
  but this checkpoint carries an explicit `head_dim=64` (q/k/v project 328 → 768), which is legal
  for Llama attention. Only that validator is disabled; the strict state-dict load still passes,
  so a real shape error would still fail loudly.
- **Effect** = mean(predicted target cells) − mean(predicted control cells), both model outputs, so
  the model's own bias cancels. 2,048 real control cells as the basal pool, 128 rows per target.
- **Truth**: H1 = the h1-test reference DE pool (`164/de164.a15.npz`, `real_sig & adj`), the exact
  object production's 0.537 is measured on. K562 = the measured `k562_essential` pseudobulk delta,
  pooled to each target's top-200 genes by `|delta|`.
- **The 20 targets are the h1-test targets the vocabulary can express** (20 of 85). H1 and K562 use
  the *same* 20 targets, so the two rows differ only in the recipient line.
- `h1_2025` was never loaded. The only H1 data used is **non-targeting control cells** from the
  2025 test split — controls, not perturbation deltas, so no oracle leakage.

---

## 4. Result

`work/scores.json`. CIs are 95% bootstrap over **targets**, the axis the metric averages over.

**HELD ROW — H1, a line the checkpoint never saw** (20 targets × 1,798 shared genes, 6,755 pool pairs):

| predictor | sign acc | 95% CI | context-mean null | excess | discrimination |
|---|---:|---|---:|---:|---:|
| **production consensus (K562+Orion)** | **0.5787** | [0.5436, 0.6195] | 0.5497 | **+0.0290** | **0.9500** |
| K562 donor alone | 0.5778 | [0.5514, 0.6023] | 0.5602 | +0.0176 | — |
| ST-HVG zeroshot, own-median norm | **0.4959** | [0.4789, 0.5126] | 0.5038 | **−0.0079** | 0.5000 |
| ST-HVG zeroshot, scale-matched | **0.5135** | [0.4865, 0.5346] | 0.5184 | **−0.0049** | 0.5711 |

**KNOWN-ANSWER ROW — K562, in-train, same 20 targets, same code path:**

| predictor | sign acc | 95% CI | context-mean null | excess | discrimination |
|---|---:|---|---:|---:|---:|
| ST-HVG zeroshot/hepg2 | **0.5962** | [0.5376, 0.6715] | 0.5200 | **+0.0762** | **0.8640** |

The contrast is the result. Identical checkpoint, identical inference code, identical targets —
**excess +0.0762 and discrimination 0.864 in-train, versus excess −0.005 and discrimination 0.50–0.57
on H1.**

**The H1 predictions are not degenerate — the model is confidently wrong, not silent.** All 20
predicted effects are distinct, with `median |delta|` 0.023 (own-median) and 0.052 (scale-matched)
against 0.027 on K562, and row norms in the same range. ST emits differentiated, full-magnitude
effects on H1 that simply carry no information about the truth. That rules out "the basal input was
out of distribution so the model output nothing" as an explanation.

**The scale-matched arm rules out the obvious confound.** H1's median library size is 54,830 against
Replogle's 13,925, so scanpy's default per-dataset normalisation puts H1 inputs ~4× higher before
`log1p` than anything the model saw. Renormalising H1 to the training line's library size moves sign
accuracy 0.4959 → 0.5135 and discrimination 0.500 → 0.571 — real, and still chance. The NO-GO is not
a preprocessing artefact.

---

## 5. Mechanism

Three things, in order of how much they close:

1. **One-hot conditioning has no out-of-vocabulary path.** The perturbation enters only as a row of
   a 2,024-way one-hot. 0/300 challenge targets and 20/85 h1-test targets are in it. This is a
   property of the released artifact, not of State the architecture, and it is why the competition
   recipe uses ESM2 perturbation features instead.
2. **What the checkpoint learned is largely line-specific.** It ranks its own targets at 0.864
   within its training distribution and 0.50–0.57 outside it. That is the same wall 159/160/161 hit
   from the other side: 161 measured 60.1% shared / 39.9% line-specific energy and capped a
   context-independent predictor at cosine 0.77 per line. A model that gets its recipient signal
   only from basal control cells has no more to work with than our donor tables do, and it has seen
   four lines where our consensus has effectively seen its donors far more deeply.
3. **The negative excess is the sharp part.** ST on H1 does not merely fail to beat our consensus —
   it fails to beat *its own across-target mean*. Whatever survives the domain shift is a common
   response axis with no target-specific content, exactly the artefact AGENTS.md §4's `fid`
   baseline (0.505–0.522) exists to expose. Arc's 0.67 held-line `de_direction_match` is published
   without such a null; our K562 row shows its own null sits at 0.5200, so a large part of a raw
   0.67 is the shared axis, not identity.

---

## 6. Consequences

- **Close public pretrained State as a donor.** Do not download the Tahoe/Parse checkpoints, do not
  blend, do not spend a slot. With 106/113 (locally-trained `state_sm`) and this (Arc's published
  multi-line weights), both ends of the State family are now measured and both are NO-GO.
- **`docs/CURRENT_STATE.md`'s §"106/113 … NO-GO" can be strengthened**: the caveat it carried —
  "treat the larger current-source `state` profile as a separate hypothesis requiring a new
  cost/benefit decision" — is now discharged. It was tested and it is dead.
- **The 0.537 sign problem is untouched and remains the one problem.** Nothing here moves it. The
  only State-shaped idea this leaves alive is *training* an ESM2-conditioned ST ourselves on our
  own donor lines, which is the competition recipe, not a pretrained donor — and 106 already
  measured that at chance on a small budget. It would need a real training budget and would be
  arguing against leg 2 above, which says the transferable part is the shared axis we already have.
- **Reusable**: this directory now holds a validated, 4-second control-only ST inference harness on
  the official code path, with two passing known-answer rows. If a future experiment wants to test
  any ST checkpoint, `prep.py`/`infer.py`/`score.py` take a new one by changing `MODEL_DIR`.

## Files

| file | what |
|---|---|
| `prep.py` | builds control-only inference inputs (H1, H1 scale-matched, K562) |
| `infer.py` | official `state tx infer` path + the `transformers` validator shim |
| `score.py` | sign accuracy vs context-mean nulls, bootstrap CIs, discrimination, known-answer rows |
| `work/scores.json` | every number in §4 |
| `work/prep_meta.json` | seeds, cell counts, normalisation target sums, gene coverage |

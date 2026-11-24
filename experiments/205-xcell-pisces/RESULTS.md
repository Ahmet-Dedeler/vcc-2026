# 205 — X-Cell and X-Atlas/Pisces: the artifacts do not exist. Gate 1 fails; gate 3 would have passed.

**Date:** 2026-09-03 · deterministic (seed `20260903`) · **no `vcc`, no `cell-eval2`, no
submission, no build, no production change.** Nothing written outside
`experiments/205-xcell-pisces/` and `data/external/xcell/`. Total network download **48 MB**
(cap 120 GB; the 20 GB reporting threshold was never approached). Peak RSS ~2.6 GB.

---

## Verdict

**NO-GO on availability, and the lead closes cleanly — but for a completely different reason than
State's, and the difference matters for what a future agent should do.**

1. **Gate 1 FAILS. Neither artifact exists.** `Xaira-Therapeutics/X-Atlas-Pisces` is
   `usedStorage: 0` with 3 files (`.gitattributes`, `LICENSE.md`, `README.md`), unchanged since
   **2026-03-25**. `Xaira-Therapeutics/X-Cell` holds **272,925 bytes** — a `.gitattributes`, a
   README, and one PNG. **There are no weights.** The GitHub repo
   `Xaira-Therapeutics/X-Cell` is a documentation stub: two commits, both 2026-03-16, no releases,
   no tags, empty `configs/`, empty `scripts/`, empty `tutorials/`, and a single 5.2 KB
   `src/xcell/model.py` in which **both `from_pretrained` and `predict` raise
   `NotImplementedError`.** `pip install xcell` installs an unrelated package by another author.
   Re-verified 2026-09-03T16:17:55Z.
2. **The paper's own availability statement is future tense for the data and false for the model.**
   §5.1: *"The following data **will be uploaded**…"* §5.2: *"we **have made** X-Cell (55M) model
   artifacts publicly available. The trained model weights are hosted on the Hugging Face Model Hub
   at `https://huggingface.com/Xaira-Therapeutics/X-Cell`"* — the repo has held one PNG for
   5.6 months, and the URL in the paper is a typo (`.com`).
3. **Gate 2 is a live question for Ahmet, not for me.** Both artifacts are CC-BY-NC-SA-4.0.
   NonCommercial + **ShareAlike**, and share-alike is a second condition State's licence does not
   carry. §2 below states it and does not resolve it.
4. **Gate 3 WOULD HAVE PASSED, and this is the finding worth carrying forward.** X-Cell does **not**
   use a one-hot over a training vocabulary. The perturbation enters as
   `C_pert ∈ R^{6×d}` — cross-attention over ESM-2, STRING, GenePT, DepMap, JUMP-Cell Painting and
   a gene-identity embedding, all **keyed by the perturbed gene**, over a vocabulary defined as
   *"19,082 protein-coding genes … the intersection of GENCODE annotations and the ESM-2
   vocabulary."* Measured here: **300/300** of the panel are HGNC protein-coding and **298/300**
   carry a STRING embedding. Context enters as **control cells**, not a context identifier, with
   explicit test-time adaptation on unlabeled controls. **X-Cell is structurally pointable at D/E/F
   and at all 300 targets.** State was not, on either count.
5. **Gate 4 cannot be run and no substitute inference is possible.** There is no checkpoint. In its
   place I measured the thing that actually decides how remaining time is spent: **what the
   promised Pisces upload would be worth if it landed.** Answer: the two Jurkat arms are the only
   panel-relevant content, they would carry **≥292/300** of the panel, and one extra effective donor
   line is worth **+21.3% cosine / +0.0062 overall** — **below §8.2's +31% single-arm bar.**

**The one-line strategic read:** the only external asset that could move our direction is real,
is designed exactly right for our problem, and **is not downloadable.** Its value if it lands is
also now priced, and it is a fifth of the gap, not the gap.

---

## 0. Known-answer rows (AGENTS.md §11.1 rule 10) — all reproduce

`gate3.py` **aborts** if any KA-1/KA-2 row misses.

| row | expected | measured | source | verdict |
|---|---:|---:|---|---|
| panel overlap `orion_hek293t` | 294 | **294** | 166 §1 | OK |
| panel overlap `orion_hct116` | 290 | **290** | 166 §1 | OK |
| panel overlap `kolf2` | 273 | **273** | 166 §1 | OK |
| panel overlap `k562_gwps` | 269 | **269** | 166 §1 | OK |
| panel overlap `h1_2025` | 20 | **20** | 166 §1 | OK |
| panel overlap `nadig_jurkat` / `rpe1` / `k562_essential` | 0 / 0 / 0 | **0 / 0 / 0** | 166 §1 | OK |
| 159's law at `k = 3.09` | surface cos 0.03814 | **0.03813** | 197 §4.2 | OK |
| 159's law at `k = 4.00` | surface cos 0.04278 | **0.04278** | 197 §4.2 | OK |
| **negative control:** the 300 panel targets must beat a random 300-gene draw on basal expression in context A (§2.1 selected on it) | ≫ | **299 vs 160.1 ± 9.5 at CPM ≥ 5 (+14.6σ)** | AGENTS.md §2.1 | OK |
| **calibration:** Xaira's own 99%-of-counts library rule applied to context A must land near Pisces' measured Jurkat library | 10,872 | **9,114 (−16.2%)** | preprint Table 1 + §4.2.1 | OK, one-signed (§4.1) |

The last row is the one that licenses §4's coverage number, and its miss is **safe**: my proxy
under-counts the library, and a larger library can only carry *more* of the panel. Contexts B
(HELA) and C (CAL33) land at 8,825 and 9,346 under the same rule, so the shortfall is my readout
axis, not Jurkat.

---

## 1. Gate 1 — AVAILABILITY. Nothing is downloadable.

### 1.1 The two Xaira artifacts, read from the APIs

| artifact | URL | files | bytes | last modified | licence |
|---|---|---:|---:|---|---|
| **X-Atlas/Pisces** (dataset) | `https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces` | 3 | **`usedStorage: 0`** | **2026-03-25T22:58:18Z** | cc-by-nc-sa-4.0 |
| **X-Cell** (model) | `https://huggingface.co/Xaira-Therapeutics/X-Cell` | 3 | **272,925** (one PNG) | **2026-03-16T23:38:37Z** | cc-by-nc-sa-4.0 |
| **X-Cell** (code) | `https://github.com/Xaira-Therapeutics/X-Cell` | 26 blobs, 1,138 KB | — | **pushed 2026-03-16T23:57:26Z** | CC BY-NC-SA 4.0 |

Pisces' entire content is `.gitattributes` (2,504 B), `LICENSE.md` (20,848 B), `README.md`
(956 B). X-Cell's is `.gitattributes` (1,575 B), `README.md` (3,629 B), `x-cell-overview.png`
(272,925 B, sha256 `094851946e94c9d8b08582c4cead82d5231c0ea03f44a69ae582295bb17b82d1`).

**The Pisces README, verbatim:**
> **(Coming Soon)** The following data will be uploaded to this dataset:
> 1. All cells from X-Atlas/Orion (cells with a valid dual-guide pair from the **HCT116** and
>    **HEK293T** screens)
> 2. Cells containing a valid dual-guide pair AND perturbations in the 200-gene validation set OR a
>    non-targeting control from the **HepG2** and **iPSC** screens
> 3. All cells with a valid dual-guide pair from the **Jurkat Resting** and **Jurkat Active** screens

**The X-Cell README, verbatim:**
> **Status: Model weights and inference code coming soon.** The Python API, model weights, and
> tutorials are under active development.

### 1.2 The GitHub repo is a stub, and its own code says so

`src/xcell/model.py` is the whole package:

```python
raise NotImplementedError(
    "Model loading is not yet implemented in this release. "
    "Full inference code is coming soon — watch the repository for updates.")
...
raise NotImplementedError("Inference is not yet implemented in this release.")
```

`configs/`, `scripts/`, `tutorials/` and `docs/` contain only `.gitkeep`. **Zero releases, zero
tags, two commits** (`6a84b66e` "X-Cell initial public release CC-BY-NC-SA-4.0" and `7195c647`
"fix repo url", both 2026-03-16). The `gh-pages` branch is the built mkdocs site of the same docs.
All **6 forks** are byte-identical at the same commit; none has pushed weights.

**Three open issues asking exactly this question have received no maintainer reply in 5.4 months:**

| issue | opened | asks | maintainer replies |
|---|---|---|---:|
| #2 Yuhang-King | 2026-03-24 | "Expected release timeline for X-Cell Mini (55M) model weights?" (bumped 2026-04-08 and 2026-06-06) | **0** |
| #3 zzheng18 | 2026-04-13 | "what is the release schedule of the full Pisces atlas across 16 contexts?" | **0** |
| #4 jw310-wjr (Harvard Med) | 2026-06-01 | the 380 held-out HepG2 split file | **0** |

### 1.3 The preprint's own availability sections, verbatim

bioRxiv **10.64898/2026.03.18.712807**, posted 2026-03-20, Europe PMC `PPR1217946`. *(The DOI
prefix is `10.64898`, **not** `10.1101` — `api.biorxiv.org/details/biorxiv/10.1101/2026.03.18.712807`
returns "no posts found", which is how this paper looks missing if you query it the old way.)*

> **5.1 Data Availability.** The following data **will be uploaded to**
> `https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces`: (1) all cells from
> X-Atlas/Orion (cells with a valid guide pair from HCT116 and HEK293T screens), (2) cells
> containing the 200-gene test perturbation set from the HepG2 and iPSC screens, and all cells with
> a valid guide pair from the Jurkat Resting and Jurkat Active screens.
>
> **5.2 Model and Code Availability.** … we **have made** X-Cell (55M) model artifacts publicly
> available. The trained model weights are hosted on the Hugging Face Model Hub at
> `https://huggingface.com/Xaira-Therapeutics/X-Cell`. Furthermore, we provide software and
> tutorials on GitHub at `https://github.com/xaira-therapeutics/x-cell`. The dataset and trained
> model weights are released under the CC BY-NC-SA 4.0 license.

No accession, no Figshare+, no Zenodo, no GEO, no S3/GCS bucket, no controlled-access route, no
"available on request". **The only named location for either artifact is a Hugging Face repo that
is empty.**

### 1.4 Everywhere else I looked, and what came back

| where | how | result |
|---|---|---|
| HF datasets `X-Atlas`, `Pisces`, `jurkat`, `perturb` (50 rows, by lastModified) | API | only `Xaira-Therapeutics/X-Atlas-{Orion,Pisces}` and `slaf-project/X-Atlas-Orion` (a re-format of Orion). **No Pisces mirror.** |
| HF models `x-cell`, `xcell`, `pisces` | API | only the empty `Xaira-Therapeutics/X-Cell`. **No re-upload.** |
| `Xaira-Therapeutics/X-Atlas-Orion` | API | `usedStorage` 128,446,247,869 B, 337 files, **lastModified 2025-09-19** — unchanged, HCT116+HEK293T only. We already hold it. |
| GitHub org `Xaira-Therapeutics` | `gh api orgs/…/repos?type=all` | **one repo**, the stub. |
| Forks (6) | `gh api …/forks` | all at the stub commit. |
| PyPI `xcell` / `x-cell` | API | `xcell` is jason-weirather's unrelated xCell package; `x-cell` is 404. **Xaira never published.** |
| GEO | E-utilities, `Xaira[All Fields]` | **0 hits.** |
| Zenodo | API (reachable today; 166 was WAF-blocked) | `Xaira` → 1 hit, an unrelated SARS-CoV-2 record. |
| CZ CELLxGENE | curation API, 388 public collections | **0** matching xaira/pisces. |
| bioRxiv supplementary material | page fetch | **no media files** — only the PDF. |
| Figshare+ (Orion's processed-data home) | paper footnotes | the paper's Figshare+ links are DepMap 24Q4 and Replogle 2022, i.e. **other people's data**. Nothing Xaira. |

### 1.5 Downloads (48 MB total)

| file | bytes | sha256 (first 16) |
|---|---:|---|
| `data/external/xcell/xcell_preprint.pdf` (111 pp) | 39,557,181 | `fbad266edbb644fc` |
| `data/external/xcell/xcell_preprint.txt` (pdftotext -layout) | 360,472 | `4d59a3047ebe5067` |
| `data/external/xcell/xcell_biorxiv_full.html` | 700,098 | `93751681a674f37c` |
| `data/external/xcell/xcell_biorxiv_full.txt` | 188,677 | `d90102a1a65a3972` |
| `data/external/xcell/hgnc_protein_coding.txt` (19,274 genes) | 9,485,022 | `33741714037a0b17` |
| `data/external/xcell/pisces_LICENSE.md` | 20,848 | `9e502ba75fd0fd30` |
| `data/external/xcell/pisces_README.md` | 956 | `9a22d02eb139312a` |
| `data/external/xcell/xcell_README.md` | 3,629 | `cb4cb69b5dd01fa6` |
| `data/external/xcell/xcell-repo/*` (13 files: README, MODEL_CARD, LICENSE, docs, `model.py`) | 28,705 | see `work/` |

---

## 2. Gate 2 — LICENCE. Stated, not resolved. **This is Ahmet's call.**

Both artifacts carry **CC BY-NC-SA 4.0**. The operative clauses, from
`data/external/xcell/pisces_LICENSE.md` (the full legal code, sha256 `9e502ba7…`):

> **§1(k) NonCommercial** means *"not primarily intended for or directed towards commercial
> advantage or **monetary compensation**."*
>
> **§2(a)(1)** grants the right to *"reproduce and Share the Licensed Material … **for
> NonCommercial purposes only**"* and to *"produce, reproduce, and Share **Adapted Material** for
> NonCommercial purposes only."*
>
> **§3(b) ShareAlike** — *"if You Share Adapted Material You produce … The Adapter's License You
> apply must be a Creative Commons license with the same License Elements, this version or later,
> or a BY-NC-SA Compatible License."*

**The open question, stated exactly as 182 stated its own and equally unresolved:**

> Is submitting a prediction derived from a CC-BY-NC-SA-4.0 model or dataset into a competition
> awarding **$100,000 / $50,000 / $25,000** a use "primarily intended for or directed towards …
> monetary compensation"? And if it is permitted, does **ShareAlike** then attach to what we
> produce?

**Three things a decision would have to weigh, and I am deliberately not weighing them:**

1. **CC's NC test is *purpose*-based and softer than Arc's State licence.** State's §1.4 reads
   *"not undertaken for direct or indirect monetary compensation"* — a use test. CC-BY-NC-SA's
   §1(k) reads *"not **primarily** intended for or directed towards commercial advantage or
   monetary compensation."* Those are different bars and 182's flag does not transfer verbatim.
2. **ShareAlike is a second, independent condition State does not have**, and it bites on
   *Adapted Material* that is *Shared*. Whether a submitted `.h5ad` of predicted counts is Adapted
   Material, and whether submitting it to a scorer is Sharing, are both genuinely arguable.
3. **The asymmetry that makes this work worth doing anyway.** Using an artifact to *measure a
   ceiling* is not the same act as *shipping its output*. Everything in §3–§4 below is measurement,
   is squarely inside any reading of NC, and would remain useful even if shipping were ruled out.

Note also the **preprint text itself** carries a different licence from the artifacts: every page
of the PDF reads *"It is made available under a **CC-BY-NC 4.0** International license"* — no
ShareAlike. So the paper is quotable on softer terms than the (unavailable) data and weights.

**No decision is recorded here and nothing was argued away. Nothing was used under this licence
beyond reading the paper and the two READMEs.**

---

## 3. Gate 3 — the structural test that killed State. **X-Cell passes it.**

182's decisive leg needed no accuracy measurement: every published Replogle State checkpoint sets
`pert_rep: onehot` over a frozen **2,024-perturbation** vocabulary, the intersection with our panel
is **0/300**, and a one-hot has no out-of-vocabulary path. I asked the same two questions of X-Cell.

**There is no checkpoint to read a config off** — that is gate 1's failure, and it means this
section is answered from the preprint's architecture specification rather than from a `config.yaml`.
That is a real weakening of the evidence and it is stated as such. But the specification is
explicit, and the parts that *can* be measured, I measured.

### 3.1 How the perturbation is represented — verbatim, preprint §A.1.2

> `h^(ℓ) ← h^(ℓ) + CrossAttn(RMSNorm(h^(ℓ)), C_pert)` … *"where the context `C_pert ∈ R^{6×d}`
> consists of six prior knowledge embeddings **associated with the perturbed gene**. These
> embeddings are derived from ESM-2 protein language model, STRING interaction networks, GenePT
> gene representations, DepMap genetic dependency profiles, JUMP-Cell Painting morphological
> features, and the gene identity embedding itself."*

And the vocabulary, preprint §4.2.11:

> *"Feature gene lists were further restricted to a set of **19,082 protein-coding genes**, defined
> by the intersection of **GENCODE annotations and the ESM-2 vocabulary**."*

with §A.1.1 giving `G ≤ 19,408 protein-coding genes` for the shared gene-identity vocabulary. The
public API takes an **HGNC symbol** (`model.predict(adata, perturbation="BRCA1")`), and the
quickstart states *"X-Cell covers all ~19,000 protein-coding genes; genes not in the vocabulary are
zero-imputed."*

**This is a feature-keyed conditioner, not a one-hot over a training vocabulary.** Five of the six
sources (ESM-2, STRING, GenePT, DepMap, JUMP) are computable for a gene the model never saw
perturbed. That is precisely the out-of-vocabulary path State's released artifacts lack.

### 3.2 The measurement (`gate3.py` §1)

| source | vocabulary size | panel covered | missing |
|---|---:|---:|---|
| HGNC protein-coding locus group | 19,274 | **299/300 direct + 1 via alias = 300/300** | **none** |
| STRING (proxy: `slim` v0.3 symbol-keyed embeddings) | 21,688 | **298/300** | `EEF1A2`, `KIF21B` |
| DepMap (on-disk `depmap26q1` build) | 100 | 0/300 | — **uninformative**, see below |
| **contrast — published pretrained State (182)** | **2,024** | **0/300** | all 300 |

- The DepMap row is **not** a finding: our on-disk build is a 100-gene local subset, not DepMap
  24Q4 Public (X-Cell's footnote 3). I did not download 24Q4 because it cannot change the verdict —
  a gene absent from **one** of six priors is degraded conditioning, not an OOV wall, and the
  gene-identity embedding covers it regardless.
- The two STRING misses are a property of *our* proxy file, not of STRING v12.

### 3.3 How the context is conditioned

The mandate's second structural question — *"zero-shot into an unseen line requires conditioning on
the basal state; if X-Cell conditions on a context identifier, it cannot be pointed at D/E/F."*

**It conditions on control cells.** §4.1.4: *"we used **32 fixed control cells** randomly selected
from NTCs in the test set to generate predictions for all test perturbations."* The public API is
`predict(control_cells_h5ad, perturbation=...)`. There is no context token. And the paper's
zero-shot protocol adds **test-time adaptation** — *"self-supervised fine-tuning on **unlabeled
non-targeting control cells** from the target domain using only the MMD loss"* — which is exactly
the material the challenge hands us for D/E/F (18,400 NTC cells per context, no perturbation
labels).

### 3.4 Verdict on gate 3

**PASS on both legs.** If the weights existed, X-Cell could be pointed at all 300 targets in
contexts D/E/F using only the released control cells. **The lead dies at gate 1, not gate 3.** That
is the operationally important difference from 182: State was closed *structurally and permanently*;
X-Cell is closed *on a missing upload*, and the upload could land at any time.

### 3.5 Two cautions about the paper's headline number, from 182's own frame

182's mechanism §3 was that a large published direction-match number *shipped without a null* is
mostly the shared response axis. The same reading applies here, and the paper's Table B7 shows it
from the inside:

| benchmark (all **held-out perturbations in a context the model trained on**) | X-Cell DE Direction Match | X-Cell **Centroid Accuracy** (specificity) |
|---|---:|---:|
| iPSC-200 | 0.9149 | **0.5525** |
| HepG2-200 | 0.9502 | **0.5473** |
| Replogle–Nadig (380 held-out HepG2 perturbations) | 0.8762 | 0.7876 |

**Direction match 0.91–0.95 with centroid accuracy 0.55.** The identity content is thin exactly
where the raw sign number is enormous — the shape 182 warned about, visible in Xaira's own table.
For orientation only (different metrics, different pools, not commensurable): our production
consensus reads discrimination **0.9500** on 182's h1-test pool.

Second: **the headline "5× Pearson Δ" is measured on held-out *perturbations*, not a held-out
*line***. The genuinely held-context rows are the melanocyte-progenitor and primary-CD4 zero-shot
panels — and those are reported in Figure 5C only, with **no numeric table anywhere in the 111-page
PDF**, so they cannot be read off. Credit where due, though: unlike Arc's held-line
`de_direction_match`, X-Cell's zero-shot benchmark **does** carry its own nulls (*"Control Mean and
Pert Mean … data-only predictions constructed using test control cells and log-fold changes from
X-Atlas/Pisces"*) — 182's missing-null criticism does not apply to this paper's design, only to the
fact that the numbers are locked in a figure.

---

## 4. Gate 4 — cannot be run. What was measured instead, and what it decides.

There is no checkpoint, so there is no in-distribution row, no held-line row, no
excess-over-context-mean-null column and no degeneracy check to compute. 182's harness
(`prep.py`/`infer.py`/`score.py`) is reusable the moment weights land and nothing here supersedes
it.

In its place I priced the decision that actually turns on this experiment: **if the promised Pisces
upload lands before Nov 5, is it worth building a donor table from?**

### 4.1 Only the two Jurkat arms can carry the panel — and they would carry ≥ 292/300

The promised upload has three parts, and 166's library-size gate prices each in one line:

| promised part | targets | panel value |
|---|---:|---|
| all Orion HCT116 + HEK293T cells | 18,924 / 18,312 | **zero new** — `data/donors/orion_*` already hold 290 and 294 of 300 |
| HepG2 + iPSC, **restricted to the 200-gene test set** | 200 | **~0–3/300**. The 200 were drawn from genes with "significant phenotype … ≥50% knockdown" across four screens; they are not the panel, and 166 §1 shows there is nothing between a genome-wide library and zero |
| **Jurkat Resting + Jurkat Active, in full** | **10,872 / 10,878** | **the only panel-relevant content** |

**Measured, not assumed.** §4.2.1 states Xaira's own library rule verbatim: *"subsequent screens
used sgRNA libraries containing perturbations for only the genes expressed in each respective cell
line as determined by bulk RNA sequencing. Expressed genes were defined as transcripts that
cumulatively account for **99% of aligned counts**."* Exp 041 identified **context A = JURKAT**, so
that rule is directly applicable to context A's own 18,400 control cells (`jurkat_library.py`):

| context | 99%-rule library | panel in library |
|---|---:|---:|
| **A (JURKAT)** | **9,114** | **292/300 (97.3%)** |
| B (HELA) | 8,825 | 275/300 |
| C (CAL33) | 9,346 | 279/300 |

Calibration: 9,114 against Pisces' measured Jurkat Resting **10,872** = **−16.2%**, and B and C land
in the same ~9k band, so the shortfall is my proxy (challenge 18,533-gene 5′ scRNA-seq readout axis
vs a full-transcriptome bulk library), not Jurkat. The bias is **one-signed** — a bigger real
library can only carry more — so **292/300 is a floor**. Independently, 299/300 panel targets are
expressed at CPM ≥ 5 in context A against **160.1 ± 9.5** for a random 300-gene draw (+14.6σ).

So a Pisces Jurkat drop would be a genuine **fifth** panel-covering line, covering the panel better
than `k562_gwps` (269/300) — and it is the *same line as a scored validation context*, which is a
second, methodological value: a held-line validation surface matched to the challenge's own panel
and platform, which we have never had.

### 4.2 What one more effective line is worth, in our units

159's law `1/cos²(k) = 66.5 + 1920/k` (KA-2 above), 197 §4.1's board↔surface anchors (×1.954 /
×2.041 / ×2.078), 197's exchange rate **+0.00365 overall per +0.01 board cosine**, and the Gaussian
sign law `q_w = ½ + arcsin(cos)/π` (187 §5 / 197 §3).

| `k_eff` | what it is | surface cos | **board cos** (a23 anchor) | [anchor range] | × ours | **`q_w`** | **Δ overall vs `rc_hybrid`** |
|---:|---|---:|---:|---|---:|---:|---:|
| 3.09 | production today | 0.03813 | 0.07782 | [0.0745, 0.0792] | 0.98× | 0.5248 | −0.0005 |
| 4.00 | best equal consensus of the four lines that exist (197 §4.2) | 0.04278 | 0.08731 | [0.0836, 0.0889] | 1.10× | 0.5278 | +0.0029 |
| **5.00** | **+ Pisces Jurkat as one effective line** | 0.04711 | **0.09616** | [0.0921, 0.0979] | **1.21×** | 0.5307 | **+0.0062** |
| 5.30 | + both Jurkat arms, treated as correlated (same background, two states) | 0.04829 | 0.09857 | [0.0944, 0.1004] | 1.24× | 0.5314 | +0.0071 |
| 6.00 | + two fully independent lines (implausible) | 0.05087 | 0.10382 | [0.0994, 0.1057] | 1.31× | 0.5331 | +0.0090 |
| 10.9 | **the board's best direction** (197: 0.13117) | 0.06247 | 0.12751 | [0.1221, 0.1298] | 1.61× | 0.5407 | +0.0176 |

*(Self-consistency: the `k = 3.09` row reproduces our measured 0.079249 to −1.8%, inside 197's own
±2% cross-family systematic.)*

**Against AGENTS.md §8.2's bar — a single direction arm must deliver +31% cosine to clear rule 13:**

| arm | cosine gain | verdict |
|---|---:|---|
| + Pisces Jurkat, one effective line | **+21.3%** | **MISSES** |
| + both Jurkat arms, correlated (k = 5.3) | +24.4% | **MISSES** |
| + two fully independent new lines (k = 6.0) | +31.0% | clears, but no such pair exists |

**So even a complete, on-time Pisces Jurkat delivery does not clear the single-arm bar on its own.**
It is worth **+0.0062 overall against a remaining gap of +0.0337** — about 18% of the way to first
place — and it would need to be composed with something else to be shippable under rule 13.

That is the number this experiment exists to produce, and it changes the answer to "should we be
waiting for Pisces?" from *maybe* to **no, but poll it cheaply**.

---

## 5. Consequences

- **Close X-Cell and X-Atlas/Pisces as usable assets for this competition cycle**, on availability.
  Do not budget work against them. Do not re-run this search.
- **Keep the 10-second poll 166 recommended, and now poll two repos, not one.** Both are one API
  call: `usedStorage` on `datasets/Xaira-Therapeutics/X-Atlas-Pisces` (0 since 2026-03-25) and on
  `models/Xaira-Therapeutics/X-Cell` (272,925 since 2026-03-16). A change in either is the trigger.
  **The trigger is now priced**: data → +0.0062 overall at best (§4.2), needs composition; weights →
  gate 3 already passed, so 182's harness can be pointed at it in hours.
- **Correct 166 §2a on one point.** It recorded Pisces' HepG2/iPSC arms as a 200-gene validation set
  and its Jurkat arms as the only genome-wide content — that is right — but it also implied Jurkat
  adds little because *"we already hold a Jurkat line"*. That is wrong in the way that matters:
  `nadig_jurkat` is a **common-essential** library at **0/300** panel coverage, whereas a Pisces
  Jurkat library would be **≥292/300**. It would be a genuine fifth panel-covering line, not a
  duplicate. It is still not enough on its own (§4.2), but for a different reason than 166 gave.
- **Record the structural asymmetry against 182.** State is closed *permanently* (a one-hot cannot
  be extended). X-Cell is closed *contingently* (an upload). Any future "should we look at a
  pretrained perturbation model?" question should ask 182's two structural questions first — they
  cost one paper-read and they separate these two cases cleanly.
- **Nothing here moves the board.** No metric was computed, no artifact was built, no slot is
  implicated.

---

## 6. What I did NOT measure, and where this is surface-limited

- **No `pds`, `mse`, `fid`, `jac`, `nmae` or `reach`; no `cell-eval2`; no submission.** Every number
  in §4 is a profile-space or arithmetic statistic, so AGENTS.md rule 9 (`n_real` in a denominator)
  is untouched by construction.
- **Gate 3 is answered from a preprint, not from a checkpoint config.** The mandate explicitly asked
  for the config, and there is none. If weights land, re-run the vocabulary check against the real
  artifact before trusting §3.
- **The DepMap prior coverage is unmeasured** (§3.2). It cannot flip the verdict but it is a gap.
  GenePT (Zenodo 10833191) and JUMP-Cell Painting coverage are likewise unmeasured.
- **§4.2's `k_eff` accounting assumes a Pisces Jurkat line behaves like our existing panel-covering
  donors in 159's law.** It might not: Jurkat is *context A itself*, so its effect vectors are
  drawn from the scored distribution rather than transferred into it. That could make it worth
  **more** than one generic line for A/B/C — and worth exactly one generic line for D/E/F, which is
  the round that pays. The table above is the D/E/F-relevant reading and is the conservative one.
- **The two Jurkat arms' independence is a judgement, not a measurement.** I bracketed it at
  k = 5.0–5.3; it cannot be measured without the data.
- **Xaira's zero-shot held-context numbers are figure-only.** If a v2 or the journal version adds a
  numeric table, that is the one thing in the paper worth re-reading.
- **The X-Cell Mini configuration table (`docs/model.md`) lists `Training: Replogle-Nadig`**, while
  the prose says the model is trained on Pisces. Read together with the paper (continual
  pre-training on four Pisces screens, then fine-tuning on Replogle–Nadig), the likeliest reading is
  that the row names the fine-tuning corpus. But it is ambiguous, and if the released Mini turns out
  to be Replogle–Nadig-only, **it carries no Pisces information at all** and is trained on data we
  already hold — which would make it a much weaker prize than the paper implies. Check this first if
  weights land.

---

## 7. Files

| file | what it does |
|---|---|
| `gate3.py` | KA-1/KA-2 gates (aborts on failure); the vocabulary/OOV measurement (§3.2); context-A basal-expression ceiling with its random-draw negative control; the `k_eff` → cosine → `q_w` → overall projection (§4.2) → `work/gate3.json` |
| `jurkat_library.py` | Xaira's own 99%-of-aligned-counts library rule reproduced on contexts A/B/C, with the Pisces Table-1 calibration row → `work/jurkat_library.json` |
| `work/gate3.json`, `work/jurkat_library.json` | every number above |
| `data/external/xcell/` | the 48 MB inventory of §1.5 — preprint PDF + text, bioRxiv HTML + text, both READMEs, the CC-BY-NC-SA legal code, the GitHub stub's 13 files, the HGNC protein-coding list |

```bash
cd ~/Code/virtual-cell-challenge-2026
uv run python experiments/205-xcell-pisces/gate3.py           # ~40 s, peak RSS ~2.6 GB
uv run python experiments/205-xcell-pisces/jurkat_library.py  # ~60 s, peak RSS ~2.6 GB
```

**Inputs read, never written:** `data/challenge/{pert_counts.csv,gene_names.csv,context_{A,B,C}.h5ad}`,
`data/donors/{orion_hek293t,orion_hct116,k562_gwps,h1_2025,nadig_jurkat,rpe1,k562_essential}/perts.txt`,
`data/external/kolf2/perts.txt`, `data/external/slim/gene_string_embeddings.v0.3.h5`,
`data/external/depmap26q1/extracted/processed/crispr_table.csv`. No production code was touched;
no file under `data/donors/` was created or modified.

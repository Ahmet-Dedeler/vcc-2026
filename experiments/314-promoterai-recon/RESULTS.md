# Experiment 314 — PromoterAI reconnaissance

**Trigger.** Team "ViRtUaL cElL builders" took first place on the Arc Virtual Cell Challenge 2026
leaderboard at **0.240185** with public model name **`PromoterAI_precomputed`**. No entry above
0.220 had previously disclosed any method, so this is the first method signal from the frontier.

**Scope.** Public web sources only. No logins, no contact with anyone, no large downloads.
**All access dates: 2026-09-07** (times UTC, ~19:16–19:40Z) unless otherwise stated.

**Method note.** Where possible I read primary artefacts (repo tree, source files, LICENSE text,
PyPI metadata, maintainer replies in the issue tracker) rather than marketing pages. Source
files quoted below were retrieved from `raw.githubusercontent.com` at the `master` branch head.

---

## 1. What is PromoterAI?

**Confirmed: yes, it is Illumina's, and yes, it is a promoter variant-effect model.** The prior
in the task brief was correct.

| Field | Value | Source |
|---|---|---|
| Origin | Illumina Artificial Intelligence Laboratory (San Diego) | press release, repo owner `Illumina` |
| Paper | Jaganathan K, Ersaro N, Novakovsky G, … Kundaje A, Montgomery SB, Ross MT, Farh KK. *Predicting expression-altering promoter mutations with deep learning.* **Science 389(6760)**, 2025 | Europe PMC |
| DOI / PMID | `10.1126/science.ads7373` / `40440429` | Europe PMC |
| Dates | Illumina announcement **2025-05-29**; Europe PMC first publication date **2025-08-07**; print issue Aug 2025. Repo created **2025-05-22** | press release, Europe PMC, GitHub API |

Notable co-authors beyond Illumina: Anshul Kundaje (Stanford), Stephen Montgomery (Stanford),
Daniel Smedley & Mark Caulfield (Genomics England), Heidi Rehm & Anne O'Donnell-Luria (Broad),
Stephan Sanders. This is a large consortium paper, not a vendor whitepaper.

**Abstract (verbatim, Europe PMC):**
> "Only a minority of patients with rare genetic diseases are presently diagnosed by exome
> sequencing… we describe PromoterAI, a deep neural network that accurately identifies noncoding
> promoter variants that dysregulate gene expression… Our estimates suggest that promoter
> variation accounts for 6% of the genetic burden associated with rare diseases."

### Input
- **DNA sequence only.** One-hot `(None, 4)` input (`promoterai/architecture.py`).
- **20,480 bp at inference** centred on a TSS (`--input_length 20480`; README and maintainer
  confirmation in issue #16). 32,768 bp windows during preprocessing/training.
- **There is no cell-type, assay, or context input.** The model conditions on nothing but sequence.

### Architecture
Dilated depthwise-convolution "metaformer" stack — 24 residual blocks, `model_dim` 1024,
kernel 5, exponentially increasing dilation, multi-scale shortcut output heads
(`promoterai/architecture.py`, `promoterai/train.py`). No attention, no pooling layers
(the maintainer cites the absence of pooling as why he expects indel calibration to hold, issue #17).

### Training targets — this matters for Q4
`data/bigwig/hg38.tsv` lists **498 human bigWig prediction targets** (ENCODE + FANTOM5):

| assay | n tracks |
|---|---|
| CAGE-seq | 188 |
| TF ChIP-seq | 161 |
| DNase-seq | 89 |
| Histone ChIP-seq | 53 |
| ATAC-seq | 7 |

135 unique target labels (TBP, NFYA, YY1, POLR2AphosphoS2/S5, BRD4, EZH2, H3K9me3, …).
Multi-species co-training on mm10 is supported (`data/bigwig/mm10.tsv`, `data/annotation/tss_mm10.tsv`).
Then **fine-tuned on GTEx expression outliers** — `data/annotation/finetune_gtex.tsv` (21 MB;
columns `chrom pos ref alt gene gene_id transcript_id strand tss_dist in_cds spliceai n p_under p_over z`).

### Output — read this carefully, it is the crux of Q4
From `promoterai/score.py`:
```python
twin_model = twin_wrap(model)
df_var['score'] = np.tanh(twin_model.predict(gen_var).round(4))
```
and from `promoterai/architecture.py`:
```python
output_ = tkl.Subtract()([output_alt, output_ref])
output_ = tkl.Lambda(lambda x: tk.backend.mean(x, axis=(1, 2)))(output_)
```

So the released score is literally:

> **score = tanh( mean over {output positions × all 498 tracks} of ( f(alt_seq) − f(ref_seq) ) )**

Three consequences, all load-bearing:
1. It is a **difference between two sequences**. There is no "score this promoter" mode in the
   released scorer — only "score this variant".
2. The mean over `axis=(1,2)` **collapses the track axis**. All 498 tracks, spanning many ENCODE
   cell types, are averaged into a single scalar *before* the tanh. The published score is
   **cell-type-blind by construction**, not by accident.
3. Range (−1, 1). Negative = predicted under-expression, positive = over-expression.
   README recommends thresholds ±0.1, ±0.2, ±0.5. Illumina Connected Annotations documents
   |score| < 0.05 as negligible.

Sources: <https://github.com/Illumina/PromoterAI>,
<https://raw.githubusercontent.com/Illumina/PromoterAI/master/promoterai/score.py>,
<https://raw.githubusercontent.com/Illumina/PromoterAI/master/promoterai/architecture.py>,
<https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1126/science.ads7373%22>,
<https://www.illumina.com/company/news-center/press-releases/2025/b9e6a150-673f-4be5-856c-6a0b246131b4.html>,
<https://help.connected.illumina.com/annotation/data-sources/promoterai>. All accessed 2026-09-07.

---

## 2. Is it publicly available?

**Summary verdict: CODE = RELEASED (restrictive licence). WEIGHTS = RELEASED BUT GATED behind a
signed licence agreement + human-in-the-loop email. PRECOMPUTED SCORES = RELEASED BUT GATED, same
gate. Nothing is anonymously downloadable.**

### 2a. Source code — RELEASED

`https://github.com/Illumina/PromoterAI` (accessed 2026-09-07 via GitHub API)

| Property | Observed value |
|---|---|
| Default branch | `master` (note: **not** `main` — raw URLs on `main` 404) |
| Created / last push | 2025-05-22 / 2026-04-03 |
| Stars | 100 |
| Repo size | 9,825 KB |
| **Releases** | **none** |
| **Tags** | **none** |
| Licence (GitHub API) | `NOASSERTION` / "Other" |

**Complete repo tree (verified — this is the whole thing):**
```
LICENSE                                       3,593 B
README.md                                     5,289 B
data/annotation/finetune_gtex.tsv        21,005,734 B
data/annotation/tss_hg38.tsv             18,226,118 B
data/annotation/tss_hg38_cage.tsv         1,140,507 B
data/annotation/tss_hg38_protein_coding.tsv 6,329,227 B
data/annotation/tss_mm10.tsv              6,449,454 B
data/benchmark/{CAGI5_saturation, GEL_RNA, GTEx_eQTL, GTEx_outlier,
                MPRA_eQTL, MPRA_saturation, UKBB_proteome}.tsv   (22–206 KB each)
data/bigwig/{hg38,mm10}.tsv               66,556 / 70,185 B
promoterai/{__init__, architecture, finetune, generator, preprocess,
            score, tfrecords, train, utils}.py                   (0–5,182 B each)
setup.py                                        786 B
```
**There are no weights anywhere in this tree.** Every blob is accounted for above.

### 2b. Licence — PolyForm Strict License 1.0.0 (READ THIS)

`LICENSE` is verbatim PolyForm Strict 1.0.0. Two clauses matter:

> "The licensor grants you a copyright license for the software to do everything you might do with
> the software that would otherwise infringe the licensor's copyright in it for any permitted
> purpose, **other than distributing the software or making changes or new works based on the software.**"

> "Any noncommercial purpose is a permitted purpose."

So: noncommercial use is fine; **redistribution and derivative works are not licensed**.
Fine-tuning the model, wrapping the code into a competition submission that is then published, or
shipping the scores inside a pipeline artefact are all in tension with "no changes or new works"
and "no distributing". PyPI metadata confirms `license: PolyForm Strict License 1.0.0`.
Commercial use: `AI_licensing@illumina.com`.

### 2c. PyPI package — code only, no weights

`https://pypi.org/pypi/promoterai/json` (accessed 2026-09-07)
- Latest version **1.0rc6**, uploaded **2025-06-22**. Still a release candidate 14 months later.
- **Wheel is 13,351 bytes.** Sdist 13,183 bytes. This is unambiguously code-only — a 13 KB
  artefact cannot contain a 24-block × 1024-dim network.
- Deps: `tensorflow>=2.13,<2.16`, `keras>=2.13,<3.0`, `numpy>=1.20`, `pandas>=2.0`,
  `pyfaidx>=0.8`, `pybigwig>=0.3`. **TF2 / Keras 2 only — no PyTorch** (issue #14 asked; closed).

### 2d. Model weights — GATED, but genuinely reaching outside users

Timeline reconstructed from the maintainer's own replies (`kishorejaganathan` = corresponding author
Kishore Jaganathan):

- **2025-06-07, issue #3:** > "At the moment, we have released the source code and precomputed
  scores, **but not the model weights.** Please check back later for updates."
- **2025-09-04, issue #10** ("Pretrained model weights release"): > "Please complete the license
  agreement at the top of the README page; the download link will be shared via email shortly
  after submission." — i.e. the gate now covers weights too.
- **2026-06-05/08, issue #16:** a third-party user is running
  `--model_folder promoterAI_v1_hg38/model/promoterAI_v1_hg38_finetune` and
  `.../promoterAI_v1_hg38_mm10_finetune`. **Weights are demonstrably in outside hands in 2026.**
- **2026-05-11, issue #15:** another user is running `promoterAI_v1_hg38_finetune` in production.

**Status: RELEASED BUT GATED.** Not "coming soon" — real users have them. But the only route is
the Adobe Sign licence agreement, after which a human emails a link. Latency is unknown and
human-dependent (issue #2, 2025-06-02: the maintainer replied and sent the link within ~2 hours;
issue #10 took **7 days** for any reply at all).

Gate URL (from README, verbatim):
`https://illumina2.na1.adobesign.com/public/esignWidget?wid=CBFCIBAA3AAABLblqZhAuRnD5FtTNwyNo-5X6njTJqQOOMu3V_0nU0MjxSi_9PLCrquWaKSRrT3e1RhHkr7w*`

### 2e. Precomputed scores — GATED, same gate

README, verbatim:
> "PromoterAI precomputed scores for all human promoter single nucleotide variants are freely
> available for academic and non-commercial research use. Please complete the license agreement;
> the download link will be shared via email shortly after submission."

- Filename is **`promoterAI_tss500.tsv.gz`** (issue #16, quoted by the user in a working `zcat`
  command — this is a real file, not a guess).
- Confirmed usable as-is: issue #2, maintainer: > "Yes, the scores can be used as is."
- **File size: not publicly observable.** I could not see it without signing.
- **No public mirror found.** Searched Zenodo/figshare/general web for `promoterAI_tss500` and
  `promoterAI_v1_hg38`; no mirror surfaced. Treat as not independently redistributable
  (consistent with PolyForm Strict).

### 2f. The one anonymous route — and why it does not scale

**SpliceAI/PromoterAI Lookup**, Broad Institute: <https://promoterailookup.broadinstitute.org/>
(accessed 2026-09-07)
- Serves "PromoterAI scores for promoters of all protein-coding transcripts in Gencode v39"
  (added 2025-06-27). hg19 and hg38.
- **Per-variant lookup only. No bulk download.** The site states it no longer uses Illumina's
  precomputed SpliceAI tables.
- Explicit rate limit: > "no more than a handful of queries per-user per-minute."
- Source code is available to self-host an API server — but self-hosting still needs the weights.

At a handful of queries/minute, harvesting even one gene's ±500 bp window (~3,000 SNVs) is a
multi-hour abuse of a shared academic service. Harvesting 300 target genes is out of the question.
**This is not a viable acquisition path and attempting it would be rude.**

### 2g. Commercial route
PromoterAI ships as a data source inside **Illumina Connected Annotations** (v3.27 documented),
GENCODE 39. Commercial licensing via `AI_licensing@illumina.com`. Price not public.

---

## 3. Coverage and semantics of the precomputed scores

| Question | Answer | Evidence |
|---|---|---|
| **Window** | TSS **±500 bp** (inferred from filename `promoterAI_tss500.tsv.gz`; not independently confirmed) | issue #16 |
| **Variant types** | **SNVs only. No indels.** | issue #17 (open, 2026-06-17). Maintainer: *"In the entire paper, we used only SNPs for **all** analyses."* Indels are scorable by the software but are **not** in the precomputed table. |
| **Genome build** | **hg38.** Models are named `promoterAI_v1_hg38*`; maintainer directs users to the UCSC hg38 FASTA. (Broad's lookup also offers hg19 — that is their own re-run, not Illumina's table.) | issues #15, #16 |
| **Annotation** | **GENCODE v39** | Illumina Connected docs; Broad lookup |
| **Granularity** | **Per variant × per transcript.** Not per-gene, not per-position. | Illumina Connected Annotations JSON output fields: `strand`, `distanceFromTss`, `geneId`, `transcriptId`, `score` |
| **Minimum columns** | `chrom`, `pos`, `ref`, `alt`, `strand` (a user's `zcat … awk` in issue #16 pulls exactly these by name) | issue #16 |
| **Score semantics** | tanh-squashed mean track-delta, (−1,1). Negative = under-expression. |±0.1/±0.2/±0.5 thresholds| README; Connected docs |

**Gene/transcript universe (measured directly from the repo's own annotation file):**
`data/annotation/tss_hg38_protein_coding.tsv` — **87,151 transcript TSS rows**,
**19,672 unique gene symbols**, 19,698 unique Ensembl gene IDs, 77,340 unique genomic positions,
25 chromosomes. So essentially the full protein-coding genome is covered. The 300 VCC target
genes will almost certainly all be present.

**Size estimate (mine, not observed):** 87,151 TSS × 1,001 bp × 3 alt alleles ≈ **262 M rows**
upper bound; deduplicating to 77,340 unique positions gives ≈ **232 M rows**. A gzipped TSV of
that shape is plausibly **2–6 GB**. Flagging clearly: **this is arithmetic, not an observation.**

**Reproducibility trap worth recording (issue #16, 2026-06-08).** The precomputed scores are the
**mean of two models** — `promoterAI_v1_hg38_finetune` and `promoterAI_v1_hg38_mm10_finetune` —
at `--input_length 20480` on the UCSC hg38 FASTA. And they are **TensorFlow-version-sensitive**:
a user got visibly wrong values on TF 2.15 and near-exact agreement only after downgrading to
**TF 2.13**. If we ever recompute anything, TF 2.13 is mandatory.

---

## 4. THE KEY QUESTION — could a promoter model help predict a CRISPRi knockdown transcriptome?

### 4.1 What the artefact literally is
For each (transcript TSS, position within ±500 bp, alt allele): one scalar Δ. That is all.
There is **no per-gene number in the file**, and **no promoter-activity number in the file** —
only derivatives with respect to point mutations.

### 4.2 The only plausible route from that file to a VCC feature
Aggregate the ~3,000 SNV scores inside each gene's promoter window into per-gene statistics:
`mean|score|`, `max|score|`, `sd(score)`, `frac(score < −0.1)`, mean of the negative tail,
a positional profile binned by `distanceFromTss`. That yields a **per-gene scalar (or short
vector) for ~19.7k genes**, joinable to the 300 VCC targets by symbol or Ensembl ID.

Since the leaderboard name says **"precomputed"**, this is almost certainly what the team did:
a static table → a groupby → a per-gene feature. **No GPU, no weights, no inference.**

**Semantic interpretation of that aggregate:** *promoter mutational sensitivity* — how tightly the
model believes this gene's regulatory output is coupled to its own promoter sequence.

### 4.3 Why that is a priori relevant to CRISPRi — the honest steelman
CRISPRi (dCas9–KRAB) silences a gene by depositing H3K9me3 and blocking Pol II **at the TSS**.
The literature on CRISPRi efficacy is unambiguous about what drives it per target:
- **Distance of the guide to the true CAGE-defined TSS** — "the FANTOM5/CAGE promoter atlas
  represents the most reliable source of TSS annotations, and proximity to the FANTOM5/CAGE-defined
  TSS predicts sgRNA functionality on a genome-wide scale."
- **Chromatin accessibility at the target site** — efficacy concentrates in open chromatin.
- **Nucleosome positioning and local sequence.**
(Radzisheuskaya et al., *NAR* 2016, <https://academic.oup.com/nar/article/44/18/e141/2468344>;
Replogle et al., *eLife* 2022, <https://elifesciences.org/articles/81856>. Accessed 2026-09-07.)

**PromoterAI was trained on exactly that feature set** — 188 CAGE, 89 DNase, 7 ATAC, 53 histone
ChIP, 161 TF ChIP tracks. Its internal representation *is* a model of promoter architecture. So the
steelman is real and not silly: a gene whose promoter the model reads as sharply CAGE-peaked,
TF-dense and accessible is plausibly a gene CRISPRi knocks down hard and cleanly.

**Where it could enter a VCC model:**
- Modulate the predicted magnitude of the **target gene's own** log-fold-change — the single most
  predictable coordinate in the whole prediction vector.
- As a **shrinkage weight**: for low-sensitivity targets, pull the predicted response toward the
  non-targeting-control mean.
- As a per-gene covariate in a calibration of response amplitude.

### 4.4 Why I think it is thin — the skeptical half, which is the important half

**(a) It is cell-type-blind by construction.** `mean(x, axis=(1,2))` in `twin_wrap` averages over
all 498 tracks before the tanh. K562, Jurkat, and a held-out line all get the same number. The VCC
task is **zero-shot generalisation to held-out cell lines** — that is the entire difficulty. A
per-gene constant contributes **exactly zero** to cross-cell-line transfer. It can only sharpen the
global, cell-line-independent ranking of target genes.

**(b) It is a Δ model, not a level model.** If what you want is "how active is this promoter in
this cell line", the precomputed table is the wrong object. The right object — the 498-track
forward pass on reference sequence, which *is* cell-type-resolved — requires the **gated weights**,
and that is emphatically not what "precomputed" means.

**(c) It is confounded with features that are free.** A per-gene `mean|PromoterAI|` will correlate
heavily with promoter CpG-island status, GC content, TF-ChIP peak density, and DNase signal — all
obtainable from ENCODE/UCSC with no licence, no form, and no email. **This is the cheap
falsification: before acquiring anything, test whether a free CpG/GC/DNase promoter feature
reproduces the effect.** If it does, PromoterAI is a expensive proxy for a free feature.

**(d) It may be redundant with our own training data.** VCC training data contains *observed*
knockdowns. For any target gene perturbed anywhere in training — even in a different cell line —
the empirical knockdown depth is directly measurable and is a **strictly better** estimate of
CRISPRi efficacy than a sequence model's guess. A promoter model only adds value for targets with
no training observation at all. That set may be small or empty.

**(e) The guide-scoring route is unavailable.** PromoterAI scores ref/alt SNVs, not 20-mer guide
binding. It does not do guide efficacy natively. And if the challenge specifies **target genes**
rather than guide sequences, there is no guide to score in the first place.

**(f) What it can never supply: the other ~17,999 coordinates.** The task is to predict a whole
transcriptome after knocking down one gene. PromoterAI has no gene–gene network, no co-expression,
no pathway structure. It speaks about the target's own promoter and is silent about every other
gene in the prediction vector. Against a metric that weights differential expression across all
genes, a per-target scalar is a **very thin** feature.

**(g) Name-reading caution.** `PromoterAI_precomputed` is a self-reported public label on a
leaderboard. It is not a method description. It may be one feature among fifty, it may be
mislabelled, it may be aspirational, or it may be deliberate misdirection. The observable is
0.240185; the attribution to PromoterAI is **the team's claim, not evidence**.

### 4.5 Net read
The mechanism is **coherent but narrow**: a cell-line-invariant, per-target-gene prior on
"how disruptable is this promoter", which at best sharpens the target gene's own predicted
knockdown depth and the amplitude calibration around it. It supplies **nothing** for the
held-out-cell-line transfer problem and **nothing** for the transcriptome-wide response — which
together are where the score actually lives. If a +0.02 leaderboard jump really came from this,
the more likely story is that a per-gene promoter-architecture prior improved response-amplitude
calibration, and **that story is testable with free ENCODE features first**.

---

## 5. Related public promoter/regulatory models — the comparison that reframes everything

All entries verified **2026-09-07**. "Downloadable" = a link that resolves without an application/review.

| Model | Weights downloadable now? | Licence | Output | Cell-type resolved? |
|---|---|---|---|---|
| **PromoterAI** | **No — signed licence agreement + emailed link** | **PolyForm Strict 1.0.0** (noncommercial; **no derivatives, no redistribution**) | one scalar Δ per variant | **No — track axis averaged away** |
| **Enformer** | **Yes.** TF-Hub `tfhub.dev/deepmind/enformer/1` now 302-redirects to <https://www.kaggle.com/models/deepmind/enformer/tensorFlow2/enformer/1>; code at <https://github.com/google-deepmind/deepmind-research/tree/master/enformer> | Apache-2.0 code, predictions CC-BY 4.0 | 5,313 human + 1,643 mouse tracks @128 bp, ~200 kb context | **Yes** |
| **Borzoi** | **Yes.** Direct GCS `.h5` links in <https://github.com/calico/borzoi> (4 folds × human/mouse); PyTorch ports <https://huggingface.co/johahi> | **Apache-2.0** (LICENSE fetched and confirmed); johahi port MIT | RNA-seq coverage + DNase/ATAC/ChIP @32 bp, many tissues | **Yes** |
| **Sei** | **Yes.** <https://zenodo.org/records/4906996> (`sei_model.tar.gz`, 3.3 GB, page live) | **Inconsistent** — GitHub says academic/research only (Princeton OTL for commercial); the Zenodo archive is tagged CC-BY 4.0 | 21,907 chromatin profiles across >1,300 cell lines/tissues + 40 sequence classes | **Yes** |
| **ChromBPNet** | **Yes.** Code <https://github.com/kundajelab/chrombpnet>; bias models <https://zenodo.org/records/7443683>; **1,512 ENCODE4 models across 408 biosamples** on the ENCODE portal + HF `kundajelab` org | **MIT** code (LICENSE fetched); ENCODE models under ENCODE data-use policy | base-resolution, bias-corrected ATAC/DNase **accessibility**, one model per biosample | **Yes — per biosample** |
| **AlphaGenome** | **Yes, as of 2026-01-28** — code <https://github.com/google-deepmind/alphagenome_research>, weights <https://huggingface.co/google/alphagenome-all-folds>. **HF-gated**: login + click-accept "Model Terms" (no review process). Local inference wants an H100 | client Apache-2.0; docs CC-BY 4.0; **weights = separate non-commercial Model Terms** | expression (RNA/CAGE), splicing, accessibility, histone ChIP, contact maps; 1 Mb; UBERON-keyed | **Yes** |
| **Evo 2** (Arc) | **Yes, fully open.** <https://huggingface.co/arcinstitute/evo2_40b> (+7b/20b/7b_262k), <https://github.com/ArcInstitute/evo2> | **Apache-2.0**, no carve-out | DNA/RNA/protein LM — likelihoods, embeddings, generation. **No tracks.** | n/a |
| **scooby** (Gagneur, Borzoi-derived) | **Yes.** <https://github.com/gagneurlab/scooby>; checkpoints <https://huggingface.co/johahi/neurips-scooby>, <https://huggingface.co/lauradmartens/onek1k-scooby>, <https://huggingface.co/lauradmartens/epicardioids-scooby>. Nat Methods 2025, <https://www.nature.com/articles/s41592-025-02854-5> | **MIT** | **single-cell-resolution** scRNA-seq coverage + scATAC, Borzoi backbone + LoRA, conditioned on a cell embedding | **Yes — per cell**, but only cell types in its training sets (PBMC, OneK1K, epicardioids) |

**Cross-cutting caveat that applies to all eight, PromoterAI included:** none of these performs
perturbation-response prediction in the VCC sense. They are all DNA-sequence-in / track-out (or
likelihood-out). The only way any of them touches "perturbation" is by substituting the input
sequence and reading the predicted difference — that is variant-effect scoring, not knockdown
modelling.

### The reframing this produces
PromoterAI is, for our purpose, **the weakest and most encumbered member of its own family**:
- It is the **only** one behind a signature-and-email gate.
- It has the **most restrictive licence** of the eight (PolyForm Strict forbids derivatives and
  redistribution outright; Borzoi is Apache-2.0, ChromBPNet and scooby are MIT).
- It is the **only one whose released output is cell-type-collapsed**. Enformer, Borzoi, Sei,
  ChromBPNet, AlphaGenome and scooby all emit **per-cell-type** tracks.

If §4's steelman is right — that a promoter-architecture prior helps calibrate CRISPRi knockdown
depth — then **ChromBPNet is the strictly better instrument**: MIT-licensed, downloadable right
now, and it emits exactly the quantity the CRISPRi-efficacy literature says matters (base-resolution
promoter accessibility), **per biosample**, for 408 biosamples. Borzoi adds predicted expression on
Apache-2.0. Both give a cell-type-resolved feature; PromoterAI's precomputed table structurally
cannot.

---

## Summary — what is worth acquiring, and what it would cost

1. **PromoterAI is real, strong science, and effectively closed.** Illumina AI Lab, *Science*
   389(6760) 2025, DOI `10.1126/science.ads7373`. Code is public under **PolyForm Strict 1.0.0**
   (noncommercial, **no derivatives, no redistribution**); the 13 KB PyPI wheel contains no weights;
   the repo tree contains no weights. Both weights (`promoterAI_v1_hg38_finetune`) and the
   precomputed table (`promoterAI_tss500.tsv.gz`, hg38, GENCODE 39, TSS ±500 bp, **SNVs only**,
   my estimate ~230–260 M rows) exist and reach outside users, but **only** by signing an Adobe Sign
   agreement and waiting for a human to email a link — observed latency 2 hours to 7 days.
2. **Cost to acquire: ~$0 and one signature, but 1–7 days of latency, a licence that arguably
   forbids what we would do with it, and a multi-GB download.** The one anonymous route (Broad's
   lookup) is per-variant with an explicit "handful of queries per minute" limit — not viable and
   not polite to abuse.
3. **The mechanism is coherent but thin.** The only route from that file to a feature is
   per-gene aggregation of ±500 bp SNV deltas into a "promoter mutational sensitivity" scalar.
   `twin_wrap` averages over `axis=(1,2)` — **all 498 tracks collapse into one number**, so the
   feature is cell-type-blind *by construction*. It contributes nothing to the held-out-cell-line
   transfer problem and nothing to the ~17,999 non-target coordinates, which is where the score lives.
4. **Test the cheap falsification before acquiring anything.** A per-gene `mean|PromoterAI|` is
   heavily confounded with promoter CpG/GC content and DNase signal — free from ENCODE/UCSC — and
   is plausibly redundant with knockdown depth already *observed* in our own training data for any
   target seen there. If a free promoter-accessibility feature reproduces the effect, PromoterAI is
   an expensive proxy for something we already have.
5. **If a promoter prior does help, acquire ChromBPNet instead, not PromoterAI.** MIT-licensed,
   downloadable today, 1,512 ENCODE4 models over 408 biosamples, base-resolution promoter
   accessibility — **per cell type**, which is exactly what PromoterAI's precomputed scores
   structurally cannot give. Borzoi (Apache-2.0) and scooby (MIT, single-cell) are the next rungs.
   Finally: `PromoterAI_precomputed` is a **self-reported leaderboard label**, not a method
   description — 0.240185 is the observable; the attribution is the team's claim, not evidence.

---

### Source index (all accessed 2026-09-07, ~19:16–19:40 UTC)

**Primary — PromoterAI**
- <https://github.com/Illumina/PromoterAI> (repo, README, issue tracker; queried via GitHub API)
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/README.md>
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/LICENSE> (PolyForm Strict 1.0.0)
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/promoterai/score.py>
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/promoterai/architecture.py>
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/data/annotation/tss_hg38_protein_coding.tsv>
- <https://raw.githubusercontent.com/Illumina/PromoterAI/master/data/bigwig/hg38.tsv>
- Issues #2, #3, #10, #15, #16, #17 — e.g. <https://github.com/Illumina/PromoterAI/issues/16>
- <https://pypi.org/pypi/promoterai/json>

**Paper & announcement**
- <https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1126/science.ads7373%22>
- <https://www.science.org/doi/10.1126/science.ads7373>
- <https://pubmed.ncbi.nlm.nih.gov/40440429/> (cookie-walled on fetch; metadata via Europe PMC)
- <https://www.illumina.com/company/news-center/press-releases/2025/b9e6a150-673f-4be5-856c-6a0b246131b4.html> (2025-05-29)

**Distribution / annotation services**
- <https://help.connected.illumina.com/annotation/data-sources/promoterai>
- <https://help.connected.illumina.com/annotation/v3.27/data-sources/promoterai>
- <https://promoterailookup.broadinstitute.org/>

**CRISPRi efficacy literature (for §4.3)**
- <https://academic.oup.com/nar/article/44/18/e141/2468344> (sgRNA position / CAGE TSS proximity)
- <https://elifesciences.org/articles/81856> (dual-sgRNA libraries, optimal effectors)

**Comparator models** — see the URL column in the §5 table; additionally
<https://github.com/google-deepmind/alphagenome_research>, <https://zenodo.org/records/4906996>,
<https://zenodo.org/records/7443683>, <https://www.nature.com/articles/s41592-025-02854-5>.

# 306 — Reconnaissance, 2026-09-07: nothing new is downloadable, and the frontier is not using data we lack

**Verdict in one line.** Between 2026-08-30 and today **no perturbation dataset that touches the 300-target
panel, the six challenge lines, or the challenge's own library became publicly downloadable** — X-Atlas/Pisces
is still an empty HuggingFace repo (`usedStorage 0`, three files, unchanged since 2026-03-25), GEO/ArrayExpress/
Zenodo/scPerturb show zero qualifying deposits after the freeze, the only genuinely new item is a cortical-neuron
CRISPRi preprint with no verified accession — and on the public board **no entry above 0.20 publishes a method
description at all** (1 of the top 30, and it is a joke), while **every rival who does describe their method names
only sources already on this disk** (Replogle K562-GWPS, X-Atlas/Orion HCT116-HEK293T, GSE314342 primary CD4,
the 2025 VCC H1 context) and scores 0.136–0.175 against our 0.2161.

**Scope and hygiene.** 2026-09-07, 15:44–17:05 UTC. Metadata and web only. No `vcc`, no `cell-eval2`, no build,
nothing under `submissions/`, no model proposed. `uv run python` throughout. Everything fetched lives in
`experiments/306-recon-oct/work/` (**1.1 MB total**, cap 200 MB); nothing was written into `data/` or any sibling
directory. Part C peak RSS **29.9 MB**, runtime 0.03 s. Every source, with the exact sentence relied on, is in
`work/sources306.md` (S1–S32); the citations below are those tags.

Live board at 15:48:55Z (`work/lb_live_306.json`, 769 entries, public unauthenticated API — S27):
leader `cqawesome`/`m0907h` **0.2393600** (submitted 15:29Z today), ours **0.2160944** rank 3, gap **+0.0232656**.

---

## 1. Part A — what became available since 2026-08-30

**State vocabulary, used strictly.** RELEASED = I observed downloadable files and their sizes.
ANNOUNCED = a page or paper says it exists; no files observed. COMING SOON = the host itself says the data is
not yet uploaded. PUBLISHED-NO-DATA = a paper exists; no deposit verified.

| item | state | what it actually contains | panel / line overlap | licence | size | what it would buy | evidence |
|---|---|---|---|---|---|---|---|
| **X-Atlas/Pisces** (HF `Xaira-Therapeutics/X-Atlas-Pisces`) | **COMING SOON** — unchanged since 2026-03-25 | Repo holds exactly `.gitattributes`, `LICENSE.md`, `README.md`. `usedStorage = 0`. README: "**(Coming Soon)** The following data will be uploaded…" listing HCT116+HEK293T (Orion re-release), HepG2+iPSC 200-gene validation subsets, and **Jurkat Resting + Jurkat Active**. Claimed atlas: 25.6 M cells, 16 contexts | Jurkat is a challenge-line identity; HepG2/HCT116/HEK293T are donor lines. **Unmeasurable** — no file | CC BY-NC-SA 4.0 | 0 B present. Scale reference: sibling **X-Atlas/Orion = 128.4 GB / 337 files** for ≈ 8 M cells, so a full Pisces ≈ 400 GB | The one item that would close F5 (more LINES) *and* give a Jurkat perturbation reference. Buys nothing today | S1, S2, S4 |
| **X-Cell** (HF model) | **COMING SOON** (no weights) | 3 files, 272,925 B (`README.md` + one PNG). `downloads = 0`, `lastModified 2026-03-16`. Preprint still **v1, 2026-03-20, never revised, never published**; its own Data Availability says the data "will be uploaded to https://huggingface…" | — | CC BY-NC-SA 4.0 | 0.27 MB | Nothing. A baseline that cannot be run | S3, S6, S7 |
| Xaira GitHub org | no change | one repo, `X-Cell`, `pushed_at 2026-03-16` | — | — | — | nothing | S5 |
| **Arc Institute — any release** | **none** | News index: last VCC post **2026-08-20**; only 2026-09-02 and 2026-09-03 posts, both unrelated genome/3D-genome stories. `cell-eval` last push 2026-07-27 (v0.8.2), `cell-eval2` 2026-08-22, `state` 2026-07-24 (v0.11.1), `cell-load` 2026-07-18. HF `arcinstitute`: 15 datasets, newest perturbation-relevant modification **2026-07-22** | — | — | — | Nothing. Timeline confirmed unchanged: **Oct 22 test release, Nov 5 deadline, 3 validation + 3 test lines, six-metric aggregate** | S10, S11, S13, S14 |
| **VCC 2026 challenge paper** — *Cell*, 2026-08-26, doi 10.1016/j.cell.2026.08.004, PMID 42648290 | **PUBLISHED-NO-DATA** | A four-sentence front-matter piece by the Arc Virtual Cell Initiative Team. Abstract read in full via PubMed efetch; body **unread** (ScienceDirect returned HTTP 403) | states "unseen cell lines"; names none | Elsevier | — | Nothing operational. No metric numbers, no baselines, no accessions in the abstract | S12 |
| **STATE in *Cell*** — 2026-08-31, doi 10.1016/j.cell.2026.07.052 | **PUBLISHED-NO-DATA (post-freeze)** | Journal version of the model already on disk (`data/external/state-pretrained`, `state-vcc`). "We further introduce Cell-Eval…" Not open access | — | Elsevier | — | Nothing new; 106/113/182 already NO-GO'd this model | S15 |
| **GEO** | **no qualifying deposit after 2026-08-30** | Six independent queries. Human Perturb-seq PDAT ≥ 2026-08-01 → 3 series, all ≤ 16 samples, all before the freeze (GSE343369 embryoid-body 225 genes; GSE344269; GSE344535 RPE-1 viral-microprotein ORF). CRISPRi/a + single-cell → 2 series. `perturb*[Title]` → count 11; `CRISPR[Title]`+single-cell → count 2; challenge-line names + CRISPR → count 14, none a transcriptomic screen at scale | ~0 | — | — | Nothing | S16–S18 |
| **ArrayExpress / BioStudies** | **none** | 22 total hits for Perturb-seq/CRISPRi; newest **2026-08-06**, bacterial *whiA*/*whiB2* CRISPRi RNA-seq | 0 | — | — | Nothing | S20 |
| **Zenodo** | **none** | 426 records created 2026-08-01→, newest 20 read: no perturbation screen | 0 | — | — | Nothing | S21 |
| **scPerturb** | **no release since 2024-08** | data record v1.4, created and updated 2024-08-26, 54 files | already inventoried by 270/292 | CC BY | — | Nothing (R3's eight files already spent by 295) | S22 |
| **A screen using the 2025 VCC library** | **not found** | No new H1-hESC CRISPRi deposit in any repository sweep; no preprint describes one. Coverage-limited: spacer sequences are not a searchable index | — | — | — | 293's open question stays open | S23 |
| **ASD cortical-neuron CRISPRi** — bioRxiv 2026-09-02, doi 10.64898/2026.09.02.748925 | **ANNOUNCED (the only post-freeze item)** | "a CRISPRi Perturb-seq screen targeting **1,408 ASD risk genes** in human embryonic stem cell-derived **immature cortical neurons**"; 215 genes with significant global dysregulation | Unknown — **no accession and no availability sentence found** in the served full-text HTML (98 KB, JS-rendered, so "not verified", not "absent") | unknown | unknown | *If* deposited: a genuinely new post-mitotic context with a 1,408-gene library → 268's stage-2 contexts and F5's "more LINES" | S24 |
| **Broad-16 preprint** — bioRxiv 2026-08-25 | already on disk | "a proof-of-concept Perturb-seq dataset targeting 100 genes in 16 diverse cancer cell lines" | 0/300 (measured, 268) | — | 78 MB on disk | Nothing new — this is `data/external/broad16line`, already used by 268/283/289 | S24 |
| **Wang et al. reliability** — bioRxiv 2026-08-12 | already cited | "Among 7,170 perturbations from 29 datasets, 65% are unreliable, 11% shared, 24% specific" | — | — | — | Already in 287's audit and 294's references | S24 |
| **IGC preprint** — bioRxiv 2026-08-24, doi 10.64898/2026.08.21.746243, **not in this repo** | PUBLISHED-NO-DATA | *Virtual-cell models compress unseen intervention geometry through a target-specific generalization bottleneck.* "same-target empirical anchoring transfers intervention identity across contexts far more effectively than increasing exposure to other interventions"; design principle: "empirically anchor intervention identity, then use models to generalize anchored effects" | — | — | — | **Reading only, no data**: an independent, mechanism-level statement that the incumbent's whole construction (transfer the *same target's* measured response) is the right axis, and that adding exposure to *other* interventions is the wrong one. Consistent with 269 (router fails), 263/264 (arbiters fail), 276 (more lines not indicated) | S24 |
| Other Aug preprints — MultiFlow, cancer-line heterogeneity, method-choice, leakage-aware benchmark, HSPC Perturb-seq, Perturb-ME | PUBLISHED-NO-DATA | methods/benchmarks | — | — | — | Nothing | S24 |
| **Any paper reporting a number on this challenge's metric** | **none exists** | Europe PMC: `"Virtual Cell Challenge"` 2026 → **hitCount 3**, all organizer/perspective pieces. `"cell-eval" OR "perturbation discrimination score" OR "pds_cosine"` 2026 → **hitCount 4**, none a 2026 entry | — | — | — | — | S25, S26 |

**Part A conclusion.** Zero items moved from ANNOUNCED/COMING SOON to RELEASED in the eight days since the freeze.
270's §4.3 finding is re-confirmed with harder evidence: `usedStorage` is *still* exactly 0 and the file list is
*still* three metadata files. The only thing that changed in the world of primary data is one preprint (ASD
neurons) whose deposit I could not verify.

---

## 2. Part B — the board's frontier, from public sources only

**B1. The frontier is silent.** The leaderboard API is public and unauthenticated (HTTP 200, 945,747 B, no
credentials — S27). Of **769 entries, 292 write a `description`** and 198 write ≥ 120 characters. But:

| band | entries | with any description |
|---|---:|---:|
| ≥ 0.220 | **2** (`cqawesome` 0.23936, `Aginglab.com` 0.22652) | **0** |
| ≥ 0.200 | 9 | **1** — and it is `Slurm Shadi`: *"Apparently the simulation has a leaderboard."* |
| top 30 | 30 | 1 (the same joke) |

The **best-scoring entry that describes a method at all is rank 42, 0.17543** — and its description is an emission
recipe ("kink mag 0.42, conc 1.0, agree 0.6 … largest-remainder exact column sums; ~+0.0015 level"), not biology.
**There is no public statement of method anywhere at 0.22–0.24.** Treat every `description` as a self-reported
claim; the fields are free text written by competitors and nothing in them is verified.

**B2. Everyone who does describe a method is standing on data we already hold.** The named sources, verbatim (S28):

- rank 92, 0.14008, `Yusif Mukhtarov`: *"Donor-pruned unit-norm fusion of public K562 genome-wide CRISPRi,
  **X-Atlas HCT116**, and **primary CD4 (Rest/Stim8/Stim48)** responses; cross-donor family averaging,
  common-response centering and Gram row-decorrelation; library-preserving coherent raw-count transport at
  footprint-calibrated per-context amplitude."*
- rank 88, 0.14200, `Ideoscopic`: *"per-target response columns from Replogle K562 gwps and from **GSE314342**
  primary CD4+ T cells … then locally decorrelated against the pool of all 296 real columns."*
- rank 109, 0.13085, `liang zuan`: *"Uses public **K562, HCT116, HEK293T and primary CD4** responses plus official
  context controls."*
- rank 78, 0.14954, `Testh`: *"… with the **VCC2025 H1 hESC context** added to the evidence pool (21/300 panel
  targets)."*

That is a near-verbatim description of the incumbent's own construction — whitened multi-source unit-norm fusion
with common-response centring and count transport — reached independently by at least four teams, all scoring
**0.136–0.150 against our 0.2161**. `GSE314342` resolves to *"Genome-scale perturb-seq in primary human CD4+ T
cells … 22 million primary CD4+ T cells from four donors under three stimulation conditions"*, PDAT 2026-02-06 —
**the same study as this repo's GWCD4i donor** (S19). **No describing team names a source that is not already on
this disk.** Two teams also report the same physics we measured: `Life Code` — *"44 percent of our delta energy
sat on genes where even the same line cannot reproduce"*; `tydymy` — *"A reference pseudobulk holds only 155 to
204 cells"*.

**B3. Public code discloses nothing above 0.08.** GitHub search returns 36 repos referencing the challenge; eight
pushed in the last ten days. Verbatim scores (S29): `biona001` publishes a dated submissions table —
**0.0184 (Aug 30), 0.0766 (Aug 31), 0.0649 (Sep 3)**, method "calibrated NB + full shrunken K562 fold changes";
`wuki-creator/magworldlatest` — *"ranked 219 with an overall scaled score of **-0.0448**"* and *"None of these
switches improved the five-fold aggregate consistently"*; `ENNCELADUS` — a negative point estimate on its own
internal metric; `saberhq/sidechain` discloses no score. **No public repository is anywhere near the frontier.**

**B4. Arc has said nothing new about October 22 or the scoring.** The only statements are the 2026-08-20
announcement (S11: *"Thursday, October 22: Final test set released"*, *"final rankings use an aggregate across a
broader panel of six metrics"*, three validation + three test lines) and the four-sentence *Cell* piece (S12).
No FAQ update, no news post, no repository push. The organizer-suggested datasets page
(`virtualcellchallenge.org/app/datasets`) is a client-rendered SPA and returned only "Loading…" — **unread**, so
I make no claim about its contents (S32).

**B5. The forum is out of reach by mandate.** The participant forum is **Discord** (REPORTED), which requires an
account; no Discourse exists. Not accessed. The one public write-up found, *"Arc Virtual Cell Challenge: Has the
Evaluation collapsed?"*, is dated **2025-11-06** and is about the **2025** challenge and its MAE-cap exploit —
already covered by `docs/COMPETITIVE-INTEL-2026-08-28.md` §3.3 (S31).

---

## 3. Part C — the falsifiable check

**The claim.** The most consequential checkable claim from Parts A/B is B2's: that a rival's addition of the 2025
VCC H1 context covers **"21/300 panel targets"** (`Testh`, `v6loco_plus_H1hESC_nb_a2`, 0.14954, read live
15:48:55Z). It is consequential because it is the only public assertion about how much of the *challenge's own
library* anyone can reach — 293's open question — and because it disagrees with this repo's own settled number.

**Known-answer row** (`work/partC_306.py`; 0.03 s; peak RSS 29.9 MB; read-only):

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 2026 panel ∩ 2025 train+val+test targets, from the two files `READ-FIRST` §6 names | **25 / 300** (READ-FIRST §6; `docs/LOG.md:4085`, 2026-09-06) | **25 / 300** | **PASS** |
| KA-1a panel size | 300 | 300 unique `target_gene` | PASS |
| KA-1b 2025 splits | 150 / 50 / 100 | 150 / 50 / 100, union exactly 300 | PASS |

Files and digests: `data/external/vcc2025/{train,validation,test}/pert_counts_*.csv`
sha256[:16] `633d202be221418b` / `9da5dec8c25aff1b` / `c4f0ba439d7c3121`; `data/donors/h1_2025/perts.txt`
`977b808e723c73ee`; `experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv` `d04e2ca81092523e`.
The 25 targets: ACLY, ADNP, AKT2, ANKZF1, BRPF1, HDAC8, HSBP1, MED13, MED15, MED25, MTA1, NFE2L1, PLAGL2, RNF2,
SHPRH, SIN3B, SLIRP, SMARCA5, STAT6, TARBP2, TRAPPC6A, TWF2, ZFP62, ZNF32, ZNF714 (13 in the 2025 train split,
4 validation, 8 test).

**Verdict on the rival's number: not reproducible, and it does not change anything.** I enumerated every
defensible reading of "the 2025 H1 context":

| reading | measured | = 21? |
|---|---:|---|
| 2025 train+val+test union | **25 / 300** | no |
| the `h1_2025` donor bank as actually built (255 perturbations that survived QC) | **20 / 300** | no |
| the 2025 library by construct id (301 genes over 348 guides) | 25 / 300 | no |
| the 2025 library restricted to spacer-resolved constructs | 25 / 300 | no |
| train split only / validation only / test only | 13 / 4 / 8 | no |

21 sits between our built bank (20) and the library ceiling (25) — consistent with a rival applying their own QC
and keeping one target more than we do. **Nobody has more of the 2025 H1 context than we do**, the whole contested
range is 20–25 of 300 targets, and the entry built on it scores **0.14954**. 293's ceiling question is not
answerable from public data because no public screen uses the 2026 library.

---

## 4. Ranked shortlist of acquisitions worth a download

**Zero items clear today.** Nothing found in Part A is both downloadable and useful; the well-supported negative
is the result. Ranked below are the two conditional acquisitions with explicit triggers. I am deliberately not
filling the third slot.

**1 — X-Atlas/Pisces, the Jurkat Resting + Jurkat Active split (and HepG2/iPSC second).**
*Open question it closes:* F5 needs more LINES and 268's stage 2 needs contexts; it is also the only prospect of a
perturbation reference in a challenge-line identity (Jurkat), and 285 already measured Jurkat's construct term
(ρ_guide 0.17) so the arrival would be immediately interpretable. *State:* COMING SOON, `usedStorage 0`, unchanged
for 165 days. *Size:* unknown; sibling Orion is **128.4 GB / 337 files** for ≈ 8 M cells, so budget 30–100 GB for a
Jurkat-only split and ~400 GB for the whole atlas — well beyond a subagent cap, so it would need a streamed,
per-file extraction to pseudobulk like the existing donor banks. *Cost to use:* a donor-bank build (compact delta
matrix), i.e. the `gse295214_h1`/`h1_2025` pattern, plus a PLOLO fold. *Trigger, and the cheap action:* poll
`https://huggingface.co/api/datasets/Xaira-Therapeutics/X-Atlas-Pisces` and act only when `usedStorage > 0` or
`siblings` exceeds 3. **This is a watch, not a download** — one HTTP request, and it is the single highest-value
recurring check between now and October 22.

**2 — The ASD cortical-neuron CRISPRi screen (bioRxiv 2026-09-02, 1,408 genes, hESC-derived neurons), conditional
on an accession appearing.** *Open question it closes:* 268's stage-2 contexts and F5's line axis, with a genuinely
new *post-mitotic* context that nothing on disk provides (Tian 2019's day-7 neurons are 24 perturbations and
restricted). *State:* ANNOUNCED; no accession or availability sentence recoverable from the served full text.
*Size:* unknown. *Cost to use:* panel-overlap count first (ASD risk genes vs the 300 panel is plausibly small —
this must be measured before any download), then a context row for the residual map. *Trigger:* a GEO/ArrayExpress
accession resolving for doi `10.64898/2026.09.02.748925`. **Do not download on the preprint alone.**

**3 — deliberately empty.** Every other candidate is either already on disk (Broad-16, Orion, CMap, GWCD4i/GSE314342,
Nadig, Feng, Jiang, Tahoe, the scPerturb R3 files — all spent or inventoried by 270/278/292/295), assessed and closed
(VIPerturb-seq by 070; GPSAdb by 270; PerturbFate by 270), or has no panel/line overlap at all (GSE343369,
GSE344535, GSE283877, the ArrayExpress bacterial CRISPRi sets). Padding this slot would misrepresent a clean
negative.

---

## 5. What this means for the programme, stated plainly

- The programme's Round-4 conclusion — *"the next investments are new data … or a new modelling idea with a stated
  mechanism"* — is **half-blocked by fact**: on 2026-09-07 there is **no new data to invest in**. The data half of
  that recommendation cannot be executed until Pisces uploads or an accession appears for the ASD screen.
- The gap to `cqawesome` (+0.0233, still widening) is **not** explained by a data advantage that is visible in
  public. Every describable public construction uses our four sources and lands 0.05–0.08 below us. Whatever the
  top two are doing, they are not doing it with data anyone has published.
- One reading worth carrying forward at no cost: the IGC preprint's *"same-target empirical anchoring transfers
  intervention identity across contexts far more effectively than increasing exposure to other interventions"*
  is an independent, mechanism-level endorsement of the incumbent's axis and an independent argument against the
  arbiter/router family this repo has already closed empirically (263, 264, 269, 276). It proposes nothing
  buildable here; it is corroboration, not a lead.
- No model is proposed. Nothing was downloaded into `data/`. No submission is discussed.

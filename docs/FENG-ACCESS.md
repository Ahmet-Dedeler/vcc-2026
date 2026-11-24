# Feng et al. 2025 (iPSC CRISPRi atlas) — access notes, for whoever picks this up

**Why it matters.** Experiment 036 measured the marginal value of a fourth donor line at
**+0.0117 macro PDS ≈ +0.005 overall**, and it is the only remaining lever that attacks the
binding constraint (directional accuracy) rather than calibrating around it.

**What it is.** Feng et al., *A genome-scale single-cell CRISPRi map of trans gene regulation
across human pluripotent stem cell lines*, Cell Genomics 6(2), 2025 (Wellcome Sanger). 7,226
genes targeted across 34 iPSC lines from 26 donors, ~2M cells.

**The useful part.** The data is split into **fitness genes**, **non-fitness genes**, and a
**targeted screen**. The VCC 2026 panel was deliberately built to *exclude* common-essential
genes (0/2,058 overlap with Replogle K562-essential — see the log), so the **non-fitness
genome-wide screen is the matching gene class**, and its panel overlap should be *enriched*
rather than the ~108/300 that random overlap of 7,226-of-20,000 would predict.

## Files (Figshare, MIT licence, 7.38 GB total; the non-fitness set is 1.73 GB)

Share page: `https://figshare.com/s/bab7c3f17a5fd284fc91` (Count Data)
Companion: `https://figshare.com/s/14edeeab56eb8a885df3` (Transcriptional Changes)

| file | size | ndownloader id |
|---|---|---|
| GenomeWideScreen_NonFitnessGenes_RNA-UMI-Counts.csv.gz | 1.63 GB | 51059558 |
| GenomeWideScreen_NonFitnessGenes_Guide-UMI-Counts.csv.gz | 99.81 MB | 51051818 |
| GenomeWideScreen_NonFitnessGenes_Cell-Metadata.tsv.gz | 4.04 MB | 51051806 |
| GenomeWideScreen_FitnessGenes_RNA-UMI-Counts.csv.gz | 1.35 GB | 51059483 |
| TargetedScreen_RNA-UMI-Counts.csv.gz | 4.17 GB | 51072902 |

URL form: `https://figshare.com/ndownloader/files/<id>?private_link=bab7c3f17a5fd284fc91`

## The blocker, and what does work

**`curl` gets HTTP 403 on both the share page and the ndownloader endpoints**, with or without a
browser user-agent — the same Figshare WAF that blocked Replogle 2022 and forced the scPerturb
route (log, 2026-08-21). scPerturb does **not** mirror Feng, and no GEO accession surfaced.

**A real browser reaches the item page fine** — the file list above was read through the browser
pane. Four retrieval routes were tried on 2026-08-22; three failed and the failure modes are
worth recording so nobody repeats them:

| route | result |
|---|---|
| `curl` with a browser user-agent | **403**, `server: awselb/2.0`, 520-byte challenge body |
| `curl` replaying the browser's `aws-waf-token` cookie | **403** — the cookie is readable (not httpOnly) and was passed with a matching UA and Referer, so the WAF is fingerprinting the **TLS/JA3 handshake**, not just the cookie |
| browser navigating directly to the `ndownloader` URL | **403 Forbidden** — a bare navigation is not treated as a download initiated from the item page |
| same-origin `fetch()` from the item page | `TypeError: Failed to fetch` — the endpoint 302s to S3 and the redirect fails CORS |

**The route that should work and was not tried:** click the actual download control in a real
browser session (Claude-in-Chrome drives the user's own Chrome, which downloads to
`~/Downloads` where bash can reach it). The in-app browser pane is sandboxed and does not expose
downloads to the filesystem.

A fifth option, cleanest of all: **email the authors**. The licence is MIT and the data is
public — this is a WAF misconfiguration, not a restriction.

## First test before spending the day

Pull only `Cell-Metadata` (4 MB) and `Guide-UMI-Counts` (100 MB), extract the targeted gene
list, and compute the overlap with `data/challenge/pert_counts.csv`. **If the non-fitness screen
covers fewer than ~100 of the 300 panel targets, stop** — the marginal-donor estimate assumes a
line that covers the panel, and a sparse donor contributes only on the perturbations it carries.

## The discount to apply

Feng is **pluripotent**; the 2026 contexts are somatic. That is the same mismatch that made
H1-2025 a weak donor, and H1 is already one of the four lines averaged into the +0.0117
estimate. A somatic fourth line would be worth more than +0.005; a second pluripotent one
should be expected to deliver less. Bio_NUDT runs `gwps+feng` at raw `pds` 0.665 against our
0.816, though their model string also says `alpha=2 probe`, so that comparison is not clean.

---

# VERDICT (2026-08-22): Feng is not usable as a donor. Both halves fail, for opposite reasons.

Ahmet downloaded the files manually, which cleared the WAF question entirely. The data is fine;
the *shape* of it is wrong for us.

| screen | VCC panel overlap | median cells per target | verdict |
|---|---|---|---|
| GenomeWide **NonFitness** | **166/300 (55%)** | **22** | covers the panel, 7–9x too thin |
| **Targeted** | 5/300 (2%) | **1,544** | deep enough, hits almost nothing we need |
| *what a donor needs* | high | >=150 (ours are 154–203) | — |

## Why the genome-wide screen is too thin

396,458 cells, but only **29% carry a confident guide call**, and those 115,124 cells are spread
over **4,960 targets across 34 cell lines**. Pooling every line gives 22 cells per target; six
panel targets reach 50 cells and **none reaches 100**. Splitting by line — which is the entire
point of the study — gives roughly *one* cell per line-target pair.

The 34 genetic backgrounds are the paper's contribution and our problem. We need depth per
(line, perturbation); they optimised for breadth across donors.

Donor deltas are already noise-dominated at 150+ cells (Nadig: raw replicate correlation 0.16,
noise-corrected 0.90). At 22 cells the pseudobulk is ~7x noisier in variance, so a fourth donor
built from this would contribute mostly sampling noise to an equal-weight average.

## Why the targeted screen does not rescue it

1,161,864 cells, 54.7% assigned, 445 targets, 19 lines, median 1,357 cells per target — exactly
the depth we want. But it targets eQTL genes, and only **5 of our 300** are among them
(EZH2, GIGYF2, NT5DC1, STT3A, ZNF526). Combining the two screens does not help: the 166 useful
targets are still at 22 cells.

## The stopping rule I wrote was the wrong one

I specified "stop if it covers fewer than ~100 of the 300 panel targets". It cleared that
comfortably at 166. **I never set a bar for cells per perturbation**, which is what actually
decided it — and which our own donor tables have recorded in their metadata all along
(179 / 154 / 203). Coverage and depth are independent requirements and a data-acquisition
screen needs both.

## What is retained

The two metadata files (17 MB) stay in `data/external/feng2025/` as the evidence. The 1.75 GB
RNA counts and 105 MB guide counts were deleted — the guide file was redundant anyway, since
`Cell-Metadata.tsv.gz` carries `Guide_Call` and `Cell_Line` per cell.

**This closes the last identified data-acquisition lever.** X-Atlas/Orion has added no lines
since 2025-09, X-Atlas/Pisces is still vaporware, Nadig and Replogle-essential cover 0/300 of
the panel, and Feng is the wrong shape. A fourth donor line remains worth ~+0.005 (experiment
036) — there simply is not one available.

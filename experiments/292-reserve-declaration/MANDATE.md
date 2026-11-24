# 292 — Declare a genuinely independent RESERVE for the write-up's final confirmations: the single-cell perturbation datasets on disk that none of the programme's maps, banks or experiments has used — inventoried from METADATA ONLY, listed with sha256 and a frozen access rule in the ledger. No response is computed. ≤ 2 hours. No build, no submission.

You are a research subagent. Work in `experiments/292-reserve-declaration/`. **Edit nothing outside it** except appending
the reserve declaration to `docs/CONFIRMATION-LEDGER.md` as instructed below. No build, no `vcc`, no `cell-eval2`,
nothing under `submissions/` or `data/`; no download. `uv run python`. Peak RSS under ~6 GB (read `.obs` / `.var` and
attribute tables only — NEVER `X`, `layers` or any expression value; open h5ad files with `h5py` and read only
`obs`, `var`, `uns`; if a file's obs is not readable without loading X, record that and skip it). **Time box: 2 hours.**
One sibling experiment (285) runs in parallel; do not read or write its directory. `docs/RESEARCH-PROGRAMME-2026-09.md`
rule 2 binds you: a reserve is data whose OUTCOMES no agent has inspected; the point of this mandate is to create one
without inspecting anything but metadata. Both existing reserves are SPENT (C1: HEK293T half-B; C2: the four Broad-16
lines) and the surfaces used by 291 (Broad-16 working lines; GWCD4i's K562 comparison table) have been looked at once.

## Why (third external review, item C.3: "Reserve a genuinely independent dataset before inspecting its outcomes")
The scientific track now has a body of claims — the 28-bank cross-line conservation map (F8) that transfers off-map
(291), the family axis (F9, 284), the reproducible per-gene remainder (284) carried by response geometry with a confirmed
SIGN (F13, 290 / 291), the reliability laws (F1, F10; 285 pending) — and every surface on disk that carries them has now
been inspected at least once. A write-up needs at least one confirmation on data no agent has looked at. The
`data/external/scperturb/` collection holds ~50 harmonised datasets; the map used a specific subset. **Which of the rest
carry genetic perturbations in human cells with enough perturbations and cells to test (i) the map's per-gene
conservation score and (ii) F13's geometry sign, and what would the frozen access rule be?**

## Inputs
`experiments/276-conservation-map/work/inventory276_sources.csv` (the 28 map sources with their paths — every dataset
named there, or whose cells feed a bank there, is USED); `docs/FINDINGS-2026-09.md` §1 and F1–F13 (which datasets each
finding touched); `experiments/270-data-reconnaissance/RESULTS.md` (the dataset map, GPSAdb coverage); `data/external/
scperturb/*.h5ad` (metadata only); `data/external/{feng-ipsc, hcasmc-igvf, igvf-telohaec, hspc_crispri500}` (check
whether used); the scPerturb paper's dataset table (Peidli 2024, on disk or cite from memory with the caveat).

## Protocol
1. **Classify every dataset on disk** as USED (in the 28 sources or in any experiment's donor bank — grep the repo for its
   file name and its key) or UNUSED, with the evidence (the path or grep hit).
2. **For each UNUSED dataset, from `obs` / `var` / `uns` only:** organism, cell line / type, perturbation modality
   (CRISPRi / CRISPRko / CRISPRa / ORF / drug / other), number of distinct genetic perturbations targeting a gene,
   cells per perturbation (median, q10), number of control cells, batch / gemgroup / guide columns present, the target
   symbols' overlap with (a) 279's 17,106 held genes (`experiments/279-znf-state-dependence/work/held279.npz` `genes`),
   (b) the 300 challenge panel targets (`experiments/276-conservation-map/work/conservation276.csv` `target`), and
   (c) 290's feature store (`experiments/290-…/work/persrc290.npz` `genes`). Record the file's sha256 and size.
3. **Rank the UNUSED human genetic-perturbation datasets by usefulness** for the two claims: n perturbations in the held
   map, cells per perturbation ≥ 50, a control set, a line NOT among the map's 21 lines (state which are new lines and
   which repeat a map line — a repeat line still adds an independent STUDY).
4. **Write the reserve declaration** (in `RESERVE.md` and appended to `docs/CONFIRMATION-LEDGER.md` under "Reserve R3"):
   the files (path, sha256, size), what each contributes, and the frozen access rule: *outcomes (responses,
   conservation, agreement) may be computed only once, for a preregistered claim written into the ledger first, by an
   experiment the coordinator names; metadata reads are free; no exploratory look.* State what the reserve CANNOT
   confirm (e.g. modalities that differ from the map's CRISPRi; tens of targets per set → pooled claims only).
5. **Do not compute any response, any cosine, any conservation value, any feature.** If a step would require reading
   `X`, stop and record why.

## Deliverable
`RESERVE.md` (the inventory table over ALL datasets on disk with USED / UNUSED and evidence; the ranked reserve; the access
rule; the sha256 list), the ledger entry, `work/inventory292.csv`, `work/sha256.txt`; one line: **RESERVE R3 DECLARED:
n datasets, n human genetic-perturbation sets, n perturbations in the held map (pooled), n new lines; NOTHING COMPUTED
BEYOND METADATA**. Then stop.

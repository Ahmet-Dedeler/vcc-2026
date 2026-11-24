# 072 — Jiang pathway-DE transfer results

Date: 2026-08-23

## Verdict

**Reject as a production donor.**

The complete archive passes structural validation and shows a modest positive
development-only transfer signal, but contains only **9 of the 300** VCC panel
targets. This fails the predeclared minimum of 30 panel targets, so the
held-panel outcomes were never opened and no panel-level model choice was
made.

## Real guarded audit

- archive: `data/external/jiang2025/DE_results_all_pathway.zip`
- bytes: `324112587`
- MD5: `f077cba680a1affc599f5153d99b0e45`
- SHA-256:
  `d2f6d920abaf5e05de6e5343acdfae6a18537d23bf9003e80aa73cce97193c80`
- DE tables: 271
- unique targets: 218
- cell lines: A549, BXPC3, HAP1, HT29, K562, MCF7
- VCC genes represented: 9,540
- non-panel proxy targets shared by all three donor banks: 171
- process-tree peak RSS: **0.56 GiB**
- Python peak RSS reported by the experiment: **0.534 GiB**
- system free-memory floor observed by the guard: **89%**
- result: `negative_metadata_gate`

The nine overlapping panel targets are `ELK1`, `FOXO4`, `IFNAR2`, `IFNGR2`,
`MED15`, `MTF1`, `SMARCA5`, `STAT6`, and `ZNF22`.

The completed artifacts are under `out-pread/`; its `results.json` SHA-256 is
`02fe71a0180945a7c787cae6eb81a0c1ea965215ce433d2323d416783fcc2d89`.

## Development-only result

On 171 non-panel perturbations shared by K562-GWPS, Orion-HCT116, and
Orion-HEK293T, the selected raw Jiang blend at fraction `0.20` improved macro
PDS by `+0.006318`:

- worst held-donor delta: `+0.005366`
- all three held donors positive
- target-cluster bootstrap 95% interval:
  `[-0.002603, +0.015079]`

The interval crosses zero, so even the development evidence does not pass the
promotion gate. The low panel coverage independently prevents a defensible
VCC-panel decision.

## Memory correction

The first two guarded attempts stopped at the experiment's internal RSS ceiling
while scanning the 2.6 GB Orion donor arrays through memory maps. The numerical
chunks were small, but macOS retained mapped file pages. `prepare_bank_effects`
now uses the existing bounded `os.pread` row reader instead. The corrected run
completed under the original 1 GiB process-tree guard.

No H5AD was downloaded, no cloud job or `cell-eval2` run was used, and no
challenge submission was built or sent.

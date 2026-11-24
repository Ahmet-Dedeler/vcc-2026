# 068 — Seed-off-target source ablation results

Date: 2026-08-23

## Decision

**REJECT. Keep experiment-064 arm D unchanged.**

Surgically deleting the source rows implicated by the public seed-effect
tables made retrieval substantially worse. The mean PDS delta across the four
predeclared target clusters was `-0.135226360`, with all four target means
negative. Across the four native experiment-064 cases, the case-weighted mean
was `-0.143344710`.

The matched random-target/source null had mean `-0.048406741` and 95th
percentile `+0.036495291`. The observed four-target mean was at the `6.64`th
null percentile; `93.36%` of matched null draws were at least as good as the
reported-target rule.

The public rows are real, but deleting the entire donor effect is too blunt:
it discards useful on-target and downstream direction along with the suspected
off-target component.

## Exact per-case rank changes

Rank is one-based within the complete held-truth gallery; lower is better.
`rank delta = ablated - fixed D`, so a negative value is an improvement.

| target | held truth / frozen geometry | removed source | gallery | fixed rank | ablated rank | rank delta | PDS delta |
|---|---|---|---:|---:|---:|---:|---:|
| `PQBP1` | HCT116 / HCT geometry | K562/PerturbMatch slot | 290 | 1 | 1 | 0 | `+0.000000000` |
| `PQBP1` | HEK293T / HEK geometry | K562/PerturbMatch slot + HCT116 | 294 | 1 | 65 | +64 | `-0.218430034` |
| `PQBP1` | KOLF / HCT geometry | K562/PerturbMatch slot | 273 | 1 | 1 | 0 | `+0.000000000` |
| `PQBP1` | KOLF / HEK geometry | K562/PerturbMatch slot + HCT116 | 273 | 1 | 137 | +136 | `-0.500000000` |
| `SHPRH` | HEK293T / HEK geometry | HCT116 | 294 | 1 | 120 | +119 | `-0.406143345` |
| `SHPRH` | KOLF / HEK geometry | HCT116 | 273 | 218 | 253 | +35 | `-0.128676471` |
| `FZD2` | HEK293T / HEK geometry | HCT116 | 294 | 42 | 27 | **-15** | **`+0.051194539`** |
| `FZD2` | KOLF / HEK geometry | HCT116 | 273 | 38 | 93 | +55 | `-0.202205882` |
| `GLYR1` | KOLF / HCT geometry | CD4 mean3 | 273 | 97 | 98 | +1 | `-0.003676471` |
| `GLYR1` | KOLF / HEK geometry | CD4 mean3 | 273 | 115 | 124 | +9 | `-0.033088235` |

`FZD2` is the only favorable native case, improving 15 ranks on held
HEK293T. It reverses by 55 ranks on independent KOLF truth, so it is not a
transferable exclusion rule.

## Target-cluster summary

| target | cases | native mean delta | secondary mean delta | all-case mean delta | matched-null mean | null p95 | one-sided `P(null >= actual)` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `PQBP1` | 4 | `-0.109215017` | `-0.250000000` | `-0.179607509` | `-0.147844685` | `+0.117710019` | `0.615385` |
| `SHPRH` | 2 | `-0.406143345` | `-0.128676471` | `-0.267409908` | `-0.043228232` | `+0.148775660` | `0.953846` |
| `FZD2` | 2 | `+0.051194539` | `-0.202205882` | `-0.075505672` | `-0.026516387` | `-0.000641500` | `1.000000` |
| `GLYR1` | 2 | not estimable | `-0.018382353` | `-0.018382353` | `+0.022977941` | `+0.023805147` | `1.000000` |

The native target-equal mean over the three evaluable targets was
`-0.154721274`. `GLYR1` has no perturbation row in either Orion bank, so its
only valid test here is the explicitly secondary KOLF held-truth extension.

The exact-availability null contained 64 controls each for `PQBP1` and
`SHPRH`, two for `FZD2` (`EEF1A2`, `REV1`), and two for `GLYR1` (`CDK8`,
`UNK`). The small `FZD2` and `GLYR1` classes are reported as a limitation, not
treated as precise p-values.

## Source-table verification

The source repository was pinned at commit
`198fe3d84263002e9cfa17e1a2178d4e8dd0a669`. Every compressed table was
streamed from that commit, checked against its Git blob SHA-1 and byte count,
and scanned in full. Each requested pair matched exactly one row.

| source table | target → off-target | target log2FC | off-target log2FC | seed length | guide Hamming distance | neighbor rank |
|---|---|---:|---:|---:|---:|---:|
| Replogle K562 | `PQBP1` → `TIMM17B` | `-3.276709901` | `-3.071498405` | 19 | 1 | 1 |
| Xaira HCT116 | `PQBP1` → `TIMM17B` | `-2.258467295` | `-1.648017421` | 19 | 1 | 9 |
| Xaira HCT116 | `SHPRH` → `CENPP` | `-2.220201274` | `-0.917026680` | 13 | 6 | 19 |
| Xaira HCT116 | `FZD2` → `P2RX7` | `-1.360193405` | `-0.707038691` | 12 | 6 | 6 |
| Zhu CD4 D1 8-hour | `GLYR1` → `UBN1` | `-0.761859837` | `-0.695266013` | 20 | 0 | 13 |

| table | compressed bytes | rows | compressed SHA-256 |
|---|---:|---:|---|
| Replogle K562 | 14,460,423 | 214,025 | `820f9e5ac81bbdbc1700e4520d7f95b0da6034f7640da1f9c0577e535454db8d` |
| Xaira HCT116 | 21,669,263 | 406,838 | `f400e357464bcd9fd1b16d743c35695b181c686fd6e65d76d8d9309f9e043e95` |
| Zhu CD4 | 15,925,175 | 366,797 | `8fd1c0ae8c9841f1171266c57f0a724ff07ef8cafab1232376432252073a2218` |

Crossref resolved DOI `10.64898/2026.03.27.714658` to *Systematic
identification of seed-driven off-target effects in Perturb-seq experiments*,
published March 28, 2026.

`TIMM17B` is itself one of the 300 VCC targets, so its direct coordinate was
already excluded by the frozen metric geometry. `CENPP`, `P2RX7`, and `UBN1`
are not panel targets and remain ordinary scored coordinates. No off-target
coordinate was manually zeroed in this experiment.

## Frozen-pipeline and leakage audit

- The low-memory implementation reconstructed experiment-064 arm-D PDS
  exactly on all 179 evaluated rows that belong to the original 221-target
  gallery: maximum absolute error `0.0`.
- Both gamma-0.65 operators and their experiment-064 arm-B preserve scalars
  were frozen.
- PerturbMatch remained outside covariance fitting.
- All 300 panel coordinates were excluded from source effects, operators, and
  retrieval.
- Held perturbation effects were used only for final ranking. They did not
  select rules, masks, controls, or hyperparameters.
- KOLF was used only as independent held truth, never as a source or operator
  input.
- The intervention was limited to source-row availability. It did not inject
  the public log-fold changes or force any reported off-target gene.
- This is still a manually targeted public-data QC list. Its challenge-rule
  eligibility would need review before production use even if it had passed;
  because the result is rejected, it must not enter a submission.

## Resource and execution audit

- Authoritative run: one Python process, one BLAS thread, no child processes.
- Authoritative peak RSS: **0.607955933 GiB**.
- Hard cap: `0.75 GiB`; internal stop threshold: `0.70 GiB`.
- A first guarded attempt stopped without result artifacts when an avoidable
  retained buffer produced `0.711258 GiB` peak RSS. That remained below the
  hard cap. The buffers were removed before the authoritative rerun.
- Maximum RSS observed across all attempts: **0.711258 GiB**.
- Largest downloaded object: **21,669,263 bytes**, below 50 MB.
- Synthetic checks: `3 passed`.
- Independent output arithmetic audit: passed.
- No `cell-eval2`, cloud launch, count emission, VCC API call, or submission.

## Files

- `run.py` — one-process streaming evaluator and provenance verifier.
- `test_run.py` — deterministic stored-NPZ, ranking, and fallback checks.
- `metrics.json` — complete machine-readable result and hashes.
- `source_evidence.csv` — exact matched public source rows.
- `per_case_results.csv` — reported and matched-control rank/PDS outcomes.
- `per_target_results.csv` — four target-cluster summaries.
- `null_controls.csv` — matched-control target deltas.
- `reconstruction.csv` — experiment-064 arm-D reproduction audit.

## Final interpretation

Seed off-target annotation is useful donor-quality information, but whole-row
deletion is not the right operation here. The contaminated rows often retain
the strongest transferable identity signal for their intended perturbation.
If this mechanism is revisited, the next defensible test is a learned,
source-wide decontamination transform trained across many annotated
source/guide pairs—not another four-target lookup or full donor-row removal.

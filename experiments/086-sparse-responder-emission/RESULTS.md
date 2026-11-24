# Experiment 086 results

Date: 2026-08-23

## Decision

The construction is mathematically valid and testable. It exactly preserves
experiment 077's integer per-gene totals while keeping fixed legal row margins,
and it handles negative effects without negative counts.

**Current launch verdict: NO-GO.** It is not yet validated on HCT116 truth. A
read-only check at `2026-08-23T09:11:18Z` found zero objects under the pinned
experiment-077 output prefix
`output/exp077-low-alpha-union/20260823T060938Z/`. More importantly, the active
exp077 runner does not produce an independently immutable binding between its
source, launch command, summarizer, and terminal result hashes. Experiment 086
now fails closed until the dedicated SSM terminal-attestation procedure in
`README.md` is completed. No cloud job and no submission were launched.

## Independent-audit repairs

- All run ids, AWS identifiers, hashes, and optional prefix values are
  validated before path construction or shell interpolation.
- The worker rejects pre-existing, noncanonical, or symlink work paths and no
  longer contains a recursive-delete path.
- The launcher rejects existing local metadata and non-empty S3 input/output
  prefixes. It uses a unique temporary archive and conditionally reserves both
  S3 prefixes before uploading.
- Exp077's local launch metadata was copied byte-for-byte and pinned together
  with the launched source archive, runner, grid, summarizer, instance, SSM
  command, and prefixes.
- Co-located output checksums are no longer trusted as authentication. The gate
  requires AWS SSM command records plus canonical terminal-attestation stdout
  that binds every terminal result hash.
- The misleading aggregate feasibility boolean was replaced by an actual
  selected-subset feasibility diagnostic. The old ratio remains only as an
  explicitly named representative-subset screening statistic.
- Fractions `(0.35, 0.50, 0.65)`, seed `8602026`, and two-scorer cap are
  enforced inside `run.py`, including direct programmatic calls.

## Bounded validation

All 19 tests passed:

- exact columns, rows, determinism, integer nonnegativity, and row caps;
- unavoidable nonresponder clipping for infeasible negative effects;
- exact preservation of nonresponder rows when responder capacity is enough;
- pseudobulk invariance and non-invariance of per-cell CPM means/jackknife;
- experiment-086 jackknife agrees with `cell-eval2` 0.16.0;
- terminal experiment-077 artifacts are SSM-attested and semantically
  recomputed, including coherent rehash and command-body tamper rejection;
- the production provenance policy fails closed while its terminal SSM command
  id is absent;
- changed fractions, seed, and scorer count are rejected;
- unsafe identifiers and reused/symlink work directories are rejected without
  deleting sentinel files;
- AWS scripts are disjoint from experiment 080, require explicit clearance,
  reject prefix reuse, do not create/start instances, contain no recursive
  work-directory deletion, and contain no submission command; and
- the local memory guard terminates an over-limit child.

Measured process-tree peak RSS for the test suite was `0.399 GiB`, below the
hard `1 GiB` limit. The in-process Python peak was `0.241 GiB`.

## Synthetic Wilcoxon result

The representative fixture uses 2,000 heterogeneous, zero-inflated control
cells over 512 genes. Synthetic truth has 40% responders and 40 planted
effect genes. The low-alpha prediction carries 40% of the true aggregate
effect. Quasi and sparse predictions have exactly equal integer columns and
bit-identical scorer pseudobulks.

| Measure | Exact-total quasi | Sparse responders |
| --- | ---: | ---: |
| Predicted significant genes | 483 | 18 |
| Real significant genes | 38 | 38 |
| Jaccard | 0.0764 | 0.4737 |
| Direction accuracy on overlap | 0.9730 | 1.0000 |
| Calls outside planted effect | 446 | 0 |
| Jackknife correction | 0.4913 | 4.1730 |

In a null-margin control with the same exact columns, quasi called 484 genes
while the unchanged control-template layout called zero. This demonstrates the
specific mechanism: flattening a skewed control distribution can create
Wilcoxon calls even when pseudobulk is unchanged.

This fixture is deliberately designed to expose that mechanism. It establishes
testability and a plausible route to better Jaccard/fidelity/reach; it does not
establish that HCT116 or the hidden challenge contexts will improve.

## Full-shape resource result

One representative `400 × 18,533` block completed in `0.126 s`.

- process-tree peak RSS: `0.646 GiB`;
- in-process peak RSS: `0.618 GiB`;
- exact column totals: yes;
- exact baseline row margins: yes;
- negative counts: none;
- nonresponder clipping: zero;
- unchanged nonresponder rows: `200/200`;
- absolute response mass in responders: `100%`;
- maximum row library: `7,647`, below `1,000,000`.

## Remaining real-data experiment

Once experiment 077 is terminal, passes its own scientific gate, and receives
the independently verified SSM terminal attestation, the guarded runner will
reconstruct its selected alpha (`0.35` or `0.40`) exactly, build fractions
`0.35`, `0.50`, and `0.65`, and score exactly the two smallest fractions that
pass all construction and all-target jackknife gates.

The decisive outputs will be:

1. exact PDS equality to experiment 077;
2. no corrected-MSE regression;
3. raw and scaled changes in fidelity, reach, Jaccard, and NMAE; and
4. per-target leakage and responder-concentration distributions.

Until those terminal artifacts and immutable provenance evidence exist, the
correct status is **NO-GO, prepared but authentication-blocked by experiment
077**, not scored or promoted.

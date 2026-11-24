# Experiment 142 — batch-disjoint donor reliability weighting

Status: **terminal NO-GO**

Date: **2026-08-31**

## Decision

Within-donor split-half retrieval reliability did not improve target-specific
weighting among the three production donor slots. Five-fold nested selection
and the complete 613-target development fit both selected exact equal weighting:

```text
fold eta: 0, 0, 0, 0, 0
final eta: 0
```

The fixed and nested out-of-fold deltas were therefore exactly zero for both
raw PDS and signed credit. Only 2 of 59 frozen promotion checks passed.

Do not rescue, retune, integrate, or submit this reliability family. The
untouched 410-target confirmation split remains sealed.

## Primary result

| endpoint | nested delta | fixed delta | bootstrap 95% |
| --- | ---: | ---: | ---: |
| raw PDS | `0.0` | `0.0` | `[0.0, 0.0]` |
| signed credit | `0.0` | `0.0` | `[0.0, 0.0]` |

Every held-fold/family cell was exactly zero because every inner fold selected
`eta=0`. Both donor-assignment and target-pairing null families also selected
identity, so their p95 values were zero and the observed candidate did not
exceed them.

## Fixed nonidentity postmortem

The reliability candidate contained a weak directional signal, but not a PDS
gain:

| eta | macro raw-PDS delta | macro signed-credit delta |
| ---: | ---: | ---: |
| `0.25` | `-0.000020` | `+0.001183` |
| `0.50` | `-0.000450` | `+0.001855` |

The PDS loss was concentrated in the K562 and RPE1 families. This is the wrong
trade for the current submission, whose main transferable strength is
retrieval geometry.

## Non-promotable support controls

Simple donor-support proxies selected `eta=0.5` and were positive, but remained
below the frozen `+0.003` practical gates:

| control | fixed raw-PDS delta | fixed signed-credit delta |
| --- | ---: | ---: |
| target cell count | `+0.00196141` | `+0.00213555` |
| target total UMI | `+0.00199606` | `+0.00208881` |
| occupied-batch count | `+0.00191654` | `+0.00209870` |
| pooled-effect norm | `0.0` | `0.0` |

The three positive controls were highly redundant and their largest PDS gains
were concentrated in Jurkat. They are evidence for a small support/precision
effect, not authorization to choose a preferred proxy or alter production.

## Integrity and access state

- The authorized 613-target development split was opened once under the
  hash-bound development authorization.
- The 410-target confirmation split was not opened.
- Reserved Xu/Cao and Tian/Kampmann expression values were not opened.
- No count generation, `cell-eval2`, challenge build, submission, commit, or
  push occurred.
- The post-run audit validates the prefreeze authorization, prepared artifact,
  development authorization, immutable result manifest, and recorded
  outcome-access lock without reopening truth arrays.
- All 45 combined Experiment 142 tests pass.

Key SHA-256 values:

- development evaluator:
  `dee877786a5cf83819ce1f954ebc82f63179a36fa31e031fb47c45f2bfcc5c1a`;
- development prefreeze authorization:
  `204d3bdccf1d11b8c1b9a21942f3888156aad6631fc860cb527e9a3b286f333a`;
- development preoutcome manifest:
  `cfe55d1e3be5d117aba5316def2cde4a2f5b240d5985fcbba319f1ea476466a4`;
- development access authorization:
  `898b035ff508aeafc8684595af2cd3f0c5d2dbf779e6ff24cdf6aefe44734e3e`;
- development result manifest:
  `a2d730451497b33aab3eeb9181a159fe72f2732cce75e62d6974d3f9ddee44e6`;
- development result:
  `46715620f4fac830be615a7bf82bfd30da16610e7816a2a94dd29fdb6e5206c4`.

## Follow-on

Experiment 097 previously reported a small positive result from
`normrel^-1.5` target-specific donor weighting. Experiment 142 removes the
old claim that exact target-level support statistics are unavailable, but its
support-control gains are too small and family-concentrated to promote
directly. Any follow-on must use a new experiment number, preregister one
latent support/precision score before outcome access, repeat nested selection,
and enforce strict family-concentration gates.

The status line in `README.md` is a historical, source-bound pre-run statement
and is intentionally not edited after the frozen evaluation. This result file
supersedes it for the terminal decision.

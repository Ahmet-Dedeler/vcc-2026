# 114 — DepMap co-essentiality for weak incumbent rows

**Verdict: terminal NO-GO after the authorized HEK293T PDS diagnostic.**

The source artifacts, role-specific a14 reconstructions, candidate matrices,
and 100 target-label permutation candidates were frozen before the held
HEK293T effect array was authenticated. The diagnostic covered the 35 frozen
eligible targets present in the HEK293T gallery and excluded all 300 VCC panel
response coordinates.

| arm | raw PDS | delta from frozen incumbent |
| --- | ---: | ---: |
| `current_fallback` / `identity` | 0.704202 | — |
| `string_neighbor` | 0.561345 | **-0.142857** |
| `depmap_neighbor` | 0.548739 | **-0.155462** |
| `string_depmap_blend_w1` | 0.584034 | **-0.120168** |

The primary blend also failed the frozen target-label null: its raw-PDS delta
was `-0.120168`, while the 100-draw null 95th percentile was `-0.119328`.
Every nonidentity arm was substantially worse than the exact incumbent.

This result closes the proposed replacement of a14's 41 weak-support rows by
leave-target-out STRING, DepMap, or blended co-essentiality neighbors. The
authorized evaluation was explicitly diagnostic-only and did not include a
directional endpoint. There is no scientific reason to open the H1 held
surface: the preregistered promotion gate required nonnegative PDS on every
clean surface, and the first authorized surface failed by a large margin.

No model was fit after truth access. No `cell-eval2`, count generation,
challenge build, submission, commit, or push occurred.

## Provenance and verification

- Frozen HEK293T candidate manifest SHA-256:
  `b27b6bfadc0c2a3e0edfafadded3430fcba2b5e64c6060ddaeb46c312f002ef2`
- Held-evaluation approval SHA-256:
  `9fd85f9df8cce2d917dcd5844eacccf9ba40fff4d70046a34d68bd018547e72c`
- Outcome access receipt SHA-256:
  `10408c9bcb2079beec70e0dd802ddd27452fc9de8d4ef87349b34aa7fdec87a2`
- Result SHA-256:
  `92586959e63c344ed0d5511851e515dc31eeb74f2efd24f34a4f70cd8100fdfb`
- The pretruth evaluator suite passed 8/8 before outcome access, as sealed in
  the approval provenance. A post-result rerun passed the seven reusable
  checks; the remaining pretruth-only test correctly refused because the
  immutable output directory now exists.

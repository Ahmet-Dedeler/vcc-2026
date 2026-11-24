# 031 — Whitening confirmed on somatic VCC-panel data; sharpening refuted

**Date:** 2026-08-22 · real `cell-eval2` on `data/harness/orion-hct116.h5ad` — **132 actual
VCC-panel targets** in HCT116 (colorectal, somatic), 18,400 control cells, predicted from
K562-gwps + Orion-HEK293T with **HCT116 dropped from the donor set**. Raw metrics rescaled with
the official 2026 anchors.

This is the first count-level measurement in the project on the right biology *and* the right
gene population. Everything before it was 40 perturbations of an embryonic line on the
common-essential-enriched 2025 panel.

## 1. Whitening transfers, and it is large

| | raw pds | scaled pds | nmae | scaled nmae | avg6 |
|---|---|---|---|---|---|
| gamma = 0 | 0.6874 | 0.414 | 0.9485 | +0.088 | 0.0816 |
| **gamma = 0.85** | **0.7633** | **0.582** | 0.9609 | +0.067 | **0.1042** |

**+0.076 raw pds, +0.168 scaled, +0.023 on the six-metric average**, measured by the real
scorer on somatic data with the actual panel. The `nmae` cost is −0.021 scaled — the price the
lower transfer cosine implies (log, 2026-08-21), and an eighth of the gain.

`fid` is **0.5220 against the official baseline of 0.5135** — positive, which no team on the
board manages. `jac` 0.0343 against a baseline of 0.029. Both are harness values rather than
predictions, but they are the first sign that the DE metrics behave differently on somatic
panel genes than on H1's essential ones.

## 2. Sharpening is refuted, and the proxy was wrong by 4x

| q | raw pds | scaled pds | raw nmae | scaled nmae | **avg6** |
|---|---|---|---|---|---|
| **1.0** | 0.7633 | 0.582 | 0.9609 | +0.067 | **0.1042** |
| 1.5 | 0.8068 | 0.678 | 1.1208 | −0.199 | 0.0803 |
| 2.0 | 0.8319 | 0.733 | 1.2596 | −0.430 | 0.0558 |

The `pds` gain is real and bigger than LODO suggested — raw `pds` 0.7633 → 0.8319 is enormous.
**And it is not worth having.** `nmae` goes to 1.26, and the six-metric average falls from
0.1042 to 0.0558.

The experiment 018 proxy put `nmae` at q = 2 at **1.05**; the real metric says **1.26**. It
understated the cost by roughly **4x**, which is more than enough to invert the decision. The
proxy was already known to be miscalibrated in alpha (it puts the `nmae` optimum at 0.75 where
the count-level run measures 1.08); this shows the same weakness on the other axis.

**This is why v7 was built and held rather than built and submitted.** The LODO evidence for
sharpening was strong, replicated on two harnesses, and mechanistically clean — and it is
wrong about the score. `pds` is one metric of six, and a transform that trades alignment for
discrimination pays for it somewhere.

## 3. alpha is flat between 0.7 and 1.08

avg6 0.1052 at alpha = 0.7 against 0.1042 at 1.08, `nmae` 0.9618 vs 0.9609 — a tie inside
noise. v6's alpha, chosen for magnitude parity with v4b, is fine.

## Caveat on `reach`

Raw `reach` here is 0.012–0.018 against an official baseline of 0.072, so it scales negative
for every configuration. Orion runs at ~154 cells per perturbation, so this reference has 150
where the competition has 400, and the reference's own DE calls are correspondingly fewer and
noisier. `reach` on this harness is not comparable to the board and should be read only within
the run, if at all.

## 4. The variance channel replicates on somatic data, modestly

| | pds | nmae | reach | fid | jac | avg6 |
|---|---|---|---|---|---|---|
| even emitter | 0.7633 | 0.9609 | 0.0124 | 0.5220 | 0.0343 | 0.1042 |
| **quasi, k=1000, spread 0.15** | 0.7624 | 0.9614 | **0.0140** | **0.5251** | 0.0342 | **0.1057** |

`reach` +13% relative, `fid` +0.0031, and `pds` / `nmae` / `jac` unmoved — which is what the
construction guarantees, since the emitter preserves each gene's panel total. The direction
matches `h1-fast` (+0.028 scaled `reach`, 3.4σ) on a different line and a different panel.

The absolute `reach` values here are not comparable to the board (see the caveat above), so
read the sign and the ratio, not the size.

## 5. The CPM gate result **reverses** on somatic data — and v6 shipped the other way

| | pds | nmae | jac | avg6 |
|---|---|---|---|---|
| CPM gate 5 (**what v6 ships**) | 0.7633 | 0.9609 | 0.0343 | 0.1042 |
| **CPM gate 20** | **0.7727** | 0.9634 | **0.0360** | **0.1074** |

Gate 20 wins by +0.0032 avg6, carried by `pds` (+0.0094 raw) and `jac`. Experiment 023 found
the opposite on `h1-fast`, six configurations out of six, and that is why v6 ships gate 5.

Both cannot be right, and the disagreement is exactly along the axis the two harnesses differ
on: `h1-fast` scores the 2025 common-essential-enriched panel in an embryonic line, this
scores the actual VCC panel in a somatic one. **The somatic harness is the one that matches the
task**, so gate 20 is probably the better setting — but the margin is small, LODO (exp 017)
had them within 0.003 of each other in the other direction, and this reference carries 150
cells per perturbation rather than 400.

Filed as unresolved, worth one clean test rather than a guess. v6's gate 5 cost at most
~0.003 avg6 on this harness, and v6 scored 0.1369 regardless.

## 6. alpha, on the right biology

| alpha | 0.7 | **1.08** | 1.5 |
|---|---|---|---|
| avg6 | 0.1052 | 0.1042 | 0.0909 |
| nmae | 0.9618 | 0.9609 | **1.0061** |

Flat between 0.7 and 1.08, and 1.5 is clearly worse — `nmae` crosses 1.0, i.e. becomes worse
than predicting nothing. Note this is the **pre-jackknife** picture: experiment 033 shows that
once emitted variance buys the `mse` credit, the alpha argument changes sign, which is what
the queued sweep tests.

## 7. gamma on the somatic harness — 0.85 confirmed

At the shipped operating point (spread 0.25, gate 5, alpha 1.08), on 132 real VCC-panel targets:

| gamma | 0.40 | 0.65 | **0.85** | 1.10 |
|---|---|---|---|---|
| avg6 | 0.1026 | 0.1080 | **0.1089** | 0.1094 |
| pds | 0.7399 | 0.7577 | 0.7620 | 0.7658 |
| nmae | 0.9576 | 0.9595 | 0.9609 | 0.9618 |

avg6 rises monotonically through 1.1 and is **flat from 0.85 onward** (0.1089 vs 0.1094, a
difference of 0.0005 — well inside noise). LODO, which computes `pds` exactly on 259 panel
perturbations, peaks at 0.85 (0.8078 against 0.8068 at 1.05 and 0.8050 at 1.2).

**0.85 sits on the plateau of the somatic harness and at the exact optimum of LODO.** No change.

Note the `nmae` gradient runs the other way — it improves as whitening weakens (0.9618 → 0.9576
from gamma 1.1 to 0.4), which is the cosine tension in §Addendum 2 showing up per-metric. The
`pds` gain dominates it by roughly 8x, so the trade stays correct, but it confirms that `nmae`
is structurally opposed to our main lever rather than merely under-tuned.

# 176 — LINCS/CMap 2020 as donor lines: **NO-GO** for both modalities — shRNA does not transfer; CRISPR-KO transfers on landmark genes only, and coverage caps it below +0.001

Date: **2026-09-02**, coordinator-run after the agent failed to authenticate (see `INVENTORY.md`).
Measurement only. Surface: h1-test reference DE pool (`164/de164.a15.npz`, 85 targets × 10,779
genes). Data: `level5/level5_beta_trt_sh_n238351x12328.gctx` (11.8 GB, MODZ z-scores), the 36,720
`CGS…` consensus-gene-signature columns (one per gene per line; multiple shRNAs collapsed, which
cancels seed effects). Scripts: `cmap176.py`; per-line matrices in `work/cgs_<line>.npz`.
Peak RSS ~3 GB (the CGS block is the first 1.8 GB of the file).

## Coverage (from `INVENTORY.md`)
9 core lines cover **98–118 of the 300 panel targets each** (union 127/300) and 50–57 of h1-test's
85. On our DE-universe axis: **867 landmark genes** (directly measured) and **7,101 inferred**.

## Gate 1 — platform agreement with our production donors, on the H1 DE pool

| pair | landmark genes | inferred genes |
|---|---:|---:|
| **K562 × HCT116** (our two donors, reference) | **0.540** | 0.529 (all) |
| CMap MCF7 × K562 | 0.531 | 0.504 |
| CMap HEPG2 × K562 | 0.525 | 0.495 |
| CMap A549 × K562 | 0.521 | 0.511 |
| CMap HT29 × K562 | 0.518 | 0.509 |
| CMap VCAP × K562 | 0.514 | 0.502 |
| CMap PC3 × K562 | 0.511 | 0.496 |
| CMap HA1E × K562 | 0.508 | 0.505 |
| CMap A375 × K562 | 0.504 | 0.498 |

CMap agrees with K562 less than HCT116 does on landmarks, and **at chance on the inferred genes**.
The inferred axis carries no signal shared with a CRISPRi screen.

## Gate 2 — sign accuracy against the H1 reference DE pool

| source | landmark acc (n) | inferred acc (n) |
|---|---:|---:|
| K562 | **0.562** (4,829) | **0.568** (28,771) |
| HCT116 | 0.527 | 0.537 |
| HEK293T | 0.537 | 0.528 |
| production trio, equal mean | 0.557 | 0.556 |
| CMap A375 | 0.499 | 0.496 |
| CMap A549 | 0.509 | 0.503 |
| CMap HA1E | 0.513 | 0.494 |
| CMap HCC515 | 0.514 | 0.500 |
| CMap HEPG2 / HT29 / MCF7 | 0.50–0.53 | 0.49–0.51 |
| CMap PC3 | 0.512 | 0.513 |
| CMap VCAP | 0.521 | 0.513 |
| **CMap 9-line mean** | **0.534** (4,189) | **0.522** (29,154) |

Every CMap line alone is between 0.48 and 0.55 on landmarks; the nine-line mean reaches 0.534 —
**a real but weak signal, well below a single CRISPRi donor (0.562).**

## Gate 3 — added to the production trio (both sources standardised), on CMap-covered pairs

| weight on CMap mean | landmark | inferred | all covered |
|---:|---:|---:|---:|
| 0 (trio alone) | 0.5582 | 0.5489 | 0.5501 |
| 0.25 | 0.5548 | 0.5499 | 0.5506 |
| 0.5 | 0.5617 | 0.5495 | 0.5510 |
| 1.0 | **0.5631** | 0.5478 | 0.5497 |
| 2.0 | 0.5493 | 0.5419 | 0.5429 |

**+0.005 on 867 landmark genes at best, nothing on the 7,101 inferred genes, ≤ +0.001 overall** — on
the ~1/3 of targets CMap covers. That is below the noise of every downstream member.

## Mechanism — the perturbations are weak

Sanity check on the sign convention doubles as the diagnosis. In each knockdown's own consensus
signature, the **target's own gene** reads:

| line | median own-`z` | fraction negative |
|---|---:|---:|
| A375 | −0.63 | 0.66 |
| A549 | −0.46 | 0.77 |
| HA1E | −0.46 | 0.75 |
| HT29 | −0.48 | 0.69 |
| MCF7 | −0.57 | 0.74 |
| PC3 | −0.53 | 0.77 |
| HEPG2 | −0.71 | 0.78 |
| VCAP | −0.85 | 0.76 |

The convention is right (negative = down), but **the knockdown itself is only 0.5–0.9 SD in the
consensus, and a quarter to a third of targets show no knockdown at all.** Our donors are CRISPRi
screens filtered at >80% knockdown. A perturbation that barely moves its own target cannot carry a
transferable downstream response, whatever the platform. Bulk L1000 at 96 h and cancer-line
biology are secondary to this.

## Verdict
**NO-GO for shRNA.** CMap's shRNA consensus signatures are not donors in the sense this project
needs. The CRISPR-KO file (`trt_xpr`, complete knockout, 82 panel targets, 41 h1 targets) is the
remaining modality and separates "knockdown depth" from "platform" — pending below.

## CRISPR-KO (`trt_xpr`) — measured (`cmap176_xpr.py`, `work/cmap176_xpr.log`)

`level5/level5_beta_trt_xpr_n142901x12328.gctx` (6.5 GB). 1,974 signatures hit h1-test targets
across 21 lines; collapsed to one mean `z` per (line, target), median 2 signatures per target. The
nine-line core library covers **41/85** h1-test targets (82/300 panel); twelve smaller lines cover 17.

**Own-target `z`** (a KO need not remove the transcript, so this is a weaker depth proxy than for
shRNA): core lines −0.26 to −0.80, 58–92% negative; the 17-target lines −0.95 to −2.5.

**Gate 1 — agreement with K562 on the H1 pool:** 0.49–0.54 landmark, 0.49–0.53 inferred (reference
K562 × HCT116: 0.540). Same as shRNA.

**Gate 2 — sign accuracy vs the H1 pool:**

| source | landmark (n) | inferred (n) |
|---|---:|---:|
| K562 (reference) | 0.562 | 0.568 |
| production trio | 0.557 | 0.556 |
| individual KO lines (41-target) | 0.506–0.545 | 0.495–0.518 |
| individual KO lines (17-target, noisy) | 0.494–0.577 | 0.471–0.570 |
| **KO 21-line mean** | **0.5550** (2,935) | **0.5196** (20,316) |

**Gate 3 — trio + w × KO-mean, on KO-covered pairs:**

| w | landmark | inferred | all covered |
|---:|---:|---:|---:|
| 0 | 0.5522 | 0.5524 | 0.5524 |
| 0.5 | 0.5635 | 0.5514 | **0.5530** |
| 1.0 | 0.5673 | 0.5488 | 0.5511 |
| 2.0 | **0.5714** | 0.5409 | 0.5448 |

**On the 867 landmark genes the 21-line KO mean adds +0.019 sign accuracy to the trio — a real,
transferable signal.** On the 7,101 inferred genes it subtracts (−0.011): the L1000 inference model
is noise for this purpose, in both modalities. Net over all covered pairs: **+0.0006.**

## Verdict, both modalities

- **shRNA: does not transfer** (0.534 vs K562 0.562 on landmarks), mechanism = weak knockdown.
- **CRISPR-KO: transfers on landmark genes only** (+0.019 on 867 genes over 41/85 targets), which is
  an existence proof that bulk L1000 landmark measurements of a knockout carry sign information
  about a single-cell CRISPRi response in hESC. **The inferred genes never do.**
- Coverage arithmetic: 82/300 panel targets × ~13% of pool pairs that are landmarks × +0.019 ≈
  **+0.0007 pooled sign accuracy** on the challenge. Below every member's noise. **NO-GO as a
  donor.** If a future arm ever needs sign help specifically on landmark genes for the 82 covered
  targets, the KO tables are in `work/xpr_<line>.npz`.

## What this closes
LINCS/CMap as a donor source for this challenge — both perturbation modalities, both gene axes. With
166 (single-cell CRISPRi exhausted) this closes public cross-line perturbation data as a route to a
better predictor at panel scale. The remaining information sources are ones nobody has: private
screens, or the D/E/F controls themselves once released (which carry no perturbation information).

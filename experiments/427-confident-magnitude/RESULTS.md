# 427 — Donor-confident genes: 87–94 % sign-accurate on H1's gates and under-sized in our rows. Raising them buys nmae −0.005 analytically; reach does not move (analytic only)

Status: 2026-09-23 05:15Z, PI, analytic on H1 (`conf427.py`). No emission and no builder.
- Base: 426's `base_s0_bps` delivered statistics (the shipped 410b construction), nmae 0.95237.
- The reach walk reproduces the delivered reach exactly (0.27474), which serves as the known answer.

**Construction.** For each h1-test target p and gene g, the six non-H1 donor tables vote:
- K562-GWPS covers 82/85 targets; K562-essential 18; RPE1 21; Nadig Jurkat 21; Orion HCT116 84; Orion HEK293T 84.
- Each vote is the donor's pseudobulk log1p-CP50k delta.
- A cell (p, g) is **confident** when ≥ n donors cover it, all agree in sign, and every |δ| ≥ τ.

| n ≥ | τ | confident cells / target | on the gate | sign accuracy on the gate | median true \|l\| | our median \|l̂\| |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.1 | 80 | 731 | 0.776 | 0.45 | 0.21 |
| 3 | 0.1 | 22 | 391 | **0.887** | 0.55 | 0.29 |
| 3 | 0.2 | 1.8 | 77 | **0.974** | 1.35 | 0.77 |
| 4 | 0.1 | 4 | 225 | **0.947** | 0.59 | — |

Our own delivered sign on those cells is already right 0.72–0.94 of the time; the magnitude is what is missing.

**nmae (Δ vs base):**
- *Replacing* our lfc by κ·mean|δ|/ln 2 loses (+0.002 … +0.008). The donors' means are smaller than our already-large
  head magnitudes on these cells.
- *Raising only*, lfc = s·max(|ours| where the signs agree, κ·mean|δ|/ln 2):

  | n ≥ | τ | κ 1 | κ 1.5 | κ 2 | κ 3 |
  |---:|---:|---:|---:|---:|---:|
  | 3 | 0.1 | −0.0006 | −0.0016 | −0.0030 | **−0.0051** |
  | 3 | 0.2 | −0.0003 | −0.0011 | −0.0024 | −0.0049 |
  | 2 | 0.1 | −0.0006 | −0.0012 | −0.0020 | −0.0023 |

  The truth on these cells is about 2–3× the donors' mean |δ|/ln 2.
- **reach:** putting the confident cells first in our order (with the donors' sign) moves reach by only +0.0002.
  They are already at the top (called head genes), and reach lives in the weak targets, where confident cells are rare
  (426 `weakreach426.py`: 406's weak tercile).

**Delivery.**
- Through the A-channel the moves are too large: up to ~1 log2 against the class structure's ~±0.5 log2 and its
  capacity.
- It would have to be a **carrier-profile** change: raise the confident cells' carrier effect, and let ballast (ratio
  mode) pin the bulk, as 416's carrier menu did.
- That makes it a carrier-design change, the class where H1 mis-signed nmae twice (392, 417c). It is priced as a
  hypothesis until an H1 emission run and a board read exist.

**Queue position:** after 426a's board read.
- If the board reads the MDE at ≥ 1× H1, run an H1 emission test (426's harness with a modified carrier matrix: 3 seeds,
  all six members).
- Otherwise shelve it with the rest of the DE-magnitude family.

## H1 emission test (05:05–05:15Z): the carrier-profile raise is CATASTROPHIC as built. Shelved
`run426.py base S --bulk ps --conf-raise κ` raises the carrier effect on the confident cells (≥ 3 donors agree,
|δ| ≥ 0.1, panel columns excluded; 22 cells/target) to κ·mean|δ|; the bulk is untouched.
- κ 2: nmae **+0.050**. κ 3: nmae **+0.082**; the other members are ~flat.
- **Mechanism: zeroed genes.**
  - Strong down-raises on low-expressed genes leave a gene with 0 counts in every carrier. Ratio-mode ballast
    (`tgt = carrier composition × bulk/DE ratio`) then assigns it 0 as well.
  - The delivered arithmetic mean is 0, so lfc = log2(1e-9/ctrl) ≈ −35: 61 cells at κ 3, against 0 non-target cells in
    the base. The base's single −34.7 cell is NGRN in its own knockdown, which is excluded from every metric.
  - One such gate cell adds ~33 to its target's nmae numerator.
- A raise would need a depth floor (e.g. effect ≥ ln 0.25) or ballast that restores dropped genes. The analytic prize
  (−0.005) comes mostly from exactly those deep down-raises, so it is shelved.
- The MDE (426) cannot create zeros: A-channel moves keep every per-gene total, so the ballast input is unchanged.

## Per-target magnitude (427b–d, `pertarget427.py`): nothing oracle-free predicts the per-target scale. Closed
- The per-target oracle scale k_p* (grid 0.5–8) is worth nmae **−0.0229** on H1, which reproduces 388's −0.023.
- k_p* follows the target's own gate sign accuracy (Spearman +0.38), which is unobservable.
- Oracle-free predictors of k_p* (Spearman):
  - our row norm −0.34;
  - donor effect norms −0.10 … −0.26;
  - donor pairwise row consistency −0.27. Consistency does not predict sign accuracy either: +0.02.
- **Rules (Δ nmae):**
  - Equalising delivered norms across targets, (median/norm)^β: +0.0000 … +0.0006 at ×1, and −0.0045 at ×1.25, which is
    only the global ×1.25.
  - Scaling by consistency tercile: −0.0018 … +0.0085.
- Do not re-propose a per-target magnitude rule without a feature that predicts the target's transfer sign accuracy.

## 427e: a plain donor average has better gate signs than our delivered LFC on H1 (`donoravg427.py`, analytic)
Donor average: l̂ = w · mean over the six non-H1 donors covering (p, g) of δ/ln 2.
- **Gate sign accuracy: donor average 0.5776 (covers 99.5 %) vs our delivered 0.5606.**
- Magnitudes are equally small: median |l̂|/|l| is 0.165 for the average and 0.163 for ours. Gate genes are
  significance-selected, so their |l| is inflated relative to any unbiased predictor, which is the MDE's point.
- corr(|l̂|, |l|): donor average +0.14, ours +0.23.

| prediction (Δ nmae vs our delivered 0.95237) | alone | + MDE add (c 0.1, cap 0.1, min 0.3) in its own sign |
|---|---:|---:|
| ours | 0 | −0.0098 |
| donor avg × 1.0 / × 1.5 | −0.0007 / −0.0090 | −0.0168 / −0.0118 |
| blend 0.5·ours + 0.5·donor × 1.0 | −0.0033 | −0.0205 |
| **blend 0.5·ours + 0.5·donor × 1.5** | −0.0120 | **−0.0245** |
| blend 0.5·ours + 0.5·donor × 2.0 | −0.0169 | −0.0229 |

**Reading.**
- On H1, the DE-channel LFC target that pays is our delivered LFC blended with the unwhitened donor average, plus
  the MDE add. That is 2.5× the plain MDE analytically (+0.0068 overall at board slopes × 1).
- **Caution:** the donor average is an unwhitened composition. 417c's γ 0 carriers were the same kind of change: H1 said
  nmae −0.0068 and the board read +0.0089, because the cancer lines' gates sit on line-specific genes.
- Delivering it means changing signs on the non-block genes, which is the A-channel flip problem (426: class over-fill)
  or a carrier-design change.
- **Queue:** after 426a's and 428a's board reads. If the MDE transfers, test the blend first on H1 through the emitter:
  carrier-profile blend (no zeros, checked with `audit428.py`) + MDE. Then price it on the board as a design change,
  with the nmae prior centred on 0, not on H1's number.

## 427f: alternative call orders for reach (`reachorder427.py`): the delivered order is best. Closed
Signs are held at our delivered lfc and only the order within the called and uncalled blocks changes. H1 reach:

| order | reach |
|---|---:|
| delivered, by p (known answer) | 0.27474 |
| by \|our lfc\| | 0.25717 |
| by \|donor avg\| | 0.24874 |
| by donor agreement, then \|avg\| | 0.25474 |
| by p, sign-disagreeing-with-donor-avg genes demoted | 0.27414 |
| same, only where \|avg\| > 0.05 | 0.27597 (+0.0012, inside noise) |

## 427e CORRECTED by 430 (2026-09-23 06:35Z): the donor-average sign advantage does not transfer to the panel
- H1's +0.017 comes from K562-essential, RPE1 and Nadig-Jurkat. They cover 18–21 of H1's 85 targets (60 % of H1's
  real-gate cells) and **0 of the 300 A/B/C targets** (the panel has no essential genes).
- With only the panel-available donors (K562-GWPS and the two Orion lines), the H1 comparison reads −0.0017.
- Cross-study on held-out cancer lines, the donor mean vs our whitened composition is +0.0012 (CI spans 0), and a blend
  +0.0018.
- The blend's increment over the plain MDE on panel-available donors is −0.0014 nmae (published: −0.0147).
- **Closed: no carrier-sign slot.** The panel-transferable part of the idea is the MDE itself (426a). See
  `experiments/430-cancer-fold-signs/RESULTS.md`.

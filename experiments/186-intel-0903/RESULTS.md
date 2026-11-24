# 186 — Competitive intel and literature refresh, 2026-09-03

**Question owned:** what mechanism gives raw `reach` ≥ 0.25 at raw `pds` ≥ 0.84?
**Surfaces used:** the live board (612 scored teams + 14 archived snapshots), every `model_name` and
`description` on it, the PertMind paper and repo, arXiv since 2026-08-01, GitHub/HF.
**Nothing was posted, no account was created, nobody was contacted, and no query or public field
carried anything about our method** (`AGENTS.md` §14).

Scripts written here, all with a known-answer row that aborts on mismatch, outputs saved as
`work_*.txt`: `reach_anatomy.py`, `natural_experiments.py`, `since_aug28.py`, `residual_rank.py`
(plus the previous run's `parse_board.py`, `clusters.py`, `diff_desc.py`, `traj.py`).

---

## 0. The answer in five lines

1. **It is not the direction.** Three teams beat us on the sign-accuracy channel (`nmae`) and lose
   0.06–0.07 raw `reach` to teams they beat. 15 of 49 high-reach × low-reach pairs inside the shrunk
   cluster are sign-inverted on `q` at **every** common magnitude ratio `k`.
2. **It is not the amplitude.** Across 29 same-team amplitude arms, `d reach = +0.0127 · d log(mse)`
   with residual sd `0.0363` — a 2.9× noise-to-signal ratio. Halving raw `mse` costs `−0.009` raw
   `reach`; five teams *gained* reach while halving it.
3. **It is not the architecture.** Three teams swapped the named model (STACK→Lingshu-Cell,
   CellFlow→Lingshu-Cell, CPA→STATE) and every one of the six members moved by less than 0.006.
4. **It is the ordering.** Six same-team arms moved `reach` by **+0.025 to +0.046** with `pds` inside
   ±0.006, `log mse` inside ±0.12 and `nmae` flat or *worse*. The best of them (`luka`, `PertMind_Pro`,
   2026-09-02) is +0.0374 raw = **+0.0070 overall** with the profile pinned.
5. **Priced for us:** on the 41 coverage-saturated teams with `pds > 0.79`, regressing raw `reach` on
   `pds + nmae + log mse + jac` (R² = 0.582) puts our residual at **−0.0411 raw reach = −0.0077
   overall, −1.92 sd, 4th-most-negative of the 41.** That is the size of the prize, and it needs no
   change to the direction, the amplitude, the profile, or the call-set size.

---

## 1. The refreshed board, and what moved since 2026-08-28

Snapshot `board/leaderboard_2026-09-03.json`, 613 rows (612 scored + 1), pulled 00:30 UTC.

### 1.1 Scorer stability — a control worth stating

Across the 15 snapshots there are **2,760 re-observations of the same `(team, submission_date,
submission_count)` row. Zero drifted on any of the six members or on `score_avg`.** The scorer, the
reference bundle and the anchors have not changed since Aug 21. Every delta below is an arm effect,
not a re-scoring. `cell-eval2` upstream is unchanged since the 0.16.0 import of 2026-08-20;
`ArcInstitute/state`'s last commit is 2026-07-24 (`normalize raw counts in model workflows`) — **no
new State release, so `experiments/182-pretrained-state/RESULTS.md` stands untouched.**

### 1.2 Field movement

| | Aug 28 13:25 | Sep 3 00:30 | Δ |
|---|---:|---:|---:|
| teams scored | 385 | 613 | +228 |
| 1st | 0.21122 | 0.21790 | **+0.0067** |
| 5th | 0.16890 | 0.19960 | +0.0307 |
| 10th | 0.15944 | 0.18904 | +0.0296 |
| 20th | 0.14853 | 0.18099 | +0.0325 |
| 50th | 0.11521 | 0.14947 | +0.0343 |
| teams over 0.170 | 4 | **30** | +26 |
| teams with raw `mse` < 1.0 | 8 | **44** | +36 |
| teams with raw `pds` > 0.82 | 6 | 21 | +15 |
| teams with raw `reach` ≥ 0.25 | 8 | 25 | +17 |

**The leader barely moved and the field closed on it.** In six days the 20th place gained five times
what 1st gained. Two waves did it:

- **The α-shrink wave.** Teams below raw `mse` 1.0 went 8 → 44. Our one unclaimed lever as of v1 is
  now the field's default. Sixteen teams did the whole 3–6 → 0.8 move inside this window.
- **A `reach` wave on 2026-09-02**, described in §3.

### 1.3 Top 22 now, against Aug 28

| # | team | now | Aug 28 | Δ | Δpds | Δmse | Δreach | Δnmae |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | cqawesome | 0.2179 | 0.1412 | +0.0767 | +0.0753 | −4.851 | **+0.0384** | −0.0435 |
| 2 | Aginglab.com | 0.2129 | 0.2112 | +0.0017 | +0.0025 | −0.007 | −0.0037 | −0.0009 |
| 3 | Andrew Liu | 0.2083 | NEW | | | | | |
| 4 | Hongjian Zhou | 0.2078 | 0.1193 | +0.0885 | +0.0560 | −8.065 | +0.0118 | −0.1066 |
| 5 | Garrrbage | 0.1996 | 0.1476 | +0.0520 | +0.0279 | −2.087 | +0.0055 | −0.0154 |
| 6 | Vivai | 0.1922 | NEW | | | | | |
| 7 | Slurm Shadi | 0.1918 | 0.1149 | +0.0769 | +0.0975 | −1.482 | **+0.0719** | −0.0249 |
| 8 | Jurassic Park | 0.1897 | 0.1594 | +0.0303 | +0.0297 | −0.396 | −0.0041 | +0.0131 |
| 9 | Patheia | 0.1894 | 0.1150 | +0.0743 | +0.0321 | −3.205 | −0.0068 | −0.1346 |
| 10 | Happy | 0.1890 | 0.1142 | +0.0749 | +0.0890 | −3.808 | +0.0322 | −0.0233 |
| 11 | Yuzhu Shi | 0.1888 | 0.1312 | +0.0576 | +0.0324 | −3.205 | +0.0192 | −0.0179 |
| 12 | GZ YAI | 0.1868 | NEW | | | | | |
| 13 | huhuhu | 0.1865 | NEW | | | | | |
| 14 | Leo | 0.1859 | NEW | | | | | |
| 15 | emira | 0.1838 | NEW | | | | | |
| 16 | Haystack | 0.1824 | 0.1472 | +0.0352 | +0.0204 | −1.592 | +0.0044 | −0.0009 |
| 17 | Dadhich | 0.1821 | 0.0588 | +0.1233 | +0.0982 | −13.873 | +0.0325 | −0.2145 |
| 18 | The_Way | 0.1816 | NEW | | | | | |
| 19 | cindy chen | 0.1814 | NEW | | | | | |
| 20 | jw_13288 | 0.1810 | −0.1467 | +0.3277 | +0.0801 | −100.868 | −0.0829 | −1.1033 |
| 21 | PumpkinChicken | 0.1808 | 0.1485 | +0.0322 | +0.0401 | −2.893 | −0.0105 | −0.0172 |
| 22 | Mountain ai | 0.1800 | NEW | | | | | |
| **24** | **Ahmet Dedeler** | **0.1740** | 0.1544 | +0.0196 | +0.0050 | −0.489 | −0.0117 | −0.0010 |

**Nine of the top 22 did not exist on the Aug-28 board**, and only 8 of the Aug-28 top 15 are still
on it at all. The compression did not come from the old leaders improving — it came from ~228 new
teams and from a cohort that was ranked 30–100 on Aug 28 arriving at 0.18–0.20. Our +0.0196 is
mid-pack among the Aug-28 top 15 survivors (Aginglab.com +0.0017, Cell-tics +0.0067, Weijing Zhao
+0.0140, GISL +0.0148, UM +0.0182, Jurassic Park +0.0303, AgingLab_BI −0.0024, We527 −0.0613).

### 1.4 The two families that reach 0.25, and the frontier

Every team with raw `reach` ≥ 0.24, on 2026-09-03:

| team | model | ovr | pds | mse | nmae | fid | **reach** | jac |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| cqawesome | `m0902b` | 0.2179 | 0.8684 | 0.760 | 0.9269 | 0.4929 | 0.2592 | 0.0274 |
| Aginglab.com | `GeroAI_v13` | 0.2129 | 0.8224 | 0.820 | 0.8883 | 0.5134 | 0.2585 | 0.0288 |
| Andrew Liu | `AL-2574` | 0.2083 | 0.8329 | 0.759 | 0.9297 | 0.5084 | 0.2490 | 0.0222 |
| Hongjian Zhou | `pm-r14d` | 0.2078 | 0.8325 | 0.760 | 0.9291 | 0.5084 | 0.2454 | 0.0224 |
| Slurm Shadi | `Cellulation` | 0.1918 | 0.8418 | 0.883 | 0.9391 | 0.5098 | 0.2584 | 0.0257 |
| Haystack | `GSA-020` | 0.1824 | 0.8014 | 0.829 | 0.9360 | 0.5067 | 0.2494 | 0.0230 |
| The_Way | `Ockham` | 0.1816 | 0.8204 | 6.665 | 0.9226 | 0.5192 | 0.2898 | 0.0283 |
| jw_13288 | `cedar_90` | 0.1810 | 0.8327 | 0.902 | 0.9370 | 0.5088 | 0.2455 | 0.0229 |
| UM | `BoomLi` | 0.1761 | 0.8241 | 3.999 | 0.9200 | 0.5100 | 0.2855 | 0.0245 |
| AgingLab_BI | `vomics-v11` | 0.1740 | 0.8014 | 6.622 | 0.8972 | 0.5167 | 0.2447 | 0.0327 |
| Mengran Li | `BoomAI` | 0.1736 | 0.8195 | 4.005 | 0.9188 | 0.5095 | 0.2796 | 0.0246 |
| luka | `PertMind_Pro` | 0.1732 | 0.8162 | 3.867 | 0.9222 | 0.5088 | **0.2908** | 0.0247 |
| pop | `quarry` | 0.1726 | 0.8416 | 2.661 | 0.9559 | 0.5139 | 0.2697 | 0.0267 |
| GISL | `GISL v15 …` | 0.1723 | 0.7958 | 5.009 | 0.9069 | 0.5173 | 0.2520 | 0.0351 |
| Weijing Zhao | `PAMI` | 0.1718 | 0.8194 | 4.136 | 0.9350 | 0.5108 | **0.2913** | 0.0246 |
| Cell-tics | `celltics-1` | 0.1692 | 0.8074 | 4.332 | 0.9264 | 0.5125 | 0.2653 | 0.0314 |
| *(21 more, all below 0.17 overall and all in the loud family)* | | | | | | | | |

Two families reach the same place from opposite amplitudes:

| family | n (reach ≥ 0.24) | median `jac` | median `fid` | median `reach` |
|---|---:|---:|---:|---:|
| **shrunk** (raw `mse` < 1.2) | 8 | 0.0243 | 0.5084 | 0.2539 |
| **loud** (raw `mse` ≥ 1.2) | 29 | 0.0283 | 0.5099 | 0.2557 |

**Identical median `reach` from raw `mse` 0.76 and from raw `mse` 4–7.** Whatever produces reach is
orthogonal to amplitude. Only **three** teams currently hold `reach` ≥ 0.25 at `pds` ≥ 0.84 —
`cqawesome` (0.8684/0.2592, mse 0.760), `Slurm Shadi` (0.8418/0.2584, mse 0.883) and `pop`
(0.8416/0.2697, mse 2.661); `Garrrbage`'s superseded `v13` was a fourth (0.8454/0.2492, mse 4.425)
before `v14` traded it away.

---

## 2. The centrepiece: ranked mechanisms for `reach` ≥ 0.25 at `pds` ≥ 0.84

Ranked by **evidence strength × plausible size**, not novelty. Every "how to test" line respects
audit 181: *any screen whose arm family already has a board point must retrodict that board point as
its known-answer row before its harness numbers are read* — for emission-law arms that means a21
(`resample_x 0.14`, reach −0.0185 scaled) and a22.

---

### M1 — The p-value ORDERING, decoupled from the profile. **Evidence: very strong. Size: +0.0077 overall.**

`reach` reads only (a) per-gene sign accuracy on the reference's DE pool and (b) the order our
emitted p-values impose on that pool. Three independent board arguments say (b) is where the gap is.

**(i) Six same-team arms move `reach` with the profile pinned.** Filter over all 1,232 consecutive
same-team submission pairs: `|Δpds| ≤ 0.006`, `|Δ log mse| ≤ 0.12`, `|Δreach| ≥ 0.020`, `pds > 0.70`.

| team | arm | Δpds | mse | Δreach | Δjac | Δfid | Δnmae |
|---|---|---:|---|---:|---:|---:|---:|
| BioFarmer | `CT-R4C Render` → `CT-R6 SafeMix` | −0.0045 | 4.36→3.90 | **+0.0455** | −0.0030 | +0.0156 | +0.0203 |
| **luka** | `PertMind_Pro` → `PertMind_Pro` (09-02) | +0.0031 | 3.93→3.87 | **+0.0374** | +0.0011 | +0.0018 | +0.0013 |
| tydymy | `xfer-orion-geneflow` → `regprop-collectri-v1` | −0.0014 | 3.97→4.42 | **+0.0353** | +0.0030 | +0.0133 | +0.0421 |
| **Weijing Zhao** | `PAMI` → `PAMI` (09-02) | +0.0024 | 3.93→4.14 | **+0.0324** | +0.0010 | +0.0036 | +0.0093 |
| VIVAI Inc | `vivai-m15` → `vivai-m43` | −0.0003 | 3.99→3.83 | **+0.0297** | −0.0041 | +0.0074 | +0.0121 |
| ggband | `TransMDM r05` → `r065 edit k5000` | +0.0058 | 4.38→4.03 | **+0.0253** | +0.0007 | +0.0023 | −0.0207 |

`luka`'s arm is the cleanest object on the board: `nmae` moves **+0.0013** (worse), so through
`nmae(k,q) = q|k−1| + (1−q)(k+1)` their sign accuracy `q` moved by **less than 0.001 in the wrong
direction at any fixed `k`** (`dq = −dnmae/2k` = −0.0009 at `k` = 0.75), while `pds` moved +0.003 and
`mse` −0.06. `reach` still moved **+0.0374 = +0.0070 overall.** Four of the six arms improved `reach`
while `nmae` got *worse*, which no sign-accuracy story can produce.

**(ii) The sign-accuracy channel is inverted across teams.** Inside the shrunk cluster (raw `mse` <
1.2, `pds` > 0.80, n = 35), invert `nmae` through the identity. Because `q` is monotone decreasing in
`nmae` at any fixed `k`, whoever has the lower `nmae` has the higher `q` at *every* common `k`:

| team | reach | nmae | q @k=0.75 | q @k=1.0 |
|---|---:|---:|---:|---:|
| Slurm Shadi | **0.2584** | 0.9391 | 0.5406 | 0.5305 |
| Haystack | **0.2494** | 0.9360 | 0.5426 | 0.5320 |
| cqawesome | **0.2592** | 0.9269 | 0.5487 | 0.5366 |
| Happy | 0.1871 | 0.9224 | **0.5517** | **0.5388** |
| Huiyan Qian | 0.1857 | 0.9189 | **0.5541** | **0.5406** |
| Mountain ai | 0.1814 | 0.9326 | 0.5449 | 0.5337 |
| *us* | 0.1773 | 0.9432 | 0.5379 | 0.5284 |

`Happy` and `Huiyan Qian` have strictly better `q` than `Slurm Shadi`, `Haystack`, `Andrew Liu`,
`jw_13288`, `Hongjian Zhou` and `cqawesome`, and lose 0.06–0.07 raw `reach` to all of them.
**15 of the 49 high-reach × low-reach pairs are sign-inverted on `q`.** A pure sign-accuracy
explanation needs zero. And the leader's own 09-02 jump (`m0830b` → `m0902a`, reach 0.2108 → 0.2606)
came with an `nmae` improvement of 0.0091, which the cluster's own `reach ~ nmae` slope
(−0.976/unit) prices at +0.0089 — **82% of cqawesome's +0.0498 is not the `q` channel.**

**(iii) Our deficit, priced.** Restrict to the 41 teams with `pds > 0.79` whose `fid` is inside
`[0.503, 0.512]` — i.e. `fid` is pure directional precision because their coverage term is saturated,
so `jac` differences are call-set differences and not precision differences:

```
reach ~ 1 + pds + nmae + log(mse) + jac        R2 = 0.582,  n = 41
        pds +0.944   nmae −0.367   log mse +0.0302   jac −2.82
corr(reach, jac)  = −0.569      <-- the strongest single correlate
corr(reach, mse)  = +0.431 (log)
corr(reach, pds)  = +0.302
corr(reach, nmae) = −0.218
our residual      = −0.0411 raw reach  =  −1.92 sd  =  −0.0077 overall
```

The residual ranking is itself the finding — the extremes are exactly the ordering arms of §M1(i):

| most negative | resid | | most positive | resid |
|---|---:|---|---|---:|
| black_handsome_egg `啊小绳头V12` | −0.0513 | | Haystack `GSA-020` | **+0.0432** |
| TR[AI]L `submission3` | −0.0467 | | Shaoqing Duan `MODEL` | +0.0271 |
| Mountain ai `m112` | −0.0414 | | **luka `PertMind_Pro`** | +0.0238 |
| **us** `We regret to inform the genes` | **−0.0411** | | **Weijing Zhao `PAMI`** | +0.0238 |

We are 4th-most-negative of 41 and second-worst among the current top 22 (behind `Mountain ai`). The
same fit on the shrunk cluster with only `(nmae, pds)` controlled gives −0.0369 raw (−0.0054
overall); the number grows when `jac` is controlled because we sit at a middling `jac` and still
under-reach. (`residual_rank.py`.)

**`Haystack` is the single best lead on this board.** It holds the largest ordering surplus
(+0.0432, +2.02 sd) at a `pds` of only 0.8014 — *below* ours — and a shrunk `mse` of 0.829, and it
converts that into rank 16 at 0.1824. It is a solo team (`SE`), it has never published a description,
and its whole history is `GSA-0NN` / `CDX-*` model names. Whatever it is doing to the emitted p-order
is worth +0.084 raw reach over us at a worse direction.

**What the mechanism physically is.** The scorer orders the reference's DE pool by
`(significant?, padj ↑, p ↑, |lfc| ↓, gene name)` — all four keys are functionals of *our emitted
counts*, and only two of them (the pseudobulk mean and the per-gene cell-to-cell spread) are things
we set. So the ordering key is a designable object that is **independent of the pseudobulk profile**
whenever the per-gene totals are held fixed. This is exactly 185's mandate and 177's `|H|` result
(production key 0.18283 → `|H|` key 0.23533 on h1-test); the board evidence above says the same
thing from the outside and puts a **board-measured** size on it, +0.025 to +0.046 raw, in a regime
177 could only forecast.

**How to test:** `debench` `Surface.reach(key=…)` on `h1-test` with the `|H|` known-answer row
(0.23533) and the production row (0.18283) both reproduced, donors restricted to
`k562_gwps + orion_hek293t` (h1 donor leakage, `AGENTS.md` §6.4 trap 1), pools scaled
*proportionally* per perturbation (trap 2). Then, before any slot: retrodict a21 and a22.

---

### M2 — A signed TF→target regulon prior mixed into the direction. **Evidence: strong (one board arm + one independent ablation). Size: +0.0066 overall as measured.**

`tydymy` (Complex Disorders Alliance) published their config verbatim:

> "reference store from X-Atlas/Orion HCT116 and HEK293T plus Replogle K562 genome-wide,
> decomposition transfer, template_weight 0, conserved_alpha 1.30, shrinkage 1.0, denoiser off,
> seed 0, **plus a CollecTRI signed TF to target correction at 25 percent of backbone L2**."

The arm that introduced it (`xfer-orion-genefill-k0-v4` → `regprop-collectri-v1`, 2026-08-26) moved
**`reach` 0.2196 → 0.2549 (+0.0353) with `pds` −0.0013 and `log mse` +0.108** (both `mse` values are
clamped to scaled 0, so amplitude was free), `fid` 0.4789 → 0.4922 (+0.0133), `jac` +0.0030,
`nmae` 0.9792 → 1.0213 (**worse** by 0.0421). In overall units: **`reach` +0.0066 and `fid` +0.0075
against `nmae` −0.0116**; the entry netted +0.0021 overall. **The regulon prior bought reach and
fidelity, and paid for it in `nmae` — i.e. it improved the ORDER and the CALLS without improving the
magnitudes.** That is the exact signature M1 predicts, from a completely independent team.

Independent, same-week confirmation from the literature: **D²R² (arXiv 2608.15288, 2026-08-15)**
reformulates perturbation prediction as regulation-guided gene-wise *ordered* generation, initialises
the ordering policy from **a gene regulatory network inferred from control cells**, and reports:
"controlled ablations holding the generator and generation budget fixed show that biological-prior
ordering improves over random ordering and is **more reliable than uncertainty-based heuristics**,
whereas reversing the biological-prior ordering degrades every metric." Evaluated on Norman19 and
**VCC-H1**. Two properties make this ours to use: the prior is built from **control cells only**,
which is all we have for A/B/C and will be all we have for D/E/F; and the paper's own ablation says
the *ordering* is where it pays, which is exactly M1's channel.

Third, weaker corroboration: `raft.bio` ships a "per-(target,gene) **permutation-null weight**
(permuteDE-style, two reference lines as replicates, ranked against 18.3k genome-wide knockdowns)"
and sits at `reach` 0.2300 — a confidence key built from a genome-wide null rather than from
`|effect|`.

**How to test:** build a signed TF→target regulon (CollecTRI via OmniPath, or a control-derived
co-expression GRN which is licence-clean and D²R²-shaped), form a per-(target, gene) regulon-agreement
score, and use it (a) as the `key` in `Surface.reach` — the cheapest possible read, no build — and
(b) as a small-weight additive correction to the direction, priced on PLOLO `pds` and h1-test `reach`
with a shuffled-regulon control, which is the control `tydymy` did not report.

---

### M3 — A small, precise call set. **Evidence: strong but the trade prices out near zero. Size: ≈ +0.0004 overall along the field's own line.**

`Testh` is an accidental controlled experiment. Its 2026-09-02 entry has **`fid` 0.2754** — far below
the 0.5 chance floor, which at ~0.51 precision means coverage `n_pred/n_real ≈ 0.54`, i.e. it calls
about half as many genes as the reference found. Same entry: `pds` 0.8112, `mse` 0.869, `nmae` 0.9322,
**`reach` 0.2301**, `jac` 0.0224. Against `Happy` (`pds` 0.8323, `mse` 0.837, better `nmae` 0.9224):
**+0.043 raw reach from a call set roughly a tenth the size.**

But the arithmetic on the field's own trade line is unkind. Using the board's fitted anchors
(reach span 0.8857, jac span 0.3668):

```
d reach / d jac  = −2.82        (coverage-saturated fit)
give up 0.005 raw jac  ->  −0.00227 overall
                        ->  +0.0141 raw reach = +0.00265 overall
net                                              +0.0004 overall
```

**Moving along the jac↔reach line is worth nothing.** The teams that are ahead are the ones sitting
*off* it. This matters because it is the obvious thing to try and it is a trap: it also matches 177's
CORRECTION, where `hy_kt8000`'s `jac` loss ate its `reach` gain to a wash. Treat call-set size as a
**constraint to hold fixed while M1/M2 move the ordering**, not as a lever.

---

### M4 — Emitted cell-to-cell dispersion, inside the current law. **Evidence: one board A/B, sign known. Size: unquantified, small.**

`SNU_LCDD` published a paired board result on exactly this knob:

> "Restores the orion-floor5-w025 configuration (**dispersion 0.15**). The **dispersion-0.35** variant
> measured +0.0047 on the 126-target held-out reference (p=0.027, 12 paired seeds) but scored
> **−0.0023 on the leaderboard: `pds` rose +0.019 as the gate predicted, while `nmae`, `fid` and
> `reach` fell more.**"

So on the board, **more emitted spread costs `reach`** — the Wilcoxon z per gene falls, the graded
part of the p-order compresses, and ties spread. This is the same axis 171 closed for us as a
*schedule*, and it is a knob, not a law change; the value here is the **sign, measured live**, and
that it points the opposite way from the local gate that team trusted. Two more teams measured the
adjacent object independently:

- **Zhongyun He**: replacing a multinomial redraw with "real control counts scaled by the fold change
  with stochastic rounding" — "the redraw was under-dispersed relative to real cells and produced
  **6,136 spurious significant genes under a zero-perturbation null (measured); the new construction
  gives 0.**" This is an outside replication of 177's zero-effect-null finding (a15: 3,883–4,312
  false calls; resample: 1–19).
- **Qigong Li** ran the *opposite* swap through the official scorer on an RPE1 fold — bootstrap
  resampling of real control cells → parametric gamma-Poisson with **per-gene dispersion fitted from
  that context's own NTC controls**, lognormal library sizes, baseline-aligned: "**pds +0.035,
  fidelity yield +0.048, overlap +0.025, nmae −0.012 (better), reach +0.003, mse flat**". **This is in
  tension with 177's h1-test forecast** (resampling worth +0.035–0.056 raw reach) and should be
  handed to 185 as a competing measurement: their emitter differs from 177's in that the dispersion is
  *fitted per gene from the context's own controls* rather than being a global constant, which is
  precisely the per-gene variance channel M1 needs.

**How to test:** this is 185's territory. The one thing to add from here is that the board's only
direct dispersion A/B says *lower*, and it cost that team overall because `pds` moved the other way —
so the arm must be priced on all six.

---

### M5 (negative) — Amplitude is not a `reach` lever. **Evidence: strong. Do not use it.**

29 same-team pairs with `|Δpds| ≤ 0.010` and an `mse` ratio ≥ 1.8 either way:

```
d reach = +0.0127 · d log(mse) + 0.0149     n=29, resid sd 0.0363, corr +0.358
d fid   = +0.0169 · d log(mse) + 0.0045     corr +0.520
d jac   = +0.0006 · d log(mse)              corr +0.227
```

The residual is **2.9× the fitted effect**. Halving raw `mse` predicts only −0.0087 raw reach, and the
realised outcomes straddle zero: `PumpkinChicken` −0.0009 at ×0.22, `Haystack` −0.0042 at ×0.35,
`xiaoyan chen` **+0.0337** at ×0.24, `Yuzhu Shi` **+0.0348** at ×0.26, `Zemin Li` **+0.0741** at ×0.36,
`Aginglab.com` **+0.1002** at ×0.53 — against `Garrrbage` −0.0382 and `Jurassic Park` −0.0391. This
directly contradicts `Alzheimer Project`'s published claim that "four of six members — discrimination,
direction fidelity, reach and set agreement — prefer higher amplitude"; that is true inside their own
narrow α sweep and false as a cross-team law. **Hand to 187/188: the field says you can take the whole
`mse` shrink and expect `reach` to be flat within ±0.036, with no systematic toll.**

---

### M6 (negative) — Published perturbation foundation models buy nothing here.

The board carries a family of auto-generated names of the form `(<model>+<6-char hash>)`. Three teams
swapped the named model between consecutive submissions:

| team | arm | Δpds | Δmse | Δreach | Δjac | Δnmae |
|---|---|---:|---:|---:|---:|---:|
| GZ YAI | `(STACK+x0hneb)` → `(lingshu+194pnb)` | +0.0012 | −0.001 | −0.0054 | −0.0012 | +0.0023 |
| huhuhu | `(cellflow+8c8vwp)` → `(lingshu+4zw6jv)` | +0.0018 | +0.003 | −0.0057 | −0.0011 | +0.0026 |
| Leo | `(CPA+pcqdq4)` → `(STATE+a6lh86)` | +0.0019 | +0.004 | −0.0052 | −0.0012 | +0.0028 |

All three land within 0.0015 of `pds` 0.820, 0.004 of `mse` 0.847 and 0.006 of `reach` 0.208. `lingshu`
is **Lingshu-Cell** (arXiv 2603.25240, Alibaba DAMO Academy, masked discrete diffusion over ~18,000
genes, weights in an HF collection, "leading performance on the Virtual Cell Challenge H1 genetic
perturbation benchmark"), and `STACK`/`CPA`/`cellflow`/`STATE` are the baselines it benchmarks
against. **Swapping the entire architecture across five published models moves every member by less
than 0.006; all of them land at `reach` ≈ 0.208, well short of 0.25.** The wrapper does the work.
Combined with `experiments/182-pretrained-state/RESULTS.md`, the whole "adopt a published virtual-cell
model" branch is closed by measurement, not by argument.

---

### M7 (negative, and it answers mandate target 1) — PertMind is not why `luka` has the board's best `reach`.

Read the paper (`web/pertmind_v1.txt`, arXiv 2608.16419, Zhejiang University + Tencent) and the repo.

**What it emits.** A query is a triplet `x = (c, d, g)` — cell line, **small-molecule** perturbation,
target gene — and the label space is `Y = {Up, Down, No}`. Discrete labels come from Tahoe-100M's
precomputed DESeq2 pseudobulk statistics:

```
Up   : padj <= 0.05  and  log2FC >=  0.585            (= log2 1.5)
Down : padj <= 0.05  and  log2FC <= -0.585
No   : |log2FC| < 0.20  and  padj >= 0.50  and  baseMean >= 10
       everything else is "ambiguous" and is DROPPED from training and evaluation
```

Each retained triplet carries a confidence weight `w ∈ [0.1, 1]` computed from `(Δ, q, µ)`. The policy
is Qwen3-4B Base → trusted-trajectory SFT → GRPO with gene-level + pathway-level + format rewards.
Held-out lines are C32, HepG2C3A, HOP62, Hs 766T, PANC-1.

**So the model outputs a ternary sign per gene, with a natural-language rationale, and no magnitude.**
It cannot produce a count matrix by itself. It is trained on **chemical** perturbations in cancer lines,
not CRISPRi knockdown. And covering the panel needs 300 × 18,533 × 3 = **16.7 M reasoning traces**.

**And the board says it is not what produced their score.** At full precision, two different teams
uploaded the *same file*:

```
UM   2026-08-28T05:52  `BoomLi`        avg 0.15793528158740663
luka 2026-08-28T13:48  `PertMind_Pro`  avg 0.15793528158740663
pds 0.811371237458194 | mse 3.7923639253137758 | nmae 0.941470762508308
fid 0.5080738615861684 | reach 0.24584781093365157 | jac 0.02482882862502385
```

Identical on all six members to 16 significant figures — one of only eight cross-team exact
six-member duplicates on the board that are *not* control-copy/baseline collisions. `luka`'s team
roster is
`['luka', 'BoomLi', 'unknown', 'Chixiang Lu']` (org `DenseAI`), and `UM`'s model name *is* `BoomLi`.
`Mengran Li`/`BoomAI`, `Weijing Zhao`/`PAMI` and `OpenPert`/`PertMind_Pro_*` sit in the same
fingerprint band (`pds` 0.813–0.819, `mse` 3.79–4.14) and all four made the **same ordering-only jump
on 2026-09-02**: luka +0.0374, Weijing Zhao +0.0324. **The board's best `reach` is one pipeline run
from ≥4 accounts, and its 09-02 gain is an emission/ordering change, not the LLM.**

**What *is* reusable from PertMind** is not the model, it is the labelling rule above: a
three-class, magnitude-free, `baseMean`-gated definition of "responds / does not respond" with an
explicit confidence weight from `(Δ, padj, baseMean)`. That is a ready-made specification for the
ordering key M1 needs, computable from any donor table we already hold, with an abstention class that
the current `|effect|` key does not have. Licence Apache-2.0. **This is the concrete thing to hand to
185.**

---

### M8 — The direction head and the profile head should be different objects. **Evidence: two independent competitor architectures. Size: unpriced.**

Two teams have independently split the model in exactly the shape `AGENTS.md` §4.2's operational
corollary predicts:

- `小桥细胞养殖基地` (Shenzhen U. of Advanced Technology), `STATE-RankTail V1`: "A **shared low-amplitude
  profile tail** optimizes perturbation discrimination/MSE, while an **adaptive-budget sparse head
  separately predicts BH-call probability, direction and conditional magnitude**."
- `望有幸`: "Separate continuous logFC and **discrete rank/detection heads**; complete real-control
  sparse background preserved."
- `HANYOUNG JEONG`: "Dense two-tier emission: full magnitude on the top-ranked calls, and the
  remaining gated genes **scaled below the rank-sum test's detection threshold instead of zeroed, so
  the pseudobulk profile stays intact while the DE call set stays precise**." Local H1-2025 blind
  benchmark with cell-eval2 0.16.0: **0.2213 vs 0.1546** for their prior sparse baseline.

Neither of the first two is scoring well yet (ranks 381 and 393), so this is architecture evidence,
not outcome evidence. `HANYOUNG JEONG`'s "scale below the detection threshold instead of zeroing" is
the cheapest expression of the idea and is directly implementable inside our emitter.

---

### The reach head is thin, and three teams measured how thin

Any mechanism has to survive this. Independent measurements of directional accuracy as a function of
rank/coverage, all published on the board:

| team | measurement |
|---|---|
| `HuYefan` | Replogle coverage bin ranks 97–104 retains **71.0% leave-one-target-out direction accuracy**; "ranks after 104 were rejected for rapidly collapsing reliability" |
| `Anagha Iyer` | raising moved-gene count 700 → 5000 multiplies recall ~4.9× while **conditional direction accuracy stays flat** (0.510→0.506 Orion, 0.476→0.492 H1) |
| `bhyi4` | KOLF2.1J sign-only patch on 7 covered targets: +0.0036 local, +0.0020 over sign-shuffled control |

Our own 177 measured `|H|`'s top-1% pool sign accuracy at 0.779. The field's independent numbers
(0.71 at the head, 0.51 everywhere else) bracket it. **Reach's 0.9 purity bar therefore cannot be met
by breadth at any depth — it can only be met by a head that is short and very pure, and by the
per-perturbation structure (`AGENTS.md`: 98% of our reach comes from the small-`n_real` tercile).**
Any ordering key must be evaluated *per perturbation*, never pooled (`AGENTS.md` §11.1 rule 8).

---

## 3. The 2026-09-02 `reach` wave, and the multi-account structure

**Not a scorer change** (§1.1). On 2026-09-02, `luka`, `Weijing Zhao`, `cqawesome`, `Happy`,
`xiaoyan chen`, `Yuzhu Shi`, `Zemin Li` and `Dadhich` all posted `reach` gains of +0.03 to +0.10.
For the first two the profile was pinned; for `cqawesome` the whole model improved.

**Coordinated multi-account submission, confirmed at full precision** (extends Aug-28 §1.7). Exact
six-member duplicates across *different* team names, excluding the ~40 control-copy/baseline
collisions:

| entry | teams |
|---|---|
| `pds 0.811371237458194 / reach 0.24584781093365157` | **UM (`BoomLi`) = luka (`PertMind_Pro`)**, 8 h apart |
| `pds 0.8286510591 / reach 0.2194286238` | **Vivai (`m110`) = MIT (`m110`)**, same minute |
| `pds 0.8196655518 / reach 0.1913804229` | **Kate P (`probe`) = Mike Chen (`vivai-m29`)** |
| `pds 0.6923188406 / reach 0.1495991823` | VCC@ZJU (`SilentOrbit-46`) = Sidechain (`Sidechain SER-3fn`) |
| `pds 0.7160981048 / reach 0.2344792311` | Aamir Fahira (`r4`) = QTEC (`r6`) |

Plus a 12-account `mNN` cluster (`zi Yan`, `雨瑄 孙`, `姜成程 姜`, `灿欣 卢`, `郭润楠 郭`, `xi cheng`,
`ptjam`, `Kai Song`, `zescy zhang`, `武加宁 武`, `Huang Chen`, `sixuan Chen`) all inside
`pds` 0.740–0.747 / `mse` 5.38–5.42 / `reach` 0.227–0.234, and a second `mNN` band at the top
(`Vivai` m110, `Ju` m114, `emira` m115, `minnesota university` m116, `Mountain ai` m112) with
different member names but the same naming scheme. **`Andrew Liu` (`AL-2574`) and `Hongjian Zhou`
(`pm-r14d`) submitted at the same minute, 2026-09-02T17:47, with `pds` 0.8329/0.8325 and `reach`
0.2490/0.2454** — ranks 3 and 4 are one pipeline on two accounts.

**Operational consequence:** the board's *effective* number of independent methods above 0.17 is
roughly **15, not 30**, and the ranks immediately above us are correlated. It also means a rank
change is a weaker signal than it looks.

---

## 4. Directly reusable published material, with licence

| artifact | what | licence | verdict |
|---|---|---|---|
| **PertMind** — `github.com/shapsider/PertMind`, weights `tzcfly/PertMind`, arXiv 2608.16419 | Qwen3-4B RL'd on Tahoe-100M chemical DE; emits `{Up, Down, No}` per (line, drug, gene) + rationale; `inference_vllm.py`, `inference_transformers.py` | **Apache-2.0** | **Model: no** (chemical, ternary, 16.7 M queries). **Label/confidence rule: yes** — the `(0.585 / 0.20 / padj / baseMean ≥ 10 / w ∈ [0.1,1])` specification is a drop-in ordering-key definition. Hand to 185; hand the Tahoe DE-table structure to 184. |
| **CollecTRI** signed TF→target regulons (`saezlab/CollecTRI`, via OmniPath / `decoupleR`) | 45,856 signed TF–gene interactions, 1,183 TFs, `source/target/mor` | **mixed** — aggregates 12 upstream resources; OmniPath separates academic vs commercial. **Check before use.** | The prior behind `tydymy`'s +0.0353 board reach. If the licence is a problem, D²R²'s control-derived co-expression GRN is licence-clean and does the same job. |
| **Lingshu-Cell** — arXiv 2603.25240, Alibaba DAMO, HF collection `Alibaba-DAMO-Academy/lingshu-cell` | masked discrete diffusion, ~18k genes, conditional perturbation simulation | weights on HF; **no public code found**; licence not stated on the collection page | **Skip** — three board arms say it scores the same as STATE/CPA/CellFlow (§M6). |
| **Sidechain** — `saberhq/sidechain` | a competitor's entire VCC-2026 pipeline, amino-acid-named model series (`GLY` baselines, `ALA` statistical shifts, `SER` cross-line delta transfer), plus its `CLAUDE.md` | **MIT** | Low method value (board `pds` 0.692) but a free read of an agentic competitor's workflow. |
| **PopPert** — arXiv 2609.01357, `github.com/whd1125/PopPert` | population-level joint distribution, low-rank Gaussian copula over genes, samples synthetic perturbed cells | licence unchecked | The only new paper that models the *emitted population* rather than the mean — adjacent to 185's emission-law question. |
| **D²R²** — arXiv 2608.15288 | GRN-ordered gene-wise generation; ablation: biological-prior ordering > uncertainty heuristics > random; reversing it degrades every metric; evaluated on VCC-H1 | code not released | The literature's direct statement of M1/M2. Reimplementable: the ordering prior is built from control cells. |
| **LGR** — arXiv 2608.01734 (Genentech) | LLM ranks candidate neighbour drugs, fixed mean aggregator; on Tahoe-100M, best gains on **unseen cell line**; "improves directional (sign) accuracy … even when magnitude-based metrics are similar" | — | Retrieval-quality-over-model-complexity, measured on the unseen-line axis that is our whole task. Relevant to 184. |
| Tahoe-100M precomputed DESeq2 pseudobulk DE release | what PertMind consumes; per (line, drug, dose, plate) with `log2FoldChange`, `padj`, `baseMean` | Arc Virtual Cell Atlas terms | 184's territory — but note PertMind's collapse rule (Appendix A.1) and cell-line-disjoint split are a ready protocol. |

**2025 VCC winners (mandate item):** Arc's wrap-up gave **two** $100k grand prizes; the **Generalist
Prize went to Altos Labs**, a **flow-matching generative model over the distribution of perturbed
cells** rather than over response averages, ranked with Altos' own PerturBench. Arc's own summary:
"winning approaches combined deep learning with classical statistical features … models are not yet
consistently outperforming naive baselines across all metrics." One other named 2025 winner,
`BM_xTVC`, was a deep-learning + classical-statistics hybrid. **Nothing in the 2025 write-ups
addresses the `reach`/ordering member, because the 2025 metric set did not contain it** — the
direction-reach member is new to `vcc2026`. This is a real and useful negative: the 2025 winners'
methods are not a source of the answer to this mandate's question.

---

## 5. What I could NOT determine

1. **What `luka`/`Weijing Zhao` actually changed on 2026-09-02.** The arm is visible (+0.037/+0.032
   raw reach, profile pinned) but there is no description, no commit, and no repo. I can bound it
   (`Δjac +0.001`, `Δfid +0.002` ⇒ the call-set size barely moved) but not name it.
2. **Any public trace of `cqawesome` (`Chenxi Qiu`).** No GitHub, HF, arXiv, or Discord result. Their
   descriptions are empty on all 25 submissions and their model names (`test_water_v*`, `test_fire_v1`,
   `mp4`, `E3a_ctrlprog`, `m0830b`, `m0902a/b`) leak only that `E3a_ctrlprog` used a **control
   programme**. The leader's method is dark.
3. **Whether the shrunk high-reach group (`cqawesome`, `Slurm Shadi`, `Andrew Liu`, `Hongjian Zhou`,
   `jw_13288`, `Haystack`) share one mechanism or six.** They differ enough on `nmae` and `jac` to be
   independent, and only `Slurm Shadi` and `Hongjian Zhou` have ever published a description.
4. **The challenge's `n_real`.** Nothing new. Every board-visible route to it runs through `fid`/`jac`,
   which are `n_real`-denominated (`AGENTS.md` §11.1 rule 9). My call-set inferences for `Testh`
   (coverage ≈ 0.54) and for the coverage-saturated cohort assume precision ≈ 0.51 and are only as
   good as that.
5. **Whether the `jac`↔`reach` anticorrelation is a call-set-size effect or an enrichment effect.**
   Both readings fit `corr = −0.569`; separating them needs `n_pred`, which the board does not expose.
   This is why M3 is ranked third and flagged as a near-zero trade.
6. **Lingshu-Cell's licence and whether any board team actually ran it.** The `(model+hash)` names may
   be cosmetic labels from a shared wrapper rather than real model swaps — the three swap arms are
   equally consistent with "the wrapper ignores the named model". Either way the conclusion in §M6
   holds, but the reason differs.
7. **Whether `Qigong Li`'s parametric-NB-beats-bootstrap result contradicts 177.** Their emitter fits
   dispersion **per gene from the context's own NTC controls**; 177's comparison used a global
   constant. The two may not be measuring the same axis. **Someone should resolve this before 185
   ships an emitter arm** — it is the only outside measurement that points the other way.
8. **Whether `AgingLab_BI`'s "activity-score re-ramp" (`vomics-v11`) and `Aginglab.com`'s
   `GeroAI_v13` share code.** Same organisation string, different member lists, very different
   `mse` (6.622 vs 0.820) and near-identical `reach`.
9. **The Discord.** Not readable without joining, and joining is a §14 violation. Not attempted.
10. **`STAT-DECODE`**, referenced by `Apex` as a "rank-23 approach" — no public trace found.
11. **What `Haystack` (`SE`, `GSA-020`) does.** It holds the board's largest ordering surplus
    (+2.02 sd) at a *worse* `pds` than ours, and has published nothing across 21 submissions. Its
    only visible history is the `GSA-0NN` / `CDX-*` name series and a hard shrink on 2026-08-30
    (`mse` 2.42 → 0.85 with `reach` flat at −0.0042, one of the cleanest amplitude-only arms on the
    board). No repo, no org, no description, one member.

---

## 6. Hand-offs

- **→ 185 (emission law).** M1 is your question and the board now prices it at **+0.0077 overall**
  for us specifically, measured on the only surface validated for `reach`. M7 gives you PertMind's
  ternary label + confidence-weight specification as a ready-made ordering key. M8 gives you
  `HANYOUNG JEONG`'s "scale gated genes below the rank-sum detection threshold instead of zeroing".
  **Item 7 above is a direct challenge to 177's emitter conclusion and should be resolved first.**
- **→ 187 / 188 (amplitude, magnitude shape).** M5: 29 board arms say `d reach/d log mse = +0.0127`
  with residual sd 0.0363. **The `mse` shrink does not have to be paid for in `reach`.** Six teams
  took the full shrink and gained reach. Also: `Alzheimer Project`'s public claim that reach prefers
  higher amplitude is not supported cross-team.
- **→ 184 (Tahoe).** PertMind consumes Tahoe's **precomputed DESeq2 pseudobulk DE release** directly
  and never touches the single-cell matrices; its Appendix A.1 dose/plate collapse rule and its
  cell-line-disjoint split protocol are reusable. Held-out lines: C32, HepG2C3A, HOP62, Hs 766T,
  PANC-1. Separately, `LGR` (2608.01734) benchmarks unseen-**cell-line** transfer on Tahoe and finds
  retrieval quality, not predictor complexity, is what carries sign accuracy.
- **→ coordinator.** The board's independent-method count above 0.17 is ~15, not 30 (§3). And the
  field gained +0.032 at rank 20 in six days while we gained +0.020 — the compression is real.

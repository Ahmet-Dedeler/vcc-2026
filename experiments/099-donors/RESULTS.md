# 099 — More donors: the inventory is four times bigger than we thought and the answer is still no

**Headline.** Three of the four questions close negative, and one of them closes on a structural
fact rather than a score.

1. **The donor pool is bigger than the brief assumed and smaller than `data/donors/` suggests.**
   Of the seven prepared donors, **three carry zero of the 300 challenge targets** —
   `rpe1`, `nadig_jurkat` and `k562_essential` are common-essential libraries and the VCC panel
   was built to exclude common-essential genes. `h1_2025` carries 20. Only the production three
   carry the panel. But **four usable donors are sitting outside `data/donors/`**:
   `data/external/kolf2` (**273/300**), and the three CD4 conditions in
   `data/external/cd4t/cloud-compact` (**283/287/284 of 300**), plus `perturbmatch` (251/300,
   *not* an independent donor) and `feng_ipsc` (182/300 on 6,151 genes). `data/external/orion`'s
   118 GB contains **only HCT116 and HEK293T** — both already prepared. Adding kolf2 and the
   three CD4 conditions takes the panel from *6 targets with one donor / 35 with two / 259 with
   three* to **0 / 2 / 298 at three-or-more, and 226 of 300 at seven donors.**

2. **Adding any of them as a voter in the mixture is negative.** The best cross-surface
   transferable multiplier on directional alignment from any added donor, at its own optimal
   weight, is **×1.026** (`cd4_rest` at weight 0.25, `gamma = 0`). Against ×1.478 for the
   confidence-weight family already built and ×1.5–2.1 for the whitening knob, that is nothing.
   At production's `gamma = 0.65` the same additions are **catastrophic** — `kappa` goes
   *negative* on the second surface — and §4 gives the mechanism, confirmed on two independently
   fitted operators.

3. **Adding donors makes the agreement KEY worse, not better, and the brief's premise is
   backwards.** The key is `donor_agreement`, sign **unanimity** across every carrying donor —
   an INTERSECTION. Its graded fraction on orion falls **0.206 → 0.104 (3 donors) → 0.029 (5) →
   0.015 (6)**. More donors do not resolve the key better; they resolve it to nearly nothing.
   Measured as a weight, every genuine new donor **loses** on the transferable surface (0.79–0.96×).
   The single exception is `perturbmatch` at **×1.012 transferable (1.073 / 1.101 / 1.012 across
   three surfaces)** — and PerturbMatch is a reanalysis of the *same Replogle K562 cells* as
   `k562_gwps`, so what it adds is our own donor's measurement reproducibility, not a fourth
   biological view.

4. **The one thing that survives is not a new donor. It is a knob on a donor already in
   production.** Adding the CD4 residual **after** the whitener — production's own
   `apply_frozen_cd4_residual` route — at roughly **twice** the shipped `CD4_FRACTION` reads
   ×1.077 transferable on direction with `pds` **+0.0065 / +0.0109** and `nmae` −0.0003 / +0.0024:
   every readable member non-negative on both surfaces, which no direction lever in this project
   has previously managed. **But the control kills most of the story**: adding back the base
   mixture's *own unwhitened effect* at the same fraction reproduces the direction curve almost
   exactly (×1.122 / ×1.123 against CD4's ×1.134 / ×1.077). The multiplier is the
   whitening-relaxation lever 097 already closed. What is CD4-specific is only the `pds` sign:
   at matched direction gain the CD4 residual is `pds`-positive on **both** surfaces where the
   whitening relaxation is **+0.0047 / −0.0114**, i.e. sign-flipping. The CD4-specific margin is
   **+0.0018 (transferable) to +0.0223**.

This experiment also built a **third measurement surface** — `data/harness/orion-hek293t-thinned.h5ad`,
132 × 150 cells, depth ratio 0.9999, PASS — because none of the new donors are readable on
h1-test (kolf2 covers 5 of its 85 targets, CD4 4–5). Without it every number above would have
been n = 1.

---

## 0. Measured / derived / assumed

**MEASURED** (the scorer's own primitives via `experiments/096-mse-floor/msekit.py` and
`decompose.py`; the exact analytic `de_wilcoxon_lfc_nmae` of `experiments/098-nmae/nmaekit.py`;
`pds_cosine` to the spec with panel-scope target exclusion):

- every `X`, `A`, `B`, `W = X²/A`, `kappa = X/√(B·ΣD)`, `pds`, `nmae` and per-perturbation cosine
  in this document, on `orion-hct116-thinned` (132 targets), `h1-test` (85) and the new
  `orion-hek293t-thinned` (132);
- every donor's coverage of the 300 challenge targets, of the 18,533-gene axis, and of each
  harness's own panel — counted from `perts.txt`/`genes.txt`, not inferred;
- every permutation null (60 draws for the mixture arms, 40 for the key), which shuffles the
  candidate's perturbation map only;
- the graded fraction of the unanimity and pairwise keys at 2–6 donors;
- the intermediate that explains the `gamma = 0.65` collapse (`diag1.py`): median ‖E‖ inside the
  CPM > 5 gate, before and after the operator.

**MEASURED, and it reproduces prior work bit-for-bit** — the check that licenses everything
else. `dkit.Ctx(h).dh_exact(effects(base2))` reproduces the shipped
`095-clean-harness` prediction's scored effect matrix at **median per-row cosine
1.0000000000000042 (orion, min 0.9999999987) and 1.0000000000000009 (h1-test, min
1.0000000000)**, with `W` agreeing to 8 significant figures and `pds`/`nmae` to all 7 printed
(orion 0.7953389 / 0.9308717; h1-test 0.8655462 / 0.9667390). And `keys.py` on h1-test
returns graded fraction **0.208 / 0.109** and `M_gene` **1.101 / 1.278** for the two- and
three-donor unanimity keys — the same four numbers `experiments/097-direction/pairwise3.py`
reports, with `rho` 0.287 / 0.348 against 097's 0.285 / 0.350.

**DERIVED**: the calibration between this document's `f` and production's `CD4_FRACTION` (§5.4).
`f` is defined self-containedly as the ratio of the added residual's median row norm to the
whitened base mixture's; production's `fraction` multiplies a differently-normalised object.
The mapping is arithmetic on two measured scalars plus an assumption about the builder's
amplitude bookkeeping, and it is the weakest number here.

**ASSUMED**: (i) that a ratio measured on these harnesses transfers to the board — the standing
load-bearing assumption, and the one that was 4x and 8x wrong in the flattering direction on a13
and a10; (ii) that `orion-hek293t-thinned` is a legitimate surface — it was built by the same
script, with the same policy, and passes the same depth gate as the incumbent, but it has never
predicted a board result; (iii) that the 031-somatic harness model's response resembles
production's, which additionally carries PerturbMatch and CD4.

**NOT MEASURED**: `fid`'s precision leg `k/n_pred` and `jac` for any arm — those need the full
emit + rank-sum DE path, which was deferred because the direction result did not justify it.
`reach` is not reported at all: it is unreadable on these harnesses (RESEARCH-METHOD §1) and no
re-seeded twin was run for it.

---

## 1. Inventory — what is actually usable

### 1.1 Prepared donors (`data/donors/`)

| donor | perts | genes | **challenge 300** | 18,533 axis | h1-test 85 | orion 132 | source |
|---|---:|---:|---:|---:|---:|---:|---|
| `k562_gwps` | 9,787 | 8,248 | **269** | 7,679 | 82 | 132 | Replogle K562 GWPS |
| `orion_hct116` | 17,767 | 38,606 | **290** | 18,106 | 84 | 132 | X-Atlas/Orion |
| `orion_hek293t` | 17,977 | 38,606 | **294** | 18,106 | 84 | 132 | X-Atlas/Orion |
| `h1_2025` | 255 | 18,080 | **20** | 18,077 | 85 | 12 | VCC2025 H1 hESC |
| `k562_essential` | 2,021 | 8,563 | **0** | 7,938 | 18 | 0 | Replogle K562 essential |
| `nadig_jurkat` | 2,270 | 8,881 | **0** | 8,284 | 21 | 0 | Nadig 2025 Jurkat |
| `rpe1` | 2,204 | 8,749 | **0** | 8,259 | 21 | 0 | Replogle RPE1 |

**Three of the seven are structurally unusable, not merely weak.** Verified directly:
`grep -cxF` of the 300 targets against each `perts.txt` returns exactly **0** for all three, and
`MYC` is in all three while `ACLY`/`ADNP`/`AGO1`/`TP53` are in `k562_gwps` only. These are
common-essential libraries and `docs/LOG.md` (2026-08-21) already records that the 2026 panel was
built to exclude common-essential genes. 097 recorded the same fact as "0/132 orion targets";
what is new here is that it is **0/300 on the challenge panel too**, so no weighting, gating or
representation trick can rescue them as *donors*. (Experiment 090 already rejected them as a
*covariance basis*, separately and on the official harness.)

### 1.2 Usable donors that are NOT in `data/donors/`

| object | path | perts | genes | **challenge 300** | h1-test 85 | orion 132 | hek293t 132 |
|---|---|---:|---:|---:|---:|---:|---:|
| `kolf2` | `data/external/kolf2` | 273 | 18,106 | **273** | 5 | 123 | 121 |
| `cd4_rest` | `data/external/cd4t/cloud-compact/cd4_compact/cd4_rest` | 283 | 10,113 | **283** | 5 | 129 | 127 |
| `cd4_stim8hr` | …`/cd4_stim8hr` | 287 | 10,113 | **287** | 4 | 129 | 128 |
| `cd4_stim48hr` | …`/cd4_stim48hr` | 284 | 10,113 | **284** | 4 | 126 | 126 |
| `perturbmatch` | `data/external/perturbmatch` | 9,471 | 9,761 | **251** | 79 | 123 | 121 |
| `feng_ipsc` | `data/external/feng-ipsc/genomewide_vcc_panel` | 182 | 18,533 (6,151 measured) | **182** | 3 | 80 | 79 |

All six load as `vcc2026.donor_transfer.Donor` (`perturbmatch` and `feng_ipsc` needed a
placeholder `b_ctrl`, written by `prep_donors.py`; neither is read, because `eb=False` and the
weighting is uniform, so `basal_weights()` is never called). `feng_ipsc`'s `lfc.npy` carries
**2,253,524 NaNs** on unmeasured genes — zeroed before use, recorded here because a silent NaN
propagates through `_renorm` and produces an all-`nan` arm table, which it did on the first run.

**`data/external/orion`'s 118 GB adds nothing.** `ls data/external/orion/data | sed
's/_Batch[0-9]*\.parquet//' | sort -u` returns exactly `HCT116` and `HEK293T`. `nadig2025`'s
14 GB holds `hepg2` (unprepared) and `jurkat` (prepared, 0/300); `scperturb`'s 37 GB was audited
by experiments 045 and 076; `vcc2025`'s 32 GB is `h1_2025`'s source.

### 1.3 What the extra donors do to panel coverage

Number of the 300 challenge targets carried by *k* donors:

| donor set | 1 | 2 | 3 | 4 | 5 | 6 | 7+ | **≥3** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **production 3** | 6 | 35 | 259 | — | — | — | — | **259** |
| + `kolf2` | 1 | 12 | 47 | 240 | — | — | — | **287** |
| + `kolf2` + CD4×3 | 0 | 2 | 4 | 5 | 16 | 47 | 226 | **298** |
| everything usable (10) | 0 | 1 | 2 | 2 | 10 | 19 | 266 | **299** |

The six targets only one production donor carries are `C5orf22`, `CDK8`, `GLYR1`, `IL10RB`,
`SPATA20`, `UNK`; `kolf2` or CD4 covers five of them. So the *coverage* premise of the brief is
correct and comfortably achievable. §3 is why it does not help.

---

## 2. Legality, per configuration

`orion-hct116-thinned` is a thinning of the `orion_hct116` donor; `h1-test` **is** `h1_2025`'s
own panel; the new `orion-hek293t-thinned` is drawn from the `orion_hek293t` donor. `dkit.Ctx`
refuses to construct an illegal pair — `set_donors` raises on the banned name rather than
warning, so an illegal number cannot be produced by accident.

| surface | reference | banned donor | legal base pair | new donors readable there |
|---|---|---|---|---|
| `orion` | Orion HCT116, 132 challenge targets | `orion_hct116` | `k562_gwps` + `orion_hek293t` | kolf2 123/132, CD4 126–129/132, PM 123/132, feng 80/132 |
| **`hek293t` (new)** | Orion HEK293T, 132 challenge targets | `orion_hek293t` | `k562_gwps` + `orion_hct116` | kolf2 121/132, CD4 126–128/132, PM 121/132, feng 79/132 |
| `h1-test` | H1 hESC, 85 **non-challenge** targets | `h1_2025` | `k562_gwps` + `orion_hek293t` | **kolf2 5/85, CD4 4–5/85** — unusable |

**The leave-one-out trap in this experiment is not the one 097 hit.** No configuration containing
`orion_hct116` is quoted on `orion`, and none containing `orion_hek293t` is quoted on `hek293t`.
The trap here is subtler: **h1-test cannot see any of the new donors at all**, because `kolf2` and
the CD4 compacts were both extracted for the VCC panel and h1-test's 85 targets are the 2025
panel. Every "two-surface" claim about kolf2 or CD4 in this document therefore rests on the two
Orion surfaces, and those two share an atlas, a protocol and a gene axis. That is weaker
independence than orion-vs-h1-test and it is stated wherever it binds.

### 2.1 The new surface

```
uv run python scripts/build_orion_harness.py --line HEK293T \
  --donors data/donors/k562_gwps data/donors/orion_hct116 \
  --n-perts 132 --cells-per-pert 150 --sampling greedy --thin-controls \
  --out data/harness/orion-hek293t-thinned.h5ad
```

132 perturbations × 150 cells + 18,400 controls, 38,200 × 18,533, median library 14,262 UMI.
`scripts/check_harness_depth.py`: **PASS, pert/control depth ratio 0.9999**, per-perturbation
median 1.0020 (IQR 0.956–1.044), 64/132 shallower than control — i.e. the depth confound that
invalidated seven experiments is absent by construction. Reference DE table built by the same
`ref_de.py` as 098: **125 of 132 rows have a ≥10-gene BH-0.05 gate, median gate 141 genes**
(orion 127/132; the challenge gates ~2,000 — rule 1's depth caveat applies here exactly as it
does to the incumbent surfaces, and harder).

A first build at the script's default `--cells-per-pert 400` filled only **30** perturbations and
was discarded; 150 matches the incumbent harness.

---

## 3. The prediction — adding a donor to the mixture

### 3.1 At `gamma = 0` (the operator out of the way), by weight

`base2 = {k562_gwps : 1, other Orion line : 1}` plus the candidate at weight `w`. Entries are
signed `M = sign(X)·√(W/W_base2)`; `MIN` is rule 1's transferable value.

| candidate | surface | w=0.1 | w=0.25 | w=0.5 | w=1.0 | w=2.0 |
|---|---|---:|---:|---:|---:|---:|
| `kolf2` | orion | 1.010 | **1.016** | 1.000 | 0.913 | 0.710 |
| | hek293t | 1.116 | 1.264 | 1.438 | 1.583 | 1.519 |
| | **MIN** | 1.010 | **1.016** | 1.000 | 0.913 | 0.710 |
| `cd4_rest` | orion | 1.014 | 1.026 | 1.025 | 0.970 | 0.796 |
| | hek293t | 1.026 | 1.049 | 1.062 | 0.989 | 0.712 |
| | **MIN** | 1.014 | **1.026** | 1.025 | 0.970 | 0.712 |
| `cd4_stim8hr` | **MIN** | 1.010 | **1.017** | 1.009 | 0.943 | 0.726 |
| `cd4_stim48hr` | **MIN** | **1.004** | 0.995 | 0.960 | 0.813 | 0.476 |
| `perturbmatch` | **MIN** | 0.988 | 0.962 | 0.910 | 0.806 | 0.655 |
| `feng_ipsc` | **MIN** | 0.954 | 0.864 | 0.689 | 0.344 | −0.171 |
| `h1_2025` | **MIN** | 1.011 | 1.024 | 1.039 | 1.054 | **1.063** |

**The number that kills it: ×1.026.** That is the largest transferable directional multiplier any
added donor achieves at any weight, and it belongs to `cd4_rest` at a 1/9 share of the mixture.
The confidence-weight family already built reads ×1.478; the whitening knob reads ×1.50–2.14.
A 2.6 % multiplier is not a route to `+0.0432`.

Every candidate is **single-peaked with the peak at a small weight and negative by w = 1**, which
is the signature of a donor whose vector is mostly orthogonal to the reference: a little
decorrelates, a lot dilutes.

`h1_2025` is the only monotone-rising row, and it is an artefact of coverage. It carries **12 of
132** orion targets and 13 of 132 hek293t targets, so the arm is `base2` on ~90 % of rows and the
whole effect sits on a tenth of the panel; on the challenge it would apply to **20 of 300**.

**Why `single:` rows are not the baseline.** On `orion`, `k562_gwps` alone reads M = 0.169
(κ = +0.0118) and `orion_hek293t` alone reads **1.380** (κ = +0.0968); on `hek293t`,
`k562_gwps` alone is **negative** (κ = −0.0052) and `orion_hct116` alone is 1.660. Each Orion
surface prefers its same-atlas relative, exactly as 097's `a3` found (tuned K562 shares of 0.60
vs 0.35 on opposite sides of 0.5). Any third vector that displaces K562 weight therefore gets a
free ride, which is why every arm here carries a permutation null.

### 3.2 The permutation null — the gain is information, and it is still too small

Shuffle the candidate's perturbation→row map. Coverage, gene content, magnitudes and the exact
mixture share are preserved; only the correspondence with our targets is destroyed. 60 draws.

**orion, `gamma = 0`, w = 0.25:**

| candidate | observed | null p50 | null p95 | null max | p |
|---|---:|---:|---:|---:|---:|
| `kolf2` | 1.0158 | 0.9868 | 1.0053 | 1.0134 | **0.000** |
| `cd4_rest` | 1.0255 | 0.9820 | 0.9976 | 1.0032 | **0.000** |
| `cd4_stim8hr` | 1.0173 | 0.9842 | 0.9990 | 1.0087 | **0.000** |
| `cd4_stim48hr` | 1.0243 | 0.9837 | 0.9981 | 0.9996 | **0.000** |
| `perturbmatch` | 0.9617 | 0.9734 | 1.0066 | 1.0259 | 0.700 |
| `feng_ipsc` | 0.9921 | 0.9880 | 1.0037 | 1.0152 | 0.333 |

The null sits **below 1** — a random third vector at w = 0.25 costs 1.3–1.8 % — so the dilution
worry is wrong in sign, and `kolf2`/CD4 genuinely carry target-specific direction (p = 0.000).
They carry about **+3.5 % against the null and +1.6 to +2.6 % against doing nothing.** On
`hek293t` at w = 1.0, `kolf2` reads 1.5835 against null p95 1.1816 (p = 0.000) — real, large, and
**not reproduced on orion** (0.913). `perturbmatch` and `feng_ipsc` are inside their nulls on
orion.

**This is RESEARCH-METHOD §10 exactly: which extra donor helps is a property of the reference
cell line, not of the model.** The two surfaces disagree by 1.7x on `kolf2` at w = 1.

### 3.3 Also measured, for the record

Donor `alone`, `gamma = 0`, signed `kappa`:

| donor | orion (HCT116) | hek293t | reading |
|---|---:|---:|---|
| `k562_gwps` | +0.0118 | **−0.0052** | the incumbent is near-useless on both Orion lines |
| same-atlas Orion line | +0.0968 | +0.0152 | dominates its own surface |
| `kolf2` | +0.0027 | +0.0035 | ~zero alone |
| `cd4_rest` | +0.0117 | −0.0031 | K562-like on one, negative on the other |
| `perturbmatch` | +0.0171 | +0.0044 | better than `k562_gwps` alone on both — consistent with 058 |
| `feng_ipsc` | −0.0099 | −0.0205 | **anti-aligned on both** |
| `h1_2025` | −0.0016 | −0.0152 | anti-aligned, 12–13 rows |

`feng_ipsc` should not be pursued further: 059 already returned NO-GO at a +0.000015 PDS gain,
and its alignment with both Orion references is negative.

---

## 4. Why `gamma = 0.65` destroys every in-mixture arm, and why that is structural

At production's whitening gamma the same arms read (orion / hek293t, `M`):

| arm | orion | hek293t |
|---|---:|---:|
| `base2` | 1.000 | 1.000 |
| `+ kolf2` | 0.511 | **−0.971** |
| `+ cd4_rest` | 0.212 | **−2.813** |
| `+ cd4×3` | 0.132 | **−3.402** |
| `+ perturbmatch` | 0.810 | 0.295 |

`diag1.py` prints the intermediate rather than inferring it (RESEARCH-METHOD §3). Median ‖E‖
**inside the CPM > 5 gate**, after the operator and the amplitude renormalisation:

| mixture | median ‖E‖ (all genes) | median ‖E‖ (gated) | pre-operator gated share of energy |
|---|---:|---:|---:|
| `base2` | 7.3624 | **6.5202** | 0.881 |
| `base2 + kolf2` | 7.3624 | **2.2213** | 0.831 |
| `base2 + cd4_rest` | 7.3624 | **1.7205** | 0.769 |

Before the operator the three are nearly identical (gated energy share 0.88 / 0.83 / 0.77). After
it, the added donors' energy has moved almost entirely **outside** the scored gate. The mechanism
is the operator's own design: `W = (C + λI)^(−γ/2)` amplifies low-variance directions, the
covariance `C` is estimated from the base donors' effects, and a new donor's novel directions are
by construction low-variance in that basis. They get multiplied by `λ^(−γ/2)`, the median-norm
renormalisation then divides everything by the amplified junk, and the signal shrinks 3x.

**This is not a stale-operator artefact.** The `hek293t` surface uses an operator fitted *from
scratch* on its own legal base pair (`k562_gwps` + `orion_hct116`, 21,787 basis rows,
17,942 columns, `sha256 4502bb4f…`) and the collapse is *worse* there — `kappa` turns negative.
Two independently fitted operators, same behaviour.

It also reproduces experiment 057's rejection of KOLF from the union-whitened model
("KOLF contributes a modest useful direction in raw space… that direction does not compose with
the covariance geometry learned from the two predictive sources") from a different statistic, and
experiment 090's monotone harm from widening the covariance basis. Three routes, one mechanism.

The honest caveat: production *would* refit the whitener if it adopted a fourth donor, and that
refit is not measured here for an adopted donor. It is not measured because §3.1 already prices
the in-mixture arm at ×1.026 with the operator out of the way entirely, which is the ceiling any
refit could approach.

---

## 5. The post-whitening residual — the only arm worth registering

### 5.1 The route

Production already adds an auxiliary source this way: `apply_frozen_cd4_residual` adds
`CD4_FRACTION × mean3(CD4)` to the whitened effect matrix, with every VCC panel coordinate zeroed
and targets missing any of the three conditions receiving nothing. `postwhiten.py` generalises
that to any candidate and sweeps the fraction. `cd4_mean3_gn_prod` reproduces both production
policies (**125/132 targets carry all three conditions; 300 panel coordinates zeroed**) and is
the only CD4 row that should be quoted, because it is the only one that is buildable.

`f` is defined here as **the ratio of the added residual's median row norm to the whitened base
mixture's median row norm, in one frame**. It is *not* `CD4_FRACTION`; see §5.4.

### 5.2 `cd4_mean3_gn_prod`, `gamma = 0.65`

| surface | stat | f=0.05 | f=0.10 | f=0.25 | f=0.50 | f=1.0 |
|---|---|---:|---:|---:|---:|---:|
| orion | `M` | 1.069 | **1.134** | 1.287 | 1.432 | 1.477 |
| | Δ`pds` | +0.00758 | **+0.01087** | +0.00763 | −0.00156 | −0.02099 |
| | Δ`nmae` | −0.00046 | **−0.00026** | +0.02146 | +0.02511 | +0.02497 |
| hek293t | `M` | 1.047 | **1.077** | 1.110 | 1.007 | 0.593 |
| | Δ`pds` | +0.00596 | **+0.00648** | +0.00509 | −0.00532 | −0.02331 |
| | Δ`nmae` | +0.00076 | **+0.00241** | +0.02257 | +0.02671 | +0.02728 |
| **MIN** | `M` | 1.047 | **1.077** | 1.110 | 1.007 | 0.593 |

`f = 0.25` and above hits a13's failure mode exactly — `nmae` jumps **+0.021 / +0.023** on both
surfaces while `pds` starts to fall. `nmae` is 34.6 % of the gap and is the member that cost a13
`−0.0203`; an arm that moves it +0.02 is not shippable regardless of what it does to `mse`.
**f ≈ 0.10 is the cliff edge**, and it is the operating point.

The nearest single-condition arm, `cd4_stim48hr`, is smaller and cleaner on orion (M 1.145,
Δpds **+0.01405**, Δnmae **−0.00439** at f = 0.5) but does not hold on hek293t (Δpds −0.0037 at
f = 0.5). Transferable it reads ×1.050 at f = 0.10 with Δpds +0.0038/+0.0067 and Δnmae ≈ 0.

`kolf2` post-whitening is the largest single number in the experiment on one surface —
×1.965 at f = 0.5 on hek293t with Δpds +0.0001 and Δnmae −0.006 — and it is **×1.037 with
Δpds −0.019 on orion.** Rule 1 kills it. Recorded because it is the kind of number that would
have been shipped.

### 5.3 The control that removes most of the story

Add back the base mixture's **own unwhitened effect** at the same `f`. This partially undoes the
operator, and 097 measured undoing the operator at ×1.50–2.14 on direction.

| added vector | surface | f=0.05 | f=0.10 | f=0.25 | f=0.50 |
|---|---|---:|---:|---:|---:|
| **SELF** (unwhitened base mix) | orion `M` | 1.064 | **1.122** | 1.263 | 1.433 |
| | orion Δ`pds` | −0.00434 | **−0.01139** | −0.02874 | −0.04545 |
| | hek293t `M` | 1.065 | **1.123** | 1.253 | 1.402 |
| | hek293t Δ`pds` | +0.00307 | **+0.00468** | +0.00382 | +0.00104 |
| `cd4_mean3_gn_prod` | orion `M` | 1.069 | **1.134** | 1.287 | 1.432 |
| | hek293t `M` | 1.047 | **1.077** | 1.110 | 1.007 |

**The direction curves are the same object.** At f = 0.10 the CD4 residual buys ×1.134 / ×1.077
and simply relaxing the whitener buys ×1.122 / ×1.123 — the second is *better* on the surface
where the new donors are weaker. So the `M` gain is not new donor information; it is the
whitening knob, and `WHITEN_GAMMA` reaches it directly (measured here: `gamma` 0.65 → 0.45 gives
M 1.169 on orion and 1.288 on hek293t).

**What is CD4-specific is the `pds` sign.** The whitening relaxation costs `pds` on orion
(−0.0114) and gains on hek293t (+0.0047) — it disagrees in sign, which is why 097 closed it. The
CD4 residual gains on **both** (+0.0109 / +0.0065). At matched direction gain the CD4-specific
margin is **+0.0223 (orion) / +0.0018 (hek293t)**, transferable **+0.0018**.

That is the honest size of the finding: a route to about ×1.08 of direction that, unlike the
`gamma` knob, does not charge `pds` for it on either surface.

The permutation null on the post-whitening arms (orion, f = 0.5, 60 draws) confirms the addition
is target-specific — every CD4 arm is at p = 0.000 with null p95 ≈ 0.96–0.98 against observed
1.09–1.48 — but note that this null cannot separate CD4 from the SELF mechanism, which §5.3 does.

### 5.4 The calibration, and why it is the weakest number here

Production's `CD4_FRACTION = 0.05` multiplies `EXP062.prepare_cd4_mean`'s output, in which each
condition was rescaled to `EFFECT_SCALE = 6.817` **before** averaging, and the average is *not*
rescaled back up. Measured here: that mean-of-three retains `r3 = 0.7026` of one condition's norm,
so `‖mean_effect‖ = 4.789`. The whitened base mixture it is added to has median norm
**4.6968 × ALPHA** in this frame, and the production builder does **not** renormalise to
`ALPHA × EFFECT_SCALE` after whitening (`apply_v9_whitening` preserves the operator-column norm
only). That puts production's effective `f` between **0.033 and 0.049** depending on the
challenge-panel mixture norm, which is not measured here.

**So the recommended `f = 0.10` is roughly 2–3× the shipped fraction, i.e. `CD4_FRACTION` in the
region of 0.10–0.15 — and that mapping must be re-solved inside the builder before anything is
built.**

**RESOLVED 2026-08-29, in the builder rather than on paper.** `make_submission_pm_cd4.py` now
carries `solve_cd4_fraction()`, which measures this document's `f` — the added residual's median
row norm over the whitened base's median row norm, full gene axis, medians over every row,
matching `postwhiten.py` exactly — and reports it in the build meta under
`effects.cd4.residual_frame` for **every** build. Because the residual is linear in the
fraction, that single measurement pins the whole line, so `VCC_CD4_TARGET_F` solves the fraction
in closed form instead of guessing at the 0.10–0.15 range above. Default unset: `CD4_FRACTION`
is used untouched and every prior build reproduces bit-for-bit. **Build a15 with
`VCC_CD4_TARGET_F=0.10`.** The range in this section is now a prediction the next build will
check, not a value to hand-pick from.

**MEASURED 2026-08-30, on the real challenge panel, and §5.4's arithmetic was wrong in BOTH
directions.** From `a15-cd4f010`'s build meta (`effects.cd4.residual_frame`):

| quantity | §5.4 predicted | **measured** |
|---|---:|---:|
| `f` per unit fraction | 0.66-0.98 | **1.0557** |
| effective `f` at the shipped `CD4_FRACTION = 0.05` | 0.033-0.049 | **0.0528** |
| `CD4_FRACTION` needed for `f = 0.10` | 0.10-0.15 | **0.09473** |

Both misses are the same miss, and it is instructive. §5.4 reasoned that the mean-of-three
retains `r3 = 0.7026` of one condition's norm against a whitened base of `4.6968 x ALPHA`; the
builder's own frame reads the residual at **4.4647** against a base of **4.2293**, so the base
is ~10% smaller than assumed and the ratio correspondingly larger. The shipped 0.05 was
therefore already slightly **above** §5.4's whole range, and the prescribed 0.10-0.15 would have
put `f` at **0.106-0.158** — the entire range past the cliff edge, with its top end most of the
way to the `f = 0.25` point where `nmae` jumps +0.021 on both surfaces and reproduces a13's
failure.

Worth stating plainly, because it cuts against the fix as much as for it: the naive
`CD4_FRACTION` 0.05 -> 0.10 doubling that `CURRENT_STATE.md` had specified would have landed
`f = 0.1056`, only 5.6% over target and **closer than anything in §5.4's own recommended
range**. It would have been roughly right by luck. The measurement is still what makes it
knowable rather than lucky, and 5.6% over a documented cliff edge is the wrong side to be on —
but the doc was not as wrong as the derivation that corrected it.

`a15-cd4f010` ships at the solved `0.09473`, realising `f = 0.09997`. Per RESEARCH-METHOD, transfer the *effect*, not the setting: the quantity to match is
**median ‖f·residual‖ / median ‖whitened base‖ ≈ 0.10**, measured in the builder's own frame. Two
codebases with a knob called `fraction` is the exact configuration that cost this project a
5x error on `gamma`.

---

## 6. The key — separated from the prediction, and worse with more donors

The prediction `dh` is frozen at each surface's legal `base2` build at `gamma = 0.65`; only the
confidence key changes. `M` is `kappa(z·α·dh)/kappa(dh)` with production's forms
(`α = rank01(unanimity)`, `z = max(rank01(donor_cos), 0.3)`), `rho` is the Spearman of the
per-perturbation key against the true per-perturbation cosine.

### 6.1 Graded fraction — the brief's premise, inverted

Unanimity is an **intersection** over donor pairs. Median share of the 18,533 panel genes graded
non-zero, on orion:

| donors | 2 | 3 (+kolf2) | 3 (+cd4) | 5 (+cd4×3) | 6 (+kolf2+cd4×3) |
|---|---:|---:|---:|---:|---:|
| unanimity | 0.206 | 0.104 | 0.091 | **0.029** | **0.015** |
| pairwise | 0.206 | 0.693 | 0.477 | 0.571 | 0.786 |

The brief's hypothesis was that "agreement over 5–7 donors is a far better-resolved statistic than
over 3". **It is a far worse-resolved statistic**: at six donors the production key grades 1.5 %
of genes. The pairwise key goes the other way and scores worse everywhere (§6.2), reproducing
097's unanimity-beats-pairwise result at 3 donors and extending it to 6.

### 6.2 As a weight — every genuine new donor loses

`M` relative to the two-donor key, per surface:

| added to the key | orion | hek293t | h1-test | **MIN** |
|---|---:|---:|---:|---:|
| `orion_hct116` (production's own third) | illegal | illegal | 1.046 | 1.046 (n = 1) |
| `perturbmatch` | 1.073 | 1.101 | 1.012 | **1.012** |
| `kolf2` | 0.958 | 1.097 | 1.005 | 0.958 |
| `cd4_stim8hr` | 0.964 | 1.016 | 1.004 | 0.964 |
| `cd4_stim48hr` | 0.981 | 0.935 | 1.004 | 0.935 |
| `cd4_rest` | 0.955 | 0.919 | 1.004 | 0.919 |
| `cd4×3` | 0.918 | 0.800 | 1.010 | 0.800 |
| `kolf2 + cd4×3` | 0.882 | 0.791 | 1.005 | **0.791** |

(The h1-test column is ~1.00 for kolf2/CD4 by construction — they cover 4–5 of its 85 targets, so
the key is unchanged. It is shown to make that visible, not as evidence.)

**Six donors in the key costs 21 % of the weight's value on the worse surface.** Every added
donor beats its own permutation null on `M` at p = 0.000 — a *random* fourth donor costs far more
(null p50 0.78–0.89 against base 1.000) — so the donors do carry information; it is just less
information than the graded fraction they destroy.

### 6.3 The one positive, and why it is not a fourth donor

`perturbmatch` is the only addition positive on all three surfaces (1.073 / 1.101 / 1.012), and it
lifts `rho` sharply on two:

| key | orion | hek293t | h1-test |
|---|---:|---:|---:|
| 2-donor | 0.295 | 0.367 | 0.287 |
| + `orion_hct116` | — | — | 0.348 |
| **+ `perturbmatch`** | **0.428** | 0.348 | **0.481** |
| + `kolf2` | 0.287 | 0.371 | 0.298 |
| + `cd4×3` | 0.216 | 0.097 | 0.285 |

Against the permutation null, `perturbmatch`'s `rho` gain is p = 0.000 on orion and h1-test and
p = 0.100 on hek293t (where it does not reproduce); every kolf2/CD4 `rho` is n.s. (p = 0.075–0.55).

**PerturbMatch is a reanalysis of the same Replogle K562 cells as `k562_gwps`** (its own
`meta.json`: "Replogle et al. (2022) K562 genome-wide CRISPRi Perturb-seq reanalysis by Yeung et
al. (2026)", limma-voom on PerturbMatch-matched controls). Sign unanimity between the two is
therefore **measurement reproducibility of one donor, not a fourth biological view** — which is
precisely why it grades 0.144–0.145 rather than the 0.09–0.11 a genuinely independent third donor
grades. Read charitably it is a useful denoiser for our noisiest donor. Read sceptically, and
this is the reading the evidence supports better: **production's prediction already blends
PerturbMatch 50:50 into the K562 slot**, so a key that votes PerturbMatch against `k562_gwps` is
scoring a component of our own prediction against itself — RESEARCH-METHOD §10's self-referential
signal, which has failed every previous time it was tried.

Either way the transferable gain is **+1.2 %**, and the shrinkage family that consumes this key
is closed: a13 shipped it and lost `−0.0203`. `nmae` under every key here moves the wrong way
(orion 0.93087 → 0.93438 at two donors, 0.93543 with PerturbMatch), which is the mechanism that
killed a13.

### 6.4 Answer to "prediction, key, or both"

**Neither, for a genuine new donor.** The prediction gains at most ×1.026 transferable; the key
*loses* 4–21 %. They are not separable in the flattering direction — the same donors that fail as
predictors also fail as voters, and for a related reason: a donor whose effect is largely
orthogonal to the reference both dilutes the mixture and vetoes genes in the intersection.

---

## 7. Re-opening exp 010's rejection of non-uniform weights, honestly

Experiment 010 rejected basal-similarity weighting on macro PDS over a 3-line LODO
(basal5 0.7243 against uniform's 0.7412 at `gene_norm`, monotonically worse with sharper
weighting) and chose uniform. 097 re-tested it on `kappa` and found `basal` a near-no-op by
construction (weights 0.531/0.469 and 0.515/0.485, 1.007 / 0.990).

**Re-opened here on a third objective — `nmae`'s direction leg and the deconfounded cosine, on
the union axis, with PerturbMatch and CD4 in the world — and uniform survives again**, with a new
and sharper reason:

- On `orion` the tuned second-voter is anything but `k562_gwps` (SWAP: kolf2 1.030, cd4_rest
  1.108, feng 1.133, h1_2025 1.442) and on `hek293t` it is also anything but `k562_gwps`
  (kolf2 2.087, perturbmatch 2.128, h1_2025 1.851). Both surfaces agree that K562 is the weak
  voter — but they disagree about *which* replacement is best, and every replacement costs
  `pds` 0.03–0.07 raw.
- The optimum weight for every added donor differs by more than 4× between the two Orion surfaces
  (`kolf2`: 0.25 on orion, 1.0 on hek293t; `perturbmatch`: 0.0 on orion, 1.0 on hek293t).

So exp 010's conclusion holds under the new objective, and 097's `a3` explanation (the tuned share
is a property of the reference cell line) is confirmed on a third, independently built surface.
**This is a lever closed three times on three objectives. It should not be re-opened again without
a fourth surface or a mechanism.**

---

## 8. Pre-registered prediction — `a15`, CD4 fraction

**Arm.** `a14` (or whatever the then-current calibration arm is) with the post-whitening CD4
residual raised from `CD4_FRACTION = 0.05` to the value that makes
`median‖residual‖ / median‖whitened base‖ = 0.10` **measured in the builder's own frame**.
Expected to be `CD4_FRACTION ≈ 0.10–0.15`. **One variable.** Nothing else changes: donor set
stays the production three, `WHITEN_GAMMA` stays 0.65, `VCC_ALPHA` unchanged, no shrinkage.

**Mechanism claimed.** The CD4 residual restores directional components the whitening operator
suppresses, and — unlike relaxing `WHITEN_GAMMA`, which reaches the same components — it does so
without the hubness cost that `pds` charges. Measured: ×1.077 transferable on `kappa` with
Δ`pds` +0.0065 / +0.0109 and Δ`nmae` −0.0003 / +0.0024, against the whitening relaxation's
Δ`pds` −0.0114 / +0.0047 at the same direction gain.

**Predictions vs the baseline arm.**

| member | predicted Δ | reasoning |
| --- | ---: | --- |
| `mse` raw | **−0.017** (score_mse +0.018) | ×1.077 on `kappa`, scaled through 097's 1.478x → −0.074 mapping |
| `pds` | **+0.002 to +0.004** | raw +0.0065 (worse surface) discounted for the harness/board `pds` gap |
| `nmae` | **0.000 ± 0.002** | measured −0.0003 / +0.0024 at f = 0.10 |
| `fid`, `reach`, `jac` | ~0 | no call-set policy change; `n_pred` should move < 2 % |
| **overall** | **≈ +0.004** | and I expect this to be optimistic |

**I expect this to under-deliver.** Two of this project's last three harness-based predictions
were 4× and 8× wrong in the flattering direction, both times because a member nobody was pricing
moved. The unpriced members here are `fid` and `reach`, and they were not measured for this arm.

**Falsifiers, in order of what each would retract:**

1. **`nmae` rises more than +0.004 raw.** That is the a13 signature and it means the f = 0.25
   cliff sits lower on the challenge's ~2,000-gene gate than on the harnesses' 89–141. Retracts
   the arm entirely and closes the CD4-fraction lever.
2. **`score_pds` falls at all.** The whole reason to prefer this over `WHITEN_GAMMA` is that it
   is `pds`-positive on both surfaces. If `pds` falls, the arm is the whitening knob in disguise
   and 097's closure of that knob applies unchanged — **stop, and do not try `f` above 0.10.**
3. **`score_mse` improves by less than +0.008** (under half the prediction). That would mean the
   `kappa` ratio does not transfer from these harnesses at the level assumed, which invalidates
   §5 and, more importantly, §3's negative would then also be unreadable — the negative and the
   positive share the assumption.
4. **`n_pred` moves more than 5 %.** Would mean the residual is changing the call set, which is
   not modelled here at all, and `fid`/`reach` would then be unpredicted rather than ~0.

**What would make me not ship it.** If a slot is scarce, this arm is worth about a tenth of the
gap to first and its mechanism is half-shared with a closed lever. The `nmae` measurement at the
challenge's own gate depth (LOG 2026-08-29 §4: our harnesses read the amplitude slope with the
*wrong sign* at BH 0.05) is the specific thing I would want before spending a slot.

---

## 9. What this closes

- **`rpe1`, `nadig_jurkat`, `k562_essential` as donors: closed on data, permanently.** 0/300
  challenge targets, verified directly. Not "weak" — inapplicable.
- **`feng_ipsc`: closed.** Anti-aligned with both Orion references (κ −0.0099 / −0.0205),
  negative at every weight on hek293t, and 059 already returned NO-GO.
- **Any donor added to the mixture: closed at ×1.026 transferable**, with a permutation null
  showing the information is real and too small, and a mechanism for why `gamma = 0.65` makes it
  much worse.
- **More donors in the agreement key: closed, and the premise inverted** — unanimity's graded
  fraction falls 0.206 → 0.015 from two to six donors and the weight loses up to 21 %.
- **Non-uniform donor weights: closed a third time**, now on the direction objective and on a
  third surface.

## 10. What this opens

- **The CD4 fraction** (§8). One knob, already in production, measured positive on all three
  readable members on both Orion surfaces at `f ≈ 0.10`.
- **`perturbmatch` as a fourth voter in the confidence key** (§6.3) — +1.2 % transferable, and
  the family that consumes it is closed, so this is a note rather than a lead.
- **`data/harness/orion-hek293t-thinned.h5ad`** is now available as a third surface and is the
  only one on which a configuration containing `orion_hct116` can be scored. Every future donor,
  operator or key question should be run on all three.
- **Unmeasured and load-bearing:** `fid`'s precision leg `k/n_pred` for the CD4-fraction arm. It
  needs the emit + rank-sum path, and it is the member most likely to contain the surprise.

## Files

| file | what it does |
|---|---|
| `dkit.py` | harness contexts, frozen operator/`L`, exact `dh`, `kappa`/`pds`/`nmae`, the key builders |
| `inventory.py`, `cover_harness.py`, `cover2.py` | §1's coverage tables |
| `prep_donors.py` | materialises `perturbmatch` and `feng_ipsc` as `Donor` directories |
| `arms.py` | §3.1/§4 — donor sets at uniform weight, both gammas |
| `sweep.py` | §3.1 — the ADD weight sweep and the SWAP design |
| `nullsweep.py` | §3.2/§3.3 — permutation nulls and donor-alone alignment |
| `keys.py`, `keys2.py` | §6 — the key with the prediction frozen, unanimity vs pairwise, nulls |
| `postwhiten.py` | §5.2 — the post-whitening residual sweep, incl. the production-shaped CD4 |
| `nullpost.py` | §5.3 — permutation null for the post-whitening arms |
| `ctrl_unwhiten.py` | **§5.3's SELF control — the one that removes most of the CD4 story** |
| `diag1.py` | §4's intermediate: gated energy before and after the operator |
| `ref_de.py` | 098's reference-DE cache, extended to the new harness |
| `work/*.json`, `work/final_tables.txt` | every number above |

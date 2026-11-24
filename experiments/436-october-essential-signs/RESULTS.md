# 436: signing the MDE add by the essential screens works on H1 (covered nmae −0.0118, leak-free with a capped dose). In October it is worth about zero, because Arc's panels have no targets these screens cover. Do not build it by default; trigger it only if the real D/E/F panel has ≥ 20 covered targets

Status: **COMPLETE** (2026-09-23 11:30Z, worker, Opus 5.5). Resumed after the 10:39Z reboot from the files on disk. Every
run whose stats/metrics were missing was re-run, and every known answer was re-verified from the stored files before any
number was read.
- H1 emission only, strictly serial behind the memory gate. Peak RSS was 6.95–7.72 GB per run
  (`/usr/bin/time -l`, `work/drive436-q*.log`), and `memory_pressure` read 71–91 % free (free+inactive 25–34 GiB) at every gate; no run was ever held.
- No git, no `vcc`, no `cell-eval2`, no builder run, no edit outside this directory.
- Production builder `2c90d9e1…` and emitter `edb5b8dc…` were checked unchanged at the start.

## 0. The answer

1. **On H1 the lever is real and clean with a capped dose.**
   - Covered targets: signing 426h's coarse MDE add by the mean of the covering donor tables buys **nmae −0.01176 ± 0.00016**
     on the 20 gated essential-covered targets (3 seeds).
   - The 38 uncovered gated targets are **bit-identical**: all 64 uncovered blocks' stats arrays are equal on every seed.
   - pds and mse are exactly 0.
   - H1 all-target: nmae −0.00405, fid +0.00077, jac +0.00006, reach −0.00011 → **+0.00155 ± 0.00004 overall** (435's
     board slopes).
2. **The plain dose (the one-factor balance on the new signs) leaks calls. The capped dose does not.**
   - Plain: in 6 of 21 covered blocks the donor signs move enough add from the down side to the up side that the down
     side's factor rises from 0.77–0.86 to 0.92–1.0. Those 6 blocks leak 200–1,850 calls each, **~6,860 extra calls per
     seed** in total (328 unsuppressible per covered block against 11 with the option off). That is 435's closed failure
     mode.
   - Capped: each side's factor is held at what the plain MDE gives it with our signs. Extra calls fall to **~190 per seed**,
     unsuppressible to **9.9 per block** (off: 11.0), and it keeps **95 %** of the covered nmae gain.
   - No refusals and no capacity failures in either dose: `mde_lambda` = 1 in all 21 covered blocks on every seed. Off, 3
     blocks had λ < 1; the donor signs relieve class pressure rather than add to it. 426b's crash does not recur.
3. **October: expected gain ≈ 0.**
   - The mandate's ~30 covered D/E/F targets per screen come from 199's synthetic bundles. Those draws carry **30–34
     K562-essential genes**, and Arc's real panel carries **0** (FINDINGS §8.1g, 238).
   - A/B/C also carries **0 of the 368 union targets that are off the K562-essential list**. A random essential-free draw
     carries ≤ 7.1 of them (P(0) ≈ 8·10⁻⁴), and 238's two realised essential-free draws carried 9 each. So Arc's filter
     also removes genes essential in RPE1/Jurkat/HepG2.
   - Point estimate: 0 covered targets. 95 % upper bound: 3.
   - Expected overall gain on D/E/F: **+0.00000** (point estimate). At the 95 % bound it is **≤ +0.00005** (all members),
     or +0.00001 at H1's off-list rate.
   - It would be +0.00015 on 238's essential-free draws (9 covered), and +0.0006–0.0008 only on a panel that breaks Arc's
     criterion.
   - These are H1-priced; the board's DE transfer is 0.5–1.0× of H1.
4. **Recommendation: do not ship it in the October stack by default, and do not build it now.**
   - The design in §5 is complete and tested on H1. On Oct 22, count the real panel's union coverage (one command, §4.3).
   - **If ≥ 20 targets are covered, build it from §5 and ship it.** 20 covered targets are worth ≥ +0.0003 H1-priced.
   - If fewer, skip it: at most +0.0003 before the board's DE discount, and most likely +0.

## 1. Known answers (all bit-equal on the stats arrays; `ka436.py`, `work/ka436.json`)

| row | reference | result |
|---|---|---|
| KA1 s0 / s1 / s2: 436's emitter copy (`rc_emit_k436.py` = `edb5b8dc` + the `mde_sign` path + pure diagnostics), option **off** | 435's `base_s{S}_bps_prodem_mde0.1cap0.1min0.3coarse` (production emitter) | **PASS ×3** |
| KA2 s0: option **on** with an empty target set (all-zero `mde_sign` vector passed on every block) | same | **PASS** |
| KA2 s1: same, seed 1 | same | **PASS** |
| KA3 s0: the emitter after the `capped` edit, plain dose, re-run of `_dsall` | the pre-edit emitter's `_dsall` output | **PASS** |
| KA-434g (analytic, `ana436.py`): 434's `ess434.py` arithmetic on 434's cell table | MDE our sign −0.0075, six-donor-mean sign −0.0242, essential-mean sign −0.0272 | −0.00746 / −0.02417 / −0.02719 (asserted) |

The 436 path is inert when unused: no key means an identical block, and an all-zero vector also means an identical block.

## 2. The measurement: H1, 3 seeds, raw deltas against the same-seed option-off run

Harness: `run436.py` (435's `run435.py` + `--mde-donor-sign/--mde-donor-dose`), `base S --bulk ps --mde-c 0.1 --mde-cap 0.1
--mde-min-count 0.3 --mde-coarse`. Covered = the target is in any of the four essential screens (K562-essential, RPE1, Nadig
Jurkat, Nadig HepG2): 21 of 85 H1 targets, 20 of the 58 with a gated nmae. Sign = sign of the mean delta over **all**
covering tables: the essential screens plus K562-GWPS / Orion HCT116 / Orion HEK293T where they cover the target
(`signs436.py`, `mode="all"`). `h1_2025` is never read (asserted). 0 (uncovered target, the target's own gene, NaN or a
mean of exactly 0) keeps the delivered sign. The donor sign redirects only the MDE add: the delivered LFC, the |z| design
and the call set are untouched. Where |delivered lfc| < add, the delivered sign flips (median 2,076 genes per covered
target).

### 2.1 Capped dose (the shippable one): `_dsallcp`, H1 overall **+0.00155 ± 0.00004**

| part | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| all targets | 0 (exact) | 0 (exact) | **−0.00405 ± 0.00006** | +0.00077 ± 0.00004 | −0.00011 ± 0.00003 | +0.00006 ± 0.00002 |
| 21 / 20 covered | 0 | · | **−0.01176 ± 0.00016** | +0.00311 ± 0.00018 | −0.00042 ± 0.00012 | +0.00024 ± 0.00009 |
| 64 / 38 uncovered | 0 | · | **0 (exact)** | 0 | 0 | 0 |

- The uncovered blocks are bit-identical (stats `sums` and `lfc`) on 64/64 per seed. They sit well inside seed noise.
- On the covered targets, nmae improves on 17 of 20. CYFIP1, EHMT2 and TSFM get worse.
- On the K562-essential list the covered gain is −0.0132 (17 targets). Off the list it is −0.0038 (3 targets: EHMT2,
  PIAS1, SET). The off-list number matters for October; see §4.2.
- The fid gain on covered targets has two sources. Delivered signs flip on called genes, and the capped add still moves
  ~+15 calls per block (`call_set_missed` falls). It is not priced in the recommendation beyond the "all members" column.
- Reach moves by −0.0004, far below H1's ~0.015 noise floor for K-class emitters.

### 2.2 Plain dose, for the record: `_dsall`, H1 overall +0.00276 ± 0.00002, **leaks calls**

| part | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| all targets | 0 | 0 | −0.00427 ± 0.00005 | +0.00224 ± 0.00003 | −0.00009 ± 0.00002 | +0.00078 ± 0.00003 |
| covered | 0 | · | −0.01237 ± 0.00016 | +0.00905 ± 0.00012 | −0.00037 ± 0.00009 | +0.00314 ± 0.00011 |
| uncovered | 0 | · | 0 (exact) | 0 | 0 | 0 |

The extra overall over the capped dose (+0.0012) is almost all fid and jac from leaked calls. 435 closed that for the
plain MDE: leaked calls are not transferable.

### 2.3 Which tables to average (H1 analytic + a held-out check)

- The mandate's rule, the mean over **all** covering tables, is what 2.1 measured.
- On H1 the analytic increment over our-sign MDE on the 20 covered targets is:
  - all −0.0149;
  - essential-only −0.0188;
  - all tables restricted to genes an essential screen measures (`allx`) −0.0145.
- The ess − all difference is **−0.0038 ± 0.0027 (1.4 se, 12/20)**. It is read on the instrument itself, so not a
  finding.
- Held out from H1 (`xcheck436.py`):
  - Truth is one essential screen at a time. The predictors use only cross-study donors: Replogle = K562-ess/RPE1/K562-GWPS,
    Nadig = Jurkat/HepG2, Orion. The genes are the truth screen's z1000 genes.
  - **Adding the production donors to the essential mean raises sign accuracy on the genes both measure in 4/4 folds**:
    +0.009 / +0.006 / +0.024 / +0.013, each ±0.001, over 1,980–2,234 targets per fold.
  - On genes only a production donor measures, the cross-study panel-only sign is 0.51–0.53 accurate, which is chance.
  - Off-list targets' donor signs are 0.64–0.66 accurate against 0.67–0.68 on-list, so they are still informative.
- The essential-only row, emitted with the capped dose (`_dsesscp`, §7), reproduces the analytic ordering: covered nmae
  −0.01511 against −0.01176. The per-target paired difference is **−0.0034 ± 0.0024 (1.4 se, ess better on 11/20)**, again
  not a finding. On the 3 off-list targets `ess` reads −0.019 against −0.004 for `all` (PIAS1 and SET carry it), with
  n = 3. `all` and `allx` differ by 0.0005 on H1.
- **Decision: keep the mandate's `all` rule, and call the table-set choice unresolved within noise.**
  - The held-out check, with thousands of targets, favours `all`.
  - The scored surface favours `ess` by 1.4 se.
  - Both rules are leak-free when capped, and the choice is second-order next to coverage (§4).
  - If the Oct 22 trigger fires, `ess` is a one-word change to `mde_sign.py`'s mode.

## 3. The capacity story (seed 0; `work/blocks436_s0.md`; seeds 1–2 agree to within a few calls)

| target | on K562-ess list | mde_lambda off / plain / capped | fu/fd off | fu/fd plain | fu/fd capped | unsuppressible off / plain / capped | delivered calls off / plain / capped | add signs reversed (block genes) |
|---|:-:|---|---|---|---|---|---|---|
| ARID3A | yes | 1 / 1 / 1 | 1.00/0.80 | 0.86/1.00 | 0.69/0.80 | 1 / **1367** / 0 | 2573 / **3959** / 2588 | 4274 (906) |
| ATF4 | yes | 1 / 1 / 1 | 1.00/0.82 | 1.00/0.72 | 1.00/0.72 | 0 / 0 / 0 | 2576 / 2587 / 2587 | 4099 (843) |
| COASY | yes | 1 / 1 / 1 | 1.00/0.78 | 1.00/0.56 | 1.00/0.56 | 0 / 0 / 0 | 2583 / 2591 / 2591 | 3721 (845) |
| CTBP2 | yes | 1 / 1 / 1 | 1.00/0.79 | 1.00/0.99 | 0.80/0.79 | 0 / **1212** / 0 | 2567 / **3794** / 2579 | 4191 (875) |
| CYCS | yes | 1 / 1 / 1 | 1.00/0.87 | 1.00/0.72 | 1.00/0.72 | 9 / 0 / 0 | 2585 / 2583 / 2583 | 3974 (773) |
| CYFIP1 | yes | 1 / 1 / 1 | 1.00/0.80 | 1.00/0.76 | 1.00/0.76 | 64 / 28 / 28 | 2664 / 2628 / 2628 | 3714 (794) |
| EHMT2 | no | 1 / 1 / 1 | 1.00/0.86 | 1.00/0.99 | 0.87/0.86 | 0 / **1104** / 0 | 2564 / **3679** / 2576 | 4000 (701) |
| HMGCS1 | yes | 0.916 / 1 / 1 | 1.00/0.88 | 1.00/0.83 | 1.00/0.83 | 68 / 74 / 74 | 2668 / 2674 / 2674 | 4180 (712) |
| HSP90B1 | yes | 1 / 1 / 1 | 1.00/0.79 | 1.00/0.66 | 1.00/0.66 | 0 / 0 / 0 | 2557 / 2582 / 2582 | 3987 (789) |
| MED19 | yes | 0.894 / 1 / 1 | 1.00/0.85 | 0.61/1.00 | 0.52/0.85 | 32 / **1038** / 85 | 2632 / **3638** / 2685 | 4227 (715) |
| NCKAP1 | yes | 1 / 1 / 1 | 1.00/0.83 | 1.00/0.92 | 0.90/0.83 | 0 / **195** / 0 | 2532 / **2737** / 2559 | 4280 (809) |
| PGAM1 | yes | 1 / 1 / 1 | 1.00/0.83 | 1.00/0.60 | 1.00/0.60 | 0 / 0 / 0 | 2548 / 2572 / 2572 | 4027 (733) |
| PHF12 | yes | 0.831 / 1 / 1 | 1.00/0.84 | 1.00/0.65 | 1.00/0.65 | 14 / 7 / 7 | 2614 / 2607 / 2607 | 3359 (619) |
| PIAS1 | no | 1 / 1 / 1 | 1.00/0.78 | 1.00/0.57 | 1.00/0.57 | 0 / 0 / 0 | 2556 / 2582 / 2582 | 3793 (780) |
| PSME1 | yes | 1 / 1 / 1 | 1.00/0.77 | 0.96/1.00 | 0.74/0.77 | 11 / **1842** / 10 | 2611 / **4442** / 2610 | 3588 (859) |
| RRM2 | yes | 1 / 1 / 1 | 1.00/0.78 | 1.00/0.76 | 1.00/0.76 | 0 / 0 / 0 | 2560 / 2575 / 2575 | 4119 (699) |
| SBDS | yes | 1 / 1 / 1 | 1.00/0.79 | 1.00/0.68 | 1.00/0.68 | 0 / 0 / 0 | 2576 / 2590 / 2590 | 3956 (846) |
| SET | no | 1 / 1 / 1 | 1.00/0.75 | 1.00/0.75 | 1.00/0.75 | 0 / 0 / 0 | 2582 / 2588 / 2588 | 4071 (892) |
| SETDB1 | yes | 1 / 1 / 1 | 1.00/0.78 | 1.00/0.79 | 1.00/0.78 | 0 / 0 / 0 | 2589 / 2591 / 2590 | 3649 (736) |
| THAP11 | yes | 1 / 1 / 1 | 1.00/0.81 | 1.00/0.72 | 1.00/0.72 | 32 / 3 / 3 | 2631 / 2601 / 2601 | 3415 (740) |
| TSFM | yes | 1 / 1 / 1 | 1.00/0.80 | 1.00/0.48 | 1.00/0.48 | 0 / 0 / 0 | 2582 / 2592 / 2592 | 3959 (808) |

**Mechanism.**
- The MDE add is mass-balanced by scaling down the heavier side (the up side at `fu`, or the down side at `fd`).
- With our signs every covered block is down-heavy (fd 0.75–0.88).
- The donor sign reverses the add on ~3,400–4,300 genes per block (620–910 of them block genes).
- In 6 blocks it moves enough add to the up side that the down side is barely heavier, or lighter. ARID3A, MED19 and PSME1
  tip over to up-heavy (fu 0.61–0.96); CTBP2, EHMT2 and NCKAP1 end at parity.
- The plain rule then gives the down side **nearly its full** add: fd 0.77–0.86 becomes 0.92–1.0.
- The larger down add on low-count complement genes pushes them past the |z| design's reach: the calls are unsuppressible.
- The capped rule holds each side at its plain-MDE factor (fu ≤ fu0, fd ≤ fd0) and scales the heavier side below its cap.
  The leak disappears in 5 of the 6 blocks. MED19 keeps 85 unsuppressible calls against 32 off.
- `coarse` removed 426b's class over-fill: λ = 1 everywhere, and no "did not meet the per-cell totals". The binding
  constraint is again the Wilcoxon design, as 435's dose ladder found.
- **The rule for any future sign source on the MDE add: use the capped per-side dose.**

## 4. October sizing (`size436.py`, `work/size436.json`)

### 4.1 Coverage

| panel | K562-essential genes on the panel | K562-ess | RPE1 | Nadig Jurkat | Nadig HepG2 | **union** | union off the K562-ess list | screens per covered target (2/3/4) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **A/B/C (Arc, real)** | **0** | 0 | 0 | 0 | 0 | **0** | **0** | — |
| 199 pessimistic (`def-bundle`) | 30 | 30 | 32 | 33 | 33 | 35 | 5 | 2 / 8 / 25 |
| 199 realistic (`def-bundle-r`) | 34 | 34 | 47 | 48 | 45 | 50 | 16 | 6 / 14 / 30 |
| 238 essential-free (`ef-bundle`, efp) | 0 | 0 | 8 | 7 | 9 | 9 | 9 | 3 / 6 / 0 |
| 238 essential-free (`ef-bundle-r`, efr) | 0 | 0 | 9 | 9 | 9 | 9 | 9 | 0 / 9 / 0 |

**Covered gene counts.**
- Each essential screen measures 7,938 (K562-ess), 8,259 (RPE1), 8,284 (Jurkat) and 9,024 (HepG2) of the 18,533 challenge
  genes.
- The production donors measure 7,679 (K562-GWPS) and 18,106 (each Orion line).
- On a covered target the `all` rule signs a median **18,085** challenge genes (pessimistic bundle; range 10,964–18,086),
  because the Orion lines cover almost every target. The `ess` rule signs a median 10,952.
- On H1's gated axis, `all` signs 99.7 % of genes on a covered target and `ess` 85.4 %.

### 4.2 Why the realistic October number is 0 (the mandate's premise does not hold)

- The screens' union is 2,389 targets, and 2,021 of them are on the K562-essential list.
- A panel built to Arc's stated criterion (non-lethal knockdown with a measurable response; FINDINGS §8.1g) excludes those
  2,021 by construction. That is why A/B/C carries 0 and why 238's essential-free draws keep only the 9 off-list genes.
- **A/B/C also carries none of the 368 off-list union genes.**
  - A random 300 from 238's essential-free pools would carry ≤ 7.1 (efp pool) or ≤ 17.2 (efr pool) of them:
    P(0) = 7.9·10⁻⁴ and 3·10⁻⁸.
  - 238's two realised draws carried 9 each.
  - So Arc's filter also removes genes that are essential in RPE1, Jurkat or HepG2. That fits a lethality filter that
    reads growth in the pilot screen.
- If Arc drew all 600 targets from one pool and split them at random, P(A/B/C gets none of k covered genes) is 0.12 for
  k = 3, 0.03 for k = 5 and 0.002 for k = 9.
- **Point estimate: 0 covered D/E/F targets. Rule-of-three 95 % upper bound: 3.**
- Any covered target will be **off** the K562-essential list. On H1 those gain less (−0.0038 per target over 3 targets,
  against −0.0132 on-list).
  - Held out (`xcheck436.py`), off-list targets' donor signs are 0.64–0.66 accurate against 0.67–0.68 on-list.
  - So the true off-list rate is probably between H1's two numbers.

### 4.3 Scaling (b)'s per-covered-target gain to the panel (capped dose; overall units via 435's board slopes; H1-priced)

- Per covered target on a 300-target panel: all members **+1.67·10⁻⁵**; nmae alone **+1.08·10⁻⁵**; nmae alone at H1's
  off-list rate **+0.35·10⁻⁵**.
  - That off-list rate rests on 3 targets and depends on the sign rule: `ess` gives +1.8·10⁻⁵ there (§7). Read the last
    column as a floor, not an estimate.
- The board's DE transfer from H1 is 0.5–1.0× on top of these (430: cross-study signs retrodict 417c's board at 0.51× H1).

| scenario | covered targets | all members | nmae only | nmae only, off-list rate |
|---|---:|---:|---:|---:|
| **Arc's criterion as observed on A/B/C (point estimate)** | **0** | **+0.00000** | **+0.00000** | **+0.00000** |
| rule-of-three 95 % upper bound | 3 | +0.00005 | +0.00003 | +0.00001 |
| 238 efp / efr (essential-free, 0 essential genes) | 9 | +0.00015 | +0.00010 | +0.00003 |
| 199 pessimistic (30 essential genes: violates the criterion) | 35 | +0.00058 | +0.00038 | +0.00012 |
| the mandate's planning number | 45 | +0.00075 | +0.00049 | +0.00016 |
| 199 realistic (34 essential genes: violates the criterion) | 50 | +0.00084 | +0.00054 | +0.00017 |

**The Oct 22 check (one command, no build).** After the D/E/F controls are downloaded, run
`uv run python experiments/436-october-essential-signs/size436.py`. Its first panel row reads `data/challenge/pert_counts.csv`;
it is labelled "abc (Arc, real)", but on Oct 22 it will read the new panel. Read that row's `union` column. At least 20
covered targets is the trigger in §0.4.

## 5. Production design (text only, not applied)

**Knob.** `VCC_DE_MDE_SIGN=essmean`, the only accepted value.
- It requires `VCC_DE_MDE` with the `coarse` field, so it refuses to run without the MDE add.
- It always uses the capped dose. The plain dose is not exposed (§3).
- Unset: no code runs, no key is added to any `lfc_target`, and no config key is written.

**Where the signs are computed: at build time, from `spec.perts`, never as a precomputed per-panel artifact.**
- §8.3's Oct 22 trap applies to anything keyed by panel. A table built for A/B/C would silently cover nothing on D/E/F.
- New module `src/vcc2026/mde_sign.py` is a port of `signs436.py`: `coverage` + `donor_signs(mode="all", min_agree=1)`.
  Its sha is pinned in the builder, the way `DE_MAG_MODULE_SHA256` is.
- It is called once after `spec` loads and before the context loop. It returns `MDE_SIGN_BY_TARGET[target]`, a float64 row
  in {−1, 0, +1} on the official 18,533-gene axis, for covered targets only, plus an info dict.
- Rules:
  - a target is covered iff it is in any essential screen's `perts.txt`;
  - the sign is the sign of the NaN-mean of `delta[target, gene]` over every table (essential + production) that has the
    target and measures the gene;
  - 0 where the mean is NaN or exactly 0, and on the target's own gene.
- Memory: 7 tables × 300 rows × 18,533 float32 ≈ 156 MB. The rows are read through `mmap`.

**Tables and pins** (sha256 of `delta.npy` / `perts.txt` / `genes.txt`; full hashes in `work/table_shas.txt`). The builder
asserts all 21. `h1_2025` is refused by name and by path.

| table | path | delta | perts | genes |
|---|---|---|---|---|
| k562_essential | `data/donors/k562_essential` | `4daf73cb…` | `b1e52795…` | `7f1617ef…` |
| rpe1 | `data/donors/rpe1` | `253cc16e…` | `4377d2dc…` | `ac7c7f7b…` |
| nadig_jurkat | `data/donors/nadig_jurkat` | `1fa16024…` | `f5a5c1bb…` | `3071f436…` |
| nadig_hepg2 | `experiments/166-donor-acquisition/donors/nadig_hepg2`: clone it (`cp -c`) to `data/donors/nadig_hepg2` first | `ca4565e8…` | `cea71d4d…` | `e808cb64…` |
| k562_gwps | `data/donors/k562_gwps` | `0c942228…` | `b784821e…` | `f499ffc3…` |
| orion_hct116 | `data/donors/orion_hct116` | `72c7a01d…` | `fd2fbbbe…` | `db8ba814…` |
| orion_hek293t | `data/donors/orion_hek293t` | `d95202d9…` | `3846edf5…` | `db8ba814…` |

**Builder wiring.** In the rc block's `lfc_target` dict (`scripts/make_submission_pm_cd4.py`, next to the
`mde_add`/`mde_coarse` entries, ~line 5589), add:

```python
**({} if (_ms := _mde_sign_row(target, rc_control.keep)) is None
   else {"mde_sign": _ms, "mde_sign_scope": "all", "mde_sign_dose": "capped"}),
```

`_mde_sign_row` returns `None` when the knob is unset, when the target is uncovered, or when the row restricted to this
context's gate (`rc_control.keep`) has no nonzero entry. So an uncovered block receives exactly today's dict.

**Emitter.** Port the `mde_sign` block and the `capped` branch from `rc_emit_k436.py` into `src/vcc2026/rc_emit_k.py`:
- the lines between `sign_rec = {}` and `move_m = …`: validation, `use_s`, `sign_rec`, and the `fu0/fd0` capped
  bisection;
- the `**sign_rec` diag keys.

The 436 diagnostics (`mde_dclass_*`, `mde_slack_min`, `mde_overfill_frac`) are pure and optional. Re-pin
`DE_SIGN_EMITTER_SHA256`.

**Config.** Written only when the knob is set, inside the `de_sign` block next to `de_mde`:

```json
"de_mde_sign": {"rule": "essmean", "dose": "capped", "scope": "all",
  "sign": "sign(nanmean delta over every covering table) on essential-screen-covered targets; 0 keeps the delivered sign",
  "module_sha256": "…", "tables": {"k562_essential": {"delta_sha256": "…", "perts_sha256": "…", "genes_sha256": "…"}, "…": {}},
  "n_covered": 0, "covered_targets": [], "sign_matrix_sha256": "…"}
```

**Tests that would prove it.**
- **T1, emitter, no key.** 436's harness (`RUN426_EMITTER=prod`) with the new emitter reproduces 435's
  `base_s{0,1,2}_bps_prodem_mde0.1cap0.1min0.3coarse` stats bit-for-bit. This is KA1.
- **T2, emitter, all-zero vector on every block.** Same result, bit-for-bit. This is KA2.
- **T3, emitter, capped donor signs.** It reproduces 436's `_dsallcp` stats for seeds 0/1/2 bit-for-bit, with the vectors
  from `mde_sign.py` on H1's 85 targets. First drop `run436.py`'s assert that refuses `--mde-donor-sign` under
  `RUN426_EMITTER=prod`.
- **T4, sign module.** On H1's targets it equals `signs436.donor_signs(mode="all")` exactly (21 covered). On A/B/C it
  covers 0. On 199's `def-bundle` it covers 35 and on 238's efp 9, matching `work/size436.json`. It raises if `h1_2025`
  appears.
- **T5, builder, knob unset.** The package is byte-identical to the incumbent recipe (the standing rule).
- **T6, builder, knob set on A/B/C.**
  - 0 covered, so no block gets the key.
  - The package's content digest equals the knob-off build.
  - The config differs from knob-off by exactly one key (`de_mde_sign`, `n_covered` 0).
  - This is the no-op proof on the real panel, and it runs today.
- **T7, rehearsal.** On 199's `def-bundle` through `scripts/def_pipeline.py`:
  - covered blocks carry `mde_sign_n > 0`, and uncovered blocks are bit-identical to knob-off;
  - every covered block's `call_set_unsuppressible` stays within +60 of its knob-off value. That is the leak guard
    (436's max was +53). Make it a postbuild check.

## 6. Recommendation

- **Not in the October stack by default, and not built now.** The expected overall gain on D/E/F is **+0.00000** (point
  estimate: 0 covered targets under the criterion Arc applied to A/B/C). The 95 % upper bound is **+0.00005** (3 covered
  targets, all members, H1-priced), and the board discount of 0.5–1.0× applies on top of it.
- The cost is real:
  - an emitter re-pin in the production path;
  - a new module;
  - 21 pinned hashes;
  - a rehearsal re-run.
- **Trigger (decided now):** on Oct 22, count the real panel's union coverage (§4.3). **If ≥ 20 targets, build it from §5
  and ship it.** 20 covered targets are worth +0.0003 H1-priced on all members (+0.0002 on nmae alone), or +0.0001–0.0003
  after the board's DE discount, and there are 14 days to the deadline. If fewer, drop it.
- What does transfer from this experiment:
  - **any sign source for the MDE add must use the capped per-side dose**. The plain one-factor balance leaks
    200–1,850 calls in the blocks whose balance it shifts (6 of 21 here);
  - the donor mean over all covering tables beats the essential-only mean on held-out screens (4/4 folds);
  - panel-only donor signs are at chance cross-study (0.51–0.53).

## 7. Secondary rows (queue 5)

**KA2 s1** (empty target set, seed 1): **PASS**, bit-equal to 435's seed-1 stats.

**Essential-only sign, capped dose (`_dsesscp`), H1 overall +0.00172 ± 0.00003**:

| part | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| all targets | 0 | 0 | −0.00521 ± 0.00005 | +0.00054 ± 0.00006 | −0.00014 ± 0.00001 | +0.00002 ± 0.00001 |
| covered | 0 | · | −0.01511 ± 0.00015 | +0.00220 ± 0.00025 | −0.00055 ± 0.00003 | +0.00006 ± 0.00003 |
| uncovered | 0 | · | 0 (exact) | 0 | 0 | 0 |

- Leak-free: 11.1 unsuppressible per covered block against 11.0 off, and ~214 extra calls per seed (max 71 in one block).
- `mde_lambda` < 1 in 1 covered block on two seeds (off: 3).
- Uncovered blocks are bit-identical, 64/64 on every seed.
- Covered nmae on-list −0.0144 (17 targets), off-list −0.0193 (3). CYFIP1, EHMT2 and SETDB1 get worse.
- Per covered target on a 300-target panel: +1.77·10⁻⁵ overall (all members), +1.39·10⁻⁵ (nmae only). At its off-list
  rate, the 95 %-bound October panel (3 covered targets) is worth +0.00005. §4's conclusion does not change.

## Files

- Code:
  - `run436.py` (harness), `rc_emit_k436.py` (emitter copy; `work/rc_emit_k436.v1.py` is the pre-`capped` version, and
    KA3 proves the two equal on plain);
  - `signs436.py` (sign tables), `ka436.py` (known answers + six-member table per variant), `ana436.py` (analytic + KA-434g);
  - `xcheck436.py` (held-out table-set check), `xscreen436.py` (cross-screen reproducibility on/off list);
  - `size436.py` (October coverage), `sum436.py` (summary + scaling);
  - `drive436.sh` (serial gated driver).
- Outputs:
  - `work/ka436.json`, `work/tab436-{dsall,dsallcp,dsesscp}.json`, `work/sum436.json`;
  - `work/ana436.json`, `work/xcheck436.json`, `work/xscreen436.json`, `work/size436.json`;
  - `work/blocks436_s0.md`, `work/table_shas.txt`, `work/h1/` (stats/metrics/diag/per-target per run);
  - `work/drive436-q{1,2,4,5}.log`, `work/logs/`.
- The reboot killed queue 2 during `base 2 … --mde-coarse`, and no partial files were left. Queue 3 never started, and its
  runs were re-queued as queue 4.

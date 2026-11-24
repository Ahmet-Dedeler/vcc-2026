# 803 — the lab's DE engine: the production emitter on the H1 surfaces, on EC2

Worker [803], 2026-10-15 22:45Z → 2026-10-16. Mandate: `MANDATE.md`. Cloud: 7 launches, $10.48 of the $60 budget (actuals from `labcloud.sh wait`).

**Bottom line.** The engine runs: one command per (surface, arms, seed). It gives per-target fid / jac / reach / nmae at
f 1 / 0.61 / 0.5 / 0.35 / cmp, the scorer's own per-target pds and mse legs, K, calls and pseudobulk sums. It runs on
EC2 from inputs that are sha-checked against the Mac.
- **The real cell-eval2 0.16.0** (on the written cells of one run) reproduces the engine's f 1 members per target: nmae to
  4e-15, fid / reach / jac / pds exactly, and mse to 5e-16.
- **Linux is not bit-equal to the Mac.** The emitted counts differ by ±1 in about 0.1–0.2 % of (target, gene) sums, and
  even identical counts give p values that differ in the last bits (≤ 3e-16). Every downstream pass (zs, UPRN, DN, the
  stats) reproduces the Mac given the same counts (b6), so the only platform-dependent step is the emission itself.
- **Every known-answer pair still reproduces the Mac to ≤ 1.4e-8 raw** (G1 zs, G2 MDE, G3 cp, the hook in both forms, DN
  on UPRN), and 619's full per-seed / per-regime table to ≤ 2.8e-6. The one exception is **the UPRN pass**: per seed
  |Δnmae| is up to 2.0e-4 (3-seed mean Δ 3.6e-5, against a seed SD of 3.9e-4 for the pair itself). That fails the
  mandate's 2e-5 rule, so every UPRN-containing arm (`upr`, `updn`, `c4hupr`, `inc`) carries a platform term (§3.3).

## 1. The driver (interface)

```
uv run python experiments/803-de-engine/de803.py run --surface {h1test,h1tv} --arms A[,B,...] --seed S
uv run python experiments/803-de-engine/de803.py ka --name N (--arm surface/arm/seed | --new stats.npz) --ref stats.npz
uv run python experiments/803-de-engine/pairs803.py <out dir> [--md F --json F]   # sanity table (no board units, no bins)
uv run python experiments/803-de-engine/kac803.py <out dir> [out json]            # KA (c) vs gate618 / score619
uv run python experiments/803-de-engine/cells803.py write|compare ...             # the real-scorer cross-check
bash experiments/803-de-engine/cloud/batch803.sh <parallel> <phase file> ...      # the single EC2 job (labcloud.sh)
```
- **How `run` works.** It takes arms that share ONE emitter configuration and runs the emitter once (`emit803.py emit`).
  For h1tv that imports `618/run618e.py`. For h1test it imports `618/run618e.py` in its no-swap `h1ka` mode when there is
  no hook, and `612/run612.py` with a hook. Both sit on `496/run496.py`, and all are read-only.
  - Then `630/rnup630.py` runs once for the union of the arms' rnup variants (`emit803.py rnup`).
  - Then the members step runs for every arm.
  - Outputs are redirected into 803: run618e's NAME argument is a relative path into `803/work/scratch/…`, and run612's
    and rnup630's WORK are re-pointed after import. **Nothing outside 803 is edited.**
  - The K arms patch `run400.CALLSET` (the harness twin of the builder's `VCC_RC_CALLSET_CPM`, 523).
- **Per arm and seed:** `work/out/<surface>/<arm>/s<S>/{members.npz, members-meta.json, manifest.json}`.
- **`members.npz`** (T targets in the reference's order):
  - `fid_r, jac_r, reach_r, nmae_r, nreal_r` for r ∈ {f1, f0.61, f0.5, f0.35, cmp}: 528's `rows_at` on 521's truth pools,
    realized order, i.e. score619's functions. f 0.61 is added for call-set pairs (662).
  - `ce_nmae, ce_fid, ce_reach, ce_jac, ce_pds`: cell-eval2's own per-target values at f 1 (evaluate323.score: the
    scorer's `dispatch_de_metrics` and `discrimination_score`, panel genes out).
  - `mse_num, mse_den`: the scorer's per-target legs of `expr_mse_unbiased_capped_norm` (096's msekit over cell-eval2's
    `mse_unbiased_capped` / `distance_unbiased` with the jackknife terms). The member is Σnum / Σden.
  - `ce_raw` (the six raw members), `K` (ballast K per target), `ncall`, `sums` (T × G per-target pseudobulk, int64),
    `jk`.
  - **pds and mse are therefore computed** from the sums against the surface's reference with the scorer's own
    functions, not a surrogate.
- **Self-checks on every members step:** mean of 521's f 1 rows == cell-eval2's raw member to 1e-9, and Σnum/Σden == raw
  mse. All 94 arm-seeds passed.
- **`manifest.json`** records:
  - argv, the arm file and its sha256, the emitter argv;
  - every step's manifest: the git blob sha1 of every imported repo module, every input opened with its sha256 (files
    < 256 MB; the big inputs are sha-checked once per instance by `cloud/setup803.sh` against
    `stage/SHA256SUMS.inputs` / `.small`), python / platform, wall and peak RSS;
  - the stats sha256 and the members.npz sha256.
- **Cost of one run on EC2 (r7i, Linux):**

  | surface | emitter run | with rnup (2 variants) | members (per arm) | peak RSS |
  |---|---|---|---|---|
  | h1test | 11–14 min | 47–49 min | | 6.0 GB |
  | h1tv | 18–20 min | 67–70 min | ≈ 1 min | 6.4 GB |

  The UP pass is single-threaded and dominates. Runs are deterministic across instances: h1test g1pre / zs s0 on an
  r7i.2xlarge (b2) and an r7i.8xlarge (b5) give identical raw members.

## 2. Arms (`arms/*.json`; pairs in `arms/pairs.json`)
All emitter configurations share:
- the production emitter `rc_emit_k` 0f56209d, the CLEAN4 q80 block prior, and `--bulk ps`;
- the surface's faithful cap15 key (323's / 618's);
- h1_2025 is never a donor.

| emit | hook | MDE | zero repair | CPM block | zs | dump | arms (stage) |
|---|---|---|---|---|---|---|---|
| E0 | none | coarse 0.1/0.1/0.3 | −5 | 2,200 | yes | yes | `g1pre` (pre), `zs`, `upr` (rnup `b_k100_q20r3p10` = UPRN), `updn` (+ `_D200` = DNMAG d_m200) |
| Ec4h | CLEAN4 τ 0.05 | coarse | −5 | 2,200 | yes | yes | `c4hpre`, `c4h`, `c4hupr`, **`inc`** |
| Eid | idall3 (612 v3, ungated) | coarse | −5 | 2,200 | yes | no | `idall3pre`, `idall3` |
| Emde | none | coarse | 0 | 2,200 | yes | no | `mde` (pre; = 618 g2mde, and on h1test = 435's run), `mdezs` |
| Enomde | none | none | 0 | 2,200 | no | no | `nomde` (= 618 g2nomde) |
| EK20 / EK16 | none | coarse | −5 | 2,000 / 1,600 | yes | no | `k20`, `k16` |
| Ecp (h1tv) | 528's call-propensity block | coarse | −5 | 2,200 | yes | no | `cp` (618's G3; WORN on h1tv) |
| Ec4hn (KA) | CLEAN4 | coarse | 0 | 2,200 | yes | no | `c4hnozr` (= 612's h1_c4kgall) |

- **`inc` is the October DE incumbent:** zs + coarse MDE + CLEAN4 hook 0.05 + UPRN + DNMAG d_m200 + zero repair −5. Every
  component has a wrapper on both surfaces.
- **The harness's DNMAG is not the production module.** It is 631's `dn631.dn_pass` through rnup630's `_D<m>` suffix,
  as 618's stacks ran it. The production module is `src/vcc2026/dn_mag.py` / `dn_mag2.py`, which 631/671/693 KA'd
  against it. 803 did not re-prove that equality.

## 3. Known answers

### 3.1 (a)/(b): Linux vs the Mac's stats, seed 0 (runs b1, b3; `work/out/ka-b1/*.json`)
**No KA is bit-equal.**

| KA | new (803) | reference (Mac) | sums cells ≠ (max) | jk ≠ | p max \|Δ\| | lfc max \|Δ\| |
|---|---|---|---|---|---|---|
| (a) run618e h1ka == 435, pre-zs | h1test/mde/s0 | 435 `stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse` | 1,508 (1) | 52 / 85 | 0.34 | 7.5e-4 |
| (a) … zs | h1test/mdezs/s0 | 435 `…coarse_dnk250zs` | 1,508 (1) | 50 | 0.34 | 7.5e-4 |
| 612 hook (no zr) | h1test/c4hnozr/s0 | 612 `h1_c4kgall/…_dnk250zs` | 1,508 (1) | 48 | 0.34 | 7.5e-4 |
| (b) 618 g1 pre-zs | h1tv/g1pre/s0 | 618 `g1/…zr-5…coarse` | 3,020 (1) | 95 / 126 | 0.39 | 6.6e-4 |
| (b) 618 g1 zs | h1tv/zs/s0 | 618 `g1/…_dnk250zs` | 3,020 (1) | 89 | 0.39 | 6.6e-4 |
| (b) 618 hk_c4h zs | h1tv/c4h/s0 | 618 `hk_c4h/…_dnk250zs` | 3,020 (1) | 94 | 0.39 | 0.06 |
| 618 g2mde | h1tv/mde/s0 | 618 `g2mde/…coarse` | 3,020 (1) | 95 | 0.39 | 6.6e-4 |
| 618 t618b (630 + DN) | h1tv/updn/s0 | 618 `stacks/stats-rnup_t618b_s0_b_k100_q20r3p10_D200` | 3,020 (1) | 88 | 0.39 | 0.73 |

- **Where it comes from.** On both surfaces the emitted counts themselves differ: 1,508 / 3,020 per-gene sums are off by
  exactly one count, the same number in every arm.
- **Identical counts still give different p.** On the Mac's own h1-test blocks, the Linux stats step already gives
  different p values (b5: rnup630's own zs KA stopped at "re-applied zs p differs"; b6 measures it, §3.4).
- **Switching off NumPy's runtime SIMD dispatch did not change it** (b3). This NumPy build dispatches only
  AVX512_ICL / AVX512_SPR, and the result stayed byte-identical to the default. The source is the platform's float
  arithmetic in the emitter / statistics stack (libm / numba / BLAS), not a dispatch path we can switch off.
- **Linux runs are deterministic across instances** (§1).

### 3.2 Member level (Linux − Mac, raw, seed 0, cell-eval2's own members)
- **Absolute members:**
  - nmae ≤ 1.9e-7, fid ≤ 1.2e-6, jac ≤ 4e-7, reach 0 on every non-rnup run;
  - pds −4.2e-4 (h1test) / +6.4e-5 (h1tv), mse 1.4e-6.
  - The pds shift is real: the bulk sums differ by those ±1 counts. It cancels in every pair (Δpds between arms is
    exactly 0 or 1e-8, as on the Mac).
- **Pair deltas reproduce the Mac to ≤ 1.4e-8 raw on every member:** G1 zs on both surfaces, the hook in both forms, and
  612's hook. DN on top of UPRN reproduces to ≤ 3e-6 on all three h1tv seeds (Mac −0.001173 / −0.000965 / −0.000986;
  Linux −0.001173 / −0.000962 / −0.000985).

### 3.3 UPRN: the one pair outside the 2e-5 rule
| h1tv upr − zs (nmae, f 1) | s0 | s1 | s2 | mean | seed SD |
|---|---|---|---|---|---|
| Mac (618 stacks `t618u` vs its zsKA) | −0.010535 | −0.009980 | −0.009774 | −0.010096 | 3.9e-4 |
| Linux (803) | −0.010738 | −0.009926 | −0.009734 | −0.010133 | 5.3e-4 |
| Δ | −2.0e-4 | +5.4e-5 | +4.0e-5 | −3.6e-5 | |

- **Mechanism (b6 diagnosis, §3.4).** The pass is deterministic across platforms on identical input, so the gap is input
  sensitivity. The UP pass is a greedy, threshold-heavy allocation over ~3,700 genes per target,
  with a call-flip guard that re-runs the pass without any gene whose call flips. Its recipients and donors are chosen by
  strict float comparisons against control values (`log1p(x / s) > cmax`, `rank_mag.cell_bounds`), so ±1-count inputs
  and last-bit p differences change which genes it moves.
- **Consequence for pricing.** UPRN-containing arms carry a platform term:
  - up to 2e-4 raw nmae per seed (2 % of the pass's effect);
  - about 4e-5 on a 3-seed mean, i.e. about 1e-5 overall at INV_nmae 1.648 / 6;
  - that is 10× below the pair's own seed SD and 25× below DE re-roll noise on the board.
- **What it does not affect.** Arms that do not touch the UP pass reproduce to ≤ 1.4e-8. A comparison between two arms
  that BOTH contain UPRN (e.g. DN on UPRN) reproduced to 3e-6.
- **Rule for the checker:** price UPRN-containing arms only against bases run on the same platform, and never mix Mac and
  Linux numbers inside a UPRN pair.

### 3.4 b6: Linux on the Mac's own blocks (h1-test s0, `630/work/blocks`; `work/out/diag/`, `work/out/ka-b6/`)
The Mac's 85 dumped pre-pass blocks of 435's construction (sha-checked, `stage/SHA256SUMS.macblocks`) were run through
the Linux stats step and the Linux rnup630 (UPRN, + DN):
- **Stats on identical counts:** sums equal. p / q / lfc / jk differ only in the last bits: max |Δp| 3.3e-16, |Δq| 7.8e-16,
  |Δlfc| 2.2e-16, |Δjk| 7e-15, 0 calls changed. This last-bit difference is what made rnup630's exact zs KA stop on b5.
- **UPRN and UPRN + DN on identical blocks:** sums bit-equal to 630's Mac stats, p / q / lfc within 8e-16, and the six
  members equal to the Mac's (nmae Δ 1.1e-16, the rest 0).
- **So the UP pass and the DN pass are platform-deterministic.** The whole Linux-vs-Mac difference enters at the
  EMISSION: the ±1-count differences in 0.1–0.2 % of the per-gene sums (§3.1). Every pass downstream reproduces given
  the same input. The UPRN gap of §3.3 is the pass's own sensitivity to its input: a different realization of the same
  emission law, i.e. re-roll class (AGENTS 11.1 rule 13). It is not a platform bug in the pass.
- The emitter's platform-dependent float step (libm / numba in `rc_emit_k` or the ballast allocation) was not
  localised; doing so would need an instrumented emitter run on both machines.

### 3.5 (c): score619's f 1 == 618's gate, on the engine's Linux re-runs (`work/out/kac803.json`, b7)
Measured on the Linux runs of all three gate pairs (G1 zs, G2 coarse MDE, G3 cp; 3 seeds each) against `gate618.json`
(f 1, seed means) and against `score619.json` (every regime × member × seed):
- **worst |Δ| = 2.8e-6 raw** (G3 jac at f 0.35, one seed): PASS under 2e-5;
- G1 nmae −0.0016278004 vs −0.0016278002;
- G2 nmae −0.014731, fid −0.000057, jac −0.000104, reach +0.000500, all equal to gate618 to ≤ 1e-6;
- G3 fid −0.004076 vs −0.004075 and jac +0.008724 vs +0.008723.

The engine's member functions are 619's.

### 3.6 Item 5: the real cell-eval2 on the engine's cells (`work/out/xcheck/`, b2)
- **How:** `cells803.py write` re-applied the production zs pass (`vcc2026.rank_mag.postpass dn_k250_zs`) to the dumped
  pre-pass blocks of h1test/zs/s0. It asserted every block's sums / p / q / lfc equal to the arm's stats file, then wrote
  the 85 × 400 cells plus h1-test.h5ad's own 18,400 controls (real file sha `06c6baea…`, on its gene axis).
- `cell-eval2 run --preset vcc2026` 0.16.0, backend scanpy, 2 h 33 min on an r7i.2xlarge.
- **Result** (85 targets; 58 nmae-eligible; mean values):

| member | cell-eval2 (real scorer) | engine `ce_*` | engine lib521 f 1 | max per-target \|Δ\| |
|---|---|---|---|---|
| nmae | 0.9408535577858443 | 0.9408535577858442 | 0.9408535577858441 | 3.8e-15 |
| fid | 0.5089291277653674 | same | same | 0 |
| reach | 0.26925697759272765 | same | same | 0 |
| jac | 0.07469853420125332 | same | same | 0 |
| pds | 0.8948179271708682 | 0.8948179271708684 | — | 0 |
| mse | 0.9253597869810147 | 0.9253597869810152 (Σnum/Σden) | — | 5e-16 on the member |

- The engine's emitter-side DE statistics (evaluate323 / exact_de323) and cell-eval2's own scanpy DE give the same calls,
  lfcs and members.
- The h5ad was written from the dumped blocks with no change to the emitter: the zs pass was re-applied, which needed
  only a wrapper.

## 4. Pair deltas (sanity table; full table `work/out/pairs803.md`, json `work/out/pairs803.json`)
3 seeds, raw Δ (arm − base), mean ± seed SD, Linux. Selected f 1 / f 0.35 columns beside 618 / 619 / 611 and the board:

| pair | surface | Δnmae f 1 | Δnmae f 0.35 | other | 618 / 619 / 611 (Mac) | board raw |
|---|---|---|---|---|---|---|
| G1 zs (zs − g1pre) | h1tv | −0.00163 ± 0.00002 | −0.00186 ± 0.00006 | — | −0.00163 / −0.00186 (619) | nmae −0.00218 |
| G1 zs | h1test | −0.00183 ± 0.00007 | −0.00185 ± 0.00003 | — | −0.00183 / −0.00185 (611) | |
| G2 coarse MDE (mde − nomde) | h1tv | −0.01473 ± 0.00008 | −0.01610 ± 0.00022 | reach f1 +0.0005, f0.35 −0.0016 | −0.01473 / −0.01610 (619) | nmae −0.01679, reach +0.0024 |
| G2 coarse MDE | h1test | −0.00953 ± 0.00021 | −0.01045 ± 0.00013 | reach f1 −0.0036 | −0.00966 / −0.01059 (611, other runs) | |
| G3 cp (cp − zs), WORN | h1tv | +0.00106 | −0.00001 | fid −0.00408 / −0.00468, jac +0.00872 / +0.00702 | fid −0.00408, jac +0.00872 (618) | fid −0.00358, jac +0.00399 |
| hook, board form (c4hpre − g1pre) | h1tv | −0.00226 ± 0.00004 | −0.00136 ± 0.00004 | reach f1 −0.0004 | | nmae −0.00158 (496a) |
| hook, 618 form (c4h − zs) | h1tv | −0.00227 ± 0.00002 | −0.00139 ± 0.00008 | reach f0.35 +0.0006 | −0.00227 / −0.00139 (618) | |
| hook, 618 form | h1test | −0.00146 ± 0.00003 | −0.00213 ± 0.00003 | reach f1 −0.0011 ± 0.0042 | (496: −0.00146 on the plain base) | |
| UPRN (upr − zs) | h1tv | −0.01013 ± 0.00053 | −0.00982 ± 0.00032 | other members 0 | Mac f1 −0.01010 (§3.3) | nmae −0.00305 (630arm) |
| UPRN | h1test | −0.00813 ± 0.00034 | −0.01048 ± 0.00066 | | | |
| DN d_m200 (updn − upr) | h1tv | −0.00104 ± 0.00012 | −0.00491 ± 0.00020 | | Mac f1 −0.00104 | nmae −0.00165 (DN) |
| DN d_m200 | h1test | −0.00046 ± 0.00018 | −0.00304 ± 0.00010 | | | |
| idall3 − c4h | h1tv | −0.00244 ± 0.00007 | −0.00074 ± 0.00006 | reach f0.35 −0.0048 | −0.00244 / −0.00074, reach −0.0048 (618) | overall −0.0005 |
| idall3 − c4h | h1test | −0.00186 ± 0.00010 | −0.00223 ± 0.00006 | reach f1 +0.0045 ± 0.0034 | | |
| K 2,000 (k20 − zs) | h1tv | +0.00048 | +0.00015 | fid f1 −0.00246, f0.61 −0.00171 | | fid −0.0017 |
| K 2,000 | h1test | +0.00086 | −0.00007 | fid f1 −0.00260, f0.61 −0.00004 | | |
| K 1,600 (k16 − zs) | h1tv | +0.00102 | −0.00053 | fid f1 −0.00947, f0.61 −0.00640 | | fid −0.0041 |
| K 1,600 | h1test | +0.00034 | −0.00218 | fid f1 −0.00997, f0.61 −0.00063 | | |
| inc − zs (stack, info) | h1tv | −0.01414 ± 0.00020 | −0.01650 ± 0.00024 | fid/jac ≤ 7e-5 | | |
| inc − zs | h1test | −0.01093 ± 0.00088 | −0.01591 ± 0.00044 | | | |

- **The Mac's numbers come back.** Every 618 / 619 value the engine re-ran is reproduced: G1, G2, G3, the hook, idall3
  incl. its f 0.35 reach −0.0048, and DN. The h1test G1 also matches 611.
- **The K pairs show 619's amendment again.** At f 1 fid falls with K (−0.0025 / −0.0095, right in sign, over-read
  1.5–2.3×). At f 0.35 fid rises (wrong sign). f 0.61 sits between, and the K20 board fid −0.0017 sits on h1tv's
  f 0.61 −0.0017.
- **pds / mse:** Δpds = 0 exactly in every pair, and Δmse ≤ 5e-8 (the jackknife term, as on the board). K, and so the
  ballast, is identical between arm and base in every pair here: no pair changes the bulk.

## 5. What is not expressible yet
- **Bulk candidates (479sp2-type).** The engine prices the DE side of a bulk change only through K / the bulk passed to
  `ps_bulk`. A bulk-dump transform needs a hook that swaps the harness's `eff_bulk` (`run496.ps_bulk`), the next item for
  check800's bulk-side DE pricing. The members step already records K, sums and the mse / pds terms per target.
- **Production DNMAG module.** `inc` uses 631's `dn631.dn_pass` (`_D200`), not `src/vcc2026/dn_mag2.py`. Swapping that in
  needs a small rnup-like wrapper (the KA would be dn_mag2 == dn631 on the dumped UPRN blocks).
- **cp on h1test.** 528's call-propensity block is wired only for h1tv (run618e's `tvcp`); run612 has no cp mode.
- **Bit-equality with the Mac.** Not achievable on Linux x86-64 with this stack (§3.1). Pricing must compare arms run on
  the same platform.
- **Call-set f 0.61** is reported for every pair, but the members step keeps 521's H1 cmp constants on both surfaces, as
  619 / price630 do.

## 6. Cloud runs and cost (labcloud.sh, all on-demand, all terminated; `labcloud.sh status` clean at the end)
| run | instance | what | wall | cost |
|---|---|---|---|---|
| 20261016T015337Z-de803-b1 | r7i.8xlarge | KA runs (5, two with rnup); KA compare not bit-equal → stopped | 1.27 h | $2.72 |
| 20261016T024954Z-de803-b2 | r7i.2xlarge | item 5: h1test E0 s0 + cells → real cell-eval2 + compare | 2.83 h | $1.56 |
| 20261016T031751Z-de803-b3 | r7i.4xlarge | NumPy SIMD-dispatch diagnosis (no effect) | 1.18 h | $1.28 |
| 20261016T031921Z-de803-b4 | r7i.8xlarge | 31 non-rnup runs (Eid, Enomde, EK20/16, Emde, Ecp) | 0.92 h | $1.97 |
| 20261016T045110Z-de803-b5 | r7i.8xlarge | 10 E0 / Ec4h runs with rnup (+ a failed Mac-block probe) | 1.20 h | $2.57 |
| 20261016T060813Z-de803-b6 | r7i.2xlarge | Linux stats + rnup on the Mac's own blocks (§3.4) | 0.70 h | $0.39 |
| 20261016T060908Z-de803-b7 | r7i.large | pairs803 + kac803 on the merged outputs | 0.04 h | $0.01 |

Staging:
- **Mac → S3, once,** 10-15 22:50–23:03Z (before the 23:30Z window): 323's `reference.npz` as zstd (755 MB), the H1
  operator (1.32 GB), and `small803.tar.zst` (392 MB: ref_tv minus its two control arrays, the tv key, the KA stats).
  keys528 (49 MB) and the Mac's h1-test s0 blocks (428 MB) went 10-16 after 01:45Z.
- `ref_tv/reference.npz` is re-spliced on EC2 from 323's file and sha-checked (`stage/splice803.py`): its `cs.npy` /
  `emit_cs.npy` (3.17 GB) are byte-identical to 323's, with the same CRC.
- Donors were copied S3→S3 from `backup-20261002` and sha-checked against the Mac. `h1-test.h5ad` was copied S3→S3 from
  `input/harness/` (sha `06c6baea…`).
- **Outputs:** `work/out/` (94 arm-seed dirs, 350 MB locally; members.npz are gitignored and mirrored at
  `s3://<your-s3-bucket>/input/lab/repo/experiments/803-de-engine/outmerged/`), plus `ka-b1/`, `xcheck/`,
  `pairs803.*`, `kac803.json`.
- **Local compute:** none heavy. A production build was running throughout, so every analysis also ran on EC2 (b7).

## 7. Follow-up 1 (lab lead, 10-16 04:20Z): bulk sharpening arms on DE-H1 — does DE-H1 see S_P2 lose to P2g?

**Answer: no.**
- On both H1 surfaces, S_P2 − P2g shows no pds loss (h1test +0.0043, h1tv −0.0001; board −0.0026).
- Its nmae cost is the board's sign but only about a tenth of its size (f 1: +0.00014 / +0.00017 against the board's
  +0.00168).
- mse agrees with the board on h1test (+0.0035 vs +0.0035) and is opposite on h1tv (−0.0075).
- DE-H1, like the bulk folds, does not reproduce 479sp2's board failure. That failure is not an emission-level effect
  this harness can see.

### Wiring (599's, generalised; `emit803.install_sharpen`, arms `sh*`, `arms/pairs_sharp.json`)
- **The run.** E0's DE configuration (production emitter, CLEAN4 q80, coarse MDE, zero repair −5, zs) on the surface's
  own bulk: `run496.ps_bulk`, the 410b key, clip 0.3 / 3. Everything is equal across arms except the ballast's bulk
  render.
- **The render.** `profile_from_effect(eff_bk[i], b_ctrl, min_cpm=10, sharpen_q=q_i, protect=panel)` from the LIVE
  builder (`scripts/make_submission_pm_cd4.py` sha `e99bbd4a…`), which is production's `profile_bulk` call. The harness
  has no effect_scale / gene_scale / shrink_center / ctx-scale terms.
- **The per-target q.** From the builder's own `sharpen_q_by_tercile(bulk, carriers, genes, perts, pt)` with bulk = the
  key-scaled rows and carriers = the cap15 design rows, as production reads the loaded dump.
  - `tercile` cuts the key ratio ||bulk|| / ||carriers|| = c·w.
  - `stercile` cuts the strength ||carriers||.
- **K follows from the sharpened profile**, as in production (`ballast_count` ↔ run496.ballast_block). Every block
  asserts that the profile run496 passed equals the unsharpened one recomputed.
- **The three arms:**
  - (a) `sh175` = scalar q 1.75 on every target, 446bj2's `VCC_BULK_SHARPEN_Q`;
  - (b) `shp2g` = `tercile:2.5,2.0,1.75`;
  - (c) `shsp2` = `stercile:2.5,2.0,1.75`.

  Each arm has `pre` and `zs` stages. The board pair's base, 446bj2, carries no zs, so **pre is the board form**.
- **Tiers** (s0 manifests, `work/out-sharp/manifests-s0/`):

  | surface | P2g LOW / MID / HIGH | S_P2 LOW / MID / HIGH | targets whose q changes P2g → S_P2 | P2g key floor → LOW |
  |---|---|---|---|---|
  | h1test | 28 / 28 / 29 | 28 / 28 / 29 | 55 of 85 | 21 |
  | h1tv | 42 / 42 / 42 | 42 / 42 / 42 | 75 of 126 | 21 |

  Mean K, s0: (a) 7.40 / 7.40; (b) 8.01 / 8.62; (c) 8.45 / 8.42 (h1test / h1tv).
- **KA:** the existing base is NOT (a). The engine's E0 runs (and 435 / 618's) render the ballast UNSHARPENED,
  because run496 never sharpens. So the known answer is the wrapper with q None (`shnoop`), which reproduces E0's
  `g1pre` / `zs` members.npz on both surfaces and all three seeds (`work/out-sharp/ka_noop803.json`):
  - every stats-derived array is bit-equal: sums, jk, K, calls, every 521 regime row, ce_pds, the mse legs;
  - the only differences are cell-eval2's own per-target `ce_nmae` / mean nmae, at ≤ 2.2e-16. That is the scorer's known
    last-bit nondeterminism on identical inputs (599's `hashseed599.py`).
  - **PASS.**

### Deltas (raw, 3 seeds, mean ± seed SD; pre-zs; full table incl. zs and f 0.61 / 0.5 / cmp in `work/out-sharp/pairs_sharp.md`)
| pair | surface | Δpds | Δmse | Δnmae f 1 | Δnmae f 0.35 | Δfid f 1 | Δjac f 1 | Δreach f 1 / f 0.35 | ΔK |
|---|---|---|---|---|---|---|---|---|---|
| **P2g − a** | h1test | +0.0033 ± 0.0019 | −0.00227 | +0.00003 ± 0.00005 | −0.00002 ± 0.00012 | +0.00001 | −0.00002 | −0.0001 / +0.0019 | +0.60 |
| | h1tv | +0.0027 ± 0.0011 | −0.00511 | −0.00008 ± 0.00025 | −0.00007 ± 0.00006 | −0.00002 | −0.00000 | +0.0018 / −0.0000 | +1.21 |
| | **board (479p2g − 446bj2)** | **+0.00159** | **−0.00374** | **+0.00012** | | +0.00000 | −0.00001 | −0.00084 | |
| **S_P2 − a** | h1test | +0.0076 ± 0.0014 | +0.00123 | +0.00017 ± 0.00010 | +0.00015 ± 0.00017 | +0.00001 | +0.00000 | −0.0004 / +0.0006 | +1.06 |
| | h1tv | +0.0027 ± 0.0006 | −0.01262 | +0.00009 ± 0.00025 | −0.00011 ± 0.00008 | −0.00007 | −0.00001 | +0.0006 / +0.0006 | +1.01 |
| | **board (479sp2 − 446bj2)** | **−0.00103** | **−0.00029** | **+0.00180** | | −0.00000 | −0.00002 | +0.00078 | |
| **S_P2 − P2g** | h1test | +0.0043 ± 0.0022 | +0.00351 | +0.00014 ± 0.00009 | +0.00017 ± 0.00008 | +0.00000 | +0.00002 | −0.0004 / −0.0013 | +0.45 |
| | h1tv | −0.0001 ± 0.0010 | −0.00750 | +0.00017 ± 0.00008 | −0.00004 ± 0.00003 | −0.00005 | −0.00000 | −0.0012 / +0.0006 | −0.20 |
| | **board (479sp2 − 479p2g)** | **−0.00262** | **+0.00345** | **+0.00168** | | −0.00000 | −0.00001 | +0.00161 | |

- mse's seed SD is ≤ 7e-5 in every row: it is fixed by the bulk sums. pds per seed:
  - h1test S_P2 − P2g: +0.0052 / +0.0018 / +0.0059;
  - h1tv: −0.0011 / +0.0001 / +0.0008.
- **The zs stage reads the same:** S_P2 − P2g Δnmae f 1 is +0.00004 (h1test) / +0.00019 (h1tv); pds and mse are
  unchanged by zs.

### Reading
1. **P2g − a agrees with the board in sign on pds and mse,** on both surfaces, with nmae ≈ 0. pds is 1.7–2.1× the board,
   mse 0.6–1.4×. This is the part of the bulk channel DE-H1 can see.
2. **S_P2's board loss is not reproduced.**
   - pds: S_P2 gains over (a) on both surfaces (+0.0076 / +0.0027 vs the board's −0.0010), and over P2g on h1test.
   - nmae: S_P2 − a is +0.00017 / +0.00009 at f 1 and S_P2 − P2g is +0.00014 / +0.00017. Both have the board's sign at
     about 1/10 of its size, and are 1–2 seed SD from zero.
   - The board's "low tercile's ballast under-priced" cost (+0.0018) does not appear when the real emitter writes the
     cells. ΔK S_P2 − P2g is only +0.45 / −0.20 here, and 599 already showed K by itself costs nothing on DE-H1.
3. **mse splits by surface** (h1test agrees with the board's S_P2 − P2g, +0.0035; h1tv is −0.0075), so it is not a
   stable read. The harness has no re-level, and the board's mse reads a re-levelled bulk.
4. **The likely gap is panel / key composition, not emission.**
   - On H1 the ps key's floor block (21 targets at w = 0.3) puts P2g's LOW tier on key floor targets.
   - On A/B/C the LODO key and the strength distribution differ, and 479's board loss sat in pds plus nmae.
   - Neither DE-H1 surface nor the bulk folds reproduce it. So 479sp2 stays a miss for the checker on both channels; the
     known-answer report should mark the shape-by-strength class as "not priceable on current instruments".
   - Next lead (not run): price S_P2 − P2g on DE-H1 with A/B/C-like key statistics (a LODO-style key with A/B/C's
     floor / cap shares) to test whether the board loss is a key-composition effect.

### Cloud
| run | instance | what | wall | cost |
|---|---|---|---|---|
| 20261016T070853Z-de803-b8 | r7i.8xlarge | 24 runs (4 sharpen emits × 2 surfaces × 3 seeds) + KA + pairs | 0.57 h | $1.21 |

Follow-up spend $1.21 of its $25. 803 total $11.69. Outputs:
- **S3:** members.npz and manifests at `s3://<your-s3-bucket>/output/lab/20261016T070853Z-de803-b8/big/out/`.
- **Local** (`work/out-sharp/`): members-meta, `pairs_sharp.{md,json}`, `ka_noop803.json`, the s0 manifests of a / P2g /
  S_P2.
- **Local compute:** only file reads and sub-second JSON aggregation. `window800.sh` was CLEAR at launch.

## 8. The bulk_transform hook and 811's c0 on DE-H1 (lab lead, 10-16 09:20Z; budget $15, spent $3.30)

**Answer.** On both H1 surfaces, the many-line consensus blend has no DE cost.
- **K moves** +3 % (h1-test) and +9–10 % (h1-trainval) at λ 0.25, so h1-trainval is above rule 6's 5 % bar.
- **fid and jac** move by ≤ 4e-5 raw.
- **reach** stays inside seed noise.
- **nmae improves slightly:** −0.0004 raw at f 1 on the zs / hook stages (2.5–9 seed SDs), and −0.0010 at f 0.35 on `inc`.
- **The bulk read on the emitted cells is mixed:**
  - mse improves on both surfaces (−0.0076 / −0.0216 raw);
  - pds falls on h1-test (−0.0094 ± 0.0011, the size of 811's W-conf Δpds −0.0088) and rises slightly on h1-trainval
    (+0.0019 ± 0.0008).
- So "ballast K moved > 5 %" does not translate into a DE side effect for this transform on DE-H1, and that question no
  longer blocks the 811 or 810 candidates. Their bulk legs still decide.

### 8.1 The hook (`de803.py` / `emit803.py`)
- An arm's `emit` may carry `bulk_transform = {module, sha256, args}`, sent to the emitter step as `DE803_BULK_TRANSFORM`.
- `emit803.install_bulk_transform` loads the module by path and asserts its sha256. It then wraps `run496.ps_bulk`, so
  `module.transform(eff, genes, targets, surface, **args) -> (rows, info)` acts on the surface's composition rows: the
  cap15 matrix, BEFORE key scaling, sharpening and the ballast render. That is where production applies a stack-dump
  transform.
- The DE carriers (`eff_de_native`), the priors, the emitter and the passes stay exactly as the arm's emit says. The
  returned rows must keep the shape and be finite.
- The manifest records the module path, sha256, git blob, args, the transform's info, and whether the rows equal the
  input.
- **It composes with the sharpening hook:** the transform is applied first, so the per-target q is read on the
  transformed bulk, as production reads the loaded dump.
- **KA: PASS** (b9 phase p2, `work/out-bt/ka_bt_identity.json`). `bt/identity803.py` on E0 (`btidpre` / `btidzs`)
  reproduces E0's `g1pre` / `zs` members.npz on both surfaces and all three seeds: every stats-derived array is
  bit-equal (sums, jk, K, calls, every 521 regime row, ce_pds, the mse legs). cell-eval2's own mean nmae differs by
  ≤ 3.3e-16 and its per-target nmae by ≤ 3.7e-15, the scorer's known last-bit nondeterminism (599).
- `kanoop803.py` now takes the arm pairs as an argument and tolerates exactly that (≤ 1e-12, cell-eval2's `ce_*` only).

### 8.2 811's c0 on the H1 surfaces (`bt/c0811h.py`, `bt/ext803.py`)
- **811's code is reused verbatim but not imported.** `base811.py` imports engine801 and therefore guard801, whose write
  guard would refuse the emitter's own outputs. So the function sources `bank_rows`, `unit`, `consensus`, `blend` and the
  `GROUPS` / `COMPOSITION` literals are taken from base811.py by ast and exec'd unchanged. base811.py's sha256
  `c29f9f97…` is pinned and asserted, and its git blob is recorded per run.
- **m811's day rule** is DAY_GROUPS = GROUPS minus KOLF2 (`m811.DAY_GROUPS`), then `day_consensus(centre=False)`, which
  is `base811.consensus` on `SimpleNamespace(n, targets, nonpanel)`. The H1 group is held out, so VCC25 H1 and Pan H1
  never enter.
- **Banks for the H1 targets.** 811's compact banks carry only 0–11 of the 211 H1 targets, because ext811 extracted the
  A/B/C + L1–L6 panels only. On those banks U would have been nearly the composition donors alone.
  - `bt/ext803.py` re-extracts the same banks for the H1 targets with 811's own reader `ext811.read_bank` (imported; it
    has no guard), on EC2.
  - **Its KA PASSES** (`bt/banks_h1_manifest.json`): every source file's sha256 equals what ext811 recorded in 811's
    banks manifest, and the same reader reproduces 811's compact rows for the first 40 A/B/C targets bit for bit, every
    bank.
  - The bank sources were copied S3→S3 from `backup-20261002`, except telohaec, altair_a549 and the L-fold CD4 compacts,
    which went Mac→S3 (≈ 0.45 GB, outside the window). KA (1) sha-checks all of them.
- **Coverage on the H1 targets** (h1-test / h1-trainval):

  | line | h1-test | h1-trainval |
  |---|---|---|
  | K562 | 85 | 126 |
  | HCT116 | 84 | 119 |
  | HEK293T | 84 | 120 |
  | DLD1 | 28 | 36 |
  | RPE1 | 21 | 36 |
  | Jurkat | 21 | 35 |
  | TeloHAEC | 15 | 26 |
  | A549 | 12 | 20 |
  | BXPC3 / HAP1 / HT29 / MCF7 | 6 each | 14 each |
  | CD4 | 5 | 11 |

  Lines per target: median 4 (2–10). Every target carries U.
- **Axis.** U is mapped from the official axis onto the cap15 gene axis by name; every cap15 gene is on the official axis.
  The blend runs over the non-panel genes and keeps each row's non-panel norm to ≤ 1.1e-14.
- **Change energy** of the composition rows: λ 0.25 0.046 / 0.050; λ 0.10 0.0067 / 0.0071.
- **Arms:** the `inc` emitter configuration (Ec4h: CLEAN4 hook, coarse MDE, zr −5, zs, UPRN + DNMAG) plus the transform.
  Stages `c0l25{c4hpre,c4h,c4hupr,inc}` and `c0l10{…}`, each against the same stage without the transform (the
  collected base runs, same platform). 3 seeds × 2 surfaces.

### 8.3 Deltas vs the untransformed arm (raw, 3 seeds, mean ± seed SD; full table `work/out-bt/pairs_bt.md`)
| arm − base | surface | ΔK rel (per seed) | Δpds | Δmse | Δnmae f 1 | Δnmae f 0.35 | Δfid f 1 | Δjac f 1 | Δreach f 1 / f 0.35 |
|---|---|---|---|---|---|---|---|---|---|
| c0_l25 + inc − inc | h1test | +2.6 / +2.9 / +3.7 % | **−0.0094 ± 0.0011** | **−0.0076** | +0.0000 ± 0.0011 | −0.0010 ± 0.0002 | +0.00004 | −0.00001 | −0.0013 ± 0.0041 / −0.0008 ± 0.0013 |
| | h1tv | **+10.3 / +9.0 / +9.0 %** | +0.0019 ± 0.0008 | **−0.0216** | −0.0004 ± 0.0005 | −0.0010 ± 0.0005 | +0.00002 | +0.00000 | −0.0006 ± 0.0016 / −0.0002 ± 0.0009 |
| c0_l25 on c4h (no UPRN/DN) | h1test | as above | same | same | −0.00044 ± 0.00018 | −0.00051 ± 0.00016 | +0.00004 | −0.00001 | −0.0014 ± 0.0042 / −0.0008 ± 0.0012 |
| | h1tv | | | | −0.00045 ± 0.00005 | −0.00047 ± 0.00004 | +0.00002 | +0.00000 | −0.0006 ± 0.0017 / −0.0002 ± 0.0009 |
| c0_l10 + inc − inc | h1test | +0.3 / +0.0 / +0.3 % | −0.0022 ± 0.0008 | −0.0045 | +0.0001 ± 0.0005 | −0.0004 ± 0.0005 | +0.00004 | −0.00001 | +0.0001 ± 0.0004 / −0.0002 ± 0.0002 |
| | h1tv | +1.6 / +1.3 / +1.4 % | +0.0005 ± 0.0008 | −0.0097 | +0.0001 ± 0.0003 | −0.0003 ± 0.0003 | −0.00002 | +0.00001 | +0.0008 ± 0.0016 / −0.0000 ± 0.0011 |

- **The DE members do not move beyond seed noise against the transform**, with one exception: nmae, which moves in the
  good direction.
  - On the stages without the UP pass (c4h, c4hpre), Δnmae is −0.00044 ± 0.00018 (h1test) and −0.00045 ± 0.00005
    (h1tv) at f 1, at λ 0.25: real, small, a GAIN.
  - On `inc`, the UPRN pass's own re-roll sensitivity (§3.3) widens the seed SD to ~0.0005–0.0011. The f 0.35 gain
    (−0.0010) stays at ~2–5 SD.
  - fid and jac stay ≤ 4e-5 raw, and reach stays inside its seed SD in every row.
  - In board units the DE side effect is ≈ +0.0001 … +0.0003 overall (nmae at INV 1.648 / 6), i.e. no cost.
- **K and DE.** The +9–10 % K rise on h1-trainval, at λ 0.25, comes with no DE cost. That agrees with 599 and with §7:
  on the real emitter, K alone does not price DE.
- **The bulk read** (pds / mse from the emitted cells, the surface's own 410b key, no re-level):
  - mse improves on both surfaces;
  - pds splits by surface. h1-test −0.0094 is a large loss, equal to 811's W-conf Δpds −0.0088; h1-trainval gains
    +0.0019.
  - These are held-out bulk reads (the H1 group is out of U), with the caveats of §7: a different key from production's
    and no re-level.
- **λ 0.10** moves K ≤ 1.6 % and the DE members not at all. Its bulk is pds −0.0022 / +0.0005, mse −0.0045 / −0.0097.

### 8.4 (b) 810's sgu1 on the H1 harness: NOT run, for these reasons
sgu1 replaces the composition's whitened donor BASE by `pre` (the unwhitened normrel'd donor mean, at the base's norm) on
the top z_src tertile, keeping the residual parts. The H1 harness stores only the final cap15 rows (`effects.npz`) and the
pre-whitening cap15 (`stage.npz` = `pre`). It does not store the base part. Building sgu1 needs three things:
1. **The whitened base part.** It has to be recomputed with 323's compose verbatim (`builder("cap15").apply_v9_whitening(
   stage cap15, op.whitener, cols)`, then KA: `apply_frozen_cd4_residual` of it == effects cap15 bit for bit), and with
   618's equivalent for h1-trainval. That needs the h1-trainval operator `618/work/tv/whitener_tv_pruned.npz` (1.31 GB),
   which is not staged.
2. **810's gate.** z_src is ka804's `zsrc_fold`: 522's `robust_strength` per source donor over the surface's targets.
   It must be re-implemented on the H1 surfaces without importing ka804, because of the same guard801 problem.
3. **A KA for the parts** (base + cd4 residual == cap15).

None of these is hard, and the hook (§8.1) takes the result as a module unchanged. It is ≈ 2 h of careful work plus a
1.3 GB upload, so I left it for a separate mandate rather than a "cheap" add-on.

### 8.5 Cloud
| run | instance | what | wall | cost |
|---|---|---|---|---|
| 20261016T091016Z-de803-b9 | r7i.8xlarge | ext803 + 6 identity-KA runs + 12 c0 runs (with UPRN/DN) + pairs | 1.54 h | $3.30 |

- **Mac → S3:** ≈ 0.45 GB of bank sources. S3→S3: the rest from `backup-20261002`.
- **Outputs:** members.npz and manifests on S3 at `output/lab/20261016T091016Z-de803-b9/big/out/`. Locally:
  `work/out-bt/` (members-meta, `pairs_bt.{md,json}`, `ka_bt_identity.json`, the s0 manifests of c0l25inc / c0l10inc) and
  `bt/banks_h1_manifest.json`.
- **Local compute:** file reads and sub-second checks only. A gate602 attempt for ext803 was withdrawn unrun, because a
  production build held the gate. `window800.sh` was CLEAR at each launch.

# 206 — The four-account cluster is ONE team (Prima Mente), it is not spreading, and the arithmetic does not say "private data"

**Date:** 2026-09-03 · live board pulled 16:09 UTC (640 teams) + 16 archived snapshots ·
**read-only intel.** No `vcc`, no build, no `cell-eval2`, no submission. Nothing was posted,
nobody contacted, no account created, no Discord joined, and **no query or public field carried
anything about our method** (`AGENTS.md` §14). Peak RSS < 0.5 GB. Everything written is inside
`experiments/206-cluster-method/`.

Scripts, each with a known-answer row that aborts on mismatch; logs in `work/`:
`mine.py`, `traj206.py`, `bursts.py`, `names.py`, `cd4.py`, `cohort.py`, `spread.py`.

---

## 0. The answer in seven lines

1. **The "four accounts" are one operator.** `jw_13288`'s first two model names are `pm-r12-yfix`
   and `pm-r14-pds3o`; `Hongjian Zhou`'s are `pm-r11g`, `pm-r12a`, `pm-r14b/d/e`. **No other team
   on the board, in any of 2,155 observations, has ever posted a `pm-` model name.** They
   interleave one run-number sequence. `pm` = **Prima Mente**: Hongjian Zhou's organisation string
   changed from `Oxford` to **`Oxford / Prima Mente`** on 2026-09-02. `Andrew Liu` and
   `MT Zeng` submit 5–6 s after the others in scripted bursts.
2. **It is not spreading — the board is shrinking its amplitude.** 197's floor is a function of
   raw `mse` alone, so it is *vacuous* for an over-amplified team no matter how good its
   direction. 47 observed arms raised the floor by ≥0.03; their **median |Δpds| is 0.0096**.
   Nobody acquired a direction this week; teams found their α optimum and made a direction they
   already held legible.
3. **The two "joins" are pure amplitude.** `jw_13288` `cedar_90→juniper_46`: **92.7% of the
   +0.02850 is `mse`, `pds` moved +0.0006.** `MT Zeng` against its own Kestrel line
   (`Kestrel4→Kestrel2`): **124% of the +0.01023 is `mse`, `pds` +0.0034, and `nmae` got worse.**
4. **Their direction advantage is real and it is older than the board can see.** It arrived
   between `pm-r12` and `pm-r14`, 2026-08-28 → 09-01, when their raw `mse` was 8.83 and the floor
   was vacuous. **The board carries no information about how they got it.**
5. **The refreshed number is bigger than the mandate's.** We ship `claim` at 100.0% of the
   `PRED_TRACE_CAP_K` cap (197 §1), so a team sitting where we sit inverts to
   `cos_us + (raw_us − raw_them)/2.074495`. That puts the cluster at **0.1345–0.1382 = 1.70–1.74×
   ours**, and `cqawesome` at **0.1392 = 1.76×**. The floors (0.1271–0.1321) remain the guaranteed
   lower bounds.
6. **Three *mutually independent* pipelines are above cosine 0.11**, not one: `cqawesome`
   (0.1321 floor), the Prima Mente cluster (0.1271–0.1312), `Garrrbage` (0.1147). Pairwise
   six-tuple z-distance inside the cluster is **0.004–0.077**; between the three parties it is
   **0.43–1.02**. Two of the three are lone individuals with no organisation string.
   **Three private atlases is a much worse explanation than one public lever.**
7. **The public lever exists, eight teams name it in their own descriptions, one of them
   board-priced it at +0.01325 overall, and we hold the data and underweight it by construction.**
   Context A is **Jurkat** (exp 041). The Marson/Zhu genome-scale CD4⁺ T-cell Perturb-seq atlas
   covers **284–291/300** of the panel and is a T-lineage match to A. We ship it at
   `cd4_target_f = 0.1`, whose **effective residual fraction is 0.0528** — 062's own operating
   point — and **062 selected that fraction by maximising LODO `pds`**, the one member 197 proved
   is least sensitive to direction. 197 §1.2
   separately measured that this residual is **the only resolvable cosine improvement in eight
   board arms, +6.8%.** *The dose was never swept against the channel it actually pays in.*

**Verdict on the mandate's discriminating question: METHOD/PUBLIC-DATA, not private data.** The
"10.9 exchangeable lines do not exist" arithmetic is correct, but it does not imply a private
atlas — it implies the leaders are **not on 159's exchangeable-line curve at all**, and the
cheapest documented way off that curve is a **lineage-matched donor**, which is public, which we
have on disk, and which we consume at 5% of its weight.

---

## 1. Known-answer gates — all pass

| gate | expected | measured | source | verdict |
|---|---:|---:|---|---|
| 197's floor on our own displayed row | 0.07221 | **0.072213** | 197 §2 a15 row | **PASS** |
| — its tightness (floor / true) | 0.911 | **0.9112** | 197 §2 | **PASS** |
| floor at raw `mse` 0.7515 | 0.13117 | **0.131170** | mandate / 197 | **PASS** |
| `k_eff` from production surface cos 0.03814 | 3.09 | **3.09** | 197 §4.2 | **PASS** |
| `k_eff` from four-line surface cos 0.04278 | 4.00 | **4.00** | 197 §4.2 | **PASS** |
| Hongjian Zhou's 09-03 00:30 row is `pm-r14d` at 0.8325 / 0.760 / 0.9291 / 0.5084 / 0.2454 / 0.0224 | — | **exact** | 186 §1.4 | **PASS** |
| `jw_13288` `cedar_90` at overall 0.1810; `Andrew Liu` `AL-2574` at 0.2083 | — | **exact** | 186 §1.3 | **PASS** |
| Andrew Liu / Hongjian Zhou co-fire at 2026-09-02T17:47 | detected | **Δt = 6 s** | 186 §3 | **PASS** |
| teams with raw `mse` < 1.0: Aug 28 / Sep 3 00:30 | 8 / 44 | **8 / 44** | 186 §1.2 | **PASS** |
| `pm-r14d` attributed to Hongjian Zhou, org `Oxford / Prima Mente` | — | **exact** | live board | **PASS** |
| GISL live row `GISL v15 …` at 0.17233 with "Marson" in the description | — | **exact** | live board | **PASS** |

**197's `k_eff` is a band, not a point.** 197 §4.1 reports it over three board↔surface anchors
(a15 ×2.078, a23 ×2.041, abc-knee ×1.954); §4.2's "board-equivalent 0.0873" uses the a23 anchor
while §4.1's leader row uses a15. Every `k_eff` below is therefore quoted as a range over all
three. This is a bookkeeping note about 197's presentation, not a correction to its result.

---

## 2. The refreshed cluster table (live, 2026-09-03 16:09 UTC, 640 teams)

Auto-detected: every live row inside tolerance on **all six** members of the cluster centroid.
It returns exactly the four accounts and nobody else.

| team | model | overall | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` | subs | **cos floor** | **cos "like us"** | ×ours | `k_eff` band from the floor |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Hongjian Zhou | `pm-r14e` | 0.21064 | 0.8331 | 0.7533 | 0.9289 | 0.5087 | 0.2524 | 0.0223 | 18 | **0.13033** | 0.13736 | 1.73 | 10.2–12.1 |
| jw_13288 | `juniper_46` | 0.20949 | 0.8333 | 0.7515 | 0.9287 | 0.5085 | 0.2441 | 0.0223 | 15 | **0.13117** | 0.13821 | 1.74 | 10.4–12.4 |
| Andrew Liu | `AL-1466` | 0.20793 | 0.8330 | 0.7599 | 0.9296 | 0.5084 | 0.2469 | 0.0222 | 12 | **0.12710** | 0.13414 | 1.69 | 9.6–11.3 |
| MT Zeng | `Kestrel2` | 0.20744 | 0.8330 | 0.7591 | 0.9308 | 0.5087 | 0.2442 | 0.0224 | 12 | **0.12749** | 0.13453 | 1.70 | 9.6–11.4 |
| **cqawesome** | `m0903c` | **0.22082** | **0.8703** | **0.7495** | 0.9279 | 0.4941 | 0.2563 | 0.0286 | 27 | **0.13211** | 0.13915 | 1.76 | 10.6–12.6 |
| Garrrbage | `v14` | 0.19960 | 0.8514 | 0.7858 | 0.9726 | 0.5104 | 0.2110 | 0.0374 | 15 | **0.11465** | 0.12168 | 1.54 | 7.3–8.6 |
| Aginglab.com | `GeroAI_v19` | 0.21288 | 0.8228 | 0.8185 | 0.8890 | 0.5135 | 0.2570 | 0.0288 | 11 | 0.09887 | 0.10590 | 1.34 | 5.1–5.9 |
| **us** | `The Mitochondria Decli…` | **0.18342** | **0.8360** | 0.8738 | 0.9442 | 0.5037 | 0.2299 | 0.0302 | 26 | 0.07221 | **0.079249** | 1.00 | 2.5–2.9 |

The last column is `k_eff` **from the floor**, so it is comparable across the table; our own
*measured* cosine 0.079249 gives 3.1–3.5, which is 197 §4.2's 3.09 at the a15 anchor.

**Two refinements to the mandate's figures.**

- `cqawesome` is now the board's **highest inverted direction**, not the cluster — it took the
  lead on `mse` today (`m0903c`, 0.7495). The mandate's table predates it.
- The floor is conservative for a reason 197 states and we can now use. Our own `claim` sits at
  **100.0%** of `PRED_TRACE_CAP_K` (197 §1's `rc_hybrid` row), which is exactly the assumption
  the floor makes. A team sitting where we sit therefore inverts to
  `0.079249 + (0.8738045 − raw)/2.074495` — the **"like us"** column. It is a point estimate, not
  a bound; if their `claim` is *below* the cap their true cosine is higher still. Both columns
  agree the answer is **1.7–1.8×**, which is 197's tightness-corrected figure arrived at from a
  different direction.

`informative rows (raw mse < 1.023611): **65 of 640**` — up from 54 of 612 on 09-03 00:30.

---

## 3. The joining diffs, member by member — and they are not what the mandate expected

### 3.1 `jw_13288`: `cedar_90` (09-02) → `juniper_46` (09-03), +0.02850 overall

| member | raw before | raw after | Δ raw | Δ scaled /6 | share of the move |
|---|---:|---:|---:|---:|---:|
| `pds` | 0.8327 | 0.8333 | **+0.0006** | +0.00016 | **0.6%** |
| `mse` | 0.9023 | 0.7515 | **−0.1508** | +0.02642 | **92.7%** |
| `nmae` | 0.9370 | 0.9287 | −0.0084 | +0.00249 | 8.7% |
| `fid` | 0.5088 | 0.5085 | −0.0003 | −0.00013 | −0.4% |
| `reach` | 0.2455 | 0.2441 | −0.0013 | −0.00020 | −0.7% |
| `jac` | 0.0229 | 0.0223 | −0.0005 | −0.00025 | −0.9% |
| **cos floor** | 0.05849 | 0.13117 | **+0.07268** | | |

### 3.2 `MT Zeng`: two readings, because `Petrel1` was a failed side-arm

`Petrel1` (#10, 09-02) is a **regression** from `Kestrel4` (#8, 09-01) — the bird names split into
a `Kestrel` line and a one-off `Petrel`. Against the immediate predecessor the diff flatters the
move; against their own Kestrel line it is the honest one.

| | vs `Petrel1` (#10) | vs `Kestrel4` (#8) |
|---|---:|---:|
| Δ overall | **+0.04627** | **+0.01023** |
| Δ `pds` | +0.0319 (25.3%) | **+0.0034 (12.1%)** |
| Δ `mse` | −0.1207 (45.7%) | **−0.0723 (124.0%)** |
| Δ `nmae` | −0.0041 (1.7%) | **+0.0129 (−36.7%, WORSE)** |
| Δ `fid` | +0.0119 (14.1%) | +0.0010 (5.4%) |
| Δ `reach` | +0.0249 (10.1%) | **−0.0044 (−7.9%)** |
| Δ `jac` | +0.0031 (3.1%) | +0.0007 (3.1%) |
| Δ cos floor | +0.05817 | **+0.03484** |

### 3.3 What the two diffs identify

**The component that moved is the amplitude, and only the amplitude.** In both cases `pds` — the
scale-invariant member — is pinned to within +0.0006 / +0.0034, `reach` and `jac` are flat or
negative, and in MT Zeng's case `nmae` got *worse*. `mse` carries 93% and 124% of the move.

This is the same signature the whole board is showing (§4), and it means **the mandate's premise
that a method spread onto two new accounts today is not what happened.** Both accounts already
held the cluster's direction — `cedar_90` was already at `pds` 0.8327 and `Kestrel4` at 0.8296,
i.e. inside the cluster on the rank channel — and what landed on 09-03 was the α calibration that
makes that direction *visible to the floor*.

**Corollary, and it is the load-bearing one:** these two arms are **not** a controlled experiment
on the cluster's direction. They are a controlled experiment on its amplitude. The direction step
is `pm-r12a` → `pm-r14b` (08-28 → 09-01: `pds` +0.0545, `nmae` −0.1200, raw `mse` 8.83 → 0.82), and
at `pm-r12a` the floor was **−3.76** — vacuous. **The board cannot see when or how they acquired
the direction, only that they had it by 09-01.**

---

## 4. "It is spreading" — measured, and refuted

197's floor `cos ≥ (1.023611 − raw)/2.074495` reads **raw `mse` and nothing else**. An
over-amplified team has a vacuous floor however good its direction. So the growth of the
high-floor cohort is ambiguous between "more teams got a better direction" and "more teams found
their amplitude optimum". `pds` separates them: it is scale-invariant, so an amplitude-only arm
moves the floor a lot and `pds` by ≈0.

| snapshot | teams | floor > 0 | > 0.05 | > 0.09 | **> 0.11** | best floor | `pds` > 0.83 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 08-21 19:06 | 84 | 3 | 0 | 0 | 0 | 0.016 | 0 |
| 08-25 19:21 | 257 | 15 | 2 | 0 | 0 | 0.074 | 1 |
| 08-28 13:25 | 385 | 14 | 4 | 1 | 0 | 0.095 | 1 |
| 08-30 15:55 | 474 | 28 | 10 | 1 | 0 | 0.095 | 3 |
| 09-02 06:30 | 581 | 50 | 32 | 8 | 1 | 0.126 | 7 |
| 09-03 00:30 | 612 | 54 | 32 | 10 | 4 | 0.127 | 8 |
| **09-03 16:09** | **640** | **65** | **42** | **12** | **6** | **0.132** | **14** |

**47 observed arms raised the cosine floor by ≥0.03 into an informative regime. Their median
|Δpds| is 0.0096; the max is 0.3054 (one team recovering from a broken submission).** For scale,
197 §3.1 measured that our own cosine span corresponds to a `pds` offset of ~0.09.

**So the wave visible on the board is the α-shrink wave 186 §1.2 already identified, continuing.**
It is revealing directions, not creating them. Six of the top rows to raise the floor this week
did it with `pds` inside ±0.01: `jw_13288` (+0.0006), `Aginglab.com` `GeroAI_v11→v12` (+0.0003),
`Ju` `m114→m120` (+0.0001), `Happy` `HP-025→HP-027` (+0.0008), `minnesota university`
`m116→m124` (−0.0000), `Garrrbage` `v13→v14` (+0.0060).

---

## 5. Attribution

### 5.1 The cluster is Prima Mente — established, with the evidence

| evidence | detail |
|---|---|
| **the `pm-` namespace, exclusive** | Across 2,155 distinct `(team, submission_count)` observations in 17 snapshots, exactly **7 rows** carry a `/^pm[-_]/i` model name. All 7 belong to two accounts: `Hongjian Zhou` (`pm-r11g`, `pm-r12a`, `pm-r14b`, `pm-r14d`, `pm-r14e`) and `jw_13288` (`pm-r12-yfix`, `pm-r14-pds3o`). **No third team, ever.** |
| **interleaved run numbers** | `pm-r11g` 08-27 06:34 (HZ) · `pm-r12-yfix` 08-27 07:52 (jw) · `pm-r14-pds3o` 08-28 02:52 (jw) · `pm-r12a` 08-28 06:23 (HZ) · `pm-r14b` 09-01 (HZ). The suffixes `-yfix` / `-pds3o` / `a` / `b` / `d` / `e` are variant tags on one run-number series. jw_13288 is the account that burns the extra slots on ablations. |
| **the organisation string** | `Hongjian Zhou`'s org read `Oxford` through 09-01 and **`Oxford / Prima Mente`** from `pm-r14d` (09-02) onward. The only other Oxford string on the board is `University of Oxford` (`PROSPECTor`, unrelated). |
| **scripted co-submission** | Within-cluster co-fires at **Δt = 5, 6, 6 s**; cross-team pairs within 60 s run **2.28× the Poisson null**. `MT Zeng`#8 / `Andrew Liu`#8 at 09-01T15:17 (5 s); `Hongjian Zhou`#17 / `Andrew Liu`#10 at 09-02T17:47 (6 s); `jw_13288`#15 / `MT Zeng`#12 at 09-03T04:11 (6 s). |
| **their own first two model names** | `effect-transfer-v1-k562gw`, `effect-transfer-v5-h1blend` (Hongjian Zhou, 08-22 / 08-25) — **donor effect transfer off K562 genome-wide plus an H1 blend.** The same family we run. `pm-r14-pds3o` additionally leaks that they were optimising `pds`. |

**Who they are.** Hongjian Zhou is a DPhil student at Oxford (IBME, supervised by David Clifton;
Clarendon Scholar, St Edmund Hall), publishing on foundation models in healthcare —
[ibme.ox.ac.uk/person/hongjian-zhou](https://ibme.ox.ac.uk/person/hongjian-zhou/),
[Scholar](https://scholar.google.com/citations?user=6JvttIYAAAAJ&hl=en). **Prima Mente** is an
AI-biology company (London / San Francisco / Dubai, CEO Ravi Solanki) that in July 2025 released
**Pleiades**, a 90M/600M/7B whole-epigenome DNA-methylation foundation model trained on 1.9T
tokens on 256 H200s, in a preprint co-authored with **University of Oxford**, Hebrew University
and Gothenburg — [primamente.com](https://www.primamente.com/),
[launch release](https://www.businesswire.com/news/home/20250716401893/en/Prima-Mente-Launches-Pleiades-a-Large-Scale-Human-Whole-Epigenome-Foundation-Model-with-Clinical-Applications-for-Alzheimers-Disease),
[NVIDIA case study](https://www.nvidia.com/en-us/case-studies/primamente/). They state they
**built a wet lab for high-throughput multiomic data generation** and describe themselves as a
company that "generates its own data", and they are currently hiring a
[Bioinformatics Engineer — Virtual Cells](https://www.builtinsf.com/job/bioinformatics-scientist-virtual-cells/8073389).

**Motive and means for private data: yes. Evidence that they used it here: none.** See §6.

### 5.2 What is NOT public

- **No Prima Mente code, model or dataset on Hugging Face** (searched `prima-mente`,
  `PrimaMente`, `pleiades` — models and datasets; zero hits attributable to the company).
- **No Prima Mente or Hongjian Zhou preprint on the VCC, on arXiv or bioRxiv.** The arXiv
  perturbation/virtual-cell feed since 2026-08-15 contains **7 items**, of which two are relevant
  and both already known to 186: `PopPert` (2609.01357, 09-01) and `D²R²` (2608.15288, 08-15);
  plus one new one, **`DeMixPert` (2608.23114, 08-24)** — decomposed response modelling with
  Gaussian mixtures for OOD single-cell perturbation prediction. Nothing matches the cluster's
  numbers.
- **No public repo.** `gh search repos "virtual cell challenge 2026"` returns two unrelated
  personal repos (`biona001/Virtual-Cell-Challenge-2026`, `forrestsheldon/…`). Our own repo is
  confirmed **private** (`gh repo view`: `isPrivate: true`) — it appeared in `gh search code`
  only because the search runs with my access.
- **All four cluster accounts have an empty `description` on every submission**, as does
  `cqawesome` on all 27 and `Garrrbage` on all 15. The top six rows on the board publish nothing.

### 5.3 `cqawesome` — nothing new, stated plainly

**No public trace found, and 186's item 2 stands unchanged.** Sole member `Chenxi Qiu`,
organisation `"."`, 27 submissions, **every description empty**. `m0903c` landed 09-03 15:23 at
0.22082 (+0.00292 over `m0902b`), taking the board's best raw `mse` (0.7495) and best `pds`
(0.8703). Their whole visible history is `test_water_v2/v3/v6` → `test_fire_v1` → `mp4` →
`E3a_ctrlprog` → `m0830b` → `m0902a/b` → `m0903c`. No GitHub, HF, arXiv or bioRxiv result. **The
leader's method remains dark.** What *is* new and useful is negative: their six-tuple z-distance
to the cluster is **0.434–0.446**, an order of magnitude larger than the cluster's internal
0.004–0.077, so **`cqawesome` and Prima Mente are not the same pipeline** (§6.1).

### 5.4 Lingshu-Cell — 186 §M6 CONFIRMED on fresh rows

The live board still separates them cleanly, and the separation is now larger, not smaller:

| team | model | overall | `pds` | `mse` | `reach` | cos floor | `k_eff` band |
|---|---|---:|---:|---:|---:|---:|---:|
| huhuhu | `(lingshu+fofj21)` | 0.19185 | 0.8219 | 0.8390 | 0.2133 | 0.08899 | 4.0–4.6 |
| Leo | `(CPA+i0bhu2)` | 0.19120 | 0.8217 | 0.8411 | 0.2129 | 0.08798 | 3.9–4.5 |
| GZ YAI | `(lingshu+194pnb)` | 0.18684 | 0.8191 | 0.8452 | 0.2089 | 0.08601 | 3.7–4.3 |
| cindy chen | `(CPA+t4gxwd)` | 0.16932 | 0.8152 | 0.8938 | 0.2013 | 0.06259 | 1.9–2.1 |

All four `(model+hash)` wrapper rows sit at cos floor **0.063–0.089**, i.e. `k_eff` 2–5, and
**0.042–0.069 below the cluster's floor**. `reach` is still pinned at **0.209–0.213** across
`lingshu`, `CPA` and `STATE` — 186 §M6's "all land at `reach` ≈ 0.208", reproduced on rows that
have since been resubmitted. Swapping the architecture still moves nothing. **Lingshu-Cell is not
the cluster's mechanism.** Refuting 186 would have required these rows to converge on the
cluster; they diverged further.

---

## 6. Is the advantage DATA or METHOD? — the discriminating evidence

### 6.1 Three mutually independent pipelines are up there, not one

Pairwise six-tuple distance, RMS of per-member gaps normalised by each member's trimmed board sd
(small = the same prediction, large = independent pipelines):

| | cqawesome | HZhou | jw_13288 | A Liu | MT Zeng | Garrrbage | Aginglab | us |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **cqawesome** | – | 0.434 | 0.446 | 0.445 | 0.445 | 0.707 | 0.241 | 0.307 |
| **Hongjian Zhou** | 0.434 | – | **0.077** | **0.052** | **0.076** | 1.023 | 0.418 | 0.537 |
| **jw_13288** | 0.446 | **0.077** | – | **0.026** | **0.004** | 0.998 | 0.433 | 0.513 |
| **Andrew Liu** | 0.445 | **0.052** | **0.026** | – | **0.025** | 1.010 | 0.432 | 0.525 |
| **MT Zeng** | 0.445 | **0.076** | **0.004** | **0.025** | – | 0.996 | 0.431 | 0.511 |
| **Garrrbage** | 0.707 | 1.023 | 0.998 | 1.010 | 0.996 | – | 0.707 | 0.490 |

The cluster is one prediction (0.004–0.077). `cqawesome`, `Garrrbage` and the cluster are three
separate objects (0.43–1.02). And two of the three are **lone individuals with no organisation
string** — `cqawesome` (`Chenxi Qiu`, org `"."`) and `Garrrbage` (`Ken`, org empty) — on
trajectories that look like hobbyists, not data programmes: `cqawesome` went `test_water_v2` →
0.22082 in thirteen days.

**A private-atlas explanation now has to be paid for three times.**

### 6.2 The `1/k` arithmetic is right, and it points somewhere else

| quantity | board cos | `k_eff` band |
|---|---:|---:|
| production (us, 197's measurement) | 0.07925 | **3.1 – 3.5** |
| best equal consensus of all FOUR public panel-covering lines (197 §4.2) | 0.08730 | 3.8 – 4.4 |
| `Garrrbage` v14 (floor) | 0.11464 | 7.3 – 8.6 |
| cluster floor (`juniper_46`) | 0.13117 | **10.4 – 12.4** |
| cluster, "like us" estimate | 0.13821 | 12.0 – 14.4 |
| `cqawesome` floor / "like us" | 0.13211 / 0.13917 | 10.6 – 12.6 / 12.3 – 14.7 |
| context-independent `k → ∞` ceiling | 0.25480 | ∞ |

166's census is unchanged and I re-confirmed its central external claim today: **X-Atlas/Pisces
is still an empty repository** — `usedStorage: 0`, three files (`.gitattributes`, `LICENSE.md`,
`README.md`), `lastModified 2026-03-25T22:58:18Z`, unchanged for five months. So the four
panel-covering lines are still four.

**The exchangeable-line reading demands 10–14 lines. `AGENTS.md` §5.3 counts ≈11 distinct human
CRISPR-perturbed cell lines *in total*, and 166 showed only four of them cover the panel** — the
panel excludes common-essentials by construction, which zeroes `nadig_jurkat`, `nadig_hepg2`,
`rpe1` and `k562_essential` at exactly 0/300. **Nobody can be at `k_eff` 12 by stacking
exchangeable public lines. Not them, not a private lab either** — Prima Mente's wet lab would
have to have run ~8 more genome-scale CRISPRi lines to buy it that way.

**Therefore the leaders are not on 159's curve.** That is the real content of the arithmetic, and
it is the opposite of the mandate's reading. 159's law caps at `cos_surface` 0.12263 as `k → ∞`
*because κ² ≈ 0.34 is the **cross-lineage** shared component*. A donor that matches the
recipient's **lineage** shares the context-specific component too, so it is not one of the `k` and
the law does not price it. **Getting off the curve is cheaper than climbing it.**

### 6.3 The lever, named in public by eight teams, board-priced by one

Exp 041 identified the contexts from control profiles alone, same top match under all eight
variants: **A = JURKAT (T-ALL), B = HELA, C = CAL33**. Two competitors state the same thing in
public without knowing it is confirmed: `AlphaBio` — *"the strongly **Jurkat-like context A**;
epithelial contexts B/C"*; `peng yu` — *"Local hold-out avg_score 0.3825 (**Jurkat** panel)"*.

The T-lineage donor that matches A is the **Marson / Zhu genome-scale CD4⁺ T-cell Perturb-seq
atlas**: CRISPRi knockdown of all expressed genes across **22 million primary CD4⁺ T cells from
4 donors**, at rest and at two post-stimulation timepoints —
[bioRxiv 2025.12.23.696273](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1),
now in [Cell (S0092-8674(26)00929-3)](https://www.cell.com/cell/fulltext/S0092-8674(26)00929-3),
CC BY-NC; data via
[CZI Virtual Cell Models](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq),
analysis code [emdann/GWT_perturbseq_analysis_2025](https://github.com/emdann/GWT_perturbseq_analysis_2025).
**166 already logged it (`GSE314342`, ~284/300 coverage) and we have it on disk.** `ElNino` and
`VRLAB` both independently report **291/300** panel coverage for it.

**Who says so on the board, in their own words:**

| team | published quote | ovr |
|---|---|---:|
| **GISL** (Columbia) | *"…each library held out as truth over the real 300 val targets, which retrodicts **Marson +0.0120 (board +0.01325)**, the six-panel ensemble +0.0007 (+0.00045) and KOLF-at-full-weight −0.0056 (−0.00245)."* | 0.17233 |
| `None` | *"context A direction = **primary CD4 T Perturb-seq (Marson 2026, Stim48hr, z≥3)**"* ; later *"Marson 2026, **w = 0.5 all contexts**"* | 0.08571 |
| `ElNino` | *"A = **cd4t 3-donor mean (291/300)**"* | 0.11263 |
| `VRLAB` | *"Context A: **CD4T Rest mean-delta overlay (291/300 perts)**"* | 0.03021 |
| `jhtao` | *"Context A retains the selected **CD4T lambda-1 augmentation** … Tests context-specific transfer from a **fifth perturbation source**."* | 0.10495 |
| `Yusif Mukhtarov` | *"…fusion of public K562 genome-wide CRISPRi, X-Atlas HCT116, and **primary CD4 (Rest/Stim8/Stim48)** responses"* | 0.14008 |
| `Zhou Li` | *"…held-out recipient lines (HCT116, HEK293T, **primary CD4 T cells**) … the sign gate now also hears a **genome-wide primary T-cell CRISPRi screen**"* | 0.12529 |
| `望有幸` | *"trained on **primary CD4 CRISPRi DE statistics** and Replogle K562"* | −0.03990 |
| `Kevin SUN` | *"A: **CD4T Rest** alpha 1.0; B/C: k562 filled2 alpha 1.5"* | 0.01526 |

**GISL's +0.01325 is a board measurement, published by them, of adding this source as a library.**
It was measured on a team whose raw `mse` is 5.0 — i.e. with `mse` clamped to scaled 0 — so the
gain arrived entirely through `pds` + `nmae` + `fid` + `reach`. **For us, in the informative `mse`
regime, the same direction gain would move `mse` as well, so +0.01325 is a floor on what it is
worth to us, not a ceiling.** At 197's exchange rate it corresponds to **+0.0363 cosine**, which
would put us at 0.1156 against the cluster's 0.1312 floor.

### 6.4 And here is the part that is ours: we selected the dose on the wrong metric

Production ships the CD4 source as a **post-whitening residual at `cd4_target_f = 0.1`**, whose
*effective* fraction is **0.0528** (`191` §: `effects.cd4.residual_frame.effective_f_at_cd4_fraction`
= 0.05278). That weight traces to **exp 062**, whose selection rule was, verbatim:

> *"selection: maximize worst HCT116/HEK293T fold **gain** … | held target role | matching baseline
> **PDS** | selected **PDS** | delta |"*

**062 chose the CD4 dose by maximising LODO `pds`.** 197 §3 then proved that `pds` is the member
*least* sensitive to the direction — we extract **26% more `pds` per unit of bulk cosine** than the
median shrunk team, so `pds` systematically under-rewards exactly the channel this source
supplies. And 197 §1.2 independently measured the consequence:

> **the CD4 residual (a14 → a15) is the ONLY resolvable cosine improvement in eight board arms:
> 0.07418 → 0.07925, +6.8%.**

So we have a measured dose-response point — **effective fraction 0.0528 buys +0.00507 cosine** —
and the dose above it has never been tested in the current (post-193, `WHITEN_LAM = 0`) geometry.
The one time a larger dose was tried, **exp 060 at fraction 0.50 pre-whitening**, it lost 0.143
PDS; 062's own conclusion was that *"source-only whitening amplified CD4 directions that did not
transfer, rather than showing that all CD4 information is harmful."* **That experiment tested the
wrong geometry and was scored on the wrong member, and it has been the reason the dose stayed at
0.05 ever since.**

---

## 7. Ranked, testable hypotheses for the cosine advantage

Priced at 197's **+0.00365 overall per +0.01 cosine**. Our gap to 1st is **+0.03740 overall**;
the whole gap to the cluster's floor is **+0.05192 cosine = +0.01895 overall**, and to the
"like us" estimate **+0.05896 cosine = +0.02152 overall**.

### H1 — The CD4⁺ T-cell dose is 5% and was chosen on `pds`. **Evidence: very strong. Size: +0.004 to +0.013 overall.**

Three independent measurements point at it: 197's +6.8% cosine at effective fraction 0.0528 (our
only resolvable direction win, ever); GISL's published board +0.01325 for shipping it as a
library; and 041's identification of context A as Jurkat, which two competitors state in public.
It is a **lineage-matched** donor, so it is off 159's exchangeable-line curve and the "no fifth
line exists" census does not bind it.

**How to test.** A dose sweep of `VCC_CD4_TARGET_F` in the *current* geometry
(`alpha 0.982, gamma 0.65, lam 0`), scored on **`mse`-channel objectives, never on `pds` alone**.
- **Known-answer row, built in and mandatory (audit 181):** the sweep must retrodict **a15 at
  effective fraction 0.0528 → cosine 0.07925 / raw `mse` 0.8738** before any other point is read.
  Both numbers are in 197 §1 and the board.
- Rig: PLOLO with the recipient restricted to a T-lineage line, plus h1-test with the donor
  restriction of `AGENTS.md` §6.4 trap 1.
- Second arm: promote CD4 from a **post-whitening residual** to a **weighted consensus member**
  (GISL's construction), and price the two constructions against each other.
- Control: the same sweep with the CD4 rows **target-shuffled**, which 062 never ran.
- **Falsifier:** if the cosine is flat or falls between effective fraction 0.05 and 0.20, H1 dies
  and the +6.8% at 0.05 was a one-point artefact.

**Value.** Linear extrapolation of the measured point to effective fraction 0.20 is +0.019 cosine
(+0.0070 overall); it will saturate, so the honest band is **+0.004 to +0.0133 overall**, the
upper end being GISL's board figure. Even the low end clears rule 13's gate.

### H2 — Lineage-matched donors for B (HeLa) and C (CAL33) too. **Evidence: moderate, by symmetry. Size: unpriced, plausibly +0.003 to +0.010.**

B and C are epithelial (`AlphaBio`: *"epithelial contexts B/C"*); today they receive the same
generic consensus as A. Orion HCT116/HEK293T are already epithelial-ish and already in, which is
why this ranks below H1, but the **weighting is uniform across contexts** and 041 gives us the
lineages. `jhtao` is running exactly this experiment in public, one context at a time, and their
KOLF arms read −0.0004 / −0.0001 — a **useful negative**: the wrong matched line buys nothing.

**How to test.** Per-context donor weights on PLOLO, with the known-answer row being the current
uniform weighting reproducing a15. Cheap: no new data.

### H3 — Prima Mente has private perturbation data. **Evidence: motive and means only. Size: unbounded but unactionable.**

They have a wet lab and say they generate their own data. **But nothing on the board supports it
and two things argue against it:** the same cosine is held by two lone individuals with no
organisation (§6.1), and even a private lab would need ~8 more genome-scale CRISPRi lines to buy
`k_eff` 12 the exchangeable way (§6.2).

**How to test.** We cannot test it directly. **It is falsified operationally if H1 + H2 reproduce
half the gap on our rig.** Run H1 first; H3 is the residual after H1/H2, not a hypothesis to
spend on.

### H4 — A better *estimator* of the shared component. **Evidence: against. Size: ≈0.**

194 (eleven estimators), 203 (hierarchical shrinkage), 193 (the eigenvalue map) and 195 (learned
transforms) all closed. And the board argues the same way from outside: the cluster's `pds`
(0.8330) is **below** ours (0.8360) at 1.7× the cosine, so their *geometry* is worse than ours and
their *raw direction* better. **We are not short of estimator; we are short of signal.** Do not
reopen.

### H5 — X-Cell / X-Atlas-Pisces. **Evidence: none yet. Size: unknown — 205's question, not mine.**

**Pisces is still empty today** (verified: 0 bytes, unchanged since 2026-03-25), so it cannot be
anyone's edge right now. `X-Cell` ([bioRxiv 2026.03.18.712807](https://www.biorxiv.org/content/10.64898/2026.03.18.712807v1.full.pdf))
does report zero-shot work on **a Jurkat CRISPRi line**, which makes it a lineage-matched
candidate for context A and therefore an H1-shaped lever rather than a new-line lever.
`kriagr` ships *"STEER **Pisces** control-weighted"* at 0.13342 — worth one look by 205 to see
whether they have access we do not.

### H6 — The cluster simply calibrated better than us. **Evidence: refuted for the gap, confirmed for the timing. Size: 0.**

187 closed calibration: the global-α optimum is a15's **0.982 to four decimals**. The cluster's
2026-09-02/03 arms *are* calibration (§3), but calibration is how their direction became
**visible**, not how it became **good** — the floor is a lower bound on cosine at any amplitude,
and at raw 0.7515 it is 0.13117 whatever they did with α. **No amount of α buys us their number.**

---

## 8. What I could NOT determine

1. **How Prima Mente's direction was actually obtained.** It arrived in the `pm-r12 → pm-r14` arm
   (2026-08-28 → 09-01) when their raw `mse` was 8.83 and the floor was **−3.76** — vacuous. The
   board records that they had it by 09-01 and nothing about where it came from. This is the
   central failure of the mandate's premise that the joining diffs are a controlled experiment on
   the method: they are a controlled experiment on the amplitude.
2. **Whether Prima Mente used private data.** Motive and means are established (wet lab, funded
   company, virtual-cell hires); use is not. No preprint, no repo, no HF artefact, no description.
3. **`cqawesome`.** Still completely dark after a second full sweep: no GitHub, HF, arXiv,
   bioRxiv or web trace for `Chenxi Qiu`; 27 empty descriptions. They now hold the board's best
   `mse`, best `pds` and best inverted cosine, and we know nothing about how.
4. **`Garrrbage` (`Ken`).** The third independent party at cos floor 0.1147. Fifteen submissions,
   **zero descriptions**, no org, no repo. Their `v13 → v14` is a clean amplitude-only arm
   (Δ`pds` +0.0060) on a direction they already had at `pds` 0.8454.
5. **Whether the CD4 dose response is monotone above effective fraction 0.05.** This is H1's whole
   question and I did not run it — it needs the rig and this was an intel mandate. The only
   adjacent datapoint, 060 at fraction 0.50, was measured pre-whitening and is not comparable.
6. **Whether GISL's "Marson" arm is the same construction as our 062 residual.** They call it a
   *library* held out *as truth* in a leave-one-library-out validator; we ship it as a
   post-whitening additive residual. The +0.01325 may be partly the construction rather than the
   source, and I cannot separate them from outside.
7. **The exact identity of `MT Zeng` and `Andrew Liu`.** Both have empty org strings. `MT Zeng`'s
   sole listed member is **`Michael Zhou`** (team name and member name disagree), which is
   suggestive next to `Hongjian Zhou` but is not evidence — the burst timing and the shared
   pipeline are.
8. **Whether Prima Mente has more than four accounts.** Four are proven by the six-tuple; a fifth
   contributing only failed arms would be invisible. The `pm-` namespace bounds it at two
   *named* accounts.
9. **Per-context decomposition of anyone's advantage.** The board reports the three-context
   average only. If the cluster's edge were concentrated in context A — which is what H1 predicts
   — the board could not show it. **This is the single measurement that would settle H1 vs H3 and
   it is not available from outside.**
10. **`X-Cell`'s actual numbers and licence** — left to 205, whose RESULTS had not landed when
    this was written.

---

## 9. Hand-offs

- **→ coordinator, the headline.** *"The method is not public"* is confirmed. *"The arithmetic says
  it is data"* is **half right and misread**: the arithmetic says the leaders are **off 159's
  exchangeable-line curve**, which no amount of data — public or private — can put them on. The
  cheapest documented exit is a **lineage-matched donor**, we hold one, and we weight it at 5%
  because 062 priced it on `pds`. **Do not stop internal direction work; redirect it at H1.**
- **→ whoever owns the next direction arm (H1).** The dose sweep of `VCC_CD4_TARGET_F` in the
  current geometry, with a15 at effective 0.0528 → cosine 0.07925 as the mandatory known-answer
  row and a target-shuffled control. This is the only direction lever left that has *three*
  independent positive measurements behind it and has never been swept.
- **→ 205.** Pisces is still 0 bytes today, verified. `kriagr` (0.13342) publicly claims a Pisces
  source — check whether that is real access or aspirational naming. `X-Cell`'s Jurkat CRISPRi
  line is an H1-shaped lever (matched to context A), not a fifth exchangeable line.
- **→ `AGENTS.md` §8.2.** Two updates. (a) The leading direction is now `cqawesome`'s **0.13211
  floor / 0.13915 point estimate**, not the cluster's; "at least 1.6× ours" should read
  **"at least 1.66×, and 1.76× on the cap-consistent point estimate."** (b) `docs/CURRENT_STATE.md`'s
  line *"more donor lines — **closed**"* is correct for *exchangeable* lines and **wrong as
  written for lineage-matched ones**; 166's census answers a different question than the one
  §8.2 uses it for.
- **→ §14 discipline.** Reconfirmed, at cost to everyone but us: this experiment reconstructed a
  competitor's donor set, validator design, per-arm board deltas and dose choices out of nine
  teams' description fields, and got Prima Mente's identity out of two model names and an
  organisation string. **The top six rows on the board publish nothing.** Neither do we.

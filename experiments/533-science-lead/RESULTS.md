# 533 — Science lead: the manuscript brought up to date as an audited record (Draft v15)

Worker 533 (Claude Code, Opus 5.5), mandate `MANDATE.md`, 2026-09-29 13:37Z → 14:20Z. Analysis and
writing only: no build, no `vcc`, no `cell-eval2`, no git, no emitter run; every process ≤ 0.3 GB RSS. Everything written
lives in `experiments/533-science-lead/`.

**Status: the complete Draft v15 exists (first draft done 2026-09-29 ~14:15Z); the 09-30 read is PENDING at 8 marked
places, with the text to insert for every outcome in §4.**

Deliverables:
- `paper/` — the updated manuscript: `PAPER.md` (Draft v15, 195.8 k chars vs v14's 156.1 k), `CLAIMS.md` (C-74 … C-93),
  `REFERENCES.md` (four new keys), `make_figs.py` + `figs/` (F-12, F-13 new; F-9 and F-10 corrected; all 13 regenerated
  from saved tables in 2.7 s), `build_html.py` (points at `paper/`).
- `paper.html` — the built page (2.8 MB, 31 contents entries, 13 figures). **Private: never post, upload or share.**
- `work/apply_v15.py` + `work/v15_blocks_{a,b,c,d}.py`, `work/apply_claims_v15.py` + `work/v15_claims.py`,
  `work/apply_refs_v15.py`: v15 is rebuilt from byte copies of the v14 files (`work/*.v14.orig.md`) with every anchor
  asserted unique, so the coordinator can review each insertion as a block and re-run after any upstream edit.
- `work/verify_board.py` + `work/verify_board.log`: the board arithmetic behind every board number v15 states.

To move it into `294-writeup`: copy `paper/{PAPER,CLAIMS,REFERENCES}.md`, `paper/figs/*`, `paper/make_figs.py`; in
`build_html.py` keep 294's `W` and take the other two edits (figure range 1–13, masthead lines). `make_figs.py`'s `REPO`
line goes back to two `..` levels in 294's directory.

---

## 0. Summary

- **Board arithmetic reproduces exactly.** All 15 captured entries used: overall = mean of the six scaled members to
  ≤ 5.6e-17. All six new BOARD-RESULT files' raw deltas reproduce to 1e-6 and their overall deltas to 1e-7 (§2.1).
- **Thirteen discrepancies found** in the record, none of which changes a decision (§2.3). The largest: the CD4 half's
  board transfer is 1.83× 438's price, not "k ≈ 1.43" (D13); a decoder with no board row counted as losing to us (D5); a
  FINDINGS line that reads a hook's increment on a leaky base as a total (D6); two stale figure captions in the manuscript
  itself (D1, D2).
- **Twenty new claims (C-74 … C-93)**, one new finding (F20, §3(l)), a new Discussion subsection (§4.7 open problems), two
  new figures. Negatives absorbed in §4.4 (516, 521, 525) and §4.5 (T9 local prices vs the board; T10, T11).
- **The central open question of §3(j) is still open**: the 09-30 entries (458d, 458g) had not been submitted when the
  draft was finished. §4 below gives the insertion text for positive / partial / null / harm (and the unbinned gap
  −0.0010 < Δfid ≤ −0.0005, which the preregistration leaves unassigned).

---

## 1. Inventory: every result since the last revision (commit `f31ad98a`, 2026-09-27 07:56Z)

Status key: **board** = read on the scored board; **H1** = local h1-test instrument only; **folds** = local PLOLO /
cancer-fold instruments only; **analytic** = computed from saved statistics; **negative** = a closure at its
construction; **pending** = preregistered, not read.

| item | finding (one line) | where | status | in v15 |
|---|---|---|---|---|
| 446bj2 (board, 09-27 00:27Z; predates f31ad98a by 7.5 h but was never in the paper) | the leave-one-donor-out bulk key transfers: +0.0022338, mse −0.014473 raw past its −0.004 … −0.012 range, pds +0.000208 below +0.001 … +0.004 | `E446/BOARD-RESULT-446bj2.md` | board | §3(k), T9, C-93 |
| 446br2 (board, 09-27; same) | zero-gene repair floor −2: nmae −0.001413 against −0.002 … −0.008, partial; −5 stays | `E446/BOARD-RESULT-446br2.md` | board | §3(i), §3(k) table, T9 |
| 453dr (board, 09-28) | the October day-of path on relabelled A/B/C: −0.0036816 vs 446b, bulk −0.0030866 inside its band: VALIDATED | `E453/BOARD-RESULT-453dr.md` | board | §3(h), F-12, C-82 |
| 446bj2z (board, 09-28) | the zs exact-rank pass transfers: nmae −0.002179, the other five members bit-identical; 1.19× H1 | `E456/BOARD-RESULT-446bj2z.md` | board | §3(i), C-74 |
| 457dr (board, 09-29) | CD4 ×1.5 on the October path: +0.0031484, all bulk, 1.83× 438's k 1 price (BOARD-RESULT says k ≈ 1.43: D13); with 453dr the refit costs −0.0005332 at the tuned dose | `E457/BOARD-RESULT-457dr.md` | board | §3(h), F-12, C-82 |
| 446bd2 (board, 09-29) | DLD-1 as a second union residual: +0.0004200, FLAT (mse −0.000471 vs −0.001 … −0.004) | `E445/BOARD-RESULT-446bd2.md` | board, negative | §3(h), T9, C-92 |
| 458d / 458g (board, 09-30) | Arc's H1 programme as every context's block prior; its design-strength gate | `E458/`, `E468/` preregistrations | **pending** | §3(j) marker, C-89, §4 |
| 513 | reach is sign-limited, not order-limited; called-gene magnitudes at H1's L1 optimum (513f) | `E513/RESULTS.md`, `work/{walk513,demote513}.json` | H1, analytic | §3(l) (C-86), §3(i) (C-75) |
| 514 | DepMap dependency does not predict cross-line response size on non-essential targets (r +0.001) | `E514/RESULTS.md` (stdout only) | folds, negative | §3(h), C-81 |
| 515 | Wiener denoising of one donor row: +0.009 pool-sign accuracy on H1 (consensus-rider class, not queued) | `E515/RESULTS.md` (stdout only) | H1 | subsumed by 516's denoising rows (§4.4) |
| 516 | target-specific public sign sources: +0.016 on H1, ≈ 0 on every cancer fold; oracle +0.030; stacking ≤ +0.002 | `E516/RESULTS.md`, `work/` | H1 + folds, negative | §3(j), T11, §4.4, §4.7, C-79 |
| 517 | the Arc-sourced uncalled-gene hook (−0.0183 H1, 12× CLEAN4's); its lower bar; the meanz_ne estimator; same-lab analogue; identity axis reverses in 3/6 lines; Flex library-size and dropout carry no sign (dropout predicts calls) | `E517/RESULTS.md` §0–§9, `work/*.json` | H1 (real emitter, 3 seeds) + public folds | §3(i), §3(j), F-13, C-76 … C-78, C-80 |
| 518 | the profile-error identity 1.032 − 2X/ΣD; the magnitude pot (+0.018 … +0.023 overall); strength key as a contingency | `E518/RESULTS.md`, `work/*.json` | analytic + folds | §3(k), §4.7, C-84, C-85 |
| 519, 526 | two cloud code reviews of the October path: four + four defects, each reproduced and fixed | `E519/RESULTS.md`, `E526/RESULTS.md` | engineering | §5.6 (Methods note) |
| 520 | no public source names D/E/F; every high-coverage bank registered; plausible lines 0–14/300 | `E520/RESULTS.md` | field fact | §3(h), C-91 |
| 521 | per-target call-set sizing closed by payoff; global K is a regime bet; board weak end thinner; "unexpressed → null" refuted | `E521/RESULTS.md`, `work/` | H1 + folds, negative | §3(l), C-87, C-88 |
| 522 | October's bulk key stays the percentile cap (P 96.647); October-like tails are real (two donors agree) | `E522/RESULTS.md` | folds (engineering) | §4.4 (a composition caveat) |
| 523 | the global call-set dose pair (K 1,600 / 2,000) preregistered; builder knob committed | `E523/PREREGISTRATION.md` | **pending** (from 10-06) | §3(l), C-88 |
| 524 | the October recipe dress rehearsal (10-05) | `E524/PREREGISTRATION.md` | pending (engineering) | not integrated |
| 525 | no public DE-emission design explains the mid-pack; the cellpred cluster is private | `E525/RESULTS.md` | field fact, negative | §3(h), §4.4, §4.7, C-90 |
| 527 | 446bj2zmz (the Arc-sourced hook on 446bj2z) built: KA-1 and KA-496 byte for byte, diff exactly `de_mde_sign` + emitter pin | `E527/RESULTS.md`, `PREREGISTRATION.md` | **pending** (10-01 s2 if 458d is not HARM, else ≥ 10-07) | §3(i), C-76 |
| 528 | call-propensity block (517 §7's lead) | `E528/RESULTS.md` | in progress (another worker) | not integrated; §3(j) mentions the lead |
| 529–531 | October pass-through bundle plan / cloud bundle / final-phase tooling | `E529/PLAN.md`, `E530`, `E531` | engineering | not integrated |
| 532 (landed 09-29 ~13:55Z) | literature and data refresh 09-15 … 09-29: nothing changes a decision; no D/E/F names; no new Arc Flex data | `E532/RESULTS.md` | field fact | §3(h) (re-check), C-91 |
| 534 (landed 09-29 ~14:00Z) | generalization audit Part A (AGENTS §11.1c) | `E534/RESULTS.md` | in progress | not read into v15 |
| 540 … 546 (committed 09-27 04:48–05:13Z, before f31ad98a, but absent from the paper) | the LODO key's cap: cap 6 +0.0019 … +0.0029 (folds); robust / kolf2 / blob-free keys negative; stack mix optimal; percentile cap for October | `E540 … E546/RESULTS.md` | folds, engineering | §3(k) (one line: "a raised cap … priced on local folds only"); §4.4 |

Board chain for context (overall): 446b 0.2700428 → 446bj2 0.2722766 → 446bj2z **0.2728260** (best); 453dr 0.2663613,
457dr 0.2695096, 446br2 0.2704006, 446bd2 0.2704629 (all recomputed, §2.1).

---

## 2. Verification (the 503 standard)

### 2.1 Board numbers → `data/board-entries/<id>.json` (`work/verify_board.py`, log `work/verify_board.log`)

- **Overall = unweighted mean of the six scaled members** for all 15 entries used (395, 410, 410b, 417e, 433a, 433h, 433k,
  446b, 446bja, 446bj2, 446br2, 453dr, 446bj2z, 457dr, 446bd2): max |difference| 5.6e-17.
- **Every raw member delta stated in the six new BOARD-RESULT files reproduces** (|d| < 6e-7), and every overall delta
  (|d| < 6e-8): 446bj2 +0.0022338, 446br2 +0.0003578, 453dr −0.0036816, 446bj2z +0.0005493, 457dr +0.0031484, 446bd2
  +0.0004200.
- **Derived statements recomputed:** 457dr − 446b = −0.0005332 (bulk +0.0000141, DE −0.0005473) ✓; 453dr − 446b bulk
  −0.0030866, DE −0.0005950 (BOARD-RESULT: −0.00309 / −0.00059) ✓; 457dr − 453dr bulk +0.0031006; 446bj2z − 446b =
  +0.0027831; 446bja − 446b = −0.0010661 ✓ (C-63).
- **Scaled members.** The captured records carry context-averaged raw members and context-averaged scaled members, not the
  per-context anchors b and r, so s = (u − b)/(r − b) cannot be recomputed from the records. Checked instead through the
  implied Δscaled / Δraw = context-weighted 1/(r − b), which must be stable across arms and consistent with AGENTS §4's
  anchor ranges: mse −1.051 … −1.054 (4 arms; r − b ≈ −0.95 ✓), reach +1.088 … +1.130 (r − b ≈ 0.87–0.91 ✓), pds
  +2.02 … +2.29 (r − b ≈ 0.43–0.48 ✓), nmae −1.51 … −1.89 (r − b ≈ −0.57 … −0.63 ✓).
- 446bj2z's five other members are bit-identical to 446bj2's (equality of the stored doubles), and 446br2's fid and jac
  bit-identical to 446b's.

### 2.2 Experiment numbers stated in v15 → primary files

| v15 number | source | check |
|---|---|---|
| reach by pool size 0.643 / 0.078 / 0.027 / 0.014, walk accuracy 0.81 / 0.51 (32 / 9 / 19 / 22 targets) | `E513/work/walk513.json` (per-target table) | **recomputed** here from the per-target rows: identical |
| demotion Δreach −0.0007 … −0.038 | `E513/work/demote513.json` | **recomputed**: demote/demote2/blocklast over B 1,000 … 4,000 give −0.0007 (demote, B 3,000) … −0.0383 (demote2, B 4,000); 513's text says −0.003 (§2.3 D8) |
| 513's known answers (fid 0.508929…, jac 0.074699…, nmae 0.942704…, reach 0.269257…) | `walk513.json` → `known_answers` | read: all `ok: True` |
| sign-accuracy decay 0.80 → 0.55; ≥ 3-donor subset 0.972 on 72 genes; scale s-scan (+0.0157 / +0.0103 / +0.0356); called 0.600 / uncalled 0.537 | `E513/RESULTS.md` | **stdout only** (`keys513.py`, `purity513.py`, `scale513.py` saved no output) |
| 516 H1 best +0.0157 [+0.0069, +0.0240], nmae −0.0067 [−0.0110, −0.0024], reach −0.030 [−0.067, −0.002] | `E516/work/best516.json` (`all85.acc.all`, `members`) | read: identical |
| 516 folds +0.0033 [+0.0004, +0.0062] / −0.0002 [−0.0027, +0.0023] / +0.0010 [−0.0030, +0.0053] | `E516/work/fold516.log` lines 42, 199, 357 | read: identical |
| 516 power +0.0148 … +0.0177 | `E516/work/power516.log` | read |
| 516 oracle per-target best source +0.030 | `explore516a.py` line 44 prints it | **stdout only** |
| 517 real runs (c4-mzg0.05: nmae −0.01832, seeds −0.018295 / −0.018257 / −0.018410; fid +0.00022; jac −0.00012; reach −0.00387; ov +0.00444, range +0.00374 … +0.00488; seed-0 CI [+0.00305, +0.00662]); CLEAN4 hook −0.001471 | `E517/work/real517.json` | **recomputed** 3-seed means here |
| attenuation curves and crossings (block gated/ungated ov 0 at r 0.237 / 0.171; fid 0 at 0.255 / 0.263; mean-z hook ov 0 at 0.282, nmae = CLEAN4 hook at 0.415; P_down hook 0.365 / 0.380) | `E517/work/{grid517,hookgrid517,mzgrid517}.json` | **recomputed** (seed means; linear interpolation); grid517's block-only rows are the `tau: None` rows (the others carry a hook) |
| same-lab analogue +0.036 / +0.081 / +0.019, +0.055 / +0.112 / +0.038; HCT116 block 0.488 / 0.523 / 0.557; complement 0.561 / 0.506 | `E517/work/labhook517.json` | read (0.5424 − 0.5063 etc.) |
| essential rows 0.500 / 0.480; meanz_ne 0.548 / 0.573 vs P_down 0.528 / 0.544 | `E517/work/progest517.json` via RESULTS §2.6 | RESULTS table consistent with the json keys |
| identity axis −0.0128 / −0.0111 / −0.0025 / −0.0071 / +0.0029 / +0.1045; HCT116→HEK293T 0.595 vs 0.541; HEK293T→HCT116 0.538 vs 0.542; Jurkat→HEK293T 0.448 | `E517/work/xline517.json` (`non-essential` block for K562/HCT116/HEK293T) | **recomputed** |
| library-size: eps×δ 0.527 / 0.523; Arc 0.640 → 0.618; CLEAN4 0.573 → 0.563; eps reliability 0.969; Spearman −0.445 (strength) | `E517/work/tech517.json` | read |
| dropout call rates 0.0328 / 0.0672 / 0.1094; capture 0.3672 vs 0.3386 | `E517/work/{dropout517,callprop517}.json` | read / averaged over the two folds |
| 518 identity: 2X/ΣD 0.40585, constant 1.03226, −0.004059 raw and +0.000711 overall per 1 % X; rivals ×1.138 / 1.133 / 1.336 | `E518/work/board518.json` | read; **arithmetic checked**: 1.032 − 0.406 = 0.626 ≈ 446bj2's raw 0.626408; 0.406 × 0.01 = 0.00406; 0.00406 × 1.052 / 6 = 0.00071; (0.462 − 0.406) × 1.052 / 6 = 0.0098 |
| magnitude pot +23.2 / +25.1 / +27.1 %, kolf2 +52.2 %, cancer mean +25.13 % | `E518/work/folds518.json` (`r_mag`) | read; −0.406 × 0.2513 = −0.102 raw; × 1.31 = −0.134; overall × 1.052 / 6 = +0.0179 / +0.0234 |
| 521: K562 thr4 r +0.601; S4 +0.551; E[u] +0.000405 / −0.003425 (8.46 : 1), break-even 0.89 (K 10) / 0.97 (K 25) | `E521/work/{pred521,purity521}.log` | read; ratio recomputed |
| 521: K 1,600 −0.01546 [−0.02511, −0.00720] (f 1), +0.00479 [+0.00090, +0.00996] (f 0.35) | `E521/work/global521.log` lines 7, 79 | read |
| 514: r +0.001 (7,051), −0.041 (8,781), within-line −0.354 / −0.062, −0.266 / −0.046 | `E514/RESULTS.md` | **stdout only** (work/ is empty) |
| 520 / 525 / 532 field facts | the documents they quote | not re-fetched (§14: nothing leaves this machine); transcribed with their own VERIFIED/UNVERIFIED labels |

### 2.3 Discrepancies (for the coordinator; none changes a decision)

| # | file:line | what it says | correct statement |
|---|---|---|---|
| D1 | `experiments/294-writeup/make_figs.py:848`, `figs/MANIFEST.md:55` | F-9 caption: DE members "within 0.0006" | within **0.0007** (largest move nmae −0.000643); 503 corrected C-61 but not the caption. **Fixed in v15.** |
| D2 | `294-writeup/make_figs.py:893`, `figs/MANIFEST.md:63`, `PAPER.md:822` | F-10's source is `labline450.json` | the figure plots `labline450_v3.json` (`make_figs.py:858`); the comment at `:855` also says "8 banks, 3 labs". **Fixed in v15.** |
| D3 | `294-writeup/PAPER.md:809` | "All eight are regenerated" | eleven in v14 (thirteen in v15). **Fixed.** |
| D4 | `294-writeup/PAPER.md:100` (§1.4) | board numbers "are confined to §3(h) and the Discussion" | stale since §3(i)–(k) (v11–v13). **Fixed.** |
| D5 | `docs/FINDINGS.md:28–29`; `E525/RESULTS.md:203` | danielchen26's decoder "loses all four DE members to ours" | no board row is mapped to that repository (`E525/RESULTS.md:127`), so its board members are unmeasured; its own shift knob measured flat |
| D6 | `docs/FINDINGS.md:38` | "the gated programme + hook at +0.0088" | +0.00876 is the **hook's increment on the gated-programme base** (real517 `g33-mzg0.05`, one seed), and its fid +0.0071 is from flipping leaked complement calls (≈ 590 unsuppressible per target on H1's g33 base; A/B/C builds leak ~1), which 517 §2.1 says to expect as ~0 on A/B/C |
| D7 | `E525/RESULTS.md:42` | "chill's 0.835 from our H1 0.9427 needs about +0.08" | the +0.08 comes from the **board** gap (446b 0.8960 − chill 0.8351 = 0.061 at 480's −0.006 … −0.008 per +0.01; `E480/RESULTS.md:40–41, 156`); from H1's 0.9427 it would be ~+0.13 … +0.18 (525's own §6 warns to keep the units apart) |
| D8 | `E513/RESULTS.md:28` | block demotion "Δreach −0.003 … −0.038" | `work/demote513.json`: **−0.0007 … −0.0383** (demote at B 3,000 −0.0007, B 4,000 −0.0013); the conclusion (every scheme loses) holds |
| D9 | `E513/RESULTS.md:19`; `docs/FINDINGS.md:53` | A/B/C's median n_real 31–45 "(AGENTS §2)" | the number is **192's board-inferred estimate** (`E192/RESULTS.md:44, 352`, f ∈ [0.35, 0.5]); AGENTS §2 does not contain it, and 521 §8.2 disputes the proportional assumption behind it |
| D10 | `E521/RESULTS.md:36, 247`; `docs/FINDINGS.md:86` | "K 1,200–1,600 … ≈ +0.005 … +0.006 at f 0.35 and in the compressed regime (CIs exclude 0)" | at f 0.35, K 1,200 is **+0.0032 [−0.0034, +0.0098]** (P ≤ 0 0.16); only K 1,600 excludes 0 there (+0.0048 [+0.0009, +0.0100]); in cmp both do (+0.0064 [+0.0031, +0.0098], +0.0061 [+0.0038, +0.0090]) (`work/global521.log:77–79, 114–115`) |
| D11 | `E517/RESULTS.md:341` (§7 table) | "Spearman −0.45 with the bank row's norm (300), −0.44 with h1-test n_real" | `work/tech517.json`: **−0.4448** and **−0.4348** (−0.44 / −0.43 at two decimals) |
| D12 | `E514/work/`, `E515/work/` (empty); `E513` (`keys513`, `purity513`, `scale513`); `E516/explore516a.py:44` | tables quoted in RESULTS with no saved output | traceability gap, not an error: re-run with `\| tee work/<name>.log` when convenient (each < 1 min, < 1 GB) |
| D13 | `experiments/457-cd4-half-def/BOARD-RESULT-457dr.md:25`; `docs/FINDINGS.md:19` | 457dr "landed … k ≈ 1.43 of 438's price" | 438's price is +0.00172 at k 1 (`E457/PREREGISTRATION.md:13, 19`); +0.0031484 / 0.00172 = **1.83** (k ≈ 1.83). 1.43 is the ratio to the k 1.3 point (+0.0022). The bin reading (≥ +0.0010, below +0.0040) is unchanged |

Scope note on the mandate's "+0.06 … +0.08 gate sign accuracy over ours": 434f's own map gives ~0.62 for the 0.87–0.88 pack
(+0.06) and ~0.65–0.66 for cqawesome/chill (+0.09 … +0.10); 480's rate gives ~+0.08 for chill; the 09-29 leader (nmae 0.804)
interpolates to ~0.69. v15 §4.7 states the range with its sources and the upper-bound caveat.

### 2.4 What could not be verified here

- Per-context scaled members (anchors not captured; §2.1).
- Stdout-only tables (D12): quoted in v15 with their source marked (C-75, C-79's oracle ceiling, C-81).
- 520 / 525 / 532's external documents: not re-fetched; their VERIFIED/UNVERIFIED labels are carried.
- The claim rows' "candidate NEW" provenance: graded by me, not by an audit (as the v12–v14 rows were).

---

## 3. What changed in v15

**Sections** (all insertions are blocks in `work/v15_blocks_{a,b,c,d}.py`):
- Top: Round 14 revision note. Abstract and plain-language summary: one added passage each (DE signs are the limit; no
  public source supplies them across cancer lines), each ending in a 09-30 marker.
- §1.4: the "confined to §3(h)" sentence corrected.
- §3(h): the pipeline refit (453dr / 457dr; C-82, F-12a); a fifth line flat (446bd2; C-92); dependency (514; C-81); field
  facts three weeks after 306 (520, 525, 532; C-90, C-91).
- §3(i): the second delivery on the board (446bj2z; C-74) and the partial repair floor (446br2); called genes at the L1
  optimum (513f; C-75); the programme-sourced uncalled-gene hook and its bar (517; C-76, C-77; F-13).
- §3(j): F18's source list; target-specific public sources negative (516; C-79); recipient-side and technical features
  (517; C-80); a same-lab programme from another line helps where the design is weakest (517; C-78); **the 09-30 test**
  (design, bins, exact decomposition; C-89; PENDING); "Why it matters" gains a PENDING clause and a public-data sentence.
- §3(k): the six-entry table and its reading (C-83); what the profile-error member reads once the bulk is free (518; C-84,
  C-85) and the transferred LODO key (C-93).
- **§3(l) new** (F20): where reach lives (C-86); per-target call-set sizing closed by payoff (C-87); the global size as a
  regime bet (C-88, 523 pending); a benchmark-design note.
- §4.1: later candidates (the hook's bar; reach and sizing; the mse identity) and a PENDING line under F18.
- §4.4: five bullets (516's and 517's scope; no competitor claim beyond opened code; call-set sizing; H1 prices are not
  board claims; the refit is about A/B/C relabelled).
- §4.5: T9 (six local prices against the board), T10 (the identity axis), T11 (516's H1 gain vs the board).
- §4.6: which statements are measured on A/B/C, updated.
- **§4.7 new**: open problems (signs on weak and uncalled responses; the public information shown insufficient — 516, 525,
  174, 184, 450h/462/466, 446bja, 517; what data would supply it; the magnitude pot; the call-set regime; the construct).
- §5.6: two known-answer lessons from 479/499, 519, 526. §5.8: Round 14's derived objects. Figure index: F-12, F-13,
  corrected F-10 source, "All thirteen".

**Claims:** C-74 … C-93 in a Round 14 addendum (`CLAIMS.md`), columns as the ledger's; C-89 is the pending 09-30 row.

**Figures:** F-12 (the October path's decomposition 446b → 453dr → 457dr, split bulk vs DE; six single-channel entries
read by the two channels; board entries only). F-13 (the hook's H1 price against attenuation r, against the block prior;
517's saved grids plus its real-emitter runs as markers). Both are reproducible by `make_figs.py` from saved tables.

**References:** VCC26-CELL, JIANG25, ATLASSHIFT, DEPMAP-CHRONOS (transcribed from 520 / 525 / 514; nothing fetched).

---

## 4. Pending items

### 4.1 The 09-30 read (458d, 458g) — the central open question of §3(j)

At the time of writing (09-29 14:17Z) neither entry existed: `experiments/459-h1-axis-def/work/submit0930s{1,2}.log` held
only their ARMED lines (10:54Z), and `data/board-entries/` had no 458d/458g capture. Read with
`experiments/468-prior-gate/read468.py` (the coordinator writes the BOARD-RESULT files); then run
`experiments/517-flex-programme/predict517.py --fid-458d <Δ> --fid-458g-vs-458d <Δ> --nmae-458d <Δ>` for the implied r.

**Markers in `paper/PAPER.md` (8):** the revision note (line 22), the abstract (91), the plain-language summary (106),
F18's header (601), the §3(j) result paragraph (816), §3(j) "Why it matters" (827), §4.1 under F18 (1002), §4.7 item 3
(1162); plus `CLAIMS.md` C-89. (`grep -n "PENDING 09-30" paper/PAPER.md`.)

**Falsifiers first** (458's and 468's bins): pds and mse must equal 446b's to ≤ 1e-6 raw; jac within ±0.0010, except the
pre-diagnosed leak band −0.001 … −0.004 (then note the leak's fid cost ≈ −0.001 … −0.002 beside the fid bin). Any other
pattern is a bug: insert nothing; write "the 09-30 read failed its falsifiers; see BOARD-RESULT" at the §3(j) marker.

Placeholders below: ⟨Δfid⟩ etc. are raw 458d − 446b; ⟨G⟩ is Δoverall 458g − 458d; ⟨S⟩ = 458g − 446b (programme on the
strong two thirds) and ⟨W⟩ = 458d − 458g (programme on the weak third), both overall; ⟨r⟩ from `predict517.py`.

**(a) §3(j) result paragraph** — replace the bold PENDING paragraph with one of:

- **POSITIVE (Δfid ≥ +0.0015).** "*The read (2026-09-30): the producer's programme carried.* 458d − 446b moved fid by
  ⟨Δfid⟩ raw, in the preregistered 'carries' bin (≥ +0.0015); nmae ⟨Δnmae⟩, reach ⟨Δreach⟩, jac ⟨Δjac⟩⟨; the jac move is the
  pre-diagnosed call leak⟩; pds and mse identical. That is ⟨Δfid/0.0198⟩× of the programme's own-line H1 fid, and on
  §3(i)'s attenuation curve it corresponds to r ≈ ⟨r⟩. The gate read ⟨G⟩ overall (⟨gate bin⟩): by the exact decomposition
  the programme was worth ⟨S⟩ on the design-strong two thirds and ⟨W⟩ on the weak third. So a programme estimated by one
  laboratory in an embryonic stem-cell line carried to three cancer lines the same laboratory screened with the same
  chemistry, where an independent screen of the same line (Jurkat, 446bja) did not. Among these lines the programme's
  laboratory-and-chemistry component transfers across lines on the scored member."
- **PARTIAL (+0.0005 ≤ Δfid < +0.0015).** As POSITIVE, with "carried in part" and "(the preregistered partial bin,
  +0.0005 … +0.0015)", ending: "The laboratory component transfers, weakly: at r ≈ ⟨r⟩ a block prior is near its break-even
  (r ≈ 0.24 gated) while the uncalled-gene hook still beats the consensus hook on nmae (§3(i)); the 0.5 blend (458e) is its
  dose partner."
- **NULL (|Δfid| < 0.0005).** "*The read (2026-09-30): the producer's programme did not carry.* 458d − 446b moved fid by
  ⟨Δfid⟩ raw, inside the preregistered null bin (|Δ| < 0.0005); nmae ⟨Δnmae⟩, reach ⟨Δreach⟩; pds and mse identical. The
  gate read ⟨G⟩ (⟨gate bin⟩; strong two thirds ⟨S⟩, weak third ⟨W⟩). Neither an independent screen of the same line
  (Jurkat, 446bja) nor the same laboratory's screen of another line (458d) supplies a challenge line's programme on the
  scored member. The H1 cross-lab transfer (Pan) stays the exception, plausibly tied to an embryonic line's strong own
  programme. A null block read corresponds to r ≈ 0.25–0.27 on §3(i)'s curve, where the uncalled-gene hook still beats the
  consensus hook on nmae (its bar is lower); the hook's board arm (446bj2zmz) is the next test."
- **HARM (Δfid ≤ −0.0010).** "*The read (2026-09-30): the producer's programme hurt.* 458d − 446b moved fid by ⟨Δfid⟩ raw
  (preregistered harm bin, ≤ −0.0010); nmae ⟨Δnmae⟩, reach ⟨Δreach⟩; pds and mse identical. The gate read ⟨G⟩ (strong two
  thirds ⟨S⟩, weak third ⟨W⟩). The H1 programme's embryonic-specific content dominates what it says about cancer lines: on
  §3(i)'s curve this is r ≥ 0.28, where the hook's six-member value is ≈ 0 even though its nmae still beats the consensus
  hook."
- **UNBINNED (−0.0010 < Δfid ≤ −0.0005).** The preregistration assigns no bin here. Report it as "a negative point
  estimate between the preregistered null and harm bins", read nmae as the secondary (≤ −0.0015 confirms a signal,
  ≥ +0.0015 contradicts fid), and use the NULL text for everything downstream (rule1001b treats anything above −0.0010 as
  not HARM).

**(b) Abstract** — replace the bold PENDING sentence with: POSITIVE "The benchmark producer's own screen of a different
line did supply it: its programme as every context's prior moved the scored sign member by ⟨Δfid⟩ (preregistered
positive bin)." · PARTIAL "…supplied part of it (⟨Δfid⟩, the preregistered partial bin)." · NULL "The benchmark producer's
own screen of a different line did not (⟨Δfid⟩, inside the preregistered null bin)." · HARM "…hurt (⟨Δfid⟩)."

**(c) Plain-language summary** — POSITIVE "The organisers' own screen of a different cell line did carry it, which
means the laboratory and its method matter as much as the cell line." · PARTIAL "…carried a little of it." · NULL "The
organisers' own screen of a different cell line did not." · HARM "…made it worse."

**(d) §3(j) "Why it matters"** — POSITIVE/PARTIAL: delete the bold PENDING clause and change "not by who did it or how"
to "or, as the 09-30 read shows, by whether the same laboratory has perturbed another line with the same chemistry".
NULL/HARM/UNBINNED: delete the bold PENDING clause; the sentence stands.

**(e) F18's header** (the "tested … by the 09-30 entries: [[PENDING 09-30]]" clause) — POSITIVE/PARTIAL: "…and the
producer's own programme from another line carried to its three validation lines (458d, C-89)." NULL/HARM: "…and the
producer's own programme from another line did not carry (458d, C-89)."

**(f) §4.1 under F18** — POSITIVE/PARTIAL: "The laboratory component transfers across lines on the challenge's own lines
(458d, C-89); caveats: one producer, three lines, one emission realisation, and fid's 0.07–0.67× attenuation history."
NULL/HARM: "The laboratory component did not transfer on the challenge's own lines (458d, C-89)."

**(g) §4.7 item 3** — keep the branch that happened and delete the other: POSITIVE/PARTIAL "The producer's own screen of
another line carried (458d): a screen's laboratory and chemistry are worth more than its line, and a benchmark producer's
own screens are the most valuable public resource for its benchmark." NULL/HARM "It did not (458d): a line's programme is
available only from a screen of that line."

**(h) C-89** (`CLAIMS.md`, status column) — "preregistered; read 2026-09-30: ⟨bin⟩ (Δfid ⟨Δfid⟩ raw; gate ⟨G⟩)"; interval
column "exact; no re-roll noise between the three entries"; add the BOARD-RESULT paths and the capture ids to the evidence
column.

**(i) Revision note** — replace the PENDING bullet with "The 09-30 entries answered F18's open question: ⟨one line⟩ (C-89)."

The gate's own bins (468; Δoverall 458g − 458d): ≥ +0.0015 the gate transfers; +0.0005 … +0.0015 weakly positive (kept,
flagged unconfirmed); −0.0005 … +0.0005 null; ≤ −0.0005 harm. H1 had put the gate's own benefit at +0.0043
[−0.0034, +0.0132] (C-71), so any gate read is the first measurement of it with no re-roll noise.

### 4.2 Other reads that will change claims

| read | when | claim(s) | how |
|---|---|---|---|
| 446bj2zmz (527, the Arc-sourced hook) | 10-01 s2 if 458d is not HARM, else ≥ 10-07 | C-76, C-77, §3(i) | bins in `E527/PREREGISTRATION.md` (Δnmae ≤ −0.0030 transfers); report the transfer factor vs H1's −0.0183 × the r implied by 09-30 |
| 496a (the CLEAN4 hook) | 10-05 | C-73, C-76 | the hook class's H1 → board factor (T = Δnmae / −0.00146) |
| 495 (LODO key cap 6) | 10-03 | C-85, C-93 | the magnitude pot's second board step |
| 523 (K 1,600 / 2,000) | from 10-06 | C-88, F20 | Δ(fid + jac) bins in `E523/PREREGISTRATION.md` |
| 524 (October recipe dress rehearsal) | 10-05/06 | C-82 (engineering) | would extend the refit validation to the full October recipe |

---

## 5. Process record

- Reading and tracing: AGENTS.md, FINDINGS §0–§2, the v14 manuscript and ledger, 503's audit, the BOARD-RESULT files,
  the preregistrations of 433-0926 / 445 / 453 / 457 / 458 / 468 / 523 / 527, and the RESULTS and work files of 513–527,
  532, 540–546 (plus 174, 184, 192, 434, 435, 480 for specific numbers).
- Computation: only `work/verify_board.py` (board JSONs), small in-line reads of saved JSON/log tables (513, 516, 517,
  518, 521), `paper/make_figs.py` (saved tables; peak 0.27 GB) and `paper/build_html.py` (pandoc; peak 0.14 GB). No raw
  matrix was opened, no statistic was computed on one, and no process came near the memory gate.
- Nothing was written outside `experiments/533-science-lead/` (and a safety copy in this session's scratchpad); no git,
  no `vcc`, no `cell-eval2`, no builds.
- At completion (2026-09-29 14:17Z) the 09-30 captures did not exist (the slots submit from 00:01Z 09-30), so the draft
  carries the 8 markers and §4.1's insertion texts, as the mandate's item 4 prescribes.

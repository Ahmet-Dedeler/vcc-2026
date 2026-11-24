# 190 — the regulatory-prior ordering key: **NO-GO. M2 is closed.**

Date **2026-09-03**. Measurement only. No `vcc`, no `cell-eval2`, no build, no submission,
nothing posted anywhere, nothing written outside this directory except the network downloads in
`data/external/regulon/` that the mandate authorised. Peak RSS **7.43 GiB** (the one-off
co-expression build), **4.61 GiB** for PLOLO, **0.70 GiB** for everything else.

---

## 0. Decision in one line

**A signed regulatory prior — CollecTRI, propagated or not, and a control-derived co-expression
GRN — adds nothing to `|H|` as an ordering key, and the reason is not coverage: on the pairs
where each network speaks, its own expected sign is LESS accurate than the sign `|H|` already
predicts.** The best key in a 13-point λ sweep over eight networks reads **`reach` 0.23670
against `|H|`'s 0.23533**, a raw **+0.00137 = +0.00026 overall**; its paired-bootstrap 95% CI is
**[−0.00140, +0.00520]**, and a **shuffled regulon selected the same way reads 0.23699 ± 0.00059
— higher than the real one**. Every network's permutation p-value is between 0.15 and 1.00.

**The hardest number in the experiment is the ceiling.** Give a *perfect oracle* control of every
(target, gene) pair the direct CollecTRI regulon covers and leave the rest at `|H|`: `reach` goes
to **0.23612, raw +0.00078 = +0.00015 overall**. That is **13× below the composed-slot-gate's
shape-only ship floor of +0.0019**, and it is an upper bound no prior built on that network can
cross. Coverage closes the direct version; *information* closes the propagated and co-expression
versions, whose ceilings are all at or above the ship floor (+0.0021 to +0.0401 overall) and whose
achieved values are zero to 2.8 % of them.

**`tydymy`'s board arm does not retrodict.** A CollecTRI correction at their published 25% of
backbone L2 reproduces **three of their five moves — `nmae` worse, `mse` worse, `pds` slightly
down — and gets the two that carry their entire gain backwards**: `reach` −0.0131 where they
measured +0.0353, `fid` −0.0004 where they measured +0.0133. And falsifier 4b shows the three
that do reproduce are reproduced **just as well by a shuffled network of identical topology and
L2**, so the cost side carries no information about the biology either. Whatever bought
`tydymy` +0.0141 gross on 2026-08-26, this experiment cannot make it be the CollecTRI prior.

`nmae` therefore stays where 187 put it — a direction problem — and 185 should not expect a
regulon-derived `rank_score` to be worth delivering.

---

## 1. Phase 0 — the instrument, and every known-answer row

| | |
|---|---|
| surface | h1-test, a15 emission (`experiments/164-callset-artefact/de164.a15.npz`, 85 × 10,779) |
| effect | `experiments/169-sharpening-six-member/work/h1test-effect.npz` (`H`, `b`, `cpm`, `gate`) |
| endpoint | `experiments/153-de-bench/debench.py:Surface`, unmodified; `experiments/172-sign-confidence-key/harness172.py` for `reach` / `topx` / `deciles` |
| pds instrument | `experiments/162-retrodicting-validator/plolo.py`, fold `orion_hct116`, held-safe operator `held_orion_hct116.npz` |
| `n_real` | never rescaled by `Surface(n_real_target=…)`; the `fid` table sweeps assumed `n_real` analytically instead (§5) |
| donors | none opened. The frozen a15 h1-test table and 169's effect are the only inputs; `data/donors/h1_2025` is never touched. |

Every script calls `harness172.check()` and aborts on any failure. Reproduced on every run:

| row | here | expected |
|---|---:|---:|
| `reach(production order)` | **0.18283** | 0.18283 ✓ |
| `reach(\|H\|)` | **0.23533** | 0.23533 ✓ |
| `reach(sig-first, then \|H\|)` | **0.21914** | 0.21914 ✓ |
| `reach(ORACLE)` | **0.64553** | 0.64553 ✓ |
| INBOX-4 §1 decile table, `\|H\|` and expression | both rows | ✓ (≤0.002) |
| `reach(midrank(\|H\|)) == reach(\|H\|)` | 0.23533 | order-preserving ✓ |
| **174 §1**: `\|r\|` median at pc0 | **0.0256** | 0.0256 ✓ |
| **174 §1**: pool sign accuracy of `−sign(r)` | **0.5158** | 0.5158 ✓ |
| **174 §1**: pool sign accuracy of `+sign(r)` | **0.4842** | 0.4842 ✓ |
| **162 HANDOFF §2**: PLOLO `library/delivered` a15 on `orion_hct116` | **0.789536** | 0.789536 ✓ |

**The 174 rows earned their keep immediately.** The first draft of this experiment stored the
co-expression matrix as `+r` and fed it to the agreement score as if it were the network's
expected knockdown direction. The sign accuracy came out at **0.4842** — exactly 174's `+sign(r)`
row — which is how the flipped convention was caught before a single conclusion was drawn from
it. Rule 10's corollary, working as designed.

For the H-oracle ceiling this experiment is measured against, 172 §4.2 established **0.61862**
(design signs perfectly ordered), not the full oracle 0.64553; 96% of the headroom is a property
of `H` itself and 5.5% of pool signs are corrupted by emission.

---

## 2. Phase 1 — the two priors, their licence, and their coverage

### 2.1 Construction

Both networks are reduced to one object so that every downstream measurement is literally the
same code:

```
N[t, g]   sign = the direction the NETWORK expects gene g to move when target t is KNOCKED DOWN
          |·|  = the network's evidence weight for that pair
          0    = the network says nothing about this pair
```

and the **per-(target, gene) regulon-agreement score**, written down in `key190.py`'s docstring
before anything was scored, is

```
R(t, g) = w(t, g) · sign(N[t, g]) · sign(H[t, g])     if N ≠ 0 and H ≠ 0, else 0
w(t, g) = |N[t, g]| / median(|N| over its nonzeros)
key_λ(t, g) = |H[t, g]| · exp(λ · R(t, g))            λ = 0 IS |H|, the known-answer row
```

`R > 0` = the network confirms the sign we already predict; `R < 0` = it contradicts it. A
rank-additive variant `midrank(|H|) + λ·R` is reported alongside because the multiplicative form
is sensitive to `|H|`'s scale; it reduces to the same λ = 0 row.

Four families of `N`:

| id | network | expected sign |
|---|---|---|
| `ct-direct` | **CollecTRI**, direct edges only — the knocked-down gene must itself be a CollecTRI TF | `−mor` |
| `ct-prop-k{5,25,∞}` | **CollecTRI propagated** through our own predicted TF expression change — `Prop[t,g] = Σ_f a[t,f]·mor(f,g)·w(f,g)` with `a` the per-target z-scored `H` over TF columns, truncated to the top-`k` TFs. **This is the shape `tydymy` published** ("signed TF to target correction"), and it covers every target, TF or not. | sign of the propagated sum |
| `coexpr-pc{0,2}`, `coexpr-pc0-top500` | **control-derived co-expression GRN**, D²R²-shaped: Pearson `r(t,g)` across h1-test's 18,400 NTC cells on log1p-CPM, optionally after projecting out the top-`k` PCs (174: 2 PCs carry the entire raw signal) | `−sign(r)` |
| `expr` | control expression alone — 189's comparator, unsigned, enters only as a key | — |

`w` for CollecTRI is `log1p(#contributing resources)`; the file carries no per-edge confidence
level (unlike DoRothEA's A–E). An unweighted variant was run throughout; at its best λ it reads
0.23538, identical to the weighted one to five decimals (`work/key190.log`).

### 2.2 Provenance and licence

| file | source | rows | sha256 |
|---|---|---:|---|
| `data/external/regulon/collectri_rescued.csv` | `https://rescued.omnipathdb.org/CollecTRI.csv` | 43,536 | `86c90b30…85bcb4e0` |
| `data/external/regulon/collectri_omnipath.tsv` | OmniPath `interactions?datasets=collectri&…&license=academic` | 64,515 | `2cb6762d…168078ed` |
| `data/external/regulon/collectri_commercial.tsv` | the same query with `license=commercial` | 64,515 | `2cb6762d…168078ed` |

CollecTRI as used here: **43,536 signed edges, 1,189 TFs, 6,838 target genes**, 37,662
activating / 5,874 repressing, `TF.category` 41,559 `dbTF` / 1,972 `coTF` / 5 `GTF`. (The
mandate's 45,856 / 1,183 is the pre-deduplication count; the difference is protein-complex
pseudo-sources and duplicate (tf, gene) rows collapsed by consensus sign.)

**Licence — checked, and it is clean.** CollecTRI's repository is GPL-3.0 and states the derived
regulons are freely available. The operative test is OmniPath's own licence filter, which exists
precisely to separate academic-only from commercially usable resources: the `license=academic` and
`license=commercial` queries return **byte-identical files (same sha256, 64,515 rows each)**, so
OmniPath classifies **no** CollecTRI edge as academic-only. The result below makes this moot, but
it is recorded so a successor does not re-litigate it.

### 2.3 Coverage — and this is the number that caps the direct version

| network | pool pairs it speaks about | % of the 65,044-pair DE pool | targets covered |
|---|---:|---:|---:|
| `ct-direct` | **116** | **0.18 %** | **17 / 85** |
| `ct-prop-k5` | 1,865 | 2.87 % | 85 / 85 |
| `ct-prop-k25` | 4,500 | 6.92 % | 85 / 85 |
| `ct-prop-∞` | 24,645 | 37.89 % | 85 / 85 |
| `coexpr-pc0-top500` | 3,966 | 6.10 % | 85 / 85 |
| `coexpr-pc0/pc2` | 65,044 | 100.00 % | 85 / 85 |

On h1-test **17 / 85 targets (20.0 %)** are CollecTRI TFs with at least one edge on the gated
axis; 463 edges are placed and the median regulon is **6 genes**. On the real panel, **41 / 300
targets (13.7 %)** are CollecTRI TFs; 1,436 edges land on the 18,533-gene axis, median regulon
**9**. So a direct-regulon arm's board effect is **multiplied by 0.137**, not by 0.20 —
per-covered-target numbers must not be quoted without it (§6.2).

### 2.4 The mechanism, stated once: the network's sign is worse than ours where it speaks

| network | its sign accuracy on the pool pairs it covers | **`\|H\|`'s accuracy on those same pairs** | `\|H\|` on the whole pool |
|---|---:|---:|---:|
| `ct-direct` | 0.5862 | **0.6379** | 0.5373 |
| `ct-prop-k25` | 0.5011 | **0.5336** | 0.5373 |
| `ct-prop-∞` | 0.5103 | **0.5364** | 0.5373 |
| `coexpr-pc0` | 0.5158 | **0.5373** | 0.5373 |
| `coexpr-pc0-top500` | 0.5229 | **0.5444** | 0.5373 |
| `coexpr-pc2` | 0.4892 | **0.5373** | 0.5373 |

CollecTRI does carry real biology — 0.5862 against a 0.5373 pool base is not nothing, and it is
the only row where the network beats chance by a useful margin. **But `|H|` is already right
0.6379 of the time on those exact 116 pairs.** A confidence key can only reorder; it cannot
import a sign. Asking a less accurate signal to re-rank a more accurate one is the whole of the
negative result, and it is why coverage is only half the story: the co-expression GRN covers
100 % of pairs and is *also* strictly worse than `|H|` everywhere.

---

## 3. Phase 2 — the key sweep against `|H|` = 0.23533

`reach(|H| · exp(λ·R))`. Full 13-point sweep and the rank-additive variant in `work/run190.log`;
the row that matters:

| network | best λ | **reach** | raw Δ | **overall Δ** | paired-bootstrap 95 % CI | P(better) |
|---|---:|---:|---:|---:|---|---:|
| `\|H\|` — **the benchmark** | — | **0.23533** | — | — | — | — |
| `ct-direct` | +1.0 | 0.23538 | +0.00004 | +0.00001 | [−0.00001, +0.00012] | 0.896 |
| `ct-prop-self` | +1.0 | 0.23538 | +0.00004 | +0.00001 | [−0.00001, +0.00012] | 0.896 |
| `ct-prop-k5` | 0 | 0.23533 | 0.00000 | 0.00000 | — | — |
| **`ct-prop-k25`** | **−0.25** | **0.23670** | **+0.00137** | **+0.00026** | **[−0.00140, +0.00520]** | 0.782 |
| `ct-prop-∞` | 0 | 0.23533 | 0.00000 | 0.00000 | — | — |
| `coexpr-pc0` | 0 | 0.23533 | 0.00000 | 0.00000 | — | — |
| `coexpr-pc2` | 0 | 0.23533 | 0.00000 | 0.00000 | — | — |
| `coexpr-pc0-top500` | −0.25 | 0.23645 | +0.00112 | +0.00021 | [−0.00038, +0.00364] | 0.829 |

Four things to read off it:

1. **Four of the eight networks select λ = 0**, i.e. the sweep's own optimum is "do not use the
   prior at all". That is the identity row and it is the honest answer for `ct-prop-k5`,
   `ct-prop-∞`, `coexpr-pc0` and `coexpr-pc2`.
2. **The two non-zero winners sit at NEGATIVE λ**, which means `reach` improves fractionally when
   you *down-weight* genes where the network agrees with your predicted sign. An anti-prior. That
   is not a mechanism anyone would ship; it is the fingerprint of taking a max over 13 λ.
3. The best magnitude is **7.3× below** the composed-slot-gate's shape-only floor of **+0.0019**,
   and the rest are indistinguishable from zero.
4. `ct-direct`'s +0.00004 is not even close to its own oracle ceiling of +0.00078 — the direct
   regulon is not "under-exploited", it is empty.

### 3.1 The ceiling — the number that actually closes the question

Give an **oracle** control of every pair the network covers (promote the truly-correct ones to the
front, demote the truly-wrong ones to the back, leave every uncovered pair at `|H|`). No prior
built on that network can beat this.

| network | ceiling `reach` | raw Δ | **overall Δ** | achieved / ceiling |
|---|---:|---:|---:|---:|
| **`ct-direct`** | **0.23612** | **+0.00078** | **+0.00015** | 5 % |
| `ct-prop-k5` | 0.24659 | +0.01126 | +0.00212 | 0 % |
| `ct-prop-k25` | 0.28403 | +0.04870 | +0.00917 | 2.8 % |
| `ct-prop-∞` | 0.44835 | +0.21301 | +0.04011 | 0 % |
| `coexpr-pc0-top500` | 0.27471 | +0.03938 | +0.00741 | 2.8 % |
| `coexpr-pc0/pc2` (100 % coverage) | 0.64553 | +0.41019 | +0.07724 | 0 % — this row is just the oracle |

**Two different failures, and they must not be conflated.** For `ct-direct` the ceiling itself is
+0.00015 overall, so coverage alone closes it — no amount of cleverness helps. For everything else
the ceiling is comfortably shippable (+0.0021 to +0.0401 overall) and the achieved value is
**zero to 2.8 % of it**, so what is missing is information, not reach.

**189 derived the same cap from the other direction and it agrees.** Its §5 mechanism: `reach`'s
purity bar is 0.90 and the median pure prefix is 3 genes (172 §7.2), so lifting front accuracy from
0.587 to 0.669 never crosses the bar, `k*` does not move, and *"reordering only inside the top
2–10 % of `|H|` is worth ≤ +0.0005 raw `reach`"*. This experiment's oracle ceiling over a network
touching 0.18 % of the pool is **+0.00078 raw**. Same phenomenon, two independent derivations: a
prior confined to a small, front-adjacent slice of pairs cannot move this metric no matter how good
it is.

### 3.2 Front precision — the objective, per 172's calibration

172 established that `reach` reads the front, and that pool rank correlation is *anti*-correlated
with it, so top-X% sign accuracy is the diagnostic. Pooled over (target, gene) pairs:

| key | 1 % | 2 % | 5 % | 10 % | 25 % | 50 % | all |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`\|H\|`** | **0.707** | **0.666** | **0.651** | **0.623** | 0.579 | 0.557 | 0.534 |
| `ct-direct` (λ=1) | 0.704 | 0.665 | 0.651 | 0.622 | 0.579 | 0.557 | 0.534 |
| `ct-prop-k25` (λ=−0.25) | 0.703 | 0.667 | 0.648 | 0.620 | 0.579 | 0.557 | 0.534 |
| `coexpr-pc0-top500` (λ=−0.25) | 0.706 | 0.667 | 0.647 | 0.621 | 0.578 | 0.557 | 0.534 |

**No key improves the front.** The macro variant moves a little for `ct-prop-k25` (0.838 vs
`|H|`'s 0.826 at 1 %, 0.787 vs 0.769 at 10 %), which is the whole of its +0.00137 `reach` — a
handful of targets, and §4.1's null covers it. Full table: `work/topx190.json`.

---

## 4. Falsifiers

### 4.1 Falsifier 1 — the shuffled-regulon control `tydymy` did not report. **This one decides it.**

Three shuffles, 20 seeds each, with **best-λ re-selected per seed exactly as for the real
network**, so the selection artefact is inside the null where it belongs.

| network | real best | shuffle = targets | shuffle = signs | shuffle = tfs |
|---|---:|---|---|---|
| `ct-direct` | 0.23538 | 0.23536 ± 0.00003 · p=0.20 | 0.23539 ± 0.00004 · p=0.60 | 0.23555 ± 0.00038 · p=0.45 |
| `ct-prop-k5` | 0.23533 | 0.23567 ± 0.00075 · p=1.00 | 0.23535 ± 0.00005 · p=1.00 | 0.23607 ± 0.00083 · p=1.00 |
| **`ct-prop-k25`** | **0.23670** | 0.23622 ± 0.00130 · p=0.25 | **0.23699 ± 0.00059 · p=0.55** | 0.23586 ± 0.00077 · p=0.15 |
| `ct-prop-∞` | 0.23533 | 0.23536 ± 0.00009 · p=1.00 | 0.23533 ± 0.00000 · p=1.00 | 0.23539 ± 0.00024 · p=1.00 |
| `coexpr-pc0` | 0.23533 | 0.23570 ± 0.00163 · p=1.00 | 0.23534 ± 0.00001 · p=1.00 | 0.23645 ± 0.00165 · p=1.00 |
| `coexpr-pc0-top500` | 0.23645 | 0.23612 ± 0.00118 · p=0.25 | 0.23668 ± 0.00027 · p=0.80 | 0.23651 ± 0.00153 · p=0.35 |

`targets` redraws each TF's gene list keeping its out-degree and its (mor, w) multiset; `signs`
permutes `mor` globally keeping topology; `tfs` reassigns whole regulons to other targets.

**Not one p-value is below 0.15, and in eleven of eighteen cells the shuffled network reads
HIGHER than the real one.** The best key in the whole experiment, `ct-prop-k25` at 0.23670, sits
*below* its own sign-shuffled null mean of 0.23699. There is no regulon signal here; there is a
max taken over a 13-point grid, and it is worth about +0.0013 raw `reach` to any signed matrix of
the right shape.

### 4.2 Falsifier 2 — coverage honesty, in overall units

Stated at §2.3 and repeated as the arithmetic: `ct-direct`'s per-covered-target `reach` delta is
**+0.00020 over its 17 targets**, which dilutes to **+0.00004 over all 85** (×0.200) and would
dilute to **×0.137** on the 300-target panel. Even at the oracle ceiling the covered-target effect
is +0.0039 per covered target → +0.00078 over 85 → **+0.00015 overall**. Quoting the
per-covered-target number would inflate the result 5–7×; it is not quoted anywhere else in this
document.

### 4.3 Falsifier 3 — does it beat `|H|` × expression?

The mandate set this comparator from 189's pilot: control expression predicts per-gene cross-line
sign reliability with a held-out **decile spread of +0.0505** (`work/pilot189.json`).
**189's own finished RESULTS.md, filed the same day, retracts that number** — it is an artefact of
the pilot's "above this line's own median |effect|" gate, and under a significance-like gate the
spread goes **+0.0505 → −0.0186 → −0.0457**. So the comparator the mandate named no longer exists
in the form it was named. It was measured here anyway, independently, and the two experiments
agree:

| key | `reach` | raw Δ vs `\|H\|` | overall Δ |
|---|---:|---:|---:|
| `\|H\|` | 0.23533 | — | — |
| expression alone (log CPM) | 0.17589 | −0.05945 | −0.01119 |
| `\|H\| · exp(+0.25·z_expr)` — best of 7 | 0.22573 | **−0.00961** | −0.00181 |
| `\|H\| · exp(−0.25·z_expr)` | 0.22453 | −0.01081 | −0.00204 |

On h1-test the expression decile spread over the DE pool is **+0.0110**, against `|H|`'s
**+0.1223** — an order of magnitude apart. **Every expression-modulated key is worse than `|H|`,
in both directions of μ**, with the best at −0.00961 raw (paired bootstrap SE 0.00567, 95 % CI
[−0.02143, +0.00092], P(better) = 0.039).

**Independent confirmation.** 189 §6 ran the same tilt on the same surface with a different prior
pipeline and reports *"every expression tilt of `|H|` loses `reach` (best −0.01001)"*. This
experiment measured **−0.00961**. Two agents, two codebases, **0.0004 apart**.

So the mandate's pre-registered framing — "the regulon must beat `|H|` × expression, not just beat
random" — resolves the other way round: **`|H|` × expression is not a bar, it is a regression, and
the only bar is `|H|` itself.** The regulon clears the expression baseline by doing nothing at
all, which is not a result; it fails against `|H|`, which is.

For scale: 189's own best key, `|H| · exp(+0.5·θ)` from a per-gene cross-line reliability prior,
reads **0.23941 (+0.00408 ± 0.00424, P = 0.828)** — three times this experiment's best regulon key
and still not separable from zero. **The learned per-gene prior, the regulatory prior and the
expression prior are now all measured on this surface and all three are null.**

### 4.4 Falsifier 4 — retrodicting `tydymy`'s board arm

Their published move (`xfer-orion-genefill-k0-v4` → `regprop-collectri-v1`, 2026-08-26), against a
CollecTRI correction at their published `c` = 0.25 of backbone L2, modelled here as
`pred_lfc' = pred_lfc + c·‖pred_lfc‖₂·unit(correction)` per target:

| member | `tydymy`, board | this model (`ct-prop-k25`, c=0.25) | instrument | verdict |
|---|---:|---:|---|---|
| `nmae` | **+0.0421 WORSE** | **+0.0326 worse** | h1-test `Surface.nmae` | ✓ sign, 77 % of magnitude |
| `mse` (log) | **+0.108 worse** | `sse` 0.999450 → 0.999483; lfc L1 error +3.8 % | PLOLO / h1-test | ✓ direction |
| `pds` | **−0.0013** | **−0.00245** (PLOLO challenge frame) | PLOLO `orion_hct116` | ✓ sign and order of magnitude |
| `reach` | **+0.0353 BETTER** | **−0.0131 worse** | h1-test `Surface.reach` | ✗ **wrong sign** |
| `fid` | **+0.0133 BETTER** | **−0.0004** | h1-test `Surface.fid` | ✗ **absent** |

**The cost side retrodicts on three of three. The benefit side — the two members that are the
entire reason to be interested — fails on two of two, and one of them fails with the wrong sign.**

**Falsifier 4b makes it worse for the hypothesis.** A *shuffled* CollecTRI of identical topology
and identical L2, same c = 0.25, 10 seeds:

| correction | Δpds | Δnmae | Δfid | Δreach |
|---|---:|---:|---:|---:|
| `ct-prop-k25` **real** | +0.00280 | **+0.03261** | −0.00039 | **−0.01309** |
| shuffle = signs | −0.00077 ± 0.00485 | **+0.03129 ± 0.00198** | −0.00036 ± 0.00014 | **−0.01308 ± 0.00462** |
| shuffle = targets | −0.00074 ± 0.00421 | +0.01532 ± 0.00158 | −0.00034 ± 0.00018 | −0.00130 ± 0.00222 |
| shuffle = tfs | −0.00529 ± 0.01078 | +0.02442 ± 0.00363 | −0.00018 ± 0.00022 | −0.01188 ± 0.00903 |

The real network's `nmae` penalty is **0.67 σ** from its sign-shuffled null, its Δ`pds` +0.00280 is
**0.74 σ** above it, and its Δ`reach` −0.01309 is **0.002 σ** from it — the sign-shuffled regulon
reproduces the real one's `reach` damage to five decimals. **The `nmae` bill `tydymy` paid is the
price of adding an off-manifold vector of that L2 to the profile, not the price of adding
CollecTRI's biology.**

Per audit-181 rule 1: **this experiment cannot reproduce the shape of the board point it was
built on, and says so.** Three possibilities remain and none is resolvable from outside their
pipeline — (i) their `reach`/`fid` gain came from something else in that arm, (ii) it is an
interaction with a backbone whose base `reach` was 0.2196 rather than ours, or (iii) the ordering
benefit is real and needs a delivery mechanism this experiment has not found. What is now
established is that **it is not obtainable from the prior as an ordering key on our backbone**,
which is the specific claim M2 made.

---

## 5. `fid` as a function of the call set the same score chooses

`F_p = k / max(n_pred, n_real)`; `n_real` is a distribution (audit 181: board median 90–110, heavy
tail), so it is swept rather than assumed. Precision is the top-`n_pred` sign accuracy over
adjudicable genes; full 5 × 5 grid in `work/comp190.log`.

| key | precision @25 | @50 | @90 | @150 | @300 |
|---|---:|---:|---:|---:|---:|
| `\|H\|` | 0.5412 | 0.5254 | 0.5255 | 0.5204 | 0.5179 |
| `ct-prop-k25` (λ=−0.25) | 0.5459 | 0.5332 | 0.5250 | 0.5182 | 0.5160 |
| `ct-direct` (λ=1) | 0.5431 | 0.5259 | 0.5244 | 0.5194 | 0.5172 |
| expression | 0.5078 | 0.5115 | 0.5102 | 0.5143 | 0.5096 |

Δ`fid` against `|H|` over the entire (n_pred, n_real) grid:

| key | Δ`fid` range | overall |
|---|---|---|
| `ct-prop-k25` | [−0.00220, **+0.00555**] | [−0.00123, +0.00311] |
| `coexpr-pc0-top500` | [−0.00125, +0.00286] | [−0.00070, +0.00160] |
| `ct-direct` | [−0.00105, +0.00067] | [−0.00059, +0.00038] |
| `ct-prop-k5` / `ct-prop-∞` / `coexpr-pc0` / `coexpr-pc2` | exactly 0 | 0 |
| expression (the 189 comparator) | [−0.01529, −0.00278] | all negative |

The one number that looks promising is `ct-prop-k25`'s +0.00555 at `n_pred` = 50, from a precision
of 0.5332 against 0.5254. **It is one sigma.** 85 targets × 50 genes = 4,250 pairs gives a
binomial SE of 0.0077 on a precision difference of +0.0078. It comes from the same key, at the
same λ, that §4.1's null cannot separate from a shuffled regulon, and `fid` was not itself
null-tested — stated as a gap, not smuggled in as a finding.

---

## 6. Phase 4 — the DIRECTION blend, priced separately so the two uses never merge

Phase 4 was conditional on Phase 2 being positive. It is not, but the mandate asked for the
measurement so that 108's use and this experiment's use are cleanly separated in the record.

### 6.1 h1-test, all five members at once (call set held fixed, so only signs and magnitudes move)

`pred_lfc' = pred_lfc + c·‖pred_lfc‖₂·unit(N)`; `jac` is invariant by construction (it reads only
the call set).

| correction | c | `pds` | `nmae` ↓ | `reach(\|H\|)` | `reach(\|pred\|)` | `fid` | sign acc |
|---|---:|---:|---:|---:|---:|---:|---:|
| **shipped a15** | — | **0.79944** | **0.96659** | **0.23533** | **0.22173** | **0.50701** | **0.5343** |
| `ct-direct` | 0.25 | 0.79930 | 0.96685 | 0.23485 | 0.22177 | 0.50700 | 0.5342 |
| `ct-prop-k25` | 0.25 | 0.80224 | **0.99920** | 0.22224 | 0.20455 | 0.50662 | 0.5330 |
| `ct-prop-∞` | 0.25 | 0.80224 | 0.99503 | 0.21035 | 0.18642 | 0.50736 | 0.5342 |
| `coexpr-pc0-top500` | 0.25 | 0.78487 | 0.98913 | 0.23512 | 0.22586 | 0.50754 | 0.5341 |
| `ct-prop-k25` | 1.00 | 0.76261 | 1.19613 | 0.21159 | 0.19120 | 0.50645 | 0.5317 |

Full c ∈ {0.05, 0.10, 0.25, 0.50, 1.00} sweep in `work/comp190.log`. Pool sign accuracy **falls
with c on every network past c = 0.10, and no setting on any network raises it by more than
0.0003** — the correction does not repair directions, it corrupts them.

### 6.2 PLOLO `pds`, the mandated instrument

Fold `orion_hct116`, held-safe operator, a15's exact shipped config, correction built from that
fold's **own** whitened effect on the 300-panel axis (41/290 fold targets are CollecTRI TFs, 1,168
TFs measurable). Known-answer row `library/delivered` a15 = **0.789536**, reproduced exactly.
162's gate map: predicted overall = Δpds/6 × 1.011.

| correction | c | Δpds (intended) | Δpds (challenge frame A) | → overall |
|---|---:|---:|---:|---:|
| `ct-direct` | 0.10 | +0.000012 | +0.000155 | +0.000026 |
| `ct-direct` | 0.25 | +0.000024 | +0.000143 | +0.000024 |
| `ct-prop-k25` | 0.25 | **−0.003878** | **−0.002446** | **−0.000412** |
| `ct-prop-k25` | 0.50 | −0.013006 | −0.008018 | −0.001351 |
| `ct-prop-∞` | 0.25 | −0.008185 | +0.001420 | +0.000239 |
| `ct-prop-∞` | 0.50 | −0.021501 | −0.001599 | −0.000269 |

`ct-direct` is flat to five decimals. `ct-prop-k25` — the configuration that matches `tydymy`'s
description — is **negative and gets worse with c**, and it is the row that reproduces their board
`pds` −0.0013. `ct-prop-∞` is the one row where the two frames **disagree in sign** (−0.0082
effect-level against +0.0014 delivered); at ±0.0002 overall that disagreement is not resolvable,
and 162's own a23 retrodiction over-predicted by 5.5× in the better of the two frames, so nothing
should be built on it. `sse` rises monotonically with c on every network, matching `tydymy`'s
`log mse` +0.108.

**This is the third independent negative on the direction use.** 108 measured it on donor lines
(sign-preserving reweight −0.00035 macro, 10 % direction blend −0.00196 / +0.00054, all four
promotion gates failed); §6.1 measures it on h1-test across five members; §6.2 measures it on the
only `pds` instrument validated against a board point. All three agree.

---

## 7. Answering INBOX 1, and the verdict in overall units

### 7.1 INBOX 1 §3 — the control-derived GRN, given equal weight

INBOX 1 landed after launch and asked for the control-derived GRN to be treated as *"the more
interesting half, not the licence-safe fallback"*, on the grounds that 160/174/155/067 closed
*"control cells cannot predict the line-specific residual"* while **using a control GRN to order or
propagate an effect we already have is a different operation on the same data, and has never been
tested here.** That is correct, it was the right thing to ask for, and it is now tested. Three
variants (`pc0`, `pc2`, `pc0-top500`), the same 13-point λ grid, the same three shuffles, 20
seeds each:

* `coexpr-pc0` and `coexpr-pc2` **select λ = 0** — the sweep's own optimum is `|H|` untouched.
* `coexpr-pc0-top500` reaches 0.23645 (+0.00112 raw, **+0.00021 overall**), at **negative λ**, with
  a target-shuffled null of 0.23612 ± 0.00118 and p = 0.25/0.80/0.35.
* Its expected sign is right **0.5229** of the time on the pairs it speaks about, against `|H|`'s
  **0.5444** on those same pairs.

**So the ordering use of the control GRN is negative for the same reason as everything else in this
document, and it is not a new closure of 174's finding — it is a genuinely different operation that
independently fails.** The one thing it does confirm is 174's diagnosis: `pc2` (cell-state axes
removed) is *worse* than `pc0` on every measurement, so what little signal the raw matrix carries is
the shared proliferation programme, not regulation. This is the answer to INBOX 1 §3's live
question; `cqawesome`'s `E3a_ctrlprog` must be doing something else with the control cells.

INBOX 1 §1's reconciliation task (AgingLab's +0.100 `reach` against 172's null) is **not addressed
here** — it is a different mechanism (a confidence head plus a shrink) and was outside this mandate.
It is still open.

### 7.2 All six members for the arm this experiment considered (rule 12)

The arm is *"replace the `rank_score` with `|H| · exp(λ·R)`"*, best case `ct-prop-k25`, λ = −0.25.

| member | Δ overall | how it is pinned |
|---|---:|---|
| `reach` | **+0.00026** | measured, CI [−0.00026, +0.00098] overall, shuffled null reads higher |
| `fid` | +0.0031 to −0.0012 across the (n_pred, n_real) grid, 1 σ | measured, not null-tested (§5) |
| `jac` | **0.000** | mechanism: a key reorders, it does not change the call set; `Surface.jac` reads only `pred_sig` |
| `pds` | **0.000** | mechanism: `pds_cosine` is a function of the profile alone and no key enters it |
| `mse` | **0.000** | same mechanism; the emitted profile is bit-identical |
| `nmae` | **0.000** | same mechanism; `nmae` is a deterministic functional of the mean profile |

**That zero-cost structure is real and is the good news buried in this negative** — it is exactly
the mandate's thesis that an ordering lever collects `tydymy`'s gross without their `nmae` bill.
The lever simply has nothing to collect. The moment the prior is moved into the *direction* to buy
coverage, all four of those zeros stop being zeros (§6), which is precisely the trade `tydymy`
made.

### 7.3 Verdict

| channel | best measured | ceiling (oracle over covered pairs) | ship floor |
|---|---:|---:|---:|
| **`reach`, as an ordering key** | **+0.00026** (CI crosses 0; shuffled null reads higher) | **+0.00015** for `ct-direct`; achieved 0–2.8 % of the larger ceilings | +0.0019 shape-only |
| **`fid`, same key** | +0.0031 at one grid cell, **1 σ**, not null-tested | — | — |
| **`pds`, direction blend** | −0.00041 at `tydymy`'s own c | — | +0.0089 amplitude-changing |
| **`nmae`, direction blend** | **−0.0090 overall** (raw 0.9666 → 0.9992, i.e. +0.0326 worse) | — | — |

**M2 is closed.** The mandate pre-authorised this outcome — *"if the honest answer is 'the regulon
adds nothing over `|H|` × expression', that is a full result"* — and the measured answer is a
degree stronger and differently shaped: the regulon adds nothing over **`|H|` alone**, and
`|H|` × expression is itself −0.0096 raw `reach` *below* `|H|`, so it was never the right bar.

**What would have to be true for a successor to reopen this.** Not more coverage — `ct-prop-∞`
covers 37.9 % of the pool with a +0.0401 overall ceiling and delivers exactly zero. The binding
constraint is §2.4: **a prior can only help if its sign is more accurate than `|H|`'s on the pairs
it speaks about**, and no version of this network is. A prior that cleared that bar — a
line-specific network with measured directions, or a per-pair confidence that is *conditionally*
better where `|H|` is weak — would be worth the same measurement; a fixed literature network is
not. 174 already showed the co-expression route dies for a related reason (2 PCs carry all of it,
and it is a shared-programme artefact).

For 185: **do not expect a regulon-derived `rank_score`.** `|H|` remains the best key measured
anywhere in this project (172's conclusion, unchanged here). For 189: your expression prior and
this regulon prior are independently negative on `reach`; §4.3 has the h1-test conversion of your
decile statistic, and it is +0.0110, not +0.0505.

---

## 8. Artifacts and reproduction

| file | what it does |
|---|---|
| `prior190.py` | CollecTRI loader, coverage report → `work/coverage190.{json,log}` |
| `coex190.py` | one-off control co-expression GRN, pc0/pc2/pc5 → `work/coex190.npz` (7.43 GiB peak) |
| `nets190.py` | the network library: `ct_direct`, `ct_prop`, `coexpr`, and all three shuffles |
| `key190.py` | `R(t,g)`, the two key forms, the first sweep → `work/key190.{json,log}` |
| `run190.py` | known-answer rows, coverage, ceilings, sweep, shuffled nulls → `work/run190.{json,log}` |
| `comp190.py` | expression comparator, `fid` grid, bootstrap, direction blend → `work/comp190.{json,log}` |
| `final190.py` | per-covered-target dilution, signed Δ`fid`, falsifier 4b → `work/final190.{json,log}` |
| `plolo190.py` | PLOLO `pds` pricing of the direction blend → `work/plolo190.{json,log}` |
| `work/topx190.json` | pooled + macro top-X% for every key |

```
uv run python prior190.py && uv run python coex190.py && uv run python key190.py
uv run python run190.py && uv run python comp190.py && uv run python final190.py
uv run python plolo190.py
```

Total wall clock ≈ 6 minutes. Every script aborts on a failed known-answer row.

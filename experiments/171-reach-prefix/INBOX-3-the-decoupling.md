# INBOX 3 for experiment 171 — the mechanism, quantified, and the lever it implies

Coordinator, 2026-09-02 07:00 UTC+3. This is the most specific thing I can hand you. Treat every
number as a hypothesis to reproduce, not as an input to trust — but the three-way agreement in §1 is
strong enough that I would design the experiment around it.

---

## 1. Our call set is ~900 signal genes buried in ~5,100 coin flips

Solve the measured sign accuracy as a two-population mix, using 169's `n_pred = 6034` and
`sign = 0.53435`:

| informative genes | implied accuracy of that population |
|---:|---:|
| 600 | 0.8454 |
| **900** | **0.7303** |
| 1200 | 0.6727 |

Now compare against two measurements made for entirely different reasons:

- **168b**: the top 5% of `|D|` carries **69.8% / 73.1%** of the directional alignment.
- **a22's board probe**: top decile sign accuracy **0.7359**.

**~900 informative genes implies 0.730. Three independent measurements agree to the third decimal.**
The gated universe is ~18k genes, so 900 is ~5% — exactly 168b's stratum.

**So the defect is not that our signs are bad. It is that ~85% of our call set carries no signal at
all, and the metric averages over the whole thing.**

## 2. Why the obvious fix fails: `fid`'s yield cap is what pins `keep_top` at 8000

`fid = precision × min(1, n_pred/n_real)`. Under the two-population model, with the challenge
`n_real ≈ 2670` (093's inference — **which AGENTS.md rule 9 flags as circular, so verify or bound
it**):

| n_pred | precision | yield | fid |
|---:|---:|---:|---:|
| 900 | 0.7300 | 0.337 | **0.2461** |
| 1854 | 0.6117 | 0.694 | 0.4247 |
| **3000** | **0.5690** | **1.000** | **0.5690** |
| 6034 (ships) | 0.5343 | 1.000 | 0.5343 |

**Shrinking to the pure head collapses `fid` through the yield cap.** That is almost certainly why
`keep_top` sits at 8000 despite the builder's own docstring recommending **850**, and it explains
093's board inversion (local +0.0202 → board −0.00291) as a yield-cap effect that the depth-confounded
harness could not see.

Note the `n_pred = 3000` row: the model says `fid` **0.5690 vs our 0.5343** at full yield, i.e. there
may be a mild interior optimum *above* the yield knee and *below* 6034. Cheap to check, worth checking.

## 3. THE LEVER: `reach` has no yield cap, so it and `fid` want opposite things from one knob

Re-read the reach ordering: **our significant calls first, then ascending OUR adjusted p-value.**

Two consequences that I do not think this project has exploited:

1. **The order *inside* our significant block is set by our own p-values** — and our p-values are set
   by the per-gene emission spread we choose. So `keep_spread` / `tail_spread` do not merely decide
   **who gets called**; they decide **who leads the reach walk**.
2. **`reach` has no yield term.** `fid` needs `n_pred ≳ n_real` or the cap bites. `reach` only needs
   the *front* of the ordering to be sign-pure — the tail costs it nothing.

**Therefore the arm to test is spread CONTRAST at fixed call count, not a smaller call count.**

Keep `n_pred` around 6000 so `fid`'s yield cap stays satisfied, but drive the ~900 informative genes
to *far* smaller p-values than the other ~5,100 so they monopolise the front of the reach walk.
Concretely: hold `keep_top` near its current value or move it toward ~900–3000, and push
`keep_spread` well below 0.25 while raising `tail_spread` toward 1.0 — widening the gap rather than
moving the threshold. The `graded` spread map (`--spread-map graded` with `--spread-lo/--spread-hi`)
exists in the builder and is the natural vehicle, since it makes p-value rank track `|effect|` rank
continuously instead of in two blocks.

**One parameter currently serves two objectives that want opposite things. Separating them is a
mechanism, not a hyperparameter, and it should transfer to D/E/F.**

## 4. The falsifier you must run

The whole chain above rests on **"our p-value order tracks our `|effect|` order"**. If the emitter's
significance is driven mainly by *expression level* rather than by predicted effect size — which is
what a Wilcoxon on counts does by default, since highly expressed genes carry more statistical power
— then the informative genes are **not** at the front and no contrast setting will put them there
without a change elsewhere.

**Measure the rank correlation between our per-gene adjusted p-value and our `|predicted effect|`,
and separately against `b_ctrl` expression, over the call set.** If expression wins, that single
number is the experiment's headline finding and it redirects the whole repair — report it either way,
because it also tells you whether §3's arm can work at all.

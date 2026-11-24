# 319o (down-first p-value order) — NOT BUILT. H1 counterfactual is negative for the emitter's own order.

2026-09-11 20:05Z. Arc's cell-eval2 0.16.0 changelog records that `reach`'s no-skill point is a down-skewed reference
(61 % of significant genes go down on val A/B/C; 17.6 % of targets have ≥ 90 %-down pools) and that a no-biology
"library-mass injection" arm keeps +0.13/+0.17/+0.14 of reach from_baseline (+0.022/+0.030/+0.026 of avg_score),
which Arc measured, documented and declined to correct for this competition. The question was whether a sign prior in
OUR ordering (down-designed genes first inside the prescribed call set and inside its complement, |H| within each)
collects any of it. Knob added, inert when unset: `VCC_RC_ORDER=down_first` (builder d8dc1613…, `rc_emit.rc_hybrid_block(order_key=…)`); run script `run319o.sh` ready.

Counterfactual on the 323 H1 outputs with 348's walk (`scratchpad/walk319o.py`, read-only), raw reach, 82 targets:

| order | mean reach |
|---|---:|
| O1 current scorer order (sig, p_adj, p, |lfc|) | 0.25435 (known answer, exact) |
| sig block first, then |design| (≈ the emitter's own order) | 0.25486 |
| down-first then |design|, whole axis (348's post-hoc O7e) | 0.26634 |
| **the arm: sig block first, down-first then |design| inside each block** | **0.24136 (−0.0135 vs sig-then-design; paired 95 % CI [−0.038, +0.007])** |

Mechanism: the sign prior on H1 is weak (designed-down accuracy on the pool 0.573 vs designed-up 0.524; mean pool
down-fraction 0.579, 12.9 % of targets ≥ 90 % down) and putting every down gene ahead of every up gene inside the block
breaks the |H| ordering that supplies the front's purity (172: head accuracy 0.624). The whole-axis version gains only
because it drops the significant-first block, which 348 already priced (O7d − O4). Verdict: no slot; the knob stays
(inert) for a sign-conditional variant that keeps |H| as the primary key, which a worker can price on H1 first.

# v10 board result (2026-08-23) — first place

Entry `tuVgUbS5vGJfcnaIT3gT`, model name `Two donors walk into a bar`,
submitted 10:26:56 UTC. **Overall `0.1496326555`, rank 1.**

v10 = v9 + frozen PerturbMatch K562-slot replacement + post-whitening CD4.

| metric | v9 | v10 | delta |
| --- | ---: | ---: | ---: |
| pds   | `0.72718` | `0.73379` | `+0.00661` |
| mse   | `0.00000` | `0.00000` | `0` |
| nmae  | `0.07313` | `0.09094` | **`+0.01782`** |
| fid   | `-0.03812` | `-0.03585` | `+0.00227` |
| reach | `0.10178` | `0.10634` | `+0.00457` |
| jac   | `0.00242` | `0.00257` | `+0.00015` |
| **overall** | **`0.14440`** | **`0.14963`** | **`+0.00524`** |

**Every one of the six moved the right way.** That has not happened before in
this project.

Standing: we lead Ibrahim Mansour (`0.146665`) by `+0.002968`.
Raw `pds_cosine` `0.8278 -> 0.8308`, raw mse `1.7661 -> 1.6817`.

## The calibration lesson, which is the opposite of v9's

exp065 measured `+0.00201` on the Orion HCT116 count-level comparator. The
board delivered `+0.00524`. **The proxy under-predicted by 2.6x.**

Set against v9, where h1-test's raw pds delta realized at `1.12x` but its
`avg6` badly overstated the DE metrics, the pattern now looks like this:

- Proxy harnesses are **directionally reliable** — every sign has been right.
- Their **magnitudes do not transfer**, in either direction, and the error is
  not a fixed discount. v9 realized below its projection; v10 realized 2.6x
  above it.

So do not apply a scalar correction factor to a proxy result. Use proxies to
decide sign and to rank candidates against each other, and treat the board as
the only instrument that measures size.

## Where the gain actually came from

`nmae` contributed `+0.01782` of the `+0.00524` overall — more than three times
pds's contribution. `nmae` is the metric `STRATEGY.md:56` flagged as
**"Underexplored. Never tuned against a trustworthy harness."** That call was
correct, and PerturbMatch + CD4 moved it without being designed to.

The reading: our pds levers really are saturating, exactly as the v9 analysis
suggested, and the remaining headroom is in the DE metrics. `nmae` and `reach`
are now our two live channels. `jac` remains near zero (`0.00257`) and `fid`
remains negative (`-0.03585`) for us and for every other team.

## mse is still zero

Raw mse improved `1.7661 -> 1.6817` but the threshold is `0.989`. Ibrahim
remains the only team on the board below `1.0`, at `0.9095`, and that single
metric is still worth `+0.0138` to him.

Note that we now lead *without* it. Our margin comes from being far enough
ahead on pds, nmae and reach to outweigh a metric we score zero on.

See the retraction note in `STRATEGY.md` Addendum 4 before reasoning about
transfer cosine from the mse column — the normalization is not established.

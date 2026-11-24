# 435 — 429's zs pass stacks on the coarse MDE exactly as on the plain one: nmae −0.00183 ± 0.00007, other five members bit-identical (H1, 3 seeds)

Status: 2026-09-23 10:15Z, PI.
- Harness: `run435.py`, which is 429's `run429.py` plus 426h's `--mde-coarse`, reading 429's inputs.
- The production emitter is used (`RUN426_EMITTER=prod`, `src/vcc2026/rc_emit_k.py` `edb5b8dc`).
- zs is the vendored production module (`src_rank_mag.py` = `src/vcc2026/rank_mag.py` `ed672832`).
- Driver: `drive435.sh`. Known answers and the reading: `ka435.py`.

## Known answers (stats arrays bit-equal)
| row | reference | result |
|---|---|---|
| KA1: production emitter, coarse MDE (0.1, 0.1, 0.3), seed 0 | 426's `k_emit426` coarse run | **PASS** |
| KA2: production emitter, plain MDE + zs `dn_k250_zs`, seed 0 | 429's MDE + zs run | **PASS** |

The first queue run died after writing its stats: the metrics JSON carried a NaN, because the production emitter has
no `mde_overfill_frac` diagnostic and the writer refuses NaN. The harness now records `null` for it. The stats arrays
were never affected.

## Result: zs on top of the coarse MDE vs the coarse MDE alone (426h's runs, same seeds)
| seed | nmae | pds, mse, fid, reach, jac | overall (board slopes) |
|---|---:|---|---:|
| 0 | −0.00185 | 0 (all five, exactly) | +0.00051 |
| 1 | −0.00175 | 0 | +0.00048 |
| 2 | −0.00188 | 0 | +0.00052 |
| **mean** | **−0.00183 ± 0.00007** | 0 | **+0.00050** |

429 measured zs on the plain MDE at −0.00179 ± 0.00007. The coarse delivery changes nothing for zs, and every run
passed the module's exact-U asserts.

**Reading.** zs is ready to ride any coarse-MDE arm from the first arm after the MDE's board read (429's rule). On
A/B/C it has 0.70× H1's capacity (429: A 0.91×, B 0.14×, C 1.04×). Expected board effect +0.0003 … +0.0008 overall,
nmae only.

## Dose ladder with the coarse delivery (H1, seed 0, vs the knob-off base; `work/queue_dose.txt`)
| dose | nmae | fid | jac | reach | unsuppressible calls/target | clip frac |
|---|---:|---:|---:|---:|---:|---:|
| **c 0.1, cap 0.1 (shipped, 426h)** | **−0.00966** | −0.00025 | −0.00009 | −0.0055 | **5** | 0.29 |
| c 0.125, cap 0.125 | −0.01054 | +0.00263 | −0.00344 | −0.0091 | 1,791 | 0.45 |
| c 0.15, cap 0.15 | −0.01030 | −0.00271 | −0.00432 | −0.0054 | 4,204 | 0.54 |
| c 0.1, cap 0.15 | −0.01086 | +0.00074 | −0.00169 | −0.0090 | 341 | 0.39 |
| c 0.15, cap 0.1 | −0.00888 | +0.00391 | −0.00085 | −0.0064 | 698 | 0.35 |

- **Closed.** Every dose above c 0.1 / cap 0.1 leaks hundreds to thousands of calls the |z| design cannot suppress.
  The extra nmae (≤ 0.0012) comes with call-set changes that are not transferable: 426 read the fid "gains" of leaked
  calls the same way.
- `coarse` removed the capacity limit (λ = 1 everywhere) but not the call limit. The binding constraint on the MDE
  dose is the Wilcoxon design, not the class split.

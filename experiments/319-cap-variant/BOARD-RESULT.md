# 319 — published result, reconstructed September 9, 2026

Entry `xuDelEik7zGWMXdNu0R4`, submitted **September 8, 05:08:50 UTC**,
published **0.21865973944346428**. Original full capture:
`data/board-entries/xuDelEik7zGWMXdNu0R4.json`. September 9 live CLI check:
published, rank 5. The earlier `RESULTS.md` is the **pre-submission build
handoff**; its “nothing was submitted” statements describe that earlier stage.

Package: `submissions/319-normrel15-cap.vcc`,
SHA-256 `57585d558f9777f066c02f5b299cacccca0dd4954b4201b8fa653bd0f6816dd2`.
Content digest under the novelty algorithm:
`359b446f3165b3e10133a260cf5f83aa39a8f5286a800277a221b9fd3729c461`.
Both reproduced in the September 9 storage audit.

| member | published scaled value | change vs 312 plain −1.5 |
|---|---:|---:|
| pds | 0.788313922524 | −0.002075105614 |
| mse | 0.181012811502 | +0.000881409059 |
| nmae | 0.107514505975 | +0.004747681025 |
| fid | −0.015759118582 | +0.000191135953 |
| reach | 0.222191666222 | +0.001932889902 |
| jac | 0.028684649020 | −0.000032882683 |
| **overall** | **0.218659739443** | **+0.000940854607** |

The preregistered ≥+0.0005 success bin fired. The measured conclusion is that
the **cap construction at this dose** beats plain on this panel. The score alone
does not prove that suppression removed noise: the cap also changes effective
dose, which donors carry each row, and downstream direction/amplitude.
`RESULTS.md` §6 already measures this coupling.

319d3 (cap −3.0) remains **unsubmitted**. It is a useful candidate for testing
whether cap shares plain's observed turnover. Neither outcome alone identifies
small-row noise as the causal mechanism. The statement in its frozen
preregistration that 318 has already excluded a magnitude explanation is too
strong; see `COORDINATOR-REVIEW-2026-09-09.md`. Do not edit the frozen bins after
the fact, and do not schedule a submission from this note.

# 319p2 board result — partial-rows dose 2 (G 0.34 → 0.5) on the 319p composition

**Entry `lqWNXX5YaYArdglu6Stn`, model "the photocopier demands a raise", 2026-09-11 slot 1** (re-armed at
19:10Z after the 2026-09-10 19:56Z reboot killed the original waiter, pid 55426, fifteen minutes after it was armed;
submitted 19:11Z, published 19:41Z). Package `submissions/319p2-cap15-partial50.vcc`, sha256 `a3f770aa9f51d56c…`.
Record: `data/board-entries/lqWNXX5YaYArdglu6Stn.json`. Reference: 319p `SvzrTOQJAAaNjawzSyrC` 0.2198190.

| | overall | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|---:|
| 319p (G 0.34) | 0.2198190 | 0.78807 | 0.18527 | 0.10813 | −0.01517 | 0.22392 | 0.02871 |
| **319p2 (G 0.5)** | **0.2201159** | 0.78708 | 0.18730 | 0.10831 | −0.01503 | 0.22435 | 0.02869 |
| Δ | **+0.000297** | −0.00099 | +0.00204 | +0.00018 | +0.00014 | +0.00043 | −0.00002 |

**Reading (bins frozen in `PREREGISTRATION-partial-dose2.md`, commit 7c88c8ee): bin 2, flat (−0.0005 … +0.0005).**
The partial-rows residual is saturated between G 0.34 and 0.5: mse carries +0.002 (raw 0.81183 → 0.80990) and pds
gives back −0.001. No third dose. 319p2 is the best-scored package (0.2201159, displayed row, rank 15 of 893 at
19:41Z) and the October composition may carry either dose; the frozen rule says ship the lower dose (0.34) since
the difference is inside the band.

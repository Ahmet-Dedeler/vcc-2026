# 433h BOARD RESULT — 0.2631120, NEW BEST (+0.0094597 vs 433a); the zero repair was worth 3× its top bin

Entry `wMoDAi1mca08FA2Bfqkd`, submitted 2026-09-25 00:01Z by `arm0925_s1.sh` ("the flow cytometer took up jazz"),
scored and captured 00:35Z, rank 7 at publication. Package `433h-…-mde010coarse-zr5-ctxlevel.vcc` sha256 `e7495644…`.
Reference 433a `FCbBcvGxF1HOFPxPISUC` 0.2536523. Bins: `PREREGISTRATION-0925.md` (written 09-23 11:25Z).
Reader: `read0925.py`.

| member | Δ raw | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.000476 | +0.001126 | +0.000188 |
| mse | **−0.013419** | +0.014167 | **+0.002361** |
| nmae | **−0.024748** | +0.041921 | **+0.006987** |
| fid | +0.000032 | +0.000111 | +0.000019 |
| reach | −0.000561 | −0.000577 | −0.000096 |
| jac | +0.000004 | +0.000009 | +0.000002 |
| **overall** | | | **+0.0094597** |

## Against the bins
- **nmae (the zero repair):** Δ −0.0247 raw is in the "repair mattered" bin (≤ −0.002), 3× beyond the preregistered
  range (0 … −0.008).
  - Mechanism: a gated gene delivered at lfc ≈ −35 costs |−35 − l| / Σ_gate |l| in its target's nmae. One such gene can add
    ~0.5 to a target with a 100-gene gate. So a handful of the 24 zero-total genes dominated whole targets on A/B/C.
  - 428's H1 price assumed the H1 gate hit rate; the A/B/C gates held the zeros more often.
- **mse (the re-level):** Δ −0.0134 raw is in the "transfers" bin (≤ −0.006). The cliff model behind 417h's factors holds.
  pds moved +0.0005 (predicted < 0.001).
- The other DE members moved inside float noise: the repair touched only the zero genes.
- **Overall +0.0095**, above the predicted +0.001 … +0.005, entirely through nmae.

## Consequences
- 433h is the incumbent. 433k (433h + q 1.5 levelled) was submitted at 00:35Z by the slot-2 rule (nmae ≤ +0.001 and
  |Δpds| ≤ 0.003 vs 433a: PASS).
- **Follow-up (queued):** repaired packages still deliver 17 gated genes with lfc < −5 and 3 with lfc > +3 (428 audit
  lines). If those are trans (not cis) genes in reference gates with modest true effects, they cost the same way the
  zeros did. Audit their identity and cis status and price a cap on the delivered |lfc|.

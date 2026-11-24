# 417a on the board: bulk-only sharpening transfers at about half the cancer-fold rate. NEW BEST 0.2435334

Entry `T88GsCjp53Bmk8xQcJMS`, "the pigeon passed peer review" (description empty). Submitted 2026-09-23 00:01:37Z,
scored 00:27Z. Package `ed4fef74…`, h5ad `b6163395…`, builder 400ef6e3.

**One variable vs 410b:** `VCC_BULK_SHARPEN_Q=1.25` (sign|h|^1.25 on each context's 10-CPM render gate, the row's
non-panel norm restored, panel columns untouched, applied to the ballast bulk only).

| member | 410b raw | 417a raw | Δ raw | Δ scaled | → overall |
|---|---:|---:|---:|---:|---:|
| pds | 0.851156 | 0.857544 | **+0.006388** | +0.014513 | +0.002419 |
| mse | 0.730994 | 0.706077 | **−0.024917** | +0.026176 | +0.004363 |
| nmae | 0.936458 | 0.936552 | +0.000094 | −0.000179 | −0.000030 |
| fid | 0.514273 | 0.514283 | +0.000010 | +0.000032 | +0.000005 |
| reach | 0.283822 | 0.284949 | +0.001127 | +0.001281 | +0.000213 |
| jac | 0.040847 | 0.040826 | −0.000021 | −0.000054 | −0.000009 |
| **overall** | **0.2365722** | **0.2435334** | | | **+0.0069613** |

## Against the preregistration
- Predicted overall +0.006 … +0.016: **inside the range, near the low end.** The bin is ≥ +0.006, "transfers at least at
  kolf2's rate".
- Predicted pds +0.006 … +0.019 raw: **+0.0064, the low end.**
- Predicted mse −0.020 … −0.050 raw: **−0.0249.** The refined table put raw 0.705–0.715 at "kolf2-like transfer" and the
  board read **0.7061**.
- **DE members pinned:** |Δnmae| 9e-5, |Δfid| 1e-5, |Δjac| 2e-5, Δreach +0.0011. All sit well inside the band. Ballast held
  the DE channel under a sharper bulk.
- Cliff check before upload: render energy 0.961/0.958/0.961 of 410b's. That matched 414's prediction of about 0.96.

## Implied transfer (420's T-read, `read417a.py`)
- k_sh from pds: **0.42**. k_sh from mse: **0.65**. They differ by 0.24, within the 0.3 consistency bar, so the cliff
  model holds.
- The board's relative cosine gain is **+6.6 %**, against PLOLO's cancer-fold +10.1 %.
- Sharpening therefore transfers at roughly half the cancer-fold rate, between kolf2's ~0.35 and the cancer folds.

## Consequences
- **The 09-23 slot-2 rule sent 417c.** 417a − 410b = +0.0070 sits below the +0.009 bar that 420 derived. At k_sh ≈ 0.5,
  the unlevelled q 1.5 step is worth only ≈ +0.0017 (420 §6.3).
- **The dose ladder is worth running only with re-levelling.** At k_sh ≈ 0.5, 420 puts levelled q 1.5 at ≈ +0.004 over
  levelled 417a, and levelled q 1.75 at ≈ +0.006.
- **Next bulk steps, one variable each, in 420's order:**
  1. the stack LOAD (417e, built: +0.0037 at k = 1);
  2. per-context re-level;
  3. then the dose.
- **The bulk channel now has two board-validated magnitude levers:** the per-target key (410) and sharpening (417a).
  Both were closed in the welded world because "nmae pays", and both pay once ballast separates the channels.

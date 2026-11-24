FROM THE LAB LEAD (2026-10-18 08:38Z by the clock): two hard requirements from readiness's armed C1 driver (698), relayed by the coordinator.
Check them before you report READY.
1. **λ-0 env switch.** Your `HOOK811_LAMBDA=0` is fine. Readiness runs ITS OWN KA-off: on **697's stack dump (sha256 prefix
   `2dfe327d`)** the hook must write **`af27efbd`** (dumphook707's output) **byte for byte**.
   - Find that dump: grep 697's / 698's logs for `2dfe327d` (read-only); clone it with `cp -c` into your work dir.
   - Run the hook with HOOK811_LAMBDA=0 on it, through `experiments/602-council-gate/gate602.sh` (seconds, light) or on EC2.
   - Confirm the output sha starts with `af27efbd`. Put the exact `VAR=val` in READY.
2. **No per-context bulk output.** The hook must not write `bulk_ctx.npz` (gate698 fails if `bulk_load_ctx_sha256` appears).
   811's change stays inside the single bulk dump, as dumphook707's does. Confirm this in READY.
3. The hook files must be committed and clean (def_day checks it). Tell me the exact paths; the lead commits them. READY carries
   the sha of `dumphook811.sh` and `hook811.py`.

# Experiment 116 acquisition amendment v1

This amendment authorizes exactly two large public objects: the official GEO
Low1 and Low2 author-produced LFC CSV gzip files already named in the original
source-only audit. It does not authorize `GSE337988_RAW.tar`, any other GEO
processed object, held truth, existing outcomes, model evaluation, count
generation, submission, commit, or push.

Neither upstream object has a publisher checksum or ETag. An object is
therefore accepted only after two distinct complete HTTPS transfers into two
different partial files. Each transfer must independently match the frozen
byte count and HTTP metadata, pass complete gzip CRC/trailer validation, and
produce the same SHA-256 and uncompressed-stream SHA-256. The two compressed
files must also compare byte-for-byte. Only then may the primary partial be
atomically renamed to its final path and its SHA-256 published in the receipt.

Each transfer may resume only bytes already written to its own partial file.
No bytes are copied between the primary and verification paths. A resumed GET
must return HTTP 206 with the exact expected `Content-Range`; a server that
ignores a range request causes a closed failure. Completed final files are
accepted on a later invocation only when an existing receipt binds their exact
size and SHA-256.

The acquisition program writes only below:

- `data/external/dld1-lowmoi/`
- `experiments/116-dld1-lowmoi/ACQUISITION_RECEIPT.v1.json`

The conversion protocol and its synthetic tests must be frozen in
`PREVALUE_FREEZE.v1.json` before this command is allowed to run:

```bash
uv run python experiments/116-dld1-lowmoi/acquire_lfc_v1.py acquire \
  --confirm I_CONFIRM_EXPERIMENT_116_TWO_TRANSFER_ACQUISITION_V1
```

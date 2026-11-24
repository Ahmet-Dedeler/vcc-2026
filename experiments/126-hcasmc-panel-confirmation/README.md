# Experiment 126 — HCASMC panel-confirmation acquisition

Status: **opaque acquisition complete and byte-verified; effect values remain
permanently unopened pending new explicit user authorization**.

This experiment owns only:

- `experiments/126-hcasmc-panel-confirmation/`
- `data/external/hcasmc-igvf/`

Its sole purpose is to acquire and byte-verify the official IGVF HCASMC source
for a possible future independent panel-response confirmation. It does not
define, fit, select, or evaluate a confirmation candidate.

The effect asset `IGVFFI0206LUDV.tsv.gz` must remain an opaque compressed
binary. No command in this experiment may decompress it, inspect its header,
read tabular rows, test its gzip stream, or derive any value-dependent
statistic until a new, hash-bound authorization is created after the
Experiment 125 candidate is frozen.

See:

- `ACQUISITION_PROTOCOL.v1.json`
- `SOURCE_MANIFEST.expected.v1.json`
- `VALUE_ACCESS_BOUNDARY.v1.json`
- `PRE_DOWNLOAD_FREEZE.v1.json`
- `ACQUISITION_RECEIPT.v1.json`

The final gzip was published only after the temporary download had exactly
1,531,381,676 bytes and MD5 `4d3067b3e934842c88b9ea52d74780f9`.
Its local compressed SHA-256 is
`b87452c2db1f508f5013bc15dff8088828ab931b35492d77ed6aea569fb5d9b9`.

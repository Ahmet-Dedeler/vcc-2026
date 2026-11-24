# 335 — local biomedical encoder passes engineering gates only

September 9, 2026. **ENGINEERING_FEASIBLE_ONLY**. Freeze **50895633**.
No response-label fit, competition score, package or submission.

The exact-revision `NeuML/pubmedbert-base-embeddings` download passed
the publisher's Git-blob/LFS checks. The model card declares Apache-2.0.
Only Safetensors weights and fixed configuration/tokenizer files were
downloaded. Encoding used CPU float32 and local-only loading, with
remote code disabled. No model API, cloud or payment.

All 1,536 description hashes reproduce 334's primary-first exact-key
NCBI/UniProt artifact. There are 1,021 available target and 507 available
response-gene rows, each encoded into 768 dimensions; the eight missing
rows remain zero. All rows finite. No text exceeds 512 tokens.

Four pure identity/pooling fixtures pass. Re-encoding the first eight
available descriptions is bit-exact. Twenty exact output-vector groups
are present, and every group has identical consumed token sequences.
These reflect identical description bodies after removing the gene-name
header, not different inputs silently mapped to one object. No
different-token/output-vector collision was found. No row was merged.

Independent array readback verifies canonical order, shape, dtype,
availability, zero missing rows and SHA. Maximum unit-norm error is
2.3842e-7. Root scoring environment remains NumPy 2.5.2 with no Torch;
335's isolated environment is 696 MB. 331's lock is unchanged.

Encoding ran 15:18:11.254–15:18:44.647 UTC, 33.39 seconds.
Process maximum RSS 1,276,166,144 bytes; sampled tree peak 1.18 GiB.
The process has exited. These are engineering measurements, not
scientific evidence about perturbation transfer.

- Model revision: `b79526d6ef3645e0df4530322e266f24c829f5ef`
- Weight SHA:
  `929ddc16369bb4ff6f8e92a6ee9f5e0748f2bef33989b5617533dd50beebe4ef`
- Input manifest SHA:
  `72acd00f168aaa993c12128e3376a2d6356756c70056381cee6c8c711a136531`
- Feature artifact SHA:
  `0a25c3ab0ad9bbdb5bf2f269b8eb7d83e4876f3f476c0565b2c4380fa67abd41`
- Download receipt SHA:
  `d3cf3be604bbc5cdaee48271d67527b0a386b7355b34ade34499eb9cf116029b`

The separately registered 336 comparison may now test these features
against 334's fixed models and controls. No biological benefit is claimed.

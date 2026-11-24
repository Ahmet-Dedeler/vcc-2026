# 333 recovery 001 — archive basename only, before deserialization

September 9, 2026. Original freeze **63e8a747** remains unchanged.
Download completed and both size and published MD5 match:

- Bytes: 574,395,233.
- MD5: `3f6ce4317e3a0091978ae5cb8fbf05a3`.
- Local SHA-256:
  `6193575dcbd7bbf214c8ca3eb518cb3d13272443f98ecd6c26402411a7ca745e`.

`unzip -l` (archive directory only, no pickle deserialization) revealed that
the actual member is
`GenePT_emebdding_v2/GenePT_gene_protein_embedding_model_3_text.pickle.`
with a **trailing period**. The release description omitted that period.
Its expanded size is 909,175,027 bytes, below the original 2 GiB cap.
The strict original selector would reject it. No conversion was attempted.

Allow exactly that observed basename instead of the originally documented
one, not an arbitrary pickle or fallback match. `convert_recovery333.py`
keeps the original acquisition code and manifest byte-for-byte, verifies
them first, and changes only the in-memory member-name constant. It pins
this addendum, its own code and the verified download receipt separately
before use.

No biological output or embedding vector has been read; no coverage,
dimension, model, resource, source-release or statistical gate changes.
The restricted unpickler stays mandatory. Conversion still waits for the
single heavy-job slot.

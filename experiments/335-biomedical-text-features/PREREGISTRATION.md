# 335 — local biomedical-text representation feasibility, no response fit

September 9, 2026. The distinct question left by 334 is representation:
word-frequency TF-IDF/SVD is not a pretrained semantic encoder.
Before another biological fit, construct an identity-bound feature bank
locally. This is engineering feasibility, not evidence that text improves
perturbation predictions.

Use `NeuML/pubmedbert-base-embeddings` at exact revision
`b79526d6ef3645e0df4530322e266f24c829f5ef`.
Its publisher's model card declares Apache-2.0 and describes
PubMed title/abstract and similar-title sentence-pair fine-tuning.
The model has 768 output dimensions, standard BERT architecture and
attention-mask mean pooling. No remote Python or pickle weights:
download only the fixed allowlist, use `model.safetensors`, verify
published LFS SHA `929ddc16369bb4ff6f8e92a6ee9f5e0748f2bef33989b5617533dd50beebe4ef`
and 437,951,328 bytes. Verify small files against Git blob IDs from the
exact-revision official tree. Do not call an embedding API or pay.

Inputs: EXACTLY 334's 1,024 target + 512 response-gene descriptions,
with its frozen primary-first, exact-symbol NCBI/UniProt fallback.
Recompute bodies and compare EVERY body SHA against 334's feature
receipt. No new gene selection, alias lookup, manual annotation or
pretrained GenePT vector access. Missing descriptions remain zeros.

Use an isolated Python environment, preserving both root scorer and
331 runtime. CPU only, torch float32, evaluation/inference mode,
six threads, batch 8, maximum length 512, deterministic algorithms.
Mean-pool token hidden states with the attention mask, then L2 normalize.
Persist canonical row names, input/body hashes, consumed token IDs'
hashes, available mask and float32 matrix.

Known-answer engineering checks:

- every description SHA matches the previous frozen direct-text artifact;
- tokenizer revision/local files fixed, truncation counts recorded;
- identity, missing-row and attention-mask pooling pure fixtures;
- finite vectors, dimension exactly 768, nonzero for every available row;
- repeat first eight available descriptions in the same batch to within
  1e-6 maximum absolute difference;
- exact duplicate output vectors may occur only when consumed token
  sequences are also identical; record both categories, never merge rows;
- all frozen inputs, runtime and downloaded-file hashes reverified.

Gates must all pass for `ENGINEERING_FEASIBLE_ONLY`. No retrieval, response
accuracy, H1 truth, reserve, count test, package or submission. A biological
comparison needs a separate registration. Failure is recorded without
switching models or token lengths until a passing result appears.

Resource budget: public download <=500 MB, local environment <=3 GiB,
45 minutes encoding, max 8 GiB process-tree RSS, repository memory launch
gate and one heavy job at a time. Start no run too late for sprint handoff.
No cloud, paid model API, production defaults, workers or external write.

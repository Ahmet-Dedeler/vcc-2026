# 333 — public GenePT feature feasibility, not a response experiment

September 9, 2026. The biological experiment 332 is still running and its
accuracy results have not been read. This independent acquisition does not
change its frozen models, features, stopping rules or outcome gates.

## Material and why this is not a renamed STRING experiment

The official GenePT repository `yiqunchen/GenePT` links Zenodo record
**10833191**, published March 18, 2024. The record declares CC-BY-4.0 and
describes pretrained embeddings of NCBI gene summaries and NCBI/UniProt
summaries. These are **text-derived features**, not new perturbation labels.
Attribution: Chen YT and Zou J, *GenePT: A Simple But Hard-to-Beat
Foundation Model for Genes and Cells Built From ChatGPT* (2024).

224 explicitly did not load GenePT; 087/123/100/104 instead used STRING,
and 132 used DepMap. A subsequent pilot would therefore test added
information, with capacity-matched STRING/source-only controls. This
acquisition alone demonstrates no response-prediction value.

Pinned released asset (exact spelling):
`GenePT_emebdding_v2.zip`, **574,395,233 bytes**, published MD5
`3f6ce4317e3a0091978ae5cb8fbf05a3`.
Official content endpoint:
`https://zenodo.org/api/records/10833191/files/GenePT_emebdding_v2.zip/content`.
Save the official metadata receipt, verify byte length and published MD5,
then record local SHA-256 before inspecting archive members.

Use the released `GenePT_gene_protein_embedding_model_3_text.pickle` member
only. Do not generate new embeddings, call a model API, send gene data to an
external service, or read gene descriptions to hand-revise predictions.
The API-model label in a precomputed artifact does not imply an API call.

## Engineering gates

1. Freeze acquisition code and this document in git before downloading.
2. Download public bytes only; no account login/payment. Preserve partial or
   failed files rather than silently replacing them.
3. Verify asset length/MD5, record SHA-256 and archive inventory.
4. Read only the exact expected pickle basename, reject duplicate matches,
   path traversal and expanded member sizes >2 GiB. Do not extract arbitrary
   archive paths.
5. Use a restricted unpickler permitting only NumPy array reconstruction
   classes/functions (and passive built-in containers). Unknown globals
   abort; do not fall back to unrestricted pickle loading.
6. Validate a string-to-finite-fixed-length-numeric-vector mapping, unique
   uppercase symbols and no zero vectors. Report row/dimension counts,
   missingness and coverage of the canonical challenge genes and targets.
   Do not drop a vector based on biological outcome.
7. Save float32 arrays with string metadata in a non-pickle NPZ plus hashes.
   Any later response experiment uses that artifact, not the original pickle.

Feasibility pass: expected artifact integrity/schema, >=285/300 target
coverage and >=70% of the 18,533 canonical response genes. Coverage is a
metadata question. No biological outcome may enter these gates.

## Scope and resource bounds

No challenge/H1/Broad/other perturbation response is read. No reserve,
source default, production model, submission or public post changes.
The small-memory download can overlap 332's CPU fits; **archive loading and
conversion wait until no other heavy job is running**.

<=1 GiB downloaded, <=2 GiB member, <=4 GiB scratch, <=8 GiB process memory.
Download <=20 minutes; conversion/coverage <=10 minutes under the repository
memory gate. Abort on a failed gate. Passing earns only design of a
separately preregistered response pilot; it is not a positive model result.

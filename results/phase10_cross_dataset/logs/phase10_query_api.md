# PHASE 10 query API validation (scvi-tools 1.3.3)

Date (UTC): 2026-09-08T06:10:55.002073+00:00
GPU: NVIDIA GeForce RTX 4060 Laptop GPU

## SCVI
- ArchesMixin: True
- prepare_query_anndata / load_query_data: present
- tiny query test: True
- trainable params after freeze: 14
- adapted model can re-encode source cells: True
- mean source latent drift under adaptation: 0.0000
- device: cuda:0

## TOTALVI
- ArchesMixin: True
- prepare_query_anndata / load_query_data: present
- tiny query test: True
- trainable params after freeze: 34
- query uses protein: True
- adapted model can re-encode source cells: True
- mean source latent drift under adaptation: 0.0018
- device: cuda:0

## Terminology
cross-dataset transfer / query generalization (scArches). Not strict zero-shot: unlabeled target RNA (and protein for totalVI) update new-batch / unfrozen parameters.

Zero-shot secondary experiment: **not added**. load_query_data expands batch embeddings for unseen categories with random padding. Embedding query cells without adaptation would use untrained batch parameters. Mapping query cells as the source batch would ignore dataset identity. Neither is a documented zero-shot mode.

Comparable SCVI vs TOTALVI query workflow: **True**

PHASE 9 was not empirically executed: totalVI 1.3.3 has no honest cell-level missing-protein mask. That is an API limitation, not a failed scientific result.

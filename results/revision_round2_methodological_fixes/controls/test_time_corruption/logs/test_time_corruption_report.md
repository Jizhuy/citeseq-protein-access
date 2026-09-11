# Test-time corruption report

- timestamp: 2026-09-11T08:43:59Z (aggregate); DONE sealed 2026-09-11
- n_rows: 1514
- branch: A (PHASE 11B adapted models)
- venue: RTX 4060 Laptop GPU / CUDA 12.6 / torch 2.14.0+cu126 / scvi 1.3.3
- seeds: Direction A 20/20, Direction B 20/20
- smoke: PASS (scVI invariance max_abs=0.0)
- integrity: param_hash_unchanged 1514/1514; rna_unchanged 1514/1514
- clean baselines: 80/80 CONSISTENT vs Phase 11B

## Clean baselines (mean over completed seeds)

- direction_A scvi_matched: macro-F1=0.6556 (n=20)
- direction_A totalvi: macro-F1=0.6524 (n=20)
- direction_B scvi_matched: macro-F1=0.6015 (n=20)
- direction_B totalvi: macro-F1=0.6540 (n=20)

## totalVI mean macro-F1 by corruption

### direction_A
- p=0.0: F1=0.6524 (Δ=0)
- p=0.1: F1=0.6499 (Δ=−0.0025)
- p=0.25: F1=0.6465 (Δ=−0.0058)
- p=0.4: F1=0.6437 (Δ=−0.0087)
- p=0.55: F1=0.6411 (Δ=−0.0113)
- p=0.7: F1=0.6370 (Δ=−0.0153)
- p=0.85: F1=0.6345 (Δ=−0.0178)

### direction_B
- p=0.0: F1=0.6540 (Δ=0)
- p=0.1: F1=0.6516 (Δ=−0.0023)
- p=0.25: F1=0.6471 (Δ=−0.0069)
- p=0.4: F1=0.6432 (Δ=−0.0108)
- p=0.55: F1=0.6366 (Δ=−0.0173)
- p=0.7: F1=0.6287 (Δ=−0.0253)
- p=0.85: F1=0.6212 (Δ=−0.0328)

## Headline

Gradual monotonic degradation; not catastrophic collapse. Dir B totalVI remains above clean scVI through p=0.85; Dir A totalVI starts slightly below clean scVI. STATE 2 retained. No P2 required by this control.

# RNA-only annotation control (PART B)

## GATE 2

- RNA-only high-confidence cells: see agreement table
- Mapping protein-independent: **YES** (PHASE 6 code audit)
- Agreement with historical labels: **1.0** by construction (same labels)

## Interpretation case

Because development labels were already RNA-only SCANVI transfers:

**CASE: label provenance does NOT explain the development totalVI advantage.**

Protein-informed *Lawlor author labels* remain a separate external concern.

## RNA-only Set A macro-F1

| Rep | F1 |
|---|---:|
| RNA PCA | 0.6665 |
| protein PCA | 0.5120 |
| concat PCA | 0.7449 |
| MOFA+ | 0.5492 |
| scVI | 0.7220 |
| totalVI | 0.7789 |

Delta totalVI−scVI = 0.0569
Delta totalVI−concat = 0.0340

# Statistical claim audit

| Claim | Experimental unit | Status |
|---|---|---|
| HC ΔF1 +0.057 CI [0.020,0.091] | Target-cell paired bootstrap n=500 | OK |
| Transfer Dir A/B mean±SD | Model seed (n=20) | OK |
| Latent AUROC mean±SD | Model seed (n=20) | OK |
| Lawlor mean±SD | Fold×seed pairs (n=50); descriptive | OK; not 50 biological replicates |
| 50/50 favor totalVI | Paired fold×seed consistency | OK with caveat |
| 19/20 seeds | Seed pairs | OK |
| Residual largest | MoM components under 5×5×2 | OK; not irreducible biology |
| AURC lower for totalVI | Integrated risk–coverage summary | OK; not conflated with single-coverage risk |
| Corruption deltas | Mask seeds (5) with fixed scVI | OK |
| ECE tied Lawlor | Fold×seed summaries | OK |

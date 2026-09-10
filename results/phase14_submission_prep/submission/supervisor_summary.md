# Supervisor summary (internal)

## Research question
How reliable are deep generative single-cell multi-omics models (especially scVI_matched vs totalVI) under protein degradation, dataset/donor shift, training stochasticity, and uncertainty-guided selective prediction?

## Datasets / models
- Development: PBMC10k / PBMC5k CITE-seq
- External: Lawlor Baseline CITE-seq (author labels; 16,175 singlets; 5,207 labeled; 12 proteins; 5 donor folds × 10 seeds)
- Models: PCA, MOFA+, scVI_matched, totalVI

## Strongest findings
1. External donor-held-out macro-F1: 0.800 → 0.888 (Δ=+0.088; 50/50)
2. Predictive uncertainty detects errors (AUROC ~0.84–0.86 internal; ~0.80–0.84 external)
3. Latent posterior uncertainty more informative for totalVI; direction replicates externally (attenuated)
4. Abstention lowers risk but creates cell-type retention imbalance (internal + external)
5. Run-level stochastic residual dominates replicated variance once separated

## Strongest negative / cautionary results
- Transfer Direction A near-null (not a robust multimodal loss)
- ECE not universally better for totalVI
- Precision-weighted PoE fusion is not a novel contribution
- Missing-modality experiment not honestly feasible in totalVI 1.3.3

## Proposed journal
Primary: **Genome Biology** (Research)  
Backup: **Bioinformatics** (Original Paper; compressed)

## Optional remaining experiments (not required)
- Stimulation shift (LPS / CD3_CD28): high-value optional
- 39-ADT sensitivity: low-value optional
- External source-composition perturbation: optional if fragility claims are foregrounded
- Non-PBMC cohort: high-value optional for breadth

## Ask for supervisor
Approve (i) Genome Biology targeting and title, (ii) whether to pursue optional stimulation later, (iii) author list / funding text, (iv) human figure polish before submission.

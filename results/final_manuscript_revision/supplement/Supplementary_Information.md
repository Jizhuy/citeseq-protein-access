# Supplementary Information — Reliability of scVI and totalVI for CITE-seq

# Supplementary Methods

This Supplementary Information expands Methods details for reproducibility. All numerical results in the main text remain frozen; no new experiments were performed for manuscript revision.

## S1. Development cohort preprocessing and annotation

PBMC10k and PBMC5k CITE-seq objects were obtained from the public scvi-tools/scverse exampledata URLs recorded in config/config.yaml (pbmc_10k_protein_v3.h5ad and pbmc_5k_protein_v3.h5ad). After QC and feature harmonization, counts were: PBMC10k 6,855 cells × 16,727 genes × 14 proteins; PBMC5k 3,994 × 16,581 × 29; combined_inner 10,849 × 15,792 × 14 shared proteins.

Raw RNA UMI counts were retained for scVI/totalVI (integer count semantics). PCA/MOFA+ used library-size normalized, log1p, scaled RNA on 2,000 HVGs. Development labels were transferred via a Seurat v4 multimodal reference pathway independent of evaluated latents (SCANVI/scArches on the reference); high-confidence cells required l2 confidence ≥0.85 (9,494/10,849).

## S2. Matched scVI configuration

scVI_matched: n_latent=20, n_hidden=256, n_layers=2, dropout_rate=0.2, gene_likelihood='nb', dispersion='gene', latent_distribution='normal'; batch_size=256; max_epochs=200; train_size=0.9; early stopping on elbo_validation (patience=45). Software: scvi-tools 1.3.3. Matched settings control major capacity parameters relative to totalVI but do not imply architectural identity.

## S3. MOFA+ preprocessing and diagnostics

MOFA+ (mofapy2 0.7.5): 20 factors; RNA view = 2,000 HVGs (normalized/log1p/scaled); protein view = CLR then z-score; Gaussian likelihoods; scale_views=True, weight_views=True. Variance-explained diagnostics are retained in project results tables.

## S4. Protein corruption generation

Corruption type: protein_nonzero_mask_to_zero. Randomly selected originally non-zero measured protein entries were set to zero; preexisting zeros unchanged. Fractions: 0, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85. Mask seeds 0–4 for nonzero levels; model_seed fixed at 0. totalVI retrained per condition; scVI_matched latent fixed. Primary endpoint: l2 high-confidence logistic-regression macro-F1 difference (totalVI − scVI_matched).

## S5. Cross-dataset scArches adaptation

Directions A (PBMC10k→PBMC5k) and B (PBMC5k→PBMC10k). Source model training followed by prepare_query_anndata + load_query_data; max_epochs=200; plan weight_decay=0.0. Post-adaptation coordinate rule: encode source and target with the adapted query model only. All manuscript transfer results use corrected coordinates.

## S6. Source-size sensitivity

PBMC10k downsampled to 3,994 cells; five source-subset seeds. Purpose: test whether Direction A/B asymmetry is explained solely by source sample size. Interpreted as partial/sensitivity explanation only.

## S7. Posterior variance extraction and Monte Carlo validation

U_latent(i)=(1/d) Σ_j Var[z_ij|x_i], d=20. scvi-tools 1.3.3 get_latent_representation(..., return_dist=True) returns qz.loc and qz.scale; variances = scale². Supplementary Monte Carlo checks confirmed variance semantics against sampled latents in project logs (phase11b posterior MC check).

## S8. Twenty-seed uncertainty analysis

n=20 model seeds; paired scVI_matched/totalVI comparisons; seed-pair bootstrap 10,000 replicates; target-cell bootstrap 500 for representation ΔF1. Metrics: predictive AUROC/AUPRC, NLL, Brier, ECE (15 bins; equal-frequency primary), AURC, risk–coverage.

## S9. Ensemble analysis

Ensemble probability averaging and disagreement diagnostics were treated as supplementary. Ensemble disagreement is not claimed as exact Bayesian mutual information. Detailed ensemble ECE behavior is retained in supplementary tables where generated.

## S10. Replicated variance-component derivation

Model: Y_smr = μ + A_s + B_m + (AB)_sm + ε_smr for 5×5×2 design per model. MoM: s2_e=MS_E; s2_AB=(MS_AB−MS_E)/n; s2_A=(MS_A−MS_AB)/(b n); s2_B=(MS_B−MS_AB)/(a n). Report raw and nonnegative truncated components. Residual = run-level residual under the scheme, not irreducible biology.

## S11. Power and identifiability notes

With two run replicates, power to exclude modest interactions is limited. Interaction estimates are descriptive. Identifiability caveats are documented in variance-replication logs.

## S12. Lawlor harmonization

Publication: Lawlor et al., Front Immunol 2021;12:636720. HCA efea6426-510a-4b60-9a19-277e52bfa815; ENA PRJEB40376/ERP124005. Baseline singlets only. Genes: Ensembl GRCh37.87 mapping; duplicate symbols summed on raw UMIs; 12,776 shared genes. Proteins (12): CD3, CD4, CD8a, CD14, CD16, CD56, CD19, CD25, CD45RA, CD45RO, PD-1↔CD279_PD1, CD127↔CD127_IL7Ra; CD15 and TIGIT absent. Author labels only (obs author_cell_type). Hao/Seurat-v4 labels were not applied.

## S13. External donor-fold details

Donors Donor01–Donor10; folds = consecutive pairs; each donor target once; 8 source donors/fold; 10 seeds; 2 models; 100 runs. Eligible classes (≥20 labeled source and ≥10 labeled target every fold): B, CD14_Mono, CD4T_Mem, CD4T_Naive, CD8T_Mem, CD8T_Naive, NK. Batch key = run_identifier; donor not used as batch.

## S14. Compute and reproducibility

Recorded environment: Python 3.11.7; scvi-tools 1.3.3; scanpy 1.10.3; anndata 0.11.4; PyTorch 2.14.0 (CUDA 12.6 on GPU runs); mofapy2 0.7.5. Scripts under scripts/ and src/experiments/. Public repository URL and Zenodo DOI remain placeholders until human release approval.

# Supplementary Results

Supplementary Results mirror main Findings with expanded per-seed, per-fold, per-protein, and variance-component tables retained in the project results tree. No obsolete mixed-coordinate transfer values are included.

# Supplementary Figure Index (S1–S14 plan)

S1 Matched scVI geometry controls. S2 MOFA+ diagnostics. S3 Additional representation metrics. S4 Per-protein recovery under corruption. S5 Transfer coordinate diagnostics. S6 Multi-seed transfer distributions. S7 Ensemble/ECE behavior. S8 Posterior variance Monte Carlo check. S9 Replicated variance heatmaps. S10 Variance power/identifiability. S11 External gene/protein harmonization. S12 External per-fold metrics. S13 Full cell-type retention curves. S14 Runtime/resource notes. Publication-facing titles omit internal phase numbers. Composite panels for S1–S14 are assembled from existing project figure artifacts; this SI document provides legends and methods sufficient for supervisor review.

# Supplementary Table Index

ST1 Dataset inventory. ST2 Full representation metrics. ST3 Corruption deltas. ST4 Corrected transfer seed table. ST5 20-seed uncertainty. ST6 Variance components. ST7 Lawlor 100-run primary. ST8 Software versions. ST9 Optional prior-art matrix for multimodal VAEs.

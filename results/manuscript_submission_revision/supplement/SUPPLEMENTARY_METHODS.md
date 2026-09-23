# Supplementary Methods

Companion to: *Evaluating protein-access and model-associated gains in CITE-seq representation learning under distribution shift*

All parameters below were verified against repository code, configs, or frozen outputs unless marked **cannot verify**.

---

## S1. Dataset provenance and QC

### Internal PBMC CITE-seq

| Item | Verified value |
| --- | --- |
| Source files | `pbmc_10k_protein_v3.h5ad`, `pbmc_5k_protein_v3.h5ad` (scvi-tools exampledata; YosefLab/scVI-data fallback) |
| Loader | Project loader (`src/data/load_citeseq.py`); does **not** use `scvi.data.pbmcs_10x_cite_seq` as primary API |
| Combined object | Inner join: **10,849** cells × **15,792** genes × **14** proteins |
| Per-dataset cells | PBMC10k **6,855**; PBMC5k **3,994** |
| Additional cell/gene QC in project | **None** (`additional_cell_filter: false`); inherited from processed source objects |
| GEO/SRA | **Not claimed** (not present in local metadata) |

### Lawlor external cohort

| Item | Verified value |
| --- | --- |
| HCA | `efea6426-510a-4b60-9a19-277e52bfa815` |
| ENA | `PRJEB40376` / `ERP124005` |
| Biological unit | Donor (*n* = 10) |
| Labels | Author labels; protein-gated |
| Batch key for deep models | Sequencing lane / `run_identifier` (donor **not** used as batch) |

---

## S2. PCA / concat preprocessing

| Step | Setting | Fit scope (development) | Fit scope (inductive transfer) |
| --- | --- | --- | --- |
| HVG | 2,000; `seurat_v3` | Train HC | Source |
| RNA norm | total count 1e4 → log1p → scale clip ±10 | Train | Source |
| RNA PCA | 10 PCs | Train | Source |
| Protein | CLR (pseudocount 1) → z-score (`ddof=1`) → 10 PCs | Train | Source |
| Concat | 10+10 PCs → `StandardScaler` | Train | Source |
| Transductive concat control | Same unsupervised steps on **source + unlabeled target** | — | Joint unlabeled features; classifier source-only |

Note: `config.yaml` lists `n_top_genes: 4000` for historical settings; the **primary PCA path hard-codes 2000**.

---

## S3. Deep model hyperparameters

### scVI_matched (primary RNA deep comparator)

`n_latent=20`, `n_hidden=256`, `n_layers=2`, `dropout_rate=0.2`, `gene_likelihood="nb"`, `dispersion="gene"`, `latent_distribution="normal"`, `use_observed_lib_size=True`; `max_epochs=200`, early stopping patience **45**, `train_size=0.9`, `batch_size=256`. Learning rate: scvi-tools TrainingPlan default (not overridden).

### scVI default (not the primary reported comparator)

`n_hidden=128`, `n_layers=1`, `gene_likelihood="zinb"`, `dropout_rate=0.1` — retained only for historical PHASE-3 runs.

### totalVI

`n_latent=20`, `n_hidden=256`, encoder layers 2 / decoder 1, dropout 0.2, `gene_likelihood="nb"`, `protein_dispersion="protein"`; `lr=4e-3`, `reduce_lr_on_plateau=True`; same epoch/early-stop/batch settings.

### scArches transfer

Source train → prepare/load query → unsupervised adaptation (`query_max_epochs=200`, `weight_decay=0`) → encode **both** cohorts with adapted model → logistic regression on post-adaptation source coordinates.

---

## S4. Classifier

`sklearn.linear_model.LogisticRegression`: `max_iter=2000`, `solver="lbfgs"`, `random_state=seed`. Defaults: `C=1.0`, `class_weight=None`. Same strategy across representations.

---

## S5. Internal split

- Stratified by `batch`, `test_fraction=0.20`, seed `0` on combined 10,849 cells (≈8,679 / 2,170).  
- Development metrics: high-confidence subset (**7,573** / **1,921**; *n*=9,494).  
- Chronology: split → train-only unsupervised fits → fixed deep embeddings → train-only classifier → held-out eval.  
- Label audit: RNA-only reannotation agreement **1.0** on 9,494 HC cells.

---

## S6. Bootstrap

- Primary: cell-level **global** paired bootstrap, *B*=10,000, seed `20260923`, percentile 2.5/97.5.  
- Sensitivity: **class-stratified** paired bootstrap (within-class resampling), seed `20260924`, same *B*.  
- Conditional on fixed prediction vectors; **not** biological CIs.

---

## S7. Degradation designs

- Entry-wise mask of originally nonzero protein counts; preexisting zeros unchanged.  
- Fractions: 0.00, 0.10, 0.25, 0.40, 0.55, 0.70, 0.85.  
- Training-time: retrain totalVI; 5 mask replicates at nonzero levels; `model_seed=0`; clean *n*=1.  
- Post-adaptation: corrupt target only after clean adaptation; classifier frozen.  
- Marker-channel dropout (Supplement): zeros entire protein channels before protein/concat PCA; **totalVI not retrained**.

---

## S8. Uncertainty and selective prediction

- \(U_{\mathrm{pred}}=1-\max_k p_k\) (downstream logistic regression).  
- Latent dispersion = mean posterior variance over latent dimensions (not epistemic uncertainty).  
- Primary selective risk uses accuracy; macro-F1 / balanced accuracy vs coverage are sensitivities.

---

## S9. Stochasticity design

5 source subsets × 5 model seeds × 2 run replicates; MoM variance components for source, seed, interaction, residual run-level variation. Residual ≠ biological variance.

---

## S10. Software (pinned project environment)

| Package | Version |
| --- | --- |
| Python | 3.11.7 |
| scvi-tools | 1.3.3 |
| scikit-learn | 1.5.2 |
| scanpy | 1.10.3 |
| anndata | 0.11.4 |
| torch | 2.14.0 |

Script: `scripts/print_environment.py`.

---

## S11. Blocked / not implemented in this revision

- totalVI retraining under marker-channel dropout.  
- Additional multimodal architecture (WNN / MultiVI / Cobolt).  
- Second external cohort with non–protein-gated labels.  
- Baseline PHASE 3–5 early-stop epoch manifests (absent from checkout).

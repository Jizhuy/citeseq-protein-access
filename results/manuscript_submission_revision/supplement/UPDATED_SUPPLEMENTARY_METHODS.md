# UPDATED Supplementary Methods

Companion to: *Evaluating protein-access and residual totalVI-associated gains in CITE-seq representation learning under distribution shift*

All parameters below were verified against repository code, configs, logs, or frozen outputs unless marked **cannot verify**.

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

| Item | Verified value | Source |
| --- | --- | --- |
| HCA / ENA | `efea6426-510a-4b60-9a19-277e52bfa815`; `PRJEB40376` / `ERP124005` | project metadata |
| Singlets | **16,175** | `data_prep_summary.json` |
| Author-labeled cells | **5,207** | same |
| Harmonized genes | **12,776** | same (`n_uniquely_mapped_symbols_in_intersection`) |
| Matched proteins | **12** | same |
| Biological donors | **10** | folds / donor tables |
| Gene mapping | Ensembl → symbol; 55 duplicated symbols aggregated by summing | `gene_report` |
| Protein aliases | alias-harmonized to internal names (e.g. PD-1 ↔ CD279) | `external_protein_overlap.csv` |
| Labels | Author labels; protein-gated | Methods |
| Batch key for deep models | Sequencing lane / `run_identifier` (donor **not** used as batch) | `data_prep_summary.json` |

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

## S3. Deep model hyperparameters and information access

### scVI_matched (primary RNA deep comparator)

`n_latent=20`, `n_hidden=256`, `n_layers=2`, `dropout_rate=0.2`, `gene_likelihood="nb"`, `dispersion="gene"`, `latent_distribution="normal"`, `use_observed_lib_size=True`; `max_epochs=200`, early stopping patience **45**, `train_size=0.9`, `batch_size=256`. Learning rate: scvi-tools TrainingPlan default (not overridden).

### scVI default (not the primary reported comparator)

`n_hidden=128`, `n_layers=1`, `gene_likelihood="zinb"`, `dropout_rate=0.1` — retained only for historical PHASE-3 runs.

### totalVI

`n_latent=20`, `n_hidden=256`, encoder layers 2 / decoder 1, dropout 0.2, `gene_likelihood="nb"`, `protein_dispersion="protein"`; `lr=4e-3`, `reduce_lr_on_plateau=True`; same epoch/early-stop/batch settings.

### Internal deep-model information access

**PRIMARY (strict train-only; matched information access).** PCA, scVI_matched, and totalVI unsupervised fits used only the **7,573** high-confidence training cells. The **1,921** held-out HC cells were encoded after fitting (frozen PCA transform; scVI/totalVI `get_latent_representation` with transferred AnnData setup; no query adaptation). Classifier: logistic regression on HC train labels only; score HC held-out only. Hyperparameters unchanged (`n_latent=20`, matched capacity, `max_epochs=200`, patience 45, `batch_size=256`, totalVI `lr=4e-3`). Artifacts: `results/manuscript_submission_revision/strict_train_only/`.

**Supplementary sensitivity (historical feature-transductive / CASE B).** Earlier repository VAEs were trained on full `pbmc_cite_combined_inner.h5ad` (**10,849** cells); `obs['split']` unused for VAE subsetting; `train_size=0.9` is an internal ELBO split only. That unequal-access evaluation yielded RNA 0.666 / concat 0.745 / scVI_matched 0.722 / totalVI 0.779; \(G_{\mathrm{assoc}}\approx+0.034\) [−0.012, 0.073]. It is **not** the primary internal result.

### scArches transfer

Source train → prepare/load query → unsupervised adaptation (`query_max_epochs=200`, `weight_decay=0`) → encode **both** cohorts with adapted model → logistic regression on post-adaptation source coordinates → target scoring. Target labels withheld until scoring. Feature-access-matched concat matches unlabeled target feature access but **not** objective, architecture, optimization, or adaptation algorithm.

---

## S4. Classifier

`sklearn.linear_model.LogisticRegression`: `max_iter=2000`, `solver="lbfgs"`, `random_state=seed`. Defaults: `C=1.0`, `class_weight=None`. Same strategy across representations.

---

## S5. Internal split chronology

1. Stratify by `batch`, `test_fraction=0.20`, seed `0` on **10,849** cells → exact **8,679** train / **2,170** test.
2. Restrict to high-confidence labels → **7,573** / **1,921** (9,494 HC).
3. **PRIMARY:** unsupervised fits for PCA, scVI_matched, and totalVI on HC train only; encode held-out after fit.
4. Classifier: HC train labels only; evaluate HC held-out.
5. **Supplementary sensitivity:** historical VAEs fit on all 10,849 features (CASE B / feature-transductive).

Label audit: RNA-only reannotation agreement **1.0** on 9,494 HC cells.

---

## S6. Bootstrap

- Primary: cell-level **global** paired bootstrap, *B*=10,000, seed `20260923`, percentile 2.5/97.5.
- Sensitivity: **class-stratified** paired bootstrap (within-class resampling), seed `20260924`, same *B*.
- Conditional on fixed prediction vectors; **not** biological CIs; **not** donor-level intervals.

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

## S9. Stochasticity design and variance model

Balanced design: 5 source subsets × 5 model seeds × 2 run replicates.

\[
Y_{ijr}^{(a)}=\mu^{(a)}+S_i^{(a)}+G_j^{(a)}+(SG)_{ij}^{(a)}+\varepsilon_{ijr}^{(a)},
\]
\[
\Delta_{ijr}=Y_{ijr}^{(\mathrm{totalVI})}-Y_{ijr}^{(\mathrm{scVI\_matched})}.
\]

Source / seed / interaction / residual are **computational** components, not biological variance. MoM expected mean squares used for variance fractions.

---

## S10. Software (verified project environment logs / `pip show`)

| Package | Version |
| --- | --- |
| Python | 3.11.7 |
| scvi-tools | 1.3.3 |
| scikit-learn | 1.5.2 |
| scanpy | 1.10.3 |
| anndata | 0.11.4 |
| PyTorch | 2.14.0 |
| numpy | 1.26.4 |
| pandas | 2.2.2 |
| scipy | 1.13.1 |
| mofapy2 | 0.7.5 |
| CUDA (local Mac log) | not available (MPS available); PHASE 10 used CUDA server builds |

Pins: `environment.yml`, `requirements.txt`. Runtime snapshot: `results/logs/environment.json`.

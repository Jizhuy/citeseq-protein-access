# Methods reproducibility audit

| Required item | Present in manuscript? | Source |
|---|---|---|
| Dataset counts PBMC10k/5k/combined | Yes | phase2_data_report |
| HC cells 9494 / threshold 0.85 | Yes | phase6 / annotation |
| Corruption nonzero definition + fractions + 5 masks | Yes | config + perturbations.py |
| scVI/totalVI hyperparameters | Yes | runners |
| scvi-tools 1.3.3 | Yes | manifests |
| MOFA+ 0.7.5 | Yes | logs |
| LogisticRegression kwargs (explicit) | Yes; unspecified C/penalty/multi_class noted as sklearn defaults | code |
| ECE 15 bins EF/EW | Yes | phase11 |
| Bootstrap 500 / 10000 | Yes | phase6/11b |
| Post-adaptation coordinate rule | Yes | phase10 erratum corrected |
| Lawlor HCA/ENA | Yes | phase12b |
| Gene/protein harmonization | Yes | phase12b |
| Donor folds / class thresholds | Yes | phase12b |
| Variance model formula | Yes | phase12a_stats |
| Missing-modality limitation | Yes | methods |
| Compute environment | Yes | manifests |
| Learning rate scVI default | Noted as not forced / totalVI 4e-3 | runners |

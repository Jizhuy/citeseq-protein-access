# Migration audit

The migration created this repository (`Jizhuy/citeseq-multimodal-reliability`) from the development repository (`Jizhuy/citeseq-protein-access`). It ran only file copying, renaming, path normalisation and documentation, plus lightweight checks: compile, import, parsing, link and ZIP-integrity checks, and figure redraws from frozen tables. No experiment, model training, query adaptation, bootstrap or statistical recomputation was performed.

Per-file actions are listed in [`MIGRATION_MANIFEST.md`](MIGRATION_MANIFEST.md); every code edit is in [`migration_edit_record.json`](migration_edit_record.json).

Historical Genome Biology-addressed draft materials were excluded from the new publication repository. Their presence in the development repository does not establish whether a formal submission occurred.

## Summary counts

| Item | Count |
| --- | --- |
| Files migrated, with a manifest row | 233 (136 copied unchanged, 46 renamed only, 50 path normalised, 1 syntax fixed) |
| Frozen result files under `results/` | 129 (124 byte-identical; 5 log/JSON files with provenance-path text normalised) |
| Copied or renamed files verified byte-identical to their source (MD5) | 182 / 182 |
| New files written for this repository | 26 (24 hand-written: root README, CITATION.cff, LICENSE, .gitignore, 11 folder READMEs, 7 `docs/` files, 2 provenance helper scripts; 2 generated: `docs/MIGRATION_MANIFEST.md`, `docs/migration_edit_record.json`) |
| Total files in repository | 259 |
| Development-repository files not migrated | 53,273 (including 17,728 cache files and 59 data/checkpoint files) |

## Migration checks (at migration time, before the first commit)

| # | Check | Result | Evidence |
| --- | --- | --- | --- |
| 1 | Primary values unchanged | PASS | `results/primary/strict_train_only/tables/*` byte-identical to source (MD5); RNA PCA 0.667, concat 0.745, scVI_matched 0.709, totalVI 0.741 |
| 2 | 314 manuscript decimals not altered | PASS | GPB package audit after the repository-URL update: 314 = 314 decimals; frozen source PDF MD5 `46da53357b431778bde3f68123ed9488` unchanged |
| 3 | Strict train-only analysis is PRIMARY | PASS | README §3–4; `docs/REPRODUCIBILITY.md` §1; `results/primary/README.md` ("COMPLETE — PRIMARY ANALYSIS") |
| 4 | Feature-transductive analysis is SECONDARY | PASS | Only under `results/supplementary/information_access_sensitivity/` ("SECONDARY INFORMATION-ACCESS SENSITIVITY") |
| 5 | G_assoc is descriptive | PASS | README §1 and §2; `docs/REPRODUCIBILITY.md` §7; no causal wording |
| 6 | Direction A unchanged | PASS | `matched_transfer_control.csv` byte-identical (0.735 / 0.636 / 0.656 / 0.652) |
| 7 | Direction B unchanged | PASS | Same file (0.734 / 0.767 / 0.601 / 0.654) |
| 8 | Lawlor n = 10 donors | PASS | `lawlor_donor_summary_stats.csv` (`n_donors` = 10; 10/10; P = 0.001953125) byte-identical |
| 9 | 38/50 described as computational | PASS | README, `results/donor_validation/README.md`, `docs/REPRODUCIBILITY.md`: fold × seed pairs, not biological replication |
| 10 | Analyses not run are stated as NOT RUN | PASS | README §4; `docs/REPRODUCIBILITY.md` §1; degradation, supplementary and donor-validation READMEs |
| 11 | Figures 1–6 migrated | PASS | `figures/Figure_1.tiff` … `Figure_6.tiff` byte-identical to the final package |
| 12 | Supplementary Figure S1 migrated | PASS | `figures/Supplementary_Figure_S1.{pdf,png,svg}`; described as a design schematic with illustrative icons |
| 13 | Redrawn Figures 2–6 used | PASS | `figures/source/Figure2…6` are the final redrawn versions; `scripts/generate_publication_figures.py` reproduces them pixel-identically from frozen tables |
| 14 | Figure 3 rounding documented | PASS | `docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md` §1 (−0.00324 → −0.003 vs 0.652 − 0.656 = −0.004); not "fixed" |
| 15 | Old drafts excluded | PASS | 11 root-level draft DOCX/PDF files and all superseded manuscript packages excluded |
| 16 | Genome Biology material excluded | PASS | 54 Genome Biology-addressed draft files excluded; two historical manuscript-draft paths cited in one frozen log (`lawlor_inconsistency_audit.md`) were replaced by a neutral description |
| 17 | No claim about a Genome Biology submission | PASS | No such claim; Genome Biology is mentioned only in the exclusion record, this audit and the author checklist |
| 18 | GPB cover letter excluded | PASS | `02_GPB_Cover_Letter.docx` not copied, not rendered |
| 19 | GPB metadata checklist excluded | PASS | `03_GPB_Metadata_Checklist.docx` not copied, not rendered |
| 20 | No personal submission information | PASS | Scan of all text, PDF and ZIP content: no postal address, phone number, ORCID, submission-system ID, funding/COI draft or reviewer suggestion. Only contact: the author email already public in the development README. |
| 21 | No raw datasets | PASS | No `.h5ad` / `.h5` / `.mtx` / `.Rds` / FASTQ files; `data/raw`, `data/processed` and `data/external` are git-ignored |
| 22 | No checkpoints | PASS | No `.pt` / `.pth` / `.ckpt` / `models/` content |
| 23 | No secrets | PASS | Scan for private keys, cloud/API tokens, password assignments and `.env`: 0 hits |
| 24 | Environment versions preserved | PASS | `environment.yml` and `requirements.txt` byte-identical (Python 3.11.7, scvi-tools 1.3.3, PyTorch 2.14.0, scikit-learn 1.5.2, …) |
| 25 | Paths portable | PASS | 0 absolute home-directory paths; code resolves paths from the repository root; the optional `CUDA_LIB_DIR` replaces a hard-coded server path |
| 26 | New URL used | PASS | README, CITATION.cff and the Supplementary Data / Information files cite `https://github.com/Jizhuy/citeseq-multimodal-reliability` |
| 27 | No DOI fabricated | PASS | No DOI anywhere; README §11 states no archival DOI has been issued |
| 28 | No license invented | PASS | No license was chosen during migration; the author later selected the MIT License (see the pre-release audit) |
| 29 | No experiment run | PASS | Only compile/import/parse/link/ZIP checks and figure redraws from frozen tables, in a scratch copy |
| 30 | No statistic regenerated | PASS | All result tables copied; none recomputed or overwritten |
| 31 | Development repository untouched | PASS | No tracked-file changes and no new commits (HEAD `56dd3b0`); FINAL DOCX, its backup and the frozen PDF MD5s unchanged. Only intended edit: the code-repository URL in the GPB build scripts and rebuilt outputs inside `results/GPB_submission_FINAL_v2/` (logged in its `logs/URL_REPLACEMENT_LOG.md`). |
| 32 | Nothing committed during migration | PASS | No `.git` existed during migration; the initial commit was made later at the author's explicit request |
| 33 | Nothing pushed during migration | PASS | No push was performed during migration |
| 34 | No release | PASS | No tag or GitHub Release created |
| 35 | No Zenodo record | PASS | None created |
| 36 | No journal submission | PASS | None made |

## Pre-release audit (first public release)

This audit was run before Git initialisation and the initial commit. Blocking checks had to PASS before any commit.

| # | Check | Result | Evidence |
| --- | --- | --- | --- |
| 1 | Primary scientific values unchanged | PASS | `results/primary/strict_train_only/tables/*` byte-identical to the development source (MD5): RNA PCA 0.667, concat 0.745, scVI_matched 0.709, totalVI 0.741; G_access +0.078 [0.042, 0.131], G_assoc −0.004 [−0.060, 0.030], Gap_total +0.074 [0.029, 0.130] |
| 2 | Strict train-only remains PRIMARY | PASS | README §3–4; `results/primary/README.md`; `docs/REPRODUCIBILITY.md` §1 |
| 3 | Historical feature-transductive analysis remains SECONDARY | PASS | Only under `results/supplementary/information_access_sensitivity/` |
| 4 | G_assoc remains descriptive | PASS | README §1–2; `docs/REPRODUCIBILITY.md` §7 |
| 5 | Transfer values unchanged | PASS | `matched_transfer_control.csv` and `transfer_replicates_*` byte-identical |
| 6 | Donor n = 10 unchanged | PASS | `lawlor_donor_summary_stats.csv` byte-identical (`n_donors` = 10; 0.690 / 0.923 / 0.942 / 0.794 / 0.884) |
| 7 | Computational vs biological replication distinguished | PASS | 38/50 described as fold × seed pairs; biological unit = 10 donors |
| 8 | Blocked analyses remain NOT RUN | PASS | README §4; Supplementary Information Table S12 lists all three as NOT RUN |
| 9 | Final figures valid | PASS | Figures 1–6 and S1 present; no unexpected or stale image files; redraws pixel-identical; figure PDFs label `scVI_matched` with no stale variants |
| 10 | Supplement valid | PASS | ZIP `testzip` OK (18 entries); PDF/ZIP text: new repository URL only, no personal data, no cover letter or journal metadata, no placeholders; blocked analyses reported as NOT RUN. The PDF refers to the archive by its journal-package file name (`05_GPB_Supplementary_Data.zip`); this is noted in the README. |
| 11 | Figure 3 rounding issue documented | PASS | `docs/KNOWN_NUMERICAL_PRESENTATION_NOTES.md` §1 |
| 12 | MIT LICENSE present | PASS | Root `LICENSE` (MIT, © 2026 Jizhu Yang); README §12 states that third-party datasets are not relicensed |
| 13 | CITATION.cff valid | PASS | `cffconvert --validate`: valid against CFF schema 1.2.0; `license: MIT`; no DOI, volume, issue, pages or dates |
| 14 | Canonical contact email is yangjizhu@outlook.com | PASS | Public contact appears in README §13 and CITATION.cff (and in this audit line); no other personal email in the repository |
| 15 | Provenance gaps accurately documented | PASS | `docs/RESULT_MANIFEST.md` §9: 2 EXACT GENERATOR ADDED, 2 FROZEN DERIVED ARTIFACT |
| 16 | No invented provenance script | PASS | Helpers were added only where output was byte-identical; none for the sensitivity table or the Lawlor conversion |
| 17 | No old manuscript drafts included | PASS | No manuscript DOCX/PDF drafts present |
| 18 | Genome Biology materials excluded | PASS | None present |
| 19 | GPB submission-only materials excluded | PASS | No cover letter, metadata checklist, submission README or manuscript DOCX |
| 20 | No secrets | PASS | 0 hits for private keys, SSH material, API/cloud tokens, password assignments; no `.env`, `*.pem` or `*.key` files |
| 21 | No personal submission metadata | PASS | No phone number, postal address, ORCID, submission-system ID, funding/COI draft, recommendation letter or application file |
| 22 | No raw public datasets committed | PASS | No `.h5ad` / `.h5` / `.mtx` / `.Rds` / FASTQ; `data/raw`, `data/processed`, `data/external` git-ignored |
| 23 | No checkpoints | PASS | No `.pt` / `.pth` / `.ckpt` / `models/` content |
| 24 | Paths portable | PASS | 0 absolute home-directory paths |
| 25 | Python syntax valid | PASS | 74 files compile |
| 26 | Imports valid | PASS | 0 unresolved `src` imports; real import of all `src` modules and import-only load of all non-figure scripts succeeded |
| 27 | YAML valid | PASS | `config/*.yaml`, `environment.yml`, `CITATION.cff` parse |
| 28 | JSON valid | PASS | All JSON files parse |
| 29 | CSV files readable | PASS | 104 CSV files parse |
| 30 | No DOI fabricated | PASS | No DOI anywhere |

## Validation results

| Check | Result |
| --- | --- |
| `py_compile` of 74 Python files | 0 errors |
| Static resolution of `src.*` imports | 0 unresolved |
| Real import of 41 `src` modules; import-only load of 30 analysis and helper scripts (scratch copy) | all succeeded |
| YAML / JSON / CSV parsing (104 CSV files) | 0 errors after the `protein_gene_map.yaml` syntax fix |
| `Supplementary_Data.zip` integrity (`testzip`) | OK, 18 entries |
| Markdown relative links | 0 broken |
| `.gitignore` tested in a scratch Git repository | no shipped file ignored; all probe paths (raw data, checkpoints, `.env`, keys, caches, logs, submission folders) ignored |
| Figure redraw from frozen tables (scratch copy) | Figures 1–6 and S1 pixel-identical |
| `__pycache__` / `.DS_Store` in repository | none |

## Changes that are not pure copies

- **Code (PATH NORMALIZED):** repository-root resolution, input/output directory strings, module renames in imports, stage-labelled output file names, usage strings, the configuration key `phase10` → `query_adaptation` (values unchanged), and an optional `CUDA_LIB_DIR`. No algorithm, parameter or seed changed.
- **Frozen logs/JSON (PATH NORMALIZED, 5 files):** development-repository path references in provenance text only. These are one absolute machine path in `data_prep_summary.json`, old directory, file and module names in `bootstrap_class_stratified_meta.json`, `class_degradation_meta.json`, `variance_model_audit.md` and `lawlor_inconsistency_audit.md`, and two historical manuscript-draft paths in the last file, replaced by a neutral description. No numerical value changed.
- **`config/protein_gene_map.yaml` (SYNTAX FIXED):** the free-text `notes` item was quoted because its trailing colon made the source file invalid YAML. No mapping entry changed.
- **`config/config.yaml`:** comments neutralised, an unused development block removed, and the `phase10` key renamed; every value used by the analyses is unchanged.

## Known, documented gaps (not fixed)

- Two of the four former provenance gaps now have exact generators: `training_time_degradation_summary.csv` and `pbmc_20seed_reliability_summary.csv`, both byte-identical on regeneration. The other two remain frozen derived artifacts with no standalone generator retained: `information_access_sensitivity_comparison.csv` and the Lawlor Matrix-Market → `.h5ad` conversion. See `docs/RESULT_MANIFEST.md` §9.
- Frozen logs and script docstrings keep historical development stage labels (for example "PHASE 11B"); `docs/REPRODUCIBILITY.md` §8 translates them. File and directory names do not use them.
- The conda environment name (`multiomics_robustness`), a logger default name and a download User-Agent string keep the development project name. They were kept to preserve the environment specification and are harmless.

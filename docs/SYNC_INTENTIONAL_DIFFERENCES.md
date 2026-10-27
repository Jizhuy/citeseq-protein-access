# Intentional differences after 2026-10-27 synchronization

SOURCE: https://github.com/Jizhuy/citeseq-multimodal-reliability  
DESTINATION: https://github.com/Jizhuy/citeseq-protein-access  

Scientific content (code, frozen results, figures, supplementary materials, config,
environment files, LICENSE) is intended to be byte-identical between the two
repositories after synchronization, except for the documentation adjustments below.

## Documentation-only intentional differences

| Path | Difference |
| --- | --- |
| `README.md` | Destination describes itself as a synchronized copy; points to the manuscript-cited companion URL; documents the archive branch/tag; repository tree label uses `citeseq-protein-access`. |
| `docs/OLD_REPOSITORY_ARCHIVE_NOTICE.md` | Rewritten for the destination repository (sync + archive pointers). The SOURCE file remains a notice prepared for the historical development repository. |
| `docs/SYNC_INTENTIONAL_DIFFERENCES.md` | Present only on DESTINATION (this file). |
| `CITATION.cff` `repository-code` | Left as `https://github.com/Jizhuy/citeseq-multimodal-reliability` so it matches the submitted manuscript Code Availability URL. |

## Non-differences / preserved separately

- SOURCE repository was not modified.
- DESTINATION pre-sync history is on `archive/legacy-pre-publication-sync-2026-10-27` / `legacy-pre-sync-2026-10-27`.
- Untracked local working files on the author's machine (submission drafts, GPB packages, etc.) are not part of either public tree.

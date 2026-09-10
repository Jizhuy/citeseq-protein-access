# Professor-Level Review

## Overall assessment
The revised manuscript is scientifically coherent and now reads as a reliability study rather than a phase diary. The Lawlor external validation is appropriately prominent, claim scope is narrowed to scVI/totalVI in PBMC CITE-seq, and key statistical wording (AURC; residual variance; fold×seed pairs) has been corrected. The paper is ready for supervisor scientific review, not yet for journal submission.

## Major strengths
1. Clear reliability framing beyond leaderboard accuracy.
2. Strong external donor-held-out replication with author labels.
3. Explicit separation of predictive vs latent uncertainty and discrimination vs calibration vs selective risk.
4. Honest negative results (Direction A near-null; retention imbalance; residual stochastic dominance).
5. Corrected post-adaptation coordinate interpretation and corruption definition.

## Major concerns
1. MAJOR — Public code/data archival identifiers still placeholders (Genome Biology blocker).
2. MAJOR — Author/funding/ethics/LLM-tool disclosures require human completion.
3. MAJOR — Live Genome Biology abstract/wording rules need human confirmation on the Springer page (secondary sources disagree on 250 vs ~350 words; current abstract is 218).
4. MAJOR — Supplementary figures S1–S14 are planned/described but not fully re-exported as a single polished SI figure set in this package.
5. MINOR/MAJOR borderline — Some NeurIPS/AAAI bibliographic page ranges remain partially verified (Wu page-range sources disagree; Naeini pages marked partial).

## Minor concerns
1. Development PBMC public landing pages are scverse/GitHub URLs rather than GEO accessions (state clearly in Data Availability).
2. Figure composites may need typography reflow for single-column print.
3. Ensemble analyses remain mostly supplementary and lightly mentioned.
4. Title is now scoped to scVI/totalVI (good); ensure Abstract/Conclusions never re-universalize.

## Statistical issues
RESOLVED: AURC no longer phrased as if reduced together with point risk incorrectly.
RESOLVED: residual not called irreducible biology.
RESOLVED: 50 pairs not treated as 50 biological replicates.
MAJOR remaining: ensure every table footnote restates experimental unit (seed vs fold×seed vs corruption mask vs cell bootstrap).

## Biological interpretation issues
RESOLVED: Lawlor independence / author labels.
RESOLVED: missing-modality vs entry-wise corruption.
MINOR: PBMC-only generalization remains a limitation and must stay visible.

## Writing issues
RESOLVED: project-diary language largely removed.
MINOR: Results still dense; supervisor may prefer moving source-size paragraph deeper into Supplement.

## Reference issues
RESOLVED: expanded from 12 to 16 verified refs; MultiVI/Cobolt/scMaui/PoE/calibration/selective/Lawlor/Hao covered.
PARTIAL: Wu NeurIPS page range conflict across databases; Naeini AAAI pages partial.

## Publication readiness
Ready for supervisor review: YES
Ready for Genome Biology submission without human completion: NO

## Exact recommended revisions before submission
1. Fill author metadata, funding, contributions, competing interests, ethics confirmation, LLM tool names.
2. Create public repo + Zenodo DOI; update Data/Code Availability.
3. Human visual QC of Figures 1–6 at single/double column.
4. Confirm live Genome Biology Research guidelines.
5. Finalize SI figure file packing for upload.
6. Resolve residual bibliographic partials (Wu/Naeini).

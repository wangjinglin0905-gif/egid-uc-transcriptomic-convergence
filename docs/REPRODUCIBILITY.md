# Reproducibility

## Public adaptation

Only local path/output-directory setup changed in the analysis copies; statistical algorithms, seeds, thresholds and frozen output bytes were retained. Provenance records original and adapted hashes. Protocol, amendment and configuration remain byte-identical so the original hash checks still hold.

## Derived-output check

Run Rscript scripts/verify_frozen_results.R . from the repository root.
The 14 checks cover UC concordance, two Hedges g estimates and approximate CIs, the 198-complete-case HC3 regression, three ORA families, two Hallmark correction hierarchies and four matched-null empirical P values.

## Full expression-level workflow

Obtain the six matrices and metadata in INPUTS.md. Set EOGE_IBD_SOURCE to the input root and EOGE_IBD_ROOT to a separate working repository copy (defaults to the working directory). Run:
1. Rscript analysis/corrective_reanalysis_v1.0.R
2. Rscript analysis/supplementary_postprocessing_v1.0.R
3. Rscript analysis/matched_null_2026-08-28/run_matched_null_v1.1.R

These scripts begin from provided or previously processed gene matrices, not FASTQ/CEL files. The core script can reuse existing null CSVs: use a clean working output directory for a new full permutation run. Never mistake stored-null reuse for a full rerun.

## Figures

Run Rscript scripts/calibrate_figures_v2_4.R from the root. Outputs go to generated_figures/main. A named figure can be supplied as a selection argument.
Run Rscript analysis/matched_null_2026-08-28/render_matched_null_figure_v1.1.R for the competitive figure under generated_figures/matched_null.

## Historical environment and limits

Core: R4.6.1, limma3.68.4, msigdbr26.1.0, sandwich3.1-2, lmtest0.9-40.
Matched-null: Matrix1.7-5, digest0.6-39, jsonlite2.0.0, ggplot2 4.0.3.
Figures additionally need patchwork, ggrepel, scales and ragg. Original sessionInfo files accompany results.

MSigDB corresponds to 2026.1 and is not bundled. Users must ensure the retrieved gene-set version matches; current-library substitution can change membership. Public accessions alone do not guarantee byte-identical historical processed matrices.

The public-copy checks run on the development computer. No fresh full expression analysis, raw-read processing or independent-machine end-to-end reproduction was claimed. No EoD full-matrix or immune-composition module is included.

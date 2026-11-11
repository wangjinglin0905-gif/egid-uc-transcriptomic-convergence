# EGID–UC transcriptomic convergence: code and derived results

Version 2.4.0 accompanies the manuscript *Transcriptomic convergence between eosinophilic colitis and ulcerative colitis with an exploratory gastric comparison*.

## Scientific question and results

Do mucosal disease effects in eosinophilic colitis (EoC) align reproducibly with active colonic ulcerative colitis (UC)? A post hoc gastric (EoG) comparison examines local overlap and selected-signature specificity.

The frozen EoC discovery result was rho=0.064 and exact two-sided P=1541/3003=0.513. Four concordant genes passed the conjunction criterion, below the minimum of 10. EoG global convergence was nonsignificant (rho=0.240; two-sided P=0.100), while a selected 209-gene score separated UC from controls. UC-control separation exceeded matched random gene sets; the within-UC histology association did not. The competitive module retains SUPPLEMENT_ONLY.

This is a secondary public-data analysis, not a clinical comorbidity cohort or systematic review. The internal plan was not publicly preregistered; earlier exploratory results preceded corrective analyses. Historical 'confirmatory' labels in unchanged outputs refer to the primary branch, not prospective preregistration. The matched-null protocol preceded its competitive results, after real-signature results were known.

## Quick check

From the repository root:

~~~sh
Rscript scripts/verify_frozen_results.R .
~~~

This checks 14 derived statistics using base R and bundled outputs; it does not require expression matrices or install packages. Results go to validation/. It does not rerun expression models or all permutations.

See [input requirements](docs/INPUTS.md), [execution and environment](docs/REPRODUCIBILITY.md), [data dictionary](docs/DATA_DICTIONARY.md) and [figure source map](docs/FIGURE_SOURCE_MAP.csv).

## Contents

- analysis/corrected_results/: gene statistics, conjunction gates, permutation distributions, pathway families, scores and models.
- analysis/matched_null_2026-08-28/: frozen protocols/configuration, scripts, 40,000 endpoint records and two RDS objects containing 10,000 matched gene sets each.
- scripts/: frozen-result checks and figure generation.
- figures/: Figure 1–6 and Supplementary Figure S1–S2, PNG and TIFF.
- provenance/: original-to-public hashes and documented deviations.
- source_metadata/: public sample-study metadata; rows include excluded samples and are not the analysis denominator.
- validation/: checks of this portable public copy.
- MANIFEST_SHA256.csv: integrity manifest.

Source expression matrices, third-party gene-set collections, article full texts, manuscript drafts, private correspondence and internal review reports are not redistributed.

## Public sources

EoC and EoG: [EGIDExpress](https://egidexpress.research.cchmc.org/data/).
UC: [GSE75214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75214), [GSE16879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16879), [GSE36807](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36807), and [GSE109142](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109142).
GSE9686 was screened but excluded from inference because disease and batch were confounded. EoD full-transcriptome data are not included.

## Authors

Yin Chen, Jinglin Wang, Yashi Peng, Fei Li, Zaili Yang, Qian Yang.
Department of Digestive Disease, Guizhou Provincial People's Hospital, Guizhou, China.
Qian Yang ORCID: https://orcid.org/0000-0002-9591-7000

CITATION.cff contains the creator list and version. Use the archived version DOI when available. Repository release does not imply peer review or journal acceptance.

## Funding

Guizhou Province High-level Innovative Talent Selection and Training Program (Contract No. gzwjrs2024-005); Science and Technology Fund Project of Guizhou Provincial Science and Technology Program (No. QKHJC-ZK [2022]-264).

## Licenses

Original project software is under MIT (LICENSE). Original derived results, figures and documentation are under CC BY 4.0 (LICENSE-DATA.md). Third-party source metadata, identifiers and annotations retain their source attribution and applicable terms; these licenses do not relicense third-party inputs. See THIRD_PARTY_NOTICES.md.

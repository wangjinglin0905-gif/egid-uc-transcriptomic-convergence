# Frozen protocol: matched-random-gene-set competitive test for the 209-gene EoG-UC signature

Protocol version: 1.0  
Freeze time: 2026-08-28T10:00:03+08:00  
Status: frozen before calculation of any matched-null result  
Scientific role: prespecified exploratory strengthening analysis; it does not replace the prespecified EoC primary endpoint.

## Question and competing explanations

The target question is whether the data-derived 209-gene EoG-UC signature transports to independent UC datasets better than comparably sized gene sets with similar technical properties. The main alternative explanation is that any sufficiently large, well-expressed, variable, or strongly UC-responsive gene set would perform similarly.

The analysis has two nested null models:

1. **Technical-property-matched null (primary competitive test):** controls signature size, up/down composition, pooled expression level, and expression variability.
2. **UC-effect-matched null (stringent sensitivity test):** additionally controls the magnitude of the UC discovery effect. This is the decisive test against the explanation that the 209 genes are merely a generic UC inflammatory/disease-activity program.

No result-dependent modification of matching variables, tolerances, endpoints, tails, covariates, or success rules is permitted. Any necessary change requires a versioned protocol amendment and explicit deviation record before rerunning.

## Authoritative inputs

- Signature definition: `analysis/corrected_results/EoG_UC_exploratory_signature_selected.csv`.
- UC discovery expression: `outputs/results_v3.1/formal_DEG_2026-08-27/expr_GSE75214_genes.csv` with the frozen sample manifest.
- UC discovery statistics: `analysis/corrected_results/GSE75214_UC_DEG_corrected.csv`.
- Cross-center disease-status evaluation: GSE36807 gene matrix and frozen sample manifest; only `sensitivity_colon_UC` and `sensitivity_colon_control` samples are eligible.
- Phenotype evaluation: GSE109142 gene matrix and frozen sample manifest; only `phenotype_validation_UC` samples are eligible.

GSE36807 and GSE109142 were not used to derive the 209-gene signature. Their disease labels and phenotype outcomes are not used to construct the matched null sets.

## Analysis universe and signature coverage

The common universe is the intersection of genes present in the GSE75214 discovery matrix, the eligible GSE36807 matrix, the eligible GSE109142 matrix, and the corrected GSE75214 differential-expression table. Genes must have finite discovery pooled mean, finite positive discovery pooled standard deviation, finite UC moderated t statistic, and non-zero UC log-fold-change direction.

All 209 signature genes are intersected with this common universe. The same covered signature and the same matched random sets are evaluated in both validation datasets. Coverage below 80% of 209 genes is a no-go condition for inferential interpretation; results would be descriptive only.

All signature genes are excluded from the random-gene candidate pool. Genes are sampled without replacement within each random set but may recur across random sets.

## Frozen matching algorithm

The signature's exact covered size and exact up/down direction counts are preserved. Direction is defined by the corrected GSE75214 UC-versus-control log-fold change, which agrees with the selected EoG direction by construction.

Matching features are calculated from the GSE75214 discovery data only and transformed to empirical percentile ranks in the common universe:

- pooled expression mean across eligible GSE75214 samples;
- log2-transformed pooled expression standard deviation;
- for the stringent null only, absolute corrected GSE75214 moderated t statistic.

For every real-signature gene slot, eligible candidates must have the same UC direction. Euclidean distance is calculated across the relevant rank-transformed features. Candidate lists are ordered by distance. For each replicate, signature slots are visited in a random order and one unused gene is sampled from the 100 nearest available candidates with weights proportional to `exp(-distance / median_distance)`. If no unused candidate remains among the nearest 100, the pool expands to the nearest 500 and then to the complete same-direction candidate pool. The scale denominator is bounded below by 1e-6. Failure to construct an exact-size, unique set is recorded; more than 1% failed replicates is a no-go condition.

Ten thousand random sets are generated per null model with seed 20260828. The full draw objects are retained as compressed RDS files.

## Match-quality gate

For each matching feature, the mean difference between every null set and the observed signature is divided by the feature standard deviation in the common universe. The match is acceptable only if, for every matched feature, the median absolute standardized difference is no greater than 0.10 and the 95th percentile is no greater than 0.25. A failed balance gate makes the corresponding competitive result descriptive and requires a prospective protocol amendment before any alternative matching is attempted.

## Frozen validation estimands

### GSE36807: cross-center disease-status effect

Each gene is z-standardized across the eligible samples. The directional signature score is the mean z score of up genes minus the mean z score of down genes. The effect estimand is Hedges g for UC versus control. The observed and null effects use the identical scoring and group definitions.

### GSE109142: histologic-severity association

Expression is transformed as `log2(x + 1)`, then each gene is z-standardized across eligible UC samples. The directional score is calculated as above. The frozen complete-case model is:

`signature score ~ eosinophil grade >=32/hpf + histologic severity + age + sex`

The competitive estimand is the HC3-robust t statistic for histologic severity, which is invariant to arbitrary score scaling. The raw coefficient, robust standard error, 95% confidence interval, two-sided model P value, and complete-case denominator are also retained. Eosinophil grade remains a specificity covariate and is not a competitive success endpoint.

## Competitive P values, multiplicity, and success rules

Positive transport was prespecified; therefore the competitive empirical P value is one-sided:

`P = (1 + number of null statistics >= observed statistic) / (10,000 + 1)`.

Two-sided empirical P values based on absolute statistics are retained as descriptive sensitivity measures. Within each null model, the two competitive P values (GSE36807 Hedges g and GSE109142 histology HC3 t) are Holm-adjusted. Both observed effects must be positive and both Holm-adjusted P values must be below 0.05 for that null model to pass.

Interpretation is frozen as follows:

- **Both null models pass:** the 209-gene signature may receive a substantive but explicitly exploratory Results subsection; evidence supports transportability beyond random technical-property-matched and equally UC-responsive gene sets.
- **Only the technical-property-matched null passes:** the signature is non-arbitrary but has not been distinguished from a generic UC-responsive program; retain only a short main-text statement or supplementary analysis.
- **Neither null model passes:** do not expand the positive narrative; report the competitive null result in the supplement and retain the EoG branch as exploratory.

Even if both models pass, the analysis does not establish causality, diagnostic performance, EoG specificity, or independence from all cell-composition effects. The EoG branch remains exploratory because the global EoG-UC convergence permutation test was not statistically significant and the 209-gene set was selected from the discovery data.

## Reproducibility outputs

- frozen configuration and protocol hash;
- common-universe and signature-coverage tables;
- match-quality diagnostics;
- full observed and null endpoint statistics;
- Holm-adjusted competitive results and frozen decision;
- compressed matched-set draws;
- exact plotted data;
- PNG, TIFF, and PDF competitive-null figure;
- R session information, console log, and SHA-256 manifest.

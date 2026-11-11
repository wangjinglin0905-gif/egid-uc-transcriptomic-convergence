# Data dictionary

CSV/TSV files have headers; missing values are blank or NA according to the frozen output.
Gene identifiers are supplied symbols/identifiers, without a new nomenclature remapping.

- DEG: gene=feature; logFC=disease-minus-comparator on modeled log scale; AveExpr=average expression; t=moderated t; P.Value=raw P; adj.P.Val=BH FDR.
- Correlation/null tables: rho or spearman_t=unitless Spearman correlation; null rows represent sample-label assignments, not independent genes.
- Conjunction tables: maximum cohort P for concordant directions, P=1 for opposite directions, then full-universe BH. The EoG selected table contains 209 genes.
- Sample scores: mean within-cohort gene z score for up genes minus mean for down genes. No clinical diagnostic calibration is implied.
- GSE109142: eos32=source eosinophil grade >=32/hpf (Y/N); histology=source ordinal severity; age=years; sex=source category.206 scored patients; 198 complete cases in the adjusted model.
- Models: estimate, robust_se, p_value, ci_low, ci_high describe HC3 regression; CIs use the recorded normal approximation.
- ORA: overlap=selected gene count; set_size=gene-set size in the common background; p_value=upper-tail hypergeometric P; q_value=BH including zero-overlap sets.
- Hallmark: cohort P/FDR and directions, candidate-level Holm validation and replicated decision;50 discovery sets.
- Matched-null: technical/uc_effect identify two models; endpoint identifies Hedges g or histology HC3 t. 10,000 gene sets per model x 2 endpoints yields 40,000 records.
- RDS files retain full matched draws; use readRDS in R.
- SUPPLEMENT_ONLY in frozen_decision.csv is the prespecified reporting decision for this post hoc competitive module, not an execution failure.

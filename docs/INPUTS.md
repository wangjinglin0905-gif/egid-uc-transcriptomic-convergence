# External inputs

Set EOGE_IBD_SOURCE to a local directory containing these relative paths. No expression matrices are bundled.

|Dataset|Relative path|
|---|---|
|EGIDExpress EoC|data_raw/EGIDExpress_EoC/EoC_RNAseq.csv|
|EGIDExpress EoG|data_raw/EGIDExpress_EoG/EGtranscriptome.csv|
|GSE75214|results_v3.1/formal_DEG_2026-08-27/expr_GSE75214_genes.csv|
|GSE16879|results_v3.1/formal_DEG_2026-08-27/expr_GSE16879_genes.csv|
|GSE36807|results_v3.1/sensitivity_2026-08-27/expr_GSE36807_genes.csv|
|GSE109142|results_v3.1/GSE109142_phenotype_2026-08-27/GSE109142_expr_genes.csv|
|Metadata|manifests/sample_manifest_v1.0.tsv|

Place the supplied source_metadata manifest under the configured source root when rerunning. It includes screened/excluded samples; its row count is not the analysis sample size. The scripts filter by role.

Historical mapping used GPL6244 gene_assignment (GSE75214) and GPL570 Gene Symbol (GSE16879/GSE36807), first gene assignment, exclusion of unmapped probes and averaging of duplicate probes. GSE109142 source files were joined by gene identifier, retaining the first duplicate and filling missing entries with zero. Transformation/filtering is in the code and manuscript. EGIDExpress inputs are not certified raw counts.

Public-study covariates originate from the source providers. No new participant information was collected. Public sources and applicable third-party terms are described in README and THIRD_PARTY_NOTICES.

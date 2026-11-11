# Prospective protocol amendment v1.1: scoreability variance gate

Amendment time: 2026-08-28T10:05:24+08:00  
Status: frozen before calculation of any matched-null endpoint statistic  
Supersedes: the common-universe eligibility paragraph in protocol v1.0 only

## Trigger

An outcome-blind scoreability preflight found 17,405 genes common to the three expression matrices. No genes had zero/non-finite pooled variance in GSE75214 or GSE36807, whereas 29 genes had zero/non-finite variance among the eligible GSE109142 UC samples. Disease-status effects, histology associations, matched-null endpoint statistics, and competitive P values were not calculated or inspected during this preflight.

## Frozen amendment

To ensure that every real and random gene contributes a defined row z score in both independent evaluation datasets, the common universe additionally requires finite positive pooled standard deviation in the eligible GSE36807 samples and finite positive pooled standard deviation in the eligible GSE109142 samples after the frozen `log2(x + 1)` transformation.

All other protocol v1.0 rules remain unchanged, including the matching features (which are still derived only from GSE75214 discovery data), random-set count, seed, endpoints, covariates, tails, multiplicity adjustment, balance gates, and decision thresholds.

This amendment is a prospective technical eligibility rule and does not use validation labels or phenotypic outcomes to select genes.

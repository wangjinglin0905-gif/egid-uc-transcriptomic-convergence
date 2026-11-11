# Supplementary post-processing for the corrected EoGE-IBD analysis
# Version 1.0; 2026-08-27

suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
})

src <- Sys.getenv("EOGE_IBD_SOURCE", unset = "")
if (!nzchar(src)) stop("Set EOGE_IBD_SOURCE to the input directory described in docs/INPUTS.md")
root <- normalizePath(Sys.getenv("EOGE_IBD_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out <- file.path(root, "analysis", "corrected_results")
cache_dir <- file.path(root, "analysis", "R_cache")
Sys.setenv(R_USER_CACHE_DIR = cache_dir)
options(stringsAsFactors = FALSE)

fit_limma <- function(mat, grp, case, ref, trend = FALSE) {
  grp <- factor(grp, levels = c(ref, case))
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  fit <- lmFit(mat, design)
  cont <- makeContrasts(contrasts = paste0(case, "-", ref), levels = design)
  fit <- contrasts.fit(fit, cont)
  fit <- eBayes(fit, trend = trend, robust = TRUE)
  tt <- topTable(fit, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt
}

rank_cor <- function(a, b) {
  g <- intersect(a$gene, b$gene)
  x <- a$t[match(g, a$gene)]
  y <- b$t[match(g, b$gene)]
  ok <- is.finite(x) & is.finite(y)
  unname(cor(x[ok], y[ok], method = "spearman"))
}

camera_table <- function(mat, grp, case, ref, pathways, trend = FALSE) {
  grp <- factor(grp, levels = c(ref, case))
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  cont <- makeContrasts(contrasts = paste0(case, "-", ref), levels = design)
  idx <- ids2indices(pathways, rownames(mat), remove.empty = TRUE)
  z <- camera(mat, idx, design, contrast = cont, inter.gene.cor = 0.01)
  z$pathway <- rownames(z)
  rownames(z) <- NULL
  z
}

cat("Computing EoC transformation sensitivity...\n")
eoc_raw <- read.csv(file.path(src, "data_raw", "EGIDExpress_EoC", "EoC_RNAseq.csv"), row.names = 1, check.names = FALSE)
grp_eoc <- factor(c(rep("Non_EoC", 8), rep("Active_EoC", 6)), levels = c("Non_EoC", "Active_EoC"))
keep_eoc <- rowSums(eoc_raw >= 1) >= 6
tt_eoc_1 <- fit_limma(log2(eoc_raw[keep_eoc, , drop = FALSE] + 1), grp_eoc, "Active_EoC", "Non_EoC", trend = TRUE)
tt_eoc_05 <- fit_limma(log2(eoc_raw[keep_eoc, , drop = FALSE] + 0.5), grp_eoc, "Active_EoC", "Non_EoC", trend = TRUE)
tt_u1 <- read.csv(file.path(out, "GSE75214_UC_DEG_corrected.csv"), check.names = FALSE)
tt_u2 <- read.csv(file.path(out, "GSE16879_UC_DEG_corrected.csv"), check.names = FALSE)

sig1 <- tt_eoc_1$gene[tt_eoc_1$adj.P.Val < 0.05]
sig05 <- tt_eoc_05$gene[tt_eoc_05$adj.P.Val < 0.05]
overlap <- intersect(sig1, sig05)
dir_ok <- sign(tt_eoc_1$logFC[match(overlap, tt_eoc_1$gene)]) ==
  sign(tt_eoc_05$logFC[match(overlap, tt_eoc_05$gene)])

transform_summary <- data.frame(
  metric = c(
    "significant_genes_log2_x_plus_1",
    "significant_genes_log2_x_plus_0.5",
    "jaccard_significant_gene_sets",
    "direction_concordance_common_significant_genes",
    "rho_with_GSE75214_log2_x_plus_1",
    "rho_with_GSE75214_log2_x_plus_0.5",
    "rho_with_GSE16879_log2_x_plus_1",
    "rho_with_GSE16879_log2_x_plus_0.5"
  ),
  value = c(
    length(sig1),
    length(sig05),
    length(overlap) / length(union(sig1, sig05)),
    mean(dir_ok),
    rank_cor(tt_eoc_1, tt_u1),
    rank_cor(tt_eoc_05, tt_u1),
    rank_cor(tt_eoc_1, tt_u2),
    rank_cor(tt_eoc_05, tt_u2)
  )
)
write.csv(transform_summary, file.path(out, "EoC_transformation_sensitivity.csv"), row.names = FALSE)

cat("Computing exploratory EoG Hallmark CAMERA results...\n")
eog_raw <- read.csv(file.path(src, "data_raw", "EGIDExpress_EoG", "EGtranscriptome.csv"), row.names = 1, check.names = FALSE)
grp_eog <- factor(c(rep("NonEG", 12), rep("EG", 9)), levels = c("NonEG", "EG"))
keep_eog <- rowSums(eog_raw >= 1) >= 9
eog_mat <- log2(eog_raw[keep_eog, , drop = FALSE] + 1)
hm <- msigdbr(species = "Homo sapiens", collection = "H")
hm_list <- split(hm$gene_symbol, hm$gs_name)
cam_eog <- camera_table(eog_mat, grp_eog, "EG", "NonEG", hm_list, trend = TRUE)
write.csv(cam_eog, file.path(out, "Hallmark_camera_EoG_exploratory.csv"), row.names = FALSE)

cam_u1 <- read.csv(file.path(out, "Hallmark_camera_GSE75214_UC.csv"), check.names = FALSE)
cam_u2 <- read.csv(file.path(out, "Hallmark_camera_GSE16879_UC.csv"), check.names = FALSE)
eog_uc <- merge(cam_eog, cam_u1, by = "pathway", suffixes = c("_EoG", "_UC_discovery"))
eog_uc <- merge(eog_uc, cam_u2, by = "pathway")
names(eog_uc)[names(eog_uc) %in% c("NGenes", "Direction", "PValue", "FDR")] <-
  paste0(names(eog_uc)[names(eog_uc) %in% c("NGenes", "Direction", "PValue", "FDR")], "_UC_validation")
eog_uc$discovery_shared <- eog_uc$FDR_EoG < 0.05 & eog_uc$FDR_UC_discovery < 0.05 &
  eog_uc$Direction_EoG == eog_uc$Direction_UC_discovery
eog_uc$validation_holm <- NA_real_
if (any(eog_uc$discovery_shared)) {
  eog_uc$validation_holm[eog_uc$discovery_shared] <-
    p.adjust(eog_uc$PValue_UC_validation[eog_uc$discovery_shared], method = "holm")
}
eog_uc$replicated <- eog_uc$discovery_shared &
  eog_uc$Direction_EoG == eog_uc$Direction_UC_validation &
  eog_uc$validation_holm < 0.05
write.csv(eog_uc, file.path(out, "Hallmark_exploratory_EoG_UC.csv"), row.names = FALSE)

cat("Supplementary post-processing complete.\n")

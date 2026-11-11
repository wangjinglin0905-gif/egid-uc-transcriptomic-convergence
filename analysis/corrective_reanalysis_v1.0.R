# Corrective reanalysis for the EoGE-IBD manuscript
# Version 1.0; 2026-08-27
# Implements the frozen v3.1 rules where the legacy scripts deviated.

suppressPackageStartupMessages({
  library(limma)
  library(msigdbr)
  library(sandwich)
  library(lmtest)
})

src <- Sys.getenv("EOGE_IBD_SOURCE", unset = "")
if (!nzchar(src)) stop("Set EOGE_IBD_SOURCE to the input directory described in docs/INPUTS.md")
root <- normalizePath(Sys.getenv("EOGE_IBD_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out <- file.path(root, "analysis", "corrected_results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(root, "analysis", "R_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = cache_dir)

options(stringsAsFactors = FALSE)
set.seed(20260824)

write_session <- function() {
  zz <- file(file.path(out, "sessionInfo.txt"), open = "wt")
  sink(zz)
  print(sessionInfo())
  sink()
  close(zz)
}

fit_limma <- function(mat, grp, case, ref, trend = FALSE) {
  grp <- factor(grp, levels = c(ref, case))
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  fit <- lmFit(mat, design)
  cont <- makeContrasts(contrasts = paste0(case, "-", ref), levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit3 <- eBayes(fit2, trend = trend, robust = TRUE)
  tt <- topTable(fit3, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  list(tt = tt, design = design, contrast = cont, fit = fit3)
}

rank_cor <- function(a, b) {
  g <- intersect(a$gene, b$gene)
  x <- a$t[match(g, a$gene)]
  y <- b$t[match(g, b$gene)]
  ok <- is.finite(x) & is.finite(y)
  list(gene = g[ok], x = x[ok], y = y[ok], rho = unname(cor(x[ok], y[ok], method = "spearman")))
}

same_direction_signature <- function(a, b) {
  g <- intersect(a$gene, b$gene)
  aa <- a[match(g, a$gene), ]
  bb <- b[match(g, b$gene), ]
  same <- sign(aa$logFC) == sign(bb$logFC) & aa$logFC != 0 & bb$logFC != 0
  pconj <- ifelse(same, pmax(aa$P.Value, bb$P.Value), 1)
  qconj <- p.adjust(pconj, method = "BH")
  data.frame(
    gene = g,
    logFC_a = aa$logFC,
    t_a = aa$t,
    p_a = aa$P.Value,
    logFC_b = bb$logFC,
    t_b = bb$t,
    p_b = bb$P.Value,
    same_direction = same,
    conjunction_p = pconj,
    conjunction_q = qconj,
    selected = same & qconj < 0.05,
    stringsAsFactors = FALSE
  )
}

zscore_rows <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  ok <- is.finite(s) & s > 0
  z <- sweep(mat[ok, , drop = FALSE], 1, mu[ok], "-")
  z <- sweep(z, 1, s[ok], "/")
  z
}

directional_score <- function(mat, up, down) {
  genes <- intersect(rownames(mat), union(up, down))
  z <- zscore_rows(mat[genes, , drop = FALSE])
  up2 <- intersect(rownames(z), up)
  dn2 <- intersect(rownames(z), down)
  if (length(up2) == 0 || length(dn2) == 0) stop("Directional score requires both up and down genes")
  colMeans(z[up2, , drop = FALSE]) - colMeans(z[dn2, , drop = FALSE])
}

hedges_g <- function(score, grp, case = "UC", ref = "CTRL") {
  x1 <- score[grp == case]
  x0 <- score[grp == ref]
  n1 <- length(x1); n0 <- length(x0)
  sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / (n1 + n0 - 2))
  d <- (mean(x1) - mean(x0)) / sp
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- J * d
  vg <- (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2))
  data.frame(
    n_case = n1,
    n_ref = n0,
    mean_case = mean(x1),
    mean_ref = mean(x0),
    hedges_g = g,
    se_g = sqrt(vg),
    ci_low = g - 1.96 * sqrt(vg),
    ci_high = g + 1.96 * sqrt(vg)
  )
}

score_label_perm <- function(score, grp, case = "UC", ref = "CTRL", nperm = 10000, seed = 20260824) {
  set.seed(seed)
  obs <- mean(score[grp == case]) - mean(score[grp == ref])
  null <- numeric(nperm)
  ncase <- sum(grp == case)
  for (i in seq_len(nperm)) {
    idx <- sample(seq_along(score), ncase)
    null[i] <- mean(score[idx]) - mean(score[-idx])
  }
  p <- (sum(abs(null) >= abs(obs)) + 1) / (nperm + 1)
  list(observed = obs, null = null, p = p)
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

run_ora <- function(sig_genes, background, collection, subcollection = NULL, prefix) {
  if (is.null(subcollection)) {
    m <- msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    m <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
  }
  sets <- split(m$gene_symbol, m$gs_name)
  sets <- lapply(sets, function(x) unique(intersect(x, background)))
  sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
  nbg <- length(background)
  nsig <- length(sig_genes)
  res <- do.call(rbind, lapply(names(sets), function(nm) {
    s <- sets[[nm]]
    ov <- intersect(sig_genes, s)
    k <- length(ov)
    p <- phyper(k - 1, length(s), nbg - length(s), nsig, lower.tail = FALSE)
    data.frame(
      term = nm,
      set_size = length(s),
      overlap = k,
      p_value = p,
      genes = paste(ov, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  res$q_value <- p.adjust(res$p_value, method = "BH")
  res <- res[order(res$q_value, res$p_value), ]
  write.csv(res, file.path(out, paste0("ORA_", prefix, "_full.csv")), row.names = FALSE)
  res
}

cat("Loading metadata and expression matrices...\n")
man <- read.delim(file.path(src, "manifests", "sample_manifest_v1.0.tsv"), check.names = FALSE)
names(man)[1] <- "dataset"  # UTF-8 BOM in the frozen TSV header

# EoC
eoc_raw <- read.csv(file.path(src, "data_raw", "EGIDExpress_EoC", "EoC_RNAseq.csv"), row.names = 1, check.names = FALSE)
stopifnot(ncol(eoc_raw) == 14)
grp_eoc <- factor(c(rep("Non_EoC", 8), rep("Active_EoC", 6)), levels = c("Non_EoC", "Active_EoC"))
keep_eoc <- rowSums(eoc_raw >= 1) >= 6
eoc_mat <- log2(eoc_raw[keep_eoc, , drop = FALSE] + 1)
fit_eoc <- fit_limma(eoc_mat, grp_eoc, "Active_EoC", "Non_EoC", trend = TRUE)
tt_eoc <- fit_eoc$tt
write.csv(tt_eoc, file.path(out, "EoC_DEG_corrected.csv"), row.names = FALSE)

# EoG exploratory
eog_raw <- read.csv(file.path(src, "data_raw", "EGIDExpress_EoG", "EGtranscriptome.csv"), row.names = 1, check.names = FALSE)
stopifnot(ncol(eog_raw) == 21)
grp_eog <- factor(c(rep("NonEG", 12), rep("EG", 9)), levels = c("NonEG", "EG"))
keep_eog <- rowSums(eog_raw >= 1) >= 9
eog_mat <- log2(eog_raw[keep_eog, , drop = FALSE] + 1)
fit_eog <- fit_limma(eog_mat, grp_eog, "EG", "NonEG", trend = TRUE)
tt_eog <- fit_eog$tt
write.csv(tt_eog, file.path(out, "EoG_DEG_corrected_exploratory.csv"), row.names = FALSE)

# UC discovery
u1_all <- read.csv(file.path(src, "results_v3.1", "formal_DEG_2026-08-27", "expr_GSE75214_genes.csv"), row.names = 1, check.names = FALSE)
m1 <- man[man$dataset == "GSE75214", ]
role1 <- m1$primary_role[match(colnames(u1_all), m1$sample_id)]
use1 <- role1 %in% c("discovery_active_colon_UC", "discovery_colon_control")
u1_mat <- u1_all[, use1, drop = FALSE]
grp_u1 <- factor(ifelse(role1[use1] == "discovery_active_colon_UC", "UC", "CTRL"), levels = c("CTRL", "UC"))
stopifnot(sum(grp_u1 == "UC") == 74, sum(grp_u1 == "CTRL") == 11)
fit_u1 <- fit_limma(u1_mat, grp_u1, "UC", "CTRL", trend = FALSE)
tt_u1 <- fit_u1$tt
write.csv(tt_u1, file.path(out, "GSE75214_UC_DEG_corrected.csv"), row.names = FALSE)

# UC within-center validation
u2_all <- read.csv(file.path(src, "results_v3.1", "formal_DEG_2026-08-27", "expr_GSE16879_genes.csv"), row.names = 1, check.names = FALSE)
u2_all <- log2(u2_all + 1)
m2 <- man[man$dataset == "GSE16879", ]
role2 <- m2$primary_role[match(colnames(u2_all), m2$sample_id)]
use2 <- role2 %in% c("validation_active_colon_UC", "validation_colon_control")
u2_mat <- u2_all[, use2, drop = FALSE]
grp_u2 <- factor(ifelse(role2[use2] == "validation_active_colon_UC", "UC", "CTRL"), levels = c("CTRL", "UC"))
stopifnot(sum(grp_u2 == "UC") == 24, sum(grp_u2 == "CTRL") == 6)
fit_u2 <- fit_limma(u2_mat, grp_u2, "UC", "CTRL", trend = FALSE)
tt_u2 <- fit_u2$tt
write.csv(tt_u2, file.path(out, "GSE16879_UC_DEG_corrected.csv"), row.names = FALSE)

# Descriptive UC reproducibility
uc_common <- intersect(tt_u1$gene, tt_u2$gene)
uc1 <- tt_u1[match(uc_common, tt_u1$gene), ]
uc2 <- tt_u2[match(uc_common, tt_u2$gene), ]
disc_sig <- uc1$adj.P.Val < 0.05
uc_repro <- data.frame(
  common_genes = length(uc_common),
  spearman_t = cor(uc1$t, uc2$t, method = "spearman"),
  discovery_significant_present = sum(disc_sig),
  direction_concordant_n = sum(sign(uc1$logFC[disc_sig]) == sign(uc2$logFC[disc_sig])),
  direction_concordant_prop = mean(sign(uc1$logFC[disc_sig]) == sign(uc2$logFC[disc_sig])),
  both_significant_same_direction = sum(disc_sig & uc2$adj.P.Val < 0.05 & sign(uc1$logFC) == sign(uc2$logFC))
)
write.csv(uc_repro, file.path(out, "UC_within_center_reproducibility.csv"), row.names = FALSE)
write.csv(data.frame(
  gene = uc_common,
  t_discovery = uc1$t,
  t_validation = uc2$t,
  logFC_discovery = uc1$logFC,
  logFC_validation = uc2$logFC,
  q_discovery = uc1$adj.P.Val,
  q_validation = uc2$adj.P.Val
), file.path(out, "UC_reproducibility_gene_level.csv"), row.names = FALSE)

# Primary EoC discovery: moderated-t ranks and exact, two-sided label enumeration
cat("Primary EoC-UC discovery permutation (3003 exact assignments)...\n")
rc_eoc_u1 <- rank_cor(tt_eoc, tt_u1)
common_eoc_u1 <- rc_eoc_u1$gene
u1_t_fixed <- rc_eoc_u1$y
disc_null_file <- file.path(out, "EoC_UC_discovery_exact_null.csv")
if (file.exists(disc_null_file)) {
  rho_eoc_u1_null <- read.csv(disc_null_file)$rho
  cat("  Reusing completed discovery null from checkpoint.\n")
} else {
  combs <- combn(ncol(eoc_mat), 6)
  rho_eoc_u1_null <- numeric(ncol(combs))
  for (i in seq_len(ncol(combs))) {
    gp <- rep("Non_EoC", ncol(eoc_mat))
    gp[combs[, i]] <- "Active_EoC"
    fp <- fit_limma(eoc_mat, factor(gp, levels = c("Non_EoC", "Active_EoC")), "Active_EoC", "Non_EoC", trend = TRUE)$tt
    xt <- fp$t[match(common_eoc_u1, fp$gene)]
    rho_eoc_u1_null[i] <- cor(xt, u1_t_fixed, method = "spearman")
  }
  write.csv(data.frame(rho = rho_eoc_u1_null), disc_null_file, row.names = FALSE)
}
p_eoc_u1 <- sum(abs(rho_eoc_u1_null) >= abs(rc_eoc_u1$rho)) / length(rho_eoc_u1_null)

# Primary EoC validation: freeze the observed EoC rank and permute UC labels
cat("Primary EoC-UC validation permutation (10000 Monte Carlo assignments)...\n")
rc_eoc_u2 <- rank_cor(tt_eoc, tt_u2)
common_eoc_u2 <- rc_eoc_u2$gene
eoc_t_fixed_u2 <- rc_eoc_u2$x
nperm <- 10000
val_null_file <- file.path(out, "EoC_UC_validation_perm_null.csv")
if (file.exists(val_null_file)) {
  rho_eoc_u2_null <- read.csv(val_null_file)$rho
  cat("  Reusing completed validation null from checkpoint.\n")
} else {
  set.seed(20260824)
  rho_eoc_u2_null <- numeric(nperm)
  for (i in seq_len(nperm)) {
    idx <- sample(seq_len(ncol(u2_mat)), sum(grp_u2 == "UC"))
    gp <- rep("CTRL", ncol(u2_mat)); gp[idx] <- "UC"
    fp <- fit_limma(u2_mat, factor(gp, levels = c("CTRL", "UC")), "UC", "CTRL", trend = FALSE)$tt
    yt <- fp$t[match(common_eoc_u2, fp$gene)]
    rho_eoc_u2_null[i] <- cor(eoc_t_fixed_u2, yt, method = "spearman")
  }
  write.csv(data.frame(rho = rho_eoc_u2_null), val_null_file, row.names = FALSE)
}
p_eoc_u2 <- (sum(abs(rho_eoc_u2_null) >= abs(rc_eoc_u2$rho)) + 1) / (nperm + 1)

# Primary signature gate
sig_eoc_tbl <- same_direction_signature(tt_eoc, tt_u1)
write.csv(sig_eoc_tbl, file.path(out, "EoC_UC_signature_gate_full.csv"), row.names = FALSE)
sig_eoc <- sig_eoc_tbl$gene[sig_eoc_tbl$selected]

# EoC transformation and leave-one-out rho sensitivity
fit_eoc_05 <- fit_limma(log2(eoc_raw[keep_eoc, , drop = FALSE] + 0.5), grp_eoc, "Active_EoC", "Non_EoC", trend = TRUE)
rho_transform_u1 <- rank_cor(fit_eoc_05$tt, tt_u1)$rho
rho_transform_u2 <- rank_cor(fit_eoc_05$tt, tt_u2)$rho
loo <- do.call(rbind, lapply(seq_len(ncol(eoc_mat)), function(j) {
  fj <- fit_limma(eoc_mat[, -j, drop = FALSE], droplevels(grp_eoc[-j]), "Active_EoC", "Non_EoC", trend = TRUE)$tt
  data.frame(
    omitted_sample = colnames(eoc_mat)[j],
    omitted_group = as.character(grp_eoc[j]),
    rho_discovery = rank_cor(fj, tt_u1)$rho,
    rho_validation = rank_cor(fj, tt_u2)$rho
  )
}))
write.csv(loo, file.path(out, "EoC_leave_one_out_primary_rho.csv"), row.names = FALSE)

primary_summary <- data.frame(
  comparison = c("EoC_vs_GSE75214_discovery", "EoC_vs_GSE16879_within_center_validation"),
  n_common = c(length(common_eoc_u1), length(common_eoc_u2)),
  rho_moderated_t = c(rc_eoc_u1$rho, rc_eoc_u2$rho),
  empirical_p_two_sided = c(p_eoc_u1, p_eoc_u2),
  permutations = c(length(rho_eoc_u1_null), length(rho_eoc_u2_null)),
  signature_genes_discovery = c(length(sig_eoc), NA)
)
write.csv(primary_summary, file.path(out, "primary_EoC_UC_summary.csv"), row.names = FALSE)

# Hallmark discovery-validation hierarchy for the primary axis
cat("Primary Hallmark camera analyses...\n")
hm <- msigdbr(species = "Homo sapiens", collection = "H")
hm_list <- split(hm$gene_symbol, hm$gs_name)
cam_eoc <- camera_table(eoc_mat, grp_eoc, "Active_EoC", "Non_EoC", hm_list, trend = TRUE)
cam_u1 <- camera_table(u1_mat, grp_u1, "UC", "CTRL", hm_list, trend = FALSE)
cam_u2 <- camera_table(u2_mat, grp_u2, "UC", "CTRL", hm_list, trend = FALSE)
write.csv(cam_eoc, file.path(out, "Hallmark_camera_EoC.csv"), row.names = FALSE)
write.csv(cam_u1, file.path(out, "Hallmark_camera_GSE75214_UC.csv"), row.names = FALSE)
write.csv(cam_u2, file.path(out, "Hallmark_camera_GSE16879_UC.csv"), row.names = FALSE)
hm_merge <- merge(cam_eoc, cam_u1, by = "pathway", suffixes = c("_EoC", "_UC_discovery"))
hm_merge <- merge(hm_merge, cam_u2, by = "pathway")
names(hm_merge)[names(hm_merge) %in% c("NGenes", "Direction", "PValue", "FDR")] <- paste0(names(hm_merge)[names(hm_merge) %in% c("NGenes", "Direction", "PValue", "FDR")], "_UC_validation")
hm_merge$discovery_candidate <- hm_merge$FDR_EoC < 0.05 & hm_merge$FDR_UC_discovery < 0.05 & hm_merge$Direction_EoC == hm_merge$Direction_UC_discovery
hm_merge$validation_holm <- NA_real_
if (any(hm_merge$discovery_candidate)) {
  hm_merge$validation_holm[hm_merge$discovery_candidate] <- p.adjust(hm_merge$PValue_UC_validation[hm_merge$discovery_candidate], method = "holm")
}
hm_merge$replicated <- hm_merge$discovery_candidate & hm_merge$Direction_EoC == hm_merge$Direction_UC_validation & hm_merge$validation_holm < 0.05
write.csv(hm_merge, file.path(out, "Hallmark_primary_discovery_validation.csv"), row.names = FALSE)

# Exploratory EoG analysis. It is never substituted for the EoC primary endpoint.
cat("Exploratory EoG-UC discovery permutation (10000 Monte Carlo assignments)...\n")
rc_eog_u1 <- rank_cor(tt_eog, tt_u1)
common_eog_u1 <- rc_eog_u1$gene
u1_t_for_eog <- rc_eog_u1$y
eog_null_file <- file.path(out, "EoG_UC_exploratory_perm_null.csv")
if (file.exists(eog_null_file)) {
  rho_eog_u1_null <- read.csv(eog_null_file)$rho
  cat("  Reusing completed exploratory EoG null from checkpoint.\n")
} else {
  set.seed(20260824)
  rho_eog_u1_null <- numeric(nperm)
  for (i in seq_len(nperm)) {
    idx <- sample(seq_len(ncol(eog_mat)), sum(grp_eog == "EG"))
    gp <- rep("NonEG", ncol(eog_mat)); gp[idx] <- "EG"
    fp <- fit_limma(eog_mat, factor(gp, levels = c("NonEG", "EG")), "EG", "NonEG", trend = TRUE)$tt
    xt <- fp$t[match(common_eog_u1, fp$gene)]
    rho_eog_u1_null[i] <- cor(xt, u1_t_for_eog, method = "spearman")
  }
  write.csv(data.frame(rho = rho_eog_u1_null), eog_null_file, row.names = FALSE)
}
p_eog_u1 <- (sum(abs(rho_eog_u1_null) >= abs(rc_eog_u1$rho)) + 1) / (nperm + 1)

rc_eog_u2 <- rank_cor(tt_eog, tt_u2)
sig_eog_tbl <- same_direction_signature(tt_eog, tt_u1)
write.csv(sig_eog_tbl, file.path(out, "EoG_UC_exploratory_signature_gate_full.csv"), row.names = FALSE)
sig_eog <- sig_eog_tbl$gene[sig_eog_tbl$selected]
write.csv(sig_eog_tbl[sig_eog_tbl$selected, ], file.path(out, "EoG_UC_exploratory_signature_selected.csv"), row.names = FALSE)

eog_summary <- data.frame(
  comparison = c("EoG_vs_GSE75214_discovery", "EoG_vs_GSE16879_within_center"),
  n_common = c(length(rc_eog_u1$gene), length(rc_eog_u2$gene)),
  rho_moderated_t = c(rc_eog_u1$rho, rc_eog_u2$rho),
  empirical_p_two_sided = c(p_eog_u1, NA),
  signature_genes_discovery = c(length(sig_eog), NA)
)

if (length(sig_eog) >= 10) {
  st <- sig_eog_tbl[sig_eog_tbl$selected, ]
  up <- st$gene[st$logFC_a > 0]
  down <- st$gene[st$logFC_a < 0]

  # Within-center sample-level validation
  score_u2 <- directional_score(u2_mat, up, down)
  hg_u2 <- hedges_g(score_u2, grp_u2)
  perm_u2 <- score_label_perm(score_u2, grp_u2, nperm = 10000, seed = 20260824)
  hg_u2$mean_difference <- perm_u2$observed
  hg_u2$permutation_p_two_sided <- perm_u2$p
  hg_u2$n_signature_total <- length(sig_eog)
  hg_u2$n_signature_scored <- length(intersect(rownames(u2_mat), sig_eog))
  write.csv(hg_u2, file.path(out, "EoG_signature_GSE16879_sample_validation.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = names(score_u2), group = as.character(grp_u2), score = as.numeric(score_u2)), file.path(out, "EoG_signature_GSE16879_scores.csv"), row.names = FALSE)
  write.csv(data.frame(mean_difference = perm_u2$null), file.path(out, "EoG_signature_GSE16879_label_perm_null.csv"), row.names = FALSE)

  # Conditional cross-center sensitivity: GSE36807 only
  e368_all <- read.csv(file.path(src, "results_v3.1", "sensitivity_2026-08-27", "expr_GSE36807_genes.csv"), row.names = 1, check.names = FALSE)
  m368 <- man[man$dataset == "GSE36807", ]
  r368 <- m368$primary_role[match(colnames(e368_all), m368$sample_id)]
  use368 <- r368 %in% c("sensitivity_colon_UC", "sensitivity_colon_control")
  e368 <- e368_all[, use368, drop = FALSE]
  g368 <- factor(ifelse(r368[use368] == "sensitivity_colon_UC", "UC", "CTRL"), levels = c("CTRL", "UC"))
  score368 <- directional_score(e368, up, down)
  hg368 <- hedges_g(score368, g368)
  perm368 <- score_label_perm(score368, g368, nperm = 10000, seed = 20260825)
  hg368$mean_difference <- perm368$observed
  hg368$permutation_p_two_sided <- perm368$p
  hg368$n_signature_total <- length(sig_eog)
  hg368$n_signature_scored <- length(intersect(rownames(e368), sig_eog))
  write.csv(hg368, file.path(out, "EoG_signature_GSE36807_conditional_sensitivity.csv"), row.names = FALSE)
  write.csv(data.frame(sample_id = names(score368), group = as.character(g368), score = as.numeric(score368)), file.path(out, "EoG_signature_GSE36807_scores.csv"), row.names = FALSE)

  # Explicitly post hoc phenotype association in UC only
  e109_all <- read.csv(file.path(src, "results_v3.1", "GSE109142_phenotype_2026-08-27", "GSE109142_expr_genes.csv"), row.names = 1, check.names = FALSE)
  sample_id_109 <- sub("_.*$", "", colnames(e109_all))
  m109 <- man[man$dataset == "GSE109142", ]
  role109 <- m109$primary_role[match(sample_id_109, m109$sample_id)]
  use109 <- role109 == "phenotype_validation_UC"
  e109 <- log2(e109_all[, use109, drop = FALSE] + 1)
  score109 <- directional_score(e109, up, down)
  sid109 <- sample_id_109[use109]
  dat109 <- data.frame(
    sample_id = sid109,
    score = as.numeric(score109),
    eos32 = factor(m109$eosinophil_grade_32[match(sid109, m109$sample_id)]),
    histology = suppressWarnings(as.numeric(m109$histology_severity[match(sid109, m109$sample_id)])),
    age = suppressWarnings(as.numeric(m109$age[match(sid109, m109$sample_id)])),
    sex = factor(m109$sex[match(sid109, m109$sample_id)])
  )
  cc <- complete.cases(dat109[, c("score", "eos32", "histology", "age", "sex")])
  fit109 <- lm(score ~ eos32 + histology + age + sex, data = dat109[cc, ])
  vc <- vcovHC(fit109, type = "HC3")
  ct <- coeftest(fit109, vcov. = vc)
  ci <- coef(fit109) + outer(sqrt(diag(vc)), c(-1.96, 1.96))
  ph <- data.frame(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    robust_se = ct[, "Std. Error"],
    ci_low = ci[, 1],
    ci_high = ci[, 2],
    p_value = ct[, "Pr(>|t|)"],
    n_complete = sum(cc),
    stringsAsFactors = FALSE
  )
  write.csv(dat109, file.path(out, "EoG_signature_GSE109142_exploratory_scores.csv"), row.names = FALSE)
  write.csv(ph, file.path(out, "EoG_signature_GSE109142_exploratory_model.csv"), row.names = FALSE)

  # Correct ORA: all eligible sets retained in the BH family, including zero-overlap sets.
  bg <- intersect(tt_eog$gene, tt_u1$gene)
  ora_h <- run_ora(sig_eog, bg, "H", NULL, "EoG_UC_Hallmark")
  ora_k <- run_ora(sig_eog, bg, "C2", "CP:KEGG_LEGACY", "EoG_UC_KEGG")
  ora_g <- run_ora(sig_eog, bg, "C5", "GO:BP", "EoG_UC_GO_BP")
}

write.csv(eog_summary, file.path(out, "exploratory_EoG_UC_summary.csv"), row.names = FALSE)

# Analysis-specific sample accounting (not the 683-row metadata inventory)
sample_accounting <- data.frame(
  dataset = c("EGIDExpress EoC", "EGIDExpress EoG", "GSE75214", "GSE16879", "GSE109142", "GSE9686", "GSE36807"),
  role = c("confirmatory disease anchor", "post hoc exploratory cross-organ anchor", "UC discovery", "within-center conditional UC validation", "conditional post hoc phenotype analysis", "QC failure; not used inferentially", "conditional sensitivity"),
  case_n = c(6, 9, 74, 24, 206, 5, 15),
  comparator_n = c(8, 12, 11, 6, 0, 8, 7),
  inferential_use = c(TRUE, TRUE, TRUE, TRUE, length(sig_eog) >= 10, FALSE, length(sig_eog) >= 10)
)
write.csv(sample_accounting, file.path(out, "analysis_specific_sample_accounting.csv"), row.names = FALSE)

writeLines(c(
  paste0("EoC discovery rho=", signif(rc_eoc_u1$rho, 5), ", exact two-sided P=", signif(p_eoc_u1, 5)),
  paste0("EoC validation rho=", signif(rc_eoc_u2$rho, 5), ", permutation two-sided P=", signif(p_eoc_u2, 5)),
  paste0("EoC gated signature genes=", length(sig_eoc)),
  paste0("EoG exploratory discovery rho=", signif(rc_eog_u1$rho, 5), ", permutation two-sided P=", signif(p_eog_u1, 5)),
  paste0("EoG gated exploratory signature genes=", length(sig_eog)),
  paste0("Primary replicated Hallmark pathways=", sum(hm_merge$replicated, na.rm = TRUE))
), file.path(out, "headline_results.txt"))

write_session()
cat("Corrective reanalysis complete. Outputs:", out, "\n")

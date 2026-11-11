# Matched-random-gene-set competitive test for the exploratory EoG-UC signature
# Protocol v1.0 plus prospective scoreability amendment v1.1
# Frozen analysis date: 2026-08-28

suppressPackageStartupMessages({
  library(jsonlite)
  library(digest)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

src <- Sys.getenv("EOGE_IBD_SOURCE", unset = "")
if (!nzchar(src)) stop("Set EOGE_IBD_SOURCE to the input directory described in docs/INPUTS.md")
root <- normalizePath(Sys.getenv("EOGE_IBD_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
base_dir <- file.path(root, "analysis", "matched_null_2026-08-28")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(base_dir, "figures")
qa_dir <- file.path(base_dir, "qa")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(base_dir, "logs"), recursive=TRUE, showWarnings=FALSE)
log_con <- file(file.path(base_dir, "logs", "run_matched_null_v1.1.log"), open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

config_file <- file.path(base_dir, "config_v1.1.json")
protocol_file <- file.path(base_dir, "frozen_protocol_v1.0_2026-08-28.md")
amendment_file <- file.path(base_dir, "protocol_amendment_v1.1_2026-08-28.md")
cfg <- fromJSON(config_file)

expected_hashes <- c(
  frozen_protocol_v1.0_2026.08.28.md = "F3D12D160B54E02A903E7C08CF2D1F6555DE776A4D5BF675585015DB67DC27A1",
  protocol_amendment_v1.1_2026.08.28.md = "1474F43DAA94BAD74974CADFDA9A63BC10CDCA730F0ECFB6BCC6FEF7FC97D1A2",
  config_v1.1.json = "BFA704AB6A705648557894AC4D3616FA70EA01DDC702ED27A8FEA219168394C1"
)
actual_hashes <- toupper(c(
  digest(file = protocol_file, algo = "sha256", serialize = FALSE),
  digest(file = amendment_file, algo = "sha256", serialize = FALSE),
  digest(file = config_file, algo = "sha256", serialize = FALSE)
))
names(actual_hashes) <- names(expected_hashes)
if (!identical(actual_hashes, expected_hashes)) {
  stop("Frozen protocol/config hash mismatch; analysis aborted.")
}
write.csv(
  data.frame(file = names(actual_hashes), sha256 = unname(actual_hashes), verified = actual_hashes == expected_hashes),
  file.path(qa_dir, "prefreeze_hash_verification.csv"), row.names = FALSE
)

stopifnot(
  identical(cfg$protocol_version, "1.1"),
  cfg$seed == 20260828,
  cfg$n_random_sets == 10000,
  cfg$nearest_k == 100,
  cfg$fallback_k == 500,
  isTRUE(cfg$require_positive_variance_in_GSE36807),
  isTRUE(cfg$require_positive_variance_in_GSE109142)
)

rank01 <- function(x) {
  (rank(x, ties.method = "average") - 0.5) / length(x)
}

row_sd <- function(mat) {
  sqrt(pmax(0, rowSums((mat - rowMeans(mat))^2) / (ncol(mat) - 1)))
}

zscore_rows_strict <- function(mat) {
  mu <- rowMeans(mat)
  s <- row_sd(mat)
  if (any(!is.finite(mu)) || any(!is.finite(s)) || any(s <= 0) || any(!is.finite(mat))) {
    stop("Non-finite values or non-positive variance reached strict z-score function.")
  }
  sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
}

hedges_g_rows <- function(score_matrix, group, case = "UC", ref = "CTRL") {
  case_mat <- score_matrix[, group == case, drop = FALSE]
  ref_mat <- score_matrix[, group == ref, drop = FALSE]
  n1 <- ncol(case_mat)
  n0 <- ncol(ref_mat)
  m1 <- rowMeans(case_mat)
  m0 <- rowMeans(ref_mat)
  v1 <- rowSums((case_mat - m1)^2) / (n1 - 1)
  v0 <- rowSums((ref_mat - m0)^2) / (n0 - 1)
  sp <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
  d <- (m1 - m0) / sp
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  J * d
}

hedges_g_summary <- function(score, group, case = "UC", ref = "CTRL") {
  x1 <- score[group == case]
  x0 <- score[group == ref]
  n1 <- length(x1)
  n0 <- length(x0)
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
    mean_difference = mean(x1) - mean(x0),
    hedges_g = g,
    se_g = sqrt(vg),
    ci_low = g - 1.96 * sqrt(vg),
    ci_high = g + 1.96 * sqrt(vg)
  )
}

prepare_hc3 <- function(covariate_data) {
  X <- model.matrix(~ eos32 + histology + age + sex, data = covariate_data)
  inv <- solve(crossprod(X))
  Hdiag <- rowSums((X %*% inv) * X)
  term_idx <- match("histology", colnames(X))
  if (is.na(term_idx)) stop("Histology term missing from frozen model matrix.")
  a <- as.numeric(X %*% inv[, term_idx])
  list(X = X, inv = inv, Hdiag = Hdiag, term_idx = term_idx, a = a, df = nrow(X) - ncol(X))
}

hc3_histology_rows <- function(score_matrix, prep) {
  Y <- t(score_matrix)
  beta <- prep$inv %*% crossprod(prep$X, Y)
  residual <- Y - prep$X %*% beta
  adjusted_residual <- residual / (1 - prep$Hdiag)
  se <- sqrt(colSums(adjusted_residual^2 * (prep$a^2)))
  estimate <- as.numeric(beta[prep$term_idx, ])
  t_value <- estimate / se
  p_value <- 2 * pt(abs(t_value), df = prep$df, lower.tail = FALSE)
  data.frame(estimate = estimate, robust_se = se, t_value = t_value, p_value_two_sided = p_value)
}

make_candidate_lists <- function(props, signature_indices, slot_direction, feature_names) {
  candidate_ok <- !props$is_signature
  out <- vector("list", length(signature_indices))
  for (i in seq_along(signature_indices)) {
    candidates <- which(candidate_ok & props$direction == slot_direction[i])
    real_vector <- as.numeric(props[signature_indices[i], feature_names, drop = TRUE])
    candidate_matrix <- as.matrix(props[candidates, feature_names, drop = FALSE])
    distance <- sqrt(rowSums((candidate_matrix - matrix(real_vector, nrow(candidate_matrix), length(real_vector), byrow = TRUE))^2))
    ord <- order(distance, props$gene[candidates])
    out[[i]] <- list(index = candidates[ord], distance = distance[ord])
  }
  out
}

draw_matched_sets <- function(candidate_lists, n_props, n_sets, nearest_k, fallback_k, seed, label) {
  set.seed(seed)
  n_slots <- length(candidate_lists)
  draws <- matrix(NA_integer_, nrow = n_slots, ncol = n_sets)
  failed <- logical(n_sets)
  used_marker <- integer(n_props)
  for (b in seq_len(n_sets)) {
    slot_order <- sample.int(n_slots)
    for (slot in slot_order) {
      idx_all <- candidate_lists[[slot]]$index
      dist_all <- candidate_lists[[slot]]$distance
      chosen_pool <- integer(0)
      for (limit in c(nearest_k, fallback_k, length(idx_all))) {
        lim <- min(limit, length(idx_all))
        probe <- seq_len(lim)
        available <- probe[used_marker[idx_all[probe]] != b]
        if (length(available) > 0) {
          chosen_pool <- available
          break
        }
      }
      if (length(chosen_pool) == 0) {
        failed[b] <- TRUE
        break
      }
      scale_distance <- median(dist_all[seq_len(min(nearest_k, length(dist_all)))])
      scale_distance <- max(scale_distance, 1e-6)
      weights <- exp(-dist_all[chosen_pool] / scale_distance)
      if (!all(is.finite(weights)) || sum(weights) <= 0) weights <- rep(1, length(chosen_pool))
      pick_position <- chosen_pool[sample.int(length(chosen_pool), 1L, prob = weights)]
      pick_index <- idx_all[pick_position]
      draws[slot, b] <- pick_index
      used_marker[pick_index] <- b
    }
    if (b %% 1000 == 0) cat(sprintf("  %s: generated %d/%d sets\n", label, b, n_sets))
  }
  list(draws = draws, failed = failed)
}

balance_diagnostics <- function(props, signature_indices, draws, feature_names, model_name) {
  universe_sd <- vapply(feature_names, function(f) sd(props[[f]]), numeric(1))
  replicate_rows <- vector("list", length(feature_names))
  summary_rows <- vector("list", length(feature_names))
  for (i in seq_along(feature_names)) {
    f <- feature_names[i]
    real_mean <- mean(props[[f]][signature_indices])
    value_matrix <- matrix(props[[f]][draws], nrow = nrow(draws), ncol = ncol(draws))
    smd <- (colMeans(value_matrix) - real_mean) / universe_sd[f]
    replicate_rows[[i]] <- data.frame(
      null_model = model_name,
      replicate = seq_along(smd),
      feature = f,
      standardized_mean_difference = smd
    )
    med <- median(abs(smd))
    p95 <- unname(quantile(abs(smd), 0.95, type = 8))
    summary_rows[[i]] <- data.frame(
      null_model = model_name,
      feature = f,
      real_mean = real_mean,
      universe_sd = universe_sd[f],
      median_abs_smd = med,
      p95_abs_smd = p95,
      median_gate = med <= cfg$balance_median_abs_smd_max,
      p95_gate = p95 <= cfg$balance_p95_abs_smd_max
    )
  }
  list(replicate = do.call(rbind, replicate_rows), summary = do.call(rbind, summary_rows))
}

score_draws <- function(draws, z_matrix, slot_weights) {
  n_sets <- ncol(draws)
  n_slots <- nrow(draws)
  W <- sparseMatrix(
    i = rep(seq_len(n_sets), each = n_slots),
    j = as.vector(draws),
    x = rep(slot_weights, n_sets),
    dims = c(n_sets, nrow(z_matrix))
  )
  as.matrix(W %*% z_matrix)
}

competitive_p <- function(null, observed) {
  c(
    p_one_sided = (1 + sum(null >= observed)) / (length(null) + 1),
    p_two_sided = (1 + sum(abs(null) >= abs(observed))) / (length(null) + 1),
    observed_percentile = (sum(null < observed) + 0.5 * sum(null == observed)) / length(null)
  )
}

cat("Loading frozen inputs...\n")
man <- read.delim(file.path(src, "manifests", "sample_manifest_v1.0.tsv"), check.names = FALSE)
names(man)[1] <- "dataset"

signature <- read.csv(file.path(root, "analysis", "corrected_results", "EoG_UC_exploratory_signature_selected.csv"), check.names = FALSE)
uc_de <- read.csv(file.path(root, "analysis", "corrected_results", "GSE75214_UC_DEG_corrected.csv"), check.names = FALSE)

u1_all <- read.csv(file.path(src, "results_v3.1", "formal_DEG_2026-08-27", "expr_GSE75214_genes.csv"), row.names = 1, check.names = FALSE)
m1 <- man[man$dataset == "GSE75214", ]
role1 <- m1$primary_role[match(colnames(u1_all), m1$sample_id)]
use1 <- role1 %in% c("discovery_active_colon_UC", "discovery_colon_control")
u1 <- as.matrix(u1_all[, use1, drop = FALSE])
storage.mode(u1) <- "double"

e368_all <- read.csv(file.path(src, "results_v3.1", "sensitivity_2026-08-27", "expr_GSE36807_genes.csv"), row.names = 1, check.names = FALSE)
m368 <- man[man$dataset == "GSE36807", ]
role368 <- m368$primary_role[match(colnames(e368_all), m368$sample_id)]
use368 <- role368 %in% c("sensitivity_colon_UC", "sensitivity_colon_control")
e368 <- as.matrix(e368_all[, use368, drop = FALSE])
storage.mode(e368) <- "double"
group368 <- factor(ifelse(role368[use368] == "sensitivity_colon_UC", "UC", "CTRL"), levels = c("CTRL", "UC"))
stopifnot(sum(group368 == "UC") == 15, sum(group368 == "CTRL") == 7)

e109_all <- read.csv(file.path(src, "results_v3.1", "GSE109142_phenotype_2026-08-27", "GSE109142_expr_genes.csv"), row.names = 1, check.names = FALSE)
sample_id109_all <- sub("_.*$", "", colnames(e109_all))
m109 <- man[man$dataset == "GSE109142", ]
role109 <- m109$primary_role[match(sample_id109_all, m109$sample_id)]
use109 <- role109 == "phenotype_validation_UC"
e109 <- log2(as.matrix(e109_all[, use109, drop = FALSE]) + 1)
storage.mode(e109) <- "double"
sample_id109 <- sample_id109_all[use109]
stopifnot(ncol(e109) == 206)

common <- Reduce(intersect, list(rownames(u1), rownames(e368), rownames(e109), uc_de$gene))
uc_match <- match(common, uc_de$gene)
u1_sd <- row_sd(u1[common, , drop = FALSE])
e368_sd <- row_sd(e368[common, , drop = FALSE])
e109_sd <- row_sd(e109[common, , drop = FALSE])
eligible <-
  is.finite(rowMeans(u1[common, , drop = FALSE])) &
  is.finite(u1_sd) & u1_sd > 0 &
  is.finite(e368_sd) & e368_sd > 0 &
  is.finite(e109_sd) & e109_sd > 0 &
  is.finite(uc_de$t[uc_match]) &
  is.finite(uc_de$logFC[uc_match]) & uc_de$logFC[uc_match] != 0
common <- common[eligible]
uc_match <- match(common, uc_de$gene)

props <- data.frame(
  gene = common,
  discovery_mean = rowMeans(u1[common, , drop = FALSE]),
  discovery_log_sd = log2(row_sd(u1[common, , drop = FALSE]) + 1e-8),
  uc_moderated_t = uc_de$t[uc_match],
  uc_logFC = uc_de$logFC[uc_match],
  direction = ifelse(uc_de$logFC[uc_match] > 0, "UP", "DOWN"),
  stringsAsFactors = FALSE
)
props$mean_rank <- rank01(props$discovery_mean)
props$log_sd_rank <- rank01(props$discovery_log_sd)
props$abs_t_rank <- rank01(abs(props$uc_moderated_t))
props$is_signature <- props$gene %in% signature$gene

covered_signature <- signature$gene[signature$gene %in% props$gene]
signature_indices <- match(covered_signature, props$gene)
slot_direction <- props$direction[signature_indices]
coverage <- length(covered_signature) / nrow(signature)

coverage_table <- data.frame(
  signature_gene = signature$gene,
  covered_common_universe = signature$gene %in% covered_signature,
  uc_direction = props$direction[match(signature$gene, props$gene)],
  stringsAsFactors = FALSE
)
write.csv(coverage_table, file.path(result_dir, "signature_coverage.csv"), row.names = FALSE)
write.csv(props, file.path(result_dir, "common_universe_properties.csv"), row.names = FALSE)

universe_summary <- data.frame(
  n_signature_original = nrow(signature),
  n_signature_covered = length(covered_signature),
  signature_coverage = coverage,
  n_signature_up = sum(slot_direction == "UP"),
  n_signature_down = sum(slot_direction == "DOWN"),
  n_common_universe = nrow(props),
  n_candidate_up = sum(!props$is_signature & props$direction == "UP"),
  n_candidate_down = sum(!props$is_signature & props$direction == "DOWN"),
  n_GSE36807_UC = sum(group368 == "UC"),
  n_GSE36807_control = sum(group368 == "CTRL"),
  n_GSE109142_UC = ncol(e109)
)
write.csv(universe_summary, file.path(result_dir, "analysis_universe_summary.csv"), row.names = FALSE)

if (coverage < cfg$minimum_signature_coverage) {
  stop(sprintf("Signature coverage %.3f is below the frozen %.3f gate.", coverage, cfg$minimum_signature_coverage))
}
if (sum(!props$is_signature & props$direction == "UP") < sum(slot_direction == "UP") ||
    sum(!props$is_signature & props$direction == "DOWN") < sum(slot_direction == "DOWN")) {
  stop("Candidate direction pool is too small for unique exact-size sets.")
}

cat(sprintf("Eligible universe: %d genes; covered signature: %d/%d (%.1f%%), %d up/%d down.\n",
            nrow(props), length(covered_signature), nrow(signature), 100 * coverage,
            sum(slot_direction == "UP"), sum(slot_direction == "DOWN")))

z368 <- zscore_rows_strict(e368[props$gene, , drop = FALSE])
z109 <- zscore_rows_strict(e109[props$gene, , drop = FALSE])
stopifnot(identical(rownames(z368), props$gene), identical(rownames(z109), props$gene))

n_up <- sum(slot_direction == "UP")
n_down <- sum(slot_direction == "DOWN")
slot_weights <- ifelse(slot_direction == "UP", 1 / n_up, -1 / n_down)
real_score368 <- as.numeric(slot_weights %*% z368[signature_indices, , drop = FALSE])
names(real_score368) <- colnames(z368)
real_score109 <- as.numeric(slot_weights %*% z109[signature_indices, , drop = FALSE])
names(real_score109) <- colnames(z109)

cov109 <- data.frame(
  sample_id = sample_id109,
  eos32 = factor(m109$eosinophil_grade_32[match(sample_id109, m109$sample_id)]),
  histology = suppressWarnings(as.numeric(m109$histology_severity[match(sample_id109, m109$sample_id)])),
  age = suppressWarnings(as.numeric(m109$age[match(sample_id109, m109$sample_id)])),
  sex = factor(m109$sex[match(sample_id109, m109$sample_id)]),
  stringsAsFactors = FALSE
)
cc109 <- complete.cases(cov109[, c("eos32", "histology", "age", "sex")])
cov109_cc <- droplevels(cov109[cc109, , drop = FALSE])
prep109 <- prepare_hc3(cov109_cc)

real_g_summary <- hedges_g_summary(real_score368, group368)
real_hc3 <- hc3_histology_rows(matrix(real_score109[cc109], nrow = 1), prep109)
real_hc3$ci_low <- real_hc3$estimate - 1.96 * real_hc3$robust_se
real_hc3$ci_high <- real_hc3$estimate + 1.96 * real_hc3$robust_se
real_hc3$n_complete <- sum(cc109)
write.csv(real_g_summary, file.path(result_dir, "GSE36807_real_signature_effect.csv"), row.names = FALSE)
write.csv(real_hc3, file.path(result_dir, "GSE109142_real_signature_histology_HC3.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = names(real_score368), group = as.character(group368), score = real_score368),
          file.path(result_dir, "GSE36807_real_signature_scores.csv"), row.names = FALSE)
write.csv(cbind(cov109, score = real_score109), file.path(result_dir, "GSE109142_real_signature_scores.csv"), row.names = FALSE)

models <- list(
  technical = c("mean_rank", "log_sd_rank"),
  uc_effect = c("mean_rank", "log_sd_rank", "abs_t_rank")
)

all_balance_replicate <- list()
all_balance_summary <- list()
all_null_stats <- list()
all_competitive <- list()
model_gate_rows <- list()
draw_objects <- list()

for (model_name in names(models)) {
  cat(sprintf("Building %s candidate lists...\n", model_name))
  feature_names <- models[[model_name]]
  candidate_lists <- make_candidate_lists(props, signature_indices, slot_direction, feature_names)
  draw_result <- draw_matched_sets(
    candidate_lists = candidate_lists,
    n_props = nrow(props),
    n_sets = cfg$n_random_sets,
    nearest_k = cfg$nearest_k,
    fallback_k = cfg$fallback_k,
    seed = cfg$seed + match(model_name, names(models)) - 1,
    label = model_name
  )
  failed_fraction <- mean(draw_result$failed)
  if (failed_fraction > cfg$maximum_failed_draw_fraction) {
    stop(sprintf("%s failed-draw fraction %.4f exceeded the frozen gate.", model_name, failed_fraction))
  }
  valid_draws <- draw_result$draws[, !draw_result$failed, drop = FALSE]
  if (ncol(valid_draws) != cfg$n_random_sets) {
    stop("Any failed draw changes the frozen Monte Carlo denominator; analysis aborted for amendment.")
  }
  draw_objects[[model_name]] <- list(
    props_gene = props$gene,
    signature_slot_gene = covered_signature,
    signature_slot_direction = slot_direction,
    feature_names = feature_names,
    draws = valid_draws
  )
  saveRDS(draw_objects[[model_name]], file.path(result_dir, paste0("matched_gene_sets_", model_name, ".rds")), compress = "xz")

  bal <- balance_diagnostics(props, signature_indices, valid_draws, feature_names, model_name)
  all_balance_replicate[[model_name]] <- bal$replicate
  all_balance_summary[[model_name]] <- bal$summary
  balance_pass <- all(bal$summary$median_gate & bal$summary$p95_gate)

  cat(sprintf("Scoring %s null sets in GSE36807 and GSE109142...\n", model_name))
  null_score368 <- score_draws(valid_draws, z368, slot_weights)
  null_g <- hedges_g_rows(null_score368, group368)
  rm(null_score368)
  gc(verbose = FALSE)

  null_score109 <- score_draws(valid_draws, z109[, cc109, drop = FALSE], slot_weights)
  null_hc3 <- hc3_histology_rows(null_score109, prep109)
  rm(null_score109)
  gc(verbose = FALSE)

  null_stats <- rbind(
    data.frame(null_model = model_name, endpoint = "GSE36807_Hedges_g", replicate = seq_along(null_g), statistic = null_g),
    data.frame(null_model = model_name, endpoint = "GSE109142_histology_HC3_t", replicate = seq_len(nrow(null_hc3)), statistic = null_hc3$t_value)
  )
  all_null_stats[[model_name]] <- null_stats

  p_g <- competitive_p(null_g, real_g_summary$hedges_g)
  p_h <- competitive_p(null_hc3$t_value, real_hc3$t_value)
  comp <- rbind(
    data.frame(null_model = model_name, endpoint = "GSE36807_Hedges_g", observed = real_g_summary$hedges_g,
               p_one_sided = p_g["p_one_sided"], p_two_sided = p_g["p_two_sided"], observed_percentile = p_g["observed_percentile"]),
    data.frame(null_model = model_name, endpoint = "GSE109142_histology_HC3_t", observed = real_hc3$t_value,
               p_one_sided = p_h["p_one_sided"], p_two_sided = p_h["p_two_sided"], observed_percentile = p_h["observed_percentile"])
  )
  comp$holm_p_one_sided <- p.adjust(comp$p_one_sided, method = "holm")
  comp$positive_observed <- comp$observed > 0
  comp$endpoint_pass <- comp$positive_observed & comp$holm_p_one_sided < 0.05
  all_competitive[[model_name]] <- comp

  model_gate_rows[[model_name]] <- data.frame(
    null_model = model_name,
    n_requested = cfg$n_random_sets,
    n_valid = ncol(valid_draws),
    failed_fraction = failed_fraction,
    coverage_pass = coverage >= cfg$minimum_signature_coverage,
    balance_pass = balance_pass,
    both_endpoints_pass = all(comp$endpoint_pass),
    model_pass = coverage >= cfg$minimum_signature_coverage & balance_pass & all(comp$endpoint_pass)
  )
}

balance_replicate <- do.call(rbind, all_balance_replicate)
balance_summary <- do.call(rbind, all_balance_summary)
null_stats_all <- do.call(rbind, all_null_stats)
competitive_results <- do.call(rbind, all_competitive)
model_gates <- do.call(rbind, model_gate_rows)

write.csv(balance_replicate, file.path(result_dir, "matching_balance_by_replicate.csv"), row.names = FALSE)
write.csv(balance_summary, file.path(result_dir, "matching_balance_summary.csv"), row.names = FALSE)
write.csv(null_stats_all, file.path(result_dir, "null_endpoint_statistics.csv"), row.names = FALSE)
write.csv(competitive_results, file.path(result_dir, "competitive_test_results.csv"), row.names = FALSE)
write.csv(model_gates, file.path(result_dir, "frozen_model_gates.csv"), row.names = FALSE)

technical_pass <- model_gates$model_pass[model_gates$null_model == "technical"]
uc_effect_pass <- model_gates$model_pass[model_gates$null_model == "uc_effect"]
decision <- if (isTRUE(technical_pass) && isTRUE(uc_effect_pass)) {
  "PROMOTE_SUBSTANTIVE_EXPLORATORY"
} else if (isTRUE(technical_pass) && !isTRUE(uc_effect_pass)) {
  "SHORT_MAIN_OR_SUPPLEMENT_GENERIC_UC_NOT_EXCLUDED"
} else if (!isTRUE(technical_pass) && !isTRUE(uc_effect_pass)) {
  "SUPPLEMENT_ONLY"
} else {
  "INCONSISTENT_REVIEW_REQUIRED"
}
decision_table <- data.frame(
  decision = decision,
  technical_model_pass = technical_pass,
  uc_effect_model_pass = uc_effect_pass,
  persistent_evidence_ceiling = "exploratory association; no causal, diagnostic, or EoG-specific claim",
  stringsAsFactors = FALSE
)
write.csv(decision_table, file.path(result_dir, "frozen_decision.csv"), row.names = FALSE)

panel_levels <- c(
  "a  Technical-property null: GSE36807 Hedges g",
  "b  Technical-property null: GSE109142 histology HC3 t",
  "c  UC-effect-matched null: GSE36807 Hedges g",
  "d  UC-effect-matched null: GSE109142 histology HC3 t"
)
plot_data <- null_stats_all
plot_data$panel <- ifelse(
  plot_data$null_model == "technical" & plot_data$endpoint == "GSE36807_Hedges_g", panel_levels[1],
  ifelse(plot_data$null_model == "technical", panel_levels[2],
         ifelse(plot_data$endpoint == "GSE36807_Hedges_g", panel_levels[3], panel_levels[4]))
)
plot_data$panel <- factor(plot_data$panel, levels = panel_levels)
write.csv(plot_data[, c("panel", "null_model", "endpoint", "replicate", "statistic")],
          file.path(result_dir, "FigureS_matched_null_plotted_data.csv"), row.names = FALSE)

annotation <- competitive_results
annotation$panel <- ifelse(
  annotation$null_model == "technical" & annotation$endpoint == "GSE36807_Hedges_g", panel_levels[1],
  ifelse(annotation$null_model == "technical", panel_levels[2],
         ifelse(annotation$endpoint == "GSE36807_Hedges_g", panel_levels[3], panel_levels[4]))
)
annotation$panel <- factor(annotation$panel, levels = panel_levels)
annotation$label <- sprintf("Observed = %.2f\ncompetitive P = %.4g\nHolm P = %.4g",
                            annotation$observed, annotation$p_one_sided, annotation$holm_p_one_sided)

p <- ggplot(plot_data, aes(x = statistic)) +
  geom_histogram(bins = 55, fill = "#9BA3AA", color = "white", linewidth = 0.15) +
  geom_vline(data = annotation, aes(xintercept = observed), inherit.aes = FALSE,
             color = "#C43C39", linewidth = 0.9, linetype = "dashed") +
  geom_text(data = annotation, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.15, size = 3.2, lineheight = 1.05, color = "#8F2323") +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(
    title = "Competitive evaluation of the 209-gene EoG-UC signature",
    subtitle = "Observed statistics (red) versus 10,000 prespecified matched random gene sets",
    x = "Null-set statistic",
    y = "Count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5),
    strip.background = element_rect(fill = "#EEF1F4", color = "#6D747A", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(face = "bold"),
    panel.spacing = grid::unit(7, "pt")
  )

ggsave(file.path(figure_dir, "FigureS_matched_null_competitive.png"), p,
       width = 190, height = 155, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "FigureS_matched_null_competitive.tiff"), p,
       width = 190, height = 155, units = "mm", dpi = 600, bg = "white", compression = "lzw")
ggsave(file.path(figure_dir, "FigureS_matched_null_competitive.pdf"), p,
       width = 190 / 25.4, height = 155 / 25.4, units = "in", device = grDevices::cairo_pdf, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(qa_dir, "sessionInfo.txt"))
writeLines(
  c(
    paste0("protocol_version=", cfg$protocol_version),
    paste0("decision=", decision),
    paste0("signature_coverage=", format(coverage, digits = 6)),
    paste0("technical_model_pass=", technical_pass),
    paste0("uc_effect_model_pass=", uc_effect_pass)
  ),
  file.path(qa_dir, "run_completion_status.txt")
)

cat("Matched-null analysis complete. Frozen decision:", decision, "\n")

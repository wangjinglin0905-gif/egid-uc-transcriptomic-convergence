# Publication figure production for the corrected EoGE-IBD manuscript
# Version 1.0; 2026-08-27

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(ragg)
})

root <- normalizePath(Sys.getenv("EOGE_IBD_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
res_dir <- file.path(root, "analysis", "corrected_results")
fig_dir <- file.path(root,"generated_figures","main")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

col <- list(
  navy = "#1F4E79",
  blue = "#4C78A8",
  cyan = "#4C9FAD",
  vermilion = "#D1495B",
  orange = "#E68600",
  green = "#2A9D8F",
  purple = "#7251A3",
  grey = "#8C8C8C",
  light_grey = "#E7E7E7",
  dark = "#242424",
  pale_blue = "#DCEAF5",
  pale_orange = "#FBE7D6",
  pale_green = "#DDEFEA",
  pale_purple = "#E9E1F2",
  pale_red = "#F6DADB"
)

theme_pub <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = col$dark),
      axis.title = element_text(size = base_size + 0.5),
      axis.text = element_text(size = base_size - 0.5, colour = col$dark),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.35),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "#555555"),
      plot.caption = element_text(size = base_size - 1, colour = "#555555", hjust = 0),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_figure <- function(name, plot, width, height) {
  target <- commandArgs(trailingOnly=TRUE)
  if (length(target) && !grepl(target[1],name,fixed=TRUE)) return(invisible(NULL))
  png_path <- file.path(fig_dir, paste0(sub("_revised$", "", name), ".png"))
  tif_path <- file.path(fig_dir, paste0(sub("_revised$", "", name), ".tiff"))
  pdf_path <- file.path(fig_dir, paste0(name, ".pdf"))

  agg_png(png_path, width = width, height = height, units = "in", res = 600,
          background = "white", scaling = 1)
  set.seed(20260827)
  print(plot)
  dev.off()

  agg_tiff(tif_path, width = width, height = height, units = "in", res = 600,
           background = "white", compression = "lzw", bitsize = 8, scaling = 1)
  set.seed(20260827)
  print(plot)
  dev.off()

  # PNG and TIFF only; both renders use identical point-jitter seed.
}

clean_hm <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  replacements <- c(
    "Tnfa Signaling Via Nfkb" = "TNF-alpha signaling via NF-kappaB",
    "G2m Checkpoint" = "G2M checkpoint",
    "E2f Targets" = "E2F targets",
    "Uv Response Dn" = "UV response down",
    "Uv Response Up" = "UV response up",
    "Il6 Jak Stat3 Signaling" = "IL6-JAK-STAT3 signaling",
    "Il2 Stat5 Signaling" = "IL2-STAT5 signaling",
    "Pi3k Akt Mtor Signaling" = "PI3K-AKT-mTOR signaling",
    "Mtorc1 Signaling" = "mTORC1 signaling",
    "Myc Targets V1" = "MYC targets V1",
    "Myc Targets V2" = "MYC targets V2",
    "Tgf Beta Signaling" = "TGF-beta signaling",
    "P53 Pathway" = "p53 pathway"
  )
  for (old in names(replacements)) x[x == old] <- replacements[[old]]
  x
}

clean_term <- function(x) {
  x <- gsub("^(HALLMARK|KEGG)_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

signed_fdr <- function(direction, fdr, cap = 15) {
  s <- ifelse(direction == "Up", 1, -1)
  s * pmin(-log10(pmax(as.numeric(fdr), .Machine$double.xmin)), cap)
}

# -----------------------------------------------------------------------------
# Figure 1: evidence architecture
# -----------------------------------------------------------------------------
nodes <- data.frame(
  x = c(1.7, 1.7, 1.7, 5.0, 5.0, 5.0, 5.0, 5.0, 8.3, 8.3, 8.3),
  y = c(8.2, 5.4, 2.7, 8.2, 6.35, 4.5, 2.65, 0.75, 8.0, 5.5, 3.0),
  w = c(2.7, 2.7, 2.7, 2.8, 2.8, 2.8, 2.8, 2.8, 2.7, 2.7, 2.7),
  h = c(1.35, 1.35, 1.35, 1.35, 1.25, 1.25, 1.25, 0.95, 1.35, 1.35, 1.35),
  label = c(
    "EoC colon\n6 active / 8 reference\nPRIMARY ANCHOR",
    "EoG stomach\n9 disease / 12 reference\nPOST HOC ANCHOR",
    "EoD duodenum\n382 preselected genes\nDESCRIPTIVE ONLY",
    "GSE75214\n74 active colonic UC / 11 controls\nDISCOVERY",
    "GSE16879\n24 active colonic UC / 6 controls\nWITHIN-CENTER VALIDATION",
    "GSE36807\n15 UC / 7 controls\nCONDITIONAL SENSITIVITY",
    "GSE109142\n206 UC; 198 complete cases\nPOST HOC PHENOTYPE",
    "GSE9686\nphenotype = batch\nEXCLUDED",
    "Genome-wide rank correlation\nmoderated t statistics\nlabel-permutation inference",
    "Local gene and pathway analyses\nSame-direction conjunction: q < 0.05\nScore only if at least 10 genes",
    "After the signature-size gate\nScore, ORA and phenotype tests\nMatched-set specificity assessment"
  ),
  fill = c(
    col$pale_blue, col$pale_orange, "#F2F2F2",
    col$pale_green, col$pale_green, col$pale_purple, col$pale_purple, col$pale_red,
    "#EDF3F6", "#EDF3F6", "#EDF3F6"
  ),
  border = c(
    col$navy, col$orange, col$grey,
    col$green, col$green, col$purple, col$purple, col$vermilion,
    col$navy, col$navy, col$navy
  )
)

edges <- data.frame(
  x = c(3.05, 3.05, 6.4, 6.4, 6.4, 6.4, 8.3, 8.3),
  y = c(8.2, 5.4, 8.2, 6.35, 4.5, 2.65, 7.3, 4.8),
  xend = c(3.60, 3.60, 6.95, 6.95, 6.95, 6.95, 8.3, 8.3),
  yend = c(8.2, 7.82, 8.0, 5.5, 3.25, 3.0, 6.2, 3.7),
  lty = c("solid", "dashed", "solid", "solid", "dashed", "dashed", "solid", "solid")
)

p1 <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = lty),
    linewidth = 0.55, colour = "#6A6A6A",
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
  ) +
  scale_linetype_identity() +
  geom_rect(
    data = nodes,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
        fill = fill, colour = border),
    linewidth = 0.7
  ) +
  scale_fill_identity() +
  scale_colour_identity() +
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            family = "Arial", size = 2.55, lineheight = 0.96, colour = col$dark) +
  annotate("text", x = 1.7, y = 9.25, label = "EGID evidence", family = "Arial",
           fontface = "bold", size = 3.2, colour = col$dark) +
  annotate("text", x = 5.0, y = 9.25, label = "UC evidence", family = "Arial",
           fontface = "bold", size = 3.2, colour = col$dark) +
  annotate("text", x = 8.3, y = 9.25, label = "Analysis gates", family = "Arial",
           fontface = "bold", size = 3.2, colour = col$dark) +
  annotate("label", x = 1.7, y = 0.55,
           label = "Preselected panel; no\nwhole-transcriptome inference",
           family = "Arial", size = 2.1, fill = "#F7F7F7",
           colour = "#666666", label.padding = grid::unit(0.14, "lines")) +
  annotate("label", x = 5.0, y = -0.42,
           label = "Primary question: EoC vs active colonic UC.\nGlobal convergence and local overlap are separate questions; EoG remains exploratory.",
           family = "Arial", size = 2.35, fill = "#FFF9E8",
           colour = col$dark, label.padding = grid::unit(0.18, "lines")) +
  coord_cartesian(xlim = c(0, 10), ylim = c(-0.75, 9.65), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(6, 6, 6, 6))

save_figure("Figure1_study_design_revised", p1, 7.2, 6.4)

# -----------------------------------------------------------------------------
# Figure 2: UC reference reproducibility
# -----------------------------------------------------------------------------
uc <- read.csv(file.path(res_dir, "UC_reproducibility_gene_level.csv"), check.names = FALSE)
uc$status <- "Not discovery-significant"
sig <- uc$q_discovery < 0.05
same <- sign(uc$logFC_discovery) == sign(uc$logFC_validation)
uc$status[sig & same] <- "Discovery-significant, concordant"
uc$status[sig & !same] <- "Discovery-significant, discordant"
uc$status <- factor(uc$status, levels = c(
  "Not discovery-significant",
  "Discovery-significant, concordant",
  "Discovery-significant, discordant"
))
lab_genes <- c("DUOX2", "S100A8", "S100A9", "LCN2", "IL1B", "CCL18", "GUCA2A", "BEST4", "GPA33")
ucl <- uc[uc$gene %in% lab_genes, ]
lim2 <- max(abs(c(uc$t_discovery, uc$t_validation)), na.rm = TRUE) * 1.04

p2 <- ggplot(uc, aes(t_discovery, t_validation)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#8E8E8E") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#8E8E8E") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = "dotted", colour = "#666666") +
  geom_point(aes(colour = status), size = 0.75, alpha = 0.42) +
  geom_text_repel(data = ucl, aes(label = gene), family = "Arial", size = 2.2,
                  min.segment.length = 0, seed = 20260827, box.padding = 0.2,
                  max.overlaps = Inf, colour = col$dark, segment.size = 0.25) +
  annotate("label", x = -0.96 * lim2, y = 0.93 * lim2,
           hjust = 0, vjust = 1,
           label = "17,730 common genes\nSpearman rho = 0.769\n8,089 / 8,956 (90.3%) discovery-significant genes\nhad the same direction in validation",
           family = "Arial", size = 2.55, lineheight = 1.05, fill = "white",
           label.size = 0.3, label.padding = grid::unit(0.14, "lines")) +
  scale_colour_manual(values = c(
    "Not discovery-significant" = "#BDBDBD",
    "Discovery-significant, concordant" = col$navy,
    "Discovery-significant, discordant" = col$vermilion
  )) +
  coord_cartesian(xlim = c(-lim2, lim2), ylim = c(-lim2, lim2)) +
  labs(
    x = "Moderated t statistic: GSE75214 discovery",
    y = "Moderated t statistic: GSE16879 validation",
    colour = NULL
  ) +
  theme_pub(8.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

save_figure("Figure2_UC_reproducibility_revised", p2, 7.2, 6.2)

# -----------------------------------------------------------------------------
# Figure 3: primary EoC convergence and sensitivity
# -----------------------------------------------------------------------------
null_disc <- read.csv(file.path(res_dir, "EoC_UC_discovery_exact_null.csv"))
null_val <- read.csv(file.path(res_dir, "EoC_UC_validation_perm_null.csv"))
primary <- read.csv(file.path(res_dir, "primary_EoC_UC_summary.csv"), check.names = FALSE)
loo <- read.csv(file.path(res_dir, "EoC_leave_one_out_primary_rho.csv"), check.names = FALSE)
obs_d <- primary$rho_moderated_t[1]
obs_v <- primary$rho_moderated_t[2]
p_d <- primary$empirical_p_two_sided[1]
p_v <- primary$empirical_p_two_sided[2]

p3a <- ggplot(null_disc, aes(rho)) +
  geom_histogram(bins = 46, fill = "#B9BDC3", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = obs_d, colour = col$vermilion, linewidth = 0.9) +
  annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
           label = sprintf("observed rho = %.3f\nexact two-sided P = %.3f", obs_d, p_d),
           family = "Arial", size = 2.5, fill = "white", label.size = 0.25) +
  labs(title = "Discovery: GSE75214", x = "Permuted Spearman rho", y = "Exact assignments") +
  theme_pub(8)

p3b <- ggplot(null_val, aes(rho)) +
  geom_histogram(bins = 46, fill = "#B9BDC3", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = obs_v, colour = col$vermilion, linewidth = 0.9) +
  annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
           label = sprintf("observed rho = %.3f\ntwo-sided P = %.3f", obs_v, p_v),
           family = "Arial", size = 2.5, fill = "white", label.size = 0.25) +
  labs(title = "Locked validation: GSE16879", x = "Permuted Spearman rho", y = "Permutations") +
  theme_pub(8)

loo_long <- rbind(
  data.frame(sample = loo$omitted_sample, group = loo$omitted_group, cohort = "Discovery", rho = loo$rho_discovery),
  data.frame(sample = loo$omitted_sample, group = loo$omitted_group, cohort = "Validation", rho = loo$rho_validation)
)
loo_long$sample <- factor(loo_long$sample, levels = rev(unique(loo$omitted_sample)))
p3c <- ggplot(loo_long, aes(rho, sample, colour = cohort, shape = group)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_point(size = 2.1, alpha = 0.9, position = position_dodge(width = 0.35)) +
  geom_vline(xintercept = obs_d, linetype = "dashed", colour = col$navy, linewidth = 0.5) +
  geom_vline(xintercept = obs_v, linetype = "dashed", colour = col$orange, linewidth = 0.5) +
  scale_colour_manual(values = c("Discovery" = col$navy, "Validation" = col$orange)) +
  scale_shape_manual(
    values = c("Non_EoC" = 16, "Active_EoC" = 17),
    labels = c("Non_EoC" = "Non-EoC", "Active_EoC" = "Active EoC")
  ) +
  labs(title = "Leave-one-out stability of the primary rank correlation",
       x = "Spearman rho after omitting one EoC sample", y = NULL, colour = NULL, shape = "Omitted group") +
  theme_pub(8) +
  theme(legend.position = "bottom")

p3 <- (p3a | p3b) / p3c +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 11)
    )
  )

save_figure("Figure3_EoC_convergence_revised", p3, 7.2, 7.0)

# -----------------------------------------------------------------------------
# Figure 4: all Hallmark evidence for the primary EoC axis
# -----------------------------------------------------------------------------
hm_primary <- read.csv(file.path(res_dir, "Hallmark_primary_discovery_validation.csv"), check.names = FALSE)
hm_long <- rbind(
  data.frame(pathway = hm_primary$pathway, cohort = "EoC", direction = hm_primary$Direction_EoC, fdr = hm_primary$FDR_EoC),
  data.frame(pathway = hm_primary$pathway, cohort = "UC discovery", direction = hm_primary$Direction_UC_discovery, fdr = hm_primary$FDR_UC_discovery),
  data.frame(pathway = hm_primary$pathway, cohort = "UC validation", direction = hm_primary$Direction_UC_validation, fdr = hm_primary$FDR_UC_validation)
)
hm_long$score <- signed_fdr(hm_long$direction, hm_long$fdr, cap = 15)
hm_long$label <- clean_hm(hm_long$pathway)
eoc_order <- hm_primary$pathway[order(hm_primary$FDR_EoC, hm_primary$pathway)]
hm_long$label <- factor(hm_long$label, levels = rev(clean_hm(eoc_order)))
hm_long$cohort <- factor(hm_long$cohort, levels = c("EoC", "UC discovery", "UC validation"),
                         labels = c("EoC", "UC\ndiscovery", "UC\nvalidation"))

p4 <- ggplot(hm_long, aes(cohort, label, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = hm_long[hm_long$fdr < 0.05, ], shape = 21, size = 1.2,
             fill = "white", colour = "#202020", stroke = 0.3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    limits = c(-15, 15), oob = squish,
    name = expression("signed " * -log[10] * "(FDR)")
  ) +
  labs(x = NULL, y = NULL,
       subtitle = "All 50 Hallmark sets are shown; circles indicate within-cohort FDR < 0.05.") +
  theme_pub(8) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6.4),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right"
  )

save_figure("Figure4_primary_Hallmark_revised", p4, 7.2, 9.4)

# -----------------------------------------------------------------------------
# Figure 5: post hoc EoG branch
# -----------------------------------------------------------------------------
null_eog <- read.csv(file.path(res_dir, "EoG_UC_exploratory_perm_null.csv"))
eog_summary <- read.csv(file.path(res_dir, "exploratory_EoG_UC_summary.csv"), check.names = FALSE)
obs_eog <- eog_summary$rho_moderated_t[1]
p_eog <- eog_summary$empirical_p_two_sided[1]

p5a <- ggplot(null_eog, aes(rho)) +
  geom_histogram(bins = 48, fill = "#B9BDC3", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = obs_eog, colour = col$vermilion, linewidth = 0.9) +
  annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
           label = sprintf("observed rho = %.3f\ntwo-sided P = %.3f", obs_eog, p_eog),
           family = "Arial", size = 2.4, fill = "white", label.size = 0.25) +
  labs(title = "Genome-wide rank correlation", x = "Permuted Spearman rho", y = "Permutations") +
  theme_pub(8)

hm_eog <- read.csv(file.path(res_dir, "Hallmark_exploratory_EoG_UC.csv"), check.names = FALSE)
hm_eog <- hm_eog[hm_eog$replicated, ]
hm_eog_long <- rbind(
  data.frame(pathway = hm_eog$pathway, cohort = "EoG", direction = hm_eog$Direction_EoG, fdr = hm_eog$FDR_EoG),
  data.frame(pathway = hm_eog$pathway, cohort = "UC discovery", direction = hm_eog$Direction_UC_discovery, fdr = hm_eog$FDR_UC_discovery),
  data.frame(pathway = hm_eog$pathway, cohort = "UC validation", direction = hm_eog$Direction_UC_validation, fdr = hm_eog$FDR_UC_validation)
)
hm_eog_long$score <- signed_fdr(hm_eog_long$direction, hm_eog_long$fdr, cap = 15)
hm_eog_long$label <- factor(clean_hm(hm_eog_long$pathway), levels = rev(clean_hm(hm_eog$pathway)))
hm_eog_long$cohort <- factor(hm_eog_long$cohort, levels = c("EoG", "UC discovery", "UC validation"),
                             labels = c("EoG", "UC\ndiscovery", "UC\nvalidation"))

p5b <- ggplot(hm_eog_long, aes(cohort, label, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
                       limits = c(-15, 15), oob = squish, name = "Signed -log10(FDR)") +
  guides(fill = guide_colourbar(direction="horizontal", title.position="top", barwidth = grid::unit(32, "mm"), barheight = grid::unit(2.5, "mm"))) +
  labs(title = "Exploratory Hallmark pathways", x = NULL, y = NULL) +
  theme_pub(7.6) +
  theme(legend.position="bottom", panel.grid = element_blank(), axis.text.y = element_text(size = 6.8), axis.text.x = element_text(size=6.8,face = "bold")) +
  scale_x_discrete(labels=c("EoG","UC\nDisc.","UC\nVal."))

v1 <- read.csv(file.path(res_dir, "EoG_signature_GSE16879_sample_validation.csv"))
v2 <- read.csv(file.path(res_dir, "EoG_signature_GSE36807_conditional_sensitivity.csv"))
vf <- data.frame(
  cohort = factor(c("GSE16879\nwithin-center", "GSE36807\ncross-center sensitivity"),
                  levels = rev(c("GSE16879\nwithin-center", "GSE36807\ncross-center sensitivity"))),
  g = c(v1$hedges_g, v2$hedges_g),
  low = c(v1$ci_low, v2$ci_low),
  high = c(v1$ci_high, v2$ci_high),
  p = c(v1$permutation_p_two_sided, v2$permutation_p_two_sided)
)
p5c <- ggplot(vf, aes(g, cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "#777777") +
  geom_segment(aes(x = low, xend = high, yend = cohort), linewidth = 0.75, colour = col$navy) +
  geom_point(size = 2.6, colour = col$navy) +
  geom_text(aes(x = high + 0.25, label = sprintf("P = %.4f", p)), hjust = 0,
            family = "Arial", size = 2.35, colour = col$dark) +
  coord_cartesian(xlim = c(0, max(vf$high) + 1.3), clip = "off") +
  labs(title = "209-gene score validation", x = "Hedges g (95% CI)", y = NULL) +
  theme_pub(8) +
  theme(plot.margin = margin(5, 20, 5, 5))

ora_h <- read.csv(file.path(res_dir, "ORA_EoG_UC_Hallmark_full.csv"), check.names = FALSE)
ora_k <- read.csv(file.path(res_dir, "ORA_EoG_UC_KEGG_full.csv"), check.names = FALSE)
ora_h$collection <- "Hallmark"; ora_k$collection <- "KEGG"
ora <- rbind(ora_h[ora_h$q_value < 0.05, ], ora_k[ora_k$q_value < 0.05, ])
ora$label <- factor(clean_term(ora$term), levels = rev(clean_term(ora$term)))
p5d <- ggplot(ora, aes(-log10(q_value), label, colour = collection, size = overlap)) +
  geom_point(alpha = 0.9) +
  scale_y_discrete(labels=function(x) vapply(x,function(s) paste(strwrap(s,width=25),collapse="\n"),character(1))) +
  scale_x_continuous(expand=expansion(mult=c(.15,.20))) +
  scale_colour_manual(values = c("Hallmark" = col$vermilion, "KEGG" = col$blue)) +
  scale_size_continuous(range = c(2.2, 4.5), breaks = sort(unique(ora$overlap))) +
  labs(title = "209-gene ORA", x = expression(-log[10] * "(BH q)"),
       y = NULL, colour = NULL, size = NULL) +
  guides(size = "none") +
  theme_pub(7.7) +
  theme(legend.position = "bottom")

p5 <- (wrap_elements(full=p5a) | wrap_elements(full=p5b)) / (wrap_elements(full=p5c) | wrap_elements(full=p5d)) +
  plot_layout(heights=c(1.25,1)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 11)
    )
  )

save_figure("Figure5_EoG_exploratory_revised", p5, 7.2, 7.8)

# -----------------------------------------------------------------------------
# Figure 6: post hoc phenotype association
# -----------------------------------------------------------------------------
ph <- read.csv(file.path(res_dir, "EoG_signature_GSE109142_exploratory_scores.csv"), check.names = FALSE)
ph <- ph[complete.cases(ph[, c("score", "eos32", "histology", "age", "sex")]), ]
ph$histology_f <- factor(ph$histology, levels = 0:4)
ph$eos32_f <- factor(ph$eos32, levels = c("N", "Y"), labels = c("<32 eosinophils/hpf", "at least 32 eosinophils/hpf"))

p6a <- ggplot(ph, aes(histology_f, score)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, fill = col$pale_blue, colour = col$navy, linewidth = 0.55) +
  geom_point(position = position_jitter(width = 0.16, height = 0, seed = 20260827), size = 0.8, alpha = 0.48, colour = col$dark) +
  labs(title = "Histologic severity", x = "Severity score", y = "Directional signature score") +
  theme_pub(8)

p6b <- ggplot(ph, aes(eos32_f, score)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "#FFF2B2", colour = "#8A6D00", linewidth = 0.55) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 20260827), size = 0.8, alpha = 0.48, colour = col$dark) +
  labs(title = "Eosinophil grade", x = NULL, y = "Directional signature score") +
  theme_pub(8) +
  theme(axis.text.x = element_text(size = 6.8))

model <- read.csv(file.path(res_dir, "EoG_signature_GSE109142_exploratory_model.csv"), check.names = FALSE)
model <- model[model$term %in% c("histology", "eos32Y"), ]
model$label <- factor(ifelse(model$term == "histology", "Histologic severity\n(per 1-point increase)",
                             "Eosinophil grade\n(at least 32 vs <32/hpf)"),
                      levels = rev(c("Histologic severity\n(per 1-point increase)", "Eosinophil grade\n(at least 32 vs <32/hpf)")))
p6c <- ggplot(model, aes(estimate, label)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#777777") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.8, colour = col$purple) +
  geom_point(size = 2.8, colour = col$purple) +
  geom_text(aes(x = ci_high + 0.035, label = sprintf("P = %.3g", p_value)), hjust = 0,
            family = "Arial", size = 2.4, colour = col$dark) +
  coord_cartesian(xlim = c(min(model$ci_low) - 0.08, max(model$ci_high) + 0.28), clip = "off") +
  labs(title = "HC3-robust adjusted model (n = 198)",
       subtitle = "Mutually adjusted for both phenotypes, age and sex", x = "Regression coefficient (95% CI)", y = NULL) +
  theme_pub(8) +
  theme(plot.margin = margin(5, 24, 5, 5))

p6 <- (p6a | p6b) / p6c +
  plot_layout(heights = c(1.15, 0.9)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 11)
    )
  )

save_figure("Figure6_phenotype_exploratory_revised", p6, 7.2, 7.0)

# -----------------------------------------------------------------------------
# Supplementary Figure S1: corrected signature construction
# -----------------------------------------------------------------------------
sg <- read.csv(file.path(res_dir, "EoG_UC_exploratory_signature_gate_full.csv"), check.names = FALSE)
sg$status <- "Not selected"
sg$status[sg$selected & sg$logFC_a > 0] <- "Selected, concordant up"
sg$status[sg$selected & sg$logFC_a < 0] <- "Selected, concordant down"
sg$status <- factor(sg$status, levels = c("Not selected", "Selected, concordant up", "Selected, concordant down"))
sel <- sg[sg$selected, ]
top_lab <- sel[order(sel$conjunction_q), ][1:min(14, nrow(sel)), ]
lims <- max(abs(c(sg$t_a, sg$t_b)), na.rm = TRUE) * 1.04
pS1 <- ggplot(sg, aes(t_a, t_b)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#888888") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#888888") +
  geom_point(aes(colour = status), size = 0.65, alpha = 0.4) +
  geom_text_repel(data = top_lab, aes(label = gene), family = "Arial", size = 2.25,
                  seed = 20260827, max.overlaps = Inf, min.segment.length = 0,
                  segment.size = 0.25, colour = col$dark) +
  annotate("label", x = -0.95 * lims, y = 0.94 * lims, hjust = 0, vjust = 1,
           label = "209 same-direction genes\n43 up / 166 down\nconjunction BH q < 0.05",
           family = "Arial", size = 2.5, fill = "white", label.size = 0.25) +
  scale_colour_manual(values = c("Not selected" = "#BDBDBD", "Selected, concordant up" = col$vermilion,
                                 "Selected, concordant down" = col$navy)) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims)) +
  labs(x = "Moderated t statistic: EoG", y = "Moderated t statistic: UC discovery", colour = NULL) +
  theme_pub(8.5) +
  theme(legend.position = "bottom")
save_figure("Supplementary_Figure_S1_EoG_signature_gate", pS1, 7.2, 6.2)

cat("Revised figures written to:", fig_dir, "\n")

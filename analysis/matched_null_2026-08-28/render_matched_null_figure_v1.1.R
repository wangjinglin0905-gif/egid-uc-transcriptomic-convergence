suppressPackageStartupMessages(library(ggplot2))

root <- normalizePath(Sys.getenv("EOGE_IBD_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
base_dir <- file.path(root, "analysis", "matched_null_2026-08-28")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(root,"generated_figures","matched_null")
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)

null_stats_all <- read.csv(file.path(result_dir, "null_endpoint_statistics.csv"), check.names = FALSE)
competitive_results <- read.csv(file.path(result_dir, "competitive_test_results.csv"), check.names = FALSE)

panel_levels <- c(
  "a  Technical-property null: GSE36807 Hedges g",
  "b  Technical-property null: GSE109142 histology HC3 t",
  "c  UC-effect-matched null: GSE36807 Hedges g",
  "d  UC-effect-matched null: GSE109142 histology HC3 t"
)

assign_panel <- function(model, endpoint) {
  ifelse(
    model == "technical" & endpoint == "GSE36807_Hedges_g", panel_levels[1],
    ifelse(model == "technical", panel_levels[2],
           ifelse(endpoint == "GSE36807_Hedges_g", panel_levels[3], panel_levels[4]))
  )
}

format_p <- function(x) {
  ifelse(x < 0.001, sprintf("%.5f", x), sprintf("%.3f", x))
}

format_percentile <- function(x) {
  ifelse(x >= 0.99995, "observed percentile > 99.99", sprintf("observed percentile = %.2f", 100 * x))
}

plot_data <- null_stats_all
plot_data$panel <- factor(assign_panel(plot_data$null_model, plot_data$endpoint), levels = panel_levels)
write.csv(plot_data[, c("panel", "null_model", "endpoint", "replicate", "statistic")],
          file.path(result_dir, "FigureS_matched_null_plotted_data.csv"), row.names = FALSE)

annotation <- competitive_results
annotation$panel <- factor(assign_panel(annotation$null_model, annotation$endpoint), levels = panel_levels)
annotation$label <- sprintf(
  "Observed = %.2f\ncompetitive P = %s\nHolm P = %s\n%s",
  annotation$observed,
  format_p(annotation$p_one_sided),
  format_p(annotation$holm_p_one_sided),
  format_percentile(annotation$observed_percentile)
)
annotation$x_position <- ifelse(annotation$endpoint == "GSE36807_Hedges_g", -Inf, Inf)
annotation$hjust_value <- ifelse(annotation$endpoint == "GSE36807_Hedges_g", -0.05, 1.05)

p <- ggplot(plot_data, aes(x = statistic)) +
  geom_histogram(bins = 55, fill = "#9BA3AA", color = "white", linewidth = 0.15) +
  geom_vline(data = annotation, aes(xintercept = observed), inherit.aes = FALSE,
             color = "#C43C39", linewidth = 0.9, linetype = "dashed") +
  geom_text(data = annotation,
            aes(x = x_position, y = Inf, label = label, hjust = hjust_value),
            inherit.aes = FALSE, vjust = 1.12, size = 3.05, lineheight = 1.04,
            color = "#8F2323") +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  labs(
    title = "Competitive evaluation of the 209-gene EoG-UC signature",
    subtitle = "Observed statistics (red) versus 10,000 matched random gene sets per null model",
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

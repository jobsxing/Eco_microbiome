#!/usr/bin/env Rscript

# R2Q20: revised ARG-class association analysis.
# Primary scale: ARGs-OAP 16S-normalized copies per cell.
# Sensitivity scales: RPKM and total-sum-scaled RPKM, all log2 transformed.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

set.seed(20260811)
argv <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/")
out_dir <- dirname(script_file)
project_dir <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
project_override <- Sys.getenv("OHMD_PROJECT_ROOT", unset = "")
if (nzchar(project_override)) project_dir <- normalizePath(project_override, winslash = "/")
metadata_file <- file.path(project_dir, "data", "metadata_final.csv")
arg_dir <- file.path(project_dir, "Figure", "Figure 3")
rpkm_file <- file.path(arg_dir, "merged_rpkm_type.tsv")
copies_file <- file.path(arg_dir, "merged_normalized16S_type.tsv")
stopifnot(all(file.exists(c(metadata_file, rpkm_file, copies_file))))

min_primary_prevalence <- 0.10

make_repeated_unit <- function(core, source) {
  tail_id <- str_match(core, "_(\\d+)$")[, 2]
  farmer_inner <- str_match(core, "(?:FFE|FFNV)(\\d+)")[, 2]
  farmer_id <- recode(
    paste(farmer_inner, tail_id, sep = "_"),
    "1_1" = "1", "1_2" = "2", "2_1" = "3", "2_2" = "4",
    .default = NA_character_
  )
  case_when(
    source %in% c("Farmer Excrement", "Farmer Nasal Vestibule") ~ paste0("Farmer_", farmer_id),
    source %in% c("Non-farmer gut", "Non-farmer nasal") ~ paste0("NonFarmer_", tail_id),
    source == "Goose Excrement" ~ paste0("Goose_", tail_id),
    source == "Ostrich Excrement" ~ paste0("Ostrich_", tail_id),
    source == "Peacock Excrement" ~ paste0("Peacock_", tail_id),
    source == "Goose Farm Soil" ~ paste0("GooseSoil_site", tail_id),
    source == "Ostrich Farm Soil" ~ paste0("OstrichSoil_site", tail_id),
    source == "Goose Paddling Pool" ~ paste0("Pond_site", tail_id),
    source == "Surrounding Rivers" ~ paste0("River_site", tail_id),
    TRUE ~ NA_character_
  )
}

meta0 <- read_csv(metadata_file, show_col_types = FALSE, name_repair = "minimal")
if (names(meta0)[1] == "" || str_detect(names(meta0)[1], "^\\.\\.\\.")) meta0 <- meta0[, -1]
source_map <- c(
  "Non-farmer gut" = "Non-farmer gut", "Farmer Excrement" = "Farmer gut",
  "Goose Excrement" = "Goose gut", "Ostrich Excrement" = "Ostrich gut",
  "Peacock Excrement" = "Peacock gut"
)
meta <- meta0 %>%
  filter(`Sample Source` %in% names(source_map)) %>%
  transmute(
    Sample_ID = Sample_core_for_mpa,
    Group = relevel(factor(unname(source_map[`Sample Source`])), ref = "Non-farmer gut"),
    Batch = factor(paste0("T", as.integer(`Sampling Time`)), levels = paste0("T", 1:14)),
    Cluster = factor(make_repeated_unit(Sample_core, `Sample Source`))
  )
if (anyDuplicated(meta$Sample_ID) || anyNA(meta$Cluster)) stop("Metadata validation failed.")

read_arg_matrix <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  names(x)[1] <- "Feature"
  x <- x %>% group_by(Feature) %>%
    summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  matrix <- x %>% column_to_rownames("Feature") %>% as.matrix()
  storage.mode(matrix) <- "double"
  matrix[!is.finite(matrix) | matrix < 0] <- 0
  t(matrix)[meta$Sample_ID, , drop = FALSE]
}

rpkm <- read_arg_matrix(rpkm_file)
copies_per_cell <- read_arg_matrix(copies_file)
stopifnot(identical(dim(rpkm), dim(copies_per_cell)), identical(colnames(rpkm), colnames(copies_per_cell)))
rpkm_tss <- rpkm / rowSums(rpkm)
rpkm_tss[!is.finite(rpkm_tss)] <- 0

# Define one testing family from the primary copies-per-cell scale and use the
# identical ARG classes for both sensitivity scales.
primary_prevalence <- colMeans(copies_per_cell > 0)
keep_features <- primary_prevalence >= min_primary_prevalence &
  colMeans(copies_per_cell) > 0
if (!any(keep_features)) stop("No ARG classes passed the primary prevalence filter.")
rpkm <- rpkm[, keep_features, drop = FALSE]
copies_per_cell <- copies_per_cell[, keep_features, drop = FALSE]
rpkm_tss <- rpkm_tss[, keep_features, drop = FALSE]

log_transform <- function(x) {
  minimum_positive <- apply(x, 2, function(z) min(z[z > 0], na.rm = TRUE))
  if (any(!is.finite(minimum_positive))) stop("A retained feature has no positive values.")
  y <- x
  for (j in seq_len(ncol(y))) y[y[, j] == 0, j] <- minimum_positive[j] / 2
  log2(y)
}

fit_matrix <- function(raw_matrix, scale_name) {
  prevalence <- colMeans(raw_matrix > 0)
  y <- log_transform(raw_matrix)
  design <- model.matrix(~ Group + Batch, data = meta)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient ARG model.")
  inverse_xtx <- solve(crossprod(design))
  coefficients <- inverse_xtx %*% crossprod(design, y)
  residuals <- y - design %*% coefficients
  leverage <- rowSums((design %*% inverse_xtx) * design)
  residuals_hc3 <- residuals / pmax(1 - leverage, 1e-8)
  clusters <- split(seq_len(nrow(meta)), meta$Cluster)
  cluster_n <- length(clusters)
  coefficient_rows <- grep("^Group", rownames(coefficients))
  comparison_levels <- levels(meta$Group)[-1]

  bind_rows(lapply(seq_along(coefficient_rows), function(k) {
    coefficient_row <- coefficient_rows[k]
    contrast_a <- inverse_xtx[coefficient_row, , drop = FALSE]
    variance <- rep(0, ncol(y))
    for (index in clusters) {
      score <- crossprod(design[index, , drop = FALSE], residuals_hc3[index, , drop = FALSE])
      projected <- as.numeric(contrast_a %*% score)
      variance <- variance + projected^2
    }
    variance <- variance * cluster_n / (cluster_n - 1)
    estimate <- coefficients[coefficient_row, ]
    standard_error <- sqrt(pmax(variance, 0))
    statistic <- estimate / standard_error
    p_value <- 2 * pt(-abs(statistic), df = cluster_n - 1)
    tibble(
      Scale = scale_name, Reference = "Non-farmer gut",
      Comparison = comparison_levels[k], ARG_class = colnames(y),
      Estimate = as.numeric(estimate), SE = standard_error,
      Df = cluster_n - 1, Statistic = statistic, P = p_value,
      Prevalence = as.numeric(prevalence[colnames(y)]),
      N_samples = nrow(meta), N_clusters = cluster_n
    )
  })) %>% mutate(Q = p.adjust(P, method = "BH"))
}

all_results <- bind_rows(
  fit_matrix(copies_per_cell, "Copies per cell"),
  fit_matrix(rpkm, "RPKM"),
  fit_matrix(rpkm_tss, "TSS-scaled RPKM")
)
primary_results <- filter(all_results, Scale == "Copies per cell") %>% arrange(Q, desc(abs(Statistic)))
write_csv(primary_results, file.path(out_dir, "sup_table_R2Q20_revised_copies_per_cell_associations.csv"))

primary <- primary_results %>%
  select(Reference, Comparison, ARG_class, Estimate_primary = Estimate,
         Statistic_primary = Statistic, Q_primary = Q)
paired <- filter(all_results, Scale != "Copies per cell") %>%
  left_join(primary, by = c("Reference", "Comparison", "ARG_class"))
sensitivity_summary <- paired %>% group_by(Scale) %>% summarise(
  N_associations = n(),
  Spearman_rho_statistic = cor(Statistic_primary, Statistic, method = "spearman"),
  Direction_concordance = mean(sign(Estimate_primary) == sign(Estimate)),
  Significant_primary = sum(Q_primary <= 0.05),
  Significant_comparator = sum(Q <= 0.05),
  Significant_both = sum(Q_primary <= 0.05 & Q <= 0.05),
  Primary_only = sum(Q_primary <= 0.05 & Q > 0.05),
  Comparator_only = sum(Q_primary > 0.05 & Q <= 0.05),
  .groups = "drop"
)
write_csv(sensitivity_summary, file.path(out_dir, "sup_table_R2Q20_scale_sensitivity.csv"))

theme_publication <- theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.tag = element_text(face = "bold", size = 14), axis.text = element_text(color = "black"))
plot_scatter <- function(scale_name, x_label) {
  dat <- filter(paired, Scale == scale_name)
  rho <- cor(dat$Statistic_primary, dat$Statistic, method = "spearman")
  ggplot(dat, aes(Statistic, Statistic_primary)) +
    geom_point(aes(color = Q_primary <= 0.05 | Q <= 0.05), size = 1.7, alpha = 0.75) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_color_manual(values = c(`TRUE` = "#9f2d2d", `FALSE` = "#a3a3a3"), guide = "none") +
    labs(x = x_label, y = "Copies-per-cell model t statistic",
         subtitle = sprintf("Spearman rho = %.3f", rho)) +
    coord_equal() + theme_publication
}

count_data <- all_results %>% group_by(Scale) %>% summarise(Significant = sum(Q <= 0.05), .groups = "drop") %>%
  mutate(Scale = factor(Scale, levels = c("Copies per cell", "RPKM", "TSS-scaled RPKM")))
panel_c <- ggplot(count_data, aes(Scale, Significant, fill = Scale)) +
  geom_col(width = 0.7) + geom_text(aes(label = Significant), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Copies per cell" = "#355d62", "RPKM" = "#3d5380", "TSS-scaled RPKM" = "#b63233"), guide = "none") +
  labs(x = NULL, y = "Significant ARG-class associations (q <= 0.05)") +
  theme_publication + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  expand_limits(y = max(count_data$Significant) * 1.12 + 1)

top_primary <- primary_results %>% filter(Q <= 0.05) %>% slice_head(n = 14) %>%
  mutate(
    Label = paste0(ARG_class, " - ", Comparison),
    Label = factor(Label, levels = rev(Label)),
    Lower = Estimate - qt(0.975, Df) * SE,
    Upper = Estimate + qt(0.975, Df) * SE
  )
panel_d <- ggplot(top_primary, aes(Estimate, Label, color = Comparison)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0) +
  geom_point(size = 2) +
  labs(x = "Log2 copies-per-cell difference", y = NULL, color = NULL) +
  theme_publication + theme(legend.position = "top", axis.text.y = element_text(size = 8))

supplementary_figure <-
  (plot_scatter("RPKM", "RPKM model t statistic") +
     plot_scatter("TSS-scaled RPKM", "TSS-scaled RPKM model t statistic")) /
  (panel_c + panel_d) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "sup_fig_R2Q20_arg_scale_sensitivity.pdf"), supplementary_figure,
       width = 12, height = 9, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "sup_fig_R2Q20_arg_scale_sensitivity.png"), supplementary_figure,
       width = 12, height = 9, units = "in", dpi = 300, bg = "white")

message("R2Q20 ARG-scale analysis completed")

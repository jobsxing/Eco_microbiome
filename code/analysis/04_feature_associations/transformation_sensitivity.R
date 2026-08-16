#!/usr/bin/env Rscript

# Transformation sensitivity for taxonomic and functional associations.
# The log2 feature-wise model is the revised primary analysis. The legacy
# arcsine-square-root branch is retained only to quantify sensitivity to the
# original transformation. Longitudinal dependence is handled using HC3
# cluster-robust covariance.

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
taxonomy_file <- file.path(project_dir, "data", "metaphlan4", "all_sample_taxonomy.tsv")
pathway_file <- file.path(project_dir, "data", "humanN3", "humann3_pathabundance_relab_unstratified.tsv")
stopifnot(all(file.exists(c(metadata_file, taxonomy_file, pathway_file))))

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
meta <- meta0 %>% transmute(
  Sample_ID = Sample_core_for_mpa,
  Sample_Source = `Sample Source`,
  Batch = factor(paste0("T", as.integer(`Sampling Time`)), levels = paste0("T", 1:14)),
  Cluster = factor(make_repeated_unit(Sample_core, `Sample Source`))
)
if (nrow(meta) != 500L || anyDuplicated(meta$Sample_ID) || anyNA(meta$Cluster) ||
    n_distinct(meta$Cluster) != 29L) stop("Metadata validation failed.")

message("Preparing terminal-genus relative-abundance matrix")
taxonomy <- read_tsv(
  taxonomy_file, show_col_types = FALSE, progress = interactive(),
  name_repair = "minimal", comment = "#"
)
genus <- taxonomy %>%
  filter(
    str_detect(ID, "(^|\\|)g__[^|]+$"),
    !str_detect(ID, regex("unclassified|unknown|uncultured|metagenome", ignore_case = TRUE))
  ) %>%
  mutate(Feature = str_remove(str_extract(ID, "g__[^|]+$"), "^g__")) %>%
  select(Feature, all_of(meta$Sample_ID)) %>%
  group_by(Feature) %>%
  summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
rm(taxonomy); invisible(gc())
genus_matrix <- genus %>% column_to_rownames("Feature") %>% as.matrix()
storage.mode(genus_matrix) <- "double"
genus_matrix[!is.finite(genus_matrix) | genus_matrix < 0] <- 0
genus_matrix <- t(genus_matrix)[meta$Sample_ID, , drop = FALSE]
genus_matrix <- genus_matrix / rowSums(genus_matrix)

message("Preparing HUMAnN3 pathway relative-abundance matrix")
pathway <- read_tsv(pathway_file, show_col_types = FALSE, name_repair = "minimal")
names(pathway)[1] <- "Feature"
pathway <- pathway %>%
  filter(!Feature %in% c("UNMAPPED", "UNINTEGRATED"), !str_detect(Feature, "\\|")) %>%
  select(Feature, all_of(meta$Sample_ID)) %>%
  group_by(Feature) %>%
  summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
pathway_matrix <- pathway %>% column_to_rownames("Feature") %>% as.matrix()
storage.mode(pathway_matrix) <- "double"
pathway_matrix[!is.finite(pathway_matrix) | pathway_matrix < 0] <- 0
pathway_matrix <- t(pathway_matrix)[meta$Sample_ID, , drop = FALSE]

maps <- list(
  Human_nasal_vs_bird = c(
    "Non-farmer nasal" = "Non-farmer nasal", "Farmer Nasal Vestibule" = "Farmer nasal",
    "Goose Excrement" = "Exotic bird", "Ostrich Excrement" = "Exotic bird",
    "Peacock Excrement" = "Exotic bird"
  ),
  Human_gut_vs_bird = c(
    "Non-farmer gut" = "Non-farmer gut", "Farmer Excrement" = "Farmer gut",
    "Goose Excrement" = "Exotic bird", "Ostrich Excrement" = "Exotic bird",
    "Peacock Excrement" = "Exotic bird"
  ),
  Goose_soil_context = c(
    "Goose Farm Soil" = "Goose-associated soil", "Ostrich Farm Soil" = "Ostrich-associated soil",
    "Goose Excrement" = "Goose gut"
  ),
  Ostrich_soil_context = c(
    "Goose Farm Soil" = "Goose-associated soil", "Ostrich Farm Soil" = "Ostrich-associated soil",
    "Ostrich Excrement" = "Ostrich gut"
  ),
  Goose_water_context = c(
    "Goose Paddling Pool" = "Goose paddling pool", "Surrounding Rivers" = "Surrounding river",
    "Goose Excrement" = "Goose gut"
  )
)
references <- c(
  Human_nasal_vs_bird = "Exotic bird", Human_gut_vs_bird = "Exotic bird",
  Goose_soil_context = "Goose gut", Ostrich_soil_context = "Ostrich gut",
  Goose_water_context = "Goose gut"
)

transform_matrix <- function(x, method) {
  if (method == "AST") return(asin(sqrt(x)))
  if (method == "LOG") {
    minimum_positive <- apply(x, 2, function(z) min(z[z > 0], na.rm = TRUE))
    if (any(!is.finite(minimum_positive))) stop("A retained feature has no positive values.")
    y <- x
    for (j in seq_len(ncol(y))) y[y[, j] == 0, j] <- minimum_positive[j] / 2
    return(log2(y))
  }
  stop("Unknown transformation")
}

fit_matrix_block <- function(raw_matrix, metadata, group_map, reference, block,
                             data_type, transformation, min_prevalence, min_mean) {
  keep_samples <- metadata$Sample_Source %in% names(group_map)
  md <- droplevels(metadata[keep_samples, , drop = FALSE])
  md$Group <- unname(group_map[md$Sample_Source])
  md$Group <- relevel(factor(md$Group), ref = reference)
  x <- raw_matrix[md$Sample_ID, , drop = FALSE]
  prevalence <- colMeans(x > 0)
  mean_abundance <- colMeans(x)
  keep_features <- prevalence >= min_prevalence & mean_abundance >= min_mean
  x <- x[, keep_features, drop = FALSE]
  prevalence <- prevalence[keep_features]
  mean_abundance <- mean_abundance[keep_features]
  if (!ncol(x)) stop("No features retained in ", block)

  y <- transform_matrix(x, transformation)
  design <- model.matrix(~ Group + Batch, data = md)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient model in ", block)
  inverse_xtx <- solve(crossprod(design))
  coefficients <- inverse_xtx %*% crossprod(design, y)
  residuals <- y - design %*% coefficients
  leverage <- rowSums((design %*% inverse_xtx) * design)
  residuals_hc3 <- residuals / pmax(1 - leverage, 1e-8)
  clusters <- split(seq_len(nrow(md)), md$Cluster)
  cluster_n <- length(clusters)
  coefficient_rows <- grep("^Group", rownames(coefficients))
  comparison_levels <- levels(md$Group)[-1]
  if (length(coefficient_rows) != length(comparison_levels)) stop("Group coefficient mapping failed.")

  results <- lapply(seq_along(coefficient_rows), function(k) {
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
      Data_Type = data_type, Block = block, Transformation = transformation,
      Reference = reference, Comparison = comparison_levels[k], Feature = colnames(y),
      Estimate = as.numeric(estimate), SE = standard_error,
      Df = cluster_n - 1, Statistic = statistic, P = p_value,
      Prevalence = as.numeric(prevalence[colnames(y)]),
      Mean_abundance = as.numeric(mean_abundance[colnames(y)]),
      N_samples = nrow(md), N_clusters = cluster_n
    )
  }) %>% bind_rows() %>% mutate(Q = p.adjust(P, method = "BH"))
  message(data_type, " / ", block, " / ", transformation, ": ",
          ncol(y), " features; ", sum(results$Q <= 0.05), " significant associations")
  results
}

run_definitions <- function(matrix, data_type, blocks, min_prevalence) {
  bind_rows(lapply(blocks, function(block) {
    bind_rows(lapply(c("LOG", "AST"), function(transformation) {
      fit_matrix_block(
        matrix, meta, maps[[block]], references[[block]], block, data_type,
        transformation, min_prevalence = min_prevalence, min_mean = 1e-5
      )
    }))
  }))
}

taxonomy_results <- run_definitions(
  genus_matrix, "Taxonomy", names(maps), min_prevalence = 0.10
)
functional_results <- run_definitions(
  pathway_matrix, "Functional pathway",
  c("Human_nasal_vs_bird", "Human_gut_vs_bird"), min_prevalence = 0.05
)
all_results <- bind_rows(taxonomy_results, functional_results)

# Keep the complete revised LOG-model output as an internal, reproducible
# figure source.  The supplementary table below remains restricted to
# significant associations for concise reporting.
write_csv(
  filter(all_results, Transformation == "LOG") %>%
    arrange(Data_Type, Block, Comparison, Q, desc(abs(Statistic))),
  file.path(out_dir, "analysis_R2Q19_all_log_associations.csv")
)

paired <- full_join(
  filter(all_results, Transformation == "LOG") %>%
    select(-Transformation) %>% rename_with(~ paste0(.x, "_LOG"),
      c(Estimate, SE, Df, Statistic, P, Q, Prevalence, Mean_abundance, N_samples, N_clusters)),
  filter(all_results, Transformation == "AST") %>%
    select(-Transformation) %>% rename_with(~ paste0(.x, "_AST"),
      c(Estimate, SE, Df, Statistic, P, Q, Prevalence, Mean_abundance, N_samples, N_clusters)),
  by = c("Data_Type", "Block", "Reference", "Comparison", "Feature")
)

summary_by_block <- paired %>% group_by(Data_Type, Block) %>% summarise(
  N_associations = n(),
  Spearman_rho_statistic = cor(Statistic_LOG, Statistic_AST, method = "spearman", use = "complete.obs"),
  Direction_concordance = mean(sign(Estimate_LOG) == sign(Estimate_AST), na.rm = TRUE),
  Significant_LOG = sum(Q_LOG <= 0.05, na.rm = TRUE),
  Significant_AST = sum(Q_AST <= 0.05, na.rm = TRUE),
  Significant_both = sum(Q_LOG <= 0.05 & Q_AST <= 0.05, na.rm = TRUE),
  LOG_only = sum(Q_LOG <= 0.05 & Q_AST > 0.05, na.rm = TRUE),
  AST_only = sum(Q_LOG > 0.05 & Q_AST <= 0.05, na.rm = TRUE),
  .groups = "drop"
)
write_csv(summary_by_block, file.path(out_dir, "sup_table_R2Q19_transformation_sensitivity.csv"))

revised_significant <- filter(all_results, Transformation == "LOG", Q <= 0.05) %>%
  arrange(Data_Type, Block, Q, desc(abs(Statistic)))
write_csv(revised_significant, file.path(out_dir, "sup_table_R2Q19_revised_log_associations.csv"))

overall_stats <- paired %>% group_by(Data_Type) %>% summarise(
  Rho = cor(Statistic_LOG, Statistic_AST, method = "spearman", use = "complete.obs"),
  .groups = "drop"
)
plot_scatter <- function(data_type, panel_label) {
  dat <- filter(paired, Data_Type == data_type)
  rho <- filter(overall_stats, Data_Type == data_type)$Rho
  ggplot(dat, aes(Statistic_AST, Statistic_LOG)) +
    stat_bin_2d(bins = 55) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_fill_viridis_c(trans = "log10", name = "Associations") +
    labs(
      x = "Arcsine model t statistic", y = "Log model t statistic",
      subtitle = sprintf("Spearman rho = %.3f", rho)
    ) + coord_equal() + theme_classic(base_family = "Arial", base_size = 11) +
    theme(plot.tag = element_text(face = "bold", size = 14))
}

count_data <- paired %>% mutate(
  Category = case_when(
    Q_LOG <= 0.05 & Q_AST <= 0.05 ~ "Both",
    Q_LOG <= 0.05 & Q_AST > 0.05 ~ "Log only",
    Q_LOG > 0.05 & Q_AST <= 0.05 ~ "Arcsine only",
    TRUE ~ "Neither"
  )
) %>% filter(Category != "Neither") %>% count(Data_Type, Category) %>%
  complete(Data_Type, Category = c("Both", "Log only", "Arcsine only"), fill = list(n = 0))
plot_counts <- function(data_type) {
  ggplot(filter(count_data, Data_Type == data_type), aes(Category, n, fill = Category)) +
    geom_col(width = 0.7) + geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Both" = "#355d62", "Log only" = "#3d5380", "Arcsine only" = "#b63233"), guide = "none") +
    labs(x = NULL, y = "Significant associations (q <= 0.05)") +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.tag = element_text(face = "bold", size = 14)) +
    expand_limits(y = max(filter(count_data, Data_Type == data_type)$n) * 1.12 + 1)
}

supplementary_figure <-
  (plot_scatter("Taxonomy") + plot_counts("Taxonomy")) /
  (plot_scatter("Functional pathway") + plot_counts("Functional pathway")) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "sup_fig_R2Q19_transformation_sensitivity.pdf"), supplementary_figure,
       width = 11, height = 9, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "sup_fig_R2Q19_transformation_sensitivity.png"), supplementary_figure,
       width = 11, height = 9, units = "in", dpi = 300, bg = "white")

message("R2Q19 transformation sensitivity analysis completed")

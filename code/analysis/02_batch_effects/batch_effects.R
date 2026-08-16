#!/usr/bin/env Rscript

# R1Q4: time-associated extraction/sequencing batch assessment.
# Sampling Time 1-14 are treated as 14 independent extraction/sequencing
# batches. Batch and biological time are therefore perfectly confounded.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lme4); library(patchwork)
  library(readr); library(stringr); library(tibble); library(vegan)
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
reads_file <- file.path(project_dir, "data", "metaphlan4", "metaphlan4_reads_summary.csv")
stopifnot(all(file.exists(c(metadata_file, taxonomy_file, reads_file))))

batch_levels <- paste0("T", 1:14)
facet_levels <- c("Human gut", "Human nasal", "Exotic bird", "Environment")
source_levels <- c(
  "Farmer Excrement", "Non-farmer gut", "Farmer Nasal Vestibule",
  "Non-farmer nasal", "Peacock Excrement", "Goose Excrement",
  "Ostrich Excrement", "Goose Farm Soil", "Ostrich Farm Soil",
  "Goose Paddling Pool", "Surrounding Rivers"
)
min_relative_abundance <- 1e-5
min_prevalence <- 0.10
permutation_n <- 999L

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
  Sample_core_for_mpa,
  Sample_Source = factor(`Sample Source`, levels = source_levels),
  Sampling_Time = as.integer(`Sampling Time`),
  Batch = factor(paste0("T", Sampling_Time), levels = batch_levels),
  Facet_Group = factor(case_when(
    `Sample Source` %in% c("Farmer Excrement", "Non-farmer gut") ~ "Human gut",
    `Sample Source` %in% c("Farmer Nasal Vestibule", "Non-farmer nasal") ~ "Human nasal",
    `Sample Source` %in% c("Goose Excrement", "Ostrich Excrement", "Peacock Excrement") ~ "Exotic bird",
    TRUE ~ "Environment"
  ), levels = facet_levels),
  Repeated_Unit_ID = factor(make_repeated_unit(Sample_core, `Sample Source`))
)
if (nrow(meta) != 500L || anyDuplicated(meta$Sample_core_for_mpa) ||
    anyNA(meta$Batch) || anyNA(meta$Repeated_Unit_ID) ||
    n_distinct(meta$Repeated_Unit_ID) != 29L) {
  stop("Metadata validation failed: expected 500 samples, 14 batches, and 29 repeated units.")
}

reads <- read_csv(reads_file, show_col_types = FALSE) %>%
  transmute(Sample_core_for_mpa = Sample, Reads_Processed = as.numeric(Reads_Processed))
meta <- left_join(meta, reads, by = "Sample_core_for_mpa")
if (anyNA(meta$Reads_Processed)) stop("Reads summary did not match all samples.")

message("Reading terminal-genus profiles")
taxonomy <- read_tsv(
  taxonomy_file, show_col_types = FALSE, progress = interactive(),
  name_repair = "minimal", comment = "#"
)
if (length(setdiff(meta$Sample_core_for_mpa, names(taxonomy)))) stop("Taxonomy table lacks metadata samples.")
genus <- taxonomy %>%
  filter(
    str_detect(ID, "(^|\\|)g__[^|]+$"),
    !str_detect(ID, regex("unclassified|unknown|uncultured|metagenome", ignore_case = TRUE))
  ) %>%
  mutate(Genus = str_remove(str_extract(ID, "g__[^|]+$"), "^g__")) %>%
  select(Genus, all_of(meta$Sample_core_for_mpa)) %>%
  group_by(Genus) %>%
  summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
rm(taxonomy); invisible(gc())

genus_matrix <- genus %>% column_to_rownames("Genus") %>% as.matrix()
storage.mode(genus_matrix) <- "double"
genus_matrix[!is.finite(genus_matrix) | genus_matrix < 0] <- 0
genus_matrix <- t(genus_matrix)[meta$Sample_core_for_mpa, , drop = FALSE]
genus_relative <- genus_matrix / rowSums(genus_matrix)
genus_relative[genus_relative < min_relative_abundance] <- 0
genus_relative <- genus_relative[, colMeans(genus_relative > 0) >= min_prevalence, drop = FALSE]
genus_relative <- genus_relative / rowSums(genus_relative)
message("Retained ", ncol(genus_relative), " genera")

alpha <- tibble(
  Sample_core_for_mpa = rownames(genus_relative),
  Richness = rowSums(genus_relative > 0),
  Shannon = diversity(genus_relative, "shannon"),
  Simpson = diversity(genus_relative, "simpson")
)
analysis_data <- meta %>%
  left_join(alpha, by = "Sample_core_for_mpa") %>%
  arrange(match(Sample_core_for_mpa, rownames(genus_relative))) %>%
  mutate(Log10_Reads = log10(Reads_Processed))
stopifnot(identical(as.character(analysis_data$Sample_core_for_mpa), rownames(genus_relative)))

# Supplementary table 1: per-batch descriptive summary.
batch_summary <- analysis_data %>%
  group_by(Batch) %>%
  summarise(
    N = n(),
    Reads_median_million = median(Reads_Processed) / 1e6,
    Reads_IQR_million = IQR(Reads_Processed) / 1e6,
    Shannon_mean = mean(Shannon), Shannon_SD = sd(Shannon),
    Richness_mean = mean(Richness), Richness_SD = sd(Richness),
    Simpson_mean = mean(Simpson), Simpson_SD = sd(Simpson),
    .groups = "drop"
  )
write_csv(batch_summary, file.path(out_dir, "sup_table_R1Q4_batch_summary.csv"))

message("Computing Bray-Curtis PCoA and PERMANOVA")
bray <- vegdist(genus_relative, method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
axis_percent <- 100 * positive_eigenvalues[1:2] / sum(positive_eigenvalues)
pcoa_data <- as.data.frame(pcoa$points) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("Sample_core_for_mpa") %>%
  left_join(analysis_data, by = "Sample_core_for_mpa")

permanova_additive <- adonis2(
  bray ~ Facet_Group + Batch, data = analysis_data,
  permutations = permutation_n, by = "margin",
  strata = analysis_data$Repeated_Unit_ID
)
permanova_interaction <- adonis2(
  bray ~ Facet_Group * Batch, data = analysis_data,
  permutations = permutation_n, by = "margin",
  strata = analysis_data$Repeated_Unit_ID
)

extract_permanova <- function(model, model_name) {
  as.data.frame(model) %>% rownames_to_column("Effect") %>%
    filter(Effect %in% c("Facet_Group", "Batch", "Facet_Group:Batch")) %>%
    transmute(
      Analysis = "PERMANOVA", Outcome = "Genus Bray-Curtis composition",
      Model = model_name, Effect, Df, Statistic = "pseudo-F",
      Statistic_value = F, R2, P = `Pr(>F)`
    )
}

mixed_model_test <- function(outcome, interaction = FALSE) {
  response <- if (outcome == "Richness") "log1p(Richness)" else outcome
  reduced <- if (interaction) {
    paste(response, "~ Facet_Group + Batch + (1|Repeated_Unit_ID)")
  } else {
    paste(response, "~ Facet_Group + (1|Repeated_Unit_ID)")
  }
  full <- if (interaction) {
    paste(response, "~ Facet_Group * Batch + (1|Repeated_Unit_ID)")
  } else {
    paste(response, "~ Facet_Group + Batch + (1|Repeated_Unit_ID)")
  }
  fit_reduced <- lmer(as.formula(reduced), analysis_data, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  fit_full <- lmer(as.formula(full), analysis_data, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  comparison <- anova(fit_reduced, fit_full)
  tibble(
    Analysis = "Linear mixed-effects model", Outcome = outcome,
    Model = if (interaction) "Facet_Group * Batch" else "Facet_Group + Batch",
    Effect = if (interaction) "Facet_Group:Batch" else "Batch",
    Df = comparison$Df[2], Statistic = "Likelihood-ratio chi-square",
    Statistic_value = comparison$Chisq[2], R2 = NA_real_, P = comparison$`Pr(>Chisq)`[2]
  )
}

global_tests <- bind_rows(
  extract_permanova(permanova_additive, "Facet_Group + Batch"),
  extract_permanova(permanova_interaction, "Facet_Group * Batch"),
  bind_rows(lapply(c("Log10_Reads", "Shannon", "Richness", "Simpson"), function(outcome) {
    bind_rows(mixed_model_test(outcome, FALSE), mixed_model_test(outcome, TRUE))
  }))
) %>% mutate(Q = p.adjust(P, method = "BH"))
write_csv(global_tests, file.path(out_dir, "sup_table_R1Q4_global_tests.csv"))

within_batch_test <- function(batch_value, grouping_variable) {
  index <- which(analysis_data$Batch == batch_value)
  metadata_subset <- droplevels(analysis_data[index, , drop = FALSE])
  distance_subset <- as.dist(as.matrix(bray)[index, index, drop = FALSE])
  result <- as.data.frame(adonis2(
    reformulate(grouping_variable, response = "distance_subset"),
    data = metadata_subset, permutations = permutation_n
  ))[1, ]
  tibble(
    Batch = batch_value, Grouping = grouping_variable, N = length(index),
    Df = result$Df, R2 = result$R2, pseudo_F = result$F, P = result$`Pr(>F)`
  )
}
within_batch <- bind_rows(lapply(batch_levels, function(batch_value) {
  bind_rows(
    within_batch_test(batch_value, "Facet_Group"),
    within_batch_test(batch_value, "Sample_Source")
  )
})) %>%
  mutate(Batch = factor(Batch, levels = batch_levels), Q = p.adjust(P, method = "BH"))
write_csv(within_batch, file.path(out_dir, "sup_table_R1Q4_within_batch_source_effects.csv"))

theme_publication <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.tag = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold")
  )
facet_palette <- c(
  "Human gut" = "#b63233", "Human nasal" = "#f4583a",
  "Exotic bird" = "#3d5380", "Environment" = "#355d62"
)
batch_palette <- setNames(hcl.colors(14, "Viridis"), batch_levels)

panel_a <- ggplot(analysis_data, aes(Batch, Reads_Processed / 1e6, fill = Batch)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.15), size = 0.7, alpha = 0.4) +
  scale_y_log10() + scale_fill_manual(values = batch_palette, guide = "none") +
  labs(x = "Batch", y = "Reads processed (million; log scale)") + theme_publication

panel_b <- ggplot(
  analysis_data,
  aes(Batch, Shannon, color = Facet_Group,
      group = interaction(Facet_Group, Repeated_Unit_ID))
) +
  geom_line(alpha = 0.16, linewidth = 0.25) +
  stat_summary(aes(group = Facet_Group), fun = mean, geom = "line", linewidth = 1) +
  stat_summary(aes(group = Facet_Group), fun = mean, geom = "point", size = 1.7) +
  scale_color_manual(values = facet_palette) +
  labs(x = "Batch", y = "Genus-level Shannon index", color = NULL) +
  theme_publication + theme(legend.position = "top")

panel_c <- ggplot(pcoa_data, aes(PCoA1, PCoA2, color = Batch)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = batch_palette) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", axis_percent[1]),
    y = sprintf("PCoA2 (%.1f%%)", axis_percent[2]), color = "Batch"
  ) + theme_publication + theme(axis.text.x = element_text(angle = 0))

panel_d <- ggplot(within_batch, aes(Batch, R2, color = Grouping, group = Grouping)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Q <= 0.05), size = 2.1) +
  scale_color_manual(
    values = c(Facet_Group = "#7f1d1d", Sample_Source = "#1f4b52"),
    labels = c(Facet_Group = "Broad group", Sample_Source = "Sample source")
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "FDR <= 0.05", `FALSE` = "FDR > 0.05")
  ) +
  labs(
    x = "Batch", y = expression("Within-batch source effect (PERMANOVA " * R^2 * ")"),
    color = NULL, shape = NULL
  ) + theme_publication + theme(legend.position = "top")

supplementary_figure <- (panel_a + panel_b) / (panel_c + panel_d) +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(out_dir, "sup_fig_R1Q4_batch_effect.pdf"), supplementary_figure,
  width = 12, height = 10, units = "in", device = cairo_pdf
)
ggsave(
  file.path(out_dir, "sup_fig_R1Q4_batch_effect.png"), supplementary_figure,
  width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)

message("R1Q4 compact supplementary outputs completed")

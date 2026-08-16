#!/usr/bin/env Rscript

# Procrustes comparison of genus-level microbiome and ARG-subtype ordinations.
# Permutations are restricted within longitudinal units. Input matrices are
# local derived files and are not distributed with this repository.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(permute)
  library(readr)
  library(stringr)
  library(tibble)
  library(vegan)
})

set.seed(20260811)
argv <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/")
out_dir <- dirname(script_file)
project_dir <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
project_override <- Sys.getenv("OHMD_PROJECT_ROOT", unset = "")
if (nzchar(project_override)) project_dir <- normalizePath(project_override, winslash = "/")

metadata_file <- file.path(project_dir, "data", "metadata_final.csv")
microbiome_file <- file.path(project_dir, "data", "derived", "genus_relative.tsv")
resistome_file <- file.path(project_dir, "data", "derived", "arg_subtype_abundance.tsv")
stopifnot(file.exists(metadata_file), file.exists(microbiome_file), file.exists(resistome_file))

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

read_sample_matrix <- function(path, feature_column = "Feature") {
  table <- read_tsv(path, show_col_types = FALSE, name_repair = "minimal")
  if (!feature_column %in% names(table)) names(table)[1] <- feature_column
  matrix <- table %>% column_to_rownames(feature_column) %>% as.matrix()
  storage.mode(matrix) <- "double"
  matrix[!is.finite(matrix) | matrix < 0] <- 0
  t(matrix)
}

metadata <- read_csv(metadata_file, show_col_types = FALSE, name_repair = "minimal")
if (names(metadata)[1] == "" || str_detect(names(metadata)[1], "^\\.\\.\\.")) metadata <- metadata[, -1]
metadata <- metadata %>% transmute(
  Sample_ID = Sample_core_for_mpa,
  Sample_Source = `Sample Source`,
  Repeated_Unit_ID = factor(make_repeated_unit(Sample_core, `Sample Source`))
)
if (anyNA(metadata$Repeated_Unit_ID)) stop("Unable to derive all longitudinal units")

microbiome <- read_sample_matrix(microbiome_file)
resistome <- read_sample_matrix(resistome_file)
samples <- Reduce(intersect, list(metadata$Sample_ID, rownames(microbiome), rownames(resistome)))
if (length(samples) < 4L) stop("Fewer than four matched samples")
metadata <- metadata[match(samples, metadata$Sample_ID), , drop = FALSE]
microbiome <- microbiome[samples, , drop = FALSE]
resistome <- resistome[samples, , drop = FALSE]

drop_empty <- function(x) {
  x <- x[rowSums(x) > 0, colSums(x) > 0, drop = FALSE]
  x / rowSums(x)
}
microbiome <- drop_empty(microbiome)
resistome <- drop_empty(resistome)
samples <- intersect(rownames(microbiome), rownames(resistome))
metadata <- metadata[match(samples, metadata$Sample_ID), , drop = FALSE]
microbiome <- microbiome[samples, , drop = FALSE]
resistome <- resistome[samples, , drop = FALSE]

microbiome_nmds <- metaMDS(vegdist(microbiome, method = "bray"), k = 2,
                           trymax = 100, autotransform = FALSE, trace = FALSE)
resistome_nmds <- metaMDS(vegdist(resistome, method = "bray"), k = 2,
                          trymax = 100, autotransform = FALSE, trace = FALSE)
fit <- procrustes(microbiome_nmds, resistome_nmds, symmetric = TRUE)
permutation_design <- how(nperm = 999, blocks = metadata$Repeated_Unit_ID)
test <- protest(microbiome_nmds, resistome_nmds, permutations = permutation_design)

micro_scores <- as.data.frame(scores(microbiome_nmds, display = "sites"))
resistome_scores <- as.data.frame(fit$Yrot)
names(micro_scores)[1:2] <- c("microbiome_x", "microbiome_y")
names(resistome_scores)[1:2] <- c("resistome_x", "resistome_y")
micro_scores$Sample_ID <- rownames(micro_scores)
resistome_scores$Sample_ID <- rownames(resistome_scores)
coordinates <- metadata %>%
  left_join(micro_scores, by = "Sample_ID") %>%
  left_join(resistome_scores, by = "Sample_ID") %>%
  mutate(Residual = sqrt((microbiome_x - resistome_x)^2 + (microbiome_y - resistome_y)^2))
write_csv(coordinates, file.path(out_dir, "procrustes_sample_coordinates.csv"))
write_csv(
  tibble(N_samples = nrow(coordinates), Procrustes_m2 = fit$ss,
         PROTEST_correlation = unname(test$t0), PROTEST_P = test$signif,
         Permutations = 999L),
  file.path(out_dir, "procrustes_summary.csv")
)

plot <- ggplot(coordinates) +
  geom_segment(aes(x = microbiome_x, y = microbiome_y,
                   xend = resistome_x, yend = resistome_y,
                   color = Sample_Source), alpha = 0.35, linewidth = 0.35) +
  geom_point(aes(microbiome_x, microbiome_y, color = Sample_Source), size = 1.4) +
  geom_point(aes(resistome_x, resistome_y, color = Sample_Source), shape = 1, size = 1.4) +
  labs(x = "Procrustes axis 1", y = "Procrustes axis 2", color = "Sample source") +
  theme_classic(base_family = "Arial", base_size = 11)
ggsave(file.path(out_dir, "microbiome_resistome_procrustes.pdf"), plot,
       width = 8, height = 7, device = cairo_pdf)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(iNEXT)
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
stopifnot(file.exists(metadata_file), file.exists(taxonomy_file))

metadata <- read_csv(metadata_file, show_col_types = FALSE, name_repair = "minimal")
if (names(metadata)[1] == "" || str_detect(names(metadata)[1], "^\\.\\.\\.")) metadata <- metadata[, -1]
stopifnot(!anyDuplicated(metadata$Sample_core_for_mpa))

taxonomy <- read_tsv(taxonomy_file, show_col_types = FALSE, name_repair = "minimal", comment = "#")
missing_samples <- setdiff(metadata$Sample_core_for_mpa, names(taxonomy))
if (length(missing_samples)) stop("Taxonomy table lacks metadata samples: ", paste(head(missing_samples), collapse = ", "))

extract_rank <- function(table, prefix) {
  token <- paste0("(^|\\|)", prefix, "__[^|]+$")
  table %>%
    filter(str_detect(ID, token),
           !str_detect(ID, regex("unclassified|unknown|uncultured|metagenome", ignore_case = TRUE))) %>%
    mutate(Feature = str_remove(str_extract(ID, paste0(prefix, "__[^|]+$")), paste0("^", prefix, "__"))) %>%
    select(Feature, all_of(metadata$Sample_core_for_mpa)) %>%
    group_by(Feature) %>%
    summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
}

build_incidence <- function(feature_table) {
  long <- feature_table %>%
    pivot_longer(-Feature, names_to = "Sample_ID", values_to = "Abundance") %>%
    left_join(metadata %>% select(Sample_ID = Sample_core_for_mpa, Sample_Source = `Sample Source`), by = "Sample_ID") %>%
    mutate(Present = as.integer(Abundance > 0))
  if (anyNA(long$Sample_Source)) stop("Sample identifiers failed to join metadata")

  sample_counts <- long %>% distinct(Sample_ID, Sample_Source) %>% count(Sample_Source, name = "T")
  incidence <- long %>%
    group_by(Sample_Source, Feature) %>%
    summarise(Frequency = sum(Present), .groups = "drop") %>%
    filter(Frequency > 0)

  setNames(lapply(sample_counts$Sample_Source, function(source) {
    c(sample_counts$T[sample_counts$Sample_Source == source],
      incidence$Frequency[incidence$Sample_Source == source])
  }), sample_counts$Sample_Source)
}

run_rank <- function(prefix, label) {
  incidence <- build_incidence(extract_rank(taxonomy, prefix))
  fit <- iNEXT(incidence, q = 0, datatype = "incidence_freq", se = TRUE, nboot = 50)
  estimates <- as_tibble(fit$AsyEst) %>% mutate(Taxonomic_rank = label, .before = 1)
  write_csv(estimates, file.path(out_dir, paste0("gamma_", tolower(label), "_estimates.csv")))

  plot <- ggiNEXT(fit, type = 1, color.var = "Assemblage") +
    labs(x = "Number of sampling units", y = "Taxonomic richness") +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(legend.title = element_blank())
  ggsave(file.path(out_dir, paste0("gamma_", tolower(label), ".pdf")), plot,
         width = 8, height = 6, device = cairo_pdf)
  estimates
}

all_estimates <- bind_rows(
  run_rank("f", "Family"),
  run_rank("g", "Genus"),
  run_rank("s", "Species")
)
write_csv(all_estimates, file.path(out_dir, "gamma_all_ranks_estimates.csv"))


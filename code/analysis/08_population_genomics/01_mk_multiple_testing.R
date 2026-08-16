#!/usr/bin/env Rscript

# R1Q2, part 1: implement the McDonald-Kreitman-style test exactly as the
# Methods describes it, apply informative-count and coverage filters, and
# correct for multiple testing.
#
# The published Fig. 5E/F counts genes with pNpS_variants > 1 & SNV_N_count >= 2
# (instrain.single.species_v2.R:333, :756). That is a threshold on a ratio, not
# the MK framework with Fisher's exact test that the Methods and the figure
# legend describe. This script implements the described test, reports the size
# of the testing family, and contrasts uncorrected with BH-corrected results.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(patchwork)
  library(readr); library(scales); library(tibble)
})
set.seed(20260811)

argv <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/")
out_dir <- dirname(script_file)

find_project_dir <- function(start) {
  current <- start
  for (i in seq_len(6)) {
    if (dir.exists(file.path(current, "data")) && dir.exists(file.path(current, "code"))) return(current)
    current <- normalizePath(file.path(current, ".."), winslash = "/")
  }
  stop("Could not locate the project root above: ", start)
}
project_dir <- find_project_dir(out_dir)
project_override <- Sys.getenv("OHMD_PROJECT_ROOT", unset = "")
if (nzchar(project_override)) project_dir <- normalizePath(project_override, winslash = "/")

focal_genome <- Sys.getenv("OHMD_FOCAL_GENOME", unset = "GF12EGFS_2_bin.26")
instrain_dir <- file.path(
  project_dir, "data", "instrain", "instrain_single_species",
  focal_genome, "combined_output"
)
gene_file <- file.path(instrain_dir, "combined_gene_info_filtered.csv")
stopifnot(file.exists(gene_file))

# Filter thresholds. The coverage/breadth pair mirrors the genome-level filter
# already used for Fig. 5D (instrain.single.species_v2.R:265), which was never
# applied to the gene-level table that feeds the selection panels.
min_gene_coverage <- 5
min_gene_breadth <- 0.5
min_axis_counts <- 3     # polymorphism and substitution axes each need >= 3
alpha_level <- 0.05

message("Reading gene-level inStrain output for ", focal_genome)
gene <- fread(gene_file, showProgress = FALSE)
stopifnot(nrow(gene) > 0)

required_columns <- c(
  "sample_id", "gene", "gene_length", "coverage", "breadth",
  "SNV_N_count", "SNV_S_count", "SNS_N_count", "SNS_S_count",
  "pNpS_variants", "dNdS_substitutions",
  "sample_source", "source_group", "sampling_time"
)
stopifnot(all(required_columns %in% names(gene)))

# ---------------------------------------------------------------- cascade ----
# Pn / Ps : nonsynonymous / synonymous polymorphism  (within-sample SNVs)
# Dn / Ds : nonsynonymous / synonymous substitution  (consensus vs reference)
# NI = (Pn/Ps) / (Dn/Ds); NI < 1 indicates an excess of nonsynonymous
# substitution relative to nonsynonymous polymorphism.
gene <- gene %>%
  rename(Pn = SNV_N_count, Ps = SNV_S_count, Dn = SNS_N_count, Ds = SNS_S_count) %>%
  mutate(across(c(Pn, Ps, Dn, Ds, coverage, breadth), as.numeric))

n_start <- nrow(gene)
step_counts <- tibble(
  Step = "gene x sample rows in inStrain output",
  Retained = n_start, Removed = NA_integer_
)

add_step <- function(table_so_far, label, retained) {
  previous <- tail(table_so_far$Retained, 1)
  bind_rows(table_so_far, tibble(
    Step = label, Retained = retained, Removed = previous - retained
  ))
}

work <- gene %>% filter(!is.na(Pn), !is.na(Ps), !is.na(Dn), !is.na(Ds))
step_counts <- add_step(step_counts, "all four MK counts present", nrow(work))

work <- work %>% filter(!is.na(coverage), !is.na(breadth),
                        coverage >= min_gene_coverage, breadth >= min_gene_breadth)
step_counts <- add_step(
  step_counts,
  sprintf("gene coverage >= %g and breadth >= %g", min_gene_coverage, min_gene_breadth),
  nrow(work)
)

work <- work %>% filter(Ps >= 1, Ds >= 1)
step_counts <- add_step(step_counts, "NI defined (Ps >= 1 and Ds >= 1)", nrow(work))

work <- work %>% filter((Pn + Ps) >= min_axis_counts, (Dn + Ds) >= min_axis_counts)
step_counts <- add_step(
  step_counts,
  sprintf("informative counts (Pn+Ps >= %d and Dn+Ds >= %d)", min_axis_counts, min_axis_counts),
  nrow(work)
)
step_counts$Step[nrow(step_counts)] <- paste0(
  step_counts$Step[nrow(step_counts)], " -- TESTING FAMILY"
)
write_csv(step_counts, file.path(out_dir, "sup_table_R1Q2_filter_cascade.csv"))
message("Testing family: ", nrow(work), " gene x sample tests across ",
        n_distinct(work$sample_id), " samples")

# ------------------------------------------------------------ MK testing ----
message("Running Fisher's exact test on ", nrow(work), " contingency tables")
work <- work %>%
  mutate(
    neutrality_index = (Pn / Ps) / (Dn / Ds),
    fisher_p = vapply(
      seq_len(n()),
      function(i) fisher.test(matrix(c(Pn[i], Ps[i], Dn[i], Ds[i]), nrow = 2, byrow = TRUE))$p.value,
      numeric(1)
    )
  )

work <- work %>%
  group_by(sample_id) %>%
  mutate(q_within_sample = p.adjust(fisher_p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q_global = p.adjust(fisher_p, method = "BH"),
    signal_raw = neutrality_index < 1 & fisher_p <= alpha_level,
    signal_within = neutrality_index < 1 & q_within_sample <= alpha_level,
    signal_global = neutrality_index < 1 & q_global <= alpha_level
  )

# The statistic actually plotted in the published Fig. 5E/F, for comparison.
published_rule <- gene %>%
  filter(!is.na(pNpS_variants), !is.na(Pn), pNpS_variants > 1, Pn >= 2) %>%
  count(sample_id, name = "published_positive_genes")

per_sample <- work %>%
  group_by(sample_id, sample_source, source_group, sampling_time) %>%
  summarise(
    tests = n(),
    mean_gene_coverage = mean(coverage),
    hits_raw = sum(signal_raw),
    hits_within = sum(signal_within),
    hits_global = sum(signal_global),
    .groups = "drop"
  ) %>%
  left_join(published_rule, by = "sample_id") %>%
  mutate(published_positive_genes = coalesce(published_positive_genes, 0L))
write_csv(per_sample, file.path(out_dir, "sup_table_R1Q2_per_sample_summary.csv"))

full_table <- work %>%
  transmute(
    sample_id, sample_source, source_group, sampling_time,
    gene, gene_length, gene_coverage = coverage, gene_breadth = breadth,
    Pn, Ps, Dn, Ds, neutrality_index,
    pNpS_variants, dNdS_substitutions,
    fisher_p, q_within_sample, q_global,
    signal_raw, signal_within, signal_global
  ) %>%
  arrange(fisher_p)
write_csv(full_table, file.path(out_dir, "sup_table_R1Q2_gene_sample_tests.csv"))

summary_table <- tibble(
  Criterion = c(
    "raw P <= 0.05 (as published)",
    "BH within sample, q <= 0.05",
    "BH across all tests, q <= 0.05"
  ),
  Gene_sample_events = c(sum(work$signal_raw), sum(work$signal_within), sum(work$signal_global)),
  Distinct_genes = c(
    n_distinct(work$gene[work$signal_raw]),
    n_distinct(work$gene[work$signal_within]),
    n_distinct(work$gene[work$signal_global])
  ),
  Genes_recurrent_in_2plus_samples = c(
    sum(table(work$gene[work$signal_raw]) >= 2),
    sum(table(work$gene[work$signal_within]) >= 2),
    sum(table(work$gene[work$signal_global]) >= 2)
  )
)
write_csv(summary_table, file.path(out_dir, "sup_table_R1Q2_correction_summary.csv"))
print(summary_table)

# The uncorrected hit rate is compared against the nominal 5% expectation. A
# ratio below 1 means the raw counts never exceeded chance to begin with.
n_tests <- nrow(work)
n_p05 <- sum(work$fisher_p <= alpha_level)
null_check <- tibble(
  Tests = n_tests,
  Tests_with_P_le_0.05 = n_p05,
  Observed_rate = n_p05 / n_tests,
  Expected_rate_uniform_null = alpha_level,
  Observed_over_expected = (n_p05 / n_tests) / alpha_level,
  Tests_with_NI_below_1 = sum(work$neutrality_index < 1),
  Median_fisher_p = median(work$fisher_p)
)
write_csv(null_check, file.path(out_dir, "sup_table_R1Q2_null_calibration.csv"))
print(null_check)

coverage_confound <- suppressWarnings(cor.test(
  per_sample$mean_gene_coverage, per_sample$published_positive_genes, method = "spearman"
))
message(sprintf(
  "Spearman(mean gene coverage, published Fig.5E statistic) = %.3f (P = %.3g)",
  coverage_confound$estimate, coverage_confound$p.value
))

# ------------------------------------------------------------------ plots ----
theme_publication <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    plot.tag = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold")
  )
group_palette <- c(
  farmer = "#b63233", `non-farmer` = "#f4583a",
  animal = "#3d5380", environment = "#355d62"
)
observed_groups <- sort(unique(per_sample$source_group))
if (!all(observed_groups %in% names(group_palette))) {
  group_palette <- setNames(hcl.colors(length(observed_groups), "Dark 3"), observed_groups)
}

cascade_plot <- step_counts %>%
  mutate(Step = factor(Step, levels = rev(Step))) %>%
  ggplot(aes(Retained, Step)) +
  geom_col(fill = "#355d62", width = 0.65) +
  geom_text(aes(label = comma(Retained)), hjust = -0.08, size = 3.2) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.28))) +
  labs(x = "Rows retained", y = NULL, title = "Filter cascade") +
  theme_publication + theme(axis.text.y = element_text(size = 8))

pvalue_plot <- ggplot(work, aes(fisher_p)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "#3d5380", color = "white", linewidth = 0.2) +
  geom_hline(yintercept = n_tests * 0.05, linetype = 2, color = "#b63233") +
  annotate("text", x = 0.5, y = n_tests * 0.05, vjust = -0.8, size = 3.1, color = "#b63233",
           label = "uniform null expectation") +
  scale_y_continuous(labels = comma) +
  labs(x = "Fisher's exact P", y = "Tests",
       title = sprintf("P-value distribution (%s tests)", comma(n_tests))) +
  theme_publication

confound_plot <- ggplot(per_sample, aes(mean_gene_coverage, published_positive_genes)) +
  geom_point(aes(color = source_group), size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = group_palette, name = NULL) +
  labs(
    x = "Mean gene coverage (log scale)",
    y = "Genes called positive, published rule (log scale)",
    title = sprintf("Published statistic tracks coverage (Spearman rho = %.2f)",
                    coverage_confound$estimate)
  ) +
  theme_publication + theme(legend.position = "top")

correction_plot <- summary_table %>%
  mutate(Criterion = factor(Criterion, levels = rev(Criterion))) %>%
  ggplot(aes(Gene_sample_events, Criterion)) +
  geom_col(fill = "#b63233", width = 0.6) +
  geom_text(aes(label = comma(Gene_sample_events)), hjust = -0.15, size = 3.4) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Gene x sample events with NI < 1 and significance", y = NULL,
       title = "Effect of multiple-testing correction") +
  theme_publication + theme(axis.text.y = element_text(size = 9))

supplementary_figure <- (cascade_plot + pvalue_plot) / (confound_plot + correction_plot) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "sup_fig_R1Q2_multiple_testing.pdf"), supplementary_figure,
       width = 13, height = 10, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "sup_fig_R1Q2_multiple_testing.png"), supplementary_figure,
       width = 13, height = 10, units = "in", dpi = 300, bg = "white")

message("R1Q2 part 1 (multiple-testing correction) completed")

#!/usr/bin/env Rscript

# R1Q2, part 3: contrast models and gene-level recurrence.
#
# Two questions not resolved by individual gene-level Fisher tests:
#
# 1. Does the MK signal DIFFER between ecological compartments or over time?
#    The manuscript claims selection intensified after cross-host spread, which
#    is a contrast, not an absolute level. It is fitted as a binomial GLMM on
#    per-sample counts with a sample random intercept and an observation-level
#    random effect for overdispersion.
#
# 2. Do the same genes recur across independent samples more often than chance?
#    Recurrence is tested against a permutation null stratified by sample, gene
#    length and coverage. This reduces, but does not eliminate, sensitivity to
#    heterogeneous gene detectability and demographic processes.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(lme4)
  library(patchwork); library(readr); library(scales); library(tibble); library(tidyr)
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
snv_file <- file.path(instrain_dir, "combined_snv_filtered.csv")
genome_file <- file.path(instrain_dir, "combined_genome_info_filtered.csv")
tests_file <- file.path(out_dir, "sup_table_R1Q2_gene_sample_tests.csv")
stopifnot(file.exists(snv_file), file.exists(genome_file), file.exists(tests_file))

min_site_coverage <- 5
min_genome_coverage <- 5
min_genome_breadth <- 0.5
permutation_n <- 2000L
# The focal lineage is recovered almost entirely from farmers and animals; only
# two environmental samples carry it. Farmer is used as the reference because it
# is the best represented compartment and the origin of the inferred chain.
reference_group <- "farmer"

# ============================================================ contrast GLMM ===
message("Building per-sample MK counts")
genome <- fread(genome_file, showProgress = FALSE)
usable_samples <- genome %>%
  filter(!is.na(coverage), !is.na(breadth),
         coverage >= min_genome_coverage, breadth >= min_genome_breadth) %>%
  pull(sample_id) %>% unique()

snv <- fread(
  snv_file, showProgress = FALSE,
  select = c("sample_id", "position_coverage", "ref_freq", "mutation_type",
             "class", "source_group", "sampling_time")
)
snv <- snv[
  sample_id %in% usable_samples &
    mutation_type %in% c("N", "S") &
    class %in% c("SNV", "con_SNV", "SNS") &
    !is.na(ref_freq) & position_coverage >= min_site_coverage
]
snv[, variant_class := fifelse(class == "SNS", "substitution", "polymorphism")]

mk_counts <- snv[
  , .N, by = .(sample_id, source_group, sampling_time, variant_class, mutation_type)
] %>%
  as_tibble() %>%
  pivot_wider(names_from = mutation_type, values_from = N, values_fill = 0) %>%
  rename(nonsynonymous = N, synonymous = S) %>%
  filter(nonsynonymous + synonymous > 0) %>%
  mutate(
    variant_class = factor(variant_class, levels = c("polymorphism", "substitution")),
    source_group = relevel(factor(source_group), ref = reference_group),
    sampling_time = as.numeric(sampling_time),
    # centred so the main effect is the signal at mean time and the model is
    # numerically identifiable
    sampling_time_c = sampling_time - mean(unique(sampling_time)),
    sample_id = factor(sample_id)
  )
write_csv(mk_counts, file.path(out_dir, "sup_table_R1Q2_per_sample_mk_counts.csv"))
message("Per-sample count rows: ", nrow(mk_counts),
        " (", n_distinct(mk_counts$sample_id), " samples)")

# The coefficient on variant_class is the log MK odds ratio; interactions test
# whether that odds ratio differs by compartment or over time.
fit_model <- function(formula) {
  tryCatch(
    glmer(formula, data = mk_counts, family = binomial,
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    error = function(e) { message("Model failed: ", conditionMessage(e)); NULL }
  )
}
mk_counts$response_success <- mk_counts$nonsynonymous
mk_counts$response_failure <- mk_counts$synonymous
# Binomial denominators here are in the tens of thousands, so residual
# overdispersion would make every standard error far too small. An
# observation-level random effect absorbs it, and is used for all inference.
mk_counts$observation_id <- factor(seq_len(nrow(mk_counts)))

model_base <- fit_model(
  cbind(response_success, response_failure) ~ variant_class +
    (1 | sample_id) + (1 | observation_id)
)
model_group <- fit_model(
  cbind(response_success, response_failure) ~ variant_class * source_group +
    (1 | sample_id) + (1 | observation_id)
)
model_time <- fit_model(
  cbind(response_success, response_failure) ~ variant_class * sampling_time_c +
    (1 | sample_id) + (1 | observation_id)
)
# Naive models without the overdispersion term, kept only for comparison.
model_base_naive <- fit_model(
  cbind(response_success, response_failure) ~ variant_class + (1 | sample_id)
)

tidy_model <- function(model, label) {
  if (is.null(model)) return(NULL)
  coefs <- summary(model)$coefficients
  tibble(
    Model = label,
    Term = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    Std_Error = coefs[, "Std. Error"],
    z = coefs[, "z value"],
    P = coefs[, "Pr(>|z|)"],
    Odds_ratio = exp(coefs[, "Estimate"]),
    CI_low = exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"]),
    CI_high = exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  )
}
model_table <- bind_rows(
  tidy_model(model_base, "variant_class only"),
  tidy_model(model_group, "variant_class x compartment"),
  tidy_model(model_time, "variant_class x sampling time"),
  tidy_model(model_base_naive, "variant_class only, no overdispersion term")
)
write_csv(model_table, file.path(out_dir, "sup_table_R1Q2_contrast_models.csv"))
print(model_table %>% filter(Term != "(Intercept)"), n = 30)

likelihood_tests <- bind_rows(
  if (!is.null(model_base) && !is.null(model_group)) {
    a <- anova(model_base, model_group)
    tibble(Comparison = "compartment interaction", Chisq = a$Chisq[2],
           Df = a$Df[2], P = a$`Pr(>Chisq)`[2])
  },
  if (!is.null(model_base) && !is.null(model_time)) {
    a <- anova(model_base, model_time)
    tibble(Comparison = "sampling-time interaction", Chisq = a$Chisq[2],
           Df = a$Df[2], P = a$`Pr(>Chisq)`[2])
  }
)
write_csv(likelihood_tests, file.path(out_dir, "sup_table_R1Q2_contrast_lrt.csv"))
print(likelihood_tests)

# Per-sample alpha, for display alongside the model.
sample_alpha <- mk_counts %>%
  select(sample_id, source_group, sampling_time, variant_class,
         nonsynonymous, synonymous) %>%
  pivot_wider(names_from = variant_class,
              values_from = c(nonsynonymous, synonymous)) %>%
  filter(
    !is.na(nonsynonymous_polymorphism), !is.na(nonsynonymous_substitution),
    synonymous_polymorphism > 0, synonymous_substitution > 0,
    nonsynonymous_substitution > 0
  ) %>%
  mutate(
    alpha = 1 - (nonsynonymous_polymorphism / synonymous_polymorphism) /
      (nonsynonymous_substitution / synonymous_substitution)
  )
write_csv(sample_alpha, file.path(out_dir, "sup_table_R1Q2_per_sample_alpha.csv"))

# ======================================================== gene recurrence ====
message("Testing gene-level recurrence against a stratified permutation null")
tests <- fread(tests_file, showProgress = FALSE)
tests <- tests[!is.na(gene) & !is.na(sample_id)]
tests[, flagged := as.logical(signal_raw)]
tests[, length_stratum := cut(gene_length, breaks = quantile(gene_length, probs = seq(0, 1, 0.25),
                                                            na.rm = TRUE),
                              include.lowest = TRUE, labels = FALSE)]
tests[, coverage_stratum := cut(gene_coverage, breaks = quantile(gene_coverage, probs = seq(0, 1, 0.25),
                                                                na.rm = TRUE),
                                include.lowest = TRUE, labels = FALSE)]
tests[, stratum := paste(sample_id, length_stratum, coverage_stratum, sep = "|")]

observed <- tests[flagged == TRUE, .(observed_samples = uniqueN(sample_id)), by = gene]
gene_levels <- sort(unique(tests$gene))
observed_full <- data.table(gene = gene_levels)[observed, on = "gene"]
observed_full <- merge(data.table(gene = gene_levels), observed, by = "gene", all.x = TRUE)
observed_full[is.na(observed_samples), observed_samples := 0L]

gene_index <- match(tests$gene, gene_levels)
stratum_index <- tests$stratum
null_counts <- matrix(0L, nrow = length(gene_levels), ncol = permutation_n)
# "Genes hit at least twice" measures how widely flags are spread, not how
# concentrated they are: if the same few genes absorb many flags, that count
# goes down. Concentration is therefore scored by pair counting, which is the
# quantity that recurrent adaptation should inflate.
null_concentration <- numeric(permutation_n)
null_max <- integer(permutation_n)
null_ge5 <- integer(permutation_n)

split_rows <- split(seq_len(nrow(tests)), stratum_index)
flag_vector <- tests$flagged
for (p in seq_len(permutation_n)) {
  permuted <- logical(length(flag_vector))
  for (rows in split_rows) {
    permuted[rows] <- sample(flag_vector[rows])
  }
  hit_genes <- gene_index[permuted]
  tabulated <- tabulate(hit_genes, nbins = length(gene_levels))
  null_counts[, p] <- tabulated
  null_concentration[p] <- sum(tabulated * (tabulated - 1) / 2)
  null_max[p] <- max(tabulated)
  null_ge5[p] <- sum(tabulated >= 5)
  if (p %% 500 == 0) message("  permutation ", p, "/", permutation_n)
}

observed_vector <- observed_full$observed_samples[match(gene_levels, observed_full$gene)]
empirical_p <- (rowSums(null_counts >= observed_vector) + 1) / (permutation_n + 1)
recurrence <- tibble(
  gene = gene_levels,
  observed_samples = observed_vector,
  null_mean = rowMeans(null_counts),
  null_sd = apply(null_counts, 1, sd),
  empirical_P = empirical_p
) %>%
  filter(observed_samples > 0) %>%
  mutate(q_value = p.adjust(empirical_P, method = "BH")) %>%
  arrange(empirical_P, desc(observed_samples))
write_csv(recurrence, file.path(out_dir, "sup_table_R1Q2_gene_recurrence.csv"))

observed_concentration <- sum(observed_vector * (observed_vector - 1) / 2)
observed_max <- max(observed_vector)
observed_ge5 <- sum(observed_vector >= 5)

global_recurrence <- tibble(
  Statistic = c(
    "concentration, sum of gene-wise sample pairs",
    "maximum samples flagging a single gene",
    "genes flagged in >= 5 samples"
  ),
  Observed = c(observed_concentration, observed_max, observed_ge5),
  Null_mean = c(mean(null_concentration), mean(null_max), mean(null_ge5)),
  Null_sd = c(sd(null_concentration), sd(null_max), sd(null_ge5)),
  Empirical_P = c(
    (sum(null_concentration >= observed_concentration) + 1) / (permutation_n + 1),
    (sum(null_max >= observed_max) + 1) / (permutation_n + 1),
    (sum(null_ge5 >= observed_ge5) + 1) / (permutation_n + 1)
  ),
  Permutations = permutation_n
)
write_csv(global_recurrence, file.path(out_dir, "sup_table_R1Q2_recurrence_global.csv"))
print(global_recurrence)
message("Genes with recurrence q <= 0.05: ", sum(recurrence$q_value <= 0.05))

# ------------------------------------------------------------------ plots ----
theme_publication <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    plot.tag = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold")
  )
observed_groups <- sort(unique(as.character(sample_alpha$source_group)))
group_palette <- setNames(hcl.colors(length(observed_groups), "Dark 3"), observed_groups)

group_sizes <- sample_alpha %>% count(source_group) %>%
  mutate(label = sprintf("%s\n(n = %d)", source_group, n))
alpha_group_plot <- sample_alpha %>%
  left_join(group_sizes, by = "source_group") %>%
  ggplot(aes(label, alpha, fill = source_group)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey40") +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.75) +
  geom_point(position = position_jitter(width = 0.15), size = 1.1, alpha = 0.6) +
  scale_fill_manual(values = group_palette, guide = "none") +
  labs(x = NULL, y = expression("Per-sample " * alpha),
       title = "Adaptive fraction by compartment") +
  theme_publication + theme(axis.text.x = element_text(angle = 20, hjust = 1))

time_effect <- model_table %>% filter(Term == "variant_classsubstitution:sampling_time_c")
alpha_time_plot <- ggplot(sample_alpha, aes(sampling_time, alpha)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey40") +
  geom_point(aes(color = source_group), size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.6) +
  scale_color_manual(values = group_palette, name = NULL) +
  scale_x_continuous(breaks = seq(1, 14, 2)) +
  labs(
    x = "Sampling time", y = expression("Per-sample " * alpha),
    title = if (nrow(time_effect)) {
      sprintf("Adaptive fraction over time (OR = %.3f/timepoint, P = %.2g)",
              time_effect$Odds_ratio[1], time_effect$P[1])
    } else "Adaptive fraction over time"
  ) +
  theme_publication + theme(legend.position = "top")

interaction_terms <- model_table %>%
  filter(grepl(":", Term)) %>%
  mutate(Term = gsub("variant_classsubstitution:", "", Term),
         Term = gsub("source_group", "", Term))
interaction_plot <- if (nrow(interaction_terms)) {
  ggplot(interaction_terms, aes(Odds_ratio, Term)) +
    geom_vline(xintercept = 1, linetype = 3, color = "grey40") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                  width = 0.18, linewidth = 0.5) +
    geom_point(size = 2.4, color = "#b63233") +
    scale_x_log10() +
    labs(x = "Interaction odds ratio (log scale)", y = NULL,
         title = "Does the MK signal differ by context?") +
    theme_publication
} else {
  ggplot() + theme_void()
}

recurrence_plot <- ggplot(tibble(null = null_concentration), aes(null)) +
  geom_histogram(bins = 40, fill = "#3d5380", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = observed_concentration, color = "#b63233", linewidth = 0.9) +
  annotate("text", x = observed_concentration, y = Inf, vjust = 1.6, hjust = 1.06,
           size = 3.2, color = "#b63233", label = "observed") +
  scale_x_continuous(labels = comma) +
  labs(x = "Concentration (sum of gene-wise sample pairs)", y = "Permutations",
       title = sprintf("Recurrence: observed %s vs null %s",
                       comma(observed_concentration), comma(round(mean(null_concentration))))) +
  theme_publication

supplementary_figure <- (alpha_group_plot + alpha_time_plot) /
  (interaction_plot + recurrence_plot) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "sup_fig_R1Q2_contrast_recurrence.pdf"), supplementary_figure,
       width = 13, height = 10, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "sup_fig_R1Q2_contrast_recurrence.png"), supplementary_figure,
       width = 13, height = 10, units = "in", dpi = 300, bg = "white")

message("R1Q2 part 3 (contrast and recurrence) completed")

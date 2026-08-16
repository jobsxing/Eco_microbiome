#!/usr/bin/env Rscript

# R1Q2, part 2: site-frequency-spectrum analysis and asymptotic MK.
#
# The gene x sample contingency tables collapse every variant to "polymorphic"
# or "fixed" and discard allele frequency. Because slightly deleterious
# nonsynonymous variants segregate at low frequency, that collapse is exactly
# what makes the standard MK test non-robust (Messer & Petrov 2013, already
# reference 32 of the manuscript). This script keeps the frequency axis:
#
#   derived allele frequency  = 1 - ref_freq
#   polymorphism (P)          = inStrain SNV and con_SNV classes
#   substitution (D)          = inStrain SNS class
#   alpha(x) = 1 - (Ds/Dn) * (Pn(x)/Ps(x))
#
# alpha(x) is estimated per derived-frequency bin and extrapolated to x -> 1,
# which removes the low-frequency deleterious contribution that biases the
# single 2x2 test downwards.

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(patchwork)
  library(readr); library(scales); library(tibble); library(tidyr)
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
stopifnot(file.exists(snv_file), file.exists(genome_file))

min_site_coverage <- 5      # per-site read depth
min_genome_coverage <- 5    # sample-level genome coverage, matching Fig. 5D
min_genome_breadth <- 0.5
bin_width <- 0.05
min_bin_variants <- 50      # bins with fewer S variants are unstable
bootstrap_n <- 2000L

message("Reading sample-level genome statistics")
genome <- fread(genome_file, showProgress = FALSE)
usable_samples <- genome %>%
  filter(!is.na(coverage), !is.na(breadth),
         coverage >= min_genome_coverage, breadth >= min_genome_breadth) %>%
  pull(sample_id) %>% unique()
message("Samples passing genome coverage/breadth: ", length(usable_samples))

message("Reading SNV table (this file is large)")
snv <- fread(
  snv_file, showProgress = FALSE,
  select = c("sample_id", "scaffold", "position", "position_coverage",
             "ref_freq", "mutation_type", "class", "gene",
             "sample_source", "source_group", "sampling_time")
)
message("SNV rows read: ", format(nrow(snv), big.mark = ","))

snv <- snv[
  sample_id %in% usable_samples &
    mutation_type %in% c("N", "S") &
    class %in% c("SNV", "con_SNV", "SNS") &
    !is.na(ref_freq) & position_coverage >= min_site_coverage
]
snv[, derived_freq := 1 - ref_freq]
snv[, variant_class := fifelse(class == "SNS", "substitution", "polymorphism")]
message("Rows retained after filtering: ", format(nrow(snv), big.mark = ","))

# --------------------------------------------------- genome-wide pooled MK ---
pooled <- snv[, .N, by = .(variant_class, mutation_type)] %>%
  as_tibble() %>%
  pivot_wider(names_from = mutation_type, values_from = N, values_fill = 0)
Pn <- pooled$N[pooled$variant_class == "polymorphism"]
Ps <- pooled$S[pooled$variant_class == "polymorphism"]
Dn <- pooled$N[pooled$variant_class == "substitution"]
Ds <- pooled$S[pooled$variant_class == "substitution"]

pooled_test <- fisher.test(matrix(c(Pn, Ps, Dn, Ds), nrow = 2, byrow = TRUE))
pooled_table <- tibble(
  Pn = Pn, Ps = Ps, Dn = Dn, Ds = Ds,
  neutrality_index = (Pn / Ps) / (Dn / Ds),
  alpha_standard_MK = 1 - (Pn / Ps) / (Dn / Ds),
  odds_ratio = as.numeric(pooled_test$estimate),
  CI_low = pooled_test$conf.int[1], CI_high = pooled_test$conf.int[2],
  fisher_p = pooled_test$p.value
)
write_csv(pooled_table, file.path(out_dir, "sup_table_R1Q2_pooled_mk.csv"))
message("Genome-wide pooled MK:")
print(pooled_table)

# ------------------------------------------------ frequency-resolved alpha ---
breaks <- seq(0, 1, by = bin_width)
poly <- snv[variant_class == "polymorphism"]
poly[, freq_bin := cut(derived_freq, breaks = breaks, include.lowest = TRUE, right = TRUE)]

spectrum <- poly[, .N, by = .(freq_bin, mutation_type)] %>%
  as_tibble() %>%
  pivot_wider(names_from = mutation_type, values_from = N, values_fill = 0) %>%
  rename(Pn_bin = N, Ps_bin = S) %>%
  filter(!is.na(freq_bin)) %>%
  mutate(
    bin_midpoint = breaks[as.integer(freq_bin)] + bin_width / 2,
    NS_ratio = Pn_bin / Ps_bin,
    alpha_x = 1 - (Ds / Dn) * (Pn_bin / Ps_bin)
  ) %>%
  arrange(bin_midpoint)
write_csv(spectrum, file.path(out_dir, "sup_table_R1Q2_site_frequency_spectrum.csv"))
print(spectrum %>% select(bin_midpoint, Pn_bin, Ps_bin, NS_ratio, alpha_x))

fit_data <- spectrum %>% filter(Ps_bin >= min_bin_variants, is.finite(alpha_x))
message("Bins used for the asymptotic fit: ", nrow(fit_data))

fit_asymptotic <- function(d) {
  if (nrow(d) < 4) return(NULL)
  tryCatch(
    nls(alpha_x ~ a + b * exp(-c * bin_midpoint), data = d,
        start = list(a = max(d$alpha_x), b = -1, c = 1),
        control = nls.control(maxiter = 500, warnOnly = TRUE)),
    error = function(e) NULL
  )
}
asymptotic_model <- fit_asymptotic(fit_data)

if (!is.null(asymptotic_model)) {
  coefs <- coef(asymptotic_model)
  alpha_asymptotic <- unname(coefs["a"] + coefs["b"] * exp(-coefs["c"]))
} else {
  alpha_asymptotic <- NA_real_
  message("Exponential fit did not converge; reporting the observed curve only")
}

# Bootstrap over samples, so the interval reflects between-sample variation
# rather than the inflated number of variant records. Counts are aggregated per
# sample first, so a replicate is a sum over resampled rows of a small matrix.
per_sample_counts <- snv[
  , .N, by = .(sample_id, variant_class, mutation_type)
] %>%
  as_tibble() %>%
  mutate(cell = paste(variant_class, mutation_type, sep = "_")) %>%
  select(sample_id, cell, N) %>%
  pivot_wider(names_from = cell, values_from = N, values_fill = 0)
for (needed in c("polymorphism_N", "polymorphism_S", "substitution_N", "substitution_S")) {
  if (!needed %in% names(per_sample_counts)) per_sample_counts[[needed]] <- 0L
}
count_matrix <- as.matrix(per_sample_counts[
  , c("polymorphism_N", "polymorphism_S", "substitution_N", "substitution_S")
])
n_samples <- nrow(count_matrix)

boot_alpha <- vapply(seq_len(bootstrap_n), function(i) {
  totals <- colSums(count_matrix[sample.int(n_samples, n_samples, replace = TRUE), , drop = FALSE])
  if (any(totals == 0)) return(NA_real_)
  1 - (totals[["polymorphism_N"]] / totals[["polymorphism_S"]]) /
    (totals[["substitution_N"]] / totals[["substitution_S"]])
}, numeric(1))
boot_alpha <- boot_alpha[is.finite(boot_alpha)]

alpha_summary <- tibble(
  Estimate = c("alpha, standard MK (pooled)", "alpha, asymptotic extrapolation"),
  Value = c(pooled_table$alpha_standard_MK, alpha_asymptotic),
  Bootstrap_CI_low = c(unname(quantile(boot_alpha, 0.025)), NA_real_),
  Bootstrap_CI_high = c(unname(quantile(boot_alpha, 0.975)), NA_real_),
  Bootstrap_replicates = c(length(boot_alpha), NA_integer_)
)
write_csv(alpha_summary, file.path(out_dir, "sup_table_R1Q2_alpha_estimates.csv"))
print(alpha_summary)

# ------------------------------------------------------- spectrum by group ---
group_spectrum <- poly[, .N, by = .(source_group, freq_bin, mutation_type)] %>%
  as_tibble() %>%
  pivot_wider(names_from = mutation_type, values_from = N, values_fill = 0) %>%
  rename(Pn_bin = N, Ps_bin = S) %>%
  filter(!is.na(freq_bin), Ps_bin >= min_bin_variants) %>%
  mutate(
    bin_midpoint = breaks[as.integer(freq_bin)] + bin_width / 2,
    NS_ratio = Pn_bin / Ps_bin
  )
write_csv(group_spectrum, file.path(out_dir, "sup_table_R1Q2_spectrum_by_group.csv"))

# ------------------------------------------------------------------ plots ----
theme_publication <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    plot.tag = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold")
  )
observed_groups <- sort(unique(group_spectrum$source_group))
group_palette <- setNames(hcl.colors(length(observed_groups), "Dark 3"), observed_groups)

substitution_ratio <- Dn / Ds
spectrum_plot <- ggplot(spectrum %>% filter(Ps_bin >= min_bin_variants),
                        aes(bin_midpoint, NS_ratio)) +
  geom_hline(yintercept = substitution_ratio, linetype = 2, color = "#b63233") +
  annotate("text", x = 0.92, y = substitution_ratio, vjust = -0.7, hjust = 1,
           size = 3.1, color = "#b63233", label = "fixed differences (Dn/Ds)") +
  geom_line(linewidth = 0.8, color = "#3d5380") +
  geom_point(aes(size = Ps_bin), color = "#3d5380") +
  scale_size_continuous(labels = comma, name = "Synonymous\nvariants") +
  labs(
    x = "Derived allele frequency", y = "Nonsynonymous / synonymous",
    title = "Site frequency spectrum of the focal lineage"
  ) + theme_publication

alpha_plot <- ggplot(fit_data, aes(bin_midpoint, alpha_x)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey40") +
  geom_line(linewidth = 0.8, color = "#355d62") +
  geom_point(size = 2, color = "#355d62") +
  {
    if (!is.null(asymptotic_model)) {
      geom_line(
        data = tibble(
          bin_midpoint = seq(min(fit_data$bin_midpoint), 1, length.out = 100)
        ) %>% mutate(alpha_x = predict(asymptotic_model, newdata = .)),
        linetype = 2, color = "#b63233", linewidth = 0.7
      )
    } else NULL
  } +
  labs(
    x = "Derived allele frequency", y = expression(alpha * "(x)"),
    title = if (is.finite(alpha_asymptotic)) {
      sprintf("Asymptotic MK: alpha -> %.3f", alpha_asymptotic)
    } else "Frequency-resolved alpha"
  ) + theme_publication

group_plot <- ggplot(group_spectrum, aes(bin_midpoint, NS_ratio, color = source_group)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  scale_color_manual(values = group_palette, name = NULL) +
  labs(x = "Derived allele frequency", y = "Nonsynonymous / synonymous",
       title = "Spectrum by ecological compartment") +
  theme_publication + theme(legend.position = "top")

boot_plot <- ggplot(tibble(alpha = boot_alpha), aes(alpha)) +
  geom_histogram(bins = 50, fill = "#3d5380", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = pooled_table$alpha_standard_MK, color = "#b63233", linewidth = 0.7) +
  labs(x = expression("Bootstrapped " * alpha * " (resampling samples)"), y = "Replicates",
       title = sprintf("alpha = %.3f [%.3f, %.3f]",
                       pooled_table$alpha_standard_MK,
                       quantile(boot_alpha, 0.025), quantile(boot_alpha, 0.975))) +
  theme_publication

supplementary_figure <- (spectrum_plot + alpha_plot) / (group_plot + boot_plot) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "sup_fig_R1Q2_site_frequency_spectrum.pdf"), supplementary_figure,
       width = 13, height = 10, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, "sup_fig_R1Q2_site_frequency_spectrum.png"), supplementary_figure,
       width = 13, height = 10, units = "in", dpi = 300, bg = "white")

message("R1Q2 part 2 (site frequency spectrum) completed")

#!/usr/bin/env Rscript

# R1Q7: repeated-measures-corrected species-level alpha-diversity sensitivity analysis.
# Inference uses HC3 cluster-robust covariance by the 29 longitudinal units.

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sandwich)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(vegan)
})

set.seed(20260813)
argv <- commandArgs(trailingOnly = FALSE)
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/")
out_dir <- dirname(script_file)
rank_arg <- commandArgs(trailingOnly = TRUE)
tax_rank <- if (length(rank_arg)) tools::toTitleCase(tolower(rank_arg[1])) else "Species"
if (!tax_rank %in% c("Family", "Species")) stop("Rank must be Family or Species.")
rank_lower <- tolower(tax_rank)
project_dir <- normalizePath(file.path(out_dir, "..", "..", ".."), winslash = "/")
project_override <- Sys.getenv("OHMD_PROJECT_ROOT", unset = "")
if (nzchar(project_override)) project_dir <- normalizePath(project_override, winslash = "/")
metadata_file <- file.path(project_dir, "data", "metadata_final.csv")
taxonomy_file <- file.path(project_dir, "data", "metaphlan4", "all_sample_taxonomy.tsv")
stopifnot(all(file.exists(c(metadata_file, taxonomy_file))))

batch_levels <- paste0("T", 1:14)
source_levels <- c(
  "Farmer Excrement", "Non-farmer gut", "Farmer Nasal Vestibule",
  "Non-farmer nasal", "Peacock Excrement", "Goose Excrement",
  "Ostrich Excrement", "Goose Farm Soil", "Ostrich Farm Soil",
  "Goose Paddling Pool", "Surrounding Rivers"
)
source_labels <- c(
  "Farmer Excrement" = "Farmer gut", "Non-farmer gut" = "Non-farmer gut",
  "Farmer Nasal Vestibule" = "Farmer nasal", "Non-farmer nasal" = "Non-farmer nasal",
  "Peacock Excrement" = "Peacock gut", "Goose Excrement" = "Goose gut",
  "Ostrich Excrement" = "Ostrich gut", "Goose Farm Soil" = "Goose-associated soil",
  "Ostrich Farm Soil" = "Ostrich-associated soil", "Goose Paddling Pool" = "Goose paddling pool",
  "Surrounding Rivers" = "Surrounding river"
)

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
  Batch = factor(paste0("T", as.integer(`Sampling Time`)), levels = batch_levels),
  Repeated_Unit_ID = factor(make_repeated_unit(Sample_core, `Sample Source`)),
  Broad_Group = factor(case_when(
    `Sample Source` %in% c("Farmer Excrement", "Non-farmer gut",
                           "Farmer Nasal Vestibule", "Non-farmer nasal") ~ "Human",
    `Sample Source` %in% c("Peacock Excrement", "Goose Excrement", "Ostrich Excrement") ~ "Exotic bird",
    TRUE ~ "Environment"
  ), levels = c("Human", "Exotic bird", "Environment"))
)
cluster_n <- n_distinct(meta$Repeated_Unit_ID)
if (nrow(meta) != 500L || anyDuplicated(meta$Sample_core_for_mpa) ||
    anyNA(meta$Repeated_Unit_ID) || cluster_n != 29L) {
  stop("Metadata validation failed: expected 500 samples and 29 repeated units.")
}

message("Reading terminal-", rank_lower, " profiles")
taxonomy <- read_tsv(
  taxonomy_file, show_col_types = FALSE, progress = interactive(),
  name_repair = "minimal", comment = "#"
)
taxon <- taxonomy %>%
  filter(
    str_detect(ID, paste0("(^|\\|)", if_else(tax_rank == "Family", "f__", "s__"), "[^|]+$")),
    !str_detect(ID, regex("unclassified|unknown|uncultured|metagenome", ignore_case = TRUE))
  ) %>%
  mutate(Taxon = str_remove(str_extract(ID, paste0(if_else(tax_rank == "Family", "f__", "s__"), "[^|]+$")), "^[fs]__")) %>%
  select(Taxon, all_of(meta$Sample_core_for_mpa)) %>%
  group_by(Taxon) %>%
  summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
rm(taxonomy); invisible(gc())

taxon_matrix <- taxon %>% column_to_rownames("Taxon") %>% as.matrix()
storage.mode(taxon_matrix) <- "double"
taxon_matrix[!is.finite(taxon_matrix) | taxon_matrix < 0] <- 0
taxon_matrix <- t(taxon_matrix)[meta$Sample_core_for_mpa, , drop = FALSE]
taxon_relative <- taxon_matrix / rowSums(taxon_matrix)
taxon_relative[taxon_relative < 1e-5] <- 0
taxon_relative <- taxon_relative[, colMeans(taxon_relative > 0) >= 0.10, drop = FALSE]
taxon_relative <- taxon_relative / rowSums(taxon_relative)
message("Retained ", ncol(taxon_relative), " ", rank_lower, " taxa")

analysis_data <- meta %>% mutate(
  Shannon = diversity(taxon_relative, "shannon"),
  Richness = rowSums(taxon_relative > 0),
  LogRichness = log1p(Richness),
  Simpson = diversity(taxon_relative, "simpson")
)

panel_sources <- list(
  Human = source_levels[1:4], Animal = source_levels[5:7],
  Environment = source_levels[8:11]
)
contrast_specs <- unlist(lapply(panel_sources, function(sources) {
  pairs <- combn(sources, 2, simplify = FALSE)
  setNames(pairs, vapply(pairs, function(pair) {
    paste(source_labels[pair[1]], "-", source_labels[pair[2]])
  }, character(1)))
}), recursive = FALSE)
contrast_vectors <- lapply(contrast_specs, function(pair) {
  x <- setNames(rep(0, length(source_levels)), source_levels)
  x[pair[1]] <- 1
  x[pair[2]] <- -1
  x
})

fit_metric <- function(metric, response) {
  model <- lm(as.formula(paste(response, "~ Sample_Source + Batch")), data = analysis_data)
  robust_vcov <- sandwich::vcovCL(
    model, cluster = analysis_data$Repeated_Unit_ID,
    type = "HC3", cadjust = TRUE, fix = TRUE
  )

  source_index <- grep("^Sample_Source", names(coef(model)))
  beta <- coef(model)[source_index]
  covariance <- robust_vcov[source_index, source_index, drop = FALSE]
  wald_chisq <- as.numeric(t(beta) %*% solve(covariance, beta))
  numerator_df <- length(source_index)
  wald_f <- wald_chisq / numerator_df
  global_p <- pf(wald_f, numerator_df, cluster_n - 1, lower.tail = FALSE)

  emm <- emmeans(
    model, ~ Sample_Source, vcov. = robust_vcov,
    df = cluster_n - 1
  )
  emm_ci <- as.data.frame(confint(emm, level = 0.95)) %>%
    transmute(
      Metric = metric, Sample_Source,
      Estimate_model_scale = emmean, SE, Df = df,
      Lower_model_scale = lower.CL, Upper_model_scale = upper.CL
    )
  if (metric == "Richness") {
    emm_ci <- emm_ci %>% mutate(
      Estimate = exp(Estimate_model_scale) - 1,
      Lower = exp(Lower_model_scale) - 1,
      Upper = exp(Upper_model_scale) - 1
    )
  } else {
    emm_ci <- emm_ci %>% mutate(
      Estimate = Estimate_model_scale,
      Lower = Lower_model_scale,
      Upper = Upper_model_scale
    )
  }

  focal <- as.data.frame(contrast(emm, method = contrast_vectors, adjust = "none")) %>%
    transmute(
      Metric = metric, Contrast = contrast, Estimate = estimate, SE,
      Df = df, Statistic = t.ratio, P = p.value
    )
  naive_focal <- as.data.frame(
    contrast(emmeans(model, ~ Sample_Source), method = contrast_vectors, adjust = "none")
  ) %>%
    transmute(
      Metric = metric, Contrast = contrast,
      SE_naive = SE, Df_naive = df, Statistic_naive = t.ratio, P_naive = p.value
    )
  focal <- left_join(focal, naive_focal, by = c("Metric", "Contrast"))
  global <- tibble(
    Metric = metric, Effect = "Sample_Source", Numerator_Df = numerator_df,
    Denominator_Df = cluster_n - 1, Wald_F = wald_f, P = global_p,
    Clusters = cluster_n, Cluster_variable = "Repeated_Unit_ID",
    Fixed_covariate = "Batch", Covariance = "HC3 cluster-robust"
  )
  list(adjusted = emm_ci, focal = focal, global = global)
}

fits <- list(
  fit_metric("Shannon", "Shannon"),
  fit_metric("Richness", "LogRichness"),
  fit_metric("Simpson", "Simpson")
)
adjusted <- bind_rows(lapply(fits, `[[`, "adjusted")) %>%
  left_join(distinct(analysis_data, Sample_Source, Broad_Group), by = "Sample_Source") %>%
  mutate(
    Display_Source = factor(source_labels[as.character(Sample_Source)],
                            levels = rev(unname(source_labels[source_levels])))
  )
focal <- bind_rows(lapply(fits, `[[`, "focal")) %>%
  group_by(Metric) %>%
  mutate(
    Q = p.adjust(P, method = "BH"), Significant = Q <= 0.05,
    Q_naive = p.adjust(P_naive, method = "BH"), Significant_naive = Q_naive <= 0.05
  ) %>% ungroup()
global <- bind_rows(lapply(fits, `[[`, "global")) %>%
  mutate(Q = p.adjust(P, method = "BH"))

write_csv(global, file.path(out_dir, paste0("sup_table_R1Q7_", rank_lower, "_global_tests.csv")))
write_csv(focal, file.path(out_dir, paste0("sup_table_R1Q7_", rank_lower, "_all_pairwise_contrasts.csv")))

if (FALSE) { # Superseded table-panel layout retained only for provenance.
group_palette <- c("Human" = "#8B2C3B", "Exotic bird" = "#315A7D", "Environment" = "#2F766D")
status_palette <- c("Higher" = "#126E82", "Lower" = "#C45A3D", "NS" = "#E8EAED")
section_palette <- c("Human" = "#F4E5E8", "Exotic bird" = "#E3EDF4", "Environment" = "#E1F0ED")

unit_summary <- analysis_data %>%
  group_by(Repeated_Unit_ID, Sample_Source, Broad_Group) %>%
  summarise(
    Shannon = mean(Shannon), Richness = mean(Richness), Simpson = mean(Simpson),
    .groups = "drop"
  ) %>%
  mutate(
    Display_Source = factor(source_labels[as.character(Sample_Source)],
                            levels = rev(unname(source_labels[source_levels])))
  )

theme_publication <- theme_minimal(base_family = "Arial", base_size = 10.5) +
  theme(
    axis.text = element_text(color = "#202124"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    panel.grid.major.y = element_line(color = "#ECEFF1", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#ECEFF1", linewidth = 0.35),
    plot.tag = element_text(face = "bold", size = 15, color = "#202124"),
    legend.position = "top", legend.title = element_blank(),
    legend.key.width = unit(12, "pt"),
    plot.margin = margin(5, 7, 4, 5)
  )

plot_metric <- function(metric, x_label, show_y = TRUE, tag) {
  dat <- filter(adjusted, Metric == metric)
  raw_dat <- unit_summary %>% select(Display_Source, Broad_Group, all_of(metric))
  names(raw_dat)[3] <- "Raw_value"
  p <- ggplot(dat, aes(Estimate, Display_Source)) +
    geom_point(
      data = raw_dat, aes(Raw_value, Display_Source), inherit.aes = FALSE,
      position = position_jitter(height = 0.11, width = 0),
      shape = 1, color = "#AEB4BA", stroke = 0.45, size = 1.45, alpha = 0.75
    ) +
    geom_errorbar(
      aes(xmin = Lower, xmax = Upper, color = Broad_Group),
      width = 0, linewidth = 0.85, orientation = "y"
    ) +
    geom_point(
      aes(fill = Broad_Group), shape = 21, color = "white", stroke = 0.8, size = 3.6
    ) +
    scale_color_manual(values = group_palette, guide = "none") +
    scale_fill_manual(values = group_palette) +
    labs(x = x_label, y = NULL, tag = tag) + theme_publication
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

contrast_display <- c(
  "Farmer nasal - farmer gut" = "Farmer nasal vs farmer gut",
  "Farmer gut - non-farmer gut" = "Farmer gut vs non-farmer gut",
  "Farmer nasal - non-farmer nasal" = "Farmer nasal vs non-farmer nasal",
  "Ostrich gut - goose gut" = "Ostrich gut vs goose gut",
  "Ostrich gut - peacock gut" = "Ostrich gut vs peacock gut",
  "Goose gut - peacock gut" = "Goose gut vs peacock gut",
  "Ostrich-associated soil - goose-associated soil" = "Ostrich-associated soil vs goose-associated soil"
)
contrast_order <- names(contrast_display)
contrast_section <- c(rep("Human", 3), rep("Exotic bird", 3), "Environment") %>%
  setNames(contrast_order)

format_q <- function(x) ifelse(x < 0.001, "q<0.001", sprintf("q=%.3f", x))

table_data <- focal %>%
  mutate(
    Comparison = factor(contrast_display[Contrast], levels = rev(unname(contrast_display))),
    Metric = factor(Metric, levels = c("Shannon", "Richness", "Simpson")),
    Status = case_when(
      Significant & Estimate > 0 ~ "Higher",
      Significant & Estimate < 0 ~ "Lower",
      TRUE ~ "NS"
    ),
    Cell = case_when(
      Status == "Higher" ~ paste0("Higher\n", format_q(Q)),
      Status == "Lower" ~ paste0("Lower\n", format_q(Q)),
      TRUE ~ paste0("NS\n", format_q(Q))
    ),
    Text_colour = if_else(Status == "NS", "#4F5459", "white")
  )

label_data <- tibble(
  Comparison = factor(unname(contrast_display), levels = rev(unname(contrast_display))),
  Section = unname(contrast_section)
)

label_plot <- ggplot(label_data, aes(1, Comparison, fill = Section)) +
  geom_tile(color = "white", linewidth = 0.85) +
  geom_text(aes(label = as.character(Comparison)), hjust = 0.03, x = 0.54,
            family = "Arial", size = 3.25, color = "#202124") +
  scale_fill_manual(values = section_palette, guide = "none") +
  scale_x_continuous(position = "top", breaks = 1, labels = "Comparison", expand = c(0, 0)) +
  labs(title = "D") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#202124",
                              margin = margin(b = 2)),
    axis.text.x.top = element_text(face = "bold", size = 10.5, color = "#202124",
                                   margin = margin(b = 6)),
    plot.margin = margin(3, 0, 4, 5)
  )

metric_table <- ggplot(table_data, aes(Metric, Comparison, fill = Status)) +
  geom_tile(color = "white", linewidth = 0.85) +
  geom_text(aes(label = Cell, color = Text_colour), family = "Arial", size = 3.15,
            lineheight = 0.88, fontface = "bold") +
  scale_fill_manual(values = status_palette, breaks = c("Higher", "Lower", "NS")) +
  scale_color_identity() +
  scale_x_discrete(position = "top", expand = c(0, 0),
                   labels = c("Shannon", "Richness\n(log scale)", "Simpson")) +
  theme_void(base_family = "Arial") +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 10.5, color = "#202124",
                                   margin = margin(b = 6)),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 9),
    legend.key.height = unit(9, "pt"), legend.key.width = unit(18, "pt"),
    plot.margin = margin(3, 5, 4, 0)
  )

top_panel <-
  plot_metric("Shannon", "Adjusted Shannon index", TRUE, "A") +
  plot_metric("Richness", "Adjusted genus richness", FALSE, "B") +
  plot_metric("Simpson", "Adjusted Simpson index", FALSE, "C") +
  plot_layout(guides = "collect", widths = c(1.43, 1, 1)) &
  theme(legend.position = "top")

table_panel <- (label_plot + metric_table + plot_layout(widths = c(1.62, 2.2))) +
  plot_annotation(
    title = "D   Prespecified pairwise contrasts after longitudinal correction",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 11.5,
                                             color = "#202124", margin = margin(b = 3)))
  )

supplementary_figure <- (top_panel / table_panel) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Genus-level alpha diversity after accounting for repeated measures",
    subtitle = paste0(
      "Open circles show longitudinal-unit means; filled points and bars show batch-adjusted estimates and 95% CIs. ",
      "Inference uses HC3 SEs clustered across ", cluster_n, " longitudinal units."
    ),
    caption = paste0(
      "In panel D, higher/lower refers to the first group relative to the second; coloured cells indicate BH-FDR q <= 0.05. ",
      "Richness contrasts were tested on log(richness + 1)."
    ),
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 17, color = "#202124"),
      plot.subtitle = element_text(family = "Arial", size = 10.5, color = "#4F5459",
                                   margin = margin(t = 4, b = 8)),
      plot.caption = element_text(family = "Arial", size = 9, color = "#5F6368",
                                  hjust = 0, margin = margin(t = 6))
    )
  )

ggsave(
  file.path(out_dir, "sup_fig_R1Q8_repeated_measures.pdf"), supplementary_figure,
  width = 13.5, height = 10.5, units = "in", device = cairo_pdf
)
ggsave(
  file.path(out_dir, "sup_fig_R1Q8_repeated_measures.png"), supplementary_figure,
  width = 13.5, height = 10.5, units = "in", dpi = 300, bg = "white"
)
}

# Final Fig. 2A-style display: raw sample distributions retain the visual
# language of the submitted figure, while diamonds, intervals, and brackets
# report the repeated-measures-corrected model results.
display_levels <- unname(source_labels[source_levels])
source_palette <- c(
  "Farmer gut" = "#770D0D", "Non-farmer gut" = "#B63233",
  "Farmer nasal" = "#AC5742", "Non-farmer nasal" = "#F4583A",
  "Peacock gut" = "#234076", "Goose gut" = "#537C87", "Ostrich gut" = "#839C97",
  "Goose-associated soil" = "#66B4BA", "Ostrich-associated soil" = "#185550",
  "Goose paddling pool" = "#237074", "Surrounding river" = "#174C53"
)

raw_long <- analysis_data %>%
  transmute(
    Sample_Source,
    Display_Source = factor(source_labels[as.character(Sample_Source)], levels = display_levels),
    Broad_Group_Display = factor(
      recode(as.character(Broad_Group), "Exotic bird" = "Animal"),
      levels = c("Human", "Animal", "Environment")
    ),
    Shannon, Richness, Simpson
  ) %>%
  pivot_longer(c(Shannon, Richness, Simpson), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Shannon", "Richness", "Simpson")))

adjusted_plot <- adjusted %>%
  transmute(
    Metric = factor(Metric, levels = c("Shannon", "Richness", "Simpson")),
    Display_Source = factor(source_labels[as.character(Sample_Source)], levels = display_levels),
    Broad_Group_Display = factor(
      recode(as.character(Broad_Group), "Exotic bird" = "Animal"),
      levels = c("Human", "Animal", "Environment")
    ),
    Estimate, Lower, Upper
  )

contrast_lookup <- tibble(
  Contrast = names(contrast_specs),
  Source_1 = vapply(contrast_specs, `[[`, character(1), 1),
  Source_2 = vapply(contrast_specs, `[[`, character(1), 2)
)
facet_ranges <- raw_long %>%
  group_by(Metric, Broad_Group_Display) %>%
  summarise(
    Y_min = min(Value, na.rm = TRUE), Y_max = max(Value, na.rm = TRUE),
    Y_range = max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE), .groups = "drop"
  )

brackets <- focal %>%
  left_join(contrast_lookup, by = "Contrast") %>%
  mutate(
    Group_1 = factor(source_labels[Source_1], levels = display_levels),
    Group_2 = factor(source_labels[Source_2], levels = display_levels)
  ) %>%
  left_join(
    analysis_data %>% distinct(Sample_Source, Broad_Group) %>%
      transmute(
        Source_1 = as.character(Sample_Source),
        Broad_Group_Display = factor(
          recode(as.character(Broad_Group), "Exotic bird" = "Animal"),
          levels = c("Human", "Animal", "Environment")
        )
      ),
    by = "Source_1"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Shannon", "Richness", "Simpson"))) %>%
  left_join(facet_ranges, by = c("Metric", "Broad_Group_Display")) %>%
  group_by(Metric, Broad_Group_Display) %>%
  mutate(
    Local_1 = match(Source_1, panel_sources[[as.character(first(Broad_Group_Display))]]),
    Local_2 = match(Source_2, panel_sources[[as.character(first(Broad_Group_Display))]]),
    Span = abs(Local_2 - Local_1)
  ) %>%
  arrange(Span, Local_1, Local_2, .by_group = TRUE) %>%
  mutate(
    Bracket_index = row_number(),
    Y = Y_max + pmax(Y_range, abs(Y_max) * 0.08) * (0.10 + 0.115 * (Bracket_index - 1)),
    Tip = pmax(Y_range, abs(Y_max) * 0.08) * 0.025,
    Label = case_when(
      Q < 0.0001 ~ "q<0.0001",
      Q < 0.001 ~ sprintf("q=%.4f", Q),
      TRUE ~ sprintf("q=%.3f", Q)
    )
  ) %>%
  ungroup()

sample_counts <- raw_long %>%
  distinct(Metric, Broad_Group_Display, Display_Source, Sample_Source) %>%
  left_join(
    analysis_data %>% count(Sample_Source, name = "N"),
    by = "Sample_Source"
  )

fig2a_style <- ggplot(raw_long, aes(Display_Source, Value)) +
  geom_rect(
    data = filter(raw_long, Broad_Group_Display == "Human") %>%
      distinct(Metric, Broad_Group_Display),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE, fill = "#F9E9EA"
  ) +
  geom_rect(
    data = filter(raw_long, Broad_Group_Display == "Animal") %>%
      distinct(Metric, Broad_Group_Display),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE, fill = "#F7E6EF"
  ) +
  geom_rect(
    data = filter(raw_long, Broad_Group_Display == "Environment") %>%
      distinct(Metric, Broad_Group_Display),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE, fill = "#E9E8F5"
  ) +
  geom_violin(
    aes(fill = Display_Source), trim = TRUE, width = 0.88, alpha = 0.88,
    linewidth = 0.28, color = "#202124"
  ) +
  geom_boxplot(
    width = 0.22, outlier.shape = NA, fill = NA,
    color = "#202124", linewidth = 0.42
  ) +
  geom_point(
    aes(color = Display_Source), position = position_jitter(width = 0.14, height = 0),
    size = 0.75, alpha = 0.42, show.legend = FALSE
  ) +
  geom_errorbar(
    data = adjusted_plot,
    aes(x = Display_Source, y = Estimate, ymin = Lower, ymax = Upper),
    inherit.aes = FALSE, width = 0.16, linewidth = 0.72, color = "#111111"
  ) +
  geom_point(
    data = adjusted_plot,
    aes(x = Display_Source, y = Estimate), inherit.aes = FALSE,
    shape = 23, size = 2.8, stroke = 0.75, fill = "white", color = "#111111"
  ) +
  geom_segment(
    data = brackets,
    aes(x = Group_1, xend = Group_2, y = Y, yend = Y),
    inherit.aes = FALSE, linewidth = 0.35, color = "#34373A"
  ) +
  geom_segment(
    data = brackets,
    aes(x = Group_1, xend = Group_1, y = Y, yend = Y - Tip),
    inherit.aes = FALSE, linewidth = 0.35, color = "#34373A"
  ) +
  geom_segment(
    data = brackets,
    aes(x = Group_2, xend = Group_2, y = Y, yend = Y - Tip),
    inherit.aes = FALSE, linewidth = 0.35, color = "#34373A"
  ) +
  geom_text(
    data = brackets,
    aes(x = Group_1, y = Y + Tip * 0.45, label = Label),
    inherit.aes = FALSE, hjust = -0.03, vjust = 0, size = 2.55,
    family = "Arial", color = "#34373A"
  ) +
  geom_text(
    data = sample_counts,
    aes(x = Display_Source, y = -Inf, label = N), inherit.aes = FALSE,
    vjust = -0.55, size = 3.0, family = "Arial", color = "#5F6368"
  ) +
  facet_grid(
    Metric ~ Broad_Group_Display, scales = "free", space = "free_x", switch = "y"
  ) +
  scale_fill_manual(values = source_palette, guide = "none", drop = FALSE) +
  scale_color_manual(values = source_palette, guide = "none", drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(position = "right", expand = expansion(mult = c(0.08, 0.06))) +
  labs(
    x = NULL, y = NULL, title = paste(tax_rank, "level"),
    caption = str_wrap(paste0(
      "Violin plots, boxes, and points show the original sample-level distributions (numbers denote samples). ",
      "White diamonds and error bars show batch-adjusted estimates and 95% CIs. All within-panel pairwise comparisons are shown; ",
      "q values are Benjamini-Hochberg adjusted within each diversity metric using HC3 SEs clustered across ",
      cluster_n, " longitudinal units."
    ), width = 175)
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    panel.border = element_rect(color = "#202124", fill = NA, linewidth = 0.65),
    strip.background = element_rect(fill = "#122526", color = "#122526"),
    strip.text = element_text(color = "white", size = 11.5),
    strip.text.y.left = element_text(angle = 90, face = "plain"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 9.2,
                               color = "#303134"),
    axis.text.y = element_text(size = 9.5, color = "#303134"),
    axis.ticks = element_line(linewidth = 0.35, color = "#303134"),
    plot.title = element_text(face = "plain", size = 12, margin = margin(b = 4)),
    plot.caption = element_text(size = 9.2, color = "#5F6368", hjust = 0,
                                margin = margin(t = 8)),
    panel.spacing.x = unit(0.08, "in"), panel.spacing.y = unit(0.10, "in"),
    plot.margin = margin(8, 8, 8, 8)
  )

# Concordance checks: raw source-level pattern versus adjusted estimates, and
# ordinary independent-sample inference versus cluster-robust inference.
raw_model_scale <- analysis_data %>%
  group_by(Sample_Source) %>%
  summarise(
    Shannon = mean(Shannon), Richness = mean(LogRichness), Simpson = mean(Simpson),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Shannon, Richness, Simpson), names_to = "Metric", values_to = "Raw_mean")

pattern_concordance <- adjusted %>%
  select(Metric, Sample_Source, Estimate_model_scale) %>%
  left_join(raw_model_scale, by = c("Metric", "Sample_Source")) %>%
  group_by(Metric) %>%
  summarise(
    Spearman_source_rank = cor(Raw_mean, Estimate_model_scale, method = "spearman"),
    Pearson_source_mean = cor(Raw_mean, Estimate_model_scale),
    .groups = "drop"
  )

raw_contrasts <- bind_rows(lapply(c("Shannon", "Richness", "Simpson"), function(metric) {
  response <- c(Shannon = "Shannon", Richness = "LogRichness", Simpson = "Simpson")[[metric]]
  bind_rows(lapply(names(contrast_specs), function(contrast_name) {
    pair <- contrast_specs[[contrast_name]]
    tibble(
      Metric = metric, Contrast = contrast_name,
      Raw_difference = mean(analysis_data[[response]][analysis_data$Sample_Source == pair[1]]) -
        mean(analysis_data[[response]][analysis_data$Sample_Source == pair[2]])
    )
  }))
}))

contrast_concordance <- focal %>%
  left_join(raw_contrasts, by = c("Metric", "Contrast")) %>%
  summarise(
    Contrasts = n(), Direction_matches = sum(sign(Estimate) == sign(Raw_difference)),
    Direction_concordance = mean(sign(Estimate) == sign(Raw_difference)),
    N_significant_naive = sum(Significant_naive),
    N_significant_clustered = sum(Significant),
    N_significant_both = sum(Significant_naive & Significant),
    Naive_only = sum(Significant_naive & !Significant),
    Clustered_only = sum(!Significant_naive & Significant)
  )

inference_concordance_by_metric <- focal %>%
  group_by(Metric) %>%
  summarise(
    Contrasts = n(), N_significant_naive = sum(Significant_naive),
    N_significant_clustered = sum(Significant),
    N_significant_both = sum(Significant_naive & Significant),
    Naive_only = sum(Significant_naive & !Significant),
    Clustered_only = sum(!Significant_naive & Significant), .groups = "drop"
  )

print(as.data.frame(pattern_concordance), digits = 8)
print(contrast_concordance)
print(inference_concordance_by_metric)

ggsave(
  file.path(out_dir, paste0("R1Q7_", rank_lower, "_all_metrics_all_pairwise_q.pdf")), fig2a_style,
  width = 12.2, height = 10.2, units = "in", device = cairo_pdf
)
ggsave(
  file.path(out_dir, paste0("R1Q7_", rank_lower, "_all_metrics_all_pairwise_q.png")), fig2a_style,
  width = 12.2, height = 10.2, units = "in", dpi = 300, bg = "white"
)

# Export manuscript-sized single-metric panels while preserving the exact
# visual grammar and statistical layers of the combined audit figure.
filter_plot_metric <- function(plot, metric) {
  # ggplot layers are reference objects; deep duplication prevents the second
  # metric export from emptying the first export's explicit layer data.
  out <- unserialize(serialize(plot, NULL))
  if (is.data.frame(out$data) && "Metric" %in% names(out$data)) {
    out$data <- dplyr::filter(out$data, as.character(Metric) == metric)
  }
  for (i in seq_along(out$layers)) {
    layer_data <- out$layers[[i]]$data
    if (is.data.frame(layer_data) && "Metric" %in% names(layer_data)) {
      out$layers[[i]]$data <- dplyr::filter(layer_data, as.character(Metric) == metric)
    }
  }
  out + labs(caption = str_wrap(paste0(
    "Violin plots, boxes, and points show sample-level distributions (numbers denote samples). ",
    "White diamonds and bars show batch-adjusted estimates and 95% CIs. All within-panel pairwise comparisons are shown; ",
    "q values are Benjamini-Hochberg adjusted within this diversity metric using HC3 SEs clustered across ",
    cluster_n, " longitudinal units."
  ), width = 165))
}

shannon_plot <- filter_plot_metric(fig2a_style, "Shannon")
richness_plot <- filter_plot_metric(fig2a_style, "Richness")
ggsave(file.path(out_dir, paste0("R1Q7_", rank_lower, "_Shannon_all_pairwise_q.pdf")), shannon_plot,
       width = 12.2, height = 4.25, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, paste0("R1Q7_", rank_lower, "_Shannon_all_pairwise_q.png")), shannon_plot,
       width = 12.2, height = 4.25, units = "in", dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0("R1Q7_", rank_lower, "_Richness_all_pairwise_q.pdf")), richness_plot,
       width = 12.2, height = 4.25, units = "in", device = cairo_pdf)
ggsave(file.path(out_dir, paste0("R1Q7_", rank_lower, "_Richness_all_pairwise_q.png")), richness_plot,
       width = 12.2, height = 4.25, units = "in", dpi = 300, bg = "white")

two_metric_plot <- (richness_plot + labs(caption = NULL)) /
  (shannon_plot + labs(title = NULL, caption = NULL)) +
  plot_layout(heights = c(1, 1))
saveRDS(two_metric_plot, file.path(out_dir, paste0("R1Q7_", rank_lower, "_two_metric_plot.rds")))
message("R1Q7 ", rank_lower, "-level cluster-robust analysis completed")

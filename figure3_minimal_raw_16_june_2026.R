# Figure 4 - minimal Tania revision
# Reads: data/MS Figure data_5.4.2026.xlsx, sheet: Fig_4_Ehsan
# Saves simple PDF/PNG figures using raw, linear ratios.

required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "readr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install missing package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

# -----------------------------
# Input / output
# -----------------------------
input_file <- file.path("data", "MS Figure data_5.4.2026.xlsx")
sheet_name <- "Fig_4_Ehsan"

if (!file.exists(input_file) && file.exists("MS Figure data_5.4.2026.xlsx")) {
  input_file <- "MS Figure data_5.4.2026.xlsx"
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found. Put 'MS Figure data_5.4.2026.xlsx' in the data/ folder.",
    call. = FALSE
  )
}

output_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read data
# -----------------------------
df <- read_excel(input_file, sheet = sheet_name)
names(df) <- trimws(names(df))

sample_col <- names(df)[grepl("^sample\\s*ID", names(df), ignore.case = TRUE)][1]
baseline_col <- "aCD3/msIgG"

if (is.na(sample_col)) stop("Could not find sample ID column.", call. = FALSE)
if (!baseline_col %in% names(df)) stop("Could not find column: ", baseline_col, call. = FALSE)

df <- df %>%
  rename(Patient = all_of(sample_col), baseline = all_of(baseline_col)) %>%
  mutate(
    Patient = as.character(Patient),
    across(-Patient, ~ suppressWarnings(as.numeric(.x)))
  )

drug_cols <- setdiff(names(df), c("Patient", "baseline"))

baseline_threshold <- 2
response_cutoff <- 1.5
teal_response <- "#008C95"

# -----------------------------
# Patient grouping
# -----------------------------
baseline_ranked <- df %>%
  filter(!is.na(baseline)) %>%
  arrange(baseline, Patient) %>%
  mutate(
    group = if_else(
      baseline < baseline_threshold,
      "Low Proliferators",
      "High Proliferators"
    ),
    group = factor(group, levels = c("Low Proliferators", "High Proliferators")),
    Patient = factor(Patient, levels = Patient)
  )

patient_groups <- baseline_ranked %>%
  select(Patient, group)

# -----------------------------
# Statistics
# -----------------------------
format_p <- function(p) {
  if (!is.finite(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", formatC(p, format = "f", digits = 3))
}

drug_response_stats <- function(data, drug_col, group_name) {
  x <- data %>%
    filter(!is.na(baseline), !is.na(.data[[drug_col]])) %>%
    mutate(group = if_else(
      baseline < baseline_threshold,
      "Low Proliferators",
      "High Proliferators"
    )) %>%
    filter(group == group_name) %>%
    pull(all_of(drug_col))

  if (length(x) < 2) {
    return(tibble(
      drug = drug_col,
      group = group_name,
      n = length(x),
      median = ifelse(length(x) == 0, NA_real_, median(x, na.rm = TRUE)),
      responders_gt_1_5 = sum(x > response_cutoff, na.rm = TRUE),
      p_wilcoxon_vs_1_greater = NA_real_
    ))
  }

  wt <- suppressWarnings(wilcox.test(x, mu = 1, alternative = "greater", exact = FALSE))

  tibble(
    drug = drug_col,
    group = group_name,
    n = length(x),
    median = median(x, na.rm = TRUE),
    responders_gt_1_5 = sum(x > response_cutoff, na.rm = TRUE),
    p_wilcoxon_vs_1_greater = wt$p.value
  )
}

stats_tbl <- bind_rows(
  drug_response_stats(df, "Entospletinib", "Low Proliferators"),
  drug_response_stats(df, "aPD1", "Low Proliferators"),
  drug_response_stats(df, "Entospletinib", "High Proliferators"),
  drug_response_stats(df, "aPD1", "High Proliferators")
)

write_csv(stats_tbl, file.path(table_dir, "figure4_minimal_statistics.csv"))

# -----------------------------
# Simple theme
# -----------------------------
theme_simple <- function() {
  theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),
      plot.tag = element_text(face = "bold", size = 14),
      legend.position = "top"
    )
}

# -----------------------------
# Panel A: violin plots, raw linear values
# -----------------------------
make_panel_a <- function(order_by = c("median", "sd", "tail90")) {
  order_by <- match.arg(order_by)

  drug_stats <- df %>%
    select(Patient, all_of(drug_cols)) %>%
    pivot_longer(-Patient, names_to = "Drug", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    summarise(
      drug_median = median(Value, na.rm = TRUE),
      drug_sd = sd(Value, na.rm = TRUE),
      drug_tail90 = as.numeric(quantile(Value, 0.90, na.rm = TRUE)),
      .by = Drug
    ) %>%
    mutate(
      order_metric = case_when(
        order_by == "median" ~ drug_median,
        order_by == "sd" ~ drug_sd,
        order_by == "tail90" ~ drug_tail90
      )
    ) %>%
    arrange(desc(order_metric), Drug)

  plot_data <- df %>%
    select(Patient, all_of(drug_cols)) %>%
    pivot_longer(-Patient, names_to = "Drug", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(Drug = factor(Drug, levels = drug_stats$Drug))

  ggplot(plot_data, aes(x = Drug, y = Value)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_violin(fill = "grey90", color = "grey30", trim = FALSE) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
    stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.5) +
    labs(
      tag = "A",
      x = NULL,
      y = "Ratio to control"
    ) +
    theme_simple()
}

# -----------------------------
# Panel B: baseline proliferation
# -----------------------------
pB <- baseline_ranked %>%
  ggplot(aes(x = baseline, y = Patient, fill = group)) +
  geom_vline(xintercept = baseline_threshold, linetype = "dashed", color = "black") +
  geom_col(width = 0.65, color = "black") +
  scale_fill_manual(values = c(
    "Low Proliferators" = "#DDEFD8",
    "High Proliferators" = "grey85"
  )) +
  labs(
    tag = "B",
    x = "Baseline stimulation (aCD3/msIgG)",
    y = "Patient",
    fill = NULL,
    caption = "Low proliferators: baseline < 2"
  ) +
  theme_simple() +
  theme(axis.text.x = element_text(angle = 0))

# -----------------------------
# Panels C/D: drug response bars, raw linear values
# -----------------------------
drug_bar_plot <- function(data, drug_col, panel_tag, title_label, group_to_show, show_legend = FALSE) {
  if (!drug_col %in% names(data)) stop("Column not found: ", drug_col, call. = FALSE)

  plot_data <- data %>%
    select(Patient, Value = all_of(drug_col)) %>%
    filter(!is.na(Value)) %>%
    left_join(patient_groups, by = "Patient") %>%
    filter(group == group_to_show) %>%
    mutate(
      Patient = factor(Patient, levels = levels(baseline_ranked$Patient)),
      response_class = if_else(
        Value > response_cutoff,
        "T cell rescue > 1.5",
        "T cell rescue <= 1.5"
      )
    )

  stat_row <- stats_tbl %>%
    filter(drug == drug_col, group == group_to_show)

  p_label <- if (nrow(stat_row) == 0) {
    "p=NA"
  } else {
    format_p(stat_row$p_wilcoxon_vs_1_greater[1])
  }

  ggplot(plot_data, aes(x = Value, y = Patient, fill = response_class)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = response_cutoff, linetype = "dotted", color = teal_response) +
    geom_col(width = 0.65, color = "black") +
    geom_text(aes(label = round(Value, 2)), hjust = -0.15, size = 3) +
    annotate(
      "text",
      x = max(plot_data$Value, na.rm = TRUE) * 0.75,
      y = Inf,
      label = p_label,
      vjust = 1.5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = c(
      "T cell rescue <= 1.5" = "grey90",
      "T cell rescue > 1.5" = teal_response
    )) +
    labs(
      tag = panel_tag,
      x = paste0(title_label, " ratio to control"),
      y = group_to_show,
      fill = "T cell rescue"
    ) +
    coord_cartesian(clip = "off") +
    theme_simple() +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.position = ifelse(show_legend, "top", "none")
    )
}

pA_median <- make_panel_a("median")
pA_sd <- make_panel_a("sd")
pA_tail90 <- make_panel_a("tail90")

pC_low <- drug_bar_plot(df, "Entospletinib", "C", "Entospletinib", "Low Proliferators", TRUE)
pD_low <- drug_bar_plot(df, "aPD1", "D", "aPD1", "Low Proliferators", FALSE)

pC_high <- drug_bar_plot(df, "Entospletinib", "C", "Entospletinib", "High Proliferators", TRUE)
pD_high <- drug_bar_plot(df, "aPD1", "D", "aPD1", "High Proliferators", FALSE)

# -----------------------------
# Main and supplementary figures
# -----------------------------
figure4_main <- pA_median / (pC_low | pB | pD_low) +
  plot_layout(heights = c(1.1, 1))

figure4_supp_high <- pC_high | pB | pD_high

save_plot_pair <- function(plot, stem, width = 14, height = 9, dpi = 300) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi)
}

save_plot_pair(figure4_main, "figure4_minimal_main", width = 15, height = 10)
save_plot_pair(figure4_supp_high, "figure4_minimal_supp_high_proliferators", width = 15, height = 5.5)

save_plot_pair(pA_median, "figure4_minimal_panelA_median", width = 12, height = 6)
save_plot_pair(pA_sd, "figure4_minimal_panelA_sd", width = 12, height = 6)
save_plot_pair(pA_tail90, "figure4_minimal_panelA_tail90", width = 12, height = 6)

message("Done. Minimal raw-scale figures saved in: ", normalizePath(output_dir))

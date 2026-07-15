# Figure 3 - Tania comments, v5 Panel A stacked bars
# Reads: data/MS Figure data_5.4.2026.xlsx, sheet: Fig_4_Ehsan
# Outputs:
#   results/figures/figure3_v5_panelA_stacked.pdf/png
#   results/figures/figure3_v5_supp_high_proliferators.pdf/png
#   results/figures/figure3_v5_panelA_stacked_only.pdf/png
#   results/tables/figure3_v5_data_source.csv
#   results/tables/figure3_v5_statistics.csv
#   results/tables/figure3_v5_panelA_stacked_bins.csv

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "scales",
  "stringr", "ggbeeswarm", "readr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the missing package(s) first: ",
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
library(scales)
library(stringr)
library(ggbeeswarm)
library(readr)

input_file <- file.path("data", "MS Figure data_5.4.2026.xlsx")
sheet_name <- "Fig_4_Ehsan"
output_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

if (!file.exists(input_file) && file.exists("MS Figure data_5.4.2026.xlsx")) {
  input_file <- "MS Figure data_5.4.2026.xlsx"
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPut 'MS Figure data_5.4.2026.xlsx' inside the project's data/ folder.",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

data_created <- suppressWarnings(system2("stat", c("-f", "%SB", shQuote(input_file)), stdout = TRUE))
if (length(data_created) == 0 || is.na(data_created[1])) {
  data_created <- NA_character_
}
data_modified <- file.info(input_file)$mtime

df <- read_excel(input_file, sheet = sheet_name)
names(df) <- trimws(names(df))

sample_col <- names(df)[str_detect(names(df), regex("^sample\\s*ID", ignore_case = TRUE))][1]
baseline_col <- "aCD3/msIgG"

if (is.na(sample_col)) {
  stop("Could not find the sample ID column. Expected a name starting with 'sample ID'.", call. = FALSE)
}
if (!baseline_col %in% names(df)) {
  stop("Could not find the baseline column: ", baseline_col, call. = FALSE)
}

df <- df |>
  rename(Patient = all_of(sample_col), baseline = all_of(baseline_col)) |>
  mutate(
    Patient = as.character(Patient),
    across(-Patient, ~ suppressWarnings(as.numeric(.x)))
  )

drug_cols <- setdiff(names(df), c("Patient", "baseline"))

baseline_threshold <- 2
response_cutoff <- 1.5
y_cap <- 5
group_gap <- 1.8
use_pseudolog_cd <- TRUE

low_bg_colour <- "#FFFFFF"
high_bg_colour <- "#FFFFFF"
speckle_colour <- "#62B45D"
teal_response <- "#008C95"
response_colours <- c("T cell rescue <= 1.5" = "#F6F6F6", "T cell rescue > 1.5" = teal_response)

baseline_ranked <- df |>
  select(Patient, baseline) |>
  filter(!is.na(baseline)) |>
  mutate(
    group = if_else(baseline < baseline_threshold, "Low Proliferators", "High Proliferators"),
    group = factor(group, levels = c("Low Proliferators", "High Proliferators"))
  ) |>
  arrange(desc(baseline), Patient) |>
  group_by(group) |>
  mutate(group_rank = row_number()) |>
  ungroup() |>
  mutate(
    n_low_total = sum(group == "Low Proliferators"),
    n_high_total = sum(group == "High Proliferators"),
    baseline_rank = row_number(),
    y_position = if_else(
      group == "High Proliferators",
      n_low_total + group_gap + n_high_total - group_rank + 1,
      n_low_total - group_rank + 1
    ),
    baseline_disp = pmin(baseline, y_cap),
    is_outlier = baseline > y_cap,
    Patient = factor(Patient, levels = Patient)
  ) |>
  select(-n_low_total, -n_high_total)

n_tot <- nrow(baseline_ranked)
n_low <- sum(baseline_ranked$group == "Low Proliferators")
n_high <- sum(baseline_ranked$group == "High Proliferators")

low_y_min <- 0.5
low_y_max <- n_low + 0.5
high_y_min <- n_low + group_gap + 0.5
high_y_max <- n_low + group_gap + n_high + 0.5
gap_y_mid <- (low_y_max + high_y_min) / 2

patient_order <- as.character(baseline_ranked$Patient)
outlier_labels <- baseline_ranked |> filter(is_outlier)

patient_groups <- baseline_ranked |>
  transmute(Patient = as.character(Patient), group, baseline_rank, y_position)

make_low_speckles <- function(xmin, xmax, n = 260) {
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin || n_low == 0) {
    return(tibble(x = numeric(), y = numeric()))
  }
  set.seed(2305)
  tibble(
    x = runif(n, xmin, xmax),
    y = runif(n, low_y_min, low_y_max)
  )
}

panel_background <- function(xmin, xmax, low_speckles = 260, show_high = TRUE) {
  low_speckles_df <- make_low_speckles(xmin, xmax, low_speckles)
  layers <- list(
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = low_y_min, ymax = low_y_max,
      fill = low_bg_colour, alpha = 1
    ),
    geom_point(
      data = low_speckles_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 16,
      size = 0.45,
      color = speckle_colour,
      alpha = 0.20
    )
  )
  if (isTRUE(show_high)) {
    layers <- c(
      list(
        annotate(
          "rect",
          xmin = -Inf, xmax = Inf,
          ymin = high_y_min, ymax = high_y_max,
          fill = high_bg_colour, alpha = 1
        )
      ),
      layers
    )
  }
  layers
}

theme_fig3 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, color = "grey20"),
      axis.text.y = element_text(color = "grey20"),
      axis.title = element_text(face = "bold", color = "grey15"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      plot.tag = element_text(size = 15, face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 8, 5, 8)
    )
}

format_p <- function(p) {
  if (!is.finite(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", formatC(p, format = "f", digits = 3))
}

drug_response_stats <- function(data, drug_col, group_name) {
  x <- data |>
    filter(!is.na(baseline), !is.na(.data[[drug_col]])) |>
    mutate(group = if_else(baseline < baseline_threshold, "Low Proliferators", "High Proliferators")) |>
    filter(group == group_name) |>
    pull(all_of(drug_col))

  if (length(x) < 2) {
    return(tibble(
      drug = drug_col,
      group = group_name,
      n = length(x),
      median = if_else(length(x) == 0, NA_real_, median(x, na.rm = TRUE)),
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

write_csv(stats_tbl, file.path(table_dir, "figure3_v5_statistics.csv"))

write_csv(
  tibble(
    input_file = normalizePath(input_file),
    sheet = sheet_name,
    file_created_macos_stat = data_created[1],
    file_modified = as.character(data_modified),
    baseline_threshold = baseline_threshold,
    response_cutoff = response_cutoff,
    n_rows = nrow(df),
    n_columns = ncol(df)
  ),
  file.path(table_dir, "figure3_v5_data_source.csv")
)

normalize_drug_name <- function(x) {
  x |>
    str_replace_all("\u00A0", " ") |>
    str_remove("\\s*\\([^)]*\\)") |>
    str_squish()
}

heatmap_single_drug_order <- function(current_drug_cols) {
  summary_file <- "heatmap_drug_waterfall_summary_tables.xlsx"
  if (!file.exists(summary_file)) {
    return(current_drug_cols)
  }

  heatmap_order <- read_excel(summary_file, sheet = "TcellProlif_Drugmedians") |>
    mutate(
      heatmap_drug = as.character(Drug),
      clean_heatmap_drug = normalize_drug_name(heatmap_drug),
      fig3_drug_key = if_else(clean_heatmap_drug == "aPD1+", "aPD1", clean_heatmap_drug)
    ) |>
    filter(clean_heatmap_drug == "aPD1+" | !str_detect(clean_heatmap_drug, "^aPD1\\+ ")) |>
    pull(fig3_drug_key)

  fig3_lookup <- tibble(
    drug = current_drug_cols,
    drug_key = normalize_drug_name(current_drug_cols)
  )

  ordered_drugs <- tibble(drug_key = heatmap_order) |>
    left_join(fig3_lookup, by = "drug_key") |>
    filter(!is.na(drug)) |>
    pull(drug)

  c(ordered_drugs, setdiff(current_drug_cols, ordered_drugs))
}

make_panel_a <- function() {
  drug_order <- heatmap_single_drug_order(drug_cols)

  log2_breaks <- c(-Inf, -2, -1, 0, 1, 2, Inf)
  log2_labels <- c("<= -2", "-2 to -1", "-1 to 0", "0 to 1", "1 to 2", "> 2")
  log2_colours <- c(
    "<= -2" = "#536F84",
    "-2 to -1" = "#8EA4B6",
    "-1 to 0" = "#D7DEE3",
    "0 to 1" = "#C9DEC2",
    "1 to 2" = "#79BB70",
    "> 2" = "#00843D"
  )

  panel_a_data <- df |>
    select(Patient, all_of(drug_cols)) |>
    pivot_longer(-Patient, names_to = "Drug", values_to = "Value") |>
    filter(!is.na(Value), Value > 0) |>
    mutate(
      Drug = factor(Drug, levels = drug_order),
      log2_ratio = log2(Value),
      log2_bin = cut(
        log2_ratio,
        breaks = log2_breaks,
        labels = log2_labels,
        include.lowest = TRUE,
        right = TRUE
      ),
      log2_bin = factor(log2_bin, levels = log2_labels)
    )

  panel_a_summary <- panel_a_data |>
    count(Drug, log2_bin, name = "n") |>
    complete(
      Drug = factor(drug_order, levels = drug_order),
      log2_bin = factor(log2_labels, levels = log2_labels),
      fill = list(n = 0)
    ) |>
    group_by(Drug) |>
    mutate(
      total_n = sum(n),
      fraction = if_else(total_n > 0, n / total_n, 0)
    ) |>
    ungroup()

  write_csv(
    panel_a_summary |>
      mutate(Drug = as.character(Drug), log2_bin = as.character(log2_bin)),
    file.path(table_dir, "figure3_v5_panelA_stacked_bins.csv")
  )

  ggplot(panel_a_summary, aes(x = Drug, y = fraction, fill = log2_bin)) +
    geom_col(width = 0.74, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = log2_colours, drop = FALSE, name = "log2 ratio") +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      tag = "A",
      x = NULL,
      y = "Patients per log2 ratio bin"
    ) +
    theme_fig3(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, size = 9, color = "grey20"),
      axis.title.y = element_text(size = 10, margin = margin(r = 8)),
      legend.position = "right"
    )
}

pA_stacked <- make_panel_a()

pB <- ggplot(baseline_ranked, aes(y = y_position)) +
  panel_background(xmin = 0, xmax = y_cap, low_speckles = 420, show_high = TRUE) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = low_y_max, ymax = high_y_min, fill = "white", alpha = 1) +
  geom_vline(xintercept = baseline_threshold, linetype = "dashed", color = "#333333", linewidth = 0.85) +
  geom_segment(
    aes(x = 0, xend = baseline_disp, yend = y_position),
    linewidth = 0.85,
    color = "grey20",
    alpha = 0.78,
    lineend = "round"
  ) +
  geom_segment(
    data = outlier_labels,
    aes(x = y_cap * 0.82, xend = y_cap * 0.90, y = y_position - 0.25, yend = y_position + 0.25),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 1.7
  ) +
  geom_segment(
    data = outlier_labels,
    aes(x = y_cap * 0.76, xend = y_cap * 0.84, y = y_position - 0.25, yend = y_position + 0.25),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 1.7
  ) +
  geom_point(
    aes(x = baseline_disp),
    shape = 21,
    size = 4.2,
    stroke = 1.1,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = ~ filter(.x, !is_outlier),
    aes(x = baseline_disp, label = round(baseline_disp, 2)),
    hjust = -0.35,
    vjust = 0.5,
    size = 2.6,
    color = "grey30"
  ) +
  geom_text(
    data = outlier_labels,
    aes(x = baseline_disp, y = y_position, label = paste0("^ ", round(baseline, 1))),
    vjust = 0.5,
    hjust = -0.25,
    size = 3.0,
    color = "#333333",
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = baseline_threshold + 0.12,
    y = gap_y_mid,
    label = "threshold = 2",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "grey20",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = y_cap * 0.85,
    y = low_y_min + n_low / 2,
    label = "LOW PROLIFERATORS\nbaseline < 2",
    color = "#2F8F35",
    fontface = "bold",
    size = 3.35,
    hjust = 0.5,
    vjust = 0.5
  ) +
  annotate(
    "text",
    x = y_cap * 0.85,
    y = high_y_min + n_high / 2,
    label = "HIGH PROLIFERATORS\nbaseline >= 2",
    color = "#596B57",
    fontface = "bold",
    size = 3.15,
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_x_continuous(breaks = 0:5, expand = expansion(mult = c(0, 0.14))) +
  scale_y_continuous(
    breaks = baseline_ranked$y_position,
    labels = as.character(baseline_ranked$Patient),
    expand = expansion(add = 0.4)
  ) +
  coord_cartesian(xlim = c(0, y_cap), clip = "off") +
  labs(
    tag = "B",
    x = "Baseline stimulation (aCD3/msIgG)",
    y = "Patient (high to low baseline proliferation)"
  ) +
  theme_fig3(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "grey20"),
    axis.text.y = element_text(size = 8, color = "grey20"),
    legend.position = "none"
  )

drug_bar_plot <- function(data, drug_col, panel_tag, title_label, group_to_show,
                          show_legend = FALSE, use_pseudolog = use_pseudolog_cd) {
  if (!drug_col %in% names(data)) {
    stop("Column not found: ", drug_col, call. = FALSE)
  }

  plot_data <- data |>
    select(Patient, Value = all_of(drug_col)) |>
    filter(!is.na(Value)) |>
    left_join(patient_groups, by = "Patient") |>
    filter(group == group_to_show, !is.na(y_position)) |>
    mutate(Patient = factor(Patient, levels = patient_order)) |>
    arrange(Patient) |>
    mutate(
      response_class = if_else(Value > response_cutoff, "T cell rescue > 1.5", "T cell rescue <= 1.5"),
      response_class = factor(response_class, levels = names(response_colours)),
      label_x = if_else(Value < 0, 0, Value)
    )

  x_min <- min(0, min(plot_data$Value, na.rm = TRUE)) * 1.08
  x_max <- max(response_cutoff * 1.15, max(plot_data$Value, na.rm = TRUE) * 1.08)
  stat_row <- stats_tbl |> filter(drug == drug_col, group == group_to_show)
  p_label <- if (nrow(stat_row) == 0) "p=NA" else format_p(stat_row$p_wilcoxon_vs_1_greater[1])
  group_note <- if_else(
    group_to_show == "Low Proliferators",
    "Low proliferators from Panel B (baseline < 2)",
    "High proliferators from Panel B (baseline >= 2)"
  )
  group_note_colour <- if_else(group_to_show == "Low Proliferators", "#2F8F35", "#596B57")

  y_breaks <- plot_data$y_position
  y_labels <- as.character(plot_data$Patient)
  y_limits <- range(plot_data$y_position, na.rm = TRUE) + c(-0.55, 1.15)

  x_scale <- if (isTRUE(use_pseudolog)) {
    scale_x_continuous(
      trans = pseudo_log_trans(sigma = 0.25, base = 10),
      breaks = c(0, 0.5, 1, response_cutoff, 2, 5, 10, 25),
      labels = c("0", "0.5", "1", "1.5", "2", "5", "10", "25"),
      expand = expansion(mult = c(0.02, 0.16))
    )
  } else {
    linear_breaks <- if (x_max > 12) {
      c(0, 1, response_cutoff, 5, 10, 25)
    } else {
      c(0, 1, response_cutoff, 2, 5, 10)
    }
    scale_x_continuous(
      breaks = linear_breaks,
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.16))
    )
  }

  ggplot(plot_data, aes(y = y_position)) +
    panel_background(xmin = x_min, xmax = x_max, low_speckles = 320, show_high = group_to_show == "High Proliferators") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_rect(
      aes(
        xmin = pmin(0, Value),
        xmax = pmax(0, Value),
        ymin = y_position - 0.35,
        ymax = y_position + 0.35,
        fill = response_class
      ),
      color = "grey20",
      linewidth = 0.25,
      alpha = 0.96
    ) +
    geom_text(
      aes(x = label_x, label = round(Value, 2)),
      hjust = -0.25,
      vjust = 0.5,
      size = 2.6,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = x_max * 0.80,
      y = max(plot_data$y_position, na.rm = TRUE) - 0.4,
      label = p_label,
      color = "grey15",
      fontface = "bold",
      size = 3.4,
      hjust = 0.5,
      vjust = 0.5
    ) +
    annotate(
      "label",
      x = x_min + (x_max - x_min) * 0.02,
      y = max(plot_data$y_position, na.rm = TRUE) + 0.78,
      label = group_note,
      hjust = 0,
      vjust = 0.5,
      size = 2.75,
      fontface = "bold",
      color = group_note_colour,
      fill = "white",
      linewidth = 0.22
    ) +
    scale_fill_manual(values = response_colours, name = "T cell rescue") +
    x_scale +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(add = 0.4)
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = y_limits, clip = "off") +
    labs(
      tag = panel_tag,
      x = paste0(
        title_label,
        " ratio to control",
        if_else(isTRUE(use_pseudolog), " (pseudo-log scale)", " (linear scale)")
      ),
      y = if_else(
        group_to_show == "Low Proliferators",
        "Low proliferators from Panel B",
        "High proliferators from Panel B"
      )
    ) +
    theme_fig3(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "grey20"),
      axis.text.y = element_text(size = 8, color = "grey20"),
      legend.position = if (show_legend) "top" else "none"
    )
}

pC_low <- drug_bar_plot(df, "Entospletinib", "C", "Entospletinib", "Low Proliferators", show_legend = TRUE)
pD_low <- drug_bar_plot(df, "aPD1", "D", "aPD1", "Low Proliferators")
pC_high <- drug_bar_plot(df, "Entospletinib", "C", "Entospletinib", "High Proliferators", show_legend = TRUE)
pD_high <- drug_bar_plot(df, "aPD1", "D", "aPD1", "High Proliferators")
pC_low_linear <- drug_bar_plot(df, "Entospletinib", "C", "Entospletinib", "Low Proliferators", show_legend = TRUE, use_pseudolog = FALSE)
pD_low_linear <- drug_bar_plot(df, "aPD1", "D", "aPD1", "Low Proliferators", use_pseudolog = FALSE)

figure3_main <- pA_stacked / (pB | pC_low | pD_low) +
  plot_layout(heights = c(1.35, 1), guides = "keep") &
  theme(plot.margin = margin(5, 8, 5, 8))

figure3_supp <- pB | pC_high | pD_high +
  plot_layout(guides = "keep") +
  plot_annotation(title = "Supplement: high proliferator drug responses from Panel B") &
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.margin = margin(5, 8, 5, 8)
  )

figure3_main_linear_cd <- pA_stacked / (pB | pC_low_linear | pD_low_linear) +
  plot_layout(heights = c(1.35, 1), guides = "keep") +
  plot_annotation(title = "Linear-scale comparison for Panels C-D") &
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.margin = margin(5, 8, 5, 8)
  )

save_plot_pair <- function(plot, stem, width = 16, height = 11, dpi = 400) {
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", device = grDevices::pdf)
  ggsave(png_file, plot = plot, width = width, height = height, units = "in", dpi = dpi)
  message("Saved: ", normalizePath(pdf_file))
  message("Saved: ", normalizePath(png_file))
}

save_plot_pair(figure3_main, "figure3_v5_panelA_stacked", width = 16, height = 11)
save_plot_pair(figure3_supp, "figure3_v5_supp_high_proliferators", width = 16, height = 6.5)
save_plot_pair(figure3_main_linear_cd, "figure3_v5_panelA_stacked_cd_linear", width = 16, height = 11)
save_plot_pair(pA_stacked, "figure3_v5_panelA_stacked_only", width = 14, height = 7.5)

message("Data source created date from macOS stat: ", data_created[1])
message("Data source modified time: ", as.character(data_modified))

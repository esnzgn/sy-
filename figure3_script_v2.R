# Figure 3 plotting script - revised aesthetic version
# Reads: data/MS Figure data_5.4.2026.xlsx, sheet: Fig_4_Ehsan
# Output: results/figures/Figure3_panels.pdf
#
# Requested updates implemented:
#   - Panel A: patient dots are nearly aligned with minimal jitter.
#   - Panel A: drugs are ordered by median signal, highest on the left and lowest on the right.
#   - Panel A: dot fill uses signal value, with zero/middle values shown in a visible gold color.
#   - Panel B: patient 19-00154 baseline value is capped for plotting only, so the 84.97 value
#              does not flatten the rest of the patients.
#   - Panel B: patients are ranked high-to-low left-to-right to match Panel A.
#   - Output is PDF only.

# -----------------------------
# 0. Packages
# -----------------------------
required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "stringr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the missing package(s) first: ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
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

# -----------------------------
# 1. User settings
# -----------------------------
input_file <- file.path("data", "MS Figure data_5.4.2026.xlsx")
sheet_name <- "Fig_4_Ehsan"
output_dir <- file.path("results", "figures")

# If the script is run while the Excel file is in the working directory instead of data/,
# this fallback keeps the script usable.
if (!file.exists(input_file) && file.exists("MS Figure data_5.4.2026.xlsx")) {
  input_file <- "MS Figure data_5.4.2026.xlsx"
}

# Pseudo-log helps keep low and high signals visible in the same panel while preserving zero.
# Set to FALSE if you want fully linear y-axes.
use_pseudolog_y <- TRUE

# The baseline outlier is capped for Panel B plotting only, not changed in the raw data.
baseline_cut_patient <- "19-00154"

# Final figure dimensions
fig_width <- 16
fig_height <- 10

# Colors for the patient groups defined from baseline aCD3/msIgG ranking.
group_colors <- c(
  "Baseline-low"  = "#2166AC",
  "Baseline-high" = "#B2182B"
)

treatment_colors <- c(
  "Entospletinib" = "#00843D",
  "aPD1" = "#C0006F"
)

# Value scale for Panel A. The midpoint/zero color is intentionally not white,
# so zero-valued points remain visible against the white background.
panel_a_value_colors <- c(
  low  = "#2C7BB6",
  mid  = "#FFD92F",
  high = "#D7191C"
)

# Minimal jitter for Panel A. Increase slightly if points overlap too strongly.
panel_a_jitter_width <- 0.035

# -----------------------------
# 2. Read and clean data
# -----------------------------
if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPut 'MS Figure data_5.4.2026.xlsx' inside the project's data/ folder, ",
    "or edit input_file in this script.",
    call. = FALSE
  )
}

df_raw <- read_excel(input_file, sheet = sheet_name)

# Trim accidental leading/trailing spaces in column names.
names(df_raw) <- trimws(names(df_raw))

sample_col <- names(df_raw)[str_detect(names(df_raw), regex("^sample\\s*ID", ignore_case = TRUE))][1]
baseline_col <- "aCD3/msIgG"

if (is.na(sample_col)) {
  stop("Could not find the sample ID column. Expected a name starting with 'sample ID'.", call. = FALSE)
}
if (!baseline_col %in% names(df_raw)) {
  stop("Could not find the baseline column: ", baseline_col, call. = FALSE)
}

# Drug columns are all treatment columns except sample ID and baseline stimulation.
drug_cols <- setdiff(names(df_raw), c(sample_col, baseline_col))

fig3_df <- df_raw %>%
  rename(
    patient_id = all_of(sample_col),
    baseline_acd3 = all_of(baseline_col)
  ) %>%
  mutate(
    patient_id = as.character(patient_id),
    across(all_of(c(drug_cols, "baseline_acd3")), ~ suppressWarnings(as.numeric(.x)))
  )

# -----------------------------
# 3. Define baseline-based patient stratification
# -----------------------------
# Ranking direction requested:
#   high values on the left, low values on the right.
baseline_ranked <- fig3_df %>%
  filter(!is.na(baseline_acd3)) %>%
  arrange(desc(baseline_acd3), patient_id) %>%
  mutate(
    baseline_rank = row_number(),
    baseline_group = if_else(
      baseline_rank <= n() / 2,
      "Baseline-high",
      "Baseline-low"
    ),
    baseline_group = factor(baseline_group, levels = names(group_colors)),
    patient_id_ranked = factor(patient_id, levels = patient_id)
  )

# Divider goes between the lower and upper half of ranked patients.
divider_x <- nrow(baseline_ranked) / 2 + 0.5

# Cap the outlier in Panel B only.
# Here the cap is the maximum baseline value among the other patients, rounded up slightly.
baseline_plot_cap <- baseline_ranked %>%
  filter(patient_id != baseline_cut_patient) %>%
  summarise(cap = max(baseline_acd3, na.rm = TRUE), .groups = "drop") %>%
  pull(cap)

if (!is.finite(baseline_plot_cap)) {
  baseline_plot_cap <- max(baseline_ranked$baseline_acd3, na.rm = TRUE)
}
baseline_plot_cap <- ceiling(baseline_plot_cap * 10) / 10

baseline_ranked <- baseline_ranked %>%
  mutate(
    baseline_rank = as.numeric(baseline_rank),
    baseline_acd3_display = pmin(baseline_acd3, baseline_plot_cap),
    baseline_was_cut = baseline_acd3 > baseline_plot_cap
  )

fig3_df <- fig3_df %>%
  left_join(
    baseline_ranked %>% select(patient_id, baseline_rank, baseline_group),
    by = "patient_id"
  ) %>%
  mutate(baseline_group = factor(baseline_group, levels = names(group_colors)))

panel_b_treatments <- fig3_df %>%
  select(patient_id, baseline_rank, baseline_group, all_of(names(treatment_colors))) %>%
  pivot_longer(
    cols = all_of(names(treatment_colors)),
    names_to = "treatment",
    values_to = "signal"
  ) %>%
  filter(!is.na(signal), !is.na(baseline_rank)) %>%
  mutate(
    treatment = factor(treatment, levels = names(treatment_colors)),
    x_position = baseline_rank + if_else(treatment == "Entospletinib", -0.17, 0.17)
  )

treatment_axis_max <- max(panel_b_treatments$signal, na.rm = TRUE)
if (!is.finite(treatment_axis_max) || treatment_axis_max <= 0) {
  treatment_axis_max <- baseline_plot_cap
}
treatment_axis_max <- ceiling(treatment_axis_max)
treatment_to_baseline_scale <- baseline_plot_cap / treatment_axis_max

panel_b_treatments <- panel_b_treatments %>%
  mutate(signal_display = signal * treatment_to_baseline_scale)

# -----------------------------
# 4. Shared plotting helpers
# -----------------------------
theme_figure3 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0, color = "grey10"),
      plot.subtitle = element_text(size = base_size, color = "grey30", margin = margin(b = 6)),
      axis.title = element_text(face = "bold", color = "grey15"),
      axis.text = element_text(color = "grey20"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      legend.box = "horizontal",
      plot.margin = margin(10, 14, 10, 10)
    )
}

maybe_pseudolog_y <- function(breaks = c(0, 0.1, 0.3, 1, 3, 10, 30, 100)) {
  if (isTRUE(use_pseudolog_y)) {
    scale_y_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = breaks,
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.03, 0.12))
    )
  } else {
    scale_y_continuous(
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.03, 0.12))
    )
  }
}

lollipop_drug_track <- function(data, drug_name, panel_title, panel_tag) {
  if (!drug_name %in% names(data)) {
    stop("Column not found: ", drug_name, call. = FALSE)
  }

  plot_data <- data %>%
    select(patient_id, baseline_group, signal = all_of(drug_name)) %>%
    filter(!is.na(signal)) %>%
    mutate(patient_id_ranked = factor(patient_id, levels = baseline_ranked$patient_id)) %>%
    arrange(patient_id_ranked)

  signal_median <- median(plot_data$signal, na.rm = TRUE)

  ggplot(plot_data, aes(x = patient_id_ranked, y = signal)) +
    annotate(
      "rect",
      xmin = 0.5,
      xmax = divider_x,
      ymin = -Inf,
      ymax = Inf,
      fill = group_colors[["Baseline-low"]],
      alpha = 0.045
    ) +
    annotate(
      "rect",
      xmin = divider_x,
      xmax = nrow(baseline_ranked) + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = group_colors[["Baseline-high"]],
      alpha = 0.045
    ) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey75") +
    geom_hline(yintercept = signal_median, linetype = "dashed", linewidth = 0.45, color = "grey35") +
    geom_segment(
      aes(
        xend = patient_id_ranked,
        y = 0,
        yend = signal,
        color = baseline_group
      ),
      linewidth = 0.75,
      alpha = 0.80,
      lineend = "round"
    ) +
    geom_point(
      aes(fill = baseline_group),
      shape = 21,
      color = "grey10",
      stroke = 0.30,
      size = 3.0
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE, name = "Baseline group") +
    scale_color_manual(values = group_colors, drop = FALSE, guide = "none") +
    maybe_pseudolog_y() +
    labs(
      tag = panel_tag,
      title = panel_title,
      subtitle = "Same patient order as Panel B; dashed line shows median drug signal",
      x = NULL,
      y = "Signal amount"
    ) +
    theme_figure3(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(7, 14, 5, 10)
    )
}

# -----------------------------
# 5. Panel A: all drugs, patient dots on each drug
# -----------------------------
panel_a_data <- fig3_df %>%
  select(patient_id, baseline_group, all_of(drug_cols)) %>%
  pivot_longer(
    cols = all_of(drug_cols),
    names_to = "drug",
    values_to = "signal"
  ) %>%
  filter(!is.na(signal))

# Order drugs by their median patient observation: highest left, lowest right.
drug_order <- panel_a_data %>%
  group_by(drug) %>%
  summarise(
    drug_median = median(signal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(drug_median), drug) %>%
  pull(drug)

panel_a_data <- panel_a_data %>%
  mutate(drug = factor(drug, levels = drug_order))

pA <- ggplot(panel_a_data, aes(x = drug, y = signal)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.45, color = "grey35") +
  geom_boxplot(
    width = 0.34,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey55",
    linewidth = 0.32,
    alpha = 0.75
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 2.7,
    fill = "black",
    color = "white",
    stroke = 0.25
  ) +
  geom_point(
    aes(fill = signal),
    shape = 21,
    color = "grey12",
    stroke = 0.18,
    size = 2.45,
    alpha = 0.96,
    position = position_jitter(width = panel_a_jitter_width, height = 0, seed = 123)
  ) +
  scale_fill_gradient2(
    low = panel_a_value_colors[["low"]],
    mid = panel_a_value_colors[["mid"]],
    high = panel_a_value_colors[["high"]],
    midpoint = 0,
    name = "Signal\nvalue",
    labels = label_number(accuracy = 0.01)
  ) +
  maybe_pseudolog_y() +
  labs(
    tag = "A",
    title = "Drug response landscape across patients",
    subtitle = "Dots are nearly aligned per drug; median signal ranks high-to-low from left to right",
    x = "Drug / stimulation condition ordered by median response",
    y = "Signal amount"
  ) +
  theme_figure3(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    legend.position = "top"
  )

# -----------------------------
# 6. Panel B: baseline aCD3/msIgG ranking and stratification divider
# -----------------------------
pB <- ggplot(baseline_ranked, aes(x = baseline_rank, y = baseline_acd3_display)) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = divider_x,
    ymin = -Inf,
    ymax = Inf,
    fill = group_colors[["Baseline-high"]],
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = divider_x,
    xmax = nrow(baseline_ranked) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = group_colors[["Baseline-low"]],
    alpha = 0.045
  ) +
  geom_hline(yintercept = median(baseline_ranked$baseline_acd3_display, na.rm = TRUE),
             linetype = "dotted", linewidth = 0.45, color = "grey35") +
  geom_col(
    width = 0.82,
    fill = "grey18",
    color = NA,
    alpha = 0.08
  ) +
  geom_col(
    aes(fill = baseline_group),
    width = 0.68,
    color = "white",
    linewidth = 0.45,
    alpha = 0.50
  ) +
  geom_col(
    aes(fill = baseline_group),
    width = 0.36,
    color = NA,
    alpha = 0.12
  ) +
  geom_point(
    aes(fill = baseline_group),
    shape = 21,
    color = "white",
    stroke = 0.40,
    size = 2.8
  ) +
  geom_segment(
    data = panel_b_treatments,
    aes(
      x = x_position,
      xend = x_position,
      y = 0,
      yend = signal_display
    ),
    inherit.aes = FALSE,
    linewidth = 2.8,
    color = "white",
    alpha = 0.95,
    lineend = "round"
  ) +
  geom_segment(
    data = panel_b_treatments,
    aes(
      x = x_position,
      xend = x_position,
      y = 0,
      yend = signal_display,
      color = treatment
    ),
    inherit.aes = FALSE,
    linewidth = 1.75,
    alpha = 1,
    lineend = "round"
  ) +
  geom_point(
    data = panel_b_treatments,
    aes(
      x = x_position,
      y = signal_display
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "white",
    stroke = 0,
    size = 5.8
  ) +
  geom_point(
    data = panel_b_treatments,
    aes(
      x = x_position,
      y = signal_display,
      color = treatment,
      shape = treatment
    ),
    inherit.aes = FALSE,
    stroke = 1.1,
    size = 4.3,
    alpha = 1
  ) +
  geom_vline(xintercept = divider_x, linetype = "dashed", linewidth = 0.55, color = "grey20") +
  geom_point(
    data = baseline_ranked %>% filter(baseline_was_cut),
    aes(x = baseline_rank, y = baseline_acd3_display),
    inherit.aes = FALSE,
    shape = 24,
    size = 3.7,
    fill = "white",
    color = "black",
    stroke = 0.45
  ) +
  geom_text(
    data = baseline_ranked %>% filter(baseline_was_cut),
    aes(
      x = baseline_rank,
      y = baseline_acd3_display,
      label = paste0("cut from ", round(baseline_acd3, 1))
    ),
    inherit.aes = FALSE,
    vjust = -0.80,
    size = 3.1,
    fontface = "bold",
    color = "grey15"
  ) +
  annotate(
    "text",
    x = divider_x,
    y = baseline_plot_cap * 0.12,
    label = "median split",
    angle = 90,
    size = 3.1,
    color = "grey25",
    fontface = "bold",
    vjust = -0.7
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = "Baseline group") +
  scale_color_manual(values = treatment_colors, drop = FALSE, name = "Treatment overlay") +
  scale_shape_manual(values = c("Entospletinib" = 16, "aPD1" = 17), drop = FALSE, name = "Treatment overlay") +
  scale_x_continuous(
    breaks = baseline_ranked$baseline_rank,
    labels = baseline_ranked$patient_id,
    expand = expansion(add = 0.5)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.20)),
    sec.axis = sec_axis(
      ~ . / treatment_to_baseline_scale,
      name = "Treatment signal overlay",
      labels = label_number(accuracy = 0.1)
    )
  ) +
  coord_cartesian(ylim = c(0, baseline_plot_cap * 1.38), clip = "off") +
  labs(
    tag = "B",
    title = "Baseline stimulation with Entospletinib and aPD1 overlays",
    subtitle = paste0(
      "Bars show baseline ranked high-to-low; colored lollipops show treatments on the right axis; ",
      baseline_cut_patient,
      " is capped at ",
      baseline_plot_cap,
      " for plotting only"
    ),
    x = "Patient ranked by baseline stimulation, high to low",
    y = "Baseline stimulation (aCD3/msIgG)"
  ) +
  theme_figure3(base_size = 10) +
  theme(
    axis.title.y.right = element_text(face = "bold", color = "grey15"),
    legend.position = "top"
  )

# -----------------------------
# 8. Combine and save Figure 3 as PDF only
# -----------------------------
figure3 <- pA / pB +
  plot_layout(heights = c(1.25, 1.15), guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, color = "grey10"),
      legend.position = "top"
    )
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(output_dir, "Figure3_panels.pdf")

ggsave(
  filename = pdf_file,
  plot = figure3,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = grDevices::pdf
)

message("Saved Figure 3 PDF to: ", normalizePath(pdf_file))

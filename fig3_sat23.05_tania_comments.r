# Figure 3 - Tania comments, Sat 23.05
# Reads: data/MS Figure data_5.4.2026.xlsx, sheet: Fig_4_Ehsan
# Outputs:
#   results/figures/fig3_sat23.05_tania_comments.pdf
#   results/figures/fig3_sat23.05_tania_comments.png

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "scales",
  "stringr", "ggbeeswarm"
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

input_file <- file.path("data", "MS Figure data_5.4.2026.xlsx")
sheet_name <- "Fig_4_Ehsan"
output_dir <- file.path("results", "figures")

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

# Green encodes proliferation. High proliferators use saturated green bars;
# low proliferators use light green bars plus a speckled background so they stand out.
group_colours <- c("High Proliferators" = "#00843D", "Low Proliferators" = "#BFE6B8")
background_colours <- c("High Proliferators" = "#D8F1D4", "Low Proliferators" = "#FFFFFF")
speckle_colour <- "#62B45D"

mid_val <- median(df$baseline, na.rm = TRUE)
y_cap <- 5

baseline_ranked <- df |>
  select(Patient, baseline) |>
  filter(!is.na(baseline)) |>
  arrange(desc(baseline), Patient) |>
  mutate(
    group = if_else(baseline >= mid_val, "High Proliferators", "Low Proliferators"),
    group = factor(group, levels = names(group_colours)),
    baseline_rank = row_number(),
    y_position = n() - baseline_rank + 1,
    baseline_disp = pmin(baseline, y_cap),
    is_outlier = baseline > y_cap,
    Patient = factor(Patient, levels = Patient)
  )

n_high <- sum(baseline_ranked$group == "High Proliferators")
n_tot <- nrow(baseline_ranked)
divider_y <- n_tot - n_high + 0.5
patient_order <- as.character(baseline_ranked$Patient)
outlier_labels <- baseline_ranked |> filter(is_outlier)

patient_groups <- baseline_ranked |>
  transmute(Patient = as.character(Patient), group, baseline_rank, y_position)

make_low_speckles <- function(xmin, xmax, n = 220) {
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin || n_high >= n_tot) {
    return(tibble(x = numeric(), y = numeric()))
  }
  set.seed(2305)
  tibble(
    x = runif(n, xmin, xmax),
    y = runif(n, 0.5, divider_y)
  )
}

group_background <- function(xmin, xmax, high_alpha = 0.55, low_speckles = 220) {
  low_speckles_df <- make_low_speckles(xmin, xmax, low_speckles)

  list(
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = divider_y,
      ymax = n_tot + 0.5,
      fill = background_colours[["High Proliferators"]],
      alpha = high_alpha
    ),
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = 0.5,
      ymax = divider_y,
      fill = background_colours[["Low Proliferators"]],
      alpha = 0.92
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

# Panel A: all drugs, ranked by median response.
drug_medians <- df |>
  select(Patient, all_of(drug_cols)) |>
  pivot_longer(-Patient, names_to = "Drug", values_to = "Value") |>
  filter(!is.na(Value), Value > 0) |>
  summarise(med = median(Value), .by = Drug) |>
  arrange(desc(med), Drug)

panel_a_data <- df |>
  select(Patient, all_of(drug_cols)) |>
  pivot_longer(-Patient, names_to = "Drug", values_to = "Value") |>
  filter(!is.na(Value), Value > 0) |>
  mutate(Drug = factor(Drug, levels = drug_medians$Drug))

pA <- ggplot(panel_a_data, aes(x = Drug, y = Value)) +
  geom_rect(
    data = drug_medians |>
      mutate(Drug = factor(Drug, levels = drug_medians$Drug), xi = as.integer(Drug)) |>
      filter(xi %% 2 == 0),
    aes(xmin = xi - 0.5, xmax = xi + 0.5, ymin = -Inf, ymax = Inf),
    fill = "grey96",
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.55) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    color = "#1A1A2E",
    fill = NA,
    linewidth = 0.55
  ) +
  geom_beeswarm(
    aes(fill = log2(Value)),
    cex = 2.2,
    size = 2.8,
    shape = 21,
    stroke = 0.35,
    color = "white",
    alpha = 0.88
  ) +
  scale_fill_gradient2(
    low = "#3182BD",
    mid = "#E8A838",
    high = "#CB181D",
    midpoint = 0,
    name = "log2 ratio",
    guide = guide_colorbar(barwidth = 0.5, barheight = 5, ticks.colour = "grey40")
  ) +
  scale_y_log10(
    breaks = c(0.1, 0.3, 0.5, 1, 2, 5, 10, 25),
    labels = c("0.1", "0.3", "0.5", "1", "2", "5", "10", "25"),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(
    tag = "A",
    x = NULL,
    y = "Ratio to control (log10 scale)"
  ) +
  theme_fig3(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9, color = "grey20"),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    legend.position = "right"
  )

# Panel B: keep the waterfall structure, rotated clockwise into horizontal bars.
pB <- ggplot(baseline_ranked, aes(y = y_position, fill = group)) +
  group_background(xmin = 0, xmax = y_cap, high_alpha = 0.58, low_speckles = 260) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = baseline_disp,
      ymin = y_position - 0.35,
      ymax = y_position + 0.35
    ),
    color = "white",
    linewidth = 0.4,
    alpha = 0.95
  ) +
  geom_vline(xintercept = mid_val, linetype = "dashed", color = "#333333", linewidth = 0.85) +
  annotate(
    "label",
    x = mid_val,
    y = 1.1,
    label = paste0("Median = ", round(mid_val, 2)),
    hjust = -0.1,
    vjust = 0.5,
    size = 3.0,
    color = "#333333",
    fill = "white",
    label.padding = unit(0.15, "lines")
  ) +
  geom_text(
    data = ~ filter(.x, !is_outlier),
    aes(x = baseline_disp, label = round(baseline_disp, 2)),
    hjust = -0.25,
    vjust = 0.5,
    size = 2.6,
    color = "grey30"
  ) +
  annotate(
    "segment",
    x = y_cap * 0.82,
    xend = y_cap * 0.88,
    y = n_tot - 0.35,
    yend = n_tot + 0.35,
    color = "white",
    linewidth = 1.8
  ) +
  annotate(
    "segment",
    x = y_cap * 0.77,
    xend = y_cap * 0.83,
    y = n_tot - 0.35,
    yend = n_tot + 0.35,
    color = "white",
    linewidth = 1.8
  ) +
  geom_text(
    data = outlier_labels,
    aes(x = baseline_disp, y = y_position, label = paste0("^ ", round(baseline, 1))),
    vjust = 0.5,
    hjust = -0.15,
    size = 3.0,
    color = "#333333",
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = y_cap * 0.88,
    y = n_tot - n_high / 2 + 0.5,
    label = "High",
    color = group_colours[["High Proliferators"]],
    fontface = "bold",
    size = 3.5,
    hjust = 0.5,
    vjust = 0.5
  ) +
  annotate(
    "text",
    x = y_cap * 0.88,
    y = (n_tot - n_high) / 2 + 0.5,
    label = "Low",
    color = "#2F8F35",
    fontface = "bold",
    size = 3.5,
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_fill_manual(values = group_colours, name = "Baseline stimulation") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_y_continuous(
    breaks = baseline_ranked$y_position,
    labels = as.character(baseline_ranked$Patient),
    expand = expansion(add = 0.4)
  ) +
  coord_cartesian(xlim = c(0, y_cap), clip = "off") +
  labs(
    tag = "B",
    x = "Baseline stimulation (aCD3/msIgG)",
    y = "Patient (ranked high to low in Panel B)"
  ) +
  theme_fig3(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "grey20"),
    axis.text.y = element_text(size = 8, color = "grey20"),
    legend.position = "top"
  )

lollipop_plot <- function(data, drug_col, panel_tag, title_label) {
  if (!drug_col %in% names(data)) {
    stop("Column not found: ", drug_col, call. = FALSE)
  }

  plot_data <- data |>
    select(Patient, Value = all_of(drug_col)) |>
    filter(!is.na(Value)) |>
    left_join(patient_groups, by = "Patient") |>
    mutate(Patient = factor(Patient, levels = patient_order)) |>
    arrange(Patient)

  x_min <- min(0, min(plot_data$Value, na.rm = TRUE)) * 1.08
  x_max <- max(plot_data$Value, na.rm = TRUE) * 1.08

  ggplot(plot_data, aes(y = y_position, x = Value)) +
    group_background(xmin = x_min, xmax = x_max, high_alpha = 0.50, low_speckles = 220) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    geom_segment(
      aes(x = 0, xend = Value, y = y_position, yend = y_position),
      linewidth = 1.05,
      alpha = 0.78,
      color = "black"
    ) +
    geom_point(size = 4.4, alpha = 0.96, color = "black") +
    geom_point(size = 2.0, color = "white", alpha = 0.9) +
    geom_text(
      aes(label = round(Value, 2)),
      hjust = -0.35,
      vjust = 0.5,
      size = 2.6,
      color = "grey25"
    ) +
    annotate(
      "text",
      x = x_max * 0.80,
      y = n_tot - n_high / 2 + 0.5,
      label = "High",
      color = group_colours[["High Proliferators"]],
      fontface = "bold",
      size = 3.2,
      hjust = 0.5,
      vjust = 0.5
    ) +
    annotate(
      "text",
      x = x_max * 0.80,
      y = (n_tot - n_high) / 2 + 0.5,
      label = "Low",
      color = "#2F8F35",
      fontface = "bold",
      size = 3.2,
      hjust = 0.5,
      vjust = 0.5
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    scale_y_continuous(
      breaks = baseline_ranked$y_position,
      labels = as.character(baseline_ranked$Patient),
      expand = expansion(add = 0.4)
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(
      tag = panel_tag,
      x = paste0(title_label, " ratio to control"),
      y = "Patient (ranked high to low in Panel B)"
    ) +
    theme_fig3(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "grey20"),
      axis.text.y = element_text(size = 8, color = "grey20")
    )
}

pC <- lollipop_plot(df, "Entospletinib", "C", "Entospletinib")
pD <- lollipop_plot(df, "aPD1", "D", "aPD1")

# Bottom row order requested: C, B, D.
fig <- pA / (pC | pB | pD) +
  plot_layout(heights = c(1.4, 1), guides = "keep") &
  theme(plot.margin = margin(5, 8, 5, 8))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(output_dir, "fig3_sat23.05_tania_comments.pdf")
png_file <- file.path(output_dir, "fig3_sat23.05_tania_comments.png")

ggsave(
  filename = pdf_file,
  plot = fig,
  width = 16,
  height = 11,
  device = grDevices::pdf
)

ggsave(
  filename = png_file,
  plot = fig,
  width = 16,
  height = 11,
  units = "in",
  dpi = 400
)

message("Saved Figure 3 PDF to: ", normalizePath(pdf_file))
message("Saved Figure 3 PNG to: ", normalizePath(png_file))

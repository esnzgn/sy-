# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# options(repos = BiocManager::repositories())
# 
# BiocManager::install("ComplexHeatmap", ask = FALSE, update = TRUE)

library(readxl)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)
library(openxlsx)
library(janitor)

file_cell_death <- "./data/Z_patient_AMLCelldeath(SMi_div_aCD3)and(SMi+aPD1_div_aCD3).xlsx"
file_tcell_prolif <- "./data/Z_patient_TcellProliferation(SMi_div_aCD3)and(SMi+aPD1_div_aCD3).xlsx"
file_genetics <- "./data/AML v2_HEATMAP_GENETICS_04.27.2026.xlsx"

read_heatmap_matrix <- function(file_path) {
  sh <- excel_sheets(file_path)[1]
  
  df <- read_excel(file_path, sheet = sh)
  
  mat <- df %>%
    as.data.frame() %>%
    rename(drug = 1) %>%
    mutate(drug = as.character(drug))
  
  rownames(mat) <- mat$drug
  mat$drug <- NULL
  
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  
  return(mat)
}

read_npm1_metadata <- function(file_path) {
  df <- read_excel(file_path, sheet = 1) %>%
    janitor::clean_names()
  
  df %>%
    transmute(
      patient_id = as.character(lab_specimen_id),
      npm1_status = case_when(
        npm1 == 100 ~ "NPM1 mutated",
        npm1 == 1   ~ "NPM1 wildtype",
        npm1 == 0   ~ "Unknown/untested",
        TRUE        ~ "Unknown/untested"
      )
    )
}

make_heatmap_drug_patient_waterfall <- function(mat,
                                                npm1_meta,
                                                title = "Heatmap",
                                                row_score_fun = function(x) median(x, na.rm = TRUE),
                                                col_score_fun = function(x) median(x, na.rm = TRUE),
                                                na_col = "grey90") {
  
  # Drug score = integrated drug effect across patients
  row_scores <- apply(mat, 1, row_score_fun)
  row_order <- order(row_scores, decreasing = TRUE)
  
  # Patient score = integrated patient response across drugs
  col_scores <- apply(mat, 2, col_score_fun)
  
  # Rank patients by median response, like drugs
  col_order <- order(col_scores, decreasing = TRUE)
  
  mat_ord <- mat[row_order, col_order, drop = FALSE]
  row_scores_ord <- row_scores[row_order]
  col_scores_ord <- col_scores[col_order]
  
  patient_ids_ord <- colnames(mat_ord)
  
  npm1_status_ord <- npm1_meta %>%
    filter(patient_id %in% patient_ids_ord) %>%
    right_join(
      tibble(patient_id = patient_ids_ord),
      by = "patient_id"
    ) %>%
    mutate(
      npm1_status = ifelse(is.na(npm1_status), "Unknown/untested", npm1_status),
      npm1_status = factor(
        npm1_status,
        levels = c("NPM1 mutated", "NPM1 wildtype", "Unknown/untested")
      )
    ) %>%
    pull(npm1_status)
  
  # Separate scales:
  # drug waterfall uses drug range
  # patient waterfall uses patient range, so patient bars are not too short
  # Use actual drug median range, so drug waterfall is scaled to drug effects
  # Symmetric drug waterfall scale around zero
  drug_max_abs <- max(abs(row_scores_ord), na.rm = TRUE)
  if (drug_max_abs == 0) drug_max_abs <- 1
  
  drug_pad <- 0.05 * drug_max_abs
  drug_ylim <- c(
    -drug_max_abs - drug_pad,
    drug_max_abs + drug_pad
  )
  
  # Symmetric patient waterfall scale around zero
  patient_max_abs <- max(abs(col_scores_ord), na.rm = TRUE)
  if (patient_max_abs == 0) patient_max_abs <- 1
  
  patient_pad <- 0.05 * patient_max_abs
  patient_ylim <- c(
    -patient_max_abs - patient_pad,
    patient_max_abs + patient_pad
  )
  
  # Heatmap colors
  col_fun <- colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#2c7fb8", "#7fcdbb", "#c7e9b4", "#fdae61", "#d7191c")
  )
  
  # Drug waterfall colors by median z-score
  drug_bar_col_fun <- colorRamp2(
    c(drug_min, 0, drug_max),
    c("#2c7fb8", "grey80", "#d7191c")
  )
  
  # Patient waterfall colors by NPM1 mutation status
  npm1_cols <- c(
    "NPM1 mutated" = "#17b3b7",
    "NPM1 wildtype" = "grey70",
    "Unknown/untested" = "white"
  )
  #  only the plotted patient bars so positive values appear above baseline
  col_scores_plot <- col_scores_ord
  
  patient_max_abs_score <- max(abs(col_scores_ord), na.rm = TRUE)
  if (patient_max_abs_score == 0) patient_max_abs_score <- 1
  
  right_ha <- rowAnnotation(
    `Drug median` = anno_barplot(
      row_scores_ord,
      gp = gpar(fill = drug_bar_col_fun(row_scores_ord), col = NA),
      border = FALSE,
      width = unit(3.5, "cm"),
      ylim = drug_ylim,
      axis_param = list(
        side = "bottom",
        at = c(-drug_max_abs, 0, drug_max_abs),
        labels = round(c(-drug_max_abs, 0, drug_max_abs), 2),
        labels_rot = 0
      ),
      bar_width = 0.85
    ),
    annotation_name_rot = 90,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  #  only the plotted patient bars so positive values appear above baseline
  col_scores_plot <- col_scores_ord
  
  # patient_min <- min(col_scores_ord, na.rm = TRUE)
  # patient_max <- max(col_scores_ord, na.rm = TRUE)
  # 
  # patient_pad <- 0.05 * (patient_max - patient_min)
  # 
  # if (patient_pad == 0) patient_pad <- 0.05
  # 
  # patient_ylim <- c(patient_min - patient_pad, patient_max + patient_pad)
  
  bottom_ha <- HeatmapAnnotation(
    `Patient median` = anno_barplot(
      col_scores_plot,
      gp = gpar(
        fill = npm1_cols[as.character(npm1_status_ord)],
        col = "black",
        lwd = 0.2
      ),
      border = FALSE,
      height = unit(2.2, "cm"),
      
      # use actual patient bar range
      ylim = patient_ylim,
      
      axis_param = list(
        side = "left",
        at = c(-patient_max_abs, 0, patient_max_abs),
        labels = round(c(-patient_max_abs, 0, patient_max_abs), 2),
        labels_rot = 0
      ),
      
      bar_width = 0.85
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    which = "column"
  )
  
  ht <- Heatmap(
    mat_ord,
    name = "Z-Score",
    col = col_fun,
    na_col = na_col,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat_ord[i, j]
      if (!is.na(val) && val > 0) {
        grid.text(
          sprintf("%.1f", val),
          x, y,
          gp = gpar(fontsize = 5.5, col = "black")
        )
      }
    },
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 7),
    
    column_names_side = "top",
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 7),
    
    right_annotation = right_ha,
    bottom_annotation = bottom_ha,
    
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "Z-Score",
      legend_height = unit(4, "cm")
    ),
    
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  npm1_legend <- Legend(
    title = "NPM1",
    labels = names(npm1_cols),
    legend_gp = gpar(fill = npm1_cols)
  )
  
  return(list(
    heatmap = ht,
    npm1_legend = npm1_legend,
    matrix_ordered = mat_ord,
    row_scores = row_scores_ord,
    col_scores = col_scores_ord,
    npm1_status = npm1_status_ord
  ))
}

mat_cell_death <- read_heatmap_matrix(file_cell_death)
mat_tcell_prolif <- read_heatmap_matrix(file_tcell_prolif)
npm1_meta <- read_npm1_metadata(file_genetics)

dim(mat_cell_death)
dim(mat_tcell_prolif)

head(mat_cell_death[, 1:5])
head(mat_tcell_prolif[, 1:5])

res_cell_death <- make_heatmap_drug_patient_waterfall(
  mat = mat_cell_death,
  npm1_meta = npm1_meta,
  title = "AML Cell Death"
)

draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_cell_death$npm1_legend)
)

pdf("AML_CellDeath_heatmap_with_drug_patient_waterfalls_NPM1.pdf", width = 10, height = 6)
draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_cell_death$npm1_legend)
)
dev.off()

png("AML_CellDeath_heatmap_with_waterfalls.png", width = 2200, height = 2600, res = 300)
draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_cell_death$npm1_legend)
)
dev.off()

res_tcell_prolif <- make_heatmap_drug_patient_waterfall(
  mat = mat_tcell_prolif,
  npm1_meta = npm1_meta,
  title = "T-cell Proliferation"
)

draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_tcell_prolif$npm1_legend)
)

pdf("Tcell_Proliferation_heatmap_with_drug_patient_waterfalls_NPM1.pdf", width = 10, height = 6)
draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_tcell_prolif$npm1_legend)
)
dev.off()

png("Tcell_Proliferation_heatmap_with_waterfalls.png", width = 2200, height = 2600, res = 300)
draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(res_tcell_prolif$npm1_legend)
)
dev.off()

write.xlsx(
  list(
    AML_CellDeath_Drugmedians = data.frame(
      Drug = rownames(res_cell_death$matrix_ordered),
      medianScore = res_cell_death$row_scores
    ),
    TcellProlif_Drugmedians = data.frame(
      Drug = rownames(res_tcell_prolif$matrix_ordered),
      medianScore = res_tcell_prolif$row_scores
    )
  ),
  file = "heatmap_drug_waterfall_summary_tables.xlsx",
  overwrite = TRUE
)


read_mutation_matrix <- function(file_path) {
  mut_df <- read_excel(file_path, sheet = 1) %>%
    janitor::clean_names()
  
  patient_ids <- mut_df$lab_specimen_id
  
  mut_mat <- mut_df %>%
    select(-lab_specimen_id) %>%
    as.data.frame()
  
  rownames(mut_mat) <- patient_ids
  
  mut_mat <- as.matrix(mut_mat)
  
  # transpose: genes as rows, patients as columns
  mut_mat <- t(mut_mat)
  
  mut_mat
}

mutation_mat <- read_mutation_matrix(file_genetics)

mutation_col_fun <- c(
  "100" = "#E64B35",  # mutation positive
  "1"   = "#4DBBD5",  # mutation negative
  "0"   = "#BDBDBD"   # unknown/untested
)

ht_mutations <- Heatmap(
  as.character(mutation_mat),
  name = "Mutation",
  col = mutation_col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 7),
  
  column_title = "AML Mutation Matrix",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  heatmap_legend_param = list(
    title = "Mutation status",
    at = c("100", "1", "0"),
    labels = c(
      "Positive",
      "Negative",
      "Unknown/untested"
    )
  ),
  
  rect_gp = gpar(col = "white", lwd = 0.5)
)

draw(ht_mutations, heatmap_legend_side = "right")

pdf("AML_mutation_matrix_all_genes.pdf", width = 10, height = 8)
draw(ht_mutations, heatmap_legend_side = "right")
dev.off()

png("AML_mutation_matrix_all_genes.png", width = 2400, height = 1800, res = 300)
draw(ht_mutations, heatmap_legend_side = "right")
dev.off()



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


file_cell_death <- "../data/Z_patient_AMLCelldeath(SMi_div_aCD3)and(SMi+aPD1_div_aCD3).xlsx"
file_tcell_prolif <- "../data/Z_patient_TcellProliferation(SMi_div_aCD3)and(SMi+aPD1_div_aCD3).xlsx"

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

make_heatmap_drug_waterfall <- function(mat,
                                        title = "Heatmap",
                                        row_score_fun = function(x) mean(x, na.rm = TRUE),
                                        na_col = "grey90") {
  
  # Row score = integrated drug effect across patients
  row_scores <- apply(mat, 1, row_score_fun)
  
  # Order rows from highest to lowest z-score summary
  row_order <- order(row_scores, decreasing = TRUE)
  
  # Keep patient order as in the input file
  col_order <- seq_len(ncol(mat))
  
  mat_ord <- mat[row_order, col_order, drop = FALSE]
  row_scores_ord <- row_scores[row_order]
  
  # Bar range
  max_abs_score <- max(abs(row_scores_ord), na.rm = TRUE)
  if (max_abs_score == 0) max_abs_score <- 1
  
  # Heatmap colors
  col_fun <- colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#2c7fb8", "#7fcdbb", "#c7e9b4", "#fdae61", "#d7191c")
  )
  
  # Drug waterfall colors
  bar_col_fun <- colorRamp2(
    c(-max_abs_score, 0, max_abs_score),
    c("#2c7fb8", "grey80", "#d7191c")
  )
  
  # Right-side drug waterfall only
  right_ha <- rowAnnotation(
    `Drug mean` = anno_barplot(
      row_scores_ord,
      gp = gpar(fill = bar_col_fun(row_scores_ord), col = NA),
      border = FALSE,
      width = unit(3.5, "cm"),
      ylim = c(-max_abs_score, max_abs_score),
      axis_param = list(
        side = "bottom",
        at = c(-max_abs_score, 0, max_abs_score),
        labels_rot = 0
      ),
      bar_width = 0.85
    ),
    annotation_name_rot = 90,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  ht <- Heatmap(
    mat_ord,
    name = "Z-Score",
    col = col_fun,
    na_col = na_col,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat_ord[i, j]
      
      if (!is.na(val) && val > 0.5) {
        grid.text(
          sprintf("%.1f", val),
          x, y,
          gp = gpar(
            fontsize = 6,
            col = ifelse(val > 1, "white", "black")
          )
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
    
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "Z-Score",
      legend_height = unit(4, "cm")
    ),
    
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  return(list(
    heatmap = ht,
    matrix_ordered = mat_ord,
    row_scores = row_scores_ord
  ))
}

mat_cell_death <- read_heatmap_matrix(file_cell_death)
mat_tcell_prolif <- read_heatmap_matrix(file_tcell_prolif)

dim(mat_cell_death)
dim(mat_tcell_prolif)

head(mat_cell_death[, 1:5])
head(mat_tcell_prolif[, 1:5])

res_cell_death <- make_heatmap_drug_waterfall(
  mat = mat_cell_death,
  title = "AML Cell Death"
)

draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

pdf("AML_CellDeath_heatmap_with_waterfalls.pdf", width = 10, height = 12)
draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

png("AML_CellDeath_heatmap_with_waterfalls.png", width = 2200, height = 2600, res = 300)
draw(
  res_cell_death$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

res_tcell_prolif <- make_heatmap_drug_waterfall(
  mat = mat_tcell_prolif,
  title = "T-cell Proliferation"
)

draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

pdf("Tcell_Proliferation_heatmap_with_waterfalls.pdf", width = 10, height = 12)
draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

png("Tcell_Proliferation_heatmap_with_waterfalls.png", width = 2200, height = 2600, res = 300)
draw(
  res_tcell_prolif$heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

write.xlsx(
  list(
    AML_CellDeath_DrugMeans = data.frame(
      Drug = rownames(res_cell_death$matrix_ordered),
      MeanScore = res_cell_death$row_scores
    ),
    TcellProlif_DrugMeans = data.frame(
      Drug = rownames(res_tcell_prolif$matrix_ordered),
      MeanScore = res_tcell_prolif$row_scores
    )
  ),
  file = "heatmap_drug_waterfall_summary_tables.xlsx",
  overwrite = TRUE
)



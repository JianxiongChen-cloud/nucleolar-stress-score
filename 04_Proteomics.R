rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Working directory and output directory
############################################################
setwd("/home/xxm_xxm/CJX_workspace/proteinomics/")
out_dir <- file.path(getwd(), "NuS_RiboSis_ssGSEA_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## Final NuStress gene set path
nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"

## RiboSis gene set path
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

############################################################
## 1. Load R packages
############################################################
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(multcomp)
library(tidyr)
library(readxl)

############################################################
## 2. Read matrix.xlsx
############################################################
expr_df <- read_xlsx("matrix.xlsx")

cat("===== Input matrix dimension =====\n")
print(dim(expr_df))
cat("===== Column names =====\n")
print(colnames(expr_df))

stopifnot(ncol(expr_df) >= 2)

############################################################
## 3. Build expression matrix
############################################################
expr_mat <- as.matrix(expr_df[, 2:ncol(expr_df), drop = FALSE])
rownames(expr_mat) <- toupper(trimws(as.character(expr_df[[1]])))
mode(expr_mat) <- "numeric"

## Remove empty gene names
valid_gene <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
expr_mat <- expr_mat[valid_gene, , drop = FALSE]

## Remove duplicates and keep the first occurrence
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

cat("===== Expression matrix dimension after cleaning =====\n")
print(dim(expr_mat))
print(colnames(expr_mat))

############################################################
## 4. Group information
############################################################
## Modify this section only if the grouping differs
group_info <- data.frame(
  Sample = colnames(expr_mat),
  Group = rep(c("Control", "OXA", "L", "L-OXA"), each = 2),
  stringsAsFactors = FALSE
)

group_info$Group <- factor(group_info$Group, levels = unique(group_info$Group))
print(group_info)

############################################################
## 5. Fill NA values with within-group means
############################################################
fill_na_by_group_mean <- function(expr_mat, group_info) {
  expr_df2 <- as.data.frame(expr_mat, check.names = FALSE)
  
  for (grp in unique(group_info$Group)) {
    grp_samples <- group_info$Sample[group_info$Group == grp]
    
    for (g in rownames(expr_df2)) {
      vals <- as.numeric(expr_df2[g, grp_samples, drop = TRUE])
      if (any(is.na(vals))) {
        grp_mean <- mean(vals, na.rm = TRUE)
        if (!is.na(grp_mean)) {
          na_idx <- which(is.na(vals))
          if (length(na_idx) > 0) {
            expr_df2[g, grp_samples[na_idx]] <- grp_mean
          }
        }
      }
    }
  }
  
  expr_out <- as.matrix(expr_df2)
  mode(expr_out) <- "numeric"
  return(expr_out)
}

expr_mat <- fill_na_by_group_mean(expr_mat, group_info)
cat("Remaining NA count after filling:", sum(is.na(expr_mat)), "\n")

## Remove genes that still contain NA values
expr_mat <- expr_mat[complete.cases(expr_mat), , drop = FALSE]

cat("===== Matrix dimension after NA removal =====\n")
print(dim(expr_mat))

############################################################
## 6. Read final NuStress gene sets
############################################################
if (!file.exists(nus_rdata)) {
  stop("NuStress gene set file not found: ", nus_rdata)
}

load(nus_rdata)

if (!exists("geneSets_final")) {
  stop("Object 'geneSets_final' not found in NuStress_geneSets_final.Rdata")
}

required_sets <- c("NuStress_UP", "NuStress_DOWN", "NuStressCore_UP", "NuStressCore_DOWN")
missing_sets <- setdiff(required_sets, names(geneSets_final))
if (length(missing_sets) > 0) {
  stop("Missing gene sets in geneSets_final: ", paste(missing_sets, collapse = ", "))
}

clean_gene_vector <- function(x) {
  x <- unique(toupper(trimws(as.character(x))))
  x <- x[!is.na(x) & x != ""]
  return(x)
}

full_up   <- clean_gene_vector(geneSets_final$NuStress_UP)
full_down <- clean_gene_vector(geneSets_final$NuStress_DOWN)
core_up   <- clean_gene_vector(geneSets_final$NuStressCore_UP)
core_down <- clean_gene_vector(geneSets_final$NuStressCore_DOWN)

cat("===== Original gene set sizes from NuStress_geneSets_final.Rdata =====\n")
cat("NuStress_UP:", length(full_up), "\n")
cat("NuStress_DOWN:", length(full_down), "\n")
cat("NuStressCore_UP:", length(core_up), "\n")
cat("NuStressCore_DOWN:", length(core_down), "\n")

## Intersect with the expression matrix
full_up_in   <- intersect(full_up, rownames(expr_mat))
full_down_in <- intersect(full_down, rownames(expr_mat))

cat("===== Full gene set overlap with expression matrix =====\n")
cat("NuStress_UP overlap:", length(full_up_in), "\n")
cat("NuStress_DOWN overlap:", length(full_down_in), "\n")

if (length(full_up_in) < 5 || length(full_down_in) < 5) {
  stop("Too few overlapping genes for Full NuStress gene set.")
}

geneSets_full <- list(
  NuStress_UP   = full_up_in,
  NuStress_DOWN = full_down_in
)

## Core
has_core <- TRUE
core_up_in   <- intersect(core_up, rownames(expr_mat))
core_down_in <- intersect(core_down, rownames(expr_mat))

cat("===== Core gene set overlap with expression matrix =====\n")
cat("NuStressCore_UP overlap:", length(core_up_in), "\n")
cat("NuStressCore_DOWN overlap:", length(core_down_in), "\n")

if (length(core_up_in) >= 5 && length(core_down_in) >= 5) {
  geneSets_core <- list(
    NuStressCore_UP   = core_up_in,
    NuStressCore_DOWN = core_down_in
  )
} else {
  warning("Too few overlapping genes for Core NuStress gene set. Core analysis skipped.")
  has_core <- FALSE
}

## Export overlap and missing gene lists for checking
write.csv(data.frame(gene_symbol = full_up_in, stringsAsFactors = FALSE),
          file.path(out_dir, "NuStress_UP_overlap_with_matrix.csv"),
          row.names = FALSE)

write.csv(data.frame(gene_symbol = full_down_in, stringsAsFactors = FALSE),
          file.path(out_dir, "NuStress_DOWN_overlap_with_matrix.csv"),
          row.names = FALSE)

write.csv(data.frame(gene_symbol = setdiff(full_up, rownames(expr_mat)), stringsAsFactors = FALSE),
          file.path(out_dir, "NuStress_UP_not_in_matrix.csv"),
          row.names = FALSE)

write.csv(data.frame(gene_symbol = setdiff(full_down, rownames(expr_mat)), stringsAsFactors = FALSE),
          file.path(out_dir, "NuStress_DOWN_not_in_matrix.csv"),
          row.names = FALSE)

if (has_core) {
  write.csv(data.frame(gene_symbol = core_up_in, stringsAsFactors = FALSE),
            file.path(out_dir, "NuStressCore_UP_overlap_with_matrix.csv"),
            row.names = FALSE)
  
  write.csv(data.frame(gene_symbol = core_down_in, stringsAsFactors = FALSE),
            file.path(out_dir, "NuStressCore_DOWN_overlap_with_matrix.csv"),
            row.names = FALSE)
  
  write.csv(data.frame(gene_symbol = setdiff(core_up, rownames(expr_mat)), stringsAsFactors = FALSE),
            file.path(out_dir, "NuStressCore_UP_not_in_matrix.csv"),
            row.names = FALSE)
  
  write.csv(data.frame(gene_symbol = setdiff(core_down, rownames(expr_mat)), stringsAsFactors = FALSE),
            file.path(out_dir, "NuStressCore_DOWN_not_in_matrix.csv"),
            row.names = FALSE)
}

############################################################
## 7. Read RiboSis gene set
############################################################
if (!file.exists(ribosis_rdata)) {
  stop("RiboSis file not found: ", ribosis_rdata)
}

load(ribosis_rdata)

if (!exists("ribosis")) {
  stop("Object 'ribosis' not found in RiboSis activity.Rdata")
}

if (is.list(ribosis)) {
  ribosis_genes <- unique(toupper(trimws(unlist(ribosis[[1]]))))
} else {
  ribosis_genes <- unique(toupper(trimws(unlist(ribosis))))
}
ribosis_genes <- ribosis_genes[!is.na(ribosis_genes) & ribosis_genes != ""]

ribosis_in <- intersect(ribosis_genes, rownames(expr_mat))

cat("===== RiboSis gene set overlap =====\n")
cat("RiboSis overlap:", length(ribosis_in), "\n")

if (length(ribosis_in) < 5) {
  stop("Too few overlapping genes for RiboSis gene set.")
}

geneSets_ribosis <- list(RiboSis = ribosis_in)

write.csv(data.frame(gene_symbol = ribosis_in, stringsAsFactors = FALSE),
          file.path(out_dir, "RiboSis_overlap_with_matrix.csv"),
          row.names = FALSE)

############################################################
## 8. ssGSEA scoring function
############################################################
run_ssgsea <- function(expr_mat, gene_sets, minSize = 5, maxSize = 500) {
  param <- GSVA::ssgseaParam(
    exprData = expr_mat,
    geneSets = gene_sets,
    minSize = minSize,
    maxSize = maxSize,
    alpha = 0.25,
    normalize = TRUE
  )
  res <- gsva(param)
  return(res)
}

############################################################
## 9. Calculate NuStress Full, Core, and RiboSis scores
############################################################
## Full
ssgsea_full <- run_ssgsea(expr_mat, geneSets_full)
NuS_Full_score <- ssgsea_full["NuStress_UP", ] - ssgsea_full["NuStress_DOWN", ]

## Core
if (has_core) {
  ssgsea_core <- run_ssgsea(expr_mat, geneSets_core)
  NuS_Core_score <- ssgsea_core["NuStressCore_UP", ] - ssgsea_core["NuStressCore_DOWN", ]
}

## RiboSis
ssgsea_ribosis <- run_ssgsea(expr_mat, geneSets_ribosis)
RiboSis_score <- ssgsea_ribosis["RiboSis", ]

############################################################
## 10. Summarize score table
############################################################
scores_list <- list(
  data.frame(
    Sample = colnames(expr_mat),
    Score = as.numeric(NuS_Full_score),
    ScoreType = "NuS_Full",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Sample = colnames(expr_mat),
    Score = as.numeric(RiboSis_score),
    ScoreType = "RiboSis",
    stringsAsFactors = FALSE
  )
)

if (has_core) {
  scores_list <- append(scores_list, list(
    data.frame(
      Sample = colnames(expr_mat),
      Score = as.numeric(NuS_Core_score),
      ScoreType = "NuS_Core",
      stringsAsFactors = FALSE
    )
  ))
}

scores_all <- dplyr::bind_rows(scores_list) %>%
  dplyr::left_join(group_info, by = "Sample")

scores_all$Group <- factor(scores_all$Group, levels = levels(group_info$Group))
scores_all$ScoreType <- factor(scores_all$ScoreType,
                               levels = c("NuS_Full", "NuS_Core", "RiboSis"))

write.csv(scores_all,
          file.path(out_dir, "All_ssGSEA_scores.csv"),
          row.names = FALSE)


############################################################
## 11. Colors and plotting functions
## Boxplot + points only, no statistical analysis
############################################################
base_fill_colors <- c(
  "Control" = "#4C72B0",
  "OXA"     = "#6A8FD7",
  "L"       = "#C44E52",
  "L-OXA"   = "#8B1E3F"
)

make_group_colors <- function(groups) {
  cols <- base_fill_colors[groups]
  cols[is.na(cols)] <- "#999999"
  return(cols)
}

plot_box_point <- function(df, score_type, out_dir) {
  plot_df <- df %>% dplyr::filter(ScoreType == score_type)
  fill_colors <- make_group_colors(levels(plot_df$Group))
  
  p <- ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(
      width = 0.42,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.65,
      alpha = 0.55
    ) +
    geom_jitter(
      aes(color = Group),
      width = 0.10,
      size = 3.2,
      shape = 16,
      alpha = 0.95
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    labs(
      title = score_type,
      x = NULL,
      y = "ssGSEA score"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  ggsave(file.path(out_dir, paste0(score_type, "_boxplot_points.pdf")),
         p, width = 3, height = 4.2, bg = "white")
  ggsave(file.path(out_dir, paste0(score_type, "_boxplot_points.png")),
         p, width = 3, height = 4.2, dpi = 300, bg = "white")
  
  return(p)
}
p1 <- plot_box_point(scores_all, "NuS_Full", out_dir)
if (has_core) p2 <- plot_box_point(scores_all, "NuS_Core", out_dir)
p3 <- plot_box_point(scores_all, "RiboSis", out_dir)


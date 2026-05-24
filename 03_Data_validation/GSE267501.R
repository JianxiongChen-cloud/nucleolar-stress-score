rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE267501"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE267501"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(readxl)
library(GSVA)
library(ggplot2)
library(dplyr)

############################################################
## 2. Read the FPKM table
############################################################
data <- read_excel(file.path(data_dir, "GSE267501_HCT116_TP53WT_fpkm.xlsx"))
data <- as.data.frame(data)

cat("===== Raw table dimension =====\n")
print(dim(data))
cat("===== Column names =====\n")
print(colnames(data))

############################################################
## 3. Select expression columns
## Consistent with the DEG pipeline:
## first 6h_rep1-3,
## then 0h_rep1-3
############################################################
expr_gene_sub <- data[, c(5:7, 2:4), drop = FALSE]
expr_gene_sub <- as.matrix(expr_gene_sub)
mode(expr_gene_sub) <- "numeric"
rownames(expr_gene_sub) <- data$gene_id

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

group_info <- data.frame(
  Sample = colnames(expr_gene_sub),
  Group = factor(
    c(rep("BMH21", 3), rep("CONTROL", 3)),
    levels = c("CONTROL", "BMH21")
  ),
  stringsAsFactors = FALSE
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## 4. Low-expression gene filtering
## FPKM > 1 in at least 3 samples
############################################################
keep <- rowSums(expr_gene_sub > 1, na.rm = TRUE) >= 3
expr_gene_sub <- expr_gene_sub[keep, , drop = FALSE]

cat("===== Filtering summary =====\n")
cat("Genes after filtering:", nrow(expr_gene_sub), "\n")

############################################################
## 5. log2 transformation
############################################################
expr_gene_sub <- log2(expr_gene_sub + 1)

cat("===== log2(FPKM+1) range =====\n")
print(range(expr_gene_sub, na.rm = TRUE))

############################################################
## 6. Remove zero-variance genes
############################################################
gene_var <- apply(expr_gene_sub, 1, var, na.rm = TRUE)
expr_gene_sub <- expr_gene_sub[!is.na(gene_var) & gene_var > 0, , drop = FALSE]

cat("===== Matrix dimension after removing zero-variance genes =====\n")
print(dim(expr_gene_sub))

############################################################
## 7. Read the full gene set
############################################################
read_gene_vector <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, 1, drop = TRUE]
  x <- unique(toupper(trimws(x)))
  x <- x[!is.na(x) & x != ""]
  return(x)
}

full_up   <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_up_genes.csv"))
full_down <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_down_genes.csv"))

############################################################
## 8. Convert gene symbols to uppercase
## and intersect with the gene sets
############################################################
rownames(expr_gene_sub) <- toupper(rownames(expr_gene_sub))

full_up_in   <- intersect(full_up, rownames(expr_gene_sub))
full_down_in <- intersect(full_down, rownames(expr_gene_sub))

cat("===== Full gene set overlap =====\n")
cat("NuS_Up:", length(full_up_in), "\n")
cat("NuS_Down:", length(full_down_in), "\n")

if (length(full_up_in) < 5 || length(full_down_in) < 5) {
  stop("Too few overlapping genes for the full gene set.")
}

geneSets_full <- list(
  NuS_Up = full_up_in,
  NuS_Down = full_down_in
)

############################################################
## 9. Calculate NuS scores using ssGSEA
############################################################
param <- ssgseaParam(
  exprData = as.matrix(expr_gene_sub),
  geneSets = geneSets_full,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_res <- gsva(param)

cat("===== ssGSEA result dimension =====\n")
print(dim(ssgsea_res))
print(rownames(ssgsea_res))

if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea_res))) {
  stop("NuS_Up or NuS_Down not found in ssGSEA result.")
}

nus_score <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]

scores_all <- data.frame(
  Sample = colnames(expr_gene_sub),
  Score = as.numeric(nus_score),
  Method = "ssGSEA",
  GeneSet = "Full",
  stringsAsFactors = FALSE
)

scores_all <- left_join(scores_all, group_info, by = "Sample")
scores_all$Group <- factor(scores_all$Group, levels = c("CONTROL", "BMH21"))

write.csv(scores_all, "GSE267501_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 10. Statistical testing: Welch's t-test
############################################################
test_res <- data.frame(
  comparison = "BMH21 vs CONTROL",
  p_welch = tryCatch(
    t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  ),
  mean_control = mean(scores_all$Score[scores_all$Group == "CONTROL"], na.rm = TRUE),
  sd_control   = sd(scores_all$Score[scores_all$Group == "CONTROL"], na.rm = TRUE),
  n_control    = sum(scores_all$Group == "CONTROL"),
  mean_treat   = mean(scores_all$Score[scores_all$Group == "BMH21"], na.rm = TRUE),
  sd_treat     = sd(scores_all$Score[scores_all$Group == "BMH21"], na.rm = TRUE),
  n_treat      = sum(scores_all$Group == "BMH21"),
  diff_mean    = mean(scores_all$Score[scores_all$Group == "BMH21"], na.rm = TRUE) -
    mean(scores_all$Score[scores_all$Group == "CONTROL"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(test_res, "GSE267501_NuS_score_statistics_ttest.csv", row.names = FALSE)
print(test_res)

############################################################
## 11. Visualization
############################################################
fill_colors <- c(
  "CONTROL" = "#C4E7C1",
  "BMH21"   = "#51B3D1"
)

p <- ggplot(scores_all, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.8, alpha = 0.85) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "NuS ssGSEA scores in GSE267501",
    x = NULL,
    y = "NuS score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("GSE267501_NuS_ssGSEA_violinplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE267501_NuS_ssGSEA_violinplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 12. Save workspace
############################################################
save(
  data,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE267501_NuS_ssGSEA_workspace.Rdata"
)
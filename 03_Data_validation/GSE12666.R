rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE12666"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE12666"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)
library(multcomp)
library(readr)

############################################################
## 2. Load the expression matrix workspace generated in the previous step
############################################################
load(file.path(data_dir, "GSE12666_step1_workspace.Rdata"))

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

cat("===== Sample names =====\n")
print(colnames(expr_gene))

############################################################
## 3. Select samples for analysis
##    Columns 13:16 and 25:28 were retained in the previous step
############################################################
expr_gene_sub <- expr_gene[, c(13:16, 25:28), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

############################################################
## 4. Group information
##    Please modify the Group labels below according to your actual sample grouping
############################################################
group_info <- data.frame(
  Sample = colnames(expr_gene_sub),
  Group = factor(
    c(
      rep("BMH-21", 4),
      rep("Control", 4)
    ),
    levels = c("Control", "BMH-21")
  ),
  stringsAsFactors = FALSE
)

print(group_info)

############################################################
## 5. Read the full gene set
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
## 6. Convert expression matrix gene symbols to uppercase
##    and intersect with the gene sets
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
## 7. Calculate NuS scores using ssGSEA
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
scores_all$Group <- factor(scores_all$Group, levels = c("Control", "BMH-21"))

write.csv(scores_all, "GSE12666_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 8. Statistical testing
##    Welch's t-test is more appropriate for two-group comparisons
############################################################
test_res <- data.frame(
  comparison = "BMH-21 vs Control",
  p_welch = tryCatch(
    t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  ),
  mean_control = mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  sd_control   = sd(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  n_control    = sum(scores_all$Group == "Control"),
  mean_treat   = mean(scores_all$Score[scores_all$Group == "BMH-21"], na.rm = TRUE),
  sd_treat     = sd(scores_all$Score[scores_all$Group == "BMH-21"], na.rm = TRUE),
  n_treat      = sum(scores_all$Group == "BMH-21"),
  diff_mean    = mean(scores_all$Score[scores_all$Group == "BMH-21"], na.rm = TRUE) -
    mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(test_res, "GSE12666_NuS_score_statistics_ttest.csv", row.names = FALSE)
print(test_res)

############################################################
## 9. Visualization
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "BMH-21" = "#51B3D1"
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
    title = "NuS ssGSEA scores in GSE12666",
    x = NULL,
    y = "NuS"
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

ggsave("GSE12666_NuS_ssGSEA_violinplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE12666_NuS_ssGSEA_violinplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 10. Save workspace
############################################################
save(
  expr_gene,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE12666_NuS_ssGSEA_workspace.Rdata"
)
rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE62593"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE62593"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. 加载R包
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)

############################################################
## 2. 读取上一步保存的表达矩阵workspace
############################################################
load(file.path(data_dir, "GSE62593_step1_workspace.Rdata"))

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

cat("===== All sample names =====\n")
print(colnames(expr_gene))

############################################################
## 3. 提取目标样本
##    columns 1:3 = CONTROL
##    columns 7:9 = ActD
############################################################
expr_gene_sub <- expr_gene[, c(1:3, 7:9), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

group_info <- data.frame(
  Sample = colnames(expr_gene_sub),
  Group = factor(
    c(rep("Control", 3), rep("ActD", 3)),
    levels = c("Control", "ActD")
  ),
  stringsAsFactors = FALSE
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## 4. 读取基因全集
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
## 5. 表达矩阵基因名统一为大写，并与基因集取交集
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
## 6. ssGSEA计算NuS score
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
scores_all$Group <- factor(scores_all$Group, levels = c("Control", "ActD"))

write.csv(scores_all, "GSE62593_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 7. 统计检验：Welch t-test
############################################################
test_res <- data.frame(
  comparison = "ActD vs Control",
  p_welch = tryCatch(
    t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  ),
  mean_control = mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  sd_control   = sd(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  n_control    = sum(scores_all$Group == "Control"),
  mean_treat   = mean(scores_all$Score[scores_all$Group == "ActD"], na.rm = TRUE),
  sd_treat     = sd(scores_all$Score[scores_all$Group == "ActD"], na.rm = TRUE),
  n_treat      = sum(scores_all$Group == "ActD"),
  diff_mean    = mean(scores_all$Score[scores_all$Group == "ActD"], na.rm = TRUE) -
    mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(test_res, "GSE62593_NuS_score_statistics_ttest.csv", row.names = FALSE)
print(test_res)

############################################################
## 8. 绘图
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "ActD"    = "#51B3D1"
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
    title = "NuS ssGSEA scores in GSE62593",
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

ggsave("GSE62593_NuS_ssGSEA_violinplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE62593_NuS_ssGSEA_violinplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 9. 保存workspace
############################################################
save(
  expr_gene,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE62593_NuS_ssGSEA_workspace.Rdata"
)
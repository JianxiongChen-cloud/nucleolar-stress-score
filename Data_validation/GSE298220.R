rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE298220"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE298220"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. 加载R包
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(tibble)

############################################################
## 2. 读取step1表达矩阵
############################################################
load(file.path(data_dir, "GSE298220_step1_output.Rdata"))

cat("===== final_expr dimension =====\n")
print(dim(final_expr))

cat("===== All sample names =====\n")
print(colnames(final_expr))

############################################################
## 3. 取9个样本：CONTROL / CX5461 / BMH21
##    sample1-3 = CONTROL
##    sample4-6 = CX5461
##    sample7-9 = BMH21
############################################################
expr_gene_sub <- final_expr[, 1:9, drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

group_info <- data.frame(
  Sample = colnames(expr_gene_sub),
  Group = factor(
    c(rep("CONTROL", 3), rep("CX5461", 3), rep("BMH21", 3)),
    levels = c("CONTROL", "CX5461", "BMH21")
  ),
  stringsAsFactors = FALSE
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## 4. 低表达过滤
##    与你的DEG流程一致：expr > 1 in at least 3 samples
############################################################
keep <- rowSums(expr_gene_sub > 1, na.rm = TRUE) >= 3
expr_gene_sub <- expr_gene_sub[keep, , drop = FALSE]

cat("===== Filtering summary =====\n")
cat("Genes after filtering:", nrow(expr_gene_sub), "\n")

############################################################
## 5. log2转换
############################################################
expr_gene_sub <- log2(expr_gene_sub + 1)

cat("===== log2(x+1) range =====\n")
print(range(expr_gene_sub, na.rm = TRUE))

############################################################
## 6. 去除零方差基因
############################################################
gene_var <- apply(expr_gene_sub, 1, var, na.rm = TRUE)
expr_gene_sub <- expr_gene_sub[!is.na(gene_var) & gene_var > 0, , drop = FALSE]

cat("===== Matrix dimension after removing zero-variance genes =====\n")
print(dim(expr_gene_sub))

############################################################
## 7. 读取基因全集
############################################################
read_gene_vector <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, 1, drop = TRUE]
  x <- unique(toupper(trimws(x)))
  x <- x[!is.na(x) & x != ""]
  x
}

full_up   <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_up_genes.csv"))
full_down <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_down_genes.csv"))

############################################################
## 8. gene symbol大写，并与gene set取交集
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
## 9. ssGSEA计算NuS score
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
scores_all$Group <- factor(scores_all$Group, levels = c("CONTROL", "CX5461", "BMH21"))

write.csv(scores_all, "GSE298220_NuS_scores_ssGSEA_all_groups.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 10. 统计检验：one-way ANOVA + Dunnett vs CONTROL
############################################################
scores_all$Group <- droplevels(scores_all$Group)
scores_all$Group <- relevel(scores_all$Group, ref = "CONTROL")

fit <- aov(Score ~ Group, data = scores_all)
anova_p <- tryCatch(
  summary(fit)[[1]][["Pr(>F)"]][1],
  error = function(e) NA_real_
)

mean_sd_df <- scores_all %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    N = sum(!is.na(Score)),
    .groups = "drop"
  )

control_row <- mean_sd_df %>% filter(Group == "CONTROL")

suppressPackageStartupMessages(library(multcomp))
dunnett_res <- tryCatch(
  summary(glht(fit, linfct = mcp(Group = "Dunnett"))),
  error = function(e) NULL
)

if (is.null(dunnett_res)) {
  test_res <- data.frame(
    comparison = NA_character_,
    treatment = NA_character_,
    p_anova = anova_p,
    p_dunnett = NA_real_,
    mean_control = control_row$Mean,
    sd_control = control_row$SD,
    n_control = control_row$N,
    mean_treat = NA_real_,
    sd_treat = NA_real_,
    n_treat = NA_real_,
    diff_mean = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  comp_names <- names(dunnett_res$test$coefficients)
  comp_diff  <- as.numeric(dunnett_res$test$coefficients)
  comp_p     <- as.numeric(dunnett_res$test$pvalues)
  
  treat_names <- sub(" - CONTROL$", "", comp_names)
  treat_df <- mean_sd_df %>% filter(Group != "CONTROL") %>% as.data.frame()
  match_idx <- match(treat_names, as.character(treat_df$Group))
  
  test_res <- data.frame(
    comparison = comp_names,
    treatment = treat_names,
    p_anova = rep(anova_p, length(comp_names)),
    p_dunnett = comp_p,
    mean_control = rep(control_row$Mean, length(comp_names)),
    sd_control = rep(control_row$SD, length(comp_names)),
    n_control = rep(control_row$N, length(comp_names)),
    mean_treat = treat_df$Mean[match_idx],
    sd_treat = treat_df$SD[match_idx],
    n_treat = treat_df$N[match_idx],
    diff_mean = comp_diff,
    stringsAsFactors = FALSE
  )
}

write.csv(test_res, "GSE298220_NuS_score_statistics_ANOVA_Dunnett.csv", row.names = FALSE)
print(test_res)

############################################################
## 11. 所有组放在一张图
############################################################
fill_colors <- c(
  "CONTROL" = "#C4E7C1",
  "CX5461"  = "#51B3D1",
  "BMH21"   = "#08589E"
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
    title = "NuS ssGSEA scores in GSE298220",
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

ggsave("GSE298220_NuS_ssGSEA_all_groups_oneplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE298220_NuS_ssGSEA_all_groups_oneplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 12. 保存workspace
############################################################
save(
  final_expr,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE298220_NuS_ssGSEA_all_groups_oneplot_workspace.Rdata"
)
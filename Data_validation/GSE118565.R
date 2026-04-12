rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE118565"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE118565"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. 加载R包
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)
library(multcomp)

############################################################
## 2. 读取 step1 保存的表达矩阵
############################################################
load(file.path(data_dir, "GSE118565_step1_workspace.Rdata"))

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

cat("===== All sample names =====\n")
print(colnames(tpm_final))

############################################################
## 3. 选择目标样本
############################################################
expr_gene_sub <- tpm_final

############################################################
## 4. 根据真实样本名建立分组信息
############################################################
group_map <- data.frame(
  Sample = c(
    "GSM3333289","GSM3333290","GSM3333291","GSM3333292","GSM3333293","GSM3333294",
    "GSM3333295","GSM3333296","GSM3333297","GSM3333298","GSM3333299","GSM3333300",
    "GSM3333301","GSM3333302","GSM3333303","GSM3333304","GSM3333305","GSM3333306"
  ),
  Group = c(
    "Veh","CX-5461","I-BET151","COMBO","Doxo","Doxo+I-BET151",
    "Veh","CX-5461","I-BET151","COMBO","Doxo","Doxo+I-BET151",
    "Veh","CX-5461","I-BET151","COMBO","Doxo","Doxo+I-BET151"
  ),
  stringsAsFactors = FALSE
)

## 只保留 group_map 中定义的样本，并按 group_map 顺序排列
expr_gene_sub <- expr_gene_sub[, group_map$Sample, drop = FALSE]

group_info <- group_map
group_info$Group <- factor(
  group_info$Group,
  levels = c("Veh", "CX-5461", "I-BET151", "COMBO", "Doxo", "Doxo+I-BET151")
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## 5. 低表达基因过滤
##    TPM > 1 in at least 3 samples
############################################################
keep <- rowSums(expr_gene_sub > 1, na.rm = TRUE) >= 3
expr_gene_sub <- expr_gene_sub[keep, , drop = FALSE]

cat("===== Matrix dimension after filtering =====\n")
print(dim(expr_gene_sub))

############################################################
## 6. log2 转换
############################################################
expr_gene_sub <- log2(expr_gene_sub + 1)

cat("===== Overall range after log2(TPM+1) =====\n")
print(range(expr_gene_sub, na.rm = TRUE))

############################################################
## 7. 去除零方差基因
############################################################
gene_var <- apply(expr_gene_sub, 1, var, na.rm = TRUE)
expr_gene_sub <- expr_gene_sub[!is.na(gene_var) & gene_var > 0, , drop = FALSE]

cat("===== Matrix dimension after removing zero-variance genes =====\n")
print(dim(expr_gene_sub))

############################################################
## 8. 读取基因全集
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
## 9. 表达矩阵基因名统一为大写，并与基因集取交集
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
## 10. ssGSEA 计算 NuS score
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
scores_all$Group <- factor(
  scores_all$Group,
  levels = c("Veh", "CX-5461", "I-BET151", "COMBO", "Doxo", "Doxo+I-BET151")
)

write.csv(scores_all, "NuS_scores_ssGSEA_6groups.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 11. 统计分析：ANOVA + Dunnett vs Veh
############################################################
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

veh_row <- mean_sd_df %>% filter(Group == "Veh")

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
    mean_control = veh_row$Mean,
    sd_control = veh_row$SD,
    n_control = veh_row$N,
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
  
  treat_names <- sub(" - Veh$", "", comp_names)
  treat_df <- mean_sd_df %>% filter(Group != "Veh") %>% as.data.frame()
  match_idx <- match(treat_names, as.character(treat_df$Group))
  
  test_res <- data.frame(
    comparison = comp_names,
    treatment = treat_names,
    p_anova = rep(anova_p, length(comp_names)),
    p_dunnett = comp_p,
    mean_control = rep(veh_row$Mean, length(comp_names)),
    sd_control = rep(veh_row$SD, length(comp_names)),
    n_control = rep(veh_row$N, length(comp_names)),
    mean_treat = treat_df$Mean[match_idx],
    sd_treat = treat_df$SD[match_idx],
    n_treat = treat_df$N[match_idx],
    diff_mean = comp_diff,
    stringsAsFactors = FALSE
  )
}

write.csv(test_res, "NuS_score_statistics_ANOVA_Dunnett_vsVeh.csv", row.names = FALSE)
print(test_res)

############################################################
## 12. 绘图：6组小提琴图
############################################################
fill_colors <- c(
  "Veh" = "#C4E7C1",
  "CX-5461" = "#93D4BC",
  "I-BET151" = "#6CC4C9",
  "COMBO" = "#51B3D1",
  "Doxo" = "#2C7FB8",
  "Doxo+I-BET151" = "#08589E"
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
    title = "NuS ssGSEA scores",
    x = NULL,
    y = "NuS score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("NuS_ssGSEA_violinplot_6groups.pdf", p, width = 6.2, height = 5, bg = "white")
ggsave("NuS_ssGSEA_violinplot_6groups.png", p, width = 6.2, height = 5, dpi = 300, bg = "white")

############################################################
## 13. 保存 workspace
############################################################
save(
  tpm_final,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "NuS_ssGSEA_6groups_workspace.Rdata"
)






#######################################只取两组##############################
rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE118565"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE118565"
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
## 2. 读取step1保存的表达矩阵
############################################################
load(file.path(data_dir, "GSE118565_step1_workspace.Rdata"))

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

cat("===== All sample names =====\n")
print(colnames(tpm_final))

############################################################
## 3. 选择目标样本并按分组顺序重排
##    original selection: 1,2,7,8,13,14
##    regrouped as: 1,7,13 = Control; 2,8,14 = CX5461
############################################################
expr_gene_sub <- tpm_final[, c(1, 2, 7, 8, 13, 14), drop = FALSE]
expr_gene_sub <- expr_gene_sub[, c(1, 3, 5, 2, 4, 6), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

group_info <- data.frame(
  Sample = colnames(expr_gene_sub),
  Group = factor(
    c(rep("Control", 3), rep("CX5461", 3)),
    levels = c("Control", "CX5461")
  ),
  stringsAsFactors = FALSE
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## 4. 低表达基因过滤
##    与前面DEG分析保持一致：TPM > 1 in at least 3 samples
############################################################
keep <- rowSums(expr_gene_sub > 1, na.rm = TRUE) >= 3
expr_gene_sub <- expr_gene_sub[keep, , drop = FALSE]

cat("===== Matrix dimension after filtering =====\n")
print(dim(expr_gene_sub))

############################################################
## 5. log2转换
############################################################
expr_gene_sub <- log2(expr_gene_sub + 1)

cat("===== Overall range after log2(TPM+1) =====\n")
print(range(expr_gene_sub, na.rm = TRUE))

############################################################
## 6. 读取基因全集
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
## 7. 表达矩阵基因名统一为大写，并与基因集取交集
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
## 8. ssGSEA计算NuS score
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
scores_all$Group <- factor(scores_all$Group, levels = c("Control", "CX5461"))

write.csv(scores_all, "GSE118565_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 9. 统计检验：Welch t-test
############################################################
test_res <- data.frame(
  comparison = "CX5461 vs Control",
  p_welch = tryCatch(
    t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  ),
  mean_control = mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  sd_control   = sd(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  n_control    = sum(scores_all$Group == "Control"),
  mean_treat   = mean(scores_all$Score[scores_all$Group == "CX5461"], na.rm = TRUE),
  sd_treat     = sd(scores_all$Score[scores_all$Group == "CX5461"], na.rm = TRUE),
  n_treat      = sum(scores_all$Group == "CX5461"),
  diff_mean    = mean(scores_all$Score[scores_all$Group == "CX5461"], na.rm = TRUE) -
    mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(test_res, "GSE118565_NuS_score_statistics_ttest.csv", row.names = FALSE)
print(test_res)

############################################################
## 10. 绘图
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "CX5461"  = "#51B3D1"
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
    title = "NuS ssGSEA scores in GSE118565",
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

ggsave("GSE118565_NuS_ssGSEA_violinplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE118565_NuS_ssGSEA_violinplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 11. 保存workspace
############################################################
save(
  tpm_final,
  expr_gene_sub,
  group_info,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE118565_NuS_ssGSEA_workspace.Rdata"
)












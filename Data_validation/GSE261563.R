rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE261563"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE261563"
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
load(file.path(data_dir, "GSE261563_step1_workspace.Rdata"))

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

cat("===== All sample names =====\n")
print(colnames(tpm_final))

############################################################
## 3. 读取基因全集
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
## 4. 通用评分函数
############################################################
run_nus_ssgsea_tpm <- function(expr_sub, group_info, prefix) {
  
  expr_sub <- as.matrix(expr_sub)
  mode(expr_sub) <- "numeric"
  
  cat("\n=============================\n")
  cat("Running:", prefix, "\n")
  cat("=============================\n")
  
  cat("===== Selected samples =====\n")
  print(colnames(expr_sub))
  
  cat("===== Group info =====\n")
  print(group_info)
  
  ##########################################################
  ## 4.1 TPM > 1 in at least 3 samples
  ##########################################################
  keep <- rowSums(expr_sub > 1, na.rm = TRUE) >= 3
  expr_sub <- expr_sub[keep, , drop = FALSE]
  
  cat("===== Matrix dimension after filtering =====\n")
  print(dim(expr_sub))
  
  ##########################################################
  ## 4.2 log2(TPM+1)
  ##########################################################
  expr_sub <- log2(expr_sub + 1)
  
  ##########################################################
  ## 4.3 去除零方差基因
  ##########################################################
  gene_var <- apply(expr_sub, 1, var, na.rm = TRUE)
  expr_sub <- expr_sub[!is.na(gene_var) & gene_var > 0, , drop = FALSE]
  
  cat("===== Matrix dimension after removing zero-variance genes =====\n")
  print(dim(expr_sub))
  
  ##########################################################
  ## 4.4 gene symbol大写，并与gene set取交集
  ##########################################################
  rownames(expr_sub) <- toupper(rownames(expr_sub))
  
  full_up_in   <- intersect(full_up, rownames(expr_sub))
  full_down_in <- intersect(full_down, rownames(expr_sub))
  
  cat("===== Full gene set overlap =====\n")
  cat("NuS_Up:", length(full_up_in), "\n")
  cat("NuS_Down:", length(full_down_in), "\n")
  
  if (length(full_up_in) < 5 || length(full_down_in) < 5) {
    stop(paste("Too few overlapping genes for", prefix))
  }
  
  geneSets_full <- list(
    NuS_Up = full_up_in,
    NuS_Down = full_down_in
  )
  
  ##########################################################
  ## 4.5 ssGSEA
  ##########################################################
  param <- ssgseaParam(
    exprData = expr_sub,
    geneSets = geneSets_full,
    minSize = 5,
    maxSize = 500,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_res <- gsva(param)
  
  if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea_res))) {
    stop(paste("NuS_Up or NuS_Down not found in", prefix))
  }
  
  nus_score <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]
  
  scores_all <- data.frame(
    Sample = colnames(expr_sub),
    Score = as.numeric(nus_score),
    Method = "ssGSEA",
    GeneSet = "Full",
    stringsAsFactors = FALSE
  )
  
  scores_all <- left_join(scores_all, group_info, by = "Sample")
  scores_all$Group <- factor(scores_all$Group, levels = c("Control", "CX5461"))
  
  ##########################################################
  ## 4.6 Welch t-test
  ##########################################################
  test_res <- data.frame(
    comparison = paste0(prefix, ": CX5461 vs Control"),
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
  
  write.csv(scores_all, paste0(prefix, "_NuS_scores_ssGSEA.csv"), row.names = FALSE)
  write.csv(test_res, paste0(prefix, "_NuS_score_statistics_ttest.csv"), row.names = FALSE)
  
  save(
    expr_sub,
    group_info,
    geneSets_full,
    scores_all,
    test_res,
    file = paste0(prefix, "_NuS_ssGSEA_workspace.Rdata")
  )
  
  return(list(
    scores_all = scores_all,
    test_res = test_res
  ))
}

############################################################
## 5. 543 cells
## columns 7:12, first 3 = Control, last 3 = CX5461
############################################################
dat_543 <- tpm_final[, 7:12, drop = FALSE]

group_info_543 <- data.frame(
  Sample = colnames(dat_543),
  Group = c(rep("Control", 3), rep("CX5461", 3)),
  Subset = "543",
  stringsAsFactors = FALSE
)

res_543 <- run_nus_ssgsea_tpm(
  expr_sub = dat_543,
  group_info = group_info_543,
  prefix = "GSE261563_543"
)

############################################################
## 6. 603 cells
## columns 13:18, first 3 = Control, last 3 = CX5461
############################################################
dat_603 <- tpm_final[, 13:18, drop = FALSE]

group_info_603 <- data.frame(
  Sample = colnames(dat_603),
  Group = c(rep("Control", 3), rep("CX5461", 3)),
  Subset = "603",
  stringsAsFactors = FALSE
)

res_603 <- run_nus_ssgsea_tpm(
  expr_sub = dat_603,
  group_info = group_info_603,
  prefix = "GSE261563_603"
)

############################################################
## 7. 合并两个子集评分结果
############################################################
scores_plot <- bind_rows(
  res_543$scores_all,
  res_603$scores_all
)

scores_plot$Subset <- factor(scores_plot$Subset, levels = c("543", "603"))

scores_plot$Group_Subset <- factor(
  paste(scores_plot$Group, scores_plot$Subset, sep = "_"),
  levels = c("Control_543", "CX5461_543", "Control_603", "CX5461_603")
)

write.csv(scores_plot, "GSE261563_NuS_scores_ssGSEA_all_groups_for_plot.csv", row.names = FALSE)

test_all <- bind_rows(
  res_543$test_res,
  res_603$test_res
)

write.csv(test_all, "GSE261563_NuS_score_statistics_all_groups_ttest.csv", row.names = FALSE)

############################################################
## 8. 所有组放在一张图
############################################################
fill_colors <- c(
  "Control_543" = "#C4E7C1",
  "CX5461_543"  = "#51B3D1",
  "Control_603" = "#B7DCE8",
  "CX5461_603"  = "#08589E"
)

p <- ggplot(scores_plot, aes(x = Group_Subset, y = Score, fill = Group_Subset)) +
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
    title = "NuS ssGSEA scores in GSE261563",
    x = NULL,
    y = "NuS score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("GSE261563_NuS_ssGSEA_all_groups_oneplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE261563_NuS_ssGSEA_all_groups_oneplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 9. 保存总workspace
############################################################
save(
  tpm_final,
  dat_543,
  dat_603,
  res_543,
  res_603,
  scores_plot,
  test_all,
  file = "GSE261563_NuS_ssGSEA_all_groups_oneplot_workspace.Rdata"
)
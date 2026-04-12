############################################################
## GSE232930
## 完整修正版：ssGSEA打分 + 统计 + 绘图
## 只分析基因全集
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE232930/")

library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(ggpubr)

geo_id <- "GSE232930"
expr_file <- "expressionmetrix_GSE.csv"
geneset_file <- "/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata"

############################################################
## Step 1. 读取表达矩阵
############################################################
expr <- read.csv(
  expr_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr.mat <- as.matrix(expr)
mode(expr.mat) <- "numeric"

cat("===== Expression matrix dimension =====\n")
print(dim(expr.mat))
cat("===== Expression range =====\n")
print(range(expr.mat, na.rm = TRUE))
cat("===== Sample names =====\n")
print(colnames(expr.mat))

## 如有 NA，去掉对应基因
expr.mat <- expr.mat[rowSums(is.na(expr.mat)) == 0, , drop = FALSE]

############################################################
## Step 2. 分组信息
############################################################
## 如果样本名已知，可直接替换成真实样本名
group_info <- data.frame(
  Sample = c("X1", "X2", "X3", "X4", "X5", "X6"),
  Group = rep(c("Control", "so", "aso"), each = 2),
  stringsAsFactors = FALSE
)

control_group <- "Control"

## 若 group_info$Sample 不是实际列名，则自动替换
if (!all(group_info$Sample %in% colnames(expr.mat))) {
  if (nrow(group_info) == ncol(expr.mat)) {
    group_info$Sample <- colnames(expr.mat)
    cat("group_info$Sample 已自动替换为表达矩阵列名。\n")
  } else {
    stop("group_info 的样本数与表达矩阵列数不一致，请检查。")
  }
}

## 保证顺序一致
expr.mat <- expr.mat[, group_info$Sample, drop = FALSE]

group_info$Group <- factor(
  group_info$Group,
  levels = unique(group_info$Group)
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## Step 3. 读取基因全集
############################################################
load(geneset_file)

geneSets_all <- list(
  NuS_Up = unique(toupper(NuS_up_all_unique)),
  NuS_Down = unique(toupper(NuS_down_all_unique))
)

rownames(expr.mat) <- toupper(trimws(rownames(expr.mat)))
expr.mat <- expr.mat[!duplicated(rownames(expr.mat)), , drop = FALSE]

geneSets_all$NuS_Up <- intersect(rownames(expr.mat), geneSets_all$NuS_Up)
geneSets_all$NuS_Down <- intersect(rownames(expr.mat), geneSets_all$NuS_Down)

cat("===== Gene set overlap =====\n")
cat("NuS_Up:", length(geneSets_all$NuS_Up), "\n")
cat("NuS_Down:", length(geneSets_all$NuS_Down), "\n")

if (length(geneSets_all$NuS_Up) < 5 || length(geneSets_all$NuS_Down) < 5) {
  stop("Too few overlapping genes for ssGSEA.")
}

############################################################
## Step 4. ssGSEA
############################################################
ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr.mat,
  geneSets = geneSets_all,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea.res <- gsva(ssgsea_param)

if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea.res))) {
  stop("ssGSEA 结果中未找到 NuS_Up 或 NuS_Down")
}

nus_score <- ssgsea.res["NuS_Up", ] - ssgsea.res["NuS_Down", ]

scores <- data.frame(
  Sample = colnames(expr.mat),
  NuS_Up = as.numeric(ssgsea.res["NuS_Up", ]),
  NuS_Down = as.numeric(ssgsea.res["NuS_Down", ]),
  NuS_Score = as.numeric(nus_score),
  stringsAsFactors = FALSE
)

scores <- merge(scores, group_info, by = "Sample")
scores$Group <- factor(scores$Group, levels = levels(group_info$Group))

plot_data <- data.frame(
  Score = scores$NuS_Score,
  Group = scores$Group
)

cat("===== NuS scores =====\n")
print(scores)

write.csv(scores, paste0(geo_id, "_NuS_all_ssGSEA_scores.csv"), row.names = FALSE)
write.csv(as.data.frame(ssgsea.res), paste0(geo_id, "_NuS_all_ssGSEA_matrix.csv"), row.names = TRUE)

############################################################
## Step 5. 统计分析
############################################################
n_groups <- length(unique(plot_data$Group))
cat("===== Number of groups =====\n")
print(n_groups)

if (n_groups == 2) {
  overall_test <- t.test(Score ~ Group, data = plot_data)
  overall_p <- overall_test$p.value
  
  stat_summary <- data.frame(
    method = "t.test",
    p.value = overall_p,
    stringsAsFactors = FALSE
  )
  
} else {
  overall_test <- aov(Score ~ Group, data = plot_data)
  overall_p <- summary(overall_test)[[1]][["Pr(>F)"]][1]
  
  tukey_res <- TukeyHSD(overall_test)
  tukey_df <- as.data.frame(tukey_res$Group)
  tukey_df$comparison <- rownames(tukey_df)
  rownames(tukey_df) <- NULL
  
  write.csv(tukey_df, paste0(geo_id, "_NuS_all_TukeyHSD.csv"), row.names = FALSE)
  
  stat_summary <- data.frame(
    method = "ANOVA",
    p.value = overall_p,
    stringsAsFactors = FALSE
  )
}

write.csv(stat_summary, paste0(geo_id, "_NuS_all_overall_stats.csv"), row.names = FALSE)

## 各处理组 vs 对照组
control_data <- plot_data$Score[plot_data$Group == control_group]
vs_control_res <- data.frame()

for (grp in unique(plot_data$Group)) {
  if (grp != control_group) {
    grp_data <- plot_data$Score[plot_data$Group == grp]
    pval <- t.test(grp_data, control_data)$p.value
    
    vs_control_res <- rbind(
      vs_control_res,
      data.frame(
        comparison = paste(grp, "vs", control_group),
        p.value = pval,
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(vs_control_res, paste0(geo_id, "_NuS_all_vsControl_stats.csv"), row.names = FALSE)

############################################################
## Step 6. 绘图颜色
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "so" = "#93D4BC",
  "aso" = "#08589E"
)

############################################################
## Step 7. 作图：总体比较（小提琴图风格）
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "so" = "#93D4BC",
  "aso" = "#08589E"
)

p1 <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(
    alpha = 0.85,
    trim = TRUE,
    scale = "width",
    width = 0.7,
    color = "black",
    linewidth = 1.0,
    adjust = 1.2
  ) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.9,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = paste0("NuS Score (ssGSEA) - ", geo_id, " - all"),
    x = NULL,
    y = "NuS Score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(
      hjust = 0.5, size = 16, face = "bold", color = "#13646B",
      margin = margin(b = 15)
    ),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

if (n_groups == 2) {
  p1 <- p1 + stat_compare_means(
    method = "t.test",
    label = "p.format",
    label.x.npc = "center",
    size = 5
  )
} else {
  p1 <- p1 + stat_compare_means(
    method = "anova",
    label = "p.format",
    label.x.npc = "center",
    size = 5
  )
}

print(p1)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.pdf"),
  plot = p1,
  width = 3,
  height = 4,
  bg = "white"
)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.png"),
  plot = p1,
  width = 3,
  height = 4,
  dpi = 600,
  bg = "white"
)
############################################################
## Step 8. 作图：相对对照组比较（小提琴图风格）
############################################################
p2 <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(
    alpha = 0.85,
    trim = TRUE,
    scale = "width",
    width = 0.7,
    color = "black",
    linewidth = 1.0,
    adjust = 1.2
  ) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.9,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = paste0("NuS Score (ssGSEA) - ", geo_id, " - all"),
    subtitle = paste0("Reference: ", control_group),
    x = NULL,
    y = "NuS Score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(
      hjust = 0.5, size = 16, face = "bold", color = "#13646B",
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  stat_compare_means(
    method = "t.test",
    ref.group = control_group,
    label = "p.signif",
    size = 6,
    vjust = 0.5
  )

print(p2)

ggsave(
  filename = paste0(geo_id, "_NuS_all_vsControl_violin.pdf"),
  plot = p2,
  width = 3,
  height = 4,
  bg = "white"
)

ggsave(
  filename = paste0(geo_id, "_NuS_all_vsControl_violin.png"),
  plot = p2,
  width = 3,
  height = 4,
  dpi = 600,
  bg = "white"
)

############################################################
## Step 9. 保存 workspace
############################################################
save(
  expr.mat,
  group_info,
  geneSets_all,
  ssgsea.res,
  scores,
  plot_data,
  stat_summary,
  vs_control_res,
  file = paste0(geo_id, "_NuS_all_workspace.Rdata")
)

cat("分析完成。\n")
cat("结果文件保存在:\n")
cat(getwd(), "\n")
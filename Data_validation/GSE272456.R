############################################################
## GSE272456
## 完整优化版：表达矩阵构建 + ssGSEA + 统计 + 绘图
## 只分析基因全集
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE272456/")

library(readr)
library(dplyr)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)

geo_id <- "GSE272456"
input_file <- "GSE272456_PN0129B_HCT116_raw_counts.txt.gz"
geneset_file <- "/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata"

############################################################
## Step 1. 读取 raw counts
############################################################
data_raw <- read_table(
  input_file,
  comment = "!",
  show_col_types = FALSE
)

cat("===== Raw table dimension =====\n")
print(dim(data_raw))
cat("===== Raw column names =====\n")
print(colnames(data_raw))

expr <- as.data.frame(data_raw, check.names = FALSE)

cat("总基因数:", nrow(expr), "\n")
cat("唯一基因名数量:", length(unique(expr$GENE)), "\n")
cat("重复基因名数量:", nrow(expr) - length(unique(expr$GENE)), "\n")

dup_genes <- expr$GENE[duplicated(expr$GENE)]
cat("重复基因名示例:\n")
print(head(dup_genes))

############################################################
## Step 2. 构建表达矩阵
## 对重复基因名添加后缀，避免行名冲突
############################################################
rownames(expr) <- make.names(expr$GENE, unique = TRUE)
expr$GENE <- NULL

expr.mat <- as.matrix(expr)
mode(expr.mat) <- "numeric"

cat("===== Expression matrix dimension =====\n")
print(dim(expr.mat))
cat("===== Expression range before transform =====\n")
print(range(expr.mat, na.rm = TRUE))
cat("===== Sample names =====\n")
print(colnames(expr.mat))

save(expr.mat, file = paste0(geo_id, "_step1_workspace.Rdata"))
write.csv(expr.mat, file = "expressionmetrix_GSE.csv")

############################################################
## Step 3. log2 转换
############################################################
expr.mat <- log2(expr.mat + 1)

cat("===== Expression range after log2(x+1) =====\n")
print(range(expr.mat, na.rm = TRUE))

## 去除含 NA 的基因
expr.mat <- expr.mat[rowSums(is.na(expr.mat)) == 0, , drop = FALSE]

############################################################
## Step 4. 分组信息
############################################################
group_info <- data.frame(
  Sample = colnames(expr.mat),
  Group = rep(c("Control", "24h_5FU", "48h_5FU"), each = 3),
  stringsAsFactors = FALSE
)

group_info$Group <- factor(
  group_info$Group,
  levels = c("Control", "24h_5FU", "48h_5FU")
)

control_group <- "Control"

cat("===== Group info =====\n")
print(group_info)

############################################################
## Step 5. 读取基因全集
############################################################
load(geneset_file)

geneSets_all <- list(
  NuS_Up = unique(toupper(NuS_up_all_unique)),
  NuS_Down = unique(toupper(NuS_down_all_unique))
)

rownames(expr.mat) <- toupper(trimws(rownames(expr.mat)))
expr.mat <- expr.mat[!duplicated(rownames(expr.mat)), , drop = FALSE]

geneSets_all$NuS_Up <- intersect(geneSets_all$NuS_Up, rownames(expr.mat))
geneSets_all$NuS_Down <- intersect(geneSets_all$NuS_Down, rownames(expr.mat))

cat("===== Gene set overlap =====\n")
cat("NuS_Up:", length(geneSets_all$NuS_Up), "\n")
cat("NuS_Down:", length(geneSets_all$NuS_Down), "\n")

if (length(geneSets_all$NuS_Up) < 5 || length(geneSets_all$NuS_Down) < 5) {
  stop("Too few overlapping genes for ssGSEA.")
}

############################################################
## Step 6. ssGSEA
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
  stop("ssGSEA result does not contain NuS_Up or NuS_Down.")
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
## Step 7. 统计分析：one-way ANOVA + Dunnett vs Control
############################################################
library(multcomp)

plot_data$Group <- factor(
  plot_data$Group,
  levels = c("Control", "24h_5FU", "48h_5FU")
)
control_group <- "Control"

cat("===== Group sizes =====\n")
print(table(plot_data$Group))

## one-way ANOVA
fit <- aov(Score ~ Group, data = plot_data)
anova_p <- tryCatch(
  summary(fit)[[1]][["Pr(>F)"]][1],
  error = function(e) NA_real_
)

cat("===== ANOVA result =====\n")
print(summary(fit))

## group summary
mean_sd_df <- plot_data %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    N = sum(!is.na(Score)),
    .groups = "drop"
  )

write.csv(
  mean_sd_df,
  paste0(geo_id, "_NuS_group_summary.csv"),
  row.names = FALSE
)

## Dunnett vs Control
dunnett_res <- tryCatch(
  summary(glht(fit, linfct = mcp(Group = "Dunnett"))),
  error = function(e) NULL
)

if (is.null(dunnett_res)) {
  dunnett_df <- data.frame(
    comparison = NA_character_,
    treatment = NA_character_,
    diff = NA_real_,
    p.value = NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  comp_names <- names(dunnett_res$test$coefficients)
  comp_diff  <- as.numeric(dunnett_res$test$coefficients)
  comp_p     <- as.numeric(dunnett_res$test$pvalues)

  dunnett_df <- data.frame(
    comparison = comp_names,
    treatment = sub(" - Control$", "", comp_names),
    diff = comp_diff,
    p.value = comp_p,
    stringsAsFactors = FALSE
  )
}

write.csv(
  data.frame(method = "one-way ANOVA", p.value = anova_p),
  paste0(geo_id, "_NuS_all_overall_ANOVA.csv"),
  row.names = FALSE
)

write.csv(
  dunnett_df,
  paste0(geo_id, "_NuS_all_Dunnett_vs_Control.csv"),
  row.names = FALSE
)

cat("===== Dunnett result =====\n")
print(dunnett_df)

############################################################
## Step 8. 小提琴图
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "24h_5FU" = "#93D4BC",
  "48h_5FU" = "#08589E"
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
  ) +
  stat_compare_means(
    method = "anova",
    label = "p.format",
    label.x.npc = "center",
    size = 5
  )

print(p1)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.pdf"),
  plot = p1,
  width = 3.2,
  height = 4,
  bg = "white"
)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.png"),
  plot = p1,
  width = 3.2,
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
  file = paste0(geo_id, "_NuS_all_workspace.Rdata")
)

cat("分析完成。\n")
cat("结果文件保存在:\n")
cat(getwd(), "\n")
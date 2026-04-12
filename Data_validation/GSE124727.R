############################################################
## GSE124727
## 完整修正版：表达矩阵构建 + ssGSEA打分 + 绘图
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE124727/")

library(readr)
library(dplyr)
library(stringr)
library(GSVA)
library(GSEABase)
library(ggplot2)

geo_id <- "GSE124727"

############################################################
## Step 1. 读取 GEO series matrix
############################################################
data_raw <- read_table(
  "GSE124727_series_matrix.txt.gz",
  comment = "!",
  show_col_types = FALSE
)

cat("===== Raw matrix dimension =====\n")
print(dim(data_raw))
cat("===== Raw column names =====\n")
print(colnames(data_raw))

############################################################
## Step 2. 读取平台注释
############################################################
ids_raw <- read_tsv(
  "GPL15207-17536.txt",
  comment = "#",
  show_col_types = FALSE
)

ids <- data.frame(
  probe_id = ids_raw$ID,
  symbol = ids_raw$`Gene Symbol`,
  stringsAsFactors = FALSE
)

ids <- ids[
  !is.na(ids$symbol) &
    ids$symbol != "" &
    !grepl("//", ids$symbol, fixed = TRUE),
]

cat("===== Valid annotation rows =====\n")
print(nrow(ids))
cat("===== Unique symbols =====\n")
print(length(unique(ids$symbol)))

############################################################
## Step 3. 构建 probe-level 表达矩阵
############################################################
data_raw <- data_raw %>%
  mutate(ID_REF = str_replace_all(`"ID_REF"`, '"', ""))

## 只保留真正的 GSM 样本列
sample_cols <- grep("^\"GSM", colnames(data_raw), value = TRUE)

expr_probe <- data_raw[, sample_cols, drop = FALSE]
expr_probe <- as.data.frame(expr_probe, check.names = FALSE)
rownames(expr_probe) <- data_raw$ID_REF
colnames(expr_probe) <- gsub('"', "", colnames(expr_probe))

expr_probe <- as.matrix(expr_probe)
mode(expr_probe) <- "numeric"

cat("===== Probe-level matrix dimension =====\n")
print(dim(expr_probe))
cat("===== Probe-level sample names =====\n")
print(colnames(expr_probe))

############################################################
## Step 4. 保留可注释 probe，并合并重复 probe
############################################################
ids <- ids[ids$probe_id %in% rownames(expr_probe), , drop = FALSE]
expr_probe <- expr_probe[ids$probe_id, , drop = FALSE]

cat("===== Matched probe matrix dimension =====\n")
print(dim(expr_probe))

## 计算每个 probe 的中位表达
ids$median_expr <- apply(expr_probe, 1, median, na.rm = TRUE)

## 按 symbol + 中位表达排序，保留每个 symbol 表达最高的 probe
ids <- ids[order(ids$symbol, -ids$median_expr), ]
ids <- ids[!duplicated(ids$symbol), ]

expr_gene <- expr_probe[ids$probe_id, , drop = FALSE]
rownames(expr_gene) <- ids$symbol

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))
cat("===== Gene-level sample names =====\n")
print(colnames(expr_gene))

save(expr_gene, ids, expr_probe, file = paste0(geo_id, "_step1_workspace.Rdata"))
write.csv(expr_gene, file = "expressionmetrix_GSE.csv")

############################################################
## Step 5. ssGSEA打分
############################################################
expr <- expr_gene
mode(expr) <- "numeric"

cat("===== Expression range before transform =====\n")
print(range(expr, na.rm = TRUE))

## log2 transform
expr <- log2(expr + 1)

cat("===== Expression range after log2(x+1) =====\n")
print(range(expr, na.rm = TRUE))

## 去除包含NA的基因
expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]

cat("===== Final expression matrix dimension =====\n")
print(dim(expr))
cat("===== Final sample names =====\n")
print(colnames(expr))

############################################################
## Step 6. 根据真实样本名建立分组信息
############################################################
group_map <- data.frame(
  Sample = c("GSM3544272", "GSM3544273", "GSM3544274", "GSM3544275"),
  Group = c("Control", "MPA", "Guanosine", "MPA + Guanosine"),
  stringsAsFactors = FALSE
)

group_info <- group_map[match(colnames(expr), group_map$Sample), , drop = FALSE]

if (any(is.na(group_info$Group))) {
  stop("Some samples in expr were not matched to group_map.")
}

group_info$Group <- factor(
  group_info$Group,
  levels = c("Control", "MPA", "Guanosine", "MPA + Guanosine")
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## Step 7. 读取基因集
############################################################
load("/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata")

geneSets_all <- list(
  NuS_Up = unique(toupper(NuS_up_all_unique)),
  NuS_Down = unique(toupper(NuS_down_all_unique))
)

geneSets_core <- list(
  NuS_Up = unique(toupper(NuS_up_core_unique)),
  NuS_Down = unique(toupper(NuS_down_core_unique))
)

## 这里选择使用 core gene set
geneSets_use <- geneSets_all
geneset_label <- "all"

############################################################
## Step 8. 与表达矩阵取交集
############################################################
rownames(expr) <- toupper(rownames(expr))

geneSets_use$NuS_Up <- intersect(geneSets_use$NuS_Up, rownames(expr))
geneSets_use$NuS_Down <- intersect(geneSets_use$NuS_Down, rownames(expr))

cat("===== Overlap with expression matrix =====\n")
cat("NuS_Up:", length(geneSets_use$NuS_Up), "\n")
cat("NuS_Down:", length(geneSets_use$NuS_Down), "\n")

if (length(geneSets_use$NuS_Up) < 5 || length(geneSets_use$NuS_Down) < 5) {
  stop("Too few overlapping genes for ssGSEA.")
}

############################################################
## Step 9. ssGSEA
############################################################
ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr,
  geneSets = geneSets_use,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_res <- gsva(ssgsea_param)

NuS_score_ssGSEA <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]

scores <- data.frame(
  Sample = colnames(expr),
  NuS_score = as.numeric(NuS_score_ssGSEA),
  stringsAsFactors = FALSE
)

scores <- merge(scores, group_info, by = "Sample")
scores$Group <- factor(scores$Group, levels = levels(group_info$Group))

cat("===== NuS scores =====\n")
print(scores)

write.csv(scores, paste0(geo_id, "_", geneset_label, "_NuS_scores.csv"), row.names = FALSE)

############################################################
## Step 10. 绘图（柱状图 + SE）
## 使用固定浅绿-浅青-蓝色风格
############################################################
summary_data <- scores %>%
  group_by(Group) %>%
  summarise(
    Mean_Score = mean(NuS_score, na.rm = TRUE),
    SE = sd(NuS_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

fill_colors <- c(
  "Control" = "#C4E7C1",
  "MPA" = "#93D4BC",
  "Guanosine" = "#51B3D1",
  "MPA + Guanosine" = "#08589E"
)

p <- ggplot(summary_data, aes(x = Group, y = Mean_Score, fill = Group)) +
  geom_col(alpha = 0.85, width = 0.65, color = "black", linewidth = 0.8) +
  geom_errorbar(
    aes(ymin = Mean_Score - SE, ymax = Mean_Score + SE),
    width = 0.2,
    color = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = round(Mean_Score, 2)),
    vjust = -1.2,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = paste0("NuS Score (ssGSEA) - ", geo_id, " - ", geneset_label),
    x = NULL,
    y = "Mean NuS Score"
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

print(p)

ggsave(
  paste0(geo_id, "_", geneset_label, "_NuS_barplot.pdf"),
  p,
  width = 4,
  height = 4,
  bg = "white"
)

ggsave(
  paste0(geo_id, "_", geneset_label, "_NuS_barplot.png"),
  p,
  width = 4,
  height = 4,
  dpi = 300,
  bg = "white"
)
############################################################
## Step 11. 保存workspace
############################################################
save(
  expr,
  expr_gene,
  ids,
  group_info,
  geneSets_use,
  ssgsea_res,
  scores,
  summary_data,
  file = paste0(geo_id, "_", geneset_label, "_NuS_workspace.Rdata")
)
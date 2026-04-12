rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 工作路径与输出路径
############################################################
setwd("/home/xxm_xxm/CJX_workspace/proteinomics/")
out_dir <- file.path(getwd(), "NuS_RiboSis_ssGSEA_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## 最终版 NuStress gene set 路径
nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"

## RiboSis 基因集路径
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

############################################################
## 1. 加载 R 包
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
## 2. 读取 matrix.xlsx
############################################################
expr_df <- read_xlsx("matrix.xlsx")

cat("===== Input matrix dimension =====\n")
print(dim(expr_df))
cat("===== Column names =====\n")
print(colnames(expr_df))

stopifnot(ncol(expr_df) >= 2)

############################################################
## 3. 构建表达矩阵
############################################################
expr_mat <- as.matrix(expr_df[, 2:ncol(expr_df), drop = FALSE])
rownames(expr_mat) <- toupper(trimws(as.character(expr_df[[1]])))
mode(expr_mat) <- "numeric"

## 去掉空基因名
valid_gene <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
expr_mat <- expr_mat[valid_gene, , drop = FALSE]

## 去重：保留首次出现
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

cat("===== Expression matrix dimension after cleaning =====\n")
print(dim(expr_mat))
print(colnames(expr_mat))

############################################################
## 4. 分组信息
############################################################
## 如分组与这里不一致，只修改这一段
group_info <- data.frame(
  Sample = colnames(expr_mat),
  Group = rep(c("Control", "OXA", "L", "L-OXA"), each = 2),
  stringsAsFactors = FALSE
)

group_info$Group <- factor(group_info$Group, levels = unique(group_info$Group))
print(group_info)

############################################################
## 5. 组内均值填充 NA
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

## 去掉仍然有 NA 的基因
expr_mat <- expr_mat[complete.cases(expr_mat), , drop = FALSE]

cat("===== Matrix dimension after NA removal =====\n")
print(dim(expr_mat))

############################################################
## 6. 读取最终版 NuStress gene sets
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

## 与表达矩阵取交集
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

## 导出 overlap / missing 方便核对
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
## 7. 读取 RiboSis gene set
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
## 8. ssGSEA 评分函数
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
## 9. 计算 NuStress Full / Core / RiboSis 评分
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
## 10. 汇总得分表
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
## 11. 统计分析函数
## four-group one-way ANOVA + TukeyHSD
############################################################
run_group_stats <- function(df) {
  df <- as.data.frame(df)
  df$Group <- droplevels(factor(df$Group))
  groups <- levels(df$Group)
  
  summary_df <- df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      sd_score   = sd(Score, na.rm = TRUE),
      n          = sum(!is.na(Score)),
      .groups = "drop"
    )
  
  if (length(groups) < 2) {
    return(data.frame(
      comparison = NA_character_,
      anova_p = NA_real_,
      diff_mean = NA_real_,
      lwr = NA_real_,
      upr = NA_real_,
      p_adj = NA_real_,
      group_high = NA_character_,
      group_low = NA_character_,
      mean_high = NA_real_,
      mean_low = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  ## only two groups
  if (length(groups) == 2) {
    g1 <- groups[1]
    g2 <- groups[2]
    
    tt <- t.test(Score ~ Group, data = df, var.equal = FALSE)
    m1 <- summary_df$mean_score[match(g1, summary_df$Group)]
    m2 <- summary_df$mean_score[match(g2, summary_df$Group)]
    
    return(data.frame(
      comparison = paste(g2, "-", g1),
      anova_p = NA_real_,
      diff_mean = m2 - m1,
      lwr = NA_real_,
      upr = NA_real_,
      p_adj = tt$p.value,
      group_high = g2,
      group_low = g1,
      mean_high = m2,
      mean_low = m1,
      stringsAsFactors = FALSE
    ))
  }
  
  ## one-way ANOVA
  fit <- aov(Score ~ Group, data = df)
  anova_p <- tryCatch(
    summary(fit)[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  ## Tukey HSD pairwise comparison
  tuk <- TukeyHSD(fit)$Group
  tuk_df <- as.data.frame(tuk, stringsAsFactors = FALSE)
  tuk_df$comparison <- rownames(tuk_df)
  rownames(tuk_df) <- NULL
  
  ## parse comparison name: A-B means mean(A) - mean(B)
  comp_split <- strsplit(tuk_df$comparison, "-")
  tuk_df$group_high <- sapply(comp_split, `[`, 1)
  tuk_df$group_low  <- sapply(comp_split, `[`, 2)
  
  tuk_df$mean_high <- summary_df$mean_score[match(tuk_df$group_high, summary_df$Group)]
  tuk_df$mean_low  <- summary_df$mean_score[match(tuk_df$group_low, summary_df$Group)]
  
  out <- data.frame(
    comparison = tuk_df$comparison,
    anova_p = anova_p,
    diff_mean = tuk_df$diff,
    lwr = tuk_df$lwr,
    upr = tuk_df$upr,
    p_adj = tuk_df$`p adj`,
    group_high = tuk_df$group_high,
    group_low = tuk_df$group_low,
    mean_high = tuk_df$mean_high,
    mean_low = tuk_df$mean_low,
    stringsAsFactors = FALSE
  )
  
  return(out)
}

stats_all <- scores_all %>%
  dplyr::group_by(ScoreType) %>%
  dplyr::group_modify(~run_group_stats(.x)) %>%
  dplyr::ungroup()

write.csv(
  stats_all,
  file.path(out_dir, "All_ssGSEA_statistics.csv"),
  row.names = FALSE
)

############################################################
## 12. 颜色与绘图函数
############################################################
base_fill_colors <- c(
  "Control" = "#C4E7C1",
  "OXA"     = "#93D4BC",
  "L"       = "#51B3D1",
  "L-OXA"   = "#08589E"
)

make_group_colors <- function(groups) {
  default_palette <- c("#C4E7C1", "#93D4BC", "#51B3D1", "#08589E",
                       "#2171B5", "#6BAED6", "#9ECAE1", "#C6DBEF")
  cols <- default_palette[seq_along(groups)]
  names(cols) <- groups
  
  overlap_names <- intersect(names(base_fill_colors), groups)
  cols[overlap_names] <- base_fill_colors[overlap_names]
  return(cols)
}

plot_violin <- function(df, score_type, out_dir) {
  plot_df <- df %>% dplyr::filter(ScoreType == score_type)
  fill_colors <- make_group_colors(levels(plot_df$Group))
  
  p <- ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
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
      title = score_type,
      x = NULL,
      y = "ssGSEA score"
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
  
  ggsave(file.path(out_dir, paste0(score_type, "_violinplot.pdf")),
         p, width = 4.8, height = 4.2, bg = "white")
  ggsave(file.path(out_dir, paste0(score_type, "_violinplot.png")),
         p, width = 4.8, height = 4.2, dpi = 300, bg = "white")
  
  return(p)
}

p1 <- plot_violin(scores_all, "NuS_Full", out_dir)
if (has_core) p2 <- plot_violin(scores_all, "NuS_Core", out_dir)
p3 <- plot_violin(scores_all, "RiboSis", out_dir)

############################################################
## 13. 合并图
############################################################
fill_colors_all <- make_group_colors(levels(scores_all$Group))

p_all <- ggplot(scores_all, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.8, alpha = 0.85) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +
  facet_wrap(~ ScoreType, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = fill_colors_all) +
  labs(
    title = "NuS and RiboSis ssGSEA scores",
    x = NULL,
    y = "ssGSEA score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(out_dir, "Combined_ssGSEA_violinplot.pdf"),
       p_all, width = 10, height = 4.5, bg = "white")
ggsave(file.path(out_dir, "Combined_ssGSEA_violinplot.png"),
       p_all, width = 10, height = 4.5, dpi = 300, bg = "white")

############################################################
## 14. NuS 与 RiboSis 相关性分析
############################################################
score_wide <- scores_all %>%
  dplyr::select(Sample, Group, ScoreType, Score) %>%
  tidyr::pivot_wider(names_from = ScoreType, values_from = Score)

## Full vs RiboSis
if (all(c("NuS_Full", "RiboSis") %in% colnames(score_wide))) {
  cor_full <- cor.test(score_wide$NuS_Full, score_wide$RiboSis, method = "pearson")
  
  p_cor_full <- ggplot(score_wide, aes(x = NuS_Full, y = RiboSis, color = Group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = make_group_colors(levels(score_wide$Group))) +
    labs(
      title = paste0("NuS_Full vs RiboSis (r = ", round(cor_full$estimate, 3), 
                     ", p = ", signif(cor_full$p.value, 3), ")"),
      x = "NuS_Full score",
      y = "RiboSis score"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")
  
  ggsave(file.path(out_dir, "Correlation_NuS_Full_vs_RiboSis.pdf"),
         p_cor_full, width = 6, height = 5, bg = "white")
}

## Core vs RiboSis
if (has_core && all(c("NuS_Core", "RiboSis") %in% colnames(score_wide))) {
  cor_core <- cor.test(score_wide$NuS_Core, score_wide$RiboSis, method = "pearson")
  
  p_cor_core <- ggplot(score_wide, aes(x = NuS_Core, y = RiboSis, color = Group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = make_group_colors(levels(score_wide$Group))) +
    labs(
      title = paste0("NuS_Core vs RiboSis (r = ", round(cor_core$estimate, 3), 
                     ", p = ", signif(cor_core$p.value, 3), ")"),
      x = "NuS_Core score",
      y = "RiboSis score"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")
  
  ggsave(file.path(out_dir, "Correlation_NuS_Core_vs_RiboSis.pdf"),
         p_cor_core, width = 6, height = 5, bg = "white")
}

############################################################
## 15. 保存 workspace
############################################################
if (has_core) {
  save(
    expr_mat,
    group_info,
    geneSets_full,
    geneSets_core,
    geneSets_ribosis,
    scores_all,
    stats_all,
    file = file.path(out_dir, "NuS_RiboSis_ssGSEA_workspace.Rdata")
  )
} else {
  save(
    expr_mat,
    group_info,
    geneSets_full,
    geneSets_ribosis,
    scores_all,
    stats_all,
    file = file.path(out_dir, "NuS_RiboSis_ssGSEA_workspace.Rdata")
  )
}

cat("\n===== Analysis completed =====\n")
cat("Results saved in:\n", out_dir, "\n")
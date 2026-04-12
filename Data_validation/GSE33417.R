rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. 路径设置
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE33417"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE33417"
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
## 2. 读取表达矩阵
############################################################
load(file.path(data_dir, "GSE33417_step1_workspace.Rdata"))

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

cat("===== All sample names =====\n")
print(colnames(expr_gene))

############################################################
## 3. 构建样本注释信息（根据GSM样本名）
############################################################
sample_info <- data.frame(
  Sample = c(
    "GSM826632","GSM826633","GSM826634","GSM826635","GSM826636","GSM826637",
    "GSM826638","GSM826639","GSM826640","GSM826641","GSM826642","GSM826643",
    "GSM826644","GSM826645","GSM826646","GSM826647","GSM826648","GSM826649",
    "GSM826650","GSM826651","GSM826652","GSM826653","GSM826654","GSM826655",
    "GSM826656","GSM826657","GSM826658","GSM826659","GSM826660","GSM826661",
    "GSM826662","GSM826663","GSM826664","GSM826665","GSM826666","GSM826667"
  ),
  CellType = c(
    rep("iPS", 18),
    rep("HFF", 18)
  ),
  Time = c(
    rep(c("0","15","30","60","120","240"), each = 3),
    rep(c("0","15","30","60","120","240"), each = 3)
  ),
  Rep = rep(c("1","2","3"), 12),
  stringsAsFactors = FALSE
)

sample_info$Time <- factor(sample_info$Time, levels = c("0","15","30","60","120","240"))
sample_info$CellType <- factor(sample_info$CellType, levels = c("iPS", "HFF"))

## 按表达矩阵顺序重排
sample_info <- sample_info[match(colnames(expr_gene), sample_info$Sample), ]

if (any(is.na(sample_info$Sample))) {
  stop("Some expression matrix sample names are not found in sample_info.")
}

cat("===== Sample info =====\n")
print(sample_info)

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
rownames(expr_gene) <- toupper(rownames(expr_gene))

full_up_in   <- intersect(full_up, rownames(expr_gene))
full_down_in <- intersect(full_down, rownames(expr_gene))

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
  exprData = as.matrix(expr_gene),
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
  Sample = colnames(expr_gene),
  Score = as.numeric(nus_score),
  Method = "ssGSEA",
  GeneSet = "Full",
  stringsAsFactors = FALSE
)

scores_all <- left_join(scores_all, sample_info, by = "Sample")
scores_all$Time <- factor(scores_all$Time, levels = c("0","15","30","60","120","240"))
scores_all$CellType <- factor(scores_all$CellType, levels = c("iPS","HFF"))

write.csv(scores_all, "GSE33417_NuS_scores_ssGSEA_all_samples.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 7. 每个细胞类型分别做统计：one-way ANOVA + Dunnett vs 0 min
############################################################
run_time_statistics <- function(dat, prefix) {
  dat$Time <- droplevels(dat$Time)
  dat$Time <- relevel(dat$Time, ref = "0")
  
  fit <- aov(Score ~ Time, data = dat)
  anova_p <- tryCatch(
    summary(fit)[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  mean_sd_df <- dat %>%
    group_by(Time) %>%
    summarise(
      Mean = mean(Score, na.rm = TRUE),
      SD = sd(Score, na.rm = TRUE),
      N = sum(!is.na(Score)),
      .groups = "drop"
    )
  
  control_row <- mean_sd_df %>% filter(Time == "0")
  
  dunnett_res <- tryCatch(
    summary(glht(fit, linfct = mcp(Time = "Dunnett"))),
    error = function(e) NULL
  )
  
  if (is.null(dunnett_res)) {
    out <- data.frame(
      comparison = NA_character_,
      treatment = NA_character_,
      p_anova = anova_p,
      p_dunnett = NA_real_,
      mean_0min = control_row$Mean,
      sd_0min = control_row$SD,
      n_0min = control_row$N,
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
    
    treat_names <- sub(" - 0$", "", comp_names)
    treat_df <- mean_sd_df %>% filter(Time != "0") %>% as.data.frame()
    match_idx <- match(treat_names, as.character(treat_df$Time))
    
    out <- data.frame(
      comparison = comp_names,
      treatment = treat_names,
      p_anova = rep(anova_p, length(comp_names)),
      p_dunnett = comp_p,
      mean_0min = rep(control_row$Mean, length(comp_names)),
      sd_0min = rep(control_row$SD, length(comp_names)),
      n_0min = rep(control_row$N, length(comp_names)),
      mean_treat = treat_df$Mean[match_idx],
      sd_treat = treat_df$SD[match_idx],
      n_treat = treat_df$N[match_idx],
      diff_mean = comp_diff,
      stringsAsFactors = FALSE
    )
  }
  
  write.csv(out, paste0(prefix, "_NuS_score_statistics_ANOVA_Dunnett.csv"), row.names = FALSE)
  return(out)
}

scores_ips <- scores_all %>% filter(CellType == "iPS")
scores_hff <- scores_all %>% filter(CellType == "HFF")

stat_ips <- run_time_statistics(scores_ips, "GSE33417_iPS")
stat_hff <- run_time_statistics(scores_hff, "GSE33417_HFF")

############################################################
## 8. 分别绘图：iPS 一张，HFF 一张
############################################################
fill_colors <- c(
  "0"   = "#C4E7C1",
  "15"  = "#B7E0C6",
  "30"  = "#93D4BC",
  "60"  = "#6EC5C6",
  "120" = "#51B3D1",
  "240" = "#08589E"
)

plot_one_celltype <- function(dat, prefix, title_text) {
  p <- ggplot(dat, aes(x = Time, y = Score, fill = Time)) +
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
      title = title_text,
      x = "Time (min)",
      y = "NuS score"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  ggsave(paste0(prefix, "_NuS_ssGSEA_violinplot.pdf"), p, width = 5.5, height = 5, bg = "white")
  ggsave(paste0(prefix, "_NuS_ssGSEA_violinplot.png"), p, width = 5.5, height = 5, dpi = 300, bg = "white")
}

plot_one_celltype(scores_ips, "GSE33417_iPS", "NuS ssGSEA scores in GSE33417 iPS")
plot_one_celltype(scores_hff, "GSE33417_HFF", "NuS ssGSEA scores in GSE33417 HFF")

############################################################
## 9. 保存workspace
############################################################
save(
  expr_gene,
  sample_info,
  geneSets_full,
  scores_all,
  scores_ips,
  scores_hff,
  stat_ips,
  stat_hff,
  file = "GSE33417_NuS_ssGSEA_workspace.Rdata"
)
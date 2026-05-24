rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE198178"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE198178"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)

############################################################
## 2. Load the expression matrix generated in Step 1
############################################################
load(file.path(data_dir, "GSE198178_step1_workspace.Rdata"))

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

cat("===== All sample names =====\n")
print(colnames(tpm_final))

############################################################
## 3. Read the full gene set
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
## 4. General scoring function:
##    return scores_all and test_res
############################################################
run_nus_ssgsea_tpm <- function(expr_sub, group_info, prefix) {
  
  expr_sub <- as.matrix(expr_sub)
  mode(expr_sub) <- "numeric"
  
  cat("\n=============================\n")
  cat("Running:", prefix, "\n")
  cat("=============================\n")
  print(colnames(expr_sub))
  print(group_info)
  
  ## TPM filter
  keep <- rowSums(expr_sub > 1, na.rm = TRUE) >= 3
  expr_sub <- expr_sub[keep, , drop = FALSE]
  
  ## log2 transform
  expr_sub <- log2(expr_sub + 1)
  
  ## gene symbol upper
  rownames(expr_sub) <- toupper(rownames(expr_sub))
  
  ## overlap
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
  
  ## ssGSEA
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
  scores_all$Group <- factor(scores_all$Group, levels = c("Control", "ACTD"))
  
  test_res <- data.frame(
    comparison = paste0(prefix, ": ACTD vs Control"),
    p_welch = tryCatch(
      t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    mean_control = mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
    sd_control   = sd(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
    n_control    = sum(scores_all$Group == "Control"),
    mean_treat   = mean(scores_all$Score[scores_all$Group == "ACTD"], na.rm = TRUE),
    sd_treat     = sd(scores_all$Score[scores_all$Group == "ACTD"], na.rm = TRUE),
    n_treat      = sum(scores_all$Group == "ACTD"),
    diff_mean    = mean(scores_all$Score[scores_all$Group == "ACTD"], na.rm = TRUE) -
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
  
  return(list(scores_all = scores_all, test_res = test_res))
}

############################################################
## 5. 1 h
############################################################
dat_1h <- tpm_final[, c(1, 5, 9, 2, 6, 10), drop = FALSE]

group_info_1h <- data.frame(
  Sample = colnames(dat_1h),
  Group = c(rep("Control", 3), rep("ACTD", 3)),
  Time = "1h",
  stringsAsFactors = FALSE
)

res_1h <- run_nus_ssgsea_tpm(
  expr_sub = dat_1h,
  group_info = group_info_1h,
  prefix = "GSE198178_1h"
)

############################################################
## 6. 2 h
############################################################
dat_2h <- tpm_final[, c(1, 5, 9, 3, 7, 11), drop = FALSE]

group_info_2h <- data.frame(
  Sample = colnames(dat_2h),
  Group = c(rep("Control", 3), rep("ACTD", 3)),
  Time = "2h",
  stringsAsFactors = FALSE
)

res_2h <- run_nus_ssgsea_tpm(
  expr_sub = dat_2h,
  group_info = group_info_2h,
  prefix = "GSE198178_2h"
)

############################################################
## 7. 4 h
############################################################
dat_4h <- tpm_final[, c(1, 5, 9, 4, 8, 12), drop = FALSE]

group_info_4h <- data.frame(
  Sample = colnames(dat_4h),
  Group = c(rep("Control", 3), rep("ACTD", 3)),
  Time = "4h",
  stringsAsFactors = FALSE
)

res_4h <- run_nus_ssgsea_tpm(
  expr_sub = dat_4h,
  group_info = group_info_4h,
  prefix = "GSE198178_4h"
)

############################################################
## 8. Merge score results from all three time points
############################################################
scores_plot <- bind_rows(
  res_1h$scores_all,
  res_2h$scores_all,
  res_4h$scores_all
)

scores_plot$Time <- factor(scores_plot$Time, levels = c("1h", "2h", "4h"))

scores_plot$Group_Time <- factor(
  paste(scores_plot$Group, scores_plot$Time, sep = "_"),
  levels = c(
    "Control_1h", "ACTD_1h",
    "Control_2h", "ACTD_2h",
    "Control_4h", "ACTD_4h"
  )
)

write.csv(scores_plot, "GSE198178_NuS_scores_ssGSEA_all_groups_for_plot.csv", row.names = FALSE)

test_all <- bind_rows(
  res_1h$test_res,
  res_2h$test_res,
  res_4h$test_res
)
write.csv(test_all, "GSE198178_NuS_score_statistics_all_groups_ttest.csv", row.names = FALSE)

############################################################
## 9. Plot all groups in a single figure
############################################################
fill_colors <- c(
  "Control_1h" = "#C4E7C1",
  "ACTD_1h"    = "#93D4BC",
  "Control_2h" = "#B7DCE8",
  "ACTD_2h"    = "#51B3D1",
  "Control_4h" = "#A9C7E8",
  "ACTD_4h"    = "#08589E"
)

p <- ggplot(scores_plot, aes(x = Group_Time, y = Score, fill = Group_Time)) +
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
    title = "NuS ssGSEA scores in GSE198178",
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

ggsave("GSE198178_NuS_ssGSEA_all_groups_oneplot.pdf", p, width = 6, height = 5, bg = "white")
ggsave("GSE198178_NuS_ssGSEA_all_groups_oneplot.png", p, width = 6, height = 5, dpi = 300, bg = "white")

############################################################
## 10. Save the combined workspace
############################################################
save(
  tpm_final,
  dat_1h,
  dat_2h,
  dat_4h,
  res_1h,
  res_2h,
  res_4h,
  scores_plot,
  test_all,
  file = "GSE198178_NuS_ssGSEA_all_groups_oneplot_workspace.Rdata"
)
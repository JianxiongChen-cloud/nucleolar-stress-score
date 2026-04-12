############################################################
## GSE254323
## Final polished version
## raw merge + gene symbol annotation + NuS ssGSEA + ANOVA + Dunnett vs DMSO
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(multcomp)

geo_id <- "GSE254323"
base_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE254323"
raw_dir <- file.path(base_dir, "GSE254323_RAW")
geneset_file <- "/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
setwd(base_dir)

############################################################
## Step 1. Merge raw files
############################################################
setwd(raw_dir)

file_list <- list.files(pattern = "\\.(csv|gz)$", full.names = TRUE)
cat("===== Found files =====\n")
print(file_list)

merged_data <- NULL

for (file in file_list) {
  temp_data <- fread(file, header = TRUE)
  
  if ("geneid" %in% colnames(temp_data)) {
    setkey(temp_data, geneid)
  } else if ("gene_id" %in% colnames(temp_data)) {
    setnames(temp_data, "gene_id", "geneid")
    setkey(temp_data, geneid)
  } else if (colnames(temp_data)[1] == "V1") {
    setnames(temp_data, "V1", "geneid")
    setkey(temp_data, geneid)
  } else {
    stop(paste("Cannot identify gene id column in:", file))
  }
  
  if (is.null(merged_data)) {
    merged_data <- temp_data
  } else {
    merged_data <- merge(merged_data, temp_data, by = "geneid", all = TRUE)
  }
  
  cat("Processed:", basename(file), "-", nrow(temp_data), "rows", ncol(temp_data), "cols\n")
}

cat("===== merged_data dimension =====\n")
print(dim(merged_data))
print(head(merged_data))

############################################################
## Step 2. Ensembl -> gene symbol
############################################################
ensembl_ids <- sub("\\..*", "", merged_data$geneid)

symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

conversion_result <- data.frame(
  geneid = merged_data$geneid,
  ensembl_id = ensembl_ids,
  gene_symbol = symbols,
  stringsAsFactors = FALSE
)

na_count <- sum(is.na(conversion_result$gene_symbol))
cat("===== Unmapped genes =====\n")
print(na_count)

############################################################
## Step 3. Build gene-level expression matrix
############################################################
expr_df <- as.data.frame(merged_data, check.names = FALSE)
expr_df$gene_symbol <- conversion_result$gene_symbol

## keep gene_symbol + expression columns
expr_df <- expr_df[, c("gene_symbol", setdiff(colnames(expr_df), "geneid")), drop = FALSE]

## remove NA symbols
expr_df <- expr_df[!is.na(expr_df$gene_symbol) & expr_df$gene_symbol != "", , drop = FALSE]

## convert expression columns to numeric
expr_cols <- setdiff(colnames(expr_df), "gene_symbol")
expr_df[, expr_cols] <- lapply(expr_df[, expr_cols, drop = FALSE], as.numeric)

## keep highest median expression for duplicated symbols
expr_df$median_expr <- apply(expr_df[, expr_cols, drop = FALSE], 1, median, na.rm = TRUE)
expr_df <- expr_df[order(expr_df$gene_symbol, -expr_df$median_expr), ]
expr_df <- expr_df[!duplicated(expr_df$gene_symbol), ]

expr_gene <- as.matrix(expr_df[, expr_cols, drop = FALSE])
mode(expr_gene) <- "numeric"
rownames(expr_gene) <- expr_df$gene_symbol

cat("===== expr_gene dimension =====\n")
print(dim(expr_gene))
print(range(expr_gene, na.rm = TRUE))

setwd(base_dir)
save(expr_gene, file = paste0(geo_id, "_step1_workspace.Rdata"))
write.csv(expr_gene, file = "expressionmetrix_GSE.csv")

############################################################
## Step 4. Load full gene set
############################################################
load(geneset_file)

geneSets_all <- list(
  NuS_Up = unique(toupper(NuS_up_all_unique)),
  NuS_Down = unique(toupper(NuS_down_all_unique))
)

############################################################
## Step 5. Analysis function
############################################################
run_nus_analysis <- function(expr_input, group_info, dataset_prefix, control_group) {
  
  expr_input <- as.matrix(expr_input)
  mode(expr_input) <- "numeric"
  
  ## align sample order
  expr_input <- expr_input[, group_info$Sample, drop = FALSE]
  
  ## log2 transform
  expr_mat <- log2(expr_input + 1)
  expr_mat <- expr_mat[rowSums(is.na(expr_mat)) == 0, , drop = FALSE]
  
  ## uppercase rownames and deduplicate
  rownames(expr_mat) <- toupper(trimws(rownames(expr_mat)))
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  
  ## gene set overlap
  geneSets_use <- geneSets_all
  geneSets_use$NuS_Up <- intersect(geneSets_use$NuS_Up, rownames(expr_mat))
  geneSets_use$NuS_Down <- intersect(geneSets_use$NuS_Down, rownames(expr_mat))
  
  cat("===== ", dataset_prefix, " overlap =====\n", sep = "")
  cat("NuS_Up:", length(geneSets_use$NuS_Up), "\n")
  cat("NuS_Down:", length(geneSets_use$NuS_Down), "\n")
  
  if (length(geneSets_use$NuS_Up) < 5 || length(geneSets_use$NuS_Down) < 5) {
    stop(paste("Too few overlapping genes for", dataset_prefix))
  }
  
  ## ssGSEA
  ssgsea_param <- GSVA::ssgseaParam(
    exprData = expr_mat,
    geneSets = geneSets_use,
    minSize = 5,
    maxSize = 500,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_res <- gsva(ssgsea_param)
  
  if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea_res))) {
    stop(paste("NuS_Up or NuS_Down not found in ssGSEA result for", dataset_prefix))
  }
  
  nus_score <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]
  
  scores <- data.frame(
    Sample = colnames(expr_mat),
    NuS_Up = as.numeric(ssgsea_res["NuS_Up", ]),
    NuS_Down = as.numeric(ssgsea_res["NuS_Down", ]),
    NuS_Score = as.numeric(nus_score),
    stringsAsFactors = FALSE
  )
  
  scores <- merge(scores, group_info, by = "Sample")
  
  ## IMPORTANT: make DMSO the reference level
  group_levels <- unique(group_info$Group)
  group_levels <- c(control_group, setdiff(group_levels, control_group))
  
  scores$Group <- factor(scores$Group, levels = group_levels)
  scores$Group <- relevel(scores$Group, ref = control_group)
  
  plot_data <- data.frame(
    Score = scores$NuS_Score,
    Group = scores$Group
  )
  
  ## ANOVA
  fit <- aov(Score ~ Group, data = plot_data)
  anova_summary <- summary(fit)
  anova_p <- anova_summary[[1]][["Pr(>F)"]][1]
  
  ## Tukey HSD
  tukey_res <- TukeyHSD(fit)
  tukey_df <- as.data.frame(tukey_res$Group)
  tukey_df$comparison <- rownames(tukey_df)
  rownames(tukey_df) <- NULL
  
  ## Dunnett vs DMSO
  dunnett_res <- summary(glht(fit, linfct = mcp(Group = "Dunnett")))
  dunnett_df <- data.frame(
    comparison = names(dunnett_res$test$coefficients),
    diff = as.numeric(dunnett_res$test$coefficients),
    p.value = as.numeric(dunnett_res$test$pvalues),
    stringsAsFactors = FALSE
  )
  dunnett_df$treatment <- sub(paste0(" - ", control_group, "$"), "", dunnett_df$comparison)
  dunnett_df <- dunnett_df[, c("comparison", "treatment", "diff", "p.value")]
  
  ## save tables
  write.csv(scores, paste0(dataset_prefix, "_NuS_scores.csv"), row.names = FALSE)
  write.csv(as.data.frame(ssgsea_res), paste0(dataset_prefix, "_NuS_ssGSEA_matrix.csv"), row.names = TRUE)
  write.csv(data.frame(method = "ANOVA", p.value = anova_p), paste0(dataset_prefix, "_NuS_ANOVA.csv"), row.names = FALSE)
  write.csv(tukey_df, paste0(dataset_prefix, "_NuS_TukeyHSD.csv"), row.names = FALSE)
  write.csv(dunnett_df, paste0(dataset_prefix, "_NuS_Dunnett_vs_", control_group, ".csv"), row.names = FALSE)
  
  ## colors
  group_levels_plot <- levels(scores$Group)
  base_colors <- c("#C4E7C1", "#93D4BC", "#51B3D1", "#08589E")
  fill_colors <- setNames(base_colors[seq_along(group_levels_plot)], group_levels_plot)
  
  ## violin plot
  p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
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
      title = paste0("NuS Score (ssGSEA) - ", dataset_prefix),
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
        margin = margin(b = 8)
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
      method = "anova",
      label = "p.format",
      label.x.npc = "center",
      size = 5
    ) +
    stat_compare_means(
      method = "t.test",
      ref.group = control_group,
      label = "p.signif",
      size = 5
    )
  
  ggsave(paste0(dataset_prefix, "_NuS_violin.pdf"), p, width = 4, height = 4, bg = "white")
  ggsave(paste0(dataset_prefix, "_NuS_violin.png"), p, width = 4, height = 4, dpi = 600, bg = "white")
  
  save(
    expr_mat, group_info, geneSets_use, ssgsea_res, scores,
    anova_summary, tukey_df, dunnett_df,
    file = paste0(dataset_prefix, "_NuS_workspace.Rdata")
  )
  
  return(list(
    scores = scores,
    anova = anova_summary,
    tukey = tukey_df,
    dunnett = dunnett_df
  ))
}

############################################################
## Step 6. DLD1 subset
############################################################
group_info_dld1 <- data.frame(
  Sample = colnames(expr_gene)[1:12],
  Group = rep(c("DLD1-Combine", "DLD1-DMSO", "DLD1-OXA", "DLD1-PRIMET"), each = 3),
  stringsAsFactors = FALSE
)

res_dld1 <- run_nus_analysis(
  expr_input = expr_gene[, 1:12, drop = FALSE],
  group_info = group_info_dld1,
  dataset_prefix = "GSE254323_DLD1",
  control_group = "DLD1-DMSO"
)

############################################################
## Step 7. HCT116 subset
############################################################
group_info_hct116 <- data.frame(
  Sample = colnames(expr_gene)[13:24],
  Group = rep(c("HCT116-Combine", "HCT116-DMSO", "HCT116-OXA", "HCT116-PRIMET"), each = 3),
  stringsAsFactors = FALSE
)

res_hct116 <- run_nus_analysis(
  expr_input = expr_gene[, 13:24, drop = FALSE],
  group_info = group_info_hct116,
  dataset_prefix = "GSE254323_HCT116",
  control_group = "HCT116-DMSO"
)

############################################################
## Step 8. Save final workspace
############################################################
save(
  expr_gene,
  res_dld1,
  res_hct116,
  file = paste0(geo_id, "_full_NuS_analysis_workspace.Rdata")
)

cat("Analysis completed.\n")
cat("Results saved in:\n")
cat(base_dir, "\n")
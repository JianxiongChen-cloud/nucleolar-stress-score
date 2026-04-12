# GSE62963
setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE62963/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(GEOquery)
library(limma)

# ==============================
# 1. Download / load GEO ExpressionSet
# ==============================
f <- "GSE62963_eSet.Rdata"
if (!file.exists(f)) {
  gset <- getGEO(
    "GSE62963",
    destdir = ".",
    AnnotGPL = FALSE,
    getGPL = FALSE
  )
  save(gset, file = f)
}
load(f)

a <- gset[[1]]
expr_probe <- exprs(a)
expr_probe <- as.matrix(expr_probe)
mode(expr_probe) <- "numeric"

cat("===== Raw matrix dimension =====\n")
print(dim(expr_probe))

cat("===== Sample names =====\n")
print(colnames(expr_probe))

pd <- pData(a)

# ==============================
# 2. Download platform annotation
# ==============================
gpl_file <- "GPL19372_platform.Rdata"
if (!file.exists(gpl_file)) {
  gpl <- getGEO(GEO = "GPL19372", destdir = ".")
  save(gpl, file = gpl_file)
}
load(gpl_file)

gpl_tab <- Table(gpl)
cat("===== GPL column names =====\n")
print(colnames(gpl_tab))

# 假设平台文件中：
# ID = probe/feature ID
# Symbol = gene symbol
ids <- data.frame(
  probe_id = as.character(gpl_tab$ID),
  symbol = as.character(gpl_tab$Symbol),
  stringsAsFactors = FALSE
)

ids <- ids[
  !is.na(ids$symbol) &
    ids$symbol != "" &
    !grepl("//", ids$symbol, fixed = TRUE),
]

ids <- unique(ids)
ids <- ids[ids$probe_id %in% rownames(expr_probe), ]

cat("===== Annotated probes =====\n")
print(nrow(ids))

# ==============================
# 3. Keep matched probes
# ==============================
expr_probe_annot <- expr_probe[ids$probe_id, , drop = FALSE]
ids <- ids[match(rownames(expr_probe_annot), ids$probe_id), ]

# ==============================
# 4. Collapse duplicated probes by highest median expression
# ==============================
ids$median_expr <- apply(expr_probe_annot, 1, median, na.rm = TRUE)
ids <- ids[order(ids$symbol, -ids$median_expr), ]
ids <- ids[!duplicated(ids$symbol), ]

expr_gene <- expr_probe_annot[ids$probe_id, , drop = FALSE]
rownames(expr_gene) <- ids$symbol

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

write.csv(expr_gene, file = "expressionmetrix_GSE62963.csv")

save(
  expr_probe,
  expr_probe_annot,
  expr_gene,
  ids,
  pd,
  file = "GSE62963_step1_workspace.Rdata"
)

# ==============================
# 5. Define sample groups by time point
#    Matrix is already log2(relative to 0h), so 0h should be ~0
# ==============================
samples_0h <- c("GSM1537251", "GSM1537256", "GSM1537261", "GSM1537266")
samples_1h <- c("GSM1537252", "GSM1537257", "GSM1537262", "GSM1537267")
samples_3h <- c("GSM1537253", "GSM1537258", "GSM1537263", "GSM1537268")
samples_6h <- c("GSM1537254", "GSM1537259", "GSM1537264", "GSM1537269")
samples_8h <- c("GSM1537255", "GSM1537260", "GSM1537265", "GSM1537270")

cat("===== 0h range =====\n")
print(range(expr_gene[, samples_0h], na.rm = TRUE))

# ==============================
# 6. Helper function:
#    analyze one time point against baseline 0
#    using:
#    (1) mean log2FC
#    (2) SD across 4 reps
#    (3) direction consistency
#    (4) moderated t-test for mean != 0
# ==============================
analyze_timepoint <- function(expr_mat, sample_vec, prefix,
                              logfc_cutoff = 1,
                              sd_cutoff = 0.5,
                              p_cutoff = 0.05,
                              use_fdr = FALSE) {
  
  expr_sub <- expr_mat[, sample_vec, drop = FALSE]
  expr_sub <- as.matrix(expr_sub)
  mode(expr_sub) <- "numeric"
  
  cat("\n=============================\n")
  cat("Analyzing:", prefix, "\n")
  cat("=============================\n")
  cat("Samples:\n")
  print(colnames(expr_sub))
  
  # Summary of replicate values
  mean_log2FC <- rowMeans(expr_sub, na.rm = TRUE)
  median_log2FC <- apply(expr_sub, 1, median, na.rm = TRUE)
  sd_log2FC <- apply(expr_sub, 1, sd, na.rm = TRUE)
  
  # Same direction across all 4 replicates
  direction_consistent <- apply(expr_sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(FALSE)
    all(x > 0) || all(x < 0)
  })
  
  # At least 3/4 in same direction (备用，较宽松)
  direction_3of4 <- apply(expr_sub, 1, function(x) {
    x <- x[!is.na(x)]
    pos_n <- sum(x > 0)
    neg_n <- sum(x < 0)
    (pos_n >= 3) || (neg_n >= 3)
  })
  
  # Moderated one-sample test: mean(log2FC) != 0
  design <- matrix(1, ncol = 1, nrow = ncol(expr_sub))
  colnames(design) <- "mean_shift"
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  
  stat_tab <- topTable(
    fit,
    coef = "mean_shift",
    number = Inf,
    sort.by = "P",
    adjust.method = "fdr"
  )
  
  stat_tab$GeneSymbol <- rownames(stat_tab)
  
  res <- data.frame(
    GeneSymbol = rownames(expr_sub),
    mean_log2FC = mean_log2FC,
    median_log2FC = median_log2FC,
    sd_log2FC = sd_log2FC,
    direction_consistent = direction_consistent,
    direction_3of4 = direction_3of4,
    stringsAsFactors = FALSE
  )
  
  res <- merge(res, stat_tab[, c("GeneSymbol", "t", "P.Value", "adj.P.Val", "B")],
               by = "GeneSymbol", all.x = TRUE, sort = FALSE)
  
  # Stable-change criteria
  if (use_fdr) {
    res$selected_strict <- with(
      res,
      abs(mean_log2FC) > logfc_cutoff &
        sd_log2FC <= sd_cutoff &
        direction_consistent &
        adj.P.Val < p_cutoff
    )
    
    res$selected_loose <- with(
      res,
      abs(mean_log2FC) > logfc_cutoff &
        sd_log2FC <= sd_cutoff &
        direction_3of4 &
        adj.P.Val < p_cutoff
    )
  } else {
    res$selected_strict <- with(
      res,
      abs(mean_log2FC) > logfc_cutoff &
        sd_log2FC <= sd_cutoff &
        direction_consistent &
        P.Value < p_cutoff
    )
    
    res$selected_loose <- with(
      res,
      abs(mean_log2FC) > logfc_cutoff &
        sd_log2FC <= sd_cutoff &
        direction_3of4 &
        P.Value < p_cutoff
    )
  }
  
  # Direction
  res$regulation <- ifelse(
    res$mean_log2FC > logfc_cutoff, "Up",
    ifelse(res$mean_log2FC < -logfc_cutoff, "Down", "NS")
  )
  
  # Output summaries
  strict_genes <- subset(res, selected_strict)
  loose_genes <- subset(res, selected_loose)
  
  strict_up <- subset(strict_genes, regulation == "Up")
  strict_down <- subset(strict_genes, regulation == "Down")
  
  loose_up <- subset(loose_genes, regulation == "Up")
  loose_down <- subset(loose_genes, regulation == "Down")
  
  cat("===== Summary:", prefix, "=====\n")
  cat("Total genes:", nrow(res), "\n")
  cat("Strict selected genes:", nrow(strict_genes), "\n")
  cat("Strict up:", nrow(strict_up), "\n")
  cat("Strict down:", nrow(strict_down), "\n")
  cat("Loose selected genes:", nrow(loose_genes), "\n")
  cat("Loose up:", nrow(loose_up), "\n")
  cat("Loose down:", nrow(loose_down), "\n")
  
  # Save outputs
  write.csv(res, paste0(prefix, "_all_results.csv"), row.names = FALSE)
  write.csv(strict_genes, paste0(prefix, "_strict_genes.csv"), row.names = FALSE)
  write.csv(strict_up, paste0(prefix, "_strict_up.csv"), row.names = FALSE)
  write.csv(strict_down, paste0(prefix, "_strict_down.csv"), row.names = FALSE)
  
  write.csv(loose_genes, paste0(prefix, "_loose_genes.csv"), row.names = FALSE)
  write.csv(loose_up, paste0(prefix, "_loose_up.csv"), row.names = FALSE)
  write.csv(loose_down, paste0(prefix, "_loose_down.csv"), row.names = FALSE)
  
  return(res)
}

# ==============================
# 7. Run analysis for each time point
#    use_fdr = FALSE because many time-course datasets are not strong enough
#    to pass adj.P.Val < 0.05 after filtering
# ==============================
res_1h <- analyze_timepoint(
  expr_mat = expr_gene,
  sample_vec = samples_1h,
  prefix = "GSE62963_1h",
  logfc_cutoff = 1,
  sd_cutoff = 0.5,
  p_cutoff = 0.05,
  use_fdr = FALSE
)

res_3h <- analyze_timepoint(
  expr_mat = expr_gene,
  sample_vec = samples_3h,
  prefix = "GSE62963_3h",
  logfc_cutoff = 1,
  sd_cutoff = 0.5,
  p_cutoff = 0.05,
  use_fdr = FALSE
)

res_6h <- analyze_timepoint(
  expr_mat = expr_gene,
  sample_vec = samples_6h,
  prefix = "GSE62963_6h",
  logfc_cutoff = 1,
  sd_cutoff = 0.5,
  p_cutoff = 0.05,
  use_fdr = FALSE
)

res_8h <- analyze_timepoint(
  expr_mat = expr_gene,
  sample_vec = samples_8h,
  prefix = "GSE62963_8h",
  logfc_cutoff = 1,
  sd_cutoff = 0.5,
  p_cutoff = 0.05,
  use_fdr = FALSE
)

# ==============================
# 8. Save full workspace
# ==============================
save(
  expr_probe,
  expr_probe_annot,
  expr_gene,
  ids,
  pd,
  res_1h,
  res_3h,
  res_6h,
  res_8h,
  file = "GSE62963_full_analysis_workspace.Rdata"
)
# ==============================
# 9. Summary table for each time point
# ==============================
get_summary_counts <- function(res, prefix) {
  strict_genes <- subset(res, selected_strict)
  strict_up <- subset(strict_genes, regulation == "Up")
  strict_down <- subset(strict_genes, regulation == "Down")
  
  loose_genes <- subset(res, selected_loose)
  loose_up <- subset(loose_genes, regulation == "Up")
  loose_down <- subset(loose_genes, regulation == "Down")
  
  data.frame(
    TimePoint = prefix,
    Strict_Total = nrow(strict_genes),
    Strict_Up = nrow(strict_up),
    Strict_Down = nrow(strict_down),
    Loose_Total = nrow(loose_genes),
    Loose_Up = nrow(loose_up),
    Loose_Down = nrow(loose_down),
    stringsAsFactors = FALSE
  )
}

summary_tab <- rbind(
  get_summary_counts(res_1h, "1h"),
  get_summary_counts(res_3h, "3h"),
  get_summary_counts(res_6h, "6h"),
  get_summary_counts(res_8h, "8h")
)

cat("===== Time-point DEG summary =====\n")
print(summary_tab)

write.csv(summary_tab, "GSE62963_timepoint_summary_counts.csv", row.names = FALSE)

# ==============================
# 10. Common DEGs across all time points
#    Use strict result first
# ==============================
strict_1h <- subset(res_1h, selected_strict)
strict_3h <- subset(res_3h, selected_strict)
strict_6h <- subset(res_6h, selected_strict)
strict_8h <- subset(res_8h, selected_strict)

common_strict_genes <- Reduce(
  intersect,
  list(strict_1h$GeneSymbol, strict_3h$GeneSymbol, strict_6h$GeneSymbol, strict_8h$GeneSymbol)
)

cat("===== Common strict DEGs across 1h/3h/6h/8h =====\n")
cat("Common strict genes:", length(common_strict_genes), "\n")

# 合并共同基因在各时间点的信息
common_strict_res <- data.frame(GeneSymbol = common_strict_genes, stringsAsFactors = FALSE)

common_strict_res <- merge(
  common_strict_res,
  strict_1h[, c("GeneSymbol", "mean_log2FC", "sd_log2FC", "P.Value", "adj.P.Val", "regulation")],
  by = "GeneSymbol",
  all.x = TRUE,
  sort = FALSE
)
colnames(common_strict_res)[2:6] <- c("mean_log2FC_1h", "sd_log2FC_1h", "P_1h", "FDR_1h", "regulation_1h")

common_strict_res <- merge(
  common_strict_res,
  strict_3h[, c("GeneSymbol", "mean_log2FC", "sd_log2FC", "P.Value", "adj.P.Val", "regulation")],
  by = "GeneSymbol",
  all.x = TRUE,
  sort = FALSE
)
colnames(common_strict_res)[7:11] <- c("mean_log2FC_3h", "sd_log2FC_3h", "P_3h", "FDR_3h", "regulation_3h")

common_strict_res <- merge(
  common_strict_res,
  strict_6h[, c("GeneSymbol", "mean_log2FC", "sd_log2FC", "P.Value", "adj.P.Val", "regulation")],
  by = "GeneSymbol",
  all.x = TRUE,
  sort = FALSE
)
colnames(common_strict_res)[12:16] <- c("mean_log2FC_6h", "sd_log2FC_6h", "P_6h", "FDR_6h", "regulation_6h")

common_strict_res <- merge(
  common_strict_res,
  strict_8h[, c("GeneSymbol", "mean_log2FC", "sd_log2FC", "P.Value", "adj.P.Val", "regulation")],
  by = "GeneSymbol",
  all.x = TRUE,
  sort = FALSE
)
colnames(common_strict_res)[17:21] <- c("mean_log2FC_8h", "sd_log2FC_8h", "P_8h", "FDR_8h", "regulation_8h")

# 判断共同基因方向是否一致
common_strict_res$direction_pattern <- apply(
  common_strict_res[, c("regulation_1h", "regulation_3h", "regulation_6h", "regulation_8h")],
  1,
  function(x) {
    if (all(x == "Up")) {
      return("Always_Up")
    } else if (all(x == "Down")) {
      return("Always_Down")
    } else {
      return("Mixed")
    }
  }
)

cat("Always Up:", sum(common_strict_res$direction_pattern == "Always_Up"), "\n")
cat("Always Down:", sum(common_strict_res$direction_pattern == "Always_Down"), "\n")
cat("Mixed direction:", sum(common_strict_res$direction_pattern == "Mixed"), "\n")

write.csv(common_strict_res, "GSE62963_common_strict_DEGs_across_all_timepoints.csv", row.names = FALSE)

# ==============================
# 11. Loose common DEGs across all time points
# ==============================
loose_1h <- subset(res_1h, selected_loose)
loose_3h <- subset(res_3h, selected_loose)
loose_6h <- subset(res_6h, selected_loose)
loose_8h <- subset(res_8h, selected_loose)

common_loose_genes <- Reduce(
  intersect,
  list(loose_1h$GeneSymbol, loose_3h$GeneSymbol, loose_6h$GeneSymbol, loose_8h$GeneSymbol)
)

cat("===== Common loose DEGs across 1h/3h/6h/8h =====\n")
cat("Common loose genes:", length(common_loose_genes), "\n")

write.csv(
  data.frame(GeneSymbol = common_loose_genes),
  "GSE62963_common_loose_DEGs_gene_list.csv",
  row.names = FALSE
)


# 提取共同 strict DEGs 的上下调基因
common_up_genes <- unique(na.omit(
  common_strict_res$GeneSymbol[common_strict_res$direction_pattern == "Always_Up"]
))

common_down_genes <- unique(na.omit(
  common_strict_res$GeneSymbol[common_strict_res$direction_pattern == "Always_Down"]
))

# 生成 long 格式
gene_list_long <- data.frame(
  dataset = c(rep("GSE62963", length(common_up_genes)),
              rep("GSE62963", length(common_down_genes))),
  group = c(rep("UP", length(common_up_genes)),
            rep("DOWN", length(common_down_genes))),
  gene_symbol = c(common_up_genes, common_down_genes)
)

write.csv(
  gene_list_long,
  "GSE62963_common_strict_DEGs_diffgenesymbol_long.csv",
  row.names = FALSE
)

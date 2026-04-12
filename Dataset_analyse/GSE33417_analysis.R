setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE33417/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(GEOquery)
library(hugene10sttranscriptcluster.db)
library(limma)
library(AnnotationDbi)

# ==============================
# 1. Download / load GEO ExpressionSet
# ==============================
f <- "GSE33417_eSet.Rdata"
if (!file.exists(f)) {
  gset <- getGEO("GSE33417", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
  save(gset, file = f)
}
load(f)

eset <- gset[[1]]
expr_probe <- exprs(eset)
expr_probe <- as.matrix(expr_probe)
mode(expr_probe) <- "numeric"

cat("===== Raw matrix dimension =====\n")
print(dim(expr_probe))
cat("===== First few sample names =====\n")
print(colnames(expr_probe)[1:min(10, ncol(expr_probe))])

# ==============================
# 2. Probe annotation
# ==============================
ids <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = rownames(expr_probe),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)

colnames(ids) <- c("probe_id", "symbol")
ids <- ids[!is.na(ids$symbol) & ids$symbol != "", ]
ids <- ids[ids$probe_id %in% rownames(expr_probe), ]
ids <- unique(ids)

cat("===== Annotated probes =====\n")
print(nrow(ids))

# ==============================
# 3. Keep matched probes
# ==============================
expr_probe_annot <- expr_probe[ids$probe_id, , drop = FALSE]
ids <- ids[match(rownames(expr_probe_annot), ids$probe_id), ]

cat("===== Matched annotated matrix dimension =====\n")
print(dim(expr_probe_annot))

# ==============================
# 4. Collapse duplicated probes
#    Strategy: keep the probe with highest median expression
# ==============================
probe_median <- apply(expr_probe_annot, 1, median, na.rm = TRUE)
ids$median_expr <- probe_median

ids <- ids[order(ids$symbol, -ids$median_expr), ]
ids <- ids[!duplicated(ids$symbol), ]

expr_gene <- expr_probe_annot[ids$probe_id, , drop = FALSE]
rownames(expr_gene) <- ids$symbol

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

write.csv(expr_gene, file = "expressionmetrix_GSE33417.csv")
save(expr_probe, expr_probe_annot, expr_gene,
     file = "GSE33417_step1_workspace.Rdata")

cat("===== All sample names =====\n")
print(colnames(expr_gene))

# ==========================================================
# Helper function: scale check + DEG
# ==========================================================
run_deg_with_check <- function(expr_sub, dataset_prefix, group_vector, contrast_name) {
  
  cat("\n=============================\n")
  cat("Running:", dataset_prefix, "\n")
  cat("=============================\n")
  
  expr_sub <- as.matrix(expr_sub)
  mode(expr_sub) <- "numeric"
  
  cat("===== Selected samples =====\n")
  print(colnames(expr_sub))
  
  # ------------------------------
  # 1. Scale check
  # ------------------------------
  cat("===== Overall range =====\n")
  print(range(expr_sub, na.rm = TRUE))
  
  cat("===== Overall quantiles =====\n")
  qx <- quantile(
    as.numeric(expr_sub),
    probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
    na.rm = TRUE
  )
  print(qx)
  
  sample_summary <- data.frame(
    Sample = colnames(expr_sub),
    Min = apply(expr_sub, 2, min, na.rm = TRUE),
    Q1 = apply(expr_sub, 2, quantile, probs = 0.25, na.rm = TRUE),
    Median = apply(expr_sub, 2, median, na.rm = TRUE),
    Mean = apply(expr_sub, 2, mean, na.rm = TRUE),
    Q3 = apply(expr_sub, 2, quantile, probs = 0.75, na.rm = TRUE),
    Max = apply(expr_sub, 2, max, na.rm = TRUE)
  )
  
  cat("===== Per-sample summary =====\n")
  print(sample_summary)
  write.csv(sample_summary,
            paste0(dataset_prefix, "_sample_summary.csv"),
            row.names = FALSE)
  
  pdf(paste0(dataset_prefix, "_boxplot_before_logcheck.pdf"), width = 10, height = 6)
  boxplot(expr_sub,
          las = 2,
          outline = FALSE,
          main = paste0(dataset_prefix, " expression before log-check"))
  dev.off()
  
  # GEO2R / limma style judgement
  qx6 <- quantile(
    as.numeric(expr_sub),
    c(0, 0.25, 0.5, 0.75, 0.99, 1.0),
    na.rm = TRUE
  )
  
  need_log2 <- (
    qx6[5] > 100 ||
      (qx6[6] - qx6[1] > 50 && qx6[2] > 0) ||
      (qx6[2] > 0 && qx6[2] < 1 && qx6[4] > 1 && qx6[4] < 2)
  )
  
  cat("===== Scale judgement =====\n")
  if (need_log2) {
    cat("This matrix is likely NOT on log2 scale. log2 transformation will be applied.\n")
    expr_sub[expr_sub <= 0] <- NA
    exp_matrix <- log2(expr_sub)
    log_status <- "log2_transformed_in_this_script"
  } else {
    cat("This matrix is likely already on log scale. No additional log2 transformation will be applied.\n")
    exp_matrix <- expr_sub
    log_status <- "already_log_scale_no_additional_log2"
  }
  
  overall_max <- max(expr_sub, na.rm = TRUE)
  overall_median <- median(expr_sub, na.rm = TRUE)
  
  write.table(
    data.frame(
      Dataset = dataset_prefix,
      Platform = "Affymetrix Human Gene 1.0 ST Array",
      Processing = "APT median normalization, no adjustment, normalized to 5th percentile",
      InputType = "processed expression matrix",
      ReNormalization = "No",
      LogDecision = log_status,
      MaxValue = overall_max,
      MedianValue = overall_median,
      DEGmethod = "limma",
      ProbeCollapse = "highest median expression"
    ),
    file = paste0(dataset_prefix, "_logcheck_decision.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  pdf(paste0(dataset_prefix, "_boxplot_used_for_DEG.pdf"), width = 10, height = 6)
  boxplot(exp_matrix,
          las = 2,
          outline = FALSE,
          main = paste0(dataset_prefix, " expression used for DEG"))
  dev.off()
  
  pdf(paste0(dataset_prefix, "_density_used_for_DEG.pdf"), width = 10, height = 6)
  plotDensities(exp_matrix,
                main = paste0(dataset_prefix, " density plot used for DEG"))
  dev.off()
  
  # ------------------------------
  # 2. DEG analysis
  # ------------------------------
  group <- factor(group_vector, levels = c("CONTROL", "ActD"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  cat("===== Design matrix =====\n")
  print(design)
  
  fit <- lmFit(exp_matrix, design)
  
  contrast_matrix <- makeContrasts(
    contrasts = contrast_name,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  all_diff <- topTable(
    fit2,
    coef = 1,
    adjust.method = "fdr",
    number = Inf,
    sort.by = "P"
  )
  
  all_diff$GeneSymbol <- rownames(all_diff)
  all_diff <- all_diff[, c("GeneSymbol", setdiff(colnames(all_diff), "GeneSymbol"))]
  
  diffSig <- subset(all_diff, abs(logFC) > 1 & adj.P.Val < 0.05)
  up_diff <- subset(all_diff, logFC > 1 & adj.P.Val < 0.05)
  down_diff <- subset(all_diff, logFC < -1 & adj.P.Val < 0.05)
  
  cat("===== DEG summary =====\n")
  cat("Total genes tested:", nrow(all_diff), "\n")
  cat("Total significant DEGs:", nrow(diffSig), "\n")
  cat("Upregulated genes:", nrow(up_diff), "\n")
  cat("Downregulated genes:", nrow(down_diff), "\n")
  
  write.csv(all_diff, paste0(dataset_prefix, "_all_DEG.csv"), row.names = FALSE)
  write.csv(diffSig, paste0(dataset_prefix, "_sig_DEG.csv"), row.names = FALSE)
  write.csv(up_diff, paste0(dataset_prefix, "_up_DEG.csv"), row.names = FALSE)
  write.csv(down_diff, paste0(dataset_prefix, "_down_DEG.csv"), row.names = FALSE)
  
# ------------------------------
  # 3. Output long-format gene symbol list
  # ------------------------------
  up_genes <- unique(na.omit(up_diff$GeneSymbol))
  down_genes <- unique(na.omit(down_diff$GeneSymbol))
  
  gene_list_long <- data.frame(
    dataset = c(rep(dataset_prefix, length(up_genes)),
                rep(dataset_prefix, length(down_genes))),
    group = c(rep("UP", length(up_genes)),
              rep("DOWN", length(down_genes))),
    gene_symbol = c(up_genes, down_genes)
  )
  
  write.csv(
    gene_list_long,
    paste0(dataset_prefix, "_diffgenesymbol_long.csv"),
    row.names = FALSE
  )

  save(expr_sub, exp_matrix, sample_summary, all_diff, diffSig, up_diff, down_diff, gene_list_long, file = paste0(dataset_prefix, "_analysis_workspace.Rdata"))
}

# ==========================================================
# Analysis 1: iPS
# original columns: 1:3 vs 16:18
# ==========================================================
expr_ips <- expr_gene[, c(1:3, 16:18), drop = FALSE]
run_deg_with_check(
  expr_sub = expr_ips,
  dataset_prefix = "GSE33417_ips",
  group_vector = c(rep("CONTROL", 3), rep("ActD", 3)),
  contrast_name = "ActD - CONTROL"
)

# ==========================================================
# Analysis 2: HFF
# original columns: 19:21 vs 34:36
# ==========================================================
expr_hff <- expr_gene[, c(19:21, 34:36), drop = FALSE]
run_deg_with_check(
  expr_sub = expr_hff,
  dataset_prefix = "GSE33417_hff",
  group_vector = c(rep("CONTROL", 3), rep("ActD", 3)),
  contrast_name = "ActD - CONTROL"
)

save(expr_probe, expr_probe_annot, expr_gene, expr_ips, expr_hff,
     file = "GSE33417_full_analysis_workspace.Rdata")
# GSE298220
rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE298220/")

library(biomaRt)
library(dplyr)
library(tibble)
library(limma)

# ==============================
# 1. Read processed matrix
# ==============================
tpm <- read.table(
  "GSE298220_Processed_Data_DMSO_CX_BMH.txt.gz",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

tpm <- as.matrix(tpm)
mode(tpm) <- "numeric"

cat("===== Raw matrix dimension =====\n")
print(dim(tpm))
cat("===== Sample names =====\n")
print(colnames(tpm))
cat("===== Raw range =====\n")
print(range(tpm, na.rm = TRUE))
cat("===== NA count =====\n")
print(sum(is.na(tpm)))

# ==============================
# 2. Map Ensembl IDs to gene symbols
# ==============================
genes_clean <- sub("\\..*$", "", rownames(tpm))

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

symbol_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(genes_clean),
  mart = ensembl
)

colnames(symbol_map) <- c("GeneID", "Symbol")
symbol_map <- symbol_map[!is.na(symbol_map$Symbol) & symbol_map$Symbol != "", ]
symbol_map <- symbol_map[!duplicated(symbol_map$GeneID), ]

cat("===== Annotated gene IDs =====\n")
print(nrow(symbol_map))

# ==============================
# 3. Merge annotation with expression matrix
# ==============================
expr_df <- as.data.frame(tpm, check.names = FALSE) %>%
  rownames_to_column(var = "GeneID_raw") %>%
  mutate(GeneID = sub("\\..*$", "", GeneID_raw))

expr_symbol <- expr_df %>%
  left_join(symbol_map, by = "GeneID") %>%
  filter(!is.na(Symbol) & Symbol != "")

cat("===== Rows after symbol mapping =====\n")
print(nrow(expr_symbol))

# sample columns
sample_cols <- setdiff(colnames(expr_symbol), c("GeneID_raw", "GeneID", "Symbol"))

# ensure numeric
expr_symbol[, sample_cols] <- lapply(expr_symbol[, sample_cols, drop = FALSE], as.numeric)

# ==============================
# 4. Collapse duplicated symbols
#    keep row with highest median expression
# ==============================
expr_symbol$median_expr <- apply(
  expr_symbol[, sample_cols, drop = FALSE],
  1,
  median,
  na.rm = TRUE
)

expr_symbol <- expr_symbol[order(expr_symbol$Symbol, -expr_symbol$median_expr), ]
expr_symbol <- expr_symbol[!duplicated(expr_symbol$Symbol), ]

final_expr <- as.matrix(expr_symbol[, sample_cols, drop = FALSE])
mode(final_expr) <- "numeric"
rownames(final_expr) <- expr_symbol$Symbol

cat("===== final_expr dimension =====\n")
print(dim(final_expr))
cat("===== final_expr range =====\n")
print(range(final_expr, na.rm = TRUE))
cat("===== final_expr NA count =====\n")
print(sum(is.na(final_expr)))

write.csv(final_expr, file = "expressionmetrix_GSE298220.csv")
save(final_expr, file = "GSE298220_step1_output.Rdata")

# ==============================
# 5. Helper function
# ==============================
run_deg_processed <- function(dat, dataset_prefix, group_vector,
                              contrast_name, logfc_cutoff = 1,
                              expr_cutoff = 1, min_samples = 3) {
  
  dat <- as.matrix(dat)
  mode(dat) <- "numeric"
  
  cat("\n=============================\n")
  cat("Running:", dataset_prefix, "\n")
  cat("=============================\n")
  cat("Selected samples:\n")
  print(colnames(dat))
  
  cat("===== Input range =====\n")
  print(range(dat, na.rm = TRUE))
  
  cat("===== Input quantiles =====\n")
  print(quantile(
    dat,
    probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
    na.rm = TRUE
  ))
  
  sample_summary_raw <- data.frame(
    Sample = colnames(dat),
    Min = apply(dat, 2, min, na.rm = TRUE),
    Q1 = apply(dat, 2, quantile, probs = 0.25, na.rm = TRUE),
    Median = apply(dat, 2, median, na.rm = TRUE),
    Mean = apply(dat, 2, mean, na.rm = TRUE),
    Q3 = apply(dat, 2, quantile, probs = 0.75, na.rm = TRUE),
    Max = apply(dat, 2, max, na.rm = TRUE)
  )
  
  write.csv(sample_summary_raw,
            paste0(dataset_prefix, "_sample_summary_raw.csv"),
            row.names = FALSE)
  
  # low-expression filtering
  keep <- rowSums(dat > expr_cutoff, na.rm = TRUE) >= min_samples
  dat_filt <- dat[keep, , drop = FALSE]
  
  cat("===== Filtering summary =====\n")
  cat("Genes before filtering:", nrow(dat), "\n")
  cat("Genes after filtering:", nrow(dat_filt), "\n")
  
  # log transform
  exp_matrix <- log2(dat_filt + 1)
  
  cat("===== log2(x+1) range =====\n")
  print(range(exp_matrix, na.rm = TRUE))
  
  sample_summary_log <- data.frame(
    Sample = colnames(exp_matrix),
    Min = apply(exp_matrix, 2, min, na.rm = TRUE),
    Q1 = apply(exp_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
    Median = apply(exp_matrix, 2, median, na.rm = TRUE),
    Mean = apply(exp_matrix, 2, mean, na.rm = TRUE),
    Q3 = apply(exp_matrix, 2, quantile, probs = 0.75, na.rm = TRUE),
    Max = apply(exp_matrix, 2, max, na.rm = TRUE)
  )
  
  write.csv(sample_summary_log,
            paste0(dataset_prefix, "_sample_summary_log.csv"),
            row.names = FALSE)
  
  pdf(paste0(dataset_prefix, "_boxplot_log_used_for_DEG.pdf"), width = 10, height = 6)
  boxplot(exp_matrix,
          las = 2,
          outline = FALSE,
          main = paste0(dataset_prefix, " log2(x+1) used for DEG"))
  dev.off()
  
  pdf(paste0(dataset_prefix, "_density_log_used_for_DEG.pdf"), width = 10, height = 6)
  plotDensities(exp_matrix, main = paste0(dataset_prefix, " density plot"))
  dev.off()
  
  # remove zero-variance genes
  gene_var <- apply(exp_matrix, 1, var, na.rm = TRUE)
  exp_matrix <- exp_matrix[gene_var > 0, , drop = FALSE]
  
  group <- factor(group_vector)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  cat("===== Design matrix =====\n")
  print(design)
  
  fit <- lmFit(exp_matrix, design)
  contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = design)
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
  
  # FDR
  diffSig_fdr <- subset(all_diff, abs(logFC) > logfc_cutoff & adj.P.Val < 0.05)
  up_fdr <- subset(all_diff, logFC > logfc_cutoff & adj.P.Val < 0.05)
  down_fdr <- subset(all_diff, logFC < -logfc_cutoff & adj.P.Val < 0.05)
  
  # nominal
  diffSig_p <- subset(all_diff, abs(logFC) > logfc_cutoff & P.Value < 0.05)
  up_p <- subset(all_diff, logFC > logfc_cutoff & P.Value < 0.05)
  down_p <- subset(all_diff, logFC < -logfc_cutoff & P.Value < 0.05)
  
  cat("===== DEG summary =====\n")
  cat("Total genes tested:", nrow(all_diff), "\n")
  cat("FDR significant DEGs:", nrow(diffSig_fdr), "\n")
  cat("FDR up:", nrow(up_fdr), "\n")
  cat("FDR down:", nrow(down_fdr), "\n")
  cat("Nominal significant DEGs:", nrow(diffSig_p), "\n")
  cat("Nominal up:", nrow(up_p), "\n")
  cat("Nominal down:", nrow(down_p), "\n")
  
  write.csv(all_diff, paste0(dataset_prefix, "_all_DEG.csv"), row.names = FALSE)
  write.csv(diffSig_fdr, paste0(dataset_prefix, "_sig_DEG_FDR.csv"), row.names = FALSE)
  write.csv(up_fdr, paste0(dataset_prefix, "_up_DEG_FDR.csv"), row.names = FALSE)
  write.csv(down_fdr, paste0(dataset_prefix, "_down_DEG_FDR.csv"), row.names = FALSE)
  write.csv(diffSig_p, paste0(dataset_prefix, "_sig_DEG_Pvalue.csv"), row.names = FALSE)
  write.csv(up_p, paste0(dataset_prefix, "_up_DEG_Pvalue.csv"), row.names = FALSE)
  write.csv(down_p, paste0(dataset_prefix, "_down_DEG_Pvalue.csv"), row.names = FALSE)
  
  # ==============================
  # output FDR long list
  # ==============================
  up_genes <- unique(na.omit(up_fdr$GeneSymbol))
  down_genes <- unique(na.omit(down_fdr$GeneSymbol))
  
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
  
  save(dat, dat_filt, exp_matrix,
       all_diff, diffSig_fdr, up_fdr, down_fdr,
       diffSig_p, up_p, down_p, gene_list_long,
       file = paste0(dataset_prefix, "_analysis_workspace.Rdata"))
  
  return(list(
    all_diff = all_diff,
    diffSig_fdr = diffSig_fdr,
    up_fdr = up_fdr,
    down_fdr = down_fdr,
    diffSig_p = diffSig_p,
    up_p = up_p,
    down_p = down_p,
    gene_list_long = gene_list_long
  ))
}

# ==============================
# 6. BMH21 vs CONTROL
# sample1-3 = CONTROL
# sample7-9 = BMH21
# ==============================
res_bmh <- run_deg_processed(
  dat = final_expr[, c(7:9, 1:3), drop = FALSE],
  dataset_prefix = "GSE298220_BMH21",
  group_vector = c(rep("BMH21", 3), rep("CONTROL", 3)),
  contrast_name = "BMH21 - CONTROL",
  logfc_cutoff = 1
)

# ==============================
# 7. CX5461 vs CONTROL
# sample4-6 = CX5461
# sample1-3 = CONTROL
# ==============================
res_cx <- run_deg_processed(
  dat = final_expr[, c(4:6, 1:3), drop = FALSE],
  dataset_prefix = "GSE298220_CX5461",
  group_vector = c(rep("CX5461", 3), rep("CONTROL", 3)),
  contrast_name = "CX5461 - CONTROL",
  logfc_cutoff = 0.5
)

save(final_expr, res_bmh, res_cx, file = "GSE298220_full_analysis_workspace.Rdata")
setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE198178/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(tibble)
library(limma)

# ==============================
# 1. Read processed TPM matrix
# ==============================
tpm <- read.table(
  "GSE198178_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

tpm <- as.matrix(tpm)
mode(tpm) <- "numeric"

cat("===== Raw TPM matrix dimension =====\n")
print(dim(tpm))
cat("===== First few sample names =====\n")
print(colnames(tpm)[1:min(10, ncol(tpm))])

# ==============================
# 2. Read gene annotation
# ==============================
probe <- fread(
  "Human.GRCh38.p13.annot.tsv.gz",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = ""
)

probe_map <- probe[, c("GeneID", "Symbol")]
probe_map <- as.data.frame(probe_map)
probe_map <- probe_map[!is.na(probe_map$GeneID) & !is.na(probe_map$Symbol) & probe_map$Symbol != "", ]
probe_map$GeneID <- as.character(probe_map$GeneID)
probe_map$Symbol <- as.character(probe_map$Symbol)
probe_map <- probe_map[!duplicated(probe_map$GeneID), ]

# ==============================
# 3. Merge symbol
# ==============================
tpm_df <- as.data.frame(tpm, check.names = FALSE) %>%
  rownames_to_column(var = "GeneID")

tpm_symbol <- tpm_df %>%
  mutate(GeneID = as.character(GeneID)) %>%
  left_join(probe_map, by = "GeneID") %>%
  filter(!is.na(Symbol) & Symbol != "")

# ==============================
# 4. Collapse duplicated symbols by mean
# ==============================
tpm_final <- tpm_symbol %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "Symbol")

tpm_final <- as.matrix(tpm_final)
mode(tpm_final) <- "numeric"

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

write.csv(tpm_final, file = "expressionmetrix_GSE198178.csv")
save(tpm_final, file = "GSE198178_step1_workspace.Rdata")

cat("===== all sample names =====\n")
print(colnames(tpm_final))

############################################################
# Helper function: TPM -> filter -> log2(TPM+1) -> limma
############################################################
run_deg_tpm <- function(dat, dataset_prefix, group_vector, contrast_name,
                        tpm_cutoff = 1, min_samples = 3) {
  
  cat("\n=============================\n")
  cat("Running:", dataset_prefix, "\n")
  cat("=============================\n")
  
  cat("Selected samples:\n")
  print(colnames(dat))
  
  cat("===== overall range (raw TPM) =====\n")
  print(range(dat, na.rm = TRUE))
  
  cat("===== overall quantiles (raw TPM) =====\n")
  print(quantile(dat,
                 probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
                 na.rm = TRUE))
  
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
            paste0(dataset_prefix, "_sample_summary_rawTPM.csv"),
            row.names = FALSE)
  
  pdf(paste0(dataset_prefix, "_boxplot_rawTPM.pdf"), width = 10, height = 6)
  boxplot(dat,
          las = 2,
          outline = FALSE,
          main = paste0(dataset_prefix, " raw TPM"))
  dev.off()
  
  # ------------------------------
  # Filter low-expression genes
  # TPM > 1 in at least 3 samples
  # ------------------------------
  keep <- rowSums(dat > tpm_cutoff, na.rm = TRUE) >= min_samples
  dat_filt <- dat[keep, , drop = FALSE]
  
  cat("===== filtering =====\n")
  cat("Genes before filtering:", nrow(dat), "\n")
  cat("Genes after filtering:", nrow(dat_filt), "\n")
  
  # ------------------------------
  # log2 transform
  # ------------------------------
  exp_matrix <- log2(dat_filt + 1)
  
  cat("===== overall range after log2(TPM+1) =====\n")
  print(range(exp_matrix, na.rm = TRUE))
  
  cat("===== overall quantiles after log2(TPM+1) =====\n")
  print(quantile(exp_matrix,
                 probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
                 na.rm = TRUE))
  
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
            paste0(dataset_prefix, "_sample_summary_logTPM.csv"),
            row.names = FALSE)
  
  write.table(
    data.frame(
      Dataset = dataset_prefix,
      InputType = "processed TPM matrix",
      ReNormalization = "No",
      Filtering = paste0("TPM > ", tpm_cutoff, " in at least ", min_samples, " samples"),
      Transformation = "log2(TPM + 1)",
      DEGmethod = "limma",
      Note = "Exploratory DEG from processed TPM"
    ),
    file = paste0(dataset_prefix, "_processing_decision.txt"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  pdf(paste0(dataset_prefix, "_boxplot_logTPM_used_for_DEG.pdf"), width = 10, height = 6)
  boxplot(exp_matrix,
          las = 2,
          outline = FALSE,
          main = paste0(dataset_prefix, " log2(TPM+1) used for DEG"))
  dev.off()
  
  pdf(paste0(dataset_prefix, "_density_logTPM_used_for_DEG.pdf"), width = 10, height = 6)
  plotDensities(exp_matrix,
                main = paste0(dataset_prefix, " density plot"))
  dev.off()
  
  # remove zero-variance genes
  gene_var <- apply(exp_matrix, 1, var, na.rm = TRUE)
  exp_matrix <- exp_matrix[gene_var > 0, , drop = FALSE]
  
  # ------------------------------
  # DEG by limma
  # ------------------------------
  group <- factor(group_vector, levels = c("Control", "ACTD"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  cat("===== design matrix =====\n")
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
  
  # strict
  diffSig_fdr <- subset(all_diff, abs(logFC) > 1 & adj.P.Val < 0.05)
  up_fdr <- subset(all_diff, logFC > 1 & adj.P.Val < 0.05)
  down_fdr <- subset(all_diff, logFC < -1 & adj.P.Val < 0.05)
  
  # exploratory
  diffSig_p <- subset(all_diff, abs(logFC) > 1 & P.Value < 0.05)
  up_p <- subset(all_diff, logFC > 1 & P.Value < 0.05)
  down_p <- subset(all_diff, logFC < -1 & P.Value < 0.05)
  
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
  
  save(dat, dat_filt, exp_matrix, all_diff,
       diffSig_fdr, up_fdr, down_fdr,
       diffSig_p, up_p, down_p,
       file = paste0(dataset_prefix, "_analysis_workspace.Rdata"))
  
  return(list(
    all_diff = all_diff,
    diffSig_fdr = diffSig_fdr,
    up_fdr = up_fdr,
    down_fdr = down_fdr,
    diffSig_p = diffSig_p,
    up_p = up_p,
    down_p = down_p
  ))
}

############################################################
# 4 h
# original columns: 1,5,9 = Control ; 4,8,12 = ACTD
############################################################
dat_4h <- tpm_final[, c(1, 5, 9, 4, 8, 12), drop = FALSE]
res_4h <- run_deg_tpm(
  dat = dat_4h,
  dataset_prefix = "GSE198178_4h",
  group_vector = c(rep("Control", 3), rep("ACTD", 3)),
  contrast_name = "ACTD - Control"
)

############################################################
# 2 h
# original columns: 1,5,9 = Control ; 3,7,11 = ACTD
############################################################
dat_2h <- tpm_final[, c(1, 5, 9, 3, 7, 11), drop = FALSE]
res_2h <- run_deg_tpm(
  dat = dat_2h,
  dataset_prefix = "GSE198178_2h",
  group_vector = c(rep("Control", 3), rep("ACTD", 3)),
  contrast_name = "ACTD - Control"
)

############################################################
# 1 h
# original columns: 1,5,9 = Control ; 2,6,10 = ACTD
############################################################
dat_1h <- tpm_final[, c(1, 5, 9, 2, 6, 10), drop = FALSE]
res_1h <- run_deg_tpm(
  dat = dat_1h,
  dataset_prefix = "GSE198178_1h",
  group_vector = c(rep("Control", 3), rep("ACTD", 3)),
  contrast_name = "ACTD - Control"
)

############################################################
# Common FDR DEGs across 4h / 2h / 1h
############################################################
up_genes <- Reduce(intersect, list(
  unique(na.omit(res_4h$up_fdr$GeneSymbol)),
  unique(na.omit(res_2h$up_fdr$GeneSymbol)),
  unique(na.omit(res_1h$up_fdr$GeneSymbol))
))

down_genes <- Reduce(intersect, list(
  unique(na.omit(res_4h$down_fdr$GeneSymbol)),
  unique(na.omit(res_2h$down_fdr$GeneSymbol)),
  unique(na.omit(res_1h$down_fdr$GeneSymbol))
))

cat("===== Common FDR DEGs across 4h / 2h / 1h =====\n")
cat("Common UP genes:", length(up_genes), "\n")
cat("Common DOWN genes:", length(down_genes), "\n")

gene_list_long <- data.frame(
  dataset = c(rep("GSE198178", length(up_genes)), rep("GSE198178", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)

write.csv(
  gene_list_long,
  "GSE198178_common_4h_2h_1h_diffgenesymbol_long.csv",
  row.names = FALSE
)

save(tpm_final, dat_1h, dat_2h, dat_4h,
     res_1h, res_2h, res_4h, gene_list_long,
     file = "GSE198178_full_analysis_workspace.Rdata")
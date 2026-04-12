setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE118565/")
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
  "GSE118565_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz",
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
# 3. Map GeneID to gene symbol
# ==============================
tpm_df <- as.data.frame(tpm, check.names = FALSE) %>%
  rownames_to_column(var = "GeneID")

tpm_symbol <- tpm_df %>%
  mutate(GeneID = as.character(GeneID)) %>%
  left_join(probe_map, by = "GeneID") %>%
  filter(!is.na(Symbol) & Symbol != "")

# ==============================
# 4. Collapse duplicated symbols
#    For RNA-seq processed TPM, averaging duplicated symbols is acceptable
# ==============================
tpm_final <- tpm_symbol %>%
  group_by(Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames(var = "Symbol")

tpm_final <- as.matrix(tpm_final)
mode(tpm_final) <- "numeric"

cat("===== Gene-level TPM matrix dimension =====\n")
print(dim(tpm_final))

write.csv(tpm_final, file = "expressionmetrix_GSE118565.csv")
save(tpm_final, file = "GSE118565_step1_workspace.Rdata")

cat("===== All sample names =====\n")
print(colnames(tpm_final))

# ==============================
# 5. Select samples
# original selection: columns 1,2,7,8,13,14
# intended regrouping:
# 1,7,13 = Control
# 2,8,14 = CX5461
# ==============================
dat <- tpm_final[, c(1, 2, 7, 8, 13, 14), drop = FALSE]
dat <- dat[, c(1, 3, 5, 2, 4, 6), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(dat))

# ==============================
# 6. Filter low-expression genes
# For TPM matrix, a practical filter is:
# TPM > 1 in at least 3 samples
# ==============================
keep <- rowSums(dat > 1, na.rm = TRUE) >= 3
dat_filt <- dat[keep, , drop = FALSE]

cat("===== Before filtering =====\n")
print(dim(dat))
cat("===== After filtering (TPM > 1 in >= 3 samples) =====\n")
print(dim(dat_filt))

# ==============================
# 7. Log2 transform TPM
# Since this is processed TPM, use log2(TPM + 1) directly
# ==============================
exp_matrix <- log2(dat_filt + 1)

cat("===== Overall range after log2(TPM+1) =====\n")
print(range(exp_matrix, na.rm = TRUE))

cat("===== Overall quantiles after log2(TPM+1) =====\n")
print(quantile(
  exp_matrix,
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
))

sample_summary <- data.frame(
  Sample = colnames(exp_matrix),
  Min = apply(exp_matrix, 2, min, na.rm = TRUE),
  Q1 = apply(exp_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(exp_matrix, 2, median, na.rm = TRUE),
  Mean = apply(exp_matrix, 2, mean, na.rm = TRUE),
  Q3 = apply(exp_matrix, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(exp_matrix, 2, max, na.rm = TRUE)
)

cat("===== Per-sample summary =====\n")
print(sample_summary)
write.csv(sample_summary, "GSE118565_sample_summary.csv", row.names = FALSE)

write.table(
  data.frame(
    Dataset = "GSE118565",
    InputType = "processed TPM matrix",
    ReNormalization = "No",
    Filtering = "TPM > 1 in at least 3 samples",
    Transformation = "log2(TPM + 1)",
    DEGmethod = "limma",
    Note = "Exploratory DEG from processed TPM, not raw-count voom analysis"
  ),
  file = "GSE118565_processing_decision.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pdf("GSE118565_boxplot_used_for_DEG.pdf", width = 10, height = 6)
boxplot(
  exp_matrix,
  las = 2,
  outline = FALSE,
  main = "GSE118565 log2(TPM+1) used for DEG"
)
dev.off()

pdf("GSE118565_density_used_for_DEG.pdf", width = 10, height = 6)
plotDensities(exp_matrix, main = "GSE118565 density plot used for DEG")
dev.off()

# ==============================
# 8. DEG analysis by limma
# ==============================
group <- factor(
  c(rep("Control", 3), rep("CX5461", 3)),
  levels = c("Control", "CX5461")
)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("===== Design matrix =====\n")
print(design)

fit <- lmFit(exp_matrix, design)

contrast_matrix <- makeContrasts(
  CX5461_vs_Control = CX5461 - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_diff <- topTable(
  fit2,
  coef = "CX5461_vs_Control",
  adjust.method = "fdr",
  number = Inf,
  sort.by = "P"
)

all_diff$GeneSymbol <- rownames(all_diff)
all_diff <- all_diff[, c("GeneSymbol", setdiff(colnames(all_diff), "GeneSymbol"))]

# FDR-based results
diffSig_fdr <- subset(all_diff, abs(logFC) > 1 & adj.P.Val < 0.05)
up_fdr <- subset(all_diff, logFC > 1 & adj.P.Val < 0.05)
down_fdr <- subset(all_diff, logFC < -1 & adj.P.Val < 0.05)

# Nominal P-value results
diffSig_p <- subset(all_diff, abs(logFC) > 1 & P.Value < 0.05)
up_p <- subset(all_diff, logFC > 1 & P.Value < 0.05)
down_p <- subset(all_diff, logFC < -1 & P.Value < 0.05)

cat("===== DEG summary =====\n")
cat("Total genes tested:", nrow(all_diff), "\n")

cat("FDR significant DEGs (|logFC| > 1 & adj.P.Val < 0.05):", nrow(diffSig_fdr), "\n")
cat("FDR upregulated genes:", nrow(up_fdr), "\n")
cat("FDR downregulated genes:", nrow(down_fdr), "\n")

cat("Nominal significant DEGs (|logFC| > 1 & P.Value < 0.05):", nrow(diffSig_p), "\n")
cat("Nominal upregulated genes:", nrow(up_p), "\n")
cat("Nominal downregulated genes:", nrow(down_p), "\n")

write.csv(all_diff, "GSE118565_all_DEG.csv", row.names = FALSE)

write.csv(diffSig_fdr, "GSE118565_sig_DEG_FDR.csv", row.names = FALSE)
write.csv(up_fdr, "GSE118565_up_DEG_FDR.csv", row.names = FALSE)
write.csv(down_fdr, "GSE118565_down_DEG_FDR.csv", row.names = FALSE)

write.csv(diffSig_p, "GSE118565_sig_DEG_Pvalue.csv", row.names = FALSE)
write.csv(up_p, "GSE118565_up_DEG_Pvalue.csv", row.names = FALSE)
write.csv(down_p, "GSE118565_down_DEG_Pvalue.csv", row.names = FALSE)

save(
  tpm,
  tpm_final,
  dat,
  dat_filt,
  exp_matrix,
  sample_summary,
  all_diff,
  diffSig_fdr,
  up_fdr,
  down_fdr,
  diffSig_p,
  up_p,
  down_p,
  file = "GSE118565_analysis_workspace.Rdata"
)



up_genes <- unique(na.omit(up_fdr$GeneSymbol))
down_genes <- unique(na.omit(down_fdr$GeneSymbol))
gene_list_long <- data.frame(
  dataset = c(rep("GSE118565", length(up_genes)), rep("GSE118565", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)
write.csv(gene_list_long, "GSE118565_diffgenesymbol_long.csv", row.names = FALSE)


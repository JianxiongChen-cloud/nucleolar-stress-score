# GSE267501
setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE267501/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(readxl)
library(limma)

# ==============================
# 1. Read FPKM table
# ==============================
data <- read_excel("GSE267501_HCT116_TP53WT_fpkm.xlsx")
data <- as.data.frame(data)

cat("===== Raw table dimension =====\n")
print(dim(data))
cat("===== Column names =====\n")
print(colnames(data))

# ==============================
# 2. Select expression columns
# original order:
# 6h_rep1-3 first, then 0h_rep1-3
# ==============================
dat <- data[, c(5:7, 2:4), drop = FALSE]
dat <- as.matrix(dat)
mode(dat) <- "numeric"
rownames(dat) <- data$gene_id

cat("===== Selected samples =====\n")
print(colnames(dat))

cat("===== Raw FPKM range =====\n")
print(range(dat, na.rm = TRUE))

cat("===== Raw FPKM quantiles =====\n")
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

write.csv(sample_summary_raw, "GSE267501_sample_summary_rawFPKM.csv", row.names = FALSE)

# ==============================
# 3. Filter low-expression genes
# FPKM > 1 in at least 3 samples
# ==============================
keep <- rowSums(dat > 1, na.rm = TRUE) >= 3
dat_filt <- dat[keep, , drop = FALSE]

cat("===== Filtering summary =====\n")
cat("Genes before filtering:", nrow(dat), "\n")
cat("Genes after filtering:", nrow(dat_filt), "\n")

# ==============================
# 4. Log2 transform
# ==============================
exp_matrix <- log2(dat_filt + 1)

cat("===== log2(FPKM+1) range =====\n")
print(range(exp_matrix, na.rm = TRUE))

cat("===== log2(FPKM+1) quantiles =====\n")
print(quantile(
  exp_matrix,
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
))

sample_summary_log <- data.frame(
  Sample = colnames(exp_matrix),
  Min = apply(exp_matrix, 2, min, na.rm = TRUE),
  Q1 = apply(exp_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(exp_matrix, 2, median, na.rm = TRUE),
  Mean = apply(exp_matrix, 2, mean, na.rm = TRUE),
  Q3 = apply(exp_matrix, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(exp_matrix, 2, max, na.rm = TRUE)
)

write.csv(sample_summary_log, "GSE267501_sample_summary_logFPKM.csv", row.names = FALSE)

pdf("GSE267501_boxplot_logFPKM_used_for_DEG.pdf", width = 10, height = 6)
boxplot(
  exp_matrix,
  las = 2,
  outline = FALSE,
  main = "GSE267501 log2(FPKM+1) used for DEG"
)
dev.off()

pdf("GSE267501_density_logFPKM_used_for_DEG.pdf", width = 10, height = 6)
plotDensities(
  exp_matrix,
  main = "GSE267501 density plot"
)
dev.off()

# remove zero-variance genes
gene_var <- apply(exp_matrix, 1, var, na.rm = TRUE)
exp_matrix <- exp_matrix[gene_var > 0, , drop = FALSE]

write.table(
  data.frame(
    Dataset = "GSE267501",
    InputType = "processed FPKM matrix",
    ReNormalization = "No",
    Filtering = "FPKM > 1 in at least 3 samples",
    Transformation = "log2(FPKM + 1)",
    DEGmethod = "limma",
    Note = "Exploratory DEG analysis based on processed FPKM"
  ),
  file = "GSE267501_processing_decision.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ==============================
# 5. DEG analysis
# first 3 = BMH21 (6h), last 3 = CONTROL (0h)
# ==============================
group <- factor(
  c(rep("BMH21", 3), rep("CONTROL", 3)),
  levels = c("CONTROL", "BMH21")
)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("===== Design matrix =====\n")
print(design)

fit <- lmFit(exp_matrix, design)

contrast_matrix <- makeContrasts(
  BMH21_vs_CONTROL = BMH21 - CONTROL,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_diff <- topTable(
  fit2,
  coef = "BMH21_vs_CONTROL",
  adjust.method = "fdr",
  number = Inf,
  sort.by = "P"
)

all_diff$GeneSymbol <- rownames(all_diff)
all_diff <- all_diff[, c("GeneSymbol", setdiff(colnames(all_diff), "GeneSymbol"))]

# strict: FDR
diffSig_fdr <- subset(all_diff, abs(logFC) > 1 & adj.P.Val < 0.05)
up_fdr <- subset(all_diff, logFC > 1 & adj.P.Val < 0.05)
down_fdr <- subset(all_diff, logFC < -1 & adj.P.Val < 0.05)

# exploratory: nominal P
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

write.csv(all_diff, "GSE267501_all_DEG.csv", row.names = FALSE)

write.csv(diffSig_fdr, "GSE267501_sig_DEG_FDR.csv", row.names = FALSE)
write.csv(up_fdr, "GSE267501_up_DEG_FDR.csv", row.names = FALSE)
write.csv(down_fdr, "GSE267501_down_DEG_FDR.csv", row.names = FALSE)

write.csv(diffSig_p, "GSE267501_sig_DEG_Pvalue.csv", row.names = FALSE)
write.csv(up_p, "GSE267501_up_DEG_Pvalue.csv", row.names = FALSE)
write.csv(down_p, "GSE267501_down_DEG_Pvalue.csv", row.names = FALSE)

save(
  data, dat, dat_filt, exp_matrix,
  sample_summary_raw, sample_summary_log,
  all_diff,
  diffSig_fdr, up_fdr, down_fdr,
  diffSig_p, up_p, down_p,
  file = "GSE267501_analysis_workspace.Rdata"
)



up_genes <- unique(na.omit(up_fdr$GeneSymbol))
down_genes <- unique(na.omit(down_fdr$GeneSymbol))
gene_list_long <- data.frame(
  dataset = c(rep("GSE267501", length(up_genes)), rep("GSE267501", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)
write.csv(gene_list_long, "GSE267501_diffgenesymbol_long.csv", row.names = FALSE)
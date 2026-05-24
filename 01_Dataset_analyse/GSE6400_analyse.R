setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE6400/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(readr)
library(hgu133plus2.db)
library(limma)
library(AnnotationDbi)

# ==============================
# 1. Read GEO series matrix
# ==============================
expr_probe <- read_table(
  "GSE6400_series_matrix.txt.gz",
  comment = "!",
  show_col_types = FALSE
)

expr_probe <- as.data.frame(expr_probe, check.names = FALSE)
rownames(expr_probe) <- expr_probe[, 1]
expr_probe <- expr_probe[, -1, drop = FALSE]

rownames(expr_probe) <- gsub('"', "", rownames(expr_probe))
colnames(expr_probe) <- gsub('"', "", colnames(expr_probe))

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
  hgu133plus2.db,
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

# Ensure that ids and expr_probe_annot are in the same order.
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

# ==============================
# 5. Select samples
#    You already checked sample order; keep 1:6 directly
# ==============================
expr_gene <- expr_gene[, 1:6, drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene))

# ==============================
# 6. Check expression scale
#    Use a more standard GEO/limma-style rule
# ==============================
cat("===== Overall range =====\n")
print(range(expr_gene, na.rm = TRUE))

cat("===== Overall quantiles =====\n")
qx <- quantile(
  as.numeric(expr_gene),
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
)
print(qx)

sample_summary <- data.frame(
  Sample = colnames(expr_gene),
  Min = apply(expr_gene, 2, min, na.rm = TRUE),
  Q1 = apply(expr_gene, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(expr_gene, 2, median, na.rm = TRUE),
  Mean = apply(expr_gene, 2, mean, na.rm = TRUE),
  Q3 = apply(expr_gene, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(expr_gene, 2, max, na.rm = TRUE)
)

cat("===== Per-sample summary =====\n")
print(sample_summary)
write.csv(sample_summary, "GSE6400_sample_summary.csv", row.names = FALSE)

pdf("GSE6400_boxplot_before_logcheck.pdf", width = 10, height = 6)
boxplot(
  expr_gene,
  las = 2,
  outline = FALSE,
  main = "GSE6400 expression before log-check"
)
dev.off()

# Commonly Used Empirical Judgments for GEO2R/limma
qx6 <- quantile(
  as.numeric(expr_gene),
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
  expr_gene[expr_gene <= 0] <- NA
  exp_matrix <- log2(expr_gene)
  log_status <- "log2_transformed_in_this_script"
} else {
  cat("This matrix is likely already on log scale. No additional log2 transformation will be applied.\n")
  exp_matrix <- expr_gene
  log_status <- "already_log_scale_no_additional_log2"
}

overall_max <- max(expr_gene, na.rm = TRUE)
overall_median <- median(expr_gene, na.rm = TRUE)

write.table(
  data.frame(
    Dataset = "GSE6400",
    Platform = "Affymetrix Human Genome U133 Plus 2.0 Array",
    Processing = "RMA normalization (GEO annotation)",
    InputType = "processed expression matrix",
    ReNormalization = "No",
    LogDecision = log_status,
    MaxValue = overall_max,
    MedianValue = overall_median,
    DEGmethod = "limma",
    ProbeCollapse = "highest median expression"
  ),
  file = "GSE6400_logcheck_decision.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pdf("GSE6400_boxplot_used_for_DEG.pdf", width = 10, height = 6)
boxplot(
  exp_matrix,
  las = 2,
  outline = FALSE,
  main = "GSE6400 expression used for DEG"
)
dev.off()

# Optional: to see if the distribution is more reasonable after logging.
pdf("GSE6400_density_used_for_DEG.pdf", width = 10, height = 6)
plotDensities(exp_matrix, main = "GSE6400 density plot used for DEG")
dev.off()

# ==============================
# 7. DEG analysis by limma
# ==============================
group <- factor(
  c(rep("ActD", 3), rep("CONTROL", 3)),
  levels = c("CONTROL", "ActD")
)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("===== Design matrix =====\n")
print(design)

fit <- lmFit(exp_matrix, design)

contrast_matrix <- makeContrasts(
  ActD_vs_CONTROL = ActD - CONTROL,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_diff <- topTable(
  fit2,
  coef = "ActD_vs_CONTROL",
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

write.csv(all_diff, "GSE6400_all_DEG.csv", row.names = FALSE)
write.csv(diffSig, "GSE6400_sig_DEG.csv", row.names = FALSE)

# ==============================
# 8. Save workspace
# ==============================
save(
  expr_probe,
  expr_probe_annot,
  expr_gene,
  exp_matrix,
  ids,
  sample_summary,
  all_diff,
  diffSig,
  file = "GSE6400_analysis_workspace.Rdata"
)

up_genes <- unique(na.omit(up_diff$GeneSymbol))
down_genes <- unique(na.omit(down_diff$GeneSymbol))
gene_list_long <- data.frame(
  dataset = c(rep("GSE6400", length(up_genes)), rep("GSE6400", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)
write.csv(gene_list_long, "GSE6400_diffgenesymbol_long.csv", row.names = FALSE)

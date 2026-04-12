setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE12666/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(GEOquery)
library(stringr)
library(hgu133a2.db)
library(limma)
library(AnnotationDbi)

# ==============================
# 1. Download / load GEO ExpressionSet
# ==============================
f <- "GSE12666_eSet.Rdata"
if (!file.exists(f)) {
  gset <- getGEO("GSE12666", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
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
# 2. Extract phenotype / group information
# ==============================
pd <- pData(eset)
group_list <- str_split(pd$title, "_", simplify = TRUE)[, 1]

cat("===== group_list table =====\n")
print(table(group_list))

# ==============================
# 3. Probe annotation
# ==============================
ids <- AnnotationDbi::select(
  hgu133a2.db,
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
# 4. Keep matched probes
# ==============================
expr_probe_annot <- expr_probe[ids$probe_id, , drop = FALSE]
ids <- ids[match(rownames(expr_probe_annot), ids$probe_id), ]

cat("===== Matched annotated matrix dimension =====\n")
print(dim(expr_probe_annot))

# ==============================
# 5. Collapse duplicated probes
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

write.csv(expr_gene, file = "expressionmetrix_GSE12666.csv")

save(
  expr_probe,
  expr_probe_annot,
  expr_gene,
  group_list,
  file = "GSE12666_step1_workspace.Rdata"
)

# ==============================
# 6. Select samples for DEG analysis
#    original selection: columns 13:16 and 25:28
# ==============================
expr_gene_sub <- expr_gene[, c(13:16, 25:28), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

# ==============================
# 7. Check expression scale
# ==============================
cat("===== Overall range =====\n")
print(range(expr_gene_sub, na.rm = TRUE))

cat("===== Overall quantiles =====\n")
qx <- quantile(
  as.numeric(expr_gene_sub),
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
)
print(qx)

sample_summary <- data.frame(
  Sample = colnames(expr_gene_sub),
  Min = apply(expr_gene_sub, 2, min, na.rm = TRUE),
  Q1 = apply(expr_gene_sub, 2, quantile, probs = 0.25, na.rm = TRUE),
  Median = apply(expr_gene_sub, 2, median, na.rm = TRUE),
  Mean = apply(expr_gene_sub, 2, mean, na.rm = TRUE),
  Q3 = apply(expr_gene_sub, 2, quantile, probs = 0.75, na.rm = TRUE),
  Max = apply(expr_gene_sub, 2, max, na.rm = TRUE)
)

cat("===== Per-sample summary =====\n")
print(sample_summary)
write.csv(sample_summary, "GSE12666_sample_summary.csv", row.names = FALSE)

pdf("GSE12666_boxplot_before_logcheck.pdf", width = 10, height = 6)
boxplot(
  expr_gene_sub,
  las = 2,
  outline = FALSE,
  main = "GSE12666 expression before log-check"
)
dev.off()

# GEO2R / limma style scale judgement
qx6 <- quantile(
  as.numeric(expr_gene_sub),
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
  expr_gene_sub[expr_gene_sub <= 0] <- NA
  exp_matrix <- log2(expr_gene_sub)
  log_status <- "log2_transformed_in_this_script"
} else {
  cat("This matrix is likely already on log scale. No additional log2 transformation will be applied.\n")
  exp_matrix <- expr_gene_sub
  log_status <- "already_log_scale_no_additional_log2"
}

overall_max <- max(expr_gene_sub, na.rm = TRUE)
overall_median <- median(expr_gene_sub, na.rm = TRUE)

write.table(
  data.frame(
    Dataset = "GSE12666",
    Platform = "Affymetrix Human Genome U133A 2.0 Array",
    Processing = "GCOS 1.4.1 + global scaling (target intensity 150)",
    InputType = "processed expression matrix",
    ReNormalization = "No",
    LogDecision = log_status,
    MaxValue = overall_max,
    MedianValue = overall_median,
    DEGmethod = "limma",
    ProbeCollapse = "highest median expression"
  ),
  file = "GSE12666_logcheck_decision.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pdf("GSE12666_boxplot_used_for_DEG.pdf", width = 10, height = 6)
boxplot(
  exp_matrix,
  las = 2,
  outline = FALSE,
  main = "GSE12666 expression used for DEG"
)
dev.off()

pdf("GSE12666_density_used_for_DEG.pdf", width = 10, height = 6)
plotDensities(exp_matrix, main = "GSE12666 density plot used for DEG")
dev.off()

# ==============================
# 8. DEG analysis
#    4 BMH21 vs 4 MOCK
# ==============================
group <- factor(
  c(rep("BMH21", 4), rep("MOCK", 4)),
  levels = c("MOCK", "BMH21")
)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cat("===== Design matrix =====\n")
print(design)

fit <- lmFit(exp_matrix, design)

contrast_matrix <- makeContrasts(
  BMH21_vs_MOCK = BMH21 - MOCK,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

all_diff <- topTable(
  fit2,
  coef = "BMH21_vs_MOCK",
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

write.csv(all_diff, "GSE12666_all_DEG.csv", row.names = FALSE)
write.csv(diffSig, "GSE12666_sig_DEG.csv", row.names = FALSE)
write.csv(up_diff, "GSE12666_up_DEG.csv", row.names = FALSE)
write.csv(down_diff, "GSE12666_down_DEG.csv", row.names = FALSE)

# ==============================
# 9. Save workspace
# ==============================
save(
  expr_probe,
  expr_probe_annot,
  expr_gene,
  expr_gene_sub,
  exp_matrix,
  ids,
  sample_summary,
  all_diff,
  diffSig,
  up_diff,
  down_diff,
  group_list,
  file = "GSE12666_analysis_workspace.Rdata"
)

up_genes <- unique(na.omit(up_diff$GeneSymbol))
down_genes <- unique(na.omit(down_diff$GeneSymbol))
gene_list_long <- data.frame(
  dataset = c(rep("GSE12666", length(up_genes)), rep("GSE12666", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)
write.csv(gene_list_long, "GSE12666_diffgenesymbol_long.csv", row.names = FALSE)
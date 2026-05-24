setwd("/home/xxm_xxm/CJX_workspace/Dataset/GSE62593/")
rm(list = ls())
options(stringsAsFactors = FALSE)

library(GEOquery)
library(limma)
library(readr)

# ==============================
# 1. Download / load GEO ExpressionSet
# ==============================
f <- "GSE62593_eSet.Rdata"
if (!file.exists(f)) {
  gset <- getGEO("GSE62593", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
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

cat("===== First few feature IDs =====\n")
print(rownames(expr_probe)[1:min(10, nrow(expr_probe))])

# ==============================
# 2. Read platform annotation file
#    Note:
#    This dataset is described as RMA at gene level.
#    So annotation is only used to map feature ID -> gene symbol.
# ==============================
anno <- read_table("GPL5175-3188.txt", comment = "#", show_col_types = FALSE)
anno <- as.data.frame(anno)
# Please first confirm in the platform's comment file:
# The ID column is the feature/probe ID
# The category column is indeed the gene symbol
ids <- data.frame(
  probe_id = as.character(anno$ID),
  symbol = as.character(anno$category),
  stringsAsFactors = FALSE
)

# remove invalid annotations
ids <- ids[
  !is.na(ids$symbol) &
    ids$symbol != "" &
    !grepl("//", ids$symbol, fixed = TRUE),
]

ids <- unique(ids)
ids <- ids[ids$probe_id %in% rownames(expr_probe), ]

cat("===== Annotated features =====\n")
print(nrow(ids))

# ==============================
# 3. Keep matched features
# ==============================
expr_annot <- expr_probe[ids$probe_id, , drop = FALSE]
ids <- ids[match(rownames(expr_annot), ids$probe_id), ]

cat("===== Matched annotated matrix dimension =====\n")
print(dim(expr_annot))

# ==============================
# 4. Collapse duplicated mapped symbols
#    Since matrix is already gene-level by processing description,
#    this step is only to handle repeated mappings in annotation.
#    Strategy: keep feature with highest median expression
# ==============================
ids$median_expr <- apply(expr_annot, 1, median, na.rm = TRUE)
ids <- ids[order(ids$symbol, -ids$median_expr), ]
ids <- ids[!duplicated(ids$symbol), ]

expr_gene <- expr_annot[ids$probe_id, , drop = FALSE]
rownames(expr_gene) <- ids$symbol

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

write.csv(expr_gene, file = "expressionmetrix_GSE62593.csv")

save(
  expr_probe,
  expr_annot,
  expr_gene,
  ids,
  file = "GSE62593_step1_workspace.Rdata"
)

# ==============================
# 5. Select samples
#    original selection: columns 1:3 and 7:9
# ==============================
expr_gene_sub <- expr_gene[, c(1:3, 7:9), drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene_sub))

# ==============================
# 6. Check expression scale
#    This dataset is expected to already be log2-scaled by RMA
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
write.csv(sample_summary, "GSE62593_sample_summary.csv", row.names = FALSE)

pdf("GSE62593_boxplot_before_DEG.pdf", width = 10, height = 6)
boxplot(
  expr_gene_sub,
  las = 2,
  outline = FALSE,
  main = "GSE62593 expression used for DEG"
)
dev.off()

pdf("GSE62593_density_used_for_DEG.pdf", width = 10, height = 6)
plotDensities(expr_gene_sub, main = "GSE62593 density plot used for DEG")
dev.off()

# For this dataset, keep original matrix directly
exp_matrix <- expr_gene_sub
log_status <- "already_log_scale_RMA_gene_level_no_additional_log2"

write.table(
  data.frame(
    Dataset = "GSE62593",
    Processing = "RMA on gene level using oligo",
    InputType = "processed gene-level expression matrix",
    ReNormalization = "No",
    LogDecision = log_status,
    MaxValue = max(expr_gene_sub, na.rm = TRUE),
    MedianValue = median(expr_gene_sub, na.rm = TRUE),
    DEGmethod = "limma",
    AnnotationHandling = "mapped to gene symbol; duplicated mappings collapsed by highest median expression"
  ),
  file = "GSE62593_logcheck_decision.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ==============================
# 7. DEG analysis
#    columns 1:3 = CONTROL
#    columns 7:9 = ActD
# ==============================
group <- factor(
  c(rep("CONTROL", 3), rep("ActD", 3)),
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

diffSig <- subset(all_diff, abs(logFC) > 1 & P.Value < 0.05)
up_diff <- subset(all_diff, logFC > 1 & P.Value < 0.05)
down_diff <- subset(all_diff, logFC < -1 & P.Value < 0.05)

cat("===== DEG summary =====\n")
cat("Total genes tested:", nrow(all_diff), "\n")
cat("Total significant DEGs:", nrow(diffSig), "\n")
cat("Upregulated genes:", nrow(up_diff), "\n")
cat("Downregulated genes:", nrow(down_diff), "\n")

write.csv(all_diff, "GSE62593_all_DEG.csv", row.names = FALSE)
write.csv(diffSig, "GSE62593_sig_DEG.csv", row.names = FALSE)
write.csv(up_diff, "GSE62593_up_DEG.csv", row.names = FALSE)
write.csv(down_diff, "GSE62593_down_DEG.csv", row.names = FALSE)

save(
  expr_probe,
  expr_annot,
  expr_gene,
  expr_gene_sub,
  exp_matrix,
  ids,
  sample_summary,
  all_diff,
  diffSig,
  up_diff,
  down_diff,
  file = "GSE62593_analysis_workspace.Rdata"
)

up_genes <- unique(na.omit(up_diff$GeneSymbol))
down_genes <- unique(na.omit(down_diff$GeneSymbol))
gene_list_long <- data.frame(
  dataset = c(rep("GSE62593", length(up_genes)), rep("GSE62593", length(down_genes))),
  group = c(rep("UP", length(up_genes)), rep("DOWN", length(down_genes))),
  gene_symbol = c(up_genes, down_genes)
)
write.csv(gene_list_long, "GSE62593_diffgenesymbol_long.csv", row.names = FALSE)

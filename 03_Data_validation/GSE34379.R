rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE34379"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"
data_file <- file.path(out_dir, "GSE34379_series_matrix.txt.gz")
annot_file <- file.path(out_dir, "GPL6244.annot.gz")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)
library(readr)

############################################################
## 2. Manually read the GSE34379 series matrix expression data
############################################################
lines <- readLines(gzfile(data_file), warn = FALSE)

begin_idx <- grep("^!series_matrix_table_begin", lines)
end_idx   <- grep("^!series_matrix_table_end", lines)

if (length(begin_idx) != 1 || length(end_idx) != 1) {
  stop("Cannot uniquely identify !series_matrix_table_begin / end in the file.")
}

expr_lines <- lines[(begin_idx + 1):(end_idx - 1)]

tmp_txt <- tempfile(fileext = ".txt")
writeLines(expr_lines, tmp_txt)

expr_df_raw <- read.delim(
  tmp_txt,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

unlink(tmp_txt)

## Remove leading/trailing quotation marks from column names
colnames(expr_df_raw) <- gsub('^"|"$', "", colnames(expr_df_raw))

cat("===== Column names after cleanup =====\n")
print(colnames(expr_df_raw))

if (!"ID_REF" %in% colnames(expr_df_raw)) {
  stop("ID_REF column not found in expression matrix after header cleanup.")
}

rownames(expr_df_raw) <- trimws(as.character(expr_df_raw$ID_REF))
expr_df_raw$ID_REF <- NULL

expr_mat <- as.matrix(expr_df_raw)
mode(expr_mat) <- "numeric"

cat("===== Expression matrix dimension =====\n")
print(dim(expr_mat))

cat("===== Sample names =====\n")
print(colnames(expr_mat))

cat("===== First few probe IDs =====\n")
print(head(rownames(expr_mat)))

############################################################
## 3. Sample group information
############################################################
sample_map <- data.frame(
  geo_accession = c(
    "GSM847867", "GSM847868", "GSM847869", "GSM847870",
    "GSM847871", "GSM847872", "GSM847873", "GSM847874"
  ),
  Sample = c(
    "CNT1", "CNT2", "CNT3", "CNT4",
    "Carbo1", "Carbo2", "Carbo3", "Carbo4"
  ),
  Group = c(
    rep("Control", 4),
    rep("Carbo", 4)
  ),
  stringsAsFactors = FALSE
)

if (!all(sample_map$geo_accession %in% colnames(expr_mat))) {
  stop("Some GEO sample IDs in sample_map were not found in expression matrix columns.")
}

expr_mat <- expr_mat[, sample_map$geo_accession, drop = FALSE]
colnames(expr_mat) <- sample_map$Sample

group_info <- data.frame(
  Sample = sample_map$Sample,
  Group = factor(sample_map$Group, levels = c("Control", "Carbo")),
  stringsAsFactors = FALSE
)

cat("===== Group information =====\n")
print(group_info)

############################################################
## 4. Manually read the GPL6244 platform annotation table
############################################################
annot_lines <- readLines(gzfile(annot_file), warn = FALSE)

begin_idx <- grep("^!platform_table_begin", annot_lines)
end_idx   <- grep("^!platform_table_end", annot_lines)

if (length(begin_idx) != 1 || length(end_idx) != 1) {
  stop("Cannot uniquely identify !platform_table_begin / end in GPL annotation file.")
}

table_lines <- annot_lines[(begin_idx + 1):(end_idx - 1)]

tmp_annot <- tempfile(fileext = ".txt")
writeLines(table_lines, tmp_annot)

annot_raw <- read.delim(
  tmp_annot,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

unlink(tmp_annot)

cat("===== Annotation dimension =====\n")
print(dim(annot_raw))

cat("===== Annotation columns =====\n")
print(colnames(annot_raw))

############################################################
## 5. Extract probe -> gene symbol mapping
##    from the GPL6244 annotation table
############################################################
if (!all(c("ID", "Gene symbol") %in% colnames(annot_raw))) {
  stop("Columns 'ID' and/or 'Gene symbol' not found in GPL6244 annotation table.")
}

gene_annot <- annot_raw[, c("ID", "Gene symbol")]
colnames(gene_annot) <- c("PROBEID", "SYMBOL")

gene_annot$PROBEID <- trimws(as.character(gene_annot$PROBEID))
gene_annot$SYMBOL  <- trimws(as.character(gene_annot$SYMBOL))

## Remove missing values
gene_annot <- gene_annot[
  !is.na(gene_annot$PROBEID) & gene_annot$PROBEID != "" &
    !is.na(gene_annot$SYMBOL) & gene_annot$SYMBOL != "",
]

## If one probe corresponds to multiple symbols,
## keep only the first symbol
gene_annot$SYMBOL <- sub("///.*$", "", gene_annot$SYMBOL)
gene_annot$SYMBOL <- toupper(trimws(gene_annot$SYMBOL))

## Keep only probes present in the expression matrix
gene_annot <- gene_annot[gene_annot$PROBEID %in% rownames(expr_mat), ]

## Remove duplicated probes
gene_annot <- gene_annot[!duplicated(gene_annot$PROBEID), ]

cat("===== Cleaned annotation =====\n")
print(dim(gene_annot))
print(head(gene_annot))

############################################################
## 6. Merge expression matrix and generate gene-level matrix
############################################################
expr_df <- data.frame(
  PROBEID = rownames(expr_mat),
  expr_mat,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expr_annot <- merge(gene_annot, expr_df, by = "PROBEID")
expr_annot <- expr_annot[!duplicated(expr_annot$PROBEID), ]

cat("===== Merged probe-level data =====\n")
print(dim(expr_annot))


## For genes with multiple probes,
## use the mean expression value
expr_gene <- expr_annot %>%
  dplyr::select(-PROBEID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

expr_gene <- as.data.frame(expr_gene)
rownames(expr_gene) <- expr_gene$SYMBOL
expr_gene <- expr_gene[, -1, drop = FALSE]
expr_gene <- as.matrix(expr_gene)
mode(expr_gene) <- "numeric"

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

cat("===== First few genes =====\n")
print(head(rownames(expr_gene)))

############################################################
## 7. Read the full nucleolar stress gene set
############################################################
read_gene_vector <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, 1, drop = TRUE]
  x <- unique(toupper(trimws(x)))
  x <- x[!is.na(x) & x != ""]
  return(x)
}

full_up   <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_up_genes.csv"))
full_down <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_down_genes.csv"))

full_up_in   <- intersect(full_up, rownames(expr_gene))
full_down_in <- intersect(full_down, rownames(expr_gene))

cat("===== Full gene set overlap =====\n")
cat("NuS_Up:", length(full_up_in), "\n")
cat("NuS_Down:", length(full_down_in), "\n")

if (length(full_up_in) < 5 || length(full_down_in) < 5) {
  stop("Too few overlapping genes for the full gene set.")
}

geneSets_full <- list(
  NuS_Up = full_up_in,
  NuS_Down = full_down_in
)

############################################################
## 8. Calculate NuS scores using ssGSEA
############################################################
param <- ssgseaParam(
  exprData = expr_gene,
  geneSets = geneSets_full,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_res <- gsva(param)

if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea_res))) {
  stop("NuS_Up or NuS_Down not found in ssGSEA result.")
}

nus_score <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]

scores_all <- data.frame(
  Sample = colnames(expr_gene),
  Score = as.numeric(nus_score),
  Method = "ssGSEA",
  GeneSet = "Full",
  stringsAsFactors = FALSE
)

scores_all <- left_join(scores_all, group_info, by = "Sample")
scores_all$Group <- factor(scores_all$Group, levels = c("Control", "Carbo"))

write.csv(scores_all, "GSE34379_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 9. Statistical testing
############################################################
test_res <- data.frame(
  comparison = "Carbo vs Control",
  p_welch = tryCatch(
    t.test(Score ~ Group, data = scores_all, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  ),
  mean_control = mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  sd_control   = sd(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  n_control    = sum(scores_all$Group == "Control"),
  mean_treat   = mean(scores_all$Score[scores_all$Group == "Carbo"], na.rm = TRUE),
  sd_treat     = sd(scores_all$Score[scores_all$Group == "Carbo"], na.rm = TRUE),
  n_treat      = sum(scores_all$Group == "Carbo"),
  diff_mean    = mean(scores_all$Score[scores_all$Group == "Carbo"], na.rm = TRUE) -
    mean(scores_all$Score[scores_all$Group == "Control"], na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(test_res, "GSE34379_NuS_score_statistics_ttest.csv", row.names = FALSE)

cat("===== Statistical test result =====\n")
print(test_res)

############################################################
## 10. Visualization
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "Carbo"   = "#51B3D1"
)

p <- ggplot(scores_all, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.8, alpha = 0.85) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "NuS ssGSEA scores in GSE34379",
    x = NULL,
    y = "NuS"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("GSE34379_NuS_ssGSEA_violinplot.pdf", p, width = 4, height = 5, bg = "white")
ggsave("GSE34379_NuS_ssGSEA_violinplot.png", p, width = 4, height = 5, dpi = 300, bg = "white")

############################################################
## 11. Save workspace
############################################################
save(
  expr_mat,
  expr_gene,
  group_info,
  gene_annot,
  geneSets_full,
  scores_all,
  test_res,
  file = "GSE34379_NuS_ssGSEA_workspace.Rdata"
)
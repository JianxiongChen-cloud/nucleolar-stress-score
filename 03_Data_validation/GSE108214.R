rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE108214"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"
data_file <- file.path(out_dir, "GSE108214_series_matrix.txt.gz")
annot_file <- file.path(out_dir, "GPL17077-17467.txt")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(GSVA)
library(ggplot2)
library(dplyr)
library(readr)
library(multcomp)

############################################################
## 2. Manually read the series matrix expression data
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

## Clean quotation marks from column names and probe IDs
colnames(expr_df_raw) <- gsub('^"|"$', "", colnames(expr_df_raw))

if (!"ID_REF" %in% colnames(expr_df_raw)) {
  stop("ID_REF column not found in expression matrix.")
}

expr_df_raw$ID_REF <- gsub('^"|"$', "", trimws(as.character(expr_df_raw$ID_REF)))
rownames(expr_df_raw) <- expr_df_raw$ID_REF
expr_df_raw$ID_REF <- NULL

expr_mat_all <- as.matrix(expr_df_raw)
mode(expr_mat_all) <- "numeric"

cat("===== Raw expression matrix dimension =====\n")
print(dim(expr_mat_all))

cat("===== Raw sample names =====\n")
print(colnames(expr_mat_all))

cat("===== First few probe IDs =====\n")
print(head(rownames(expr_mat_all)))

############################################################
## 3. Sample group information
##    grouped according to biological replicates
############################################################
sample_map <- data.frame(
  geo_accession = c(
    "GSM2892607","GSM2892608","GSM2892609","GSM2892610","GSM2892611",
    "GSM2892612","GSM2892613","GSM2892614","GSM2892615","GSM2892616",
    "GSM2892617","GSM2892618","GSM2892619","GSM2892620","GSM2892621",
    "GSM2892622","GSM2892623","GSM2892624","GSM2892625","GSM2892626",
    "GSM2892627","GSM2892628"
  ),
  Sample = c(
    "A549sens_ctrl_1","A549sens_11uM_1","A549res_ctrl_1","A549res_11uM_1","A549res_34uM_1",
    "A549sens_ctrl_2","A549sens_11uM_2","A549res_ctrl_2","A549sens_ctrl_3","A549sens_11uM_3",
    "A549res_ctrl_3","A549res_11uM_3","A549res_34uM_3","A549res_11uM_2","A549res_34uM_2",
    "A549sens_ctrl_4","A549res_ctrl_4","A549res_11uM_4","A549res_34uM_4","A549res_ctrl_5",
    "A549res_11uM_5","A549res_34uM_5"
  ),
  Group = c(
    "A549sens_ctrl","A549sens_11uM","A549res_ctrl","A549res_11uM","A549res_34uM",
    "A549sens_ctrl","A549sens_11uM","A549res_ctrl","A549sens_ctrl","A549sens_11uM",
    "A549res_ctrl","A549res_11uM","A549res_34uM","A549res_11uM","A549res_34uM",
    "A549sens_ctrl","A549res_ctrl","A549res_11uM","A549res_34uM","A549res_ctrl",
    "A549res_11uM","A549res_34uM"
  ),
  stringsAsFactors = FALSE
)

if (!all(sample_map$geo_accession %in% colnames(expr_mat_all))) {
  stop("Some GEO sample IDs in sample_map were not found in expression matrix columns.")
}

expr_mat <- expr_mat_all[, sample_map$geo_accession, drop = FALSE]
colnames(expr_mat) <- sample_map$Sample

group_info <- data.frame(
  Sample = sample_map$Sample,
  Group = factor(
    sample_map$Group,
    levels = c(
      "A549sens_ctrl",
      "A549sens_11uM",
      "A549res_ctrl",
      "A549res_11uM",
      "A549res_34uM"
    )
  ),
  stringsAsFactors = FALSE
)

cat("===== Selected samples =====\n")
print(colnames(expr_mat))

cat("===== Group information =====\n")
print(group_info)

############################################################
## 4. Read the GPL17077-17467 annotation table
##    (skip the first 16 description lines)
############################################################
annot_raw <- read.delim(
  annot_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  skip = 16
)

cat("===== Annotation dimension =====\n")
print(dim(annot_raw))

cat("===== Annotation columns =====\n")
print(colnames(annot_raw))

cat("===== First 6 rows =====\n")
print(head(annot_raw[, 1:min(10, ncol(annot_raw))]))

############################################################
## 5. Extract probe -> gene symbol annotation
############################################################
if (!all(c("ID", "GENE_SYMBOL") %in% colnames(annot_raw))) {
  stop("Columns 'ID' and/or 'GENE_SYMBOL' not found in annotation file.")
}

gene_annot <- annot_raw[, c("ID", "GENE_SYMBOL")]
colnames(gene_annot) <- c("PROBEID", "SYMBOL")

gene_annot$PROBEID <- gsub('^"|"$', "", trimws(as.character(gene_annot$PROBEID)))
gene_annot$SYMBOL  <- gsub('^"|"$', "", trimws(as.character(gene_annot$SYMBOL)))

gene_annot <- gene_annot[
  !is.na(gene_annot$PROBEID) & gene_annot$PROBEID != "" &
    !is.na(gene_annot$SYMBOL) & gene_annot$SYMBOL != "",
]

## When multiple gene annotations exist,
## keep only the first symbol
gene_annot$SYMBOL <- sub("///.*$", "", gene_annot$SYMBOL)
gene_annot$SYMBOL <- sub("//.*$", "", gene_annot$SYMBOL)
gene_annot$SYMBOL <- toupper(trimws(gene_annot$SYMBOL))

## Keep only probes present in the expression matrix
gene_annot <- gene_annot[gene_annot$PROBEID %in% rownames(expr_mat), ]
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
scores_all$Group <- factor(
  scores_all$Group,
  levels = c(
    "A549sens_ctrl",
    "A549sens_11uM",
    "A549res_ctrl",
    "A549res_11uM",
    "A549res_34uM"
  )
)

write.csv(scores_all, "GSE108214_NuS_scores_ssGSEA.csv", row.names = FALSE)

cat("===== NuS scores =====\n")
print(scores_all)

############################################################
## 9. Statistical testing:
##    one-way ANOVA + Tukey multiple comparison
############################################################
anova_fit <- aov(Score ~ Group, data = scores_all)
anova_res <- summary(anova_fit)
print(anova_res)

anova_table <- as.data.frame(anova_res[[1]])
anova_table$term <- rownames(anova_table)
rownames(anova_table) <- NULL
write.csv(anova_table, "GSE108214_NuS_score_ANOVA.csv", row.names = FALSE)

tukey_res <- TukeyHSD(anova_fit, "Group")
tukey_df <- as.data.frame(tukey_res$Group)
tukey_df$comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL
colnames(tukey_df) <- c("diff", "lwr", "upr", "p_adj", "comparison")
write.csv(tukey_df, "GSE108214_NuS_score_TukeyHSD.csv", row.names = FALSE)

cat("===== Tukey multiple comparison =====\n")
print(tukey_df)

############################################################
## 10. Visualization
############################################################
fill_colors <- c(
  "A549sens_ctrl" = "#C4E7C1",
  "A549sens_11uM" = "#93D4BC",
  "A549res_ctrl"  = "#51B3D1",
  "A549res_11uM"  = "#2C7FB8",
  "A549res_34uM"  = "#08589E"
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
    title = "NuS ssGSEA scores in GSE108214",
    x = NULL,
    y = "NuS"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("GSE108214_NuS_ssGSEA_violinplot.pdf", p, width = 5.5, height = 5, bg = "white")
ggsave("GSE108214_NuS_ssGSEA_violinplot.png", p, width = 5.5, height = 5, dpi = 300, bg = "white")

############################################################
## 11. Save workspace
############################################################
save(
  expr_mat_all,
  expr_mat,
  expr_gene,
  group_info,
  gene_annot,
  geneSets_full,
  scores_all,
  anova_fit,
  tukey_df,
  file = "GSE108214_NuS_ssGSEA_workspace.Rdata"
)
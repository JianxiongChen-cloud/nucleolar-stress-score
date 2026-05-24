###############################################################################
## Task 1: Pan-cancer analysis of RiboSis and NuS ########Step3
## Focus: robust GSVA scoring with clear progress messages
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 1e6)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
  library(GSVA)
  library(matrixStats)
  library(forcats)
  library(scales)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"
processed_dir <- file.path(base_dir, "processed_data")

task_dir  <- file.path(base_dir, "task1_RiboSis_NuS_analysis")
plot_dir  <- file.path(task_dir, "plots")
table_dir <- file.path(task_dir, "tables")
rds_dir   <- file.path(task_dir, "rds")
qc_dir    <- file.path(task_dir, "qc")
log_dir   <- file.path(task_dir, "logs")

dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

expr_file    <- file.path(processed_dir, "expr_main_tcga_tumor_gtex_normal.rds")
annot_file   <- file.path(processed_dir, "annot_main_tcga_tumor_gtex_normal_with_clinical.rds")
nus_rdata    <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

stopifnot(file.exists(expr_file))
stopifnot(file.exists(annot_file))
stopifnot(file.exists(nus_rdata))
stopifnot(file.exists(ribosis_rdata))

###############################################################################
## 1. Helper functions
###############################################################################
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

msg <- function(...) {
  cat(sprintf("[%s] ", timestamp()), ..., "\n", sep = "")
  flush.console()
}

save_csv <- function(df, filename) {
  write.csv(df, file.path(table_dir, filename), row.names = FALSE)
}

save_rds <- function(obj, filename) {
  saveRDS(obj, file.path(rds_dir, filename))
}

save_qc_csv <- function(df, filename) {
  write.csv(df, file.path(qc_dir, filename), row.names = FALSE)
}

load_rdata_to_list <- function(rdata_path) {
  e <- new.env()
  load(rdata_path, envir = e)
  as.list(e)
}

assert_numeric_matrix <- function(mat, object_name = "matrix") {
  if (!is.matrix(mat)) stop(object_name, " is not a matrix.")
  if (!is.numeric(mat)) stop(object_name, " is not numeric.")
  invisible(TRUE)
}

###############################################################################
## 2. Load data
###############################################################################
msg("Step 1/7: Loading expression and annotation ...")

expr_mat <- readRDS(expr_file)
annot <- readRDS(annot_file)

expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"
assert_numeric_matrix(expr_mat, "expr_mat")

stopifnot(all(colnames(expr_mat) %in% annot$sample))
annot <- annot %>% arrange(match(sample, colnames(expr_mat)))
stopifnot(identical(colnames(expr_mat), annot$sample))

msg("Expression matrix loaded: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
msg("Annotation loaded: ", nrow(annot), " rows x ", ncol(annot), " columns")

###############################################################################
## 3. Load gene sets
###############################################################################
msg("Step 2/7: Loading NuS and RiboSis gene sets ...")

nus_obj <- load_rdata_to_list(nus_rdata)
ribo_obj <- load_rdata_to_list(ribosis_rdata)

if (!"geneSets_final" %in% names(nus_obj)) {
  stop("Object 'geneSets_final' not found in NuStress_geneSets_final.Rdata")
}
if (!all(c("NuStress_UP", "NuStress_DOWN") %in% names(nus_obj$geneSets_final))) {
  stop("NuStress_UP / NuStress_DOWN not found in geneSets_final")
}
if (!"ribosis" %in% names(ribo_obj)) {
  stop("Object 'ribosis' not found in RiboSis activity.Rdata")
}
if (!is.list(ribo_obj$ribosis) || !"set" %in% names(ribo_obj$ribosis)) {
  stop("Cannot find ribosis$set in RiboSis activity.Rdata")
}

nus_sets_raw <- list(
  NuStress_UP   = unique(as.character(nus_obj$geneSets_final$NuStress_UP)),
  NuStress_DOWN = unique(as.character(nus_obj$geneSets_final$NuStress_DOWN))
)

ribo_sets_raw <- list(
  RiboSis = unique(as.character(ribo_obj$ribosis$set))
)

msg("Raw gene sets loaded: NuStress_UP = ", length(nus_sets_raw$NuStress_UP),
    ", NuStress_DOWN = ", length(nus_sets_raw$NuStress_DOWN),
    ", RiboSis = ", length(ribo_sets_raw$RiboSis))

###############################################################################
## 4. Convert Ensembl IDs to gene symbols
###############################################################################
msg("Step 3/7: Converting Ensembl IDs to gene symbols ...")

# IMPORTANT: remove version suffix correctly
ensembl_ids <- sub("\\..*$", "", rownames(expr_mat))

gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

gene_annot <- gene_annot %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

map_idx <- match(ensembl_ids, gene_annot$ENSEMBL)
mapped_symbol <- gene_annot$SYMBOL[map_idx]

keep_mapped <- !is.na(mapped_symbol) & mapped_symbol != ""
expr_mapped <- expr_mat[keep_mapped, , drop = FALSE]
mapped_symbol <- mapped_symbol[keep_mapped]

msg("Mapped rows retained: ", nrow(expr_mapped), " / ", nrow(expr_mat))

# collapse duplicate SYMBOLs by highest median expression
row_med <- matrixStats::rowMedians(expr_mapped, na.rm = TRUE)

collapse_df <- data.frame(
  row_index = seq_len(nrow(expr_mapped)),
  SYMBOL = mapped_symbol,
  median_expr = row_med,
  stringsAsFactors = FALSE
) %>%
  arrange(SYMBOL, desc(median_expr)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

expr_mat_symbol <- expr_mapped[collapse_df$row_index, , drop = FALSE]
rownames(expr_mat_symbol) <- collapse_df$SYMBOL
storage.mode(expr_mat_symbol) <- "numeric"
assert_numeric_matrix(expr_mat_symbol, "expr_mat_symbol")

msg("Expression matrix after SYMBOL conversion: ",
    nrow(expr_mat_symbol), " genes x ", ncol(expr_mat_symbol), " samples")

save_rds(expr_mat_symbol, "expr_mat_symbol_collapsed.rds")

mapping_table <- data.frame(
  ENSEMBL = ensembl_ids[keep_mapped][collapse_df$row_index],
  SYMBOL  = collapse_df$SYMBOL,
  stringsAsFactors = FALSE
)
save_csv(mapping_table, "ensembl_to_symbol_mapping_table.csv")

###############################################################################
## 5. Match gene sets to expression matrix
###############################################################################
msg("Step 4/7: Matching gene sets to SYMBOL-space expression matrix ...")

clean_gene_set <- function(gs, expr_genes) {
  gs <- unique(gs)
  gs <- gs[!is.na(gs) & gs != ""]
  intersect(gs, expr_genes)
}

nus_sets <- lapply(nus_sets_raw, clean_gene_set, expr_genes = rownames(expr_mat_symbol))
ribo_sets <- lapply(ribo_sets_raw, clean_gene_set, expr_genes = rownames(expr_mat_symbol))

msg("Matched genes: NuStress_UP = ", length(nus_sets$NuStress_UP),
    ", NuStress_DOWN = ", length(nus_sets$NuStress_DOWN),
    ", RiboSis = ", length(ribo_sets$RiboSis))

if (length(nus_sets$NuStress_UP) < 5) stop("Too few NuStress_UP genes matched.")
if (length(nus_sets$NuStress_DOWN) < 5) stop("Too few NuStress_DOWN genes matched.")
if (length(ribo_sets$RiboSis) < 5) stop("Too few RiboSis genes matched.")

gene_set_qc <- tibble(
  GeneSet = c("NuStress_UP", "NuStress_DOWN", "RiboSis"),
  N_genes_raw = c(length(nus_sets_raw$NuStress_UP),
                  length(nus_sets_raw$NuStress_DOWN),
                  length(ribo_sets_raw$RiboSis)),
  N_genes_matched = c(length(nus_sets$NuStress_UP),
                      length(nus_sets$NuStress_DOWN),
                      length(ribo_sets$RiboSis))
)

save_csv(gene_set_qc, "gene_set_qc_summary.csv")
save_rds(nus_sets, "nus_gene_sets_used.rds")
save_rds(ribo_sets, "ribosis_gene_set_used.rds")

###############################################################################
## 6. Clean expression matrix before GSVA
###############################################################################
msg("Step 5/7: Cleaning expression matrix before GSVA ...")

expr_mat_symbol_clean <- expr_mat_symbol
expr_mat_symbol_clean[!is.finite(expr_mat_symbol_clean)] <- NA_real_

msg("Any NA in matrix? ", any(is.na(expr_mat_symbol_clean)))
msg("Total NA count: ", sum(is.na(expr_mat_symbol_clean)))

row_non_na <- rowSums(!is.na(expr_mat_symbol_clean))

msg("Genes with all NA: ", sum(row_non_na == 0))
msg("Genes with < 2 non-NA values: ", sum(row_non_na < 2))

expr_mat_symbol_clean <- expr_mat_symbol_clean[row_non_na >= 2, , drop = FALSE]

msg("After removing low-information genes: ",
    nrow(expr_mat_symbol_clean), " genes x ", ncol(expr_mat_symbol_clean), " samples")

row_meds <- matrixStats::rowMedians(expr_mat_symbol_clean, na.rm = TRUE)

na_idx <- which(is.na(expr_mat_symbol_clean), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  expr_mat_symbol_clean[na_idx] <- row_meds[na_idx[, 1]]
}

msg("Remaining NA after imputation: ", sum(is.na(expr_mat_symbol_clean)))

row_sd <- matrixStats::rowSds(expr_mat_symbol_clean, na.rm = TRUE)

msg("Genes with SD > 0: ", sum(row_sd > 0, na.rm = TRUE))
msg("Genes with SD == 0: ", sum(row_sd == 0, na.rm = TRUE))

expr_mat_symbol_var <- expr_mat_symbol_clean[row_sd > 0 & !is.na(row_sd), , drop = FALSE]
storage.mode(expr_mat_symbol_var) <- "numeric"

msg("Final matrix for GSVA: ",
    nrow(expr_mat_symbol_var), " genes x ", ncol(expr_mat_symbol_var), " samples")

saveRDS(expr_mat_symbol_clean, file.path(rds_dir, "expr_mat_symbol_clean.rds"))
saveRDS(expr_mat_symbol_var, file.path(rds_dir, "expr_mat_symbol_var_for_GSVA.rds"))

###############################################################################
## 7. Rematch gene sets after cleaning
###############################################################################
msg("Step 6/7: Rematching gene sets after cleaning ...")

nus_sets <- lapply(nus_sets, function(gs) {
  gs <- gs[!is.na(gs) & gs != ""]
  intersect(gs, rownames(expr_mat_symbol_var))
})

ribo_sets <- lapply(ribo_sets, function(gs) {
  gs <- gs[!is.na(gs) & gs != ""]
  intersect(gs, rownames(expr_mat_symbol_var))
})

msg("NuStress_UP genes after cleaning: ", length(nus_sets$NuStress_UP))
msg("NuStress_DOWN genes after cleaning: ", length(nus_sets$NuStress_DOWN))
msg("RiboSis genes after cleaning: ", length(ribo_sets$RiboSis))

gs_list <- list(
  NuStress_UP   = nus_sets$NuStress_UP,
  NuStress_DOWN = nus_sets$NuStress_DOWN,
  RiboSis       = ribo_sets$RiboSis
)

stopifnot(length(gs_list$NuStress_UP) >= 5)
stopifnot(length(gs_list$NuStress_DOWN) >= 5)
stopifnot(length(gs_list$RiboSis) >= 5)

###############################################################################
## 8. ssGSEA scoring sample-by-sample with progress messages
###############################################################################
msg("Step 7/7: Running ssGSEA scoring sample-by-sample ...")
msg("This step may take a while on 17,599 samples.")

dir.create(file.path(rds_dir, "ssgsea_by_sample"), recursive = TRUE, showWarnings = FALSE)

sample_ids <- colnames(expr_mat_symbol_var)
n_samples <- length(sample_ids)

msg("Total samples: ", n_samples)
msg("Gene sets to score: ", paste(names(gs_list), collapse = ", "))

run_single_sample_ssgsea <- function(expr_mat_full, sample_id, gs_list) {
  expr_one <- expr_mat_full[, sample_id, drop = FALSE]
  storage.mode(expr_one) <- "numeric"
  
  sample_scores <- numeric(length(gs_list))
  names(sample_scores) <- names(gs_list)
  
  for (set_name in names(gs_list)) {
    gene_set <- gs_list[[set_name]]
    
    par_obj <- GSVA::ssgseaParam(
      exprData = expr_one,
      geneSets = list(tmp_set = gene_set),
      alpha = 0.25,
      normalize = FALSE
    )
    
    score <- GSVA::gsva(par_obj, verbose = FALSE)
    sample_scores[set_name] <- as.numeric(score[1, 1])
  }
  
  return(sample_scores)
}

score_mat <- matrix(
  NA_real_,
  nrow = length(gs_list),
  ncol = n_samples,
  dimnames = list(names(gs_list), sample_ids)
)

pb <- txtProgressBar(min = 0, max = n_samples, style = 3)

t_all_start <- Sys.time()

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  t0 <- Sys.time()
  
  msg("------------------------------------------------------------")
  msg("Sample ", i, "/", n_samples, ": ", sid)
  
  out_file <- file.path(
    rds_dir, "ssgsea_by_sample",
    paste0("ssgsea_", sid, ".rds")
  )
  
  if (file.exists(out_file)) {
    msg("  -> Cached result found, loading.")
    sample_scores <- readRDS(out_file)
  } else {
    sample_scores <- run_single_sample_ssgsea(
      expr_mat_full = expr_mat_symbol_var,
      sample_id = sid,
      gs_list = gs_list
    )
    saveRDS(sample_scores, out_file)
    msg("  -> Result saved: ", basename(out_file))
  }
  
  score_mat[, sid] <- sample_scores
  
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  elapsed_min <- round(as.numeric(difftime(Sys.time(), t_all_start, units = "mins")), 2)
  avg_sec <- round(as.numeric(difftime(Sys.time(), t_all_start, units = "secs")) / i, 2)
  remain_min <- round(avg_sec * (n_samples - i) / 60, 2)
  
  msg("  -> Sample finished in ", dt, " sec")
  msg("  -> Elapsed: ", elapsed_min, " min; Estimated remaining: ", remain_min, " min")
  msg("  -> Scores: ",
      paste(names(sample_scores), round(sample_scores, 4), sep = "=", collapse = "; "))
  
  setTxtProgressBar(pb, i)
}
close(pb)

ssgsea_scores <- score_mat[, sample_ids, drop = FALSE]

msg("ssGSEA finished successfully.")
msg("Final ssGSEA score matrix: ", nrow(ssgsea_scores), " x ", ncol(ssgsea_scores))

saveRDS(ssgsea_scores, file.path(rds_dir, "ssgsea_scores_raw_by_sample.rds"))

###############################################################################
## 9. Construct final scores and merge into annotation
###############################################################################
msg("Step 8/8: Constructing final NuStress / RiboSis scores ...")

score_df <- data.frame(
  sample = colnames(ssgsea_scores),
  NuStress_UP = as.numeric(ssgsea_scores["NuStress_UP", ]),
  NuStress_DOWN = as.numeric(ssgsea_scores["NuStress_DOWN", ]),
  RiboSis = as.numeric(ssgsea_scores["RiboSis", ]),
  stringsAsFactors = FALSE
)

## Standardize each score across all samples before constructing composite NuStress
score_df <- score_df %>%
  mutate(
    NuStress_UP_z = as.numeric(scale(NuStress_UP)),
    NuStress_DOWN_z = as.numeric(scale(NuStress_DOWN)),
    RiboSis_z = as.numeric(scale(RiboSis)),
    NuStress = NuStress_UP_z - NuStress_DOWN_z,
    NuStress_z = as.numeric(scale(NuStress))
  )

msg("Score table constructed: ", nrow(score_df), " samples")
msg("NuStress summary: min = ", round(min(score_df$NuStress, na.rm = TRUE), 4),
    ", median = ", round(median(score_df$NuStress, na.rm = TRUE), 4),
    ", max = ", round(max(score_df$NuStress, na.rm = TRUE), 4))
msg("RiboSis_z summary: min = ", round(min(score_df$RiboSis_z, na.rm = TRUE), 4),
    ", median = ", round(median(score_df$RiboSis_z, na.rm = TRUE), 4),
    ", max = ", round(max(score_df$RiboSis_z, na.rm = TRUE), 4))

annot2 <- annot %>%
  left_join(score_df, by = "sample")

msg("Merged annotation table dimension: ", nrow(annot2), " x ", ncol(annot2))

save_csv(score_df, "ssgsea_scores_per_sample.csv")
saveRDS(score_df, file.path(rds_dir, "ssgsea_scores_per_sample.rds"))
saveRDS(annot2, file.path(rds_dir, "annot_with_ssgsea_scores.rds"))

msg("All done.")


dim(score_df)
head(score_df)

summary(score_df[, c("NuStress_UP", "NuStress_DOWN", "RiboSis",
                     "NuStress_UP_z", "NuStress_DOWN_z", "RiboSis_z",
                     "NuStress", "NuStress_z")])

cor(score_df[, c("NuStress_UP_z", "NuStress_DOWN_z", "RiboSis_z", "NuStress")],
    use = "pairwise.complete.obs", method = "spearman")



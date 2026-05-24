### Final version; use this script as the reference.
rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 6000)

###############################################################################
## 0. Paths
###############################################################################
workdir <- "/home/xxm_xxm/CJX_workspace/scRNAseq/"
setwd(workdir)

clus_file <- "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
meta_file <- "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"

nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

outdir  <- file.path(workdir, "Final_NuStress_scRNA_optimized")
plotdir <- file.path(outdir, "plots")
tabdir  <- file.path(outdir, "tables")
rdsdir  <- file.path(outdir, "rdata")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabdir,  showWarnings = FALSE, recursive = TRUE)
dir.create(rdsdir,  showWarnings = FALSE, recursive = TRUE)

###############################################################################
## 1. Libraries
###############################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggpubr)
  library(viridis)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(broom)
  library(scales)
  library(tidyr)
})

theme_set(
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

###############################################################################
## 2. Helper functions
###############################################################################
fig_counter <- 0

save_pdf <- function(plot_obj, filename_stub, width = 7, height = 5) {
  fig_counter <<- fig_counter + 1
  filename <- sprintf("Figure_%02d_%s.pdf", fig_counter, filename_stub)
  
  ggsave(
    filename = file.path(plotdir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = "pdf"
  )
  
  message("Saved: ", filename)
}

save_csv <- function(df, filename) {
  write.csv(df, file.path(tabdir, filename), row.names = FALSE)
}

upper_unique <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

extract_nus_sets <- function(env_obj) {
  obj_names <- ls(env_obj)
  
  if ("geneSets_all" %in% obj_names) {
    gs <- get("geneSets_all", envir = env_obj)
  } else if ("geneSets_final" %in% obj_names) {
    gs <- get("geneSets_final", envir = env_obj)
  } else if ("NuStress_geneSets_final" %in% obj_names) {
    gs <- get("NuStress_geneSets_final", envir = env_obj)
  } else {
    stop("Cannot find a valid NuStress gene set object in nus_rdata. Objects found: ",
         paste(obj_names, collapse = ", "))
  }
  
  if (!is.list(gs)) stop("NuStress object is not a list.")
  
  if (all(c("NuS_Up", "NuS_Down") %in% names(gs))) {
    up_genes   <- upper_unique(gs$NuS_Up)
    down_genes <- upper_unique(gs$NuS_Down)
  } else if (all(c("NuStress_UP", "NuStress_DOWN") %in% names(gs))) {
    up_genes   <- upper_unique(gs$NuStress_UP)
    down_genes <- upper_unique(gs$NuStress_DOWN)
  } else {
    stop("NuStress object found, but required elements are missing.")
  }
  
  list(NuS_Up = up_genes, NuS_Down = down_genes)
}

extract_ribosis_set <- function(env_obj) {
  obj_names <- ls(env_obj)
  
  if ("ribosis" %in% obj_names) {
    rb <- get("ribosis", envir = env_obj)
  } else {
    stop("Cannot find ribosis object in ribosis_rdata.")
  }
  
  if (is.list(rb) && "set" %in% names(rb)) {
    upper_unique(rb$set)
  } else {
    stop("ribosis object does not contain $set.")
  }
}

calc_module_score_simple <- function(obj, features, score_name, min_genes = 5) {
  DefaultAssay(obj) <- "RNA"
  
  features <- intersect(upper_unique(features), rownames(obj))
  
  if (length(features) < min_genes) {
    warning(score_name, ": available genes < ", min_genes)
    obj[[score_name]] <- NA_real_
    return(list(obj = obj, available = features))
  }
  
  need_norm <- FALSE
  tryCatch({
    dat <- GetAssayData(obj, assay = "RNA", slot = "data")
    if (nrow(dat) == 0 || ncol(dat) == 0) need_norm <- TRUE
  }, error = function(e) {
    need_norm <<- TRUE
  })
  
  if (need_norm) {
    obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  }
  
  obj <- AddModuleScore(
    object = obj,
    features = list(features),
    name = paste0(score_name, "_tmp"),
    assay = "RNA"
  )
  
  tmp_col <- paste0(score_name, "_tmp1")
  obj[[score_name]] <- obj@meta.data[[tmp_col]]
  obj@meta.data[[tmp_col]] <- NULL
  
  list(obj = obj, available = features)
}

make_cor_df <- function(meta, x, y) {
  df <- meta[, c(x, y), drop = FALSE]
  colnames(df) <- c("x", "y")
  df <- df[complete.cases(df), , drop = FALSE]
  df
}

run_spearman_test <- function(df) {
  if (nrow(df) < 3) {
    return(data.frame(rho = NA, p = NA, n = nrow(df)))
  }
  ct <- suppressWarnings(cor.test(df$x, df$y, method = "spearman"))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = nrow(df))
}

run_phase_specific_condition_test <- function(df, score_col = "NuStress",
                                              condition_col = "specimen_type",
                                              phase_col = "Phase") {
  phases <- unique(as.character(df[[phase_col]]))
  res <- lapply(phases, function(ph) {
    sub <- df[df[[phase_col]] == ph, , drop = FALSE]
    sub <- sub[!is.na(sub[[score_col]]) & !is.na(sub[[condition_col]]), , drop = FALSE]
    
    if (length(unique(sub[[condition_col]])) < 2) {
      return(data.frame(
        Phase = ph, group1 = NA, group2 = NA, n = nrow(sub),
        p_value = NA, median_group1 = NA, median_group2 = NA
      ))
    }
    
    lv <- sort(unique(as.character(sub[[condition_col]])))
    wt <- suppressWarnings(wilcox.test(sub[[score_col]] ~ sub[[condition_col]]))
    
    data.frame(
      Phase = ph,
      group1 = lv[1],
      group2 = lv[2],
      n = nrow(sub),
      p_value = wt$p.value,
      median_group1 = median(sub[sub[[condition_col]] == lv[1], score_col], na.rm = TRUE),
      median_group2 = median(sub[sub[[condition_col]] == lv[2], score_col], na.rm = TRUE)
    )
  })
  do.call(rbind, res)
}

plot_feature_umap <- function(obj, feature, title,
                              reduction_use = "umap",
                              pt_size = 3,
                              low_quantile = 0.01,
                              high_quantile = 0.99) {
  vals <- obj@meta.data[[feature]]
  vals <- vals[is.finite(vals)]
  
  vmin <- quantile(vals, low_quantile, na.rm = TRUE)
  vmax <- quantile(vals, high_quantile, na.rm = TRUE)
  
  p <- FeaturePlot(
    obj,
    features = feature,
    reduction = reduction_use,
    pt.size = pt_size,
    order = TRUE
  ) +
    scale_color_viridis_c(
      option = "magma",
      limits = c(vmin, vmax),
      oob = scales::squish
    ) +
    labs(title = title) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  p
}

###############################################################################
## 3. Load saved object and metadata
###############################################################################
cat("=== Step 1: Load saved object ===\n")

load("precessed_obj.Rdata")
clus <- fread(clus_file, data.table = FALSE)
meta <- fread(meta_file, data.table = FALSE)

if (!exists("seurat_obj")) stop("seurat_obj not found in precessed_obj.Rdata")

cat("Cells in Seurat object:", ncol(seurat_obj), "\n")
cat("Rows in cluster table:", nrow(clus), "\n")
cat("Rows in meta table:", nrow(meta), "\n")

###############################################################################
## 4. Attach metadata
###############################################################################
cat("=== Step 2: Attach metadata by row order ===\n")

if (ncol(seurat_obj) != nrow(clus)) {
  stop("Cell number in Seurat object does not match rows in clus.")
}
if (ncol(seurat_obj) != nrow(meta)) {
  stop("Cell number in Seurat object does not match rows in meta.")
}

if ("clTopLevel" %in% colnames(clus)) {
  seurat_obj$cluster <- clus$clTopLevel
} else {
  stop("Column clTopLevel not found in clus.")
}

if ("cl295v11SubFull" %in% colnames(clus)) {
  seurat_obj$epi_subcluster_full <- clus$cl295v11SubFull
} else {
  seurat_obj$epi_subcluster_full <- NA
  warning("Column cl295v11SubFull not found in clus.")
}

if ("SPECIMEN_TYPE" %in% colnames(meta)) {
  seurat_obj$specimen_type <- meta$SPECIMEN_TYPE
} else {
  stop("Column SPECIMEN_TYPE not found in meta.")
}

seurat_obj$clus_row_index <- seq_len(nrow(clus))
seurat_obj$meta_row_index <- seq_len(nrow(meta))

cat("Cluster distribution:\n")
print(table(seurat_obj$cluster, useNA = "ifany"))

cat("Specimen type distribution:\n")
print(table(seurat_obj$specimen_type, useNA = "ifany"))

###############################################################################
## 5. Basic QC field
###############################################################################
if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}

###############################################################################
## 6. Prepare expression data
###############################################################################
cat("=== Step 3: Prepare expression data ===\n")

rownames(seurat_obj) <- toupper(rownames(seurat_obj))

need_normalize <- FALSE
tryCatch({
  x <- GetAssayData(seurat_obj, slot = "data")
  if (nrow(x) == 0 || ncol(x) == 0) need_normalize <- TRUE
}, error = function(e) {
  need_normalize <<- TRUE
})

if (need_normalize) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
}

if (!"pca" %in% names(seurat_obj@reductions)) {
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
}
if (!"umap" %in% names(seurat_obj@reductions)) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
}

###############################################################################
## 7. Load gene sets
###############################################################################
cat("=== Step 4: Load NEW gene sets ===\n")

if (!file.exists(nus_rdata)) {
  stop("NuStress gene set file not found: ", nus_rdata)
}
if (!file.exists(ribosis_rdata)) {
  stop("RiboSis gene set file not found: ", ribosis_rdata)
}

nus_env <- new.env()
load(nus_rdata, envir = nus_env)
nus_sets <- extract_nus_sets(nus_env)

ribo_env <- new.env()
load(ribosis_rdata, envir = ribo_env)
ribosis_set <- extract_ribosis_set(ribo_env)

cat("NuS_Up genes:", length(nus_sets$NuS_Up), "\n")
cat("NuS_Down genes:", length(nus_sets$NuS_Down), "\n")
cat("RiboSis genes:", length(ribosis_set), "\n")

###############################################################################
## 8. Recalculate NuStress and RiboSis
###############################################################################
cat("=== Step 5: Recalculate NuStress and RiboSis ===\n")

old_cols <- c("NuS_Up1", "NuS_Up_score", "NuS_Down1", "NuS_Down_score",
              "NuS_combined_score", "NuStress", "RiboSis1", "RiboSis")
old_cols <- old_cols[old_cols %in% colnames(seurat_obj@meta.data)]
if (length(old_cols) > 0) {
  seurat_obj@meta.data[, old_cols] <- list(NULL)
}

sc1 <- calc_module_score_simple(seurat_obj, nus_sets$NuS_Up, "NuS_Up_score")
seurat_obj <- sc1$obj
available_Up <- sc1$available

sc2 <- calc_module_score_simple(seurat_obj, nus_sets$NuS_Down, "NuS_Down_score")
seurat_obj <- sc2$obj
available_Down <- sc2$available

sc3 <- calc_module_score_simple(seurat_obj, ribosis_set, "RiboSis")
seurat_obj <- sc3$obj
available_RiboSis <- sc3$available

seurat_obj$NuStress <- seurat_obj$NuS_Up_score - seurat_obj$NuS_Down_score
seurat_obj$NuS_combined_score <- seurat_obj$NuStress

gene_coverage_df <- data.frame(
  GeneSet = c("NuS_Up", "NuS_Down", "RiboSis"),
  Total_Genes = c(length(nus_sets$NuS_Up), length(nus_sets$NuS_Down), length(ribosis_set)),
  Available_Genes = c(length(available_Up), length(available_Down), length(available_RiboSis))
)
save_csv(gene_coverage_df, "GeneSet_Coverage.csv")

###############################################################################
## 9. Whole dataset plots
###############################################################################
cat("=== Step 6: Whole dataset plots ===\n")

custom_colors <- c(
  "B" = "#1f77b4",
  "Epi" = "#ff7f0e",
  "Mast" = "#2ca02c",
  "Myeloid" = "#d62728",
  "Plasma" = "#9467bd",
  "Strom" = "#8c564b",
  "TNKILC" = "#e377c2"
)

p_umap_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "cluster",
  cols = custom_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 3
) + labs(title = "UMAP of major cell types")
save_pdf(p_umap_cluster, "wholedata_umap_celltypes", width = 7.5, height = 6)

p_umap_nus <- plot_feature_umap(
  seurat_obj,
  feature = "NuStress",
  title = "UMAP of all cells: NuStress"
)
save_pdf(p_umap_nus, "wholedata_umap_NuStress", width = 7, height = 6)

p_umap_ribo <- plot_feature_umap(
  seurat_obj,
  feature = "RiboSis",
  title = "UMAP of all cells: RiboSis"
)
save_pdf(p_umap_ribo, "wholedata_umap_RiboSis", width = 7, height = 6)

comp_df <- seurat_obj@meta.data %>%
  dplyr::count(specimen_type, cluster) %>%
  group_by(specimen_type) %>%
  mutate(Fraction = n / sum(n))

p_comp_bar <- ggplot(comp_df, aes(x = specimen_type, y = Fraction, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Cell-type composition by specimen type",
    x = "Specimen type",
    y = "Fraction of cells"
  )
save_pdf(p_comp_bar, "wholedata_bar_celltype_composition_N_vs_T", width = 6.5, height = 5.5)

score_summary_whole <- seurat_obj@meta.data %>%
  group_by(cluster, specimen_type) %>%
  summarise(
    NuStress_mean = mean(NuStress, na.rm = TRUE),
    NuStress_se   = sd(NuStress, na.rm = TRUE) / sqrt(sum(!is.na(NuStress))),
    RiboSis_mean  = mean(RiboSis, na.rm = TRUE),
    RiboSis_se    = sd(RiboSis, na.rm = TRUE) / sqrt(sum(!is.na(RiboSis))),
    .groups = "drop"
  )

p_bar_nus_whole <- ggplot(score_summary_whole,
                          aes(x = cluster, y = NuStress_mean, fill = specimen_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = NuStress_mean - NuStress_se,
                    ymax = NuStress_mean + NuStress_se),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  labs(title = "Mean NuStress by cell type and specimen type",
       x = "Cell type", y = "Mean NuStress") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_bar_nus_whole, "wholedata_bar_mean_NuStress_by_celltype", width = 8, height = 5)

p_bar_ribo_whole <- ggplot(score_summary_whole,
                           aes(x = cluster, y = RiboSis_mean, fill = specimen_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = RiboSis_mean - RiboSis_se,
                    ymax = RiboSis_mean + RiboSis_se),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  labs(title = "Mean RiboSis by cell type and specimen type",
       x = "Cell type", y = "Mean RiboSis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_bar_ribo_whole, "wholedata_bar_mean_RiboSis_by_celltype", width = 8, height = 5)

p_box_nus_whole <- ggplot(
  seurat_obj@meta.data,
  aes(x = cluster, y = NuStress, fill = specimen_type)
) +
  geom_boxplot(
    width = 0.55,
    position = position_dodge(0.75),
    outlier.size = 0.25,
    color = "black"
  ) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  stat_compare_means(
    aes(group = specimen_type),
    method = "wilcox.test",
    label = "p.signif",
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "NuStress by cell type and specimen type",
    x = "Cell type",
    y = "NuStress"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_box_nus_whole, "wholedata_boxplot_NuStress_by_celltype", width = 9, height = 5.5)

p_box_ribo_whole <- ggplot(
  seurat_obj@meta.data,
  aes(x = cluster, y = RiboSis, fill = specimen_type)
) +
  geom_boxplot(
    width = 0.55,
    position = position_dodge(0.75),
    outlier.size = 0.25,
    color = "black"
  ) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  stat_compare_means(
    aes(group = specimen_type),
    method = "wilcox.test",
    label = "p.signif",
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "RiboSis by cell type and specimen type",
    x = "Cell type",
    y = "RiboSis"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_box_ribo_whole, "wholedata_boxplot_RiboSis_by_celltype", width = 9, height = 5.5)

###############################################################################
## 10. Extract epithelial cells
###############################################################################
cat("=== Step 7: Epithelial subset ===\n")

epi_cells <- colnames(seurat_obj)[seurat_obj$cluster == "Epi"]
cat("Number of Epi cells:", length(epi_cells), "\n")
if (length(epi_cells) == 0) stop("No Epi cells found in seurat_obj.")

epi_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[, epi_cells, drop = FALSE]
rownames(epi_counts) <- toupper(rownames(epi_counts))

dup_genes_epi <- duplicated(rownames(epi_counts))
if (sum(dup_genes_epi) > 0) {
  epi_counts <- epi_counts[!dup_genes_epi, , drop = FALSE]
}

subset_epi <- CreateSeuratObject(
  counts = epi_counts,
  project = "GSE178341_Epi"
)
DefaultAssay(subset_epi) <- "RNA"

epi_meta <- seurat_obj@meta.data[epi_cells, , drop = FALSE]
epi_meta <- epi_meta[colnames(subset_epi), , drop = FALSE]
subset_epi <- AddMetaData(subset_epi, metadata = epi_meta)

if ("epi_subcluster_full" %in% colnames(subset_epi@meta.data)) {
  subset_epi$epi_subcluster <- as.character(subset_epi$epi_subcluster_full)
} else {
  subset_epi$epi_subcluster <- "Unknown"
}

subset_epi <- NormalizeData(subset_epi, normalization.method = "LogNormalize")

sc1_epi <- calc_module_score_simple(subset_epi, nus_sets$NuS_Up, "NuS_Up_score")
subset_epi <- sc1_epi$obj
sc2_epi <- calc_module_score_simple(subset_epi, nus_sets$NuS_Down, "NuS_Down_score")
subset_epi <- sc2_epi$obj
sc3_epi <- calc_module_score_simple(subset_epi, ribosis_set, "RiboSis")
subset_epi <- sc3_epi$obj

subset_epi$NuStress <- subset_epi$NuS_Up_score - subset_epi$NuS_Down_score
subset_epi$NuS_combined_score <- subset_epi$NuStress

subset_epi <- FindVariableFeatures(subset_epi, selection.method = "vst", nfeatures = 2000)
subset_epi <- ScaleData(subset_epi, features = rownames(subset_epi))
subset_epi <- RunPCA(subset_epi, features = VariableFeatures(subset_epi))
subset_epi <- FindNeighbors(subset_epi, dims = 1:30)
subset_epi <- FindClusters(subset_epi, resolution = 0.8)
subset_epi <- RunUMAP(subset_epi, dims = 1:30)

epi_colors <- c(
  "cE01 (Stem/TA-like)" = "#1f77b4",
  "cE02 (Stem/TA-like/Immature Goblet)" = "#ff7f0e",
  "cE03 (Stem/TA-like prolif)" = "#2ca02c",
  "cE04 (Enterocyte 1)" = "#d62728",
  "cE05 (Enterocyte 2)" = "#9467bd",
  "cE06 (Immature Goblet)" = "#8c564b",
  "cE07 (Goblet/Enterocyte)" = "#e377c2",
  "cE08 (Goblet)" = "#7f7f7f",
  "cE09 (Best4)" = "#bcbd22",
  "cE10 (Tuft)" = "#17becf",
  "cE11 (Enteroendocrine)" = "#ff9896"
)

present_epi_levels <- unique(as.character(subset_epi$epi_subcluster))
epi_colors_use <- epi_colors[names(epi_colors) %in% present_epi_levels]

p_epi_umap <- DimPlot(
  subset_epi,
  reduction = "umap",
  group.by = "epi_subcluster",
  cols = epi_colors_use,
  label = TRUE,
  repel = TRUE,
  pt.size = 3
) + labs(title = "UMAP of epithelial subclusters")
save_pdf(p_epi_umap, "epi_umap_subclusters", width = 8, height = 5.5)

p_epi_nus <- plot_feature_umap(
  subset_epi,
  feature = "NuStress",
  title = "UMAP of epithelial cells: NuStress"
)
save_pdf(p_epi_nus, "epi_umap_NuStress", width = 5.5, height = 5)

p_epi_ribo <- plot_feature_umap(
  subset_epi,
  feature = "RiboSis",
  title = "UMAP of epithelial cells: RiboSis"
)
save_pdf(p_epi_ribo, "epi_umap_RiboSis", width = 5.5, height = 5)

###############################################################################
## 11. Tumor vs normal epithelial comparison
###############################################################################
cat("=== Step 8: Tumor vs normal epithelial comparison ===\n")

p_box_epi_nus <- ggplot(
  subset_epi@meta.data,
  aes(x = epi_subcluster, y = NuStress, fill = specimen_type)
) +
  geom_boxplot(width = 0.55, position = position_dodge(0.75),
               outlier.size = 0.25, color = "black") +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  stat_compare_means(
    aes(group = specimen_type),
    method = "wilcox.test",
    label = "p.signif",
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "NuStress in epithelial subclusters: normal vs tumor",
    x = "Epithelial subcluster",
    y = "NuStress"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_box_epi_nus, "epi_boxplot_NuStress_N_vs_T", width = 9, height = 5.5)

p_box_epi_ribo <- ggplot(
  subset_epi@meta.data,
  aes(x = epi_subcluster, y = RiboSis, fill = specimen_type)
) +
  geom_boxplot(width = 0.55, position = position_dodge(0.75),
               outlier.size = 0.25, color = "black") +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  stat_compare_means(
    aes(group = specimen_type),
    method = "wilcox.test",
    label = "p.signif",
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "RiboSis in epithelial subclusters: normal vs tumor",
    x = "Epithelial subcluster",
    y = "RiboSis"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_box_epi_ribo, "epi_boxplot_RiboSis_N_vs_T", width = 9, height = 5.5)

epi_comp_df <- subset_epi@meta.data %>%
  dplyr::count(specimen_type, epi_subcluster) %>%
  group_by(specimen_type) %>%
  mutate(Fraction = n / sum(n))

p_epi_comp_bar <- ggplot(epi_comp_df, aes(x = specimen_type, y = Fraction, fill = epi_subcluster)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = epi_colors_use) +
  labs(
    title = "Composition of epithelial subclusters by specimen type",
    x = "Specimen type",
    y = "Fraction of epithelial cells"
  )
save_pdf(p_epi_comp_bar, "epi_bar_subcluster_composition_N_vs_T", width = 7, height = 5.5)

epi_score_summary <- subset_epi@meta.data %>%
  group_by(epi_subcluster, specimen_type) %>%
  summarise(
    NuStress_mean = mean(NuStress, na.rm = TRUE),
    NuStress_se   = sd(NuStress, na.rm = TRUE) / sqrt(sum(!is.na(NuStress))),
    RiboSis_mean  = mean(RiboSis, na.rm = TRUE),
    RiboSis_se    = sd(RiboSis, na.rm = TRUE) / sqrt(sum(!is.na(RiboSis))),
    .groups = "drop"
  )

p_epi_bar_nus <- ggplot(epi_score_summary,
                        aes(x = epi_subcluster, y = NuStress_mean, fill = specimen_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = NuStress_mean - NuStress_se,
                    ymax = NuStress_mean + NuStress_se),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  labs(title = "Mean NuStress in epithelial subclusters",
       x = "Epithelial subcluster", y = "Mean NuStress") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_epi_bar_nus, "epi_bar_mean_NuStress_by_subcluster", width = 9, height = 5.5)

p_epi_bar_ribo <- ggplot(epi_score_summary,
                         aes(x = epi_subcluster, y = RiboSis_mean, fill = specimen_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = RiboSis_mean - RiboSis_se,
                    ymax = RiboSis_mean + RiboSis_se),
                position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  labs(title = "Mean RiboSis in epithelial subclusters",
       x = "Epithelial subcluster", y = "Mean RiboSis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_epi_bar_ribo, "epi_bar_mean_RiboSis_by_subcluster", width = 9, height = 5.5)

nu_min <- min(subset_epi$NuStress, na.rm = TRUE)
nu_max <- max(subset_epi$NuStress, na.rm = TRUE)

p_epi_nus_split <- FeaturePlot(
  subset_epi,
  features = "NuStress",
  reduction = "umap",
  split.by = "specimen_type",
  pt.size = 3,
  order = TRUE,
  keep.scale = "feature",
  min.cutoff = nu_min,
  max.cutoff = nu_max
) &
  scale_color_viridis_c(option = "magma", limits = c(nu_min, nu_max)) &
  theme(legend.position = "right")
save_pdf(p_epi_nus_split, "epi_umap_split_NuStress_by_specimen", width = 10, height = 5)

ribo_min <- min(subset_epi$RiboSis, na.rm = TRUE)
ribo_max <- max(subset_epi$RiboSis, na.rm = TRUE)

p_epi_ribo_split <- FeaturePlot(
  subset_epi,
  features = "RiboSis",
  reduction = "umap",
  split.by = "specimen_type",
  pt.size = 3,
  order = TRUE,
  keep.scale = "feature",
  min.cutoff = ribo_min,
  max.cutoff = ribo_max
) &
  scale_color_viridis_c(option = "magma", limits = c(ribo_min, ribo_max)) &
  theme(legend.position = "right")
save_pdf(p_epi_ribo_split, "epi_umap_split_RiboSis_by_specimen", width = 10, height = 5)

###############################################################################
## 12. Cell cycle confounding control
###############################################################################
cat("=== Step 9: Cell cycle confounding control ===\n")

s_genes   <- intersect(toupper(cc.genes.updated.2019$s.genes), rownames(subset_epi))
g2m_genes <- intersect(toupper(cc.genes.updated.2019$g2m.genes), rownames(subset_epi))

cat("Available S genes:", length(s_genes), "\n")
cat("Available G2M genes:", length(g2m_genes), "\n")

subset_epi <- CellCycleScoring(
  object = subset_epi,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

p_phase_umap <- DimPlot(
  subset_epi,
  reduction = "umap",
  group.by = "Phase",
  cols = c("G1" = "#4DBBD5", "S" = "#E64B35", "G2M" = "#3C5488"),
  pt.size = 3
) + labs(title = "Epithelial UMAP by cell-cycle phase")
save_pdf(p_phase_umap, "epi_umap_cellcycle_phase", width = 5.5, height = 5)

cor_s_df   <- make_cor_df(subset_epi@meta.data, "NuStress", "S.Score")
cor_g2m_df <- make_cor_df(subset_epi@meta.data, "NuStress", "G2M.Score")

cor_s_res   <- run_spearman_test(cor_s_df)
cor_g2m_res <- run_spearman_test(cor_g2m_df)

cor_summary <- bind_rows(
  data.frame(Comparison = "NuStress_vs_S.Score",   cor_s_res),
  data.frame(Comparison = "NuStress_vs_G2M.Score", cor_g2m_res)
)
save_csv(cor_summary, "NuStress_cellcycle_correlation_summary.csv")

p_cor_s <- ggplot(cor_s_df, aes(x = x, y = y)) +
  geom_point(size = 1, alpha = 0.35, color = "#E64B35") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  labs(
    title = paste0("NuStress vs S.Score (rho = ", round(cor_s_res$rho, 3),
                   ", p = ", signif(cor_s_res$p, 3), ")"),
    x = "NuStress", y = "S.Score"
  )
save_pdf(p_cor_s, "epi_scatter_NuStress_vs_SScore", width = 5, height = 5)

p_cor_g2m <- ggplot(cor_g2m_df, aes(x = x, y = y)) +
  geom_point(size = 1, alpha = 0.35, color = "#3C5488") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  labs(
    title = paste0("NuStress vs G2M.Score (rho = ", round(cor_g2m_res$rho, 3),
                   ", p = ", signif(cor_g2m_res$p, 3), ")"),
    x = "NuStress", y = "G2M.Score"
  )
save_pdf(p_cor_g2m, "epi_scatter_NuStress_vs_G2MScore", width = 5, height = 5)

phase_df <- subset_epi@meta.data %>%
  dplyr::select(NuStress, RiboSis, Phase, specimen_type, epi_subcluster, S.Score, G2M.Score) %>%
  filter(!is.na(NuStress), !is.na(RiboSis), !is.na(Phase))

phase_kw_nus  <- kruskal.test(NuStress ~ Phase, data = phase_df)
phase_kw_ribo <- kruskal.test(RiboSis ~ Phase, data = phase_df)

phase_kw_res <- data.frame(
  Score = c("NuStress", "RiboSis"),
  p_value = c(phase_kw_nus$p.value, phase_kw_ribo$p.value)
)
save_csv(phase_kw_res, "cellcycle_kruskal_summary.csv")

p_phase_box_nus <- ggplot(phase_df, aes(x = Phase, y = NuStress, fill = Phase)) +
  geom_boxplot(outlier.size = 0.2, color = "black") +
  scale_fill_manual(values = c("G1" = "#4DBBD5", "S" = "#E64B35", "G2M" = "#3C5488")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(phase_df$NuStress, na.rm = TRUE) * 1.05) +
  labs(title = "NuStress across cell-cycle phases", x = "Phase", y = "NuStress")
save_pdf(p_phase_box_nus, "epi_boxplot_NuStress_by_cellcycle_phase", width = 5.5, height = 5)

p_phase_box_ribo <- ggplot(phase_df, aes(x = Phase, y = RiboSis, fill = Phase)) +
  geom_boxplot(outlier.size = 0.2, color = "black") +
  scale_fill_manual(values = c("G1" = "#4DBBD5", "S" = "#E64B35", "G2M" = "#3C5488")) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(phase_df$RiboSis, na.rm = TRUE) * 1.05) +
  labs(title = "RiboSis across cell-cycle phases", x = "Phase", y = "RiboSis")
save_pdf(p_phase_box_ribo, "epi_boxplot_RiboSis_by_cellcycle_phase", width = 5.5, height = 5)

###############################################################################
## 12.1 Linear models controlling for S/G2M
###############################################################################
lm_df <- phase_df
lm_df$Phase <- factor(lm_df$Phase, levels = c("G1", "S", "G2M"))
lm_df$specimen_type <- factor(lm_df$specimen_type, levels = c("N", "T"))

fit1 <- lm(NuStress ~ specimen_type + S.Score + G2M.Score, data = lm_df)
fit2 <- lm(NuStress ~ specimen_type + S.Score + G2M.Score + epi_subcluster, data = lm_df)

save_csv(broom::tidy(fit1), "LM_NuStress_specimenType_S_G2M.csv")
save_csv(broom::tidy(fit2), "LM_NuStress_specimenType_S_G2M_epiSubcluster.csv")

###############################################################################
## 12.2 Phase-stratified N vs T comparison
###############################################################################
phase_condition_res <- run_phase_specific_condition_test(
  df = lm_df,
  score_col = "NuStress",
  condition_col = "specimen_type",
  phase_col = "Phase"
)
save_csv(phase_condition_res, "NuStress_phase_stratified_N_vs_T.csv")

phase_stat_plot <- phase_condition_res

phase_stat_plot$p.signif <- dplyr::case_when(
  is.na(phase_stat_plot$p_value) ~ "NA",
  phase_stat_plot$p_value < 0.0001 ~ "****",
  phase_stat_plot$p_value < 0.001  ~ "***",
  phase_stat_plot$p_value < 0.01   ~ "**",
  phase_stat_plot$p_value < 0.05   ~ "*",
  TRUE ~ "ns"
)

y_pos_df <- lm_df %>%
  dplyr::group_by(Phase) %>%
  dplyr::summarise(
    y.position = max(NuStress, na.rm = TRUE) + 0.04,
    .groups = "drop"
  )

phase_stat_plot <- phase_stat_plot %>%
  dplyr::left_join(y_pos_df, by = "Phase") %>%
  dplyr::mutate(
    group1 = "N",
    group2 = "T",
    xmin = c(0.8, 1.8, 2.8)[match(Phase, c("G1", "S", "G2M"))],
    xmax = c(1.2, 2.2, 3.2)[match(Phase, c("G1", "S", "G2M"))]
  )

p_phase_cond <- ggplot(lm_df, aes(x = Phase, y = NuStress, fill = specimen_type)) +
  geom_boxplot(
    width = 0.65,
    outlier.size = 0.2,
    position = position_dodge(0.75),
    color = "black"
  ) +
  scale_fill_manual(values = c("N" = "#1f77b4", "T" = "#d62728")) +
  ggpubr::stat_pvalue_manual(
    phase_stat_plot,
    label = "p.signif",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4,
    inherit.aes = FALSE
  ) +
  labs(
    title = "NuStress by phase and specimen type",
    x = "Phase",
    y = "NuStress"
  )
save_pdf(p_phase_cond, "boxplot_NuStress_phase_by_specimenType", width = 6, height = 5)

###############################################################################
## 12.3 Recompute UMAP after regressing out cell-cycle signals
###############################################################################
subset_epi_reg <- subset_epi
subset_epi_reg <- ScaleData(subset_epi_reg, vars.to.regress = c("S.Score", "G2M.Score"))
subset_epi_reg <- RunPCA(subset_epi_reg, features = VariableFeatures(subset_epi_reg))
subset_epi_reg <- RunUMAP(subset_epi_reg, dims = 1:30, reduction.name = "umap_ccreg")

p_umap_ccreg_nus <- plot_feature_umap(
  subset_epi_reg,
  feature = "NuStress",
  title = "NuStress after regressing out S.Score and G2M.Score",
  reduction_use = "umap_ccreg"
)
save_pdf(p_umap_ccreg_nus, "umap_ccreg_NuStress", width = 5.5, height = 5)

###############################################################################
## 13. Correlation between NuStress and RiboSis before/after cell-cycle regression
###############################################################################
cat("=== Step 10: NuStress-RiboSis correlation before/after cell-cycle regression ===\n")

nr_df <- subset_epi@meta.data %>%
  dplyr::select(NuStress, RiboSis, S.Score, G2M.Score, specimen_type, epi_subcluster, Phase) %>%
  dplyr::filter(
    !is.na(NuStress),
    !is.na(RiboSis),
    !is.na(S.Score),
    !is.na(G2M.Score)
  )

raw_cor_res <- suppressWarnings(cor.test(
  nr_df$NuStress, nr_df$RiboSis,
  method = "spearman"
))

fit_nus_cc  <- lm(NuStress ~ S.Score + G2M.Score, data = nr_df)
fit_ribo_cc <- lm(RiboSis  ~ S.Score + G2M.Score, data = nr_df)

nr_df$NuStress_ccresid <- resid(fit_nus_cc)
nr_df$RiboSis_ccresid  <- resid(fit_ribo_cc)

ccreg_cor_res <- suppressWarnings(cor.test(
  nr_df$NuStress_ccresid, nr_df$RiboSis_ccresid,
  method = "spearman"
))

cor_nr_summary <- data.frame(
  Analysis = c("Raw", "CellCycle_regressed"),
  rho = c(unname(raw_cor_res$estimate), unname(ccreg_cor_res$estimate)),
  p_value = c(raw_cor_res$p.value, ccreg_cor_res$p.value),
  n = c(nrow(nr_df), nrow(nr_df))
)
save_csv(cor_nr_summary, "NuStress_RiboSis_correlation_before_after_cellcycle_regression.csv")
save_csv(broom::tidy(fit_nus_cc),  "LM_NuStress_regress_cellcycle.csv")
save_csv(broom::tidy(fit_ribo_cc), "LM_RiboSis_regress_cellcycle.csv")

p_cor_nr_raw <- ggplot(nr_df, aes(x = RiboSis, y = NuStress)) +
  geom_point(size = 0.8, alpha = 0.35, color = "#6D28D9") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  labs(
    title = paste0(
      "NuStress vs RiboSis before cell-cycle regression\n",
      "Spearman rho = ", round(raw_cor_res$estimate, 3),
      ", p = ", signif(raw_cor_res$p.value, 3)
    ),
    x = "RiboSis",
    y = "NuStress"
  ) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
save_pdf(p_cor_nr_raw, "epi_scatter_NuStress_vs_RiboSis_raw", width = 5.5, height = 5)

p_cor_nr_ccreg <- ggplot(nr_df, aes(x = RiboSis_ccresid, y = NuStress_ccresid)) +
  geom_point(size = 0.8, alpha = 0.35, color = "#A855F7") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  labs(
    title = paste0(
      "NuStress vs RiboSis after regressing out S.Score and G2M.Score\n",
      "Spearman rho = ", round(ccreg_cor_res$estimate, 3),
      ", p = ", signif(ccreg_cor_res$p.value, 3)
    ),
    x = "RiboSis residual",
    y = "NuStress residual"
  ) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
save_pdf(p_cor_nr_ccreg, "epi_scatter_NuStress_vs_RiboSis_cellcycle_regressed", width = 5.8, height = 5)

cor_nr_plot_df <- cor_nr_summary %>%
  dplyr::mutate(
    Analysis = factor(Analysis, levels = c("Raw", "CellCycle_regressed"))
  )

p_cor_nr_bar <- ggplot(cor_nr_plot_df, aes(x = Analysis, y = rho, fill = Analysis)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(
    aes(label = paste0("rho = ", round(rho, 3))),
    vjust = ifelse(cor_nr_plot_df$rho >= 0, -0.3, 1.2),
    size = 4
  ) +
  scale_fill_manual(values = c("Raw" = "#8B5CF6", "CellCycle_regressed" = "#C084FC")) +
  labs(
    title = "NuStress-RiboSis correlation before and after cell-cycle regression",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )
save_pdf(p_cor_nr_bar, "epi_bar_NuStress_RiboSis_correlation_before_after_ccreg", width = 5.5, height = 5)

###############################################################################
## 14. Compare NuStress distribution before vs after cell-cycle regression
###############################################################################
cat("=== Step 11: Compare NuStress before vs after cell-cycle regression ===\n")

nus_compare_df <- subset_epi@meta.data %>%
  dplyr::select(
    NuStress, RiboSis, S.Score, G2M.Score,
    specimen_type, epi_subcluster, Phase
  ) %>%
  dplyr::filter(
    !is.na(NuStress),
    !is.na(S.Score),
    !is.na(G2M.Score),
    !is.na(epi_subcluster),
    !is.na(specimen_type)
  )

nus_compare_df$specimen_type <- factor(nus_compare_df$specimen_type, levels = c("N", "T"))
nus_compare_df$epi_subcluster <- factor(
  nus_compare_df$epi_subcluster,
  levels = unique(as.character(subset_epi$epi_subcluster))
)

fit_nus_cc_simple <- lm(NuStress ~ S.Score + G2M.Score, data = nus_compare_df)
nus_compare_df$NuStress_ccresid <- resid(fit_nus_cc_simple)
nus_compare_df$NuStress_ccadj <- nus_compare_df$NuStress_ccresid +
  mean(nus_compare_df$NuStress, na.rm = TRUE)

save_csv(
  broom::tidy(fit_nus_cc_simple),
  "LM_NuStress_regress_S_G2M_for_subcluster_comparison.csv"
)

nus_long_df <- nus_compare_df %>%
  dplyr::select(specimen_type, epi_subcluster, Phase, NuStress, NuStress_ccadj) %>%
  tidyr::pivot_longer(
    cols = c(NuStress, NuStress_ccadj),
    names_to = "Analysis",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Analysis = dplyr::recode(
      Analysis,
      "NuStress" = "Raw_NuStress",
      "NuStress_ccadj" = "CellCycle_adjusted_NuStress"
    ),
    Analysis = factor(
      Analysis,
      levels = c("Raw_NuStress", "CellCycle_adjusted_NuStress")
    )
  )

save_csv(nus_long_df, "NuStress_before_after_cellcycle_regression_longformat.csv")

p_nus_before_after_box <- ggplot(
  nus_long_df,
  aes(x = epi_subcluster, y = Score, fill = Analysis)
) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.2,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c(
    "Raw_NuStress" = "#8B5CF6",
    "CellCycle_adjusted_NuStress" = "#C084FC"
  )) +
  labs(
    title = "NuStress before and after cell-cycle regression across epithelial subclusters",
    x = "Epithelial subcluster",
    y = "NuStress"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(
  p_nus_before_after_box,
  "epi_boxplot_NuStress_before_after_cellcycle_regression",
  width = 10,
  height = 5.5
)

p_nus_before_after_by_specimen <- ggplot(
  nus_long_df,
  aes(x = epi_subcluster, y = Score, fill = Analysis)
) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.15,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~ specimen_type, nrow = 1) +
  scale_fill_manual(values = c(
    "Raw_NuStress" = "#8B5CF6",
    "CellCycle_adjusted_NuStress" = "#C084FC"
  )) +
  labs(
    title = "NuStress before and after cell-cycle regression in normal and tumor epithelial cells",
    x = "Epithelial subcluster",
    y = "NuStress"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(
  p_nus_before_after_by_specimen,
  "epi_boxplot_NuStress_before_after_cellcycle_regression_by_specimen",
  width = 11,
  height = 5.5
)

nus_mean_compare_df <- nus_long_df %>%
  group_by(epi_subcluster, Analysis) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score   = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
    .groups = "drop"
  )

p_nus_before_after_mean <- ggplot(
  nus_mean_compare_df,
  aes(x = epi_subcluster, y = mean_score, fill = Analysis)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
    position = position_dodge(width = 0.75),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    "Raw_NuStress" = "#8B5CF6",
    "CellCycle_adjusted_NuStress" = "#C084FC"
  )) +
  labs(
    title = "Mean NuStress before and after cell-cycle regression",
    x = "Epithelial subcluster",
    y = "Mean NuStress"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(
  p_nus_before_after_mean,
  "epi_bar_mean_NuStress_before_after_cellcycle_regression",
  width = 10,
  height = 5.5
)

before_after_stats <- nus_long_df %>%
  dplyr::group_by(epi_subcluster, specimen_type) %>%
  dplyr::summarise(
    p_value = tryCatch(
      wilcox.test(
        Score[Analysis == "Raw_NuStress"],
        Score[Analysis == "CellCycle_adjusted_NuStress"]
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p.signif = dplyr::case_when(
      is.na(p_adj)   ~ "NA",
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )
save_csv(
  before_after_stats,
  "NuStress_before_after_cellcycle_regression_stats_by_subcluster_and_specimen.csv"
)

###############################################################################
## 15. Tumor epithelial only
###############################################################################
cat("=== Step 12: Tumor epithelial subset ===\n")

subset_epi_tumor <- subset(subset_epi, subset = specimen_type == "T")

p_tumor_nus <- ggplot(
  subset_epi_tumor@meta.data,
  aes(x = epi_subcluster, y = NuStress, fill = epi_subcluster)
) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.2) +
  scale_fill_manual(values = epi_colors_use) +
  coord_flip() +
  labs(title = "Tumor epithelial cells: NuStress across subclusters", x = "", y = "NuStress") +
  theme(legend.position = "none")
save_pdf(p_tumor_nus, "tumor_epi_subclusters_NuStress", width = 6, height = 5.5)

p_tumor_ribo <- ggplot(
  subset_epi_tumor@meta.data,
  aes(x = epi_subcluster, y = RiboSis, fill = epi_subcluster)
) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.2) +
  scale_fill_manual(values = epi_colors_use) +
  coord_flip() +
  labs(title = "Tumor epithelial cells: RiboSis across subclusters", x = "", y = "RiboSis") +
  theme(legend.position = "none")
save_pdf(p_tumor_ribo, "tumor_epi_subclusters_RiboSis", width = 6, height = 5.5)

tumor_nus_kw  <- kruskal.test(NuStress ~ epi_subcluster, data = subset_epi_tumor@meta.data)
tumor_ribo_kw <- kruskal.test(RiboSis ~ epi_subcluster, data = subset_epi_tumor@meta.data)

tumor_kw_res <- data.frame(
  Score = c("NuStress", "RiboSis"),
  p_value = c(tumor_nus_kw$p.value, tumor_ribo_kw$p.value)
)
save_csv(tumor_kw_res, "tumor_epi_subcluster_kruskal_summary.csv")

###############################################################################
## 15.1 Tumor epithelial only: before vs after cell-cycle regression
###############################################################################
nus_long_tumor_df <- nus_long_df %>%
  dplyr::filter(specimen_type == "T")

p_tumor_nus_before_after <- ggplot(
  nus_long_tumor_df,
  aes(x = epi_subcluster, y = Score, fill = Analysis)
) +
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.2,
    color = "black",
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c(
    "Raw_NuStress" = "#8B5CF6",
    "CellCycle_adjusted_NuStress" = "#C084FC"
  )) +
  coord_flip() +
  labs(
    title = "Tumor epithelial cells: NuStress before and after cell-cycle regression",
    x = "",
    y = "NuStress"
  )
save_pdf(
  p_tumor_nus_before_after,
  "tumor_epi_subclusters_NuStress_before_after_cellcycle_regression",
  width = 7,
  height = 6
)

###############################################################################
## 16. Save objects
###############################################################################
saveRDS(seurat_obj,      file = file.path(rdsdir, "seurat_obj_with_NuStress.rds"))
saveRDS(subset_epi,      file = file.path(rdsdir, "subset_epi_with_NuStress.rds"))
saveRDS(subset_epi_reg,  file = file.path(rdsdir, "subset_epi_regressed_cellcycle.rds"))

cat("=== Analysis completed successfully ===\n")

###############################################################################
## 17. Differential expression and pathway analysis based on NuStress/RiboSis quadrants
## Purpose:
## - Define four quadrants in epithelial cells using NuStress and RiboSis quartiles.
## - Identify abnormal expression features and potential new signals across quadrants.
## - Key comparisons:
##   1. NuS_high / Ribo_high vs NuS_low / Ribo_low: overall abnormal state.
##   2. NuS_low / Ribo_high vs NuS_low / Ribo_low: elevated ribosome biogenesis program.
##   3. NuS_high / Ribo_low vs NuS_low / Ribo_low: elevated nucleolar stress program.
###############################################################################
cat("=== Step 13: Quadrant grouping and differential expression analysis ===\n")

DefaultAssay(subset_epi) <- "RNA"

###############################################################################
## 17.1 Define quartile cutoffs
###############################################################################
nus_q25  <- quantile(subset_epi$NuStress, 0.25, na.rm = TRUE)
nus_q75  <- quantile(subset_epi$NuStress, 0.75, na.rm = TRUE)
ribo_q25 <- quantile(subset_epi$RiboSis, 0.25, na.rm = TRUE)
ribo_q75 <- quantile(subset_epi$RiboSis, 0.75, na.rm = TRUE)

cutoff_df <- data.frame(
  Metric = c("NuStress_Q25", "NuStress_Q75", "RiboSis_Q25", "RiboSis_Q75"),
  Value  = c(nus_q25, nus_q75, ribo_q25, ribo_q75)
)
save_csv(cutoff_df, "Quadrant_cutoffs_NuStress_RiboSis.csv")

###############################################################################
## 17.2 Assign four quadrants
###############################################################################
subset_epi$NuS_Ribo_quadrant <- NA_character_

subset_epi$NuS_Ribo_quadrant[
  subset_epi$NuStress <= nus_q25 & subset_epi$RiboSis <= ribo_q25
] <- "NuS_low_Ribo_low"

subset_epi$NuS_Ribo_quadrant[
  subset_epi$NuStress <= nus_q25 & subset_epi$RiboSis >= ribo_q75
] <- "NuS_low_Ribo_high"

subset_epi$NuS_Ribo_quadrant[
  subset_epi$NuStress >= nus_q75 & subset_epi$RiboSis <= ribo_q25
] <- "NuS_high_Ribo_low"

subset_epi$NuS_Ribo_quadrant[
  subset_epi$NuStress >= nus_q75 & subset_epi$RiboSis >= ribo_q75
] <- "NuS_high_Ribo_high"

subset_epi$NuS_Ribo_quadrant <- factor(
  subset_epi$NuS_Ribo_quadrant,
  levels = c(
    "NuS_low_Ribo_low",
    "NuS_low_Ribo_high",
    "NuS_high_Ribo_low",
    "NuS_high_Ribo_high"
  )
)

quad_count_df <- as.data.frame(table(subset_epi$NuS_Ribo_quadrant, useNA = "ifany"))
colnames(quad_count_df) <- c("Quadrant", "Cell_count")
save_csv(quad_count_df, "Quadrant_cell_counts.csv")

cat("Quadrant distribution:\n")
print(table(subset_epi$NuS_Ribo_quadrant, useNA = "ifany"))

###############################################################################
## 17.3 Subset cells in four quadrants only
###############################################################################
cells_quad <- colnames(subset_epi)[!is.na(subset_epi$NuS_Ribo_quadrant)]
cat("Cells retained in four quadrants:", length(cells_quad), "\n")

if (length(cells_quad) == 0) {
  stop("No cells assigned to NuS/RiboSis quadrants.")
}

subset_epi_quad <- subset(subset_epi, cells = cells_quad)

subset_epi_quad$NuS_Ribo_quadrant <- factor(
  subset_epi_quad$NuS_Ribo_quadrant,
  levels = c(
    "NuS_low_Ribo_low",
    "NuS_low_Ribo_high",
    "NuS_high_Ribo_low",
    "NuS_high_Ribo_high"
  )
)

Idents(subset_epi_quad) <- "NuS_Ribo_quadrant"

print(table(subset_epi_quad$NuS_Ribo_quadrant, useNA = "ifany"))

###############################################################################
## 17.4 Basic visualization of quadrant grouping
###############################################################################
quad_cols <- c(
  "NuS_low_Ribo_low"   = "#3C5488",
  "NuS_low_Ribo_high"  = "#00A087",
  "NuS_high_Ribo_low"  = "#E64B35",
  "NuS_high_Ribo_high" = "#7E6148"
)

## (1) Scatter plot with cutoff lines
scatter_df <- subset_epi@meta.data %>%
  dplyr::select(NuStress, RiboSis, NuS_Ribo_quadrant, specimen_type, epi_subcluster) %>%
  dplyr::filter(!is.na(NuStress), !is.na(RiboSis))

p_quad_scatter <- ggplot(scatter_df, aes(x = RiboSis, y = NuStress, color = NuS_Ribo_quadrant)) +
  geom_point(size = 0.35, alpha = 0.45, na.rm = TRUE) +
  geom_vline(xintercept = c(ribo_q25, ribo_q75), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(nus_q25, nus_q75), linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = quad_cols, na.value = "grey80") +
  labs(
    title = "NuStress-RiboSis quadrants in epithelial cells",
    x = "RiboSis",
    y = "NuStress",
    color = "Quadrant"
  )
save_pdf(p_quad_scatter, "epi_scatter_NuStress_RiboSis_quadrants", width = 6.5, height = 5.8)

## (2) UMAP of quadrants
p_quad_umap <- DimPlot(
  subset_epi,
  reduction = "umap",
  group.by = "NuS_Ribo_quadrant",
  cols = quad_cols,
  pt.size = 3,
  na.value = "grey85"
) + labs(title = "UMAP of NuStress-RiboSis quadrants in epithelial cells")
save_pdf(p_quad_umap, "epi_umap_quadrants_NuStress_RiboSis", width = 7, height = 5.5)

## (3) Boxplots to confirm quadrant separation
p_box_quad_nus <- ggplot(
  subset_epi_quad@meta.data,
  aes(x = NuS_Ribo_quadrant, y = NuStress, fill = NuS_Ribo_quadrant)
) +
  geom_boxplot(width = 0.6, outlier.size = 0.2, color = "black") +
  scale_fill_manual(values = quad_cols) +
  labs(title = "NuStress across NuStress-RiboSis quadrants", x = "Quadrant", y = "NuStress") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_pdf(p_box_quad_nus, "epi_boxplot_NuStress_by_quadrant", width = 7.2, height = 5)

p_box_quad_ribo <- ggplot(
  subset_epi_quad@meta.data,
  aes(x = NuS_Ribo_quadrant, y = RiboSis, fill = NuS_Ribo_quadrant)
) +
  geom_boxplot(width = 0.6, outlier.size = 0.2, color = "black") +
  scale_fill_manual(values = quad_cols) +
  labs(title = "RiboSis across NuStress-RiboSis quadrants", x = "Quadrant", y = "RiboSis") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_pdf(p_box_quad_ribo, "epi_boxplot_RiboSis_by_quadrant", width = 7.2, height = 5)

## (4) Quadrant composition by specimen type
quad_spec_df <- subset_epi@meta.data %>%
  dplyr::filter(!is.na(NuS_Ribo_quadrant)) %>%
  dplyr::count(specimen_type, NuS_Ribo_quadrant) %>%
  dplyr::group_by(specimen_type) %>%
  dplyr::mutate(Fraction = n / sum(n))

p_quad_spec <- ggplot(quad_spec_df, aes(x = specimen_type, y = Fraction, fill = NuS_Ribo_quadrant)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = quad_cols) +
  labs(
    title = "Quadrant composition by specimen type",
    x = "Specimen type",
    y = "Fraction"
  )
save_pdf(p_quad_spec, "epi_bar_quadrant_composition_by_specimen", width = 6.2, height = 5)

## (5) Quadrant composition across epithelial subclusters
quad_subcluster_df <- subset_epi@meta.data %>%
  dplyr::filter(!is.na(NuS_Ribo_quadrant), !is.na(epi_subcluster)) %>%
  dplyr::count(epi_subcluster, NuS_Ribo_quadrant) %>%
  dplyr::group_by(epi_subcluster) %>%
  dplyr::mutate(Fraction = n / sum(n))

p_quad_subcluster <- ggplot(
  quad_subcluster_df,
  aes(x = epi_subcluster, y = Fraction, fill = NuS_Ribo_quadrant)
) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = quad_cols) +
  labs(
    title = "Quadrant composition across epithelial subclusters",
    x = "Epithelial subcluster",
    y = "Fraction"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_quad_subcluster, "epi_bar_quadrant_composition_by_subcluster", width = 10, height = 5.5)

###############################################################################
## 17.5 Differential expression analysis
###############################################################################
perform_de_analysis <- function(seurat_obj, group_col, ident_1, ident_2,
                                comparison_name,
                                logfc_threshold = 0.25,
                                min_pct = 0.1,
                                sig_fc = 0.5,
                                sig_fdr = 0.05) {
  Idents(seurat_obj) <- group_col
  
  de_genes <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    test.use = "wilcox"
  )
  
  de_genes <- de_genes %>%
    rownames_to_column("gene") %>%
    mutate(
      significant = case_when(
        p_val_adj < sig_fdr & avg_log2FC >  sig_fc ~ "Up",
        p_val_adj < sig_fdr & avg_log2FC < -sig_fc ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  de_genes_sig <- de_genes %>%
    filter(p_val_adj < sig_fdr, abs(avg_log2FC) > sig_fc)
  
  save_csv(de_genes, paste0("DEG_", comparison_name, "_all.csv"))
  save_csv(de_genes_sig, paste0("DEG_", comparison_name, "_sig.csv"))
  
  p_volcano <- ggplot(
    de_genes,
    aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300), color = significant)
  ) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = c(
      "Up" = "#E64B35",
      "Down" = "#3C5488",
      "Not Sig" = "grey75"
    )) +
    geom_hline(yintercept = -log10(sig_fdr), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-sig_fc, sig_fc), linetype = "dashed", color = "red") +
    labs(
      title = comparison_name,
      x = "avg_log2FC",
      y = "-log10(adj.P)"
    )
  
  save_pdf(p_volcano, paste0("volcano_", comparison_name), width = 7, height = 5)
  
  list(
    de_all = de_genes,
    de_sig = de_genes_sig,
    volcano = p_volcano
  )
}


cat("=== Step 14: Differential expression analyses ===\n")

de_hh_vs_ll <- perform_de_analysis(
  seurat_obj = subset_epi_quad,
  group_col = "NuS_Ribo_quadrant",
  ident_1 = "NuS_high_Ribo_high",
  ident_2 = "NuS_low_Ribo_low",
  comparison_name = "Epi_NuSHigh_RiboHigh_vs_NuSLow_RiboLow",
  logfc_threshold = 0.25,
  min_pct = 0.1,
  sig_fc = 0.5,
  sig_fdr = 0.05
)

de_lh_vs_ll <- perform_de_analysis(
  seurat_obj = subset_epi_quad,
  group_col = "NuS_Ribo_quadrant",
  ident_1 = "NuS_low_Ribo_high",
  ident_2 = "NuS_low_Ribo_low",
  comparison_name = "Epi_NuSLow_RiboHigh_vs_NuSLow_RiboLow",
  logfc_threshold = 0.25,
  min_pct = 0.1,
  sig_fc = 0.5,
  sig_fdr = 0.05
)

de_hl_vs_ll <- perform_de_analysis(
  seurat_obj = subset_epi_quad,
  group_col = "NuS_Ribo_quadrant",
  ident_1 = "NuS_high_Ribo_low",
  ident_2 = "NuS_low_Ribo_low",
  comparison_name = "Epi_NuSHigh_RiboLow_vs_NuSLow_RiboLow",
  logfc_threshold = 0.25,
  min_pct = 0.1,
  sig_fc = 0.5,
  sig_fdr = 0.05
)

deg_summary_df <- data.frame(
  Comparison = c(
    "NuS_high_Ribo_high_vs_NuS_low_Ribo_low",
    "NuS_low_Ribo_high_vs_NuS_low_Ribo_low",
    "NuS_high_Ribo_low_vs_NuS_low_Ribo_low"
  ),
  Total_sig_DEG = c(
    nrow(de_hh_vs_ll$de_sig),
    nrow(de_lh_vs_ll$de_sig),
    nrow(de_hl_vs_ll$de_sig)
  ),
  Up_in_ident1 = c(
    sum(de_hh_vs_ll$de_sig$avg_log2FC > 0, na.rm = TRUE),
    sum(de_lh_vs_ll$de_sig$avg_log2FC > 0, na.rm = TRUE),
    sum(de_hl_vs_ll$de_sig$avg_log2FC > 0, na.rm = TRUE)
  ),
  Down_in_ident1 = c(
    sum(de_hh_vs_ll$de_sig$avg_log2FC < 0, na.rm = TRUE),
    sum(de_lh_vs_ll$de_sig$avg_log2FC < 0, na.rm = TRUE),
    sum(de_hl_vs_ll$de_sig$avg_log2FC < 0, na.rm = TRUE)
  )
)
save_csv(deg_summary_df, "Quadrant_DEG_summary.csv")

###############################################################################
## 17.6 Enrichment analysis
###############################################################################
cat("=== Step 15: Enrichment analyses ===\n")
perform_enrichment_analysis <- function(sig_df, comparison_name) {
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    warning("No significant DEG available for enrichment: ", comparison_name)
    return(NULL)
  }
  
  up_genes   <- sig_df$gene[sig_df$avg_log2FC > 0]
  down_genes <- sig_df$gene[sig_df$avg_log2FC < 0]
  
  perform_go <- function(genes, suffix) {
    genes <- unique(genes)
    if (length(genes) < 5) {
      warning("Too few genes for GO enrichment: ", comparison_name, " - ", suffix)
      return(NULL)
    }
    
    ego <- tryCatch({
      enrichGO(
        gene = genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      )
    }, error = function(e) {
      warning("GO enrichment failed for ", comparison_name, " - ", suffix, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      return(NULL)
    }
    
    ego2 <- tryCatch({
      simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }, error = function(e) ego)
    
    res_df <- as.data.frame(ego2)
    save_csv(res_df, paste0("GO_", comparison_name, "_", suffix, ".csv"))
    
    p <- dotplot(ego2, showCategory = 15) +
      scale_color_viridis_c() +
      labs(title = paste0(comparison_name, " GO ", suffix))
    save_pdf(p, paste0("go_", comparison_name, "_", suffix), width = 8, height = 9)
    
    list(result = ego2, plot = p)
  }
  
  perform_kegg <- function(genes, suffix) {
    genes <- unique(genes)
    if (length(genes) < 5) {
      warning("Too few genes for KEGG enrichment: ", comparison_name, " - ", suffix)
      return(NULL)
    }
    
    entrez <- tryCatch({
      bitr(
        genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
    }, error = function(e) {
      warning("SYMBOL-to-ENTREZ conversion failed for ", comparison_name, " - ", suffix, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(entrez) || nrow(entrez) < 5) {
      return(NULL)
    }
    
    kk <- tryCatch({
      enrichKEGG(
        gene = unique(entrez$ENTREZID),
        organism = "hsa",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      )
    }, error = function(e) {
      warning("KEGG enrichment failed for ", comparison_name, " - ", suffix, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
      return(NULL)
    }
    
    res_df <- as.data.frame(kk)
    save_csv(res_df, paste0("KEGG_", comparison_name, "_", suffix, ".csv"))
    
    p <- dotplot(kk, showCategory = 15) +
      scale_color_viridis_c() +
      labs(title = paste0(comparison_name, " KEGG ", suffix))
    save_pdf(p, paste0("kegg_", comparison_name, "_", suffix), width = 8, height = 9)
    
    list(result = kk, plot = p)
  }
  
  list(
    GO_Up = perform_go(up_genes, "Up"),
    GO_Down = perform_go(down_genes, "Down"),
    KEGG_Up = perform_kegg(up_genes, "Up"),
    KEGG_Down = perform_kegg(down_genes, "Down")
  )
}
enrich_hh_vs_ll <- perform_enrichment_analysis(
  sig_df = de_hh_vs_ll$de_sig,
  comparison_name = "Epi_NuSHigh_RiboHigh_vs_NuSLow_RiboLow"
)

enrich_lh_vs_ll <- perform_enrichment_analysis(
  sig_df = de_lh_vs_ll$de_sig,
  comparison_name = "Epi_NuSLow_RiboHigh_vs_NuSLow_RiboLow"
)

enrich_hl_vs_ll <- perform_enrichment_analysis(
  sig_df = de_hl_vs_ll$de_sig,
  comparison_name = "Epi_NuSHigh_RiboLow_vs_NuSLow_RiboLow"
)

###############################################################################
## 17.7 Marker-style visualization for top genes
###############################################################################
plot_top_markers_dotplot <- function(seurat_obj, de_df, group_col, groups_use,
                                     comparison_name, top_n = 10) {
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
  
  de_up <- de_df %>%
    dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = top_n)
  
  if (nrow(de_up) == 0) return(NULL)
  
  features_use <- unique(de_up$gene)
  obj_sub <- subset(seurat_obj, subset = !!as.name(group_col) %in% groups_use)
  Idents(obj_sub) <- group_col
  
  p <- DotPlot(
    obj_sub,
    features = rev(features_use),
    group.by = group_col,
    cols = c("lightgrey", "#E64B35")
  ) +
    RotatedAxis() +
    labs(title = paste0("Top upregulated genes: ", comparison_name))
  
  save_pdf(p, paste0("dotplot_topmarkers_", comparison_name), width = 8, height = 5.5)
  return(p)
}

plot_top_markers_dotplot(
  seurat_obj = subset_epi_quad,
  de_df = de_hh_vs_ll$de_sig,
  group_col = "NuS_Ribo_quadrant",
  groups_use = c("NuS_high_Ribo_high", "NuS_low_Ribo_low"),
  comparison_name = "Epi_NuSHigh_RiboHigh_vs_NuSLow_RiboLow",
  top_n = 12
)

plot_top_markers_dotplot(
  seurat_obj = subset_epi_quad,
  de_df = de_lh_vs_ll$de_sig,
  group_col = "NuS_Ribo_quadrant",
  groups_use = c("NuS_low_Ribo_high", "NuS_low_Ribo_low"),
  comparison_name = "Epi_NuSLow_RiboHigh_vs_NuSLow_RiboLow",
  top_n = 12
)

plot_top_markers_dotplot(
  seurat_obj = subset_epi_quad,
  de_df = de_hl_vs_ll$de_sig,
  group_col = "NuS_Ribo_quadrant",
  groups_use = c("NuS_high_Ribo_low", "NuS_low_Ribo_low"),
  comparison_name = "Epi_NuSHigh_RiboLow_vs_NuSLow_RiboLow",
  top_n = 12
)

###############################################################################
## 17.8 KEGG / GO barplots for key comparisons
## Color gradient represents -log10(adjusted P).
###############################################################################
cat("=== Step 16: KEGG / GO barplots for key comparisons ===\n")

plot_enrich_bar_gradient <- function(enrich_obj, title, gradient_colors, top_n = 10,
                                     xlab = expression(-log[10]("adjusted P")),
                                     ylab = "Pathway / term") {
  if (is.null(enrich_obj) || is.null(enrich_obj$result)) {
    warning("No enrichment result available for: ", title)
    return(NULL)
  }
  
  enrich_df <- as.data.frame(enrich_obj$result)
  if (nrow(enrich_df) == 0) {
    warning("Empty enrichment result for: ", title)
    return(NULL)
  }
  
  enrich_df <- enrich_df %>%
    dplyr::arrange(p.adjust, desc(Count)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      minus_log10_FDR = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))
    )
  
  grad_fun <- colorRampPalette(gradient_colors)
  
  p <- ggplot(enrich_df, aes(x = minus_log10_FDR, y = Description, fill = minus_log10_FDR)) +
    geom_col(width = 0.72) +
    scale_fill_gradientn(
      colours = grad_fun(100),
      name = expression(-log[10]("adj.P"))
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme(
      axis.text.y = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "right"
    )
  
  return(p)
}

purple_gradient <- c("#EE9E78", "#BB4D6E", "#8F2F6D", "#622360")

## KEGG
p_kegg_lh <- plot_enrich_bar_gradient(
  enrich_obj = enrich_lh_vs_ll$KEGG_Up,
  title = "KEGG: NuS low / RiboSis high vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "KEGG pathway"
)
if (!is.null(p_kegg_lh)) {
  save_pdf(
    p_kegg_lh,
    "kegg_bar_gradient_NuSLow_RiboHigh_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

p_kegg_hl <- plot_enrich_bar_gradient(
  enrich_obj = enrich_hl_vs_ll$KEGG_Up,
  title = "KEGG: NuS high / RiboSis low vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "KEGG pathway"
)
if (!is.null(p_kegg_hl)) {
  save_pdf(
    p_kegg_hl,
    "kegg_bar_gradient_NuSHigh_RiboLow_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

p_kegg_hh <- plot_enrich_bar_gradient(
  enrich_obj = enrich_hh_vs_ll$KEGG_Up,
  title = "KEGG: NuS high / RiboSis high vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "KEGG pathway"
)
if (!is.null(p_kegg_hh)) {
  save_pdf(
    p_kegg_hh,
    "kegg_bar_gradient_NuSHigh_RiboHigh_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

## GO
p_go_lh <- plot_enrich_bar_gradient(
  enrich_obj = enrich_lh_vs_ll$GO_Up,
  title = "GO: NuS low / RiboSis high vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "GO biological process"
)
if (!is.null(p_go_lh)) {
  save_pdf(
    p_go_lh,
    "go_bar_gradient_NuSLow_RiboHigh_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

p_go_hl <- plot_enrich_bar_gradient(
  enrich_obj = enrich_hl_vs_ll$GO_Up,
  title = "GO: NuS high / RiboSis low vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "GO biological process"
)
if (!is.null(p_go_hl)) {
  save_pdf(
    p_go_hl,
    "go_bar_gradient_NuSHigh_RiboLow_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

p_go_hh <- plot_enrich_bar_gradient(
  enrich_obj = enrich_hh_vs_ll$GO_Up,
  title = "GO: NuS high / RiboSis high vs NuS low / RiboSis low",
  gradient_colors = purple_gradient,
  top_n = 10,
  ylab = "GO biological process"
)
if (!is.null(p_go_hh)) {
  save_pdf(
    p_go_hh,
    "go_bar_gradient_NuSHigh_RiboHigh_vs_NuSLow_RiboLow",
    width = 8.2,
    height = 5.4
  )
}

###############################################################################
## 17.9 Save quadrant objects
###############################################################################
saveRDS(subset_epi_quad, file = file.path(rdsdir, "subset_epi_quadrant_grouped.rds"))

cat("=== Quadrant-based DEG and pathway analysis completed successfully ===\n")

# Step5_pre/post-treatment and primary/metastasis comparison + untreated DE/GO (final version)
# ======================================================================
# 1. Create environment
# ======================================================================
getwd()
setwd("/home/xxm_xxm/CJX_workspace/SC-ST/")
rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 1e6)

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ======================================================================
# 2. Parameter settings
# ======================================================================
colon_samples <- c("ST-colon1", "ST-colon2", "ST-colon3", "ST-colon4")
liver_samples <- c("ST-liver1", "ST-liver2", "ST-liver3", "ST-liver4")

untreated_colon_samples <- c("ST-colon1", "ST-colon2")
treated_colon_samples   <- c("ST-colon3", "ST-colon4")

untreated_liver_samples <- c("ST-liver1", "ST-liver2")
treated_liver_samples   <- c("ST-liver3", "ST-liver4")

# Unified color palette
treatment_colors <- c(
  "Untreated" = "#EE9E78",
  "Treated" = "#8F2F6D"
)

tissue_colors <- c(
  "Colon Cancer" = "#EE9E78",
  "Liver Metastasis" = "#8F2F6D"
)

# ======================================================================
# 3. Utility functions
# ======================================================================

# ------------------------------
# 3.1 Safely read data frames from *_data.rda files
# Compatible with previously saved files where the object name is data4
# ------------------------------
load_score_data <- function(sample_name, tissue_type, treatment) {
  f <- paste0(sample_name, "_data.rda")
  if (!file.exists(f)) stop("File not found: ", f)

  env_tmp <- new.env()
  load(f, envir = env_tmp)

  obj_names <- ls(env_tmp)
  if (length(obj_names) == 0) stop("No object found in file: ", f)

  # Prefer data4; otherwise use the first data.frame
  if ("data4" %in% obj_names) {
    df <- get("data4", envir = env_tmp)
  } else {
    candidate <- obj_names[sapply(obj_names, function(x) is.data.frame(get(x, envir = env_tmp)))]
    if (length(candidate) == 0) stop("No data.frame found in file: ", f)
    df <- get(candidate[1], envir = env_tmp)
  }

  df <- as.data.frame(df)
  df$sample_id <- sample_name
  df$tissue_type <- tissue_type
  df$treatment <- treatment
  df$group <- paste(gsub(" ", "_", tissue_type), treatment, sep = "_")

  return(df)
}

# ------------------------------
# 3.2 Safely read Seurat objects
# ------------------------------
load_seurat_object <- function(sample_name, object_name) {
  f <- paste0(sample_name, ".rda")
  if (!file.exists(f)) stop("File not found: ", f)

  env_tmp <- new.env()
  load(f, envir = env_tmp)

  if (!exists(object_name, envir = env_tmp)) {
    stop("File ", f, " does not contain object ", object_name)
  }

  obj <- get(object_name, envir = env_tmp)
  return(obj)
}

# ------------------------------
# 3.3 Add grouping information to Seurat objects
# ------------------------------
annotate_seurat_obj <- function(obj, sample_name, tissue_type, treatment) {
  obj$sample_id <- sample_name
  obj$tissue_type <- tissue_type
  obj$treatment <- treatment
  obj$group <- paste(gsub(" ", "_", tissue_type), treatment, sep = "_")
  return(obj)
}

# ------------------------------
# 3.4 Boxplot/violin plot function for continuous variables
# ------------------------------
plot_violin_box <- function(df, xvar, yvar, fillvar,
                            facetvar = NULL,
                            colors = NULL,
                            title = NULL,
                            ylab = NULL,
                            xlab = NULL,
                            compare_method = "wilcox.test",
                            save_file = NULL,
                            width = 5,
                            height = 4) {

  p <- ggplot(df, aes_string(x = xvar, y = yvar, fill = fillvar)) +
    geom_violin(alpha = 0.7, trim = TRUE, color = NA) +
    geom_boxplot(width = 0.2, alpha = 0.85, outlier.shape = NA,
                 color = "black", linewidth = 0.5)

  if (!is.null(facetvar)) {
    p <- p + facet_wrap(as.formula(paste("~", facetvar)), scales = "free_x")
  }

  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }

  p <- p +
    stat_compare_means(
      method = compare_method,
      label = "p.signif"
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab,
      fill = fillvar
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(fill = "lightgray", color = "black"),
      strip.text = element_text(face = "bold", size = 11)
    )

  if (!is.null(save_file)) {
    ggsave(save_file, plot = p, width = width, height = height)
  }

  return(p)
}

# ------------------------------
# 3.5 Save statistical summary tables
# ------------------------------
save_wilcox_summary <- function(df, group_var, score_var, split_var = NULL, outfile) {
  if (is.null(split_var)) {
    res <- df %>%
      wilcox_test(as.formula(paste(score_var, "~", group_var))) %>%
      add_significance()
  } else {
    res <- df %>%
      group_by(.data[[split_var]]) %>%
      wilcox_test(as.formula(paste(score_var, "~", group_var))) %>%
      add_significance()
  }
  write.csv(res, outfile, row.names = FALSE)
  return(res)
}

# ======================================================================
# 4. Load all score data
# ======================================================================
score_list <- list(
  load_score_data("ST-colon1", "Colon Cancer", "Untreated"),
  load_score_data("ST-colon2", "Colon Cancer", "Untreated"),
  load_score_data("ST-colon3", "Colon Cancer", "Treated"),
  load_score_data("ST-colon4", "Colon Cancer", "Treated"),
  load_score_data("ST-liver1", "Liver Metastasis", "Untreated"),
  load_score_data("ST-liver2", "Liver Metastasis", "Untreated"),
  load_score_data("ST-liver3", "Liver Metastasis", "Treated"),
  load_score_data("ST-liver4", "Liver Metastasis", "Treated")
)

combined_data <- bind_rows(score_list)

cat("Basic information for combined_data:\n")
print(table(combined_data$tissue_type, combined_data$treatment))
print(table(combined_data$name))

# Tumor cells only
tumor_data <- combined_data %>%
  filter(name == "Tumor")

cat("Basic information for tumor_data:\n")
print(table(tumor_data$tissue_type, tumor_data$treatment))

# Pretreatment data
untreated_tumor_data <- tumor_data %>%
  filter(treatment == "Untreated")

# ======================================================================
# 5. Global visualization: all cell types
# ======================================================================
p_ribo_all <- ggplot(combined_data, aes(x = name, y = Ribosis, fill = name)) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA,
               color = "black", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "RiboSis Score Distribution by Cell Type",
    x = "Cell Type",
    y = "RiboSis Score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )
ggsave("RiboSis_all_celltypes.pdf", p_ribo_all, width = 6, height = 4)

p_nus_all <- ggplot(combined_data, aes(x = name, y = NuS, fill = name)) +
  geom_violin(alpha = 0.8, scale = "width", trim = TRUE, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA,
               color = "black", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "NuS Score Distribution by Cell Type",
    x = "Cell Type",
    y = "NuS Score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )
ggsave("NuS_all_celltypes.pdf", p_nus_all, width = 6, height = 4)



# ======================================================================
# 6. Analysis 1: NuS / RiboSis changes in tumor cells before and after chemotherapy
# ======================================================================
# Keep this section: compare Untreated vs Treated within each tissue_type

p_treat_ribo <- plot_violin_box(
  df = tumor_data,
  xvar = "treatment",
  yvar = "Ribosis",
  fillvar = "treatment",
  facetvar = "tissue_type",
  colors = treatment_colors,
  title = "RiboSis Score: Treatment Effect in Tumor Cells",
  xlab = "Treatment",
  ylab = "RiboSis Score",
  save_file = "Tumor_RiboSis_treatment_by_tissue.pdf",
  width = 5,
  height = 4
)

p_treat_nus <- plot_violin_box(
  df = tumor_data,
  xvar = "treatment",
  yvar = "NuS",
  fillvar = "treatment",
  facetvar = "tissue_type",
  colors = treatment_colors,
  title = "NuS Score: Treatment Effect in Tumor Cells",
  xlab = "Treatment",
  ylab = "NuS Score",
  save_file = "Tumor_NuS_treatment_by_tissue.pdf",
  width = 5,
  height = 4
)

treat_ribo_stats <- save_wilcox_summary(
  df = tumor_data,
  group_var = "treatment",
  score_var = "Ribosis",
  split_var = "tissue_type",
  outfile = "Tumor_RiboSis_treatment_stats.csv"
)

treat_nus_stats <- save_wilcox_summary(
  df = tumor_data,
  group_var = "treatment",
  score_var = "NuS",
  split_var = "tissue_type",
  outfile = "Tumor_NuS_treatment_stats.csv"
)

# ======================================================================
# 7. Analysis 2: NuS / RiboSis changes in primary vs metastatic tumor cells
# ======================================================================
# Change: add separate comparisons before and after treatment
# Specifically:
#   7.1 untreated: Colon Cancer vs Liver Metastasis
#   7.2 treated:   Colon Cancer vs Liver Metastasis

untreated_tumor_data <- tumor_data %>%
  filter(treatment == "Untreated")

treated_tumor_data <- tumor_data %>%
  filter(treatment == "Treated")

# 7.1 untreated
p_untreated_site_nus <- plot_violin_box(
  df = untreated_tumor_data,
  xvar = "tissue_type",
  yvar = "NuS",
  fillvar = "tissue_type",
  colors = tissue_colors,
  title = "Untreated Tumor Cells: NuS in Primary vs Metastatic Lesions",
  xlab = "Tissue Type",
  ylab = "NuS Score",
  save_file = "Untreated_Tumor_NuS_primary_vs_metastasis.pdf",
  width = 4,
  height = 4
)

p_untreated_site_ribo <- plot_violin_box(
  df = untreated_tumor_data,
  xvar = "tissue_type",
  yvar = "Ribosis",
  fillvar = "tissue_type",
  colors = tissue_colors,
  title = "Untreated Tumor Cells: RiboSis in Primary vs Metastatic Lesions",
  xlab = "Tissue Type",
  ylab = "RiboSis Score",
  save_file = "Untreated_Tumor_RiboSis_primary_vs_metastasis.pdf",
  width = 4,
  height = 4
)

untreated_site_nus_stats <- save_wilcox_summary(
  df = untreated_tumor_data,
  group_var = "tissue_type",
  score_var = "NuS",
  outfile = "Untreated_Tumor_NuS_primary_vs_metastasis_stats.csv"
)

untreated_site_ribo_stats <- save_wilcox_summary(
  df = untreated_tumor_data,
  group_var = "tissue_type",
  score_var = "Ribosis",
  outfile = "Untreated_Tumor_RiboSis_primary_vs_metastasis_stats.csv"
)

# 7.2 treated
p_treated_site_nus <- plot_violin_box(
  df = treated_tumor_data,
  xvar = "tissue_type",
  yvar = "NuS",
  fillvar = "tissue_type",
  colors = tissue_colors,
  title = "Treated Tumor Cells: NuS in Primary vs Metastatic Lesions",
  xlab = "Tissue Type",
  ylab = "NuS Score",
  save_file = "Treated_Tumor_NuS_primary_vs_metastasis.pdf",
  width = 4,
  height = 4
)

p_treated_site_ribo <- plot_violin_box(
  df = treated_tumor_data,
  xvar = "tissue_type",
  yvar = "Ribosis",
  fillvar = "tissue_type",
  colors = tissue_colors,
  title = "Treated Tumor Cells: RiboSis in Primary vs Metastatic Lesions",
  xlab = "Tissue Type",
  ylab = "RiboSis Score",
  save_file = "Treated_Tumor_RiboSis_primary_vs_metastasis.pdf",
  width = 4,
  height = 4
)

treated_site_nus_stats <- save_wilcox_summary(
  df = treated_tumor_data,
  group_var = "tissue_type",
  score_var = "NuS",
  outfile = "Treated_Tumor_NuS_primary_vs_metastasis_stats.csv"
)

treated_site_ribo_stats <- save_wilcox_summary(
  df = treated_tumor_data,
  group_var = "tissue_type",
  score_var = "Ribosis",
  outfile = "Treated_Tumor_RiboSis_primary_vs_metastasis_stats.csv"
)

# ======================================================================
# 8. Analysis 3: Differential genes in primary vs metastatic tumor cells
# ======================================================================
# Change: perform DE + GO separately for untreated and treated samples

# ------------------------------
# 8.1 Read all complete Seurat objects
# ------------------------------
obj_colon1 <- load_seurat_object("ST-colon1", "obj_colon")
obj_colon2 <- load_seurat_object("ST-colon2", "obj_colon")
obj_colon3 <- load_seurat_object("ST-colon3", "obj_colon")
obj_colon4 <- load_seurat_object("ST-colon4", "obj_colon")

obj_liver1 <- load_seurat_object("ST-liver1", "obj_liver")
obj_liver2 <- load_seurat_object("ST-liver2", "obj_liver")
obj_liver3 <- load_seurat_object("ST-liver3", "obj_liver")
obj_liver4 <- load_seurat_object("ST-liver4", "obj_liver")

obj_colon1 <- annotate_seurat_obj(obj_colon1, "ST-colon1", "Colon Cancer", "Untreated")
obj_colon2 <- annotate_seurat_obj(obj_colon2, "ST-colon2", "Colon Cancer", "Untreated")
obj_colon3 <- annotate_seurat_obj(obj_colon3, "ST-colon3", "Colon Cancer", "Treated")
obj_colon4 <- annotate_seurat_obj(obj_colon4, "ST-colon4", "Colon Cancer", "Treated")

obj_liver1 <- annotate_seurat_obj(obj_liver1, "ST-liver1", "Liver Metastasis", "Untreated")
obj_liver2 <- annotate_seurat_obj(obj_liver2, "ST-liver2", "Liver Metastasis", "Untreated")
obj_liver3 <- annotate_seurat_obj(obj_liver3, "ST-liver3", "Liver Metastasis", "Treated")
obj_liver4 <- annotate_seurat_obj(obj_liver4, "ST-liver4", "Liver Metastasis", "Treated")

# ------------------------------
# 8.2 Build untreated tumor object
# ------------------------------
untreated_colon_tumor_list <- list(
  subset(obj_colon1, subset = manual.annot == "Tumor"),
  subset(obj_colon2, subset = manual.annot == "Tumor")
)

untreated_liver_tumor_list <- list(
  subset(obj_liver1, subset = manual.annot == "Tumor"),
  subset(obj_liver2, subset = manual.annot == "Tumor")
)

untreated_tumor_obj <- merge(
  x = untreated_colon_tumor_list[[1]],
  y = c(
    untreated_colon_tumor_list[[2]],
    untreated_liver_tumor_list[[1]],
    untreated_liver_tumor_list[[2]]
  ),
  add.cell.ids = c("colon1", "colon2", "liver1", "liver2"),
  project = "Untreated_Primary_vs_Metastasis_Tumor"
)

Idents(untreated_tumor_obj) <- "tissue_type"

cat("Identity groups for untreated_tumor_obj:\n")
print(table(Idents(untreated_tumor_obj)))

de_untreated_primary_vs_met <- FindMarkers(
  untreated_tumor_obj,
  ident.1 = "Liver Metastasis",
  ident.2 = "Colon Cancer",
  assay = DefaultAssay(untreated_tumor_obj),
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox"
)

de_untreated_primary_vs_met$gene <- rownames(de_untreated_primary_vs_met)
de_untreated_primary_vs_met <- de_untreated_primary_vs_met %>%
  arrange(p_val_adj, desc(avg_log2FC))

write.csv(
  de_untreated_primary_vs_met,
  "Untreated_Tumor_LiverMet_vs_Primary_DEG.csv",
  row.names = FALSE
)

deg_sig_untreated <- de_untreated_primary_vs_met %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

write.csv(
  deg_sig_untreated,
  "Untreated_Tumor_LiverMet_vs_Primary_DEG_sig.csv",
  row.names = FALSE
)

deg_up_untreated <- deg_sig_untreated %>%
  filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  unique()

deg_down_untreated <- deg_sig_untreated %>%
  filter(avg_log2FC < 0) %>%
  pull(gene) %>%
  unique()

write.csv(
  data.frame(gene = deg_up_untreated),
  "Untreated_Tumor_LiverMet_vs_Primary_up_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = deg_down_untreated),
  "Untreated_Tumor_LiverMet_vs_Primary_down_genes.csv",
  row.names = FALSE
)

# ------------------------------
# 8.3 Build treated tumor object
# ------------------------------
treated_colon_tumor_list <- list(
  subset(obj_colon3, subset = manual.annot == "Tumor"),
  subset(obj_colon4, subset = manual.annot == "Tumor")
)

treated_liver_tumor_list <- list(
  subset(obj_liver3, subset = manual.annot == "Tumor"),
  subset(obj_liver4, subset = manual.annot == "Tumor")
)

treated_tumor_obj <- merge(
  x = treated_colon_tumor_list[[1]],
  y = c(
    treated_colon_tumor_list[[2]],
    treated_liver_tumor_list[[1]],
    treated_liver_tumor_list[[2]]
  ),
  add.cell.ids = c("colon3", "colon4", "liver3", "liver4"),
  project = "Treated_Primary_vs_Metastasis_Tumor"
)

Idents(treated_tumor_obj) <- "tissue_type"

cat("Identity groups for treated_tumor_obj:\n")
print(table(Idents(treated_tumor_obj)))

de_treated_primary_vs_met <- FindMarkers(
  treated_tumor_obj,
  ident.1 = "Liver Metastasis",
  ident.2 = "Colon Cancer",
  assay = DefaultAssay(treated_tumor_obj),
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox"
)

de_treated_primary_vs_met$gene <- rownames(de_treated_primary_vs_met)
de_treated_primary_vs_met <- de_treated_primary_vs_met %>%
  arrange(p_val_adj, desc(avg_log2FC))

write.csv(
  de_treated_primary_vs_met,
  "Treated_Tumor_LiverMet_vs_Primary_DEG.csv",
  row.names = FALSE
)

deg_sig_treated <- de_treated_primary_vs_met %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

write.csv(
  deg_sig_treated,
  "Treated_Tumor_LiverMet_vs_Primary_DEG_sig.csv",
  row.names = FALSE
)

deg_up_treated <- deg_sig_treated %>%
  filter(avg_log2FC > 0) %>%
  pull(gene) %>%
  unique()

deg_down_treated <- deg_sig_treated %>%
  filter(avg_log2FC < 0) %>%
  pull(gene) %>%
  unique()

write.csv(
  data.frame(gene = deg_up_treated),
  "Treated_Tumor_LiverMet_vs_Primary_up_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = deg_down_treated),
  "Treated_Tumor_LiverMet_vs_Primary_down_genes.csv",
  row.names = FALSE
)

# ======================================================================
# 9. GO enrichment analysis
# ======================================================================
run_go_bp <- function(gene_vec, out_csv) {
  if (length(gene_vec) == 0) return(NULL)

  entrez_df <- bitr(
    gene_vec,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (is.null(entrez_df) || nrow(entrez_df) == 0) return(NULL)

  ego <- enrichGO(
    gene = unique(entrez_df$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), out_csv, row.names = FALSE)
  }

  return(ego)
}

# untreated GO
ego_up_untreated <- run_go_bp(
  deg_up_untreated,
  "Untreated_Tumor_LiverMet_vs_Primary_GO_BP_up.csv"
)

ego_down_untreated <- run_go_bp(
  deg_down_untreated,
  "Untreated_Tumor_LiverMet_vs_Primary_GO_BP_down.csv"
)

if (!is.null(ego_up_untreated) && nrow(as.data.frame(ego_up_untreated)) > 0) {
  p_go_up_untreated <- dotplot(ego_up_untreated, showCategory = 15) +
    ggtitle("GO BP enrichment: Up in Liver Metastasis vs Primary (Untreated Tumor)")
  ggsave(
    "Untreated_Tumor_LiverMet_vs_Primary_GO_BP_up_dotplot.pdf",
    p_go_up_untreated,
    width = 8,
    height = 6
  )
}

if (!is.null(ego_down_untreated) && nrow(as.data.frame(ego_down_untreated)) > 0) {
  p_go_down_untreated <- dotplot(ego_down_untreated, showCategory = 15) +
    ggtitle("GO BP enrichment: Down in Liver Metastasis vs Primary (Untreated Tumor)")
  ggsave(
    "Untreated_Tumor_LiverMet_vs_Primary_GO_BP_down_dotplot.pdf",
    p_go_down_untreated,
    width = 8,
    height = 6
  )
}

# treated GO
ego_up_treated <- run_go_bp(
  deg_up_treated,
  "Treated_Tumor_LiverMet_vs_Primary_GO_BP_up.csv"
)

ego_down_treated <- run_go_bp(
  deg_down_treated,
  "Treated_Tumor_LiverMet_vs_Primary_GO_BP_down.csv"
)

if (!is.null(ego_up_treated) && nrow(as.data.frame(ego_up_treated)) > 0) {
  p_go_up_treated <- dotplot(ego_up_treated, showCategory = 15) +
    ggtitle("GO BP enrichment: Up in Liver Metastasis vs Primary (Treated Tumor)")
  ggsave(
    "Treated_Tumor_LiverMet_vs_Primary_GO_BP_up_dotplot.pdf",
    p_go_up_treated,
    width = 8,
    height = 6
  )
}

if (!is.null(ego_down_treated) && nrow(as.data.frame(ego_down_treated)) > 0) {
  p_go_down_treated <- dotplot(ego_down_treated, showCategory = 15) +
    ggtitle("GO BP enrichment: Down in Liver Metastasis vs Primary (Treated Tumor)")
  ggsave(
    "Treated_Tumor_LiverMet_vs_Primary_GO_BP_down_dotplot.pdf",
    p_go_down_treated,
    width = 8,
    height = 6
  )
}

# ======================================================================
# 10. Save summary information
# ======================================================================
summary_info <- list(
  combined_data_dim = dim(combined_data),
  tumor_data_dim = dim(tumor_data),
  untreated_tumor_data_dim = dim(untreated_tumor_data),
  treated_tumor_data_dim = dim(treated_tumor_data),
  untreated_tumor_obj_cells = ncol(untreated_tumor_obj),
  treated_tumor_obj_cells = ncol(treated_tumor_obj),
  deg_total_untreated = nrow(de_untreated_primary_vs_met),
  deg_sig_untreated = nrow(deg_sig_untreated),
  deg_up_untreated = length(deg_up_untreated),
  deg_down_untreated = length(deg_down_untreated),
  deg_total_treated = nrow(de_treated_primary_vs_met),
  deg_sig_treated = nrow(deg_sig_treated),
  deg_up_treated = length(deg_up_treated),
  deg_down_treated = length(deg_down_treated)
)

save(summary_info, file = "Step5_summary_info.rda")

cat("\n=============================\n")
cat("Step5 analysis completed\n")
cat("=============================\n")
print(summary_info)

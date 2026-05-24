rm(list = ls())
options(stringsAsFactors = FALSE)

#===============================
# 0. Initialization
#===============================
setwd("/home/xxm_xxm/CJX_workspace/geneset")

library(Seurat)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(readr)

# Path settings: keep the style consistent with Step 2 / Step 3
input_dir <- file.path(getwd(), "final_nucleolar_gene_sets")
out_dir <- file.path(getwd(), "step4_reference_overlap")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

overlap_gene_dir <- file.path(out_dir, "overlap_gene_lists")
if (!dir.exists(overlap_gene_dir)) dir.create(overlap_gene_dir, recursive = TRUE)

cat("Input directory:", input_dir, "\n")
cat("Output directory:", out_dir, "\n\n")

#===============================
# 1. Read NuStress / NuStressCore gene sets
#===============================
read_gene_symbols <- function(file_path) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  x <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("gene_symbol" %in% colnames(x))) {
    stop(paste("Missing 'gene_symbol' column in:", file_path))
  }
  
  g <- unique(trimws(toupper(x$gene_symbol)))
  g <- g[!is.na(g) & g != ""]
  return(sort(g))
}

NuStress_UP <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_up_genes.csv"))
NuStress_DOWN <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_down_genes.csv"))
NuStressCore_UP <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_core_up_genes.csv"))
NuStressCore_DOWN <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_core_down_genes.csv"))

NuStress_ALL <- sort(unique(c(NuStress_UP, NuStress_DOWN)))
NuStressCore_ALL <- sort(unique(c(NuStressCore_UP, NuStressCore_DOWN)))

cat("NuStress UP:", length(NuStress_UP), "\n")
cat("NuStress DOWN:", length(NuStress_DOWN), "\n")
cat("NuStress ALL:", length(NuStress_ALL), "\n\n")

cat("NuStressCore UP:", length(NuStressCore_UP), "\n")
cat("NuStressCore DOWN:", length(NuStressCore_DOWN), "\n")
cat("NuStressCore ALL:", length(NuStressCore_ALL), "\n\n")

#===============================
# 2. Retrieve reference gene sets
#===============================
## 2.1 Seurat cell cycle genes
cc <- Seurat::cc.genes.updated.2019
Seurat_S_genes <- sort(unique(toupper(trimws(cc$s.genes))))
Seurat_G2M_genes <- sort(unique(toupper(trimws(cc$g2m.genes))))

## 2.2 MSigDB gene sets
msig_all <- msigdbr(species = "Homo sapiens")

get_msig_set <- function(msig_df, set_name) {
  x <- msig_df %>%
    filter(gs_name == set_name) %>%
    pull(gene_symbol) %>%
    unique() %>%
    toupper() %>%
    trimws()
  
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    warning(paste("Gene set not found or empty:", set_name))
  }
  return(sort(unique(x)))
}

reference_sets <- list(
  Seurat_S_genes = Seurat_S_genes,
  Seurat_G2M_genes = Seurat_G2M_genes,
  KEGG_CELL_CYCLE = get_msig_set(msig_all, "KEGG_CELL_CYCLE"),
  HALLMARK_G2M_CHECKPOINT = get_msig_set(msig_all, "HALLMARK_G2M_CHECKPOINT"),
  HALLMARK_E2F_TARGETS = get_msig_set(msig_all, "HALLMARK_E2F_TARGETS"),
  GO_DNA_DAMAGE_RESPONSE = get_msig_set(msig_all, "GO_DNA_DAMAGE_RESPONSE"),
  HALLMARK_UNFOLDED_PROTEIN_RESPONSE = get_msig_set(msig_all, "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
  HALLMARK_P53_PATHWAY = get_msig_set(msig_all, "HALLMARK_P53_PATHWAY")
)

reference_meta <- data.frame(
  reference_set = names(reference_sets),
  category = c(
    "Cell cycle",
    "Cell cycle",
    "Cell cycle",
    "Cell cycle",
    "Cell cycle",
    "DNA damage",
    "ER stress",
    "p53 signaling"
  ),
  stringsAsFactors = FALSE
)

# Save reference gene sets for direct reuse in downstream analyses
for (nm in names(reference_sets)) {
  write.csv(
    data.frame(gene_symbol = reference_sets[[nm]], stringsAsFactors = FALSE),
    file.path(out_dir, paste0(nm, ".csv")),
    row.names = FALSE
  )
}

#===============================
# 3. Overlap calculation functions
#===============================
calc_overlap_one <- function(query_genes, query_name, ref_genes, ref_name) {
  overlap_genes <- sort(intersect(query_genes, ref_genes))
  
  data.frame(
    query_set = query_name,
    reference_set = ref_name,
    query_size = length(query_genes),
    reference_size = length(ref_genes),
    overlap_n = length(overlap_genes),
    overlap_ratio_in_query = length(overlap_genes) / length(query_genes),
    overlap_percent_in_query = length(overlap_genes) / length(query_genes) * 100,
    overlap_ratio_in_reference = ifelse(length(ref_genes) > 0,
                                        length(overlap_genes) / length(ref_genes),
                                        NA),
    overlap_percent_in_reference = ifelse(length(ref_genes) > 0,
                                          length(overlap_genes) / length(ref_genes) * 100,
                                          NA),
    overlap_genes = paste(overlap_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

calc_overlap_all <- function(query_genes, query_name, reference_sets) {
  do.call(rbind, lapply(names(reference_sets), function(nm) {
    calc_overlap_one(
      query_genes = query_genes,
      query_name = query_name,
      ref_genes = reference_sets[[nm]],
      ref_name = nm
    )
  }))
}

#===============================
# 4. Calculate overlaps between NuStress / NuStressCore and reference gene sets
#===============================
overlap_nustress <- calc_overlap_all(NuStress_ALL, "NuStress_ALL", reference_sets)
overlap_nustresscore <- calc_overlap_all(NuStressCore_ALL, "NuStressCore_ALL", reference_sets)

overlap_summary <- rbind(overlap_nustress, overlap_nustresscore)
overlap_summary <- merge(overlap_summary, reference_meta, by = "reference_set", all.x = TRUE)

overlap_summary <- overlap_summary[, c(
  "query_set", "category", "reference_set",
  "query_size", "reference_size",
  "overlap_n",
  "overlap_ratio_in_query", "overlap_percent_in_query",
  "overlap_ratio_in_reference", "overlap_percent_in_reference",
  "overlap_genes"
)]

overlap_summary <- overlap_summary[order(
  overlap_summary$query_set,
  factor(overlap_summary$category, levels = c("Cell cycle", "DNA damage", "ER stress", "p53 signaling")),
  overlap_summary$reference_set
), ]

write.csv(
  overlap_summary,
  file.path(out_dir, "NuStress_vs_reference_sets_overlap_summary.csv"),
  row.names = FALSE
)

#===============================
# 5. Export overlap gene lists for each comparison
#===============================
for (qname in c("NuStress_ALL", "NuStressCore_ALL")) {
  qgenes <- if (qname == "NuStress_ALL") NuStress_ALL else NuStressCore_ALL
  
  for (rname in names(reference_sets)) {
    overlap_genes <- sort(intersect(qgenes, reference_sets[[rname]]))
    
    out_df <- data.frame(
      gene_symbol = overlap_genes,
      stringsAsFactors = FALSE
    )
    
    write.csv(
      out_df,
      file.path(overlap_gene_dir, paste0(qname, "_overlap_", rname, ".csv")),
      row.names = FALSE
    )
  }
}

#===============================
# 6. Visualization
# Recommended for manuscript figures: horizontal grouped barplot
# Display overlap_percent_in_query
#===============================
my_cols <- c("#C4E7C1", "#93D4BC", "#51B3D1", "#08589E")

plot_df <- overlap_summary
plot_df$query_set <- factor(plot_df$query_set, levels = c("NuStress_ALL", "NuStressCore_ALL"))
plot_df$reference_set <- factor(
  plot_df$reference_set,
  levels = rev(c(
    "Seurat_S_genes",
    "Seurat_G2M_genes",
    "KEGG_CELL_CYCLE",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS",
    "GO_DNA_DAMAGE_RESPONSE",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_P53_PATHWAY"
  ))
)

## 6.1 All reference gene sets
p1 <- ggplot(plot_df, aes(x = overlap_percent_in_query, y = reference_set, fill = query_set)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", overlap_percent_in_query)),
    position = position_dodge(width = 0.75),
    hjust = -0.15,
    size = 3
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(my_cols[3], my_cols[4])) +
  labs(
    title = "Overlap of nucleolar stress gene sets with reference pathways",
    x = "Overlapping genes as % of NuStress / NuStressCore",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(
    xlim = c(0, max(plot_df$overlap_percent_in_query, na.rm = TRUE) * 1.22)
  )

pdf(file.path(out_dir, "Figure_step4_overlap_percent_in_query.pdf"),
    width = 8.5, height = 6.5)
print(p1)
dev.off()

## 6.2 Cell cycle-related gene sets only
plot_df_cell <- subset(plot_df, category == "Cell cycle")

p2 <- ggplot(plot_df_cell, aes(x = overlap_percent_in_query, y = reference_set, fill = query_set)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", overlap_percent_in_query)),
    position = position_dodge(width = 0.75),
    hjust = -0.15,
    size = 3
  ) +
  scale_fill_manual(values = c(my_cols[3], my_cols[4])) +
  labs(
    title = "Overlap of nucleolar stress gene sets with cell cycle-related gene sets",
    x = "Overlapping genes as % of NuStress / NuStressCore",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  coord_cartesian(
    xlim = c(0, max(plot_df_cell$overlap_percent_in_query, na.rm = TRUE) * 1.25)
  )

pdf(file.path(out_dir, "Figure_step4_cell_cycle_overlap_percent_in_query.pdf"),
    width = 8, height = 4.8)
print(p2)
dev.off()

#===============================
# 7. Save intermediate Step 4 objects
# This allows direct reuse in downstream analyses
#===============================
save(
  NuStress_UP, NuStress_DOWN, NuStress_ALL,
  NuStressCore_UP, NuStressCore_DOWN, NuStressCore_ALL,
  reference_sets, reference_meta,
  overlap_nustress, overlap_nustresscore, overlap_summary,
  file = file.path(out_dir, "step4_reference_overlap_objects.RData")
)

#===============================
# 8. Summary statistics
#===============================
summary_df <- overlap_summary[, c(
  "query_set", "category", "reference_set",
  "query_size", "reference_size",
  "overlap_n", "overlap_percent_in_query", "overlap_percent_in_reference"
)]

write.csv(
  summary_df,
  file.path(out_dir, "step4_summary_statistics.csv"),
  row.names = FALSE
)

#===============================
# 9. Console output
#===============================
cat("Step 4 finished.\n\n")

cat("Main result table:\n")
cat(file.path(out_dir, "NuStress_vs_reference_sets_overlap_summary.csv"), "\n\n")

cat("Figures:\n")
cat(file.path(out_dir, "Figure_step4_overlap_percent_in_query.pdf"), "\n")
cat(file.path(out_dir, "Figure_step4_cell_cycle_overlap_percent_in_query.pdf"), "\n\n")

cat("Saved RData:\n")
cat(file.path(out_dir, "step4_reference_overlap_objects.RData"), "\n\n")

print(summary_df)
###############################################################################
## Task 1: Pan-cancer analysis of RiboSis and NuS ######## Step7 final
## CMS classification in combined COAD + READ using CMScaller
## Input matrix: expr_mat_symbol_clean.rds
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 1e6)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(matrixStats)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"
task_dir <- file.path(base_dir, "task1_RiboSis_NuS_analysis")

plot_dir  <- file.path(task_dir, "plots_step7_CMScaller")
table_dir <- file.path(task_dir, "tables_step7_CMScaller")
rds_dir   <- file.path(task_dir, "rds_step7_CMScaller")
qc_dir    <- file.path(task_dir, "qc_step7_CMScaller")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

annot4_file <- file.path(task_dir, "rds_step4", "step4_annot_main_tumor_normal_untreated.rds")
expr_symbol_file <- file.path(task_dir, "rds", "expr_mat_symbol_clean.rds")

stopifnot(file.exists(annot4_file))
stopifnot(file.exists(expr_symbol_file))

###############################################################################
## 1. Helper functions
###############################################################################
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

msg <- function(...) {
  cat(sprintf("[%s] ", timestamp()), ..., "\n", sep = "")
  flush.console()
}

save_pdf <- function(p, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width, height = height, units = "in",
    device = cairo_pdf, dpi = 300
  )
}

save_png <- function(p, filename, width = 10, height = 7, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width, height = height, units = "in",
    dpi = dpi
  )
}

save_csv <- function(df, filename) {
  write.csv(df, file.path(table_dir, filename), row.names = FALSE)
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1)
    )
}

cms_pal <- c(
  "CMS1" = "#d73027",
  "CMS2" = "#4575b4",
  "CMS3" = "#66a61e",
  "CMS4" = "#984ea3"
)

normalize_cms_label <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- toupper(x)
  
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x %in% c("CMS1", "CMS 1") ~ "CMS1",
    x %in% c("CMS2", "CMS 2") ~ "CMS2",
    x %in% c("CMS3", "CMS 3") ~ "CMS3",
    x %in% c("CMS4", "CMS 4") ~ "CMS4",
    x %in% c("NOLBL", "UNCLASSIFIED", "UNASSIGNED", "", "NA", "NONE") ~ NA_character_,
    TRUE ~ x
  )
}

safe_kw <- function(df, value_col, group_col, min_groups = 2, min_total = 10) {
  sub <- df %>%
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  
  if (nrow(sub) < min_total || length(unique(sub[[group_col]])) < min_groups) {
    return(data.frame(
      value = value_col,
      n = nrow(sub),
      n_groups = length(unique(sub[[group_col]])),
      p_value = NA_real_
    ))
  }
  
  kw <- kruskal.test(sub[[value_col]] ~ sub[[group_col]])
  data.frame(
    value = value_col,
    n = nrow(sub),
    n_groups = length(unique(sub[[group_col]])),
    p_value = kw$p.value
  )
}

safe_spearman <- function(df, x, y, min_n = 5) {
  sub <- df %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  
  if (nrow(sub) < min_n) {
    return(data.frame(
      x = x, y = y, n = nrow(sub),
      rho = NA_real_, p_value = NA_real_
    ))
  }
  
  ct <- suppressWarnings(cor.test(sub[[x]], sub[[y]], method = "spearman"))
  data.frame(
    x = x, y = y, n = nrow(sub),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

plot_corr_scatter <- function(df, x, y, color_var = NULL,
                              title = NULL, xlab = NULL, ylab = NULL,
                              point_size = 0.8, alpha = 0.5,
                              palette = NULL) {
  p <- ggplot(df, aes_string(x = x, y = y, color = color_var)) +
    geom_point(size = point_size, alpha = alpha) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    labs(title = title, x = xlab, y = ylab, color = NULL) +
    theme_pub(base_size = 12)
  
  if (!is.null(palette) && !is.null(color_var)) {
    present_levels <- intersect(names(palette), unique(as.character(df[[color_var]])))
    p <- p + scale_color_manual(values = palette[present_levels], drop = FALSE)
  }
  
  p
}

###############################################################################
## 2. Load data
###############################################################################
msg("Step 1/8: Loading annot4_main and cleaned expression matrix ...")

annot4_main <- readRDS(annot4_file)
expr_mat_symbol <- readRDS(expr_symbol_file)

expr_mat_symbol <- as.matrix(expr_mat_symbol)
storage.mode(expr_mat_symbol) <- "numeric"

msg("annot4_main dim: ", nrow(annot4_main), " x ", ncol(annot4_main))
msg("expr_mat_symbol_clean dim: ", nrow(expr_mat_symbol), " x ", ncol(expr_mat_symbol))
msg("NA count in expr_mat_symbol_clean: ", sum(is.na(expr_mat_symbol)))
msg("NaN count in expr_mat_symbol_clean: ", sum(is.nan(expr_mat_symbol)))

###############################################################################
## 3. Prepare combined COAD + READ tumor samples
###############################################################################
msg("Step 2/8: Preparing combined COAD + READ tumor samples ...")

crc_annot <- annot4_main %>%
  dplyr::filter(
    tumor_normal_group == "Tumor",
    cancer_abbr %in% c("COAD", "READ")
  ) %>%
  dplyr::mutate(
    cancer_abbr = as.character(cancer_abbr)
  )

msg("CRC tumor sample count: ", nrow(crc_annot))
print(table(crc_annot$cancer_abbr, useNA = "ifany"))

stopifnot(all(crc_annot$sample %in% colnames(expr_mat_symbol)))

expr_crc <- expr_mat_symbol[, crc_annot$sample, drop = FALSE]
expr_crc <- expr_crc[, crc_annot$sample, drop = FALSE]

save_csv(
  data.frame(sample = crc_annot$sample, cancer_abbr = crc_annot$cancer_abbr),
  "step7_crc_sample_manifest_COAD_READ_combined.csv"
)

###############################################################################
## 4. Convert SYMBOL to Entrez for CMScaller
###############################################################################
msg("Step 3/8: Converting gene symbols to Entrez IDs ...")

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(expr_crc),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

gene_map_clean <- gene_map %>%
  dplyr::filter(!is.na(ENTREZID), ENTREZID != "") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

map_idx <- match(rownames(expr_crc), gene_map_clean$SYMBOL)
entrez_ids <- gene_map_clean$ENTREZID[map_idx]

keep_idx <- which(!is.na(entrez_ids) & entrez_ids != "")
expr_crc_entrez <- expr_crc[keep_idx, , drop = FALSE]
entrez_ids_keep <- entrez_ids[keep_idx]

row_med <- matrixStats::rowMedians(expr_crc_entrez, na.rm = TRUE)

collapse_df <- data.frame(
  row_index = seq_len(nrow(expr_crc_entrez)),
  ENTREZID = entrez_ids_keep,
  median_expr = row_med,
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(ENTREZID, dplyr::desc(median_expr)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

expr_crc_entrez2 <- expr_crc_entrez[collapse_df$row_index, , drop = FALSE]
rownames(expr_crc_entrez2) <- collapse_df$ENTREZID
storage.mode(expr_crc_entrez2) <- "numeric"

msg("Expression matrix for CMScaller: ", nrow(expr_crc_entrez2), " genes x ", ncol(expr_crc_entrez2), " samples")
msg("Remaining NA count after SYMBOL->ENTREZ: ", sum(is.na(expr_crc_entrez2)))

saveRDS(expr_crc_entrez2, file.path(rds_dir, "step7_expr_crc_entrez_for_CMScaller.rds"))

qc_map <- data.frame(
  n_symbol_input = nrow(expr_crc),
  n_symbol_mapped = length(unique(gene_map_clean$SYMBOL)),
  n_entrez_after_filter = length(keep_idx),
  n_entrez_after_collapse = nrow(expr_crc_entrez2),
  na_count_after_collapse = sum(is.na(expr_crc_entrez2))
)
save_csv(qc_map, "step7_qc_symbol_to_entrez_mapping_for_CMScaller.csv")

###############################################################################
## 5. CMS classification using CMScaller
###############################################################################
msg("Step 4/8: Running CMScaller ...")

if (!requireNamespace("CMScaller", quietly = TRUE)) {
  stop(
    "CMScaller is not installed.\n",
    "Install with:\n",
    "remotes::install_github('peterawe/CMScaller')"
  )
}

cms_pred <- CMScaller::CMScaller(
  emat = expr_crc_entrez2,
  RNAseq = TRUE,
  FDR = 0.05,
  doPlot = FALSE,
  verbose = TRUE
)

saveRDS(cms_pred, file.path(rds_dir, "step7_raw_CMScaller_output.rds"))
msg("CMScaller output columns: ", paste(colnames(cms_pred), collapse = ", "))

cms_df <- data.frame(
  sample = colnames(expr_crc_entrez2),
  stringsAsFactors = FALSE
)

if ("prediction" %in% colnames(cms_pred)) {
  cms_df$CMS_predicted <- as.character(cms_pred$prediction)
} else if ("class" %in% colnames(cms_pred)) {
  cms_df$CMS_predicted <- as.character(cms_pred$class)
} else if ("CMS" %in% colnames(cms_pred)) {
  cms_df$CMS_predicted <- as.character(cms_pred$CMS)
} else {
  stop("Cannot find CMS class column in CMScaller output.")
}

for (nm in c("FDR", "p.value", "corP", "min.dist", "prediction")) {
  if (nm %in% colnames(cms_pred)) {
    cms_df[[nm]] <- cms_pred[[nm]]
  }
}

cms_df <- cms_df %>%
  dplyr::mutate(
    CMS_predicted = normalize_cms_label(CMS_predicted),
    CMS = factor(CMS_predicted, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
  )

crc_cms <- crc_annot %>%
  dplyr::left_join(cms_df, by = "sample") %>%
  dplyr::mutate(
    CMS = factor(as.character(CMS), levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
  )

msg("Final CMS counts (combined COAD + READ):")
print(table(crc_cms$CMS, useNA = "ifany"))

msg("CMScaller raw predicted counts:")
print(table(cms_df$CMS_predicted, useNA = "ifany"))

saveRDS(crc_cms, file.path(rds_dir, "step7_crc_cmscaller_dataset.rds"))
save_csv(crc_cms, "step7_crc_cmscaller_annotation_merged.csv")
save_csv(as.data.frame(table(crc_cms$CMS, useNA = "ifany")),
         "step7_crc_cmscaller_counts_final.csv")
save_csv(as.data.frame(table(cms_df$CMS_predicted, useNA = "ifany")),
         "step7_crc_cmscaller_predicted_raw_counts.csv")

###############################################################################
## 6. CMS association with NuStress and RiboSis
###############################################################################
msg("Step 5/8: Testing CMS association with NuStress and RiboSis ...")

crc_cms_main <- crc_cms %>%
  dplyr::filter(CMS %in% c("CMS1", "CMS2", "CMS3", "CMS4")) %>%
  dplyr::mutate(
    CMS = factor(CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
  )

kw_nus <- safe_kw(crc_cms_main, "NuStress_z", "CMS")
kw_ribo <- safe_kw(crc_cms_main, "RiboSis_z", "CMS")

kw_summary <- dplyr::bind_rows(kw_nus, kw_ribo) %>%
  dplyr::mutate(
    p_label = dplyr::case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

save_csv(kw_summary, "step7_kw_CMScaller_vs_NuStress_RiboSis.csv")

pairwise_nus <- ggpubr::compare_means(
  NuStress_z ~ CMS,
  data = crc_cms_main,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

pairwise_ribo <- ggpubr::compare_means(
  RiboSis_z ~ CMS,
  data = crc_cms_main,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

save_csv(pairwise_nus, "step7_pairwise_wilcox_NuStress_by_CMScaller.csv")
save_csv(pairwise_ribo, "step7_pairwise_wilcox_RiboSis_by_CMScaller.csv")

###############################################################################
## 7. Visualization
###############################################################################
msg("Step 6/8: Drawing CMScaller figures ...")

p_cms_count <- ggplot(crc_cms, aes(x = CMS, fill = CMS)) +
  geom_bar(color = "black", linewidth = 0.4) +
  scale_fill_manual(values = cms_pal, drop = FALSE) +
  labs(
    title = "CMS subtype distribution in combined COAD + READ",
    subtitle = "CMScaller-based classification",
    x = NULL,
    y = "Sample count"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

save_pdf(p_cms_count, "step7_CMScaller_counts_combined_COAD_READ.pdf", width = 6.0, height = 5.2)
save_png(p_cms_count, "step7_CMScaller_counts_combined_COAD_READ.png", width = 6.0, height = 5.2)

p_nus_cms <- ggplot(crc_cms_main, aes(x = CMS, y = NuStress_z, fill = CMS)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_fill_manual(values = cms_pal[c("CMS1", "CMS2", "CMS3", "CMS4")]) +
  labs(
    title = "NuStress across CMS subtypes",
    subtitle = paste0(
      "Combined COAD + READ tumors; Kruskal-Wallis P = ",
      kw_summary$p_label[kw_summary$value == "NuStress_z"]
    ),
    x = NULL,
    y = "NuStress_z"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

p_ribo_cms <- ggplot(crc_cms_main, aes(x = CMS, y = RiboSis_z, fill = CMS)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_fill_manual(values = cms_pal[c("CMS1", "CMS2", "CMS3", "CMS4")]) +
  labs(
    title = "RiboSis across CMS subtypes",
    subtitle = paste0(
      "Combined COAD + READ tumors; Kruskal-Wallis P = ",
      kw_summary$p_label[kw_summary$value == "RiboSis_z"]
    ),
    x = NULL,
    y = "RiboSis_z"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

save_pdf(p_nus_cms, "step7_CMScaller_vs_NuStress_violin_combined_COAD_READ.pdf", width = 6.4, height = 5.8)
save_png(p_nus_cms, "step7_CMScaller_vs_NuStress_violin_combined_COAD_READ.png", width = 6.4, height = 5.8)

save_pdf(p_ribo_cms, "step7_CMScaller_vs_RiboSis_violin_combined_COAD_READ.pdf", width = 6.4, height = 5.8)
save_png(p_ribo_cms, "step7_CMScaller_vs_RiboSis_violin_combined_COAD_READ.png", width = 6.4, height = 5.8)

p_corr_cms <- plot_corr_scatter(
  df = crc_cms_main,
  x = "RiboSis_z",
  y = "NuStress_z",
  color_var = "CMS",
  palette = cms_pal,
  title = "Relationship between NuStress and RiboSis across CMS subtypes",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

save_pdf(p_corr_cms, "step7_CMScaller_colored_correlation_NuStress_vs_RiboSis.pdf", width = 6.8, height = 5.8)
save_png(p_corr_cms, "step7_CMScaller_colored_correlation_NuStress_vs_RiboSis.png", width = 6.8, height = 5.8)

cms_corr_summary <- split(crc_cms_main, crc_cms_main$CMS) %>%
  lapply(function(dd) {
    res <- safe_spearman(dd, "NuStress_z", "RiboSis_z")
    res$CMS <- unique(as.character(dd$CMS))[1]
    res
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::select(CMS, n, rho, p_value)

save_csv(cms_corr_summary, "step7_CMScaller_specific_correlation_summary.csv")

label_df <- cms_corr_summary %>%
  dplyr::mutate(
    label = paste0(
      "rho = ", round(rho, 2),
      "\nP = ",
      ifelse(p_value < 2.2e-16, "< 2.2e-16", signif(p_value, 2))
    )
  )

p_corr_facet <- ggplot(crc_cms_main, aes(x = RiboSis_z, y = NuStress_z, color = CMS)) +
  geom_point(size = 0.9, alpha = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  facet_wrap(~ CMS, nrow = 2) +
  scale_color_manual(values = cms_pal[c("CMS1", "CMS2", "CMS3", "CMS4")]) +
  geom_text(
    data = label_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.1, size = 3.8
  ) +
  labs(
    title = "CMS-specific relationships between NuStress and RiboSis",
    subtitle = "Combined COAD + READ tumors",
    x = "RiboSis_z",
    y = "NuStress_z"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "none")

save_pdf(p_corr_facet, "step7_CMScaller_facet_correlation_NuStress_vs_RiboSis.pdf", width = 8.5, height = 7.2)
save_png(p_corr_facet, "step7_CMScaller_facet_correlation_NuStress_vs_RiboSis.png", width = 8.5, height = 7.2)

crc_comp <- crc_cms %>%
  dplyr::filter(CMS %in% c("CMS1", "CMS2", "CMS3", "CMS4")) %>%
  dplyr::group_by(cancer_abbr, CMS) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup()

save_csv(crc_comp, "step7_COAD_READ_composition_by_CMScaller.csv")

p_comp <- ggplot(crc_comp, aes(x = cancer_abbr, y = freq, fill = CMS)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.35) +
  scale_fill_manual(values = cms_pal[c("CMS1", "CMS2", "CMS3", "CMS4")]) +
  labs(
    title = "CMS composition in COAD and READ",
    subtitle = "Classification performed on combined CRC cohort",
    x = NULL,
    y = "Proportion",
    fill = "CMS"
  ) +
  theme_pub(base_size = 12)

save_pdf(p_comp, "step7_COAD_READ_composition_by_CMScaller.pdf", width = 5.8, height = 5.4)
save_png(p_comp, "step7_COAD_READ_composition_by_CMScaller.png", width = 5.8, height = 5.4)

###############################################################################
## 8. Save summary
###############################################################################
msg("Step 7/8: Saving outputs ...")

saveRDS(kw_summary, file.path(rds_dir, "step7_kw_summary_CMScaller.rds"))
saveRDS(pairwise_nus, file.path(rds_dir, "step7_pairwise_nus_CMScaller.rds"))
saveRDS(pairwise_ribo, file.path(rds_dir, "step7_pairwise_ribo_CMScaller.rds"))
saveRDS(cms_corr_summary, file.path(rds_dir, "step7_cms_corr_summary_CMScaller.rds"))

summary_df <- data.frame(
  metric = c(
    "n_total_crc_tumor",
    "n_classified_main",
    "n_CMS1",
    "n_CMS2",
    "n_CMS3",
    "n_CMS4",
    "n_NA_unclassified"
  ),
  value = c(
    nrow(crc_cms),
    sum(!is.na(crc_cms$CMS)),
    sum(crc_cms$CMS == "CMS1", na.rm = TRUE),
    sum(crc_cms$CMS == "CMS2", na.rm = TRUE),
    sum(crc_cms$CMS == "CMS3", na.rm = TRUE),
    sum(crc_cms$CMS == "CMS4", na.rm = TRUE),
    sum(is.na(crc_cms$CMS))
  )
)

save_csv(summary_df, "step7_summary_metrics_CMScaller.csv")

msg("CMScaller analysis finished.")
print(kw_summary)
print(cms_corr_summary)
print(summary_df)
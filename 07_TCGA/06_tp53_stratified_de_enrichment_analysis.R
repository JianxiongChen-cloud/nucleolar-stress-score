###############################################################################
## TP53-stratified NuStress analysis
## Final integrated script for publication-ready analysis and figures
##
## Goals:
## 1) Test whether NuStress remains associated with RiboSis in TP53-mut tumors
## 2) Assess whether this association is consistent across cancer types
## 3) Compare RiboSis between NuStress-high and NuStress-low TP53-mut tumors
## 4) Perform DE and GO analyses in TP53-mut tumors stratified by NuStress
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 1e6)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(forcats)
  library(broom)
  library(purrr)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"
task_dir <- file.path(base_dir, "task1_RiboSis_NuS_analysis")

plot_dir  <- file.path(task_dir, "plots_step6_tp53_final")
table_dir <- file.path(task_dir, "tables_step6_tp53_final")
rds_dir   <- file.path(task_dir, "rds_step6_tp53_final")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

annot4_file <- file.path(task_dir, "rds_step4", "step4_annot_main_tumor_normal_untreated.rds")
expr_file   <- file.path(base_dir, "processed_data", "expr_tcga_tumor_only_final.rds")

stopifnot(file.exists(annot4_file))
stopifnot(file.exists(expr_file))

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

save_pdf <- function(p, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    dpi = 300
  )
}

save_png <- function(p, filename, width = 7, height = 5, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

parse_tp53_status <- function(df) {
  df %>%
    dplyr::mutate(
      TP53_status = dplyr::case_when(
        TP53_mut %in% c("WT", "Wildtype", "Wild-Type", "Wild Type") ~ "WT",
        TP53_mut %in% c("Mut", "MUT", "Mutant", "Mutation", "Altered") ~ "Mut",
        TP53_mut %in% c(TRUE, "TRUE") ~ "Mut",
        TP53_mut %in% c(FALSE, "FALSE") ~ "WT",
        TP53 %in% c("WT", "Wildtype", "Wild-Type", "Wild Type") ~ "WT",
        TP53 %in% c("Mut", "MUT", "Mutant", "Mutation", "Altered") ~ "Mut",
        TRUE ~ NA_character_
      )
    )
}

safe_spearman <- function(df, x, y, min_n = 10) {
  sub <- df %>% dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  if (nrow(sub) < min_n) {
    return(data.frame(
      n = nrow(sub),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ct <- suppressWarnings(cor.test(sub[[x]], sub[[y]], method = "spearman"))
  data.frame(
    n = nrow(sub),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

add_p_label <- function(df, p_col = "p_value") {
  df %>%
    mutate(
      p_label = case_when(
        is.na(.data[[p_col]]) ~ "NA",
        .data[[p_col]] < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", .data[[p_col]])
      )
    )
}

format_p_label <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("< 2.2e-16")
  if (p < 0.001) return("< 0.001")
  paste0("= ", signif(p, 2))
}

plot_corr_scatter <- function(df, x, y, color = "#d95f02",
                              title = NULL, xlab = NULL, ylab = NULL,
                              point_size = 0.35, alpha = 0.18) {
  sub <- df %>% dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  ct <- suppressWarnings(cor.test(sub[[x]], sub[[y]], method = "spearman"))
  rho <- unname(ct$estimate)
  pval <- ct$p.value
  
  label_txt <- paste0(
    "rho = ", round(rho, 2),
    "\nP ", format_p_label(pval)
  )
  
  x_pos <- quantile(sub[[x]], 0.05, na.rm = TRUE)
  y_pos <- quantile(sub[[y]], 0.93, na.rm = TRUE)
  
  ggplot(sub, aes_string(x = x, y = y)) +
    geom_point(size = point_size, alpha = alpha, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    annotate(
      "text",
      x = x_pos, y = y_pos,
      label = label_txt,
      hjust = 0, vjust = 1, size = 4
    ) +
    labs(title = title, x = xlab, y = ylab) +
    theme_pub(base_size = 12) +
    theme(plot.title = element_text(size = 14, face = "bold"))
}

run_go_bp <- function(gene_vec, out_csv, map_csv = NULL) {
  gene_vec <- unique(as.character(gene_vec))
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  if (length(gene_vec) == 0) {
    message("Input gene vector is empty.")
    return(NULL)
  }
  
  gene_vec_clean <- sub("\\..*$", "", gene_vec)
  
  entrez_df <- tryCatch({
    bitr(
      gene_vec_clean,
      fromType = "ENSEMBL",
      toType = c("ENTREZID", "SYMBOL"),
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) {
    message("ID conversion failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(entrez_df) || nrow(entrez_df) == 0) {
    message("No valid Ensembl IDs could be mapped for GO enrichment.")
    return(NULL)
  }
  
  if (!is.null(map_csv)) {
    write.csv(entrez_df, map_csv, row.names = FALSE)
  }
  
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
  
  ego_df <- as.data.frame(ego)
  if (nrow(ego_df) > 0) {
    write.csv(ego_df, out_csv, row.names = FALSE)
  } else {
    message("GO enrichment returned no significant terms.")
  }
  
  return(ego)
}

make_volcano_plot <- function(de_df,
                              title = "TP53-mut: NuStress high vs low",
                              p_cut = 0.05,
                              fc_cut = 0.5,
                              n_label = 15) {
  plot_df <- de_df %>%
    dplyr::mutate(
      neglog10FDR = -log10(pmax(adj.P.Val, 1e-300)),
      sig_group = dplyr::case_when(
        adj.P.Val < p_cut & logFC > fc_cut  ~ "Up",
        adj.P.Val < p_cut & logFC < -fc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  label_df <- plot_df %>%
    dplyr::filter(sig_group != "NS") %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = n_label)
  
  ggplot(plot_df, aes(x = logFC, y = neglog10FDR, color = sig_group)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = -log10(p_cut), linetype = 2, linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3
    ) +
    scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey75")) +
    labs(
      title = title,
      x = "log2 fold change",
      y = expression(-log[10]("FDR")),
      color = NULL
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

###############################################################################
## 2. Load data
###############################################################################
msg("Step 1/10: Loading annotation and expression matrix ...")

annot4_main <- readRDS(annot4_file)
expr_mat <- readRDS(expr_file)
expr_mat <- as.matrix(expr_mat)

msg("Expression matrix loaded: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")

###############################################################################
## 3. Prepare tumor annotation with TP53 status
###############################################################################
msg("Step 2/10: Preparing tumor-only annotation ...")

tumor_df <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Tumor") %>%
  parse_tp53_status() %>%
  dplyr::filter(TP53_status %in% c("WT", "Mut")) %>%
  dplyr::mutate(
    TP53_status = factor(TP53_status, levels = c("WT", "Mut"))
  )

stopifnot("sample" %in% colnames(tumor_df))

saveRDS(tumor_df, file.path(rds_dir, "tumor_df_with_tp53.rds"))
save_csv(tumor_df, "tumor_df_with_tp53.csv")

###############################################################################
## 4. Match samples between annotation and expression matrix
###############################################################################
msg("Step 3/10: Matching samples ...")

common_samples <- intersect(colnames(expr_mat), tumor_df$sample)
msg("Matched samples: ", length(common_samples))

expr_mat <- expr_mat[, common_samples, drop = FALSE]
tumor_df <- tumor_df %>%
  dplyr::filter(sample %in% common_samples) %>%
  dplyr::arrange(match(sample, common_samples))

stopifnot(identical(colnames(expr_mat), tumor_df$sample))

###############################################################################
## 5. TP53-stratified correlations
###############################################################################
msg("Step 4/10: TP53-stratified correlations ...")

tumor_wt <- tumor_df %>%
  dplyr::filter(
    TP53_status == "WT",
    !is.na(NuStress_z),
    !is.na(RiboSis_z),
    !is.na(cancer_abbr),
    cancer_abbr != ""
  )

tumor_mut <- tumor_df %>%
  dplyr::filter(
    TP53_status == "Mut",
    !is.na(NuStress_z),
    !is.na(RiboSis_z),
    !is.na(cancer_abbr),
    cancer_abbr != ""
  )

corr_wt_all <- safe_spearman(tumor_wt, "NuStress_z", "RiboSis_z") %>%
  mutate(group = "WT")
corr_mut_all <- safe_spearman(tumor_mut, "NuStress_z", "RiboSis_z") %>%
  mutate(group = "Mut")

corr_tp53_summary <- bind_rows(corr_wt_all, corr_mut_all)
save_csv(corr_tp53_summary, "tp53_stratified_overall_correlation_NuStress_vs_RiboSis.csv")
saveRDS(corr_tp53_summary, file.path(rds_dir, "tp53_stratified_overall_corr.rds"))

p_corr_wt <- plot_corr_scatter(
  tumor_wt,
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#9ecae1",
  title = "TP53 WT tumors",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

p_corr_mut <- plot_corr_scatter(
  tumor_mut,
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#f4a582",
  title = "TP53-mut tumors",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

save_pdf(p_corr_wt, "tp53_wt_correlation_NuStress_vs_RiboSis.pdf", width = 6, height = 5.3)
save_pdf(p_corr_mut, "tp53_mut_correlation_NuStress_vs_RiboSis.pdf", width = 6, height = 5.3)

###############################################################################
## 6. Per-cancer correlation in TP53-mut tumors
###############################################################################
msg("Step 5/10: Per-cancer correlations in TP53-mut tumors ...")

save_csv(
  tumor_mut %>%
    dplyr::count(cancer_abbr, name = "n_mut") %>%
    dplyr::arrange(dplyr::desc(n_mut)),
  "tp53_mut_sample_counts_by_cancer.csv"
)
mut_count_by_cancer <- tumor_mut %>%
  dplyr::count(cancer_abbr, name = "n_mut") %>%
  dplyr::filter(n_mut >= 20)

cancers_use <- mut_count_by_cancer$cancer_abbr

corr_mut_pc <- tumor_mut %>%
  filter(cancer_abbr %in% cancers_use) %>%
  split(.$cancer_abbr) %>%
  lapply(function(dd) {
    safe_spearman(dd, "NuStress_z", "RiboSis_z") %>%
      mutate(cancer_abbr = unique(dd$cancer_abbr)[1])
  }) %>%
  bind_rows() %>%
  relocate(cancer_abbr) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_p_label("p_adj") %>%
  arrange(rho)

save_csv(corr_mut_pc, "tp53_mut_per_cancer_correlation_NuStress_vs_RiboSis.csv")
saveRDS(corr_mut_pc, file.path(rds_dir, "tp53_mut_per_cancer_corr.rds"))

corr_mut_pc_plot <- corr_mut_pc %>%
  mutate(
    cancer_abbr = factor(cancer_abbr, levels = cancer_abbr),
    label_y = ifelse(rho < 0, rho + 0.03, rho - 0.03)
  )

p_corr_mut_pc <- ggplot(corr_mut_pc_plot, aes(x = cancer_abbr, y = rho)) +
  geom_col(fill = "#f4a582", color = "black", width = 0.72) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = p_label), size = 3) +
  coord_flip() +
  labs(
    title = "Per-cancer correlation in TP53-mut tumors",
    x = "Cancer type",
    y = "Spearman rho"
  ) +
  theme_pub(base_size = 12) +
  theme(plot.title = element_text(size = 14, face = "bold"))

save_pdf(p_corr_mut_pc, "tp53_mut_per_cancer_correlation_barplot.pdf", width = 7.2, height = 7.5)
save_png(p_corr_mut_pc, "tp53_mut_per_cancer_correlation_barplot.png", width = 7.2, height = 7.5)

###############################################################################
## 7. NuStress-high vs low in TP53-mut tumors
###############################################################################
msg("Step 6/10: NuStress-high vs low grouping within TP53-mut tumors ...")

mut_meta <- tumor_mut %>%
  dplyr::mutate(
    NuStress_group = ifelse(
      NuStress_z >= median(NuStress_z, na.rm = TRUE),
      "High", "Low"
    ),
    NuStress_group = factor(NuStress_group, levels = c("Low", "High"))
  )

save_csv(mut_meta, "tp53_mut_meta.csv")
saveRDS(mut_meta, file.path(rds_dir, "tp53_mut_meta.rds"))

save_csv(
  mut_meta %>%
    dplyr::count(NuStress_group, name = "n"),
  "tp53_mut_NuStress_high_low_counts.csv"
)

wilcox_mut_ribo <- wilcox.test(RiboSis_z ~ NuStress_group, data = mut_meta)

mut_ribo_compare <- data.frame(
  comparison = "TP53_mut_RiboSis_by_NuStress_group",
  n_low = sum(mut_meta$NuStress_group == "Low"),
  n_high = sum(mut_meta$NuStress_group == "High"),
  p_value = wilcox_mut_ribo$p.value
) %>% add_p_label("p_value")

save_csv(mut_ribo_compare, "tp53_mut_RiboSis_by_NuStress_group_wilcox.csv")

p_mut_group_ribo <- ggplot(mut_meta, aes(x = NuStress_group, y = RiboSis_z, fill = NuStress_group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_fill_manual(values = c("Low" = "#fddbc7", "High" = "#d6604d")) +
  labs(
    title = "RiboSis in TP53-mut tumors\nstratified by NuStress",
    x = "NuStress group",
    y = "RiboSis_z"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(mut_meta$RiboSis_z, na.rm = TRUE) * 1.05,
    label = paste0("Wilcoxon P ", format_p_label(mut_ribo_compare$p_value)),
    size = 4
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

save_pdf(p_mut_group_ribo, "tp53_mut_RiboSis_by_NuStress_group.pdf", width = 4.4, height = 5.2)
save_png(p_mut_group_ribo, "tp53_mut_RiboSis_by_NuStress_group.png", width = 4.4, height = 5.2)

###############################################################################
## 8. Simple adjusted models
###############################################################################
msg("Step 7/10: Simplified adjusted models ...")

model_df <- tumor_df %>%
  filter(
    !is.na(NuStress_z),
    !is.na(RiboSis_z),
    !is.na(TP53_status),
    !is.na(cancer_abbr),
    cancer_abbr != ""
  ) %>%
  mutate(
    TP53_status = factor(TP53_status, levels = c("WT", "Mut")),
    cancer_abbr = factor(cancer_abbr)
  )

fit_main <- lm(RiboSis_z ~ NuStress_z + TP53_status + cancer_abbr, data = model_df)
fit_mut_only <- lm(RiboSis_z ~ NuStress_z + cancer_abbr, data = model_df %>% filter(TP53_status == "Mut"))

coef_main <- broom::tidy(fit_main) %>% mutate(model = "all_tumors_adjusted")
coef_mut_only <- broom::tidy(fit_mut_only) %>% mutate(model = "TP53_mut_only")

coef_key <- bind_rows(coef_main, coef_mut_only) %>%
  filter(term == "NuStress_z") %>%
  mutate(
    term_label = c(
      "NuStress effect in all tumors (adjusted for TP53 and cancer type)",
      "NuStress effect in TP53-mut tumors (adjusted for cancer type)"
    )
  )

save_csv(coef_main, "lm_RiboSis_by_NuStress_TP53_cancer_main.csv")
save_csv(coef_mut_only, "lm_RiboSis_by_NuStress_cancer_TP53mut_only.csv")
save_csv(coef_key, "lm_key_NuStress_effect_summary.csv")

###############################################################################
## 9. DE in TP53-mut tumors by NuStress group
###############################################################################
msg("Step 8/10: Differential expression in TP53-mut tumors ...")

expr_mut <- expr_mat[, mut_meta$sample, drop = FALSE]
mut_meta <- mut_meta[match(colnames(expr_mut), mut_meta$sample), , drop = FALSE]

stopifnot(identical(colnames(expr_mut), mut_meta$sample))

design <- model.matrix(~ NuStress_group + cancer_abbr, data = mut_meta)

fit <- lmFit(expr_mut, design)
fit <- eBayes(fit)

coef_name <- "NuStress_groupHigh"
if (!(coef_name %in% colnames(coef(fit)))) {
  stop("Coefficient ", coef_name, " not found. Available coefficients: ",
       paste(colnames(coef(fit)), collapse = ", "))
}

de_res <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "P"
)

de_res$gene <- rownames(de_res)
de_res <- de_res %>%
  dplyr::select(gene, everything()) %>%
  dplyr::arrange(adj.P.Val, desc(abs(logFC)))

save_csv(de_res, "tp53_mut_DE_NuStressHigh_vs_Low_limma.csv")
saveRDS(de_res, file.path(rds_dir, "tp53_mut_DE_NuStressHigh_vs_Low_limma.rds"))

deg_sig <- de_res %>%
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 0.5)

deg_up <- deg_sig %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::arrange(adj.P.Val)

deg_down <- deg_sig %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::arrange(adj.P.Val)

save_csv(deg_sig, "tp53_mut_DE_NuStressHigh_vs_Low_limma_sig.csv")
save_csv(deg_up, "tp53_mut_DE_NuStressHigh_vs_Low_up.csv")
save_csv(deg_down, "tp53_mut_DE_NuStressHigh_vs_Low_down.csv")

p_volcano <- make_volcano_plot(
  de_res,
  title = "TP53-mut tumors: NuStress-high vs NuStress-low",
  p_cut = 0.05,
  fc_cut = 0.5,
  n_label = 20
)

save_pdf(p_volcano, "tp53_mut_DE_volcano.pdf", width = 6.5, height = 5.5)
save_png(p_volcano, "tp53_mut_DE_volcano.png", width = 6.5, height = 5.5)

###############################################################################
## 10. GO BP enrichment and final figure
###############################################################################
msg("Step 9/10: GO BP enrichment ...")

genes_up <- unique(deg_up$gene)
genes_down <- unique(deg_down$gene)

ego_up <- run_go_bp(
  genes_up,
  file.path(table_dir, "tp53_mut_GO_BP_up.csv"),
  file.path(table_dir, "tp53_mut_GO_BP_up_gene_mapping.csv")
)

ego_down <- run_go_bp(
  genes_down,
  file.path(table_dir, "tp53_mut_GO_BP_down.csv"),
  file.path(table_dir, "tp53_mut_GO_BP_down_gene_mapping.csv")
)

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  p_go_up <- dotplot(ego_up, showCategory = 15) +
    ggtitle("GO BP: up in TP53-mut / NuStress-high") +
    theme_pub()
  save_pdf(p_go_up, "tp53_mut_GO_BP_up_dotplot.pdf", width = 8, height = 6)
  save_png(p_go_up, "tp53_mut_GO_BP_up_dotplot.png", width = 8, height = 6)
}

if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  p_go_down <- dotplot(ego_down, showCategory = 15) +
    ggtitle("GO BP: down in TP53-mut / NuStress-high") +
    theme_pub()
  save_pdf(p_go_down, "tp53_mut_GO_BP_down_dotplot.pdf", width = 8, height = 6)
  save_png(p_go_down, "tp53_mut_GO_BP_down_dotplot.png", width = 8, height = 6)
}

msg("Step 10/10: Assembling final publication-style figure ...")

p_top <- p_corr_wt + p_corr_mut + plot_layout(ncol = 2)
p_bottom <- p_mut_group_ribo + p_corr_mut_pc + plot_layout(ncol = 2, widths = c(1, 1.25))

p_final_tp53 <- p_top / p_bottom +
  plot_annotation(
    title = "NuStress remains associated with lower RiboSis in TP53-mut tumors",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

save_pdf(p_final_tp53, "Figure_TP53_NuStress_RiboSis_final.pdf", width = 12.5, height = 10)
save_png(p_final_tp53, "Figure_TP53_NuStress_RiboSis_final.png", width = 12.5, height = 10)
saveRDS(p_final_tp53, file.path(rds_dir, "Figure_TP53_NuStress_RiboSis_final.rds"))

summary_df <- data.frame(
  metric = c(
    "n_genes",
    "n_samples_matched",
    "n_tp53_mut_samples",
    "n_low_NuStress",
    "n_high_NuStress",
    "n_sig_deg",
    "n_up_deg",
    "n_down_deg"
  ),
  value = c(
    nrow(expr_mat),
    ncol(expr_mat),
    nrow(mut_meta),
    sum(mut_meta$NuStress_group == "Low"),
    sum(mut_meta$NuStress_group == "High"),
    nrow(deg_sig),
    nrow(deg_up),
    nrow(deg_down)
  )
)

save_csv(summary_df, "tp53_mut_final_summary.csv")
saveRDS(summary_df, file.path(rds_dir, "tp53_mut_final_summary.rds"))

msg("All TP53 final analyses finished successfully.")
print(corr_tp53_summary)
print(mut_ribo_compare)
print(coef_key)
print(summary_df)

run_kegg <- function(gene_vec, out_csv, map_csv = NULL) {
  gene_vec <- unique(as.character(gene_vec))
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  if (length(gene_vec) == 0) {
    message("Input gene vector is empty.")
    return(NULL)
  }
  
  gene_vec_clean <- sub("\\..*$", "", gene_vec)
  
  entrez_df <- tryCatch({
    bitr(
      gene_vec_clean,
      fromType = "ENSEMBL",
      toType = c("ENTREZID", "SYMBOL"),
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) {
    message("ID conversion failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(entrez_df) || nrow(entrez_df) == 0) {
    message("No valid Ensembl IDs could be mapped for KEGG enrichment.")
    return(NULL)
  }
  
  if (!is.null(map_csv)) {
    write.csv(entrez_df, map_csv, row.names = FALSE)
  }
  
  ekegg <- enrichKEGG(
    gene         = unique(entrez_df$ENTREZID),
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  ekegg_df <- as.data.frame(ekegg)
  if (nrow(ekegg_df) > 0) {
    write.csv(ekegg_df, out_csv, row.names = FALSE)
  } else {
    message("KEGG enrichment returned no significant terms.")
  }
  
  return(ekegg)
}

###############################################################################
## 10. GO BP + KEGG enrichment
###############################################################################
msg("Step 9/10: GO BP + KEGG enrichment ...")

genes_up <- unique(deg_up$gene)
genes_down <- unique(deg_down$gene)

## GO
ego_up <- run_go_bp(
  genes_up,
  file.path(table_dir, "tp53_mut_GO_BP_up.csv"),
  file.path(table_dir, "tp53_mut_GO_BP_up_gene_mapping.csv")
)

ego_down <- run_go_bp(
  genes_down,
  file.path(table_dir, "tp53_mut_GO_BP_down.csv"),
  file.path(table_dir, "tp53_mut_GO_BP_down_gene_mapping.csv")
)

## KEGG
ekegg_up <- run_kegg(
  genes_up,
  file.path(table_dir, "tp53_mut_KEGG_up.csv"),
  file.path(table_dir, "tp53_mut_KEGG_up_gene_mapping.csv")
)

ekegg_down <- run_kegg(
  genes_down,
  file.path(table_dir, "tp53_mut_KEGG_down.csv"),
  file.path(table_dir, "tp53_mut_KEGG_down_gene_mapping.csv")
)

## KEGG plots
if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) {
  p_kegg_up <- dotplot(ekegg_up, showCategory = 15) +
    ggtitle("KEGG: up in TP53-mut / NuStress-high") +
    theme_pub()
  
  save_pdf(p_kegg_up, "tp53_mut_KEGG_up_dotplot.pdf", width = 8, height = 6)
  save_png(p_kegg_up, "tp53_mut_KEGG_up_dotplot.png", width = 8, height = 6)
}

if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) {
  p_kegg_down <- dotplot(ekegg_down, showCategory = 15) +
    ggtitle("KEGG: down in TP53-mut / NuStress-high") +
    theme_pub()
  
  save_pdf(p_kegg_down, "tp53_mut_KEGG_down_dotplot.pdf", width = 8, height = 6)
  save_png(p_kegg_down, "tp53_mut_KEGG_down_dotplot.png", width = 8, height = 6)
}

step6_tp53_final_pvalue_table_list <- data.frame(
  Table = c(
    "tp53_stratified_overall_correlation_NuStress_vs_RiboSis.csv",
    "tp53_mut_per_cancer_correlation_NuStress_vs_RiboSis.csv",
    "tp53_mut_RiboSis_by_NuStress_group_wilcox.csv",
    "lm_RiboSis_by_NuStress_TP53_cancer_main.csv",
    "lm_RiboSis_by_NuStress_cancer_TP53mut_only.csv",
    "lm_key_NuStress_effect_summary.csv",
    "tp53_mut_DE_NuStressHigh_vs_Low_limma.csv",
    "tp53_mut_GO_BP_up.csv",
    "tp53_mut_GO_BP_down.csv",
    "tp53_mut_KEGG_up.csv",
    "tp53_mut_KEGG_down.csv"
  ),
  Analysis = c(
    "Overall Spearman correlation between NuStress_z and RiboSis_z stratified by TP53 status",
    "Per-cancer Spearman correlation between NuStress_z and RiboSis_z in TP53-mut tumors",
    "Wilcoxon test of RiboSis_z between NuStress-low and NuStress-high TP53-mut tumors",
    "Linear model coefficients for RiboSis_z by NuStress_z, TP53 status, and cancer type",
    "Linear model coefficients for RiboSis_z by NuStress_z and cancer type in TP53-mut tumors",
    "Key NuStress_z coefficient summary from linear models",
    "limma differential expression results for NuStress-high vs NuStress-low TP53-mut tumors",
    "GO BP enrichment for genes up in NuStress-high TP53-mut tumors",
    "GO BP enrichment for genes down in NuStress-high TP53-mut tumors",
    "KEGG enrichment for genes up in NuStress-high TP53-mut tumors",
    "KEGG enrichment for genes down in NuStress-high TP53-mut tumors"
  ),
  stringsAsFactors = FALSE
)
save_csv(step6_tp53_final_pvalue_table_list, "step6_tp53_final_pvalue_table_file_list.csv")

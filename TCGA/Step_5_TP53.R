###############################################################################
## Task 1: Pan-cancer analysis of RiboSis and NuS ######## Step5 final
## TP53 association analysis based on Step4 output
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
  library(forcats)
  library(scales)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"
task_dir <- file.path(base_dir, "task1_RiboSis_NuS_analysis")

plot_dir  <- file.path(task_dir, "plots_step5")
table_dir <- file.path(task_dir, "tables_step5")
rds_dir   <- file.path(task_dir, "rds_step5")
qc_dir    <- file.path(task_dir, "qc_step5")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

annot4_file <- file.path(task_dir, "rds_step4", "step4_annot_main_tumor_normal_untreated.rds")
stopifnot(file.exists(annot4_file))

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

safe_wilcox <- function(df, value_col, group_col, group1, group2, min_n_each = 3) {
  sub <- df %>%
    dplyr::filter(.data[[group_col]] %in% c(group1, group2)) %>%
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  
  n1 <- sum(sub[[group_col]] == group1)
  n2 <- sum(sub[[group_col]] == group2)
  
  if (n1 < min_n_each || n2 < min_n_each) {
    return(data.frame(
      value = value_col,
      group1 = group1,
      group2 = group2,
      n1 = n1,
      n2 = n2,
      p_value = NA_real_
    ))
  }
  
  wt <- wilcox.test(sub[[value_col]] ~ sub[[group_col]])
  data.frame(
    value = value_col,
    group1 = group1,
    group2 = group2,
    n1 = n1,
    n2 = n2,
    p_value = wt$p.value
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

add_p_label <- function(df, p_col = "p_adj") {
  df %>%
    dplyr::mutate(
      p_label = dplyr::case_when(
        is.na(.data[[p_col]]) ~ "NA",
        .data[[p_col]] < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", .data[[p_col]])
      )
    )
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

theme_violin <- function(base_size = 12) {
  theme_pub(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

fill_tp53 <- c(
  "WT" = "#9ecae1",
  "Mut" = "#f4a582"
)

plot_corr_scatter <- function(df, x, y, color = "#d95f02",
                              title = NULL, xlab = NULL, ylab = NULL,
                              cor_method = "spearman",
                              point_size = 0.35, alpha = 0.18) {
  sub <- df %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  
  ct <- suppressWarnings(cor.test(sub[[x]], sub[[y]], method = cor_method))
  rho <- unname(ct$estimate)
  pval <- ct$p.value
  
  label_txt <- paste0(
    "rho = ", round(rho, 2),
    ", p = ",
    ifelse(pval < 2.2e-16, "< 2.2e-16", signif(pval, 2))
  )
  
  x_pos <- quantile(sub[[x]], 0.05, na.rm = TRUE)
  y_pos <- quantile(sub[[y]], 0.95, na.rm = TRUE)
  
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
    theme_pub(base_size = 12)
}

###############################################################################
## 2. Load Step4 result
###############################################################################
msg("Step 1/7: Loading Step4 result ...")

annot4_main <- readRDS(annot4_file)

msg("Loaded annot4_main: ", nrow(annot4_main), " rows x ", ncol(annot4_main), " columns")

###############################################################################
## 3. Keep tumor samples and define TP53 status
###############################################################################
msg("Step 2/7: Preparing tumor-only TP53 dataset ...")

tumor_tp53 <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Tumor") %>%
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
  ) %>%
  dplyr::filter(TP53_status %in% c("WT", "Mut")) %>%
  dplyr::mutate(
    TP53_status = factor(TP53_status, levels = c("WT", "Mut"))
  )

msg("Tumor TP53 status counts:")
print(table(tumor_tp53$TP53_status, useNA = "ifany"))

qc_tp53 <- tumor_tp53 %>%
  dplyr::count(TP53_status, name = "n")

save_csv(qc_tp53, "step5_tp53_status_counts_overall.csv")
saveRDS(tumor_tp53, file.path(rds_dir, "step5_tumor_tp53_dataset.rds"))

###############################################################################
## 4. Overall TP53 WT vs Mut
###############################################################################
msg("Step 3/7: Overall TP53 WT vs Mut comparison ...")

overall_tp53_nus <- safe_wilcox(
  tumor_tp53,
  value_col = "NuStress_z",
  group_col = "TP53_status",
  group1 = "WT",
  group2 = "Mut"
)

overall_tp53_ribo <- safe_wilcox(
  tumor_tp53,
  value_col = "RiboSis_z",
  group_col = "TP53_status",
  group1 = "WT",
  group2 = "Mut"
)

overall_tp53_tests <- dplyr::bind_rows(overall_tp53_nus, overall_tp53_ribo) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_p_label("p_adj")

save_csv(overall_tp53_tests, "overall_tp53_wt_vs_mut_wilcox.csv")

p_tp53_overall_nus <- ggplot(
  tumor_tp53,
  aes(x = TP53_status, y = NuStress_z, fill = TP53_status)
) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_fill_manual(values = fill_tp53) +
  labs(
    title = "Overall TP53 association with NuStress",
    subtitle = "Tumor samples only",
    x = NULL,
    y = "NuStress_z"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(tumor_tp53$NuStress_z, na.rm = TRUE) * 1.05,
    label = paste0(
      "BH-adjusted P = ",
      overall_tp53_tests$p_label[overall_tp53_tests$value == "NuStress_z"]
    ),
    size = 4
  ) +
  theme_violin() +
  theme(legend.position = "none")

p_tp53_overall_ribo <- ggplot(
  tumor_tp53,
  aes(x = TP53_status, y = RiboSis_z, fill = TP53_status)
) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.5) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_fill_manual(values = fill_tp53) +
  labs(
    title = "Overall TP53 association with RiboSis",
    subtitle = "Tumor samples only",
    x = NULL,
    y = "RiboSis_z"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(tumor_tp53$RiboSis_z, na.rm = TRUE) * 1.05,
    label = paste0(
      "BH-adjusted P = ",
      overall_tp53_tests$p_label[overall_tp53_tests$value == "RiboSis_z"]
    ),
    size = 4
  ) +
  theme_violin() +
  theme(legend.position = "none")

save_pdf(p_tp53_overall_nus, "overall_tp53_wt_vs_mut_NuStress_z.pdf", width = 3, height = 5.8)

save_pdf(p_tp53_overall_ribo, "overall_tp53_wt_vs_mut_RiboSis_z.pdf", width = 3, height = 5.8)

###############################################################################
## 5. Per-cancer TP53 WT vs Mut
###############################################################################
msg("Step 4/7: Per-cancer TP53 WT vs Mut comparison ...")

tp53_count_df <- tumor_tp53 %>%
  dplyr::filter(!is.na(cancer_abbr), cancer_abbr != "") %>%
  dplyr::count(cancer_abbr, TP53_status, name = "n") %>%
  tidyr::pivot_wider(
    names_from = TP53_status,
    values_from = n,
    values_fill = 0
  )

if (!"WT" %in% colnames(tp53_count_df)) tp53_count_df$WT <- 0
if (!"Mut" %in% colnames(tp53_count_df)) tp53_count_df$Mut <- 0

tp53_count_df <- tp53_count_df %>%
  dplyr::mutate(has_both = WT > 0 & Mut > 0) %>%
  dplyr::arrange(dplyr::desc(has_both), cancer_abbr)

save_csv(tp53_count_df, "per_cancer_tp53_status_counts.csv")

cancers_tp53_both <- tp53_count_df %>%
  dplyr::filter(has_both) %>%
  dplyr::pull(cancer_abbr)

msg("Cancer types with both TP53 WT and Mut samples: ", length(cancers_tp53_both))
print(cancers_tp53_both)

tumor_tp53_pc <- tumor_tp53 %>%
  dplyr::filter(cancer_abbr %in% cancers_tp53_both)

per_cancer_tp53_test <- function(df, value_col) {
  split(df, df$cancer_abbr) %>%
    lapply(function(dd) {
      res <- safe_wilcox(
        dd,
        value_col = value_col,
        group_col = "TP53_status",
        group1 = "WT",
        group2 = "Mut",
        min_n_each = 3
      )
      res$cancer_abbr <- unique(dd$cancer_abbr)[1]
      res
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::relocate(cancer_abbr) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    add_p_label("p_adj") %>%
    dplyr::arrange(p_adj, p_value)
}

per_cancer_tp53_nus <- per_cancer_tp53_test(tumor_tp53_pc, "NuStress_z")
per_cancer_tp53_ribo <- per_cancer_tp53_test(tumor_tp53_pc, "RiboSis_z")

save_csv(per_cancer_tp53_nus, "per_cancer_tp53_wt_vs_mut_NuStress_z_wilcox.csv")
save_csv(per_cancer_tp53_ribo, "per_cancer_tp53_wt_vs_mut_RiboSis_z_wilcox.csv")

cancer_order_nus <- tumor_tp53_pc %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(med = median(NuStress_z[TP53_status == "Mut"], na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>%
  dplyr::pull(cancer_abbr)

cancer_order_ribo <- tumor_tp53_pc %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(med = median(RiboSis_z[TP53_status == "Mut"], na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>%
  dplyr::pull(cancer_abbr)

tumor_tp53_pc_nus <- tumor_tp53_pc %>%
  dplyr::mutate(cancer_abbr = factor(cancer_abbr, levels = cancer_order_nus))

tumor_tp53_pc_ribo <- tumor_tp53_pc %>%
  dplyr::mutate(cancer_abbr = factor(cancer_abbr, levels = cancer_order_ribo))

label_df_nus <- tumor_tp53_pc_nus %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(y = max(NuStress_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
  dplyr::left_join(
    per_cancer_tp53_nus %>% dplyr::select(cancer_abbr, p_label),
    by = "cancer_abbr"
  ) %>%
  dplyr::filter(!is.na(p_label), p_label != "NA")

label_df_ribo <- tumor_tp53_pc_ribo %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(y = max(RiboSis_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
  dplyr::left_join(
    per_cancer_tp53_ribo %>% dplyr::select(cancer_abbr, p_label),
    by = "cancer_abbr"
  ) %>%
  dplyr::filter(!is.na(p_label), p_label != "NA")

###############################################################################
## 6. Pan-cancer violin plots for TP53
###############################################################################
msg("Step 5/7: Drawing pan-cancer TP53 violin plots ...")

p_pan_tp53_nus <- ggplot(
  tumor_tp53_pc_nus,
  aes(x = cancer_abbr, y = NuStress_z, fill = TP53_status)
) +
  geom_violin(
    position = position_dodge(width = 0.85),
    trim = FALSE, scale = "width",
    color = "black", linewidth = 0.45, alpha = 0.9
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.14, outlier.shape = NA,
    fill = "white", color = "black", linewidth = 0.45
  ) +
  stat_summary(
    fun = median, geom = "point",
    position = position_dodge(width = 0.85),
    shape = 95, size = 4.8, color = "black"
  ) +
  geom_text(
    data = label_df_nus,
    aes(x = cancer_abbr, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  scale_fill_manual(values = fill_tp53) +
  labs(
    title = "Pan-cancer TP53 association with NuStress",
    subtitle = "Tumor samples only",
    x = "Cancer type",
    y = "NuStress_z",
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_violin()

p_pan_tp53_ribo <- ggplot(
  tumor_tp53_pc_ribo,
  aes(x = cancer_abbr, y = RiboSis_z, fill = TP53_status)
) +
  geom_violin(
    position = position_dodge(width = 0.85),
    trim = FALSE, scale = "width",
    color = "black", linewidth = 0.45, alpha = 0.9
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.14, outlier.shape = NA,
    fill = "white", color = "black", linewidth = 0.45
  ) +
  stat_summary(
    fun = median, geom = "point",
    position = position_dodge(width = 0.85),
    shape = 95, size = 4.8, color = "black"
  ) +
  geom_text(
    data = label_df_ribo,
    aes(x = cancer_abbr, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  scale_fill_manual(values = fill_tp53) +
  labs(
    title = "Pan-cancer TP53 association with RiboSis",
    subtitle = "Tumor samples only",
    x = "Cancer type",
    y = "RiboSis_z",
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_violin()

save_pdf(p_pan_tp53_nus, "pan_cancer_tp53_wt_vs_mut_NuStress_z.pdf", width = 14.5, height = 6.8)

save_pdf(p_pan_tp53_ribo, "pan_cancer_tp53_wt_vs_mut_RiboSis_z.pdf", width = 14.5, height = 6.8)

###############################################################################
## 7. Correlation between NuStress_z and RiboSis_z stratified by TP53 status
###############################################################################
msg("Step 6/7: Correlation analysis stratified by TP53 status ...")

corr_wt <- tumor_tp53 %>%
  dplyr::filter(TP53_status == "WT") %>%
  safe_spearman("NuStress_z", "RiboSis_z") %>%
  dplyr::mutate(group = "WT")

corr_mut <- tumor_tp53 %>%
  dplyr::filter(TP53_status == "Mut") %>%
  safe_spearman("NuStress_z", "RiboSis_z") %>%
  dplyr::mutate(group = "Mut")

corr_summary_tp53 <- dplyr::bind_rows(corr_wt, corr_mut) %>%
  dplyr::select(group, n, rho, p_value)

save_csv(corr_summary_tp53, "tp53_stratified_correlation_NuStress_vs_RiboSis.csv")

p_corr_wt <- plot_corr_scatter(
  df = tumor_tp53 %>% dplyr::filter(TP53_status == "WT"),
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#9ecae1",
  title = "Correlation in TP53 WT tumors",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

p_corr_mut <- plot_corr_scatter(
  df = tumor_tp53 %>% dplyr::filter(TP53_status == "Mut"),
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#f4a582",
  title = "Correlation in TP53 mutant tumors",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

save_pdf(p_corr_wt, "correlation_tp53_wt_NuStress_z_vs_RiboSis_z.pdf", width = 6.2, height = 5.6)

save_pdf(p_corr_mut, "correlation_tp53_mut_NuStress_z_vs_RiboSis_z.pdf", width = 6.2, height = 5.6)

###############################################################################
## 8. Save final objects and print concise summary
###############################################################################
msg("Step 7/7: Saving final Step5 objects ...")

saveRDS(overall_tp53_tests, file.path(rds_dir, "step5_overall_tp53_tests.rds"))
saveRDS(per_cancer_tp53_nus, file.path(rds_dir, "step5_per_cancer_tp53_nus_tests.rds"))
saveRDS(per_cancer_tp53_ribo, file.path(rds_dir, "step5_per_cancer_tp53_ribo_tests.rds"))
saveRDS(corr_summary_tp53, file.path(rds_dir, "step5_corr_summary_tp53.rds"))

msg("Overall TP53 WT vs Mut:")
print(overall_tp53_tests)

msg("TP53-stratified correlation summary:")
print(corr_summary_tp53)

msg("All Step5 analyses finished.")

###############################################################################
## 9. Quick checks
###############################################################################
table(tumor_tp53$TP53_status, useNA = "ifany")
table(tumor_tp53_pc$cancer_abbr, tumor_tp53_pc$TP53_status)
head(per_cancer_tp53_nus)
head(per_cancer_tp53_ribo)
corr_summary_tp53
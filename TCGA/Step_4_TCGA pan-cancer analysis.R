###############################################################################
## Task 1: Pan-cancer analysis of RiboSis and NuS ######## Step4 final polished
## Based on Step3 output: annot_with_ssgsea_scores.rds
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

plot_dir  <- file.path(task_dir, "plots_step4")
table_dir <- file.path(task_dir, "tables_step4")
rds_dir   <- file.path(task_dir, "rds_step4")
qc_dir    <- file.path(task_dir, "qc_step4")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

annot2_file <- file.path(task_dir, "rds", "annot_with_ssgsea_scores.rds")
stopifnot(file.exists(annot2_file))

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

fill_pal_2 <- c(
  "Normal" = "#9ecae1",
  "Tumor" = "#f4a582"
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
## 2. Load Step3 result
###############################################################################
msg("Step 1/8: Loading Step3 result ...")

annot2 <- readRDS(annot2_file)

msg("Loaded annot2: ", nrow(annot2), " rows x ", ncol(annot2), " columns")

###############################################################################
## 3. Remove treated samples (remove only TRUE; keep FALSE and NA)
###############################################################################
msg("Step 2/8: Removing treated samples (treated_flag_sample == TRUE) ...")

qc_treated_summary <- data.frame(
  stage = c("before_filter", "after_filter"),
  n_total = c(
    nrow(annot2),
    sum(is.na(annot2$treated_flag_sample) | annot2$treated_flag_sample == FALSE)
  ),
  n_FALSE = c(
    sum(annot2$treated_flag_sample == FALSE, na.rm = TRUE),
    sum(annot2$treated_flag_sample == FALSE, na.rm = TRUE)
  ),
  n_TRUE = c(
    sum(annot2$treated_flag_sample == TRUE, na.rm = TRUE),
    0
  ),
  n_NA = c(
    sum(is.na(annot2$treated_flag_sample)),
    sum(is.na(annot2$treated_flag_sample))
  )
)

save_csv(qc_treated_summary, "step4_treated_filter_qc.csv")

annot4 <- annot2 %>%
  dplyr::filter(is.na(treated_flag_sample) | treated_flag_sample == FALSE)

msg("Samples before filtering: ", nrow(annot2))
msg("Samples after filtering: ", nrow(annot4))

###############################################################################
## 3A. Visualize treated-sample removal
###############################################################################
msg("Step 2A/8: Visualizing treated-sample composition and filtering process ...")

annot2_plot <- annot2 %>%
  dplyr::mutate(
    treated_status_plot = dplyr::case_when(
      treated_flag_sample == TRUE  ~ "Treated",
      treated_flag_sample == FALSE ~ "Untreated",
      is.na(treated_flag_sample)   ~ "Unknown",
      TRUE ~ "Unknown"
    )
  )

treated_comp_before <- annot2_plot %>%
  dplyr::count(treated_status_plot, name = "n") %>%
  dplyr::mutate(
    treated_status_plot = factor(
      treated_status_plot,
      levels = c("Untreated", "Treated", "Unknown")
    )
  ) %>%
  dplyr::arrange(treated_status_plot) %>%
  dplyr::mutate(
    proportion = n / sum(n),
    percent = scales::percent(proportion, accuracy = 0.1)
  )

treated_filter_process <- tibble::tibble(
  stage = c("Before filter", "After filter"),
  n = c(nrow(annot2), nrow(annot4))
) %>%
  dplyr::mutate(
    stage = factor(stage, levels = c("Before filter", "After filter"))
  )

treated_filter_summary_for_plot <- treated_comp_before %>%
  dplyr::select(treated_status_plot, n, proportion, percent) %>%
  dplyr::rename(treated_status = treated_status_plot)

save_csv(treated_filter_summary_for_plot, "treated_filter_summary_for_plot.csv")

treated_fill_pal <- c(
  "Untreated" = "#80b1d3",
  "Treated"   = "#fb8072",
  "Unknown"   = "#bdbdbd"
)

stage_fill_pal <- c(
  "Before filter" = "#9ecae1",
  "After filter"  = "#f4a582"
)

###############################################################################
## Plot 1. Composition of treated status before filtering (pie chart)
###############################################################################
treated_comp_before <- treated_comp_before %>%
  dplyr::arrange(desc(treated_status_plot)) %>%
  dplyr::mutate(
    ymax = cumsum(proportion),
    ymin = dplyr::lag(ymax, default = 0),
    label_pos = (ymin + ymax) / 2,
    pie_label = paste0(as.character(treated_status_plot), "\n", percent, "\n(n=", n, ")")
  )

p_treated_comp <- ggplot(
  treated_comp_before,
  aes(x = "", y = proportion, fill = treated_status_plot)
) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y") +
  geom_text(
    aes(y = label_pos, label = pie_label),
    size = 4,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = treated_fill_pal, drop = FALSE) +
  labs(
    title = "Composition of treatment status before filtering",
    subtitle = "Proportion of treated, untreated, and treatment-unknown samples",
    fill = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

save_pdf(p_treated_comp, "treated_status_composition_before_filter_pie.pdf", width = 6.0, height = 5.8)


n_removed <- nrow(annot2) - nrow(annot4)
removed_label <- paste0("Removed treated\nn = ", scales::comma(n_removed))

p_filter_process <- ggplot(
  treated_filter_process,
  aes(x = stage, y = n, fill = stage)
) +
  geom_col(width = 0.62, color = "black", linewidth = 0.5) +
  geom_text(
    aes(label = paste0("n = ", scales::comma(n))),
    vjust = -0.45,
    size = 4.2
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = max(treated_filter_process$n) * 1.08,
    yend = max(treated_filter_process$n) * 1.08,
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = max(treated_filter_process$n) * 1.06,
    yend = max(treated_filter_process$n) * 1.08,
    linewidth = 0.7
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = max(treated_filter_process$n) * 1.06,
    yend = max(treated_filter_process$n) * 1.08,
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(treated_filter_process$n) * 1.11,
    label = removed_label,
    size = 4
  ) +
  scale_fill_manual(values = stage_fill_pal) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    labels = scales::comma_format()
  ) +
  labs(
    title = "Filtering process for treated samples",
    subtitle = "Samples with treated_flag_sample == TRUE were removed",
    x = NULL,
    y = "Number of samples"
  ) +
  theme_pub(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

save_pdf(p_filter_process, "treated_filter_before_after_counts.pdf", width = 5.4, height = 5.6)
save_png(p_filter_process, "treated_filter_before_after_counts.png", width = 5.4, height = 5.6, dpi = 300)

###############################################################################
## 4. Define tumor vs normal groups
###############################################################################
msg("Step 3/8: Defining tumor/normal groups ...")

annot4 <- annot4 %>%
  dplyr::mutate(
    tumor_normal_group = dplyr::case_when(
      sample_type_final %in% c("Tumor", "Primary Tumor", "Tumour", "Cancer", "TCGA_Tumor") ~ "Tumor",
      sample_type_final %in% c("Normal", "Solid Tissue Normal", "GTEx_Normal", "Adjacent Normal", "TCGA_Normal") ~ "Normal",
      stringr::str_detect(tolower(sample_type_final), "tumor|tumour|cancer") ~ "Tumor",
      stringr::str_detect(tolower(sample_type_final), "normal") ~ "Normal",
      TRUE ~ NA_character_
    )
  )

msg("Tumor/Normal counts:")
print(table(annot4$tumor_normal_group, useNA = "ifany"))

annot4_main <- annot4 %>%
  dplyr::filter(tumor_normal_group %in% c("Tumor", "Normal")) %>%
  dplyr::mutate(
    tumor_normal_group = factor(tumor_normal_group, levels = c("Normal", "Tumor"))
  )

saveRDS(annot4_main, file.path(rds_dir, "annot_step4_main_tumor_normal.rds"))

###############################################################################
## 5. Overall tumor vs normal comparison
###############################################################################
msg("Step 4/8: Overall tumor vs normal comparison ...")

overall_n <- annot4_main %>%
  dplyr::group_by(tumor_normal_group) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

overall_test_nus <- safe_wilcox(
  annot4_main,
  value_col = "NuStress_z",
  group_col = "tumor_normal_group",
  group1 = "Normal",
  group2 = "Tumor"
)

overall_test_ribo <- safe_wilcox(
  annot4_main,
  value_col = "RiboSis_z",
  group_col = "tumor_normal_group",
  group1 = "Normal",
  group2 = "Tumor"
)

overall_tests <- dplyr::bind_rows(overall_test_nus, overall_test_ribo) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_p_label("p_adj")

save_csv(overall_n, "overall_tumor_normal_counts.csv")
save_csv(overall_tests, "overall_tumor_normal_wilcox.csv")

p_overall_nus <- ggplot(
  annot4_main,
  aes(x = tumor_normal_group, y = NuStress_z, fill = tumor_normal_group)
) +
  geom_violin(
    trim = FALSE, scale = "width",
    color = "black", linewidth = 0.5, alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA,
    fill = "white", color = "black", linewidth = 0.5
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 95, size = 6, color = "black"
  ) +
  scale_fill_manual(values = fill_pal_2) +
  labs(
    title = "Overall comparison of NuStress",
    subtitle = "Untreated and treatment-unknown samples",
    x = NULL,
    y = "NuStress_z"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(annot4_main$NuStress_z, na.rm = TRUE) * 1.05,
    label = paste0(
      "BH-adjusted P = ",
      overall_tests$p_label[overall_tests$value == "NuStress_z"]
    ),
    size = 4
  ) +
  theme_violin() +
  theme(legend.position = "none")

p_overall_ribo <- ggplot(
  annot4_main,
  aes(x = tumor_normal_group, y = RiboSis_z, fill = tumor_normal_group)
) +
  geom_violin(
    trim = FALSE, scale = "width",
    color = "black", linewidth = 0.5, alpha = 0.9
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA,
    fill = "white", color = "black", linewidth = 0.5
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 95, size = 6, color = "black"
  ) +
  scale_fill_manual(values = fill_pal_2) +
  labs(
    title = "Overall comparison of RiboSis",
    subtitle = "Untreated and treatment-unknown samples",
    x = NULL,
    y = "RiboSis_z"
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(annot4_main$RiboSis_z, na.rm = TRUE) * 1.05,
    label = paste0(
      "BH-adjusted P = ",
      overall_tests$p_label[overall_tests$value == "RiboSis_z"]
    ),
    size = 4
  ) +
  theme_violin() +
  theme(legend.position = "none")

save_pdf(p_overall_nus, "overall_tumor_vs_normal_NuStress_z.pdf", width = 3, height = 5.8)
save_pdf(p_overall_ribo, "overall_tumor_vs_normal_RiboSis_z.pdf", width = 3, height = 5.8)

###############################################################################
## 6. Per-cancer tumor vs matched normal comparison
###############################################################################
msg("Step 5/8: Per-cancer tumor vs matched normal comparison ...")

tcga_gtex_map <- tibble::tribble(
  ~cancer_abbr, ~normal_tissue,
  "ACC",  "Adrenal Gland",
  "BLCA", "Bladder",
  "BRCA", "Breast",
  "CESC", "Cervix Uteri",
  "CHOL", "Liver",
  "COAD", "Colon",
  "DLBC", NA_character_,
  "ESCA", "Esophagus",
  "GBM",  "Brain",
  "HNSC", "Salivary Gland",
  "KICH", "Kidney",
  "KIRC", "Kidney",
  "KIRP", "Kidney",
  "LAML", "Blood",
  "LGG",  "Brain",
  "LIHC", "Liver",
  "LUAD", "Lung",
  "LUSC", "Lung",
  "MESO", NA_character_,
  "OV",   "Ovary",
  "PAAD", "Pancreas",
  "PCPG", "Adrenal Gland",
  "PRAD", "Prostate",
  "READ", "Colon",
  "SARC", "Muscle",
  "SKCM", "Skin",
  "STAD", "Stomach",
  "TGCT", "Testis",
  "THCA", "Thyroid",
  "THYM", "Thymus",
  "UCEC", "Uterus",
  "UCS",  "Uterus",
  "UVM",  NA_character_
)

tumor_df <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Tumor", !is.na(cancer_abbr), cancer_abbr != "") %>%
  dplyr::left_join(tcga_gtex_map, by = "cancer_abbr")

normal_df <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Normal", !is.na(tissue), tissue != "")

cancer_count_df <- tumor_df %>%
  dplyr::filter(!is.na(normal_tissue), normal_tissue != "") %>%
  dplyr::group_by(cancer_abbr, normal_tissue) %>%
  dplyr::summarise(Tumor = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(
    normal_df %>%
      dplyr::group_by(tissue) %>%
      dplyr::summarise(Normal = dplyr::n(), .groups = "drop"),
    by = c("normal_tissue" = "tissue")
  ) %>%
  dplyr::mutate(
    Normal = ifelse(is.na(Normal), 0L, Normal),
    has_both = Tumor > 0 & Normal > 0
  ) %>%
  dplyr::arrange(desc(has_both), cancer_abbr)

save_csv(cancer_count_df, "per_cancer_tumor_normal_counts.csv")

cancers_with_both <- cancer_count_df %>%
  dplyr::filter(has_both) %>%
  dplyr::pull(cancer_abbr)

msg("Cancer types with both tumor and matched normal samples: ", length(cancers_with_both))
print(cancers_with_both)

annot4_match_tumor <- tumor_df %>%
  dplyr::filter(cancer_abbr %in% cancers_with_both) %>%
  dplyr::mutate(match_group = cancer_abbr)

annot4_match_normal <- normal_df %>%
  dplyr::inner_join(
    tcga_gtex_map %>%
      dplyr::filter(cancer_abbr %in% cancers_with_both, !is.na(normal_tissue)) %>%
      dplyr::transmute(
        tissue = normal_tissue,
        match_group = cancer_abbr
      ),
    by = "tissue",
    relationship = "many-to-many"
  )

annot4_match <- dplyr::bind_rows(
  annot4_match_tumor,
  annot4_match_normal
) %>%
  dplyr::mutate(
    match_group = factor(match_group, levels = unique(cancers_with_both)),
    tumor_normal_group = factor(tumor_normal_group, levels = c("Normal", "Tumor"))
  )

saveRDS(annot4_match, file.path(rds_dir, "step4_per_cancer_matched_dataset.rds"))

per_cancer_test <- function(df, value_col) {
  split(df, df$match_group) %>%
    lapply(function(dd) {
      res <- safe_wilcox(
        dd,
        value_col = value_col,
        group_col = "tumor_normal_group",
        group1 = "Normal",
        group2 = "Tumor",
        min_n_each = 3
      )
      res$cancer_abbr <- unique(as.character(dd$match_group))[1]
      res
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::relocate(cancer_abbr) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    add_p_label("p_adj") %>%
    dplyr::arrange(p_adj, p_value)
}

per_cancer_nus <- per_cancer_test(annot4_match, "NuStress_z")
per_cancer_ribo <- per_cancer_test(annot4_match, "RiboSis_z")

save_csv(per_cancer_nus, "per_cancer_tumor_normal_NuStress_z_wilcox.csv")
save_csv(per_cancer_ribo, "per_cancer_tumor_normal_RiboSis_z_wilcox.csv")

cancer_order_nus <- annot4_match %>%
  dplyr::filter(tumor_normal_group == "Tumor") %>%
  dplyr::group_by(match_group) %>%
  dplyr::summarise(med = median(NuStress_z, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(med)) %>%
  dplyr::pull(match_group) %>%
  as.character()

cancer_order_ribo <- annot4_match %>%
  dplyr::filter(tumor_normal_group == "Tumor") %>%
  dplyr::group_by(match_group) %>%
  dplyr::summarise(med = median(RiboSis_z, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(med)) %>%
  dplyr::pull(match_group) %>%
  as.character()

annot4_match_nus <- annot4_match %>%
  dplyr::mutate(match_group = factor(as.character(match_group), levels = cancer_order_nus))

annot4_match_ribo <- annot4_match %>%
  dplyr::mutate(match_group = factor(as.character(match_group), levels = cancer_order_ribo))

label_df_nus <- annot4_match_nus %>%
  dplyr::group_by(match_group) %>%
  dplyr::summarise(y = max(NuStress_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
  dplyr::left_join(
    per_cancer_nus %>% dplyr::select(cancer_abbr, p_label),
    by = c("match_group" = "cancer_abbr")
  ) %>%
  dplyr::filter(!is.na(p_label), p_label != "NA")

label_df_ribo <- annot4_match_ribo %>%
  dplyr::group_by(match_group) %>%
  dplyr::summarise(y = max(RiboSis_z, na.rm = TRUE) + 0.25, .groups = "drop") %>%
  dplyr::left_join(
    per_cancer_ribo %>% dplyr::select(cancer_abbr, p_label),
    by = c("match_group" = "cancer_abbr")
  ) %>%
  dplyr::filter(!is.na(p_label), p_label != "NA")

###############################################################################
## 7. Pan-cancer violin plots
###############################################################################
msg("Step 6/8: Drawing pan-cancer violin plots ...")

p_pan_nus <- ggplot(
  annot4_match_nus,
  aes(x = match_group, y = NuStress_z, fill = tumor_normal_group)
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
    aes(x = match_group, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  scale_fill_manual(values = fill_pal_2) +
  labs(
    title = "Pan-cancer comparison of NuStress",
    subtitle = "Untreated and treatment-unknown samples; TCGA tumors matched to GTEx normal tissues",
    x = "Cancer type",
    y = "NuStress_z",
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_violin()

p_pan_ribo <- ggplot(
  annot4_match_ribo,
  aes(x = match_group, y = RiboSis_z, fill = tumor_normal_group)
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
    aes(x = match_group, y = y, label = p_label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  scale_fill_manual(values = fill_pal_2) +
  labs(
    title = "Pan-cancer comparison of RiboSis",
    subtitle = "Untreated and treatment-unknown samples; TCGA tumors matched to GTEx normal tissues",
    x = "Cancer type",
    y = "RiboSis_z",
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_violin()

save_pdf(p_pan_nus, "pan_cancer_violin_NuStress_z_tumor_vs_normal.pdf", width = 14.5, height = 5)
save_pdf(p_pan_ribo, "pan_cancer_violin_RiboSis_z_tumor_vs_normal.pdf", width = 14.5, height = 5)

###############################################################################
## 8. Correlation between NuStress_z and RiboSis_z in all normals and all tumors
###############################################################################
msg("Step 7/8: Correlation analysis of NuStress_z vs RiboSis_z ...")

corr_normal <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Normal") %>%
  safe_spearman("NuStress_z", "RiboSis_z") %>%
  dplyr::mutate(group = "Normal")

corr_tumor <- annot4_main %>%
  dplyr::filter(tumor_normal_group == "Tumor") %>%
  safe_spearman("NuStress_z", "RiboSis_z") %>%
  dplyr::mutate(group = "Tumor")

corr_summary <- dplyr::bind_rows(corr_normal, corr_tumor) %>%
  dplyr::select(group, n, rho, p_value)

save_csv(corr_summary, "overall_normal_vs_tumor_correlation_NuStress_vs_RiboSis.csv")

p_corr_normal <- plot_corr_scatter(
  df = annot4_main %>% dplyr::filter(tumor_normal_group == "Normal"),
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#4c72b0",
  title = "Correlation in all normal samples",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

p_corr_tumor <- plot_corr_scatter(
  df = annot4_main %>% dplyr::filter(tumor_normal_group == "Tumor"),
  x = "RiboSis_z",
  y = "NuStress_z",
  color = "#ef8a62",
  title = "Correlation in all tumor samples",
  xlab = "RiboSis_z",
  ylab = "NuStress_z"
)

save_pdf(p_corr_normal, "correlation_normal_NuStress_z_vs_RiboSis_z.pdf", width = 5, height = 5)
save_pdf(p_corr_tumor, "correlation_tumor_NuStress_z_vs_RiboSis_z.pdf", width = 5, height = 5)

###############################################################################
## 9. Save final objects and print concise summary
###############################################################################
msg("Step 8/8: Saving final Step4 objects ...")

saveRDS(overall_tests, file.path(rds_dir, "step4_overall_tests.rds"))
saveRDS(per_cancer_nus, file.path(rds_dir, "step4_per_cancer_nus_tests.rds"))
saveRDS(per_cancer_ribo, file.path(rds_dir, "step4_per_cancer_ribo_tests.rds"))
saveRDS(corr_summary, file.path(rds_dir, "step4_corr_summary.rds"))
saveRDS(annot4_main, file.path(rds_dir, "step4_annot_main_tumor_normal_untreated.rds"))

msg("Overall tumor vs normal:")
print(overall_tests)

msg("Normal vs tumor correlation summary:")
print(corr_summary)

msg("All Step4 analyses finished.")

###############################################################################
## 10. Quick checks
###############################################################################
table(annot4_main$tumor_normal_group, useNA = "ifany")
table(annot4_match$match_group, annot4_match$tumor_normal_group)
head(per_cancer_nus)
head(per_cancer_ribo)
corr_summary
###############################################################################
## Task 1: Pan-cancer analysis of RiboSis and NuS
## Step6 extended final
## Nucleolar functional stratification + joint survival + single-marker survival
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
  library(survival)
  library(survminer)
  library(broom)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"
task_dir <- file.path(base_dir, "task1_RiboSis_NuS_analysis")

plot_dir  <- file.path(task_dir, "plots_step6")
table_dir <- file.path(task_dir, "tables_step6")
rds_dir   <- file.path(task_dir, "rds_step6")
qc_dir    <- file.path(task_dir, "qc_step6")

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

nufunc_pal <- c(
  "H-H" = "#DB2462",
  "H-L" = "#47957F",
  "L-H" = "#7DC18C",
  "L-L" = "#F0C6DD"
)

single_pal <- c(
  "High" = "#DB2462",
  "Low"  = "#47957F"
)

tp53_pal <- c(
  "WT" = "#9ecae1",
  "Mut" = "#f4a582"
)

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
    ", p = ", ifelse(pval < 2.2e-16, "< 2.2e-16", signif(pval, 2))
  )
  
  x_pos <- quantile(sub[[x]], 0.05, na.rm = TRUE)
  y_pos <- quantile(sub[[y]], 0.95, na.rm = TRUE)
  
  ggplot(sub, aes_string(x = x, y = y)) +
    geom_point(size = point_size, alpha = alpha, color = color) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    annotate("text", x = x_pos, y = y_pos, label = label_txt,
             hjust = 0, vjust = 1, size = 4) +
    labs(title = title, x = xlab, y = ylab) +
    theme_pub(base_size = 12)
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 2.2e-16, "< 2.2e-16", as.character(signif(p, 3))))
}

extract_hr_from_cox <- function(sub_data, time_col, event_col, group_col, ref_level = "Low") {
  sub2 <- sub_data %>%
    dplyr::filter(!is.na(.data[[time_col]]),
                  !is.na(.data[[event_col]]),
                  !is.na(.data[[group_col]])) %>%
    dplyr::mutate(
      .group = factor(.data[[group_col]], levels = c(ref_level, setdiff(unique(.data[[group_col]]), ref_level)))
    )
  
  if (nrow(sub2) < 10 || length(unique(sub2$.group)) < 2) {
    return(data.frame(
      HR = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      cox_p = NA_real_
    ))
  }
  
  fit <- tryCatch(
    coxph(Surv(sub2[[time_col]], sub2[[event_col]]) ~ .group, data = sub2),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      HR = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      cox_p = NA_real_
    ))
  }
  
  tb <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  if (nrow(tb) < 1) {
    return(data.frame(
      HR = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      cox_p = NA_real_
    ))
  }
  
  data.frame(
    HR = tb$estimate[1],
    lower95 = tb$conf.low[1],
    upper95 = tb$conf.high[1],
    cox_p = tb$p.value[1]
  )
}

run_survival_by_group <- function(
    data,
    group_col,
    analysis_name,
    legend_title,
    palette_vec,
    out_prefix,
    min_n_total = 20,
    min_n_each = 3,
    time_col = "os_time",
    event_col = "os",
    cancer_col = "cancer_abbr",
    plot_only_sig = TRUE,
    sig_cutoff = 0.05
) {
  msg("Running survival analysis: ", analysis_name)
  
  surv_df <- data %>%
    dplyr::filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(.data[[cancer_col]]),
      .data[[cancer_col]] != "",
      !is.na(.data[[group_col]])
    ) %>%
    dplyr::mutate(
      cancer_abbr = trimws(as.character(.data[[cancer_col]])),
      .group = as.character(.data[[group_col]])
    )
  
  cancer_count_check <- surv_df %>%
    dplyr::count(cancer_abbr, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  
  save_csv(cancer_count_check, paste0(out_prefix, "_survival_input_cancer_counts.csv"))
  
  cancer_ids <- unique(surv_df$cancer_abbr)
  results_list <- list()
  plots_list <- list()
  
  for (ct in cancer_ids) {
    sub_data <- surv_df %>%
      dplyr::filter(cancer_abbr == .env$ct) %>%
      dplyr::filter(!is.na(.group))
    
    msg("Checking ", analysis_name, " | ", ct, " | n = ", nrow(sub_data))
    
    if (nrow(sub_data) < min_n_total) {
      msg("Skipping ", ct, ": insufficient total sample size")
      next
    }
    
    group_counts <- table(sub_data$.group)
    if (length(group_counts) < 2) {
      msg("Skipping ", ct, ": only one group present")
      next
    }
    
    if (any(group_counts < min_n_each)) {
      msg("Skipping ", ct, ": some groups have < ", min_n_each, " samples")
      next
    }
    
    sub_data$.group <- factor(sub_data$.group, levels = names(group_counts))
    
    ## Explicitly write time/event into data.frame to avoid losing the surv_obj environment
    sub_data <- sub_data %>%
      dplyr::mutate(
        .time  = .data[[time_col]],
        .event = .data[[event_col]]
      )
    
    fit <- tryCatch(
      survival::survfit(
        survival::Surv(.time, .event) ~ .group,
        data = sub_data
      ),
      error = function(e) {
        msg("survfit failed for ", ct, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(fit)) next
    
    logrank_test <- tryCatch(
      survival::survdiff(
        survival::Surv(.time, .event) ~ .group,
        data = sub_data
      ),
      error = function(e) {
        msg("survdiff failed for ", ct, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(logrank_test)) next
    
    p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
    
    res_row <- data.frame(
      analysis = analysis_name,
      cancer_abbr = ct,
      n = nrow(sub_data),
      n_groups = length(group_counts),
      p_value = p_value,
      significant = ifelse(p_value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    count_df <- as.data.frame(group_counts, stringsAsFactors = FALSE)
    colnames(count_df) <- c("group", "count")
    for (i in seq_len(nrow(count_df))) {
      nm <- paste0("count_", make.names(count_df$group[i]))
      res_row[[nm]] <- count_df$count[i]
    }
    
    ## Output Cox HR for binary classifications
    if (length(group_counts) == 2) {
      hr_df <- extract_hr_from_cox(
        sub_data = sub_data,
        time_col = ".time",
        event_col = ".event",
        group_col = ".group",
        ref_level = sort(unique(as.character(sub_data$.group)))[1]
      )
      res_row <- dplyr::bind_cols(res_row, hr_df)
    } else {
      res_row <- dplyr::bind_cols(
        res_row,
        data.frame(HR = NA_real_, lower95 = NA_real_, upper95 = NA_real_, cox_p = NA_real_)
      )
    }
    
    results_list[[ct]] <- res_row
    
    do_plot <- if (plot_only_sig) {
      !is.na(p_value) && p_value < sig_cutoff
    } else {
      !is.na(p_value)
    }
    
    if (do_plot) {
      present_levels <- levels(droplevels(sub_data$.group))
      pal_use <- palette_vec[present_levels]
      pal_use <- pal_use[!is.na(pal_use)]
      
      p <- tryCatch(
        survminer::ggsurvplot(
          fit = fit,
          data = sub_data,
          pval = TRUE,
          conf.int = FALSE,
          risk.table = TRUE,
          risk.table.height = 0.25,
          palette = unname(pal_use),
          legend.title = legend_title,
          legend.labs = present_levels,
          xlab = "Overall survival time (days)",
          ylab = "Overall survival probability",
          title = paste0(ct, ": OS by ", analysis_name),
          ggtheme = theme_classic(base_size = 11)
        ),
        error = function(e) {
          msg("ggsurvplot failed for ", ct, ": ", e$message)
          return(NULL)
        }
      )
      
      if (!is.null(p)) {
        plots_list[[ct]] <- p
      }
    }
    
    msg("Completed ", analysis_name, " | ", ct,
        " | p = ", signif(p_value, 3))
  }
  
  if (length(results_list) == 0) {
    msg("No cancer type passed filters for ", analysis_name)
    return(list(
      results = data.frame(),
      plots = list()
    ))
  }
  
  results_summary <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      significant_bh = ifelse(p_adj < 0.05, "Yes", "No"),
      log10_p = -log10(p_value),
      log10_p_adj = -log10(p_adj)
    ) %>%
    dplyr::arrange(p_value)
  
  save_csv(results_summary, paste0(out_prefix, "_survival_results_by_cancer.csv"))
  saveRDS(results_summary, file.path(rds_dir, paste0(out_prefix, "_survival_results_by_cancer.rds")))
  
  sig_results <- results_summary %>%
    dplyr::filter(p_value < 0.05)
  save_csv(sig_results, paste0(out_prefix, "_survival_results_significant_p_lt_0.05.csv"))
  
  km_dir <- file.path(plot_dir, paste0("km_plots_", out_prefix, "_significant"))
  dir.create(km_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(plots_list) > 0) {
    for (ct in names(plots_list)) {
      p <- plots_list[[ct]]
      pdf(file.path(km_dir, paste0(ct, "_OS_", out_prefix, ".pdf")), width = 8.5, height = 7)
      print(p)
      dev.off()
    }
  }
  
  bubble_data <- results_summary %>%
    dplyr::mutate(
      cancer_abbr = factor(cancer_abbr, levels = cancer_abbr[order(p_value)])
    )
  
  p_bubble <- ggplot(bubble_data, aes(x = cancer_abbr, y = log10_p)) +
    geom_point(aes(size = n, fill = significant_bh), shape = 21, color = "black", alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("Yes" = "#DB2462", "No" = "#47957F")) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = paste0("Cancer-specific OS association: ", analysis_name),
      subtitle = "Bubble size indicates sample size; fill indicates BH significance",
      x = "Cancer type",
      y = "-log10(raw P value)",
      fill = "BH < 0.05",
      size = "Sample size"
    ) +
    theme_pub(base_size = 12)
  
  save_pdf(p_bubble, paste0(out_prefix, "_survival_bubbleplot_by_cancer_BH.pdf"), width = 12.5, height = 5.4)
  save_png(p_bubble, paste0(out_prefix, "_survival_bubbleplot_by_cancer_BH.png"), width = 12.5, height = 5.4)
  
  return(list(
    results = results_summary,
    plots = plots_list
  ))
}

###############################################################################
## 2. Load data
###############################################################################
msg("Step 1/9: Loading Step4 result ...")

annot4_main <- readRDS(annot4_file)

msg("Loaded annot4_main: ", nrow(annot4_main), " rows x ", ncol(annot4_main), " columns")

###############################################################################
## 3. Prepare untreated tumor dataset and TP53 status
###############################################################################
msg("Step 2/9: Preparing tumor dataset ...")

df_tumor <- annot4_main %>%
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
  )

msg("Tumor sample count: ", nrow(df_tumor))
saveRDS(df_tumor, file.path(rds_dir, "step6_tumor_dataset.rds"))

###############################################################################
## 4. Nucleolar functional stratification
###############################################################################
msg("Step 3/9: Defining nucleolar functional groups ...")

## 4.1 Global tumor-level grouping
df_tumor <- df_tumor %>%
  dplyr::mutate(
    Ribo_group_global = ifelse(RiboSis_z > median(RiboSis_z, na.rm = TRUE), "H", "L"),
    NuS_group_global  = ifelse(NuStress_z > median(NuStress_z, na.rm = TRUE), "H", "L"),
    NuFunc_group_global = paste(Ribo_group_global, NuS_group_global, sep = "-")
  )

## 4.2 Cancer-specific grouping for survival analysis
df_tumor <- df_tumor %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::mutate(
    Ribo_group_cancer = ifelse(RiboSis_z > median(RiboSis_z, na.rm = TRUE), "High", "Low"),
    NuS_group_cancer  = ifelse(NuStress_z > median(NuStress_z, na.rm = TRUE), "High", "Low"),
    NuFunc_group_cancer = paste(
      ifelse(Ribo_group_cancer == "High", "H", "L"),
      ifelse(NuS_group_cancer == "High", "H", "L"),
      sep = "-"
    )
  ) %>%
  dplyr::ungroup()

saveRDS(df_tumor, file.path(rds_dir, "step6_tumor_dataset_with_nufunc.rds"))

## 4.3 Tumor-only correlation between NuStress and RiboSis
corr_tumor <- safe_spearman(df_tumor, "NuStress_z", "RiboSis_z")
save_csv(corr_tumor, "step6_tumor_correlation_NuStress_vs_RiboSis.csv")

p_corr_tumor <- plot_corr_scatter(
  df = df_tumor,
  x = "NuStress_z",
  y = "RiboSis_z",
  color = "#4c72b0",
  title = "Tumor-wide relationship between NuStress and RiboSis",
  xlab = "NuStress_z",
  ylab = "RiboSis_z"
)

save_pdf(p_corr_tumor, "step6_tumor_correlation_NuStress_vs_RiboSis.pdf", width = 6.4, height = 5.8)
save_png(p_corr_tumor, "step6_tumor_correlation_NuStress_vs_RiboSis.png", width = 6.4, height = 5.8)

## 4.4 Cancer mean-level functional map
cancer_means <- df_tumor %>%
  dplyr::filter(!is.na(cancer_abbr), cancer_abbr != "") %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(
    mean_RiboSis = mean(RiboSis_z, na.rm = TRUE),
    mean_NuStress = mean(NuStress_z, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

median_ribo_mean <- median(cancer_means$mean_RiboSis, na.rm = TRUE)
median_nus_mean  <- median(cancer_means$mean_NuStress, na.rm = TRUE)

cancer_means <- cancer_means %>%
  dplyr::mutate(
    Ribo_group = ifelse(mean_RiboSis > median_ribo_mean, "H", "L"),
    NuS_group  = ifelse(mean_NuStress > median_nus_mean, "H", "L"),
    NuFunc_group = paste(Ribo_group, NuS_group, sep = "-")
  )

save_csv(cancer_means, "step6_cancer_mean_nucleolar_function.csv")

p_cancer_map <- ggplot(cancer_means, aes(x = mean_RiboSis, y = mean_NuStress, color = NuFunc_group)) +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_hline(yintercept = median_nus_mean, linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = median_ribo_mean, linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(label = cancer_abbr), vjust = -0.7, size = 3.2) +
  scale_color_manual(values = nufunc_pal) +
  labs(
    title = "Cancer-level nucleolar functional stratification",
    subtitle = "Based on mean NuStress_z and mean RiboSis_z in untreated tumors",
    x = "Mean RiboSis_z",
    y = "Mean NuStress_z",
    color = "NuFunc group",
    size = "Sample size"
  ) +
  theme_pub(base_size = 12)

save_pdf(p_cancer_map, "step6_cancer_mean_nucleolar_function_map.pdf", width = 6.6, height = 5.8)
save_png(p_cancer_map, "step6_cancer_mean_nucleolar_function_map.png", width = 6.6, height = 5.8)

## 4.5 Composition of nucleolar functional groups by cancer
df_summary <- df_tumor %>%
  dplyr::filter(!is.na(cancer_abbr), cancer_abbr != "") %>%
  dplyr::group_by(cancer_abbr, NuFunc_group_global) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup()

save_csv(df_summary, "step6_nufunc_group_proportion_by_cancer.csv")

p_stack <- ggplot(df_summary, aes(x = cancer_abbr, y = freq, fill = NuFunc_group_global)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = nufunc_pal) +
  labs(
    title = "Distribution of nucleolar functional groups across cancers",
    x = "Cancer type",
    y = "Proportion",
    fill = "NuFunc group"
  ) +
  theme_pub(base_size = 12)

save_pdf(p_stack, "step6_nufunc_group_proportion_by_cancer.pdf", width = 12.5, height = 5.2)
save_png(p_stack, "step6_nufunc_group_proportion_by_cancer.png", width = 12.5, height = 5.2)

###############################################################################
## 5. Survival analysis input
###############################################################################
msg("Step 4/9: Preparing survival input ...")

surv_df <- df_tumor %>%
  dplyr::filter(!is.na(os), !is.na(os_time), !is.na(cancer_abbr), cancer_abbr != "") %>%
  dplyr::mutate(
    cancer_abbr = trimws(as.character(cancer_abbr)),
    NuFunc_group_cancer = as.character(NuFunc_group_cancer),
    Ribo_group_cancer   = as.character(Ribo_group_cancer),
    NuS_group_cancer    = as.character(NuS_group_cancer)
  )

saveRDS(surv_df, file.path(rds_dir, "step6_survival_input_dataset.rds"))

###############################################################################
## 6. Joint survival: RiboSis + NuS
###############################################################################
msg("Step 5/9: Joint survival analysis (RiboSis + NuS) ...")

joint_res <- run_survival_by_group(
  data = surv_df,
  group_col = "NuFunc_group_cancer",
  analysis_name = "joint nucleolar functional group",
  legend_title = "NuFunc group",
  palette_vec = nufunc_pal,
  out_prefix = "step6_joint_nufunc",
  min_n_total = 20,
  min_n_each = 3,
  plot_only_sig = TRUE,
  sig_cutoff = 0.05
)

###############################################################################
## 7. Single survival: RiboSis alone
###############################################################################
msg("Step 6/9: Single-marker survival analysis for RiboSis ...")

ribo_res <- run_survival_by_group(
  data = surv_df,
  group_col = "Ribo_group_cancer",
  analysis_name = "RiboSis group",
  legend_title = "RiboSis",
  palette_vec = single_pal,
  out_prefix = "step6_single_RiboSis",
  min_n_total = 20,
  min_n_each = 3,
  plot_only_sig = TRUE,
  sig_cutoff = 0.05
)

###############################################################################
## 8. Single survival: NuS alone
###############################################################################
msg("Step 7/9: Single-marker survival analysis for NuS ...")

nus_res <- run_survival_by_group(
  data = surv_df,
  group_col = "NuS_group_cancer",
  analysis_name = "NuS group",
  legend_title = "NuS",
  palette_vec = single_pal,
  out_prefix = "step6_single_NuS",
  min_n_total = 20,
  min_n_each = 3,
  plot_only_sig = TRUE,
  sig_cutoff = 0.05
)

###############################################################################
## 9. Merge all results and export summary
###############################################################################
msg("Step 8/9: Merging survival results ...")

joint_tbl <- joint_res$results %>% dplyr::mutate(model = "Joint")
ribo_tbl  <- ribo_res$results  %>% dplyr::mutate(model = "RiboSis")
nus_tbl   <- nus_res$results   %>% dplyr::mutate(model = "NuS")

all_surv_results <- bind_rows(joint_tbl, ribo_tbl, nus_tbl) %>%
  dplyr::arrange(model, p_value)

save_csv(all_surv_results, "step6_all_survival_results_combined.csv")
saveRDS(all_surv_results, file.path(rds_dir, "step6_all_survival_results_combined.rds"))

## Keep only significant result summaries
all_sig_results <- all_surv_results %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::mutate(
    p_value_label = format_p(p_value),
    p_adj_label   = format_p(p_adj),
    HR_label = ifelse(
      is.na(HR),
      NA_character_,
      paste0(sprintf("%.2f", HR), " (",
             sprintf("%.2f", lower95), "-",
             sprintf("%.2f", upper95), ")")
    )
  ) %>%
  dplyr::select(
    model, analysis, cancer_abbr, n, n_groups,
    p_value, p_value_label, p_adj, p_adj_label,
    significant, significant_bh,
    HR, lower95, upper95, HR_label, cox_p,
    dplyr::everything()
  )

save_csv(all_sig_results, "step6_all_survival_results_significant_only.csv")

###############################################################################
## Summary heatmap-like dotplot
###############################################################################
plot_sig_df <- all_surv_results %>%
  dplyr::mutate(
    model = factor(model, levels = c("Joint", "RiboSis", "NuS")),
    cancer_abbr = factor(cancer_abbr, levels = unique(cancer_abbr[order(p_value)])),
    log10_p = -log10(p_value)
  )

p_compare <- ggplot(plot_sig_df, aes(x = model, y = cancer_abbr)) +
  geom_point(aes(size = n, fill = log10_p), shape = 21, color = "black", alpha = 0.9) +
  scale_size_continuous(range = c(2.5, 8)) +
  scale_fill_gradient(low = "#d9f0d3", high = "#7a0177") +
  labs(
    title = "Comparison of survival associations across three models",
    subtitle = "Joint = RiboSis + NuS; single-marker models are cancer-specific median splits",
    x = "Model",
    y = "Cancer type",
    size = "Sample size",
    fill = "-log10(P)"
  ) +
  theme_pub(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

save_pdf(p_compare, "step6_survival_model_comparison_dotplot.pdf", width = 8.2, height = 10.5)
save_png(p_compare, "step6_survival_model_comparison_dotplot.png", width = 8.2, height = 10.5)

###############################################################################
## 10. Final message and session info
###############################################################################
msg("Step 9/9: Finished.")

msg("Joint significant cancers (raw p < 0.05):")
print(joint_tbl %>% dplyr::filter(p_value < 0.05))

msg("RiboSis-alone significant cancers (raw p < 0.05):")
print(ribo_tbl %>% dplyr::filter(p_value < 0.05))

msg("NuS-alone significant cancers (raw p < 0.05):")
print(nus_tbl %>% dplyr::filter(p_value < 0.05))

sessionInfo()


###############################################################################
## Summary comparison plot (horizontal) 
## significant cancers highlighted by color
###############################################################################
plot_sig_df <- all_surv_results %>%
  dplyr::mutate(
    model = factor(model, levels = c("Joint", "RiboSis", "NuS")),
    sig_flag = ifelse(p_value < 0.05, "P < 0.05", "NS"),
    sig_bh_flag = ifelse(p_adj < 0.05, "BH < 0.05", "BH >= 0.05")
  )

## Sort cancer types by minimum P value
cancer_order <- plot_sig_df %>%
  dplyr::group_by(cancer_abbr) %>%
  dplyr::summarise(min_p = min(p_value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(min_p) %>%
  dplyr::pull(cancer_abbr)

plot_sig_df <- plot_sig_df %>%
  dplyr::mutate(
    cancer_abbr = factor(cancer_abbr, levels = cancer_order),
    sig_flag = factor(sig_flag, levels = c("NS", "P < 0.05")),
    sig_bh_flag = factor(sig_bh_flag, levels = c("BH >= 0.05", "BH < 0.05"))
  )

## Version A: mark significance by raw P < 0.05
p_compare_horizontal_raw <- ggplot(plot_sig_df, aes(x = cancer_abbr, y = model)) +
  geom_point(aes(size = n, fill = sig_flag), shape = 21, color = "black", alpha = 0.95) +
  geom_text(
    data = subset(plot_sig_df, p_value < 0.05),
    aes(label = "*"),
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("NS" = "grey85", "P < 0.05" = "#DB2462")) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "Comparison of survival associations across three models",
    subtitle = "Red indicates statistically significant association (raw P < 0.05); * marks significant results",
    x = "Cancer type",
    y = "Model",
    fill = "Significance",
    size = "Sample size"
  ) +
  theme_pub(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold")
  )

save_pdf(p_compare_horizontal_raw, "step6_survival_model_comparison_horizontal_rawP.pdf",
         width = 13, height = 4.8)
save_png(p_compare_horizontal_raw, "step6_survival_model_comparison_horizontal_rawP.png",
         width = 13, height = 4.8)

## Version B: mark significance by BH < 0.05 (more stringent; recommended for the main text)
p_compare_horizontal_bh <- ggplot(plot_sig_df, aes(x = cancer_abbr, y = model)) +
  geom_point(aes(size = n, fill = sig_bh_flag), shape = 21, color = "black", alpha = 0.95) +
  geom_text(
    data = subset(plot_sig_df, p_adj < 0.05),
    aes(label = "*"),
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("BH >= 0.05" = "grey85", "BH < 0.05" = "#6A1B9A")) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "Comparison of survival associations across three models",
    subtitle = "Purple indicates statistically significant association after BH correction; * marks BH-significant results",
    x = "Cancer type",
    y = "Model",
    fill = "BH significance",
    size = "Sample size"
  ) +
  theme_pub(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold")
  )

save_pdf(p_compare_horizontal_bh, "step6_survival_model_comparison_horizontal_BH.pdf",
         width = 13, height = 4.8)
save_png(p_compare_horizontal_bh, "step6_survival_model_comparison_horizontal_BH.png",
         width = 13, height = 4.8)

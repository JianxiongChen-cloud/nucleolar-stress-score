rm(list = ls())
options(stringsAsFactors = FALSE)

#===============================
# 0. 初始化
#===============================
setwd("/home/xxm_xxm/CJX_workspace/geneset")

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

gene_set_dir <- file.path(getwd(), "final_nucleolar_gene_sets")
depmap_file <- file.path(getwd(), "CRISPRGeneEffect.csv")

step6_dir <- file.path(getwd(), "step6_DepMap_dependency")
output_dir <- file.path(step6_dir, "output")

if (!dir.exists(step6_dir)) dir.create(step6_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Gene set directory:", gene_set_dir, "\n")
cat("DepMap file:", depmap_file, "\n")
cat("Output directory:", output_dir, "\n\n")

if (!file.exists(depmap_file)) {
  stop("Cannot find DepMap file: ", depmap_file)
}

#===============================
# 1. 读取 NuStress_UP / NuStressCore_UP
#===============================
read_gene_symbols <- function(file_path) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  x <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("gene_symbol" %in% colnames(x))) {
    stop("Missing 'gene_symbol' column in: ", file_path)
  }
  
  g <- unique(trimws(toupper(x$gene_symbol)))
  g <- g[!is.na(g) & g != ""]
  sort(g)
}

# 若你的文件名不同，在这里改
NuStress_UP <- read_gene_symbols(file.path(gene_set_dir, "NuStress_UP.csv"))
NuStressCore_UP <- read_gene_symbols(file.path(gene_set_dir, "NuStressCore_UP.csv"))

cat("NuStress_UP genes:", length(NuStress_UP), "\n")
cat("NuStressCore_UP genes:", length(NuStressCore_UP), "\n\n")

#===============================
# 2. 读取 DepMap 矩阵
#===============================
depmap_dt <- fread(depmap_file, data.table = FALSE)

cat("DepMap matrix dimension:", nrow(depmap_dt), "x", ncol(depmap_dt), "\n")
cat("First columns:\n")
print(colnames(depmap_dt)[1:min(8, ncol(depmap_dt))])
cat("\n")

#===============================
# 3. 清洗 DepMap 列名
# 你的格式是: A1BG (1), A1CF (29974), ...
#===============================
extract_gene_symbol <- function(x) {
  x <- gsub("\\s*\\([^\\)]*\\)$", "", x)   # 去掉末尾括号ID
  x <- trimws(toupper(x))
  x
}

all_colnames <- colnames(depmap_dt)
id_col <- all_colnames[1]

gene_cols_raw <- all_colnames[-1]
gene_cols_clean <- extract_gene_symbol(gene_cols_raw)

gene_map <- data.frame(
  raw_col = gene_cols_raw,
  gene_symbol = gene_cols_clean,
  stringsAsFactors = FALSE
)

# 如有重复 gene symbol，只保留第一次出现
gene_map_unique <- gene_map[!duplicated(gene_map$gene_symbol), ]

cat("Unique DepMap genes available:", nrow(gene_map_unique), "\n\n")

#===============================
# 4. 提取目标基因 dependency matrix
#===============================
get_depmap_subset <- function(gene_set, gene_set_name, depmap_df, gene_map_unique, id_col) {
  genes_found <- intersect(gene_set, gene_map_unique$gene_symbol)
  genes_missing <- setdiff(gene_set, gene_map_unique$gene_symbol)
  
  gene_map_sub <- gene_map_unique[match(genes_found, gene_map_unique$gene_symbol), , drop = FALSE]
  gene_map_sub <- gene_map_sub[!is.na(gene_map_sub$raw_col), , drop = FALSE]
  
  if (nrow(gene_map_sub) == 0) {
    sub_df <- depmap_df[, id_col, drop = FALSE]
  } else {
    sub_df <- depmap_df[, c(id_col, gene_map_sub$raw_col), drop = FALSE]
    colnames(sub_df) <- c(id_col, gene_map_sub$gene_symbol)
  }
  
  list(
    data = sub_df,
    genes_found = sort(genes_found),
    genes_missing = sort(genes_missing),
    gene_map = gene_map_sub,
    gene_set_name = gene_set_name
  )
}

dep_full <- get_depmap_subset(NuStress_UP, "Full Set", depmap_dt, gene_map_unique, id_col)
dep_core <- get_depmap_subset(NuStressCore_UP, "Core Set", depmap_dt, gene_map_unique, id_col)

cat("Full Set found in DepMap:", length(dep_full$genes_found), "\n")
cat("Full Set missing in DepMap:", length(dep_full$genes_missing), "\n")
cat("Core Set found in DepMap:", length(dep_core$genes_found), "\n")
cat("Core Set missing in DepMap:", length(dep_core$genes_missing), "\n\n")

write.csv(
  data.frame(gene_symbol = dep_full$genes_missing, stringsAsFactors = FALSE),
  file.path(output_dir, "FullSet_missing_in_DepMap.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = dep_core$genes_missing, stringsAsFactors = FALSE),
  file.path(output_dir, "CoreSet_missing_in_DepMap.csv"),
  row.names = FALSE
)

#===============================
# 5. 每个基因聚合 dependency score
# 默认用 median；若想改成 mean，把 aggregate_fun 改为 "mean"
#===============================
aggregate_fun <- "median"   # "median" or "mean"
essential_cutoff <- -1

summarize_dependency <- function(dep_obj, aggregate_fun = "median", essential_cutoff = -1) {
  df <- dep_obj$data
  gene_names <- colnames(df)[-1]
  
  if (length(gene_names) == 0) {
    stop(
      paste0(
        "No matched genes found in DepMap for gene set: ",
        dep_obj$gene_set_name,
        ". Please check gene symbols."
      )
    )
  }
  
  gene_scores_list <- lapply(gene_names, function(g) {
    x <- suppressWarnings(as.numeric(df[[g]]))
    
    agg_score <- if (aggregate_fun == "median") {
      median(x, na.rm = TRUE)
    } else if (aggregate_fun == "mean") {
      mean(x, na.rm = TRUE)
    } else {
      stop("aggregate_fun must be 'median' or 'mean'")
    }
    
    data.frame(
      gene_symbol = g,
      n_cell_lines = sum(!is.na(x)),
      aggregated_dependency = agg_score,
      is_essential = agg_score < essential_cutoff,
      stringsAsFactors = FALSE
    )
  })
  
  gene_scores <- do.call(rbind, gene_scores_list)
  gene_scores <- as.data.frame(gene_scores, stringsAsFactors = FALSE)
  
  gene_scores$n_cell_lines <- as.numeric(gene_scores$n_cell_lines)
  gene_scores$aggregated_dependency <- as.numeric(gene_scores$aggregated_dependency)
  gene_scores$is_essential <- as.logical(gene_scores$is_essential)
  gene_scores$gene_set <- dep_obj$gene_set_name
  
  gene_scores <- gene_scores[order(gene_scores$aggregated_dependency), , drop = FALSE]
  
  essential_n <- sum(gene_scores$is_essential, na.rm = TRUE)
  total_n <- nrow(gene_scores)
  essential_percent <- essential_n / total_n * 100
  
  summary_row <- data.frame(
    gene_set = dep_obj$gene_set_name,
    aggregate_fun = aggregate_fun,
    total_genes_in_input = length(unique(c(dep_obj$genes_found, dep_obj$genes_missing))),
    genes_found_in_DepMap = length(dep_obj$genes_found),
    genes_missing_in_DepMap = length(dep_obj$genes_missing),
    analyzed_genes = total_n,
    essential_cutoff = essential_cutoff,
    essential_genes_n = essential_n,
    essential_genes_percent = essential_percent,
    stringsAsFactors = FALSE
  )
  
  list(
    gene_level = gene_scores,
    summary = summary_row
  )
}

res_full <- summarize_dependency(dep_full, aggregate_fun, essential_cutoff)
res_core <- summarize_dependency(dep_core, aggregate_fun, essential_cutoff)

gene_level_df <- rbind(res_full$gene_level, res_core$gene_level)
summary_df <- rbind(res_full$summary, res_core$summary)

write.csv(
  gene_level_df,
  file.path(output_dir, "DepMap_dependency_gene_level_summary.csv"),
  row.names = FALSE
)

write.csv(
  summary_df,
  file.path(output_dir, "DepMap_dependency_set_level_summary.csv"),
  row.names = FALSE
)

write.csv(
  subset(gene_level_df, gene_set == "Full Set" & is_essential),
  file.path(output_dir, "FullSet_broadly_essential_genes.csv"),
  row.names = FALSE
)

write.csv(
  subset(gene_level_df, gene_set == "Core Set" & is_essential),
  file.path(output_dir, "CoreSet_broadly_essential_genes.csv"),
  row.names = FALSE
)

#===============================
# 6. 绘图函数：尽量贴近原图
#===============================
make_dep_density_plot <- function(plot_df_one, essential_percent, set_title,
                                  essential_cutoff = -1,
                                  xlim_use = c(-3, 0.5),
                                  hist_fill = "#C8D6A3",
                                  density_line_col = "#D97B93",
                                  density_fill_col = "#F4C7D3",
                                  line_col = "#5A84B3") {
  
  dens_obj <- density(plot_df_one$aggregated_dependency, na.rm = TRUE, adjust = 1)
  dens_df <- data.frame(
    x = dens_obj$x,
    y = dens_obj$y
  )
  
  # 同时考虑 histogram 和 density 的最大高度，避免柱子或曲线被截断
  hist_obj <- hist(plot_df_one$aggregated_dependency,
                   breaks = 28,
                   plot = FALSE)
  hist_density_max <- ifelse(
    sum(hist_obj$counts) > 0,
    max(hist_obj$density, na.rm = TRUE),
    0
  )
  
  dens_max <- max(dens_df$y, na.rm = TRUE)
  ymax_raw <- max(dens_max, hist_density_max, na.rm = TRUE)
  
  # 顶部留白更充分，确保完整显示
  ymax <- ymax_raw * 1.30
  
  arrow_y <- ymax * 0.72
  text_y <- ymax * 0.84
  
  ggplot(plot_df_one, aes(x = aggregated_dependency)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 28,
      fill = hist_fill,
      color = hist_fill,
      alpha = 0.9
    ) +
    geom_area(
      data = dens_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      fill = density_fill_col,
      alpha = 0.45
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = density_line_col,
      linewidth = 0.9
    ) +
    geom_vline(
      xintercept = essential_cutoff,
      linetype = "dashed",
      color = line_col,
      linewidth = 0.8
    ) +
    annotate("segment",
             x = xlim_use[1] + 0.03,
             xend = essential_cutoff - 0.03,
             y = arrow_y,
             yend = arrow_y,
             colour = line_col,
             linewidth = 0.65,
             arrow = arrow(length = unit(0.18, "cm"), ends = "both")) +
    annotate("text",
             x = (xlim_use[1] + essential_cutoff) / 2,
             y = text_y,
             label = paste0(sprintf("%.0f", essential_percent), "%"),
             colour = line_col,
             size = 4.1,
             fontface = "bold") +
    labs(
      title = set_title,
      subtitle = "Distribution of Gene Dependency Scores",
      x = "Dependency Score",
      y = "Density"
    ) +
    coord_cartesian(
      xlim = xlim_use,
      ylim = c(0, ymax),
      expand = TRUE,
      clip = "off"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_classic(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.6),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 10.5),
      axis.text = element_text(size = 10),
      plot.margin = margin(8, 10, 8, 8)
    )
}
x_min <- min(gene_level_df$aggregated_dependency, na.rm = TRUE)
x_max <- max(gene_level_df$aggregated_dependency, na.rm = TRUE)
xlim_use <- c(min(-3, floor(x_min * 2) / 2), max(0.5, ceiling(x_max * 2) / 2))

full_percent <- summary_df$essential_genes_percent[summary_df$gene_set == "Full Set"]
core_percent <- summary_df$essential_genes_percent[summary_df$gene_set == "Core Set"]

p_full <- make_dep_density_plot(
  plot_df_one = subset(gene_level_df, gene_set == "Full Set"),
  essential_percent = full_percent,
  set_title = "Full Set",
  essential_cutoff = essential_cutoff,
  xlim_use = xlim_use
)

p_core <- make_dep_density_plot(
  plot_df_one = subset(gene_level_df, gene_set == "Core Set"),
  essential_percent = core_percent,
  set_title = "Core Set",
  essential_cutoff = essential_cutoff,
  xlim_use = xlim_use
)

pdf(file.path(output_dir, "Figure_step6_DepMap_dependency_density_full.pdf"),
    width = 8, height = 3)
print(p_full)
dev.off()

pdf(file.path(output_dir, "Figure_step6_DepMap_dependency_density_core.pdf"),
    width = 8, height = 3)
print(p_core)
dev.off()

#===============================
# 7. 辅助图：essential gene proportion
#===============================
my_cols <- c("#C4E7C1", "#93D4BC", "#51B3D1", "#08589E")

p_bar <- ggplot(summary_df, aes(x = gene_set, y = essential_genes_percent, fill = gene_set)) +
  geom_col(width = 0.58, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", essential_genes_percent)),
            vjust = -0.35, size = 3.8) +
  scale_fill_manual(values = c(my_cols[3], my_cols[4])) +
  labs(
    title = "Broadly essential genes in nucleolar stress gene sets",
    subtitle = paste0("Threshold: aggregated dependency score < ", essential_cutoff),
    x = NULL,
    y = "Essential genes (%)",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  coord_cartesian(ylim = c(0, max(summary_df$essential_genes_percent, na.rm = TRUE) * 1.18))

pdf(file.path(output_dir, "Figure_step6_DepMap_dependency_essential_percent.pdf"),
    width = 5.0, height = 4.5)
print(p_bar)
dev.off()

#===============================
# 8. 保存对象
#===============================
save(
  NuStress_UP, NuStressCore_UP,
  dep_full, dep_core,
  res_full, res_core,
  gene_level_df, summary_df,
  aggregate_fun, essential_cutoff,
  file = file.path(output_dir, "step6_DepMap_dependency_objects.RData")
)

#===============================
# 9. 输出
#===============================
cat("Step 6 finished.\n\n")

cat("Main tables:\n")
cat(file.path(output_dir, "DepMap_dependency_gene_level_summary.csv"), "\n")
cat(file.path(output_dir, "DepMap_dependency_set_level_summary.csv"), "\n\n")

cat("Figures:\n")
cat(file.path(output_dir, "Figure_step6_DepMap_dependency_density_full.pdf"), "\n")
cat(file.path(output_dir, "Figure_step6_DepMap_dependency_density_core.pdf"), "\n")
cat(file.path(output_dir, "Figure_step6_DepMap_dependency_essential_percent.pdf"), "\n\n")

cat("Saved RData:\n")
cat(file.path(output_dir, "step6_DepMap_dependency_objects.RData"), "\n\n")

print(summary_df)